
// testDlg.cpp : implementation file
//

#include "stdafx.h"
#include "test.h"
#include "testDlg.h"
#include "afxdialogex.h"
#include "plot.h"

#include <iostream>
#include <array>
#include <numeric>
#include <atomic>
#include <thread>
#include <ctime>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#define WM_UPDATEDATA_MESSAGE (WM_USER + 1000)

namespace
{

    phd::simple_list_plot frequency_bias_last_plot,
                          frequency_bias_avg_plot;

    std::atomic<bool> working;
}

CtestDlg::CtestDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CtestDlg::IDD, pParent)
    , signal_magnitude_from(1)
    , signal_magnitude_to(3)
    , signal_frequency_from(0.02)
    , signal_frequency_to(0.28)
    , snr_percents_from(0)
    , snr_percents_to(100)
    , snr_percents_step_count(100)
    , acm_order(3)
    , number_of_experiments(1000)
    , sample_count(128)
    , sampling_period(1)
    , progress(0)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
    
    frequency_bias_last_plot_ctrl.plot_layer.with(
        (plot::plot_builder() << frequency_bias_last_plot)
        .with_ticks(plot::palette::pen(RGB(150, 150, 0)))
        .with_x_ticks(0, 10, 0)
        .with_y_ticks(0, 5, 5)
        .build()
    );
    frequency_bias_avg_plot_ctrl.plot_layer.with(
        (plot::plot_builder() << frequency_bias_avg_plot)
        .with_ticks(plot::palette::pen(RGB(150, 150, 0)))
        .with_x_ticks(0, 10, 0)
        .with_y_ticks(0, 5, 5)
        .build()
    );
}

void CtestDlg::DoDataExchange(CDataExchange* pDX)
{
    CDialog::DoDataExchange(pDX);
    DDX_Text(pDX, IDC_SIGNAL_MAGNITUDE_FROM, signal_magnitude_from);
    DDX_Text(pDX, IDC_SIGNAL_MAGNITUDE_TO, signal_magnitude_to);
    DDX_Text(pDX, IDC_SIGNAL_FREQUENCY_FROM, signal_frequency_from);
    DDX_Text(pDX, IDC_SIGNAL_FREQUENCY_TO, signal_frequency_to);
    DDX_Text(pDX, IDC_SNR_FROM, snr_percents_from);
    DDX_Text(pDX, IDC_SNR_TO, snr_percents_to);
    DDX_Text(pDX, IDC_SNR_STEPS, snr_percents_step_count);
    DDX_Text(pDX, IDC_AUTOCORRELATION_MATRIX_ORDER_2, acm_order);
    DDX_Text(pDX, IDC_NUMBER_OF_EXPERIMENTS, number_of_experiments);
    DDX_Control(pDX, IDC_LAST_BIASES_PLOT, frequency_bias_last_plot_ctrl);
    DDX_Control(pDX, IDC_AVERAGE_BIASES_PLOT, frequency_bias_avg_plot_ctrl);
    DDX_Text(pDX, IDC_SAMPLE_COUNT_2, sample_count);
    DDX_Text(pDX, IDC_SAMPLING_PERIOD_2, sampling_period);
    DDX_Text(pDX, IDC_PROGRESS, progress);
}

BEGIN_MESSAGE_MAP(CtestDlg, CDialog)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
    ON_BN_CLICKED(IDC_MANUAL, &CtestDlg::OnBnClickedManual)
    ON_BN_CLICKED(IDC_START_PAUSE, &CtestDlg::OnBnClickedStartPause)
    ON_BN_CLICKED(IDC_STOP, &CtestDlg::OnBnClickedStop)
    ON_MESSAGE(WM_UPDATEDATA_MESSAGE, &CtestDlg::OnUpdateDataMessage)
END_MESSAGE_MAP()


// CtestDlg message handlers

BOOL CtestDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

    UpdateData(FALSE);
    UpdateData(TRUE);

    manual_view_dialog.Create(ManualViewDialog::IDD, this);

	// TODO: Add extra initialization here

	return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CtestDlg::OnPaint()
{
	CDialog::OnPaint();
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CtestDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CtestDlg::OnBnClickedManual()
{
    manual_view_dialog.ShowWindow(SW_SHOW);
    manual_view_dialog.UpdateWindow();

    return;
}


void CtestDlg::OnBnClickedStartPause()
{
    bool expected = false;
    if (!working.compare_exchange_strong(expected, true))
    {
        return;
    }
    std::thread worker(&CtestDlg::DoWork, this);
    worker.detach();
}


void CtestDlg::OnBnClickedStop()
{
    working = false;
}


void CtestDlg::DoWork()
{
    DoWork0();
}


void CtestDlg::DoWork0()
{
    RequestUpdateData(TRUE);
    
    srand(clock());

    // setup periods and bounds

    double snr_period = (snr_percents_to - snr_percents_from) / snr_percents_step_count;

    phd::sampled_t frequency_bias_last_sampled = phd::allocate_sampled(snr_percents_step_count + 1, snr_period);
    phd::sampled_t frequency_bias_avg_sampled  = phd::allocate_sampled(snr_percents_step_count + 1, snr_period);
    phd::sampled_t signal_sampled              = phd::allocate_sampled(sample_count,                sampling_period);
    phd::sampled_t noise_sampled               = phd::allocate_sampled(sample_count,                sampling_period);
    phd::sampled_t autocorrelation_sampled     = phd::allocate_sampled(acm_order,                   sampling_period);
    double *singular_u            = new double[acm_order * acm_order];
    double *singular_v            = new double[acm_order * acm_order];
    double *singular_sigma        = new double[acm_order];
    double *succeeded_experiments = new double[snr_percents_step_count + 1]();
    double roots_r[2];
    double roots_i[2];

    phd::map(frequency_bias_avg_sampled, [] (size_t, double d) { return 0; });

    for (size_t i = 0; i < number_of_experiments; i++)
    {
        double magnitude = phd::random(signal_magnitude_from, signal_magnitude_to);
        double frequency = phd::random(signal_frequency_from, signal_frequency_to);

        for (size_t j = 0; j < snr_percents_step_count + 1; j++)
        {
            phd::sin_signal_noised(magnitude, frequency,
                                   snr_period * j + snr_percents_from,
                                   signal_sampled, noise_sampled);
            phd::autocorrelation(signal_sampled, autocorrelation_sampled);
            int status = svd_toepliz(acm_order, acm_order,
                                     autocorrelation_sampled.samples,
                                     singular_u, singular_v, singular_sigma);
            // filter SVD fails
            if (status < 0)
            {
                frequency_bias_last_sampled.samples[j] = 0;
                continue;
            }

            assert(acm_order == 3);
            phd::solve_quad(singular_u[2 * 3 + 2],
                            singular_u[1 * 3 + 2],
                            singular_u[0 * 3 + 2],
                            roots_r, roots_i);

            double bias = abs(abs(acos(roots_r[0])) / (2 * M_PI * sampling_period) - frequency);

            // filter any fails of either PHD or SVD algorithms
            if ((bias > 1. / (2 * sampling_period))
                || (bias < -1. / (2 * sampling_period))
                || (bias == std::numeric_limits<double>::infinity())
                || (bias != bias))
            {
                frequency_bias_last_sampled.samples[j] = 0;
                continue;
            }

            ++succeeded_experiments[j];

            frequency_bias_last_sampled.samples[j] = bias;
            frequency_bias_avg_sampled.samples[j] += bias;

            if (!working)
            {
                return;
            }
        }

        phd::setup(frequency_bias_avg_plot, frequency_bias_avg_sampled, 0,
                   [succeeded_experiments] (size_t idx, double s) {
                       return s / succeeded_experiments[idx];
                   }
        );

        phd::setup(frequency_bias_last_plot, frequency_bias_last_sampled);

        frequency_bias_avg_plot_ctrl.RedrawWindow();
        frequency_bias_last_plot_ctrl.RedrawWindow();

        progress = i + 1;
        RequestUpdateData(FALSE);
    }

    phd::free_sampled(autocorrelation_sampled);
    phd::free_sampled(signal_sampled);
    phd::free_sampled(noise_sampled);
    delete[] singular_u;
    delete[] singular_v;
    delete[] singular_sigma;
    delete[] succeeded_experiments;

    working = false;
}


LRESULT CtestDlg::OnUpdateDataMessage(WPARAM wpD, LPARAM lpD)
{
    UpdateData(wpD == TRUE);
    return 0;
}

void CtestDlg::RequestUpdateData(BOOL saveAndValidate)
{
    SendMessage(WM_UPDATEDATA_MESSAGE, saveAndValidate);
}
