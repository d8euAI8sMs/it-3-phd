// ManualViewDialog.cpp : implementation file
//

#include "stdafx.h"
#include "test.h"
#include "ManualViewDialog.h"
#include "afxdialogex.h"


// ManualViewDialog dialog

IMPLEMENT_DYNAMIC(ManualViewDialog, CDialogEx)

namespace
{

    phd::simple_list_plot signal_plot,
                          autocorrelation_plot,
                          singular_plot,
                          estimated_frequency_plot  (0, { -1.3, 1.3, -1.3, 1.3 }, RGB(255, 000, 000), 5, false),
                          psd_plot, one_to_psd_plot;
}

ManualViewDialog::ManualViewDialog(CWnd* pParent /*=NULL*/)
	: CDialogEx(ManualViewDialog::IDD, pParent)
    , signal_magnitude(1)
    , signal_frequency(0.2)
    , signal_to_noise_ratio_percents(0)
    , acm_order(3)
    , sample_count(256)
    , sampling_period(1)
    , output_frequency(0)
    , output_bias(0)
    , method(FEM_PHD)
{
    signal_plot_ctrl.plot_layer.with(
        (plot::plot_builder() << signal_plot)
        .with_ticks(plot::palette::pen(RGB(150, 150, 0)))
        .with_x_ticks(0, 10, 0)
        .with_y_ticks(0, 5, 2)
        .build()
    );
    autocorrelation_plot_ctrl.plot_layer.with(
        (plot::plot_builder() << autocorrelation_plot)
        .with_ticks(plot::palette::pen(RGB(150, 150, 0)))
        .with_x_ticks(0, 10, 0)
        .with_y_ticks(0, 5, 2)
        .build()
    );
    singular_plot_ctrl.plot_layer.with(
        (plot::plot_builder() << singular_plot)
        .with_ticks(plot::palette::pen(RGB(150, 150, 0)))
        .with_x_ticks(0, 10, 1)
        .with_y_ticks(0, 5, 2)
        .build()
    );
    estimated_frequency_plot_ctrl.symmetric = true;
    estimated_frequency_plot_ctrl.plot_layer.with(
        (plot::plot_builder() << estimated_frequency_plot)
        .with_custom([this] (CDC &dc, const plot::viewport &bounds)
        {
            auto pen = plot::palette::pen(RGB(150, 150, 0));
            dc.SelectObject(*pen);
            // draw circle e^{i a}
            plot::point<int> left_top     = bounds.world_to_screen().xy({ -1,  1 });
            plot::point<int> bottom_right = bounds.world_to_screen().xy({  1, -1 });
            dc.Arc(left_top.x, left_top.y, bottom_right.x, bottom_right.y,
                   left_top.x, left_top.y, left_top.x,     left_top.y);
            // draw e^{i 2 pi +-f0 +- pi}
            double c = bounds.world.width()  * std::cos(2 * M_PI * signal_frequency * sampling_period);
            double s = bounds.world.height() * std::sin(2 * M_PI * signal_frequency * sampling_period);
            dc.MoveTo(bounds.world_to_screen().xy({  0,  0 }));
            dc.LineTo(bounds.world_to_screen().xy({  c,  s }));
            dc.MoveTo(bounds.world_to_screen().xy({  0,  0 }));
            dc.LineTo(bounds.world_to_screen().xy({  c, -s }));
        })
        .build()
    );
    psd_plot_ctrl.plot_layer.with(
        (plot::plot_builder() << psd_plot)
        .with_ticks(plot::palette::pen(RGB(150, 150, 0)))
        .with_x_ticks(0, 10, 1)
        .with_y_ticks(0, 5, 2)
        .build()
    );
    one_to_psd_plot_ctrl.plot_layer.with(
        (plot::plot_builder() << one_to_psd_plot)
        .with_ticks(plot::palette::pen(RGB(150, 150, 0)))
        .with_x_ticks(0, 10, 1)
        .with_y_ticks(0, 5, 2)
        .build()
    );
}

ManualViewDialog::~ManualViewDialog()
{
}

void ManualViewDialog::DoDataExchange(CDataExchange* pDX)
{
    CDialogEx::DoDataExchange(pDX);
    DDX_Control(pDX, IDC_SIGNAL_PLOT, signal_plot_ctrl);
    DDX_Control(pDX, IDC_AUTOCORRELATION_PLOT, autocorrelation_plot_ctrl);
    DDX_Control(pDX, IDC_SINGULAR_VALUE_PLOT, singular_plot_ctrl);
    DDX_Control(pDX, IDC_POLYNOMIAL_ROOTS_PLOT, estimated_frequency_plot_ctrl);
    DDX_Text(pDX, IDC_SIGNAL_MAGNITUDE, signal_magnitude);
    DDX_Text(pDX, IDC_SIGNAL_FREQUENCY, signal_frequency);
    DDX_Text(pDX, IDC_SIGNAL_NOISE_PERCENT, signal_to_noise_ratio_percents);
    DDX_Text(pDX, IDC_AUTOCORRELATION_ORDER, acm_order);
    DDX_Text(pDX, IDC_SAMPLE_COUNT, sample_count);
    DDX_Text(pDX, IDC_SAMPLING_PERIOD, sampling_period);
    DDX_Text(pDX, IDC_OUTPUT_FREQUENCY, output_frequency);
    DDX_Text(pDX, IDC_OUTPUT_BIAS, output_bias);
    DDX_Control(pDX, IDC_PSD_PLOT, psd_plot_ctrl);
    DDX_Control(pDX, IDC_ONE_TO_PSD_PLOT, one_to_psd_plot_ctrl);
    DDX_Radio(pDX, IDC_RADIO_MUSIC, method);
    DDX_Control(pDX, IDC_AUTOCORRELATION_ORDER, acm_order_edit_ctrl);
}


BEGIN_MESSAGE_MAP(ManualViewDialog, CDialogEx)
    ON_BN_CLICKED(IDC_APPLY, &ManualViewDialog::OnBnClickedApply)
    ON_CONTROL_RANGE(BN_CLICKED, IDC_RADIO_MUSIC, IDC_RADIO_MUSIC + 2, &ManualViewDialog::OnBnsClickedMethod)
END_MESSAGE_MAP()


// ManualViewDialog message handlers


void ManualViewDialog::OnBnClickedApply()
{
    UpdateData(TRUE);

    // hide the acm_order class field
    size_t acm_order = ((method == FEM_PHD) ? 3 : this->acm_order);

    size_t psd_sample_count    = 1024;
    size_t psd_sample_count_x2 = 1024 * 2 + 1;
    double psd_sampling_period = 1. / (sampling_period * psd_sample_count_x2);

    phd::sampled_t signal_sampled          = phd::allocate_sampled(sample_count,        sampling_period);
    phd::sampled_t noise_sampled           = phd::allocate_sampled(sample_count,        sampling_period);
    phd::sampled_t autocorrelation_sampled = phd::allocate_sampled(acm_order,           sampling_period);
    phd::sampled_t singular_sampled        = phd::allocate_sampled(acm_order,           1);
    phd::sampled_t psd_sampled             = phd::allocate_sampled(psd_sample_count_x2, psd_sampling_period);
    phd::sampled_t one_to_psd_sampled      = phd::allocate_sampled(psd_sample_count_x2, psd_sampling_period);

    double *singular_u     = new double[acm_order * acm_order];
    double *singular_v     = new double[acm_order * acm_order];
    double *psd_polynomial = new double[acm_order];
    double roots_r[2];
    double roots_i[2];

    // produce and setup sine signal

    phd::sin_signal_noised(signal_magnitude, signal_frequency,
                           signal_to_noise_ratio_percents,
                           signal_sampled, noise_sampled);
    phd::setup(signal_plot, signal_sampled, 0, phd::identity_un_op(), true, true);

    // produce and setup autocorrelation

    phd::autocorrelation(signal_sampled, autocorrelation_sampled);
    phd::setup(autocorrelation_plot, autocorrelation_sampled, 0, phd::identity_un_op(), true, true);

    // get eigenvalues and eigenvectors of the autocorrelation matrix, setup eigenvalues

    svd_toepliz(acm_order, acm_order,
                autocorrelation_sampled.samples,
                singular_u, singular_v, singular_sampled.samples);
    phd::setup(singular_plot, singular_sampled, 1);

    // MUSIC method implementation
    // weighing noise eigenvectors and obtaining a_i coefficients:
    // PSD ~ 1 / { sum_i w_i |e.sv_i|^2 } = 1 / { sum_i a_i cos(2pi f T i) }

    double *eigenvalues = singular_sampled.samples;
    auto weight_func = [&] (size_t eigenvect) {
        switch (method)
        {
        case FEM_MUSIC:
        case FEM_PHD:
        default:
            return 1.;
        case FEM_EV:
            return 1. / eigenvalues[eigenvect];
        }
    };

    for (size_t p = 0; p < acm_order; p++)
    {
        double c = 0;
        for (size_t i = 2; i < acm_order; i++) // for each noise eigenvector
        {
            for (size_t j = 0; j < acm_order - p; j++)
            {
                c += weight_func(i) * singular_u[j * (acm_order) + i] * singular_u[(j + p) * (acm_order) + i];
            }
        }
        psd_polynomial[p] = ((p == 0) ? c : 2 * c);
    }

    phd::continuous_t one_to_psd = [&] (double f)
    {
        double psd = 0;
        for (size_t p = 0; p < acm_order; p++)
        {
            psd += psd_polynomial[p] * std::cos(2 * M_PI * f * sampling_period * p);
        }
        return psd;
    };

    // PSD calculation and discrete-freq max search
    // to estimate an interval for iterative max search

    double one_to_psd_min_freq, one_to_psd_min = (std::numeric_limits<double>::max)();
    phd::map(one_to_psd_sampled, [&] (size_t i, double d) {
        double f = (double) (int(i) - int(psd_sample_count)) / (psd_sample_count_x2 * sampling_period);
        double val = one_to_psd(f);
        if (one_to_psd_min > val)
        {
            one_to_psd_min = val;
            one_to_psd_min_freq = f;
        }
        return val;
    });

    phd::map(psd_sampled, one_to_psd_sampled, [] (size_t, double, double s) {
        return 1./ s;
    });

    phd::setup(psd_plot,        psd_sampled, - int(psd_sample_count));
    phd::setup(one_to_psd_plot, one_to_psd_sampled, - int(psd_sample_count));

    // PSD iterative max search / PHD method

    if (method == FEM_PHD)
    {
        assert(acm_order == 3);

        // Solve an equation az^2 + bz + c == 0
        phd::solve_quad(singular_u[2 * 3 + 2],
                        singular_u[1 * 3 + 2],
                        singular_u[0 * 3 + 2],
                        roots_r, roots_i);

        output_frequency = abs(acos(roots_r[0])) / (2 * M_PI * sampling_period);
        output_bias      = abs(signal_frequency - output_frequency);
    }
    else
    {
        one_to_psd_min_freq = phd::find_min(one_to_psd,
                                            one_to_psd_min_freq - one_to_psd_sampled.period,
                                            one_to_psd_min_freq + one_to_psd_sampled.period);

        output_frequency = abs(one_to_psd_min_freq);
        output_bias = abs(signal_frequency - output_frequency);

        roots_r[0] = roots_r[1] = std::cos(2 * M_PI * output_frequency * sampling_period);
        roots_i[0] = -(roots_i[1] = std::sin(2 * M_PI * output_frequency * sampling_period));
    }

    phd::setup(estimated_frequency_plot, roots_r, roots_i, 2, false);

    UpdateData(FALSE);

    phd::free_sampled(signal_sampled);
    phd::free_sampled(noise_sampled);
    phd::free_sampled(autocorrelation_sampled);
    phd::free_sampled(singular_sampled);
    phd::free_sampled(psd_sampled);
    phd::free_sampled(one_to_psd_sampled);
    delete[] singular_u;
    delete[] singular_v;
    delete[] psd_polynomial;

    // redraw

    signal_plot_ctrl.Invalidate();
    autocorrelation_plot_ctrl.Invalidate();
    singular_plot_ctrl.Invalidate();
    estimated_frequency_plot_ctrl.Invalidate();
    psd_plot_ctrl.Invalidate();
    one_to_psd_plot_ctrl.Invalidate();

    UpdateData(FALSE);
}


BOOL ManualViewDialog::OnInitDialog()
{
    CDialogEx::OnInitDialog();

    OnBnsClickedMethod(0);

    //UpdateData(FALSE);
    //UpdateData(TRUE);

    return TRUE;  // return TRUE unless you set the focus to a control
    // EXCEPTION: OCX Property Pages should return FALSE
}


void ManualViewDialog::OnBnsClickedMethod(UINT nID)
{
    UpdateData(TRUE);
    switch (method)
    {
    case FEM_MUSIC:
    case FEM_EV:
    default:
        acm_order_edit_ctrl.SetReadOnly(FALSE);
        break;
    case FEM_PHD:
        acm_order_edit_ctrl.SetReadOnly(TRUE);
        break;
    }
}
