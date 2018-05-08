#pragma once

#include "PlotStatic.h"
#include "plot.h"

#include "common.h"
#include "afxwin.h"


// ManualViewDialog dialog

enum frequency_estimation_method
{
    FEM_MUSIC = 0,
    FEM_EV,
    FEM_PHD
};

class ManualViewDialog : public CDialogEx
{
	DECLARE_DYNAMIC(ManualViewDialog)

public:
	ManualViewDialog(CWnd* pParent = NULL);   // standard constructor
	virtual ~ManualViewDialog();

// Dialog Data
	enum { IDD = IDD_MANUALVIEWDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
    PlotStatic signal_plot_ctrl;
    PlotStatic autocorrelation_plot_ctrl;
    PlotStatic singular_plot_ctrl;
    PlotStatic estimated_frequency_plot_ctrl;
    PlotStatic psd_plot_ctrl;
    PlotStatic one_to_psd_plot_ctrl;
    afx_msg void OnBnClickedApply();
    afx_msg void OnBnsClickedMethod(UINT nID);
    double signal_magnitude;
    double signal_frequency;
    double signal_to_noise_ratio_percents;
    size_t acm_order;
    size_t sample_count;
    double sampling_period;
    virtual BOOL OnInitDialog();
    int phd_weight;
    double output_frequency;
    double output_bias;
    int method;
    CEdit acm_order_edit_ctrl;
};
