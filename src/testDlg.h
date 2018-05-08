
// testDlg.h : header file
//

#pragma once
#include "afxwin.h"
#include "afxcmn.h"

#include "common.h"

#include "ManualViewDialog.h"
#include "PlotStatic.h"


// CtestDlg dialog
class CtestDlg : public CDialog
{
// Construction
public:
	CtestDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	enum { IDD = IDD_TEST_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	CWnd plot;
    afx_msg void OnBnClickedManual();
    ManualViewDialog manual_view_dialog;
    double signal_magnitude_from;
    double signal_magnitude_to;
    double signal_frequency_from;
    double signal_frequency_to;
    double snr_percents_from;
    double snr_percents_to;
    size_t snr_percents_step_count;
    size_t acm_order;
    size_t number_of_experiments;
    afx_msg void OnBnClickedStartPause();
    afx_msg void OnBnClickedStop();
    PlotStatic frequency_bias_last_plot_ctrl;
    PlotStatic frequency_bias_avg_plot_ctrl;
    size_t sample_count;
    double sampling_period;
    void DoWork();
    LRESULT OnUpdateDataMessage(WPARAM wpD, LPARAM lpD);
    void RequestUpdateData(BOOL saveAndValidate);
    void DoWork0();
    double progress;
};
