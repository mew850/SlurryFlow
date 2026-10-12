#ifndef CB_FINAL_CONFIG_H
#define CB_FINAL_CONFIG_H
/* Demonstration parameters, NOT a material calibration. */
#define CB_FINAL_NBIN 15
#define CB_FINAL_PHASE 2
#define CB_FINAL_GRAPHITE 1
/* Verify against Fluent's actual Bin-i order before accepting results. */
#define CB_FINAL_BIN0_LARGEST 1
#define CB_FINAL_DF 2.3
#define CB_FINAL_DP 1.0e-7
#define CB_FINAL_REF_DV 8.483428982e-7
/* A positive override is S_ref from a measured reference PSD. */
#define CB_FINAL_SREF_OVERRIDE 0.0
#define CB_FINAL_FEEDBACK 1
#define CB_FINAL_LAMBDA 0.3
#define CB_FINAL_PROBE_GAMMA 10.0
#define CB_FINAL_ES 0.5
#define CB_FINAL_EB 1.0
#define CB_FINAL_T 298.15
#define CB_FINAL_BREAK_C 0.01
#define CB_FINAL_BREAK_DREF 1.0e-6
#define CB_FINAL_BREAK_GREF 100.0
#define CB_FINAL_BREAK_P 1.0
#define CB_FINAL_BREAK_Q 1.0
/* 1: effective Model C viscosity in Brownian kernel; 0: measured carrier. */
#define CB_FINAL_BROWNIAN_EFFECTIVE 1
#define CB_FINAL_CARRIER_MU 0.01
static const double cb_final_dv_ascending[CB_FINAL_NBIN] = {
  1.000000000e-7,1.533598213e-7,2.351952636e-7,
  3.606749764e-7,5.531681198e-7,8.483428982e-7,
  1.301180000e-6,1.995262315e-6,3.059949687e-6,
  4.692762460e-6,7.196856730e-6,1.104089514e-5,
  1.693813981e-5,2.598526445e-5,3.986071706e-5
};
#endif
