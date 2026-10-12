#ifndef CB_FINAL_MATH_H
#define CB_FINAL_MATH_H
#include <math.h>
#include "cb_final_config.h"
#include "model_c_table.h"
static double limit(double x, double a, double b) { return MAX(a, MIN(x,b)); }
static double dv_bin(int i) {
    return cb_final_dv_ascending[CB_FINAL_BIN0_LARGEST ? CB_FINAL_NBIN-1-i : i];
}
static double envelope(double d) {
    return CB_FINAL_DP * pow(d / CB_FINAL_DP, 3.0 / CB_FINAL_DF);
}
static double structure_one(double d) {
    return pow(d / CB_FINAL_DP, 9.0 / CB_FINAL_DF - 3.0);
}
static double reference_structure(void) {
    return CB_FINAL_SREF_OVERRIDE > 0.0 ? CB_FINAL_SREF_OVERRIDE : structure_one(CB_FINAL_REF_DV);
}
static int bracket(const real *x, int n, double q) {
    int lo=0, hi=n-1, mid;
    while(hi-lo>1) { mid=(hi+lo)/2; if(x[mid]<=q) lo=mid; else hi=mid; }
    return lo;
}
static double eta_table(double phi, double gamma) {
    int i,j;
    double p,g,a,b,z;
    p=limit(phi,MODEL_C_PHI[0],MODEL_C_PHI[MODEL_C_N_PHI-1]);
    g=log10(MAX(gamma,pow(10.0,MODEL_C_LOG_GAMMA[0])));
    g=limit(g,MODEL_C_LOG_GAMMA[0],MODEL_C_LOG_GAMMA[MODEL_C_N_GAMMA-1]);
    i=bracket(MODEL_C_PHI,MODEL_C_N_PHI,p);
    j=bracket(MODEL_C_LOG_GAMMA,MODEL_C_N_GAMMA,g);
    a=(p-MODEL_C_PHI[i])/(MODEL_C_PHI[i+1]-MODEL_C_PHI[i]);
    b=(g-MODEL_C_LOG_GAMMA[j])/(MODEL_C_LOG_GAMMA[j+1]-MODEL_C_LOG_GAMMA[j]);
    z=(1-a)*(1-b)*MODEL_C_LOG_ETA[i*MODEL_C_N_GAMMA+j]
      +a*(1-b)*MODEL_C_LOG_ETA[(i+1)*MODEL_C_N_GAMMA+j]
      +(1-a)*b*MODEL_C_LOG_ETA[i*MODEL_C_N_GAMMA+j+1]
      +a*b*MODEL_C_LOG_ETA[(i+1)*MODEL_C_N_GAMMA+j+1];
    return limit(pow(10.0,z),MODEL_C_ETA_MIN,MODEL_C_ETA_MAX);
}
#endif
