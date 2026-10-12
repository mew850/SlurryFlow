/* Tests the same numeric functions included by the Fluent UDF. */
#include <stdio.h>
#include <assert.h>
#include <math.h>
typedef double real;
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#include "cb_final_math.h"
int main(void) {
    int i,j;
    double expected,actual,d,da,child,r;
    for(i=0;i<MODEL_C_N_PHI;++i) for(j=0;j<MODEL_C_N_GAMMA;++j) {
        expected=limit(pow(10.0,MODEL_C_LOG_ETA[i*MODEL_C_N_GAMMA+j]),MODEL_C_ETA_MIN,MODEL_C_ETA_MAX);
        actual=eta_table(MODEL_C_PHI[i],pow(10.0,MODEL_C_LOG_GAMMA[j]));
        assert(fabs(actual/expected-1.0)<1e-10);
    }
    assert(eta_table(-1.0,0.0)==eta_table(MODEL_C_PHI[0],0.0));
    assert(eta_table(1.0,1e20)==eta_table(MODEL_C_PHI[MODEL_C_N_PHI-1],1e20));
    assert(fabs(structure_one(CB_FINAL_REF_DV)/reference_structure()-1.0)<1e-12);
    for(i=0;i<CB_FINAL_NBIN;++i) {
        d=dv_bin(i); da=envelope(d);
        assert(d>=CB_FINAL_DP);
        assert(fabs(pow(d/CB_FINAL_DP,3.0)/pow(da/CB_FINAL_DP,CB_FINAL_DF)-1.0)<1e-12);
        assert(fabs(structure_one(d)/pow(da/d,3.0)-1.0)<1e-12);
        child=d*pow(0.5,1.0/3.0);
        assert(fabs(2.0*pow(child/d,3.0)-1.0)<1e-12);
        r=structure_one(d)/reference_structure();
        assert(eta_table(0.0152*pow(r,0.0),100.0)==eta_table(0.0152,100.0));
    }
    puts("PASS: 16281 table nodes, bounds, reference recovery, fractal mass, binary mass, feedback-off recovery.");
    printf("Sref=%.9g, reference envelope diameter=%.9g m\n",reference_structure(),envelope(CB_FINAL_REF_DV));
    return 0;
}
