/* Mixture / homogeneous Discrete CB PBM. Fixed fractal dimension.
 * Compile ONLY this C file for the final model (replace old hooks).
 * Headers: udf.h, sg_mphase.h, sg_pb.h, cb_final_config.h, model_c_table.h.
 * Allocate 24 mixture UDMs globally. Three non-granular phases, no slip.
 * PBM coordinate: solid-volume-equivalent diameter, not envelope diameter.
 * Binary equal-volume daughter PDF is selected in Fluent, not here.
 */
#include "udf.h"
#include "sg_mphase.h"
#include "sg_pb.h"
#include "cb_final_config.h"
#include "model_c_table.h"
#include <math.h>

#include "cb_final_math.h"

#if !RP_HOST
static int phases(Thread *tm, Thread **cb, Thread **gr, Thread **carrier) {
    Thread *tp;
    int k,n=0;
    if(tm==NULL || THREAD_SUB_THREADS(tm)==NULL) return 0;
    sub_thread_loop(tp,tm,k) { ++n; }
    if(n!=3 || CB_FINAL_PHASE<1 || CB_FINAL_PHASE>=n ||
       CB_FINAL_GRAPHITE<1 || CB_FINAL_GRAPHITE>=n || CB_FINAL_PHASE==CB_FINAL_GRAPHITE) return 0;
    *carrier=THREAD_SUB_THREAD(tm,0);
    *cb=THREAD_SUB_THREAD(tm,CB_FINAL_PHASE);
    *gr=THREAD_SUB_THREAD(tm,CB_FINAL_GRAPHITE);
    return *cb!=NULL && *gr!=NULL && *carrier!=NULL;
}

/* Pure evaluation of current solution fields: no iteration-dependent
 * integration, no reference recapture, no dependence on previous UDMs.
 * Status: 1 OK; 2 no CB; 0 initialization; -1 phases; -2 storage;
 * -3 invalid PSD; -4 inconsistent raw sum; -5 invalid phase fractions. */
static double evaluate(cell_t c, Thread *tm, int output) {
    Thread *cb,*gr,*carrier;
    double raw[CB_FINAL_NBIN],sum=0.0,neg=0.0,f,d,da;
    double s=0.0,sref=reference_structure(),ratio=1.0,phi=0.0;
    double ag=0.0,ap=1.0,gamma=0.0,eff,eta,d32v=0.0,d32a=0.0;
    double mv2=0.0,ma2=0.0,ma3=0.0,n0=0.0,largest=0.0,smallest=0.0;
    double expected=0.0,sumerr=0.0;
    int i,status=1,pflag,gflag;
    if(!phases(tm,&cb,&gr,&carrier)) status=-1;
    else if(!Data_Valid_P()) status=0;
    else {
        phi=C_VOF(c,cb); ag=C_VOF(c,gr); ap=C_VOF(c,carrier);
        gamma=C_STRAIN_RATE_MAG(c,tm);
        if(!isfinite(gamma) || gamma<0.0) gamma=0.0;
        if(!isfinite(phi) || !isfinite(ag) || !isfinite(ap) || phi<0.0 || phi>1.0) status=-5;
        else if(phi<1e-12) status=2;
        else {
            for(i=0;i<CB_FINAL_NBIN;++i) {
                if(THREAD_STORAGE(cb,SV_PB_DISC_I(i))==NULL) { status=-2; break; }
                raw[i]=C_PB_DISCI(c,cb,i);
                if(!isfinite(raw[i])) { status=-3; break; }
                sum+=raw[i]; neg+=MAX(-raw[i],0.0);
            }
            if(status==1) {
                /* Accept either globally scaled or phase-normalized storage;
                 * retain the raw sum and reject large conservation errors. */
                expected=fabs(sum-1.0)<fabs(sum-phi)/MAX(phi,1e-12) ? 1.0 : phi;
                sumerr=fabs(sum-expected)/MAX(expected,1e-12);
                if(sum<=0.0 || neg>1e-8*MAX(expected,1e-12)) status=-3;
                else if(sumerr>0.01) status=-4;
                else {
                    for(i=0;i<CB_FINAL_NBIN;++i) {
                        f=MAX(raw[i],0.0)/(sum+neg); d=dv_bin(i); da=envelope(d);
                        s+=f*structure_one(d);
                        mv2+=f/d; ma2+=f*da*da/(d*d*d); ma3+=f*da*da*da/(d*d*d);
                        n0+=phi*f/(0.5235987755982989*d*d*d);
                        if(i==(CB_FINAL_BIN0_LARGEST?0:CB_FINAL_NBIN-1)) largest=f;
                        if(i==(CB_FINAL_BIN0_LARGEST?CB_FINAL_NBIN-1:0)) smallest=f;
                    }
                    ratio=s/sref; d32v=1.0/mv2; d32a=ma3/ma2;
                }
            }
        }
    }
    if(status!=1) { s=sref; ratio=1.0; }
    if(!isfinite(phi) || phi<0.0) phi=0.0;
    eff=phi*pow(ratio,CB_FINAL_FEEDBACK ? CB_FINAL_LAMBDA : 0.0);
    eta=eta_table(eff,gamma);
    pflag=eff<MODEL_C_PHI[0]?-1:(eff>MODEL_C_PHI[MODEL_C_N_PHI-1]?1:0);
    gflag=gamma<pow(10.0,MODEL_C_LOG_GAMMA[0])?-1:(gamma>pow(10.0,MODEL_C_LOG_GAMMA[MODEL_C_N_GAMMA-1])?1:0);
    if(output && n_udm>=24 && THREAD_STORAGE(tm,SV_UDM_I)!=NULL) {
        C_UDMI(c,tm,0)=phi; C_UDMI(c,tm,1)=ag; C_UDMI(c,tm,2)=eff;
        C_UDMI(c,tm,3)=eta; C_UDMI(c,tm,4)=status; C_UDMI(c,tm,5)=gamma;
        C_UDMI(c,tm,6)=pflag; C_UDMI(c,tm,7)=gflag;
        C_UDMI(c,tm,8)=ap; C_UDMI(c,tm,9)=phi+ag+ap;
        C_UDMI(c,tm,10)=sum; C_UDMI(c,tm,11)=s; C_UDMI(c,tm,12)=sref;
        C_UDMI(c,tm,13)=ratio; C_UDMI(c,tm,14)=eff;
        C_UDMI(c,tm,15)=limit(eff,MODEL_C_PHI[0],MODEL_C_PHI[MODEL_C_N_PHI-1]);
        C_UDMI(c,tm,16)=d32v; C_UDMI(c,tm,17)=d32a;
        C_UDMI(c,tm,18)=largest; C_UDMI(c,tm,19)=pflag!=0;
        C_UDMI(c,tm,20)=eta_table(eff,CB_FINAL_PROBE_GAMMA);
        C_UDMI(c,tm,21)=n0; C_UDMI(c,tm,22)=sumerr; C_UDMI(c,tm,23)=smallest;
    }
    return eta;
}
static void update_all(Domain *domain) {
    Thread *tm; cell_t c;
    thread_loop_c(tm,domain) {
        if(!FLUID_THREAD_P(tm)) continue;
        begin_c_loop(c,tm) { evaluate(c,tm,1); } end_c_loop(c,tm)
    }
}
#endif

DEFINE_PROPERTY(cb_final_viscosity, c, phase_thread)
{
#if !RP_HOST
    Thread *tm=THREAD_SUPER_THREAD(phase_thread);
    if(tm==NULL) {
        Error("CB final viscosity must be hooked to a multiphase primary material.\n");
        return 1.0;
    }
    return evaluate(c,tm,1);
#else
    return 1.0;
#endif
}
DEFINE_ADJUST(cb_final_adjust, domain)
{
#if !RP_HOST
    update_all(domain);
#endif
}
DEFINE_EXECUTE_AT_END(cb_final_end)
{
#if !RP_HOST
    update_all(Get_Domain(1));
#endif
}
DEFINE_ON_DEMAND(cb_final_check)
{
    int i;
    Message0("\nCB final: 24 UDMs, 3 phases, 15 CB bins, no slip.\n");
    Message0("feedback=%d lambda=%g df=%g dp=%g Sref=%g\n",CB_FINAL_FEEDBACK,CB_FINAL_LAMBDA,CB_FINAL_DF,CB_FINAL_DP,reference_structure());
    for(i=0;i<CB_FINAL_NBIN;++i) Message0("assumed Bin-%d: dv=%g m, da=%g m\n",i,dv_bin(i),envelope(dv_bin(i)));
#if !RP_HOST
    update_all(Get_Domain(1));
#endif
}
DEFINE_PB_COALESCENCE_RATE(cb_final_aggregation, c, t1, d1, t2, d2)
{
#if !RP_HOST
    Thread *tm=THREAD_SUPER_THREAD(t1);
    double a,b,g,mu;
    if(tm==NULL || !Data_Valid_P() || d1<CB_FINAL_DP || d2<CB_FINAL_DP) return 0.0;
    a=envelope(d1); b=envelope(d2);
    g=MAX(C_STRAIN_RATE_MAG(c,tm),0.0);
    mu=CB_FINAL_BROWNIAN_EFFECTIVE ? evaluate(c,tm,0) : CB_FINAL_CARRIER_MU;
    mu=MAX(mu,1e-6);
    return CB_FINAL_ES*g*(a+b)*(a+b)*(a+b)/6.0
      +CB_FINAL_EB*2.0*1.380649e-23*CB_FINAL_T/(3.0*mu)*(a+b)*(1.0/a+1.0/b);
#else
    return 0.0;
#endif
}
DEFINE_PB_BREAK_UP_RATE_FREQ(cb_final_breakage, c, t, d)
{
#if !RP_HOST
    Thread *tm=THREAD_SUPER_THREAD(t);
    double g;
    if(tm==NULL || !Data_Valid_P() || d<pow(2.0,1.0/3.0)*cb_final_dv_ascending[0]) return 0.0;
    g=MAX(C_STRAIN_RATE_MAG(c,tm),0.0);
    return CB_FINAL_BREAK_C*CB_FINAL_BREAK_GREF*pow(g/CB_FINAL_BREAK_GREF,CB_FINAL_BREAK_P)
      *MAX(pow(d/CB_FINAL_BREAK_DREF,CB_FINAL_BREAK_Q)-1.0,0.0);
#else
    return 0.0;
#endif
}
