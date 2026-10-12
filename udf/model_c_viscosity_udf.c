/* Fluent property UDF for the tabulated hybrid microrheology Model C.
 *
 * UDM contract:
 *   UDM 0 = local nominal CB volume fraction
 *   UDM 1 = dimensionless PBM structure ratio F_struct/F_struct,ref
 *   UDM 2 = effective CB fraction (diagnostic, written here)
 *   UDM 3 = viscosity in Pa.s (diagnostic, written here)
 *
 * The structure ratio should be one when PBM structure equals the rheometer
 * reference state.  Populate UDM 0/1 in a case-specific DEFINE_ADJUST UDF.
 */
#include "udf.h"
#include "generated/model_c_table.h"

#define UDM_PHI_CB 0
#define UDM_STRUCTURE_RATIO 1
#define UDM_PHI_EFFECTIVE 2
#define UDM_VISCOSITY 3
#define PHI_PACKING_LIMIT 0.599999

static real clamp_real(real value, real lower, real upper)
{
    return MAX(lower, MIN(value, upper));
}

static int lower_index(const real *grid, int n, real value)
{
    int lo = 0, hi = n - 1, mid;
    if (value <= grid[0]) return 0;
    if (value >= grid[n - 1]) return n - 2;
    while (hi - lo > 1) {
        mid = (lo + hi) / 2;
        if (grid[mid] <= value) lo = mid;
        else hi = mid;
    }
    return lo;
}

static real model_c_lookup(real phi_cb, real shear_rate)
{
    int ip, ig;
    real p, lg, wp, wg, z00, z10, z01, z11, log_eta;
    p = clamp_real(phi_cb, MODEL_C_PHI[0], MODEL_C_PHI[MODEL_C_N_PHI - 1]);
    lg = log10(clamp_real(shear_rate,
        pow(10.0, MODEL_C_LOG_GAMMA[0]),
        pow(10.0, MODEL_C_LOG_GAMMA[MODEL_C_N_GAMMA - 1])));
    ip = lower_index(MODEL_C_PHI, MODEL_C_N_PHI, p);
    ig = lower_index(MODEL_C_LOG_GAMMA, MODEL_C_N_GAMMA, lg);
    wp = (p - MODEL_C_PHI[ip]) / (MODEL_C_PHI[ip + 1] - MODEL_C_PHI[ip]);
    wg = (lg - MODEL_C_LOG_GAMMA[ig]) /
         (MODEL_C_LOG_GAMMA[ig + 1] - MODEL_C_LOG_GAMMA[ig]);
    z00 = MODEL_C_LOG_ETA[ip * MODEL_C_N_GAMMA + ig];
    z10 = MODEL_C_LOG_ETA[(ip + 1) * MODEL_C_N_GAMMA + ig];
    z01 = MODEL_C_LOG_ETA[ip * MODEL_C_N_GAMMA + ig + 1];
    z11 = MODEL_C_LOG_ETA[(ip + 1) * MODEL_C_N_GAMMA + ig + 1];
    log_eta = (1.0-wp)*(1.0-wg)*z00 + wp*(1.0-wg)*z10
            + (1.0-wp)*wg*z01 + wp*wg*z11;
    return clamp_real(pow(10.0, log_eta), MODEL_C_ETA_MIN, MODEL_C_ETA_MAX);
}

DEFINE_PROPERTY(model_c_viscosity, cell, thread)
{
    real phi_cb = MAX(C_UDMI(cell, thread, UDM_PHI_CB), 0.0);
    real structure = C_UDMI(cell, thread, UDM_STRUCTURE_RATIO);
    real phi_effective, viscosity;
    if (structure != structure || structure <= 0.0) structure = 1.0;
    phi_effective = clamp_real(phi_cb * structure, 0.0, PHI_PACKING_LIMIT);
    viscosity = model_c_lookup(phi_effective, C_STRAIN_RATE_MAG(cell, thread));
    C_UDMI(cell, thread, UDM_PHI_EFFECTIVE) = phi_effective;
    C_UDMI(cell, thread, UDM_VISCOSITY) = viscosity;
    return viscosity;
}
