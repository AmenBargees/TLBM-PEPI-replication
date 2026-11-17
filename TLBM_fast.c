#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))
#define clip(x, low, high) (max(low, min(x, high)))

//-------------------------------------------------------------------------
// Minimal NumPy .npy writer for 3-D arrays.
// dtype: "<f4" for float32, "<f8" for float64. Little-endian host assumed.
//-------------------------------------------------------------------------
static void write_npy_3d(const char *path, const void *data,
                         int n0, int n1, int n2, const char *dtype) {
    FILE *fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "write_npy_3d: cannot open %s\n", path); return; }
    char hdr[256];
    int  hlen_body = snprintf(hdr, sizeof(hdr),
        "{'descr': '%s', 'fortran_order': False, 'shape': (%d, %d, %d), }",
        dtype, n0, n1, n2);
    int  preamble  = 10;                              // 6 magic + 2 ver + 2 hlen
    int  total     = preamble + hlen_body + 1;        // +1 trailing newline
    int  pad       = (64 - (total % 64)) % 64;
    unsigned short hlen    = (unsigned short)(hlen_body + pad + 1);
    unsigned char  magic[6] = {0x93, 'N', 'U', 'M', 'P', 'Y'};
    unsigned char  ver[2]   = {0x01, 0x00};
    fwrite(magic, 1, 6, fp);
    fwrite(ver,   1, 2, fp);
    fputc( hlen        & 0xFF, fp);
    fputc((hlen >> 8)  & 0xFF, fp);
    fwrite(hdr, 1, hlen_body, fp);
    for (int i = 0; i < pad; i++) fputc(' ', fp);
    fputc('\n', fp);
    size_t elem = (dtype[strlen(dtype)-1] == '8') ? 8 : 4;
    fwrite(data, elem, (size_t)n0 * (size_t)n1 * (size_t)n2, fp);
    fclose(fp);
}

// float version

#define FLOAT float
#define POW   powf
#define EXP   expf
#define SQRT  sqrtf


void TLBM_fast32(
    int na, int D, int nx, int ny, int nz, int mx, int my, int mz, int xedge, int yedge, int zedge,
    FLOAT w[na], int c[na][D], bool shifted_distribution, bool overlap, bool interior,
    FLOAT rho[ny][nz][nx], FLOAT eta[ny][nz][nx], FLOAT u[ny][nz][D][nx],
    FLOAT f_last[ny][nz][na][nx], FLOAT f[ny][nz][na][nx], FLOAT g_last[ny][nz][na][nx], FLOAT g[ny][nz][na][nx],
    int x_xa[na][nx], int y_ya[na][ny], int z_za[na][nz],
    FLOAT nu_f[nz], FLOAT tau_g[nz], FLOAT tau_h[nz], FLOAT beta_g[nz],
    char *temperature, char *rheology, char *nu_formula,
    FLOAT nu_min, FLOAT nu_min_cap, FLOAT nu_max_cap, FLOAT gamma_T, FLOAT Tau_y[nz], FLOAT N[nz], FLOAT Delta_g[na],
    FLOAT nu_f_last[ny][nz][nx], FLOAT tau_f_last[ny][nz][nx], FLOAT nu_total[ny][nz][nx],
    FLOAT nu_star,
    FLOAT P1_thresh, FLOAT P2_thresh,
    FLOAT M1_pd, FLOAT M2_pd, FLOAT M3_pd, FLOAT M4_pd,
    FLOAT M1_bs, FLOAT M2_bs, FLOAT M3_bs, FLOAT M4_bs,
    FLOAT rho0_pd, FLOAT rho_melt_pd, FLOAT M_ext_pd[ny][nz][nx],
    FLOAT rho0_bs, FLOAT rho_melt_bs, FLOAT M_ext_bs[ny][nz][nx],
    FLOAT M_ext_bs_LP_TTG[ny][nz][nx], FLOAT M_ext_bs_MP_TTG[ny][nz][nx], FLOAT M_ext_bs_HP_TTG[ny][nz][nx],
    bool mask_bs[ny][nz][nx], bool mask_pd[ny][nz][nx],
    FLOAT g_l,
    FLOAT T_nd[ny][nz][nx],
    FLOAT cumulative_extracted_pd[ny][nz][nx], FLOAT cumulative_extracted_bs[ny][nz][nx],
    FLOAT cumulative_extracted_LP[ny][nz][nx], FLOAT cumulative_extracted_MP[ny][nz][nx], FLOAT cumulative_extracted_HP[ny][nz][nx],
    FLOAT P0[ny][nz][nx], FLOAT rho_M[ny][nz][nx], FLOAT L_ed[ny][nz][nx],
    FLOAT h_last[ny][nz][na][nx], FLOAT h[ny][nz][na][nx],
    FLOAT DR_2, FLOAT T0_p, FLOAT dT_p, FLOAT z_nd[ny][nz][nx],
    FLOAT rho_solid[ny][nz][nx], FLOAT rho_eff[ny][nz][nx],
    bool boundary_l, int xl, int naL, int nbL, int aL[na], int aiL[na],
    bool boundary_r, int xr, int naR, int nbR, int aR[na], int aiR[na],
    bool boundary_b, int yb, int naB, int nbB, int aB[na], int aiB[na],
    bool boundary_f, int yf, int naF, int nbF, int aF[na], int aiF[na],
    bool boundary_u, int zu, int naU, int nbU, int aU[na], int aiU[na],
    bool boundary_d, int zd, int naD, int nbD, int aD[na], int aiD[na],
    FLOAT eta_u[ny][nx], FLOAT eta_d[ny][nx],
    bool  pre_streaming_BC,
    int rank
)

{
#include "TLBM_fast_kernel.c"
}

// double version

#undef FLOAT
#undef POW
#undef EXP
#undef SQRT

#define FLOAT double
#define POW   pow
#define EXP   exp
#define SQRT  sqrt

void TLBM_fast64(
    int na, int D, int nx, int ny, int nz, int mx, int my, int mz, int xedge, int yedge, int zedge,
    FLOAT w[na], int c[na][D], bool shifted_distribution, bool overlap, bool interior,
    FLOAT rho[ny][nz][nx], FLOAT eta[ny][nz][nx], FLOAT u[ny][nz][D][nx],
    FLOAT f_last[ny][nz][na][nx], FLOAT f[ny][nz][na][nx], FLOAT g_last[ny][nz][na][nx], FLOAT g[ny][nz][na][nx],
    int x_xa[na][nx], int y_ya[na][ny], int z_za[na][nz],
    FLOAT nu_f[nz], FLOAT tau_g[nz], FLOAT tau_h[nz], FLOAT beta_g[nz],
    char *temperature, char *rheology, char *nu_formula,
    FLOAT nu_min, FLOAT nu_min_cap, FLOAT nu_max_cap, FLOAT gamma_T, FLOAT Tau_y[nz], FLOAT N[nz], FLOAT Delta_g[na],
    FLOAT nu_f_last[ny][nz][nx], FLOAT tau_f_last[ny][nz][nx], FLOAT nu_total[ny][nz][nx],
    FLOAT nu_star,
    FLOAT P1_thresh, FLOAT P2_thresh,
    FLOAT M1_pd, FLOAT M2_pd, FLOAT M3_pd, FLOAT M4_pd,
    FLOAT M1_bs, FLOAT M2_bs, FLOAT M3_bs, FLOAT M4_bs,
    FLOAT rho0_pd, FLOAT rho_melt_pd, FLOAT M_ext_pd[ny][nz][nx],
    FLOAT rho0_bs, FLOAT rho_melt_bs, FLOAT M_ext_bs[ny][nz][nx],
    FLOAT M_ext_bs_LP_TTG[ny][nz][nx], FLOAT M_ext_bs_MP_TTG[ny][nz][nx], FLOAT M_ext_bs_HP_TTG[ny][nz][nx],
    bool mask_bs[ny][nz][nx], bool mask_pd[ny][nz][nx],
    FLOAT g_l,
    FLOAT T_nd[ny][nz][nx],
    FLOAT cumulative_extracted_pd[ny][nz][nx], FLOAT cumulative_extracted_bs[ny][nz][nx],
    FLOAT cumulative_extracted_LP[ny][nz][nx], FLOAT cumulative_extracted_MP[ny][nz][nx], FLOAT cumulative_extracted_HP[ny][nz][nx],
    FLOAT P0[ny][nz][nx], FLOAT rho_M[ny][nz][nx], FLOAT L_ed[ny][nz][nx],
    FLOAT h_last[ny][nz][na][nx], FLOAT h[ny][nz][na][nx],
    FLOAT DR_2, FLOAT T0_p, FLOAT dT_p, FLOAT z_nd[ny][nz][nx],
    FLOAT rho_solid[ny][nz][nx], FLOAT rho_eff[ny][nz][nx],
    bool boundary_l, int xl, int naL, int nbL, int aL[na], int aiL[na],
    bool boundary_r, int xr, int naR, int nbR, int aR[na], int aiR[na],
    bool boundary_b, int yb, int naB, int nbB, int aB[na], int aiB[na],
    bool boundary_f, int yf, int naF, int nbF, int aF[na], int aiF[na],
    bool boundary_u, int zu, int naU, int nbU, int aU[na], int aiU[na],
    bool boundary_d, int zd, int naD, int nbD, int aD[na], int aiD[na],
    FLOAT eta_u[ny][nx], FLOAT eta_d[ny][nx],
    bool  pre_streaming_BC,
    int rank
)

{
#include "TLBM_fast_kernel.c"
}

// ============================================================================
//                RKTLBM strmac32/64 and TLBM_rest32/64
// ============================================================================

// float version

#undef FLOAT
#undef POW
#undef EXP
#undef SQRT

#define FLOAT float
#define POW   powf
#define EXP   expf
#define SQRT  sqrtf


void RKTLBM_strmac32(
     int   na, int D, int nx, int ny, int nz, int xedge, int yedge, int zedge,
     int   mx, int my, int mz,
     FLOAT gg[ny][nz][na][nx], FLOAT g[ny][nz][na][nx],
     int   x_xa[na][nx], int y_ya[na][ny], int z_za[na][nz],
     FLOAT rho_r[ny][nz][nx], FLOAT rho_b[ny][nz][nx],
     FLOAT ff_r[ny][nz][na][nx], FLOAT f_r[ny][nz][na][nx],
     FLOAT ff_b[ny][nz][na][nx], FLOAT f_b[ny][nz][na][nx],
     char  *rheology,
     FLOAT hh[ny][nz][na][nx], FLOAT h[ny][nz][na][nx],
     FLOAT hhd[ny][nz][na][nx], FLOAT hd[ny][nz][na][nx],
     FLOAT hhb[ny][nz][na][nx], FLOAT hb[ny][nz][na][nx],
     bool  boundary_l, int xl, int naL, int nbL, int aL[na], int aiL[na],
     bool  boundary_r, int xr, int naR, int nbR, int aR[na], int aiR[na],
     bool  boundary_b, int yb, int naB, int nbB, int aB[na], int aiB[na],
     bool  boundary_f, int yf, int naF, int nbF, int aF[na], int aiF[na],
     bool  boundary_u, int zu, int naU, int nbU, int aU[na], int aiU[na],
     bool  boundary_d, int zd, int naD, int nbD, int aD[na], int aiD[na],
     FLOAT BC_r_wetting_u, FLOAT BC_b_wetting_u, FLOAT BC_r_wetting_d, FLOAT BC_b_wetting_d,
     bool  shifted_distribution,
     FLOAT rho_r_init[ny][nz][nx], FLOAT rho_b_init[ny][nz][nx],
     int   rank
)
{
#include "TLBM_strmac_kernel.c"
}

// double version

#undef FLOAT
#undef POW
#undef EXP
#undef SQRT

#define FLOAT double
#define POW   pow
#define EXP   exp
#define SQRT  sqrt

void RKTLBM_strmac64(
     int   na, int D, int nx, int ny, int nz, int xedge, int yedge, int zedge,
     int   mx, int my, int mz,
     FLOAT gg[ny][nz][na][nx], FLOAT g[ny][nz][na][nx],
     int   x_xa[na][nx], int y_ya[na][ny], int z_za[na][nz],
     FLOAT rho_r[ny][nz][nx], FLOAT rho_b[ny][nz][nx],
     FLOAT ff_r[ny][nz][na][nx], FLOAT f_r[ny][nz][na][nx],
     FLOAT ff_b[ny][nz][na][nx], FLOAT f_b[ny][nz][na][nx],
     char  *rheology,
     FLOAT hh[ny][nz][na][nx], FLOAT h[ny][nz][na][nx],
     FLOAT hhd[ny][nz][na][nx], FLOAT hd[ny][nz][na][nx],
     FLOAT hhb[ny][nz][na][nx], FLOAT hb[ny][nz][na][nx],
     bool  boundary_l, int xl, int naL, int nbL, int aL[na], int aiL[na],
     bool  boundary_r, int xr, int naR, int nbR, int aR[na], int aiR[na],
     bool  boundary_b, int yb, int naB, int nbB, int aB[na], int aiB[na],
     bool  boundary_f, int yf, int naF, int nbF, int aF[na], int aiF[na],
     bool  boundary_u, int zu, int naU, int nbU, int aU[na], int aiU[na],
     bool  boundary_d, int zd, int naD, int nbD, int aD[na], int aiD[na],
     FLOAT BC_r_wetting_u, FLOAT BC_b_wetting_u, FLOAT BC_r_wetting_d, FLOAT BC_b_wetting_d,
     bool  shifted_distribution,
     FLOAT rho_r_init[ny][nz][nx], FLOAT rho_b_init[ny][nz][nx],
     int   rank
)
{
#include "TLBM_strmac_kernel.c"
}

// float version

#undef FLOAT
#undef POW
#undef EXP
#undef SQRT

#define FLOAT float
#define POW   powf
#define EXP   expf
#define SQRT  sqrtf


void RKTLBM_rest32(
     int   na, int D, int nx, int ny, int nz, int xedge, int yedge, int zedge,
     int   mx, int my, int mz,
     FLOAT w[na], int c[na][D],
     FLOAT eta[ny][nz][nx],
     FLOAT u[ny][nz][D][nx],
     FLOAT g[ny][nz][na][nx],
     int   x_xa[na][nx], int y_ya[na][ny], int z_za[na][nz],
     FLOAT beta_g,
     FLOAT tau_f_r, FLOAT tau_f_b, FLOAT tau_g[ny][nz][nx],
     FLOAT rho_r[ny][nz][nx], FLOAT rho_b[ny][nz][nx],
     FLOAT f_r[ny][nz][na][nx],
     FLOAT f_b[ny][nz][na][nx],
     FLOAT C_r[na], FLOAT C_b[na],
     FLOAT rho1_p, FLOAT rho2_p, FLOAT rho3_p, FLOAT beta_g_factor, FLOAT g_comp_factor,
     FLOAT FF[D][ny][nz][nx],
     FLOAT A, FLOAT B[na], FLOAT interface_thickness, FLOAT g_comp,
     FLOAT BC_r_wetting_u, FLOAT BC_b_wetting_u, FLOAT BC_r_wetting_d, FLOAT BC_b_wetting_d,
     char *rheology, char *temperature, char *nu_formula, 
     FLOAT tau_h[nz],
     FLOAT nu_f_r[nz], FLOAT nu_f_b[nz], FLOAT nu_star,
     FLOAT nu_min, FLOAT nu_min_cap, FLOAT nu_max_cap, FLOAT gamma_T, FLOAT Tau_y[nz], FLOAT N[nz], FLOAT Delta_g[na],
     FLOAT nu_f_last[ny][nz][nx], FLOAT tau_f_last[ny][nz][nx],
     FLOAT T_nd[ny][nz][nx],
     FLOAT act_E[nz], FLOAT act_V[nz], FLOAT Delta_eta_array[nz], FLOAT eta_p, FLOAT DR_2, FLOAT T0_p, FLOAT dT_p, FLOAT z_nd[ny][nz][nx],
     FLOAT h[ny][nz][na][nx], FLOAT L_ed[ny][nz][nx],
     FLOAT hd[ny][nz][na][nx], FLOAT D_ed[ny][nz][nx],
     FLOAT hb[ny][nz][na][nx], FLOAT D_eb[ny][nz][nx],
     FLOAT P0[ny][nz][nx],FLOAT cumulative_extracted_pd[ny][nz][nx], FLOAT cumulative_extracted_bs[ny][nz][nx],
     FLOAT cumulative_extracted_LP[ny][nz][nx], FLOAT cumulative_extracted_MP[ny][nz][nx], FLOAT cumulative_extracted_HP[ny][nz][nx], 
     FLOAT P1_thresh, FLOAT P2_thresh, 
     FLOAT M1_pd, FLOAT M2_pd, FLOAT M4_pd,
     FLOAT M1_bs, FLOAT M2_bs, FLOAT M4_bs, 
     bool collision, bool mask_bs[ny][nz][nx], bool mask_pd[ny][nz][nx],
     FLOAT M_ext_bs_LP_TTG[ny][nz][nx], FLOAT M_ext_bs_MP_TTG[ny][nz][nx], FLOAT M_ext_bs_HP_TTG[ny][nz][nx], 
     FLOAT rho0_pd, FLOAT rho_melt_pd, FLOAT M_ext_pd[ny][nz][nx],
     FLOAT rho0_bs, FLOAT rho_melt_bs,
     FLOAT rho_solid[ny][nz][nx], FLOAT rho_eff[ny][nz][nx],
     FLOAT g_l, FLOAT rho_M[ny][nz][nx], FLOAT rho_reff[ny][nz][nx],
     bool  melt_extraction, char *melt_position,
     bool  boundary_l, int xl, int naL, int nbL, int aL[na], int aiL[na],
     bool  boundary_r, int xr, int naR, int nbR, int aR[na], int aiR[na],
     bool  boundary_b, int yb, int naB, int nbB, int aB[na], int aiB[na],
     bool  boundary_f, int yf, int naF, int nbF, int aF[na], int aiF[na],
     bool  boundary_u, int zu, int naU, int nbU, int aU[na], int aiU[na],
     bool  boundary_d, int zd, int naD, int nbD, int aD[na], int aiD[na],
     FLOAT eta_u[ny][nx], FLOAT eta_d[ny][nx],
     bool  shifted_distribution,
     FLOAT rho_r_init[ny][nz][nx], FLOAT rho_b_init[ny][nz][nx],
     int rank,
     FLOAT lp_gate, FLOAT mp_gate, FLOAT hp_gate
)
{
#include "TLBM_rest_kernel.c"
}

// double version
#undef FLOAT
#undef POW
#undef EXP
#undef SQRT

#define FLOAT double
#define POW   pow
#define EXP   exp
#define SQRT  sqrt

void RKTLBM_rest64(
     int   na, int D, int nx, int ny, int nz, int xedge, int yedge, int zedge,
     int   mx, int my, int mz,
     FLOAT w[na], int c[na][D],
     FLOAT eta[ny][nz][nx],
     FLOAT u[ny][nz][D][nx],
     FLOAT g[ny][nz][na][nx],
     int   x_xa[na][nx], int y_ya[na][ny], int z_za[na][nz],
     FLOAT beta_g,
     FLOAT tau_f_r, FLOAT tau_f_b, FLOAT tau_g[ny][nz][nx],
     FLOAT rho_r[ny][nz][nx], FLOAT rho_b[ny][nz][nx],
     FLOAT f_r[ny][nz][na][nx],
     FLOAT f_b[ny][nz][na][nx],
     FLOAT C_r[na], FLOAT C_b[na],
     FLOAT rho1_p, FLOAT rho2_p, FLOAT rho3_p, FLOAT beta_g_factor, FLOAT g_comp_factor,
     FLOAT FF[D][ny][nz][nx],
     FLOAT A, FLOAT B[na], FLOAT interface_thickness, FLOAT g_comp,
     FLOAT BC_r_wetting_u, FLOAT BC_b_wetting_u, FLOAT BC_r_wetting_d, FLOAT BC_b_wetting_d,
     char *rheology, char *temperature, char *nu_formula,
     FLOAT tau_h[nz],
     FLOAT nu_f_r[nz], FLOAT nu_f_b[nz], FLOAT nu_star,
     FLOAT nu_min, FLOAT nu_min_cap, FLOAT nu_max_cap, FLOAT gamma_T, FLOAT Tau_y[nz], FLOAT N[nz], FLOAT Delta_g[na],
     FLOAT nu_f_last[ny][nz][nx], FLOAT tau_f_last[ny][nz][nx],
     FLOAT T_nd[ny][nz][nx],
     FLOAT act_E[nz], FLOAT act_V[nz], FLOAT Delta_eta_array[nz], FLOAT eta_p, FLOAT DR_2, FLOAT T0_p, FLOAT dT_p, FLOAT z_nd[ny][nz][nx],
     FLOAT h[ny][nz][na][nx], FLOAT L_ed[ny][nz][nx],
     FLOAT hd[ny][nz][na][nx], FLOAT D_ed[ny][nz][nx],
     FLOAT hb[ny][nz][na][nx], FLOAT D_eb[ny][nz][nx],
     FLOAT P0[ny][nz][nx],FLOAT cumulative_extracted_pd[ny][nz][nx], FLOAT cumulative_extracted_bs[ny][nz][nx],
     FLOAT cumulative_extracted_LP[ny][nz][nx], FLOAT cumulative_extracted_MP[ny][nz][nx], FLOAT cumulative_extracted_HP[ny][nz][nx],
     FLOAT P1_thresh, FLOAT P2_thresh,
     FLOAT M1_pd, FLOAT M2_pd, FLOAT M4_pd,
     FLOAT M1_bs, FLOAT M2_bs, FLOAT M4_bs,
     bool collision, bool mask_bs[ny][nz][nx], bool mask_pd[ny][nz][nx],
     FLOAT M_ext_bs_LP_TTG[ny][nz][nx], FLOAT M_ext_bs_MP_TTG[ny][nz][nx], FLOAT M_ext_bs_HP_TTG[ny][nz][nx],
     FLOAT rho0_pd, FLOAT rho_melt_pd, FLOAT M_ext_pd[ny][nz][nx],
     FLOAT rho0_bs, FLOAT rho_melt_bs,
     FLOAT rho_solid[ny][nz][nx], FLOAT rho_eff[ny][nz][nx],
     FLOAT g_l, FLOAT rho_M[ny][nz][nx], FLOAT rho_reff[ny][nz][nx],
     bool  melt_extraction, char *melt_position,
     bool  boundary_l, int xl, int naL, int nbL, int aL[na], int aiL[na],
     bool  boundary_r, int xr, int naR, int nbR, int aR[na], int aiR[na],
     bool  boundary_b, int yb, int naB, int nbB, int aB[na], int aiB[na],
     bool  boundary_f, int yf, int naF, int nbF, int aF[na], int aiF[na],
     bool  boundary_u, int zu, int naU, int nbU, int aU[na], int aiU[na],
     bool  boundary_d, int zd, int naD, int nbD, int aD[na], int aiD[na],
     FLOAT eta_u[ny][nx], FLOAT eta_d[ny][nx],
     bool  shifted_distribution,
     FLOAT rho_r_init[ny][nz][nx], FLOAT rho_b_init[ny][nz][nx],
     int rank,
     FLOAT lp_gate, FLOAT mp_gate, FLOAT hp_gate
)
{
#include "TLBM_rest_kernel.c"
}

