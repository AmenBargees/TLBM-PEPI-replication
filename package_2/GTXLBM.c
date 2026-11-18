#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

// float version

#define FLOAT float
#define POW   powf
#define EXP   expf
#define SQRT  sqrtf


void GTXLBM_32bit(
    char *state, int na, int D, int nx, int ny, int nz, int mx, int my, int mz, int xedge, int yedge, int zedge,
    FLOAT w[na], int c[na][D], bool shifted_distribution, bool overlap, bool interior,
    FLOAT rho[ny][nz][nx], FLOAT eta[ny][nz][nx], FLOAT u[ny][nz][D][nx],
    FLOAT f_last[ny][nz][na][nx], FLOAT f[ny][nz][na][nx], FLOAT g_last[ny][nz][na][nx], FLOAT g[ny][nz][na][nx],
    int x_xa[na][nx], int y_ya[na][ny], int z_za[na][nz],
    FLOAT nu_f[nz], FLOAT tau_g[nz], FLOAT beta_g[nz],
    char *Tdep_viscosity, char *rheology,
    FLOAT nu_min, FLOAT gamma_T, FLOAT Tau_y[nz], FLOAT N[nz], FLOAT Delta_g[na],
    FLOAT nu_f_last[ny][nz][nx], FLOAT tau_f_last[ny][nz][nx],
    bool boundary_l, int xl, int naL, int nbL, int aL[na], int aiL[na],
    bool boundary_r, int xr, int naR, int nbR, int aR[na], int aiR[na],
    bool boundary_b, int yb, int naB, int nbB, int aB[na], int aiB[na],
    bool boundary_f, int yf, int naF, int nbF, int aF[na], int aiF[na],
    bool boundary_u, int zu, int naU, int nbU, int aU[na], int aiU[na],
    bool boundary_d, int zd, int naD, int nbD, int aD[na], int aiD[na],
    FLOAT eta_u[ny][nx], FLOAT eta_d[ny][nx], bool TTG_partial_melting,
    char *geometry, FLOAT cylindrical_factor[na][nz], FLOAT dx[ny][nz][nx], int rank
)
// ooo

{
#include "GTXLBM_kernel.c"
}

// double version

#define FLOAT double
#define POW   pow
#define EXP   exp
#define SQRT  sqrt

void GTXLBM_64bit(
    char *state, int na, int D, int nx, int ny, int nz, int mx, int my, int mz, int xedge, int yedge, int zedge,
    FLOAT w[na], int c[na][D], bool shifted_distribution, bool overlap, bool interior,
    FLOAT rho[ny][nz][nx], FLOAT eta[ny][nz][nx], FLOAT u[ny][nz][D][nx],
    FLOAT f_last[ny][nz][na][nx], FLOAT f[ny][nz][na][nx], FLOAT g_last[ny][nz][na][nx], FLOAT g[ny][nz][na][nx],
    int x_xa[na][nx], int y_ya[na][ny], int z_za[na][nz],
    FLOAT nu_f[nz], FLOAT tau_g[nz], FLOAT beta_g[nz],
    char *Tdep_viscosity, char *rheology,
    FLOAT nu_min, FLOAT gamma_T, FLOAT Tau_y[nz], FLOAT N[nz], FLOAT Delta_g[na],
    FLOAT nu_f_last[ny][nz][nx], FLOAT tau_f_last[ny][nz][nx],
    bool boundary_l, int xl, int naL, int nbL, int aL[na], int aiL[na],
    bool boundary_r, int xr, int naR, int nbR, int aR[na], int aiR[na],
    bool boundary_b, int yb, int naB, int nbB, int aB[na], int aiB[na],
    bool boundary_f, int yf, int naF, int nbF, int aF[na], int aiF[na],
    bool boundary_u, int zu, int naU, int nbU, int aU[na], int aiU[na],
    bool boundary_d, int zd, int naD, int nbD, int aD[na], int aiD[na],
    FLOAT eta_u[ny][nx], FLOAT eta_d[ny][nx], bool TTG_partial_melting,
    char *geometry, FLOAT cylindrical_factor[na][nz], FLOAT dx[ny][nz][nx], int rank
)
// ooo
{
#include "GTXLBM_kernel.c"
}
