//-------------------------------------------------------------------------
//  C code for TLBM_rest
//-------------------------------------------------------------------------
//  Archean2 sub-block toggles (1=enabled, 0=disabled). Production = all 1.
//  Companion toggle A2_HSTRM lives in TLBM_strmac_kernel.c.
#define A2_LED_MACRO 1   // L_ed = sum h in macroscopic block
#define A2_HEQ       1   // h_eq, h_neq in equilibrium block
#define A2_HCOLL     1   // h relaxation in collision (needs A2_HEQ)
#define A2_MELT      1   // melt-extraction block (needs A2_LED_MACRO)
#define A2_BUOY      1   // archean2 density branch in buoyancy (needs A2_MELT)
#define TAU_HD     0.52
#ifndef RECOLOR_HD
#define RECOLOR_HD 0   // DISABLED: LK recoloring spinodally decomposes the continuous depletion field (smoke A/B 2026-11-09); FCT-style limiter is the open design item
#ifndef BETA_HD
#define BETA_HD 0.05   // segregation strength for hd/hb, sized against TAU_HD diffusivity
#endif
#endif  // relaxation time of the advected-depletion scalars hd/hb.
                         // tau_h ~ 0.5 is dispersive for the spiky melt sources and
                         // lets the depletion moment oscillate negative, which holds
                         // the melt gate open and runs away; 0.56 adds just enough
                         // diffusion to damp the oscillation.
//-------------------------------------------------------------------------
    int   incx, startx, endx, starty, endy, startz, endz;
    int   x, y, z, a, d, i, j, xx_xa, yy_ya, zz_za;
    bool  edgey, edgez, x_loop, boundary_bb, boundary_ff, boundary_uu, boundary_dd;
    FLOAT tmp_r;
    FLOAT cu[na][nx], u2[nx];
    FLOAT gi[nx];
    FLOAT feq_r[na][nx],feq_b[na][nx], f_eq[na][nx], f_neq[na][nx];
    FLOAT frac_r[nx], frac_b[nx];
    FLOAT g_eq[na][nx];
    FLOAT rho[nx];
    FLOAT g_r[nx], g_b[nx];
    FLOAT tau[ny][nz][nx];
    FLOAT tau_fi[nx], tau_gi[nz];
    FLOAT fm[nx], ttmp[nx], tempp, ff_sum, feq, cm[na];
    FLOAT F[D][nx], cos_phi[na][nx], cos_lambda[na][nx];
    FLOAT zero = 1.e-10;
    FLOAT rho_comp [nx];
    FLOAT delta_rho[nx];
//-------------------------------------------------------------------------
//  Archean2 variables
//-------------------------------------------------------------------------
    FLOAT weight_r[nx], weight_b[nx];
    FLOAT f_r_neq[na][nx], f_b_neq[na][nx];
    FLOAT h_eq[na][nx],  h_neq[na][nx], g_neq[na][nx];
    FLOAT hd_eq[na][nx], hd_neq[na][nx], hb_eq[na][nx], hb_neq[na][nx];
    double strain_rate[nx], S[D][D][nx], S_factor[nx], Tau[nx];
    double nu_f[nx], nu_linear[nx], nu_plastic[nx], nu_total[nx];
    double P[nx], T[nx], T_K[nx];
    double gas_C, ref_K, max_val, arg[nx];
    double eqn1_sol_bs, eqn2_sol_bs, eqn3_sol_bs;
    double eqn1_liq_bs, eqn2_liq_bs;
    double T_solidus_bs, T_liquidus_bs, T_solidus_pd, T_liquidus_pd;
    double M_pd[nx], M_bs[nx];
    FLOAT M_ext_bs[ny][nz][nx];
    double M0_bs, M0_pd;
    double L_den[nx];
    double P_garnet_in[nx], P_plag_out[nx];
    double coef_phase_bs[nx], coef_res_bs [nx], coef_phase_pd[nx], coef_res_pd[nx];
    FLOAT rho_solid_pd [ny][nz][nx], rho_solid_bs[ny][nz][nx], rho_ext_pd[ny][nz][nx];
    bool  T1_cond, P1_cond, mask_lpttg, mask_mpttg, LP_TTG, MP_TTG;
    bool  T2_cond, P2_cond, mask_hpttg, HP_TTG;
    FLOAT tmp;
    int   itmp;     

//-------------------------------------------------------------------------
//  Determine loop ranges for x, y and z
//-------------------------------------------------------------------------
    startx = xedge  ; endx = mx-xedge;
    starty = yedge  ; endy = my-yedge;
    startz = zedge  ; endz = mz-zedge;
//---------------------------------------------------------------
//  Start y and z loops and edge bool for y edges or z edges
//---------------------------------------------------------------
    if(D==2) {
        for (a=1;a<na;a++) {
	    cm[a] = sqrt(c[a][0] * c[a][0] + c[a][1] * c[a][1]);
        }
    }
    else if(D==3) {
        for (a=1;a<na;a++) {
	    cm[a] = sqrt(c[a][0] * c[a][0] + c[a][1] * c[a][1] + c[a][2] * c[a][2]);
        }
    }
//---------------------------------------------------------------
//  Wetting BC: write halo macroscopic densities. The wetting effect
//  propagates through the color gradient F (which uses halo rho via
//  the FD stencil) into the recoloring step's cos(lambda) bias.
//  Standard halfway bounce-back in strmac handles f reflection;
//  halo f is NOT written here.
//
//    rho_r_halo = (rho_r[zu] + rho_b[zu]) * BC_r_wetting   ((1+cos θ)/2)
//    rho_b_halo = (rho_r[zu] + rho_b[zu]) * BC_b_wetting   ((1-cos θ)/2)
//
//  The reference scaling is the LOCAL cell total (rho_r + rho_b) —
//  matches the original V3 formulation. Mora paper uses a constant
//  rho_t = rho_2 (initial majority density), but for incompressible
//  LBM the local total stays close to 1 and the two formulations
//  agree in the bulk. Using local-total keeps consistency with the
//  cell's actual mass at the wall.
//
//  GUARD: each branch requires the corresponding halo (zedge=1 for
//  upper/lower, xedge=1 for left/right). Without the halo the BC
//  would read/write out-of-bounds; under that condition it's a no-op.
//---------------------------------------------------------------
    for(y=0;y<ny;y++) {
        if (boundary_u && zedge == 1) {
            for(x=0;x<nx;x++) {
                /* BC_r_wetting < 0 signals neutral (zero-gradient) ghosts:
                   the fixed 50/50 wetting split (theta = 90) acts as a color
                   SOURCE wherever the local red fraction differs from 1/2
                   (red/ce_bs grew 113 -> 124% over 110 My in M1_hpc). */
                if (BC_r_wetting_u < 0) {
                    rho_r[y][zu-1][x] = rho_r[y][zu][x];
                    rho_b[y][zu-1][x] = rho_b[y][zu][x];
                } else {
                    tmp               = rho_r[y][zu][x] + rho_b[y][zu][x];
                    rho_r[y][zu-1][x] = tmp * BC_r_wetting_u;
                    rho_b[y][zu-1][x] = tmp * BC_b_wetting_u;
                }
            }
        }
        if (boundary_d && zedge == 1) {
            for(x=0;x<nx;x++) {
                if (BC_r_wetting_d < 0) {
                    rho_r[y][zd+1][x] = rho_r[y][zd][x];
                    rho_b[y][zd+1][x] = rho_b[y][zd][x];
                } else {
                    tmp               = rho_r[y][zd][x] + rho_b[y][zd][x];
                    rho_r[y][zd+1][x] = tmp * BC_r_wetting_d;
                    rho_b[y][zd+1][x] = tmp * BC_b_wetting_d;
                }
            }
        }
        if (boundary_l && xedge == 1) {
            for(z=0;z<nz;z++) {
                tmp               = rho_r[y][z][xl] + rho_b[y][z][xl];
                rho_r[y][z][xl-1] = tmp * BC_r_wetting_u;
                rho_b[y][z][xl-1] = tmp * BC_b_wetting_u;
            }
        }
        if (boundary_r && xedge == 1) {
            for(z=0;z<nz;z++) {
                tmp               = rho_r[y][z][xr] + rho_b[y][z][xr];
                rho_r[y][z][xr+1] = tmp * BC_r_wetting_d;
                rho_b[y][z][xr+1] = tmp * BC_b_wetting_d;
            }
        }
    }
//---------------------------------------------------------------
//  Start of loops to do calculations
//---------------------------------------------------------------
//---------------------------------------------------------------
//  Advected depletion D_ed/D_eb: van-Leer-limited upwind advection
//  (2D). Replaces the hd/hb population transport. A stationary
//  parcel holds its depletion EXACTLY (zero flux at u = 0); a
//  moving parcel carries it monotonically with sharp fronts.
//  Lattice units dx = dt = 1; u is the macroscopic velocity.
//  Wall faces carry zero flux; x ghosts come via the rho exchange.
//---------------------------------------------------------------
    if (D == 2 && strcmp(rheology,"archean2") == 0 && melt_extraction) {
        int nmax = (nx > nz) ? nx : nz;
        FLOAT Dold[nmax], slp[nmax];
        for (int fld = 0; fld < 2; fld++) {
            for(y=0;y<ny;y++) {
                for(z=zu;z<=zd;z++) {              /* x sweep */
                    for(x=0;x<nx;x++) Dold[x] = (fld==0)? D_ed[y][z][x] : D_eb[y][z][x];
                    for(x=1;x<nx-1;x++) {
                        FLOAT a1 = Dold[x+1]-Dold[x], b1 = Dold[x]-Dold[x-1];
                        slp[x] = (a1*b1 > 0) ? 2*a1*b1/(a1+b1) : (FLOAT)0.0;
                    }
                    slp[0] = 0; slp[nx-1] = 0;
                    for(x=1;x<nx-1;x++) {
                        FLOAT ufm = (x==1)    ? u[y][z][0][x] : (FLOAT)0.5*(u[y][z][0][x-1]+u[y][z][0][x]);
                        FLOAT ufp = (x==nx-2) ? u[y][z][0][x] : (FLOAT)0.5*(u[y][z][0][x]+u[y][z][0][x+1]);
                        FLOAT Dm = (ufm > 0) ? Dold[x-1]+(FLOAT)0.5*slp[x-1]*(1-ufm)
                                             : Dold[x]  -(FLOAT)0.5*slp[x]  *(1+ufm);
                        FLOAT Dp = (ufp > 0) ? Dold[x]  +(FLOAT)0.5*slp[x]  *(1-ufp)
                                             : Dold[x+1]-(FLOAT)0.5*slp[x+1]*(1+ufp);
                        FLOAT val = Dold[x] - (ufp*Dp - ufm*Dm) + Dold[x]*(ufp-ufm);
                        if (fld==0) D_ed[y][z][x] = val; else D_eb[y][z][x] = val;
                    }
                }
                for(x=1;x<nx-1;x++) {              /* z sweep, zero-flux walls */
                    for(z=zu;z<=zd;z++) Dold[z] = (fld==0)? D_ed[y][z][x] : D_eb[y][z][x];
                    for(z=zu+1;z<zd;z++) {
                        FLOAT a1 = Dold[z+1]-Dold[z], b1 = Dold[z]-Dold[z-1];
                        slp[z] = (a1*b1 > 0) ? 2*a1*b1/(a1+b1) : (FLOAT)0.0;
                    }
                    slp[zu] = 0; slp[zd] = 0;
                    for(z=zu;z<=zd;z++) {
                        FLOAT ufm = (z==zu) ? (FLOAT)0.0 : (FLOAT)0.5*(u[y][z-1][D-1][x]+u[y][z][D-1][x]);
                        FLOAT ufp = (z==zd) ? (FLOAT)0.0 : (FLOAT)0.5*(u[y][z][D-1][x]+u[y][z+1][D-1][x]);
                        FLOAT Dm = (z==zu || ufm==0) ? (FLOAT)0.0 : ((ufm > 0) ? Dold[z-1]+(FLOAT)0.5*slp[z-1]*(1-ufm)
                                                                               : Dold[z]  -(FLOAT)0.5*slp[z]  *(1+ufm));
                        FLOAT Dp = (z==zd || ufp==0) ? (FLOAT)0.0 : ((ufp > 0) ? Dold[z]  +(FLOAT)0.5*slp[z]  *(1-ufp)
                                                                               : Dold[z+1]-(FLOAT)0.5*slp[z+1]*(1+ufp));
                        FLOAT val = Dold[z] - (ufp*Dp - ufm*Dm) + Dold[z]*(ufp-ufm);
                        if (fld==0) D_ed[y][z][x] = val; else D_eb[y][z][x] = val;
                    }
                }
            }
        }
    }

    for(y=starty;y<endy;y++) {
        edgey = (yedge==1&&(y==1||y==my-2)); // Set whether on y edge or not
        if(D==3) {
	    if(boundary_b&&y==yb) { boundary_bb = true;  }
	    else                  { boundary_bb = false; }
	    if(boundary_f&&y==yf) { boundary_ff = true;  }
	    else                  { boundary_ff = false; }
	}
	else {
	    boundary_bb = false;
	    boundary_ff = false;
	}
        for(z=startz;z<endz;z++) {
            edgez = (zedge==1&&(z==1||z==mz-2)); // Set whether on z edge or not
	    if(boundary_u&&z==zu) { boundary_uu = true;  }
	    else                  { boundary_uu = false; }
	    if(boundary_d&&z==zd) { boundary_dd = true;  }
	    else                  { boundary_dd = false; }

	    incx = 1;   x_loop = true;

            if(x_loop) {
//---------------------------------------------------------------
//              Macroscopic properties. Backards a loop to
//              optimize precision by summing small values first
//---------------------------------------------------------------

               for(x=startx;x<endx;x+=incx) {
                    eta[y][z][x] = g[y][z][na-1][x];
                     if (A2_LED_MACRO && strcmp(rheology,"archean2") == 0) {
                         L_ed[y][z][x] = h[y][z][na-1][x];
                         /* D_ed/D_eb are prognostic scalars advected by the
                            TVD pass above; they are no longer hd/hb moments. */
                     }
               }
               for(a=na-2;a>=0;a--) {
                    for(x=startx;x<endx;x+=incx) {
                        eta[y][z][x] += g[y][z][a][x];
                        if (A2_LED_MACRO && strcmp(rheology,"archean2") == 0) {
                            L_ed[y][z][x] += h[y][z][a][x];
                        }
                    }
               }

               for(x=startx;x<endx;x+=incx) {
                   L_ed[y][z][x]  = max(0.0, min(1, L_ed[y][z][x]));
                   D_ed[y][z][x]  = max(0.0, D_ed[y][z][x]);
                   D_eb[y][z][x]  = max(0.0, D_eb[y][z][x]);
                   if (A2_LED_MACRO && strcmp(rheology,"archean2") == 0) {
                       /* Clamp advected depletion to the extraction cap and
                          rescale the populations so the moment stays consistent */
                       if (D_ed[y][z][x] > M4_pd) {D_ed[y][z][x] = M4_pd;}
                       if (D_eb[y][z][x] > M4_bs) {D_eb[y][z][x] = M4_bs;}
                   }
                   rho       [x]  = rho_r[y][z][x] + rho_b[y][z][x];
                   frac_r    [x]  = rho_r[y][z][x]/rho[x]; // Calculate red  fraction
                   frac_b    [x]  = rho_b[y][z][x]/rho[x]; // Calculate blue fraction
                   eta [y][z][x] /= rho[x];
               }

//             u calculation

               for(d=0;d<D;d++) {
                    for(x=startx;x<endx;x+=incx) {
                        u[y][z][d][x] = 0;
                    }
                    for(a=na-1;a>0;a--) {
                        for(x=startx;x<endx;x+=incx) {
                            u[y][z][d][x] += ((f_r[y][z][a][x] + f_b[y][z][a][x]) * c[a][d]);
                        }
                    }
                    for(x=startx;x<endx;x+=incx) {
                        u[y][z][d][x] /= rho[x];
                    }
               }
//---------------------------------------------------------------
//              Thermal boundary conditions
//---------------------------------------------------------------

                if(boundary_uu) { // Upper thermal BC
                    for(x=startx;x<endx;x+=incx) { eta[y][z][x] = eta_u[y][x]; }
                }
                if(boundary_dd) { // Lower thermal BC
                    for(x=startx;x<endx;x+=incx) { eta[y][z][x] = eta_d[y][x]; }
                } 
//---------------------------------------------------------------
//              Archean2 melt extraction
//---------------------------------------------------------------
                if (A2_MELT && strcmp(rheology,"archean2") == 0 && melt_extraction) {
                    for(x=startx;x<endx;x+=incx) {
                        T[x]                   =  eta[y][z][x]/DR_2 + T0_p + dT_p/2;
                        T[x]                   =  max(T0_p, min(T0_p + dT_p, T[x]));
                        P[x]                   =  P0[y][z][x] * rho[x];
                        P[x]                   =  max(0.0, min(137, P[x])); 
                        // Basalt Solidii and Liquidii
                        eqn1_sol_bs            =  699.85 - (70400/(P[x]*1000 + 354)) + (77800000/POW(P[x]*1000 + 354,2));
                        eqn2_sol_bs            =  661.85 + 3.5*P[x] + 6.2*POW(P[x],2);
                        eqn3_sol_bs            =(-1.0116e-12) * POW(P[x], 7) + (8.9986e-10) * POW(P[x], 6) +
                                                (-2.9466e-7 ) * POW(P[x], 5) + (4.781e-5  ) * POW(P[x], 4) +
                                                (-0.0039836)  * POW(P[x], 3) + (0.0072596 ) * POW(P[x], 2) + 36.75 * P[x] + 984.75;
                        if ((P[x] <  P1_thresh))           {T_solidus_bs = eqn1_sol_bs;}
                        // else if( P1_thresh <= P[x] < 5) {T_solidus_bs = eqn2_sol_bs;}
                        else                               {T_solidus_bs = eqn2_sol_bs;}
                        eqn1_liq_bs            = 1149.85 + 105*P[x];
                        eqn2_liq_bs            = (1.3728e-10) * POW(P[x], 6) - (3.7739e-8) * POW(P[x], 5) - (5.0861e-7) * POW(P[x], 4) +
                                                  0.0011277   * POW(P[x], 3) -  0.15346    * POW(P[x], 2) +  23.869 * P[x] + 2580.85;
                        if ((P[x] < 5)) {T_liquidus_bs = eqn1_liq_bs;}
                        else {T_liquidus_bs = eqn2_liq_bs;}
                        // Peridotite Solidii and Liquidii from Denis Andrault Solidus and liquidus profiles of chondritic mantle
                        // Solidus: Sizova et al. (2015) quadratic below 13 GPa, Andrault et al. (2011)
                        // Simon-Glatzel fit above; the two curves join at ~13 GPa (1988 vs 1991 C).
                        // Liquidus: Andrault at all P (1667 C at P=0, 4474 C at 130 GPa); the Sizova
                        // linear 1799.85+114P is only calibrated shallow and diverges at depth.
                        if (P[x] < 13.0) {T_solidus_pd = -5.1 * POW(P[x],2) + 133*P[x] + 1120.85;}
                        else             {T_solidus_pd = (2045 * POW(1.0 + P[x]/92.0,(1/1.3))) - 273.15;}
                        T_liquidus_pd          = (1940 * POW(1.0 + P[x]/29.0,(1/1.9))) - 273.15;
                        // T/P conditions for LP-TTG  
                        T1_cond                = (T[x] > (760 - (60 * POW((P[x] - 1), 2)))) && (T[x] < (1000 - (150 * POW(((P[x] - 1.2) / 1.2), 2))));
                        P1_cond                = (P[x] > (-0.5 * ((T[x] - 870) / 220))) && (P[x] < (1.5 + (0.7 * ((T[x] - 700) / 200))));
                        mask_lpttg             = P[x]  < 1;
                        LP_TTG                 = T1_cond && P1_cond && mask_lpttg;
                        // T/P conditions for MP-TTG
                        mask_mpttg             = (P[x] >= 1);
                        MP_TTG                 = T1_cond && P1_cond && mask_mpttg;
                        // T/P conditions for HP-TTG
                        T2_cond                = (T[x] > 1000) && (T[x] < (1100 + (50 * POW(((P[x] - 3.5) / 3.5), 2))));
                        P2_cond                = (P[x] > (2.35 + 0.15 * ((T[x] - 1000) / 100))) && (P[x] < 5); // from Charitra Jain
                        HP_TTG                 = T2_cond && P2_cond;
                        if (L_ed[y][z][x] <= 1 ) {L_den[x] = L_ed[y][z][x];}
                        else                     {L_den[x] = 1;}
                        if (L_ed[y][z][x] > 0.1) {L_den[x] = L_ed[y][z][x];}
                        else                     {L_den[x] = 0;}

                        // Melt calculation Basalt
                        if ((T[x] > T_solidus_bs) && (T[x] < T_liquidus_bs) && (T_liquidus_bs > T_solidus_bs)){
                             M0_bs = ((T[x] - T_solidus_bs)/(T_liquidus_bs - T_solidus_bs));}
                        else if (T[x] > T_liquidus_bs) {M0_bs = 1;}
                        else {M0_bs = 0;}
                        if ((M0_bs - D_eb[y][z][x]) > 0){M_bs[x] = M0_bs - D_eb[y][z][x];}
                        else {M_bs[x] = 0;}
                        if  ((M_bs[x] > M2_bs) && ((M_bs[x] - M1_bs) < M4_bs) && ((M_bs[x] - M1_bs) > 0)) {
                              M_ext_bs[y][z][x] = ( M_bs[x] - M1_bs) * L_den[x] * mask_bs[y][z][x];}
                        else {M_ext_bs[y][z][x] = 0; }
                        // never extract past the remaining advected capacity
                        tmp = M4_bs - D_eb[y][z][x];
                        if (tmp < 0) {tmp = 0;}
                        if (M_ext_bs[y][z][x] > tmp) {M_ext_bs[y][z][x] = tmp;}

                        // Keep track of melt for LP, MP and HP TTG

                        if   (LP_TTG && !MP_TTG && !HP_TTG) {M_ext_bs_LP_TTG[y][z][x] = M_ext_bs[y][z][x];}
                        else {M_ext_bs_LP_TTG[y][z][x] = 0;}
                        if   (MP_TTG && !LP_TTG && !HP_TTG) {M_ext_bs_MP_TTG[y][z][x] = M_ext_bs[y][z][x];}
                        else {M_ext_bs_MP_TTG[y][z][x] = 0;}
                        if   (HP_TTG && !LP_TTG && !MP_TTG) {M_ext_bs_HP_TTG[y][z][x] = M_ext_bs[y][z][x];}
                        else {M_ext_bs_HP_TTG[y][z][x] = 0;}
                        M_ext_bs_LP_TTG[y][z][x] *= lp_gate;   /* window share caps */
                        M_ext_bs_MP_TTG[y][z][x] *= mp_gate;
                        M_ext_bs_HP_TTG[y][z][x] *= hp_gate;
                        M_ext_bs[y][z][x] = M_ext_bs_LP_TTG[y][z][x] + M_ext_bs_MP_TTG[y][z][x] + M_ext_bs_HP_TTG[y][z][x];
                        
                        if ((D_eb[y][z][x] +  M_ext_bs[y][z][x]) < M4_bs){
                             cumulative_extracted_bs[y][z][x] += M_ext_bs[y][z][x];
                             cumulative_extracted_LP[y][z][x] += M_ext_bs_LP_TTG[y][z][x];
                             cumulative_extracted_MP[y][z][x] += M_ext_bs_MP_TTG[y][z][x];
                             cumulative_extracted_HP[y][z][x] += M_ext_bs_HP_TTG[y][z][x];
                             }
                        else {M_ext_bs[y][z][x] = 0;
                              M_ext_bs_LP_TTG[y][z][x] = 0;
                              M_ext_bs_MP_TTG[y][z][x] = 0;
                              M_ext_bs_HP_TTG[y][z][x] = 0;}
                        if (cumulative_extracted_bs[y][z][x] < 0) {cumulative_extracted_bs[y][z][x] = 0;} /* unclamped: global DI/window gates limit production */
                        D_eb[y][z][x] += M_ext_bs[y][z][x];   /* source into the advected scalar */
                        // Cap the per-type counters like the total: without this the
                        // LP/MP/HP counters run past M4 when the advected D_eb drains
                        // away from a persistent melting site (ce_MP reached 406%).
                        if (cumulative_extracted_LP[y][z][x] < 0) {cumulative_extracted_LP[y][z][x] = 0;}
                        if (cumulative_extracted_MP[y][z][x] < 0) {cumulative_extracted_MP[y][z][x] = 0;}
                        if (cumulative_extracted_HP[y][z][x] < 0) {cumulative_extracted_HP[y][z][x] = 0;}

                        // Melt calculation Peridotite
                        if ((T[x] > T_solidus_pd) && (T[x] < T_liquidus_pd) && (T_liquidus_pd > T_solidus_pd)){
                             M0_pd  = ((T[x] - T_solidus_pd)/(T_liquidus_pd - T_solidus_pd));}
                        else if (T[x] > T_liquidus_pd) {M0_pd = 1.0;}
                        else {M0_pd = 0;}
                        if ((M0_pd - D_ed[y][z][x]) > 0){M_pd[x] = M0_pd - D_ed[y][z][x];}
                        else {M_pd[x] = 0;}
                        if ((M_pd[x] > M2_pd) && ((M_pd[x] - M1_pd) < M4_pd) && ((M_pd[x] - M1_pd) > 0)) {
                             M_ext_pd[y][z][x] = ( M_pd[x] - M1_pd) * (1 - L_den[x]) * mask_pd[y][z][x];}
                        else {M_ext_pd[y][z][x] = 0;}
                        if (M_ext_pd[y][z][x] < 0) {M_ext_pd[y][z][x] = 0;}
                        // never extract past the remaining advected capacity
                        tmp = M4_pd - D_ed[y][z][x];
                        if (tmp < 0) {tmp = 0;}
                        if (M_ext_pd[y][z][x] > tmp) {M_ext_pd[y][z][x] = tmp;}
                        if ((D_ed[y][z][x] +  M_ext_pd[y][z][x]) < M4_pd) {
                                 cumulative_extracted_pd[y][z][x] += M_ext_pd[y][z][x];}
                        else {M_ext_pd[y][z][x] = 0;}
                        if (cumulative_extracted_pd[y][z][x] < 0) {cumulative_extracted_pd[y][z][x] = 0;}
                        D_ed[y][z][x] += M_ext_pd[y][z][x];   /* source into the advected scalar */
                        // Peridotite melt enriches the local basalt fraction. The
                        // source is injected into the h POPULATIONS in the collision
                        // (h[a] += w[a]*M_ext), NOT by incrementing the moment here:
                        // a moment increment enters h_eq and the BGK relaxation then
                        // amplifies it by 1/tau_h (~2x at tau_h ~ 0.5), which inflated
                        // L by ~26 points in the first 5 My of the smF smoke.
                        // Adding TTG to rho_r — insitu deposit branch
                        if (strcmp(melt_position,"insitu") == 0) {
                            rho_r[y][z][x] += M_ext_bs[y][z][x];
                            rho_b[y][z][x] -= M_ext_bs[y][z][x];
                            frac_r     [x]  = rho_r[y][z][x]/rho[x];
                            frac_b     [x]  = rho_b[y][z][x]/rho[x];
                        }
                        // Surface deposit as a COMPOSITION CHANGE at z = zu.
                        // Each cell carries rho = rho_r + rho_b ≈ 1
                        // (incompressibility), so we do not transport mass
                        // between cells — we just *convert* basalt (rho_b) to
                        // TTG (rho_r) at the surface. The depth cell where the
                        // melt was generated keeps its full rho; the extraction
                        // is already bookkept in cumulative_extracted_bs. The
                        // f_b/f_r distributions at zu are updated alongside so
                        // the conversion survives the next strmac rho=sum(f)
                        // recomputation. The amount converted is capped by the
                        // basalt available at zu (take ≤ rho_b[y][zu][x]) so
                        // rho_b never goes negative; once the surface saturates
                        // (rho_b[zu] → 0) further melt is silently dropped.
                        // Only fires on the rank that owns the upper boundary.
                        // Deposit at z = zu+1 (one cell below the wall) to keep
                        // the wall BC's bounce-back from amplifying the deposit.
                        // (Even at zu+1 mass eventually streams up to zu and the
                        // BC amplifies it, but more slowly. Long-term, a proper
                        // half-way bounce-back BC fix would be needed in strmac.)
                        if (strcmp(melt_position,"surface") == 0 && boundary_u) {
                            /* Deposit walks DOWN from zu+1 until it finds basalt to
                               convert (crustal thickening) instead of silently
                               dropping melt on a saturated surface. Keeps the red
                               budget equal to the bookkept extraction. */
                            FLOAT rem = M_ext_bs[y][z][x];
                            int   zd_dep;
                            for (zd_dep = zu + 1; zd_dep <= zu + 40 && rem > 0; zd_dep++) {
                                FLOAT take = rem;
                                if (take > rho_b[y][zd_dep][x]) take = rho_b[y][zd_dep][x];
                                if (take < 0)                    take = 0;
                                rho_b[y][zd_dep][x] -= take;
                                rho_r[y][zd_dep][x] += take;
                                for (int _aa = 0; _aa < na; _aa++) {
                                    f_b[y][zd_dep][_aa][x] -= take * w[_aa];
                                    f_r[y][zd_dep][_aa][x] += take * w[_aa];
                                }
                                rem -= take;
                            }
                        }
                    }
                }

//---------------------------------------------------------------
//              Calculate total base viscosity from red/blue fractions
//---------------------------------------------------------------

                if(strcmp(temperature,"reynolds") == 0 && gamma_T > 1.e-10)  {                       // Reynolds    temperature dependence
                    for(x=startx;x<endx;x+=incx) {
                        nu_linear[x]   = (nu_f_r[z] * frac_r[x] + nu_f_b[z] * frac_b[x]) * EXP(-gamma_T*eta[y][z][x]);
                    }
                }
                else if (strcmp(temperature,"arrhenius") == 0)               {                       // Arrhenius   temperature dependence
                    // eta_p is the physical reference viscosity (Pa.s) — used to calculate nu_f_b and nu_f_r (set in Python).
                    // nu_f_b[z] already incorporates the Delta_eta depth profile (set in Python).
                    // arg is a dimensionless multiplier: exp(0)=1 at T_ref=ref_K, P=0.
                    gas_C = 8.314;    //  Gas constant              units: J/(mol*K)
                    ref_K = 1600;     //  Reference temperature     units: K
                    for(x=startx;x<endx;x+=incx) {
                        T        [x]      = eta[y][z][x]/DR_2 + T0_p + dT_p/2;
                        T        [x]      = max(T0_p, min(T0_p + dT_p, T[x]));
                        T_K      [x]      = T[x] + 273.15;
                        // P[x] set by archean2 melt block only when melt_extraction==true;
                        // recompute here so Arrhenius works regardless of melt_extraction.
                        P        [x]      = P0[y][z][x] * rho[x];
                        P        [x]      = max(0.0, min(137.0, P[x]));
                        double exponent   = (act_E[z] + P[x]*1e9 * act_V[z])/(gas_C * T_K[x])
                                          - act_E[z]/(ref_K * gas_C);
                        exponent          = min(exponent, (FLOAT)40.0);           // prevent EXP overflow -> NaN
                        arg      [x]      = exp(exponent);                        // dimensionless multiplier
                        nu_linear[x]      = (nu_f_r[z] * frac_r[x] + nu_f_b[z] * frac_b[x]) * arg[x];
                    }
                }
                else {                                                                               // Unspecified temperature dependence
                    for(x=startx;x<endx;x+=incx) {
                        nu_linear [x]  = nu_f_r[z] * frac_r[x] + nu_f_b[z] * frac_b[x];
                    }
                }
//---------------------------------------------------------------
//              Equilibrium and non-equilibrium distributions
//---------------------------------------------------------------

                if(D==2) {
                   for(x=startx;x<endx;x+=incx) {
                       u2[x]    = 1.5 * (u[y][z][0][x] * u[y][z][0][x] + u[y][z][1][x] * u[y][z][1][x]);
                       cu[0][x] = -w[0] * u2[x];
                   }
                   for(a=1;a<na;a++) {
                       for(x=startx;x<endx;x+=incx) {
                           cu[a][x] = 3 *(c[a][0] * u[y][z][0][x] + c[a][1] * u[y][z][1][x]);
                           cu[a][x] = w[a] * (cu[a][x] + 0.5 * cu[a][x] * cu[a][x] - u2[x]);
                       }
                   }
                }
                else if(D==3) {
                   for(x=startx;x<endx;x+=incx) {
                       u2[x] = 1.5 * (u[y][z][0][x] * u[y][z][0][x] + u[y][z][1][x] * u[y][z][1][x] + u[y][z][2][x] * u[y][z][2][x]);
                       cu[0][x] = -w[0] *  u2[x];
                   }
                   for(a=1;a<na;a++) {
                       for(x=startx;x<endx;x+=incx) {
                           cu[a][x] = 3 * (c[a][0] * u[y][z][0][x] + c[a][1] * u[y][z][1][x] + c[a][2] * u[y][z][2][x]);
                           cu[a][x] = w[a] * (cu[a][x] + 0.5 * cu[a][x] * cu[a][x] - u2[x]);
                       }
                   }
                }
                for(a=0;a<na;a++) {
                    for(x=startx;x<endx;x+=incx) {
                        feq_r  [a][x]  = rho_r[y][z][x] * (C_r[a] + cu[a][x]);
                        feq_b  [a][x]  = rho_b[y][z][x] * (C_b[a] + cu[a][x]);
                        // Shifted distribution: equilibrium is feq̃ = feq - rho_init*C[a]
                        // so that f̃_neq = f̃ - feq̃ = f_neq (collision unchanged)
                        if (shifted_distribution) {
                            feq_r[a][x] -= rho_r_init[y][z][x] * C_r[a];
                            feq_b[a][x] -= rho_b_init[y][z][x] * C_b[a];
                        }
                        g_eq   [a][x]  = rho[x] * eta[y][z][x]   * (w[a] + cu[a][x]);
                        f_r_neq[a][x]  = f_r[y][z][a][x] - feq_r[a][x];
                        f_b_neq[a][x]  = f_b[y][z][a][x] - feq_b[a][x];
                        g_neq  [a][x]  = g  [y][z][a][x] - (g_eq[a][x] + Delta_g[a]);
                        if (A2_HEQ && strcmp(rheology,"archean2") == 0) {
                            h_eq   [a][x]    = L_ed[y][z][x] * (w[a] + cu[a][x]);
                            h_neq  [a][x]    = h   [y][z][a][x] - h_eq[a][x];

                        }
                    }
                }
//---------------------------------------------------------------
//              Calculate nu total: allow for nu_plastic if needed
//---------------------------------------------------------------

                if(strcmp(rheology,"linear") == 0)  {  // nu linear if linear rheology
                    for(x=startx;x<endx;x+=incx) {
                        nu_total[x]        = max(nu_linear[x],nu_min);
                        // Physical viscosity clamp: min = eta_min_p, max = eta_max_p
                        if (nu_total[x] < nu_min_cap) nu_total[x] = nu_min_cap;
                        if (nu_total[x] > nu_max_cap) nu_total[x] = nu_max_cap;
                        nu_f_last[y][z][x] = nu_total[x];
                        tau_fi  [x]        = nu_total[x] * 3 + 0.5;
                        tau_gi  [z]        = tau_g[y][z][x];
                    }
                }
                else {
//                  Calculate S_ij and (strain rate)^2
                    for(i=0;i<D;i++) {
                        for(j=0;j<D;j++) {
                            for(x=startx;x<endx;x+=incx) {
                                S[i][j][x] = 0.0;
                            }
                        }
                    }
                    for(i=0;i<D;i++) {
                        for(j=0;j<D;j++) {
                            for(a=1;a<na;a++) {
                                itmp = c[a][i] * c[a][j];
                                if(itmp != 0) {
                                    for(x=startx;x<endx;x+=incx) {
                                        S[i][j][x] += ((f_r_neq[a][x] + f_b_neq[a][x]) * itmp);
                                    }
                                }
                            }
                        }
                    }
                    for(x=startx;x<endx;x+=incx) {
                        strain_rate[x] = 0.0;
                        S_factor   [x] = 1.5/(rho[x] * tau_f_last[y][z][x]);
                    }
                    for(i=0;i<D;i++) {
                        for(j=0;j<D;j++) {
                            for(x=startx;x<endx;x+=incx) {
                                S    [i][j][x] *= S_factor[x];
                                strain_rate[x] += (S[i][j][x] * S[i][j][x]);
                            }
                        }
                    }
//                  Calculate shear stress
                    for(x=startx;x<endx;x+=incx) {
                        strain_rate[x] = SQRT(0.5 * strain_rate[x]);
                        Tau        [x] = rho[x] * nu_f_last[y][z][x] * strain_rate[x];
                    }
//                  Plastic rheology
                    if (strcmp(nu_formula,"Tau") == 0) {            // nu_plastic formula Tau_y/Tau
                        for(x=startx;x<endx;x+=incx) {
                            nu_plastic[x] = (nu_f_r[z] * frac_r[x] + nu_f_b[z] * frac_b[x]) * pow((double)Tau_y[z]/(Tau[x]+1.e-5*(double)Tau_y[z]),(double)N[z]);
                        }
                    }
                    if (strcmp(nu_formula,"strain_rate") == 0) {    // nu_plastic formula Tau_y/(2*strain_rate)
                        for(x=startx;x<endx;x+=incx) {
                            nu_plastic[x] = nu_star + pow((double)Tau_y[z]/((2 * strain_rate[x]) + (1.e-5*(double)Tau_y[z])),(double)N[z]);
                        }
                    }
//                  Calculate total viscosity due to liner & plastic terms, and relaxation time

                    for(x=startx;x<endx;x+=incx) {
                        nu_total        [x] = max(1/(1/nu_linear[x] + 1/nu_plastic[x]),nu_min);
                        // Physical viscosity clamp: min = eta_min_p, max = eta_max_p
                        if (nu_total[x] < nu_min_cap) nu_total[x] = nu_min_cap;
                        if (nu_total[x] > nu_max_cap) nu_total[x] = nu_max_cap;
                        nu_f_last [y][z][x] = nu_total[x];
                        tau_fi          [x] = nu_total[x] * 3 + 0.5;
                        tau_f_last[y][z][x] = tau_fi[x];
                        tau_gi          [z] = tau_g[y][z][x];
                    }
                }

//---------------------------------------------------------------
//              Add collision terms
//---------------------------------------------------------------
                for(a=0;a<na;a++) {
                    for(x=startx;x<endx;x+=incx) {
                        f_r[y][z][a][x] -= (f_r_neq[a][x]/tau_fi[x]);
                        f_b[y][z][a][x] -= (f_b_neq[a][x]/tau_fi[x]);
                        g  [y][z][a][x] -= (g_neq  [a][x]/tau_gi[z]);
                        if (A2_HCOLL && strcmp(rheology,"archean2") == 0) {
                            h[y][z][a][x]   -= (h_neq[a][x]/tau_h[z]);
                            /* depletion transport is handled by the TVD pass;
                               melt sources enter D_ed/D_eb directly in the melt block */
                            if (melt_extraction) {
                                h[y][z][a][x] += w[a] * M_ext_pd[y][z][x];   /* exact L source */
                            }
                        }
                    }
                }
//---------------------------------------------------------------
                if(collision) {
                   for(d=0;d<D;d++) {
                       for(x=startx;x<endx;x+=incx) {
                         F[d][x] = 0;
                       }
                       for (a=1;a<na;a++) {
                            FLOAT b, w = 0.298, W = 2 + 4*w;
                            if(a < 5) {b = 1/W;}
                            else {b = w/W;}
                            yy_ya = y_ya[a][y];
                            zz_za = z_za[a][z];
                            for(x=startx;x<endx;x+=incx) {
                                xx_xa    = x_xa[a][x];
                                F[d][x] -= b * c[a][d] * (rho_r[yy_ya][zz_za][xx_xa] - rho_b[yy_ya][zz_za][xx_xa]);
                            }
                       }
                   }
                   if(D==2) {
                      for(x=startx;x<endx;x+=incx) {
                          fm[x] = sqrt(F[0][x] * F[0][x] + F[1][x] * F[1][x]);
                      }
                      for(a=0;a<na;a++){
                          for(x=startx;x<endx;x+=incx) {
                              ttmp[x]          = c[a][0] * F[0][x] + c[a][1] * F[1][x];
                              cos_lambda[a][x] = ttmp[x]/(fm[x] * cm[a] + zero);
                              cos_phi   [a][x] = ttmp[x]/(fm[x] + zero);
                          }
                      }
                   }  
                   else if(D==3) {
                      for(x=startx;x<endx;x+=incx) {
                          fm[x] = sqrt(F[0][x] * F[0][x] + F[1][x] * F[1][x] + F[2][x] * F[2][x]);
                      }
                      for(a=0;a<na;a++){
                          for(x=startx;x<endx;x+=incx) {
                              ttmp[x]   = c[a][0]*F[0][x] + c[a][1]*F[1][x] + c[a][2]*F[2][x];
                              cos_lambda[a][x] = ttmp[x]/(fm[x] * cm[a] + zero);
                              cos_phi   [a][x] = ttmp[x]/(fm[x] + zero);
                          }  
                      }
                   }
//---------------------------------------------------------------
//              Save color gradient
//---------------------------------------------------------------
                   for(d=0;d<D;d++) {
                       for(x=startx;x<endx;x+=incx) {
                           FF[d][y][z][x] = F[d][x];
                       }
                   }
//---------------------------------------------------------------
//                 Collision term 2
//---------------------------------------------------------------
                   if(A > 0) {
                      for(a=0;a<na;a++) {
                           for(x=startx;x<endx;x+=incx) {
                               f_r[y][z][a][x] += A * fm[x]*(w[a] * cos_phi[a][x] * cos_phi[a][x] - B[a]);
                               f_b[y][z][a][x] += A * fm[x]*(w[a] * cos_phi[a][x] * cos_phi[a][x] - B[a]);
                           }
                      }
                   }
//---------------------------------------------------------------
//                 Recoloring step
//
//  Shifted distribution: stored f̃ = f - rho_init*w, so the true
//  total is f_r+f_b = f̃_r+f̃_b + (rho_r_init+rho_b_init)*w.
//  ff_sum must use the unshifted total, and the output must remain
//  in shifted form (subtract rho_init*w after writing back).
//---------------------------------------------------------------
                   for(x=startx;x<endx;x+=incx) {
                       // Recolor wherever any colour gradient exists. The old cutoff
                       // (fm > 0.002) let colour that diffused beyond the interface
                       // escape recapture: at low Mach / high nu_l the colour field
                       // then diffuses at the momentum viscosity and red mass "grows"
                       // out of blue regions. A ~zero threshold keeps the LKR
                       // anti-diffusion active on every mixed cell.
                       if(fm[x]>1.e-9) {
                          for(a=0;a<na;a++) {
                              ff_sum          = (f_r[y][z][a][x] + f_b[y][z][a][x])/rho[x];
                              if (shifted_distribution)
                                  ff_sum += (rho_r_init[y][z][x] + rho_b_init[y][z][x]) * w[a] / rho[x];
                              feq             =  w[a] * rho[x];
                              tempp           = (interface_thickness * rho_b[y][z][x] * rho_r[y][z][x])/(rho[x] * rho[x]) * cos_lambda[a][x] * feq;
                              f_r[y][z][a][x] =  rho_r[y][z][x] * ff_sum + tempp;
                              f_b[y][z][a][x] =  rho_b[y][z][x] * ff_sum - tempp;
                              if (shifted_distribution) {
                                  f_r[y][z][a][x] -= rho_r_init[y][z][x] * w[a];
                                  f_b[y][z][a][x] -= rho_b_init[y][z][x] * w[a];
                              }
                          }
                       }
                   }
                }
//---------------------------------------------------------------
//              LKR recoloring for the advected depletion scalars:
//              treat (D, 1-D) as an immiscible pair on a unit carrier
//              and segregate along the depletion gradient with the
//              same interface_thickness as the phases. This keeps
//              depletion fronts sharp so material parcels do not lose
//              their depletion value by numerical diffusion, while
//              TAU_HD keeps the bulk diffusion near zero. Away from
//              fronts (|grad D| ~ 0) the BGK result is kept.
//---------------------------------------------------------------
                if (RECOLOR_HD && collision && strcmp(rheology,"archean2") == 0 && melt_extraction) {
                    for (int sfield = 0; sfield < 2; sfield++) {
                        for(d=0;d<D;d++) {
                            for(x=startx;x<endx;x+=incx) F[d][x] = 0;
                            for (a=1;a<na;a++) {
                                FLOAT bb, ww = 0.298, WW = 2 + 4*ww;
                                if(a < 5) {bb = 1/WW;} else {bb = ww/WW;}
                                yy_ya = y_ya[a][y];
                                zz_za = z_za[a][z];
                                for(x=startx;x<endx;x+=incx) {
                                    xx_xa    = x_xa[a][x];
                                    F[d][x] -= bb * c[a][d] * (sfield==0 ? D_ed[yy_ya][zz_za][xx_xa]
                                                                         : D_eb[yy_ya][zz_za][xx_xa]);
                                }
                            }
                        }
                        if (D==2) {
                            for(x=startx;x<endx;x+=incx) fm[x] = sqrt(F[0][x]*F[0][x] + F[1][x]*F[1][x]);
                        } else {
                            for(x=startx;x<endx;x+=incx) fm[x] = sqrt(F[0][x]*F[0][x] + F[1][x]*F[1][x] + F[2][x]*F[2][x]);
                        }
                        for(a=0;a<na;a++){
                            for(x=startx;x<endx;x+=incx) {
                                ttmp[x] = c[a][0]*F[0][x] + c[a][1]*F[1][x];
                                if (D==3) ttmp[x] += c[a][2]*F[2][x];
                                cos_lambda[a][x] = ttmp[x]/(fm[x] * cm[a] + zero);
                            }
                        }
                        // Additive anti-diffusion only: the TAU_HD-relaxed populations
                        // are KEPT (overwriting them with carrier equilibrium imposed a
                        // tau=1 diffusivity and homogenized the field) and a small
                        // counter-gradient segregation flux is added. BETA_HD is sized
                        // to balance the TAU_HD diffusivity (0.007), not the phase
                        // recolorer's 1/6; the phase value 0.75 over-sharpens ~25x and
                        // produces lattice-scale herringbone oscillations.
                        for(x=startx;x<endx;x+=incx) {
                            if (fm[x] > 1.e-9) {
                                FLOAT Dm   = (sfield==0) ? D_ed[y][z][x]     : D_eb[y][z][x];
                                FLOAT Dsrc = (sfield==0) ? M_ext_pd[y][z][x] : M_ext_bs[y][z][x];
                                FLOAT Dtot = Dm + Dsrc;
                                if (Dtot < 0) Dtot = 0;
                                for(a=0;a<na;a++) {
                                    FLOAT segr = BETA_HD * Dtot * (1 - Dtot) * cos_lambda[a][x] * w[a];
                                    if (sfield==0) hd[y][z][a][x] += segr;
                                    else           hb[y][z][a][x] += segr;
                                }
                            }
                        }
                    }
                }
//---------------------------------------------------------------
//              Add Boussinesq bouyancy force (assumes that d = D-1 is the z-axis)
//---------------------------------------------------------------
                for(x=startx;x<endx;x+=incx) {
                    gi[x]   = beta_g_factor * rho[x] * beta_g * eta[y][z][x]; // Boussinesq thermal bouyancy term
                    g_r[x]  = frac_r[x] * gi[x];
                    g_b[x]  = frac_b[x] * gi[x];
                    if (A2_BUOY && strcmp(rheology,"archean2") == 0 && melt_extraction ) {
                        P_garnet_in[x]    = 0.014 * T[x]  - 5.4;
                        P_plag_out [x]    = 0.024 * T[x]  + 4;
                        // Transition lines are in kbar: compare P*10 (P is GPa).
                        // Above plagioclase-out the basalt is full eclogite and the
                        // coefficient SATURATES at 0.16 (Jain et al. 2019; deep basalt
                        // stays denser than olivine to the lowermost mantle).
                        if      (P[x]*10 <  P_garnet_in[x]) {coef_phase_bs[x] = 0.0;}
                        else if (P[x]*10 <= P_plag_out[x])  {coef_phase_bs[x] = 0.16 * (P[x]*10 - P_garnet_in[x]) /(P_plag_out[x] - P_garnet_in[x]);}
                        else                                {coef_phase_bs[x] = 0.16;}
                        if (M_ext_bs[y][z][x] > 0) {coef_res_bs[x]   =  0.4  * M_ext_bs[y][z][x];}
                        else {coef_res_bs[x] = 0;}
                        if (M_ext_pd[y][z][x] > 0) {coef_res_pd[x]   = -0.04 * M_ext_pd[y][z][x];}
                        else {coef_res_pd[x] = 0;}
                        coef_phase_pd[x] = 0.0;
                        // Change of Peridotite density due to melt extraction
                        rho_solid_pd[y][z][x]  = rho_b[y][z][x] * (1- L_den[x]) * (1 - M_ext_pd[y][z][x]) * rho2_p * (1 + coef_res_pd[x] + coef_phase_pd[x]);
                        // Extracted basalt from Peridotite
                        rho_ext_pd  [y][z][x]  = rho_b[y][z][x] * (1- L_den[x]) *      M_ext_pd[y][z][x]  * rho3_p;
                        // Change of Basalt density due to melt extraction
                        rho_solid_bs[y][z][x]  = rho_b[y][z][x] *     L_den[x]  * (1 - M_ext_bs[y][z][x]) * rho3_p * (1 + coef_res_bs[x] + coef_phase_bs[x]);
                        rho_solid   [y][z][x]  = rho_solid_pd[y][z][x] + rho_ext_pd[y][z][x] + rho_solid_bs[y][z][x];
                        rho_eff     [y][z][x]  = rho_solid[y][z][x];
                        // Divide by the ORIGINAL peridotite reference (rho2_p, constant),
                        // not by the per-cell peridotite mass rho_solid_pd which can
                        // collapse to ~0 in basalt-rich cells (L_den -> 1) and produce
                        // a runaway delta_rho. The depletion physics ((1 + coef_res_pd))
                        // is preserved in the rho_solid_pd / rho_solid_bs build above,
                        // so depleted peridotite still gets a lighter rho_solid and a
                        // negative delta_rho (rises). The (rho_solid_pd > 0) guard is
                        // no longer needed because rho2_p is a constant.
                        rho_comp [x]  = (rho_r[y][z][x] * rho1_p + rho_solid[y][z][x]) / rho2_p;
                        delta_rho[x]  = rho_comp[x] - rho[x];
                    }
                    else {
                        rho_comp [x] = rho_r[y][z][x] * rho1_p/rho2_p \
                                     + rho_b[y][z][x] * ((1 - L_ed[y][z][x]) * rho2_p + L_ed[y][z][x] * rho3_p) / rho2_p;
                        delta_rho[x] = rho_comp[x] - rho[x];
                    }
                    gi [x]      += g_comp_factor * 3 * delta_rho[x] * g_l;  // Add compositional buoyancy to gi
                    g_r[x]       = frac_r[x] * gi[x];
                    g_b[x]       = frac_b[x] * gi[x];
                }
                if(false&&((rank==0&&z<startz+1)||(rank==31&&z>endz-2))) {
                   float eta_min=1.e10,eta_max=-1.e10;
                   float g_r_min=1.e10,g_r_max=-1.e10;
                   float g_b_min=1.e10,g_b_max=-1.e10;
                   for(x=startx;x<endx;x++) {
                       eta_min = min(eta_min,eta[0][z][x]);
                       eta_max = max(eta_max,eta[0][z][x]);
                       g_r_min = min(g_r_min,g_r[x]);
                       g_r_max = max(g_r_max,g_r[x]);
                       g_b_min = min(g_b_min,g_b[x]);
                       g_b_max = max(g_b_max,g_b[x]);
                    }
                    printf("------------\n");
                    printf("rank = %d, z = %d, beta_g = %8.6g, g_r E [%8.3g,%8.3g], g_b E [%8.3g,%8.3g], eta E [%8.3g,%8.3g]\n",rank,z,beta_g,g_r_min,g_r_max,g_b_min,g_b_max,eta_min,eta_max);
                }
                for(a=1;a<na;a++) {
                    if(c[a][D-1]>0) {
                        for(x=startx;x<endx;x+=incx) {
                            f_r[y][z][a][x] +=  (w[a] * g_r[x]);
                            f_b[y][z][a][x] +=  (w[a] * g_b[x]);
                        }
                    }
                    else if(c[a][D-1]<0) {
                        for(x=startx;x<endx;x+=incx) {
                            f_r[y][z][a][x] -=  (w[a] * g_r[x]);
                            f_b[y][z][a][x] -=  (w[a] * g_b[x]);
    		        }
                    }
    	        }
//---------------------------------------------------------------
            }
        }
    }
