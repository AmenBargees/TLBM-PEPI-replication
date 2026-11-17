//-------------------------------------------------------------------------
//  C code for TLBM_fast with overlapping option of calcuations and
//  communications (overlap bool)
//  Allows for interior or edges to be calculated according to whether
//  interior = true or false
//-------------------------------------------------------------------------

// I need to speak to Henry (T) to understand whether pragma block loop factor can
// get better performance, which it seems is should. Namely, from the description,
// pragma block_loop factor directives should be able to be applied to the loops
// where the loop levels for the TLBM_fast structure would be 1=x,2=a,3=z,4=y.
// So by changing the block loop factor well for each loop level, one should be
// able to make cache blocing more optimal on a given chip and to suite it's L3
// cache size. However, when I use the pragma block_loop directive, it makes no
// difference to performance. And if I use large block loop factors where the memory
// required for f,f_last,g and g_last on a domain would exceed the L3 cache size,
// the speed remains unchanged! Maybe my data structure and loops are already optimal???
//  Maybe by putting a and x into one loop will speed up calculations and streaming???

// The performance tests below were on an Intel Max 9480 HBM chip with 112.5MB
// L3 cache, and used cache block loop factors of 256(1),19(2),128(3),128(4)
//
// 3D 512^3           (nxi_min=nxi_max=256)    : 8.8  MFlups (16  cores)
//                    (nxi_min=nxi_max=512)    : 6.8  MFlups (112 cores)
// 2D small 512x63    (nxi_min=nxi_max=256)    : 27.5 MFlups (16  cores), z-domain=8
// 2D small 512x63    (nxi_min=128=nxi_max=128): 5.5  MFlups (112 cores), z-domain=2
// 2D large 8192x2048 (nxi_min=nxi_max=512)    : 13.7 MFlups (112 cores)
// MFLups for this large 2D run dropped to ~11 Mflups after 8 mins and was still dropping???

// The pragma block_loop directives make no differnce no matter what values were used
// for each level. As such, they are commented out for now. Need to consult Henry T!

//#pragma block_loop factor(256) level(1) 
//#pragma block_loop factor (19) level(2)
//#pragma block_loop factor(128) level(3)
//#pragma block_loop factor(128) level(4)

    int incx, startx, endx, starty, endy, startz, endz;
    int x, y, z, a, d, i, j, xx_xa, yy_ya, zz_za, itmp;
    FLOAT tmp;
    bool edgey, edgez, x_loop, boundary_bb, boundary_ff, boundary_uu, boundary_dd;

    FLOAT f_eq[na][nx], f_neq[na][nx], g_neq[na][nx];
    FLOAT drho[nx], cu[na][nx], u2[nx];
    FLOAT strain_rate[nx], S[D][D][nx], S_factor[nx], Tau[nx];
    FLOAT nu_linear[nx], nu_plastic[nx], tau_f[nx], gi[nx];


//  TTG variables
    FLOAT h_eq[na][nx], h_neq[na][nx];
    FLOAT P[ny][nz][nx], PKbar[ny][nz][nx], T[ny][nz][nx], arg[ny][nz][nx] ;
    FLOAT eqn1_sol_bs, eqn2_sol_bs;
    FLOAT T_solidus_bs, T_liquidus_bs, T_solidus_pd, T_liquidus_pd;
    FLOAT M_pd[y][z][x], M_bs[y][z][x];
    FLOAT M0_bs, M0_pd;
    FLOAT L_den[ny][nz][nx];
    FLOAT P_garnet_in, P_plag_out;
    FLOAT coef_phase_bs, coef_res_bs, coef_phase_pd, coef_res_pd;
    FLOAT rho_solid_pd[ny][nz][nx], rho_solid_bs[ny][nz][nx];
    bool T1_cond, P1_cond, mask_lpttg, mask_mpttg, LP_TTG, MP_TTG;
    bool T2_cond, P2_cond, mask_hpttg, HP_TTG;

//---------------------------------------------------------------
//  Set the x loop specifications as required for both the
//  case of no overlapping and overlapping which calculates
//  either on interior regions or on edges of domains
//---------------------------------------------------------------

//  DBG print throttle: only print every DBG_INTERVAL kernel calls.
#define DBG_INTERVAL 100

//  DBG: step counter and single-cell probe location (rank 0, middle of domain)
#if 0
    static long long _dbg_call_fast = -1;
    _dbg_call_fast++;
    int _pry = (my > 1) ? my/2 : 0;
    int _prz = mz/2;
    int _prx = mx/2;
#endif

    startx = xedge  ; endx = mx-xedge;
    starty = yedge  ; endy = my-yedge;
    startz = zedge  ; endz = mz-zedge;

//  Re-set loop limits for interior regions of physical domains for the
//  case of overlapping, ie. without physical domain edges on this core

    if(overlap&&interior) {
        startx += xedge; endx -= xedge;
        starty += yedge; endy -= yedge;
        startz += zedge; endz -= zedge;
    }

//---------------------------------------------------------------
//  Start y and z loops and edge bool for y edges or z edges
//---------------------------------------------------------------

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

	    // For overlaping/calculate edges case: calculate only at x
	    // edges if not at a y or z edge, otherwise calculate all x
            if(overlap&&!interior) {
                if(!edgez&&!edgey) { incx = endx-startx-1; }
	        else               { incx = 1;             }
//              If not at a yedge or zedge & there is 1 x-domain, no x loop is needed
                if(xedge==0&&(!edgez&&!edgey)) { x_loop = false; }
                else                           { x_loop = true;  }
            }
	    else                   { incx = 1;   x_loop = true;  }

//          Do the calculations if required. Not required for
//          edge calculations for overlapping case if not at
//          a y edge or a z edge if there is just one x domain

            if(x_loop) {

//---------------------------------------------------------------
//              Streaming
//---------------------------------------------------------------

                for(a=0;a<na;a++) {
                    yy_ya = y_ya[a][y];
                    zz_za = z_za[a][z];
                    for(x=startx;x<endx;x+=incx) {
                        xx_xa         = x_xa[a][x];
                        f[y][z][a][x] = f_last[yy_ya][zz_za][a][xx_xa];
                        g[y][z][a][x] = g_last[yy_ya][zz_za][a][xx_xa];
                        if (strcmp(rheology,"archean2") == 0) {
                            h[y][z][a][x] = h_last[yy_ya][zz_za][a][xx_xa];}
                    }
                }
#if 0
                // DBG probe: after streaming
                if (rank == 0 && y == _pry && z == _prz && (_dbg_call_fast % DBG_INTERVAL == 0)) {
                    double _sum_f = 0, _sum_g = 0;
                    for (int _aa = 0; _aa < na; _aa++) {
                        _sum_f += f[y][z][_aa][_prx];
                        _sum_g += g[y][z][_aa][_prx];
                    }
                    printf("DBG fast   call=%lld stage=STREAM   y=%d z=%d x=%d sum_f=%.9e sum_g=%.9e\n",
                           _dbg_call_fast, y, z, _prx, _sum_f, _sum_g);
                }
#endif

//---------------------------------------------------------------
//              Boundary conditions
//
//  pre_streaming_BC=true  : read incoming wall distributions from
//    f_last/g_last/h_last (the pre-streaming values at the wall node),
//    matching the RKTLBM strmac convention — correct halfway bounce-back.
//  pre_streaming_BC=false : read from f/g/h (post-streaming from the
//    neighbouring cell) — legacy behaviour.
//  Applies only to the pure no-slip loop (a < naX).  The halfway loop
//  (naX <= a < nbX) uses post-streaming values in both conventions.
//---------------------------------------------------------------
#define BC_F(y,z,a,x) (pre_streaming_BC ? f_last[y][z][a][x] : f[y][z][a][x])
#define BC_G(y,z,a,x) (pre_streaming_BC ? g_last[y][z][a][x] : g[y][z][a][x])
#define BC_H(y,z,a,x) (pre_streaming_BC ? h_last[y][z][a][x] : h[y][z][a][x])

                // Left BC
                if(boundary_l) {
                    for(a=0;a<naR;a++) {
                        f[y][z][aR[a]][xl] = BC_F(y,z,aiR[a],xl);
                        g[y][z][aR[a]][xl] = BC_G(y,z,aiR[a],xl);
                        if (strcmp(rheology,"archean2") == 0){
                            h[y][z][aR[a]][xl] = BC_H(y,z,aiR[a],xl);}
                    }
                    for(a=naR;a<nbR;a++) {
                        tmp = 0.5*(f[y][z][aR[a]][xl] + f[y][z][aiR[a]][xl]);
                        f[y][z][aR [a]][xl] = tmp;
                        f[y][z][aiR[a]][xl] = tmp;
                        tmp = 0.5*(g[y][z][aR[a]][xl] + g[y][z][aiR[a]][xl]);
                        g[y][z][aR [a]][xl] = tmp;
                        g[y][z][aiR[a]][xl] = tmp;
                        if (strcmp(rheology,"archean2") == 0) {
                            tmp = 0.5*(h[y][z][aR[a]][xl] + h[y][z][aiR[a]][xl]);
                            h[y][z][aR [a]][xl] = tmp;
                            h[y][z][aiR[a]][xl] = tmp;}
                   }
                }
                // Right BC
                if(boundary_r) {
                    for(a=0;a<naL;a++) {
                        f[y][z][aL[a]][xr] = BC_F(y,z,aiL[a],xr);
                        g[y][z][aL[a]][xr] = BC_G(y,z,aiL[a],xr);
                        if (strcmp(rheology,"archean2") == 0){
                            h[y][z][aL[a]][xr] = BC_H(y,z,aiL[a],xr);}
                    }
                    for(a=naL;a<nbL;a++) {
                        tmp = 0.5*(f[y][z][aL[a]][xr] + f[y][z][aiL[a]][xr]);
                        f[y][z][aL [a]][xr] = tmp;
                        f[y][z][aiL[a]][xr] = tmp;
                        tmp = 0.5*(g[y][z][aL[a]][xr] + g[y][z][aiL[a]][xr]);
                        g[y][z][aL [a]][xr] = tmp;
                        g[y][z][aiL[a]][xr] = tmp;
                        if (strcmp(rheology,"archean2") == 0) {
                            tmp = 0.5*(h[y][z][aL[a]][xr] + h[y][z][aiL[a]][xr]);
                            h[y][z][aL [a]][xr] = tmp;
                            h[y][z][aiL[a]][xr] = tmp;}
                    }
                 }
                // Upper BC
                if(boundary_uu) {
                    for(a=0;a<naD;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][aD[a]][x] = BC_F(y,z,aiD[a],x);
                            g[y][z][aD[a]][x] = BC_G(y,z,aiD[a],x);
                            if (strcmp(rheology,"archean2") == 0) {
                                h[y][z][aD[a]][x] = BC_H(y,z,aiD[a],x);}
                        }
                    }
                    for(a=naD;a<nbD;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            tmp = 0.5*(f[y][z][aD[a]][x] + f[y][z][aiD[a]][x]);
                            f[y][z][aD [a]][x] = tmp;
                            f[y][z][aiD[a]][x] = tmp;
                            tmp = 0.5*(g[y][z][aD[a]][x] + g[y][z][aiD[a]][x]);
                            g[y][z][aD [a]][x] = tmp;
                            g[y][z][aiD[a]][x] = tmp;
                            if (strcmp(rheology,"archean2") == 0) {
                                tmp = 0.5*(h[y][z][aD[a]][x] + h[y][z][aiD[a]][x]);
                                h[y][z][aD [a]][x] = tmp;
                                h[y][z][aiD[a]][x] = tmp;}
                        }
                    }
                }
                // Lower BC
                if(boundary_dd) {
                    for(a=0;a<naU;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][aU[a]][x] = BC_F(y,z,aiU[a],x);
                            g[y][z][aU[a]][x] = BC_G(y,z,aiU[a],x);
                            if (strcmp(rheology,"archean2") == 0) {
                                h[y][z][aU[a]][x] = BC_H(y,z,aiU[a],x);}
                        }
                    }
                    for(a=naU;a<nbU;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            tmp = 0.5*(f[y][z][aU[a]][x] + f[y][z][aiU[a]][x]);
                            f[y][z][aU [a]][x] = tmp;
                            f[y][z][aiU[a]][x] = tmp;
                            tmp = 0.5*(g[y][z][aU[a]][x] + g[y][z][aiU[a]][x]);
                            g[y][z][aU [a]][x] = tmp;
                            g[y][z][aiU[a]][x] = tmp;
                            if (strcmp(rheology,"archean2") == 0) {
                                tmp = 0.5*(h[y][z][aU[a]][x] + h[y][z][aiU[a]][x]);
                                h[y][z][aU [a]][x] = tmp;
                                h[y][z][aiU[a]][x] = tmp;}
                        }
                    }
                }
                // Backward BC
                if(boundary_bb) {
                    for(a=0;a<naF;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][aF[a]][x] = BC_F(y,z,aiF[a],x);
                            g[y][z][aF[a]][x] = BC_G(y,z,aiF[a],x);
                            if (strcmp(rheology,"archean2") == 0) {
                                h[y][z][aF[a]][x] = BC_H(y,z,aiF[a],x);}
                        }
                    }
                    for(a=naF;a<nbF;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            tmp = 0.5*(f[y][z][aF[a]][x] + f[y][z][aiF[a]][x]);
                            f[y][z][aF [a]][x] = tmp;
                            f[y][z][aiF[a]][x] = tmp;
                            tmp = 0.5*(g[y][z][aF[a]][x] + g[y][z][aiF[a]][x]);
                            g[y][z][aF [a]][x] = tmp;
                            g[y][z][aiF[a]][x] = tmp;
                            if (strcmp(rheology,"archean2") == 0) {
                                tmp = 0.5*(h[y][z][aF[a]][x] + h[y][z][aiF[a]][x]);
                                h[y][z][aF [a]][x] = tmp;
                                h[y][z][aiF[a]][x] = tmp;}
                        }
                    }
                }
                // Forward BC
                if(boundary_ff) {
                    for(a=0;a<naB;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][aB[a]][x] = BC_F(y,z,aiB[a],x);
                            g[y][z][aB[a]][x] = BC_G(y,z,aiB[a],x);
                            if (strcmp(rheology,"archean2") == 0) {
                                h[y][z][aB[a]][x] = BC_H(y,z,aiB[a],x);}
                        }
                    }
                    for(a=naB;a<nbB;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            tmp = 0.5*(f[y][z][aB[a]][x] + f[y][z][aiB[a]][x]);
                            f[y][z][aB [a]][x] = tmp;
                            f[y][z][aiB[a]][x] = tmp;
                            tmp = 0.5*(g[y][z][aB[a]][x] + g[y][z][aiB[a]][x]);
                            g[y][z][aB [a]][x] = tmp;
                            g[y][z][aiB[a]][x] = tmp;
                            if (strcmp(rheology,"archean2") == 0) {
                                tmp = 0.5*(h[y][z][aB[a]][x] + h[y][z][aiB[a]][x]);
                                h[y][z][aB [a]][x] = tmp;
                                h[y][z][aiB[a]][x] = tmp;}
                        }
                    }
                }

#undef BC_F
#undef BC_G
#undef BC_H

//---------------------------------------------------------------
//              Macroscopic properties. Backards a loop to
//              optimize precision by summing small values first
//---------------------------------------------------------------

#if 0
                // DBG probe: after all BCs, before macroscopic
                if (rank == 0 && y == _pry && z == _prz && (_dbg_call_fast % DBG_INTERVAL == 0)) {
                    double _sum_f = 0, _sum_g = 0;
                    for (int _aa = 0; _aa < na; _aa++) {
                        _sum_f += f[y][z][_aa][_prx];
                        _sum_g += g[y][z][_aa][_prx];
                    }
                    printf("DBG fast   call=%lld stage=BC       y=%d z=%d x=%d sum_f=%.9e sum_g=%.9e\n",
                           _dbg_call_fast, y, z, _prx, _sum_f, _sum_g);
                }
#endif

                for(x=startx;x<endx;x+=incx) {
                    rho[y][z][x] = f[y][z][na-1][x];
                    eta[y][z][x] = g[y][z][na-1][x];
                    if (strcmp(rheology,"archean2") == 0) {
                        L_ed[y][z][x] = h[y][z][na-1][x];}
                }
                for(a=na-2;a>=0;a--) {
                    for(x=startx;x<endx;x+=incx) {
                        rho[y][z][x] += f[y][z][a][x];
                        eta[y][z][x] += g[y][z][a][x];
                        if (strcmp(rheology,"archean2") == 0) {
                            L_ed[y][z][x] += h[y][z][a][x];}
                    }
                }
                if(shifted_distribution) {
                    for(x=startx;x<endx;x+=incx) {
                        drho[x]       = rho[y][z][x];
                        rho[y][z][x] += 1;
                    }
                }
                for(x=startx;x<endx;x+=incx) {
                    eta[y][z][x] /= rho[y][z][x];
                }
                for(d=0;d<D;d++) {
                    for(x=startx;x<endx;x+=incx)     { u[y][z][d][x] = 0; }
                    for(a=na-1;a>0;a--) {
                        for(x=startx;x<endx;x+=incx) { u[y][z][d][x] += (f[y][z][a][x] * c[a][d]); }
                    }
                    for(x=startx;x<endx;x+=incx) {
                        u[y][z][d][x] /= rho[y][z][x];
                    }
                }
#if 0
                // DBG probe: after macroscopic (rho, eta, u)
                if (rank == 0 && y == _pry && z == _prz && (_dbg_call_fast % DBG_INTERVAL == 0)) {
                    printf("DBG fast   call=%lld stage=MACRO    y=%d z=%d x=%d rho=%.9e eta=%.9e u_x=%.9e u_z=%.9e\n",
                           _dbg_call_fast, y, z, _prx,
                           rho[y][z][_prx], eta[y][z][_prx],
                           u[y][z][0][_prx], u[y][z][D-1][_prx]);
                }
#endif

//---------------------------------------------------------------
//              Archean2 melt extraction
//---------------------------------------------------------------
                if (strcmp(rheology,"archean2") == 0) {
                    for(x=startx;x<endx;x+=incx) {
                        T[y][z][x]             = (eta[y][z][x]/DR_2 + T0_p + dT_p/2);
                        P[y][z][x]             =  P0[y][z][x] * rho[y][z][x];
                        eqn1_sol_bs            =  699.85 - (70400/(P[y][z][x]*1000 + 354)) + (77800000/POW(P[y][z][x]*1000 + 354,2));
                        eqn2_sol_bs            =  661.85 + 3.5*P[y][z][x] + 6.2*POW(P[y][z][x],2);
                        if ((P[y][z][x] <  P1_thresh))    {T_solidus_bs = eqn1_sol_bs;}
                        else if((P[y][z][x] >= P1_thresh) && (P[y][z][x] < P2_thresh)) {T_solidus_bs = eqn2_sol_bs;}
                        else                              {T_solidus_bs = 0;}
                        T_liquidus_bs          =  1149.85 + 105*P[y][z][x];
                        T_solidus_pd           = -5.1 * POW(P[y][z][x],2) + 133*P[y][z][x] + 1120.85;
                        T_liquidus_pd          =  1799.85 + 114 * P[y][z][x];

                        // T/P conditions for LP/MP-TTG  
                        PKbar[y][z][x]         = P[y][z][x] ;
                        T1_cond                = (T[y][z][x] >= (760 - (60 * POW((PKbar[y][z][x] - 1), 2)))) && (T[y][z][x] <= (1000 - (150 * POW(((PKbar[y][z][x] - 1.2) / 1.2), 2))));
                        P1_cond                = (PKbar[y][z][x] >= (-0.5 * ((T[y][z][x] - 870) / 220)));
                        mask_lpttg             = PKbar[y][z][x] < 1;
                        LP_TTG                 = T1_cond && P1_cond && mask_lpttg;
                        mask_mpttg             = (PKbar[y][z][x] >= 1) && (PKbar[y][z][x] < (1.5 + (0.7 * ((T[y][z][x] - 700) / 200))));
                        MP_TTG                 = T1_cond && P1_cond && mask_mpttg;
                           // HP-TTG
//                        T2_cond                = (T[y][z][x] >= 1000) && (T[y][z][x] <= 1100 + 50 * POW(((PKbar[y][z][x] - 3.5) / 3.5), 2));
//                        P2_cond                = (PKbar[y][z][x] >= (2.35 + 0.15 * ((T[y][z][x] - 1000) / 100))) && (PKbar[y][z][x] <= 5);
//                        HP_TTG                 = T2_cond && P2_cond;

                        // melt calculation
                        if ((T[y][z][x] > T_solidus_bs) && (T[y][z][x] < T_liquidus_bs) && (T_liquidus_bs > T_solidus_bs) && (P[y][z][x] < P2_thresh)){
                             M0_bs = ((T[y][z][x] - T_solidus_bs)/(T_liquidus_bs - T_solidus_bs));}
                        else if (T[y][z][x] > T_liquidus_bs) {M0_bs = 1.0;}
                        else {M0_bs = 0;}
                        if ((T[y][z][x] > T_solidus_pd) && (T[y][z][x] < T_liquidus_pd) && (T_liquidus_pd > T_solidus_pd)){
                             M0_pd  = ((T[y][z][x] - T_solidus_pd)/(T_liquidus_pd - T_solidus_pd));}
                        else if (T[y][z][x] > T_liquidus_pd) {M0_pd = 1.0;}
                        else {M0_pd = 0;}
                        if ((M0_bs - cumulative_extracted_bs[y][z][x]) > 0){M_bs[y][z][x] = M0_bs - cumulative_extracted_bs[y][z][x];}
                        else {M_bs[y][z][x] = 0;}
                        if ((M0_pd - cumulative_extracted_pd[y][z][x]) > 0){M_pd[y][z][x] = M0_pd - cumulative_extracted_pd[y][z][x];}
                        else {M_pd[y][z][x] = 0;}
                        if (L_ed[y][z][x] > 0.1) {L_den[y][z][x] = L_ed[y][z][x];}
                        else                     {L_den[y][z][x] = 0;}
                        // Basalt
                        if  ((M_bs[y][z][x] > M2_bs) && ((M_bs[y][z][x] - M1_bs) < M4_bs) && ((M_bs[y][z][x] - M1_bs) > 0) && LP_TTG) {
                              M_ext_bs_LP_TTG[y][z][x] = (M_bs[y][z][x] - M1_bs) * L_den[y][z][x] * mask_bs[y][z][x]; }
                        else {M_ext_bs_LP_TTG[y][z][x] = 0; }
                        if  ((M_bs[y][z][x] > M2_bs) && ((M_bs[y][z][x] - M1_bs) < M4_bs) && ((M_bs[y][z][x] - M1_bs) > 0) && MP_TTG) {
                              M_ext_bs_MP_TTG[y][z][x] = (M_bs[y][z][x] - M1_bs) * L_den[y][z][x] * mask_bs[y][z][x]; }
                        else {M_ext_bs_MP_TTG[y][z][x] = 0; }
//                        if  ((M_bs[y][z][x] > M2_bs) && ((M_bs[y][z][x] - M1_bs) < M4_bs) && ((M_bs[y][z][x] - M1_bs) > 0) && HP_TTG) {
//                              M_ext_bs_HP_TTG[y][z][x] = (M_bs[y][z][x] - M1_bs) * L_den[y][z][x] * mask_bs[y][z][x]; }
//                        else {M_ext_bs_HP_TTG[y][z][x] = 0; }
                        M_ext_bs[y][z][x] = M_ext_bs_LP_TTG[y][z][x] + M_ext_bs_MP_TTG[y][z][x]; // + M_ext_bs_HP_TTG[y][z][x];
                        // Peridotite
                        if ((M_pd[y][z][x] > M2_pd) && ((M_pd[y][z][x] - M1_pd) < M4_pd) && ((M_pd[y][z][x] - M1_pd) > 0)) {
                             M_ext_pd[y][z][x] = (M_pd[y][z][x] - M1_pd) * (1 - L_den[y][z][x]) * mask_pd[y][z][x];}
                        else {M_ext_pd[y][z][x] = 0;}
                        if  ((cumulative_extracted_LP[y][z][x] + M_ext_bs_LP_TTG[y][z][x]) < M4_bs) {cumulative_extracted_LP[y][z][x] += M_ext_bs_LP_TTG[y][z][x];}
                        if  ((cumulative_extracted_MP[y][z][x] + M_ext_bs_MP_TTG[y][z][x]) < M4_bs) {cumulative_extracted_MP[y][z][x] += M_ext_bs_MP_TTG[y][z][x];}
//                        if  ((cumulative_extracted_HP[y][z][x] + M_ext_bs_HP_TTG[y][z][x]) < M4_bs) {cumulative_extracted_HP[y][z][x] += M_ext_bs_HP_TTG[y][z][x];}
                        if ((cumulative_extracted_pd[y][z][x] + M_ext_pd[y][z][x]) < M4_pd) {cumulative_extracted_pd[y][z][x] += M_ext_pd[y][z][x];}
                        if ((cumulative_extracted_bs[y][z][x] + M_ext_bs[y][z][x]) < M4_bs) {cumulative_extracted_bs[y][z][x] += M_ext_bs[y][z][x];}
                    }
                }

//---------------------------------------------------------------
//              Thermal boundary conditions and temperature
//---------------------------------------------------------------

                if(boundary_uu) { // Upper thermal BC
                    for(x=startx;x<endx;x+=incx) { eta[y][z][x] = eta_u[y][x]; }
                }
                if(boundary_dd) { // Lower thermal BC
                    for(x=startx;x<endx;x+=incx) { eta[y][z][x] = eta_d[y][x]; }
                }

//---------------------------------------------------------------
//              Equilibrium and non-equilibrium distributions
//---------------------------------------------------------------

                if(D==2) {
                    if(shifted_distribution) {
                        for(x=startx;x<endx;x+=incx) {
                            u2[x]    = 1.5 * (u[y][z][0][x] * u[y][z][0][x] + u[y][z][1][x] * u[y][z][1][x]);
                            cu[0][x] = - u2[x];
                        }
                        for(a=1;a<na;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                cu[a][x] = 3 *(c[a][0] * u[y][z][0][x] + c[a][1] * u[y][z][1][x]);
                                cu[a][x] = cu[a][x] + 0.5 * cu[a][x] * cu[a][x] - u2[x];
                            }
                        }
                    }
                    else {
                        for(x=startx;x<endx;x+=incx) {
                            u2[x]    = 1.5 * (u[y][z][0][x] * u[y][z][0][x] + u[y][z][1][x] * u[y][z][1][x]);
                            cu[0][x] = 1 - u2[x];
                        }
                        for(a=1;a<na;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                cu[a][x] = 3 *(c[a][0] * u[y][z][0][x] + c[a][1] * u[y][z][1][x]);
                                cu[a][x] = 1 + cu[a][x] + 0.5 * cu[a][x] * cu[a][x] - u2[x];
                            }
                        }
                    }
                }
                else if(D==3) {
                    if(shifted_distribution) {
                        for(x=startx;x<endx;x+=incx) {
                            u2[x]    = 1.5 * (u[y][z][0][x] * u[y][z][0][x] + u[y][z][1][x] * u[y][z][1][x] + u[y][z][2][x] * u[y][z][2][x]);
                            cu[0][x] = - u2[x];
                        }
                        for(a=1;a<na;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                cu[a][x] = 3 * (c[a][0] * u[y][z][0][x] + c[a][1] * u[y][z][1][x] + c[a][2] * u[y][z][2][x]);
                                cu[a][x] = cu[a][x] + 0.5 * cu[a][x] * cu[a][x] - u2[x];
                            }
                        }
                    }
                    else {
                        for(x=startx;x<endx;x+=incx) {
                            u2[x] = 1.5 * (u[y][z][0][x] * u[y][z][0][x] + u[y][z][1][x] * u[y][z][1][x] + u[y][z][2][x] * u[y][z][2][x]);
                            cu[0][x] = 1 - u2[x];
                        }
                        for(a=1;a<na;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                cu[a][x] = 3 * (c[a][0] * u[y][z][0][x] + c[a][1] * u[y][z][1][x] + c[a][2] * u[y][z][2][x]);
                                cu[a][x] = 1 + cu[a][x] + 0.5 * cu[a][x] * cu[a][x] - u2[x];
                            }
                        }
                    }
                }
                if(shifted_distribution) {
                    for(a=0;a<na;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f_eq [a][x] = w[a] * (rho[y][z][x] * cu[a][x] + drho[x]);
                            f_neq[a][x] = f[y][z][a][x] - f_eq[a][x];
                            g_neq[a][x] = g[y][z][a][x] - (eta [y][z][x] * (f_eq[a][x] + w[a]) + Delta_g[a]);
                            if (strcmp(rheology,"archean2") == 0) {
                                // h is NOT shifted; equilibrium uses unshifted cu = cu[a]+1
                                h_eq [a][x] = L_ed[y][z][x] * w[a] * (cu[a][x] + 1);
                                h_neq[a][x] = h[y][z][a][x] - h_eq[a][x];
                            }
                        }
                    }
                }
                else {
                    for(a=0;a<na;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f_eq [a][x] = rho [y][z][x] * w[a] * cu[a][x];
                            f_neq[a][x] = f[y][z][a][x] - f_eq[a][x];
                            g_neq[a][x] = g[y][z][a][x] - (eta [y][z][x] * f_eq[a][x] + Delta_g[a]);
                            if (strcmp(rheology,"archean2") == 0) {
                                h_eq [a][x] = L_ed[y][z][x] * w[a] * cu[a][x];
                                h_neq[a][x] = h[y][z][a][x] - h_eq[a][x];}
                        }
                    }
                }
#if 0
                // DBG probe: after equilibrium / non-eq
                if (rank == 0 && y == _pry && z == _prz && (_dbg_call_fast % DBG_INTERVAL == 0)) {
                    double _sum_feq = 0, _sum_fneq = 0, _sum_gneq = 0;
                    for (int _aa = 0; _aa < na; _aa++) {
                        _sum_feq  += f_eq [_aa][_prx];
                        _sum_fneq += f_neq[_aa][_prx];
                        _sum_gneq += g_neq[_aa][_prx];
                    }
                    printf("DBG fast   call=%lld stage=EQUIL    y=%d z=%d x=%d sum_f_eq=%.9e sum_f_neq=%.9e sum_g_neq=%.9e\n",
                           _dbg_call_fast, y, z, _prx, _sum_feq, _sum_fneq, _sum_gneq);
                }
#endif

//---------------------------------------------------------------
//	        Calculate nu_linear, so any temperture dependence
//---------------------------------------------------------------

                if(strcmp(temperature,"reynolds") == 0 && gamma_T > 1.e-10)  { // Reynolds temperature dependence
                    for(x=startx;x<endx;x+=incx) {
                        nu_linear[x] = nu_f[z] * EXP(-gamma_T*eta[y][z][x]);
                    }
                }
                else if (strcmp(temperature,"arrhenius") == 0) {               // Arrhenius temperature dependence
                    for(x=startx;x<endx;x+=incx) {
                        T_nd[y][z][x] = (T[y][z][x] - T0_p)/dT_p;
                        arg[y][z][x]  = 1;
                        nu_linear[x]  =  nu_f[z] * EXP(arg[y][z][x]);
                    }
                }
		else {                                                         // No       temperature dependence
                    for(x=startx;x<endx;x+=incx) {
                        nu_linear[x]     = nu_f[z];
                    }
                }

//---------------------------------------------------------------
//	        Calculate nu total: allow for nu_plastic if needed
//---------------------------------------------------------------

                if(strcmp(rheology,"linear") == 0)  {  // nu linear if linear rheology
                    for(x=startx;x<endx;x+=incx) {
                        nu_total [y][z][x] = max(nu_linear[x],nu_min);
                        if (nu_total[y][z][x] < nu_min_cap) nu_total[y][z][x] = nu_min_cap;
                        if (nu_total[y][z][x] > nu_max_cap) nu_total[y][z][x] = nu_max_cap;
                        nu_f_last[y][z][x] = nu_total [y][z][x];
                        tau_f          [x] = nu_total [y][z][x] * 3 + 0.5;
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
                                        S[i][j][x] += f_neq[a][x] * itmp;
                                    }
                                }
                            }
                        }
                    }
                    for(x=startx;x<endx;x+=incx) {
                        strain_rate[x] = 0.0;
                        S_factor   [x] = 1.5/(rho[y][z][x] * tau_f_last[y][z][x]);
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
                        strain_rate[x] = SQRT(strain_rate[x]);
                        Tau        [x] = rho[y][z][x] * nu_f_last[y][z][x] * strain_rate[x];
                    }

//                  Plastic rheology

                    if (strcmp(nu_formula,"Tau") == 0) {            // nu_plastic formula Tau_y/Tau
                        for(x=startx;x<endx;x+=incx) {
                            nu_plastic[x] = nu_f[z] * POW(Tau_y[z]/(Tau[x]+1.e-5*Tau_y[z]),N[z]);
                        }
                    }
                    
                    if (strcmp(nu_formula,"strain_rate") == 0) {    // nu_plastic formula Tau_y/(2*strain_rate)
                        for(x=startx;x<endx;x+=incx) {
                            nu_plastic[x] = nu_star + POW(Tau_y[z]/((2 * strain_rate[x]) + (1.e-5*Tau_y[z])),N[z]);
                        }
                    }

//                  Calculate total viscosity due to liner & plastic terms, and relaxation time

                    for(x=startx;x<endx;x+=incx) {
                        nu_total  [y][z][x] = max(1/(1/nu_linear[x] + 1/nu_plastic[x]),nu_min);
                        if (nu_total[y][z][x] < nu_min_cap) nu_total[y][z][x] = nu_min_cap;
                        if (nu_total[y][z][x] > nu_max_cap) nu_total[y][z][x] = nu_max_cap;
                        nu_f_last [y][z][x] = nu_total[y][z][x];
                        tau_f           [x] = nu_total[y][z][x] * 3 + 0.5;
                        tau_f_last[y][z][x] = tau_f[x];
                    }
                }

//---------------------------------------------------------------
//              Add collision terms
//---------------------------------------------------------------

#if 0
                // DBG probe: after nu_total / tau_f
                if (rank == 0 && y == _pry && z == _prz && (_dbg_call_fast % DBG_INTERVAL == 0)) {
                    printf("DBG fast   call=%lld stage=NU_TOT   y=%d z=%d x=%d nu_total=%.9e tau_f=%.9e tau_g=%.9e\n",
                           _dbg_call_fast, y, z, _prx,
                           nu_total[y][z][_prx], tau_f[_prx], tau_g[z]);
                }
#endif

                for(a=0;a<na;a++) {
                    for(x=startx;x<endx;x+=incx) {
                        f[y][z][a][x] -= f_neq[a][x]/tau_f[x];
                        g[y][z][a][x] -= g_neq[a][x]/tau_g[z];
                        if (strcmp(rheology,"archean2") == 0){
                            h[y][z][a][x] -= (h_neq[a][x]/tau_h[z]);}
                    }
                }
#if 0
                // DBG probe: after BGK collision (pre-buoyancy)
                if (rank == 0 && y == _pry && z == _prz && (_dbg_call_fast % DBG_INTERVAL == 0)) {
                    double _sum_f = 0, _sum_g = 0;
                    for (int _aa = 0; _aa < na; _aa++) {
                        _sum_f += f[y][z][_aa][_prx];
                        _sum_g += g[y][z][_aa][_prx];
                    }
                    printf("DBG fast   call=%lld stage=COLLIDE  y=%d z=%d x=%d sum_f=%.9e sum_g=%.9e\n",
                           _dbg_call_fast, y, z, _prx, _sum_f, _sum_g);
                }
#endif

//---------------------------------------------------------------
//              Add Boussinesq bouyancy force (assumes that d = D-1 is the z-axis)
//---------------------------------------------------------------

                for(x=startx;x<endx;x+=incx) {
                    gi[x] = rho[y][z][x] * beta_g * eta[y][z][x];
                    if (strcmp(rheology,"archean2") == 0) {
                        P_garnet_in    = 0.014 * T[y][z][x]  - 5.4;
                        P_plag_out     = 0.020 * T[y][z][x]  + 4;
                        if ((M_ext_bs[y][z][x] > 0) && ((P[y][z][x]*10) >= P_garnet_in) && ((P[y][z][x]*10) <= P_plag_out)) {
                             coef_phase_bs  = 0.16 * (P[y][z][x]*10 - P_garnet_in) /(P_plag_out - P_garnet_in);}
                        else{coef_phase_bs  = 0.0;}
                        if ((M_ext_bs[y][z][x] > 0)) {coef_res_bs =  0.4  * cumulative_extracted_bs[y][z][x];}
                        else{coef_res_bs    = 0.0;}
                        if ((M_ext_pd[y][z][x] > 0)) {coef_res_pd = -0.04 * cumulative_extracted_pd[y][z][x];}
                        else{coef_res_pd    = 0.0;}
                        if ((M_pd[y][z][x] >= M4_pd)) {coef_res_pd =  0.4;}
                        coef_phase_pd         = 0.0;
                        rho_solid_pd[y][z][x] = rho0_pd  * (1 + coef_res_pd + coef_phase_pd);
                        rho_solid_bs[y][z][x] = rho0_bs  * (1 + coef_res_bs + coef_phase_bs);
                        rho_solid[y][z][x]    =    L_den[y][z][x]    * rho_solid_bs[y][z][x] + (1 - L_den[y][z][x]) * rho_solid_pd[y][z][x];
                        rho_eff[y][z][x]      =    L_den[y][z][x]    *(rho_solid_bs[y][z][x] * (1 - M_bs[y][z][x])  + rho_melt_bs * M_bs[y][z][x]) +
                                              (1 - L_den[y][z][x])   *(rho_solid_pd[y][z][x] * (1 - M_pd[y][z][x])  + rho_melt_pd * M_pd[y][z][x]);
                        gi[x]                += (((rho_eff[y][z][x]-rho_M[y][z][x])/rho_M[y][z][x])*g_l);
                    }
                }

                for(a=1;a<na;a++) {
                    if(c[a][D-1]>0) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][a][x] +=  (w[a] * gi[x]);
                        }
                    }
                    else if(c[a][D-1]<0) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][a][x] -=  (w[a] * gi[x]);
    		        }
                    }
    	        }
#if 0
                // DBG probe: after buoyancy force applied to f
                if (rank == 0 && y == _pry && z == _prz && (_dbg_call_fast % DBG_INTERVAL == 0)) {
                    double _sum_f = 0;
                    for (int _aa = 0; _aa < na; _aa++) _sum_f += f[y][z][_aa][_prx];
                    printf("DBG fast   call=%lld stage=BUOY     y=%d z=%d x=%d gi=%.9e sum_f=%.9e beta_g=%.9e\n",
                           _dbg_call_fast, y, z, _prx, gi[_prx], _sum_f, beta_g);
                }
#endif
//---------------------------------------------------------------
//---------------------------------------------------------------
            }
        }
    }

