//-------------------------------------------------------------------------
//  C code for GTXLBM with overlapping option of calcuations and
//  communications (overlap bool)
//  Allows for interior or edges to be calculated according to whether
//  interior = true or false
//-------------------------------------------------------------------------

// I need to speak to Henry (T) to understand whether pragma block loop factor can
// get better performance, which it seems is should. Namely, from the description,
// pragma block_loop factor directives should be able to be applied to the loops
// where the loop levels for the GTXLBM structure would be 1=x,2=a,3=z,4=y.
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

    FLOAT f_eq[na][nx], f_neq[na][nx], g_eq[na][nx], g_neq[na][nx];
    FLOAT drho[nx], cu[na][nx], u2[nx];
    FLOAT strain_rate[nx], S[D][D][nx], S_factor[nx], Tau[nx];
    FLOAT nu_linear[nx], nu_plastic[nx], nu_total[nx], tau_f[nx], gi[nx];

//---------------------------------------------------------------
//  Set the x loop specifications as required for both the
//  case of no overlapping and overlapping which calculates
//  either on interior regions or on edges of domains
//---------------------------------------------------------------

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
                    }
                }

//---------------------------------------------------------------
//              Boundary conditions
//---------------------------------------------------------------

                // Left BC

                if(boundary_l) {
                    for(a=0;a<naR;a++) {
                        f[y][z][aR[a]][xl] = f[y][z][aiR[a]][xl];
                        g[y][z][aR[a]][xl] = g[y][z][aiR[a]][xl];
                    }
                    for(a=naR;a<nbR;a++) {
                        tmp = 0.5*(f[y][z][aR[a]][xl] + f[y][z][aiR[a]][xl]);
                        f[y][z][aR [a]][xl] = tmp;
                        f[y][z][aiR[a]][xl] = tmp;
                        tmp = 0.5*(g[y][z][aR[a]][xl] + g[y][z][aiR[a]][xl]);
                        g[y][z][aR [a]][xl] = tmp;
                        g[y][z][aiR[a]][xl] = tmp;
                   }
                }

                // Right BC

                if(boundary_r) {
                    for(a=0;a<naL;a++) {
                        f[y][z][aL[a]][xr] = f[y][z][aiL[a]][xr];
                        g[y][z][aL[a]][xr] = g[y][z][aiL[a]][xr];
                    }
                    for(a=naL;a<nbL;a++) {
                        tmp = 0.5*(f[y][z][aL[a]][xr] + f[y][z][aiL[a]][xr]);
                        f[y][z][aL [a]][xr] = tmp;
                        f[y][z][aiL[a]][xr] = tmp;
                        tmp = 0.5*(g[y][z][aL[a]][xr] + g[y][z][aiL[a]][xr]);
                        g[y][z][aL [a]][xr] = tmp;
                        g[y][z][aiL[a]][xr] = tmp;
                    }
                 }

                // Upper BC

                if(boundary_uu) {
                    for(a=0;a<naD;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][aD[a]][x] = f[y][z][aiD[a]][x];
                        }
                    }
                    for(a=naD;a<nbD;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            tmp = 0.5*(f[y][z][aD[a]][x] + f[y][z][aiD[a]][x]);
                            f[y][z][aD [a]][x] = tmp;
                            f[y][z][aiD[a]][x] = tmp;
                        }
                    }

                    // Set eta = eta_u on the upper boundary for state = solid => static_BC

                    if(strcmp(state,"solid")==0) {
                        for(a=0;a<na;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                g[y][z][a][x] = eta_u[y][x]*w[a];
                            }
                        }
                    }
                    else {
                        for(a=0;a<naD;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                g[y][z][aD[a]][x] = g[y][z][aiD[a]][x];
                            }
                        }
                        for(a=naD;a<nbD;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                tmp = 0.5*(g[y][z][aD[a]][x] + g[y][z][aiD[a]][x]);
                                g[y][z][aD [a]][x] = tmp;
                                g[y][z][aiD[a]][x] = tmp;
                            }
                        }
                    }
                }

                // Lower BC

                if(boundary_dd) {
                    for(a=0;a<naU;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][aU[a]][x] = f[y][z][aiU[a]][x];
                        }
                    }
                    for(a=naU;a<nbU;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            tmp = 0.5*(f[y][z][aU[a]][x] + f[y][z][aiU[a]][x]);
                            f[y][z][aU [a]][x] = tmp;
                            f[y][z][aiU[a]][x] = tmp;
                        }
                    }

                    // Set eta = eta_d on the lower boundary for state = solid => static_BC

                    if(strcmp(state,"solid")==0&&boundary_dd) {
                        for(a=0;a<na;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                g[y][z][a][x] = eta_d[y][x]*w[a];
                            }
                        }
                    }
                    else {
                       for(a=0;a<naU;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                g[y][z][aU[a]][x] = g[y][z][aiU[a]][x];
                            }
                        }
                        for(a=naU;a<nbU;a++) {
                            for(x=startx;x<endx;x+=incx) { 
                                tmp = 0.5*(g[y][z][aU[a]][x] + g[y][z][aiU[a]][x]);
                                g[y][z][aU [a]][x] = tmp;
                                g[y][z][aiU[a]][x] = tmp;
                            }
                        }
                    }
                }

                // Backward BC

                if(boundary_bb) {
                    for(a=0;a<naF;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][aF[a]][x] = f[y][z][aiF[a]][x];
                            g[y][z][aF[a]][x] = g[y][z][aiF[a]][x];
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
                        }
                    }
                }

                // Forward BC

                if(boundary_ff) {
                    for(a=0;a<naB;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][aB[a]][x] = f[y][z][aiB[a]][x];
                            g[y][z][aB[a]][x] = g[y][z][aiB[a]][x];
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
                        }
                    }
                }

//---------------------------------------------------------------
//              Macroscopic properties. Backards a loop to
//              optimize precision by summing small values first
//---------------------------------------------------------------

                for(x=startx;x<endx;x+=incx) {
                    rho[y][z][x] = f[y][z][na-1][x];
                    eta[y][z][x] = g[y][z][na-1][x];
                }
                for(a=na-2;a>=0;a--) {
                    for(x=startx;x<endx;x+=incx) {
                        rho[y][z][x] += f[y][z][a][x];
                        eta[y][z][x] += g[y][z][a][x];
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
//              Equilibrium and non-equilibrium distributions
//---------------------------------------------------------------

                if(D==2) {
                    if(shifted_distribution) {
//                        FLOAT factor1=1/SQRT(dx[y][z][startx]), factor2=1/dx[y][z][startx];
                        FLOAT factor1=1, factor2=1;
                        for(x=startx;x<endx;x+=incx) {
                            u2[x]    = 1.5 * factor2 * (u[y][z][0][x] * u[y][z][0][x] + u[y][z][1][x] * u[y][z][1][x]); // ooo
                            cu[0][x] = - u2[x];
                        }
                        for(a=1;a<na;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                cu[a][x] = 3 * factor1 * (c[a][0] * u[y][z][0][x] + c[a][1] * u[y][z][1][x]); // ooo
                                cu[a][x] = cu[a][x] + 0.5 * cu[a][x] * cu[a][x] - u2[x];
                            }
                        }
                    }
                    else {
                        for(x=startx;x<endx;x+=incx) {
                            u2[x]    = 1.5 * (u[y][z][0][x] * u[y][z][0][x] + u[y][z][1][x] * u[y][z][1][x]);
                            cu[0][x] = 1 - u2[x]; // ooo
                        }
                        for(a=1;a<na;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                cu[a][x] = 3 * (c[a][0] * u[y][z][0][x] + c[a][1] * u[y][z][1][x]);
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
                    if(strcmp(state,"magma")==0||!(boundary_uu||boundary_dd)) {
                        for(a=0;a<na;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                f_eq [a][x] = w[a] * (rho[y][z][x] * cu[a][x] + drho[x]);
                                f_neq[a][x] = f[y][z][a][x] - f_eq[a][x];
                                g_neq[a][x] = g[y][z][a][x] - (eta [y][z][x] * (f_eq[a][x] + w[a]));
                            }
                        }
                    }

                    // Set upper and lower boundaries to thermal equilibrium for state = solid

                    else if(strcmp(state,"solid")==0) {
                        if(boundary_uu) {
                            for(a=0;a<na;a++) {
                                for(x=startx;x<endx;x+=incx) {
                                    f_eq   [a][x] = w[a] * (rho[y][z][x] * cu[a][x] + drho[x]);
                                    f_neq  [a][x] = f[y][z][a][x] - f_eq[a][x];
                                    g[y][z][a][x] = eta_u[y][x] * (f_eq[a][x] + w[a]);
                                }
                            }
                        }
                        else if(boundary_dd) {
                            for(a=0;a<na;a++) {
                                for(x=startx;x<endx;x+=incx) {
                                    f_eq   [a][x] = w[a] * (rho[y][z][x] * cu[a][x] + drho[x]);
                                    f_neq  [a][x] = f[y][z][a][x] - f_eq[a][x];
                                    g[y][z][a][x] = eta_d[y][x] * (f_eq[a][x] + w[a]);
                                }
                            }
                        }
                    }
                }
                else {
                    if(strcmp(state,"magma")==0||!(boundary_uu||boundary_dd)) {
                        for(a=0;a<na;a++) {
                            for(x=startx;x<endx;x+=incx) {
                                f_eq [a][x] = rho [y][z][x] * w[a] * cu[a][x];
                                f_neq[a][x] = f[y][z][a][x] - f_eq[a][x];
                                g_neq[a][x] = g[y][z][a][x] - (eta [y][z][x] * f_eq[a][x]);
                            }
                        }
                    }

                    // Set upper and lower boundaries to thermal equilibrium for state = solid

                    else if(strcmp(state,"solid")==0) {
                        if(boundary_uu) {
                            for(a=0;a<na;a++) {
                                for(x=startx;x<endx;x+=incx) {
                                    f_eq   [a][x] = rho [y][z][x] * w[a] * cu[a][x];
                                    f_neq  [a][x] = f[y][z][a][x] - f_eq[a][x];
                                    g[y][z][a][x] = eta_u  [y][x] * f_eq[a][x];
                                }
                            }
                        }
                        else if(boundary_dd) {
                            for(a=0;a<na;a++) {
                                for(x=startx;x<endx;x+=incx) {
                                    f_eq   [a][x] = rho [y][z][x] * w[a] * cu[a][x];
                                    f_neq  [a][x] = f[y][z][a][x] - f_eq[a][x];
                                    g[y][z][a][x] = eta_d  [y][x] * f_eq[a][x];
                                }
                            }
                        }
                    }
                }

//---------------------------------------------------------------
//	        Calculate nu_linear, so any temperture dependence
//---------------------------------------------------------------

                if(strcmp(Tdep_viscosity,"reynolds") == 0 && gamma_T > 1.e-10)  { // Reynolds temperature dependence
                    for(x=startx;x<endx;x+=incx) {
                        nu_linear[x] = nu_f[z] * EXP(-gamma_T*eta[y][z][x]);
                    }
                }
		else {                                                            // No       temperature dependence
                    for(x=startx;x<endx;x+=incx) {
                        nu_linear[x] = nu_f[z];
                    }
                }

//---------------------------------------------------------------
//	        Calculate nu total: allow for nu_plastic if needed
//---------------------------------------------------------------

                if(strcmp(rheology,"linear") == 0)  {  // Linear rheology
                    for(x=startx;x<endx;x+=incx) {
                        nu_total       [x] = max(nu_linear[x],nu_min);
                        nu_f_last[y][z][x] = nu_total [x];
                        tau_f          [x] = nu_total [x] * 3 + 0.5;
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
                        strain_rate[x] = SQRT(0.5*strain_rate[x]);
                        Tau        [x] = rho[y][z][x] * nu_f_last[y][z][x] * strain_rate[x];
                    }

//                  Plastic rheology: (yield stress and strain rate dependent power law rheology: nu_plastic = nu_0 . (Tau_y/Tau)^N where Tau = rho.nu_0.strain_rate

                    if (strcmp(rheology,"plastic") == 0 || strcmp(rheology,"plastic_switch") == 0) {
                        for(x=startx;x<endx;x+=incx) {
                            nu_plastic[x] = nu_f[z] * POW(Tau_y[z]/(Tau[x]+1.e-5*Tau_y[z]),N[z]);
                        }
                    }

//                  Calculate total viscosity due to linear & plastic terms, and relaxation time

                    if (strcmp(rheology,"plastic") == 0) {            // Plastic rheology: nu_total = (1/nu_linear+1/nu_plastic)^(-1)
                        for(x=startx;x<endx;x+=incx) {
                            nu_total        [x] = max(1/(1/nu_linear[x] + 1/nu_plastic[x]),nu_min);
                            nu_f_last [y][z][x] = nu_total[x];
                            tau_f           [x] = nu_total[x] * 3 + 0.5;
                            tau_f_last[y][z][x] = tau_f[x];
                        }
                    }
                    else if (strcmp(rheology,"plastic_switch") == 0) { // Plastic switch rheology: nu_total = nu_linear (strain_rate<=Tau_y), nu_total = nu_plastic (strain_rate>Tau_y)
                        for(x=startx;x<endx;x+=incx) {
                           if(strain_rate[x]<=Tau_y[z]) { nu_total[x] = max(nu_linear [x],nu_min); }
                           else                         { nu_total[x] = max(nu_plastic[x],nu_min); }
                        }
                   }
                }

//---------------------------------------------------------------
//              Add collision terms and reduce energy
//              due to radiogenic decay
//---------------------------------------------------------------

                if(strcmp(state,"magma")==0||!(boundary_uu||boundary_dd)) {
                    for(a=0;a<na;a++) {
//                        FLOAT factor = cylindrical_factor[a][z]; ooo
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][a][x] -=  f_neq[a][x]/tau_f[x];
                            g[y][z][a][x] -= (g_neq[a][x]/tau_g[z] + Delta_g[a]);
                        }
                    }
                }

               // On the upper and lower boundaries for state = solid, only update f as g was already set to equilibrium

                else {
                    for(a=0;a<na;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][a][x] -=  f_neq[a][x]/tau_f[x];
                        }
                    }
                }

//---------------------------------------------------------------
//              Add Boussinesq bouyancy (d = D-1 is the z-axis)
//---------------------------------------------------------------

                for(x=startx;x<endx;x+=incx) {
                    gi[x] = beta_g[z] * eta[y][z][x];
                }
                for(a=1;a<na;a++) {
                    if(c[a][D-1]>0) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][a][x] += (w[a] * gi[x]);
                        }
                    }
                    else if(c[a][D-1]<0) {
                        for(x=startx;x<endx;x+=incx) {
                            f[y][z][a][x] -= (w[a] * gi[x]);
    		        }
                    }
    	        }
            }
        }
    }
