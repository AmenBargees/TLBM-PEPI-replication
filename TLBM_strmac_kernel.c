//-------------------------------------------------------------------------
//  C code for strmac
//-------------------------------------------------------------------------
//  Archean2 sub-block toggle (1=enabled, 0=disabled). Starts OFF so the
//  archean2 input file runs without any archean2 work in C. Flip 0 -> 1
//  to enable this single block in isolation. See TLBM_rest_kernel.c for
//  the other A2_* toggles. NOTE: with A2_HSTRM=0 the h-distribution does
//  not stream and BCs are not applied to it, so h stays at Python init
//  values throughout the run.
#define A2_HSTRM 1   // h streaming and h boundary conditions
//-------------------------------------------------------------------------
    int incx, startx, endx, starty, endy, startz, endz;
    int x, y, z, a, d, i, j, xx_xa, yy_ya, zz_za;
    FLOAT tmp;
    bool edgey, edgez, x_loop, boundary_bb, boundary_ff, boundary_uu, boundary_dd;
//-------------------------------------------------------------------------
    startx = xedge  ; endx = mx-xedge;
    starty = yedge  ; endy = my-yedge;
    startz = zedge  ; endz = mz-zedge;
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

            incx = 1;
 
//---------------------------------------------------------------
//              Streaming
//---------------------------------------------------------------

                for(a=0;a<na;a++) {
                    yy_ya = y_ya[a][y];
                    zz_za = z_za[a][z];
                    for(x=startx;x<endx;x++) {
                        xx_xa                = x_xa[a][x];
                        f_r[y][z][a][x]      = ff_r[yy_ya][zz_za][a][xx_xa];
                        f_b[y][z][a][x]      = ff_b[yy_ya][zz_za][a][xx_xa];
                        g  [y][z][a][x]      = gg  [yy_ya][zz_za][a][xx_xa];
                        if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                            h  [y][z][a][x]  = hh  [yy_ya][zz_za][a][xx_xa];
                            }
                    }
                }

//---------------------------------------------------------------
//              Boundary conditions
//---------------------------------------------------------------
                
                // Left BC — halfway bounce-back from pre-streaming ff_r[aiR]
                if(boundary_l) {
                    for(a=0;a<naR;a++) {
                        f_r[y][z][aR[a]][xl] = ff_r[y][z][aiR[a]][xl];
                        f_b[y][z][aR[a]][xl] = ff_b[y][z][aiR[a]][xl];
                        g  [y][z][aR[a]][xl] = gg  [y][z][aiR[a]][xl];
                        if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                            h  [y][z][aR[a]][xl] = hh  [y][z][aiR[a]][xl];
                        }
                    }
                    for(a=naR;a<nbR;a++) {
                        tmp = 0.5*(f_r[y][z][aR[a]][xl] + f_r[y][z][aiR[a]][xl]);
                        f_r[y][z][aR [a]][xl] = tmp;
                        f_r[y][z][aiR[a]][xl] = tmp;
                        tmp = 0.5*(f_b[y][z][aR[a]][xl] + f_b[y][z][aiR[a]][xl]);
                        f_b[y][z][aR [a]][xl] = tmp;
                        f_b[y][z][aiR[a]][xl] = tmp;
                        tmp = 0.5*(g[y][z][aR[a]][xl] + g[y][z][aiR[a]][xl]);
                        g[y][z][aR [a]][xl] = tmp;
                        g[y][z][aiR[a]][xl] = tmp;
                        if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                            tmp = 0.5*(h[y][z][aR[a]][xl] + h[y][z][aiR[a]][xl]);
                            h[y][z][aR [a]][xl] = tmp;
                            h[y][z][aiR[a]][xl] = tmp;
                        }
                   }
                }
                // Right BC — halfway bounce-back from pre-streaming ff_r[aiL]
                if(boundary_r) {
                    for(a=0;a<naL;a++) {
                        f_r[y][z][aL[a]][xr] = ff_r[y][z][aiL[a]][xr];
                        f_b[y][z][aL[a]][xr] = ff_b[y][z][aiL[a]][xr];
                        g  [y][z][aL[a]][xr] = gg  [y][z][aiL[a]][xr];
                        if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                            h  [y][z][aL[a]][xr] = hh  [y][z][aiL[a]][xr];
                        }
                    }
                    for(a=naL;a<nbL;a++) {
                        tmp = 0.5*(f_r[y][z][aL[a]][xr] + f_r[y][z][aiL[a]][xr]);
                        f_r[y][z][aL [a]][xr] = tmp;
                        f_r[y][z][aiL[a]][xr] = tmp;
                        tmp = 0.5*(f_b[y][z][aL[a]][xr] + f_b[y][z][aiL[a]][xr]);
                        f_b[y][z][aL [a]][xr] = tmp;
                        f_b[y][z][aiL[a]][xr] = tmp;
                        tmp = 0.5*(g[y][z][aL[a]][xr]   + g  [y][z][aiL[a]][xr]);
                        g[y][z][aL [a]][xr] = tmp;
                        g[y][z][aiL[a]][xr] = tmp;
                        if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                            tmp = 0.5*(h[y][z][aL[a]][xr] + h[y][z][aiL[a]][xr]);
                            h[y][z][aL [a]][xr] = tmp;
                            h[y][z][aiL[a]][xr] = tmp;
                        }
                    }
                 }
                // Upper BC. Halfway bounce-back: reflect the pre-streaming f
                // that was going INTO the wall (ff_r[aiD]) back as f_r[aD].
                // Reading from ff_r (pre-streaming) is critical — using post-
                // streaming f_r[aiD] destroys the halo-reflected mass and
                // causes massive per-color leak.
                if(boundary_uu) {
                    for(a=0;a<naD;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f_r[y][z][aD[a]][x] = ff_r[y][z][aiD[a]][x];
                            f_b[y][z][aD[a]][x] = ff_b[y][z][aiD[a]][x];
                            g  [y][z][aD[a]][x] = gg  [y][z][aiD[a]][x];
                            if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                                h  [y][z][aD[a]][x] = hh  [y][z][aiD[a]][x];
                            }
                        }
                    }
                    for(a=naD;a<nbD;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            tmp = 0.5*(f_r[y][z][aD[a]][x] + f_r[y][z][aiD[a]][x]);
                            f_r[y][z][aD [a]][x] = tmp;
                            f_r[y][z][aiD[a]][x] = tmp;
                            tmp = 0.5*(f_b[y][z][aD[a]][x] + f_b[y][z][aiD[a]][x]);
                            f_b[y][z][aD [a]][x] = tmp;
                            f_b[y][z][aiD[a]][x] = tmp;
                            tmp = 0.5*(g[y][z][aD[a]][x]   + g  [y][z][aiD[a]][x]);
                            g[y][z][aD [a]][x] = tmp;
                            g[y][z][aiD[a]][x] = tmp;
                            if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                                tmp = 0.5*(h[y][z][aD[a]][x] + h[y][z][aiD[a]][x]);
                                h[y][z][aD [a]][x] = tmp;
                                h[y][z][aiD[a]][x] = tmp;
                            }
                        }
                    }
                }
                // Lower BC
                if(boundary_dd) {
                    for(a=0;a<naU;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f_r[y][z][aU[a]][x] = ff_r[y][z][aiU[a]][x];
                            f_b[y][z][aU[a]][x] = ff_b[y][z][aiU[a]][x];
                            g  [y][z][aU[a]][x] = gg  [y][z][aiU[a]][x];
                            if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                                h  [y][z][aU[a]][x] = hh  [y][z][aiU[a]][x];
                            }
                        }
                    }
                    for(a=naU;a<nbU;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            tmp = 0.5*(f_r[y][z][aU[a]][x] + f_r[y][z][aiU[a]][x]);
                            f_r[y][z][aU [a]][x] = tmp;
                            f_r[y][z][aiU[a]][x] = tmp;
                            tmp = 0.5*(f_b[y][z][aU[a]][x] + f_b[y][z][aiU[a]][x]);
                            f_b[y][z][aU [a]][x] = tmp;
                            f_b[y][z][aiU[a]][x] = tmp;
                            tmp = 0.5*(g[y][z][aU[a]][x]   + g  [y][z][aiU[a]][x]);
                            g[y][z][aU [a]][x] = tmp;
                            g[y][z][aiU[a]][x] = tmp;
                            if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                                tmp = 0.5*(h[y][z][aU[a]][x] + h[y][z][aiU[a]][x]);
                                h[y][z][aU [a]][x] = tmp;
                                h[y][z][aiU[a]][x] = tmp;
                            }
                        }
                    }
                }
                // Backward BC — halfway bounce-back from pre-streaming ff_r[aiF]
                if(boundary_bb) {
                    for(a=0;a<naF;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f_r[y][z][aF[a]][x] = ff_r[y][z][aiF[a]][x];
                            f_b[y][z][aF[a]][x] = ff_b[y][z][aiF[a]][x];
                            g  [y][z][aF[a]][x] = gg  [y][z][aiF[a]][x];
                            if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                                h  [y][z][aF[a]][x] = hh  [y][z][aiF[a]][x];
                            }
                        }
                    }
                    for(a=naF;a<nbF;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            tmp = 0.5*(f_r[y][z][aF[a]][x] + f_r[y][z][aiF[a]][x]);
                            f_r[y][z][aF [a]][x] = tmp;
                            f_r[y][z][aiF[a]][x] = tmp;
                            tmp = 0.5*(f_b[y][z][aF[a]][x] + f_b[y][z][aiF[a]][x]);
                            f_b[y][z][aF [a]][x] = tmp;
                            f_b[y][z][aiF[a]][x] = tmp;
                            tmp = 0.5*(g[y][z][aF[a]][x]   + g  [y][z][aiF[a]][x]);
                            g[y][z][aF [a]][x] = tmp;
                            g[y][z][aiF[a]][x] = tmp;
                            if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                                tmp = 0.5*(h[y][z][aF[a]][x] + h[y][z][aiF[a]][x]);
                                h[y][z][aF [a]][x] = tmp;
                                h[y][z][aiF[a]][x] = tmp;
                            }
                        }
                    }
                }
                // Forward BC — halfway bounce-back from pre-streaming ff_r[aiB]
                if(boundary_ff) {
                    for(a=0;a<naB;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            f_r[y][z][aB[a]][x] = ff_r[y][z][aiB[a]][x];
                            f_b[y][z][aB[a]][x] = ff_b[y][z][aiB[a]][x];
                            g  [y][z][aB[a]][x] = gg  [y][z][aiB[a]][x];
                            if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                                h  [y][z][aB[a]][x] = hh  [y][z][aiB[a]][x];
                            }
                        }
                    }
                    for(a=naB;a<nbB;a++) {
                        for(x=startx;x<endx;x+=incx) {
                            tmp = 0.5*(f_r[y][z][aB[a]][x] + f_r[y][z][aiB[a]][x]);
                            f_r[y][z][aB [a]][x] = tmp;
                            f_r[y][z][aiB[a]][x] = tmp;
                            tmp = 0.5*(f_b[y][z][aB[a]][x] + f_b[y][z][aiB[a]][x]);
                            f_b[y][z][aB [a]][x] = tmp;
                            f_b[y][z][aiB[a]][x] = tmp;
                            tmp = 0.5*(g[y][z][aB[a]][x]   + g  [y][z][aiB[a]][x]);
                            g[y][z][aB [a]][x] = tmp;
                            g[y][z][aiB[a]][x] = tmp;
                            if (A2_HSTRM && strcmp(rheology,"archean2") == 0) {
                                tmp = 0.5*(h[y][z][aB[a]][x] + h[y][z][aiB[a]][x]);
                                h[y][z][aB [a]][x] = tmp;
                                h[y][z][aiB[a]][x] = tmp;
                            }
                        }
                    }
                }
//---------------------------------------------------------------
//              Macroscopic properties. Backards a loop to
//              optimize precision by summing small values first.
//              With shifted_distribution the stored f̃_r/f̃_b are
//              deviations from rho_r_init/rho_b_init, so add the
//              reference back: rho_r = rho_r_init + Σ f̃_r.
//---------------------------------------------------------------

                for(x=startx;x<endx;x+=incx) {
                    rho_r[y][z][x] = f_r[y][z][na-1][x];
                    rho_b[y][z][x] = f_b[y][z][na-1][x];
                }
                for(a=na-2;a>=0;a--) {
                    for(x=startx;x<endx;x+=incx) {
                        rho_r[y][z][x] += f_r[y][z][a][x];
                        rho_b[y][z][x] += f_b[y][z][a][x];
                    }
                }
                if (shifted_distribution) {
                    for(x=startx;x<endx;x+=incx) {
                        rho_r[y][z][x] += rho_r_init[y][z][x];
                        rho_b[y][z][x] += rho_b_init[y][z][x];
                    }
                }
//---------------------------------------------------------------
        }
    }
