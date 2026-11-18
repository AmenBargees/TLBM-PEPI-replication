#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int cmpfunc (const void * a, const void * b) {
   return ( *(int*)b - *(int*)a );
}

void main(int argc, char *argv[] )
{
    if(argc!= 8) {
	printf("7 arguments are required, namely: D np size nx ny nz option\n");
	printf("option = 1 ... same aspect ratio\n");
	printf("option = 2 ... variable aspect ratio best fitting power 2 sizes\n");
	exit(-1);
    }
    int D    = atoi(argv[1]);
    int np   = atoi(argv[2]);
    int size = atoi(argv[3]);
    int nx   = atoi(argv[4]);
    int ny   = atoi(argv[5]);
    int nz   = atoi(argv[6]);
    int option = atoi(argv[7]);

    if(D==2) {
      if(option == 2) {
	int i = -1, nxout=nx, nzout=nz;
        while(nxout*nzout<nx*nz*np) {
            i ++;
            if     (i%2 == 0) { nzout *= 2; }
            else if(i%2 == 1) { nxout *= 2; }
        }
        if(nxout*nzout>nx*nz*np) {
            if     (i%2 == 0) { nzout /= 2; }
            else if(i%2 == 1) { nxout /= 2; }
	}
        if(nxout*nzout<nx*nz*np) {
            if(nxout>nzout) { nzout *= 3; nzout /= 2; }
            else            { nxout *= 3; nxout /= 2; }
	}
	nx = nxout * size;
	nz = nzout * size;
	int nxnynz[2];
	nxnynz[0] = nx;
	nxnynz[1] = nz;
	qsort(nxnynz,2,sizeof(int),cmpfunc);
	nx = nxnynz[0];
	nz = nxnynz[1];
      }
      else {

        float factor = sqrt(np) * size;
        float n = np * size * size * nx * nz;
        nx = (nx * factor);
        nz = (nz * factor);
	int i = -1;
	while(nx * nz < n) {
	    i ++;
	    if(i%2 == 0) { nx += 1; }
	    else         { nz += 1; }
	}
        float error0 = abs(n - (nx * nz));
	if(i%2 == 0) { nx -= 1; }
	else         { nz -= 1; }
        float error = n - (nx * nz);
	if(error > error0) {
	    if(i%2 == 0) { nx += 1; }
	    else         { nz += 1; }
	}
      }
      ny = 1;
    }
    else if(D == 3) {
      if(option == 2) {
	int i = -1, nxout=nx, nyout = ny, nzout=nz;
        while(nxout*nyout*nzout<nx*ny*nz*np) {
            i ++;
            if     (i%3 == 0) { nzout *= 2; }
            else if(i%3 == 1) { nyout *= 2; }
            else if(i%3 == 2) { nxout *= 2; }
        }
        if(nxout*nyout*nzout>nx*ny*nz*np) {
            if     (i%3 == 0) { nzout /= 2; }
            else if(i%3 == 1) { nyout /= 2; }
            else if(i%3 == 2) { nxout /= 2; }
	}
        if(nxout*nyout*nzout<nx*ny*nz*np) {
            if     (nxout>nzout) { nzout *= 3; nzout /= 2; }
            else if(nxout>nyout) { nyout *= 3; nyout /= 2; }
            else                 { nxout *= 3; nxout /= 2; }
	}
	nx = nxout * size;
	ny = nyout * size;
	nz = nzout * size;
	int nxnynz[3];
	nxnynz[0] = nx;
	nxnynz[1] = ny;
	nxnynz[2] = nz;
	qsort(nxnynz,3,sizeof(int),cmpfunc);
	nx = nxnynz[0];
	ny = nxnynz[1];
	nz = nxnynz[2];
      }
      else {
        float factor = pow(np,0.33333333) * size;
        float n = factor * factor * factor * nx * ny * nz;
        nx = (nx * factor);
        ny = (ny * factor);
        nz = (nz * factor);
	int i = -1;
	while(nx * ny * nz < n) {
	    i ++;
	    if(i%3   == 0)    { nx += 1; }
	    else if(i%3 == 1) { ny += 1; }
	    else if(i%3 == 2) { nz += 1; }
	}
        float error0 = abs(n - (nx * ny *  nz));
	if(i%3 == 0) { nx -= 1; }
	else if(i%3 == 1) { ny -= 1; }
	else if(i%3 == 2) { nz -= 1; }
        float error = n - (nx * ny * nz);
	if(error > error0) {
	    if(i%3 == 0)      { nx += 1; }
    	    else if(i%3 == 1) { ny += 1; }
	    else if(i%3 == 2) { nz += 1; }
	}
      }
    }
    printf("%d %d %d",nx,ny,nz);
}
