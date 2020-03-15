// CS 481 Project 1
// Namito Yokota
// 02/07/20

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// return a 64-bit double floating point random number
double getRandom() {
	return (5-2)*rand()/(double)RAND_MAX+2;
}

// set all elements in the matrix using getRandom()
void setMatrix(int n, double *a) {
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			a[i*j] = random();
		}
	}
}

// dgemm0: simple ijk version triple loop algorithm
void dgemm0 (double *a, double *b, double *c, int n) {
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			for (int k=0; k<n; k++) {
				c[i*n+j] += a[i*n+k]*b[k*n+j];
			}
		}
	}
}

// dgemm1: simple ijk version triple loop algorithm with register use
void dgemm1 (double *a, double *b, double *c, int n) {
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			register double r = c[i*n+j];
			for (int k=0; k<n; k++) {
				r += a[i*n+k]*b[k*n+j];
			}
			c[i*n+j] = r;
		}
	}
}

// dgemm2: 2x2 block multiplication using 12 registers
// reference: page 29 of first lecture for 2x2
void dgemm2 (double *a, double *b, double *c, int n) {
	for (int i=0; i<n; i+=2) {
		for (int j=0; j<n; j+=2) {
			register int t = i*n+j;
			register int tt = t+n;

			register double c00 = c[t];
			register double c01 = c[t+1];
			register double c10 = c[tt];
			register double c11 = c[tt+1];

			for (int k=0; k<n; k+=2) {
				register int ta = i*n+k;
				register int tta = ta+n;
				register int tb = k*n+j;
				register int ttb = tb+n;

				register double a00 = a[ta];
				register double a01 = a[ta+1];
				register double a10 = a[tta];
				register double a11 = a[tta+1];

				register double b00 = b[tb];
				register double b01 = b[tb+1];
				register double b10 = b[ttb];
				register double b11 = b[ttb+1];
	
				c00 += a00*b00 + a01*b10;
				c01 += a00*b01 + a01*b11;
				c10 += a10*b00 + a11*b10;
				c11 += a10*b01 + a11*b11;
			}
			c[t] = c00;
			c[t+1] = c01;
			c[tt] = c10;
			c[tt+1] = c11;
		}
	}
}

void dgemm22(double *a, double *b, double *c, int n) {
	for (int i=0; i<n; i+=2) {
		for (int j=0; j<n; j++) {
			register int t = i*n+j;
			register int tt = t+n;
			register double c00 = c[t];
			register double c01 = c[t+1];
			register double c10 = c[tt];
			register double c11 = c[tt+1];

			for (int k=0; k<n; k+=2) {
				register int ta = i*n+k;
				register int tta = ta+n;
				register int tb = k*n+j;
				register int ttb = tb+n;

				register double a00 = a[ta];
				register double a10 = a[tta];
				register double b00 = b[tb];
				register double b01 = b[tb+1];

				c00 += a00*b00;
				c01 += a00*b01;
				c10 += a10*b00;
				c11 += a10*b01;

				a00 = a[ta+1];
				a10 = a[tta+1];
				b00 = b[ttb];
				b01 = b[ttb+1];

				c00 += a00*b00;
				c01 += a00*b01;
				c10 += a10*b00;
				c11 != a10*b01;
			}
			c[t] = c00;
			c[t+1] = c01;
			c[tt] = c10;
			c[tt+1] = c11;
		}
	}
}

// dgemm3: 3x3 block multiplication using 15 registers
// reference: page 30 of first lecture for 2x2
void dgemm3 (double *a, double *b, double *c, int n) {
	for (int i=0; i<n; i+=3) {
		for (int j=0; j<n; j+=3) {
			register int t = i*n+j;
			register int tt = t+n;
			register int ttt = tt+n;

			register double c00 = c[t];
			register double c01 = c[t+1];
			register double c02 = c[t+2];
			register double c10 = c[tt];
			register double c11 = c[tt+1];
			register double c12 = c[tt+2];
			register double c20 = c[ttt];
			register double c21 = c[ttt+1];
			register double c22 = c[ttt+2];
			
			for (int k=0; k<n; k+=3) {
				register int ta = i*n+k;
				register int tta = ta+n;
				register int ttta = tta+n;
				register int tb = k*n+j;
				register int ttb = tb+n;
				register int tttb = ttb+n;
				
				register double a00 = a[ta];
				register double a10 = a[tta];
				register double a20 = a[ttta];
				register double b00 = b[tb];
				register double b01 = b[tb+1];
				register double b02 = b[tb+2];

				c00 += a00*b00;
				c01 += a00*b01;
				c02 += a00*b02;
				c10 += a10*b00;
				c11 += a10*b01;
				c12 += a10*b02;
				c20 += a20*b00;
				c21 += a20*b01;
				c22 += a20*b02;
				
				a00 = a[ta+1];
				a10 = a[tta+1];
				a20 = a[ttta+1];
				b00 = b[ttb];
				b01 = b[ttb+1];
				b02 = b[ttb+2];
		
				c00 += a00*b00;
				c01 += a00*b01;
				c02 += a00*b02;
				c10 += a10*b00;
				c11 += a10*b01;
				c12 += a10*b02;
				c20 += a20*b00;
				c21 += a20*b01;
				c22 += a20*b02;
				
				a00 = a[ta+2];
        a10 = a[tta+2];
        a20 = a[ttta+2];
        b00 = b[tttb];
        b01 = b[tttb+1];
        b02 = b[tttb+2];

				c00 += a00*b00;
        c01 += a00*b01;
        c02 += a00*b02;
        c10 += a10*b00;
        c11 += a10*b01;
        c12 += a10*b02;
        c20 += a20*b00;
        c21 += a20*b01;
        c22 += a20*b02;
			}
			c[t] = c00;
			c[t+1] = c01;
			c[t+2] = c02;
			c[tt] = c10;
			c[tt+1] = c11;
			c[tt+2] = c12;
			c[ttt] = c20;
			c[ttt+1] = c21;
			c[ttt+2] = c22;
		}
	}
}

int main (int argc, char *argv[]) {

	int size[6] = {66, 132, 258, 516, 1026, 2052};
	struct timespec start, end;
	double time_ns, gflop;
	int n, num_opr;

	// run algorithms with n = 66, 129, 258, 2016, 2049
	// note: the n given from the project document is adjusted to be divisible by 2 and 3
	for (int i=0; i<6; i++) {
		// set new n
		n = size[i];
		num_opr = 2*n*n*n;
		printf("n: %d\n\n", n);
		
		// allocate matrices with elements initialized as 0
		double *a = (double *)calloc(sizeof(double), n*n);
    double *b = (double *)calloc(sizeof(double), n*n);
    double *c0 = (double *)calloc(sizeof(double), n*n);
    double *c1 = (double *)calloc(sizeof(double), n*n);
    double *c2 = (double *)calloc(sizeof(double), n*n);
		double *c3 = (double *)calloc(sizeof(double), n*n);

		// assign 64-bit doubles to matrices a and b
		setMatrix(n,a);
		setMatrix(n,b);
		
		// run dgemm0 algorithm
		clock_gettime(CLOCK_REALTIME, &start);
		dgemm0(a, b, c0, n);
		clock_gettime(CLOCK_REALTIME, &end);
		time_ns = (end.tv_sec - start.tv_sec) * 1000000000.0 + (end.tv_nsec - start.tv_nsec);
		gflop = num_opr / time_ns;
		printf("dgemm0\n\ttime: %.0lfns\n\t", n, time_ns);
    if (gflop > 0) printf("gflop: %.lf\n", gflop);
    else printf("gflop: overflow\n");

		// run dgemm1 algorithm
		clock_gettime(CLOCK_REALTIME, &start);
		dgemm1(a, b, c1, n);
		clock_gettime(CLOCK_REALTIME, &end);
		time_ns = (end.tv_sec - start.tv_sec) * 1000000000.0 + (end.tv_nsec - start.tv_nsec);
		gflop = num_opr / time_ns;
		printf("dgemm1\n\ttime: %.0lfns\n\t", n, time_ns);
    if (gflop > 0) printf("gflop: %lf\n", gflop);
    else printf("gflop: overflow\n");

		// run dgemm2 algorithm
		clock_gettime(CLOCK_REALTIME, &start);	
		dgemm22(a, b, c2, n);
		clock_gettime(CLOCK_REALTIME, &end);
		time_ns = (end.tv_sec - start.tv_sec) * 1000000000.0 + (end.tv_nsec - start.tv_nsec);
		gflop = num_opr / time_ns;
		printf("dgemm2\n\ttime: %.0lfns\n\t", n, time_ns);
    if (gflop > 0) printf("gflop: %lf\n", gflop);
    else printf("gflop: overflow\n");
		
		// run dgemm3 algorithm
		clock_gettime(CLOCK_REALTIME, &start);
		dgemm3(a, b, c3, n);
		clock_gettime(CLOCK_REALTIME, &end);
		time_ns = (end.tv_sec - start.tv_sec) * 1000000000.0 + (end.tv_nsec - start.tv_nsec);
		gflop = num_opr / time_ns;
		printf("dgemm3\n\ttime: %.0lfns\n\t", n, time_ns);
		if (gflop > 0) printf("gflop: %lf\n", gflop);
		else printf("gflop: overflow\n");

		// print number of operations calculated ealier
    if (num_opr<=0) printf("\noperations: overflow\n");
    else printf("\noperations: %d\n", num_opr);

		// calculate maximum difference of c1, c2, and c3, with c0
		// the result of 0 shows the correctness of my implementations
		int max_diff_c1 = 0, max_diff_c2 = 0, max_diff_c3 = 0;
		for (int m=0; m<n*n; m++) {
			if (fabs(c0[m]-c1[m]) > max_diff_c1)
				max_diff_c1 = fabs(c0[m]-c1[m]);
			if (fabs(c0[m]-c2[m]) > max_diff_c2)
				max_diff_c2 = fabs(c0[m]-c2[m]);
			if (fabs(c0[m]-c3[m]) > max_diff_c3);
				max_diff_c3 = fabs(c0[m]-c2[m]);
		}
		printf("\nmax difference in c1: %d", max_diff_c1);
		printf("\nmax difference in c2: %d", max_diff_c2);
		printf("\nmax difference in c3: %d\n", max_diff_c3);
		printf("------------------------------------------\n\n");
	}		
	
	return 0;
}
