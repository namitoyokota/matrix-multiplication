#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double getRandom() {
        return (5-2)*rand()/(double)RAND_MAX+2;
}

void setMatrix(int n, double *a) {
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			a[i*j] = getRandom();
		}
	}
}

void simple_ijk (double *a, double *b, double *c, int n) {
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			register double r = c[i*n+j];
			for (int k=0; k<n; k++)
				r += a[i*n+k] * b[k*n+j];
			c[i*n+j] = r;
		}
	}
}

void simple_ikj (double *a, double *b, double *c, int n) {
        for (int i=0; i<n; i++) {
                for (int k=0; k<n; k++) {
                        register double r = a[i*n+k];
                        for (int j=0; j<n; j++)
				c[i*n+j] += r * b[k*n+j];
                }
        }
}

void simple_jik (double *a, double *b, double *c, int n) {
        for (int j=0; j<n; j++) {
                for (int i=0; i<n; i++) {
                        register double r = c[i*n+j];
                        for (int k=0; k<n; k++)
                                r += a[i*n+k] * b[k*n+j];
                        c[i*n+j] = r;
                }
        }
}

void simple_jki (double *a, double *b, double *c, int n) {
        for (int j=0; j<n; j++) {
                for (int k=0; k<n; k++) {
                        register double r = b[k*n+j];
                        for (int i=0; i<n; i++)
				c[i*n+j] += a[i*n+k] * r;
                }
        }
}

void simple_kij (double *a, double *b, double *c, int n) {
        for (int k=0; k<n; k++) {
                for (int i=0; i<n; i++) {
                        register double r = a[i*n+k];
                        for (int j=0; j<n; j++)
                                c[i*n+j] += r * b[k*n+j];
                }
        }
}

void simple_kji (double *a, double *b, double *c, int n) {
        for (int k=0; k<n; k++) {
                for (int j=0; j<n; j++) {
                        register double r = b[k*n+j];
                        for (int i=0; i<n; i++)
                                c[i*n+j] += a[i*n+k] * r;
                }
        }
}

void blocked_ijk (double *a, double *b, double *c, int n, int B) {
	for (int i=0; i<n; i+=B) {
		for (int j=0; j<n; j+=B) {
			for (int k=0; k<n; k+=B) {
				for (int i1=i; i1<i+B; i1++) {
					for (int j1=j; j1<j+B; j1++) {
						for (int k1=k; k1<k+B; k1++)
							c[i1*n+j1] += a[i1*n+k1]*b[k1*n+j1];
					}
				}
			}
		}
	}
}

void blocked_ikj (double *a, double *b, double *c, int n, int B) {
        for (int i=0; i<n; i+=B) {
                for (int k=0; k<n; k+=B) {
                        for (int j=0; j<n; j+=B) {
                                for (int i1=i; i1<B+i; i1++) {
                                        for (int k1=k; k1<B+k; k1++) {
                                                for (int j1=j; j1<B+j; j1++) {
                                                        c[i1*n+j1] += a[i1*n+k1] * b[k1*n+j1];
                                                }
                                        }
                                }
                        }
                }
        }
}

void blocked_jik (double *a, double *b, double *c, int n, int B) {
        for (int j=0; j<n; j+=B) {
                for (int i=0; i<n; i+=B) {
                        for (int k=0; k<n; k+=B) {
                                for (int i1=j; i1<B+j; i1++) {
                                        for (int j1=i; j1<B+i; j1++) {
                                                for (int k1=k; k1<B+k; k1++) {
                                                        c[j1*n+i1] += a[j1*n+k1] * b[k1*n+i1];
                                                }
                                        }
                                }
                        }
                }
        }
}

void blocked_jki (double *a, double *b, double *c, int n,  int B) {
        for (int j=0; j<n; j+=B) {
                for (int k=0; k<n; k+=B) {
                        for (int i=0; i<n; i+=B) {
                                for (int i1=j; i1<B+j; i1++) {
                                        for (int j1=k; j1<B+k; j1++) {
                                                for (int k1=i; k1<B+i; k1++) {
                                                        c[k1*n+i1] += a[k1*n+j1] * b[j1*n+i1];
                                                }
                                        }
                                }
                        }
                }
        }
}

void blocked_kij (double *a, double *b, double *c, int n, int B) {
        for (int k=0; k<n; k+=B) {
                for (int i=0; i<n; i+=B) {
                        for (int j=0; j<n; j+=B) {
                                for (int i1=k; i1<B+k; i1++) {
                                        for (int j1=i; j1<B+i; j1++) {
                                                for (int k1=j; k1<B+j; k1++) {
                                                        c[j1*n+k1] += a[j1*n+i1] * b[i1*n+k1];
                                                }
                                        }
                                }
                        }
                }
        }
}

void blocked_kji (double *a, double *b, double *c, int n, int B) { 
        for (int k=0; k < n; k+=B) {
                for (int j=0; j<n; j+=B) {
                        for (int i=0; i<n; i+=B) {
                                for (int i1=k; i1<B+k; i1++) {
                                        for (int j1=j; j1<B+j; j1++) {
                                                for (int k1=i; k1<B+i; k1++) {
                                                        c[k1*n+j1] += a[k1*n+i1] * b[i1*n+j1];
                                                }
                                        }
                                }
                        }
                }
        }
}

void hybrid_ijk (double *a, double *b, double *c, int n, int B) {
        for (int i=0; i<n; i+=B) {
                for (int j=0; j<n; j+=B) {
                        for (int k=0; k<n; k+=B) {
                                for (int i1=i; i1<B+i; i1+=2) {
                                        for (int j1=j; j1<B+j; j1+=2) {
                                                register int t = i1*n+j1;
                        			register int tt = t+n;

                       				register double c00 = c[t];
                        			register double c01 = c[t+1];
                        			register double c10 = c[tt];
                        			register double c11 = c[tt+1];

                        			for (int k1=k; k1<B+k; k1+=2) {

                                			register int ta = i1*n+k1;
                                			register int tta = ta+n;
                                			register int tb = k1*n+j1;
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
                                			c11 += a10*b01;
                        			}

                        			c[t] = c00;
                        			c[t+1] = c01;
                        			c[tt] = c10;
                        			c[tt+1] = c11;
                                        }
                                }
                        }
                }
        }
}

void hybrid_ikj (double *a, double *b, double *c, int n, int B) {
        for (int i=0; i<n; i+=B) {
                for (int k=0; k<n; k+=B) {
                        for (int j=0; j<n; j+=B) {
                                for (int i1=i; i1<B+i; i1+=2) {
                                        for (int k1=k; k1<B+k; k1+=2) {
                                                register int t = i1*n+k1;
                        			register int tt = t+n;

                       				register double c00 = c[t];
                        			register double c01 = c[t+1];
                        			register double c10 = c[tt];
                        			register double c11 = c[tt+1];

                        			for (int j1=j; j1<B+j; j1+=2) {

                                			register int ta = i1*n+j1;
                                			register int tta = ta+n;
                                			register int tb = j1*n+k1;
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
                                			c11 += a10*b01;
                        			}

                        			c[t] = c00;
                        			c[t+1] = c01;
                        			c[tt] = c10;
                        			c[tt+1] = c11;
                                        }
                                }
                        }
                }
        }
}

void hybrid_jik (double *a, double *b, double *c, int n, int B) {
        for (int j=0; j<n; j+=B) {
                for (int i=0; i<n; i+=B) {
                        for (int k=0; k<n; k+=B) {
                                for (int j1=j; j1<B+j; j1+=2) {
                                        for (int i1=i; i1<B+i; i1+=2) {
                                                register int t = j1*n+i1;
                        			register int tt = t+n;

                       				register double c00 = c[t];
                        			register double c01 = c[t+1];
                        			register double c10 = c[tt];
                        			register double c11 = c[tt+1];

                        			for (int k1=k; k1<B+k; k1+=2) {

                                			register int ta = j1*n+k1;
                                			register int tta = ta+n;
                                			register int tb = k1*n+i1;
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
                                			c11 += a10*b01;
                        			}

                        			c[t] = c00;
                        			c[t+1] = c01;
                        			c[tt] = c10;
                        			c[tt+1] = c11;
                                        }
                                }
                        }
                }
        }
}

void hybrid_jki (double *a, double *b, double *c, int n, int B) {
        for (int j=0; j<n; j+=B) {
                for (int k=0; k<n; k+=B) {
                        for (int i=0; i<n; i+=B) {
                                for (int j1=j; j1<B+j; j1+=2) {
                                        for (int k1=k; k1<B+k; k1+=2) {
                                                register int t = j1*n+k1;
                        			register int tt = t+n;

                       				register double c00 = c[t];
                        			register double c01 = c[t+1];
                        			register double c10 = c[tt];
                        			register double c11 = c[tt+1];

                        			for (int i1=i; i1<B+i; i1+=2) {

                                			register int ta = j1*n+i1;
                                			register int tta = ta+n;
                                			register int tb = i1*n+k1;
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
                                			c11 += a10*b01;
                        			}

                        			c[t] = c00;
                        			c[t+1] = c01;
                        			c[tt] = c10;
                        			c[tt+1] = c11;
                                        }
                                }
                        }
                }
        }
}

void hybrid_kij (double *a, double *b, double *c, int n, int B) {
        for (int k=0; k<n; k+=B) {
                for (int i=0; i<n; i+=B) {
                        for (int j=0; j<n; j+=B) {
                                for (int k1=k; k1<B+k; k1+=2) {
                                        for (int i1=i; i1<B+i; i1+=2) {
                                                register int t = k1*n+i1;
                        			register int tt = t+n;

                       				register double c00 = c[t];
                        			register double c01 = c[t+1];
                        			register double c10 = c[tt];
                        			register double c11 = c[tt+1];

                        			for (int j1=j; j1<B+j; j1+=2) {

                                			register int ta = k1*n+j1;
                                			register int tta = ta+n;
                                			register int tb = j1*n+i1;
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
                                			c11 += a10*b01;
                        			}

                        			c[t] = c00;
                        			c[t+1] = c01;
                        			c[tt] = c10;
                        			c[tt+1] = c11;
                                        }
                                }
                        }
                }
        }
}

void hybrid_kji (double *a, double *b, double *c, int n, int B) {
        for (int k=0; k<n; k+=B) {
                for (int j=0; j<n; j+=B) {
                        for (int i=0; i<n; i+=B) {
                                for (int k1=k; k1<B+k; k1+=2) {
                                        for (int j1=j; j1<B+j; j1+=2) {
                                                register int t = k1*n+j1;
                        			register int tt = t+n;

                       				register double c00 = c[t];
                        			register double c01 = c[t+1];
                        			register double c10 = c[tt];
                        			register double c11 = c[tt+1];

                        			for (int i1=i; i1<B+i; i1+=2) {

                                			register int ta = k1*n+i1;
                                			register int tta = ta+n;
                                			register int tb = i1*n+j1;
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
                                			c11 += a10*b01;
                        			}

                        			c[t] = c00;
                        			c[t+1] = c01;
                        			c[tt] = c10;
                        			c[tt+1] = c11;
                                        }
                                }
                        }
                }
        }
}

int main (int agrc, char *argv[]) {
        // matrix 
	int n = 2048;
        int block[6] = {16, 32, 64, 128, 256, 512};

        // max difference
        int max_diff_block = 0, max_diff_hybrid = 0;
        
        // functions
        char *functions[6] = {"ijk", "ikj", "jik", "jki", "kij", "kji"};
        void (*simple[6])(double*, double*, double*, int) = {simple_ijk, simple_ikj, simple_jik, simple_jki, simple_kij, simple_kji};
	void (*blocked[6])(double*, double*, double*, int, int) = {blocked_ijk, blocked_ikj, blocked_jik, blocked_jki, blocked_kij, blocked_kji};
	void (*hybrid[6])(double*, double*, double*, int, int) = {hybrid_ijk, hybrid_ikj, hybrid_jik, hybrid_jki, hybrid_kij, hybrid_kji};

        // time
	struct timespec start, end;
	double t;

        // a and b arrays
      	double *a = (double *)calloc(sizeof(double), n*n);
        double *b = (double *)calloc(sizeof(double), n*n);
        setMatrix(n, a);
        setMatrix(n, b);

        // loop through all functions: ijk, ikj, jik, jki, kij, kji
        for (int i=0; i<6; i++) {
                // reset max difference for each algorithm
                max_diff_block = 0; max_diff_hybrid = 0;

                // run single register reuse
                printf("\n-----------------------------\n");
                printf("Simple (%s)\n\n", functions[i]);
                double *c1 = (double *)calloc(sizeof(double), n*n);
                clock_gettime(CLOCK_REALTIME, &start);
                simple[i](a, b, c1, n);
                clock_gettime(CLOCK_REALTIME, &end);
                t = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;
                printf("\t%s: %lf seconds\n", functions[i], t);

                // run cache blocking
                printf("\n-----------------------------\n");
                printf("Blocked (%s)\n\n", functions[i]);
                // look through all block numbers: 10, 16, 32, 64, 128, 256, 512
                for (int j=0; j<6; j++) {
                        double *c2 = (double *)calloc(sizeof(double), n*n);
                        clock_gettime(CLOCK_REALTIME, &start);
                        blocked[i](a, b, c2, n, block[j]);
                        clock_gettime(CLOCK_REALTIME, &end);
                        t = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;
                        printf("\t%d: %lf seconds\n", block[j], t);
                        // calculate max difference
                        for (int m=0; m<n*n; m++) {
			        if (fabs(c1[m]-c2[m]) > max_diff_block)
				        max_diff_block = fabs(c1[m]-c2[m]);
		        }
                        printf("\tmax diff: %d\n", max_diff_block);
                }

                // run hybrid (register + cache blocking)
                printf("\n-----------------------------\n");
                printf("Hybrid (%s)\n\n", functions[i]);
                // loop through all block numbers: 16, 32, 64, 128, 256, 512
                for (int j=0; j<6; j++) {
                        double *c3 = (double *)calloc(sizeof(double), n*n);
                        clock_gettime(CLOCK_REALTIME, &start);
                        hybrid[i](a, b, c3, n, block[j]);
                        clock_gettime(CLOCK_REALTIME, &end);
                        t = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;
                        printf("\t%d: %lf seconds\n", block[j], t);
                        // calculate max difference
                        for (int m=0; m<n*n; m++) {
			        if (fabs(c1[m]-c3[m]) > max_diff_hybrid)
				        max_diff_hybrid = fabs(c1[m]-c3[m]);
		        }
                        printf("\tmax diff: %d\n", max_diff_hybrid);
                }
        }
        printf("\n-----------------------------\n\n");

        return 0;
}
