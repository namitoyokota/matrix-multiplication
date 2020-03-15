#include "helper.h"
#include "simple.h"
#include "blocked.h"
#include "hybrid.h"

int main(int agrc, char *argv[])
{
	// matrix
	int n = 10;
	int block[6] = {16, 32, 64, 128, 256, 512};

	// max difference
	int max_diff_block = 0, max_diff_hybrid = 0;

	// functions
	char *functions[6] = {"ijk", "ikj", "jik", "jki", "kij", "kji"};
	void (*simple[6])(double *, double *, double *, int) = {simple_ijk, simple_ikj, simple_jik, simple_jki, simple_kij, simple_kji};
	void (*blocked[6])(double *, double *, double *, int, int) = {blocked_ijk, blocked_ikj, blocked_jik, blocked_jki, blocked_kij, blocked_kji};
	void (*hybrid[6])(double *, double *, double *, int, int) = {hybrid_ijk, hybrid_ikj, hybrid_jik, hybrid_jki, hybrid_kij, hybrid_kji};

	// time
	struct timespec start, end;
	double t;

	// a and b arrays
	double *a = (double *)calloc(sizeof(double), n * n);
	double *b = (double *)calloc(sizeof(double), n * n);
	setMatrix(n, a);
	setMatrix(n, b);

	// loop through all functions: ijk, ikj, jik, jki, kij, kji
	for (int i = 0; i < 6; i++)
	{
		// reset max difference for each algorithm
		max_diff_block = 0;
		max_diff_hybrid = 0;

		// run single register reuse
		printf("\n-----------------------------\n");
		printf("Simple (%s)\n\n", functions[i]);
		double *c1 = (double *)calloc(sizeof(double), n * n);
		clock_gettime(CLOCK_REALTIME, &start);
		simple[i](a, b, c1, n);
		clock_gettime(CLOCK_REALTIME, &end);
		t = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf("\t%s: %lf seconds\n", functions[i], t);

		// run cache blocking
		printf("\n-----------------------------\n");
		printf("Blocked (%s)\n\n", functions[i]);
		// look through all block numbers: 10, 16, 32, 64, 128, 256, 512
		for (int j = 0; j < 6; j++)
		{
			double *c2 = (double *)calloc(sizeof(double), n * n);
			clock_gettime(CLOCK_REALTIME, &start);
			blocked[i](a, b, c2, n, block[j]);
			clock_gettime(CLOCK_REALTIME, &end);
			t = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;
			printf("\t%d: %lf seconds\n", block[j], t);
			// calculate max difference
			for (int m = 0; m < n * n; m++)
			{
				if (fabs(c1[m] - c2[m]) > max_diff_block)
					max_diff_block = fabs(c1[m] - c2[m]);
			}
			printf("\tmax diff: %d\n", max_diff_block);
		}

		// run hybrid (register + cache blocking)
		printf("\n-----------------------------\n");
		printf("Hybrid (%s)\n\n", functions[i]);
		// loop through all block numbers: 16, 32, 64, 128, 256, 512
		for (int j = 0; j < 6; j++)
		{
			double *c3 = (double *)calloc(sizeof(double), n * n);
			clock_gettime(CLOCK_REALTIME, &start);
			hybrid[i](a, b, c3, n, block[j]);
			clock_gettime(CLOCK_REALTIME, &end);
			t = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;
			printf("\t%d: %lf seconds\n", block[j], t);
			// calculate max difference
			for (int m = 0; m < n * n; m++)
			{
				if (fabs(c1[m] - c3[m]) > max_diff_hybrid)
					max_diff_hybrid = fabs(c1[m] - c3[m]);
			}
			printf("\tmax diff: %d\n", max_diff_hybrid);
		}
	}
	printf("\n-----------------------------\n\n");

	return 0;
}
