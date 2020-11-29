#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <upc.h>
#include <upc_collective.h>

#define SYNC_MODE (UPC_IN_ALLSYNC | UPC_OUT_ALLSYNC)

struct pair {
    int first;
    double second;
};

shared struct pair values[THREADS];
shared struct pair max;

shared struct pair *cmp_max(shared struct pair *p1, shared struct pair *p2) {
    return (p2->second > p1->second) ? p2 : p1;
}

void reduce(shared struct pair *dst, shared struct pair *src, int src_tam,
            shared struct pair *(*f)(shared struct pair *,
                                     shared struct pair *)) {
    *dst = src[0];
    for (size_t i = 1; i < THREADS; i++) {
        *dst = *f(dst, &src[i]);
    }
}

int LUPDecompose(shared[] double *A, int n, double Tol) {
    int i, j, k, imax;
    double maxA, absA;

    for (i = 0; i < n; i++) {
        values[MYTHREAD].first = i;
        values[MYTHREAD].second = 0.0;

        upc_forall(k = i; k < n; k++; k) {
            if ((absA = fabs(A[k * n + i])) > values[MYTHREAD].second) {
                values[MYTHREAD].second = absA;
                values[MYTHREAD].first = k;
            }
        }

        upc_barrier;

        if (MYTHREAD == 0) {
            reduce(&max, values, THREADS, cmp_max);
        }
        upc_barrier;

        if (max.first != i) {
            // pivoting rows of A
            upc_forall(j = 0; j < n; j++; j) {
                double aux = A[i * n + j];
                A[i * n + j] = A[max.first * n + j];
                A[max.first * n + j] = aux;
            }
        }

        upc_barrier;

        // if (maxA <= Tol) return 0;  // failure, matrix is degenerate

        for (j = i + 1; j < n; j++) {
            if (MYTHREAD == 0) A[j * n + i] /= A[i * n + i];
            upc_barrier;

            upc_forall(k = i + 1; k < n; k++; k) {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];
            }

            upc_barrier;
        }
    }

    return 1;  // decomposition done
}

void help_menu(const char *prog_name) {
    printf("Usage: %s [flags]\n", prog_name);
    printf("    -h               prints this usage guide\n");
    printf(
        "    -v               prints the input and output matrix after the "
        "process\n");
    printf("    -o <file_name>   sets an output file\n");
    printf(
        "    -n <number> generate random a matrix of size number x number\n");
}

double now() {
    const double ONE_BILLION = 1000000000.0;
    struct timespec current_time;

    clock_gettime(CLOCK_REALTIME, &current_time);

    return current_time.tv_sec + (current_time.tv_nsec / ONE_BILLION);
}

int main(int argc, char **argv) {
    int N = 4, verbose = 0, ch;
    FILE *original_matrix_out = fopen("original-matrix-par.out", "w+"),
         *output = stdout;

    while ((ch = getopt(argc, argv, "o:n:hv")) != -1) {
        switch (ch) {
            case 'o':
                output = fopen(optarg, "w+");
                break;

            case 'n':
                N = atoi(optarg);
                break;

            case 'h':
                help_menu(argv[0]);
                exit(EXIT_SUCCESS);

            case 'v':
                verbose = 1;
                break;

            default:
                help_menu(argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    shared[] double *A =
        (shared[] double *)upc_all_alloc(THREADS, N * N * sizeof(double));

    if (MYTHREAD == 0) {
        srand(time(NULL));
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i * N + j] = rand() % 100;
            }
        }
    }

    upc_barrier;
    if (MYTHREAD == 0) {
        printf("\nGenerating inputs...\n");

        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                fprintf(original_matrix_out, "%.2lf ", A[i * N + j]);
            }
            fprintf(original_matrix_out, "\n");
        }

        printf("Original matrix written to original-matrix-par.out\n\n");
    }

    if (MYTHREAD == 0) printf("Calculating...\n");
    double start_time = now();
    LUPDecompose(A, N, 0);
    double end_time = now();
    if (MYTHREAD == 0) printf("Done!\n");

    if (MYTHREAD == 0) {
        if (verbose) {
            printf("Writing output...\n");
            fprintf(output, "\nLU decomposed matrix:\n");
            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                    fprintf(output, "% 6.2lf ", A[i * N + j]);
                }
                fprintf(output, "\n");
            }
        }

        printf("Time elapsed: %lf\n", end_time - start_time);
    }

    return 0;
}
