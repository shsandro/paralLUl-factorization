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

        if (maxA <= Tol) return 0;  // failure, matrix is degenerate

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
    printf("    -i <file_name>   sets an input file\n");
    printf("    -o <file_name>   sets an output file\n");
    printf(
        "    -n <number> generate random a matrix of size number x number\n");

    printf(
        "\nIt's important to note that -o option will only work if -v option "
        "is set\n");
}

int main(int argc, char const *argv[]) {
    int ch, N, n_set = 0, verbose = 0;
    FILE *output = stdout, *input = NULL;

    while ((ch = getopt(argc, argv, "i:o:n:hv")) != -1) {
        switch (ch) {
            case 'i':
                input = fopen(optarg, "r");
                break;

            case 'o':
                output = fopen(optarg, "w+");
                break;

            case 'n':
                N = atoi(optarg);
                n_set = 1;
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
        if (input != NULL) {
            printf("Reading inputs...\n");
            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                    fscanf(input, "%lf", &A[i * N + j]);
                }
            }
        } else {
            printf("Generating inputs...\n");
            srand(time(NULL));
            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                    A[i * N + j] = rand() % 100;
                }
            }
        }
    }

    upc_barrier;
    if (MYTHREAD == 0) {
        if (verbose) {
            fprintf(output, "\nOriginal matrix:\n");
            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                    fprintf(output, "% 6.2lf ", A[i * N + j]);
                }
                fprintf(output, "\n");
            }
            printf("\n");
        }
    }

    if (LUPDecompose(A, N, 0.0001) == 0) {
        printf("Erro matriz degenerada\n");
        upc_global_exit(1);
    }

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
    }

    return 0;
}
