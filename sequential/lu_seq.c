#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int LUPDecompose(double **A, int n, double Tol) {
    int i, j, k, imax;
    double maxA, *ptr, absA;

    for (i = 0; i < n; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < n; k++)
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0;  // failure, matrix is degenerate

        if (imax != i) {
            // pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;
        }

        for (j = i + 1; j < n; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < n; k++) A[j][k] -= A[j][i] * A[i][k];
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

int main(int argc, char **argv) {
    int ch, N, n_set = 0, verbose = 0;
    double **A;
    FILE *input = NULL, *output = stdout;

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

    if (n_set) {
        double **A = calloc(N, sizeof(double *));

        for (size_t i = 0; i < N; i++) A[i] = calloc(N, sizeof(double));

        if (input != NULL) {
            printf("Reading inputs...\n");
            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                    fscanf(input, "%lf", &A[i][j]);
                }
            }
        } else {
            printf("Generating inputs...\n");
            srand(time(NULL));
            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                    A[i][j] = rand() % 100;
                }
            }
        }

        if (verbose) {
            fprintf(output, "\nOriginal matrix:\n");
            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                    fprintf(output, "% 6.2lf ", A[i][j]);
                }
                fprintf(output, "\n");
            }
            printf("\n");
        }

        printf("Calculating...\n");
        LUPDecompose(A, N, 0);
        printf("Done!\n");

        if (verbose) {
            printf("Writing output...\n");
            fprintf(output, "\nLU decomposed matrix:\n");
            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < N; j++) {
                    fprintf(output, "% 6.2lf ", A[i][j]);
                }
                fprintf(output, "\n");
            }
        }
    } else {
        printf("Required flag '-N' is not set! Aborting...\n\n");
        help_menu(argv[0]);
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
