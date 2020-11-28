#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 4

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

        printf("Matrix A after Swap:\n");
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                printf("%.2lf ", A[i][j]);
            }
            printf("\n");
        }

        for (j = i + 1; j < n; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < n; k++) A[j][k] -= A[j][i] * A[i][k];
        }

        printf("Matrix A after Muls:\n");
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                printf("%.2lf ", A[i][j]);
            }
            printf("\n");
        }
    }

    return 1;  // decomposition done
}

int main(int argc, char const *argv[]) {
    double **A = calloc(N, sizeof(double *));

    for (size_t i = 0; i < N; i++) A[i] = calloc(N, sizeof(double));

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            // A[i][j] = rand() % 10;
            scanf("%lf", &A[i][j]);
        }
    }

    printf("Matrix A before:\n");
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            printf("%.2lf ", A[i][j]);
        }
        printf("\n");
    }

    LUPDecompose(A, N, 0);

    printf("\nMatrix A after:\n");
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            printf("%.2lf ", A[i][j]);
        }
        printf("\n");
    }

    return 0;
}
