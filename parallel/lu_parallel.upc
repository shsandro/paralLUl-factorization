#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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

int main(int argc, char const *argv[]) {
    FILE *out = fopen("teste.txt", "w+");

    srand(time(NULL));
    int N = 4;
    shared[] double *A =
        (shared[] double *)upc_all_alloc(THREADS, N * N * sizeof(double));

    if (MYTHREAD == 0) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                // scanf("%lf", &A[i * N + j]);
                A[i * N + j] = rand() % 10;
            }
        }
    }

    upc_barrier;
    if (MYTHREAD == 0) {
        printf("Matrix A before:\n");
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                fprintf(out, "%.2lf ", A[i * N + j]);
            }
            fprintf(out, "\n");
        }
    }

    // if (MYTHREAD == 0) {
    if (LUPDecompose(A, N, 0.0001) == 0) {
        printf("Erro matriz degenerada\n");
        upc_global_exit(1);
    }
    // }

    if (MYTHREAD == 0) {
        printf("Matrix A in LU:\n");
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                printf("%.2lf ", A[i * N + j]);
            }
            printf("\n");
        }
    }

    return 0;
}
