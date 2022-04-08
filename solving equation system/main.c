

#include <stdio.h>
#include "gsl/gsl_math.h"
#include "gsl/gsl_linalg.h"

void reading_matrix_from_file(gsl_matrix* Matrix, FILE* file_pointer, int rows_and_colums, double buffer) {
    if ((file_pointer = fopen("macierz4.txt", "rt")) == NULL) {
        printf("file couldn't be opened");
    }
    else {
        for (int i = 0; i < rows_and_colums; i++) {
            for (int j = 0; j < rows_and_colums; j++) {
                fscanf(file_pointer, "%lf", &buffer);
                gsl_matrix_set(Matrix, i, j, buffer);
            }
        }
        fclose(file_pointer);
    }
}

void square_matrix_print(gsl_matrix* Matrix, int rows_and_columns) {
    printf("\n");
    for (int i = 0; i < rows_and_columns; i++) {
        for (int j = 0; j < rows_and_columns; j++) {
            printf("%.2f", gsl_matrix_get(Matrix, i, j));
            printf(" ");
        }
        printf("\n");
    }
}

void vector_print(gsl_vector* vector, int lenght) {
    printf("\n");
    for (int i = 0; i < lenght; i++) {
        printf("%.2f ", gsl_vector_get(vector, i));
    }
    printf("\n");
}

int main() {

    int n = 4;
    double val = 0;
    gsl_matrix* A = gsl_matrix_calloc(n, n);
    FILE* fp = NULL;
    reading_matrix_from_file(A, fp, n, val);
    square_matrix_print(A, n);

    gsl_vector* b = gsl_vector_calloc(n);
    for (int i = 0; i < n; i++) {
        gsl_vector_set(b, i, i);
    }

    gsl_permutation* p = gsl_permutation_calloc(n);
    int s;
    gsl_linalg_LU_decomp(A, p, &s);

    gsl_vector* b1 = gsl_vector_calloc(n);
    for (int i = 0; i < n; i++) {
        gsl_vector_set(b1, i, gsl_vector_get(b, gsl_permutation_get(p, i)));
    }

    gsl_vector* y = gsl_vector_calloc(n);

    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i-1; j++) {
            if (j < i) {
                sum += gsl_matrix_get(A, i, j);
                sum *= gsl_vector_get(y, j);
            }
        }
        gsl_vector_set(y, i, (gsl_vector_get(b1, i) - sum)/gsl_matrix_get(A, i, i));
    }

    
    vector_print(y, n);

    gsl_vector* x = gsl_vector_calloc(n);

    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            if (j >= i) {
                sum += gsl_matrix_get(A, i, j);
                sum *= gsl_vector_get(x, j);
            }
        }
        gsl_vector_set(x, i, (gsl_vector_get(y, i) - sum) / gsl_matrix_get(A, i, i));
    }

    vector_print(x, n);

    gsl_linalg_LU_solve(A, p, b, x);
    vector_print(x, n);

    //cleaning
    gsl_vector_free(b);
    gsl_vector_free(b1);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_matrix_free(A);
    gsl_permutation_free(p);

    return 0;

}

//projekt 1
//rozwiazywanie ukladu rownan liniowych
// PA = LU - decompost
// Ax = b / *P
// PAx = Pb zmieniamy kolejnosc wierszy, Pb = b'
// LUx = b'

// 1) Ly = b -> y
// 2) Ux = y -> x

// dla macierzy n = 10
// w sprawozdaniu lda n = 4



