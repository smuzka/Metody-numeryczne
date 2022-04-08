

#include <iostream>
#include "gsl/gsl_math.h"
#include "gsl/gsl_linalg.h"
#include <string>
#include <cstring>
#include <fstream>
#include <iomanip>

#define ROWS_COLS 5
#define NUMBER_OF_ITERATIONS 50


gsl_matrix* reading_matrix_from_file(std::string name, int size){
    gsl_matrix* Matrix = gsl_matrix_calloc(size, size);
    double buffer;
    double denominator;
    char division_extractor;

    std::ifstream file(name);

        for(int i = 0; i < size; i++){
            for(int j = 0; j < size; j++){
                file >> buffer;             
                if(file.peek() == '/'){
                    file >> division_extractor;
                    file >> denominator;
                    gsl_matrix_set(Matrix, i, j, buffer/denominator);  
                }
                else{
                    gsl_matrix_set(Matrix, i, j, buffer);  
                }             
            }
        }
    file.close();

    return Matrix;
}

gsl_vector* reading_vector_from_file(std::string name, int size){

    gsl_vector* Vector = gsl_vector_calloc(size);
    double buffer;
    std::ifstream file(name);

    for(int i = 0; i < size; i++){
        file >> buffer;
        gsl_vector_set(Vector, i, buffer);
    }
    file.close();

    return Vector;
}

void square_matrix_print(gsl_matrix* Matrix, int size){
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            std::cout << std::setprecision(6) << gsl_matrix_get(Matrix,i,j) << " ";
        }
        std::cout << std::endl;
    }
}

void vector_print(gsl_vector* vector, int size) {
    std::cout << std::endl;
    for (int i = 0; i < size; i++) {
        std::cout << std::setprecision(6) << gsl_vector_get(vector, i) << " ";
    }
    std::cout << std::endl;
}

bool is_matrix_dominant_diagonally(gsl_matrix* Matrix, int size){
    double diagonal, sum;

    for(int i = 0; i < size ; i++){
        sum = 0;
        diagonal = fabs(gsl_matrix_get(Matrix, i, i));
        for(int j = 0; j < size; j++){
            if(i != j){
                sum += fabs(gsl_matrix_get(Matrix, i, j));
            }
        }
        if(sum > diagonal){
            std::cout << "Matrix is not diagonally dominant" << std::endl;
            return false;
        }
    }

    std::cout << "Matrix is diagonally dominant" << std::endl;
    return true;

}

void jacobi_solver(gsl_matrix* Matrix, gsl_vector* vector, int size, gsl_vector* x){
    gsl_vector* xp1 = gsl_vector_calloc(size);

    for(int k = 0; k < NUMBER_OF_ITERATIONS; k++){
        if(k % 5 == 0){
            std::cout << "iteration: " << k;
            for(int i = 0; i < size; i++){
                std::cout << ", x[" << i << "]=" << std::setprecision(6) << gsl_vector_get(x, i); 
            }
            std::cout << std::endl;
        }

        for(int i = 0; i < size; i++){
            double sum = 0;
            for(int j = 0; j < size; j++){
                if(i != j){
                    sum += gsl_matrix_get(Matrix, i, j) * gsl_vector_get(x, j);
                }
            }
            gsl_vector_set(xp1, i, (gsl_vector_get(vector, i) - sum) / gsl_matrix_get(Matrix, i, i));
        }

        for(int i = 0; i < size; i++){
            gsl_vector_set(x, i, gsl_vector_get(xp1, i));
        }
    }
}

gsl_vector* jacobi_solver_without_printing(gsl_matrix* Matrix, gsl_vector* vector, int size){
    gsl_vector* x = gsl_vector_calloc(size);
    gsl_vector* xp1 = gsl_vector_calloc(size);

    for(int k = 0; k < NUMBER_OF_ITERATIONS; k++){
        for(int i = 0; i < size; i++){
            double sum = 0;
            for(int j = 0; j < size; j++){
                if(i != j){
                    sum += gsl_matrix_get(Matrix, i, j) * gsl_vector_get(x, j);
                }
            }
            gsl_vector_set(xp1, i, (gsl_vector_get(vector, i) - sum) / gsl_matrix_get(Matrix, i, i));
        }

        for(int i = 0; i < size; i++){
            gsl_vector_set(x, i, gsl_vector_get(xp1, i));
        }
    }

    return x;
}


gsl_matrix* matrix_inversion(gsl_matrix* Matrix, int size){

    gsl_matrix* inverted_matrix = gsl_matrix_calloc(size, size);
    gsl_vector* temp_vector = gsl_vector_calloc(size);

    for(int i = 0; i < size; i++){
        gsl_vector_set_all(temp_vector, 0);
        gsl_vector_set(temp_vector, i, 1);
        gsl_matrix_set_col(inverted_matrix, i, jacobi_solver_without_printing(Matrix, temp_vector, size));
    }

    return inverted_matrix;
}


int main(){


    gsl_matrix* first_matrix = gsl_matrix_calloc(ROWS_COLS, ROWS_COLS);

    first_matrix = reading_matrix_from_file("macierz" + std::to_string(ROWS_COLS) + ".txt", ROWS_COLS);
    std::cout << "Matrix " << ROWS_COLS << "x" << ROWS_COLS << ":" << std::endl;
    square_matrix_print(first_matrix, ROWS_COLS);
    if(!is_matrix_dominant_diagonally(first_matrix, ROWS_COLS)){
        return 1;
    }

    gsl_vector* first_vector = gsl_vector_calloc(ROWS_COLS);
    first_vector = reading_vector_from_file("wektor" + std::to_string(ROWS_COLS) + ".txt", ROWS_COLS);
    std::cout << "Vector:";
    vector_print(first_vector, ROWS_COLS);

    // solution taken from online equation system solver https://matrixcalc.org
    std::cout << "expected solution: x[0]=" << 7763./60924 << ", x[1]=" << 7781./60924 << ", x[2]=" << 3583./20308 
            << ", x[3]=" << 11809./60924 << ", x[4]=" << 33445./60924 << std::endl;

    std::cout << std::endl;
    gsl_vector* first_x = gsl_vector_calloc(ROWS_COLS);
    jacobi_solver(first_matrix, first_vector, ROWS_COLS, first_x);

    std::cout << "result:";
    vector_print(first_x, ROWS_COLS);


    std::cout << std::endl;
    #undef ROWS_COLS
    #define ROWS_COLS 10

    std::cout << ROWS_COLS << std::endl;
    gsl_matrix* second_matrix = gsl_matrix_calloc(ROWS_COLS, ROWS_COLS);
    second_matrix = reading_matrix_from_file("macierz" + std::to_string(ROWS_COLS) + ".txt", ROWS_COLS);
    std::cout << "Matrix " << ROWS_COLS << "x" << ROWS_COLS << ":" << std::endl;
    square_matrix_print(second_matrix, ROWS_COLS);
        if(!is_matrix_dominant_diagonally(second_matrix, ROWS_COLS)){
        return 1;
    }

    gsl_vector* second_vector = gsl_vector_calloc(ROWS_COLS);
    second_vector = reading_vector_from_file("wektor" + std::to_string(ROWS_COLS) + ".txt", ROWS_COLS);
    std::cout << "Vector:";
    vector_print(second_vector, ROWS_COLS);

    // solution taken from online equation system solver https://matrixcalc.org
    std:: cout << "expected solution: x[0]=" << 2908327292459./225249041913918 << ", x[1]=" << -392862387934./10238592814269
            << ", x[2]=" << 20897385315557./225249041913918 << ", x[3]=" << 6625654370648./112624520956959 << ", x[4]="
            << 33308105940101./225249041913918 << ", x[5]=" << 38838301167053./225249041913918 << ", x[6]="
            << 1897711064035./20477185628538 << ", x[7]=" << 22354217271387./75083013971306 << ", x[8]="
            << 35867487982519./112624520956959 << ", x[9]=" << 2045253258869./6825728542846 << std::endl;

    std::cout << std::endl;
    gsl_vector* second_x = gsl_vector_calloc(ROWS_COLS);
    jacobi_solver(second_matrix, second_vector, ROWS_COLS, second_x);

    std::cout << "result:";
    vector_print(second_x, ROWS_COLS);


    // 3)
    std::cout << std::endl;
    #undef ROWS_COLS
    #define ROWS_COLS 5

    std::cout << "Matrix " << ROWS_COLS << "x" << ROWS_COLS << ":" << std::endl;
    square_matrix_print(first_matrix, ROWS_COLS);
    std::cout << "Vector:";
    vector_print(first_vector, ROWS_COLS);
    std::cout << std::endl;
    gsl_vector* x = gsl_vector_calloc(ROWS_COLS);
    for(int i = 0; i < ROWS_COLS; i++){
        gsl_vector_set(x, i, i*10);
    }
    std::cout << "starting vector x:";
    vector_print(x, ROWS_COLS);
    jacobi_solver(first_matrix, first_vector, ROWS_COLS, x);
    std::cout << "result:";
    vector_print(x, ROWS_COLS);

    for(int i = 0; i < ROWS_COLS; i++){
        gsl_vector_set(x, i, i*100);
    }
    std::cout << "starting vector x:";
    vector_print(x, ROWS_COLS);
    jacobi_solver(first_matrix, first_vector, ROWS_COLS, x);
    std::cout << "result:";
    vector_print(x, ROWS_COLS);

    for(int i = 0; i < ROWS_COLS; i++){
        gsl_vector_set(x, i, i*1000);
    }
    std::cout << "starting vector x:";
    vector_print(x, ROWS_COLS);
    jacobi_solver(first_matrix, first_vector, ROWS_COLS, x);
    std::cout << "result:";
    vector_print(x, ROWS_COLS);


    // 4)
        std::cout << std::endl;
    #undef ROWS_COLS
    #define ROWS_COLS 6

    gsl_matrix* new_matrix = gsl_matrix_calloc(ROWS_COLS, ROWS_COLS);

    new_matrix = reading_matrix_from_file("macierz" + std::to_string(ROWS_COLS) + ".txt", ROWS_COLS);
    std::cout << "Matrix " << ROWS_COLS << "x" << ROWS_COLS << ":" << std::endl;
    square_matrix_print(new_matrix, ROWS_COLS);

    gsl_vector* e = gsl_vector_calloc(ROWS_COLS);
    gsl_matrix* inverted_matrix = matrix_inversion(new_matrix, ROWS_COLS);
    std::cout << std::endl;
    std::cout << "inverted matrix:" << std::endl;
    square_matrix_print(inverted_matrix, ROWS_COLS);

    std::cout << std::endl;
    std::cout << "result:" << std::endl;
    gsl_matrix_mul_elements(new_matrix, inverted_matrix);
    square_matrix_print(new_matrix, ROWS_COLS);

    return 0;
}