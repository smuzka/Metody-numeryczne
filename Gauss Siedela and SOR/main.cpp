

#include "gsl/gsl_math.h"
#include "gsl/gsl_linalg.h"
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <cstdlib>


#define ROWS_COLS 4
#define MAX_ITERATIONS 200

gsl_matrix* create_diagonally_dominant_matrix(int size){

    gsl_matrix* Matrix = gsl_matrix_calloc(size, size);
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            if(i != j){
                gsl_matrix_set(Matrix, i, j, rand() % 5);
            }
            else{
                gsl_matrix_set(Matrix, i, j, rand() % 100 + 100);
            }
        }
    }

    return Matrix;
}


void export_matrix_to_file(std::string name, int size, gsl_matrix *matrix){
    
    std::ofstream file(name);

    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){    
            file << gsl_matrix_get(matrix, i, j);
            file << " ";
        }
        file.seekp(-1 + file.tellp());
        file << std::endl;
    }
    file.close();
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

double vector_norm(gsl_vector *first, int size){

    double buffer = 0;

    for (int i = 0; i < size; i++){
        buffer += gsl_vector_get(first, i) * gsl_vector_get(first, i);
    }

    return sqrt(buffer);
}

double vector_norm(gsl_vector *first, gsl_vector *second, int size){

    gsl_vector* temp_vector = gsl_vector_calloc(size);

    for (int i = 0; i < size; i++){
        gsl_vector_set(temp_vector, i, gsl_vector_get(second, i) - gsl_vector_get(first, i));
    }

    return vector_norm(temp_vector, size);
}

gsl_vector* gauss_seidlea_solver(gsl_matrix* matrix, gsl_vector* vector, int size, double reltolerr){
    gsl_vector* first_x = gsl_vector_calloc(size);
    gsl_vector* previous_x = gsl_vector_calloc(size);

    for(int i = 0; i < MAX_ITERATIONS && vector_norm(previous_x, first_x, size) >= reltolerr * vector_norm(first_x, size); i++){

        for(int j = 0; j < size; j++){
            gsl_vector_set(previous_x, j, gsl_vector_get(first_x, j));
        }

        for(int j = 0; j < size; j++){
            double sum = 0;
            for(int k = 0; k < size; k++){
                if(j != k){
                    sum += gsl_matrix_get(matrix, j, k) * gsl_vector_get(first_x, k);
                }
            }
            gsl_vector_set(first_x, j, (gsl_vector_get(vector, j) - sum) / gsl_matrix_get(matrix, j, j)); 
        }

        if(true){
            std::cout << "iteracja: " << i;
            vector_print(first_x, size); 
        }
    }

    return first_x;
}

gsl_vector* SOR_method_solver(gsl_matrix *matrix, gsl_vector* vector, int size, double reltolerr){
    gsl_vector* first_x = gsl_vector_calloc(size);
    gsl_vector* previous_x = gsl_vector_calloc(size);

    double w = 0.8;

    for(int i = 0; i < MAX_ITERATIONS && vector_norm(previous_x, first_x, size) >= reltolerr * vector_norm(first_x, size); i++){

        for(int j = 0; j < size; j++){
            gsl_vector_set(previous_x, j, gsl_vector_get(first_x, j));
        }

        for(int j = 0; j < size; j++){
            double sum1 = 0;
            double sum2 = 0;
            for(int k = 0; k < size; k++){
                if(k < j){
                    sum1 += gsl_matrix_get(matrix, j, k) * gsl_vector_get(first_x, k);
                }
                if(k > j){
                    sum2 += gsl_matrix_get(matrix, j, k) * gsl_vector_get(first_x, k);
                }
            double sum2 = 0;
            }

        gsl_vector_set(first_x, j, w / gsl_matrix_get(matrix, j, j) * (gsl_vector_get(vector, j) - sum1 - sum2) + (1-w) * gsl_vector_get(first_x, j));
        }

        if(true){
            std::cout << "iteracja: " << i;
            vector_print(first_x, size); 
        }

    }

    return first_x;

}

int SOR_method_iteration_counter(gsl_matrix *matrix, gsl_vector* vector, int size, double reltolerr, double w){
    gsl_vector* first_x = gsl_vector_calloc(size);
    gsl_vector* previous_x = gsl_vector_calloc(size);

    int i;

    for(i = 0; i < MAX_ITERATIONS && vector_norm(previous_x, first_x, size) >= reltolerr * vector_norm(first_x, size); i++){

        for(int j = 0; j < size; j++){
            gsl_vector_set(previous_x, j, gsl_vector_get(first_x, j));
        }

        for(int j = 0; j < size; j++){
            double sum1 = 0;
            double sum2 = 0;
            for(int k = 0; k < size; k++){
                if(k < j){
                    sum1 += gsl_matrix_get(matrix, j, k) * gsl_vector_get(first_x, k);
                }
                if(k > j){
                    sum2 += gsl_matrix_get(matrix, j, k) * gsl_vector_get(first_x, k);
                }
            double sum2 = 0;
            }

        gsl_vector_set(first_x, j, w / gsl_matrix_get(matrix, j, j) * (gsl_vector_get(vector, j) - sum1 - sum2) + (1-w) * gsl_vector_get(first_x, j));
        }
    }

    return i;

}


int main(){

    // 4x4 Matrix

    gsl_matrix* matrix_from_exercise = gsl_matrix_calloc(4, 4);
    matrix_from_exercise = reading_matrix_from_file("macierz4.txt", 4);
    std::cout << "Matrix 4x4:" << std::endl;
    square_matrix_print(matrix_from_exercise, 4);
    if(!is_matrix_dominant_diagonally(matrix_from_exercise, 4)){
        return 1;
    }

    gsl_vector* vector_from_exercise = gsl_vector_calloc(4);
    vector_from_exercise = reading_vector_from_file("wektor4.txt", 4);
    std::cout << "Vector:";
    vector_print(vector_from_exercise, 4);


    // testing for different reltolerr

    std::cout << std::endl;
    gsl_vector *x_from_exercise = gsl_vector_calloc(4);
    std::cout << "reltolerr = 0.1" << std::endl;
    x_from_exercise = gauss_seidlea_solver(matrix_from_exercise, vector_from_exercise, 4, 0.1);
    std::cout << "\nreltolerr = 0.01" << std::endl;
    x_from_exercise = gauss_seidlea_solver(matrix_from_exercise, vector_from_exercise, 4, 0.01);
    std::cout << "\nreltolerr = 0.01" << std::endl;
    x_from_exercise = gauss_seidlea_solver(matrix_from_exercise, vector_from_exercise, 4, 0.001);
    std::cout << "\nreltolerr = 0.01" << std::endl;
    x_from_exercise = gauss_seidlea_solver(matrix_from_exercise, vector_from_exercise, 4, 0.0001);



    // 10x10 Matrix

    gsl_matrix *first_matrix = gsl_matrix_calloc(10, 10);
    first_matrix = reading_matrix_from_file("macierz10.txt", 10);
    std::cout << std::endl;
    std::cout << "Matrix 10x10:" << std::endl;
    square_matrix_print(first_matrix, 10);
    if(!is_matrix_dominant_diagonally(first_matrix, 10)){
        return 1;
    }

    gsl_vector* first_x = gsl_vector_calloc(10);
    for(int i = 0; i < 10; i++){
        gsl_vector_set(first_x, i, i);
    }
    std::cout << std::endl;
    std::cout << "x vector:";
    vector_print(first_x, 10);

    gsl_vector* first_solution_vector = gsl_vector_calloc(10);
    gsl_blas_dgemv(CblasNoTrans, 1.0, first_matrix, first_x, 0.0, first_solution_vector);
    std::cout << std::endl;
    std::cout << "first_solution_vector:";
    vector_print(first_solution_vector, 10);

    std::cout << "\nGauss Seidlea:" << std::endl;
    first_x = gauss_seidlea_solver(first_matrix, first_solution_vector, 10, 0.0001);
    std::cout << "\nSOR method:" << std::endl;
    first_x = SOR_method_solver(first_matrix, first_solution_vector, 10, 0.0001);



    // 20x20 matrix

    gsl_matrix *second_matrix = gsl_matrix_calloc(20, 20);
    second_matrix = reading_matrix_from_file("macierz20.txt", 20);
    std::cout << std::endl;
    std::cout << "Matrix 20x20:" << std::endl;
    square_matrix_print(second_matrix, 20);
    if(!is_matrix_dominant_diagonally(second_matrix, 20)){
        return 1;
    }

    gsl_vector* second_x = gsl_vector_calloc(20);
    for(int i = 0; i < 20; i++){
        gsl_vector_set(second_x, i, i);
    }
    std::cout << std::endl;
    std::cout << "Solution vector:";
    vector_print(second_x, 20);

    gsl_vector* second_solution_vector = gsl_vector_calloc(20);
    gsl_blas_dgemv(CblasNoTrans, 1.0, second_matrix, second_x, 0.0, second_solution_vector);
    std::cout << std::endl;
    std::cout << "solution_vector:";
    vector_print(second_solution_vector, 20);

    std::cout << "\nGauss Seidlea:" << std::endl;
    second_x = gauss_seidlea_solver(second_matrix, second_solution_vector, 20, 0.0001);
    std::cout << "\nSOR method:" << std::endl;
    second_x = SOR_method_solver(second_matrix, second_solution_vector, 20, 0.0001);



    // finding right omega for SOR method

    std::ofstream file("0.01iterationnumber");
    std::ofstream file2("0.01omega");
    std::cout << "\nreltolerr = 0.01";
    std::cout << std::endl;
    for(double omega = 0.2; omega < 1.9; omega += 0.05){
        std::cout << "relaxation parameter, omega = " << omega << ", iterations: " << SOR_method_iteration_counter(second_matrix, second_solution_vector, 20, 0.01, omega) << std::endl;
        file << SOR_method_iteration_counter(second_matrix, second_solution_vector, 20, 0.01, omega);
        file << ", ";
        file2 << omega;
        file2 << ", ";
    }
    file.close();
    file2.close();


    std::ofstream file3("0.0001iterationnumber");
    std::ofstream file4("0.0001omega");
    std::cout << "\nreltolerr = 0.0001";
    std::cout << std::endl;
    for(double omega = 0.2; omega < 1.9; omega += 0.05){
        std::cout << "relaxation parameter, omega = " << omega << ", iterations: " << SOR_method_iteration_counter(second_matrix, second_solution_vector, 20, 0.0001, omega) << std::endl;
        file3 << SOR_method_iteration_counter(second_matrix, second_solution_vector, 20, 0.0001, omega);
        file3 << ", ";
        file4 << omega;
        file4 << ", ";
    }
    file3.close();
    file4.close();


    return 0;
}