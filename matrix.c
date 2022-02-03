//
// Created by Никита Иксарица on 30.01.2022.
//

#include "matrix.h"

size_t get_dynamic_matrix_size(int rows, int cols, size_t size_of_type, size_t size_of_type_pointer) {
    return rows * size_of_type_pointer + rows * cols * size_of_type;
}

void **alloc_matrix(int rows, int cols, size_t size_of_type, size_t size_of_type_pointer) {
    size_t size = get_dynamic_matrix_size(rows, cols, size_of_type, size_of_type_pointer);
    void **matrix = malloc(size);

    if (matrix != NULL) {
        // ptr is now pointing to the first element of 2D array.
        void *ptr = (void *) ((size_t) matrix + rows * size_of_type_pointer);

        // for-loop to point rows pointer to appropriate location in 2D array.
        for (int i = 0; i < rows; i++) {
            matrix[i] = (void *) ((size_t) ptr + i * cols * size_of_type);
        }
    }

    return matrix;
}

float **alloc_matrix_f(int rows, int cols) {
    return (float **) alloc_matrix(rows, cols, sizeof(float), sizeof(float *));
}

float **copy_matrix_f(int rows, int cols, float **matrix) {
    float **new_matrix = alloc_matrix_f(rows, cols);

    if (new_matrix != NULL) {
        size_t size = get_dynamic_matrix_size(rows, cols, sizeof(float), sizeof(float *));
        size_t offset = rows * sizeof(float *);
        memcpy((float **) ((size_t) new_matrix + offset),
               (float **) ((size_t) matrix + offset),
               size - offset);
    }

    return new_matrix;
}

void print_matrix_f(int rows, int cols, float **matrix) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%*.*f", 10, 2, matrix[i][j]);
        }

        printf("\n");
    }
}

float **random_matrix_init_f(int rows, int cols, float **matrix) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = (rand() % 201) - 100;
        }
    }

    return matrix;
}

float **set_every_matrix_el_to_val_f(int rows, int cols, float **matrix, float val) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = val;
        }
    }

    return matrix;
}
