//
// Created by Никита Иксарица on 30.01.2022.
//

#ifndef LUDECOMPOSITION_MATRIX_H
#define LUDECOMPOSITION_MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

size_t get_dynamic_matrix_size(int rows, int cols, size_t size_of_type, size_t size_of_type_pointer);

/**
 * Выделяет память для матрицы любого типа, которую затем можно освободить одним вызовом free().
 *
 * @param rows Количество строк.
 * @param cols Количество столбцов.
 * @param size_of_type Размер типа матрицы.
 * @param size_of_type_pointer Размера указателя на тип матрицы.
 *
 * @return Возвращает указатель на выделенную область памяти.
 */
void **alloc_matrix(int rows, int cols, size_t size_of_type, size_t size_of_type_pointer);

float **alloc_matrix_f(int rows, int cols);

float **copy_matrix_f(int rows, int cols, float **matrix);

void print_matrix_f(int rows, int cols, float **matrix);

/**
 * Заполнаяет матрицу случайными числами в диапазоне от -100 до 100.
 * Ядро рандома должно быть задано извне.
 *
 * @param rows Количество строк.
 * @param cols Количество столбцов.
 * @param[in, out] matrix Предварительно выделенная матрица.
 *
 * @return Возвращает указатель на matrix, инициализированную случайными числами.
 * Новую матрицу не создаёт.
 */
float **random_matrix_init_f(int rows, int cols, float **matrix);

float **set_every_matrix_el_to_val_f(int rows, int cols, float **matrix, float val);

#endif //LUDECOMPOSITION_MATRIX_H
