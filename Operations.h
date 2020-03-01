//
// Created by alex on 20.02.2020.
//

#ifndef LAB1_OPERATIONS_H
#define LAB1_OPERATIONS_H

#include <stdbool.h>

void fill_matr(const int N, int const rank, int const size, double **matr);

void preset_solution(double **matr, double **answ, double **vector, int const N, int const rank, int const size);

double * make_new_vec(int const N);

double *
mul_vec_by_matr(double const *const vector, double const *const matrix, int const N, int const rank, int const size);

void subtract_vectors(double * const minuend_vec, double const * const sub_vec, int const N);

void mul_vec_by_numb(double * const vector, double const number, int const N);

void approximate(double *vector, double const *const matrix, double const *const answers, int const N, double const t,
                 int const rank, int const size);

double calc_norm(double const *const vec, int const N, int const size);

bool is_answer_correct(double const *const vec, double const *const matr, double const *const answ, int const N,
                       int const rank, int const size);

#endif //LAB1_OPERATIONS_H
