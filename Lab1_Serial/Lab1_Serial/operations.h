#pragma once

void fill_matr(const int N, double **matr);

void preset_solution(double **matr, double **answ, double **vector, int const N);

double * make_new_vec(int const N);

double *
mul_vec_by_matr(double const *const vector, double const *const matrix, int const N, int const iteration);

void subtract_vectors(double * const minuend_vec, double const * const sub_vec, int const N);

void mul_vec_by_numb(double * const vector, double const number, int const N);

double squares_sum(double const * const vec, int const N);

double approximate(double *vector, double const *const matrix, double const *const answers,
 int const N,

	double const t,
	int const iteration);

double calc_norm(double * vect, int N);

bool is_answer_correct(double const sq_sum,

	double const *const answ,


	int const N,



	int const iteration);

void print_vec(char const * const name, double const * const vec, int const rank, int const iteration, char const * const mark, int const N, int const size);

