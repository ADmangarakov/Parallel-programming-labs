#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <malloc.h>
#include <iostream>

#include "operations.h"

void fill_matr(const int N, double **matr)
{
	*matr = make_new_vec(N * N);

	//FILE *f;
	//char filename[32];
	//snprintf(filename, sizeof(char) * 32, "matr%i.txt", rank);
	//f = fopen(filename, "w");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i + N == j) {
				(*matr)[i * N + j] = 1000;
			}
			else {
				(*matr)[i * N + j] = 1;
			}
			//fprintf(f, "%f  ", (*matr)[i * N * size + j]);
		}
		//fprintf(f, "\n");
	}
	//fclose(f);
}

void preset_solution(double **matr, double **answ, double **vector, int const N)
{
	*vector = make_new_vec(N);
	for (int i = 0; i < N; i++) {
		(*vector)[i] = 0;
	}

	fill_matr(N, matr);
	*answ = make_new_vec(N);
	//FILE *f;
	//char filename[32];
	//snprintf(filename, sizeof(char) * 32, "leftPart%i.txt", rank);
	//f = fopen(filename, "w");
	for (int i = 0; i < N; i++) {
		(*answ)[i] = N + 1;
		//fprintf(f, "%f \n", (*answ)[i]);
	}
	//fclose(f);
}

double *make_new_vec(int const N)
{
	return (double *)malloc(N * sizeof(double));
}

double *mul_vec_by_matr(double const *const vector, double const *const matrix, int const N,
	int const iteration)
{
	double *new_vector = make_new_vec(N);
	for (int i = 0; i < N; i++) {
		double sum = 0;
		for (int j = 0; j < N; j++) {
			sum = sum + (matrix[i * N + j] * vector[j]);
		}
		new_vector[i] = sum;
	}
	return new_vector;
}

void subtract_vectors(double *const minuend_vec, double const *const sub_vec, int const N)
{
	for (int i = 0; i < N; i++) {
		minuend_vec[i] -= sub_vec[i];
	}
}

void mul_vec_by_numb(double *const vector, double const number, int const N)
{
	for (int i = 0; i < N; i++) {
		vector[i] *= number;
	}
}

double squares_sum(double const * const vec, int const N)
{
	double norm = 0;
	for (int i = 0; i < N; i++) {
		norm += vec[i] * vec[i];
	}
	return norm;
}

double approximate(double *vector, double const *const matrix, double const *const answers,
	int const N,
	double const t,
	int const iteration)
{
	double *buff;
	
	buff = mul_vec_by_matr(vector, matrix, N, iteration);
	
	subtract_vectors(buff, answers, N);
	
	double sq_sum = squares_sum(buff, N);
	mul_vec_by_numb(buff, t, N);
	subtract_vectors(vector, buff, N);
	
	free(buff);
	return sq_sum;
}

double calc_norm(double const * const vect, int const N)
{
	double norm = 0;
	for (int i = 0; i < N; i++) {
		norm += (vect[i] * vect[i]);
	}
	return sqrt(norm);
}

bool is_answer_correct(double const sq_sum,	double const *const answ, int const N, int const iteration)
{
	double const static E = 1e-2;
	double answ_norm = 0;
	for (int i = 0; i < N; i++) {
		answ_norm += answ[0] * answ[0];
	}
	answ_norm = sqrt(answ_norm);
	double norm = sqrt(sq_sum) / answ_norm;

	//FILE *f;
	//char filename[32];
	//snprintf(filename, sizeof(char) * 32, "deviation%i.txt", rank);
	//f = fopen(filename, "a");
	//fprintf(f, "deviation#%d is : %f \n", iteration, norm);
	//fprintf(f, "\n");
	//fclose(f);

	bool is_approx_correct = norm < E;
	return is_approx_correct;
}

void print_vec(char const * const name, double const * const vec, int const rank, int const iteration,
	char const * const mark, int const N, int const size)
{
	FILE *f;
	char filename[32];
	snprintf(filename, sizeof(char) * 32, "%s%i.txt", name, rank);
	f = fopen(filename, "a");
	fprintf(f, "iteration_%d %s", iteration, mark);
	for (int i = 0; i < N*size; i++) {
		fprintf(f, "%f \n", vec[i]);
	}
	fprintf(f, "-----------------------------\n");
	fclose(f);
}