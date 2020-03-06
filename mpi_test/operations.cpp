#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <malloc.h>
#include <mpi.h>
#include <iostream>

#include "operations.h"

void fill_matr(const int N, int const rank, int const size, double **matr)
{
	*matr = make_new_vec(N * (N * size));

	FILE *f;
	char filename[32];
	snprintf(filename, sizeof(char) * 32, "matr%i.txt", rank);
	f = fopen(filename, "w");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N * size; j++) {
			if (i + (rank * N) == j) {
				(*matr)[i * N + j] = 1000;
			}
			else {
				(*matr)[i * N + j] = 1;
			}
			fprintf(f, "%f  ", (*matr)[i * N + j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void preset_solution(double **matr, double **answ, double **vector, int const N, int const rank, int const size)
{
	*vector = make_new_vec(N);
	for (int i = 0; i < N; i++) {
		(*vector)[i] = 0;
	}

	fill_matr(N, rank, size, matr);
	*answ = make_new_vec(N);
	FILE *f;
	char filename[32];
	snprintf(filename, sizeof(char) * 32, "leftPart%i.txt", rank);
	f = fopen(filename, "w");
	for (int i = 0; i < N; i++) {
		(*answ)[i] = N * size + 1;
		fprintf(f, "%f \n", (*answ)[i]);
	}
	fclose(f);
}

double *make_new_vec(int const N)
{
	return (double *)malloc(N * sizeof(double));
}

void collect_whole_vec(double *const whole_vec, double const *const part_of_vec, int const N, int const rank, int const size, char const * const mark, int const iteration)
{
	FILE *g;
	char filename2[32];
	snprintf(filename2, sizeof(char) * 32, "imm_part_vec%i.txt", rank);
	g = fopen(filename2, "a");
	fprintf(g, "%s: %d", mark, iteration);
	for (int i = 0; i < N; i++) {
		fprintf(g, "%f \n", part_of_vec[i]);
	}
	fprintf(g, "-----------------------------\n");
	fclose(g);

	MPI_Allgather(part_of_vec, N, MPI_DOUBLE, whole_vec, N, MPI_DOUBLE, MPI_COMM_WORLD);

	FILE *f;
	char filename[32];
	snprintf(filename, sizeof(char) * 32, "imm_whole_vec%i.txt", rank);
	f = fopen(filename, "a");
	fprintf(f, "%s: %d", mark, iteration);
	for (int i = 0; i < N*size; i++) {
		fprintf(f, "%f \n", whole_vec[i]);
	}
	fprintf(f, "-----------------------------\n");
	fclose(f);
}

double *
mul_vec_by_matr(double const *const vector, double const *const matrix, int const N, int const rank, int const size, int const iteration)
{
	double *new_vector = make_new_vec(N);
	double *whole_vec = make_new_vec(N*size);
	collect_whole_vec(whole_vec, vector, N, rank, size, "---------MULT-----------\n", iteration);
	for (int i = 0; i < N; i++) {
		double sum = 0;
		for (int j = 0; j < N * size; j++) {
			sum += matrix[i * N + j] * whole_vec[j];
		}
		new_vector[i] = sum;
	}
	free(whole_vec);
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

void approximate(double *vector, double const *const matrix, double const *const answers, int const N,
	double const t,
	int const rank, int const size, int const iteration)
{
	double *buff;
	buff = mul_vec_by_matr(vector, matrix, N, rank, size, iteration);
	subtract_vectors(buff, answers, N);
	mul_vec_by_numb(buff, t, N);
	subtract_vectors(vector, buff, N);
	free(buff);
}

double calc_norm(double const *const vec, int const N, int const rank, int const size, int const iteration)
{
	double norm = 0;
	double * const whole_vec = make_new_vec(N*size);

	collect_whole_vec(whole_vec, vec, N, rank, size, "---------NORM--------\n", iteration);

	for (int i = 0; i < N*size; i++) {
		norm += (whole_vec[i] * whole_vec[i]);
	}
	free(whole_vec);

	return sqrt(norm);
}

bool is_answer_correct(double const *const vec, double const *const matr, double const *const answ,
	int const N,
	int const rank, int const size, int const iteration)
{
	double const static E = 1e-1;
	double answ_norm = calc_norm(answ, N, rank, size, iteration);
	double *const buff = mul_vec_by_matr(vec, matr, N, rank, size, iteration);
	subtract_vectors(buff, answ, N);
	double norm = calc_norm(buff, N, rank, size, iteration) / answ_norm;
	FILE *f;
	char filename[32];
	snprintf(filename, sizeof(char) * 32, "deviation%i.txt", rank);
	f = fopen(filename, "a");
	fprintf(f, "deviation is : %f \n", norm);
	fprintf(f, "\n");
	fclose(f);
	bool is_approx_correct = norm < E;
	free(buff);
	return is_approx_correct;
}