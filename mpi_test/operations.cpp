#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <malloc.h>
#include <mpi.h>
#include <iostream>

#include "operations.h"

void fill_matr(const int N, int const rank, int const size, double **matr)
{
	*matr = make_new_vec(N * (N * size));

	//FILE *f;
	//char filename[32];
	//snprintf(filename, sizeof(char) * 32, "matr%i.txt", rank);
	//f = fopen(filename, "w");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N * size; j++) {
			if (i + (rank * N) == j) {
				(*matr)[i * N * size + j] = 1000;
			}
			else {
				(*matr)[i * N * size + j] = 1;
			}
			//fprintf(f, "%f  ", (*matr)[i * N * size + j]);
		}
		//fprintf(f, "\n");
	}
	//fclose(f);
}

void preset_solution(double **matr, double **answ, double **vector, int const N, int const rank, int const size)
{
	*vector = make_new_vec(N);
	for (int i = 0; i < N; i++) {
		(*vector)[i] = 0;
	}

	fill_matr(N, rank, size, matr);
	*answ = make_new_vec(N);
	//FILE *f;
	//char filename[32];
	//snprintf(filename, sizeof(char) * 32, "leftPart%i.txt", rank);
	//f = fopen(filename, "w");
	for (int i = 0; i < N; i++) {
		(*answ)[i] = N * size + 1;
		//fprintf(f, "%f \n", (*answ)[i]);
	}
	//fclose(f);
}

double *make_new_vec(int const N)
{
	return (double *)malloc(N * sizeof(double));
}

void collect_whole_vec(double *const whole_vec, double const *const part_of_vec, int const N, int const rank,
	int const size, char const * const mark, int const iteration)
{
	//print_vec("imm_part_vec", part_of_vec, rank, iteration, mark, N, 1);
	MPI_Allgather(part_of_vec, N, MPI_DOUBLE, whole_vec, N, MPI_DOUBLE, MPI_COMM_WORLD);

}

double *
mul_vec_by_matr(double const *const vector, double const *const matrix, int const N, int const rank,
	int const size, int const iteration)
{
	double *new_vector = make_new_vec(N);
	double *whole_vec = make_new_vec(N*size);
	collect_whole_vec(whole_vec, vector, N, rank, size, "---------MULT-----------\n", iteration);
	//print_vec("imm_whole_vec", whole_vec, rank, iteration, "---------MULT-----------\n", N, size);
	for (int i = 0; i < N; i++) {
		double sum = 0;
		//print_vec("matr_while_mult", matrix + i * N * size, rank, iteration, "---------MULT-----------\n", N, size);
		for (int j = 0; j < N * size; j++) {
			sum = sum + (matrix[i * N * size + j] * whole_vec[j]);
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

double squares_sum(double const * const vec, int const N)
{
	double norm = 0;
	for (int i = 0; i < N; i++) {
		norm += vec[i] * vec[i];
	}
	return norm;
}

double approximate(double *vector, double const *const matrix, double const *const answers, int const N,
	double const t,
	int const rank, int const size, int const iteration)
{
	double *buff;
	//print_vec("input_vec", vector, rank, iteration, "-------NEXT------\n", N, 1);//size = 1
	buff = mul_vec_by_matr(vector, matrix, N, rank, size, iteration);
	//print_vec("after_mult_matr", buff, rank, iteration, "-------NEXT------\n", N, 1);//size = 1
	subtract_vectors(buff, answers, N);
	//print_vec("after_sub_ans", buff, rank, iteration, "-------NEXT------\n", N, 1);//size = 1
	double sq_sum = squares_sum(buff, N);
	mul_vec_by_numb(buff, t, N);
	subtract_vectors(vector, buff, N);
	//print_vec("after_sub_buf", vector, rank, iteration, "-------NEXT------\n", N, 1);//size = 1
	free(buff);
	return sq_sum;
}

double calc_norm(double const sq_sum)
{
	double norm = 0;
	MPI_Allreduce(&sq_sum, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sqrt(norm);
}

bool is_answer_correct(double const sq_sum,
	double const *const answ,
	int const N, int const rank, int const size, int const iteration)
{
	double const static E = 1e-2;

	double answ_norm = 0;
	for (int i = 0; i < N*size; i++) {
		answ_norm += answ[0] * answ[0];
	}
	answ_norm = sqrt(answ_norm);
	double norm = calc_norm(sq_sum) / answ_norm;

	FILE *f;
	char filename[32];
	snprintf(filename, sizeof(char) * 32, "deviation%i.txt", rank);
	f = fopen(filename, "a");
	fprintf(f, "deviation#%d is : %f \n", iteration, norm);
	fprintf(f, "\n");
	fclose(f);

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