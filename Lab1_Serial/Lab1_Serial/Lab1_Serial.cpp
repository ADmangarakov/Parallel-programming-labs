// mpi_test.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include "Operations.h"

int main(int argc, char *argv[])
{
	int const N = 12;
	//gettimeofday(&tv1, NULL);
	double *matr = NULL, *answ = NULL, *vector = NULL;
	preset_solution(&matr, &answ, &vector, N);
	double const T = 1e-5;
	int count = 0;
	double sq_norm;
	do {
		sq_norm = approximate(vector, matr, answ, N, T, count);
	} while (!is_answer_correct(sq_norm, answ, N, count++));

	//gettimeofday(&tv2, NULL);
	//double dt_sec = (tv2.tv_sec - tv1.tv_sec);
	//double dt_usec = (tv2.tv_usec - tv1.tv_usec);
	//double dt = dt_sec + 1e-6 * dt_usec;

	print_vec("answerName", vector, 1, 1, "Answer\n", N, 1);

	free(matr);
	free(answ);
	free(vector);
	return 0;
}
