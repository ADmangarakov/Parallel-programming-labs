// mpi_test.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <malloc.h>
#include "Operations.h"

int main(int argc, char *argv[])
{
	//struct timeval tv1, tv2;

	int size, rank;
	MPI_Init(&argc, &argv); // Инициализация MPI
	MPI_Comm_size(MPI_COMM_WORLD, &size); // Получение числа процессов
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Получение номера процесса

	int const PROCESS_NUMBER = size;
	int const N = 4 / PROCESS_NUMBER;

	//gettimeofday(&tv1, NULL);
	double *matr = NULL, *answ = NULL, *vector = NULL;
	preset_solution(&matr, &answ, &vector, N, rank, size);
	double const T = 1e-4;
	int count = 0;
	do {
		approximate(vector, matr, answ, N, T, rank, size, count);
	} while (!is_answer_correct(vector, matr, answ, N, rank, size, count++));

	//gettimeofday(&tv2, NULL);
	//double dt_sec = (tv2.tv_sec - tv1.tv_sec);
	//double dt_usec = (tv2.tv_usec - tv1.tv_usec);
	//double dt = dt_sec + 1e-6 * dt_usec;



	FILE *f = NULL;
	char filename[50];
	snprintf(filename, sizeof(char) * 50, "out_file%i.txt", rank);
	f = fopen(filename, "w");

	//fprintf(f, "time diff %e \n", dt);
	for (int i = 0; i < N; i++) {
		fprintf(f, "%f\n", vector[i]);
	}
	fclose(f);
	free(matr);
	free(answ);
	free(vector);
	MPI_Finalize();
	return 0;
}
