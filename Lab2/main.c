#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include <sys/time.h>

#include "operations.h"

int main(int argc, char *argv[]) {
    int const N = 100000;
    struct timeval tv1, tv2;
    double *matr = NULL, *answ = NULL, *vector = NULL;
    preset_solution(&matr, &answ, &vector, N);
    double const T = 1e-5;
    int count = 0;
    double sq_norm;
    gettimeofday(&tv1, NULL);
    do {
        sq_norm = approximate(vector, matr, answ, N, T, count);
        print_vec("normParallel", &sq_norm, 1, count, "\n\r", 1, 1 );
    } while (!is_answer_correct(sq_norm, answ, N, count++));
    gettimeofday(&tv2, NULL);
    double dt_sec = (tv2.tv_sec - tv1.tv_sec);
    double dt_usec = (tv2.tv_usec - tv1.tv_usec);
    double dt = dt_sec + 1e-6 * dt_usec;

    print_vec("answerPragmaMod", vector, 1, 1, "Answer\n\r", N, 1);
    print_vec("TimeWithModPragma", &dt, 1, 1, "\n\r", 1, 1);
//    printf("%f", dt);
    free(matr);
    free(answ);
    free(vector);
    return 0;
}
