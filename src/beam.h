#ifndef BEAM_H_
#define BEAM_H_

#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_cdf.h>

/*!
 * Расчет крена балки по заданным параметрам
 * Input:
 * deviation - стандартное отклонение случайной величины жесткости основания;
 * springStiffnessMean - среднее значение случайной величины жесткости основания;
 * lBeam - длина балки;
 * corRadius- радиус корреляции случайной величины жесткости основания.
 * Output:
 * tilt - вычесленное значение крена.
 */
void tilt_solver(const double deviation, const double springStiffnessMean, const double lBeam, const double corRadius, double* result, uint8_t tests);



#endif /* BEAM_H_ */