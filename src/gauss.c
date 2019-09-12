#include <stdio.h>
#include <glib.h>
#include <gsl/gsl_matrix.h>

#define MATRIX_SIZE 3

static double d_matrix[MATRIX_SIZE][MATRIX_SIZE] = {
  {1., 4., 2.},
  {-1., -2., 1.},
  {3., 20., 19.},
};
static double d_vector[3] = {8, 3, 71};

int main(void) {
  printf("Gauss.\n\n");

  // Alloc matrix`s
  gsl_matrix *matrix = gsl_matrix_calloc(MATRIX_SIZE, MATRIX_SIZE);
  gsl_matrix *vector = gsl_matrix_calloc(1, MATRIX_SIZE);
  // gsl_matrix *result = gsl_matrix_calloc(1, MATRIX_SIZE);

  // Init matrix
  for (size_t i = 0; i < MATRIX_SIZE; i++) {
    for (size_t j = 0; j < MATRIX_SIZE; j++) {
      gsl_matrix_set(matrix, i, j, d_matrix[i][j]);
    }
  }

  // Init vector
  for (size_t j = 0; j < MATRIX_SIZE; j++) {
    gsl_matrix_set(vector, 0, j, d_vector[j]);
  }

  printf("First matrix row:\n");

  for (size_t j = 0; j < MATRIX_SIZE; j++) {
    double el = gsl_matrix_get(matrix, 0, j);
    printf("%zu - %f\n", j, el);
  }
  printf("\n");


  // Gauss

  // Forward

  // Steps by equations
  for (size_t step = 0; step < MATRIX_SIZE - 1; step++) {

    // Walk by matrix rows
    for (size_t row = step + 1; row < MATRIX_SIZE; row++) {
      double multiplier = d_matrix[row][step] / d_matrix[step][step];

      // Update vector value
      d_vector[row] = d_vector[row] - multiplier * d_vector[step];

      // Walk by matrix cells
      for (size_t col = step + 1, col0 = step; col < MATRIX_SIZE; col++, col0++) {
        double *cell = &(d_matrix[row][col]);

        *cell = *cell - multiplier * d_matrix[step][col];

        // Zero processed cells
        if (col0 < row) {
          d_matrix[row][col0] = 0;
        }
      }
    }
  }

  // /Forward

  // Check det

  // Back

  //for (size_t) {

  //}

  // /Back

  // /Gauss

/*
  printf("Result:\n");

  for (size_t j = 0; j < MATRIX_SIZE; j++) {
    double el = gsl_matrix_get(result, 0, j);
    printf("%zu - %f\n", j, el);
  }
*/
  return 0;
}
