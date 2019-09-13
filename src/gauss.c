#include <stdio.h>
#include <glib.h>
#include <gsl/gsl_matrix.h>

#define MATRIX_SIZE 3

static double d_matrix[MATRIX_SIZE][MATRIX_SIZE] = {
  {1., 1., 1.},
  {2., 3., 1.},
  {1., -1., -1.},
};
static double d_orig_matrix[MATRIX_SIZE][MATRIX_SIZE] = {
  {1., 1., 1.},
  {2., 3., 1.},
  {1., -1., -1.},
};
static double d_vector[] = {4., 9., -2.};
static double d_vector_orig[] = {4., 9., -2.};
static double r_vector[] = {0, 0, 0};

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
      d_matrix[row][step] = 0;

      // Update vector value
      d_vector[row] = d_vector[row] - multiplier * d_vector[step];

      // Walk by matrix cells
      for (size_t col = step + 1; col < MATRIX_SIZE; col++) {
        double *cell = &(d_matrix[row][col]);

        *cell = *cell - multiplier * d_matrix[step][col];
      }
    }
  }

  // /Forward

  // Check det

  // Back

  for (ssize_t eqIdx = MATRIX_SIZE - 1; eqIdx >= 0; eqIdx--) { //1
    double sum = 0;

    for (size_t col = eqIdx + 1; col <= MATRIX_SIZE - 1; col++) { // 2
      sum += d_matrix[eqIdx][col] * r_vector[col];
    }

    r_vector[eqIdx] = (d_vector[eqIdx] - sum) / d_matrix[eqIdx][eqIdx];
  }

  // /Back

  // /Gauss


  printf("Result:\n");

  for (size_t i = 0; i < MATRIX_SIZE; i++) {
    // double el = gsl_matrix_get(r_vector, 0, j);
    printf("%zu - %f\n", i, r_vector[i]);
  }

  printf("\nCheck:\n");

  for (size_t row = 0; row < MATRIX_SIZE; row++) {
    double sum = 0;

    for (size_t col = 0; col < MATRIX_SIZE; col++) {
      sum += d_orig_matrix[row][col] * r_vector[col];
    }

    // double el = gsl_matrix_get(r_vector, 0, j);
    printf("%f = %f\n", sum, d_vector_orig[row]);
  }

  return 0;
}
