#include <stdio.h>
#include <glib.h>
#include <gsl/gsl_matrix.h>

#define MATRIX_SIZE 3

static double d_matrix[MATRIX_SIZE][MATRIX_SIZE] = {
  {1., 1., 1.},
  {2., 3., 1.},
  {1., -1., -1.},
};
static double d_vector[] = {4., 9., -2.};

int main(void) {
  printf("Gauss.\n\n");

  // Alloc matrix`s
  gsl_matrix *matrix = gsl_matrix_calloc(MATRIX_SIZE, MATRIX_SIZE);
  gsl_matrix *orig_matrix = gsl_matrix_calloc(MATRIX_SIZE, MATRIX_SIZE);

  gsl_vector *vector = gsl_vector_calloc(MATRIX_SIZE);
  gsl_vector *orig_vector = gsl_vector_calloc(MATRIX_SIZE);

  gsl_vector *result = gsl_vector_calloc(MATRIX_SIZE);

  // Init matrix
  for (size_t i = 0; i < MATRIX_SIZE; i++) {
    for (size_t j = 0; j < MATRIX_SIZE; j++) {
      gsl_matrix_set(matrix, i, j, d_matrix[i][j]);
      gsl_matrix_set(orig_matrix, i, j, d_matrix[i][j]);
    }
  }

  // Init vector
  for (size_t j = 0; j < MATRIX_SIZE; j++) {
    gsl_vector_set(vector, j, d_vector[j]);
    gsl_vector_set(orig_vector, j, d_vector[j]);
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
    for (size_t eq_idx = step + 1; eq_idx < MATRIX_SIZE; eq_idx++) {
      double multiplier = gsl_matrix_get(matrix, eq_idx, step) / gsl_matrix_get(matrix, step, step);

      gsl_matrix_set(matrix, eq_idx, step, 0);

      // Update vector value
      double vector_val = gsl_vector_get(vector, eq_idx) - multiplier * gsl_vector_get(vector, step);
      gsl_vector_set(vector, eq_idx, vector_val);

      // Walk by eq cells
      for (size_t col = step + 1; col < MATRIX_SIZE; col++) {
        double cell_val = gsl_matrix_get(matrix, eq_idx, col) - multiplier * gsl_matrix_get(matrix, step, col);

        gsl_matrix_set(matrix, eq_idx, col, cell_val);
      }
    }
  }

  // /Forward

  // Check det

  // Back

  for (ssize_t eqIdx = MATRIX_SIZE - 1; eqIdx >= 0; eqIdx--) { //1
    double sum = 0;

    for (size_t col = eqIdx + 1; col <= MATRIX_SIZE - 1; col++) { // 2
      sum += gsl_matrix_get(matrix, eqIdx, col) * gsl_vector_get(vector, col);
    }

    gsl_vector_set(result, eqIdx, (gsl_vector_get(vector, eqIdx) - sum) / gsl_matrix_get(matrix, eqIdx, eqIdx));
  }

  // /Back

  // /Gauss


  printf("Result:\n");

  for (size_t i = 0; i < MATRIX_SIZE; i++) {
    // double el = gsl_matrix_get(r_vector, 0, j);
    printf("%zu - %f\n", i, gsl_vector_get(result, i));
  }

  printf("\nCheck:\n");

  for (size_t row = 0; row < MATRIX_SIZE; row++) {
    double sum = 0;

    for (size_t col = 0; col < MATRIX_SIZE; col++) {
      sum += gsl_matrix_get(orig_matrix, row, col) * gsl_vector_get(result, col);
    }

    // double el = gsl_matrix_get(r_vector, 0, j);
    printf("%f = %f\n", sum, gsl_vector_get(orig_vector, row));
  }

  return 0;
}
