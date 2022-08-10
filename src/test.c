#include <check.h>
#include "s21_matrix.h"

START_TEST(matrix_test1) {
  matrix_t matrix1, matrix2;
  int ret = s21_create_matrix(0, 1, &matrix1);
  ck_assert_int_eq(ret, 1);
  ret = s21_create_matrix(3, 3, &matrix1);
  ck_assert_int_eq(ret, 0);
  s21_remove_matrix(&matrix1);
  s21_create_matrix(3, 3, &matrix1);
  s21_create_matrix(3, 3, &matrix2);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      matrix2.matrix[i][j] = k++;
    }
  }
  ck_assert_int_eq(s21_eq_matrix(&matrix1, &matrix2), 0);
  k = 0;
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      matrix2.matrix[i][j] = k++;
    }
  }
  ck_assert_int_eq(s21_eq_matrix(&matrix1, &matrix2), 1);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
}
END_TEST

START_TEST(matrix_add_sub_test) {
  matrix_t matrix1, matrix2, result;
  s21_create_matrix(3, 3, &matrix1);
  s21_create_matrix(3, 3, &matrix2);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  k = 0;
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      matrix2.matrix[i][j] = k++;
    }
  }
  s21_sum_matrix(&matrix2, &matrix1, &result);
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      ck_assert_int_eq(result.matrix[i][j],
                       matrix2.matrix[i][j] + matrix1.matrix[i][j]);
    }
  }
  s21_remove_matrix(&result);
  s21_sub_matrix(&matrix2, &matrix1, &result);
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j],
                              matrix2.matrix[i][j] - matrix1.matrix[i][j],
                              1e-6);
    }
  }
  s21_remove_matrix(&result);

  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      srand(i + j);
      matrix1.matrix[i][j] = (double)rand() / (double)rand();
      srand((j + i) * 2);
      matrix2.matrix[i][j] = (double)rand() / (double)rand();
    }
  }

  s21_sum_matrix(&matrix2, &matrix1, &result);
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j],
                              matrix2.matrix[i][j] + matrix1.matrix[i][j],
                              1e-6);
    }
  }
  s21_remove_matrix(&result);

  s21_sub_matrix(&matrix2, &matrix1, &result);
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j],
                              matrix2.matrix[i][j] - matrix1.matrix[i][j],
                              1e-6);
    }
  }

  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
  s21_remove_matrix(&result);
  s21_create_matrix(2, 3, &matrix1);
  s21_create_matrix(3, 3, &matrix2);
  s21_create_matrix(3, 3, &result);
  int ret = s21_sum_matrix(&matrix1, &matrix2, &result);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
  s21_remove_matrix(&result);
  s21_create_matrix(2, 3, &matrix1);
  s21_create_matrix(3, 3, &matrix2);
  s21_create_matrix(3, 3, &result);
  ret = s21_sub_matrix(&matrix1, &matrix2, &result);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
  s21_remove_matrix(&result);

  matrix1.columns = 0;
  matrix1.rows = 3;
  matrix2.columns = 3;
  matrix2.rows = 3;

  ret = s21_sub_matrix(&matrix1, &matrix2, &result);
  ck_assert_int_eq(ret, 1);

  ret = s21_sum_matrix(&matrix1, &matrix2, &result);
  ck_assert_int_eq(ret, 1);
}

START_TEST(matrix_mul_test) {
  matrix_t matrix1, matrix2, result;
  s21_create_matrix(3, 3, &matrix1);
  s21_create_matrix(3, 3, &matrix2);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  k = 0;
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      matrix2.matrix[i][j] = k++;
    }
  }
  s21_mult_number(&matrix2, 2, &result);
  for (int i = 0; i < matrix2.rows; i++) {
    for (int j = 0; j < matrix2.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j], matrix2.matrix[i][j] * 2,
                              1e-6);
    }
  }
  s21_remove_matrix(&result);
  int ret = s21_mult_matrix(&matrix1, &matrix2, &result);
  ck_assert_int_eq(ret, 0);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
  s21_remove_matrix(&result);

  s21_create_matrix(2, 3, &matrix1);
  s21_create_matrix(2, 3, &matrix2);
  s21_create_matrix(3, 3, &result);
  ret = s21_mult_matrix(&matrix1, &matrix2, &result);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
  s21_remove_matrix(&result);

  matrix1.columns = 0;
  matrix1.rows = 3;

  ret = s21_mult_number(&matrix1, 2, &result);
  ck_assert_int_eq(ret, 1);

  matrix_t matrix;
  matrix.matrix = NULL;
  ret = s21_mult_matrix(&matrix1, &matrix, &result);
  ck_assert_int_eq(ret, 1);
}
END_TEST

START_TEST(crazy_s21_mul_test) {
  matrix_t matrix1;
  matrix_t matrix2;

  s21_create_matrix(3, 3, &matrix1);
  s21_create_matrix(3, 3, &matrix2);

  matrix1.matrix[0][0] = 5.474748;
  matrix1.matrix[0][1] = -23.365346;
  matrix1.matrix[0][2] = 7.464645;
  matrix1.matrix[1][0] = 13.235363;
  matrix1.matrix[1][1] = -17.326774;
  matrix1.matrix[1][2] = -0.000034;
  matrix1.matrix[2][0] = -12.235567;
  matrix1.matrix[2][1] = 11.124526;
  matrix1.matrix[2][2] = 5.325634;

  matrix2.matrix[0][0] = 15.352534;
  matrix2.matrix[0][1] = 143.532636;
  matrix2.matrix[0][2] = 345.35356;
  matrix2.matrix[1][0] = 124.523552;
  matrix2.matrix[1][1] = -654.234562;
  matrix2.matrix[1][2] = 123.353578;
  matrix2.matrix[2][0] = -245.636465;
  matrix2.matrix[2][1] = 6324.235523;
  matrix2.matrix[2][2] = -2353.632542;
  matrix_t result;
  int ret = s21_mult_matrix(&matrix1, &matrix2, &result);
  ck_assert_int_eq(ret, 0);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(transpose_test) {
  matrix_t matrix1, matrix2, result;
  s21_create_matrix(3, 3, &matrix1);
  s21_create_matrix(3, 3, &matrix2);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  matrix2.matrix[0][0] = 0;
  matrix2.matrix[0][1] = 3;
  matrix2.matrix[0][2] = 6;
  matrix2.matrix[1][0] = 1;
  matrix2.matrix[1][1] = 4;
  matrix2.matrix[1][2] = 7;
  matrix2.matrix[2][0] = 2;
  matrix2.matrix[2][1] = 5;
  matrix2.matrix[2][2] = 8;
  int ret = s21_transpose(&matrix1, &result);
  for (int i = 0; i < result.rows; i++) {
    for (int j = 0; j < result.columns; j++) {
      ck_assert_double_eq_tol(result.matrix[i][j], matrix2.matrix[i][j], 1e-6);
    }
  }
  ck_assert_int_eq(ret, 0);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&matrix2);
  s21_remove_matrix(&result);

  matrix_t matrix;
  matrix.matrix = NULL;
  ret = s21_transpose(&matrix, &result);
  ck_assert_int_eq(ret, 1);
}
END_TEST

START_TEST(calc_complements_test) {
  matrix_t matrix1, result;
  s21_create_matrix(3, 2, &matrix1);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  int ret = s21_calc_complements(&matrix1, &result);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);

  s21_create_matrix(3, 3, &matrix1);
  k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  ret = s21_calc_complements(&matrix1, &result);
  ck_assert_int_eq(ret, 0);
  s21_remove_matrix(&matrix1);

  s21_remove_matrix(&result);
  s21_create_matrix(3, 3, &matrix1);
  matrix1.matrix[0][0] = 1;
  matrix1.matrix[0][1] = 2;
  matrix1.matrix[0][2] = 3;

  matrix1.matrix[1][0] = 0;
  matrix1.matrix[1][1] = 4;
  matrix1.matrix[1][2] = 2;

  matrix1.matrix[2][0] = 5;
  matrix1.matrix[2][1] = 2;
  matrix1.matrix[2][2] = 1;
  ret = s21_calc_complements(&matrix1, &result);
  matrix_t expected;
  s21_create_matrix(3, 3, &expected);
  expected.matrix[0][0] = 0;
  expected.matrix[0][1] = 10;
  expected.matrix[0][2] = -20;

  expected.matrix[1][0] = 4;
  expected.matrix[1][1] = -14;
  expected.matrix[1][2] = 8;

  expected.matrix[2][0] = -8;
  expected.matrix[2][1] = -2;
  expected.matrix[2][2] = 4;
  for (int i = 0; i < result.rows; i++) {
    for (int j = 0; j < result.columns; j++) {
      ck_assert_double_eq(result.matrix[i][j], expected.matrix[i][j]);
    }
  }
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&result);
  s21_remove_matrix(&expected);

  matrix_t matrix;
  matrix.matrix = NULL;
  ret = s21_calc_complements(&matrix, &result);
  ck_assert_int_eq(ret, 1);
}
END_TEST

START_TEST(determinant_test) {
  matrix_t matrix1;
  s21_create_matrix(3, 2, &matrix1);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  double res;
  int ret = s21_determinant(&matrix1, &res);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);

  s21_create_matrix(3, 3, &matrix1);
  k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  ret = s21_determinant(&matrix1, &res);
  ck_assert_int_eq(ret, 0);
  ck_assert_double_eq(res, 0);
  s21_remove_matrix(&matrix1);

  s21_create_matrix(3, 3, &matrix1);
  matrix1.matrix[0][0] = 0;
  matrix1.matrix[0][1] = 9;
  matrix1.matrix[0][2] = 5;

  matrix1.matrix[1][0] = 4;
  matrix1.matrix[1][1] = 3;
  matrix1.matrix[1][2] = -5;

  matrix1.matrix[2][0] = -1;
  matrix1.matrix[2][1] = 6;
  matrix1.matrix[2][2] = -4;
  ret = s21_determinant(&matrix1, &res);
  ck_assert_int_eq(ret, 0);
  ck_assert_double_eq(res, 324);
  s21_remove_matrix(&matrix1);

  matrix_t matrix;
  matrix.matrix = NULL;
  ret = s21_determinant(&matrix, &res);
  ck_assert_int_eq(ret, 1);

  s21_create_matrix(1, 1, &matrix1);
  matrix1.matrix[0][0] = 12;
  ret = s21_determinant(&matrix1, &res);
  ck_assert_int_eq(ret, 0);
  ck_assert_double_eq(res, 12);
  s21_remove_matrix(&matrix1);
}
END_TEST

START_TEST(inverse_test) {
  matrix_t A, result, matrix1;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 2;
  A.matrix[0][1] = 5;
  A.matrix[0][2] = 7;

  A.matrix[1][0] = 6;
  A.matrix[1][1] = 3;
  A.matrix[1][2] = 4;

  A.matrix[2][0] = 5;
  A.matrix[2][1] = -2;
  A.matrix[2][2] = -3;

  matrix_t expected_result;
  s21_create_matrix(3, 3, &expected_result);

  expected_result.matrix[0][0] = 1;
  expected_result.matrix[0][1] = -1;
  expected_result.matrix[0][2] = 1;

  expected_result.matrix[1][0] = -38;
  expected_result.matrix[1][1] = 41;
  expected_result.matrix[1][2] = -34;

  expected_result.matrix[2][0] = 27;
  expected_result.matrix[2][1] = -29;
  expected_result.matrix[2][2] = 24;

  matrix_t actual_result;
  s21_inverse_matrix(&A, &actual_result);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ck_assert_double_eq(expected_result.matrix[i][j],
                          actual_result.matrix[i][j]);
    }
  }
  s21_remove_matrix(&expected_result);
  s21_remove_matrix(&actual_result);
  s21_remove_matrix(&A);

  s21_create_matrix(3, 3, &matrix1);
  s21_create_matrix(3, 3, &result);
  int k = 0;
  for (int i = 0; i < matrix1.rows; i++) {
    for (int j = 0; j < matrix1.columns; j++) {
      matrix1.matrix[i][j] = k++;
    }
  }
  int ret = s21_inverse_matrix(&matrix1, &result);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);
  s21_remove_matrix(&result);

  s21_create_matrix(2, 3, &matrix1);
  ret = s21_inverse_matrix(&matrix1, &result);
  ck_assert_int_eq(ret, 2);
  s21_remove_matrix(&matrix1);

  matrix_t matrix;
  matrix.matrix = NULL;
  ret = s21_inverse_matrix(&matrix, &result);
  ck_assert_int_eq(ret, 1);
}
END_TEST

int main() {
  Suite *s1 = suite_create("Core");
  TCase *s21_test = tcase_create("Test");
  SRunner *sr = srunner_create(s1);
  suite_add_tcase(s1, s21_test);

  tcase_add_test(s21_test, matrix_test1);
  tcase_add_test(s21_test, matrix_add_sub_test);
  tcase_add_test(s21_test, matrix_mul_test);
  tcase_add_test(s21_test, transpose_test);
  tcase_add_test(s21_test, calc_complements_test);
  tcase_add_test(s21_test, determinant_test);
  tcase_add_test(s21_test, inverse_test);
  tcase_add_test(s21_test, crazy_s21_mul_test);

  srunner_run_all(sr, CK_VERBOSE);
  int errors = srunner_ntests_failed(sr);
  srunner_free(sr);

  return errors == 0 ? 0 : 1;
}
