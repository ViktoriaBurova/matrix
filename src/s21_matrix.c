#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
    int res = 0;
    result->columns = columns;
    result->rows = rows;
    if (result->columns > 0 && result->rows > 0) {
        result->matrix = (double **)malloc(result->rows*sizeof(double *));
        for (int i = 0; i < result->rows; i++) {
            result->matrix[i] = (double *)malloc(result->columns*sizeof(double));
        }
    } else {
        res = 1;
    }
    return res;
}

void s21_remove_matrix(matrix_t *A) {
    for (int i = 0; i < A->rows; i++) {
        free(A->matrix[i]);
    }
    free(A->matrix);
    A->columns = 0;
    A->rows = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
    int res = FAILURE;
    int count = 0;
    if (A->columns == B->columns && A->rows == B->rows && A->columns > 0
    && B->columns > 0 && A->rows > 0 && B->rows > 0) {
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                if (fabs(A->matrix[i][j] - B->matrix[i][j]) <= 1e-6) {
                    count++;
                }
            }
        }
        if (count == (A->columns * A->rows)) {
            res = SUCCESS;
        }
    }
    return res;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int res = 0;
    if (A->matrix != NULL && B->matrix != NULL && A->columns > 0
    && B->columns > 0 && A->rows > 0 && B->rows > 0) {
        if (A->columns == B->columns && A->rows == B->rows) {
            s21_create_matrix(A->rows, A->columns, result);
            for (int i = 0; i < A->rows; i++) {
                for (int j = 0; j < A->columns; j++) {
                    result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
                }
            }
        } else {
            res = 2;
        }
    } else {
        res = 1;
    }
    return res;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int res = 0;
    if (A->matrix != NULL && B->matrix != NULL && A->columns > 0
    && B->columns > 0 && A->rows > 0 && B->rows > 0) {
        if (A->columns == B->columns && A->rows == B->rows) {
            s21_create_matrix(A->rows, A->columns, result);
            for (int i = 0; i < A->rows; i++) {
                for (int j = 0; j < A->columns; j++) {
                    result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
                }
            }
        } else {
            res = 2;
        }
    } else {
        res = 1;
    }
    return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
    int res = 0;
    if (A->matrix != NULL && A->columns > 0 && A->rows > 0) {
        s21_create_matrix(A->rows, A->columns, result);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                result->matrix[i][j] = number * A->matrix[i][j];
            }
        }
    } else {
        res = 1;
    }
    return res;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int res = 0;
    if (A->matrix != NULL && B->matrix != NULL && A->columns > 0
    && B->columns > 0 && A->rows > 0 && B->rows > 0) {
        if (A->columns == B->rows) {
            s21_create_matrix(A->rows, B->columns, result);
            for (int i = 0; i < A->rows; i++) {
                for (int j = 0; j < B->columns; j++) {
                    int sum = 0;
                    for (int m = 0; m < B->rows; m++) {
                        sum = sum + A->matrix[i][m] * B->matrix[m][j];
                    }
                    result->matrix[i][j] = sum;
                }
            }
        } else {
            res = 2;
        }
    } else {
        res = 1;
    }
    return res;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
    int res = 0;
    if (A->matrix != NULL && A->columns > 0 && A->rows > 0) {
        s21_create_matrix(A->columns, A->rows, result);
        for (int i = 0; i < result->rows; i++) {
            for (int j = 0; j < result->columns; j++) {
                result->matrix[i][j] = A->matrix[j][i];
            }
        }
    } else {
        res = 1;
    }
    return res;
}

int s21_determinant(matrix_t *A, double *result) {
    int res = 0;
    if (A->matrix != NULL && A->rows > 0 && A->columns > 0) {
        if (A->columns == A->rows) {
            if (A->columns == 2) {
                *result = A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
            }
            if (A->columns == 1) {
                *result = A->matrix[0][0];
            }
            if (A->columns > 2) {
                matrix_t B;
                s21_create_matrix(A->rows - 1, A->rows - 1, &B);
                double result1 = 0;
                *result = 0;
                for (int j = 0; j < A->rows; j++) {
                    for (int n = 0; n < A->rows - 1; n++) {
                        for (int m = 0; m < A->rows - 1; m++) {
                            if (m >= j) {
                                B.matrix[n][m] = A->matrix[n + 1][m + 1];
                            } else {
                                B.matrix[n][m] = A->matrix[n + 1][m];
                            }
                        }
                    }
                    s21_determinant(&B, &result1);
                    *result = *result + A->matrix[0][j] * result1 * pow(-1, j);
                }
                s21_remove_matrix(&B);
            }
        } else {
            res = 2;
        }
    } else {
        res = 1;
    }
    return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
    int res = 0;
    if (A->matrix != NULL && A->rows > 0 && A->columns > 0) {
        if (A->columns == A->rows) {
            s21_create_matrix(A->rows, A->columns, result);
            matrix_t B;
            s21_create_matrix(A->rows-1, A->rows-1, &B);
            double result1 = 0;
            for (int i = 0; i < A->rows; i++) {
                for (int j = 0; j < A->rows; j++) {
                    for (int n = 0; n < A->rows - 1; n++) {
                        for (int m = 0; m < A->rows - 1; m++) {
                            if (n < i && m < j) {
                                B.matrix[n][m] = A->matrix[n][m];
                            } else if (n < i && m >= j) {
                                B.matrix[n][m] = A->matrix[n][m + 1];
                            } else if (n >= i && m < j) {
                                B.matrix[n][m] = A->matrix[n + 1][m];
                            } else if (n >= i && m >= j) {
                                B.matrix[n][m] = A->matrix[n + 1][m + 1];
                            }
                        }
                    }
                    s21_determinant(&B, &result1);
                    result->matrix[i][j] = result1 * pow(-1, j + i);
                }
            }
            s21_remove_matrix(&B);
        } else {
            res = 2;
        }
    } else {
        res = 1;
    }
    return res;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
    int res = 0;
    if (A->matrix != NULL && A->rows > 0 && A->columns > 0) {
        if (A->columns == A->rows) {
            double result1 = 0;
            s21_determinant(A, &result1);
            if (fabs(result1) > 1e-7) {
                matrix_t calc_complements;
                matrix_t transpose;
                s21_calc_complements(A, &calc_complements);
                s21_transpose(&calc_complements, &transpose);
                s21_mult_number(&transpose, 1 / result1, result);
                s21_remove_matrix(&transpose);
                s21_remove_matrix(&calc_complements);
            } else {
                res = 2;
            }
        } else {
            res = 2;
        }
    } else {
        res = 1;
    }
    return res;
}
