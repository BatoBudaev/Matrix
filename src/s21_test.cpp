#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

TEST(S21MatrixTest, DefaultConstructor) {
  S21Matrix matrix;
  EXPECT_EQ(matrix.getRows(), 0);
  EXPECT_EQ(matrix.getCols(), 0);
  EXPECT_EQ(matrix.getMatrix(), nullptr);
}

TEST(S21MatrixTest, ConstructorWithSize) {
  S21Matrix matrix(3, 4);
  EXPECT_EQ(matrix.getRows(), 3);
  EXPECT_EQ(matrix.getCols(), 4);
  double** matrixData = matrix.getMatrix();
  EXPECT_NE(matrixData, nullptr);

  for (int i = 0; i < matrix.getRows(); i++) {
    for (int j = 0; j < matrix.getCols(); j++) {
      EXPECT_EQ(matrixData[i][j], 0.0);
    }
  }
}

TEST(S21MatrixTest, DecreaseRowsAndCols) {
  S21Matrix matrix1(3, 3);
  S21Matrix expected(2, 2);
  double value1 = 1.0;

  for (int i = 0; i < matrix1.getRows(); ++i) {
    for (int j = 0; j < matrix1.getCols(); ++j) {
      matrix1.getMatrix()[i][j] = value1;
      value1 += 1.0;
    }
  }

  expected(0, 0) = 1;
  expected(0, 1) = 2;
  expected(1, 0) = 4;
  expected(1, 1) = 5;

  matrix1.setRows(2);
  matrix1.setCols(2);

  EXPECT_TRUE(matrix1 == expected);
}

TEST(S21MatrixTest, IncreaseRowsAndCols) {
  S21Matrix matrix1(3, 3);
  S21Matrix expected(4, 4);
  double value1 = 1.0;

  for (int i = 0; i < matrix1.getRows(); ++i) {
    for (int j = 0; j < matrix1.getCols(); ++j) {
      matrix1.getMatrix()[i][j] = value1;
      value1 += 1.0;
    }
  }

  value1 = 1.0;

  for (int i = 0; i < matrix1.getRows(); ++i) {
    for (int j = 0; j < matrix1.getCols(); ++j) {
      expected.getMatrix()[i][j] = value1;
      value1 += 1.0;
    }
  }

  matrix1.setRows(4);
  matrix1.setCols(4);

  EXPECT_TRUE(matrix1 == expected);
}

TEST(S21MatrixTest, MatrixEquality) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(3, 3);

  matrix1(0, 0) = 0.123456;
  matrix2(0, 0) = 0.123456;

  EXPECT_TRUE(matrix1 == matrix2);
}

TEST(S21MatrixTest, MatrixInequality) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(3, 3);

  matrix1(0, 0) = 1;
  matrix2(0, 0) = 2;

  EXPECT_FALSE(matrix1 == matrix2);
}

TEST(S21MatrixTest, MatrixSum) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(3, 3);
  S21Matrix expected(3, 3);

  matrix1(0, 0) = 1;
  matrix1(0, 1) = 2;
  matrix1(0, 2) = 3;
  matrix1(1, 0) = 0;
  matrix1(1, 1) = 4;
  matrix1(1, 2) = 5;
  matrix1(2, 0) = 0;
  matrix1(2, 1) = 0;
  matrix1(2, 2) = 6;

  matrix2(0, 0) = 1;
  matrix2(0, 1) = 0;
  matrix2(0, 2) = 0;
  matrix2(1, 0) = 2;
  matrix2(1, 1) = 0;
  matrix2(1, 2) = 0;
  matrix2(2, 0) = 3;
  matrix2(2, 1) = 4;
  matrix2(2, 2) = 1;

  expected(0, 0) = 2;
  expected(0, 1) = 2;
  expected(0, 2) = 3;
  expected(1, 0) = 2;
  expected(1, 1) = 4;
  expected(1, 2) = 5;
  expected(2, 0) = 3;
  expected(2, 1) = 4;
  expected(2, 2) = 7;

  S21Matrix result(3, 3);

  result = matrix1 + matrix2;

  EXPECT_TRUE(result == expected);
}

TEST(S21MatrixTest, MatrixSubtract) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(3, 3);
  S21Matrix expected(3, 3);

  matrix1(0, 0) = 1;
  matrix1(0, 1) = 2;
  matrix1(0, 2) = 3;
  matrix1(1, 0) = 0;
  matrix1(1, 1) = 4;
  matrix1(1, 2) = 5;
  matrix1(2, 0) = 0;
  matrix1(2, 1) = 0;
  matrix1(2, 2) = 6;

  matrix2(0, 0) = 1;
  matrix2(0, 1) = 0;
  matrix2(0, 2) = 0;
  matrix2(1, 0) = 2;
  matrix2(1, 1) = 0;
  matrix2(1, 2) = 0;
  matrix2(2, 0) = 3;
  matrix2(2, 1) = 4;
  matrix2(2, 2) = 1;

  expected(0, 0) = 0;
  expected(0, 1) = 2;
  expected(0, 2) = 3;
  expected(1, 0) = -2;
  expected(1, 1) = 4;
  expected(1, 2) = 5;
  expected(2, 0) = -3;
  expected(2, 1) = -4;
  expected(2, 2) = 5;

  S21Matrix result(3, 3);

  result = matrix1 - matrix2;

  EXPECT_TRUE(result == expected);
}

TEST(S21MatrixTest, MatrixMultiplyNumber) {
  S21Matrix matrix1(3, 3);
  S21Matrix expected(3, 3);
  double number = 2;

  matrix1(0, 0) = 1;
  matrix1(0, 1) = 2;
  matrix1(0, 2) = 3;
  matrix1(1, 0) = 0;
  matrix1(1, 1) = 4;
  matrix1(1, 2) = 2;
  matrix1(2, 0) = 2;
  matrix1(2, 1) = 3;
  matrix1(2, 2) = 4;

  expected(0, 0) = 2;
  expected(0, 1) = 4;
  expected(0, 2) = 6;
  expected(1, 0) = 0;
  expected(1, 1) = 8;
  expected(1, 2) = 4;
  expected(2, 0) = 4;
  expected(2, 1) = 6;
  expected(2, 2) = 8;

  S21Matrix result(3, 3);

  result = matrix1 * number;

  EXPECT_TRUE(result == expected);
}

TEST(S21MatrixTest, MatrixMultiply) {
  S21Matrix matrix1(3, 2);
  S21Matrix matrix2(2, 3);
  S21Matrix expected(3, 3);

  matrix1(0, 0) = 1;
  matrix1(0, 1) = 4;
  matrix1(1, 0) = 2;
  matrix1(1, 1) = 5;
  matrix1(2, 0) = 3;
  matrix1(2, 1) = 6;

  matrix2(0, 0) = 1;
  matrix2(0, 1) = -1;
  matrix2(0, 2) = 1;
  matrix2(1, 0) = 2;
  matrix2(1, 1) = 3;
  matrix2(1, 2) = 4;

  expected(0, 0) = 9;
  expected(0, 1) = 11;
  expected(0, 2) = 17;
  expected(1, 0) = 12;
  expected(1, 1) = 13;
  expected(1, 2) = 22;
  expected(2, 0) = 15;
  expected(2, 1) = 15;
  expected(2, 2) = 27;

  S21Matrix result(3, 3);

  result = matrix1 * matrix2;

  EXPECT_TRUE(result == expected);
}

TEST(S21MatrixTest, MatrixTranspose) {
  S21Matrix matrix1(3, 2);
  S21Matrix expected(2, 3);

  matrix1(0, 0) = 1;
  matrix1(0, 1) = 4;
  matrix1(1, 0) = 2;
  matrix1(1, 1) = 5;
  matrix1(2, 0) = 3;
  matrix1(2, 1) = 6;

  expected(0, 0) = 1;
  expected(0, 1) = 2;
  expected(0, 2) = 3;
  expected(1, 0) = 4;
  expected(1, 1) = 5;
  expected(1, 2) = 6;

  S21Matrix result(2, 3);

  result = matrix1.Transpose();

  EXPECT_TRUE(result == expected);
}

TEST(S21MatrixTest, MatrixComplements) {
  S21Matrix matrix1(3, 3);
  S21Matrix expected(3, 3);

  matrix1(0, 0) = 1;
  matrix1(0, 1) = 2;
  matrix1(0, 2) = 3;
  matrix1(1, 0) = 0;
  matrix1(1, 1) = 4;
  matrix1(1, 2) = 2;
  matrix1(2, 0) = 5;
  matrix1(2, 1) = 2;
  matrix1(2, 2) = 1;

  expected(0, 0) = 0;
  expected(0, 1) = 10;
  expected(0, 2) = -20;
  expected(1, 0) = 4;
  expected(1, 1) = -14;
  expected(1, 2) = 8;
  expected(2, 0) = -8;
  expected(2, 1) = -2;
  expected(2, 2) = 4;

  S21Matrix result(3, 3);

  result = matrix1.CalcComplements();

  EXPECT_TRUE(result == expected);
}

TEST(S21MatrixTest, MatrixDeterminant) {
  S21Matrix matrix1(3, 3);
  double expected = -230;

  matrix1(0, 0) = -1;
  matrix1(0, 1) = 2;
  matrix1(0, 2) = 5;
  matrix1(1, 0) = 7;
  matrix1(1, 1) = -4;
  matrix1(1, 2) = 3;
  matrix1(2, 0) = -5;
  matrix1(2, 1) = 0;
  matrix1(2, 2) = 10;

  double result = matrix1.Determinant();

  EXPECT_DOUBLE_EQ(result, expected);
}

TEST(S21MatrixTest, MatrixInverse) {
  S21Matrix matrix1(3, 3);
  S21Matrix expected(3, 3);

  matrix1(0, 0) = 2;
  matrix1(0, 1) = 5;
  matrix1(0, 2) = 7;
  matrix1(1, 0) = 6;
  matrix1(1, 1) = 3;
  matrix1(1, 2) = 4;
  matrix1(2, 0) = 5;
  matrix1(2, 1) = -2;
  matrix1(2, 2) = -3;

  expected(0, 0) = 1;
  expected(0, 1) = -1;
  expected(0, 2) = 1;
  expected(1, 0) = -38;
  expected(1, 1) = 41;
  expected(1, 2) = -34;
  expected(2, 0) = 27;
  expected(2, 1) = -29;
  expected(2, 2) = 24;

  S21Matrix result(3, 3);

  result = matrix1.InverseMatrix();

  EXPECT_TRUE(result == expected);
}

TEST(S21MatrixTest, MatrixSumWithEmptyMatrix) {
  S21Matrix matrix1(3, 3);
  S21Matrix emptyMatrix;
  S21Matrix result(3, 3);

  matrix1(0, 0) = 2;
  matrix1(0, 1) = 5;
  matrix1(0, 2) = 7;
  matrix1(1, 0) = 6;
  matrix1(1, 1) = 3;
  matrix1(1, 2) = 4;
  matrix1(2, 0) = 5;
  matrix1(2, 1) = -2;
  matrix1(2, 2) = -3;

  EXPECT_THROW(result = matrix1 + emptyMatrix, std::invalid_argument);
}

TEST(S21MatrixTest, MatrixOperatorPlusEqual) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(3, 3);

  matrix1(0, 0) = 1;
  matrix1(0, 1) = 2;
  matrix1(0, 2) = 7;
  matrix1(1, 0) = 0;
  matrix1(1, 1) = 4;
  matrix1(1, 2) = 2;
  matrix1(2, 0) = 5;
  matrix1(2, 1) = 2;
  matrix1(2, 2) = 1;

  double value1 = 1.0;

  for (int i = 0; i < matrix2.getRows(); ++i) {
    for (int j = 0; j < matrix2.getCols(); ++j) {
      matrix2.getMatrix()[i][j] = value1;
      value1 += 1.0;
    }
  }

  matrix1 += matrix2;

  EXPECT_EQ(matrix1(0, 0), 2.0);
  EXPECT_EQ(matrix1(0, 2), 10.0);
}

TEST(S21MatrixTest, MatrixOperatorMinusEqual) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(3, 3);

  matrix1(0, 0) = 1;
  matrix1(0, 1) = 2;
  matrix1(0, 2) = 7;
  matrix1(1, 0) = 0;
  matrix1(1, 1) = 4;
  matrix1(1, 2) = 2;
  matrix1(2, 0) = 5;
  matrix1(2, 1) = 2;
  matrix1(2, 2) = 1;

  matrix2(0, 0) = 1;
  matrix2(0, 1) = 2;
  matrix2(0, 2) = 1;
  matrix2(1, 0) = 0;
  matrix2(1, 1) = 1;
  matrix2(1, 2) = 2;
  matrix2(2, 0) = 2;
  matrix2(2, 1) = 2;
  matrix2(2, 2) = 2;

  matrix1 -= matrix2;

  EXPECT_EQ(matrix1(0, 0), 0.0);
  EXPECT_EQ(matrix1(0, 1), 0.0);
  EXPECT_EQ(matrix1(0, 2), 6.0);
}

TEST(S21MatrixTest, MatrixOperatorMultiplyEqual) {
  S21Matrix matrix1(3, 3);
  S21Matrix matrix2(3, 3);
  S21Matrix expected(3, 3);

  matrix1(0, 0) = 1;
  matrix1(0, 1) = 4;
  matrix1(1, 0) = 2;
  matrix1(1, 1) = 5;
  matrix1(2, 0) = 3;
  matrix1(2, 1) = 6;

  matrix2(0, 0) = 1;
  matrix2(0, 1) = -1;
  matrix2(0, 2) = 1;
  matrix2(1, 0) = 2;
  matrix2(1, 1) = 3;
  matrix2(1, 2) = 4;

  expected(0, 0) = 9;
  expected(0, 1) = 11;
  expected(0, 2) = 17;
  expected(1, 0) = 12;
  expected(1, 1) = 13;
  expected(1, 2) = 22;
  expected(2, 0) = 15;
  expected(2, 1) = 15;
  expected(2, 2) = 27;

  S21Matrix result(3, 3);

  matrix1 *= matrix2;

  EXPECT_TRUE(matrix1 == expected);
}

TEST(S21MatrixTest, MatrixOperatorMultiplyEqualNumber) {
  S21Matrix matrix(3, 3);
  double number = 2.0;

  matrix(0, 0) = 1;
  matrix(0, 1) = 2;
  matrix(0, 2) = 3;
  matrix(1, 0) = 4;
  matrix(1, 1) = 5;
  matrix(1, 2) = 6;
  matrix(2, 0) = 7;
  matrix(2, 1) = 8;
  matrix(2, 2) = 9;

  matrix *= number;

  EXPECT_EQ(matrix(0, 0), 2.0);
  EXPECT_EQ(matrix(0, 1), 4.0);
  EXPECT_EQ(matrix(0, 2), 6.0);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
