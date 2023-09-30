#ifndef CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H
#define CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H
#define EPSILON 0.000001

#include <cmath>
#include <iostream>

class S21Matrix {
 private:
  int rows_, cols_;
  double** matrix_;

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();

  int getRows() const;
  int getCols() const;
  double** getMatrix() const;
  void setRows(int rows);
  void setCols(int cols);
  void setMatrix(double** matrix);
  void printMatrix() const;
  bool EqMatrix(const S21Matrix& other) const;
  static bool comparator(double a, double b);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  friend S21Matrix operator+(S21Matrix& matrix1, S21Matrix& matrix2);
  friend S21Matrix operator-(S21Matrix& matrix1, S21Matrix& matrix2);
  friend S21Matrix operator*(S21Matrix& matrix1, S21Matrix& matrix2);
  friend S21Matrix operator*(double number, S21Matrix& matrix);
  friend S21Matrix operator*(S21Matrix& matrix, double number);
  friend bool operator==(S21Matrix& matrix1, S21Matrix& matrix2);
  S21Matrix& operator=(const S21Matrix& matrix);
  friend S21Matrix operator+=(S21Matrix& matrix1, S21Matrix& matrix2);
  friend S21Matrix operator-=(S21Matrix& matrix1, S21Matrix& matrix2);
  friend S21Matrix operator*=(S21Matrix& matrix1, S21Matrix& matrix2);
  friend S21Matrix operator*=(S21Matrix& matrix, double number);
  friend S21Matrix operator*=(double number, S21Matrix& matrix);
  double& operator()(int row, int col);
};

#endif  // CPP1_S21_MATRIXPLUS_1_S21_MATRIX_OOP_H
