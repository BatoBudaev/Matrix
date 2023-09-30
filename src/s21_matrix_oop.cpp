#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr){};

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < 1 || cols < 1) {
    throw std::invalid_argument{"Invalid matrix"};
  }

  matrix_ = new double *[rows];

  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];

    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = 0.0;
    }
  }
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];

    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }

  delete[] matrix_;
}

void S21Matrix::printMatrix() const {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      std::cout << matrix_[i][j] << " ";
    }

    std::cout << std::endl;
  }
}

int S21Matrix::getRows() const { return rows_; }

void S21Matrix::setRows(int rows) {
  if (rows == rows_) {
    return;
  }

  double **newMatrix = new double *[rows];

  for (int i = 0; i < rows; i++) {
    newMatrix[i] = new double[cols_];

    for (int j = 0; j < cols_; j++) {
      if (i < rows_) {
        newMatrix[i][j] = matrix_[i][j];
      } else {
        newMatrix[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }

  delete[] matrix_;

  matrix_ = newMatrix;
  rows_ = rows;
}

int S21Matrix::getCols() const { return cols_; }

void S21Matrix::setCols(int cols) {
  if (cols == cols_) {
    return;
  }

  double **newMatrix = new double *[rows_];

  for (int i = 0; i < rows_; i++) {
    newMatrix[i] = new double[cols];

    for (int j = 0; j < cols; j++) {
      if (j < cols_) {
        newMatrix[i][j] = matrix_[i][j];
      } else {
        newMatrix[i][j] = 0.0;
      }
    }
  }

  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;

  matrix_ = newMatrix;
  cols_ = cols;
}

double **S21Matrix::getMatrix() const { return matrix_; }

void S21Matrix::setMatrix(double **matrix) { matrix_ = matrix; }

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (!comparator(matrix_[i][j], other.matrix_[i][j])) {
        return false;
      }
    }
  }

  return true;
}

bool S21Matrix::comparator(double a, double b) { return fabs(a - b) < EPSILON; }

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument{"Matrix sizes are not equal"};
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument{"Matrix sizes are not equal"};
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_) {
    throw std::invalid_argument(
        "The number of columns of the first matrix is not equal to the "
        "number of rows of the second matrix");
  }

  S21Matrix result(rows_, other.cols_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int k = 0; k < cols_; k++) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      matrix_[i][j] = result.matrix_[i][j];
    }
  }

  cols_ = other.cols_;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_, rows_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      result.matrix_[j][i] = matrix_[i][j];
    }
  }

  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::invalid_argument("The matrix must be square");
  }

  S21Matrix result(rows_, cols_);
  S21Matrix minor(rows_ - 1, cols_ - 1);

  for (int k = 0; k < rows_; k++) {
    for (int l = 0; l < cols_; l++) {
      int minorRow = 0;
      int minorCol = 0;
      double determinant = 0.0;

      for (int i = 0; i < rows_; i++) {
        if (i == k) {
          continue;
        }

        minorCol = 0;

        for (int j = 0; j < cols_; j++) {
          if (j == l) {
            continue;
          }

          minor.matrix_[minorRow][minorCol] = matrix_[i][j];
          minorCol++;
        }

        minorRow++;
      }

      determinant = minor.Determinant();
      result.matrix_[k][l] = pow(-1, (k + l)) * determinant;
    }
  }

  return result;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::invalid_argument("The matrix must be square");
  }

  double determinant = 0.0;
  int n = rows_;

  if (n == 1) {
    determinant = matrix_[0][0];
  } else if (n == 2) {
    determinant = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    S21Matrix tempMatrix(n - 1, n - 1);

    for (int k = 0; k < n; k++) {
      for (int i = 1; i < n; i++) {
        int t = 0;

        for (int j = 0; j < n; j++) {
          if (j == k) {
            continue;
          }

          tempMatrix.matrix_[i - 1][t] = matrix_[i][j];
          t++;
        }
      }

      determinant += pow(-1, k + 2) * matrix_[0][k] * tempMatrix.Determinant();
    }
  }

  return determinant;
}

S21Matrix S21Matrix::InverseMatrix() {
  double determinant = S21Matrix::Determinant();

  if (comparator(determinant, 0)) {
    throw std::invalid_argument("The determinant of the matrix is ​​0");
  }

  S21Matrix calcCompMatrix = S21Matrix::CalcComplements();
  S21Matrix transposeMatrix = calcCompMatrix.Transpose();
  transposeMatrix.MulNumber(1 / determinant);

  return transposeMatrix;
}

S21Matrix operator+(S21Matrix &matrix1, S21Matrix &matrix2) {
  if (matrix1.rows_ != matrix2.rows_ || matrix1.cols_ != matrix2.cols_) {
    throw std::invalid_argument{"Matrix sizes are not equal"};
  }

  S21Matrix result(matrix1.rows_, matrix1.cols_);

  for (int i = 0; i < matrix1.rows_; i++) {
    for (int j = 0; j < matrix1.cols_; j++) {
      result.matrix_[i][j] = matrix1.matrix_[i][j] + matrix2.matrix_[i][j];
    }
  }

  return result;
}

S21Matrix operator-(S21Matrix &matrix1, S21Matrix &matrix2) {
  if (matrix1.rows_ != matrix2.rows_ || matrix1.cols_ != matrix2.cols_) {
    throw std::invalid_argument{"Matrix sizes are not equal"};
  }

  S21Matrix result(matrix1.rows_, matrix1.cols_);

  for (int i = 0; i < matrix1.rows_; i++) {
    for (int j = 0; j < matrix1.cols_; j++) {
      result.matrix_[i][j] = matrix1.matrix_[i][j] - matrix2.matrix_[i][j];
    }
  }

  return result;
}

S21Matrix operator*(S21Matrix &matrix1, S21Matrix &matrix2) {
  if (matrix1.cols_ != matrix2.rows_) {
    throw std::invalid_argument(
        "The number of columns of the first matrix is not equal to the "
        "number of rows of the second matrix");
  }

  S21Matrix result(matrix1.rows_, matrix2.cols_);

  for (int i = 0; i < matrix1.rows_; i++) {
    for (int j = 0; j < matrix2.cols_; j++) {
      for (int k = 0; k < matrix1.cols_; k++) {
        result.matrix_[i][j] += matrix1.matrix_[i][k] * matrix2.matrix_[k][j];
      }
    }
  }

  return result;
}

S21Matrix operator*(double number, S21Matrix &matrix) {
  S21Matrix result(matrix.rows_, matrix.cols_);

  for (int i = 0; i < matrix.rows_; i++) {
    for (int j = 0; j < matrix.cols_; j++) {
      result.matrix_[i][j] = number * matrix.matrix_[i][j];
    }
  }

  return result;
}

S21Matrix operator*(S21Matrix &matrix, double number) {
  S21Matrix result(matrix.rows_, matrix.cols_);

  for (int i = 0; i < matrix.rows_; i++) {
    for (int j = 0; j < matrix.cols_; j++) {
      result.matrix_[i][j] = number * matrix.matrix_[i][j];
    }
  }

  return result;
}

bool operator==(S21Matrix &matrix1, S21Matrix &matrix2) {
  bool result = matrix1.EqMatrix(matrix2);

  return result;
}

S21Matrix operator+=(S21Matrix &matrix1, S21Matrix &matrix2) {
  matrix1.SumMatrix(matrix2);

  return matrix1;
}

S21Matrix operator-=(S21Matrix &matrix1, S21Matrix &matrix2) {
  matrix1.SubMatrix(matrix2);

  return matrix1;
}

S21Matrix operator*=(S21Matrix &matrix1, S21Matrix &matrix2) {
  matrix1.MulMatrix(matrix2);

  return matrix1;
}

S21Matrix operator*=(S21Matrix &matrix, double number) {
  matrix.MulNumber(number);

  return matrix;
}

S21Matrix operator*=(double number, S21Matrix &matrix) {
  matrix.MulNumber(number);

  return matrix;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &matrix) {
  if (this == &matrix) {
    return *this;
  }

  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
  }

  rows_ = matrix.rows_;
  cols_ = matrix.cols_;

  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix.matrix_[i][j];
    }
  }

  return *this;
}

double &S21Matrix::operator()(int row, int col) {
  if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
    throw std::out_of_range("Index is out of range");
  }

  return matrix_[row][col];
}
