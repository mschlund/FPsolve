#pragma once

#include <assert.h>
#include <initializer_list>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

template <typename SR>
class Matrix {
  private:
    std::size_t columns;
    std::size_t rows;
    std::vector<SR> elements;

    bool Sanity() const {
      return columns * rows == elements.size();
    }

    // this is a naive implementation which creates lots of matrices
    // maybe we can work directly on the elements (but then multiplication etc. are harder)
    static Matrix recursive_star(Matrix matrix) {
      assert(matrix.rows == matrix.columns);
      if (matrix.rows == 1) {
        // just a scalar in a matrix
        matrix.elements[0] = matrix.elements[0].star(); // use semiring-star
        return matrix;
      }
      // peel mode if n%2 != 0, split in middle otherwise
      std::size_t split = matrix.rows%2 == 0 ? matrix.columns/2 : matrix.columns-1;

      // TODO: this might be slightly inefficient
      Matrix a_11 = matrix.submatrix(0,split,0,split);
      Matrix a_12 = matrix.submatrix(split,matrix.columns,0,split);
      Matrix a_21 = matrix.submatrix(0,split,split,matrix.rows);
      Matrix a_22 = matrix.submatrix(split,matrix.columns,split,matrix.rows);
      Matrix as_11 = recursive_star(a_11);
      Matrix as_22 = recursive_star(a_22);
      Matrix A_11 = recursive_star(a_11 + a_12 * as_22 * a_21);
      Matrix A_22 = recursive_star(a_22 + a_21 * as_11 * a_12);
      Matrix A_12 = as_11 * a_12 * A_22;
      Matrix A_21 = as_22 * a_21 * A_11;
      return block_matrix(A_11,A_12,A_21,A_22);
    }

    static Matrix block_matrix(Matrix a_11, Matrix a_12, Matrix a_21, Matrix a_22) {
      std::vector<SR> ret;
      assert(a_11.rows == a_12.rows && a_21.rows == a_22.rows);
      assert(a_11.columns == a_21.columns && a_12.columns == a_22.columns);
      for (int r=0; r < a_11.rows; r++) {
        ret.insert(     ret.end(),
            a_11.elements.begin()+(r*a_11.columns),
            a_11.elements.begin()+((r+1)*a_11.columns));
        ret.insert(     ret.end(),
            a_12.elements.begin()+(r*a_12.columns),
            a_12.elements.begin()+((r+1)*a_12.columns));
      }
      for (int r=0; r < a_21.rows; r++) {
        ret.insert(     ret.end(),
            a_21.elements.begin()+(r*a_21.columns),
            a_21.elements.begin()+((r+1)*a_21.columns));
        ret.insert(     ret.end(),
            a_22.elements.begin()+(r*a_22.columns),
            a_22.elements.begin()+((r+1)*a_22.columns));
      }
      return Matrix(a_11.columns+a_12.columns, a_11.rows+a_21.rows, ret);
    }

    // get the submatrix starting from colum cs,...
    Matrix submatrix(std::size_t cs, std::size_t ce, std::size_t rs, std::size_t re) {
      assert(cs >= 0 && cs < columns && ce <= columns && ce > cs);
      assert(rs >= 0 && rs < rows && re <= rows && re > rs);
      std::size_t nc = ce-cs; // new column count
      std::size_t nr = re-rs; // new row count
      std::vector<SR> ret;
      for (std::size_t r = rs; r < re; r++) {
        for (std::size_t c = cs; c < ce; c++) {
          // copy the needed values from elements to ret
          //ret[nc*(r-rs)+cs-c] = elements[columns*r+c];
          ret.push_back(elements[columns*r+c]);
        }
      }
      return Matrix(nc, nr, std::move(ret));
    }


  public:
    Matrix(const Matrix &m) = default;
    Matrix(Matrix &&m) = default;

    Matrix(std::size_t c, std::size_t r, const std::vector<SR> &es)
        : columns(c), rows(r), elements(es) {
      assert(Sanity());
    }

    Matrix(std::size_t c, std::size_t r, std::vector<SR> &&es)
      : columns(c), rows(r), elements(std::move(es)) {
      assert(Sanity());
    }

    Matrix(std::size_t c, std::size_t r, std::initializer_list<SR> es)
        : columns(c), rows(r), elements(es) {
      assert(Sanity());
    }

    Matrix(std::size_t c, std::size_t r)
        : columns(c), rows(r), elements(columns * rows, SR::null()) {}

    Matrix(std::size_t c, std::size_t r, const SR &elem)
        : columns(c), rows(r), elements(columns * rows, elem) {}

    Matrix& operator=(const Matrix &matrix) = default;
    Matrix& operator=(Matrix &&matrix) = default;

    Matrix operator+(const Matrix &mat) const {
      assert(rows == mat.rows && columns == mat.columns &&
             elements.size() == mat.elements.size());
      std::vector<SR> result;
      result.reserve(elements.size());
      for (std::size_t i = 0; i < columns * rows; ++i) {
        result.emplace_back(elements[i] + mat.elements[i]);
      }

      return Matrix(columns, rows, std::move(result));
    };

    Matrix operator*(const Matrix &rhs) const {
      assert(columns == rhs.rows);
      // TODO: naive implementation, tune this
      std::vector<SR> result(rows * rhs.columns, SR::null());
      for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < rhs.columns; ++c) {
          SR tmp = SR::null();
          for (std::size_t i = 0; i < columns; ++i) {
            assert(i * rhs.columns + c < rhs.elements.size());
            assert(r * columns + i < elements.size());
            tmp += elements[r * columns + i] * rhs.elements[i * rhs.columns + c];
          }
          assert(r * rhs.columns + c < result.size());
          result[r * rhs.columns + c] = tmp;
        }
      }
      return Matrix(rhs.columns, rows, std::move(result));
    };

    bool operator==(const Matrix &rhs) {
      assert(rows == rhs.rows && columns == rhs.columns &&
             elements.size() == rhs.elements.size());
      return elements == rhs.elements;
    }

    Matrix star() const {
      assert(columns == rows);
      return recursive_star(*this);
    };

    Matrix transpose() const {
      std::vector<SR> result;
      result.reserve(elements.size());
      for (std::size_t c = 0; c < columns; ++c) {
        for (std::size_t r = 0; r < rows; ++r) {
          result.push_back(elements[r * columns + c]);
        }
      }
      return Matrix(rows, columns, std::move(result));
    };

    std::size_t getRows() const {
      return rows;
    };

    std::size_t getColumns() const {
      return columns;
    };

    const std::vector<SR>& getElements() const {
      return elements;
    };

    std::string string() const {
      std::stringstream ss;
      for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < columns; ++c) {
          ss << elements[r * columns + c] << " ";
        }
        ss << std::endl;
      }
      return ss.str();
    }

    static Matrix<SR> const null(std::size_t size) {
      return Matrix(size, size, SR::null());
    }

    static Matrix<SR> const one(std::size_t size) {
      std::vector<SR> result;
      result.reserve(size * size);
      for (std::size_t i = 0; i < size * size; ++i) {
        /* Diagonal entry, [0, size + 1, 2 * (size + 1), ...]. */
        if (i % (size + 1) == 0)
          result.emplace_back(SR::one());
        else
          result.emplace_back(SR::null());
      }
      return Matrix(size, size, std::move(result));
    }
};

template <typename SR>
Matrix<SR> operator*(SR elem, const Matrix<SR> &mat) {
  std::vector<SR> ret;
  for (std::size_t i = 0; i < mat.rows * mat.columns; ++i) {
    ret.push_back(elem * mat[i]); // semiring multiplication
  }
}

template <typename SR>
std::ostream& operator<<(std::ostream &os, const Matrix<SR> &matrix) {
  return os << matrix.string();
}
