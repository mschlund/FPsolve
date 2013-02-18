#pragma once

#include <cassert>
#include <initializer_list>
#include <sstream>
#include <string>
#include <vector>

template <typename SR>
class Matrix {
  public:
    Matrix(const Matrix &m) = default;
    Matrix(Matrix &&m) = default;

    Matrix(std::size_t r, const std::vector<SR> &es)
        : rows_(r), columns_(es.size() / r), elements_(es) {
      assert(Sanity());
    }

    Matrix(std::size_t r, std::vector<SR> &&es)
      : rows_(r), columns_(es.size() / r), elements_(std::move(es)) {
      assert(Sanity());
    }

    Matrix(std::size_t r, std::initializer_list<SR> es)
        : rows_(r), columns_(es.size() / r), elements_(es) {
      assert(Sanity());
    }

    Matrix(std::size_t r, std::size_t c)
        : rows_(r), columns_(c), elements_(rows_ * columns_, SR::null()) {}

    Matrix(std::size_t r, std::size_t c, const SR &elem)
        : rows_(r), columns_(c), elements_(rows_ * columns_, elem) {}

    inline std::size_t GetIndex(std::size_t r, std::size_t c) const {
      return r * columns_ + c;
    }

    inline SR& At(std::size_t r, std::size_t c) {
      assert(GetIndex(r, c) < elements_.size());
      return elements_[GetIndex(r, c)];
    }

    inline const SR& At(std::size_t r, std::size_t c) const {
      assert(GetIndex(r, c) < elements_.size());
      return elements_[GetIndex(r, c)];
    }

    Matrix& operator=(const Matrix &matrix) = default;
    Matrix& operator=(Matrix &&matrix) = default;

    Matrix operator+(const Matrix &mat) const {
      assert(rows_ == mat.rows_ && columns_ == mat.columns_ &&
             elements_.size() == mat.elements_.size());
      std::vector<SR> result;
      result.reserve(elements_.size());
      for (std::size_t i = 0; i < columns_ * rows_; ++i) {
        result.emplace_back(elements_[i] + mat.elements_[i]);
      }

      return Matrix{rows_, std::move(result)};
    };

    Matrix operator*(const Matrix &rhs) const {
      assert(columns_ == rhs.rows_);
      Matrix result{rows_, rhs.columns_, SR::null()};
      for (std::size_t r = 0; r < rows_; ++r) {
        for (std::size_t c = 0; c < rhs.columns_; ++c) {
          for (std::size_t i = 0; i < columns_; ++i) {
            assert(GetIndex(r, i) < elements_.size());
            assert(rhs.GetIndex(i, c) < rhs.elements_.size());
            assert(result.GetIndex(r, c) < result.elements_.size());
            result.At(r, c) += At(r, i) * rhs.At(i, c);
          }
        }
      }
      return result;
    };

    bool operator==(const Matrix &rhs) {
      assert(rows_ == rhs.rows_ && columns_ == rhs.columns_ &&
             elements_.size() == rhs.elements_.size());
      return elements_ == rhs.elements_;
    }

    Matrix FloydWarshall() const {
      assert(columns_ == rows_);
      Matrix result = *this;
      for (std::size_t k = 0; k < rows_; ++k) {
        for (std::size_t i = 0; i < rows_; ++i) {
          for (std::size_t j = 0; j < rows_; ++j) {
            result.At(i, j) +=
              result.At(i, k) * result.At(k, k).star() * result.At(k, j);
          }
        }
      }
      /* Add element 1, i.e., Floyd-Warshall will give us A+ matrix, so
       *   1 + A+ = A*
       * Reusing one() and operator+ would be much slower (additional
       * allocations and unnecessary traersals). */
      for (std::size_t i = 0; i < rows_; ++i) {
        result.At(i, i) += SR::one();
      }
      return result;
    }

    Matrix star2() const {
      assert(columns_ == rows_);
      return FloydWarshall();
    }

    Matrix star() const {
      assert(columns_ == rows_);
      return recursive_star(*this);
    }

    std::size_t getRows() const {
      return rows_;
    };

    std::size_t getColumns() const {
      return columns_;
    };

    const std::vector<SR>& getElements() const {
      return elements_;
    };

    std::string string() const {
      std::stringstream ss;
      for (std::size_t r = 0; r < rows_; ++r) {
        for (std::size_t c = 0; c < columns_; ++c) {
          ss << At(r, c) << " | ";
        }
        ss << std::endl;
      }
      return ss.str();
    }

    static Matrix const null(std::size_t size) {
      return Matrix{size, size, SR::null()};
    }

    static Matrix const one(std::size_t size) {
      std::vector<SR> result;
      result.reserve(size * size);
      for (std::size_t i = 0; i < size * size; ++i) {
        /* Diagonal entry, [0, size + 1, 2 * (size + 1), ...]. */
        if (i % (size + 1) == 0)
          result.emplace_back(SR::one());
        else
          result.emplace_back(SR::null());
      }
      return Matrix{size, std::move(result)};
    }

  private:
    std::size_t rows_;
    std::size_t columns_;
    std::vector<SR> elements_;

    bool Sanity() const {
      return columns_ * rows_ == elements_.size();
    }

    // this is a naive implementation which creates lots of matrices
    // maybe we can work directly on the elements_ (but then multiplication etc. are harder)
    static Matrix recursive_star(Matrix matrix) {
      assert(matrix.rows_ == matrix.columns_);
      if (matrix.rows_ == 1) {
        // just a scalar in a matrix
        matrix.elements_[0] = matrix.elements_[0].star(); // use semiring-star
        return matrix;
      }
      // peel mode if n%2 != 0, split in middle otherwise
      std::size_t split = matrix.rows_%2 == 0 ? matrix.columns_/2 : matrix.columns_-1;

      // TODO: this might be slightly inefficient
      Matrix a_11 = matrix.submatrix(0,split,0,split);
      Matrix a_12 = matrix.submatrix(split,matrix.columns_,0,split);
      Matrix a_21 = matrix.submatrix(0,split,split,matrix.rows_);
      Matrix a_22 = matrix.submatrix(split,matrix.columns_,split,matrix.rows_);
      Matrix as_11 = recursive_star(a_11);
      Matrix as_22 = recursive_star(a_22);
      Matrix A_11 = recursive_star(a_11 + a_12 * as_22 * a_21);
      Matrix A_22 = recursive_star(a_22 + a_21 * as_11 * a_12);
      Matrix A_12 = as_11 * a_12 * A_22;
      Matrix A_21 = as_22 * a_21 * A_11;
      // FIXME: should be:
      // Matrix A_12 = as_11 * a_12 * as_22;
      // Matrix A_21 = as_22 * a_21 * as_11;
      return block_matrix(A_11,A_12,A_21,A_22);
    }

    static Matrix block_matrix(const Matrix &a_11, const Matrix &a_12,
                               const Matrix &a_21, const Matrix &a_22) {
      std::vector<SR> ret;
      assert(a_11.rows_ == a_12.rows_ && a_21.rows_ == a_22.rows_);
      assert(a_11.columns_ == a_21.columns_ && a_12.columns_ == a_22.columns_);
      for (int r=0; r < a_11.rows_; r++) {
        ret.insert(     ret.end(),
            a_11.elements_.begin()+(r*a_11.columns_),
            a_11.elements_.begin()+((r+1)*a_11.columns_));
        ret.insert(     ret.end(),
            a_12.elements_.begin()+(r*a_12.columns_),
            a_12.elements_.begin()+((r+1)*a_12.columns_));
      }
      for (int r=0; r < a_21.rows_; r++) {
        ret.insert(     ret.end(),
            a_21.elements_.begin()+(r*a_21.columns_),
            a_21.elements_.begin()+((r+1)*a_21.columns_));
        ret.insert(     ret.end(),
            a_22.elements_.begin()+(r*a_22.columns_),
            a_22.elements_.begin()+((r+1)*a_22.columns_));
      }
      return Matrix{a_11.rows_+a_21.rows_,
                    std::move(ret)};
    }

    // get the submatrix starting from colum cs,...
    Matrix submatrix(std::size_t cs, std::size_t ce,
                     std::size_t rs, std::size_t re) const {
      assert(cs >= 0 && cs < columns_ && ce <= columns_ && ce > cs);
      assert(rs >= 0 && rs < rows_ && re <= rows_ && re > rs);
      // std::size_t nc = ce-cs; // new column count
      std::size_t nr = re-rs; // new row count
      std::vector<SR> ret;
      for (std::size_t r = rs; r < re; r++) {
        for (std::size_t c = cs; c < ce; c++) {
          // copy the needed values from elements_ to ret
          //ret[nc*(r-rs)+cs-c] = elements_[columns_*r+c];
          ret.push_back(elements_[columns_*r+c]);
        }
      }
      return Matrix{nr, std::move(ret)};
    }


};

template <typename SR>
Matrix<SR> operator*(SR elem, const Matrix<SR> &mat) {
  std::vector<SR> ret;
  for (std::size_t i = 0; i < mat.rows_ * mat.columns_; ++i) {
    ret.push_back(elem * mat[i]); // semiring multiplication
  }
}

template <typename SR>
std::ostream& operator<<(std::ostream &os, const Matrix<SR> &matrix) {
  return os << matrix.string();
}
