#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <sstream>
#include <vector>
#include <assert.h>
#include <initializer_list>
#include <memory>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace ub = boost::numeric::ublas;

template <typename SR>
class Matrix
{
private:
	ub::matrix<SR> m;

	inline static ub::matrix<SR> topLeft(ub::matrix<SR>& matrix, int split) {return ub::subrange(matrix, 0,split, 0,split);}
	inline static ub::matrix<SR> topRight(ub::matrix<SR>& matrix, int split) {return ub::subrange(matrix, 0,split, split,matrix.size2());}
	inline static ub::matrix<SR> bottomLeft(ub::matrix<SR>& matrix, int split) {return ub::subrange(matrix, split,matrix.size1(), 0,split);}
	inline static ub::matrix<SR> bottomRight(ub::matrix<SR>& matrix, int split) {return ub::subrange(matrix, split,matrix.size1(), split,matrix.size2());}

	static ub::matrix<SR> recursive_star(ub::matrix<SR> matrix)
	{
		assert(matrix.size1() == matrix.size2());
		if(matrix.size1() == 1) // just a scalar
		{
			auto ret = ub::matrix<SR>(1,1);
			ret(0,0) = matrix(0,0).star();
			return ret;
		}
		// peel mode if n%2 != 0, split in middle otherwise
		int split = matrix.size1()%2 == 0 ? matrix.size2()/2 : matrix.size2()-1;

		ub::matrix_range<ub::matrix<SR> > a_11 (matrix, ub::range(0, split), ub::range(0, split));
		ub::matrix_range<ub::matrix<SR> > a_12 (matrix, ub::range(0, split), ub::range(split, matrix.size2()));
		ub::matrix_range<ub::matrix<SR> > a_21 (matrix, ub::range(split, matrix.size1()), ub::range(0, split));
		ub::matrix_range<ub::matrix<SR> > a_22 (matrix, ub::range(split, matrix.size1()), ub::range(split, matrix.size2()));
		ub::matrix<SR> as_11 = recursive_star(a_11);
		ub::matrix<SR> as_22 = recursive_star(a_22);

		// create result matrix and fill it
		ub::matrix<SR> ret(matrix.size1(), matrix.size2());
		// TODO: implement topLeft etc. so that the result can be used on the left hand side
		// at the moment, using them as lhs seem to work on a copy. Try to return a reference somehow

		//topLeft(ret, split) = recursive_star(a_11 + ub::prod(ub::matrix<SR>(ub::prod(a_12, as_22)), a_21));
		ub::subrange(ret, 0,split, 0,split) = recursive_star(a_11 + ub::prod(ub::matrix<SR>(ub::prod(a_12, as_22)), a_21));
		//bottomRight(ret, split) = recursive_star(a_22 + ub::prod(ub::matrix<SR>(ub::prod(a_21, as_11)), a_12));
		ub::subrange(ret, split,matrix.size1(), split,matrix.size2()) = recursive_star(a_22 + ub::prod(ub::matrix<SR>(ub::prod(a_21, as_11)), a_12));
		//topRight(ret, split) = ub::prod(as_11, ub::prod<ub::matrix<SR>>(a_12, bottomRight(ret,split)));
		ub::subrange(ret, 0,split, split,matrix.size2()) = ub::prod(as_11, ub::prod<ub::matrix<SR>>(a_12, bottomRight(ret,split)));
		//bottomLeft(ret, split) = ub::prod(as_22, ub::prod<ub::matrix<SR>>(a_21, topLeft(ret,split)));
		ub::subrange(ret, split,matrix.size1(), 0,split) = ub::prod(as_22, ub::prod<ub::matrix<SR>>(a_21, topLeft(ret,split)));
		return ret;
	}

public:
	Matrix(const Matrix& matrix)
	{
		this->m = matrix.m;
	}

	Matrix(int c, int r, std::initializer_list<SR> elements)
	{
		this->m = ub::matrix<SR>(r,c);
		int i_c = 0;
		int i_r = 0;
		for(auto element_it = elements.begin(); element_it != elements.end(); ++element_it)
		{
			this->m(i_r,i_c) = *element_it;
			if(++i_c >= c)
			{
				i_c = 0;
				++i_r;
			}
			if(i_r >= r) assert(false);
		}
	}

	Matrix& operator=(const Matrix& matrix)
	{
		this->m = matrix.m;
		return (*this);
	};

	 // initialize c x r matrix with the null element
	Matrix(int c, int r){
		this->m = ub::matrix<SR>(r,c);
		for(int i_r = 0; i_r < r; ++i_r)
		{
			for(int i_c = 0; i_c < c; ++i_c)
			{
				this->m(i_r,i_c) = SR::null();
			}
		}

	}
	Matrix(int c, int r, const SR& elem){ // initialize c x r matrix with elem
		this->m = ub::matrix<SR>(r,c);
		for(int i_r = 0; i_r < r; ++i_r)
		{
			for(int i_c = 0; i_c < c; ++i_c)
			{
				this->m(i_r,i_c) = elem;
			}
		}
	}
	Matrix(unsigned int c, unsigned int r, const std::vector<SR>& elements)
	{
		assert(c*r == elements.size());
		this->m = ub::matrix<SR>(r,c);
		for(int i_r = 0; i_r < r; ++i_r)
		{
			for(int i_c = 0; i_c < c; ++i_c)
			{
				this->m(i_r,i_c) = elements[i_r*c+i_c];
			}
		}
	}
	Matrix operator += (const Matrix& mat)
	{
		this->m = this->m + mat.m;
		return *this;
	};
	Matrix operator + (const Matrix& mat) const
	{
		Matrix result = *this;
		result += mat;
		return result;
	}
	Matrix operator *= (const Matrix& mat)
	{
		this->m = ub::prod(this->m,mat.m);
		return *this;
	};
	Matrix operator * (const Matrix& mat) const
	{
		Matrix result = *this;
		result *= mat;
		return result;
	}

	bool operator== (const Matrix& mat)
	{
		unsigned int rows1 = this->m.size1();
		unsigned int cols1 = this->m.size2();
		unsigned int rows2 = mat.m.size1();
		unsigned int cols2 = mat.m.size2();

		if(rows1 != rows2 || cols1 != cols2)
			return false;
		else {
			// if the dimensions match, do element-wise comparison
			for(unsigned int i=0; i<rows1; ++i) {
				for(unsigned int j=0; j<cols2; ++j)	{
					if(not(this->m(i,j) == mat.m(i,j)))
						return false;
				}
			}
		}

		//no difference found -> matrices are equal
		return true;
	}

	Matrix star ()
	{
		assert(this->m.size1() == this->m.size2());
		Matrix ret(*this);
		ret.m = recursive_star(this->m);
		return ret;
	};

	Matrix transpose() const
	{
		Matrix result = *this;
		result->m = ub::trans(this->m);
		return result;
	};

	int getRows() const
	{
		return this->m.size1();
	};

	int getColumns() const
	{
		return this->m.size2();
	};

	std::vector<SR> getElements() const
	{
		std::vector<SR> result;
		for(auto it = this->m.data().begin(); it != this->m.data().end(); ++it)
		{
			result.push_back(*it);
		}
		return result;
	};

	std::string string() const
	{
		std::stringstream ss;
		int r,c;
		int rows = this->getRows();
		int columns = this->getColumns();
		std::vector<SR> elements = this->getElements();
		for(r = 0; r < rows; r++)
		{
			for(c = 0; c < columns; c++)
			{
				ss << elements.at(r*columns+c) << " ";
			}
			ss << std::endl;
		}
		// TODO: use operator<< of ublas::matrix (in io.hpp)
		return ss.str();
	}

	static std::shared_ptr<Matrix<SR>> elem_null;
	static std::shared_ptr<Matrix<SR>> elem_one;

	static Matrix<SR> const null(int size)
	{
		if(!Matrix::elem_null)
			Matrix::elem_null = std::shared_ptr<Matrix<SR>>(new Matrix(size, size, SR::null()));
		return *Matrix::elem_null;
	}

	static Matrix<SR> const one(int size)
	{
		if(!Matrix::elem_one)
		{
			std::vector<SR> ret;
			for(int i = 0; i < size*size; ++i)
			{
				if(i%(size+1) == 0) // diagonal entry, [0,size+1,2*(size+1),...]
					ret.push_back(SR::one());
				else
					ret.push_back(SR::null());
			}
			Matrix::elem_one = std::shared_ptr<Matrix<SR>>(new Matrix(size, size, ret));
		}
		return *Matrix::elem_one;
	}
};

template <typename SR> std::shared_ptr<Matrix<SR>> Matrix<SR>::elem_null;
template <typename SR> std::shared_ptr<Matrix<SR>> Matrix<SR>::elem_one;

template <typename SR>
std::ostream& operator<<(std::ostream& os, const Matrix<SR>& matrix)
{
	return os << matrix.string();
}
#endif
