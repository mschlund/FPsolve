#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <sstream>
#include <vector>
#include <assert.h>

template <typename SR>
class Matrix
{
private:
	std::vector<SR> elements;
	int columns;
	int rows;

	// this is a naive implementation which creates lots of matrices
	// maybe we can work directly on the elements (but then multiplication etc. are harder)
	static Matrix recursive_star(Matrix matrix)
	{
		assert(matrix.rows == matrix.columns);
		if(matrix.rows == 1) // just a scalar in a matrix
		{
			matrix.elements[0] = matrix.elements[0].star(); // use semiring-star
			return matrix;
		}
		// peel mode if n%2 != 0, split in middle otherwise
		int split = matrix.rows%2 == 0 ? matrix.columns/2 : matrix.columns-1;

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

	static Matrix block_matrix(Matrix a_11, Matrix a_12, Matrix a_21, Matrix a_22)
	{
		std::vector<SR> ret;
		assert(a_11.rows == a_12.rows && a_21.rows == a_22.rows);
		assert(a_11.columns == a_21.columns && a_12.columns == a_22.columns);
		for(int r=0; r < a_11.rows; r++)
		{
			ret.insert(	ret.end(),
					a_11.elements.begin()+(r*a_11.columns),
					a_11.elements.begin()+((r+1)*a_11.columns));
			ret.insert(	ret.end(),
					a_12.elements.begin()+(r*a_12.columns),
					a_12.elements.begin()+((r+1)*a_12.columns));
		}
		for(int r=0; r < a_21.rows; r++)
		{
			ret.insert(	ret.end(),
					a_21.elements.begin()+(r*a_21.columns),
					a_21.elements.begin()+((r+1)*a_21.columns));
			ret.insert(	ret.end(),
					a_22.elements.begin()+(r*a_22.columns),
					a_22.elements.begin()+((r+1)*a_22.columns));
		}
		return Matrix(a_11.columns+a_12.columns, a_11.rows+a_21.rows, ret);
	}
public:
	Matrix(const Matrix& matrix)
	{
		this->columns = matrix.columns;
		this->rows = matrix.rows;
		this->elements = matrix.elements;
	}

	Matrix& operator=(const Matrix& matrix)
	{
		this->columns = matrix.columns;
		this->rows = matrix.rows;
		this->elements = matrix.elements;
		return (*this);
	};

	 // initialize c x r matrix with the null element
	Matrix(int c, int r){
		this->columns = c;
		this->rows = r;
		elements.assign(this->rows*this->columns,SR::null());
	}
	Matrix(int c, int r, const SR& elem){ // initialize c x r matrix with elem
		this->columns = c;
		this->rows = r;
		elements.assign(this->rows*this->columns,elem);
	}
	Matrix(int c, int r, const std::vector<SR>& elements)
	{
		assert(c*r == elements.size());
		this->columns = c;
		this->rows = r;
		this->elements = elements;
	}
	// get the submatrix starting from colum cs,...
	Matrix submatrix(int cs, int ce, int rs, int re)
	{
		assert(cs>=0 && cs<this->columns && ce <= this->columns && ce > cs);
		assert(rs>=0 && rs<this->rows && re <= this->rows && re > rs);
		int nc = ce-cs; // new column count
		int nr = re-rs; // new row count
		std::vector<SR> ret;
		for(int r=rs; r<re; r++)
		{
			for(int c=cs; c<ce; c++)
			{
				// copy the needed values from elements to ret
				//ret[nc*(r-rs)+cs-c] = this->elements[this->columns*r+c];
				ret.push_back(this->elements[this->columns*r+c]);
			}
		}
		return Matrix(nc, nr, ret);
	}
	Matrix operator + (const Matrix& mat)
	{
		assert(this->rows == mat.rows && this->columns == mat.columns);
		std::vector<SR> ret;
		for(int i=0; i<this->columns * this->rows; i++)
		{
			ret.push_back(this->elements.at(i)+mat.elements.at(i));
		}

		return Matrix(this->rows, this->columns, ret);
	};
	Matrix operator * (const Matrix& mat)
	{
		assert(this->columns == mat.rows);
		// TODO: naive implementation, tune this
		std::vector<SR> ret;
		ret.assign(this->rows*mat.columns,SR::null());
		int i,r,c;
		for(r = 0; r<this->rows; r++)
		{
			for(c = 0; c<mat.columns; c++)
			{
				SR tmp = SR::null();
				for(i = 0; i<this->columns; i++)
				{
					assert(i*mat.columns+c < mat.elements.size());
					assert(r*this->columns+i < this->elements.size());
					tmp = tmp + this->elements.at(r*this->columns+i) * mat.elements.at(i*mat.columns+c);
				}
				assert(r*mat.columns+c < ret.size());
				ret.at(r*mat.columns+c) = tmp;
			}
		}
		return Matrix(mat.columns, this->rows, ret);
	};

	Matrix star ()
	{
		assert(this->columns == this->rows);
		Matrix ret(*this);
		return recursive_star(ret);
	};

	Matrix transpose() const
	{
		std::vector<SR> ret;
		for(int c = 0; c<this->columns; c++)
		{
			for(int r = 0; r<this->rows; r++)
			{
				ret.push_back(this->elements.at(r*this->columns+c));
			}
		}
		return Matrix(this->rows, this->columns, ret);
	};

	int getRows() const
	{
		return this->rows;
	};

	int getColumns() const
	{
		return this->columns;
	};

	std::vector<SR> getElements() const
	{
		std::vector<SR> ret = this->elements;
		return ret;
	};

	std::string string() const
	{
		std::stringstream ss;
		int r;
		int c;
		for(r = 0; r < this->rows; r++)
		{
			for(c = 0; c < this->columns; c++)
			{
				ss << this->elements.at(r*this->columns+c) << " ";
			}
			ss << std::endl;
		}
		return ss.str();
	}
};

// friend method of Matrix
template <typename SR>
Matrix<SR> operator * (SR elem, const Matrix<SR>& mat)
{
	std::vector<SR> ret;
	for(int i=0; i<mat.rows*mat.columns; i++)
	{
		ret.push_back(elem*mat[i]); // semiring multiplication
	}
}

template <typename SR>
std::ostream& operator<<(std::ostream& os, const Matrix<SR>& matrix)
{
	return os << matrix.string();
}
#endif
