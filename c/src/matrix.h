#ifndef MATRIX_H
#define MATRIX

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
public:
	 // initialize c x r matrix with the null element
	Matrix(int c, int r){
		this->columns = c;
		this->rows = r;
		elements.assign(this->rows*this->columns,SR::null());
	}
	Matrix(int c, int r, SR elem){ // initialize c x r matrix with elem
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
					tmp = tmp + this->elements.at(r*this->columns+i) * mat.elements.at(i*mat.rows+c);
				}
				ret.at(r*mat.columns+c) = tmp;
			}
		}
		return Matrix(this->rows, mat.columns, ret);
	};
	Matrix star ()
	{return *this;};
	std::string string()
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

template <typename SR>
std::ostream& operator<<(std::ostream& os, Matrix<SR>& matrix)
{
	return os << matrix.string();
}
#endif
