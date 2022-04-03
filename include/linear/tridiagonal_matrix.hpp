#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "vector.hpp"

template <class T>
class TDMatrix {
private:
	std::vector<std::vector<T>> _data;
	size_t _size;

public:
	TDMatrix(size_t n) : _data(n, std::vector<T>(3)), _size(n) { }

	std::vector<T>& operator[](size_t index) {
		return _data[index];
	}

	const std::vector<T>& operator[](size_t index) const {
		return _data[index];
	}

	size_t Size() const {
		return _size;
	}

	Vector<T> Solve(const Vector<T>& b) {
		if (b.Size() != _size) {
			throw std::runtime_error("Dimension mismatch");
		}

		Vector<T> p(_size-1), q(_size-1), x(_size);
		p[0] = - _data[0][2] / _data[0][1];
		q[0] = b[0] / _data[0][1];

		for (size_t i = 1; i < _size-1; ++i) {
			p[i] = - _data[i][2] / (_data[i][1] + _data[i][0] * p[i-1]);
			q[i] = (b[i] - _data[i][0] * q[i-1]) / (_data[i][1] + _data[i][0] * p[i-1]);
		}

		x[_size-1] = (b[_size-1] - _data[_size-1][0] * q[_size-2]) / (_data[_size-1][1] + _data[_size-1][0] * p[_size-2]);
		for (int i = _size - 2; i >= 0; --i) {
			x[i] = p[i] * x[i+1] + q[i];
		}

		return x;
	}

	template <class U>
	friend std::ostream& operator<<(std::ostream&, const TDMatrix<U>&);
	template <class U>
	friend std::istream& operator>>(std::istream&, TDMatrix<U>&);
};

template <class T>
std::ostream& operator<<(std::ostream& os, const TDMatrix<T>& matrix) {
	os << matrix[0][1] << ' ' << matrix[0][2] << '\n';
	for (size_t i = 1; i < matrix.Size()-1; ++i) {
		for (size_t j = 0; j < 3; ++j) {
			os << std::setw(8) << matrix[i][j] << ' ';
		}
		os << '\n';
	}
	os << matrix[matrix.Size()-1][0] << ' ' << matrix[matrix.Size()-1][1] << '\n';
	return os;
}

template <class T>
std::istream& operator>>(std::istream& is, TDMatrix<T>& matrix) {
	is >> matrix[0][1] >> matrix[0][2];
	for (size_t i = 1; i < matrix.Size()-1; ++i) {
		for (size_t j = 0; j < 3; ++j) {
			is >> matrix[i][j];
		}
	}
	is >> matrix[matrix.Size()-1][0] >> matrix[matrix.Size()-1][1];
	return is;
}
