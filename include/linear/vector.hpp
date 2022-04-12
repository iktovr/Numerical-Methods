#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <initializer_list>
#include <stdexcept>

template <class T>
class Vector {
private:
	std::vector<T> _data;
	size_t _size;

	Vector() : _data(), _size(0) { }

public:
	Vector(size_t n) : _data(n), _size(n) { }

	Vector(size_t n, T value) : _data(n, value), _size(n) { }

	Vector(std::initializer_list<T> list) : _data(list), _size(_data.size()) { }

	T& operator[](size_t index) {
		return _data[index];
	}

	const T& operator[](size_t index) const {
		return _data[index];
	}

	size_t Size() const {
		return _size;
	}

	T Norm(int p = 0) {
		T norm = 0;
		if (p > 0) {
			for (const T& a: _data) {
				norm += std::pow(a, p);
			}
			norm = std::pow(norm, 1.0 / p);
		} else {
			for (const T& a: _data) {
				norm = std::max(norm, std::abs(a));
			}
		}
		return norm;
	}

	Vector<T>& operator+=(const Vector<T>& other) {
		if (_size != other.Size()) {
			throw std::runtime_error("Dimension mismatch");
		}
		for (size_t i = 0; i < _size; ++i) {
			_data[i] += other[i];
		}
		return *this;
	}

	Vector<T>& operator*=(const T& value) {
		for (size_t i = 0; i < _size; ++i) {
			_data[i] *= value;
		}
		return *this;
	}

	template <class U>
	friend std::ostream& operator<<(std::ostream&, const Vector<U>&);

	template <class U>
	friend std::istream& operator>>(std::istream&, Vector<U>&);

	template<class U>
	friend Vector<U> operator-(const Vector<U>&, const Vector<U>&);

	template<class U>
	friend Vector<U> operator+(const Vector<U>&, const Vector<U>&);

	template<class U>
	friend Vector<U> operator-(const Vector<U>&);

	template<class U>
	friend U operator*(const Vector<U>&, const Vector<U>&);

	template<class U>
	friend Vector<U> operator*(const Vector<U>&, const U&);
};


template <class T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& vec) {
	for (size_t i = 0; i < vec.Size(); ++i) {
		os << std::setw(8) << vec[i] << ' ';
	}
	return os;
}

template <class T>
std::istream& operator>>(std::istream& is, Vector<T>& vec) {
	for (size_t i = 0; i < vec.Size(); ++i) {
		is >> vec[i];
	}
	return is;
}

template<class T>
Vector<T> operator-(const Vector<T>& a, const Vector<T>& b) {
	if (a.Size() != b.Size()) {
		throw std::runtime_error("Dimension mismatch");
	}
	Vector<T> res(a);
	for (size_t i = 0; i < b.Size(); ++i) {
		res[i] -= b[i];
	}
	return res;
}

template<class T>
Vector<T> operator+(const Vector<T>& a, const Vector<T>& b) {
	if (a.Size() != b.Size()) {
		throw std::runtime_error("Dimension mismatch");
	}
	Vector<T> res(a);
	for (size_t i = 0; i < b.Size(); ++i) {
		res[i] += b[i];
	}
	return res;
}

template<class T>
Vector<T> operator-(const Vector<T>& a) {
	Vector<T> res(a);
	for (size_t i = 0; i < a.Size(); ++i) {
		res[i] *= -1;
	}
	return res;
}

template<class T>
T operator*(const Vector<T>& a, const Vector<T>& b) {
	if (a.Size() != b.Size()) {
		throw std::runtime_error("Dimension mismatch");
	}
	T res = 0;
	for (size_t i = 0; i < b.Size(); ++i) {
		res += a[i] * b[i];
	}
	return res;
}

template<class T>
Vector<T> operator*(const Vector<T>& a, const T& b) {
	Vector<T> res(a);
	for (size_t i = 0; i < a.Size(); ++i) {
		res[i] *= b;
	}
	return res;
}