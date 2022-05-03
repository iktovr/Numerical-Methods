#pragma once

#include <vector>
#include <iostream>

#include "polynomial.hpp"
#include "../linear/tridiagonal_matrix.hpp"
#include "../linear/vector.hpp"

template <class T>
class CubicSpline {
public:
	std::vector<T> x;
	std::vector<T> y;
	std::vector<Polynomial<T>> segment;

	CubicSpline(const std::vector<T>& x, const std::vector<T>& y) : x(x), y(y), segment(x.size() - 1) {
		if (x.size() != y.size()) {
			throw std::runtime_error("Incompatible arrays");
		}

		std::vector<T> h(x.size() - 1);
		for (size_t i = 1; i < x.size(); ++i) {
			h[i-1] = x[i] - x[i-1];
		}

		size_t n = h.size();
		TDMatrix<T> matrix(n - 1);
		Vector<T> vec(n - 1);

		matrix[0][1] = 2 * (h[0] + h[1]);
		matrix[0][2] = h[1];
		vec[0] = 3 * ((y[2] - y[1]) / h[1] - (y[1] - y[0]) / h[0]);

		for (size_t i = 2; i < n - 1; ++i) {
			matrix[i-1][0] = h[i-2];
			matrix[i-1][1] = 2 * (h[i-2] + h[i-1]);
			matrix[i-1][2] = h[i-1];
			vec[i-1] = 3 * ((y[i+1] - y[i]) / h[i-1] - (y[i] - y[i-1]) / h[i-2]);
		}

		matrix[n-2][0] = h[n-2];
		matrix[n-2][1] = 2 * (h[n-2] + h[n-1]);
		vec[n-2] = 3 * ((y[n] - y[n-1]) / h[n-1] - (y[n-1] - y[n-2]) / h[n-2]);

		Vector<T> c2 = matrix.Solve(vec);
		std::vector<T> a(n), b(n), c(n), d(n);
		c[0] = 0;
		for (size_t i = 0; i < n-1; ++i) {
			a[i] = y[i];
			c[i+1] = c2[i];
			b[i] = (y[i+1] - y[i]) / h[i] - h[i] * (c[i+1] + 2 * c[i]) / 3.;
			d[i] = (c[i+1] - c[i]) / h[i] / 3.;
		}
		a[n-1] = y[n-1];
		b[n-1] = (y[n] - y[n-1]) / h[n-1] - 2. / 3. * h[n-1] * c[n-1];
		d[n-1] = - c[n-1] / (3 * h[n-1]);

		for (size_t i = 0; i < n; ++i) {
			Polynomial<T> xb{-x[i], 1}, xc = xb * Polynomial<T>{-x[i], 1}, xd = xc * Polynomial<T>{-x[i], 1};
			segment[i] = d[i] * xd + c[i] * xc + b[i] * xb + a[i];
		}
	}

	size_t Size() const {
		return segment.size();
	}

	const Polynomial<T>& operator[](size_t index) {
		return segment[index];
	}

	T operator()(const T& value) {
		if (value < x[0] || value > x.back()) {
			throw std::runtime_error("Out of range");
		}

		for (size_t i = 0; i < x.size()-1; ++i) {
			if (value >= x[i] && value <= x[i+1] ) {
				return segment[i](value);
			}
		}

		return 0; // just for silence the warning
	}
};

template <class T>
std::ostream& operator<<(std::ostream& os, const CubicSpline<T>& spline) {
	for (size_t i = 0; i < spline.x.size() - 1; ++i) {
		os << '[' << spline.x[i] << ':' << spline.x[i+1] << "] " << spline.segment[i] << '\n';
	}
	return os;
}