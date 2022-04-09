#include <functional>
#include <iostream>
#include <cmath>

#include "../include/function/polynomial.hpp"
#include "../include/function/interpolation_polynomial.hpp"

int main() {
	std::cout.precision(4);
	std::cout << std::fixed;

	std::function<double(double)> f = [](double x){ return std::cos(x); };
	Polynomial<double> p;

	int n;
	std::cin >> n;

	std::vector<double> x(n), y(n);
	for (int i = 0; i < n; ++i) {
		std::cin >> x[i];
		y[i] = f(x[i]);
	}

	double x_test;
	std::cin >> x_test;

	std::cout << "Многочлен в форме Лагранжа:\n";
	p = LagrangePolynomial(x, y);
	std::cout << p << '\n';

	std::cout << "Многочлен в форме Ньютона:\n";
	p = NewtonPolynomial(x, y);
	std::cout << p << '\n';

	std::cout << "Погрешность интерполяции:\n" << std::abs(f(x_test) - p(x_test)) << '\n';
}