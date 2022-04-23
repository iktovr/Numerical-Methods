#include <functional>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstring>

#include "../include/function/polynomial.hpp"
#include "../include/function/interpolation_polynomial.hpp"

int main(int argc, char* argv[]) {
	bool plot = false;
	if (argc == 2 && std::strcmp(argv[1], "plot") == 0) {
		plot = true;
	}

	std::cout.precision(4);
	std::cout << std::fixed;

	std::function<double(double)> f = [](double x){ return std::cos(x); };
	Polynomial<double> p1, p2;

	int n;
	std::cin >> n;

	std::vector<double> x(n), y(n);
	for (int i = 0; i < n; ++i) {
		std::cin >> x[i];
		y[i] = f(x[i]);
	}

	double x_test;
	std::cin >> x_test;

	p1 = LagrangePolynomial(x, y);
	p2 = NewtonPolynomial(x, y);

	if (!plot) {
		std::cout << "Многочлен в форме Лагранжа:\n" << p1 << '\n'
				  << "Многочлен в форме Ньютона:\n" << p2 << '\n'
				  << "Погрешность интерполяции:\n" << std::abs(f(x_test) - p1(x_test)) << '\n';

	} else {
		std::cout << "set size ratio -1\n";
		const auto [xmin, xmax] = std::minmax_element(x.begin(), x.end());
		const auto [ymin, ymax] = std::minmax_element(y.begin(), y.end());
		double dx = (*xmax - *xmin) * 0.1, dy = (*ymax - *ymin) * 0.1;
		std::cout << "set xrange [" << (*xmin-dx) << ':' << (*xmax+dx) << "]\nset yrange [" << (*ymin-dy) << ':' << (*ymax+dy) << "]\n"
				  << "$points << EOD\n\n";
		for (size_t i = 0; i < x.size(); ++i) {
			std::cout << x[i] << ' ' << y[i] << '\n';
		}
		std::cout << "EOD\np(x) = " << p1 << "\nf(x) = cos(x)\n"
				  << "plot $points with points ps 1.5 title 'points', f(x) title 'cos(x)' dt 2 lw 2, p(x) title 'polynomial' lw 2\n";
	}
}