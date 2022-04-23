#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>

#include "../include/function/derivation.hpp"
#include "../include/function/interpolation_polynomial.hpp"

int main(int argc, char* argv[]) {
	bool plot = false;
	if (argc == 2 && std::strcmp(argv[1], "plot") == 0) {
		plot = true;
	}

	std::cout.precision(4);
	std::cout << std::fixed;

	std::vector<double> x, y;
	double x0, a, b;
	std::cin >> x0;
	while (std::cin >> a >> b) {
		x.push_back(a);
		y.push_back(b);
	}

	double first_derivative = TableFirstDerivative(x, y, x0);
	double second_derivative = TableSecondDerivative(x, y, x0);

	if (!plot) {
		std::cout << "Первая производная: " << first_derivative << '\n';
		std::cout << "Вторая производная: " << second_derivative << '\n';

	} else {
		std::cout << "set size ratio -1\n";
		const auto [xmin, xmax] = std::minmax_element(x.begin(), x.end());
		const auto [ymin, ymax] = std::minmax_element(y.begin(), y.end());
		double dx = (*xmax - *xmin) * 0.2, dy = (*ymax - *ymin) * 0.5;
		std::cout << "set xrange [" << (*xmin-dx) << ':' << (*xmax+dx) << "]\nset yrange [" << (*ymin-dy) << ':' << (*ymax+dy) << "]\n"
				  << "$points << EOD\n\n";
		for (size_t i = 0; i < x.size(); ++i) {
			std::cout << x[i] << ' ' << y[i] << '\n';
		}
		std::cout << "EOD\n";
		size_t i = 0;
		double min = x0 - x[0];
		for (size_t j = 1; j < x.size(); ++j) {
			if (min > std::abs(x0 - x[j])) {
				min = std::abs(x0 - x[j]);
				i = j - 1;
			}
		}
		if (i == x.size() - 2) {
			--i;
		}
		Polynomial<double> p = NewtonPolynomial(std::vector<double>{x[i], x[i+1], x[i+2]}, std::vector<double>{y[i], y[i+1], y[i+2]});
		std::cout << "plot sample [" << x[i] << ':' << x[i+2] << "] " << p << " title 'polynomial' lw 1.5, $points title 'points', ";
		std::cout << '[' << *xmin << ':' << *xmax << "] " << p(x0) << '+' << first_derivative << "*(x-" << x0 << ") title 'tangent' lw 1.5 dt 2";
	}
}