#include <functional>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>

#include "../include/function/lsm.hpp"
#include "../include/function/polynomial.hpp"

int main(int argc, char* argv[]) {
	bool plot = false;
	if (argc == 2 && std::strcmp(argv[1], "plot") == 0) {
		plot = true;
	}

	std::cout.precision(4);
	std::cout << std::fixed;

	std::vector<double> x, y;
	double a, b;
	while (std::cin >> a >> b) {
		x.push_back(a);
		y.push_back(b);
	}

	std::vector<std::function<double(double)>> basis{
		[](double) { return 1; },
		[](double x) { return x; }
	};
	Polynomial<double> p1(LSM(x, y, basis));

	basis.push_back(
		[](double x) { return x*x; }
	);
	Polynomial<double> p2(LSM(x, y, basis));

	basis.push_back(
		[](double x) { return x*x*x; }
	);
	Polynomial<double> p3(LSM(x, y, basis));

	if (!plot) {
		std::cout << "Приближающий многочлен 1-ой степени:\n" << p1 << '\n'
				  << "Ошибка: " << SquareError(std::function<double(double)>{p1}, x, y) << "\n\n"

				  << "Приближающий многочлен 2-ой степени:\n" << p2 << '\n'
				  << "Ошибка: " << SquareError(std::function<double(double)>{p2}, x, y) << "\n\n"

				  << "Приближающий многочлен 3-ой степени:\n" << p3 << '\n'
				  << "Ошибка: " << SquareError(std::function<double(double)>{p3}, x, y) << '\n';

	} else {
		std::cout << "set size ratio -1\n";
		const auto [xmin, xmax] = std::minmax_element(x.begin(), x.end());
		const auto [ymin, ymax] = std::minmax_element(y.begin(), y.end());
		double dx = (*xmax - *xmin) * 0.5, dy = (*ymax - *ymin);
		std::cout << "set xrange [" << (*xmin-dx) << ':' << (*xmax+dx) << "]\nset yrange [" << (*ymin-dy) << ':' << (*ymax+dy) << "]\n"
				  << "$points << EOD\n\n";
		for (size_t i = 0; i < x.size(); ++i) {
			std::cout << x[i] << ' ' << y[i] << '\n';
		}
		std::cout << "EOD\np1(x) = " << p1 << "\np2(x) = " << p2 << "\np3(x) = " << p3 << '\n'
				  << "plot $points ps 1.5 with points title 'points', p1(x) lw 1.5 title 'first degree', p2(x) lw 1.5 title 'second degree', p3(x) lw 1.5 title 'third degree'\n";
	}
}