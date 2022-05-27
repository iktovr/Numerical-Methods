/*
xy'' + 2y' - xy=0
y(1) = 1 / e
y(2) = 0,5 / e^2
x in [1, 2]

y' = z
z' = -2z / x + y

alpha1 = 0, beta1 = 1, gamma1 = 1 / e
alpha2 = 0, beta2 = 1, gamma2 = 0,5 / e^2

y'' + 2 / x * y' - y + 0 = 0
p(x) = 2 / x
q(x) = -1
r(x) = 0
*/

#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>

const double e = std::exp(1);

#include "../include/differential/boundary_value_problem.hpp"

int main(int argc, char* argv[]) {
	bool plot = false;
	if (argc == 2 && std::strcmp(argv[1], "plot") == 0) {
		plot = true;
	}

	function_t<double> f = [](double, double, double z){ return z; };
	function_t<double> g = [](double x, double y, double z){ return -2 * z / x + y; };

	std::function<double(double)> p = [](double x){ return 2 / x; };
	std::function<double(double)> q = [](double){ return -1; };
	std::function<double(double)> r = [](double){ return 0; };

	double start = 1, end = 2, alpha1 = 0, beta1 = 1, gamma1 = 1 / e, alpha2 = 0, beta2 = 1, gamma2 = 0.5 / (e * e), h;
	double y1 = 1 / e, y2 = 0.5 / (e * e);
	std::cin >> h;

	std::vector<std::vector<double>> shooting = ShootingMethod(f, g, y1, y2, start, end, h);
	std::vector<std::vector<double>> finite_difference = FiniteDifferenceMethod(p, q, r, alpha1, beta1, gamma1, alpha2, beta2, gamma2, start, end, h);

	if (!plot) {
		std::cout << "Метод стрельбы:\n";
		for (size_t i = 0; i < shooting[0].size(); ++i) {
			std::cout << shooting[1][i] << ' ';
		}
		std::cout << "\nПогрешность: " << RungeRombergError(shooting[1], ShootingMethod(f, g, y1, y2, start, end, 2 * h)[1], 2) << "\n\n";

		std::cout << "Метод конечных разностей:\n";
		for (size_t i = 0; i < finite_difference[0].size(); ++i) {
			std::cout << finite_difference[1][i] << ' ';
		}
		std::cout << "\nПогрешность: " << RungeRombergError(finite_difference[1], FiniteDifferenceMethod(p, q, r, alpha1, beta1, gamma1, alpha2, beta2, gamma2, start, end, 2 * h)[1], 2) << "\n";

	} else {
		std::cout << "set size ratio -1\nset key off\n";
		const auto [xmin, xmax] = std::minmax_element(shooting[0].begin(), shooting[0].end());
		const auto [ymin, ymax] = std::minmax_element(shooting[1].begin(), shooting[1].end());
		double dx = (*xmax - *xmin) * 0.1, dy = (*ymax - *ymin) * 0.5;
		std::cout << "set xrange [" << (*xmin-dx) << ':' << (*xmax+dx) << "]\nset yrange [" << (*ymin-dy) << ':' << (*ymax+dy) << "]\n";
		
		std::cout << "$shooting << EOD\n\n";
		for (size_t i = 0; i < shooting[0].size(); ++i) {
			std::cout << shooting[0][i] << ' ' << shooting[1][i] << '\n';
		}
		std::cout << "EOD\n\n";

		std::cout << "$finite_difference << EOD\n\n";
		for (size_t i = 0; i < finite_difference[0].size(); ++i) {
			std::cout << finite_difference[0][i] << ' ' << finite_difference[1][i] << '\n';
		}
		std::cout << "EOD\n\n";

		std::cout << "f(x) = exp(-x) / x\n";

		std::cout << "\nset title 'Shooting method'\n"
		          << "plot f(x) with lines ls 4 lw 1.5, $shooting with " << ((shooting[0].size() <= 11) ? "linespoints ls 22 pt 22 " : "lines ls 22 ") << "lw 1.5 dt 3\n";

		std::cout << "\n# set title 'Finite difference method'\n"
		          << "# plot f(x) with lines ls 4 lw 1.5, $finite_difference with " << ((finite_difference[0].size() <= 11) ? "linespoints ls 22 pt 22 " : "lines ls 22 ") << "lw 1.5 dt 3\n";
	}
}