/*
y'' + y - 2cosx = 0
y(0) = 1
y'(0) = 0
x in [0, 1]
h = 0.1

y' = z
z' = 2cosx - y
*/

#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>

#include "../include/differential/cauchy_problem.hpp"

int main(int argc, char* argv[]) {
	bool plot = false;
	if (argc == 2 && std::strcmp(argv[1], "plot") == 0) {
		plot = true;
	}

	function_t<double> f = [](double, double, double z){ return z; };
	function_t<double> g = [](double x, double y, double){ return 2 * std::cos(x) - y; };

	double start, end, y0, z0, h;
	std::cin >> start >> end >> y0 >> z0 >> h;

	std::vector<std::vector<double>> euler = EulerMethod(f, g, y0, z0, start, end, h);
	std::vector<std::vector<double>> runge = RungeKuttaMethod(f, g, y0, z0, start, end, h);
	std::vector<std::vector<double>> adams = AdamsMethod(f, g, y0, z0, start, end, h);

	if (!plot) {
		std::cout << "Метод Эйлера:\n";
		for (size_t i = 0; i < euler[0].size(); ++i) {
			std::cout << euler[1][i] << ' ';
		}
		std::cout << "\nПогрешность: " << RungeRombergError(euler[1], EulerMethod(f, g, y0, z0, start, end, 2 * h)[1], 2) << "\n\n";

		std::cout << "Метод Рунге-Кутты:\n";
		for (size_t i = 0; i < runge[0].size(); ++i) {
			std::cout << runge[1][i] << ' ';
		}
		std::cout << "\nПогрешность: " << RungeRombergError(runge[1], RungeKuttaMethod(f, g, y0, z0, start, end, 2 * h)[1], 4) << "\n\n";

		std::cout << "Метод Адамса:\n";
		for (size_t i = 0; i < adams[0].size(); ++i) {
			std::cout << adams[1][i] << ' ';
		}
		std::cout << "\nПогрешность: " << RungeRombergError(adams[1], AdamsMethod(f, g, y0, z0, start, end, 2 * h)[1], 4) << "\n";

	} else {
		std::cout << "set size ratio -1\nset key off\n";
		const auto [xmin, xmax] = std::minmax_element(euler[0].begin(), euler[0].end());
		const auto [ymin, ymax] = std::minmax_element(euler[1].begin(), euler[1].end());
		double dx = (*xmax - *xmin) * 0.1, dy = (*ymax - *ymin) * 0.5;
		std::cout << "set xrange [" << (*xmin-dx) << ':' << (*xmax+dx) << "]\nset yrange [" << (*ymin-dy) << ':' << (*ymax+dy) << "]\n";
		
		std::cout << "$euler << EOD\n\n";
		for (size_t i = 0; i < euler[0].size(); ++i) {
			std::cout << euler[0][i] << ' ' << euler[1][i] << '\n';
		}
		std::cout << "EOD\n\n";

		std::cout << "$runge << EOD\n\n";
		for (size_t i = 0; i < runge[0].size(); ++i) {
			std::cout << runge[0][i] << ' ' << runge[1][i] << '\n';
		}
		std::cout << "EOD\n\n";

		std::cout << "$adams << EOD\n\n";
		for (size_t i = 0; i < adams[0].size(); ++i) {
			std::cout << adams[0][i] << ' ' << adams[1][i] << '\n';
		}
		std::cout << "EOD\n\n";

		std::cout << "f(x) = x * sin(x) + cos(x)\n";

		std::cout << "\nset title 'Euler method'\n"
		          << "plot f(x) with lines ls 4 lw 1.5, $euler with " << ((euler[0].size() <= 11) ? "linespoints ls 22 pt 22 " : "lines ls 22 ") << "lw 1.5 dt 3\n";

		std::cout << "\n# set title 'Runge-Kutta method'\n"
		          << "# plot f(x) with lines ls 4 lw 1.5, $runge with " << ((runge[0].size() <= 11) ? "linespoints ls 22 pt 22 " : "lines ls 22 ") << "lw 1.5 dt 3\n";

		std::cout << "\n# set title 'Adams method'\n"
		          << "# plot f(x) with lines ls 4 lw 1.5, $adams with " << ((adams[0].size() <= 11) ? "linespoints ls 22 pt 22 " : "lines ls 22 ") << "lw 1.5 dt 3\n";
	}
}