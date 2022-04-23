#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>

#include "../include/function/spline.hpp"

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

	CubicSpline<double> spline(x, y);

	if (!plot) {
		std::cout << "Значение в точке " << x0 << ": " << spline(x0) << '\n';
		std::cout << "Сегменты кубического сплайна:\n" << spline;

	} else {
		std::cout << "set size ratio -1\n";
		const auto [xmin, xmax] = std::minmax_element(spline.x.begin(), spline.x.end());
		const auto [ymin, ymax] = std::minmax_element(spline.y.begin(), spline.y.end());
		double dx = (*xmax - *xmin) * 0.2, dy = (*ymax - *ymin) * 0.5;
		std::cout << "set xrange [" << (*xmin-dx) << ':' << (*xmax+dx) << "]\nset yrange [" << (*ymin-dy) << ':' << (*ymax+dy) << "]\n";
		std::cout << "plot sample ";
		for (size_t i = 0; i < spline.Size(); ++i) {
			std::cout << '[' << spline.x[i] << ':' << spline.x[i+1] << "] " << spline[i] << " title 'segment " << i << "' lw 1.5";
			if (i != spline.Size() - 1) {
				std::cout << "\nreplot ";
			}
		}
	}
}