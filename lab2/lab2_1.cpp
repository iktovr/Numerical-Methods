#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <vector>
#include <utility>

#include "../include/nonlinear/solvers.hpp"

int main() {
	std::cout.precision(8);
	std::cout << std::fixed;

	std::function<double(double)> f = [](double x){ return std::log(x + 2) - x * x; };
	std::function<double(double)> df = [](double x){ return 1.0 / (x + 2) - 2 * x; };
	std::vector<interval_t<double>> intervals = {{-1.9, 0}, {0.5, 2}};

	double eps;
	std::cin >> eps;

	double x0;

	std::cout << "Метод простой итерации:\n";
	for (auto interval: intervals) {
		x0 = (interval.first + interval.second) / 2.0;
		auto [x, iter_count] = IterationMethod(x0, interval, f, df, eps);

		std::cout << "Корень: " << x << "\nКоличество итераций: " << iter_count << '\n';
	}

	std::cout << '\n';
	std::cout << "Метод Ньютона:\n";
	for (auto [a, b]: intervals) {
		x0 = (a + b) / 2.0;
		auto [x, iter_count] = NewtonMethod(x0, f, df, eps);

		std::cout << "Корень: " << x << "\nКоличество итераций: " << iter_count << '\n';
	}
}