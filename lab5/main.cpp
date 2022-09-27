#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstring>

#include "../include/partial_differential/parabolic_pde.hpp"

const double PI = 2 * std::acos(0.0);

int main(int argc, char* argv[]) {
	bool plot = false;
	if (argc == 2 && std::strcmp(argv[1], "plot") == 0) {
		plot = true;
	}

	double a = 1, b = 0, c = 0;
	auto f = [](double, double){ return 0; };
	double start = 0, end = 1, t_end;
	auto psi = [](double x){ return std::sin(PI * x); };
	auto phi_start = [](double){ return 0; };
	auto phi_end = [](double){ return 0; };

	int n;
	double sigma;
	std::cin >> t_end >> n >> sigma;

	std::vector<double> t, x;
	grid_t<double> u;
	std::cout << std::setprecision(12) << std::fixed;
	std::tie(t, x, u) = ParabolicPDESolver<double>(a, b, c, f, start, end, t_end, psi, phi_start, phi_end, n, sigma);

	if (!plot) {
		for (size_t k = 0; k < t.size(); ++k) {
			for (size_t i = 0; i < x.size(); ++i) {
				std::cout << u[k][i] << ' ';
			}
			std::cout << '\n';
		}
	} else {
		std::cout << "set size ratio -1\nset key off\nset xlabel \"t\"\nset ylabel \"x\"\nset zlabel \"u(x, t)\"\n$data << EOD\n";
		for (size_t k = 0; k < t.size(); ++k) {
			for (size_t i = 0; i < x.size(); ++i) {
				std::cout << t[k] << ' ' << x[i] << ' ' << u[k][i] << '\n';
			}
			std::cout << '\n';
		}
		std::cout << "EOD\nsplot ";
		for (size_t k = 0; k < t.size()-1; ++k) {
			std::cout << "$data every :::" << k << "::" << k << " with lines, \\\n";
		}
		std::cout << "$data every :::" << t.size()-1 << "::" << t.size()-1 << " with lines";
	}
}
