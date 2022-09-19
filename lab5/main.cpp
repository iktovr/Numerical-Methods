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

	int h_count, tau_count;
	std::cin >> t_end >> h_count >> tau_count;

	grid_t<double> u = ParabolicPDESolver<double>(a, b, c, f, start, end, t_end, psi, phi_start, phi_end, h_count, tau_count);

	if (!plot) {
		std::cout << std::setprecision(2) << std::fixed;
		for (int k = 0; k <= tau_count; ++k) {
			for (int i = 0; i <= h_count; ++i) {
				std::cout << u[k][i] << ' ';
			}
			std::cout << '\n';
		}
	} else {
		double h = (end - start) / h_count;
		double tau = t_end / tau_count;

		std::cout << "set size ratio -1\nset key off\nset xlabel \"t\"\nset ylabel \"x\"\nset zlabel \"u(x, t)\"\n$data << EOD\n";
		for (int k = 0; k <= tau_count; ++k) {
			for (int i = 0; i <= h_count; ++i) {
				std::cout << tau * k << ' ' << start + h * i << ' ' << u[k][i] << '\n';
			}
			std::cout << '\n';
		}
		std::cout << "EOD\nsplot ";
		for (int k = 0; k < tau_count; ++k) {
			std::cout << "$data every :::" << k << "::" << k << " with lines, \\\n";
		}
		std::cout << "$data every :::" << tau_count << "::" << tau_count << " with lines";
	}
}