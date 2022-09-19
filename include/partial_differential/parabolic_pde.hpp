#pragma once

#include <functional>
#include <vector>

template <class T>
using grid_t = std::vector<std::vector<T>>;

template <class T>
grid_t<T> ParabolicPDESolver(T a, T b, T c, const std::function<T(T, T)>& f, T start, T end, T t_end, const std::function<T(T)>& psi, 
                             const std::function<T(T)>& phi_start, const std::function<T(T)>& phi_end, int h_count, int tau_count) {
	// u_t = a * u_xx + b * u_x + c * u + f
	// x in [start, end]
	// t in [0, t_end]
	// u(x, 0) = psi(x)
	// u(start, t) = phi_start(t)
	// u(end, t) = phi_end(t)
	// h = (end - start) / h_count
	// tau = t_end / tau_count

	grid_t<T> grid(tau_count + 1, std::vector<T>(h_count + 1));
	std::vector<T> x(h_count + 1), t(tau_count + 1);
	T h = (end - start) / h_count;
	T tau = t_end / tau_count;

	for (int i = 0; i <= h_count; ++i) {
		x[i] = start + h * i;
	}
	for (int k = 0; k <= tau_count; ++k) {
		t[k] = tau * k;
	}

	for (int i = 0; i <= h_count; ++i) {
		grid[0][i] = psi(x[i]);
	}
	for (int k = 1; k <= tau_count; ++k) {
		grid[k][0] = phi_start(t[k]);
		grid[k][h_count] = phi_end(t[k]);
	}

	for (int k = 0; k < tau_count; ++k) {
		for (int i = 1; i < h_count; ++i) {
			grid[k+1][i] = ((grid[k][i-1] - 2 * grid[k][i] + grid[k][i+1]) * a / (h * h) + (grid[k][i+1] - grid[k][i-1]) * b / (2 * h) +
			               c * grid[k][i] + f(x[i], t[k])) * tau + grid[k][i];
		}
	}
	return grid;
}