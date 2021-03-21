#include "RandomAlghorithms.h"

// Generate global variable
std::random_device rd;
std::mt19937 gen(rd());


/// <summary>
/// Makes Euclidean modulo
/// </summary>
/// <param name="a">first value</param>
/// <param name="b">second value</param>
/// <returns>a%b</returns>
int myModuloEuclidean(int a, int b)
{
	int m = a % b;
	if (m < 0) {
		// m += (b < 0) ? -b : b; // avoid this form: it is UB when b == INT_MIN
		m = (b < 0) ? m - b : m + b;
	}
	return m;
}

/// <summary>
/// Produces value from 0.0 to 1.0 in uniform distribution
/// </summary>
/// <returns>value from 0.0 to 1.0 in uniform distribution</returns>
double randZero_One_Uni()
{
	std::uniform_real_distribution<> dist(0.0, 1.0);
	return dist(gen);
}

int rand0_N_Uni(int L) {
	std::uniform_int_distribution<> dist(0, L);
	return dist(gen);
}
