#pragma once
#include <random>


constexpr double PI = 3.14159265359;
constexpr double PI_half = PI / 2.0;
constexpr double TWO_PI = PI * 2;

// Euclidean modulo function denoting also the negative sign
int myModuloEuclidean(int a, int b);
// Function that gives number from 0 to 1 from Uniform Distribution
double randZero_One_Uni();
// Function that gives number from 0 to N from Uniform Distribution
int rand0_N_Uni(int L);


