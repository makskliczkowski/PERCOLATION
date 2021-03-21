#include "General.h"
#include "RandomAlghorithms.h"
#include "Percolation.h"


using namespace std;


int main(int argc, char* argv[])
{

	int processors = omp_get_num_procs();
	int mcsteps_per_proc = 5000;
	int Lx = 10;
	int Ly = 10;
	double rho = 0.4;
	//mkl_disable_fast_mm();
	//omp_set_max_active_levels(1);
	//std::cout << "Number of processors aviable = " << omp_get_num_procs() << "\n\n";
	/* PARAMETERS */
	/*
#pragma omp parallel for
	for (int pnum = 0; pnum < 16; pnum++) {
		double rho = 0.52 + pnum * 0.005;
		double prob = percolation::makeSingleSimulation(rho, Lx, Ly, mcsteps_per_proc);
#pragma omp critical
		cout << "rho = " << rho << ", probability=" <<  prob / mcsteps_per_proc << endl;
	}*/
	general::square_lattice lattice(Lx, Ly);
	percolation::SitePercolationMonteCarlo sim(rho, lattice);
	auto a = sim.getClusters();
	return 0;
}

