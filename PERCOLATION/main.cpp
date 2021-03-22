#include "General.h"
#include "RandomAlghorithms.h"
#include "Percolation.h"


using namespace std;


int main(int argc, char* argv[])
{

	int processors = omp_get_num_procs();
	percolation::makeSingleSimulation();
	return 0;
}

