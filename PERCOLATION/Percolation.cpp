#include "Percolation.h"

/// <summary>
/// If we ever enconquer negative index in mass function we need to recursively get into it to find the real cluster
/// </summary>
/// <param name="badk">negative</param>
/// <param name="massDist">current mass distribution in hoshen kopelman algorithm</param>
/// <returns></returns>
int percolation::SitePercolationMonteCarlo::recursiveClusterNumberDetection(int k, const std::vector<int>& massDist)
{
	if (massDist[k] > 0) return k;
	else return recursiveClusterNumberDetection(abs(massDist[k]) , massDist);
}

percolation::SitePercolationMonteCarlo::SitePercolationMonteCarlo()
{
	this->probability = 0;
	this->isPercolating = 0;
	this->maxClusterSize = 0;
	this->shortestPath = 0;

}

percolation::SitePercolationMonteCarlo::SitePercolationMonteCarlo(double p, general::lattice2D& lattice)
{
	this->probability = p;
	this->lattice = lattice.move_clone();
	this->isPercolating = 0;
	this->maxClusterSize = 0;
	this->shortestPath = 0;
	/* initialize */
	this->sites = std::vector<std::vector<int>>(this->lattice->get_Lx(), std::vector<int>(this->lattice->get_Ly()));
	this->setFields();
}

percolation::SitePercolationMonteCarlo::SitePercolationMonteCarlo(const SitePercolationMonteCarlo& A)
{
	this->probability = A.probability;
	this->lattice = A.lattice->move_clone();
	this->sites = A.sites;
	this->setFields();
	this->burned = A.burned;
	this->clusters = A.burned;
	this->maxClusterSize = A.maxClusterSize;
	this->shortestPath = A.shortestPath;

}

percolation::SitePercolationMonteCarlo::SitePercolationMonteCarlo(SitePercolationMonteCarlo&& A) noexcept
{
	this->probability = A.probability;
	this->lattice = A.lattice->move_clone();
	this->sites = A.sites;
	this->setFields();
	this->burned = A.burned;
	this->clusters = A.burned;
	this->maxClusterSize = A.maxClusterSize;
	this->shortestPath = A.shortestPath;
}

void percolation::SitePercolationMonteCarlo::setFields()
{
	double p = this->probability;
	this->isPercolating = 0;
	for_each(this->sites.begin(), this->sites.end(), [p](std::vector<int>& a) mutable
		{
			for_each(a.begin(), a.end(), [a,p](int& b) mutable {b = (randZero_One_Uni()) < p ? 1 : 0; });
		}
	);
	this->burned = sites;
	this->clusters = sites;
	this->sizeDistribution = std::vector<int>(this->lattice->get_Ns(), 0); // biggest cluster possible shall be Ns
	this->massArray = std::vector<int>(2, 0);
	/* calculate again */
	this->burningMethod();
	this->hoshen_kopelman();
}

void percolation::SitePercolationMonteCarlo::burningMethod()
{
	/* --------------> X
	|
	|
	|
	v
	y
	*/

	int t = 2;
	// set all occupied top to 2
	std::vector<std::tuple<int,int>> currentVector;
	for (int i = 0; i < this->sites[0].size(); i++) 
	{
		if (this->burned[0][i] == 1)
		{
			currentVector.push_back(std::make_tuple(i, 0));
			burned[0][i] = t;
		}
	}
	bool haveChangedSomething = true;
	while (haveChangedSomething) {
		std::vector<std::tuple<int, int>> tempVector;
		haveChangedSomething = false;
		for (auto tuple : currentVector) {
			int n_num = this->lattice->get_nn_number(std::get<0>(tuple), std::get<1>(tuple));
			// iterate all nearest neighbors 
			for (int i = 0; i < n_num; i++) {
				int xneighbor = std::get<0>(this->lattice->get_nn(std::get<0>(tuple), std::get<1>(tuple), i));
				int yneighbor = std::get<1>(this->lattice->get_nn(std::get<0>(tuple), std::get<1>(tuple), i));
				// neighbors are outside allowed space -> no boundary conditions
				if (yneighbor == -1 && std::get<1>(tuple) != 0) {
					// got till the end!
					this->isPercolating = 1;
					this->shortestPath = t;
					break;
				}
				else if (yneighbor == -1 || xneighbor == -1) continue;
				else if(this->burned[yneighbor][xneighbor] == 1) {
					this->burned[yneighbor][xneighbor] = t+1;
					haveChangedSomething = true;
					// add to next iteration
					tempVector.push_back(std::make_tuple(xneighbor, yneighbor));
				}
				else {
					break;
				}
			}
		}
		currentVector = tempVector;
		t++;
	}

	this->shortestPath = t - 1;
	return;
}

void percolation::SitePercolationMonteCarlo::hoshen_kopelman()
{ 
	/* 
		LEFT AND TOP ARE FIRST TWO
	*/
	int k = 2; 

	/* SWEEP ALL SITES */
	for (int y = 0; y < this->lattice->get_Ly(); y++) {
		for (int x = 0; x < this->lattice->get_Lx(); x++) {
			// start building cluster
			if (this->clusters[y][x] == 1) {
				int xtop = std::get<0>(this->lattice->get_nn(x, y, 1));
				int ytop = std::get<1>(this->lattice->get_nn(x, y, 1));
				int xleft = std::get<0>(this->lattice->get_nn(x, y, 0));
				int yleft = std::get<1>(this->lattice->get_nn(x, y, 0));

				// check top
				int ktop = (xtop == -1 || ytop == -1) ? 0 : this->clusters[ytop][xtop];
				int kleft = (xleft == -1 || yleft == -1) ? 0 : this->clusters[yleft][xleft];

				// both empty
				if (kleft + ktop == 0) {
					// found new cluster
					this->massArray.push_back(1);
					// set new element to k
					this->clusters[y][x] = k;
					k++;
				}
				// one of them belongs to a cluster already
				else if (kleft * ktop == 0) {
					int which_k = kleft > 0 ? this->recursiveClusterNumberDetection(kleft, massArray) : this->recursiveClusterNumberDetection(ktop, massArray);
					// we mark by negative value of the clusters combined
					// so we need to make sure that this isn't the one combining to another cluster
					// else simply addd 1
					massArray[which_k]++;
					// set new element to that cluster
					this->clusters[y][x] = which_k;
				}
				else{
					// if they are the same do only one thing
					if (kleft == ktop) {
						int which = this->recursiveClusterNumberDetection(kleft, massArray);
						massArray[which]++;
						this->clusters[y][x] = which;
					}
					// if they are different but should belong to the same cluster we take the left one!
					else if (kleft < ktop) {
						// change the left k mass function to M1 + M2 + 1 but remember both can be negative => search recursively again
						int ksmaller = this->recursiveClusterNumberDetection(kleft, massArray);
						int kbigger = this->recursiveClusterNumberDetection(ktop, massArray);
						massArray[ksmaller] += massArray[kbigger] + 1;
						// set this array to -k_smaller
						massArray[kbigger] = -ksmaller;
						// change the smaller k mass function to M1 + M2 + 1 but remember both can be negative => search recursively again
						this->clusters[y][x] = ksmaller;
					}
					else {
						// change the smaller k mass function to M1 + M2 + 1 but remember both can be negative => search recursively again
						int ksmaller = this->recursiveClusterNumberDetection(ktop, massArray);
						int kbigger = this->recursiveClusterNumberDetection(kleft, massArray);
						massArray[ksmaller] += massArray[kbigger] + 1;
						// set this array to -k_smaller
						massArray[kbigger] = -ksmaller;
						// change the smaller k mass function to M1 + M2 + 1 but remember both can be negative => search recursively again
						this->clusters[y][x] = ksmaller;
					}
				}
			}
			continue;
		}
	}
	/* CREATE DISTRIBUTION */
	for (int elem : massArray) {
		if (this->maxClusterSize < elem) maxClusterSize = elem;
		if(elem > 0) this->sizeDistribution[elem]++;
	}
}

int percolation::SitePercolationMonteCarlo::getIfPercolating() const
{
	return double(this->isPercolating);
}

std::vector<std::vector<int>> percolation::SitePercolationMonteCarlo::getBurned() const
{
	return this->burned;
}

std::vector<std::vector<int>> percolation::SitePercolationMonteCarlo::getClusters()
{
	for (int y = 0; y < this->lattice->get_Ly(); y++) {
		for (int x = 0; x < this->lattice->get_Lx(); x++) {
			// change again to see how does it look;
			int trueClusterNumber = this->clusters[y][x] > 0 ? this->recursiveClusterNumberDetection(this->clusters[y][x], this->massArray) : 0;
			this->clusters[y][x] =  trueClusterNumber;
			std::cout << this->clusters[y][x] << "\t";
		}
		std::cout << std::endl;
	}
	return this->clusters;
}

void percolation::SitePercolationMonteCarlo::setProbability(double p)
{
	this->probability = p;
	this->isPercolating = 0;
	/* calculate again */
	this->setFields();
}





/* SIMULATION FUNCTIONS */

double percolation::makeSingleSimulation(double rho, int Lx, int Ly, int mcs) {
	double prob = 0;
	general::square_lattice lattice(Lx, Ly);
	percolation::SitePercolationMonteCarlo sim(rho, lattice);
	for (int i = 0; i < mcs; i++) {
		sim.setFields();
		//cout << sim.getIfPercolating() << endl;
		prob += sim.getIfPercolating();
	}
	return prob;
}