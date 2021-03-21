#pragma once
#ifndef PERCOLATION_H
#define PERCOLATION_H
#include "General.h"
#include "RandomAlghorithms.h"

namespace percolation{

	class SitePercolationMonteCarlo {
	private:
		double probability;
		double isPercolating;
		int shortestPath;
		int maxClusterSize;
		// lattice for the model
		std::unique_ptr<general::lattice2D> lattice;
		// all the sites
		std::vector<std::vector<int>> sites;
		std::vector<std::vector<int>> burned;
		std::vector<std::vector<int>> clusters;
		// size distribution for the cluster
		std::vector<int> massArray;
		std::vector<int> sizeDistribution;
	protected:
		int recursiveClusterNumberDetection(int k, const std::vector<int>& massDist);
	public:
		/* General */
		SitePercolationMonteCarlo();
		SitePercolationMonteCarlo(double p, general::lattice2D& lattice);
		~SitePercolationMonteCarlo() = default;
		SitePercolationMonteCarlo(const SitePercolationMonteCarlo& A);
		SitePercolationMonteCarlo(SitePercolationMonteCarlo&& A) noexcept;

		/* SIMULATION */
		void setFields();
		void burningMethod();
		void hoshen_kopelman();

		/* GETTERS */
		int getIfPercolating() const;
		std::vector<std::vector<int>> getBurned()const;
		std::vector<std::vector<int>> getClusters();

		/* PRINTERS */
		void printBurned(std::ostream& output, std::string separator = "\t");
		void printClusters(std::ostream& output, std::string separator = "\t");

		/* SETTERS */
		void setProbability(double p);




	};
	/* SIMULATION FUNCTIONS */

	double makeSingleSimulation(double rho, int Lx, int Ly, int mcs);


}



#endif