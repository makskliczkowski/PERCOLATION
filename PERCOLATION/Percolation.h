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
		general::lattice2D* lattice;
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
		SitePercolationMonteCarlo(double p, general::lattice2D* lattice);
		~SitePercolationMonteCarlo() = default;
		SitePercolationMonteCarlo(const SitePercolationMonteCarlo& A);
		SitePercolationMonteCarlo(SitePercolationMonteCarlo&& A) noexcept;

		/* SIMULATION */
		void setFields();
		void burningMethod();
		void hoshen_kopelman();

		/* GETTERS */
		int getIfPercolating() const;
		int getSmax() const;
		std::vector<std::vector<int>> getBurned()const;
		std::vector<std::vector<int>> getClusters();
		std::vector<int> getDistribution() const;
		int getDistributionElem(int elem);

		/* PRINTERS */
		void printBurned(std::string directory, std::string separator, int numberOfPicture = 0);
		void printClusters(std::string directory, std::string separator, int numberOfPicture = 0);

		/* SETTERS */
		void setProbability(double p);




	};
	/* SIMULATION FUNCTIONS */
	
	enum parsers {
		T,
		L,
		l,
		p0,
		pk,
		dp,
		d
	};



	void exit_with_help();


	void makeSingleSimulation(std::string param_file = "input.txt");
	void singleFor(int mcsteps, int L, double rho,general::lattice2D* lattice, std::string lat_type_str, std::string sv_dir);
	// print flow probability and maximal cluster size
	void printProbability(std::string directory, int L, int mcsteps, double rho, double av_pflow, double av_smax, std::string, std::string separator = "  ");
	// print distribution of the clusters
	void printDistribution(std::string directory, int L, int mcsteps, double rho, std::string lat_type, const std::vector<double>& ns, std::string separator = "  ");
}


#endif