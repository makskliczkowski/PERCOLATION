#include "Percolation.h"

/// <summary>
/// If we ever enconquer negative index in mass function we need to recursively get into it to find the real cluster
/// </summary>
/// <param name="badk">negative</param>
/// <param name="massDist">current mass distribution in hoshen kopelman algorithm</param>
/// <returns></returns>
int percolation::SitePercolationMonteCarlo::recursiveClusterNumberDetection(int k, const std::vector<int>& massDist)
{
	if (massDist[k] >= 0) return k;
	else return recursiveClusterNumberDetection(abs(massDist[k]) , massDist);
}
/// <summary>
/// default constructor
/// </summary>
percolation::SitePercolationMonteCarlo::SitePercolationMonteCarlo()
{
	this->probability = 0;
	this->isPercolating = 0;
	this->maxClusterSize = 0;
	this->shortestPath = 0;

}
/// <summary>
/// used constructor
/// </summary>
/// <param name="p"></param>
/// <param name="lattice"></param>
percolation::SitePercolationMonteCarlo::SitePercolationMonteCarlo(double p,  general::lattice2D* lattice)
{
	this->probability = p;
	this->lattice = lattice;
	this->isPercolating = 0;
	this->maxClusterSize = 0;
	this->shortestPath = 0;
	/* initialize */
	this->sites = std::vector<std::vector<int>>(this->lattice->get_Ly(), std::vector<int>(this->lattice->get_Lx()));
	// note that we got reversed order of lx and ly here because of neighbors convinience!
	this->setFields();
}
/// <summary>
/// copy constructor
/// </summary>
/// <param name="A"></param>
percolation::SitePercolationMonteCarlo::SitePercolationMonteCarlo(const SitePercolationMonteCarlo& A)
{
	this->probability = A.probability;
	this->lattice = A.lattice;
	this->sites = A.sites;
	this->setFields();
	this->burned = A.burned;
	this->clusters = A.burned;
	this->maxClusterSize = A.maxClusterSize;
	this->shortestPath = A.shortestPath;

}
/// <summary>
/// move constructor
/// </summary>
/// <param name="A"></param>
/// <returns></returns>
percolation::SitePercolationMonteCarlo::SitePercolationMonteCarlo(SitePercolationMonteCarlo&& A) noexcept
{
	this->probability = A.probability;
	this->lattice = A.lattice;
	this->sites = A.sites;
	this->setFields();
	this->burned = A.burned;
	this->clusters = A.burned;
	this->maxClusterSize = A.maxClusterSize;
	this->shortestPath = A.shortestPath;
}
/// <summary>
/// Setting fields - > resets actually
/// </summary>
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
/// <summary>
/// Burning method for finding minimal cluster span
/// </summary>
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
					return;
				}
				else if (yneighbor == -1 || xneighbor == -1) continue;
				else if(this->burned[yneighbor][xneighbor] == 1) {
					this->burned[yneighbor][xneighbor] = t+1;
					haveChangedSomething = true;
					// add to next iteration
					tempVector.push_back(std::make_tuple(xneighbor, yneighbor));
				}
				else {
					continue;
				}
			}
		}
		currentVector = tempVector;
		t++;
	}

	this->shortestPath = t - 1;
	return;
}
/// <summary>
/// Hoshen-Kopelman algorithm for clusters recognition
/// </summary>
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
				int ktop = (xtop == -1 || ytop == -1) ? 0 : this->recursiveClusterNumberDetection(this->clusters[ytop][xtop],this->massArray);
				int kleft = (xleft == -1 || yleft == -1) ? 0 : this->recursiveClusterNumberDetection(this->clusters[yleft][xleft],this->massArray);

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
					int which_k = kleft > 0 ? kleft:ktop;
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
						int which = kleft;
						massArray[which]++;
						this->clusters[y][x] = which;
					}
					// if they are different but should belong to the same cluster we take the left one!
					else if (kleft < ktop) {
						// change the left k mass function to M1 + M2 + 1 but remember both can be negative => search recursively again
						int ksmaller = kleft;
						int kbigger = ktop;
						massArray[ksmaller] += massArray[kbigger] + 1;
						// set this array to -k_smaller
						massArray[kbigger] = -ksmaller;
						// change the smaller k mass function to M1 + M2 + 1 but remember both can be negative => search recursively again
						this->clusters[y][x] = ksmaller;
					}
					else {
						// change the smaller k mass function to M1 + M2 + 1 but remember both can be negative => search recursively again
						int ksmaller = ktop;
						int kbigger = kleft;
						massArray[ksmaller] += massArray[kbigger] + 1;
						// set this array to -k_smaller
						massArray[kbigger] = -ksmaller;
						// change the smaller k mass function to M1 + M2 + 1 but remember both can be negative => search recursively again
						this->clusters[y][x] = ksmaller;
					}
				}
			}
			else {
				continue;
			}
		}
	}
	/* CREATE DISTRIBUTION */
	for (int elem : massArray) {
		if (this->maxClusterSize < elem) maxClusterSize = elem;
		if(elem > 0) this->sizeDistribution[elem]++;
	}
}
/// <summary>
/// Tells us if we have percolating cluster or not
/// </summary>
/// <returns>True or false</returns>
int percolation::SitePercolationMonteCarlo::getIfPercolating() const
{
	return double(this->isPercolating);
}
/// <summary>
/// Maximum cluster size
/// </summary>
/// <returns></returns>
int percolation::SitePercolationMonteCarlo::getSmax() const
{
	return this->maxClusterSize;
}
/// <summary>
/// returning burned
/// </summary>
/// <returns>burned matrix</returns>
std::vector<std::vector<int>> percolation::SitePercolationMonteCarlo::getBurned() const
{
	return this->burned;
}
/// <summary>
/// returning clusters
/// </summary>
/// <returns>clusters matrix</returns>
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
/// <summary>
/// get distribution matrix
/// </summary>
/// <returns></returns>
std::vector<int> percolation::SitePercolationMonteCarlo::getDistribution() const
{
	return this->sizeDistribution;
}
/// <summary>
/// get element from distribution matrix
/// </summary>
/// <param name="elem">elem number</param>
/// <returns>element from distribution matrix</returns>
int percolation::SitePercolationMonteCarlo::getDistributionElem(int elem)
{
	return this->sizeDistribution[elem];
}
/// <summary>
/// Printing the picture of burned elements
/// </summary>
/// <param name="directory">directory to be saved into</param>
/// <param name="separator">separator between values</param>
/// <param name="numberOfPicture">picture number if we want more of them</param>
void percolation::SitePercolationMonteCarlo::printBurned(std::string directory, std::string separator, int numberOfPicture)
{
	using namespace std;
	std::string lat = this->lattice->getString_type() + std::string(kPathSeparator);
	/* create folder if doesn't exist */
	string folder = directory + lat + "BurnedPictures" +  string(kPathSeparator) + "L" + std::to_string(this->lattice->get_Lx()) + std::string(kPathSeparator);
	std::filesystem::create_directories(folder);

	/* open file */
	string filename = folder + "Burned_L" + to_string(this->lattice->get_Lx()) + "p" + to_string(this->probability) + "_" + to_string(numberOfPicture) + ".txt";
	std::ofstream file;
	file.open(filename, ios::out);
	if (!file.is_open()) {
		cout << "Couldn't open a file : " + filename << endl;
		exit(-1);
	}
	cout << "Saving data to : " + filename << endl;
	// printing matrix to file
	for (auto row : this->burned) {
		for (auto elem : row) {
			file << elem << separator;
		}
		file << endl;
	}
	file.close();
}
/// <summary>
/// Printing clusters as a matrix
/// </summary>
/// <param name="directory">directory to be saved into</param>
/// <param name="separator">separator between values</param>
/// <param name="numberOfPicture">picture number if we want more of them</param>
void percolation::SitePercolationMonteCarlo::printClusters(std::string directory, std::string separator, int numberOfPicture)
{
	using namespace std;
	std::string lat = this->lattice->getString_type() + std::string(kPathSeparator);
	/* create folder if doesn't exist */
	string folder = directory + lat + "ClusterPictures" + string(kPathSeparator) + "L" + std::to_string(this->lattice->get_Lx()) + std::string(kPathSeparator);
	std::filesystem::create_directories(folder);

	/* open file */
	string filename = folder + "Cluster_L" + to_string(this->lattice->get_Lx()) + "p" + general::to_string_prec(this->probability) + "_" + to_string(numberOfPicture) + ".txt";
	std::ofstream file;
	file.open(filename, ios::out);
	if (!file.is_open()) {
		cout << "Couldn't open a file : " + filename << endl;
		exit(-1);
	}
	cout << "Saving data to : " + filename << endl;
	// printing matrix to file
	for (auto row : this->clusters) {
		for (auto elem : row) {
			file << this->recursiveClusterNumberDetection(elem,this->massArray) << separator;
		}
		file << endl;
	}
	file.close();
}
/// <summary>
/// Setting new probability -> restarts system
/// </summary>
/// <param name="p">occupation probability</param>
void percolation::SitePercolationMonteCarlo::setProbability(double p)
{
	this->probability = p;
	this->isPercolating = 0;
	/* calculate again */
	this->setFields();
}





/* SIMULATION FUNCTIONS */

/// <summary>
/// To give message
/// </summary>
void percolation::exit_with_help()
{
	printf(
		" LOOKING FOR FILLE IN CURRENT DIRECTORY -> input.txt\n"
		"Usage: value % [options] \n"
		"options:\n"
		" T monte carlo steps : bigger than 0 (default 1000)\n"
		" L lattice size >0 (default 100)\n"
		" l lattice type : (default square):\n"
		"   square, triangle, hexagonal \n"
		" p0 starting occupation probability : p0>=0 (default 0)\n"
		" pk ending occupation probability : 0<pk<=1 (default 1)\n"
		" dp probability step : (default 0.1)\n"
		" d directory to save files (default current directory)\n"
	);
	exit(1);
}

/// <summary>
/// Creates a single simulation
/// </summary>
/// <param name="rho">occupation probability</param>
/// <param name="Lx">x size</param>
/// <param name="Ly">y size</param>
/// <param name="mcs">Monte Carlo steps</param>
/// <returns></returns>
void percolation::makeSingleSimulation(std::string param_file) {
	std::unordered_map <std::string, percolation::parsers> const parse_table = {
	{"T",T},
	{"d",d},
	{"l",l},
	{"L",L},
	{"p0",p0},
	{"pk",pk},
	{"dp",dp},
	};
	/* SETTING DEFAULT */
	general::lattice2D* lattice;
	int mcsteps = 1000;
	int Lx = 100;
	double p_0 = 0;
	double p_k = 1;
	double d_p = 0.1;
	general::lattice2D::lattice_types lat_type = general::lattice2D::lattice_types::square;
	std::string lat_type_str = "square";
	std::string sv_dir = "." + std::string(kPathSeparator);
	/* PARSE FILE */
	std::fstream file;
	file.open(param_file, std::ios::in);
	if (!file.is_open()) {
		std::cout << "Couldn't open input file ;c\n";
		percolation::exit_with_help();
		exit(-1);
	}
	std::string line;
	bool changedSmthing = false;
	while (std::getline(file, line)) {
		std::vector<std::string> split = general::split_str(line, "%");
		if (split.size() >= 2) {
			// our argument to switch
			percolation::parsers enum_arg;
			// find command in parser
			auto it = parse_table.find(split[1]);
			if (it != parse_table.end()) {
				enum_arg = it->second;
				changedSmthing = true;
			}
			// if not found tell us!
			else {
				std::cout << (split[1]) <<  " NOT FOUND " << std::endl;
				//cout << "Setting default training model parameters\n";
				//exit_with_help();
				continue;
			}
			switch (enum_arg) {
			default:
				std::cout << "unrecognized argument\n";
				exit_with_help();
				break;
			case T:
				mcsteps = std::stoi(split[0]);
				if (mcsteps < 0) {
					std::cout << "Mcsteps can't be negative! Setting default!\n";
					mcsteps = 1000;
				}
				else changedSmthing = true;
				break;
			case L:
				Lx = std::stoi(split[0]);
				if (L < 0) {
					std::cout << "L can't be negative! Setting default!\n";
					Lx = 100;
				}
				else changedSmthing = true;
				break;
			case p0:
				p_0 = std::stod(split[0]);
				if (p_0 < 0 || p_0 >= 1) {
					std::cout << "p0 can't be like that! Setting default!\n";
					p_0 = 0.0;
				}
				else changedSmthing = true;
				break;
			case pk:
				p_k = std::stod(split[0]);
				if (p_k < 0 || p_k > 1 || p_k < p_0) {
					std::cout << "pk can't be like that! Setting default!\n";
					p_k = 0.0;
				}
				else changedSmthing = true;
				break;
			case dp:
				d_p = std::stod(split[0]);
				if (d_p < 0 || d_p > 1) {
					std::cout << "dp can't be like that! Setting default!\n";
					p_k = 0.0;
				}
				else changedSmthing = true;
				break;
			case l:
				// our argument to switch
				general::lattice2D::lattice_types enum_arg;
				// find lattice type in parser
				auto it = general::lattice_parse_table.find(split[0]);
				if (it != general::lattice_parse_table.end()) {
					enum_arg = it->second;
					changedSmthing = true;
					lat_type = enum_arg;
					lat_type_str = split[0];
				}
				else {
					std::cout << "I don't have that lattice : " + split[0] + "! Setting default!\n";
				}
				break;
			}
		}
	}

	// get lattice type
	switch (lat_type) {
	case general::lattice2D::lattice_types::square:
		lattice = new general::square_lattice(Lx, Lx);
		lat_type_str = "square";
		break;
	case general::lattice2D::lattice_types::triangle:
		lattice = new general::triangle_lattice(Lx, Lx);
		lat_type_str = "triangle";
		break;
	default:
		std::cout << "Sorry, must add them\n";
		lattice = new general::square_lattice(Lx, Lx);
	}
	/* SIMULATION PART */

	int pnum = int( (p_k - p_0) / d_p);

	// FOR OTHER VALUES
#pragma omp parallel for
	for (int i = 0; i <= pnum; i++) {
		double rho = p_0 + i * d_p;
		percolation::singleFor(mcsteps, Lx, rho, lattice, lat_type_str, sv_dir);
	}
	file.close();

}
void percolation::singleFor(int mcsteps, int L, double rho, general::lattice2D* lattice,std::string lat_type_str, std::string sv_dir)
{
	double tol = 1.0e-6;
	double pflow = 0.0;
	double smax = 0.0;
	std::vector<double> avDist(lattice->get_Ns(), 0.0);
	/* make sim */
	percolation::SitePercolationMonteCarlo sim(rho, lattice);
	for (int step = 0; step < mcsteps; step++) {
		pflow += double(sim.getIfPercolating());
		smax += double(sim.getSmax());
		for (int elem = 0; elem < lattice->get_Ns(); elem++) {
			avDist[elem] += sim.getDistributionElem(elem);
		}
		sim.setFields();
	}
	/* normalise */
	pflow = 1.0 * pflow / (1.0 * mcsteps);
	smax = 1.0 * smax / (1.0 * mcsteps * lattice->get_Ns()); // relative size
	for (int elem = 0; elem < lattice->get_Ns(); elem++) {
		avDist[elem] = avDist[elem] / (1.0 * mcsteps);
	}
	/* print */

#pragma omp critical
	percolation::printProbability(sv_dir, L, mcsteps, rho, pflow, smax, lat_type_str, "  ");
#pragma omp critical
	if ((rho > 0.2-tol && rho < 0.2 + tol) || (rho > 0.3 - tol && rho < 0.3 + tol) || (rho > 0.4 - tol && rho < 0.4 + tol) || (rho > 0.5 - tol && rho < 0.5 + tol) || (rho > 0.6 - tol && rho < 0.6 + tol) || (rho > 0.7 - tol && rho < 0.7 + tol) || (rho > 0.8 - tol && rho < 0.8 + tol) || rho == 0.592746) {
		sim.printBurned(sv_dir, "  ");
		sim.printClusters(sv_dir, "  ");
		percolation::printDistribution(sv_dir, L, mcsteps, rho, lat_type_str, avDist, "  ");
	}
}
void percolation::printProbability(std::string directory,int L, int mcsteps,double rho, double av_pflow, double av_smax, std::string lat_type, std::string separator)
{
	using namespace std;
	std::string lat = lat_type + std::string(kPathSeparator);
	/* create folder if doesn't exist */
	string folder = directory + lat ;
	std::filesystem::create_directories(folder);

	/* open file */
	string filename = folder + "Ave_L" + to_string(L) + "T" + to_string(mcsteps) + ".txt";
	std::ofstream file;
	file.open(filename, ios::app | ios::out);
	if (!file.is_open()) {
		cout << "Couldn't open a file : " + filename << endl;
		exit(-1);
	}
	cout << "Saving data to : " + filename << endl;
	file << to_string(rho) << separator << to_string(double(av_pflow)) << separator << to_string(av_smax) << endl;
	file.close();
}

void percolation::printDistribution(std::string directory, int L, int mcsteps, double rho, std::string lat_type, const std::vector<double>& ns, std::string separator)
{
	using namespace std;
	std::string lat = lat_type + std::string(kPathSeparator);

	/* create folder if doesn't exist */
	string folder = directory + lat + "Distributions" + string(kPathSeparator) + "L" + std::to_string(L) + std::string(kPathSeparator) ;
	std::filesystem::create_directories(folder);

	/* open file */
	string filename = folder + "'Dist_p" + to_string(rho) + "L" + to_string(L) +  "T" + to_string(mcsteps) + ".txt";
	std::ofstream file;
	file.open(filename, ios::out);
	if (!file.is_open()) {
		cout << "Couldn't open a file : " + filename << endl;
		exit(-1);
	}
	cout << "Saving data to : " + filename << endl;

	/* saving to file */
	for (int s = 0; s < ns.size(); s++) {
//if (ns[s] > 0.0) {
			file << s << separator << ns[s] << endl;
//	}
		//else continue;
		
	}
	file.close();
}

