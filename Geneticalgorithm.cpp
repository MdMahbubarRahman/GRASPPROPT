#include "Geneticalgorithm.h"



//default constructor
Chromosome::Chromosome() {
	fitness					  = 0.0;
	sepInt					  = 0;
	isFeasible				  = false;
	maxRouteCapacity		  = 10;
}

//copy constructor
Chromosome::Chromosome(const Chromosome& chrom) {
	chromosome				  = chrom.chromosome;
	fitness					  = chrom.fitness;
	sepInt					  = chrom.sepInt;
	isFeasible				  = chrom.isFeasible;
	maxRouteCapacity		  = chrom.maxRouteCapacity;
}

//constructs chromosome from elements
Chromosome::Chromosome(std::vector<int> chrom, double fitNess, int spInt, bool isFeb, int maxCap) {
	chromosome = chrom;
	fitness = fitNess;
	sepInt = spInt;
	isFeasible = isFeb;
	maxRouteCapacity = maxCap;
}


//returns chromosome
std::vector<int> Chromosome::getChromosomeRepresentation() {
	return chromosome;
}

//returns fitness value
double Chromosome::getFitness() {
	return fitness;
}

//returns seperator integer value
int Chromosome::getSepInt() {
	return sepInt;
}

//returns feasibility status
bool Chromosome::getFeasibilityStatus() {
	return isFeasible;
}

//returns size of the chromosome
int Chromosome::getSize() {
	return chromosome.size();
}

//returns route capacity
int Chromosome::getMaxRouteCapacity() {
	return maxRouteCapacity;
}

//prints the chromosome
void Chromosome::showChromosome() {
	std::cout << "\nThe Chromosome representation is : " << std::endl;
	for (auto it : chromosome) {
		std::cout << it << " ";
	}
	std::cout << ";" << std::endl;
	std::cout << "The fitness value is : " << fitness << ";" << std::endl;
	std::cout << "The seperator integer is : " << sepInt << ";" << std::endl;
	std::cout << "The feasibility status is : " << isFeasible << ";" << std::endl;
	std::cout << "The max route capacity is : " << maxRouteCapacity << ";" << std::endl;
}

//updates to feasible chromosome
void Chromosome::updateToFeasibleChromosome(std::vector<int> demand, std::vector<std::vector<double>> distance) {
	if (isFeasible) {
		std::cout << "\nThe chromosome representation is already feasible!" << std::endl;
	}
	else {
		std::cout << "\nThe chromosome is not feasible! Feasibility restoration is in progress!!!" << std::endl;
		std::vector<int> infSol = chromosome;
		std::vector<int> tempStorage;
		chromosome.clear(); 
		int capacity = 0;
		int numSepInt = 0;
		for (auto it = infSol.begin(); it != infSol.end(); ++it) {
			if (*it != sepInt) {
				capacity += demand.at(*it);
				if (capacity <= maxRouteCapacity) {
					tempStorage.push_back(*it);
				}
				else {
					for (auto iter : tempStorage) {
						chromosome.push_back(iter);
					}
					if (numSepInt >= 1) {
						for (int i = 0; i < numSepInt; ++i) {
							chromosome.push_back(sepInt);
						}
						numSepInt = 0;
					}
					else {
						chromosome.push_back(sepInt);
						numSepInt -= 1;
					}					
					tempStorage.clear();
					capacity = demand.at(*it);
					tempStorage.push_back(*it);
				}
			}
			else {
				numSepInt += 1;
			}
		}
		if (!tempStorage.empty()) {
			for (auto it : tempStorage) {
				chromosome.push_back(it);
			}
		}
		if (numSepInt>=1) {
			for (int i = 0; i < numSepInt; ++i)
				chromosome.push_back(sepInt);
		}
		else {
			chromosome.push_back(sepInt);
		}
		// update fitness
		double cost = 0.0;
		int preval = 0;
		int val = 0;
		for (auto it = chromosome.begin(); it != chromosome.end(); ++it) {
			if (*it == sepInt) {
				val = 0;
			}
			else {
				val = *it;
			}
			cost += distance[preval][val];
			preval = val;
		}
		cost += distance[preval][0];
		fitness = cost;
		//update feasibility
		isFeasible = true;
	}
}

//default constructor
CrossOver::CrossOver() {
	std::cout << "\nThe default constructor for cross over class is called!" << std::endl;
}

//copy constructor
CrossOver::CrossOver(const CrossOver & crossr) {
	std::cout << "\nThe copy constructor of crossover class is called!" << std::endl;
	parent1		= crossr.parent1;
	parent1		= crossr.parent2;
	offspring	= crossr.offspring;
}

//crossover class construction with parents
CrossOver::CrossOver(Chromosome par1, Chromosome par2) {
	parent1 = par1;
	parent2 = par2;
}

//performs partially mapped crossover operation
void CrossOver::performPertiallyMappedCrossover(std::vector<int> demand, std::vector<std::vector<double>> distance) {
	std::vector<int> par1 = parent1.getChromosomeRepresentation();
	std::vector<int> par2 = parent2.getChromosomeRepresentation();
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distr(int(par1.size()/2), (par1.size()-1));
	std::map<int, int> dict;
	int randVar = distr(gen);
	int sepInt = parent1.getSepInt();
	//dictionary	
	for (int i = 0; i < par2.size(); ++i) {
		if (par2.at(i) != sepInt) {
			dict.insert(std::pair<int, int>(par2.at(i), i));
		}
	}
	//partially mapped crossover
	for (int i = randVar; i < par1.size(); ++i) {
		int val1 = par1.at(i);
		int val2 = par2.at(i);
		if (val1 != val2 && val1 != sepInt && val2 != sepInt) {
			int indx = dict[val1];
			par2[indx] = val2;
			par2[i] = val1;
			dict.erase(val1);
			dict.erase(val2);
			dict.insert(std::pair<int, int>(val2, indx));
			dict.insert(std::pair<int, int>(val1, i));
		}
		else if (val1 != val2 && val1 == sepInt) {
			par2[i] = val1;
			for (int j = 0; j < par2.size(); ++j) {
				if (par2.at(j) == 10) {
					dict.erase(val2);
					dict.insert(std::pair<int, int>(val2, j));
					par2[j] = val2;
					break;
				}
			}
		}
		else if (val1 != val2 && val2 == sepInt) {
			int indx = dict[val1];
			par2[i] = val1;
			dict.erase(val1);
			dict.insert(std::pair<int, int>(val1, i));
			par2[indx] = val2;
		}
		else
			continue;
	}
	//fitness
	double cost = 0.0;
	int preval = 0;
	int val = 0;
	for (auto it = par2.begin(); it != par2.end(); ++it) {
		if (*it == sepInt) {
			val = 0;
		}
		else {
			val = *it;
		}
		cost += distance[preval][val];
		preval = val;
	}
	cost += distance[preval][0];
	//check feasibility
	int capacity = 0;
	bool feasible = true;
	for (auto it : par2) {
		if (it != sepInt) {
			capacity += demand[it];
			if (capacity >= parent1.getMaxRouteCapacity()) {
				feasible = false;
				break;
			}
		}
		else {
			capacity = 0;
		}
	}
	//offspring after crossover
	Chromosome child(par2, cost, sepInt, feasible, parent1.getMaxRouteCapacity());
	if(!child.getFeasibilityStatus()){
		child.updateToFeasibleChromosome(demand, distance);
	}
	offspring = child;
}

//prints parent 1
void CrossOver::showParent1() {
	std::cout << "The first parent is : " << std::endl;//need to catch in buffer
	parent1.showChromosome();
}

//prints parent 2
void CrossOver::showParent2() {
	std::cout << "The second parent is : " << std::endl;
	parent2.showChromosome();
}

//prints offspring
void CrossOver::showOffspring() {
	std::cout << "The offspring is : " << std::endl;
	offspring.showChromosome();
}

//returns parent 1
Chromosome CrossOver::getParent1() {
	return parent1;
}

//returns parent 2
Chromosome CrossOver::getParent2() {
	return parent2;
}

//returns offspring
Chromosome CrossOver::getOffspring() {
	return offspring;
}

//default constructor
Mutation::Mutation() {
	std::cout << "The default mutation constructor is called!" << std::endl;
}
	
//copy constructor
Mutation::Mutation(const Mutation & mutatn){
	std::cout << "The copy constructor of the mutation class is called!" << std::endl;
	originalChrom = mutatn.originalChrom;
	mutatedChrom = mutatn.mutatedChrom;
}

//mutation constructs with original chromosome
Mutation::Mutation(Chromosome origChrom) {
	originalChrom = origChrom;
}

//returns premutation chromosome
Chromosome Mutation::getOriginalChromosome() {
	return originalChrom;
}

//returns mutated chromosome
Chromosome Mutation::getMutatedOffspring() {
	return mutatedChrom;
}

//prints premutation chromosome
void Mutation::showOriginalChromosome() {
	std::cout << "The premuated chromosome is : " << std::endl;
	originalChrom.showChromosome();
}

//prints mutated chromosome
void Mutation::showMutatedChromosome() {
	std::cout << "The mutated chromosome is : " << std::endl;
	mutatedChrom.showChromosome();
}

//performs chain mutation
void Mutation::performChainMutation(std::vector<int> demand, std::vector<std::vector<double>> distance) {
	std::vector<int> chromToMutat = getOriginalChromosome().getChromosomeRepresentation();
	int size = getOriginalChromosome().getSize();
	int sepInt = getOriginalChromosome().getSepInt();
	int capacityLimit = getOriginalChromosome().getMaxRouteCapacity();
	int val = ceil(sqrt(size));
	//no. of genes
	std::random_device rd;
	std::mt19937 gen(rd());
	int numGene = 2;
	if (val>2) {
		std::uniform_int_distribution<> distr(2, val);
		numGene = distr(gen);
	}
	//choose the genes
	std::vector<int> geneIndxVec;
	bool flag = false;
	while(geneIndxVec.size() != numGene) {
		std::uniform_int_distribution<> distr(0, (size-1));
		int geneIndx = distr(gen);
		//std::cout << "the indx is : " << geneIndx << std::endl;
		if (chromToMutat[geneIndx] != sepInt && chromToMutat[geneIndx] != 0) {
			if (geneIndxVec.size() > 0) {
				for (int i = 0; i < geneIndxVec.size(); ++i) {
					if (geneIndxVec.at(i) == geneIndx) {
						flag = true;
						break;
					}
				}
			}
			if (flag == false) {
				geneIndxVec.push_back(geneIndx);
			}
			flag = false;
		}
	}
	//perform mutation
	int prev = chromToMutat[geneIndxVec.at(0)];
	int current = 0;
	for (int i = 1; i < geneIndxVec.size(); ++i) {
		current = chromToMutat[geneIndxVec.at(i)];
		chromToMutat[geneIndxVec.at(i)] = prev;
		prev = current;
	}
	chromToMutat[geneIndxVec.at(0)] = prev;
	//calculate fitness/cost
	double cost = 0;
	int preVal = 0;
	int vall = 0;
	for (auto& it : chromToMutat) {
		if (it != sepInt) {
			vall = it;
		}
		else {
			vall = 0;
		}
		cost = cost + distance[preVal][vall];
		//std::cout << preVal << " " << val << " " << distance[preVal][vall] << " " << std::endl;
		preVal = vall;
	}
	cost = cost + distance[preVal][0];
	//check feasibility
	int capacity = 0;
	bool feasible = true;
	for (auto it : chromToMutat) {
		if (it != sepInt) {
			capacity += demand[it];
			if (capacity > capacityLimit) {
				feasible = false;
				break;
			}
		}
		else {
			capacity = 0;
		}
	}
	//get the mutated child
	Chromosome child(chromToMutat, cost, sepInt, feasible, capacityLimit);
	if (!child.getFeasibilityStatus()) {
		child.updateToFeasibleChromosome(demand, distance);
	}
	mutatedChrom = child;
}

//comparator
bool Comparator::operator()(Chromosome &a, Chromosome &b){
	return (a.getFitness() > b.getFitness());
}


//default constructor
Population::Population() {
	std::cout << "The default constructor of the population class is called!" << std::endl;
	populationSize = 50;
	numberOfNodes  = 0;
	depotNode      = 0;
	capacityLimit  = 0;
	kChainLength   = 0;
	sWapLength     = 0;
}

//copy constructor
Population::Population(const Population & gentn) {
	std::cout << "The copy constructor of the population class is called!" << std::endl;
	populationSize = gentn.populationSize;
	numberOfNodes  = gentn.numberOfNodes;
	depotNode      = gentn.depotNode;
	capacityLimit  = gentn.capacityLimit;
	kChainLength   = gentn.kChainLength;
	sWapLength     = gentn.sWapLength;
	population     = gentn.population;
	demand		   = gentn.demand;
	distance	   = gentn.distance;
	offspring	   = gentn.offspring;
	populationBest = gentn.populationBest;
	crossOverChild = gentn.crossOverChild;
	mutationChild  = gentn.mutationChild;
}

//construct population with parameters
Population::Population(int popSize, int numNodes, int depNode, int capLimit, int kChain, int sWap, std::vector<int> dem, std::vector<std::vector<double>> dist) {
	populationSize = popSize;
	numberOfNodes  = numNodes;
	depotNode      = depNode;
	capacityLimit  = capLimit;
	kChainLength   = kChain;
	sWapLength     = sWap;
	demand         = dem;
	distance       = dist;
}


//return population
std::priority_queue<Chromosome, std::vector<Chromosome>, Comparator> Population::getPopulation() {
	return population;
}

//print population
void Population::showPopulation() {
	std::priority_queue<Chromosome, std::vector<Chromosome>, Comparator> abc = population;
	std::cout << "The chromosomes of the current population are : " << std::endl;
	int i = 1;
	while (!abc.empty())
	{
		std::cout << "\nThe chromosome no. is : " << i << std::endl;
		Chromosome crm = abc.top();
		crm.showChromosome();
		abc.pop();
		i++;
	}
}

//updates population with the offspring
void Population::updatePopulationWithOffspring() {
	population.push(offspring);
}

//returns offspring after genetic operations 
Chromosome Population::getOffspring() {
	return offspring;
}

//removes clones with random Chromosome
void Population::manageClones() {
	//needs implementation
}

//returns generation best chromosome
Chromosome Population::getBestChromosome() {
	populationBest = population.top();
	return populationBest;
}

//returns offspring after cross over
Chromosome Population::getCrossOverChild() {
	return crossOverChild;
}

//returns offspring after mutation
Chromosome Population::getMutationChild() {
	return mutationChild;
}

//generates a random chromosome
Chromosome Population::generateRandomChromosome() {
	std::random_device rd;
	std::mt19937 gen(rd());
	//populate chromosome representation
	int totalDemand = 0;
	for (int i = 0; i < demand.size(); ++i) {
		totalDemand += demand.at(i);
	}
	int numDepotNodes = ceil(double(totalDemand) / double(capacityLimit));
	std::list<int> nodeBox;
	for (int i = 0; i < numberOfNodes; ++i)
		nodeBox.push_back(i);
	for (int j = 0; j < (4 * numDepotNodes - 1); ++j)
		nodeBox.push_back(depotNode);
	std::vector<int> chrom_container;
	while (!nodeBox.empty()) {
		auto it = nodeBox.begin();
		std::uniform_int_distribution<> distr(0, (nodeBox.size() - 1));
		int num = distr(gen);
		std::advance(it, num);
		chrom_container.push_back(*it);
		nodeBox.erase(it);
	}
	chrom_container.push_back(depotNode);
	//calculate fitness/cost
	double cost = 0;
	int preVal = 0;
	int val = 0;
	for (auto& it : chrom_container) {
		if (it != depotNode) {
			val = it;
		}
		else {
			val = 0;
		}
		cost = cost + distance[preVal][val];
		//std::cout << preVal << " " << val << " " << distance[preVal][val] << " " << std::endl;
		preVal = val;
	}
	cost = cost + distance[preVal][0];
	//check feasibility
	int capacity = 0;
	bool feasible = true;
	for (auto it : chrom_container) {
		if (it != depotNode) {
			capacity += demand[it];
			if (capacity > capacityLimit) {
				feasible = false;
				break;
			}
		}
		else {
			capacity = 0;
		}
	}
	Chromosome crm(chrom_container, cost, depotNode, feasible, capacityLimit);
	if (!crm.getFeasibilityStatus()) {
		crm.updateToFeasibleChromosome(demand, distance);
	}
	return crm;
}

//populates first generation of population
void Population::populatePopulation() {
	for (int i = 0; i < populationSize; ++i) {
		Chromosome crm = generateRandomChromosome();
		population.push(crm);
	}
}

//perfoms crossover on the current population
void Population::performCrossOver() {
	Chromosome parent1 = population.top();
	population.pop();
	Chromosome parent2 = population.top();
	population.push(parent1);
	CrossOver cross(parent1, parent2);
	cross.performPertiallyMappedCrossover(demand, distance);
	crossOverChild = cross.getOffspring();
}

//performs mutation on the current population
void Population::performMutation() {
	Mutation mut(crossOverChild);
	mut.performChainMutation(demand, distance);
	mutationChild = mut.getMutatedOffspring();
	offspring = mutationChild;//needs work
}


//default constructor
Geneticalgorithm::Geneticalgorithm() {
	std::cout << "The default constructor of the genetic algorithm has been called!" << std::endl;
	maxIterations			= 50;
	populationSize			= 0;
	numberOfNodes			= 0;
	depotNode				= 0;
	capacityLimit			= 0;
	kChainLength			= 0;
	sWapLength				= 0;
}

//copy constructor	
Geneticalgorithm::Geneticalgorithm(const Geneticalgorithm & ga) {
	std::cout << "The default constructor of the genetic algorithm has been called!" << std::endl;
	maxIterations			= ga.maxIterations;
	populationSize			= ga.populationSize;
	numberOfNodes			= ga.numberOfNodes;
	depotNode				= ga.depotNode;
	capacityLimit			= ga.capacityLimit;
	kChainLength		    = ga.kChainLength;
	sWapLength				= ga.sWapLength;
	ppl						= ga.ppl;
	initialSolution			= ga.incumbentSolution;
	incumbentSolution		= ga.incumbentSolution;
	bestSolution			= ga.bestSolution;
	generationBestSolutions = ga.generationBestSolutions;
	generationalOffsprings	= ga.generationalOffsprings;
}

//construct ga with initial values
Geneticalgorithm::Geneticalgorithm(int popSize, int numNodes, int depNode, int capLimit, int kChain, int sWap, std::vector<int> dem, std::vector<std::vector<double>> dist) {
	maxIterations			= 50;
	populationSize			= popSize;
	numberOfNodes		    = numNodes;
	depotNode				= depNode;
	capacityLimit			= capLimit;
	kChainLength			= kChain;
	sWapLength				= sWap;
	demand					= dem;
	distance				= dist;
}


//populate initial generation
void Geneticalgorithm::populateInitialGeneration() {
	Population population(populationSize, numberOfNodes, depotNode, capacityLimit, kChainLength, sWapLength, demand, distance);
	ppl = population;
	ppl.populatePopulation();
}

//returns best solutions of every generation
std::vector<Chromosome> Geneticalgorithm::getGenerationBestSolutions() {
	return generationBestSolutions;
}

//returns offsprings of every generation
std::vector<FeasibleSolution> Geneticalgorithm::getGenerationalOffsprings() {
	return generationalOffsprings;
}

//returns ga solution
Chromosome Geneticalgorithm::getGASolution() {
	return bestSolution;
}

//runs ga algorithm
void Geneticalgorithm::runGeneticAlgorithm() {
	//needs implementation
}

//prints ga solution
void Geneticalgorithm::showGASolution() {
	bestSolution.showChromosome();
}

//converst a chromosome to a feasible solution 
FeasibleSolution Geneticalgorithm::getFeasibleSolutionFromChromosome(Chromosome chrom) {
	//needs implementation 
	FeasibleSolution febSol;
	return febSol;
}

//converts feasible solution to a chromosome
Chromosome Geneticalgorithm::getChromosomeFromFeasibleSolution(FeasibleSolution febSol) {
	//needs implementation
	Chromosome crm;
	return crm;
}

//prints current generations chromosomes
void Geneticalgorithm::showCurrentGeneration() {
	ppl.showPopulation();
}

