#ifndef _BIHARMONIC_SOLVER_INCLUDE
#define _BIHARMONIC_SOLVER_INCLUDE


#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "mesh.h"
#include "ScalarField.h"

enum Normal 	{gradient, sampling};
enum Smoothness {singleDimension, twoDimension};
enum Solver {BiCGSTAB, ConjugateGradient, LeastSquaresConjugateGradient};

static const char* NORMAL_STRING [] = {"Gradient", "Sampling"};
static const char* SMOOTHNESS_STRING [] = {"1D", "2D"};
static const char* SOLVER_STRING [] = {"BiCGSTAB", "ConjugateGradient", "LeastSquaresConjugateGradient"};
static const char* RESULT_STRING [] = {"No", "Yes"};

struct SolverData
{
	int 		resolution;
	float		gaussianNoise;
	std::string normalAlgorithm;
	std::string smoothingAlgorithm;
	std::string solverMethod;
	std::string isMultigrid;
	int			multigridIterations;
	int			numberThreads;
	time_t 		systemBuildTime;
	time_t 		matrixBuildTime = 0;
	time_t 		solverResolutionTime;
	uint32_t 	iterations;
	bool		isSolved = 1;
	double		error;	
};

class BiharmonicSolver
{
public:
	BiharmonicSolver();

	void setWeights(double pointEqWeight, double gradientEqWeight, double smoothnessEqWeight);
	SolverData computeWith(const data_representation::Mesh &cloud, ScalarField &field, const int &normal_type, const int &smoothness_type, const int &solver_method, const int &numThreads, const bool &printLogs);
	SolverData computeMultigrid(const data_representation::Mesh &cloud, ScalarField &field, int iterations, const int &normal_type, const int &smoothness_type, const int &solver_method, const int &numThreads, const bool &printLogs);

	static Eigen::VectorXd computeConjugateGradient(const Eigen::SparseMatrix<double> &AtA, const Eigen::VectorXd &Atb, const bool &printLogs, SolverData &outData);
	static Eigen::VectorXd computeBiconjugateGradient(const Eigen::SparseMatrix<double> &AtA, const Eigen::VectorXd &Atb, const bool &printLogs, SolverData &outData);
	static Eigen::VectorXd computeLeastSquares(const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b, const bool &printLogs, SolverData &outData);

private:
	void addPointEquation(unsigned int eqIndex, const glm::vec2 &P, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b, float value = 0.0f);
	void addGradientEquations(unsigned int eqIndex, const glm::vec2 &P, const glm::vec2 &N, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addPointAndGradientEquations(unsigned int eqIndex, const glm::vec2 &P, const glm::vec2 &N, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addSmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addBilaplacianSmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);

	void addHorizontalBoundarySmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addVerticalBoundarySmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	
	unsigned int unknownIndex(const ScalarField &field, unsigned int i, unsigned int j) const;

	void addPointAndGradientEquations(const data_representation::Mesh &cloud, ScalarField &field, unsigned int &eqIndex, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addSamplingEquations(const data_representation::Mesh &cloud, ScalarField &field, unsigned int &eqIndex, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addOneDimensionalSmoothnessEquations(unsigned int &eqIndex, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addTwoDimensionalSmoothnessEquations(unsigned int &eqIndex, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);

	Eigen::VectorXd computeConjugateGradientWithGuess(const Eigen::SparseMatrix<double> &AtA, const Eigen::VectorXd &Atb, const Eigen::VectorXd &guess, const bool &printLogs, SolverData &outData);
	Eigen::VectorXd computeBiconjugateGradientWithGuess(const Eigen::SparseMatrix<double> &AtA, const Eigen::VectorXd &Atb, const Eigen::VectorXd &guess, const bool &printLogs, SolverData &outData);
	Eigen::VectorXd computeLeastSquaresWithGuess(const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b, const Eigen::VectorXd &guess, const bool &printLogs, SolverData &outData);


	Eigen::VectorXd divideField(const ScalarField &field);
private:
	double pW, gW, sW;

};

#endif 