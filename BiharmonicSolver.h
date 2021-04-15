#ifndef _BIHARMONIC_SOLVER_INCLUDE
#define _BIHARMONIC_SOLVER_INCLUDE


#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "mesh.h"
#include "ScalarField.h"

enum Normal 	{gradient, sampling};
enum Smoothness {singleDimension, twoDimension};

class BiharmonicSolver
{
public:
	BiharmonicSolver();

	void setWeights(double pointEqWeight, double gradientEqWeight, double smoothnessEqWeight);
	void compute(const data_representation::Mesh &cloud, ScalarField &field);
	void computeComponentWise(const data_representation::Mesh &cloud, ScalarField &field);
	void computeBilaplacian(const data_representation::Mesh &cloud, ScalarField &field);
	void computeNoGradient(const data_representation::Mesh &cloud, ScalarField &field);
	void computeWith(const data_representation::Mesh &cloud, ScalarField &field, const int &normal_type, const int &smoothness_type);

private:
	void addPointEquation(unsigned int eqIndex, const glm::vec2 &P, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b, float value = 0.0f);
	void addGradientEquations(unsigned int eqIndex, const glm::vec2 &P, const glm::vec2 &N, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addSmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addBilaplacianSmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);

	void addHorizontalBoundarySmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addVerticalBoundarySmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	
	unsigned int unknownIndex(const ScalarField &field, unsigned int i, unsigned int j) const;

	void addGradientEquations(const data_representation::Mesh &cloud, ScalarField &field, unsigned int &eqIndex, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addSamplingEquations(const data_representation::Mesh &cloud, ScalarField &field, unsigned int &eqIndex, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addOneDimensionalSmoothnessEquations(unsigned int &eqIndex, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);
	void addTwoDimensionalSmoothnessEquations(unsigned int &eqIndex, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b);

private:
	double pW, gW, sW;

};


#endif 


