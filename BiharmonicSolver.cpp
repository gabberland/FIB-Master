#include <iostream>
#include <cmath>
#include <omp.h>

#include "BiharmonicSolver.h"
#include "timing.h"


BiharmonicSolver::BiharmonicSolver()
{
	pW = gW = sW = 1.0f;
}


void BiharmonicSolver::setWeights(double pointEqWeight, double gradientEqWeight, double smoothnessEqWeight)
{
	pW = pointEqWeight;
	gW = gradientEqWeight;
	sW = smoothnessEqWeight;
}

void BiharmonicSolver::compute(const data_representation::Mesh &cloud, ScalarField &field)
{
	unsigned int nEquations, nUnknowns;
	
	cout << "Preparing the system" << endl;
	long lastTime = getTimeMilliseconds();
	
	nEquations = 3 * cloud.vertices_.size() + (field.width() - 2) * (field.height() - 2);
	nEquations += 2 * field.width() + 2 * field.height() - 8;
	nUnknowns = field.width() * field.height();
	
	cout << nEquations << " equations and " << nUnknowns << " unknowns" << endl;

	Eigen::SparseMatrix<double> A(nEquations, nUnknowns), AtA;
	Eigen::VectorXd b(nEquations), Atb, x;
	vector<Eigen::Triplet<double>> triplets;
	unsigned int eqIndex = 0;

	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), field, triplets, b);
		eqIndex++;
	}
	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		addPointAndGradientEquations(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), field, triplets, b);
		eqIndex += 2;
	}
	for(unsigned int j=1; j<(field.height()-1); j++)
		for(unsigned int i=1; i<(field.width()-1); i++)
		{
			addSmoothnessEquation(eqIndex, field, i, j, triplets, b);
			eqIndex++;
		}

	
	for(unsigned int i=1; i<(field.width()-1); i++)
	{
		addHorizontalBoundarySmoothnessEquation(eqIndex, field, i, 0, triplets, b);
		eqIndex++;
		addHorizontalBoundarySmoothnessEquation(eqIndex, field, i, field.height()-1, triplets, b);
		eqIndex++;
	}
	for(unsigned int j=1; j<(field.height()-1); j++)
	{
		addVerticalBoundarySmoothnessEquation(eqIndex, field, 0, j, triplets, b);
		eqIndex++;
		addVerticalBoundarySmoothnessEquation(eqIndex, field, field.width()-1, j, triplets, b);
		eqIndex++;
	}
	
	
	cout << "Total equations added: " << eqIndex << endl;
	
	A.setFromTriplets(triplets.begin(), triplets.end());

	cout << "Building A & b in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "Computing normal equation" << endl;
	lastTime = getTimeMilliseconds();
	
	AtA = A.transpose() * A;
	Atb = A.transpose() * b;
	
	cout << "Computed AtA & Atb in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "Solving least squares" << endl;
	lastTime = getTimeMilliseconds();

	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> cg;
	cg.setTolerance(1e-5);
	cg.compute(AtA);
	Eigen::ComputationInfo info = cg.info();
	if(info!=Eigen::Success)
		cout << "Decomposition failed!!!" << endl;	
	x = cg.solve(Atb);
	info = cg.info();
	if(info!=Eigen::Success)
		cout << "Solving failed!!!" << endl;

	cout << "Least squares solved in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;
	
	double relative_error = (A*x - b).norm() / b.norm();
	cout << "The relative error is: " << relative_error << endl;
	cout << "Error: " << cg.error() << endl;
	cout << "Iterations: " << cg.iterations() << endl;
		
	for(unsigned int j=0, pos=0; j<field.height(); j++)
		for(unsigned int i=0; i<field.width(); i++, pos++)
			field(i, j) = x(pos);
}

void BiharmonicSolver::computeComponentWise(const data_representation::Mesh &cloud, ScalarField &field)
{
	unsigned int nEquations, nUnknowns;
	
	cout << "Preparing the system" << endl;
	long lastTime = getTimeMilliseconds();
	
	nEquations = 3 * cloud.vertices_.size() + 2 * field.width() * field.height() - 2 * field.width() - 2 * field.height();
	nUnknowns = field.width() * field.height();
	
	cout << nEquations << " equations and " << nUnknowns << " unknowns" << endl;

	Eigen::SparseMatrix<double> A(nEquations, nUnknowns), AtA;
	Eigen::VectorXd b(nEquations), Atb, x;
	vector<Eigen::Triplet<double>> triplets;
	unsigned int eqIndex = 0;

	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), field, triplets, b);
		eqIndex++;
	}
	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		addPointAndGradientEquations(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), field, triplets, b);
		eqIndex += 2;
	}
	for(unsigned int j=0; j<field.height(); j++)
		for(unsigned int i=1; i<(field.width()-1); i++)
	{
		addHorizontalBoundarySmoothnessEquation(eqIndex, field, i, j, triplets, b);
		eqIndex++;
	}
	for(unsigned int j=1; j<(field.height()-1); j++)
		for(unsigned int i=0; i<field.width(); i++)
	{
		addVerticalBoundarySmoothnessEquation(eqIndex, field, i, j, triplets, b);
		eqIndex++;
	}
	
	cout << "Total equations added: " << eqIndex << endl;
	
	A.setFromTriplets(triplets.begin(), triplets.end());

	cout << "Building A & b in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "Computing normal equation" << endl;
	lastTime = getTimeMilliseconds();
	
	AtA = A.transpose() * A;
	Atb = A.transpose() * b;
	
	cout << "Computed AtA & Atb in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "Solving least squares" << endl;
	lastTime = getTimeMilliseconds();

	
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> cg;
	cg.setTolerance(1e-5);
	cg.compute(AtA);
	Eigen::ComputationInfo info = cg.info();
	if(info!=Eigen::Success)
		cout << "Decomposition failed!!!" << endl;	
	x = cg.solve(Atb);
	info = cg.info();
	if(info!=Eigen::Success)
		cout << "Solving failed!!!" << endl;
	
	cout << "Least squares solved in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;
	
	double relative_error = (A*x - b).norm() / b.norm();
	cout << "The relative error is: " << relative_error << endl;
	cout << "Error: " << cg.error() << endl;
	cout << "Iterations: " << cg.iterations() << endl;
		
	for(unsigned int j=0, pos=0; j<field.height(); j++)
		for(unsigned int i=0; i<field.width(); i++, pos++)
			field(i, j) = x(pos);
}

void BiharmonicSolver::computeBilaplacian(const data_representation::Mesh &cloud, ScalarField &field)
{
	unsigned int nEquations, nUnknowns;
	
	cout << "Preparing the system" << endl;
	long lastTime = getTimeMilliseconds();
	
	nEquations = 3 * cloud.vertices_.size() + (field.width() - 2) * (field.height() - 2) + (field.width() - 4) * (field.height() - 4);
	nUnknowns = field.width() * field.height();
	
	cout << nEquations << " equations and " << nUnknowns << " unknowns" << endl;

	Eigen::SparseMatrix<double> A(nEquations, nUnknowns), AtA;
	Eigen::VectorXd b(nEquations), Atb, x;
	vector<Eigen::Triplet<double>> triplets;
	unsigned int eqIndex = 0;

	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), field, triplets, b);
		eqIndex++;
	}
	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		addPointAndGradientEquations(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), field, triplets, b);
		eqIndex += 2;
	}
	for(unsigned int j=1; j<(field.height()-1); j++)
		for(unsigned int i=1; i<(field.width()-1); i++)
		{
			addSmoothnessEquation(eqIndex, field, i, j, triplets, b);
			eqIndex++;
		}
	for(unsigned int j=2; j<(field.height()-2); j++)
		for(unsigned int i=2; i<(field.width()-2); i++)
		{
			addBilaplacianSmoothnessEquation(eqIndex, field, i, j, triplets, b);
			eqIndex++;
		}
	
	cout << "Total equations added: " << eqIndex << endl;
	
	A.setFromTriplets(triplets.begin(), triplets.end());

	cout << "Building A & b in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "Computing normal equation" << endl;
	lastTime = getTimeMilliseconds();
	
	AtA = A.transpose() * A;
	Atb = A.transpose() * b;
	
	cout << "Computed AtA & Atb in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "Solving least squares" << endl;
	lastTime = getTimeMilliseconds();

	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> cg;
	cg.setTolerance(1e-5);
	cg.compute(AtA);
	Eigen::ComputationInfo info = cg.info();
	if(info!=Eigen::Success)
		cout << "Decomposition failed!!!" << endl;	
	x = cg.solve(Atb);
	info = cg.info();
	if(info!=Eigen::Success)
		cout << "Solving failed!!!" << endl;

	cout << "Least squares solved in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;
	
	double relative_error = (A*x - b).norm() / b.norm();
	cout << "The relative error is: " << relative_error << endl;
	cout << "Error: " << cg.error() << endl;
	cout << "Iterations: " << cg.iterations() << endl;
		
	for(unsigned int j=0, pos=0; j<field.height(); j++)
		for(unsigned int i=0; i<field.width(); i++, pos++)
			field(i, j) = x(pos);
}

void BiharmonicSolver::computeNoGradient(const data_representation::Mesh &cloud, ScalarField &field)
{
	unsigned int nEquations, nUnknowns;
	
	cout << "Preparing the system" << endl;
	long lastTime = getTimeMilliseconds();
	
	nEquations = 3 * cloud.vertices_.size() + 2 * field.width() * field.height() - 2 * field.width() - 2 * field.height();
	nUnknowns = field.width() * field.height();
	
	cout << nEquations << " equations and " << nUnknowns << " unknowns" << endl;

	Eigen::SparseMatrix<double> A(nEquations, nUnknowns), AtA;
	Eigen::VectorXd b(nEquations), Atb, x;
	vector<Eigen::Triplet<double>> triplets;
	unsigned int eqIndex = 0;

	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), field, triplets, b);
		eqIndex++;
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()) + (1.0f / field.width()) * glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), field, triplets, b, 1.0f);
		eqIndex++;
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()) - (1.0f / field.height()) * glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), field, triplets, b, -1.0f);
		eqIndex++;
	}
	for(unsigned int j=0; j<field.height(); j++)
		for(unsigned int i=1; i<(field.width()-1); i++)
	{
		addHorizontalBoundarySmoothnessEquation(eqIndex, field, i, j, triplets, b);
		eqIndex++;
	}
	for(unsigned int j=1; j<(field.height()-1); j++)
		for(unsigned int i=0; i<field.width(); i++)
	{
		addVerticalBoundarySmoothnessEquation(eqIndex, field, i, j, triplets, b);
		eqIndex++;
	}
	
	cout << "Total equations added: " << eqIndex << endl;
	
	A.setFromTriplets(triplets.begin(), triplets.end());

	cout << "Building A & b in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "Computing normal equation" << endl;
	lastTime = getTimeMilliseconds();
	
	AtA = A.transpose() * A;
	Atb = A.transpose() * b;
	
	cout << "Computed AtA & Atb in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "Solving least squares" << endl;
	lastTime = getTimeMilliseconds();

	
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> cg;
	cg.setTolerance(1e-5);
	cg.compute(AtA);
	Eigen::ComputationInfo info = cg.info();
	if(info!=Eigen::Success)
		cout << "Decomposition failed!!!" << endl;	
	x = cg.solve(Atb);
	info = cg.info();
	if(info!=Eigen::Success)
		cout << "Solving failed!!!" << endl;
	
	cout << "Least squares solved in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;
	
	double relative_error = (A*x - b).norm() / b.norm();
	cout << "The relative error is: " << relative_error << endl;
	cout << "Error: " << cg.error() << endl;
	cout << "Iterations: " << cg.iterations() << endl;
		
	for(unsigned int j=0, pos=0; j<field.height(); j++)
		for(unsigned int i=0; i<field.width(); i++, pos++)
			field(i, j) = x(pos);
}

void BiharmonicSolver::addPointEquation(unsigned int eqIndex, const glm::vec2 &P, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b, float value)
{
	unsigned int i, j;
	float x, y;
	
	i = (unsigned int)floor(P.x * (field.width() - 1));
	i = glm::max(0u, glm::min(i, field.width()-2));
	j = (unsigned int)floor(P.y * (field.height() - 1));
	j = glm::max(0u, glm::min(j, field.height()-2));
	
	x = P.x * (field.width() - 1) - float(i);
	y = P.y * (field.height() - 1) - float(j);
	
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i+1, j+1), pW*x*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j+1), pW*(1.0f-x)*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i+1, j), pW*x*(1.0f-y)));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j), pW*(1.0f-x)*(1.0f-y)));
	
	b(eqIndex) = value;
}

void BiharmonicSolver::addPointAndGradientEquations(unsigned int eqIndex, const glm::vec2 &P, const glm::vec2 &N, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b)
{
	unsigned int i, j;
	float x, y;
	
	i = (unsigned int)floor(P.x * (field.width() - 1));
	i = glm::max(0u, glm::min(i, field.width()-2));
	j = (unsigned int)floor(P.y * (field.height() - 1));
	j = glm::max(0u, glm::min(j, field.height()-2));
	
	x = P.x * (field.width() - 1) - float(i);
	y = P.y * (field.height() - 1) - float(j);
	
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i+1, j+1), gW*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j+1), -gW*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i+1, j), gW*(1.0f-y)));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j), -gW*(1.0f-y)));
	
	b(eqIndex) = gW* 64 / (field.width()-1)*N.x;

	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, unknownIndex(field, i+1, j+1), gW*x));
	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, unknownIndex(field, i, j+1), gW*(1.0f-x)));
	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, unknownIndex(field, i+1, j), -gW*x));
	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, unknownIndex(field, i, j), -gW*(1.0f-x)));
	
	b(eqIndex+1) = gW* 64 / (field.height()-1) * N.y;
}

void BiharmonicSolver::addSmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b)
{
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j), -sW*4.0f));

	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i+1, j), sW*1.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i-1, j), sW*1.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j+1), sW*1.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j-1), sW*1.0f));
	
	b(eqIndex) = 0.0f;
}

void BiharmonicSolver::addBilaplacianSmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b)
{
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j), sW*20.0f/16.0f));

	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i+1, j), -sW*8.0f/16.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i-1, j), -sW*8.0f/16.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j+1), -sW*8.0f/16.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j-1), -sW*8.0f/16.0f));

	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i+1, j+1), sW*2.0f/16.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i+1, j-1), sW*2.0f/16.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i-1, j+1), sW*2.0f/16.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i-1, j-1), sW*2.0f/16.0f));

	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i+2, j), sW*1.0f/16.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i-2, j), sW*1.0f/16.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j+2), sW*1.0f/16.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j-2), sW*1.0f/16.0f));
	
	b(eqIndex) = 0.0f;
}

void BiharmonicSolver::addHorizontalBoundarySmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b)
{
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j), -sW * (field.width()-1)/64 * 2.0f));

	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i+1, j), sW * (field.width()-1)/64 * 1.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i-1, j), sW * (field.width()-1)/64 * 1.0f));
	
	b(eqIndex) = 0.0f;
}

void BiharmonicSolver::addVerticalBoundarySmoothnessEquation(unsigned int eqIndex, const ScalarField &field, unsigned int i, unsigned int j, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b)
{
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j), -sW * (field.height()-1)/64 * 2.0f));

	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j+1), sW * (field.height()-1)/64 * 1.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j-1), sW * (field.height()-1)/64 * 1.0f));
	
	b(eqIndex) = 0.0f;
}

unsigned int BiharmonicSolver::unknownIndex(const ScalarField &field, unsigned int i, unsigned int j) const
{
	return j * field.width() + i;
}

void BiharmonicSolver::addPointAndGradientEquations(const data_representation::Mesh &cloud, ScalarField &field, unsigned int &eqIndex, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b)
{
	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), field, triplets, b);
		eqIndex++;
	}
	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		addPointAndGradientEquations(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), field, triplets, b);
		eqIndex += 2;
	}
}

void BiharmonicSolver::addSamplingEquations(const data_representation::Mesh &cloud, ScalarField &field, unsigned int &eqIndex, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b)
{
	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), field, triplets, b);
		eqIndex++;
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()) + (1.0f / 64) * glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), field, triplets, b, 1.0f);
		eqIndex++;
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()) - (1.0f / 64) * glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), field, triplets, b, -1.0f);
		eqIndex++;
	}
}

void BiharmonicSolver::addOneDimensionalSmoothnessEquations(unsigned int &eqIndex, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b)
{
	for(unsigned int j=0; j<field.height(); j++)
		for(unsigned int i=1; i<(field.width()-1); i++)
		{
			addHorizontalBoundarySmoothnessEquation(eqIndex, field, i, j, triplets, b);
			eqIndex++;
		}

	for(unsigned int j=1; j<(field.height()-1); j++)
		for(unsigned int i=0; i<field.width(); i++)
		{
			addVerticalBoundarySmoothnessEquation(eqIndex, field, i, j, triplets, b);
			eqIndex++;
		}

}

void BiharmonicSolver::addTwoDimensionalSmoothnessEquations(unsigned int &eqIndex, const ScalarField &field, vector<Eigen::Triplet<double>> &triplets, Eigen::VectorXd &b)
{	
	for(unsigned int j=1; j<(field.height()-1); j++)
		for(unsigned int i=1; i<(field.width()-1); i++)
		{
			addSmoothnessEquation(eqIndex, field, i, j, triplets, b);
			eqIndex++;
		}	

	for(unsigned int i=1; i<(field.width()-1); i++)
	{
		addHorizontalBoundarySmoothnessEquation(eqIndex, field, i, 0, triplets, b);
		eqIndex++;
		addHorizontalBoundarySmoothnessEquation(eqIndex, field, i, field.height()-1, triplets, b);
		eqIndex++;
	}
	for(unsigned int j=1; j<(field.height()-1); j++)
	{
		addVerticalBoundarySmoothnessEquation(eqIndex, field, 0, j, triplets, b);
		eqIndex++;
		addVerticalBoundarySmoothnessEquation(eqIndex, field, field.width()-1, j, triplets, b);
		eqIndex++;
	}
}

SolverData BiharmonicSolver::computeWith(const data_representation::Mesh &cloud, ScalarField &field, const int &normal_type, const int &smoothness_type, const int &solver_method, const int &numThreads, const bool &printLogs)
{
	// Init data
	//
	unsigned int nEquations, nUnknowns;
	SolverData outData;
	
	if(printLogs) cout << "[MESSAGE] Preparing the system with " << normal_type << " normal equations and " << smoothness_type << " smoothing equation ..." << endl;
	long lastTime = getTimeMilliseconds();
	
	// Compute number of equations depending on the smoothing algorithm
	//
	if(smoothness_type == Smoothness::singleDimension)
	{
		nEquations = 3 * cloud.vertices_.size() + 2 * field.width() * field.height() - 2 * field.width() - 2 * field.height();

	}
	else if (smoothness_type == Smoothness::twoDimension)
	{
		nEquations = 3 * cloud.vertices_.size() + (field.width() - 2) * (field.height() - 2);
		nEquations += 2 * field.width() + 2 * field.height() - 8;
	}	

	// Compute number of unkowns
	//
	nUnknowns = field.width() * field.height();
	
	if(printLogs) cout << "[MESSAGE] Added " << nEquations << " equations and " << nUnknowns << " unknowns" << endl;

	// Init data containers
	//
	Eigen::SparseMatrix<double> A(nEquations, nUnknowns), AtA;
	Eigen::VectorXd b(nEquations), Atb, x;
	vector<Eigen::Triplet<double>> triplets;
	unsigned int eqIndex = 0;

	// Add equations with the choosen configuration
	//
	if (normal_type == Normal::gradient)
	{
		addPointAndGradientEquations(cloud, field, eqIndex, triplets, b);
	}

	else if(normal_type == Normal::sampling)
	{
		addSamplingEquations(cloud, field, eqIndex, triplets, b);
	}

	if(smoothness_type == Smoothness::singleDimension)
	{
		addOneDimensionalSmoothnessEquations(eqIndex, field, triplets, b);
	} 
	
	else if(smoothness_type == Smoothness::twoDimension)
	{
		addTwoDimensionalSmoothnessEquations(eqIndex, field, triplets, b);
	}
	
	if(printLogs) cout << "[DATA] Total equations added: " << eqIndex << endl;
	
	A.setFromTriplets(triplets.begin(), triplets.end());

	// Compute matrix build time
	//
	time_t sysTime = getTimeMilliseconds() - lastTime;
	outData.systemBuildTime = sysTime;
	if(printLogs) cout << "[TIME] Building A & b in " << sysTime << " ms" << endl;

	if(printLogs) cout << "[MESSAGE] Computing normal equation..." << endl;
	lastTime = getTimeMilliseconds();
	
	// Compute transposed matrices
	//
	AtA = A.transpose() * A;
	Atb = A.transpose() * b;
	
	time_t buildTime = getTimeMilliseconds() - lastTime;
	if(solver_method == Solver::LeastSquaresConjugateGradient)
	{
		buildTime = 0;
	}
	outData.matrixBuildTime = buildTime;
	if(printLogs) cout << "[TIME]  AtA & Atb in " << buildTime << " ms" << endl;

	if(printLogs) cout << "[MESSAGE] Solving System ..." << endl;

	int n = Eigen::nbThreads();
	if(printLogs) std::cout << "[INFO] Initial number Threads: " << n << std::endl;

	Eigen::initParallel();
	omp_set_num_threads(numThreads);
	Eigen::setNbThreads(numThreads);
	n = Eigen::nbThreads();
	if(printLogs) std::cout << "[INFO] Number Threads used: " << n << std::endl;
		
	if(solver_method == Solver::ConjugateGradient)
	{
		x = computeConjugateGradient(AtA, Atb, printLogs, outData);
	}
	else if(solver_method == Solver::BiCGSTAB)
	{
		x = computeBiconjugateGradient(AtA, Atb, printLogs, outData);
	}
	else if(solver_method == Solver::LeastSquaresConjugateGradient)
	{
		x = computeLeastSquares(A, b, printLogs, outData);
	}
	double relative_error = (A*x - b).norm() / b.norm();
	if(printLogs) cout << "[DATA] The relative error is: " << relative_error << endl;
		
	for(unsigned int j=0, pos=0; j<field.height(); j++)
		for(unsigned int i=0; i<field.width(); i++, pos++)
			field(i, j) = x(pos);

	return outData;
}

SolverData BiharmonicSolver::computeMultigrid(const data_representation::Mesh &cloud, ScalarField &field, int iterations, const int &normal_type, const int &smoothness_type, const int &solver_method, const int &numThreads, const bool &printLogs)
{

	if(printLogs) cout << "[MESSAGE] First will generate a base solution of field size: " << field.width() << endl;

	// Solve first grid
	computeWith(cloud, field, normal_type, smoothness_type, solver_method, numThreads, 0);

	if(printLogs) cout << "[MESSAGE] Base solution generated!" << endl;

	SolverData outData;
	ScalarField actualField = field;
	
	while (iterations > 0)
	{
		// Compute new field size & guess
		size_t newFieldSize = std::max(actualField.height()-1, actualField.width()-1) * 2 + 1;
		if(printLogs) cout << endl << "[MESSAGE] Dividing new field with size: " << newFieldSize << endl;
		Eigen::VectorXd guess = divideField(actualField);

		if(printLogs) cout << "[MESSAGE] New Guess field generated " << endl;

		// Initialize new grid field
		//
		field.init(newFieldSize, newFieldSize);

		unsigned int nEquations, nUnknowns;
		
		if(printLogs) cout << "[MESSAGE] Preparing the system with " << normal_type << " normal equations and " << smoothness_type << " smoothing equation ..." << endl;
		long lastTime = getTimeMilliseconds();
		
		if(smoothness_type == Smoothness::singleDimension)
		{
			nEquations = 3 * cloud.vertices_.size() + 2 * field.width() * field.height() - 2 * field.width() - 2 * field.height();

		}
		else if (smoothness_type == Smoothness::twoDimension)
		{
			nEquations = 3 * cloud.vertices_.size() + (field.width() - 2) * (field.height() - 2);
			nEquations += 2 * field.width() + 2 * field.height() - 8;
		}	

		nUnknowns = field.width() * field.height();
		
		if(printLogs) cout << "[MESSAGE] Added " << nEquations << " equations and " << nUnknowns << " unknowns" << endl;

		Eigen::SparseMatrix<double> A(nEquations, nUnknowns), AtA;
		Eigen::VectorXd b(nEquations), Atb, x;
		vector<Eigen::Triplet<double>> triplets;
		unsigned int eqIndex = 0;

		if (normal_type == Normal::gradient)
		{
			addPointAndGradientEquations(cloud, field, eqIndex, triplets, b);
		}

		else if(normal_type == Normal::sampling)
		{
			addSamplingEquations(cloud, field, eqIndex, triplets, b);
		}

		if(smoothness_type == Smoothness::singleDimension)
		{
			addOneDimensionalSmoothnessEquations(eqIndex, field, triplets, b);
		} 
		
		else if(smoothness_type == Smoothness::twoDimension)
		{
			addTwoDimensionalSmoothnessEquations(eqIndex, field, triplets, b);
		}
		
		if(printLogs) cout << "[DATA] Total equations added: " << eqIndex << endl;
		
		A.setFromTriplets(triplets.begin(), triplets.end());

		time_t sysTime = getTimeMilliseconds() - lastTime;
		outData.systemBuildTime = sysTime;
		if(printLogs) cout << "[TIME] Building A & b in " << sysTime << " ms" << endl;

		if(printLogs) cout << "[MESSAGE] Computing normal equation..." << endl;
		lastTime = getTimeMilliseconds();
		
		AtA = A.transpose() * A;
		Atb = A.transpose() * b;
		
		time_t buildTime = getTimeMilliseconds() - lastTime;
		outData.matrixBuildTime = buildTime;
		if(printLogs) cout << "[TIME]  AtA & Atb in " << buildTime << " ms" << endl;

		if(printLogs) cout << "[MESSAGE] Solving least squares ..." << endl;
		lastTime = getTimeMilliseconds();

		int n = Eigen::nbThreads();
		if(printLogs) std::cout << "[INFO] Initial number Threads: " << n << std::endl;

		Eigen::initParallel();
		omp_set_num_threads(numThreads);
		Eigen::setNbThreads(numThreads);
		n = Eigen::nbThreads();
		if(printLogs) std::cout << "[INFO] Number Threads used: " << n << std::endl;
		
		if(solver_method == Solver::ConjugateGradient)
		{
			x = computeBiconjugateGradientWithGuess(AtA, Atb, guess, printLogs, outData);
		}
		else if(solver_method == Solver::BiCGSTAB)
		{
			x = computeBiconjugateGradientWithGuess(AtA, Atb, guess, printLogs, outData);
		}
		else if(solver_method == Solver::LeastSquaresConjugateGradient)
		{
			x = computeLeastSquaresWithGuess(A, b, guess, printLogs, outData);
		}
		double relative_error = (A*x - b).norm() / b.norm();
		if(printLogs) cout << "[DATA] The relative error is: " << relative_error << endl;
			
		for(unsigned int j=0, pos=0; j<field.height(); j++)
			for(unsigned int i=0; i<field.width(); i++, pos++)
				field(i, j) = x(pos);

		--iterations;
		actualField = field;
	}
	return outData;
}

Eigen::VectorXd BiharmonicSolver::divideField(const ScalarField &field)
{
	ScalarField dividedField;
	size_t newFieldSize = std::max(field.height()-1, field.width()-1) * 2 + 1;
	dividedField.init(newFieldSize, newFieldSize);

	for(unsigned int j=0, pos=0; j < newFieldSize; j++)
		for(unsigned int i=0; i < newFieldSize; i++)
		{
			//If we move through the original row divisions, we will assign itselves
			if(j % 2 == 0)
			{
				// If we move through the original column divisions, we will assign itselves
				if(i % 2 == 0)
				{
					dividedField(i,j) = field((size_t)i/2, (size_t)j/2);
				}
				// If we move through a new column division, we will interpolate its 2 neighbor values
				else
				{
					float leftValue  = field((size_t)i/2,    (size_t)j/2);
					float rightValue = field((size_t)i/2 + 1, (size_t)j/2);
					dividedField(i,j) = (float)(leftValue+rightValue) / 2;
				}
			}

			else
			{
				// If we move through a new column division, we will floaterpolate its 2 neighbor values
				if(i % 2 == 0)
				{
					float topValue    = field((size_t)i/2, (size_t)j/2);
					float bottomValue = field((size_t)i/2, (size_t)j/2 + 1);
					dividedField(i,j) = (float)(topValue+bottomValue) / 2;
				}
				// If we move through a new row & column subdivision, we will floaterpolate the value of its 4 original neighbors
				else
				{					
					float topLeftValue     = field((size_t)i/2,     (size_t)j/2);
					float topRightValue    = field((size_t)i/2 + 1, (size_t)j/2);
					float bottomLeftValue  = field((size_t)i/2,     (size_t)j/2 + 1);
					float bottomRightValue = field((size_t)i/2 + 1, (size_t)j/2 + 1);
					
					// New field value is the average of the 4 neighbors
					dividedField(i,j) = (float)(topLeftValue+topRightValue+bottomLeftValue+bottomRightValue) / 4;
				}
			}
		}

	Eigen::VectorXd out(newFieldSize*newFieldSize);
	for(unsigned int j=0, pos=0; j<dividedField.height(); j++)
		for(unsigned int i=0; i<dividedField.width(); i++, pos++)
			out(pos) = dividedField(i, j);

	return out;
}

Eigen::VectorXd BiharmonicSolver::computeConjugateGradient(const Eigen::SparseMatrix<double> &AtA, const Eigen::VectorXd &Atb, const bool &printLogs, SolverData &outData)
{
	long lastTime = getTimeMilliseconds();

	// Initialize Conjugate Gradient
	//
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
	
	// Set error tolerance of the solver
	//
	cg.setTolerance(1e-5);

	// Store a reference to the matrix A
	//
	cg.compute(AtA);
	Eigen::ComputationInfo info = cg.info();

	// Check if was succeed
	//
	outData.isSolved = 1;
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Decomposition failed!!!" << endl;	
	}

	// Solve the equation system
	//
	Eigen::VectorXd x = cg.solve(Atb);
	info = cg.info();

	// Check if was succeed
	//
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Solving failed!!!" << endl;
	}
	
	// Compute solving time
	//
	time_t solveTime =(getTimeMilliseconds() - lastTime);
	if(printLogs) cout << "[TIME] Least squares solved in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;
	outData.solverResolutionTime = solveTime;
	
	// Obtain Error & Iterations
	//
	outData.error = cg.error();
	outData.iterations = cg.iterations();
	if(printLogs) cout << "[DATA] Iterations: " << cg.iterations() << endl;
	if(printLogs) cout << "[DATA] Error: " << cg.error() << endl;

	return x;
}

Eigen::VectorXd BiharmonicSolver::computeBiconjugateGradient(const Eigen::SparseMatrix<double> &AtA, const Eigen::VectorXd &Atb, const bool &printLogs, SolverData &outData)
{
	long lastTime = getTimeMilliseconds();

	// Initialize Biconjugate Gradient
	//
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> cg;	
	
	// Set error tolerance of the solver
	//
	cg.setTolerance(1e-5);

	// Store a reference to the matrix A
	//
	cg.compute(AtA);
	Eigen::ComputationInfo info = cg.info();

	// Check if was succeed
	//
	outData.isSolved = 1;
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Decomposition failed!!!" << endl;	
	}

	// Solve the equation system
	//
	Eigen::VectorXd x = cg.solve(Atb);
	info = cg.info();

	// Check if was succeed
	//
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Solving failed!!!" << endl;
	}
	
	// Compute solving time
	//
	time_t solveTime =(getTimeMilliseconds() - lastTime);
	if(printLogs) cout << "[TIME] Least squares solved in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;
	outData.solverResolutionTime = solveTime;
	
	// Obtain Error & Iterations
	//
	outData.error = cg.error();
	outData.iterations = cg.iterations();
	if(printLogs) cout << "[DATA] Iterations: " << cg.iterations() << endl;
	if(printLogs) cout << "[DATA] Error: " << cg.error() << endl;

	return x;
}
Eigen::VectorXd BiharmonicSolver::computeLeastSquares(const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b, const bool &printLogs, SolverData &outData)
{
	long lastTime = getTimeMilliseconds();

	// Initialize Biconjugate Gradient
	//
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>> cg;	
	
	// Set error tolerance of the solver
	//
	cg.setTolerance(1e-5);

	// Store a reference to the matrix A
	//
	cg.compute(A);
	Eigen::ComputationInfo info = cg.info();

	// Check if was succeed
	//
	outData.isSolved = 1;
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Decomposition failed!!!" << endl;	
	}

	// Solve the equation system
	//
	Eigen::VectorXd x = cg.solve(b);
	info = cg.info();

	// Check if was succeed
	//
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Solving failed!!!" << endl;
	}
	
	// Compute solving time
	//
	time_t solveTime =(getTimeMilliseconds() - lastTime);
	if(printLogs) cout << "[TIME] Least squares solved in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;
	outData.solverResolutionTime = solveTime;
	
	// Obtain Error & Iterations
	//
	outData.error = cg.error();
	outData.iterations = cg.iterations();
	if(printLogs) cout << "[DATA] Iterations: " << cg.iterations() << endl;
	if(printLogs) cout << "[DATA] Error: " << cg.error() << endl;

	return x;
}

Eigen::VectorXd BiharmonicSolver::computeConjugateGradientWithGuess(const Eigen::SparseMatrix<double> &AtA, const Eigen::VectorXd &Atb, const Eigen::VectorXd &guess, const bool &printLogs, SolverData &outData)
{
	long lastTime = getTimeMilliseconds();

	// Initialize Conjugate Gradient
	//
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
	
	// Set error tolerance of the solver
	//
	cg.setTolerance(1e-5);

	// Store a reference to the matrix A
	//
	cg.compute(AtA);
	Eigen::ComputationInfo info = cg.info();

	// Check if was succeed
	//
	outData.isSolved = 1;
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Decomposition failed!!!" << endl;	
	}

	// Solve the equation system
	//
	Eigen::VectorXd x = cg.solveWithGuess(Atb, guess);
	info = cg.info();

	// Check if was succeed
	//
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Solving failed!!!" << endl;
	}
	
	// Compute solving time
	//
	time_t solveTime =(getTimeMilliseconds() - lastTime);
	if(printLogs) cout << "[TIME] Least squares solved in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;
	outData.solverResolutionTime = solveTime;
	
	// Obtain Error & Iterations
	//
	outData.error = cg.error();
	outData.iterations = cg.iterations();
	if(printLogs) cout << "[DATA] Iterations: " << cg.iterations() << endl;
	if(printLogs) cout << "[DATA] Error: " << cg.error() << endl;

	return x;
}

Eigen::VectorXd BiharmonicSolver::computeBiconjugateGradientWithGuess(const Eigen::SparseMatrix<double> &AtA, const Eigen::VectorXd &Atb, const Eigen::VectorXd &guess, const bool &printLogs, SolverData &outData)
{
	long lastTime = getTimeMilliseconds();

	// Initialize Biconjugate Gradient
	//
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> cg;	
	
	// Set error tolerance of the solver
	//
	cg.setTolerance(1e-5);

	// Store a reference to the matrix A
	//
	cg.compute(AtA);
	Eigen::ComputationInfo info = cg.info();

	// Check if was succeed
	//
	outData.isSolved = 1;
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Decomposition failed!!!" << endl;	
	}

	// Solve the equation system
	//
	Eigen::VectorXd x = cg.solveWithGuess(Atb, guess);
	info = cg.info();

	// Check if was succeed
	//
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Solving failed!!!" << endl;
	}
	
	// Compute solving time
	//
	time_t solveTime =(getTimeMilliseconds() - lastTime);
	if(printLogs) cout << "[TIME] Least squares solved in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;
	outData.solverResolutionTime = solveTime;
	
	// Obtain Error & Iterations
	//
	outData.error = cg.error();
	outData.iterations = cg.iterations();
	if(printLogs) cout << "[DATA] Iterations: " << cg.iterations() << endl;
	if(printLogs) cout << "[DATA] Error: " << cg.error() << endl;

	return x;
}

Eigen::VectorXd BiharmonicSolver::computeLeastSquaresWithGuess(const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b, const Eigen::VectorXd &guess, const bool &printLogs, SolverData &outData)
{
	long lastTime = getTimeMilliseconds();

	// Initialize Biconjugate Gradient
	//
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>> cg;	
	
	// Set error tolerance of the solver
	//
	cg.setTolerance(1e-5);

	// Store a reference to the matrix A
	//
	cg.compute(A);
	Eigen::ComputationInfo info = cg.info();

	// Check if was succeed
	//
	outData.isSolved = 1;
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Decomposition failed!!!" << endl;	
	}

	// Solve the equation system
	//
	Eigen::VectorXd x = cg.solveWithGuess(b, guess);
	info = cg.info();

	// Check if was succeed
	//
	if(info!=Eigen::Success)
	{
		outData.isSolved = 0;
		if(printLogs) cout << "[ERROR] Solving failed!!!" << endl;
	}
	
	// Compute solving time
	//
	time_t solveTime =(getTimeMilliseconds() - lastTime);
	if(printLogs) cout << "[TIME] Least squares solved in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;
	outData.solverResolutionTime = solveTime;
	
	// Obtain Error & Iterations
	//
	outData.error = cg.error();
	outData.iterations = cg.iterations();
	if(printLogs) cout << "[DATA] Iterations: " << cg.iterations() << endl;
	if(printLogs) cout << "[DATA] Error: " << cg.error() << endl;

	return x;
}

