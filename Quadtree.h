#ifndef _QUADTREE_INCLUDE
#define _QUADTREE_INCLUDE


#include <vector>
#include <map>
#include <glm/glm.hpp>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "ScalarField.h"
#include "mesh.h"

using namespace std;


class Quadtree;


struct QuadtreeNode
{
	~QuadtreeNode();
	
	void subdivide(int levels);
	
	bool isLeaf() const;
	QuadtreeNode *pointToCell(const glm::vec2 &P);
	void collectCorners(Quadtree *qtree);

	void draw(Image &image) const;

public:
	glm::vec2 minCoords, maxCoords;
	vector<glm::vec2> points, normals;
	QuadtreeNode *children[4];
	// Unknown id for each of the node's corners (0-> (minx, miny), 1-> (maxx, miny), 2-> (minx, maxy), 3-> (maxx, maxy))
	int cornerUnknowns[4];
	
};


class Quadtree
{
public:
	Quadtree();
	~Quadtree();
	
	void setWeights(double pointEqWeight, double gradientEqWeight, double smoothnessEqWeight);
	void compute(const data_representation::Mesh &cloud, unsigned int levels, ScalarField &field);
	
	void draw(Image &image);
	
private:
	QuadtreeNode *pointToCell(const glm::vec2 &P);
	unsigned int pointToInteger(const glm::vec2 &P) const;

	void addPointEquation(unsigned int eqIndex, const glm::vec2 &P, vector<Eigen::Triplet<double>> &triplets, vector<float> &bCoeffs);
	void addGradientEquations(unsigned int eqIndex, const glm::vec2 &P, const glm::vec2 &N, vector<Eigen::Triplet<double>> &triplets, vector<float> &bCoeffs);

private:
	unsigned int nLevels, nUnknowns;
	QuadtreeNode *root;
	// Maps unique corner id to unknown id in the linear system
	map<int, int> cornerToUnknown;
	double pW, gW, sW;
	
	friend class QuadtreeNode;

};


#endif


