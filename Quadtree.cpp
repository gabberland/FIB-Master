#include <iostream>
#include "Quadtree.h"
#include "timing.h"


QuadtreeNode::~QuadtreeNode()
{
	for(unsigned int i=0; i<4; i++)
		delete children[i];
}

void QuadtreeNode::subdivide(int levels)
{
	if((levels > 0) && (points.size() > 0))
	{
		float sx, sy;
		unsigned int node;
		
		sx = (minCoords.x + maxCoords.x) / 2.0f;
		sy = (minCoords.y + maxCoords.y) / 2.0f;
		for(unsigned int i=0; i<4; i++)
		{
			children[i] = new QuadtreeNode();
			for(unsigned int j=0; j<4; j++)
				children[i]->children[j] = NULL;
		}
		for(unsigned int j=0; j<points.size(); j++)
		{
			node = 0;
			node += (points[j].x < sx)?0:1;
			node += (points[j].y < sy)?0:2;
			children[node]->points.push_back(points[j]);
			children[node]->normals.push_back(normals[j]);
		}
		children[0]->minCoords = minCoords;
		children[0]->maxCoords = glm::vec2(sx, sy);
		children[1]->minCoords = glm::vec2(sx, minCoords.y);
		children[1]->maxCoords = glm::vec2(maxCoords.x, sy);
		children[2]->minCoords = glm::vec2(minCoords.x, sy);
		children[2]->maxCoords = glm::vec2(sx, maxCoords.y);
		children[3]->minCoords = glm::vec2(sx, sy);
		children[3]->maxCoords = maxCoords;

		for(unsigned int i=0; i<4; i++)
			children[i]->subdivide(levels-1);
	}
}

bool QuadtreeNode::isLeaf() const
{
	return (children[0] == NULL) && (children[1] == NULL) && (children[2] == NULL) && (children[3] == NULL);
}

QuadtreeNode *QuadtreeNode::pointToCell(const glm::vec2 &P)
{
	if(isLeaf())
		return this;
	else
	{
		unsigned int node = 0;
		node += (P.x < children[0]->maxCoords.x)?0:1;
		node += (P.y < children[0]->maxCoords.y)?0:2;
		return children[node]->pointToCell(P);
	}
}

void QuadtreeNode::collectCorners(Quadtree *qtree)
{
	glm::vec2 corner;
	unsigned int cornerIndex;
	map<int, int>:: iterator it;
	
	if(isLeaf())
	{
		corner = minCoords;
		cornerIndex = qtree->pointToInteger(corner);
		it = qtree->cornerToUnknown.find(cornerIndex);
		if(it == qtree->cornerToUnknown.end())
		{
			qtree->cornerToUnknown[cornerIndex] = qtree->nUnknowns;
			cornerUnknowns[0] = qtree->nUnknowns++;
		}
		else
			cornerUnknowns[0] = qtree->cornerToUnknown[cornerIndex];

		corner = glm::vec2(maxCoords.x, minCoords.y);
		cornerIndex = qtree->pointToInteger(corner);
		it = qtree->cornerToUnknown.find(cornerIndex);
		if(it == qtree->cornerToUnknown.end())
		{
			qtree->cornerToUnknown[cornerIndex] = qtree->nUnknowns;
			cornerUnknowns[1] = qtree->nUnknowns++;
		}
		else
			cornerUnknowns[1] = qtree->cornerToUnknown[cornerIndex];

		corner = glm::vec2(minCoords.x, maxCoords.y);
		cornerIndex = qtree->pointToInteger(corner);
		it = qtree->cornerToUnknown.find(cornerIndex);
		if(it == qtree->cornerToUnknown.end())
		{
			qtree->cornerToUnknown[cornerIndex] = qtree->nUnknowns;
			cornerUnknowns[2] = qtree->nUnknowns++;
		}
		else
			cornerUnknowns[2] = qtree->cornerToUnknown[cornerIndex];

		corner = maxCoords;
		cornerIndex = qtree->pointToInteger(corner);
		it = qtree->cornerToUnknown.find(cornerIndex);
		if(it == qtree->cornerToUnknown.end())
		{
			qtree->cornerToUnknown[cornerIndex] = qtree->nUnknowns;
			cornerUnknowns[3] = qtree->nUnknowns++;
		}
		else
			cornerUnknowns[3] = qtree->cornerToUnknown[cornerIndex];
	}
	else
	{
		for(unsigned int i=0; i<4; i++)
			children[i]->collectCorners(qtree);
	}
}

void QuadtreeNode::draw(Image &image) const
{
	if(isLeaf())
	{
		glm::ivec2 P, Q;
		
		P = glm::ivec2(minCoords * glm::vec2(image.width(), image.height()));
		Q = glm::ivec2(maxCoords * glm::vec2(image.width(), image.height()));
		P = glm::max(P, glm::ivec2(0, 0));
		Q = glm::min(Q, glm::ivec2(image.width()-1, image.height()-1));
		image.drawRectangle(P.x, P.y, Q.x, Q.y, glm::vec3(1.0f, 1.0f, 1.0f));
		
		for(unsigned int i=0; i<points.size(); i++)
			image.drawFilledCircle((unsigned int)(image.width() * points[i].x), (unsigned int)(image.height() * points[i].y), 4, glm::vec3(0.75f, 0.75f, 0.75f));
	}
	else
	{
		for(unsigned int i=0; i<4; i++)
			children[i]->draw(image);
	}
}


Quadtree::Quadtree()
{
	root = NULL;
	nUnknowns = 0;
	pW = gW = sW = 1.0f;
}

Quadtree::~Quadtree()
{
	if(root != NULL)
		delete root;
}

void Quadtree::setWeights(double pointEqWeight, double gradientEqWeight, double smoothnessEqWeight)
{
	pW = pointEqWeight;
	gW = gradientEqWeight;
	sW = smoothnessEqWeight;
}

void Quadtree::compute(const data_representation::Mesh &cloud, unsigned int levels, ScalarField &field)
{
	nUnknowns = 0;
	nLevels = levels;
	field.init((1 << levels) + 1, (1 << levels) + 1);

	size_t cloudSize = cloud.vertices_.size();

	// Create quatree
	cout << "[INFO] Generating Quadtree" << endl;
	root = new QuadtreeNode();
	for(unsigned int i=0; i<4; i++)
		root->children[i] = NULL;
	root->points.resize(cloudSize);
	for(unsigned int i=0; i<cloudSize; i++)
		root->points[i] = glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y());
	root->normals.resize(cloudSize);
	for(unsigned int i=0; i<cloudSize; i++)
		root->normals[i] = glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y());
	
	root->minCoords = glm::vec2(0.0f, 0.0f);
	root->maxCoords = glm::vec2(1.0f, 1.0f);
	
	root->subdivide(levels);
	
	// Collect corners
	root->collectCorners(this);
	
	cout << "[DATA] # Unknowns = " << nUnknowns << endl;
	
	// Prepare linear system
	unsigned int nEquations;
	
	cout << "[INFO] Preparing Quadtree the system" << endl;
	long lastTime = getTimeMilliseconds();
	
	nEquations = 0;

	vector<Eigen::Triplet<double>> triplets;
	vector<float> bCoeffs;
	
	for(unsigned int i=0; i<cloudSize; i++)
	{
		addPointEquation(nEquations, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), triplets, bCoeffs);
		nEquations++;
	}
	for(unsigned int i=0; i<cloudSize; i++)
	{
		addGradientEquations(nEquations, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), 
							glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), triplets, bCoeffs);
		nEquations += 2;
	}
	/*
	for(unsigned int j=0; j<field.height(); j++)
		for(unsigned int i=1; i<(field.width()-1); i++)
	{
		addHorizontalBoundarySmoothnessEquation(eqIndex, field, i, j, triplets, b);
		nEquations++;
	}
	for(unsigned int j=1; j<(field.height()-1); j++)
		for(unsigned int i=0; i<field.width(); i++)
	{
		addVerticalBoundarySmoothnessEquation(eqIndex, field, i, j, triplets, b);
		nEquations++;
	}
	*/
	cout << nEquations << " equations and " << nUnknowns << " unknowns" << endl;
	
	Eigen::SparseMatrix<double> A(nEquations, nUnknowns), AtA;
	Eigen::VectorXd b(nEquations), Atb, x;

	A.setFromTriplets(triplets.begin(), triplets.end());
	for(unsigned int i=0; i<bCoeffs.size(); i++)
		b(i) = bCoeffs[i];

	cout << "Building A & b in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "Computing normal equation" << endl;
	lastTime = getTimeMilliseconds();
	
	AtA = A.transpose() * A;
	Atb = A.transpose() * b;
	
	cout << "Computed AtA & Atb in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "Solving least squares" << endl;
	lastTime = getTimeMilliseconds();
	
	// TODO
}

void Quadtree::draw(Image &image)
{
	if(image.width() == 0)
		return;
	image.fill(glm::vec3(0.0f, 0.0f, 0.0f));
	root->draw(image);
	/*
	for(map<int, int>::iterator it=cornerToUnknown.begin(); it!=cornerToUnknown.end(); it++)
	{
		int cornerIndex = it->first;
		float y = float(cornerIndex / ((1 << nLevels) + 1)) / (1 << nLevels);
		float x = float(cornerIndex % ((1 << nLevels) + 1)) / (1 << nLevels);
		image.drawFilledCircle(int(image.width() * x), int(image.height() * y), 4, glm::vec3(0.75f, 0.75f, 0.75f));
	}
	*/
}

// From: Point in sampling space
// To:   Quadtree node that contains it

QuadtreeNode *Quadtree::pointToCell(const glm::vec2 &P)
{
	return root->pointToCell(P);
}

// From: Point on a corner of the limit uniform grid
// To:   Unique integer id

unsigned int Quadtree::pointToInteger(const glm::vec2 &P) const
{
	return (unsigned int)round((1 << nLevels) * P.y) * ((1 << nLevels) + 1) + (unsigned int)round((1 << nLevels) * P.x);
}

void Quadtree::addPointEquation(unsigned int eqIndex, const glm::vec2 &P, vector<Eigen::Triplet<double>> &triplets, vector<float> &bCoeffs, const float &value)
{
	QuadtreeNode *node;
	float x, y;
	
	node = pointToCell(P);
	x = (P.x - node->minCoords.x) / (node->maxCoords.x - node->minCoords.x);
	y = (P.y - node->minCoords.y) / (node->maxCoords.y - node->minCoords.y);
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[3], pW*x*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[2], pW*(1.0f-x)*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[1], pW*x*(1.0f-y)));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[0], pW*(1.0f-x)*(1.0f-y)));
	bCoeffs.push_back(value);
}

void Quadtree::addGradientEquations(unsigned int eqIndex, const glm::vec2 &P, const glm::vec2 &N, vector<Eigen::Triplet<double>> &triplets, vector<float> &bCoeffs)
{
	QuadtreeNode *node;
	unsigned int i, j;
	float x, y;
	
	node = pointToCell(P);
	uint nodeWidth  = node->maxCoords.x - node->minCoords.x;
	uint nodeHeight = node->maxCoords.y - node->minCoords.y;
	i = (unsigned int)floor(P.x * (nodeWidth - 1));
	i = glm::max(0u, glm::min(i, nodeWidth-2));
	j = (unsigned int)floor(P.y * (nodeHeight - 1));
	j = glm::max(0u, glm::min(j, nodeHeight-2));
	
	x = P.x * (nodeWidth - 1) - float(i);
	y = P.y * (nodeHeight - 1) - float(j);
	
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[3], gW*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[2], -gW*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[1], gW*(1.0f-y)));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[0], -gW*(1.0f-y)));
	
	bCoeffs.push_back(gW*N.x);

	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, node->cornerUnknowns[3], gW*x));
	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, node->cornerUnknowns[2], gW*(1.0f-x)));
	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, node->cornerUnknowns[1], -gW*x));
	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, node->cornerUnknowns[0], -gW*(1.0f-x)));
	
	bCoeffs.push_back(gW*N.y);
}
void Quadtree::addSamplingEquations(const data_representation::Mesh &cloud, unsigned int eqIndex, const glm::vec2 &P, const glm::vec2 &N, vector<Eigen::Triplet<double>> &triplets, vector<float> &bCoeffs)
{
	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		QuadtreeNode *node = pointToCell(P);
		uint nodeWidth  = node->maxCoords.x - node->minCoords.x;
		uint nodeHeight = node->maxCoords.y - node->minCoords.y;
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), triplets, bCoeffs);
		eqIndex++;
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()) + (1.0f /  (1 << nLevels)) * glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), triplets, bCoeffs, 1.0f);
		eqIndex++;
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()) - (1.0f / (1 << nLevels)) * glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), triplets, bCoeffs, -1.0f);
		eqIndex++;
	}
}

int Quadtree::addHorizontalBoundarySmoothnessEquation(unsigned int eqIndex, const glm::vec2 &P, vector<Eigen::Triplet<double>> &triplets, vector<float> &bCoeffs, const float &value)
{
	QuadtreeNode *node = pointToCell(P);
	

	//Construir estructura unknowns index by position (menstres construim octree setCorners() )

	// Comprobar que realment es pot fer la equacio
	// no pot haver ni x min ni x max
	// si true retorn 0

	// si false calcula eq i retorna 1
		//

/*
	//Casi identica
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i, j), -sW*2.0f));

	//Comprovar si son corners
	//si no ho son, traduir en quin node quadtree es troba -> interpolacio similar point eq

	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i+1, j), sW*1.0f));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, unknownIndex(field, i-1, j), sW*1.0f));
*/	
	bCoeffs.push_back(0.0f);

	return 0;
}





