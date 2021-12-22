#include <iostream>
#include <set>
#include "Quadtree.h"
#include "timing.h"


QuadtreeNode::QuadtreeNode()
{
	for(unsigned int i=0; i<4; i++)
		children[i] = NULL; 
}
QuadtreeNode::~QuadtreeNode()
{
	for(unsigned int i=0; i<4; i++)
		delete children[i];
}

void QuadtreeNode::subdivide(int levels, const bool &fullGridSubdivision)
{
	if((levels > 0) && (points.size() > 0 || fullGridSubdivision))
	{
		float sx, sy;
		unsigned int node;
		int newNodeCornerX, newNodeCornerY;
		
		sx = (minCoords.x + maxCoords.x) / 2.0f;
		sy = (minCoords.y + maxCoords.y) / 2.0f;
		
		newNodeCornerX = (minNodeCorner.x + maxNodeCorner.x) / 2;
		newNodeCornerY = (minNodeCorner.y + maxNodeCorner.y) / 2;

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

		children[0]->minNodeCorner = minNodeCorner;
		children[0]->maxNodeCorner = glm::vec2(newNodeCornerX, newNodeCornerY);
		children[1]->minNodeCorner = glm::vec2(newNodeCornerX, minNodeCorner.y);
		children[1]->maxNodeCorner = glm::vec2(maxNodeCorner.x, newNodeCornerY);
		children[2]->minNodeCorner = glm::vec2(minNodeCorner.x, newNodeCornerY);
		children[2]->maxNodeCorner = glm::vec2(newNodeCornerX, maxNodeCorner.y);
		children[3]->minNodeCorner = glm::vec2(newNodeCornerX, newNodeCornerY);
		children[3]->maxNodeCorner = maxNodeCorner;

		//std::cout << "MAX NODE CORNER (" << maxNodeCorner.x << ", " << maxNodeCorner.y << ")" << std::endl; 

		for(unsigned int i=0; i<4; i++)
			children[i]->subdivide(levels-1, fullGridSubdivision);
	}
}

bool QuadtreeNode::isLeaf() const
{
	return (children[0] == NULL) && (children[1] == NULL) && (children[2] == NULL) && (children[3] == NULL);
}

QuadtreeNode *QuadtreeNode::pointToCell(const glm::vec2 &P)
{
	assert(P.x >= minCoords.x && P.y >= minCoords.y && P.x <= maxCoords.x && P.y <= maxCoords.y);

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

void QuadtreeNode::collectCorners(Quadtree *qtree, std::set<std::pair<int,int>> & cornersArray)
{
	glm::ivec2 corner;
	unsigned int cornerIndex;
	map<int, int>:: iterator it;
	
	corner = maxNodeCorner;
	//std::cout << "INTERCORNER INSERTED: (" << corner.x << ", " << corner.y << ")" << std::endl;
	//std::cout << (maxNodeCorner.x - minNodeCorner.x) << std::endl ;
	//std::cout << (maxNodeCorner.y - minNodeCorner.y) << std::endl ;

	if(isLeaf())
	{
		//No hauria de ser corner vec2 sino ivec2
		corner = minNodeCorner;

		cornerIndex = qtree->cornerToInteger(corner);
		qtree->indexToCorner[cornerIndex] = corner;

		it = qtree->cornerIndexToUnknown.find(cornerIndex);
		if(it == qtree->cornerIndexToUnknown.end())
		{
			qtree->cornerIndexToUnknown[cornerIndex] = qtree->nUnknowns;
			qtree->unknownToCornerIndex[qtree->nUnknowns] = cornerIndex;
			cornersArray.insert(std::pair<int,int>(corner.x, corner.y));
			//std::cout << "CORNER INSERTED: (" << corner.x << ", " << corner.y << ")" << std::endl;
			cornerUnknowns[0] = qtree->nUnknowns++;
		}
		else
			cornerUnknowns[0] = qtree->cornerIndexToUnknown[cornerIndex];

		corner = glm::vec2(maxNodeCorner.x, minNodeCorner.y);
		cornerIndex = qtree->cornerToInteger(corner);
		it = qtree->cornerIndexToUnknown.find(cornerIndex);
		if(it == qtree->cornerIndexToUnknown.end())
		{
			qtree->cornerIndexToUnknown[cornerIndex] = qtree->nUnknowns;
			qtree->unknownToCornerIndex[qtree->nUnknowns] = cornerIndex;
			cornersArray.insert(std::pair<int,int>(corner.x, corner.y));
			//std::cout << "CORNER INSERTED: (" << corner.x << ", " << corner.y << ")" << std::endl;
			cornerUnknowns[1] = qtree->nUnknowns++;
		}
		else
			cornerUnknowns[1] = qtree->cornerIndexToUnknown[cornerIndex];

		corner = glm::vec2(minNodeCorner.x, maxNodeCorner.y);
		cornerIndex = qtree->cornerToInteger(corner);
		it = qtree->cornerIndexToUnknown.find(cornerIndex);
		if(it == qtree->cornerIndexToUnknown.end())
		{
			qtree->cornerIndexToUnknown[cornerIndex] = qtree->nUnknowns;
			qtree->unknownToCornerIndex[qtree->nUnknowns] = cornerIndex;
			cornersArray.insert(std::pair<int,int>(corner.x, corner.y));
			//std::cout << "CORNER INSERTED: (" << corner.x << ", " << corner.y << ")" << std::endl;
			cornerUnknowns[2] = qtree->nUnknowns++;
		}
		else
			cornerUnknowns[2] = qtree->cornerIndexToUnknown[cornerIndex];

		corner = maxNodeCorner;
		cornerIndex = qtree->cornerToInteger(corner);
		it = qtree->cornerIndexToUnknown.find(cornerIndex);
		if(it == qtree->cornerIndexToUnknown.end())
		{
			qtree->cornerIndexToUnknown[cornerIndex] = qtree->nUnknowns;
			qtree->unknownToCornerIndex[qtree->nUnknowns] = cornerIndex;
			cornersArray.insert(std::pair<int,int>(corner.x, corner.y));
			//std::cout << "CORNER INSERTED: (" << corner.x << ", " << corner.y << ")" << std::endl;
			cornerUnknowns[3] = qtree->nUnknowns++;
		}
		else
			cornerUnknowns[3] = qtree->cornerIndexToUnknown[cornerIndex];
	}
	else
	{
		for(unsigned int i=0; i<4; i++)
			children[i]->collectCorners(qtree, cornersArray);
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

SolverData Quadtree::compute(const data_representation::Mesh &cloud, unsigned int levels, ScalarField &field, const int &normal_type, const int &smoothness_type, const int &solver_method, const int &numThreads, const bool &fullGridSubdivision, const bool &printLogs)
{
    if(printLogs) cout << "[MESSAGE] Preparing the system with " << NORMAL_STRING[normal_type] << " Normal equations and " << SMOOTHNESS_STRING[smoothness_type] << " Smoothing equations and " << SOLVER_STRING[solver_method] << "Solver ..." << endl;
        
	SolverData outData;
	nUnknowns = 0;
	nLevels = levels;

	field.init(((1 << levels) + 1), ((1 << levels) + 1));

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
	
	// Init root min max
	//
	root->minCoords = glm::vec2(0.0f, 0.0f);
	root->maxCoords = glm::vec2(1.0f, 1.0f);
	
	// Integer ids for the grid nodes
	//
	root->minNodeCorner = glm::ivec2(0, 0);
	root->maxNodeCorner = glm::ivec2(std::pow(2, levels), std::pow(2, levels));

	// Subdivide levels
	//
	root->subdivide(levels, fullGridSubdivision);
	
	// Collect corners
	//
	std::set<std::pair<int,int>> cornersArray;
	root->collectCorners(this, cornersArray);
	
	cout << "[DATA] # Unknowns = " << nUnknowns << endl;
	
	// Prepare linear system
	unsigned int nEquations;
	
	cout << "[INFO] Preparing Quadtree the system..." << endl;
	long lastTime = getTimeMilliseconds();
	
	nEquations = 0;

	vector<Eigen::Triplet<double>> triplets;
	vector<float> bCoeffs;	

	for(unsigned int i=0; i<cloudSize; i++)
	{
		addPointEquation(nEquations, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), triplets, bCoeffs);
		nEquations++;
	}
	cout << "[INFO] Point Eq Added" << endl;
    if (normal_type == Normal::gradient)
    {
        for(unsigned int i=0; i<cloudSize; i++)
        {
            addGradientEquations(nEquations, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()), 
                                glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), triplets, bCoeffs);
            nEquations += 2;
        }
    }
    if (normal_type == Normal::sampling)
    {
        for(unsigned int i=0; i<cloudSize; i++)
        {
            addPointEquation(nEquations, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()) + (1.0f /  (1 << nLevels)) * glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), triplets, bCoeffs, 1.0f);
            nEquations++;
            addPointEquation(nEquations, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()) - (1.0f / (1 << nLevels)) * glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), triplets, bCoeffs, -1.0f);
            nEquations++;
        }
    }
	cout << "[INFO] Normal Eq Added" << endl;
	std::set<std::pair<int, int>>::iterator it;
	for(it = cornersArray.begin(); it != cornersArray.end(); ++it)
	{
		bool valid = addHorizontalBoundarySmoothnessEquation(nEquations, glm::ivec2(it->first, it->second), triplets, bCoeffs);
		if(valid)
			nEquations++;
	}
	cout << "[INFO] Horizontal Eq Added" << endl;
	for(it = cornersArray.begin(); it != cornersArray.end(); ++it)
	{
		bool valid = addVerticalBoundarySmoothnessEquation(nEquations, glm::ivec2(it->first, it->second), triplets, bCoeffs);
		if(valid)
			nEquations++;
	}
	
	cout << "[INFO] Vertical Eq Added" << endl;

	cout << "[INFO]" << nEquations << " equations and " << nUnknowns << " unknowns" << endl;
	//std::cout << "CORNER INDEX NUMBER: " << cornerIndexToUnknown.size() << std::endl;

	Eigen::SparseMatrix<double> A(nEquations, nUnknowns), AtA;
	Eigen::VectorXd b(nEquations), Atb, x;

	A.setFromTriplets(triplets.begin(), triplets.end());
	for(unsigned int i=0; i<bCoeffs.size(); i++)
		b(i) = bCoeffs[i];

	time_t sysTime = getTimeMilliseconds() - lastTime;
	outData.systemBuildTime = sysTime;

	cout << "[INFO] Building A & b in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "[MESSAGE] Computing normal equation" << endl;
	lastTime = getTimeMilliseconds();
	
	// Compute transposed matrices
	//
	AtA = A.transpose() * A;
	Atb = A.transpose() * b;
	
	cout << "[INFO] Computed AtA & Atb in " << (getTimeMilliseconds() - lastTime) << " ms" << endl;

	cout << "[INFO] Solving least squares.." << endl;
	//lastTime = getTimeMilliseconds();
	
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
		x = BiharmonicSolver::computeConjugateGradient(AtA, Atb, printLogs, outData);
	}
	else if(solver_method == Solver::BiCGSTAB)
	{
		x = BiharmonicSolver::computeBiconjugateGradient(AtA, Atb, printLogs, outData);
	}
	else if(solver_method == Solver::LeastSquaresConjugateGradient)
	{
		x = BiharmonicSolver::computeLeastSquares(A, b, printLogs, outData);
	}
	double relative_error = (A*x - b).norm() / b.norm();
	if(printLogs) cout << "[DATA] The relative error is: " << relative_error << endl;
	
	std::cout << "[DATA] Field Height: " << field.height() << " Field Width: " << field.width() << std::endl;
	std::cout << "[DATA] Solver Rows : " << x.rows()  << std::endl;

	//NO ES RECORRE PER FILES!

	if(fullGridSubdivision)
	{
		for(size_t pos = 0; pos < x.size(); ++pos)
		{
			int cornerIndex = unknownToCornerIndex[pos];		
			glm::vec2 corner = indexToCorner[cornerIndex];

			field(corner.x, corner.y) = x(pos);
		}
	}
	
	//cas que no sigui full sibdivision:
	else
	{
		fillField(root, field, x);
	}
	return outData;

}

//- Recorrer totes les fulles del quadtree:
		//- Saps 4 corners:
			// - obtens valors dels index de cada corner del vector x
			// - per tota la fulla interpolar els 4 valors del pas anterior
				// - Fer un cas general, no intentar fer per casos en particular
void Quadtree::fillField(QuadtreeNode *node, ScalarField &field, const Eigen::VectorXd &x)
{
	if(node->isLeaf())
	{
		int cornerIndexBottomLeft  = cornerToInteger(node->minNodeCorner);
		int cornerIndexTopRight    = cornerToInteger(node->maxNodeCorner);		
		int cornerIndexTopLeft     = cornerToInteger(glm::vec2(node->minNodeCorner.x, node->maxNodeCorner.y));
		int cornerIndexBottomRight = cornerToInteger(glm::vec2(node->maxNodeCorner.x, node->minNodeCorner.y));

		int unknownBottomLeft  = cornerIndexToUnknown[cornerIndexBottomLeft];
		int unknownTopRight    = cornerIndexToUnknown[cornerIndexTopRight];
		int unknownTopLeft     = cornerIndexToUnknown[cornerIndexTopLeft];
		int unknownBottomRight = cornerIndexToUnknown[cornerIndexBottomRight];

		float value_bottom_left = x(unknownBottomLeft);
		float value_bottom_right = x(unknownBottomRight);
		float value_top_left = x(unknownTopLeft);
		float value_top_right = x(unknownTopRight);

		glm::ivec2 texel;
		glm::vec2 node_coords;
		for(texel.y = node->minNodeCorner.y; texel.y <= node->maxNodeCorner.y; texel.y++)
		{
			for(texel.x = node->minNodeCorner.x; texel.x <= node->maxNodeCorner.x; texel.x++)
			{
				node_coords.x = (float)(texel.x - node->minNodeCorner.x)/(node->maxNodeCorner.x - node->minNodeCorner.x);
				node_coords.y = (float)(texel.y - node->minNodeCorner.y)/(node->maxNodeCorner.y - node->minNodeCorner.y);

				field(texel.x, texel.y) = (1-node_coords.x)*(1-node_coords.y)*value_bottom_left + 
						node_coords.x*(1-node_coords.y)*value_bottom_right +
						(1 - node_coords.x)*node_coords.y*value_top_left +
						node_coords.x*node_coords.y*value_top_right;												
			}
		}
		//Doble bucle recorrer totes les posicions iterant desde min.x a max.x i min.y a max.y
		//Ens dius la posicio on es troba
	}
	else
	{
		for(unsigned int i=0; i<4; i++)
			fillField(node->children[i], field, x);
	}
}

void Quadtree::draw(Image &image)
{
	if(image.width() == 0)
		return;
	image.fill(glm::vec3(0.0f, 0.0f, 0.0f));
	root->draw(image);
	
	for(map<int, int>::iterator it=cornerIndexToUnknown.begin(); it!=cornerIndexToUnknown.end(); it++)
	{
		int cornerIndex = it->first;
		float y = float(cornerIndex / ((1 << nLevels) + 1)) / (1 << nLevels);
		float x = float(cornerIndex % ((1 << nLevels) + 1)) / (1 << nLevels);
		image.drawFilledCircle(int(image.width() * x), int(image.height() * y), 4, glm::vec3(0.75f, 0.75f, 0.75f));
	}
	
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

//NOU CORNER TO INTEGER
unsigned int Quadtree::cornerToInteger(const glm::ivec2 &C) const
{
	return C.y * ((1 << nLevels) + 1) + C.x;
}

//NOU INTEGER TO CORNER


unsigned int Quadtree::nodeToInteger(const glm::ivec2 &P) const
{
	return (((1 << nLevels)+1) * P.y) + P.x;
}


void Quadtree::addPointEquation(unsigned int eqIndex, const glm::vec2 &P, vector<Eigen::Triplet<double>> &triplets, vector<float> &bCoeffs, const float &value)
{
	unsigned int i, j;
	QuadtreeNode *node;
	float x, y;
	
	node = pointToCell(P);

	unsigned int fieldSize = std::pow(2, nLevels)+1;

	i = (unsigned int)floor(P.x * (fieldSize - 1));
	i = std::max(0u, glm::min(i, fieldSize-2));
	j = (unsigned int)floor(P.y * (fieldSize - 1));
	j = glm::max(0u, glm::min(j, fieldSize - 2));
	
	x = P.x * (fieldSize - 1) - float(i);
	y = P.y * (fieldSize - 1) - float(j);
	
	//std::cout << "QUAD FIELD EQ I: " << node->minNodeCorner.x << " J: " << node->minNodeCorner.y << std::endl;
	//std::cout << "QUAD POINT EQ I: " << i << " J: " << j << std::endl;
	//std::cout << "QUAD POINT EQ X: " << x << " Y: " << y << std::endl << std::endl;

	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[3], pW*x*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[2], pW*(1.0f-x)*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[1], pW*x*(1.0f-y)));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[0], pW*(1.0f-x)*(1.0f-y)));
	
	bCoeffs.push_back(value);
}

void Quadtree::addGradientEquations(unsigned int eqIndex, const glm::vec2 &P, const glm::vec2 &N, vector<Eigen::Triplet<double>> &triplets, vector<float> &bCoeffs)
{
	unsigned int i, j;
	QuadtreeNode *node;
	float x, y;
	
	node = pointToCell(P);

	unsigned int fieldSize = std::pow(2, nLevels)+1;

	i = (unsigned int)floor(P.x * (fieldSize - 1));
	i = std::max(0u, glm::min(i, fieldSize-2));
	j = (unsigned int)floor(P.y * (fieldSize - 1));
	j = glm::max(0u, glm::min(j, fieldSize - 2));
	
	x = P.x * (fieldSize - 1) - float(i);
	y = P.y * (fieldSize - 1) - float(j);
	
	//std::cout << "QUAD FIELD    EQ I: " << node->minNodeCorner.x << " J: " << node->minNodeCorner.y << std::endl;
	//std::cout << "QUAD GRADIENT EQ I: " << i << " J: " << j << std::endl;
	//std::cout << "QUAD GRADIENT EQ X: " << x << " Y: " << y << std::endl << std::endl;

	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[3], gW*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[2], -gW*y));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[1], gW*(1.0f-y)));
	triplets.push_back(Eigen::Triplet<double>(eqIndex, node->cornerUnknowns[0], -gW*(1.0f-y)));
	
	bCoeffs.push_back( gW* 64 / (fieldSize-1)*N.x);

	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, node->cornerUnknowns[3], gW*x));
	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, node->cornerUnknowns[2], gW*(1.0f-x)));
	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, node->cornerUnknowns[1], -gW*x));
	triplets.push_back(Eigen::Triplet<double>(eqIndex+1, node->cornerUnknowns[0], -gW*(1.0f-x)));
	
	bCoeffs.push_back(gW* 64 / (fieldSize-1) * N.y);
}
void Quadtree::addSamplingEquations(const data_representation::Mesh &cloud, unsigned int eqIndex, const glm::vec2 &P, const glm::vec2 &N, vector<Eigen::Triplet<double>> &triplets, vector<float> &bCoeffs)
{
	for(unsigned int i=0; i<cloud.vertices_.size(); i++)
	{
		QuadtreeNode *node = pointToCell(P);
		uint nodeWidth  = node->maxCoords.x - node->minCoords.x;
		uint nodeHeight = node->maxCoords.y - node->minCoords.y;

		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()) + (1.0f /  (1 << nLevels)) * glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), triplets, bCoeffs, 1.0f);
		eqIndex++;
		addPointEquation(eqIndex, glm::vec2(cloud.vertices_[i].x(), cloud.vertices_[i].y()) - (1.0f / (1 << nLevels)) * glm::vec2(cloud.normals_[i].x(), cloud.normals_[i].y()), triplets, bCoeffs, -1.0f);
		eqIndex++;
	}
}

int Quadtree::addHorizontalBoundarySmoothnessEquation(unsigned int eqIndex, const glm::ivec2 &P, vector<Eigen::Triplet<double>> &triplets, vector<float> &bCoeffs, const float &value)
{	
	glm::ivec2 left = P;
	--left.x;

	glm::ivec2 right = P;
	++right.x;

	//std::cout << "LOOKING FOR POINT P(" << P.x << ", " << P.y << ")-> " << "NEIGHBORS RIGHT: (" << right.x << ", " << P.y << ")  LEFT: (" << left.x << ", " << P.y << ")... ";

	if(right.x <= pow(2, nLevels) && left.x >= 0)
	{		
		std::map<int, float> myMap;

		// First add the point itself
		//
		accumulateToMap(myMap, cornerIndexToUnknown[cornerToInteger(P)], -2.0f);

		// Check if exists right & left neighbor nodes
		//
		std::map<int, int>::iterator it_right = cornerIndexToUnknown.find(cornerToInteger(right));
		std::map<int, int>::iterator it_left = cornerIndexToUnknown.find(cornerToInteger(left));

		// Case right
		//
		if(it_right != cornerIndexToUnknown.end())
		{
			accumulateToMap(myMap, cornerIndexToUnknown[cornerToInteger(right)], 1.0f);
		}

		// Case no right
		//
		else
		{
			// Get the point cell
			glm::vec2 normalized_right = glm::vec2(right.x/(pow(2, nLevels)+1), right.y/(pow(2, nLevels)+1));
			//std::cout << "NO RIGHT FOUND (" << normalized_right.x << ", " << normalized_right.y << ")" << std::endl;
			QuadtreeNode *point_cell = pointToCell(normalized_right);
			glm::vec2 top_left     = point_cell->minCoords;
			glm::vec2 bottom_right = point_cell->maxCoords;
			
			// v(x0, y0) = top-left coner (min-max)
			// v(x1, y1) = bottom-right corner (max-min)
			// v(x', y') = ((x-x0)/(x1,x0), (y-y0)/(y1-y0))
			
			float x = (normalized_right.x-top_left.x)/(bottom_right.x - top_left.x);
			float y = (normalized_right.y-top_left.y)/(bottom_right.y - top_left.y);

			//std::cout << "NO RIGHT FOUND (" << normalized_right.x << ", " << normalized_right.y << ").. NODE TOP LEFT (" << top_left.x << ", " << top_left.y << ") NODE BOTTOM RIGHT (" << bottom_right.x << ", " << bottom_right.y << ").. X = " <<  x << " Y = " << y << std::endl;

			// v(x, y) = (1-y')*(1-x')*v00 + (1-y')*x'*v10 + y'*(1-x')*v01 + y'*x'*v11			
			accumulateToMap(myMap, point_cell->cornerUnknowns[0], (1-y)*(1-x)); //top-left v00
			accumulateToMap(myMap, point_cell->cornerUnknowns[1], (1-y)*x);	 	//top-right v10
			accumulateToMap(myMap, point_cell->cornerUnknowns[2], y*(1-x)); 	//bottom-left v10
			accumulateToMap(myMap, point_cell->cornerUnknowns[3], y*x); 		//bottom-right v11
		}
		
		// Case left
		//
		if(it_left != cornerIndexToUnknown.end())
		{
			accumulateToMap(myMap, cornerIndexToUnknown[cornerToInteger(left)], 1.0f);
		}

		// Case no left
		//
		else
		{
			// Get the point cell
			glm::vec2 normalized_left = glm::vec2(left.x/(pow(2, nLevels)+1), left.y/(pow(2, nLevels)+1));
			QuadtreeNode *point_cell = pointToCell(normalized_left);
			glm::vec2 top_left     = point_cell->minCoords;
			glm::vec2 bottom_right = point_cell->maxCoords;
			
			// v(x0, y0) = top-left coner (min-max)
			// v(x1, y1) = bottom-right corner (max-min)
			// v(x', y') = ((x-x0)/(x1,x0), (y-y0)/(y1-y0))

			float x = (normalized_left.x-top_left.x)/(bottom_right.x - top_left.x);
			float y = (normalized_left.y-top_left.y)/(bottom_right.y - top_left.y);

			// v(x, y) = (1-y')*(1-x')*v00 + (1-y')*x'*v10 + y'*(1-x')*v01 + y'*x'*v11			
			accumulateToMap(myMap, point_cell->cornerUnknowns[0], (1-y)*(1-x)); // v00
			accumulateToMap(myMap, point_cell->cornerUnknowns[1], (1-y)*x);	 	// v10
			accumulateToMap(myMap, point_cell->cornerUnknowns[2], y*(1-x)); 	// v10
			accumulateToMap(myMap, point_cell->cornerUnknowns[3], y*x); 		// v11
		}

		// Read all map and add Eq
		//
		for(std::map<int, float>::iterator it = myMap.begin(); it != myMap.end(); it++)
		{
			triplets.push_back(Eigen::Triplet<double>(eqIndex, it->first, sW * it->second));
		}

		bCoeffs.push_back(0.0f);
		return 1;
	}
	else 
	{
		//std::cout << "\tONE OF THE NEIGHBORS IS OUT OF BOUNDS" << std::endl;
		return 0;
	}
	std::cout << endl;
}

int Quadtree::addVerticalBoundarySmoothnessEquation(unsigned int eqIndex, const glm::ivec2 &P, vector<Eigen::Triplet<double>> &triplets, vector<float> &bCoeffs, const float &value)
{
	
	glm::ivec2 bottom = P;
	--bottom.y;

	glm::ivec2 top = P;
	++top.y;

	//std::cout << "LOOKING FOR POINT P(" << P.x << ", " << P.y << ")-> " << "NEIGHBORS TOP: (" << P.x << ", " << top.y << ")  BOTTOM: (" << P.x << ", " << bottom.y << ")... ";

	if(top.y <= pow(2, nLevels) && bottom.y >= 0)
	{
		std::map<int, float> myMap;

		// First add the point itself
		//
		accumulateToMap(myMap, cornerIndexToUnknown[cornerToInteger(P)], -2.0f);

		// Check if exists right & left neighbor nodes
		//
		std::map<int, int>::iterator it_top = cornerIndexToUnknown.find(cornerToInteger(top));
		std::map<int, int>::iterator it_bottom = cornerIndexToUnknown.find(cornerToInteger(bottom));

		// Case right
		//
		if(it_top != cornerIndexToUnknown.end())
		{
			accumulateToMap(myMap, cornerIndexToUnknown[cornerToInteger(top)], 1.0f);
		}

		// Case no right
		//
		else
		{
			// Get the point cell
			glm::vec2 normalized_top = glm::vec2(top.x/(pow(2, nLevels)+1), top.y/(pow(2, nLevels)+1));
			QuadtreeNode *point_cell = pointToCell(normalized_top);
			glm::vec2 top_left     = point_cell->minCoords;
			glm::vec2 bottom_right = point_cell->maxCoords;
			
			// v(x0, y0) = top-left coner (min-max)
			// v(x1, y1) = bottom-right corner (max-min)
			// v(x', y') = ((x-x0)/(x1,x0), (y-y0)/(y1-y0))

			float x = (normalized_top.x-top_left.x)/(bottom_right.x - top_left.x);
			float y = (normalized_top.y-top_left.y)/(bottom_right.y - top_left.y);

			// v(x, y) = (1-y')*(1-x')*v00 + (1-y')*x'*v10 + y'*(1-x')*v01 + y'*x'*v11			
			accumulateToMap(myMap, point_cell->cornerUnknowns[0], (1-y)*(1-x)); //top-left v00
			accumulateToMap(myMap, point_cell->cornerUnknowns[1], (1-y)*x);	 	//top-right v10
			accumulateToMap(myMap, point_cell->cornerUnknowns[2], y*(1-x)); 	//bottom-left v10
			accumulateToMap(myMap, point_cell->cornerUnknowns[3], y*x); 		//bottom-right v11
		}
		
		// Case bottom
		//
		if(it_bottom != cornerIndexToUnknown.end())
		{
			accumulateToMap(myMap, cornerIndexToUnknown[cornerToInteger(bottom)], 1.0f);
		}

		// Case no bottom
		//
		else
		{
			// Get the point cell
			glm::vec2 normalized_bottom = glm::vec2(bottom.x/(pow(2, nLevels)+1), bottom.y/(pow(2, nLevels)+1));
			QuadtreeNode *point_cell = pointToCell(normalized_bottom);
			glm::vec2 top_left     = point_cell->minCoords;
			glm::vec2 bottom_right = point_cell->maxCoords;
			
			// v(x0, y0) = top-left coner (min-max)
			// v(x1, y1) = bottom-right corner (max-min)
			// v(x', y') = ((x-x0)/(x1,x0), (y-y0)/(y1-y0))
			float x = (normalized_bottom.x-top_left.x)/(bottom_right.x - top_left.x);
			float y = (normalized_bottom.y-top_left.y)/(bottom_right.y - top_left.y);

			// v(x, y) = (1-y')*(1-x')*v00 + (1-y')*x'*v10 + y'*(1-x')*v01 + y'*x'*v11			
			accumulateToMap(myMap, point_cell->cornerUnknowns[0], (1-y)*(1-x)); //top-left v00
			accumulateToMap(myMap, point_cell->cornerUnknowns[1], (1-y)*x);	 	//top-right v10
			accumulateToMap(myMap, point_cell->cornerUnknowns[2], y*(1-x)); 	//bottom-left v10
			accumulateToMap(myMap, point_cell->cornerUnknowns[3], y*x); 		//bottom-right v11
		}

		// Read all map and add Eq
		//
		for(std::map<int, float>::iterator it = myMap.begin(); it != myMap.end(); it++)
		{
			triplets.push_back(Eigen::Triplet<double>(eqIndex, it->first, sW * it->second));
		}

		bCoeffs.push_back(0.0f);
		return 1;
	}
	else 
	{
		//std::cout << "\tONE OF THE NEIGHBORS IS OUT OF BOUNDS" << std::endl;
		return 0;
	}

	std::cout << endl;
}

void Quadtree::accumulateToMap(std::map<int, float> &map, const int &key, const float &value)
{
	if(map.find(key) == map.end())
	{
		map[key] = value;
	}
	else
	{
		map[key] += value;
	}
}
