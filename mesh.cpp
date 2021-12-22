#include "mesh.h"

#include <algorithm>
#include <limits>
#include <random>
#include <GL/glew.h>
#include <GL/glut.h>

namespace data_representation {

Mesh::Mesh() { Clear(); }

void Mesh::Clear() {
  vertices_.clear();
  normals_.clear();

  min_ = Eigen::Vector2f(std::numeric_limits<float>::max(),
                         std::numeric_limits<float>::max());
  max_ = Eigen::Vector2f(std::numeric_limits<float>::lowest(),
                         std::numeric_limits<float>::lowest());
}

glm::vec2 Mesh::point(const int & i)
{
  return glm::vec2(vertices_[i].x(), vertices_[i].y());
}

glm::vec2 Mesh::normal(const int & i)
{
  return glm::vec2(normals_[i].x(), normals_[i].y());
}

void Mesh::writeTxtFile(const std::string &filename)
{
  // Create and open a text file
  std::ofstream MyFile(filename);

  // Write to the file
  MyFile << vertices_.size() << std::endl;
  
  for(size_t i = 0; i < vertices_.size(); ++i)
  {
    MyFile << vertices_[i].x() << " " << vertices_[i].y() << " " << normals_[i].x() << " " << normals_[i].y() << std::endl;
  }

  // Close the file
  MyFile.close();
}

void Mesh::normalizePointValues()
{
/*
  //Compute Bounding Box Center
  float bboxX = (max_.x() - min_.x());
  float bboxY = (max_.y() - min_.y());
  
  //Compute Longest edge in order to make it  squared
  float longestEdge = std::max(max_.x() - min_.x(), max_.y() - min_.y() );

  float max = std::max(max_[0], max_[1]);
  float offset = std::abs(max_[0]-max_[1]);
*/
  for(size_t i = 0; i < vertices_.size(); ++i)
  {
    vertices_[i].x() = (float)(vertices_[i].x() * 0.8f / max_[0] + 0.1);
    vertices_[i].y() = (float)(vertices_[i].y() * 0.8f / max_[1] + 0.1);
  }
}

void Mesh::computeBoundingBox() 
{
  const size_t kVertices = vertices_.size();
  for (size_t i = 0; i < kVertices; ++i) {
    min_[0] = std::min(min_[0], vertices_[i][0]);
    min_[1] = std::min(min_[1], vertices_[i][1]);

    max_[0] = std::max(max_[0], vertices_[i][0]);
    max_[1] = std::max(max_[1], vertices_[i][1]);
  }
}

void Mesh::addNoise(const double &mean, const double &standardDev)
{
      // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device rd; 

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(rd()); 
    
    int i;
    double sample;
    for(i = 0; i < vertices_.size(); ++i)
    {
        // instance of class std::normal_distribution with specific mean and stddev
        std::normal_distribution<double> d(mean, standardDev); 

        // get random number with normal distribution using gen as random source
        sample = d(gen); 
        Coordinate2d displ = normals_[i] * sample;
        vertices_[i] += displ;
    }
}


}  // namespace data_representation
