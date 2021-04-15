#include "mesh.h"

#include <algorithm>
#include <limits>
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
  float max = std::max(max_[0], max_[1]);

  for(size_t i = 0; i < vertices_.size(); ++i)
  {
    vertices_[i].x() = (float)(vertices_[i].x() * 0.8f / max + 0.1);
    vertices_[i].y() = (float)(vertices_[i].y() * 0.8f / max + 0.1);
  }
}
}  // namespace data_representation
