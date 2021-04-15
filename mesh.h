#ifndef MESH_H_
#define MESH_H_

#include <iostream>
#include <fstream>
#include <vector>

#include <eigen3/Eigen/Geometry>
#include <glm/glm.hpp>

namespace data_representation {

typedef Eigen::Vector2f Coordinate2d;

class Mesh {
 public:
  /**
   * @brief Mesh Constructor of the class. Calls clear.
   */
  Mesh();

  /**
   * @brief ~Mesh Destructor of the class.
   */
  ~Mesh() {}

  /**
   * @brief Clear Empties the data arrays and resets the bounding box vertices.
   */
  void Clear();

 public:

  std::vector<Coordinate2d> vertices_;
  std::vector<Coordinate2d> normals_;

  glm::vec2 point(const int & i);

  glm::vec2 normal(const int & i);

  void writeTxtFile(const std::string &filename);

  void normalizePointValues();

  /**
   * @brief min The minimum point of the bounding box.
   */
  Eigen::Vector2f min_;

  /**
   * @brief max The maximum point of the bounding box.
   */
  Eigen::Vector2f max_;

  private:
    /*  Render data  */
};

}  // namespace data_representation

#endif  //  MESH_H_