#ifndef MESH_IO_H_
#define MESH_IO_H_

#include "mesh.h"

#include <string>

namespace data_representation {

/**
 * @brief ReadFromTXT Read the mesh stored in TXT format at the path filename
 * and stores the corresponding Mesh representation
 * @param filename The path to the TXT mesh.
 * @param mesh The resulting representation with computed per-vertex normals.
 * @return Whether it was able to read the file.
 */
bool ReadFromTXT(const std::string &filename, Mesh *mesh);

/**
 * @brief ReadFromSVG Read the mesh stored in SVG format at the path filename
 * and stores the corresponding Mesh representation
 * @param filename The path to the SVG mesh.
 * @param mesh The resulting representation with computed per-vertex normals.
 * @return Whether it was able to read the file.
 */
bool ReadFromSVG(const std::string &filepathName, Mesh *mesh);

/**
 * @brief WriteToTXT Stores the mesh representation in TXT format at the path
 * filename.
 * @param filename The path where the mesh will be stored.
 * @param mesh The mesh to be stored.
 * @return Whether it was able to store the file.
 */
bool WriteToTXT(const std::string &filename, const Mesh &mesh);

/**
 * @brief given a string array which represents two dimesional coordinates separed by comas,
 * provides the i-essim coordiate
 * @param array The array of string coordinates
 * @param i the i-essim value to return the coordinate
 * @return An Eigen two dimensional coordinate
 **/
Coordinate2d readCoordinate(std::vector<std::string> array, const int &i);

}  // namespace data_representation

#endif  // MESH_IO_H_