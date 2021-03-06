#include "mesh_io.h"
#include "points_generator.h"

#include "pugixml.hpp"
#include "stdlib.h"
#include <assert.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <boost/tokenizer.hpp>

namespace data_representation {

void AddItems(float i1, float i2, size_t index, std::vector<Eigen::Vector2f> *vector) 
{
  (*vector)[index] = Eigen::Vector2f(i1, i2);
}

Coordinate2d readCoordinate(std::vector<std::string> array, const int &i)
{
  std::stringstream ss (array[i]);
  std::string token;
  
  std::getline(ss, token, ',');
  
  float vx = std::stof(token);
  std::getline(ss, token, ',');
  float vy = std::stof(token);
  return Coordinate2d(vx,vy);
}

Coordinate2d readHorizontalCoordinate(std::vector<std::string> array, const int &i)
{
  std::stringstream ss (array[i]);
  std::string token;
  
  std::getline(ss, token, ',');
  
  float vx = std::stof(token);
  return Coordinate2d(vx,0);
}

Coordinate2d readVerticalCoordinate(std::vector<std::string> array, const int &i)
{
  std::stringstream ss (array[i]);
  std::string token;
  
  std::getline(ss, token, ',');

  float vy = std::stof(token);
  return Coordinate2d(0,vy);
}

bool ReadTxtHeader(std::ifstream *fin, int &vertices) 
{

  char line[100];

  fin->getline(line, 100);
  vertices = atoi(&line[0]);
  //std::cout << "Loading point cloud mesh" << std::endl;
  //std::cout << "\tVertices = " << vertices << std::endl;

  return true;
}

void ReadTxtVertices(std::ifstream *fin, Mesh *mesh) 
{
  char line[100];
  const size_t kVertices = mesh->vertices_.size();
  float z = 0;
  for (size_t i = 0; i < kVertices; ++i) {    

    fin->getline(line, 100);

    std::string s = line;
    boost::char_separator<char> sep(" ");
    boost::tokenizer<boost::char_separator<char>> tokens(s, sep);
    std::vector<std::string> array;
    for (const auto& t : tokens)
    {
        array.push_back(t);
    }

    float vx, vy, nx, ny;

    vx = std::stof(array[0]);
    vy = std::stof(array[1]);
    AddItems(vx, vy, i , &(mesh->vertices_));

    nx = std::stof(array[2]);
    ny = std::stof(array[3]);
    AddItems(nx, ny, i, &(mesh->normals_));
  }
}

void ReadTxtVerticesWithDensity(std::ifstream *fin, Mesh *mesh, const float& density) 
{
  char line[100];
  const size_t kVertices = mesh->vertices_.size();
  mesh->vertices_.clear();
  mesh->normals_.clear();

  int spacing = 0;
  for (size_t i = 0; i < kVertices; ++i) {    

    fin->getline(line, 100);
    if(spacing <= 0)
    {
      int rand = std::rand() % 100;
      if(rand > density)
      {
        spacing = std::rand() % (int)(density*0.1) + (int)(density*0.05);
      }
    }
    if(spacing <= 0)
    {
      std::string s = line;
      boost::char_separator<char> sep(" ");
      boost::tokenizer<boost::char_separator<char>> tokens(s, sep);
      std::vector<std::string> array;
      for (const auto& t : tokens)
      {
          array.push_back(t);
      }

      float vx, vy, nx, ny;

      vx = std::stof(array[0]);
      vy = std::stof(array[1]);
      mesh->vertices_.push_back(Eigen::Vector2f(vx, vy));

      nx = std::stof(array[2]);
      ny = std::stof(array[3]);
      mesh->normals_.push_back(Eigen::Vector2f(nx, ny));
    }
    else
    {
      --spacing;
    }    
  }
}

void ComputeVertexNormals(Mesh *mesh) 
{
  mesh->normals_.resize(mesh->vertices_.size());
  const size_t kVertices = mesh->vertices_.size();
  std::vector<Coordinate2d>normals(kVertices);
  for (size_t i = 0; i < kVertices; ++i) {
    size_t j = (i+1) % kVertices;
    Coordinate2d edge_vec = mesh->vertices_[j] - mesh->vertices_[i];
    normals[i] = Coordinate2d(-edge_vec.y(), edge_vec.x()).normalized();
  }

  for (size_t i = 0; i < kVertices; ++i) {
    size_t j = (i+1) % kVertices;
    mesh->normals_[i] = (normals[i]+normals[j]).normalized();
  }
}

std::vector<Coordinate2d> GiveVertexNormals(const std::vector<Coordinate2d> &vertices) 
{
  const size_t kVertices = vertices.size();
  std::vector<Coordinate2d> n(kVertices);
  std::vector<Coordinate2d>normals(kVertices);
  for (size_t i = 0; i < kVertices; ++i) {
    size_t j = (i+1) % kVertices;
    Coordinate2d edge_vec = vertices[j] - vertices[i];
    normals[i] = Coordinate2d(-edge_vec.y(), edge_vec.x()).normalized();
  }

  for (size_t i = 0; i < kVertices; ++i) {
    size_t j = (i+1) % kVertices;
    n[i] = (normals[i]+normals[j]).normalized();
  }
  return n;
}

std::vector<data_representation::Coordinate2d> ComputeVertexNormals(const std::vector<data_representation::Coordinate2d> &v) 
{
  const size_t kVertices = v.size();
  std::vector<data_representation::Coordinate2d>normals(kVertices);
  for (size_t i = 0; i < kVertices; ++i) {
    size_t j = (i+1) % kVertices;
    data_representation::Coordinate2d edge_vec = v[j] - v[i];
    normals[i] = data_representation::Coordinate2d(-edge_vec.y(), edge_vec.x()).normalized();
  }
  std::vector<data_representation::Coordinate2d> r ;
  for (size_t i = 0; i < kVertices; ++i) {
    size_t j = (i+1) % kVertices;
    r.push_back((normals[i]+normals[j]).normalized());
  }
  return r;
}

bool ReadFromTXT(const std::string &filename, Mesh *mesh, const float &noise, const float &density) 
{
    std::ifstream fin;

    fin.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!fin.is_open() || !fin.good()) return false;

    int vertices = 0;
    if (!ReadTxtHeader(&fin, vertices)) {
      fin.close();
      return false;
    }

    mesh->vertices_.resize(static_cast<size_t>(vertices));
    mesh->normals_.resize(static_cast<size_t>(vertices));

    ReadTxtVerticesWithDensity(&fin, mesh, density);

    fin.close();
    mesh->addNoise(0, noise);
    mesh->computeBoundingBox();

    return true;
}

bool ReadFromSVG(const std::string &filepathName, Mesh *mesh) 
{

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filepathName.c_str());
    if (!result)
    {
      return false;
    }
      
    mesh->vertices_.clear();
    mesh->normals_.clear();

    for (pugi::xml_node tool: doc.child("svg").child("g").children("path"))
    {
        std::string attr = tool.attribute("d").value();
        std::stringstream ss(attr);
        std::vector<std::string> array;
        std::string value;
        while(std::getline(ss, value, ' '))
        {
          array.push_back(value);
        }

        std::vector<data_representation::Coordinate2d> vertices;
        std::vector<data_representation::Coordinate2d> normals;
        Coordinate2d currPathPoint = Coordinate2d(-1,-1);
        size_t i = 0;
        while (i < array.size() && array[i] != "z" && array[i] != "Z")
        {
            if (array[i] == "M" || array[i] == "L")
            {
              ++i;

              while (array[i] != "M" && array[i] != "L" && array[i] != "m" && array[i] != "l" && array[i] != "z" && array[i] != "Z" && array[i] != "v" && array[i] != "V" && array[i] != "h" && array[i] != "H" && i < array.size())
              {
                Coordinate2d coord;
                
                coord = readCoordinate(array, i);

                vertices.push_back(coord);

                currPathPoint = coord;

                ++i;

                std::cout << "ADDED " << coord << std::endl;
              }
              
            }

            else if (array[i] == "m" || "l")
            {
              ++i;
              
              while (array[i] != "M" && array[i] != "m" && array[i] != "L" && array[i] != "l" && array[i] != "z" && array[i] != "Z" && array[i] != "v" && array[i] != "V" && array[i] != "h" && array[i] != "H" && i < array.size())
              {
                Coordinate2d coord = readCoordinate(array, i);
                Coordinate2d newcoord = Coordinate2d(currPathPoint[0]+coord[0],currPathPoint[1]+coord[1]);
                vertices.push_back(newcoord);
                currPathPoint = newcoord;
                
                ++i;

                std::cout << "ADDED " << coord << std::endl;

              }             
            }

            else if(array[i] == "v" || array[i] != "V")
            {
              while (array[i] != "M" && array[i] != "m" && array[i] != "L" && array[i] != "l" && array[i] != "z" && array[i] != "Z" && array[i] != "v" && array[i] != "V" && array[i] != "h" && array[i] != "H" && i < array.size())
              {
                Coordinate2d coord = readVerticalCoordinate(array, i);
                Coordinate2d newcoord = Coordinate2d(currPathPoint[0]+coord[0],currPathPoint[1]+coord[1]);
                vertices.push_back(newcoord);
                currPathPoint = newcoord;
                
                ++i;

                std::cout << "ADDED " << coord << std::endl;

              }  
            }
            else if(array[i] == "h" || array[i] != "H")
            {
              while (array[i] != "M" && array[i] != "m" && array[i] != "L" && array[i] != "l" && array[i] != "z" && array[i] != "Z" && array[i] != "v" && array[i] != "V" && array[i] != "h" && array[i] != "H" && i < array.size())
              {
                Coordinate2d coord = readHorizontalCoordinate(array, i);
                Coordinate2d newcoord = Coordinate2d(currPathPoint[0]+coord[0],currPathPoint[1]+coord[1]);
                vertices.push_back(newcoord);
                currPathPoint = newcoord;
                
                ++i;

                std::cout << "ADDED " << coord << std::endl;
              }
            }  
            else if(array[i] == "z" || array[i] == "Z")
            {
              break;
            }
        }

      normals = GiveVertexNormals(vertices);
      mesh->vertices_.insert(mesh->vertices_.end(), vertices.begin(), vertices.end());
      mesh->normals_.insert(mesh->normals_.end(), normals.begin(), normals.end());
      mesh->computeBoundingBox(); 

    }
      mesh->normalizePointValues();
      mesh->computeBoundingBox();

    
    mesh->writeTxtFile("out.txt");
    return true;
}



}  // namespace data_representation
