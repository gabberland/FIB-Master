#ifndef POINTS_GENERATOR_H_
#define POINTS_GENERATOR_H_

#include "mesh.h"

namespace data_representation 
{

class pointsGenerator 
{

public:

    bool generateRandomPoints(Mesh *mesh);

    /**
     * @brief norm Computes the norm between the two input values
     * @param x The first value.
     * @param y The second value.
     * @return The norm between the two values.
     */
    float norm(const float &x, const float &y);

        /**
     * @brief normal Computes the normalization between the two input values
     * @param x The first value.
     * @param y The second value.
     * @return The normal between the two values.
     */
    Eigen::Vector2f normalize(const float &x, const float &y);

    /**
     * @brief distance Computes the distance between the two input values
     * @param x The first value.
     * @param y The second value.
     * @return The distance between the two values.
     */
    float distance(const float &x, const float &y);

    /**
     * @brief gradient Computes the gradient between the two input values
     * @param x The first value.
     * @param y The second value.
     * @return The gradient between the two values.
     */
    Eigen::Vector2f gradient(const float &x, const float &y);

    /**
     * @brief project Computes the projection between the two input values
     * @param x The first value.
     * @param y The second value.
     * @return The projection between the two values.
     */
    Eigen::Vector2f project(const float &x, const float &y);    

    /**
     * @brief randomDisplacement Computes a random displacement 
     * @return a random displacement.
     */
    Eigen::Vector2f randomDisplacement();

};

}  // namespace data_representation

#endif  // POINTS_GENERATOR_H_