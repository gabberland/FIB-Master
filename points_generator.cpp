#include "points_generator.h"

#include <math.h> 
#include <algorithm>
#include <stdlib.h>
#include <iostream>

namespace data_representation 
{
    #define RANDOM_BETWEEN_0_1 static_cast <float> (rand()) / static_cast <float> (RAND_MAX)
    float pointsGenerator::norm(const float &x, const float &y)
    {
        return sqrt(x*x + y*y);
    }

    Eigen::Vector2f pointsGenerator::normalize(const float &x, const float &y)
    {
        float n = norm(x, y);
	    return Eigen::Vector2f(x / n, y / n);
    }

    float pointsGenerator::distance(const float &x, const float &y)
    {
        float cx = std::fmax(0.25, std::fmin(x, 0.75));
        float cy = 0.5;
        return norm(cx-x, cy-y) - 0.125;
    }

    Eigen::Vector2f pointsGenerator::gradient(const float &x, const float &y)
    {
        float prec = 0.001;
        float gX = (distance(x + prec, y) - distance(x - prec, y)) / prec / 2;
        float gY = (distance(x, y + prec) - distance(x, y - prec)) / prec / 2;
        return Eigen::Vector2f(gX, gY);
    }

    Eigen::Vector2f pointsGenerator::project(const float &x, const float &y)
    {
        size_t iter = 0;
        float cX = x;
        float cY = y;
        float d = distance(cX,cY);
        while(abs(d) > 0.001 && iter < 10)
        {
            Eigen::Vector2f g = gradient(cX,cY);
            cX -= d * g.x();
            cY -= d * g.y();
            d = distance(cX,cY);
            iter += 1;
        }
        return Eigen::Vector2f(cX,cY);
    }

    Eigen::Vector2f pointsGenerator::randomDisplacement()
    {
        float alpha = 2 * M_PI * RANDOM_BETWEEN_0_1;
        float c = RANDOM_BETWEEN_0_1;
        float d = c * c;
        return Eigen::Vector2f(d * cos(alpha), d * sin(alpha));
    }

    bool pointsGenerator::generateRandomPoints(Mesh *mesh)
    {
        mesh->vertices_.clear();
        mesh->normals_.clear();

        float noiseSize = 0.05;
        int nPoints = 100;

        for( size_t i = 0; i < nPoints; ++i)
        {
            Eigen::Vector2f d = randomDisplacement();
            float dx = noiseSize * d.x();
            float dy = noiseSize * d.y();
            float x = RANDOM_BETWEEN_0_1;
            float y = RANDOM_BETWEEN_0_1;
            Eigen::Vector2f p = project(x, y);
            x = p.x() + dx;
            y = p.y() + dy;
            Eigen::Vector2f g = gradient(x,y);
            Eigen::Vector2f n = normalize(g.x(),g.y());
            
            mesh->vertices_.push_back(data_representation::Coordinate2d(x,y));
            mesh->normals_.push_back(data_representation::Coordinate2d(n.x(), n.y()));
        }

        return true;
    }
} //namespace data_representation