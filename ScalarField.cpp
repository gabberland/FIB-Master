#include <cmath>
#include <glm/glm.hpp>
#include "ScalarField.h"


glm::vec3 palette[7] = {
	glm::vec3(0.0f, 0.0f, 1.0f), 
	glm::vec3(0.0f, 0.5f, 1.0f), 
	glm::vec3(0.0f, 1.0f, 1.0f), 
	glm::vec3(1.0f, 1.0f, 1.0f), 
	glm::vec3(1.0f, 1.0f, 0.0f), 
	glm::vec3(1.0f, 0.5f, 0.0f), 
	glm::vec3(1.0f, 0.0f, 0.0f)
};


void ScalarField::init(unsigned int width, unsigned int height)
{
	w = width;
	h = height;
	values.resize(w*h);
}
	
unsigned int ScalarField::width() const
{
	return w;
}
	
unsigned int ScalarField::height() const
{
	return h;
}
	
float &ScalarField::operator()(unsigned int xIndex, unsigned int yIndex)
{
	return values[yIndex * w + xIndex];
}
	
const float &ScalarField::operator()(unsigned int xIndex, unsigned int yIndex) const
{
	return values[yIndex * w + xIndex];
}

void ScalarField::sampleFill()
{
	glm::vec2 P;

	for(unsigned int y=0; y<h; y++)
		for(unsigned int x=0; x<w; x++)
		{
			P = glm::vec2(float(x) / w - 0.5f, float(y) / h - 0.5f);
			values[y * w + x] = glm::length(P) - 0.25f;
		}
}

Image *ScalarField::toImage(float maxValue, float zeroThickness) const
{
	Image *img;
	float scaledValue, lambda;
	unsigned int colorIndex;
	
	img = new Image();
	img->init(w, h);
	for(unsigned int y=0; y<h; y++)
		for(unsigned int x=0; x<w; x++)
		{
			//scaledValue = fabs(values[y * w + x]) / maxValue;
			scaledValue = values[y * w + x] / maxValue;
			scaledValue = (scaledValue + 1.0f) / 2.0f;
			scaledValue = max(0.0f, min(scaledValue, 1.0f));
			colorIndex = min(int(6.0f * scaledValue), 6);
			lambda = 6.0f * (scaledValue - colorIndex / 6.0f);
			(*img)(x, y) = glm::mix(palette[colorIndex], palette[colorIndex+1], lambda);
			if(zeroThickness != 0.0f && fabs(values[y * w + x]) < zeroThickness)
				(*img)(x, y) = glm::vec3(0.0f);
			/*
			if(values[y * w + x] > 0.0f)
				(*img)(x, y) = glm::vec3(1.0f, 0.0f, 0.0f);
			else
				(*img)(x, y) = glm::vec3(0.0f, 0.0f, 1.0f);
			*/
		}
	
	return img;
}


