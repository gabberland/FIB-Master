#ifndef _SCALAR_FIELD_INCLUDE
#define _SCALAR_FIELD_INCLUDE


#include <vector>
#include "Image.h"


using namespace std;


class ScalarField
{
public:
	void init(unsigned int width, unsigned int height);
	
	unsigned int width() const;
	unsigned int height() const;
	float &operator()(unsigned int xIndex, unsigned int yIndex);
	const float &operator()(unsigned int xIndex, unsigned int yIndex) const;
	
	void sampleFill();
	
	Image *toImage(float maxValue = 1.f, float zeroThickness = 0.005f) const;

private:
	unsigned int w, h;
	vector<float> values;

};


#endif 


