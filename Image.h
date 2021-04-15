#ifndef _IMAGE_INCLUDE
#define _IMAGE_INCLUDE


#include <glm/glm.hpp>
#include <vector>
#include <string>


using namespace std;


class Image
{
public:
	Image();

	void init(unsigned int width, unsigned int height);
	
	unsigned int width() const;
	unsigned int height() const;
	glm::vec3 &operator()(unsigned int x, unsigned int y);
	const glm::vec3 &operator()(unsigned int x, unsigned int y) const;
	
	void fill(const glm::vec3 &color);
	void drawHorizontalLine(unsigned int x0, unsigned int x1, unsigned int y, const glm::vec3 &color);
	void drawVerticalLine(unsigned int x, unsigned int y0, unsigned int y1, const glm::vec3 &color);
	void drawRectangle(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, const glm::vec3 &color);
	void drawFilledCircle(unsigned int x, unsigned int y, float radius, const glm::vec3 &color);
	
	bool savePNG(const string &filename) const;
	
private:
	unsigned int w, h;
	vector<glm::vec3> pixels;

};


#endif


