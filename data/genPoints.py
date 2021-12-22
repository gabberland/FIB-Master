import math, random
from PIL import Image, ImageDraw 


def norm(x, y):
	return math.sqrt(x*x + y*y)

def normalize(x, y):
	n = norm(x, y)
	return x / n, y / n

#def distance(x, y):
#	return norm(x - 0.5, y - 0.5) - 0.25

def distance(x, y):
	cx, cy = max(0.25, min(x, 0.75)), 0.5
	return norm(cx-x, cy-y) - 0.125

def gradient(x, y):
	prec = 0.001
	gX = (distance(x + prec, y) - distance(x - prec, y)) / prec / 2
	gY = (distance(x, y + prec) - distance(x, y - prec)) / prec / 2
	return gX, gY

def project(x, y):
	iter = 0
	cX = x
	cY = y
	d = distance(cX, cY)
	while abs(d) > 0.001 and iter < 10:
		gX, gY = gradient(cX, cY)
		cX -= d * gX
		cY -= d * gY
		d = distance(cX, cY)
		iter += 1
	return cX, cY

def randomDisplacement():
	alpha = 2 * math.pi * random.random()
	c = random.random()
	d = c ** 2
	return d * math.cos(alpha), d * math.sin(alpha)

# print("Distance: " + str(distance(0.5, 0)))
# print("Gradient: " + str(gradient(0, 0.5)))
# gx, gy = gradient(0.8, 0.8)
# print("Gradient: " + str(normalize(gx, gy)))
# print("Project: " + str(project(0.2, 0.3)))

noiseSize = 0.05
w, h = 800, 800
img = Image.new("RGB", (w, h)) 
iDraw = ImageDraw.Draw(img) 
npoints = 1000
f = open("points.txt", "w")
f.write(str(npoints) + "\n")
for i in range(npoints):
	dx, dy = randomDisplacement()
	dx = noiseSize * dx
	dy = noiseSize * dy
	x = random.random()
	y = random.random()
	x, y = project(x, y)
	x = x + dx
	y = y + dy
	nX, nY = gradient(x, y)
	nX, nY = normalize(nX, nY)
	# print("P = (" + str(x) + ", " + str(y) + ")")
	# print("N = (" + str(nX) + ", " + str(nY) + ")")
	# print("")
	#if abs(nX) < 0.3:
	f.write(str(x) + " " + str(y) + " " + str(nX) + " " + str(nY) + "\n")
	iDraw.ellipse([w*x-1, h*y-1, w*x+1, h*y+1], "white", "white")
	iDraw.line([w*x, h*y, w*(x + 0.03*nX), h*(y + 0.03*nY)], "yellow")
	
f.close()
img.show()
img.save("points.png", "PNG")





