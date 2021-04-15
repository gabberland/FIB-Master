#include <GL/glut.h>
#include <iostream> 
#include <memory>
#include <unistd.h>

#include "mesh_io.h"
#include "points_generator.h"
#include "mesh.h"
#include "ScalarField.h"
#include "Image.h"
#include "BiharmonicSolver.h"

std::shared_ptr<data_representation::Mesh> mesh_;

static const int 	PIXEL_MARGINS = 1;
static const string GRADIENT = "gradient";
static const string SMOOTHNESS = "2D";
int					FIELD_RESOLUTION = 64;
float 				NORMAL_SIZE = .2;
string				RELATIVE_OUT_PATH_NAME = "out.png";
string 				RELATIVE_IN_PATH_NAME;
bool				MODEL_GENERATION_MODE = 0;
bool 				MODEL_GENERATED = 0;
Normal				NORMAL_ALGORITHM = Normal::sampling;
Smoothness			SMOOTHNESS_ALGORTITHM = Smoothness::singleDimension;

void printUsage()
{
	std::cout << "Usage: reconstruction.exe [ -f | -g ] [ relative path | model generation mode ] [ -n | -r | -a | -s ] [ normal size | resolution | normal algorithm | smoothing algorithm ]" << std::endl;
}
void readFlagArguments(const int &argc, char **argv)
{
	int index;
	int c;
	opterr = 0;

	if(argc == 1)
	{
		printUsage();
		abort();
	}
	while ((c = getopt (argc, argv, "ghf:n:r:a:s:")) != -1)
	{
		switch (c)
		{
			case 'g':
				MODEL_GENERATION_MODE = 1;
				break;
						
			case 'h':
				printUsage();
				break;
			
			case 'f':
				RELATIVE_IN_PATH_NAME = optarg;
				break;
			
			case 'n':
				NORMAL_SIZE = stof(optarg);
				break;
			case 'r':
				FIELD_RESOLUTION = stoi(optarg);
				break;

			case 'a':
				if (GRADIENT == optarg)
					NORMAL_ALGORITHM = Normal::gradient;
				break;

			case 's':
				if (SMOOTHNESS == optarg)
					SMOOTHNESS_ALGORTITHM = Smoothness::twoDimension;
				break;

			case '?':
				if (optopt == 'f')
					fprintf (stderr, "Option -f requires a filepath as argument.\n");
				
				else if (optopt == 'n')
					fprintf (stderr, "Option -n requires a float value as argument.\n");				
				
				else if (optopt == 'r')
					fprintf (stderr, "Option -r requires an integer value as argument.\n");	

				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				
				else
				{
					fprintf (stderr,
							"Unknown option character `\\x%x'.\n",
							optopt);
					printUsage();
				}
			default:
				abort ();
		}
	}

	for (index = optind; index < argc; index++)
	{
		printf ("Non-option argument %s\n", argv[index]);
	}
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0, 0.0, 0.0);

	glBegin(GL_POINTS);
	
	for(size_t i = 0; i < mesh_->vertices_.size();  ++i)
	{
			glVertex2f(mesh_->vertices_[i][0], mesh_->vertices_[i][1]);
	}
	glEnd();

	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	for(size_t i = 0; i < mesh_->normals_.size();  ++i)
	{
			glVertex2f(mesh_->vertices_[i][0], mesh_->vertices_[i][1]);

			glVertex2f(mesh_->vertices_[i][0] + (NORMAL_SIZE *mesh_->normals_[i][0]), mesh_->vertices_[i][1] + (NORMAL_SIZE * mesh_->normals_[i][1]));
	}
	glEnd();

	glFlush();

}

void myinit() {
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(5.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(mesh_->min_.x() - PIXEL_MARGINS, mesh_->max_.x() + PIXEL_MARGINS, mesh_->min_.y() - PIXEL_MARGINS, mesh_->max_.y() + PIXEL_MARGINS);
}

void OnReshape(int w, int h)
{
	GLfloat formato;

	if(h == 0) h = 1;
		
	glViewport(0, 0, w, h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	formato = (GLfloat)w / (GLfloat)h;
	if (w <= h) glOrtho (-100.0f, 100.0f, -100.0f / formato, 100.0f / formato, 1.0f, -1.0f);
	else glOrtho (-100.0f * formato, 100.0f * formato, -100.0f, 100.0f, 1.0f, -1.0f);
}

int main(int argc, char** argv) {

	readFlagArguments(argc, argv);

	// Read file file | create random datapoints
	if(MODEL_GENERATION_MODE && !RELATIVE_IN_PATH_NAME.empty())
	{
		std::cout << "[ERROR] Please use -f | -g one by one" << std::endl;
	}

	if(!RELATIVE_IN_PATH_NAME.empty())
	{
		std::string file = (std::string)argv[2];
		size_t pos = file.find_last_of(".");
		std::string type = file.substr(pos + 1);

		mesh_ = std::make_shared<data_representation::Mesh>();

		bool res = false;
		if (type.compare("txt") == 0) 
		{
			res = data_representation::ReadFromTXT(file, mesh_.get());
			MODEL_GENERATED = true;
		}

		else if (type.compare("svg") == 0)
		{
			res = data_representation::ReadFromSVG(file, mesh_.get());
			std::cout << "[DATA] New file .txt generated, please run the script again with the new file" ;
		}

		else
		{
			std::cout << "[ERROR] " << type << " is not a valid .txt file type." << std::endl;
			return 0;
		}
	}

	else if(MODEL_GENERATION_MODE)
	{
		mesh_ = std::make_shared<data_representation::Mesh>();
		data_representation::pointsGenerator pg;
		pg.generateRandomPoints(mesh_.get());
		MODEL_GENERATED = true;
	}

	if(MODEL_GENERATED)
	{
		ScalarField field;
		BiharmonicSolver solver;
		Image *img;

		//field.init(256, 256);
		field.init(FIELD_RESOLUTION, FIELD_RESOLUTION);

		solver.setWeights(1.0f, 1.0f, 1.0f);
		//solver.compute(*mesh_.get(), field);
		//solver.computeBilaplacian(*mesh_.get(), field);
		//solver.computeComponentWise(*mesh_.get(), field);
		//solver.computeNoGradient(*mesh_.get(), field);
		solver.computeWith(*mesh_.get(), field, NORMAL_ALGORITHM, SMOOTHNESS_ALGORTITHM);


		img = field.toImage(16.0f, 0.0f);
		if(!img->savePNG(RELATIVE_OUT_PATH_NAME))
		{
			cout << "[ERROR] Could not save file!" << endl;
			delete img;
			return -1;
		}
		delete img;

		// Initialize Viewer
		//
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(100, 100);
		glutInitWindowPosition(0, 0);
		glutCreateWindow("Point Cloud Viewer");


		// Display Function
		glutDisplayFunc(display);

		// set the function to handle changes in screen size
		//glutReshapeFunc(OnReshape);

		// Properties Init
		myinit();
		glutMainLoop();
	}
	return 0;
}



