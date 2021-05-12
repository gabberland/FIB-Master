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
#include "Quadtree.h"

std::shared_ptr<data_representation::Mesh> mesh_;

static const int 	PIXEL_MARGINS = 1;
static const string GRADIENT = "gradient";
static const string SMOOTHNESS = "2D";

int					Field_Resolution_ = 64;
int					Multigrid_Iterations_ = 1;
float 				Normal_Size_ = .2;
string				Relative_Out_Path_Name_ = "out.png";
string 				Relative_In_Path_Name_;
bool				Model_Generation_Mode_ = 0;
bool 				Model_Generated_ = 0;
Normal				Normal_Algorithm = Normal::sampling;
Smoothness			Smoothness_Algorithm_ = Smoothness::singleDimension;
Solver				Solver_Method_ = Solver::BiCGSTAB;
bool				Print_Logs_ = 1;

void printUsage()
{
	std::cout << "Usage: " << std::endl << "  reconstruction.exe [ -f | -g ] [ relative path | model generation mode ]" << std::endl << 
	"Optional Parameters:" << std::endl << " \t[ -n | normal size ]" << std::endl << "\t[ -r | resolution ]" << std::endl << "\t[ -a | normal algorithm (gradient or sampling) ]" << std::endl << "\t[ -s | smoothing algorithm (1D or 2D) ]" << 
	std::endl << "\t[ -x | solver method (BiCGSTAB, ConjugateGradient or LeastSquaresConjugateGradient) ] " << std::endl << "\t[ -i | Multigrid Iterations ] " << std::endl << "\t[ -p | print logs, 0 or 1 ]" << std::endl;
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
	while ((c = getopt (argc, argv, "ghf:n:r:a:s:x:i:p:")) != -1)
	{
		switch (c)
		{
			case 'g':
				Model_Generation_Mode_ = 1;
				break;
						
			case 'h':
				printUsage();
				break;
			
			case 'f':
				Relative_In_Path_Name_ = optarg;
				break;
			
			case 'n':
				Normal_Size_ = stof(optarg);
				break;
			case 'r':
				Field_Resolution_ = stoi(optarg);
				break;

			case 'a':
				if (GRADIENT == optarg)
					Normal_Algorithm = Normal::gradient;
				break;

			case 's':
				if (SMOOTHNESS == optarg)
					Smoothness_Algorithm_ = Smoothness::twoDimension;
				break;

			case 'x':
				if(optarg == "LeastSquaresConjugateGradient")
					Solver_Method_ = Solver::LeastSquaresConjugateGradient;

				else if(optarg == "ConjugateGradient")
					Solver_Method_ = Solver::ConjugateGradient;
				break;

			case 'i':
				Multigrid_Iterations_ = stoi(optarg);
				break;

			case 'p':
				Print_Logs_ = stoi(optarg);
				break;

			case '?':
				if (optopt == 'f')
					fprintf (stderr, "Option -f requires a filepath as argument.\n");
				
				else if (optopt == 'n')
					fprintf (stderr, "Option -n requires a float value as argument.\n");				
				
				else if (optopt == 'r')
					fprintf (stderr, "Option -r requires an integer value as argument.\n");	

				else if (optopt == 'a')
					fprintf (stderr, "Option -a requires a normal type as argument.\n");	

				else if (optopt == 's')
					fprintf (stderr, "Option -s requires a smoothing type as argument.\n");	

				else if (optopt == 'x')
					fprintf (stderr, "Option -x requires a smoothing type as argument.\n");	

				else if (optopt == 'p')
					fprintf (stderr, "Option -p requires a bool value as argument.\n");	

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
			glVertex2f( 1 - mesh_->vertices_[i][0], 1 - mesh_->vertices_[i][1]);
	}
	glEnd();

	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	for(size_t i = 0; i < mesh_->normals_.size();  ++i)
	{
			glVertex2f(mesh_->vertices_[i][0], mesh_->vertices_[i][1]);

			glVertex2f(mesh_->vertices_[i][0] + (Normal_Size_ *mesh_->normals_[i][0]), mesh_->vertices_[i][1] + (Normal_Size_ * mesh_->normals_[i][1]));
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
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, 1.0 * (GLfloat) w / (GLfloat) h, 1.0, 30.0);
    glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv) {

	readFlagArguments(argc, argv);

	// Read file file | create random datapoints
	if(Model_Generation_Mode_ && !Relative_In_Path_Name_.empty())
	{
		std::cout << "[ERROR] Please use -f | -g one by one" << std::endl;
	}

	if(!Relative_In_Path_Name_.empty())
	{
		std::string file = (std::string)argv[2];
		size_t pos = file.find_last_of(".");
		std::string type = file.substr(pos + 1);

		mesh_ = std::make_shared<data_representation::Mesh>();

		bool res = false;
		if (type.compare("txt") == 0) 
		{
			res = data_representation::ReadFromTXT(file, mesh_.get());
			Model_Generated_ = true;
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

	else if(Model_Generation_Mode_)
	{
		mesh_ = std::make_shared<data_representation::Mesh>();
		data_representation::pointsGenerator pg;
		pg.generateRandomPoints(mesh_.get());
		Model_Generated_ = true;
	}

	if(Model_Generated_)
	{
		ScalarField field;
		BiharmonicSolver solver;
		Eigen::VectorXd guess;
		Image *img;

		//field.init(256, 256);
		field.init(Field_Resolution_, Field_Resolution_);

		solver.setWeights(1.0f, 1.0f, 1.0f);
		//solver.compute(*mesh_.get(), field);
		//solver.computeBilaplacian(*mesh_.get(), field);
		//solver.computeComponentWise(*mesh_.get(), field);
		//solver.computeNoGradient(*mesh_.get(), field);
		SolverData s = solver.computeWith(*mesh_.get(), field, Normal_Algorithm, Smoothness_Algorithm_, Solver_Method_, Print_Logs_);
		//SolverData s = solver.computeMultigrid(*mesh_.get(), field, Multigrid_Iterations_, Normal_Algorithm, Smoothness_Algorithm_, Solver_Method_, Print_Logs_);
		std::cout.precision(5);
		std::cout << Relative_In_Path_Name_ << ", " << Field_Resolution_ << ", " << NORMAL_STRING[Normal_Algorithm] << ", " << SMOOTHNESS_STRING[Smoothness_Algorithm_] << ", " << SOLVER_STRING[Solver_Method_] << ", " << RESULT_STRING[s.isSolved] << ", " << s.systemBuildTime << ", " <<s.matrixBuildTime << ", " <<s.solverResolutionTime << ", " << s.iterations << ", " << s.error << std::endl;

		img = field.toImage(16.0f, 0.0f);
		if(!img->savePNG(Relative_Out_Path_Name_))
		{
			cout << "[ERROR] Could not save file!" << endl;
			delete img;
			return -1;
		}
		delete img;
		
		Quadtree qtree;
		qtree.compute(*mesh_.get(), 8, field);
		Image qtreeImg;
		qtreeImg.init(1024, 1024);
		qtree.draw(qtreeImg);
		qtreeImg.savePNG("quadtree.png");
		
		// Initialize Viewer
		//
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(100, 100);
		glutInitWindowPosition(0, 0);
		glutCreateWindow("Point Cloud Viewer");
		//glutFullScreen();

		// Display Function
		glutDisplayFunc(display);

		// set the function to handle changes in screen size
		//glutReshapeFunc(OnReshape);
		int retCode = system("./png_visualizer.exe quadtree.png");

		// Properties Init
		myinit();
		glutMainLoop();
	}
	return 0;
}



