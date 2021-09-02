#include <GL/glut.h>
#include <iostream>
#include <iomanip> 
#include <memory>
#include <unistd.h>
#include <boost/algorithm/string.hpp>
#include "stb_image_write.h"

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
static const string LEAST_SQUARES = "LeastSquares";
static const string CONJUGATE = "Conjugate";

int					Field_Resolution_ = 64;
int					Multigrid_Iterations_ = 1;
int					Num_Threads_ = 8;
int					QuadTree_Levels_ = 8;
float 				Normal_Size_ = .2;
float 				Gaussian_Noise_ = 0;
string				Relative_Out_Path_Name_ = "out.png";
string				Relative_Quadtree_Out_Path_Name_;
string 				Relative_In_Path_Name_;
bool				Model_Generation_Mode_ = 0;
bool 				Model_Generated_ = 0;
bool 				Experimentation_Mode_ = 0;
bool				Fixed_Normal_Algorithm = 0;
bool				Fixed_Smoothness_Algorithm = 0;
bool 				Fixed_Solver_Method = 0;
bool				Fixed_Multigrid_Mode = 0;
Normal				Normal_Algorithm_ = Normal::sampling;
Smoothness			Smoothness_Algorithm_ = Smoothness::singleDimension;
Solver				Solver_Method_ = Solver::BiCGSTAB;
bool				Multigrid_ = 0;
bool				Print_Logs_ = 1;

void printUsage()
{
	std::cout << "Usage: " << std::endl << "  reconstruction.exe [ -f | -g ] [ Single Analysis | Model Generation | In Depth Study ]" << std::endl << 
	"Optional Parameters:" << std::endl <<  " \t[ -g | experiment mode (0 or 1) ]" << std::endl << " \t[ -n | normal's size ]" << std::endl << " \t[ -z | Gaussian noise (between 0-0.025 recommended) ]" << std::endl << "\t[ -r | resolution ]" << std::endl << "\t[ -a | normal algorithm (gradient or sampling) ]" << std::endl << "\t[ -s | smoothing algorithm (1D or 2D) ]" << 
	std::endl << "\t[ -x | solver method (Conjugate, Biconjugate or LeastSquares) ] " << std::endl << "\t[ -m | Multigrid Solving  Mode on (0 or 1) ] " << std::endl << "\t[ -i | Multigrid Iterations ] " << std::endl << "\t[ -t | Number Threads (default 8) ]" << std::endl << "\t[ -p | print logs (0 or 1) ]" << std::endl;
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
	while ((c = getopt (argc, argv, "ghf:n:z:r:a:s:x:m:i:l:t:p:q:")) != -1)
	{
		switch (c)
		{
			case 'g':
				Model_Generation_Mode_ = 1;
				break;

			case 'q':
				Experimentation_Mode_ = 1;
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

			case 'z':
				Gaussian_Noise_ = stof(optarg);
				break;

			case 'r':
				Field_Resolution_ = stoi(optarg);
				break;

			case 'a':
				Fixed_Normal_Algorithm = true;
				if (GRADIENT == optarg)
					Normal_Algorithm_ = Normal::gradient;
				break;

			case 's':
				Fixed_Smoothness_Algorithm = true;
				if (SMOOTHNESS == optarg)
					Smoothness_Algorithm_ = Smoothness::twoDimension;
				break;

			case 'x':
				Fixed_Solver_Method = true;
				if(optarg == LEAST_SQUARES)
					Solver_Method_ = Solver::LeastSquaresConjugateGradient;

				else if(optarg == CONJUGATE)
					Solver_Method_ = Solver::ConjugateGradient;
				
				break;

			case 'm':
				Fixed_Multigrid_Mode = true;
				Multigrid_ = stoi(optarg);
				break;

			case 'i':
				Multigrid_Iterations_ = stoi(optarg);
				break;

			case 'l':
				QuadTree_Levels_ = stoi(optarg);
				break;

			case 't':
				Num_Threads_ = stoi(optarg);
				break;

			case 'p':
				Print_Logs_ = stoi(optarg);
				break;

			case '?':
				if (optopt == 'f')
					fprintf (stderr, "Option -f requires a filepath as argument.\n");
				
				else if (optopt == 'n')
					fprintf (stderr, "Option -n requires a float value as argument.\n");				
				
				else if (optopt == 'z')
					fprintf (stderr, "Option -z requires a float value as argument.\n");	
				
				else if (optopt == 'r')
					fprintf (stderr, "Option -r requires an integer value as argument.\n");	

				else if (optopt == 'a')
					fprintf (stderr, "Option -a requires a normal type as argument.\n");	

				else if (optopt == 's')
					fprintf (stderr, "Option -s requires a smoothing type as argument.\n");	
				
				else if (optopt == 'm')
					fprintf (stderr, "Option -m requires a bool value as argument.\n");	

				else if (optopt == 'i')
					fprintf (stderr, "Option -i requires an integer value as argument.\n");	

				else if (optopt == 'l')
					fprintf (stderr, "Option -l requires an integer value as argument.\n");	

				else if (optopt == 't')
					fprintf (stderr, "Option -t requires an integer value as argument.\n");	

				else if (optopt == 'x')
					fprintf (stderr, "Option -x requires a solver type as argument.\n");	

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
			glVertex2f(1- mesh_->vertices_[i][0], 1- mesh_->vertices_[i][1]);

			glVertex2f(1 - mesh_->vertices_[i][0] - (Normal_Size_ *mesh_->normals_[i][0]), 1 - mesh_->vertices_[i][1] - (Normal_Size_ * mesh_->normals_[i][1]));
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

void takeScreenshot(const char *fname)
{
	GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    int x = viewport[0];
    int y = viewport[1];
    int width = viewport[2];
    int height = viewport[3];

    char *data = (char*) malloc((size_t) (width * height * 3)); // 3 components (R, G, B)

    if (!data)
        return;

    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, data);

    int saved = stbi_write_png(fname, width, height, 3, data, 0);

    free(data);
}

void computeSingleReconstruction(int &argc, char** argv)
{		
		ScalarField field;
		BiharmonicSolver solver;
		Eigen::VectorXd guess;
		Image *img;

		// Init & Set Weights
		//
		field.init(Field_Resolution_, Field_Resolution_);
		solver.setWeights(1.0f, 1.0f, 1.0f);

		// Reconstruct
		//
		SolverData s;		
		if(Multigrid_)
		{
			s = solver.computeMultigrid(*mesh_.get(), field, Multigrid_Iterations_, Normal_Algorithm_, Smoothness_Algorithm_, Solver_Method_, Num_Threads_, Print_Logs_);
		}
		else
		{
			s = solver.computeWith(*mesh_.get(), field, Normal_Algorithm_, Smoothness_Algorithm_, Solver_Method_, Num_Threads_, Print_Logs_);
		}
		// Create & Manage Directory Hierarchy
		//
		std::vector<std::string> results;
		std::vector<std::string> endResults;
		boost::algorithm::split(results, Relative_In_Path_Name_, boost::is_any_of("/"));
		boost::algorithm::split(endResults, results[results.size()-1], boost::is_any_of("."));

		std::string directoryName = "../singleAnalysis";
		std::string directoryCommand = "mkdir -p " + directoryName;
		system(directoryCommand.c_str());
		std::string subDirName = directoryName + "/analysis-" + endResults[0] + NORMAL_STRING[Normal_Algorithm_] + "-" + SMOOTHNESS_STRING[Smoothness_Algorithm_] + "-" + SOLVER_STRING[Solver_Method_] + "-" + RESULT_STRING[Multigrid_] + "Multigrid" + std::to_string(Num_Threads_) + "threads";
		std::string subDirCommand = "mkdir -p " + subDirName;
		system(subDirCommand.c_str());
		Relative_Quadtree_Out_Path_Name_ = subDirName + "/quadtree-" + Relative_Out_Path_Name_;
		Relative_Out_Path_Name_ = subDirName + "/" + Relative_Out_Path_Name_;
		img = field.toImage(16.0f, 0.0f);

		// Save PNG Image
		if(!img->savePNG(Relative_Out_Path_Name_))
		{
			cout << "[ERROR] Could not save file!" << endl;
			delete img;
			return;
		}
		delete img;

		// Save Text Reports
		ofstream myReport;
		std::string strReport = subDirName + "/report.txt";
		myReport.open (strReport);
		myReport << std::string(204, '-') << std::endl;
		myReport << "|" << setw(25) << "Relative In Path Name" << "|" << setw(10) << "Resolution" << "|" << setw(5) << "Noise" <<  "|" << setw(9) << "Normal" << "|" << setw(6) << "Smooth" << "|" << setw(29) << "Solver Method" << "|" << setw(10) << "Multigrid?" << "|" << setw(16) << "Multi Iterations" << "|" << setw(7) << "Threads" << "|" << setw(7) << "Solved?" << "|" << setw(8) << "Sys Time" << "|" << setw(8) << "Mat Time" << "|" << setw(13) << "Solv Res Time" << "|" << setw(13) << "Total Time" << "|" << setw(10) << "Iterations" << "|" << setw(12) << "Error" << "|" << std::endl;
		myReport << std::string(204, '-') << std::endl;
		time_t totalTime = s.systemBuildTime + s.matrixBuildTime + s.solverResolutionTime;
		myReport << "|" << setw(25) << Relative_In_Path_Name_ << "|" << setw(10) << Field_Resolution_ << "|" << setw(5) << Gaussian_Noise_ <<  "|" << setw(9) << NORMAL_STRING[Normal_Algorithm_] << "|" << setw(6) << SMOOTHNESS_STRING[Smoothness_Algorithm_] << "|" << setw(29) << SOLVER_STRING[Solver_Method_] << "|" << setw(10) << RESULT_STRING[Multigrid_] << "|" << setw(16) << Multigrid_Iterations_ << "|" << setw(7) << Num_Threads_ << "|" << setw(7) << RESULT_STRING[s.isSolved] << "|" << setw(8) << s.systemBuildTime << "|" << setw(8) <<s.matrixBuildTime << "|" << setw(13) <<s.solverResolutionTime << "|" << setw(13) << totalTime << "|" << setw(10) << s.iterations << "|" << setw(12) << s.error << "|" << std::endl;
		myReport << std::string(204, '-') << std::endl;
		myReport.close();

		ofstream excel;
		std::string strExcel = subDirName + "/report.csv";
		excel.open(strExcel);
		excel << "Relative In Path Name" << ", " << "Resolution" << ", " << "Noise" <<  ", " << "Normal" << ", " << "Smooth" << ", " << "Solver Method" << ", " << "Multigrid?" << ", " << "Multi Iterations" << ", " << "Threads" << ", " << "Is Solved?" << ", " << "Syst Time" << ", " << "Mat Time" << ", " << "Solver Res Time" << ", " << "Total Time" << ", " << "Iterations" << ", " << "Error" << std::endl;
		std::cout.precision(5);
		excel << Relative_In_Path_Name_ << ", " << Field_Resolution_ << ", " << Gaussian_Noise_ <<  ", " << NORMAL_STRING[Normal_Algorithm_] << ", " << SMOOTHNESS_STRING[Smoothness_Algorithm_] << ", " << SOLVER_STRING[Solver_Method_] << ", " << RESULT_STRING[Multigrid_] << ", " << Multigrid_Iterations_ << ", " << Num_Threads_ << ", " << RESULT_STRING[s.isSolved] << ", " << s.systemBuildTime << ", " <<s.matrixBuildTime << ", " <<s.solverResolutionTime << ", " << totalTime << ", " << s.iterations << ", " << s.error << std::endl;
		excel.close();

		std::cout.precision(5);
		std::cout << "[RESULT] " << Relative_In_Path_Name_ << ", " << Field_Resolution_ << ", " << Gaussian_Noise_ <<  ", " << NORMAL_STRING[Normal_Algorithm_] << ", " << SMOOTHNESS_STRING[Smoothness_Algorithm_] << ", " << SOLVER_STRING[Solver_Method_] << ", " << RESULT_STRING[Multigrid_] << ", " << Multigrid_Iterations_ << ", " << Num_Threads_ << ", " << RESULT_STRING[s.isSolved] << ", " << s.systemBuildTime << ", " <<s.matrixBuildTime << ", " <<s.solverResolutionTime << ", " << s.iterations << ", " << s.error << std::endl;
		
		// Now Compute the Quadree Solution
		//

		Quadtree qtree;
		qtree.compute(*mesh_.get(), QuadTree_Levels_, field, Normal_Algorithm_, Smoothness_Algorithm_, Solver_Method_, Num_Threads_, Print_Logs_);
		
		img = field.toImage(16.0f, 0.0f);
		if(!img->savePNG(Relative_Quadtree_Out_Path_Name_))
		{
			cout << "[ERROR] Could not save file!" << endl;
			delete img;
			return;
		}
		delete img;

		Image qtreeImg;
		qtreeImg.init(1024, 1024);
		qtree.draw(qtreeImg);
		std::string strQuadtree = subDirName + "/quadtree.png";
		qtreeImg.savePNG(strQuadtree);
		
		// Initialize Viewer
		//
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(1000, 1000);
		glutInitWindowPosition(0, 0);
		glutCreateWindow("Point Cloud Viewer");
		//glutFullScreen();

		// Display Function
		glutDisplayFunc(display);

		// set the function to handle changes in screen size
		//glutReshapeFunc(OnReshape);
		std::string exePath = "./png_visualizer " + Relative_Out_Path_Name_ + " " + std::to_string(Field_Resolution_) + " " + std::to_string(Field_Resolution_) + " &";
		int retCode = system(exePath.c_str());
		std::string command = "./png_visualizer " + strQuadtree + " 1024 1024 &";
		int retCode_quadtree = system(command.c_str());

		// Take screenshoot
		std::string strScreenshot = subDirName + "/screenshot.png";
		takeScreenshot(strScreenshot.c_str());

		// Properties Init
		myinit();
		glutMainLoop();
		std::string kill = "pkill -9 -f" + exePath;
		system(kill.c_str());
}

void computeMultipleReconstruction(int &argc, char** argv)
{
	// Create & Manage Directory Hierarchy
	//
	std::vector<std::string> results;
	std::vector<std::string> endResults;
	boost::algorithm::split(results, Relative_In_Path_Name_, boost::is_any_of("/"));
	boost::algorithm::split(endResults, results[results.size()-1], boost::is_any_of("."));

	std::string directoryName = "../multipleAnalysis";
	std::string directoryCommand = "mkdir -p " + directoryName;
	system(directoryCommand.c_str());
	
	
	directoryName = directoryName + "/multiple-analysis-" + endResults[0];
	directoryCommand = "mkdir -p " + directoryName;
	system(directoryCommand.c_str());
	
	std::vector<SolverData> solverOuts;

	for(size_t isMultiIt = 0; isMultiIt < 2; ++isMultiIt)
	{	
		for(size_t normalAlgIt = 0; normalAlgIt < 2; ++normalAlgIt)
		{
			for(size_t solverAlgIt = 0; solverAlgIt < 3; ++solverAlgIt)
			{
				for(size_t smoothAlgIt = 0; smoothAlgIt < 2; ++smoothAlgIt)	
				{
					for(float numberThreadsIt = Num_Threads_; numberThreadsIt > 1; numberThreadsIt/=2)
					{					
						ScalarField field;
						BiharmonicSolver solver;
						Eigen::VectorXd guess;
						Image *img;

						
						// Init & Set Weights
						//

						if(isMultiIt)	 
							field.init((Field_Resolution_/(Multigrid_Iterations_+1)+1), Field_Resolution_/((Multigrid_Iterations_+1)+1));
						else
							field.init(Field_Resolution_, Field_Resolution_);
						solver.setWeights(1.0f, 1.0f, 1.0f);

						// Reconstruct
						//
						SolverData s;		
						if(isMultiIt)
						{
							s = solver.computeMultigrid(*mesh_.get(), field, Multigrid_Iterations_, normalAlgIt, smoothAlgIt, solverAlgIt, numberThreadsIt, Print_Logs_);
						}
						else
						{
							s = solver.computeWith(*mesh_.get(), field, normalAlgIt, smoothAlgIt, solverAlgIt, numberThreadsIt, Print_Logs_);
						}

						s.resolution = Field_Resolution_;
						s.gaussianNoise = Gaussian_Noise_;
						s.normalAlgorithm = NORMAL_STRING[normalAlgIt];
						s.smoothingAlgorithm = SMOOTHNESS_STRING[smoothAlgIt];
						s.solverMethod = SOLVER_STRING[solverAlgIt];
						s.isMultigrid = RESULT_STRING[isMultiIt];
						s.multigridIterations = Multigrid_Iterations_;
						s.numberThreads = numberThreadsIt;

						// Create & Manage Directory Hierarchy
						//
						std::string singleSubDirCommand = "mkdir -p " + directoryName + "/singleAnalysis";
						system(singleSubDirCommand.c_str());
						std::string subDirName = directoryName + "/singleAnalysis/analysis-" + endResults[0] + s.normalAlgorithm + "-" + s.smoothingAlgorithm + "-" + s.solverMethod + "-" + RESULT_STRING[isMultiIt] + "Multigrid-" + std::to_string(s.numberThreads) + "threads";
						std::string subDirCommand = "mkdir -p " + subDirName;
						system(subDirCommand.c_str());
						Relative_Quadtree_Out_Path_Name_ = subDirName + "/quadtree-" + Relative_Out_Path_Name_;
						std::string relative_Out_Path_Name = subDirName + "/" + Relative_Out_Path_Name_;
						img = field.toImage(16.0f, 0.0f);

						// Save PNG Image
						if(!img->savePNG(relative_Out_Path_Name))
						{
							cout << "[ERROR] Could not save file!" << endl;
							delete img;
							return;
						}
						delete img;

						// Save Text Reports
						ofstream myReport;
						std::string strReport = subDirName + "/report.txt";
						myReport.open (strReport);
						myReport << std::string(204, '-') << std::endl;
						myReport << "|" << setw(25) << "Relative In Path Name" << "|" << setw(10) << "Resolution" << "|" << setw(5) << "Noise" <<  "|" << setw(9) << "Normal" << "|" << setw(6) << "Smooth" << "|" << setw(29) << "Solver Method" << "|" << setw(10) << "Multigrid?" << "|" << setw(16) << "Multi Iterations" << "|" << setw(7) << "Threads" << "|" << setw(7) << "Solved?" << "|" << setw(8) << "Sys Time" << "|" << setw(8) << "Mat Time" << "|" << setw(13) << "Solv Res Time" << "|" << setw(13) << "Total Time" << "|" << setw(10) << "Iterations" << "|" << setw(12) << "Error" << "|" << std::endl;
						myReport << std::string(204, '-') << std::endl;
						time_t totalTime = s.systemBuildTime + s.matrixBuildTime + s.solverResolutionTime;
						myReport << "|" << setw(25) << Relative_In_Path_Name_ << "|" << setw(10) << Field_Resolution_ << "|" << setw(5) << Gaussian_Noise_ <<  "|" << setw(9) << NORMAL_STRING[Normal_Algorithm_] << "|" << setw(6) << SMOOTHNESS_STRING[Smoothness_Algorithm_] << "|" << setw(29) << SOLVER_STRING[Solver_Method_] << "|" << setw(10) << RESULT_STRING[Multigrid_] << "|" << setw(16) << Multigrid_Iterations_ << "|" << setw(7) << Num_Threads_ << "|" << setw(7) << RESULT_STRING[s.isSolved] << "|" << setw(8) << s.systemBuildTime << "|" << setw(8) <<s.matrixBuildTime << "|" << setw(13) <<s.solverResolutionTime << "|" << setw(13) << totalTime << "|" << setw(10) << s.iterations << "|" << setw(12) << s.error << "|" << std::endl;
						myReport << std::string(204, '-') << std::endl;
						myReport.close();

						ofstream excel;
						std::string strExcel = subDirName + "/report.csv";
						excel.open(strExcel);
						excel << "Relative In Path Name" << ", " << "Resolution" << ", " << "Noise" <<  ", " << "Normal" << ", " << "Smooth" << ", " << "Solver Method" << ", " << "Multigrid?" << ", " << "Multi Iterations" << ", " << "Threads" << ", " << "Is Solved?" << ", " << "Syst Time" << ", " << "Mat Time" << ", " << "Solver Res Time" << ", " << "Total Time" << " ," << "Iterations" << ", " << "Error" << std::endl;
						std::cout.precision(5);
						excel << Relative_In_Path_Name_ << ", " << s.resolution << ", " << s.gaussianNoise <<  ", " << s.normalAlgorithm << ", " << s.smoothingAlgorithm << ", " << s.solverMethod << ", " << s.isMultigrid << ", " << s.multigridIterations << ", " << s.numberThreads << ", " << RESULT_STRING[s.isSolved] << ", " << s.systemBuildTime << ", " << s.matrixBuildTime << ", " << s.solverResolutionTime << ", " << totalTime << s.iterations << ", " << s.error << std::endl;
						excel.close();

						std::cout.precision(5);
						std::cout << "[RESULT] " << Relative_In_Path_Name_ << ", " << s.resolution << ", " << s.gaussianNoise <<  ", " << s.normalAlgorithm << ", " << s.smoothingAlgorithm << ", " << s.solverMethod << ", " << s.isMultigrid << ", " << s.multigridIterations << ", " << s.numberThreads << ", " << RESULT_STRING[s.isSolved] << ", " << s.systemBuildTime << ", " << s.matrixBuildTime << ", " << s.solverResolutionTime << ", " << s.iterations << ", " << s.error << std::endl;

						// Now Compute the Quadree Solution
						//

						Quadtree qtree;
						qtree.compute(*mesh_.get(), QuadTree_Levels_, field, normalAlgIt, smoothAlgIt, solverAlgIt, numberThreadsIt, Print_Logs_);
						
						img = field.toImage(16.0f, 0.0f);
						if(!img->savePNG(Relative_Quadtree_Out_Path_Name_))
						{
							cout << "[ERROR] Could not save file!" << endl;
							delete img;
							return;
						}
						delete img;

						Image qtreeImg;
						qtreeImg.init(1024, 1024);
						qtree.draw(qtreeImg);
						std::string strQuadtree = subDirName + "/quadtree.png";
						qtreeImg.savePNG(strQuadtree);

						std::cout << endl;

						solverOuts.push_back(s);
					}
				}
			}
		}
	}
			// Save Text Reports
		ofstream myReport;
		std::string strReport = directoryName + "/report.txt";
		myReport.open (strReport);

		ofstream excel;
		std::string strExcel = directoryName + "/report.csv";
		excel.open(strExcel);

		myReport << std::string(204, '-') << std::endl;
		myReport << "|" << setw(25) << "Relative In Path Name" << "|" << setw(10) << "Resolution" << "|" << setw(5) << "Noise" <<  "|" << setw(9) << "Normal" << "|" << setw(6) << "Smooth" << "|" << setw(29) << "Solver Method" << "|" << setw(10) << "Multigrid?" << "|" << setw(16) << "Multi Iterations" << "|" << setw(7) << "Threads" << "|" << setw(7) << "Solved?" << "|" << setw(8) << "Sys Time" << "|" << setw(8) << "Mat Time" << "|" << setw(13) << "Solv Res Time" << "|" << setw(13) << "Total Time" << "|" << setw(10) << "Iterations" << "|" << setw(12) << "Error" << "|" << std::endl;

		excel << "Relative In Path Name" << ", " << "Resolution" << ", " << "Noise" <<  ", " << "Normal" << ", " << "Smooth" << ", " << "Solver Method" << ", " << "Multigrid?" << ", " << "Multi Iterations" << ", " << "Threads" << ", " << "Is Solved?" << ", " << "Syst Time" << ", " << "Mat Time" << ", " << "Solver Res Time" << ", " << "Total Time" << ", " << "Iterations" << ", " << "Error" << std::endl;
		std::cout.precision(5);

		SolverData bestCombi;
		bestCombi.solverResolutionTime = 9999999999999999;

		for(size_t i = 0; i < solverOuts.size(); ++i)
		{
			myReport << std::string(204, '-') << std::endl;
			time_t totalTime = solverOuts[i].systemBuildTime + solverOuts[i].matrixBuildTime + solverOuts[i].solverResolutionTime;
			myReport << "|" << setw(25) << Relative_In_Path_Name_ << "|" << setw(10) << solverOuts[i].resolution << "|" << setw(5) << solverOuts[i].gaussianNoise <<  "|" << setw(9) << solverOuts[i].normalAlgorithm << "|" << setw(6) << solverOuts[i].smoothingAlgorithm << "|" << setw(29) << solverOuts[i].solverMethod << "|" << setw(10) << solverOuts[i].isMultigrid << "|" << setw(16) << Multigrid_Iterations_ << "|" << setw(7) << solverOuts[i].numberThreads << "|" << setw(7) << RESULT_STRING[solverOuts[i].isSolved] << "|" << setw(8) << solverOuts[i].systemBuildTime << "|" << setw(8) << solverOuts[i].matrixBuildTime << "|" << setw(13) << solverOuts[i].solverResolutionTime << "|" << setw(13) << totalTime <<  "|" << setw(10) <<  solverOuts[i].iterations << "|" << setw(12) << solverOuts[i].error << "|" << std::endl;
			excel << Relative_In_Path_Name_ << ", " << solverOuts[i].resolution << ", " << solverOuts[i].gaussianNoise <<  ", " << solverOuts[i].normalAlgorithm << ", " << solverOuts[i].smoothingAlgorithm << ", " << solverOuts[i].solverMethod << ", " << solverOuts[i].isMultigrid << ", " << totalTime << ", " << solverOuts[i].multigridIterations << ", " << solverOuts[i].numberThreads << ", " << RESULT_STRING[solverOuts[i].isSolved] << ", " << solverOuts[i].systemBuildTime << ", " << solverOuts[i].matrixBuildTime << ", " << solverOuts[i].solverResolutionTime << ", " << solverOuts[i].iterations << ", " << solverOuts[i].error << std::endl;

			if(totalTime < bestCombi.systemBuildTime + bestCombi.matrixBuildTime + bestCombi.solverResolutionTime && solverOuts[i].isSolved == true)
				bestCombi = solverOuts[i];
		}

		myReport << std::string(204, '-') << std::endl << endl;
		myReport << "Fastest Combination With Solution: " << endl;
		myReport << "\t"  << "- Normal Method:     " << bestCombi.normalAlgorithm << endl; 
		myReport << "\t"  << "- Smoothness Method: " << bestCombi.smoothingAlgorithm << endl; 
		myReport << "\t"  << "- Solver System:     " << bestCombi.solverMethod << endl; 
		myReport << "\t"  << "- MultiGrid Needed?  " << bestCombi.isMultigrid << endl; 
		/*if(bestCombi.isMultigrid == "Yes")
			myReport << "\t"  << "- MultiGrid Iter:    " << bestCombi.isMultigrid << endl; */
		myReport << "\t"  << "- Number Threads:    " << bestCombi.numberThreads << endl; 
/*		time_t totalTime = bestCombi.systemBuildTime + bestCombi.matrixBuildTime + bestCombi.solverResolutionTime;
		myReport << "\t"  << "- Total Time (ms.):  " << totalTime << endl << endl;*/

		myReport.close();
		excel.close();

		// Compress and Remove subdirectory
		std::string compressDirectoryCommand = "tar -zcvf " + directoryName + "/singleAnalysis.tgz " + directoryName + "/singleAnalysis";
		std::cout << compressDirectoryCommand << endl;
		system(compressDirectoryCommand.c_str());
		std::string removeDirectoryCommand = "rm -rf " + directoryName + "/singleAnalysis/";
		system(removeDirectoryCommand.c_str());

	// Initialize Viewer
	//
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(1000, 1000);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Point Cloud Viewer");
	//glutFullScreen();

	// Display Function
	glutDisplayFunc(display);
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
			res = data_representation::ReadFromTXT(file, mesh_.get(), Gaussian_Noise_);
			Model_Generated_ = true;
		}

		else if (type.compare("svg") == 0)
		{
			res = data_representation::ReadFromSVG(file, mesh_.get());
			std::cout << "[DATA] New file .txt generated, please run the script again with the new file" << std::endl;
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

	if(Model_Generated_ && !Experimentation_Mode_)
	{
		computeSingleReconstruction(argc, argv);
	}
	else if(Model_Generated_ && Experimentation_Mode_)
	{
		computeMultipleReconstruction(argc, argv);
	}

	return 0;
}



