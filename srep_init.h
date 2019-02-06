#ifndef __SREP_INIT_H
#define __SREP_INIT_H

#include <string>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <igl/opengl/glfw/Viewer.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

class srep_init
{
public:
  srep_init(double d, double smooth, int max);
  int set_mesh(Eigen::MatrixXd V, Eigen::MatrixXi F);
  int step_forwardflow();
  int update_viewer(igl::opengl::glfw::Viewer *);
  int fit_ellipsoid(vtkSmartPointer<vtkPolyData> polydata_smooth, int);
  int write_ellipsoid();
  int generate_ellipsoid_srep();
  int write_header(int);
  int backward_flow();

  double dt;
  double smoothAmount;
  int max_iter;
  int iter = 0;
  int nRows = 5;
  int nCols = 5;
  double q=1;
  double tol = 0.01;
  std::string output_folder = "";
  std::string input_mesh;
  Eigen::VectorXd radii;    //Radii of best fitting ellipsoid
  Eigen::MatrixXd rotation;
  Eigen::MatrixXd center;
  Eigen::MatrixXd V;        //Initial mesh vertices
  Eigen::MatrixXd U;        //Latest mesh vertices
  Eigen::MatrixXi F;        //Initial mesh faces
  Eigen::SparseMatrix<double> L;   //Discrete laplacian cotangent stiffness matrix
  Eigen::MatrixXd ell_U;        //ellipsoid mesh vertices
  Eigen::MatrixXi ell_F;        //ellipsoid mesh faces
  vtkSmartPointer<vtkPolyData> best_fitting_ellipsoid_polydata = vtkSmartPointer<vtkPolyData>::New();
};

#endif
