#ifndef __SREP_INIT_H
#define __SREP_INIT_H

#include <string>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <igl/opengl/glfw/Viewer.h>

class srep_init
{
public:
  srep_init(double d, double smooth, int max);
  int set_mesh(Eigen::MatrixXd V, Eigen::MatrixXi F);
  int step_forwardflow();
  int update_viewer(igl::opengl::glfw::Viewer *);

  double dt = 0.1f;
  double smoothAmount = 0.1f;
  int maxIter = 10;
  int iter = 0;
  double q=1;
  std::string inputMesh = "../test_data/bunny.off";
  Eigen::MatrixXd V;        //Initial mesh vertices
  Eigen::MatrixXd U;        //Latest mesh vertices
  Eigen::MatrixXi F;        //Initial mesh faces
  Eigen::SparseMatrix<double> L;   //Discrete laplacian cotangent stiffness matrix

};

#endif
