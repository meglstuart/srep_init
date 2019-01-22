#include "srep_init.h"

#include <cstdio>
#include <iostream>


srep_init::srep_init(double d, double smooth, int max)
{
  this->dt = d;
  this->smoothAmount = smooth;
  this->maxIter = max;
}

int srep_init::process_inputMesh()
{
  // size_t last_dot = this->inputMesh.rfind('.');
  // if (last_dot == std::string::npos)
  // {
  //   std::cerr<<"Error: No file extension found in "<<
  //   this->inputMesh<<std::endl;
  //   return false;
  // }
  // std::string extension = this->inputMesh.substr(last_dot+1);
  //
  // if (extension == "off" || extension =="OFF")
  // {
  //   Eigen::MatrixXd V;
  //   Eigen::MatrixXi F;
  //   if (!igl::readOFF(this->inputMesh, V, F))
  //   {
  //     return false;
  //   }
  //   data().set_mesh(V,F);
  // }
  // else if (extension == "obj" || extension =="OBJ")
  // {
  //   Eigen::MatrixXd corner_normals;
  //   Eigen::MatrixXi fNormIndices;
  //
  //   Eigen::MatrixXd UV_V;
  //   Eigen::MatrixXi UV_F;
  //   Eigen::MatrixXd V;
  //   Eigen::MatrixXi F;
  //
  //   if (!(
  //     igl::readOBJ(
  //       this->inputMesh,
  //       V, UV_V, corner_normals, F, UV_F, fNormIndices)))
  //       {
  //         return false;
  //       }
  //
  //       data().set_mesh(V,F);
  //       data().set_uv(UV_V,UV_F);
  //
  //     }
  //     else
  //     {
  //       // unrecognized file type
  //       printf("Error: %s is not a recognized file type.\n",extension.c_str());
  //       return false;
  //     }

}
