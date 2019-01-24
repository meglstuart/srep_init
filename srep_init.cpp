#include "srep_init.h"

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/barycenter.h>
#include <igl/writeOFF.h>


#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkMassProperties.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCell.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkIdList.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunctionSource.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkCleanPolyData.h>

srep_init::srep_init(double d, double smooth, int max)
{
  this->dt = d;
  this->smoothAmount = smooth;
  this->maxIter = max;
};

int srep_init::set_mesh(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
  this->V = V;
  this->F = F;
  this->U = V;
  igl::cotmatrix(V,F,this->L);

};

int srep_init::update_viewer(igl::opengl::glfw::Viewer *viewer)
{
  viewer->data().set_vertices(U);
  viewer->data().compute_normals();
  viewer->core.align_camera_center(U,F);
};

int srep_init::step_forwardflow()
{
  Eigen::MatrixXd rotation;
  Eigen::VectorXd radii;
  char temp[128];
  // compute mean curvature flow
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(U,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
  // Solve (M-delta*L) U = M*U
  const auto & S = (M - dt*L);
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
  assert(solver.info() == Eigen::Success);
  U = solver.solve(M*U).eval();
  // Compute centroid and subtract (also important for numerics)
  Eigen::VectorXd dblA;
  igl::doublearea(U,F,dblA);
  double area = 0.5*dblA.sum();
  Eigen::MatrixXd BC;
  igl::barycenter(U,F,BC);
  Eigen::RowVector3d centroid(0,0,0);
  for(int i = 0;i<BC.rows();i++)
  {
      centroid += 0.5*dblA(i)/area*BC.row(i);
  }
  U.rowwise() -= centroid;
  // Normalize to unit surface area (important for numerics)
  U.array() /= sqrt(area);
  sprintf(temp,"%d", ++iter);
  std::string prefix = temp;

  vtkSmartPointer<vtkPoints> points =
      vtkSmartPointer<vtkPoints>::New();
  for(int i = 0; i < U.rows(); ++i) {
      points->InsertNextPoint(U(i,0), U(i,1), U(i,2));
  }

  // construct polys from F
  vtkSmartPointer<vtkCellArray> polys =
      vtkSmartPointer<vtkCellArray>::New();
  for(int i = 0; i < F.rows(); ++i) {
      vtkSmartPointer<vtkIdList> ids =
          vtkSmartPointer<vtkIdList>::New();
      ids->InsertNextId(F(i,0));
      ids->InsertNextId(F(i,1));
      ids->InsertNextId(F(i,2));
      polys->InsertNextCell(ids);
  }

  vtkSmartPointer<vtkPolyData> polydata =
      vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPolys(polys);
  polydata->SetPoints(points);

  // polydata writer
  vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();

  // smoother
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
      vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  smoother->SetNumberOfIterations(20);
  smoother->BoundarySmoothingOff();
  smoother->FeatureEdgeSmoothingOff();
  smoother->SetPassBand(smoothAmount);
  smoother->NonManifoldSmoothingOn();
  smoother->NormalizeCoordinatesOn();



  // smooth polydata
  smoother->SetInputData(polydata);
  smoother->Update();
  vtkSmartPointer<vtkPolyData> polydata_smooth = smoother->GetOutput();

  // mass filter
  vtkSmartPointer<vtkMassProperties> mass = vtkSmartPointer<vtkMassProperties>::New();
  mass->SetInputData(polydata_smooth);
  mass->Update();
  double current_volume = mass->GetVolume();

  // set U to smoothed points
  for(int i = 0; i < U.rows(); ++i) {
      double p[3];
      polydata_smooth->GetPoint(i,p);
      U(i,0) = p[0];
      U(i,1) = p[1];
      U(i,2) = p[2];
  }


  // std::sort(distance_vector.begin(), distance_vector.end());
  // std::string off_filename = "temp_off/temp_" + prefix;
  // off_filename += ".off";
  // igl::writeOFF(off_filename, U, F);
  // int quantile_idx = static_cast<int>(U.rows() * 0.95);
  //
  // q = distance_vector[quantile_idx];
  // cout << "iter " << iter << ": " << q << endl;
  // std::string vtk_filename = "ell/ell_" + prefix;
  // vtk_filename += ".vtk";
  // writer->SetFileName(vtk_filename.c_str());
  // writer->SetInputData(best_fitting_ellipsoid_polydata);
  // writer->Update();

  std::string vtk_filename = "temp_vtk/" + prefix;
  vtk_filename += ".vtk";
  writer->SetFileName(vtk_filename.c_str());
  writer->SetInputData(polydata_smooth);
  writer->Update();
  // std::string cleaned_vtk_filename = "cleaned/cleaned_" + prefix;
  // cleaned_vtk_filename+=".vtk";
  // writer->SetFileName(cleaned_vtk_filename.c_str());
  // writer->SetInputData(cleaned_polydata);
  // writer->Update();
  char buffer[128];
  // sprintf(buffer, "%f, %f, %f\n", radii(0)*volume_factor, radii(1) * volume_factor, radii(2) * volume_factor);
  // radii_file << buffer;


  return true;
};

// srep_init::get_best_fitting_ellipsoid()
// {
//   vtkSmartPointer<vtkCleanPolyData> cleaner =
//       vtkSmartPointer<vtkCleanPolyData>::New();
//   cleaner->SetTolerance(0.05);
//   cleaner->SetInputData(polydata_smooth);
//   cleaner->Update();
//   vtkSmartPointer<vtkPolyData> cleaned_polydata = cleaner->GetOutput();
//   Eigen::MatrixXd U_temp(cleaned_polydata->GetNumberOfPoints(),3);
//   // set U to smoothed points
//   for(int i = 0; i < U_temp.rows(); ++i) {
//       double p[3];
//       cleaned_polydata->GetPoint(i,p);
//       U_temp(i,0) = p[0];
//       U_temp(i,1) = p[1];
//       U_temp(i,2) = p[2];
//   }
//   // filter out clustered points
//   Eigen::MatrixXd cog = U_temp.colwise().mean();
//   Eigen::MatrixXd U_centered = U_temp - cog.replicate(U_temp.rows(), 1); // N by 3
//   Eigen::MatrixXd U_transposed = U_centered.transpose();
//   Eigen::Matrix3d U_second_moment = U_transposed * U_centered;
//
//   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(U_second_moment);
//   rotation = es.eigenvectors(); // 3 by 3 rotation matrix
//   radii = es.eigenvalues();
//   radii(0) = sqrt(radii(0));
//   radii(1) = sqrt(radii(1));
//   radii(2) = sqrt(radii(2));
//
//   double ellipsoid_volume = 4 / 3.0 * M_PI * radii(0) * radii(1) * radii(2);
//   double volume_factor = pow(current_volume / ellipsoid_volume, 1.0 / 3.0);
//   radii(0) *= volume_factor;
//   radii(1) *= volume_factor;
//   radii(2) *= volume_factor;
//   // obtain the best fitting ellipsoid from the second moment matrix
//   vtkSmartPointer<vtkParametricEllipsoid> ellipsoid =
//       vtkSmartPointer<vtkParametricEllipsoid>::New();
//   ellipsoid->SetXRadius(radii(0));
//   ellipsoid->SetYRadius(radii(1));
//   ellipsoid->SetZRadius(radii(2));
//
//   vtkSmartPointer<vtkParametricFunctionSource> parametric_function =
//       vtkSmartPointer<vtkParametricFunctionSource>::New();
//   parametric_function->SetParametricFunction(ellipsoid);
//   parametric_function->SetUResolution(30);
//   parametric_function->SetVResolution(30);
//   parametric_function->Update();
//   vtkSmartPointer<vtkPolyData> ellipsoid_polydata = parametric_function->GetOutput();
//
//   // Get ellipsoid points into the matrix
//   Eigen::MatrixXd ellipsoid_points_matrix(ellipsoid_polydata->GetNumberOfPoints(), 3);
//   for(int i = 0; i < ellipsoid_polydata->GetNumberOfPoints(); ++i) {
//       double p[3];
//       ellipsoid_polydata->GetPoint(i,p);
//       ellipsoid_points_matrix(i,0) = p[0];
//       ellipsoid_points_matrix(i,1) = p[1];
//       ellipsoid_points_matrix(i,2) = p[2];
//   }
//   // rotate the points
//   Eigen::MatrixXd rotated_ellipsoid_points = rotation * (ellipsoid_points_matrix.transpose());
//   rotated_ellipsoid_points.transposeInPlace(); // n x 3
//   // translate the points
//   Eigen::MatrixXd translated_points = rotated_ellipsoid_points + cog.replicate(rotated_ellipsoid_points.rows(),1);
//
//   // convert eigen matrix to vtk polydata
//   vtkSmartPointer<vtkPolyData> best_fitting_ellipsoid_polydata =
//       vtkSmartPointer<vtkPolyData>::New();
//   vtkSmartPointer<vtkPoints> best_fitting_ellipsoid_points =
//       vtkSmartPointer<vtkPoints>::New();
//   for(int i = 0; i < translated_points.rows(); ++i) {
//       double p[3] = {translated_points(i,0), translated_points(i,1), translated_points(i,2)};
//       best_fitting_ellipsoid_points->InsertNextPoint(p);
//   }
//   best_fitting_ellipsoid_polydata->SetPoints(best_fitting_ellipsoid_points);
//   best_fitting_ellipsoid_polydata->SetPolys(ellipsoid_polydata->GetPolys());
//   best_fitting_ellipsoid_polydata->Modified();
//
//   vtkSmartPointer<vtkImplicitPolyDataDistance> implicit_distance_filter =
//       vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
//   implicit_distance_filter->SetInput(best_fitting_ellipsoid_polydata);
//   std::vector<double> distance_vector;
//
//   for(int i = 0; i < U.rows(); ++i) {
//       double p[3];
//       polydata->GetPoint(i,p);
//       double current_distance = abs(implicit_distance_filter->EvaluateFunction(p));
//       distance_vector.push_back(current_distance);
//   }
//
// }
