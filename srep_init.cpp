#include "srep_init.h"

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

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
#include <vtkTriangleFilter.h>
#include <vtkCellArray.h>
#include <vtkCell.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkIdList.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunctionSource.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkCleanPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkVector.h>
#include <vtkLine.h>
#include <vtkQuad.h>
#include <vtkCellArray.h>
#include <vtkObjectBase.h>
#include <vtksys/SystemTools.hxx>

srep_init::srep_init(double d, double smooth, int max)
{
  this->dt = d;
  this->smoothAmount = smooth;
  this->max_iter = max;
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
  viewer->data().clear();
  viewer->data().set_mesh(U,F);
  viewer->data().compute_normals();
  viewer->core.align_camera_center(U,F);
};

int srep_init::step_forwardflow()
{

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

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for(int i = 0; i < U.rows(); ++i) {
    points->InsertNextPoint(U(i,0), U(i,1), U(i,2));
  }

  // construct polys from F
  vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
  for(int i = 0; i < F.rows(); ++i) {
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    ids->InsertNextId(F(i,0));
    ids->InsertNextId(F(i,1));
    ids->InsertNextId(F(i,2));
    polys->InsertNextCell(ids);
  }

  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPolys(polys);
  polydata->SetPoints(points);


  // polydata writer
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();

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

  // set U to smoothed points
  for(int i = 0; i < U.rows(); ++i) {
    double p[3];
    polydata_smooth->GetPoint(i,p);
    U(i,0) = p[0];
    U(i,1) = p[1];
    U(i,2) = p[2];
  }

  std::string vtk_filename = output_folder + "forward/" + prefix;
  vtk_filename += ".vtk";
  writer->SetFileName(vtk_filename.c_str());
  writer->SetInputData(polydata_smooth);
  writer->Update();

  this->fit_ellipsoid(polydata_smooth);
  return true;
};

int srep_init::fit_ellipsoid(vtkSmartPointer<vtkPolyData> polydata_smooth)
{
  // Eigen::MatrixXd rotation;
  // mass filter
  vtkSmartPointer<vtkMassProperties> mass = vtkSmartPointer<vtkMassProperties>::New();
  mass->SetInputData(polydata_smooth);
  mass->Update();
  double current_volume = mass->GetVolume();

  // filter out clustered points for smaller dimension matrix operations
  vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetTolerance(0.02);
  cleaner->SetInputData(polydata_smooth);
  cleaner->Update();
  vtkSmartPointer<vtkPolyData> cleaned_polydata = cleaner->GetOutput();
  Eigen::MatrixXd U_temp(cleaned_polydata->GetNumberOfPoints(),3);
  // Eigen::MatrixXd U_temp(polydata_smooth->GetNumberOfPoints(),3);

  for(int i = 0; i < U_temp.rows(); ++i) {
      double p[3];
      cleaned_polydata->GetPoint(i,p);
      // polydata_smooth->GetPoint(i,p);
      U_temp(i,0) = p[0];
      U_temp(i,1) = p[1];
      U_temp(i,2) = p[2];
  }
  Eigen::MatrixXd cog = U_temp.colwise().mean();
  Eigen::MatrixXd U_centered = U_temp - cog.replicate(U_temp.rows(), 1); // N by 3
  Eigen::MatrixXd U_transposed = U_centered.transpose();
  Eigen::Matrix3d U_second_moment = U_transposed * U_centered;
  center = cog;

  // obtain the best fitting ellipsoid from the second moment matrix
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(U_second_moment);
  rotation = es.eigenvectors(); // 3 by 3 rotation matrix
  radii = es.eigenvalues();
  radii(0) = sqrt(radii(0));
  radii(1) = sqrt(radii(1));
  radii(2) = sqrt(radii(2));

  double ellipsoid_volume = 4 / 3.0 * M_PI * radii(0) * radii(1) * radii(2);
  double volume_factor = pow(current_volume / ellipsoid_volume, 1.0 / 3.0);
  radii(0) *= volume_factor;
  radii(1) *= volume_factor;
  radii(2) *= volume_factor;

  vtkSmartPointer<vtkParametricEllipsoid> ellipsoid = vtkSmartPointer<vtkParametricEllipsoid>::New();
  ellipsoid->SetXRadius(radii(0));
  ellipsoid->SetYRadius(radii(1));
  ellipsoid->SetZRadius(radii(2));

  vtkSmartPointer<vtkParametricFunctionSource> parametric_function = vtkSmartPointer<vtkParametricFunctionSource>::New();
  parametric_function->SetParametricFunction(ellipsoid);
  parametric_function->SetUResolution(30);
  parametric_function->SetVResolution(30);
  parametric_function->Update();
  vtkSmartPointer<vtkPolyData> ellipsoid_polydata = parametric_function->GetOutput();

  //ensure polys are a triangle mesh
  vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
  triangleFilter->SetInputData(ellipsoid_polydata);
  triangleFilter->Update();
  vtkSmartPointer<vtkPolyData> triangle_ellipsoid = triangleFilter->GetOutput();


  // Get ellipsoid points into the matrix in order to rotate and translate
  Eigen::MatrixXd ellipsoid_points_matrix(triangle_ellipsoid->GetNumberOfPoints(), 3);
  for(int i = 0; i < triangle_ellipsoid->GetNumberOfPoints(); ++i) {
      double p[3];
      triangle_ellipsoid->GetPoint(i,p);
      ellipsoid_points_matrix(i,0) = p[0];
      ellipsoid_points_matrix(i,1) = p[1];
      ellipsoid_points_matrix(i,2) = p[2];
  }
  //Get ellipsoid faces
  Eigen::MatrixXi ellipsoid_faces(triangle_ellipsoid->GetNumberOfPolys(), 3);
  for(int i = 0; i < triangle_ellipsoid->GetNumberOfPolys(); ++i) {
      vtkNew<vtkIdList> idL;
      triangle_ellipsoid->GetCellPoints(i,idL);
      ellipsoid_faces(i,0) = idL->GetId(0);
      ellipsoid_faces(i,1) = idL->GetId(1);
      ellipsoid_faces(i,2) = idL->GetId(2);
  }



  // rotate the points
  Eigen::MatrixXd rotated_ellipsoid_points = rotation * (ellipsoid_points_matrix.transpose());
  rotated_ellipsoid_points.transposeInPlace(); // n x 3
  // translate the points
  Eigen::MatrixXd translated_points = rotated_ellipsoid_points + cog.replicate(rotated_ellipsoid_points.rows(),1);

  // convert eigen matrix to vtk polydata
  vtkSmartPointer<vtkPoints> best_fitting_ellipsoid_points = vtkSmartPointer<vtkPoints>::New();
  for(int i = 0; i < translated_points.rows(); ++i) {
      double p[3] = {translated_points(i,0), translated_points(i,1), translated_points(i,2)};
      best_fitting_ellipsoid_points->InsertNextPoint(p);
  }
  best_fitting_ellipsoid_polydata->SetPoints(best_fitting_ellipsoid_points);
  best_fitting_ellipsoid_polydata->SetPolys(triangle_ellipsoid->GetPolys());
  best_fitting_ellipsoid_polydata->Modified();

  //compute per-point average error between mesh and best fitting ellipsoid
  vtkSmartPointer<vtkImplicitPolyDataDistance> implicit_distance_filter = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
  implicit_distance_filter->SetInput(best_fitting_ellipsoid_polydata);
  double distance = 0.0;
  for(int i = 0; i < U.rows(); ++i) {
    distance += abs(implicit_distance_filter->EvaluateFunction(U.row(i).data()));
  }
  distance/=U.rows();
  q = distance;
  ell_U = translated_points;
  ell_F = ellipsoid_faces;

 return 0;
};

int srep_init::write_ellipsoid()
{
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  std::string vtk_filename = output_folder + "ellipsoid.vtk";
  writer->SetFileName(vtk_filename.c_str());
  writer->SetInputData(this->best_fitting_ellipsoid_polydata);
  writer->Update();
};

int srep_init::generate_ellipsoid_srep()
{
  const double ELLIPSE_SCALE = 0.9;
  const double EPS = 0.0001;
  double shift = 0.02;
  std::string srep_folder(output_folder);
  srep_folder += "model/";

  //copy final ellipsoid to forward flow folder
  std::string ellSurfaceFile(output_folder);
  ellSurfaceFile += "ellipsoid.vtk";
  std::string  newEllSurfaceFile(output_folder);
  newEllSurfaceFile = newEllSurfaceFile + "forward/" + std::to_string(iter + 1) + ".vtk";
  std::ifstream  src(ellSurfaceFile, std::ios::binary);
  std::ofstream  dst(newEllSurfaceFile,   std::ios::binary);
  dst << src.rdbuf();

  double rz = radii(0);
  double ry = radii(1);
  double rx = radii(2);

  double mrx_o = (rx*rx-rz*rz)/rx;
  double mry_o = (ry*ry-rz*rz)/ry;
  double mrb = mry_o * ELLIPSE_SCALE;
  double mra = mrx_o * ELLIPSE_SCALE;

  // 1. compute the skeletal points
  int nCrestPoints = nRows*2 + (nCols-2)*2;
  double deltaTheta = 2*vtkMath::Pi()/nCrestPoints;
  Eigen::MatrixXd skeletal_points_x(nRows, nCols);
  Eigen::MatrixXd skeletal_points_y(nRows, nCols);
  //MatrixXd skeletal_points_z(nRows, nCols);
  int r = 0, c = 0;
  for(int i = 0; i < nCrestPoints; ++i)
  {
      double theta = vtkMath::Pi() - deltaTheta * floor(nRows/2) - deltaTheta*i;
      double x = mra * cos(theta);
      double y = mrb * sin(theta);

      // these crest points have no inward points (side or corner of the s-rep)
      skeletal_points_x(r, c) = x;
      skeletal_points_y(r, c) = y;
      //skeletal_points_z(r, c) = z;
      if(i < nCols - 1)
      {
          // top row of crest points
          c += 1;
      }
      else if(i < nCols - 1 + nRows - 1)
      {
          // right side col of crest points ( if the top-left point is the origin)
          r = r + 1;
      }
      else if(i < nCols - 1 + nRows - 1 + nCols - 1)
      {
          // bottom row of crest points
          c = c - 1;
      }
      else
      {
          // left side col of crest points
          r = r - 1;
      }
      if((i < nCols - 1 && i > 0) || (i > nCols + nRows - 2 && i < 2*nCols + nRows - 3))
      {
          // compute skeletal points inward
          double mx_ = (mra * mra - mrb * mrb) * cos(theta) / mra; // this is the middle line
          double my_ = .0;
          double dx_ = x - mx_;
          double dy_ = y - my_;
          int numSteps = floor(nRows/2); // steps from crest point to the skeletal point
          double stepSize = 1.0 / double(numSteps); // step size on the half side of srep
          for(int j = 0; j <= numSteps; ++j)
          {
              double tempX_ = mx_ + stepSize * j * dx_;
              double tempY_ = my_ + stepSize * j * dy_;
              if(i < nCols - 1)
              {
                  // step from medial to top at current iteration on the top line
                  int currR = numSteps - j;
                  skeletal_points_x(currR, c-1) = tempX_;
                  skeletal_points_y(currR, c-1) = tempY_;
              }
              else
              {
                  int currR = j + numSteps;
                  skeletal_points_x(currR, c+1) = tempX_;
                  skeletal_points_y(currR, c+1) = tempY_;
              }

          }

      }
  }
  // 2. compute the head points of spokes
  Eigen::MatrixXd skeletal_points(nRows*nCols, 3);
  Eigen::MatrixXd bdry_points_up(nRows*nCols, 3);
  Eigen::MatrixXd bdry_points_down(nRows*nCols, 3);
  Eigen::MatrixXd bdry_points_crest(nCrestPoints, 3);
  Eigen::MatrixXd skeletal_points_crest(nCrestPoints, 3);
  int id_pt = 0; int id_crest = 0;
  Eigen::MatrixXd shift_dir(nCrestPoints, 3); // shift direction for every crest point
  for(int i = 0; i < nRows; ++i)
  {
      for(int j = 0; j < nCols; ++j)
      {
          double mx = skeletal_points_x(i,j);
          double my = skeletal_points_y(i,j);
          double sB = my * mrx_o;
          double cB = mx * mry_o;
          double l = sqrt(sB*sB + cB*cB);
          double sB_n, cB_n;
          if(l == 0)
          {
              sB_n = sB;
              cB_n = cB;
          }
          else
          {
              sB_n = sB / l;
              cB_n = cB / l;
          }
          double cA = l / (mrx_o * mry_o);
          double sA = sqrt(1 - cA*cA);
          double sx = rx * cA * cB_n - mx;
          double sy = ry * cA * sB_n - my;
          double sz = rz * sA;

          double bx = (sx + mx);
          double by = (sy + my);
          double bz = (sz);

          skeletal_points.row(id_pt) << mx, my, 0.0;
          bdry_points_up.row(id_pt) << bx, by, bz;
          bdry_points_down.row(id_pt) << bx, by, -bz;
          id_pt++;
          // fold curve
          if(i == 0 || i == nRows - 1 || j == 0 || j == nCols - 1)
          {
              double cx = rx * cB_n - mx;
              double cy = ry * sB_n - my;
              double cz = 0;
              Eigen::Vector3d v, v2, v3;
              v << cx, cy, cz;
              v2 << sx, sy, 0.0;
              double v_n = v.norm();
              v2.normalize(); // v2 is the unit vector pointing out to norm dir
              v3 = v_n * v2;
              double bx = (v3(0) + mx);
              double by = (v3(1) + my);
              double bz = v3(2);
              bdry_points_crest.row(id_crest) << bx, by, bz;
              skeletal_points_crest.row(id_crest) << mx, my, 0.0;
              shift_dir.row(id_crest) << v2(0), v2(1), v2(2);
              id_crest++;
          }
      }
  }

  // 3. transform the s-rep
  Eigen::MatrixXd transpose_srep = skeletal_points.transpose(); // 3xn
  Eigen::Matrix3d srep_secondMoment = transpose_srep * skeletal_points; // 3x3
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_srep(srep_secondMoment);

  Eigen::Matrix3d rot_srep;
  rot_srep = es_srep.eigenvectors().transpose();
  rotation = rotation * rot_srep;

  // all skeletal points
  Eigen::MatrixXd trans_srep = (rotation * transpose_srep).transpose();
  Eigen::MatrixXd transformed_skeletal_points = trans_srep+
          center.replicate(trans_srep.rows(), 1);

  // up spoke head point on the bdry
  Eigen::MatrixXd transpose_up_pdm = bdry_points_up.transpose();
  Eigen::MatrixXd trans_up_pdm = (rotation * transpose_up_pdm).transpose();
  Eigen::MatrixXd transformed_up_pdm =  trans_up_pdm +
          center.replicate(trans_up_pdm.rows(), 1);

  // down spoke head point on the bdry
  Eigen::MatrixXd transpose_down_pdm = bdry_points_down.transpose();
  Eigen::MatrixXd trans_down_pdm = (rotation * transpose_down_pdm).transpose();
  Eigen::MatrixXd transformed_down_pdm = trans_down_pdm +
          center.replicate(trans_down_pdm.rows(), 1);

  // crest head point on the bdry
  Eigen::MatrixXd transpose_crest_pdm = bdry_points_crest.transpose();
  Eigen::MatrixXd trans_crest_pdm = (rotation * transpose_crest_pdm).transpose();
  Eigen::MatrixXd transformed_crest_pdm = trans_crest_pdm + center.replicate(trans_crest_pdm.rows(), 1);

  // crest base point on the skeletal sheet
  Eigen::MatrixXd transpose_crest_base = skeletal_points_crest.transpose();
  Eigen::MatrixXd trans_crest_base = (rotation * transpose_crest_base).transpose();
  Eigen::MatrixXd transformed_crest_base = trans_crest_base + center.replicate(trans_crest_base.rows(), 1);

  // 4. transfer points to polydata
  // srep_poly is supposed to form a mesh grid connecting skeletal points
  vtkSmartPointer<vtkPolyData>  srep_poly       = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>    skeletal_sheet  = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> skeletal_mesh   = vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkPolyData>  upSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>    upSpokes_pts       = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> upSpokes_lines     = vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkPolyData>  downSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>    downSpokes_pts       = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> downSpokes_lines     = vtkSmartPointer<vtkCellArray>::New();

  // TODO:crest spokes should be a little off the inner spokes
  vtkSmartPointer<vtkPolyData>  crestSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>    crestSpokes_pts       = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> crestSpokes_lines     = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkPolyData> foldCurve_poly         = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>    foldCurve_pts         = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> fold_curve            = vtkSmartPointer<vtkCellArray>::New();

  skeletal_sheet->SetDataTypeToDouble();
  upSpokes_pts->SetDataTypeToDouble();
  downSpokes_pts->SetDataTypeToDouble();
  crestSpokes_pts->SetDataTypeToDouble();

  vtkSmartPointer<vtkDoubleArray> upSpokeLengths = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> downSpokeLengths = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> crestSpokeLengths = vtkSmartPointer<vtkDoubleArray>::New();
  upSpokeLengths->SetNumberOfComponents(1);
  downSpokeLengths->SetNumberOfComponents(1);
  crestSpokeLengths->SetNumberOfComponents(1);

  upSpokeLengths->SetName("spokeLength");
  downSpokeLengths->SetName("spokeLength");
  crestSpokeLengths->SetName("spokeLength");

  vtkSmartPointer<vtkDoubleArray> upSpokeDirs = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> downSpokeDirs = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> crestSpokeDirs = vtkSmartPointer<vtkDoubleArray>::New();

  upSpokeDirs->SetNumberOfComponents(3);
  downSpokeDirs->SetNumberOfComponents(3);
  crestSpokeDirs->SetNumberOfComponents(3);

  upSpokeDirs->SetName("spokeDirection");
  downSpokeDirs->SetName("spokeDirection");
  crestSpokeDirs->SetName("spokeDirection");


  for(int i = 0; i < nRows * nCols; ++i)
  {
      // skeletal points
      double mx = transformed_skeletal_points(i,0);
      double my = transformed_skeletal_points(i,1);
      double mz = transformed_skeletal_points(i,2);
      int id0 = upSpokes_pts->InsertNextPoint(mx,my, mz);

      double bx_up = transformed_up_pdm(i, 0);
      double by_up = transformed_up_pdm(i, 1);
      double bz_up = transformed_up_pdm(i, 2);
      int id1 = upSpokes_pts->InsertNextPoint(bx_up, by_up, bz_up);

      // spoke length and dir
      vtkVector3d upSpoke(bx_up-mx, by_up-my, bz_up-mz);
      double upSpokeLength = upSpoke.Normalize();
      upSpokeLengths->InsertNextTuple1(upSpokeLength);
      upSpokeDirs->InsertNextTuple3(upSpoke.GetX(), upSpoke.GetY(), upSpoke.GetZ());

      // form up spokes
      vtkSmartPointer<vtkLine> up_arrow = vtkSmartPointer<vtkLine>::New();
      up_arrow->GetPointIds()->SetId(0, id0);
      up_arrow->GetPointIds()->SetId(1, id1);
      upSpokes_lines->InsertNextCell(up_arrow);

      // form down spokes
      int id2 = downSpokes_pts->InsertNextPoint(mx, my, mz);
      double bx_down = transformed_down_pdm(i,0);
      double by_down = transformed_down_pdm(i,1);
      double bz_down = transformed_down_pdm(i,2);
      int id3 = downSpokes_pts->InsertNextPoint(bx_down,by_down,bz_down);

      // spoke length and dir
      vtkVector3d downSpoke(bx_down-mx, by_down-my, bz_down-mz);
      double downSpokeLength = downSpoke.Normalize();
      downSpokeLengths->InsertNextTuple1(downSpokeLength);
      downSpokeDirs->InsertNextTuple3(downSpoke.GetX(), downSpoke.GetY(), downSpoke.GetZ());

      vtkSmartPointer<vtkLine> down_arrow = vtkSmartPointer<vtkLine>::New();
      down_arrow->GetPointIds()->SetId(0, id2);
      down_arrow->GetPointIds()->SetId(1, id3);
      downSpokes_lines->InsertNextCell(down_arrow);

  }

  std::string model_prefix(srep_folder);
  model_prefix+= std::to_string(iter);
  if (!vtksys::SystemTools::FileExists(model_prefix, false))
  {
    if (!vtksys::SystemTools::MakeDirectory(model_prefix))
    {
      std::cout << "Failed to create folder : " << model_prefix << std::endl;
    }
  }

  // display up spokes
  upSpokes_poly->SetPoints(upSpokes_pts);
  upSpokes_poly->SetLines(upSpokes_lines);

  upSpokes_poly->GetPointData()->AddArray(upSpokeDirs);
  upSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
  upSpokes_poly->GetPointData()->AddArray(upSpokeLengths);
  upSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");

  // write to file
  vtkSmartPointer<vtkPolyDataWriter> upSpokeWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  std::string upFileName(model_prefix);
  std::string downFileName(model_prefix);
  std::string crestFileName(model_prefix);
  std::string meshFileName(model_prefix);
  std::string curveFileName(model_prefix);


  upFileName  += "/up.vtk";
  upSpokeWriter->SetFileName(upFileName.c_str());
  upSpokeWriter->SetInputData(upSpokes_poly);
  upSpokeWriter->Update();

  // display down spokes
  downSpokes_poly->SetPoints(downSpokes_pts);
  downSpokes_poly->SetLines(downSpokes_lines);

  downSpokes_poly->GetPointData()->AddArray(downSpokeDirs);
  downSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
  downSpokes_poly->GetPointData()->AddArray(downSpokeLengths);
  downSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");

  downFileName +="/down.vtk";
  vtkSmartPointer<vtkPolyDataWriter> downSpokeWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  downSpokeWriter->SetFileName(downFileName.c_str());
  downSpokeWriter->SetInputData(downSpokes_poly);
  downSpokeWriter->Update();

  // deal with skeletal mesh
  for(int i = 0; i < nRows * nCols; ++i)
  {
      double mx = transformed_skeletal_points(i, 0);
      double my = transformed_skeletal_points(i, 1);
      double mz = transformed_skeletal_points(i, 2);
      int current_id = skeletal_sheet->InsertNextPoint(mx, my, mz);
      int current_row = floor(i / nRows);
      int current_col = i - current_row * nRows;
      if(current_col >= 0 && current_row >= 0
              && current_row < nRows-1 && current_col < nCols - 1)
      {
          vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
          quad->GetPointIds()->SetId(0, current_id);
          quad->GetPointIds()->SetId(1, current_id + nCols);
          quad->GetPointIds()->SetId(2, current_id + nCols + 1);
          quad->GetPointIds()->SetId(3, current_id + 1);
          skeletal_mesh->InsertNextCell(quad);
      }
  }
  srep_poly->SetPoints(skeletal_sheet);
  srep_poly->SetPolys(skeletal_mesh);

  meshFileName += "/mesh.vtk";
  vtkSmartPointer<vtkPolyDataWriter> meshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  meshWriter->SetFileName(meshFileName.c_str());
  meshWriter->SetInputData(srep_poly);
  meshWriter->Update();

  // deal with crest spokes
  for(int i = 0; i < nCrestPoints; ++i)
  {
    // tail point
    double cx_t = transformed_crest_base(i, 0);
    double cy_t = transformed_crest_base(i, 1);
    double cz_t = transformed_crest_base(i, 2);
    // head point (_b means boundary)
    double cx_b = transformed_crest_pdm(i, 0);
    double cy_b = transformed_crest_pdm(i, 1);
    double cz_b = transformed_crest_pdm(i, 2);

    if(shift > 0)
    {
      double shift_x = (cx_b - cx_t) * shift;
      double shift_y = (cy_b - cy_t) * shift;
      double shift_z = (cz_b - cz_t) * shift;

      cx_t += shift_x;
      cy_t += shift_y;
      cz_t += shift_z;
    }

    int id0 = crestSpokes_pts->InsertNextPoint(cx_t, cy_t, cz_t);
    int id1 = crestSpokes_pts->InsertNextPoint(cx_b, cy_b, cz_b);

    vtkSmartPointer<vtkLine> crest_arrow = vtkSmartPointer<vtkLine>::New();
    crest_arrow->GetPointIds()->SetId(0, id0);
    crest_arrow->GetPointIds()->SetId(1, id1);
    crestSpokes_lines->InsertNextCell(crest_arrow);

    vtkVector3d crestSpoke(cx_b-cx_t, cy_b-cy_t, cz_b-cz_t);
    double crestSpokeLength = crestSpoke.Normalize();

    crestSpokeLengths->InsertNextTuple1(crestSpokeLength);
    crestSpokeDirs->InsertNextTuple3(crestSpoke.GetX(), crestSpoke.GetY(), crestSpoke.GetZ());
  }
  crestSpokes_poly->SetPoints(crestSpokes_pts);
  crestSpokes_poly->SetLines(crestSpokes_lines);

  crestSpokes_poly->GetPointData()->AddArray(crestSpokeDirs);
  crestSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
  crestSpokes_poly->GetPointData()->AddArray(crestSpokeLengths);
  crestSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");

  crestFileName += "/crest.vtk";
  vtkSmartPointer<vtkPolyDataWriter> crestSpokeWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  crestSpokeWriter->SetFileName(crestFileName.c_str());
  crestSpokeWriter->SetInputData(crestSpokes_poly);
  crestSpokeWriter->Update();

  // deal with fold curve
  for(int i = 0; i < nCrestPoints; ++i)
  {
      double cx_t = transformed_crest_base(i, 0);
      double cy_t = transformed_crest_base(i, 1);
      double cz_t = transformed_crest_base(i, 2);
      double cx_b = transformed_crest_pdm(i, 0);
      double cy_b = transformed_crest_pdm(i, 1);
      double cz_b = transformed_crest_pdm(i, 2);

      if(shift > 0)
      {
          double shift_x = (cx_b - cx_t) * shift;
          double shift_y = (cy_b - cy_t) * shift;
          double shift_z = (cz_b - cz_t) * shift;

          cx_t += shift_x;
          cy_t += shift_y;
          cz_t += shift_z;
      }
      int id0 = foldCurve_pts->InsertNextPoint(cx_t, cy_t, cz_t);

      if(id0 > 0 && i < nCols)
      {
          // first row
          vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
          fold_seg->GetPointIds()->SetId(0, id0-1);
          fold_seg->GetPointIds()->SetId(1, id0);
          fold_curve->InsertNextCell(fold_seg);
      }

      if(i > nCols && i < nCols + 2*(nRows-2) + 1 && (i-nCols) % 2 == 1)
      {
          // right side of crest
          vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
          fold_seg->GetPointIds()->SetId(0, id0-2);
          fold_seg->GetPointIds()->SetId(1, id0);
          fold_curve->InsertNextCell(fold_seg);
      }
      if(i > nCols && i < nCols + 2*(nRows-2) + 1 && (i-nCols) % 2 == 0)
      {
          // part of left side
          vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
          fold_seg->GetPointIds()->SetId(0, id0-2);
          fold_seg->GetPointIds()->SetId(1, id0);
          fold_curve->InsertNextCell(fold_seg);
      }

      if(i == nCols)
      {
          // remaining part of left side
          vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
          fold_seg->GetPointIds()->SetId(0, 0);
          fold_seg->GetPointIds()->SetId(1, id0);
          fold_curve->InsertNextCell(fold_seg);
      }
      if(i > nCols + 2*(nRows-2))
      {
          //bottom side
          vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
          fold_seg->GetPointIds()->SetId(0, id0-1);
          fold_seg->GetPointIds()->SetId(1, id0);
          fold_curve->InsertNextCell(fold_seg);
      }
      if(i == nCrestPoints - 1)
      {
          // bottome right
          vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
          fold_seg->GetPointIds()->SetId(0, id0-nCols);
          fold_seg->GetPointIds()->SetId(1, id0);
          fold_curve->InsertNextCell(fold_seg);
      }
  }
  foldCurve_poly->SetPoints(foldCurve_pts);
  foldCurve_poly->SetLines(fold_curve);

  curveFileName += "/curve.vtk";
  vtkSmartPointer<vtkPolyDataWriter> curveWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  curveWriter->SetFileName(curveFileName.c_str());
  curveWriter->SetInputData(foldCurve_poly);
  curveWriter->Update();
  return 0;
};
