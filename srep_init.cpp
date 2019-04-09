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
#include <igl/principal_curvature.h>
#include <igl/decimate.h>


#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkMassProperties.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
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
#include <vtkThinPlateSplineTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPointLocator.h>
#include <vtkPointSet.h>
#include <vtkTransformFilter.h>

srep_init::srep_init(double d, double smooth, int max)
{
  this->dt = d;
  this->smoothAmount = smooth;
  this->max_iter = max;
  if (output_folder != "" && !vtksys::SystemTools::FileExists(output_folder, false))
  {
    if (!vtksys::SystemTools::MakeDirectory(output_folder))
    {
      std::cout << "Failed to create folder : " << output_folder << std::endl;
    }
  }

  std::string forward_folder(output_folder);
  forward_folder += "forward/";
  if (!vtksys::SystemTools::FileExists(forward_folder, false))
  {
    if (!vtksys::SystemTools::MakeDirectory(forward_folder))
    {
      std::cout << "Failed to create folder : " << forward_folder << std::endl;
    }
  }
};

srep_init::srep_init(std::string inMesh, std::string outFolder, int nRows, int nCols, double d, double smooth, double tolerance, int samplingDensity, int max)
{
  this->input_mesh = inMesh;
  this->output_folder = outFolder;
  this->nRows = nRows;
  this->nCols = nCols;
  this->dt = d;
  this->smoothAmount = smooth;
  this->tol = tolerance;
  this->sampling_density = samplingDensity;
  this->max_iter = max;
  if (output_folder != "" && !vtksys::SystemTools::FileExists(output_folder, false))
  {
    if (!vtksys::SystemTools::MakeDirectory(output_folder))
    {
      std::cout << "Failed to create folder : " << output_folder << std::endl;
    }
  }

  std::string forward_folder(output_folder);
  forward_folder += "forward/";
  if (!vtksys::SystemTools::FileExists(forward_folder, false))
  {
    if (!vtksys::SystemTools::MakeDirectory(forward_folder))
    {
      std::cout << "Failed to create folder : " << forward_folder << std::endl;
    }
  }
};

int srep_init::set_mesh(Eigen::MatrixXd V, Eigen::MatrixXi F)
{
  // Eigen::VectorXi J;
  // Eigen::MatrixXd W;
  // Eigen::MatrixXi G;
  // igl::decimate(V, F, 2000, W, G, J);
  // this->V = W;
  // this->F = G;
  // this->U = W;
  // igl::cotmatrix(W,G,this->L);

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

  if (iter == 0)
  {
    std::string init_file = output_folder + "forward/0.vtk";
    writer->SetFileName(init_file.c_str());
    writer->SetInputData(polydata);
    writer->Update();
  }

  sprintf(temp,"%d", ++iter);
  std::string prefix = temp;


  // // smoother
  // vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
  // vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  // smoother->SetNumberOfIterations(20);
  // smoother->BoundarySmoothingOff();
  // smoother->FeatureEdgeSmoothingOff();
  // smoother->SetPassBand(smoothAmount);
  // smoother->NonManifoldSmoothingOn();
  // smoother->NormalizeCoordinatesOn();
  //
  //
  // // smooth polydata
  // smoother->SetInputData(polydata);
  // smoother->Update();
  // vtkSmartPointer<vtkPolyData> polydata_smooth = smoother->GetOutput();
  //
  // // set U to smoothed points
  // for(int i = 0; i < U.rows(); ++i) {
  //   double p[3];
  //   polydata_smooth->GetPoint(i,p);
  //   U(i,0) = p[0];
  //   U(i,1) = p[1];
  //   U(i,2) = p[2];
  // }

  std::string vtk_filename = output_folder + "forward/" + prefix + ".vtk";
  writer->SetFileName(vtk_filename.c_str());
  // writer->SetInputData(polydata_smooth);
  writer->SetInputData(polydata);
  writer->Update();

  this->fit_ellipsoid(polydata, floor(sqrt(polydata->GetNumberOfPoints())));
  // this->fit_ellipsoid(polydata_smooth, floor(sqrt(polydata->GetNumberOfPoints())));
  return true;
};

int srep_init::fit_ellipsoid(vtkSmartPointer<vtkPolyData> polydata_smooth, int resolution)
{
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

  for(int i = 0; i < U_temp.rows(); ++i) {
      double p[3];
      cleaned_polydata->GetPoint(i,p);
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
  // resolution should be set to ensure the ellipsoid and the input mesh have ~ the same number of points.
  // this makes landmark sampling work better
  parametric_function->SetUResolution(resolution);
  parametric_function->SetVResolution(resolution);
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

  //compute per-point average error between mesh and nearest point in best fitting ellipsoid
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
  double shift = 0.2;
  std::string srep_folder(output_folder);
  srep_folder += "model/";

  // create folder for srep model. vtksys tool is OS independant
  std::string model_prefix(srep_folder);
  model_prefix+= std::to_string(iter+1);
  if (!vtksys::SystemTools::FileExists(model_prefix, false))
  {
    if (!vtksys::SystemTools::MakeDirectory(model_prefix))
    {
      std::cout << "Failed to create folder : " << model_prefix << std::endl;
    }
  }

  // write header to file
  write_header(iter+1);

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
  // In addition to using these to create the srep, we will save them so we can successively deform them during backflow
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
              if(i==0) id_crest = j;
              else if (j == (nCols - 1)) id_crest = nCols-1 + i;
              else if (i == (nRows - 1)) id_crest = 2*(nCols-1) + nRows - 1 -j;
              else id_crest = 2*(nCols-1) + 2*(nRows - 1) -i;
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
  vtkSmartPointer<vtkPoints>    skeletal_sheet  = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> skeletal_mesh   = vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkPoints>    foldCurve_pts         = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> fold_curve            = vtkSmartPointer<vtkCellArray>::New();

  //unneccesary?
  skeletal_sheet->SetDataTypeToDouble();
  foldCurve_pts->SetDataTypeToDouble();

  vtkSmartPointer<vtkPolyData>  upSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData>  downSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData>  crestSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();

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


  // This loop defines all neccesary data for up & down spokes
  for(int i = 0; i < nRows * nCols; ++i)
  {
      // skeletal points
      double mx = transformed_skeletal_points(i,0);
      double my = transformed_skeletal_points(i,1);
      double mz = transformed_skeletal_points(i,2);

      // up boundary points
      double bx_up = transformed_up_pdm(i, 0);
      double by_up = transformed_up_pdm(i, 1);
      double bz_up = transformed_up_pdm(i, 2);

      // down boundary points
      double bx_down = transformed_down_pdm(i,0);
      double by_down = transformed_down_pdm(i,1);
      double bz_down = transformed_down_pdm(i,2);

      backflow_upPoints->InsertNextPoint(mx, my, mz);
      backflow_upPoints->InsertNextPoint(bx_up, by_up, bz_up);
      backflow_downPoints->InsertNextPoint(mx, my, mz);
      backflow_downPoints->InsertNextPoint(bx_down, by_down, bz_down);


      // up spoke length and dir
      vtkVector3d upSpoke(bx_up-mx, by_up-my, bz_up-mz);
      double upSpokeLength = upSpoke.Normalize();
      upSpokeLengths->InsertNextTuple1(upSpokeLength);
      upSpokeDirs->InsertNextTuple3(upSpoke.GetX(), upSpoke.GetY(), upSpoke.GetZ());

      // down spoke length and dir
      vtkVector3d downSpoke(bx_down-mx, by_down-my, bz_down-mz);
      double downSpokeLength = downSpoke.Normalize();
      downSpokeLengths->InsertNextTuple1(downSpokeLength);
      downSpokeDirs->InsertNextTuple3(downSpoke.GetX(), downSpoke.GetY(), downSpoke.GetZ());

      // add spoke hub to skeletal sheet
      int current_id = skeletal_sheet->InsertNextPoint(mx, my, mz);

      // define quad mesh on spoke hubs
      // TODO: this could be changed to tri-mesh? I think vizualization code would not have to change, but I should verify
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

  // this loop defines all necessary data for crest spokes
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

    backflow_crestPoints->InsertNextPoint(cx_t, cy_t, cz_t);
    backflow_crestPoints->InsertNextPoint(cx_b, cy_b, cz_b);

    int id0 = foldCurve_pts->InsertNextPoint(cx_t, cy_t, cz_t);


    // if(id0 > 0)
    // {
    //     // first row
    //     vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
    //     fold_seg->GetPointIds()->SetId(0, id0-1);
    //     fold_seg->GetPointIds()->SetId(1, id0);
    //     fold_curve->InsertNextCell(fold_seg);
    // }
    // vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
    // fold_seg->GetPointIds()->SetId(0, id0);
    // fold_seg->GetPointIds()->SetId(1, 0);
    // fold_curve->InsertNextCell(fold_seg);

    // if(i > nCols && i < nCols + 2*(nRows-2) + 1 && (i-nCols) % 2 == 1)
    // {
    //     // right side of crest
    //     vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
    //     fold_seg->GetPointIds()->SetId(0, id0-2);
    //     fold_seg->GetPointIds()->SetId(1, id0);
    //     fold_curve->InsertNextCell(fold_seg);
    // }
    // if(i > nCols && i < nCols + 2*(nRows-2) + 1 && (i-nCols) % 2 == 0)
    // {
    //     // part of left side
    //     vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
    //     fold_seg->GetPointIds()->SetId(0, id0-2);
    //     fold_seg->GetPointIds()->SetId(1, id0);
    //     fold_curve->InsertNextCell(fold_seg);
    // }
    //
    // if(i == nCols)
    // {
    //     // remaining part of left side
    //     vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
    //     fold_seg->GetPointIds()->SetId(0, 0);
    //     fold_seg->GetPointIds()->SetId(1, id0);
    //     fold_curve->InsertNextCell(fold_seg);
    // }
    // if(i > nCols + 2*(nRows-2))
    // {
    //     //bottom side
    //     vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
    //     fold_seg->GetPointIds()->SetId(0, id0-1);
    //     fold_seg->GetPointIds()->SetId(1, id0);
    //     fold_curve->InsertNextCell(fold_seg);
    // }
    // if(i == nCrestPoints - 1)
    // {
    //     // bottome right
    //     vtkSmartPointer<vtkLine> fold_seg = vtkSmartPointer<vtkLine>::New();
    //     fold_seg->GetPointIds()->SetId(0, id0-nCols);
    //     fold_seg->GetPointIds()->SetId(1, id0);
    //     fold_curve->InsertNextCell(fold_seg);
    // }

    vtkVector3d crestSpoke(cx_b-cx_t, cy_b-cy_t, cz_b-cz_t);
    double crestSpokeLength = crestSpoke.Normalize();

    crestSpokeLengths->InsertNextTuple1(crestSpokeLength);
    crestSpokeDirs->InsertNextTuple3(crestSpoke.GetX(), crestSpoke.GetY(), crestSpoke.GetZ());
  }

  // Create crest curve
  vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
  polyLine->GetPointIds()->SetNumberOfIds(nCrestPoints+1);
  for (int i = 0; i < nCrestPoints; i++)
  {
    polyLine->GetPointIds()->SetId(i, i);
  }
  polyLine->GetPointIds()->SetId(nCrestPoints,0);
  fold_curve->InsertNextCell(polyLine);

  std::string upFileName(model_prefix);
  upFileName  += "/up.vtp";
  std::string downFileName(model_prefix);
  downFileName +="/down.vtp";
  std::string crestFileName(model_prefix);
  crestFileName += "/crest.vtp";


  // Add hubs, connectivity, spoke lengths, and directions to the upspoke poly
  upSpokes_poly->SetPoints(skeletal_sheet);
  upSpokes_poly->SetPolys(skeletal_mesh);

  upSpokes_poly->GetPointData()->AddArray(upSpokeDirs);
  upSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
  upSpokes_poly->GetPointData()->AddArray(upSpokeLengths);
  upSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");

  // then write to file
  vtkSmartPointer<vtkXMLPolyDataWriter> upSpokeWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  upSpokeWriter->SetDataModeToAscii();
  upSpokeWriter->SetFileName(upFileName.c_str());
  upSpokeWriter->SetInputData(upSpokes_poly);
  upSpokeWriter->Update();

  // Add hubs, connectivity, spoke lengths, and directions to the downspoke poly
  downSpokes_poly->SetPoints(skeletal_sheet);
  downSpokes_poly->SetPolys(skeletal_mesh);

  downSpokes_poly->GetPointData()->AddArray(downSpokeDirs);
  downSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
  downSpokes_poly->GetPointData()->AddArray(downSpokeLengths);
  downSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");

  // then write to file
  vtkSmartPointer<vtkXMLPolyDataWriter> downSpokeWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  downSpokeWriter->SetDataModeToAscii();
  downSpokeWriter->SetFileName(downFileName.c_str());
  downSpokeWriter->SetInputData(downSpokes_poly);
  downSpokeWriter->Update();

  // Add hubs, connectivity, spoke lengths, and directions to the crest spoke poly
  crestSpokes_poly->SetPoints(foldCurve_pts);
  crestSpokes_poly->SetLines(fold_curve); // Visualizer must have lines instead of polys here, which is just so silly. Could be changed easily.

  crestSpokes_poly->GetPointData()->AddArray(crestSpokeDirs);
  crestSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
  crestSpokes_poly->GetPointData()->AddArray(crestSpokeLengths);
  crestSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");

  // then write to file
  vtkSmartPointer<vtkXMLPolyDataWriter> crestSpokeWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  crestSpokeWriter->SetDataModeToAscii();
  crestSpokeWriter->SetFileName(crestFileName.c_str());
  crestSpokeWriter->SetInputData(crestSpokes_poly);
  crestSpokeWriter->Update();

  // save srep polydata for backflow
  this->up_mesh->DeepCopy(skeletal_mesh);
  this->down_mesh->DeepCopy(skeletal_mesh);
  this->crest_mesh->DeepCopy(fold_curve);
  return 0;
};

int srep_init::write_header(int prefix)
{
  std::stringstream output;

  output<<"<s-rep>"<<std::endl;
  output<<"  <nRows>"<<nRows<<"</nRows>"<<std::endl;
  output<<"  <nCols>"<<nCols<<"</nCols>"<<std::endl;
  output<<"  <color>"<<std::endl;
  output<<"    <red>0</red>"<<std::endl;
  output<<"    <green>0.5</green>"<<std::endl;
  output<<"    <blue>0</blue>"<<std::endl;
  output<<"  </color>"<<std::endl;
  output<<"  <isMean>False</isMean>"<<std::endl;
  output<<"  <meanStatPath/>"<<std::endl;
  output<<"  <upSpoke>up.vtp</upSpoke>"<<std::endl;
  output<<"  <downSpoke>down.vtp</downSpoke>"<<std::endl;
  output<<"  <crestSpoke>crest.vtp</crestSpoke>"<<std::endl;
  output<<"</s-rep>"<<std::endl;

  std::string header_file(output_folder);
  header_file = header_file + "model/" + std::to_string(prefix) + "/header.xml";
  std::ofstream out_file;
  out_file.open(header_file);
  out_file << output.rdbuf();
  out_file.close();
}

int srep_init::write_srep(int prefix)
{
  std::string srep_folder(output_folder);
  srep_folder += "model/";

  // create folder for srep model. vtksys tool is OS independant
  std::string model_prefix(srep_folder);
  model_prefix+= std::to_string(prefix);
  if (!vtksys::SystemTools::FileExists(model_prefix, false))
  {
    if (!vtksys::SystemTools::MakeDirectory(model_prefix))
    {
      std::cout << "Failed to create folder : " << model_prefix << std::endl;
    }
  }

  write_header(prefix);

  std::string upFileName(model_prefix);
  upFileName  += "/up.vtp";
  std::string downFileName(model_prefix);
  downFileName +="/down.vtp";
  std::string crestFileName(model_prefix);
  crestFileName += "/crest.vtp";

  vtkSmartPointer<vtkPoints> old_upPoints = this->backflow_upPoints;
  vtkSmartPointer<vtkPoints> old_downPoints = this->backflow_downPoints;
  vtkSmartPointer<vtkPoints> old_crestPoints = this->backflow_crestPoints;

  vtkSmartPointer<vtkPolyData>  upSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData>  downSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData>  crestSpokes_poly      = vtkSmartPointer<vtkPolyData>::New();

  vtkSmartPointer<vtkPoints> upPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> downPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> crestPoints = vtkSmartPointer<vtkPoints>::New();

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
      // up spoke hubs
      double m_up[3];
      old_upPoints->GetPoint(2*i,m_up);
      // up spoke surface points
      double b_up[3];
      old_upPoints->GetPoint(2*i +1, b_up);

      // down spoke hubs
      double m_down[3];
      old_downPoints->GetPoint(2*i,m_down);
      // down spoke surface points
      double b_down[3];
      old_downPoints->GetPoint(2*i +1, b_down);

      upPoints->InsertNextPoint(m_up);
      downPoints->InsertNextPoint(m_down);

      // up spoke length and dir
      vtkVector3d upSpoke(b_up[0]-m_up[0], b_up[1]-m_up[1], b_up[2]-m_up[2]);
      double upSpokeLength = upSpoke.Normalize();
      upSpokeLengths->InsertNextTuple1(upSpokeLength);
      upSpokeDirs->InsertNextTuple3(upSpoke.GetX(), upSpoke.GetY(), upSpoke.GetZ());

      // down spoke length and dir
      vtkVector3d downSpoke(b_down[0]-m_down[0], b_down[1]-m_down[1], b_down[2]-m_down[2]);
      double downSpokeLength = downSpoke.Normalize();
      downSpokeLengths->InsertNextTuple1(downSpokeLength);
      downSpokeDirs->InsertNextTuple3(downSpoke.GetX(), downSpoke.GetY(), downSpoke.GetZ());

  }

  for (int i = 0; i < old_crestPoints->GetNumberOfPoints()/2; i++)
  {
    // crest spoke hubs
    double m_crest[3];
    old_crestPoints->GetPoint(2*i,m_crest);
    // crest spoke surface points
    double b_crest[3];
    old_crestPoints->GetPoint(2*i +1, b_crest);

    crestPoints->InsertNextPoint(m_crest);

    // crest spoke length and dir
    vtkVector3d crestSpoke(b_crest[0]-m_crest[0], b_crest[1]-m_crest[1], b_crest[2]-m_crest[2]);
    double crestSpokeLength = crestSpoke.Normalize();
    crestSpokeLengths->InsertNextTuple1(crestSpokeLength);
    crestSpokeDirs->InsertNextTuple3(crestSpoke.GetX(), crestSpoke.GetY(), crestSpoke.GetZ());
  }

  // Add hubs, connectivity, spoke lengths, and directions to the upspoke poly
  upSpokes_poly->SetPolys(this->up_mesh);
  upSpokes_poly->SetPoints(upPoints);

  upSpokes_poly->GetPointData()->AddArray(upSpokeDirs);
  upSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
  upSpokes_poly->GetPointData()->AddArray(upSpokeLengths);
  upSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");

  // then write to file
  vtkSmartPointer<vtkXMLPolyDataWriter> upSpokeWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  upSpokeWriter->SetDataModeToAscii();
  upSpokeWriter->SetFileName(upFileName.c_str());
  upSpokeWriter->SetInputData(upSpokes_poly);
  upSpokeWriter->Update();

  // Add hubs, connectivity, spoke lengths, and directions to the downspoke poly
  downSpokes_poly->SetPolys(this->down_mesh);
  downSpokes_poly->SetPoints(downPoints);

  downSpokes_poly->GetPointData()->AddArray(downSpokeDirs);
  downSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
  downSpokes_poly->GetPointData()->AddArray(downSpokeLengths);
  downSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");

  // then write to file
  vtkSmartPointer<vtkXMLPolyDataWriter> downSpokeWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  downSpokeWriter->SetDataModeToAscii();
  downSpokeWriter->SetFileName(downFileName.c_str());
  downSpokeWriter->SetInputData(downSpokes_poly);
  downSpokeWriter->Update();

  // Add hubs, connectivity, spoke lengths, and directions to the crest spoke poly
  crestSpokes_poly->SetLines(this->crest_mesh); // Visualizer must have lines instead of polys here, which is just so silly. Could be changed easily.
  crestSpokes_poly->SetPoints(crestPoints);

  crestSpokes_poly->GetPointData()->AddArray(crestSpokeDirs);
  crestSpokes_poly->GetPointData()->SetActiveVectors("spokeDirection");
  crestSpokes_poly->GetPointData()->AddArray(crestSpokeLengths);
  crestSpokes_poly->GetPointData()->SetActiveScalars("spokeLength");

  // then write to file
  vtkSmartPointer<vtkXMLPolyDataWriter> crestSpokeWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  crestSpokeWriter->SetDataModeToAscii();
  crestSpokeWriter->SetFileName(crestFileName.c_str());
  crestSpokeWriter->SetInputData(crestSpokes_poly);
  crestSpokeWriter->Update();
}

int srep_init::backward_flow()
{
  //Folder to output flowed mesh at each step for sanity check
  std::string intermediate_steps_folder(output_folder);
  intermediate_steps_folder += "backward/";
  if (!vtksys::SystemTools::FileExists(intermediate_steps_folder, false))
  {
    if (!vtksys::SystemTools::MakeDirectory(intermediate_steps_folder))
    {
      std::cout << "Failed to create folder : " << intermediate_steps_folder << std::endl;
    }
  }


  vtkSmartPointer<vtkPolyDataReader> targetSurfaceReader = vtkSmartPointer<vtkPolyDataReader>::New();


  std::string firstMeshFile(output_folder);
  firstMeshFile = firstMeshFile + "forward/" + std::to_string(iter) + ".vtk";
  targetSurfaceReader->SetFileName(firstMeshFile.c_str());
  targetSurfaceReader->Update();

  vtkSmartPointer<vtkPolyData> polyData_target = targetSurfaceReader->GetOutput();
  vtkSmartPointer<vtkPolyData> polyData_source = this->best_fitting_ellipsoid_polydata;

  vtkSmartPointer<vtkPoints> source_landmarks = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> target_landmarks = vtkSmartPointer<vtkPoints>::New();

  // First deal with the best_fitting_ellipsoid -> last deformed object transformation.
  // This is a special case because we do not have point to point correspondance
  // between the generated ellipsoid mesh and the deformed object.
  // I use only surface landmarks, not spoke hub landmarks here. I wonder if it makes much of a difference.
  vtkSmartPointer<vtkPointLocator> closestPointFinder = vtkSmartPointer<vtkPointLocator>::New();
  closestPointFinder->SetDataSet(polyData_source);
  closestPointFinder->AutomaticOn();
  closestPointFinder->SetNumberOfPointsPerBucket(2);
  closestPointFinder->BuildLocator();

  for(int i = 0; i < polyData_target->GetNumberOfPoints(); i+=sampling_density)
  {
    double p[3];
    double s[3];
    polyData_target->GetPoint(i, p);
    target_landmarks->InsertNextPoint(p);
    polyData_source->GetPoint(closestPointFinder->FindClosestPoint(p),s);
    source_landmarks->InsertNextPoint(s);
  }

  vtkSmartPointer<vtkThinPlateSplineTransform> tps = vtkSmartPointer<vtkThinPlateSplineTransform>::New();
  tps->SetSourceLandmarks(source_landmarks);
  tps->SetTargetLandmarks(target_landmarks);
  tps->SetSigma(this->elasticity);
  tps->SetBasisToR();

  // transform ellipsoid mesh
  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformFilter->SetInputData(polyData_source);
  transformFilter->SetTransform(tps);
  transformFilter->Update();

  //We will successively deform the ellipsoid mesh, so we can see how the backflow performs on the whole surface.
  vtkSmartPointer<vtkPolyData> deformed_ellipsoid = transformFilter->GetOutput();


  std::string deformed_mesh_file(intermediate_steps_folder);
  deformed_mesh_file = deformed_mesh_file + std::to_string(iter) + ".vtk";
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName(deformed_mesh_file.c_str());
  writer->SetInputConnection(transformFilter->GetOutputPort());
  writer->Update();

  // now transform srep polydata
  vtkSmartPointer<vtkTransformPolyDataFilter> srep_transformFilter =
  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  srep_transformFilter->SetTransform(tps);

  vtkSmartPointer<vtkPolyData> temp_upPoly = vtkSmartPointer<vtkPolyData>::New();
  temp_upPoly->SetPoints(this->backflow_upPoints);

  srep_transformFilter->SetInputData(temp_upPoly);
  srep_transformFilter->Update();
  this->backflow_upPoints->DeepCopy(srep_transformFilter->GetOutput()->GetPoints());

  vtkSmartPointer<vtkPolyData> temp_downPoly = vtkSmartPointer<vtkPolyData>::New();
  temp_downPoly->SetPoints(this->backflow_downPoints);

  srep_transformFilter->SetInputData(temp_downPoly);
  srep_transformFilter->Update();
  this->backflow_downPoints->DeepCopy(srep_transformFilter->GetOutput()->GetPoints());

  vtkSmartPointer<vtkPolyData> temp_crestPoly = vtkSmartPointer<vtkPolyData>::New();
  temp_crestPoly->SetPoints(this->backflow_crestPoints);

  srep_transformFilter->SetInputData(temp_crestPoly);
  srep_transformFilter->Update();
  this->backflow_crestPoints->DeepCopy(srep_transformFilter->GetOutput()->GetPoints());

  // then write to file
  this->write_srep(iter);

  // Now deal with the rest of the steps
  // We use both surface landmarks and spoke hub landmarks

  //first, the old target becomes the new source
  polyData_source = polyData_target;
  source_landmarks = vtkSmartPointer<vtkPoints>::New();

  // Eigen matrices of surface mesh points for curvature calculations
  // Topology does not change, so we can use the same face matrix (F) as we have been all along
  Eigen::MatrixXd source_V(polyData_source->GetNumberOfPoints(),3);
  // Per-vertex normals
  Eigen::MatrixXd source_N;

  for(int i = 0; i < source_V.rows(); ++i) {
    double p[3];
    polyData_source->GetPoint(i,p);
    source_V(i,0) = p[0];
    source_V(i,1) = p[1];
    source_V(i,2) = p[2];
  }
  igl::per_vertex_normals(source_V,this->F, source_N);

  // now get principal curvatures for each point via quadratic fitting
  // (more robust than descrete curvature calculation)
  // and compute mean curvature
  Eigen::VectorXd PV1,PV2;   // Curvature
  Eigen::MatrixXd PD1,PD2;   // Curvature direction
  Eigen::VectorXd H;         // Mean Curvature
  igl::principal_curvature(source_V,this->F,PD1,PD2,PV1,PV2);
  H = 0.5*(PV1+PV2);


  for(unsigned int i = 0; i < polyData_source->GetNumberOfPoints(); i += sampling_density)
  {
    //resample the surface mesh
    double p[3];
    polyData_source->GetPoint(i,p);
    source_landmarks->InsertNextPoint(p);

  }

  int crestPoint_indices[backflow_crestPoints->GetNumberOfPoints()/2];

  closestPointFinder = vtkSmartPointer<vtkPointLocator>::New();
  closestPointFinder->SetDataSet(polyData_source);
  closestPointFinder->AutomaticOn();
  closestPointFinder->SetNumberOfPointsPerBucket(2);
  closestPointFinder->BuildLocator();
  for(int i = 0; i < backflow_crestPoints->GetNumberOfPoints()/2; i++)
  {
    double p[3];
    double pt[3];
    //Get crest boundary point
    backflow_crestPoints->GetPoint(2*i+1,p);
    //find nearest point on mesh, store index for use in later iterations
    crestPoint_indices[i] = closestPointFinder->FindClosestPoint(p);
    polyData_source->GetPoint(crestPoint_indices[i],pt);

    //then compute where the crest spoke hubs should be
    Eigen::Vector3d v(pt[0], pt[1], pt[2]);
    Eigen::Vector3d s;
    s = v + (1/H(crestPoint_indices[i]))*source_N.row(crestPoint_indices[i]).transpose();
    double sp[3];
    sp[0] = s[0];
    sp[1] = s[1];
    sp[2] = s[2];

    source_landmarks->InsertNextPoint(sp);
  }



  for (int stepNum = iter-1; stepNum >= 0; stepNum-- )
  {
    std::cout<<"Backflow step: "<<stepNum<<std::endl;
    vtkSmartPointer<vtkPoints> target_landmarks = vtkSmartPointer<vtkPoints>::New();

    std::string nextMeshFile(output_folder);
    nextMeshFile = nextMeshFile + "forward/" + std::to_string(stepNum) + ".vtk";
    targetSurfaceReader->SetFileName(nextMeshFile.c_str());
    targetSurfaceReader->Update();

    // next surface mesh to flow to
    vtkSmartPointer<vtkPolyData> polyData_target = targetSurfaceReader->GetOutput();

    Eigen::MatrixXd target_V(polyData_target->GetNumberOfPoints(),3);
    // Per-vertex normals
    Eigen::MatrixXd target_N;

    for(int i = 0; i < target_V.rows(); ++i) {
      double p[3];
      polyData_target->GetPoint(i,p);
      target_V(i,0) = p[0];
      target_V(i,1) = p[1];
      target_V(i,2) = p[2];
    }
    igl::per_vertex_normals(target_V,this->F, target_N);

    igl::principal_curvature(target_V,this->F,PD1,PD2,PV1,PV2);
    // mean curvature
    H = 0.5*(PV1+PV2);

    //Number of points in source and target should be equal

    for(unsigned int i = 0; i < polyData_target->GetNumberOfPoints(); i += sampling_density)
    {
        double p[3];
        polyData_target->GetPoint(i,p);
        target_landmarks->InsertNextPoint(p);

    }

    for (int i = 0; i < backflow_crestPoints->GetNumberOfPoints()/2; i++)
    {
      double pt[3];
      polyData_target->GetPoint(crestPoint_indices[i],pt);

      //then compute where the crest spoke hubs should be
      Eigen::Vector3d v(pt[0], pt[1], pt[2]);
      Eigen::Vector3d s;
      s = v + (1/H(crestPoint_indices[i]))*source_N.row(crestPoint_indices[i]).transpose();
      double sp[3];
      sp[0] = s[0];
      sp[1] = s[1];
      sp[2] = s[2];

      target_landmarks->InsertNextPoint(sp);
    }

    vtkSmartPointer<vtkThinPlateSplineTransform> tps = vtkSmartPointer<vtkThinPlateSplineTransform>::New();
    tps->SetSourceLandmarks(source_landmarks);
    tps->SetTargetLandmarks(target_landmarks);
    tps->SetSigma(0.05);
    tps->SetBasisToR();

    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter->SetInputData(deformed_ellipsoid);
    transformFilter->SetTransform(tps);
    transformFilter->Update();

    deformed_ellipsoid->DeepCopy(transformFilter->GetOutput());

    std::string deformed_mesh_file(intermediate_steps_folder);
    deformed_mesh_file = deformed_mesh_file + std::to_string(stepNum) + ".vtk";
    writer->SetFileName(deformed_mesh_file.c_str());
    writer->SetInputData(deformed_ellipsoid);
    writer->Update();


    // now transform srep polydata
    vtkSmartPointer<vtkTransformPolyDataFilter> srep_transformFilter =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    srep_transformFilter->SetTransform(tps);


    vtkSmartPointer<vtkPolyData> temp_upPoly = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> temp_upPoints = vtkSmartPointer<vtkPoints>::New();
    temp_upPoints->DeepCopy(this->backflow_upPoints);

    temp_upPoly->SetPoints(temp_upPoints);

    srep_transformFilter->SetInputData(temp_upPoly);
    srep_transformFilter->Update();
    // this->backflow_upPoints = vtkSmartPointer<vtkPoints>::New();
    this->backflow_upPoints->DeepCopy(srep_transformFilter->GetOutput()->GetPoints());

    vtkSmartPointer<vtkPolyData> temp_downPoly = vtkSmartPointer<vtkPolyData>::New();
    temp_downPoly->SetPoints(this->backflow_downPoints);

    srep_transformFilter->SetInputData(temp_downPoly);
    srep_transformFilter->Update();
    this->backflow_downPoints->DeepCopy(srep_transformFilter->GetOutput()->GetPoints());

    vtkSmartPointer<vtkPolyData> temp_crestPoly = vtkSmartPointer<vtkPolyData>::New();
    temp_crestPoly->SetPoints(this->backflow_crestPoints);

    srep_transformFilter->SetInputData(temp_crestPoly);
    srep_transformFilter->Update();
    this->backflow_crestPoints->DeepCopy(srep_transformFilter->GetOutput()->GetPoints());

    // then write to file
    this->write_srep(stepNum);

    //Deforming the source landmark points for the next iteration
    vtkSmartPointer<vtkTransformPolyDataFilter> landmark_transformFilter =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    landmark_transformFilter->SetTransform(tps);

    vtkSmartPointer<vtkPoints> temp_points = vtkSmartPointer<vtkPoints>::New();
    temp_points->DeepCopy(source_landmarks);
    vtkSmartPointer<vtkPolyData> landmark_poly = vtkSmartPointer<vtkPolyData>::New();
    landmark_poly->SetPoints(temp_points);

    landmark_transformFilter->SetInputData(landmark_poly);
    landmark_transformFilter->Update();
    source_landmarks->DeepCopy(landmark_transformFilter->GetOutput()->GetPoints());
  }

  return 0;
}
