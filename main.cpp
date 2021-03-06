#include <igl/remove_duplicate_vertices.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/readSTL.h>
#include <imgui/imgui.h>

#include <iostream>
#include <string>

#include <cxxopts.hpp>
#include "srep_init.h"

int main(int argc, char *argv[])
{
  try
  {
    cxxopts::Options options("srep_init", "Initializes srep structure from input mesh");

    options.positional_help("file").show_positional_help();
    options.add_options("Configuration")
    ("f,file","Input mesh file", cxxopts::value<std::string>()->default_value(""))
    ("d,dt","Length of time object is allowed to flow at each forward flow step",cxxopts::value<double>()->default_value("0.01"))
    ("s,smoothAmount","Amount to smooth the mesh at each forward flow step",cxxopts::value<double>()->default_value("0.001"))
    ("m,maxIter", "Maximum number of forward flow steps", cxxopts::value<int>()->default_value("10"))
    ("p,samplingDensity", "Sample mesh every arg points",cxxopts::value<int>()->default_value("50"))
    ("e,elasticity","Elasticity parameter of thin plate spline, lower = more elastic",cxxopts::value<double>()->default_value("0.05"))
    ("o,outputFolder","Folder to output the generated sreps and intermediate results. If not specified, current directory will be used",cxxopts::value<std::string>()->default_value(""))
    ("r,nRows","Number of rows in medial mesh of produced srep",cxxopts::value<int>()->default_value("9"))
    ("c,nCols","Number of columns in medial mesh of produced srep",cxxopts::value<int>()->default_value("9"))
    ("t, tolerance", "foward flow will stop when it reaches this average point-to-point error from an ellipsoid, or reaches the maximum number of allowed iterations", cxxopts::value<double>()->default_value("0.01"))
    ("g,gui", "Use gui. All other arguements will be ignored",cxxopts::value<bool>()->default_value("false"))
    ("h,help", "Print help")
    ;

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    bool gui = result["gui"].as<bool>();

    if(result.count("help") || (result.count("file")!=1 && !gui))
    {
      std::cout << options.help({"Configuration"})<<std::endl;
      exit(0);
    }
    srep_init *init = new srep_init(result["file"].as<std::string>(), result["outputFolder"].as<std::string>(), result["nRows"].as<int>(), result["nCols"].as<int>(), result["dt"].as<double>(), result["smoothAmount"].as<double>(), result["tolerance"].as<double>(), result["elasticity"].as<double>(), result["samplingDensity"].as<int>(), result["maxIter"].as<int>());

    if(!gui)
    {
      size_t last_dot = init->input_mesh.rfind('.');
      if (last_dot == std::string::npos)
      {
        std::cerr<<"Error: No file extension found in "<<
        init->input_mesh<<std::endl;
        return false;
      }

      Eigen::MatrixXd V, N;
      Eigen::MatrixXi F;
      Eigen::MatrixXd newV;
      Eigen::MatrixXi newF, SVI, SVJ;
      double eps = 0;

      std::string extension = init->input_mesh.substr(last_dot+1);

      if (extension == "stl" || extension =="STL")
      {
        if (!igl::readSTL(init->input_mesh, V, F, N))
        {
          std::cerr<<"failed to read stl"<<std::endl;
          return false;
        }
      }else if (extension == "off" || extension =="OFF")
      {
        if (!igl::readOFF(init->input_mesh, V, F))
        {
          std::cerr<<"failed to read off"<<std::endl;
          return false;
        }
      }else if (extension == "obj" || extension =="OBJ")
      {
        Eigen::MatrixXd corner_normals;
        Eigen::MatrixXi fNormIndices;

        Eigen::MatrixXd UV_V;
        Eigen::MatrixXi UV_F;

        if (!(igl::readOBJ(init->input_mesh,V, UV_V, corner_normals, F, UV_F, fNormIndices)))
        {
          std::cerr<<"failed to read obj"<<std::endl;
          return false;
        }

      }else
      {
        std::cerr<<"Unrecognized file format"<<std::endl;
        return false;
      }

      igl::remove_duplicate_vertices(V,F,eps,newV,SVI,SVJ,newF);
      init->set_mesh(newV, newF);

      while(init->q > init->tol && init->iter < init->max_iter)
      {
        init->step_forwardflow();
        std::cout<<"Iteration "<<init->iter<<": error = "<<init->q<<std::endl;
      }
      init->write_ellipsoid();
      init->generate_ellipsoid_srep();
      init->backward_flow();

    }else
    {

      // Init the viewer
      igl::opengl::glfw::Viewer viewer;

      // Attach a menu plugin
      igl::opengl::glfw::imgui::ImGuiMenu menu;
      viewer.plugins.push_back(&menu);

      // Variables to be defined using the menu
      double dt = 0.01f;
      double smoothAmount = 0.001f;
      int max_iter = 10;

      // Add content to the default menu window
      menu.callback_draw_viewer_menu = [&]()
      {
        // Draw parent menu content
        // menu.draw_viewer_menu();

        if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen))
        {
          float w = ImGui::GetContentRegionAvailWidth();
          float p = ImGui::GetStyle().FramePadding.x;
          if (ImGui::Button("Load##Mesh", ImVec2((w-p)/2.f, 0)))
          {
            init->input_mesh = igl::file_dialog_open();
            if(init->input_mesh.length() != 0)
            {
              Eigen::MatrixXd newV;
              Eigen::MatrixXi newF, SVI, SVJ;
              double eps = 0;
              size_t last_dot = init->input_mesh.rfind('.');
              if (last_dot == std::string::npos)
              {
                std::cerr<<"Error: No file extension found in "<<
                init->input_mesh<<std::endl;
                return false;
              }

              std::string extension = init->input_mesh.substr(last_dot+1);

              if (extension == "stl" || extension =="STL")
              {
                Eigen::MatrixXd V, N;
                Eigen::MatrixXi F;
                if (!igl::readSTL(init->input_mesh, V, F, N))
                {
                  std::cerr<<"failed to read stl"<<std::endl;
                  return false;
                }
                viewer.data().set_mesh(V,F);
                viewer.data().compute_normals();
                viewer.core.align_camera_center(V,F);
              }
              else
              {
                viewer.load_mesh_from_file(init->input_mesh.c_str());
              }
              igl::remove_duplicate_vertices(viewer.data().V,viewer.data().F,eps,newV,SVI,SVJ,newF);
              init->set_mesh(newV, newF);
              viewer.data().clear();
              viewer.data().set_mesh(newV, newF);
            }
          }

          ImGui::SameLine(0, p);
          if (ImGui::Button("Save##Mesh", ImVec2((w-p)/2.f, 0)))
          {
            viewer.open_dialog_save_mesh();
          }
        }

        // Viewing options
        if (ImGui::CollapsingHeader("Viewing Options"))
        {
          if (ImGui::Button("Center object", ImVec2(-1, 0)))
          {
            viewer.core.align_camera_center(viewer.data().V, viewer.data().F);
          }
          if (ImGui::Button("Snap canonical view", ImVec2(-1, 0)))
          {
            viewer.snap_to_canonical_quaternion();
          }

          // Zoom
          ImGui::PushItemWidth(80 * menu.menu_scaling());
          ImGui::DragFloat("Zoom", &(viewer.core.camera_zoom), 0.05f, 0.1f, 20.0f);

          // Select rotation type
          int rotation_type = static_cast<int>(viewer.core.rotation_type);
          static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
          static bool orthographic = true;
          if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\0002D Mode\0\0"))
          {
            using RT = igl::opengl::ViewerCore::RotationType;
            auto new_type = static_cast<RT>(rotation_type);
            if (new_type != viewer.core.rotation_type)
            {
              if (new_type == RT::ROTATION_TYPE_NO_ROTATION)
              {
                trackball_angle = viewer.core.trackball_angle;
                orthographic = viewer.core.orthographic;
                viewer.core.trackball_angle = Eigen::Quaternionf::Identity();
                viewer.core.orthographic = true;
              }
              else if (viewer.core.rotation_type == RT::ROTATION_TYPE_NO_ROTATION)
              {
                viewer.core.trackball_angle = trackball_angle;
                viewer.core.orthographic = orthographic;
              }
              viewer.core.set_rotation_type(new_type);
            }
          }

          // Orthographic view
          ImGui::Checkbox("Orthographic view", &(viewer.core.orthographic));
          ImGui::PopItemWidth();
        }

        // Draw options
        if (ImGui::CollapsingHeader("Draw Options"))
        {
          if (ImGui::Checkbox("Face-based", &(viewer.data().face_based)))
          {
            viewer.data().dirty =  igl::opengl::MeshGL::DIRTY_ALL;
          }
          ImGui::Checkbox("Show texture", &(viewer.data().show_texture));
          if (ImGui::Checkbox("Invert normals", &(viewer.data().invert_normals)))
          {
            viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_NORMAL;
          }
          ImGui::Checkbox("Show overlay", &(viewer.data().show_overlay));
          ImGui::Checkbox("Show overlay depth", &(viewer.data().show_overlay_depth));
          ImGui::ColorEdit4("Background", viewer.core.background_color.data(),
          ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
          ImGui::ColorEdit4("Line color", viewer.data().line_color.data(),
          ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
          ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.3f);
          ImGui::DragFloat("Shininess", &(viewer.data().shininess), 0.05f, 0.0f, 100.0f);
          ImGui::PopItemWidth();
        }

        // Overlays
        if (ImGui::CollapsingHeader("Overlays"))
        {
          ImGui::Checkbox("Wireframe", &(viewer.data().show_lines));
          ImGui::Checkbox("Fill", &(viewer.data().show_faces));
          ImGui::Checkbox("Show vertex labels", &(viewer.data().show_vertid));
          ImGui::Checkbox("Show faces labels", &(viewer.data().show_faceid));
        }



        // Add new group
        if (ImGui::CollapsingHeader("Forward Flow", ImGuiTreeNodeFlags_DefaultOpen))
        {
          float w = ImGui::GetContentRegionAvailWidth();
          float p = ImGui::GetStyle().FramePadding.x;
          // Expose variable directly ...
          ImGui::InputDouble("dt", &(init->dt), 0, 0, "%.4f");
          ImGui::InputDouble("smoothAmount", &(init->smoothAmount), 0, 0, "%.4f");
          ImGui::InputDouble("Per-vertex best-fitting-ellipsoid tolerance", &(init->tol), 0, 0, "%.4f");
          ImGui::InputDouble("Rigidity of backward flow deformations", &(init->elasticity), 0, 0, "%.4f");
          ImGui::InputInt("max Iterations", &(init->max_iter));
          ImGui::InputInt("TPS sampling density", &(init->sampling_density));

          // Add Step Button
          if (ImGui::Button("Reset Mesh", ImVec2(-1,0)))
          {
            init->U = init->V;
            init->iter = 0;
            init->q = 1;
            init->update_viewer(&viewer);
          }
          if (ImGui::Button("Step Forward Flow Once", ImVec2(-1,0)))
          {
            init->step_forwardflow();
            init->write_ellipsoid();
            init->update_viewer(&viewer);
          }
          if (ImGui::Button("Run Foward Flow", ImVec2(-1,0)))
          {
            while(init->q > init->tol && init->iter < init->max_iter)
            {
              init->step_forwardflow();
              std::cout<<"Iteration "<<init->iter<<": error = "<<init->q<<std::endl;
            }
            init->update_viewer(&viewer);
            init->write_ellipsoid();
          }
          if (ImGui::Button("Run Backward Flow", ImVec2(-1,0)))
          {
            init->generate_ellipsoid_srep();
            init->backward_flow();
          }
          if (ImGui::Button("Show Ellipsoid", ImVec2((w-p)/2.f, 0)))
          {
            viewer.data().clear();
            viewer.data().set_mesh(init->ell_U, init->ell_F);
            viewer.data().compute_normals();
            viewer.core.align_camera_center(init->U,init->F);
          }
          ImGui::SameLine(0, p);
          if (ImGui::Button("Show Object", ImVec2((w-p)/2.f, 0)))
          {
            init->update_viewer(&viewer);
          }
          if (ImGui::Button("Generate Ellipsoid S-Rep", ImVec2(-1,0)))
          {
            init->generate_ellipsoid_srep();
          }
        }
      };

      // Plot the mesh
      // viewer.data().set_mesh(V, F);
      viewer.launch();
    }
  }catch (const cxxopts::OptionException& e)
  {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}
