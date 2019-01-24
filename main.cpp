#include <igl/remove_duplicate_vertices.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include <iostream>
#include <string>

#include "srep_init.h"

int main(int argc, char *argv[])
{


  // Load a mesh in OFF format

  // Init the viewer
  igl::opengl::glfw::Viewer viewer;

  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  // Variables to be defined using the menu
  double dt = 0.1f;
  double smoothAmount = 0.1f;
  int maxIter = 10;
  std::string inputMesh = "../test_data/bunny.off";
  srep_init *init = new srep_init(dt, smoothAmount, maxIter);

  // igl::readOFF(inputMesh, V, F);

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
        init->inputMesh = igl::file_dialog_open();
        if(init->inputMesh.length() != 0)
        {
          Eigen::MatrixXd newV;
          Eigen::MatrixXi newF, SVI, SVJ;
          double eps = 0;
          viewer.load_mesh_from_file(init->inputMesh.c_str());
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
      // Expose variable directly ...
      ImGui::InputDouble("dt", &(init->dt), 0, 0, "%.4f");
      ImGui::InputDouble("smoothAmount", &(init->smoothAmount), 0, 0, "%.4f");
      ImGui::InputInt("maxIter", &(init->maxIter));

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
        init->update_viewer(&viewer);
      }
    }
  };

  // Plot the mesh
  // viewer.data().set_mesh(V, F);
  viewer.launch();
}
