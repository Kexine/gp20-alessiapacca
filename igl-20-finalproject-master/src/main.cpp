#include <iostream>
#include "Lasso.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>



std::unique_ptr<Lasso> lasso;
using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;
Viewer viewer;

int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}


bool load_mesh(string filename)
{
    igl::read_triangle_mesh(filename,V,F);
    viewer.data().clear();
    viewer.data().set_mesh(V, F);

    viewer.core.align_camera_center(V);
    V_cp = V;
    handle_id.setConstant(V.rows(), 1, -1);
    // Initialize selector
    lasso = std::unique_ptr<Lasso>(new Lasso(V, F, viewer));

    selected_v.resize(0,1);

    return true;
}

bool callback_mouse_down(Viewer& viewer, int button, int modifier)
{
    if (button == (int) Viewer::MouseButton::Right)
        return false;

    down_mouse_x = viewer.current_mouse_x;
    down_mouse_y = viewer.current_mouse_y;

    if (mouse_mode == SELECT)
    {
        if (lasso->strokeAdd(viewer.current_mouse_x, viewer.current_mouse_y) >=0)
            doit = true;
        else
            lasso->strokeReset();
    }
    else if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
    {
        int vi = lasso->pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
        if(vi>=0 && handle_id[vi]>=0)  //if a region was found, mark it for translation/rotation
        {
            moving_handle = handle_id[vi];
            get_new_handle_locations();
            doit = true;
        }
    }
    return doit;
}

bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y)
{
    if (!doit)
        return false;
    if (mouse_mode == SELECT)
    {
        lasso->strokeAdd(mouse_x, mouse_y);
        return true;
    }
    if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
    {
        if (mouse_mode == TRANSLATE) {
            translation = computeTranslation(viewer,
                                             mouse_x,
                                             down_mouse_x,
                                             mouse_y,
                                             down_mouse_y,
                                             handle_centroids.row(moving_handle));
        }
        else {
            rotation = computeRotation(viewer,
                                       mouse_x,
                                       down_mouse_x,
                                       mouse_y,
                                       down_mouse_y,
                                       handle_centroids.row(moving_handle));
        }
        get_new_handle_locations();
#ifndef UPDATE_ONLY_ON_UP
        solve(viewer);
        down_mouse_x = mouse_x;
        down_mouse_y = mouse_y;
#endif
        return true;

    }
    return false;
}

bool callback_mouse_up(Viewer& viewer, int button, int modifier)
{
    if (!doit)
        return false;
    doit = false;
    if (mouse_mode == SELECT)
    {
        selected_v.resize(0,1);
        lasso->strokeFinish(selected_v);
        return true;
    }

    if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
    {
#ifdef UPDATE_ONLY_ON_UP
        if(moving_handle>=0)
      solve(viewer);
#endif
        translation.setZero();
        rotation.setZero(); rotation[3] = 1.;
        moving_handle = -1;

        compute_handle_centroids();

        return true;
    }

    return false;
};

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers)
{
    bool handled = false;

    //TODO fill if needed
    return handled;
}