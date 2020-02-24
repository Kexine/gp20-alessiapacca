#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
/*** insert any libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/facet_components.h>
#include <igl/jet.h>


#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>



using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V; //V is a Matrix numerox3, dove numero è il numero dei vertici
// Face array, #F x3
Eigen::MatrixXi F; //F is a Matrix numerox3, che contiene per ogni faccia i 3 vertici
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations. Vector of vector of ints (list, 3*#F). Contains the indexes of adjacent faces on each vertex
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to face relations here;
        // store in VF,VFi.
        //VF and VFi have already been declared at the start
        //n = number of vertices, so rows of V
        int n = V.rows(); //remember, rows() gives us the number of the rows of a matrix in eigen
        igl::vertex_triangle_adjacency(n,F,VF, VFi);
        //now I have VF and VFi updated
        for(int i = 0; i < VF.size(); i++) //for every vertex, we print the faces. VF is not a matrix, but a vector of vector, therefore it has no ROWS()
        {
            cout << "The vertex " << i << " is adjacent to the following faces: ";
            for(int j = 0; j < VF[i].size(); j++)
            {
                cout << VF[i][j] << endl;
            }
        }
        //cout << "CIAO";
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V,F);
        // Add your code for computing vertex to vertex relations here:
        // store in VV.

        igl::adjacency_list(F,VV);
        for(int i = 0; i < VV.size(); i++) //for every vertex, we print the vertices between which it is connected
        {
            cout << "The vertex " << i << " is connected to the following vertices: ";
            for(int j = 0; j < VV[i].size(); j++)
            {
                cout << VV[i][j] << endl;
            }
        }
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        FN.setZero(F.rows(),3);
        // Add your code for computing per-face normals here: store in FN.
        igl::per_face_normals(V,F,FN);
        //now we have #Fx3 eigen Matrix of mesh face (triangle) 3D normals (normale per ogni faccia della mesh)
        // Set the viewer normals.
        viewer.data().set_normals(FN); //already written, we want the normals to be in FN then
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-vertex normals here: store in VN.
        igl::per_vertex_normals(V,F,VN);
        //now we have #Fx3 eigen Matrix of normals of the vertices


        // Set the viewer normals.
        viewer.data().set_normals(VN);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-corner normals here: store in CN.

        int threshold = 0;
        igl::per_corner_normals(V,F,threshold,CN);
        //Set the viewer normals
        viewer.data().set_normals(CN);

    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        component_colors_per_face.setZero(F.rows(),3);
        // Add your code for computing per-face connected components here:
        // store the component labels in cid.
        igl::facet_components(F, cid); //cid was pre defined, it is a list of connected component

        int number = cid.maxCoeff() + 1; //indice massimo + 1
        cout << "The number of connected components is: " << number << endl;

        //now we have to count the number of elements for each component. we do that by increasing the size for every component found in cid.
        vector<int> sizeComponent(number,0); //a vector of dimension number [number of components], initialized to zero so then we increase
        for (int i = 0; i < F.rows(); i++) {
            sizeComponent[cid(i)]++; //cid(i) mi da, per ogni input, l'output a cui appartiene (la componente). quindi aumento la casella della componente
        }

        for (int j = 0; j < number; j++) {
            cout << "The size of the component " << j << "is: " << sizeComponent[j] << endl; //here we print the number of faces for each component found.
        }
        
        // Compute colors for the faces based on components, storing them in
        // component_colors_per_face.
        //we use the function jet as suggested in the assignment.
        igl::jet(cid, true, component_colors_per_face);

        // Set the viewer colors
        viewer.data().set_colors(component_colors_per_face);
    }

    if (key == '7') {
        Eigen::MatrixXd Vout=V;
        Eigen::MatrixXi Fout=F;
        // Add your code for sqrt(3) subdivision here.

        // Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
    }

    return true;
}

bool load_mesh(Viewer& viewer,string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
  igl::readOFF(filename,V,F);
  viewer.data().clear();
  viewer.data().set_mesh(V,F);
  viewer.data().compute_normals();
  viewer.core.align_camera_center(V, F);
  return true;
}

int main(int argc, char *argv[]) {
    // Show the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    
    std::string filename;
    if (argc == 2) {
        filename = std::string(argv[1]);
    }
    else {
        filename = std::string("../data/honda.off");
    }
    load_mesh(viewer,filename,V,F);

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
