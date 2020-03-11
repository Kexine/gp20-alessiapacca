#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
/*** insert any libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/adjacency_list.h>
#include <igl/edge_topology.h>
#include <igl/barycenter.h>
#include <igl/facet_components.h>
#include <igl/jet.h>


//#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
//#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>




#include <igl/vertex_triangle_adjacency.h>



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
//   For ex 6
//   EV  #Ex2 matrix storing the edge description as pair of indices to vertices
//   FE  #Fx3 matrix storing the Triangle-Edge relation
//   EF  #Ex2 matrix storing the Edge-Triangle relation
Eigen::MatrixXi EV;
Eigen::MatrixXi FE;
Eigen::MatrixXi EF;

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
                cout << VF[i][j]<< ",";
            }
            cout << endl;
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
                cout << VV[i][j] << ",";
            }
            cout << endl;
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

        int threshold = 20;
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
        igl::facet_components(F, cid); //cid was pre defined, it is a list of connected component. so it gives me the component at which every face belongs to.

        int number = cid.maxCoeff() + 1; //maximum index of the component found in cid [so number of connected components, required].
        cout << "The number of connected components is: " << number << endl;

        //now we have to count the number of elements for each component. we do that by increasing the size for every component found in cid.
        vector<int> sizeComponent(number); //a vector of dimension number [number of components], initialized to zero so then we increase
        for(int i = 0; i < number; i++){
            sizeComponent[i] = 0;
        }
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

        //barycenter modifies BC, which is a matrix! we can see it from the code of barycenter
        Eigen::MatrixXd BC;
        igl::barycenter(V, F, BC);
        igl::adjacency_list(F, VV); //VV contains the adjacent vertices of every vertex (so #V rows)
        //now we have to append the newly creaed vertices to Vout, but in rows as specified in the assignment.
        //we put the matrix to zero, and then we merge the matrices

        //2. Modify position of old vertices
        double a;
        int n;
        Eigen::MatrixXd p(V.rows(),V.cols()); //positions for vertices ,    Eigen::MatrixXd P, P.setZero(V.rows(),V.cols());
        p.setZero();

        for (int i = 0; i < V.rows(); i++) { //for every vertex
            n = VV[i].size(); // number of adjacent vertices for every vertex
            a = (4 - 2 * cos(2 * M_PI/double(n))) / 9; //just as defined in the slides
            for (int j = 0; j < n; j++) {
                p.row(i) = p.row(i) + V.row(VV[i][j]); //VV gives us the index of the vertex we want
            }
            p.row(i) = p.row(i) /double(n);
            p.row(i) = a * p.row(i) + (1 - a) * V.row(i);
        }



        //1.1. Update Vout
        Eigen::MatrixXd C(p.rows()+ BC.rows(), p.cols());
        C.setZero(); //we put the matrix to zero, and then we merge the matrices //Vout.resize(P.rows()+BC.rows(), P.cols());Vout << P, BC;
        C << p, BC;  //(https://stackoverflow.com/questions/21496157/eigen-how-to-concatenate-matrix-along-a-specific-dimension)
        Vout = C; //V'' updated

        // we can use both triangle triangle adjacency of edge topology
        // initialize Fout
        Eigen::MatrixXi U(F.rows()*3, F.cols());
        U.setZero();
        Fout = U;

        //   EV  #Ex2 matrix storing the edge description as pair of indices to
        //       vertices (INDICES!!!)
        //   FE  #Fx3 matrix storing the Triangle-Edge relation
        //   EF  #Ex2 matrix storing the Edge-Triangle relation (index of barycenter of face)
        igl::edge_topology(V,F, EV, FE, EF);

        //debug prints
        //cout << BC;
        //cout << endl;
        //cout << BC.rows();
        //cout << endl;
        //cout << BC.cols();
        //cout << endl << "Max coeff of triangle: ";
        //cout << EF.maxCoeff();
        //cout << EF;






        //1.2. Update Fout
        int counter = 0;
        int numberEdges = EV.rows();
        int lastIndexVertices = V.rows();

        //for every edge, we want to add the new face/faces, depending on the fact that we are or not at the border of the mesh
        for(int i = 0; i < numberEdges; i++){
            //the two vertices of the edge
            int tr1 = (EF(i,0));
            int tr2 = (EF(i,1));
            int A = EV(i,0);
            int B = EV(i,1);

            if((tr1 != -1)&&(tr2 != -1)){ //in this case we are not at the border, cause both of the values are different than -1
                //in this case, we have to add 2 triangles for that edge
                //the indices of the barycenters of the two triangles are the indices of the face + V.rows()
                //We "flip" it here by adding the two new faces with the correct new vertices
                int b1 = lastIndexVertices + tr1;
                int b2 = lastIndexVertices + tr2;
                Eigen::Vector3i newFace1(b2, b1, A);//always counter clock wise
                Fout.row(counter) = newFace1;
                counter++;
                Eigen::Vector3i newFace2(b1, b2, B);
                Fout.row(counter) =  newFace2;
            }
            else{ //in this case, we are at the border: we have to add only one new face
                Eigen::Vector3i newFace1;
                if(tr1 == -1){ //then we are missing the left face, we only want to add the right one
                    int b1 = lastIndexVertices + tr2;
                    Eigen::Vector3i newFace1(b1, B, A);
                    Fout.row(counter) = newFace1;
                }
                else if(tr2 == -1){ //then we are missing the right face
                    int b1 = lastIndexVertices + tr1;
                    Eigen::Vector3i newFace1(b1, A, B);
                    Fout.row(counter) = newFace1;
                }
            }
            counter++;
        }





        cout << "faces before: " << F.rows() << endl;
        cout << "faces after: " << Fout.rows();



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

    //igl::opengl::glfw::imgui::ImGuiMenu menu;
    //viewer.plugins.push_back(&menu);

    std::string filename;
    if (argc == 2) {
        filename = std::string(argv[1]);
    }
    else {
        filename = std::string("../data/bumpy_cube.off");
    }
    load_mesh(viewer,filename,V,F);

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
