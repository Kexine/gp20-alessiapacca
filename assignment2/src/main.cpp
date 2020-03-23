#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>
#include <cmath>
#include <igl/slice.h>

#define MAX 100000

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
int polyDegree = 0;
double diag;
double stepRate = 0.1;
double scale = 1.2;
double grid_step_x, grid_step_y, grid_step_z;
double offset_view = 1.1;
double epsConst = 0.05;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;
double wendlandRadiusDiag;

// Parameter: grid resolution
int resolutionX = 20;
int resolutionY = 20;
int resolutionZ = 20;
bool N_CONSTRAINT = false;

// Intermediate result: grid points, at which the implicit function will be evaluated, #G x3
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

//new structure for grid
std::vector<std::vector<int>> newGrid;
int size_newGrid_x, size_newGrid_y, size_newGrid_z;
double newGrid_step;
bool PCA = false;

double minimumDistance = 1000000.0;
Eigen::VectorXd saveConstrValues;
Eigen::MatrixXd closepoints;
Eigen::VectorXi neighbors_points;
string shapeName;
Eigen::MatrixXd tempP, tempN;
Eigen::Matrix3d temp;
// Global variable for grid bounds
Eigen::RowVector3d  bb_min, bb_max, dim;
bool singleton_PCA = true;
std::vector<double> distanceVector;

// Functions
void createGrid();
void evaluateImplicitFunc();
void getLines();
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);
void neighbors(Eigen::RowVector3d p);

void init_global_variable(){

    tempP = P;
    tempN = N;

    if(PCA && singleton_PCA) {
        Eigen::MatrixXd centerPoints = P.rowwise() - P.colwise().mean();
        Eigen::MatrixXd cov = centerPoints.adjoint() * centerPoints;
        //cov = cov / (P.rows() - 1);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(cov);
        Eigen::MatrixXd eigenVectors = eigen.eigenvectors();
        eigenVectors = eigenVectors.rightCols(3);
        tempP = (P * eigenVectors);
        tempN = (N * eigenVectors);
        singleton_PCA = false;
    }

    P = tempP;
    N = tempN;

    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();
    dim = (bb_max-bb_min)*scale;
    diag = dim.norm();
    wendlandRadiusDiag = wendlandRadius*diag;
}
int newGrid_coordinates_to_idx(int x, int y, int z)
{
    return x + size_newGrid_x * y + size_newGrid_x * size_newGrid_y * z ;
}
int grid_coordinates_to_idx(int x, int y, int z)
{
    return x + resolutionX*y + resolutionX*resolutionY*z;
}

//returns the wendLand function
Eigen::VectorXd wendLand(std::vector<double> r){
    Eigen::VectorXd dist = Eigen::VectorXd::Map(&r[0], r.size());
    Eigen::VectorXd result = (1 - (dist.array()/(wendlandRadiusDiag))).pow(4) * (4 * dist.array()/(wendlandRadiusDiag) + 1);
    return result;
}

void createGrid() {
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines.resize(0, 6);
    grid_values.resize(0);
    V.resize(0, 3);
    F.resize(0, 3);
    FN.resize(0, 3);


    //Global variable
    grid_step_x = dim[0] / (double)(resolutionX - 1);
    grid_step_y = dim[1] / (double)(resolutionY - 1);
    grid_step_z = dim[2] / (double)(resolutionZ - 1);

    //enlarge the grid
    Eigen::RowVector3d diff(3);
    diff << offset_view * grid_step_x, offset_view * grid_step_y, offset_view * grid_step_z;

    grid_points.resize(resolutionX*resolutionY*resolutionZ, 3);
    for (int i = 0; i < resolutionX; ++i) {
        for (int j = 0; j < resolutionY; ++j) {
            for (int k = 0; k < resolutionZ; ++k) {
                grid_points.row(grid_coordinates_to_idx(i, j, k)) = bb_min -diff + Eigen::RowVector3d(i * grid_step_x, j * grid_step_y, k * grid_step_z);
            }
        }
    }

}

int newGrid_index(Eigen::RowVector3d p){
    Eigen::RowVector3d pDist = (p - bb_min)/newGrid_step;
    int X = floor(pDist.x());
    int Y = floor(pDist.y());
    int Z = floor((pDist.z()));
    return newGrid_coordinates_to_idx(X, Y, Z);
}

void createNewGrid(){
    newGrid.clear();
    newGrid_step = stepRate * diag;

    size_newGrid_x = ceil(dim.x() / newGrid_step);
    size_newGrid_y = ceil(dim.y() / newGrid_step);
    size_newGrid_z = ceil(dim.z() / newGrid_step);
    int newGrid_size = size_newGrid_x * size_newGrid_y * size_newGrid_z;
    newGrid.resize(newGrid_size);

    for (int i = 0; i < P.rows(); ++i) {
        int indexI = newGrid_index(P.row(i));
        newGrid[indexI].push_back(i);
    }
}

std::vector<int> p_to_coordinates_newGrid(Eigen::RowVector3d p){
    Eigen::RowVector3d pDist = (p - bb_min)/newGrid_step;
    int X = floor(pDist.x());
    int Y = floor(pDist.y());
    int Z = floor(pDist.z());
    std::vector<int> result{X,Y,Z};
    return result;
}

int closest_point(Eigen::RowVector3d p){
    minimumDistance = 100000000.0;
    std::vector<int> xyz = p_to_coordinates_newGrid(p);
    int X = xyz[0];
    int Y = xyz[1];
    int Z = xyz[2];
    int checkDistance = 1;
    int minBlockOffsetX = max(0, X-checkDistance);
    int maxBlockOffsetX = min(X + 1 + checkDistance, size_newGrid_x);
    int minBlockOffsetY = max(0, Y-checkDistance);
    int maxBlockOffsetY = min(Y + 1 + checkDistance, size_newGrid_y);
    int minBlockOffsetZ = max(0, Z-checkDistance);
    int maxBlockOffsetZ = min(Z + 1 + checkDistance, size_newGrid_z);
    int minDistanceIndex = -1;

    for (int i = minBlockOffsetX; i < maxBlockOffsetX ; ++i) {
        for (int j = minBlockOffsetY; j < maxBlockOffsetY; ++j) {
            for (int k = minBlockOffsetZ; k < maxBlockOffsetZ ; ++k) {
                int indexI = newGrid_coordinates_to_idx(i, j, k);
                for (int u = 0; u < newGrid[indexI].size(); ++u) {
                    Eigen::RowVector3d dist = (p - P.row(newGrid[indexI][u]));
                    if (dist.norm() < minimumDistance){
                        minDistanceIndex = newGrid[indexI][u];
                        minimumDistance = dist.norm();
                    }
                }
            }
        }
    }
    minimumDistance = 100000000.0;
    return minDistanceIndex;
}

void neighbors(Eigen::RowVector3d p){
    saveConstrValues.setZero(constrained_points.rows(), 1);
    neighbors_points.setZero(constrained_points.rows());
    distanceVector.clear();
    std::vector<int> xyz = p_to_coordinates_newGrid(p);
    int X = xyz[0];
    int Y = xyz[1];
    int Z = xyz[2];
    int checkDistance = floor(wendlandRadiusDiag/newGrid_step)+1;
    int minBlockOffsetX = max(0, X-checkDistance);
    int maxBlockOffsetX = min(X + 1 + checkDistance, size_newGrid_x);
    int minBlockOffsetY = max(0, Y-checkDistance);
    int maxBlockOffsetY = min(Y + 1 + checkDistance, size_newGrid_y);
    int minBlockOffsetZ = max(0, Z-checkDistance);
    int maxBlockOffsetZ = min(Z + 1 + checkDistance, size_newGrid_z);

    int cont = 0;
    for (int i = minBlockOffsetX; i < maxBlockOffsetX ; ++i) {
        for (int j = minBlockOffsetY; j < maxBlockOffsetY; ++j) {
            for (int k = minBlockOffsetZ; k < maxBlockOffsetZ ; ++k) {
                int indexI = newGrid_coordinates_to_idx(i, j, k);
                for (int u = 0; u < newGrid[indexI].size(); ++u) {
                    double dist = (p - constrained_points.row(newGrid[indexI][u])).norm();
                    if (dist < wendlandRadiusDiag){
                        distanceVector.push_back(dist);
                        neighbors_points(cont) = newGrid[indexI][u];
                        saveConstrValues.row(cont) = constrained_values.row(newGrid[indexI][u]);
                        cont++;
                    }
                    dist = (p - constrained_points.row(P.rows() + newGrid[indexI][u])).norm();
                    if (dist < wendlandRadiusDiag){
                        distanceVector.push_back(dist);
                        neighbors_points(cont) = P.rows() + newGrid[indexI][u];
                        saveConstrValues.row(cont) = constrained_values.row(P.rows() + newGrid[indexI][u]);
                        cont++;
                    }
                    dist = (p - constrained_points.row(2*P.rows() + newGrid[indexI][u])).norm();
                    if (dist < wendlandRadiusDiag){
                        distanceVector.push_back(dist);
                        neighbors_points(cont) = 2*P.rows() + newGrid[indexI][u];
                        saveConstrValues.row(cont) = constrained_values.row(2*P.rows() + newGrid[indexI][u]);
                        cont++;
                    }
                }
            }
        }
    }
    neighbors_points.conservativeResize(cont);
    saveConstrValues.conservativeResize(cont, 1);
}

void evaluateImplicitFunc(){

    grid_values.resize(resolutionX*resolutionY*resolutionZ);
    grid_values.setZero(resolutionX*resolutionY*resolutionZ);

    for (int i = 0; i < resolutionX; ++i) {
        for (int j = 0; j < resolutionY; ++j) {
            for (int k = 0; k < resolutionZ; ++k) {
                int grid_index = grid_coordinates_to_idx(i,j,k);

                Eigen::MatrixXd b;
                Eigen::VectorXd weightVec;


                neighbors(grid_points.row(grid_index));

                int closepointsSize = neighbors_points.size();
                if(closepointsSize == 0)
                    grid_values[grid_index] = MAX;
                else{
                    Eigen::VectorXd finalDot(1);
                    weightVec = wendLand(distanceVector);
                    if(polyDegree == 0){
                        b.resize(closepointsSize, 1);
                        for (int i = 0; i < closepointsSize; i++) {
                            b.row(i) << 1;
                        }
                        finalDot << 1;
                    }
                    else if(polyDegree == 1){
                        b.resize(closepointsSize, 4);
                        for (int i = 0; i < closepointsSize; i++) {
                            b.row(i) <<  1, constrained_points(neighbors_points(i), 0), constrained_points(neighbors_points(i), 1),constrained_points(neighbors_points(i), 2);
                        }
                        finalDot.resize(4);
                        finalDot << 1, grid_points(grid_index,0), grid_points(grid_index,1), grid_points(grid_index,2);
                    }
                    else if(polyDegree == 2){
                        b.resize(closepointsSize, 10);
                        for (int i = 0; i < closepointsSize; i++) {
                            b.row(i) << 1, constrained_points(neighbors_points(i), 0), constrained_points(neighbors_points(i), 1), constrained_points(neighbors_points(i), 2),
                                    pow(constrained_points(neighbors_points(i), 0),2), pow(constrained_points(neighbors_points(i), 1),2), pow(constrained_points(neighbors_points(i), 2),2),
                                    constrained_points(neighbors_points(i), 0)*constrained_points(neighbors_points(i), 1), constrained_points(neighbors_points(i), 1)*constrained_points(neighbors_points(i), 2),
                                    constrained_points(neighbors_points(i), 0) *constrained_points(neighbors_points(i), 2);
                        }
                        finalDot.resize(10);
                        finalDot << 1, grid_points(grid_index,0), grid_points(grid_index,1), grid_points(grid_index,2),
                                pow(grid_points(grid_index,0),2),  pow(grid_points(grid_index,1),2),  pow(grid_points(grid_index,2), 2),
                                grid_points(grid_index,0)* grid_points(grid_index,1),  grid_points(grid_index,1)* grid_points(grid_index,2),  grid_points(grid_index,0)* grid_points(grid_index,2);
                    }

                    //TODO CHECK HAT VECTORXI IS OK
                    if(N_CONSTRAINT) {
                        for(int i=0; i<saveConstrValues.size(); i++) {
                            saveConstrValues(i) += (grid_points.row(grid_index) - constrained_points.row(neighbors_points(i))).dot(N.row(neighbors_points(i)%N.rows()));
                        }
                    }

                    Eigen::MatrixXd A = weightVec.asDiagonal()*b;
                    Eigen::VectorXd c = A.colPivHouseholderQr().solve((weightVec).asDiagonal()*saveConstrValues);
                    grid_values[grid_index] = finalDot.dot(c);


                    //Eigen::VectorXd c = (b.transpose() * weightVec.asDiagonal() * weightVec.asDiagonal() * b).ldlt().solve(A.transpose() * weightVec.asDiagonal() * saveConstrValues);
                }
            }
        }
    }
}


// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines() {
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x<resolutionX; ++x) {
        for (unsigned int y = 0; y < resolutionY; ++y) {
            for (unsigned int z = 0; z < resolutionZ; ++z) {
                int index = x + resolutionX * (y + resolutionY * z);
                if (x < resolutionX - 1) {
                    int index1 = (x + 1) + y * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolutionY - 1) {
                    int index1 = x + (y + 1) * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolutionZ - 1) {
                    int index1 = x + y * resolutionX + (z + 1) * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }
    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        // Show imported points
        init_global_variable();

        viewer.data().clear();
        viewer.core.align_camera_center(P);
        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 0));
    }
    if (key == '2') {
        // Show all constraints
        init_global_variable();
        viewer.data().clear();
        viewer.core.align_camera_center(P);
        // Add your code for computing auxiliary constraint points here


        //acc data structure created by me, one array that contains all the points
        createNewGrid();
        // normalize normal of every point
        N.rowwise().normalize();
        int sizeP = P.rows();
        //constrained values contains the value of every constrain, so it is numberPoints*3
        constrained_values.resize(sizeP * 3);
        //constrained points contains the coordinated of the points, so it is numberPoints*3,3
        constrained_points.resize(sizeP * 3, 3);

        //just like in create grid

        //define epsilon
        double eps = epsConst * diag;
        for (int i = 0; i < sizeP; ++i) {
            constrained_points.row(i) = P.row(i);
            constrained_values(i) = 0;
            constrained_points.row(i+sizeP) = P.row(i) + eps*N.row(i);
            while (closest_point(P.row(i) + eps*N.row(i)) != i){
                eps *= 0.5;
            }
            constrained_points.row(i+sizeP) = P.row(i) + eps*N.row(i);
            constrained_values(i+sizeP) = eps;
            eps = epsConst * diag;
            constrained_points.row(i+2*sizeP) = P.row(i) - eps*N.row(i);
            while (closest_point(P.row(i) - eps*N.row(i)) != i){
                eps *= 0.5;
            }
            constrained_points.row(i+2*sizeP) = P.row(i) - eps*N.row(i);
            constrained_values(i+2*sizeP) = -eps;
        }
        // Add code for displaying all points, as above
        viewer.data().clear();
        viewer.core.align_camera_center(constrained_points);
        viewer.data().point_size = 6;
        //blue
        viewer.data().add_points(constrained_points.block(0, 0, sizeP, 3), Eigen::RowVector3d(0,0,1));
        //red
        viewer.data().add_points(constrained_points.block(sizeP, 0, sizeP, 3), Eigen::RowVector3d(1,0,0));
        //green
        viewer.data().add_points(constrained_points.block(sizeP * 2, 0, sizeP, 3), Eigen::RowVector3d(0,1,0));
    }
    if (key == '3') {
        // Show grid points with colored nodes and connected with lines
        if (constrained_points.rows() == 0)
            callback_key_down(viewer, '2', 0);
        viewer.data().clear();
        viewer.core.align_camera_center(P);
        // Add code for creating a grid
        // Add your code for evaluating the implicit function at the grid points
        // Add code for displaying points and lines
        // You can use the following example:

        //cout << constrained_values << endl;

        /*** begin: sphere example, replace (at least partially) with your code ***/
        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc();
        //cout << "ok finito la prima stampa" << saveConstrValues << endl;

        // get grid lines
        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i) {
            double value = grid_values(i);
            if (value < 0) {
                grid_colors(i, 1) = 1;
            }
            else {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }
        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
                                grid_lines.block(0, 3, grid_lines.rows(), 3),
                                Eigen::RowVector3d(0.8, 0.8, 0.8));
        /*** end: sphere example ***/

    }

    if (key == '4') {
        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0)) {
            /* cerr << "Not enough data for Marching Cubes !" << endl;
             return true;*/
            callback_key_down(viewer, '3', 0);
            viewer.data().clear();
        }


        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolutionX, resolutionY, resolutionZ, V, F);
        if (V.rows() == 0) {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }


        igl::per_face_normals(V, F, FN);
        viewer.data().set_mesh(V, F);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);
    }

    //Save
    if (key == '0') {
        igl::writeOFF("../results/res.off", V, F);
        cout << "Saved reconstructed " << shapeName << endl;
    }
    return true;
}

bool callback_load_mesh(Viewer& viewer,string filename)
{
    igl::readOFF(filename,P,F,N);
    callback_key_down(viewer,'1',0);
    return true;
}


int main(int argc, char *argv[]){
    if (argc != 2) {
        cout << "Usage ex2_bin <mesh.off>" << endl;
        igl::readOFF("../data/sphere.off",P,F,N);
    }
    else
    {
        // Read points and normals
        igl::readOFF(argv[1],P,F,N);
    }

    Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    viewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Expose variable directly ...
            ImGui::InputInt("Resolution X", &resolutionX, 0, 0);
            ImGui::InputInt("Resolution Y", &resolutionY, 0, 0);
            ImGui::InputInt("Resolution Z", &resolutionZ, 0, 0);
            ImGui::InputInt("polyDegree", &polyDegree, 0, 0);
            ImGui::InputDouble("Wendland Radius Default", &wendlandRadius, 0, 0);
            ImGui::Checkbox("Use PCA", &PCA);
            ImGui::Checkbox("Use Normal Constrain", &N_CONSTRAINT);
            ImGui::InputDouble("offset view", &offset_view, 0, 0);
            ImGui::InputDouble("Epsilon", &epsConst, 0, 0);

            if (ImGui::Button("Reset Grid", ImVec2(-1,0)))
            {
                std::cout << "ResetGrid\n";
                // Recreate the grid
                createGrid();
                // Switch view to show the grid
                callback_key_down(viewer,'3',0);
            }

            // TODO: Add more parameters to tweak here...
        }

    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}