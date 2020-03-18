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
double stepRate = 0.01;
double scale = 1.3;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;

// Parameter: grid resolution
int resolution = 20;

// Intermediate result: grid points, at which the imlicit function will be evaluated, #G x3
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
int ngx, ngy, ngz;
double minimumDistance = 1000000.0;
Eigen::VectorXd saveConstrValues;
Eigen::MatrixXd closepoints;
string shapeName;


// Functions
void createGrid();
void evaluateImplicitFunc();
void getLines();
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);
void neighbors(Eigen::RowVector3d p);

// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid() {
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines. resize(0, 6);
    grid_values.resize(0);
    V. resize(0, 3);
    F. resize(0, 3);
    FN.resize(0, 3);



    // Grid bounds: axis-aligned bounding box
    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();
    Eigen::RowVector3d s;
    s << 1.1, 1.1, 1.1;

    // Bounding box dimensions
    Eigen::RowVector3d dim = (bb_max - bb_min)*scale;


    // Grid spacing
    const double dx = dim[0] / (double)(resolution - 1);
    const double dy = dim[1] / (double)(resolution - 1);
    const double dz = dim[2] / (double)(resolution - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolution * resolution * resolution, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min * scale + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }
}


void createNewGrid(){
    newGrid.clear();
    //follow the implementation of createGrid,
    Eigen::RowVector3d dim = P.colwise().maxCoeff() - P.colwise().minCoeff();
    //diagonal of the bbox
    double diag = dim.norm();
    //diag of cube
    double step = stepRate * diag;

    //number of steps in x
    int pX = ceil(dim.x() / step);
    ngx = pX;
    //number of steps in y
    int pY = ceil(dim.y() / step);
    ngy = pY;
    //number of steps in z
    int pZ = ceil(dim.z() / step);
    ngz = pY;

    //total size of the grid
    int sizeGrid = pX*pY*pZ;
    int numberPoints = P.rows();
    newGrid.clear();
    newGrid.resize((sizeGrid));
    int sizeCreatedGrid =  newGrid.size();

    for (int i = 0; i < numberPoints; i++) //For every point in there
    {
        // point i
        Eigen::RowVector3d p = P.row(i);
        //find distance between p and min bbox, expressed in number of steps
        Eigen::RowVector3d pDist = (p - P.colwise().minCoeff())/step;

        int X = floor(pDist.x());
        int Y = floor(pDist.y());
        int Z = floor(pDist.z());
        //cout << "\nX of " << i << " = " << X << endl;
        //cout << "\nY of " << i << " = " << Y << endl;
        //cout << "\nZ of " << i << " = " << Z << endl;

        //one long vector that has size pX*pY*pZ. every point will have the indices corresponding in P, that lend in that cube
        int indexI = X + (Y * pX) + (Z * pX * pY);
        //cout << "\nIndex in the newGrid of " << i << " = " << indexI << endl;
        newGrid[indexI].push_back(i);
        //for(int k = 0; k < newGrid[indexI].size(); k++)
        //cout << "at index" << indexI << "we have points: " << newGrid[indexI][k] << endl;
    }
    //cout << "\nyey I created the new grid!!" << endl;
}


int closest_point(Eigen::RowVector3d p){
    Eigen::RowVector3d dim = P.colwise().maxCoeff() - P.colwise().minCoeff();
    double diag = dim.norm();
    double step = stepRate * diag; //diag of cube

    Eigen::RowVector3d pDist = (p - P.colwise().minCoeff()) / step;
    //TODO ask if needs to be positive?
    int X = fabs(floor(pDist.x())); //needs to be positive?
    int Y = fabs(floor(pDist.y()));
    int Z = fabs(floor(pDist.z()));
    //number of steps in x
    int pX = ceil(dim.x() / step);
    //number of steps in y
    int pY = ceil(dim.y() / step);
    //number of steps in z
    int pZ = ceil(dim.z() / step);
    /*cout << "\nX of the point = " << X << endl;
    cout << "\nY of the point = " << Y << endl;
    cout << "\nZ of the point = " << Z << endl;
    cout << "\nSize of the grid X= " << pX << endl;
    cout << "\nSize of the grid Y= " << pY << endl;
    cout << "\nSize of the grid Z= " << pZ << endl;*/

    //in reality we just use 1 for the closest point.
    int checkDistance = floor(wendlandRadius*diag/step) + 1;;

    //FOR THE LOOP we have to make sure we are not accessing wrong indices
    int minBlockOffsetX;
    int minBlockOffsetY;
    int minBlockOffsetZ;
    int maxBlockOffsetX;
    int maxBlockOffsetY;
    int maxBlockOffsetZ;

    minBlockOffsetX = max(0, X-checkDistance);
    maxBlockOffsetX = min(X + 1 + checkDistance, pX);
    minBlockOffsetY = max(0, Y-checkDistance);
    maxBlockOffsetY = min(Y + 1 + checkDistance, pY);
    minBlockOffsetZ = max(0, Z-checkDistance);
    maxBlockOffsetZ = min(Z + 1 + checkDistance, pZ);

    int minDistanceIndex = -1;
    //minimumDistance = diag;
    //here we have to access all the cubes around the one i am, at checkDistance distance
    //for (int i = minBlockOffsetX; i < maxBlockOffsetX; i++){
    for (int i = minBlockOffsetX; i < maxBlockOffsetX; i++){
        //cout << "\ni =   " << i << endl;
        //cout << "\nminBlockOffsetX =   " << minBlockOffsetX << endl;
        //cout << "\nmaxBlockOffsetX =   " << maxBlockOffsetX << endl;
        for(int j = minBlockOffsetY; j < maxBlockOffsetY; j++){
            //cout << "\nj =   " << j << endl;
            //cout << "\nminBlockOffsetY =   " << minBlockOffsetY << endl;
            //cout << "\nmaxBlockOffsetY =   " << maxBlockOffsetY << endl;
            for(int t = minBlockOffsetZ; t < maxBlockOffsetZ; t++){
                //cout << "\nt =   " << t << endl;
                //cout << "\nminBlockOffsetZ =   " << minBlockOffsetZ << endl;
                //cout << "\nmaxBlockOffsetZ =   " << maxBlockOffsetZ << endl;

                int indexI = (i) + (j) * pX + (t) * pY * pX;
                //add neighbours
                for (int u = 0; u < newGrid[indexI].size(); u++) //we access every index inside that grid point
                {
                    //cout << "\n Index point =  " << newGrid[indexI][u] << endl;
                    Eigen::RowVector3d dist = p - P.row(newGrid[indexI][u]);
                    if (dist.norm() < minimumDistance){
                        minDistanceIndex = (newGrid[indexI][u]);
                        minimumDistance = dist.norm();
                        //cout << "\nCurrent minimum distance =  " << minimumDistance << endl;
                    }
                }
            }
        }
    }
    minimumDistance = 1000000000;
    return minDistanceIndex;
}




void neighbors(Eigen::RowVector3d p) {
    saveConstrValues.setZero(constrained_points.rows(), 1);
    Eigen::RowVector3d dim = P.colwise().maxCoeff() - P.colwise().minCoeff();
    double diag = dim.norm();
    double step = stepRate * diag; //diag of cube

    Eigen::RowVector3d pDist = (p - P.colwise().minCoeff()) / step;
    int X = floor(pDist.x());
    int Y = floor(pDist.y());
    int Z = floor(pDist.z());
    int pX = ceil(dim.x() / step); //numero di cubi in x
    int pY = ceil(dim.y() / step); //numero di cubi in y
    int pZ = ceil(dim.z() / step); //numero di cubi in z

    //we want to have a look at the kernel size wendlandRadius
    int checkDistance = floor(wendlandRadius*diag/step) + 1;
    int minBlockOffsetX;
    int minBlockOffsetY;
    int minBlockOffsetZ;
    int maxBlockOffsetX;
    int maxBlockOffsetY;
    int maxBlockOffsetZ;

    minBlockOffsetX = max(0, X-checkDistance);
    maxBlockOffsetX = min(X + checkDistance + 1, pX);
    minBlockOffsetY = max(0, Y-checkDistance);
    maxBlockOffsetY = min(Y + 1 + checkDistance, pY);
    minBlockOffsetZ = max(0, Z-checkDistance);
    maxBlockOffsetZ = min(Z + 1 + checkDistance, pZ);


    //std::vector<int> neighbours;
    Eigen::MatrixXd neighbours;
    neighbours.setZero(constrained_points.rows(), 3);
    int cont = 0;
    int constrained_pointsSize = constrained_points.rows();

    for (int i = minBlockOffsetX; i < maxBlockOffsetX; i++){
        for(int j = minBlockOffsetY; j < maxBlockOffsetY; j++){
            for(int t = minBlockOffsetZ; t < maxBlockOffsetZ; t++){
                int indexI = (i) + (j) * pX + (t) * pY * pX;
                //int indexI = (i) + (j) * ngx + (t) * ngy * ngx;
                for (int u = 0; u < newGrid[indexI].size(); u++) //we access every index inside that grid point
                {
                    double dist = (p - constrained_points.row(newGrid[indexI][u])).norm();
                    //cout << "dist = " << dist << endl;
                    //cout << "newGrid at index " << indexI << " at u " << u << " = " << newGrid[indexI][u] << endl;
                    if (dist <= wendlandRadius*diag){
                        //neighbours.push_back(newGrid[indexI][u]);
                        neighbours(cont) = constrained_points(newGrid[indexI][u]);
                        saveConstrValues.row(cont) = constrained_values.row(newGrid[indexI][u]);
                        cont++;
                    }

                    dist = (p - constrained_points.row(P.rows()+newGrid[indexI][u])).norm();
                    if (dist <= wendlandRadius*diag) {
                        //neighbours.push_back(P.rows()+newGrid[indexI][u]);
                        neighbours(cont) = constrained_points(P.rows()+newGrid[indexI][u]);
                        saveConstrValues.row(cont) = constrained_values.row(P.rows()+newGrid[indexI][u]);
                        cont++;
                    }
                    dist = (p - constrained_points.row(P.rows()*2+newGrid[indexI][u])).norm();
                    if (dist <= wendlandRadius*diag) {
                        //neighbours.push_back(P.rows()*2+newGrid[indexI][u]);
                        neighbours(cont) = constrained_points(P.rows()*2+newGrid[indexI][u]);
                        saveConstrValues.row(cont) = constrained_values.row(P.rows()*2 + newGrid[indexI][u]);
                        cont++;
                    }
                }
            }
        }
    }

    saveConstrValues.conservativeResize(cont, 1);
    neighbours.conservativeResize(cont,3);
    closepoints = neighbours;
}


//returns the wendLand function
double wendLand(double r){
    double result = pow(1 - (r/(wendlandRadius*diag)),(4)) * (4 * r/(wendlandRadius*diag) + 1);
    return result;
}


// Function for explicitly evaluating the implicit function for a sphere of
// radius r centered at c : f(p) = ||p-c|| - r, where p = (x,y,z).
// This will NOT produce valid results for any mesh other than the given
// sphere.
// Replace this with your own function for evaluating the implicit function
// values at the grid points using MLS
void evaluateImplicitFunc() {
    // Sphere center
    auto bb_min = grid_points.colwise().minCoeff().eval();
    auto bb_max = grid_points.colwise().maxCoeff().eval();
    //Eigen::RowVector3d center = 0.5 * (bb_min + bb_max);

    //double radius = 0.5 * (bb_max - bb_min).minCoeff();

    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolution * resolution * resolution);
    grid_values.setZero(resolution * resolution * resolution);

    createNewGrid(); //now we have newgrid
    int size = newGrid.size();
    double radius = wendlandRadius * (bb_max - bb_min).minCoeff();

    int prova= 0;
    // Evaluate sphere's signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);


                Eigen::MatrixXd f;
                Eigen::MatrixXd b;
                Eigen::VectorXd weightVec;

                //cout << grid_points.row(index) <<endl;
                neighbors(grid_points.row(index));
                int closepointsSize = closepoints.rows();


                //if no constraint points are within wendlandRadius, asign a large positive (outside) value to the grid point
                if(closepointsSize == 0)
                    grid_values[index] = MAX;
                else {
                    if(polyDegree == 0) { //a0 . b has size (1)
                        weightVec.setZero(closepointsSize, 1);
                        b.resize(closepointsSize, 1);
                        for (int i = 0; i < closepointsSize; i++) {
                            weightVec(i) = wendLand( (grid_points.row(index) - closepoints.row(i)).norm());
                            b.row(i) << 1;
                        }
                        Eigen::MatrixXd A = weightVec.asDiagonal()*b;
                        int sizeW = weightVec.rows();
                        int sizeSave = saveConstrValues.rows();
                        Eigen::VectorXd c = A.colPivHouseholderQr().solve(weightVec.asDiagonal()* saveConstrValues);
                        //cout << "product" << weightVec.asDiagonal()* saveConstrValues << endl;

                        //cout << c;

                        Eigen::VectorXd finalDot(1);
                        finalDot << 1;
                        grid_values[index] = finalDot.dot(c);

                    }
                    else if(polyDegree == 1) {
                        //a0 + a1x + a2y + a3z (4), b has size 4
                        weightVec.setZero(closepointsSize, 1);
                        for (int i = 0; i < closepointsSize; i++) {
                            weightVec(i) = wendLand((grid_points.row(index) - closepoints.row(i)).norm());
                            b.resize(closepointsSize, 4);
                            f.resize(closepointsSize,1);
                            b.row(i) << 1, closepoints(i,0), closepoints(i,1), closepoints(i,2);
                        }
                        Eigen::MatrixXd A = weightVec.asDiagonal()*b;
                        Eigen::VectorXd c = A.colPivHouseholderQr().solve(weightVec.asDiagonal()*saveConstrValues);
                        Eigen::VectorXd finalDot(4);
                        finalDot << 1, grid_points(index,0), grid_points(index,1), grid_points(index,2);
                        grid_values[index] = finalDot.dot(c);
                    }
                    else if(polyDegree == 2){
                        //a0 + a1x + a2y + a3z + a4x^2 + a5y^2 + a6 z^2 + a7 xy + a8 yz + a9 xz (10), b has size 10
                        weightVec.setZero(closepointsSize, 1);
                        for (int i = 0; i < closepointsSize; i++) {
                            weightVec(i) = wendLand((grid_points.row(index) - closepoints.row(i)).norm());
                            b.resize(closepointsSize, 10);
                            f.resize(closepointsSize,1);
                            b.row(i) << 1, closepoints(i,0), closepoints(i,1), closepoints(i,2), pow(closepoints(i,0),2), pow(closepoints(i,1),2), pow(closepoints(i,2), 2),(closepoints(i,0)) * closepoints(i,1), (closepoints(i,1))*closepoints(i,2), closepoints(i,0)*closepoints(i,2);
                        }
                        Eigen::MatrixXd A = weightVec.asDiagonal()*b;
                        Eigen::VectorXd c = A.colPivHouseholderQr().solve(weightVec.asDiagonal()*saveConstrValues);
                        Eigen::VectorXd finalDot(10);
                        finalDot << 1, grid_points(index,0), grid_points(index,1), grid_points(index,2),
                                pow(grid_points(index,0),2),  pow(grid_points(index,1),2),  pow(grid_points(index,0), 2),
                                grid_points(index,0)* grid_points(index,1),  grid_points(index,1)* grid_points(index,2),  grid_points(index,0)* grid_points(index,2);
                        grid_values(index) = finalDot.dot(c);
                    }

                }

                // Value at (x,y,z) = implicit function for the sphere
                //grid_values[index] = (grid_points.row(index) - center).norm() - radius;
            }
        }
    }

    //cout << "\ndebug print" << endl;
}






// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines() {
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x<resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {
                int index = x + resolution * (y + resolution * z);
                if (x < resolution - 1) {
                    int index1 = (x + 1) + y * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolution - 1) {
                    int index1 = x + (y + 1) * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolution - 1) {
                    int index1 = x + y * resolution + (z + 1) * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}






/*int closest_point(Eigen::RowVector3d p){
    double dmin = 10000000;
    double d;
    int minindex = 0;
    for(int i = 0; i < P.rows(); i++){
        d = (P.row(i) - p).norm();
        if(d < dmin)
            minindex = i;
    }
    return minindex;
}*/




bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        // Show imported points
        Eigen::RowVector3d bb_min, bb_max;
        bb_min = P.colwise().minCoeff();
        bb_max = P.colwise().maxCoeff();
        Eigen::RowVector3d dim = bb_max - bb_min;
        diag = dim.norm();

        viewer.data().clear();
        viewer.core.align_camera_center(P);
        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0,0,0));
    }

    if (key == '2') {
        // Show all constraints
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
        Eigen::RowVector3d bb_min, bb_max;
        bb_min = P.colwise().minCoeff();
        bb_max = P.colwise().maxCoeff();
        // Bounding box dimensions
        Eigen::RowVector3d dim = bb_max - bb_min;
        diag = dim.norm();
        //define epsilon
        double eps = 0.01 * diag;
        //cout << "size of P: " << P.rows()  <<endl;


        //now we have to add the constraints, 2 for every point
        for(int i = 0; i < sizeP; i++){
            //we add the point we are considering
            constrained_points.row(i) = P.row(i);
            constrained_values(i) = 0;



            //the new constrained point is in i+number points, and it's the point + epsilon*n, just as defined in the slides
            eps = stepRate * diag;
            constrained_points.row(i + sizeP) = P.row(i) + eps * N.row(i);
            //we have to check that epsilon is sufficiently small, so that p is the closest point to him
            //int count = 0;
            while (closest_point(P.row(i) + eps * N.row(i)) != i) {
                //halve epsilon and recompute p until it is the case
                eps *= 0.5;
                //cout << "epsilon " << eps <<endl;
                //cout << "calling the distance function for the " << count << "time" <<endl;
                //count++;

            }
            constrained_points.row(i + sizeP) = P.row(i) + eps * N.row(i);
            //cout << "alright, finished" <<endl;
            //add his value that is epsilon
            constrained_values(i + sizeP) = eps;



            eps = stepRate * diag;
            //the second new constrained point is in i+2*number points, and it's the point + epsilon*n, just as defined in the slides
            constrained_points.row(i + (2 * sizeP)) = P.row(i) - eps * N.row(i);
            //we have to check that epsilon is sufficiently small, so that p is the closest point to him
            //count = 0;
            minimumDistance = 1000000.0;
            while (closest_point(P.row(i) - eps * N.row(i)) != i) {
                //halve epsilon and recompute p until it is the case
                eps *= 0.5;
                //cout << "epsilon " << eps <<endl;
                //count++;
            }
            constrained_points.row(i + (2 * sizeP)) = P.row(i) - eps * N.row(i);
            //add his value that is epsilon
            constrained_values(i + 2*sizeP) = -eps;
            eps = stepRate * diag;
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
        viewer.data().point_size = 7;
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
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }


        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolution, resolution, resolution, V, F);
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

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage ex2_bin <mesh.off>" << endl;
        igl::readOFF("../data/cat.off",P,F,N);
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
            ImGui::InputInt("Resolution", &resolution, 0, 0);
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
