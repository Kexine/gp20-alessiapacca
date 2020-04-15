#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/local_basis.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>


/*** insert any necessary libigl headers here ***/
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/cat.h>
#include <igl/dijkstra.h>
#include <igl/adjacency_matrix.h>

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

bool showingUV = false;
bool freeBoundary = false;
double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;
bool lengthPreservation = false;
bool anglePreservation = false;
bool areaPreservation = false;
Eigen::MatrixXd colors;
bool five = false;
bool cotagentLaplacianMethod = false;
bool notInitialized = true;


void Redraw()
{
	viewer.data().clear();

	if (!showingUV)
	{
		viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(false);

    if(UV.size() != 0)
    {
      viewer.data().set_uv(TextureResolution*UV);
      viewer.data().show_texture = true;
    }
	}
	else
	{
		viewer.data().show_texture = false;
		viewer.data().set_mesh(UV, F);
	}

	if(five){
	    viewer.data().show_texture = true;
	    viewer.data().set_colors(colors);
	}
}

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
	return false;
}

static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}
static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}

void ConvertConstraintsToMatrixForm(VectorXi indices, MatrixXd positions, Eigen::SparseMatrix<double> &C, VectorXd &d)
{
	// Convert the list of fixed indices and their fixed positions to a linear system
	// Hint: The matrix C should contain only one non-zero element per row and d should contain the positions in the correct order.

	//first, resize the matrices
    int numberVertices = V.rows();
    int rowsIndices = indices.rows();
    d.setZero(2 * rowsIndices);
    C.resize(rowsIndices * 2, numberVertices * 2);

    std::vector<Eigen::Triplet<double>> _triplet;
    //fill C with a 1 in those particular indices, where we want to fix the values
    for (int i = 0; i < rowsIndices; i++) {
        _triplet.push_back(Eigen::Triplet<double>(i, indices(i),1));
        _triplet.push_back(Eigen::Triplet<double>(rowsIndices+i,numberVertices+indices(i),1));
    }
    C.setFromTriplets(_triplet.begin(), _triplet.end());

    //fill d
    for (int i = 0; i < rowsIndices; i++)
    {
        d(i) = positions(i, 0);
        d(i + rowsIndices) = positions(i, 1);
    }

}


void computeA_double_half(Eigen::SparseMatrix<double> & A_double, VectorXd & doubleArea, Eigen::SparseMatrix<double> & Dx, Eigen::SparseMatrix<double> & Dy){
    //compute the gradients
    //compute double area
    //copy into A_Resized (#f x #f)
    computeSurfaceGradientMatrix(Dx, Dy);
    igl::doublearea(V, F, doubleArea);
    doubleArea = doubleArea/2;
    for(int i = 0; i < doubleArea.size(); i++) {
        A_double.insert(i, i) = (doubleArea(i));
    }
}

void computeA_double(Eigen::SparseMatrix<double> & A_double, VectorXd & doubleArea, Eigen::SparseMatrix<double> & Dx, Eigen::SparseMatrix<double> & Dy){
    //compute the gradients
    //compute double area
    //copy into A_Resized (#f x #f)
    computeSurfaceGradientMatrix(Dx, Dy);
    igl::doublearea(V, F, doubleArea);
    for(int i = 0; i < doubleArea.size(); i++) {
        A_double.insert(i, i) = (doubleArea(i));
    }
}


//A LSCM
void fill_A_LSCM(Eigen::SparseMatrix<double> & A, const Eigen::SparseMatrix<double> & A_double, const Eigen::SparseMatrix<double> & Dx, const Eigen::SparseMatrix<double> & Dy){
    Eigen::SparseMatrix<double> el1, el2, el3, zeros;
    SparseMatrix<double> res1, res2;
    // see class derivation, fill the matrix A
    el1 = (Dx.transpose() * A_double * Dx) + (Dy.transpose() * A_double * Dy);
    el2 = (-Dx.transpose() * A_double * Dy) + (Dy.transpose() * A_double * Dx);
    el3 = (-Dy.transpose() * A_double * Dx) + (Dx.transpose() * A_double * Dy);

    zeros.resize(el1.rows(), el1.cols());
    zeros.setZero();

    //concatenate matrices
    igl::cat(2, el1, el2, res1);
    igl::cat(2, el3, el1, res2);
    igl::cat(1, res1, res2, A);
}



void checkDeterminant(Matrix2d & U_vi, Matrix2d & Vtr, Matrix2d & sign){
    if(signbit((U_vi * Vtr).determinant() == 1))
        sign << 1, 0, 0, -1;
    else
        sign << 1, 0, 0, 1;
}



void computeRotation(MatrixXd & R, const Eigen::SparseMatrix<double> & Dx, const Eigen::SparseMatrix<double> & Dy){
    for (int i = 0; i < F.rows(); i++)
    {
        Matrix2d U_vi,V_vi,R_vi,S_vi,UV_vi;
        Matrix2d J_vi;
        MatrixXd D_vi(2, V.rows());
        //jacobian matrix construction
        D_vi.row(0) = Dx.row(i);
        D_vi.row(1) = Dy.row(i);
        J_vi = D_vi * UV;
        //SVD with the function SSVD2x2
        SSVD2x2(J_vi, U_vi, S_vi, V_vi);

        Matrix2d sign;
        Matrix2d Vtr = V_vi.transpose();
        checkDeterminant(U_vi, Vtr, sign);
        UV_vi = U_vi * sign * V_vi.transpose();
        R_vi = UV_vi.transpose();

        R(i,0) = R_vi(0,0);
        R(i + F.rows(), 0) = R_vi(0,1);
        R(i + F.rows() * 2, 0) = R_vi(1,0);
        R(i + F.rows() * 3, 0) = R_vi(1,1);
    }
}


void compute_Dijkstra(VectorXi & fixed_UV_indices, MatrixXd & fixed_UV_positions){
    vector<vector<int>> VV;
    igl::adjacency_list(F, VV);
    VectorXd min_distance(V.rows());
    VectorXi previous(V.rows());
    set<int> targets;
    int indexMaxDist;
    int indexMaxGlobalDist1, indexMaxGlobalDist2;
    double maxDist = 0;
    double maxGlobalDist = 0;

    vector<int> boundaryIndices;

    for (int i = 0; i < V.rows(); i++)
    {
        igl::dijkstra(i, targets, VV, min_distance, previous);


        maxDist = min_distance.maxCoeff(&indexMaxDist);
        if (maxDist > maxGlobalDist)
        {
            maxGlobalDist = maxDist;
            indexMaxGlobalDist1 = indexMaxDist;
            indexMaxGlobalDist2 = i;
        }
    }

    fixed_UV_indices.setZero(2);
    fixed_UV_indices << indexMaxGlobalDist2, indexMaxGlobalDist1;

    //cout << "indices with dikstra are: " << indexMaxGlobalDist2 << "and" << indexMaxGlobalDist1 << endl;
    //cout << "max distance is" << maxGlobalDist << endl;
    igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);

    //cout << "positions map to circle: " << fixed_UV_positions << endl;
}


void computeParameterization(int type)
{
	VectorXi fixed_UV_indices;
	MatrixXd fixed_UV_positions;

	SparseMatrix<double> A;
	VectorXd b;
	Eigen::SparseMatrix<double> C;
	VectorXd d;
	// Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices
	if (!freeBoundary)
	{
		// The boundary vertices should be fixed to positions on the unit disc. Find these position and
		// save them in the #V x 2 matrix fixed_UV_position.

		// F is the list of mesh faces
		// fixed_UV_indices is the L ordered list of boundary vertices
        igl::boundary_loop(F, fixed_UV_indices);
        igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
        //cout << fixed_UV_positions;
	}
	else
	{
		// Fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.

        //we could also use the function dijkstra
        //INPUTS: source(index of source vertex), targets (target vector set), VV (adjacency list)
        //OUTPUTS: min_distance (#V by 1 list of the minimum distances from source to all vertices),previous(#V by 1 list of the previous visited vertices)
        fixed_UV_indices.resize(2);
        fixed_UV_positions.resize(2,2);
		compute_Dijkstra(fixed_UV_indices, fixed_UV_positions);
		//we will only have 2 fixed vertices now
    }

	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);

	// Find the linear system for the parameterization (1- Tutte, 2- Harmonic, 3- LSCM, 4- ARAP)
	// and put it in the matrix A.
	// The dimensions of A should be 2#V x 2#V.
	if (type == '1') {
		// Add your code for computing uniform Laplacian for Tutte parameterization
		// Hint: use the adjacency matrix of the mesh

        notInitialized = false;

		// size of b: 2 * number of vertices, indeed we multiply it with A
        b.setZero(2 * V.rows(), 1);

        //just as described in the documentation
        Eigen::SparseMatrix<double> AA;
        SparseVector<double> Asum;
        Eigen::SparseMatrix<double> Adiag;
        SparseMatrix<double> U;

        igl::adjacency_matrix(F, AA);
        // sum each row
        igl::sum(AA,1,Asum);

        // Convert row sums into diagonal of sparse matrix
        igl::diag(Asum,Adiag);
        // Build uniform laplacian
        U = AA-Adiag;
        // we need it for u and for v, so we have to repeat it twice
        igl::repdiag(U, 2, A);
    }

	// ha i constraints
	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~

        notInitialized = false;

        b.setZero(2 * V.rows(), 1);
        Eigen::SparseMatrix<double> L;
        Eigen::SparseMatrix<double> Dx, Dy;
        VectorXd doubleArea;
        SparseMatrix<double> A_double (F.rows(), F.rows());

        //simple call to igl method
        if(!cotagentLaplacianMethod){
           igl::cotmatrix(V, F, L);
        }
        // cot stiffness matrix
        else{
            computeA_double_half(A_double, doubleArea, Dx, Dy);
            L = (Dx.transpose() * A_double * Dx + Dy.transpose() * A_double * Dy);
        }

        //repeat the cot matrix along the diagonal, L, 0, 0, L
        igl::repdiag(L, 2, A);
	}

	if (type == '3') {
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!

        notInitialized = false;
        b.setZero(2 * V.rows(), 1);

		Eigen::SparseMatrix<double> Dx, Dy;
        SparseMatrix<double> A_double (F.rows(), F.rows());
        VectorXd doubleArea;

        computeA_double(A_double, doubleArea, Dx, Dy);
        fill_A_LSCM(A, A_double, Dx, Dy);
	}

	if (type == '4') {
		// Add your code for computing ARAP system and right-hand side
		// Implement a function that computes the local step first
		// Then construct the matrix with the given rotation matrices

		//first, initialize
		if(notInitialized){
                computeParameterization('3');
                notInitialized = false;
		}

        SparseMatrix<double> Dx, Dy;
        SparseMatrix<double> L(V.rows(), V.rows());
        VectorXd doubleArea;
        SparseMatrix<double> A_corsive(F.rows() * 4, V.rows() * 2);
        SparseMatrix<double> A_corsive_tr(F.rows() * 4, V.rows() * 2);
        SparseMatrix<double> zeros(V.rows(), V.rows());
        SparseMatrix<double> result1, result2, result3, result4, result5, result6, result7, result8, result9, result10;
        zeros.setZero();

        //rotation matrix, 2x2 for every face (4x1 vec)
        MatrixXd R(4 * F.rows(), 1);
        computeSurfaceGradientMatrix(Dx, Dy);


        SparseMatrix<double> A_double (F.rows(), F.rows());

        //second, for every face compute the jacobian and find the closest rotation
        computeRotation(R, Dx, Dy);

        //third, minimize the energies.
        // double of the faces areas
        computeA_double_half(A_double, doubleArea, Dx, Dy);

        //A_corsive_transpose * A_corsive
        L = (Dx.transpose() * A_double * Dx) + (Dy.transpose() * A_double * Dy);
        igl::cat(2, L, zeros, result1);
        igl::cat(2, zeros, L, result2);
        igl::cat(1, result1, result2, A);

        //now we want to find C_transpose, so then we can find b

        zeros.resize(F.rows(), V.rows());
        result3 = A_double * Dx;
        result4 = A_double * Dy;
        igl::cat(2, result3, zeros, result5);
        igl::cat(2, result4, zeros, result6);
        igl::cat(2, zeros, result3, result7);
        igl::cat(2, zeros, result4, result8);
        igl::cat(1, result5, result6, result9);
        igl::cat(1, result7, result8, result10);
        igl::cat(1, result9, result10, A_corsive);

        A_corsive_tr = A_corsive.transpose();
        b = A_corsive_tr * R;
	}

	// Solve the linear system.
	// Construct the system as discussed in class and the assignment sheet
	// Use igl::cat to concatenate matrices
	// Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail
	Eigen::SparseMatrix<double> first_block;
    Eigen::SparseLU <Eigen::SparseMatrix<double>> solver;
    Eigen::SparseMatrix<double> result_1;
    Eigen::SparseMatrix<double> result_2;
    Eigen::SparseMatrix<double> Ctranspose = C.transpose();
    Eigen::SparseMatrix<double> zeros_vec(C.rows(), C.rows());
    igl::cat(2, A, Ctranspose, result_1); //[A C']
    igl::cat(2, C, zeros_vec, result_2); //[C 0]
    igl::cat(1, result_1, result_2, first_block);
    // Compute the ordering permutation vector from the structural pattern of A


    VectorXd x_vector, b_vector;
    b_vector.resize(b.rows() + d.rows()); //b, d in column
    //fill b_vector
    b_vector.segment(0, b.rows()) << b;
    b_vector.segment(b.rows(), d.rows()) << d;

    solver.analyzePattern(first_block);
    solver.factorize(first_block);
    x_vector = solver.solve(b_vector);


	// The solver will output a vector
	UV.resize(V.rows(), 2);

	//the first col starts at element 0 and contains the first V.rows elements (u)
	UV.col(0) = x_vector.segment(0, V.rows());
    //the first col starts at element V.rows() and contains the second V.rows elements (v)
	UV.col(1) = x_vector.segment(V.rows(), V.rows());
}

void setColorsMatrix(Eigen::VectorXd & J, const Eigen::VectorXd & ones){
    double eigSmallWhite, eigBigRed;
    Eigen::VectorXd smallVec = Eigen::VectorXd(F.rows());

    colors.setZero();
    eigSmallWhite = J.minCoeff();
    eigBigRed = J.maxCoeff();
    smallVec.setConstant(eigSmallWhite);
    J = J - smallVec;
    J = J / eigBigRed;
    colors.col(0) << ones;
    colors.col(1) << ones - J;
    colors.col(2) << ones - J;
//if J(i) is 0, it means that we have a smallVec   1 1 1  white
// if J(i) is 1, it means that we have a bigVec    1 0 0  red
}


void calculateDistortion(){
    Eigen::Matrix2d identityMatrix = Eigen::MatrixXd::Identity(2,2);
    //color each face
    colors.resize(F.rows(), 3);
    Eigen::VectorXd J = Eigen::VectorXd(F.rows());
    Eigen::VectorXd ones = Eigen::VectorXd(F.rows());
    Eigen::SparseMatrix<double> Dx, Dy;
    ones.setOnes();
    computeSurfaceGradientMatrix(Dx, Dy);
    MatrixXd J_vi(2, 2);
    double distorsion;
    MatrixXd D_vi(2, V.rows());


    //for every face
    for(int i = 0; i < F.rows(); i++){
        D_vi.row(0) = Dx.row(i);
        D_vi.row(1) = Dy.row(i);
        J_vi = D_vi * UV;


        //conformal parametrization LCSM
        if(anglePreservation){
            J_vi = J_vi + J_vi.trace() * identityMatrix;
            distorsion = J_vi.norm();
        }
        //ARAP
        else if(lengthPreservation){
            Matrix2d U_vi,V_vi,R_vi,S_vi,UV_vi;
            //SVD with the function SSVD2x2
            SSVD2x2(J_vi, U_vi, S_vi, V_vi);

            Matrix2d sign;
            Matrix2d Vtr = V_vi.transpose();
            checkDeterminant(U_vi, Vtr, sign);
            UV_vi = U_vi * sign * V_vi.transpose();
            R_vi = UV_vi.transpose();
            J_vi = J_vi - R_vi;
            distorsion = J_vi.norm();
        }
        J(i) = distorsion;
    }
    setColorsMatrix(J, ones);
}

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
	switch (key) {
	case '1':
	case '2':
	case '3':
	case '4':
		computeParameterization(key);
		five = false;
		break;
	case '5':
			// Add your code for detecting and displaying flipped triangles in the
			// UV domain here
			calculateDistortion();
			five = true;
		break;
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
    if(showingUV)
    {
      temp2D = viewer.core;
      viewer.core = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        temp3D = viewer.core;
        viewer.core = temp2D;
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	Redraw();

	return true;
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  Redraw();
  viewer.core.align_camera_center(V);
  showingUV = false;

  return true;
}

bool callback_init(Viewer &viewer)
{
	temp3D = viewer.core;
	temp2D = viewer.core;
	temp2D.orthographic = true;

	return false;
}

int main(int argc,char *argv[]) {
  if(argc != 2) {
    cout << "Usage ex4_bin <mesh.off/obj>" << endl;
    load_mesh("../data/octo_cut2.obj");
  }
  else
  {
    // Read points and normals
    load_mesh(argv[1]);
  }

	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Parmaterization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::Checkbox("Free boundary", &freeBoundary);
			ImGui::Checkbox("Angle Preservation ", &anglePreservation);
            ImGui::Checkbox("Length Preservation ", &lengthPreservation);
            ImGui::Checkbox("Cotangent Laplacian calculation", &cotagentLaplacianMethod);
		}
	};

  viewer.callback_key_pressed = callback_key_pressed;
  viewer.callback_mouse_move = callback_mouse_move;
	viewer.callback_init = callback_init;

  viewer.launch();
}
