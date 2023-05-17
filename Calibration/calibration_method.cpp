/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"


using namespace easy3d;

bool Calibration::calibration(
        const std::vector<Vector3D> &points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D> &points_2d, /// input: An array of 2D image points.
        double &fx,    /// output: focal length (i.e., K[0][0], which is equal to 'alpha' in our slides).
        double &fy,    /// output: focal length (i.e., K[1][1], which is equal to 'beta/sin(theta)' in our slides).
        double &cx,    /// output: x component of the principal point (i.e., K[0][2], which is 'u0' in our slides).
        double &cy,    /// output: y component of the principal point (i.e., K[1][2], which is 'v0' in our slides).
        double &skew,  /// output: skew factor (i.e., K[0][1], which is equal to '-alpha * cot(theta)' in our slides).
        Matrix33 &R,   /// output: the 3x3 rotation matrix encoding camera rotation.
        Vector3D &t)   /// output：a 3D vector encoding camera translation.
{
//    std::cout << "\nTODO: I am going to implement the calibration() function in the following file:\n"
//                 "\t    - calibration_method.cpp\n\n";
//
//    std::cout << "[Liangliang]:\n"
//                 "\tCamera calibration requires computing the SVD and inverse of matrices.\n"
//                 "\tIn this assignment, I provide you with a 'Matrix' and a 'Vector' data structures for storing and\n"
//                 "\tmanipulating matrices and vectors of arbitrary sizes. I also wrote some code to show you how to:\n"
//                 "\t    - compute the SVD of a matrix;\n"
//                 "\t    - compute the inverse of a matrix;\n"
//                 "\t    - compute the transpose of a matrix.\n\n"
//                 "\tFeel free to use any of the provided data structures and functions. The commonly used linear algebra\n"
//                 "\tfunctions are provided in the following files:\n"
//                 "\t    - Calibration/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
//                 "\t    - Calibration/vector.h  Vectors of arbitrary dimensions and related functions.\n"
//                 "\t    - Calibration/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
//                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
//                 "\tIn your final submission, please\n"
//                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
//                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
//                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
//                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;
//
//    std::cout << "\n[Liangliang]:\n"
//                 "\tThe input parameters of this function are:\n"
//                 "\t\t- points_3d: An array of 3D points (input to this function)\n"
//                 "\t\t- points_2d: An array of 2D image points (input to this function)\n"
//                 "\tThis function must return either 'true' on success or 'false' otherwise. On success, the camera\n"
//                 "\tparameters are returned by the following variables:\n"
//                 "\t\t- fx and fy: the focal lengths (in our slides, we use 'alpha' and 'beta')\n"
//                 "\t\t- cx and cy: the principal point (in our slides, we use 'u0' and 'v0')\n"
//                 "\t\t- skew:      the skew factor ('-alpha * cot_theta')\n"
//                 "\t\t- R:         the 3x3 rotation matrix encoding camera orientation\n"
//                 "\t\t- t:         a 3D vector encoding camera location.\n"
//                 "\tIMPORTANT: don't forget to write your recovered parameters to the above variables." << std::endl;

    // check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    if (points_3d.size() != points_2d.size() || points_3d.size() < 6) {
        std::cout << "Invalid input!" << std::endl;
        return false;
    }

    // Construct the P matrix and populate it with the points from the input file.
    Matrix P(points_3d.size() * 2, 12, 0.0);
    for (int i = 0; i < points_3d.size(); ++i) {
        P.set_row(i * 2, {points_3d[i].x(), points_3d[i].y(), points_3d[i].z(), 1.0,
                          0.0, 0.0, 0.0, 0.0,
                          -points_2d[i].x() * points_3d[i].x(), -points_2d[i].x() * points_3d[i].y(),
                          -points_2d[i].x() * points_3d[i].z(), -points_2d[i].x()});
        P.set_row(i * 2 + 1, {0.0, 0.0, 0.0, 0.0,
                              points_3d[i].x(), points_3d[i].y(), points_3d[i].z(), 1.0,
                              -points_2d[i].y() * points_3d[i].x(), -points_2d[i].y() * points_3d[i].y(),
                              -points_2d[i].y() * points_3d[i].z(), -points_2d[i].y()});
    }
    std::cout << P << std::endl;

    // Solve for projection matrix M
    Matrix U(points_3d.size() * 2, points_3d.size() * 2, 0.0);   // initialized with 0s
    Matrix S(points_3d.size() * 2, 12, 0.0);   // initialized with 0s
    Matrix V(12, 12, 0.0);   // initialized with 0s

    // Compute the SVD decomposition of A and use it to get the projection matrix.
    svd_decompose(P, U, S, V);
    std::cout << U << std::endl;
    std::cout << S << std::endl;
    std::cout << V << std::endl;
    Vector m = V.get_column(V.cols() - 1);
    Matrix M(3, 4, m.data());
    std::cout << M << std::endl;

    // extract intrinsic parameters from M.
    Vector a1(std::vector<double>{M(0, 0), M(0, 1), M(0, 2)});
    Vector a2(std::vector<double>{M(1, 0), M(1, 1), M(1, 2)});
    Vector a3(std::vector<double>{M(2, 0), M(2, 1), M(2, 2)});
    double b1 = M(0, 3);
    double b2 = M(1, 3);
    double b3 = M(2, 3);

    // calculate the value of rho; sign of rho is checked later since below parameters are the same regardless of sign
    double rho = 1.0 / length(a3);
    cx = pow(rho, 2) * dot(a1, a3);
    cy = pow(rho, 2) * dot(a2, a3);
    double theta = acos(dot(cross(a1, a3), cross(a2, a3)) /
                        (length(cross(a1, a3)) * length(cross(a2, a3))));
    std::cout << "theta:" << theta << std::endl;
    double alpha = pow(rho, 2) * length(cross(a1, a3)) * sin(theta);
    double beta = pow(rho, 2) * length(cross(a2, a3)) * sin(theta);

    // calculate K matrix values
    fx = alpha;
    fy = beta / sin(theta);
    skew = -alpha * (cos(theta) / sin(theta));
    Matrix33 K(fx, skew, cx, 0, fy, cy, 0, 0, 1);

    // extract extrinsic parameters from M.
    Matrix inverseK;
    inverse(K, inverseK);
    Vector b_vector(std::vector<double>{b1, b2, b3});
    t = rho * inverseK * b_vector;

    // check if rho should change signs to negative
    if (t.z() < 0.0) {
        rho = -rho;
        t = rho * inverseK * b_vector;
    }

    // Calculate extrinsic parameters and assign them to matrix
    Vector r1 = cross(a2, a3) / length(cross(a2, a3));
    Vector r3 = rho * a3;
    Vector r2 = cross(r3, r1);

    R.set_row(0, r1);
    R.set_row(1, r2);
    R.set_row(2, r3);

    // Compute the accuracy of the calibration
    std::cout << "Computing accuracy of calibration in pixel coordinates..." << std::endl;
    Matrix34 test_Rt(r1[0], r1[1], r1[2], t.x(),
                    r2[0], r2[1], r2[2], t.y(),
                    r3[0], r3[1], r3[2], t.z());

    double total_error = 0.0;
    for (int i = 0; i < points_3d.size(); ++i){
        Vector P_test(std::vector<double>{points_3d[i].x(), points_3d[i].y(), points_3d[i].z(), 1.0});
        Vector3D p_test = K * (test_Rt * P_test);
        //go from homogeneous coords to image coords
        double x_test = p_test.x() / p_test.z();
        double y_test = p_test.y() / p_test.z();
        double error_x = x_test - points_2d[i].x();
        double error_y = y_test - points_2d[i].y();
        total_error += sqrt(pow(error_x, 2) + pow(error_y, 2));
    }
    double avg_error = total_error / points_3d.size();
    std::cout << "The average error is " << avg_error << " pixels." << std::endl;

    if (avg_error < 3.0){
        std::cout << "Calibration successful" << std::endl;
        return true;
    }

    return false;
}


















