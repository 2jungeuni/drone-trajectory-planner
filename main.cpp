#include "Eigen/LU"
#include <iostream>
#include <valarray>
#include <algorithm>
#include <fstream>
#include <cmath>

using namespace std;
using namespace Eigen;

const double deg2rad = M_PI / 180.0;

int sgn(double v) {
    if (v < 0) return -1;
    if (v > 0) return 1;
    return 0;
}

VectorXf get_poly(int order, int der, float t) {
    VectorXf poly(order);
    VectorXf D(order);

    // Initialization
    fill_n(poly.begin(), order, 1);

    int count = 0;
    for (int k = order-1; k >= 0; k--) {
        D[count] = k;
        count ++;
    }

    for (int i = 0; i < order; i++) {
        for (int j = 0; j < der; j++) {
            poly[i] = poly[i] * D[i];
            D[i] = D[i] - 1;
            if (D[i] == -1) {
                D[i] = 0;
            }
        }
    }

    for (int i = 0; i < order; i++) {
        poly[i] = poly[i] * pow(t, D[i]);
    }

    return poly;
}

VectorXf minTraj(VectorXf waypoints, VectorXf times, int order, int time_step) {
    int n = waypoints.size() - 1;

    // Initialize A, and B matrix
    MatrixXf A(order*n, order*n);
    VectorXf B(order*n);

    for (int i = 0; i < order*n; i++) {
        B[i] = 0;
        for (int j = 0; j < order*n; j++) {
            A(i, j) = 0;
        }
    }

    // B matrix
    for (int i = 0; i < n; i++) {
        B[i] = waypoints[i];
        B[i + n] = waypoints[i + 1];
    }

    // Constraint 1 - Starting position for every segment
    VectorXf poly_sp = get_poly(order, 0, 0); // polynomial at staring point
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < poly_sp.size(); j++) {
            A(i, order*i+j) = poly_sp[j];
        }
    }

    // Constraint 2 - Ending position for every segment
    for (int i = 0; i < n; i++) {
        VectorXf poly_ep = get_poly(order, 0, times[i]); // polynomial at ending point
        for (int j = 0; j <poly_ep.size(); j++) {
            A(i+n, order*i+j) = poly_ep[j];
        }
    }

    // Constraint 3 - Starting position derivatives (up to order) are null
    int half_order = order / 2;
    for (int k = 1; k < half_order; k++) {
        VectorXf poly_dev_sp = get_poly(order, k, 0); // polynomial derivatives at starting point
        for (int l = 0; l < order; l++) {
            A(2*n+k-1, l) = poly_dev_sp[l];
        }
    }

    // Constraint 4 - Ending position derivatives (up to order) are null
    for (int k = 1; k < half_order; k++) {
        VectorXf poly_dev_ep = get_poly(order, k, times[time_step-1]);
        for (int l = 0; l < order; l++) {
            A(2*n+(half_order-1)+k-1, order*n-l-1) = poly_dev_ep[order-l-1];
        }
    }

    // Constant 5 - All derivatives are continuous at each waypoint transition
    for (int i = 0; i < n-1; i++) {
        for (int k = 1; k < order-1; k++) {
            VectorXf get_poly_smooth_1 = get_poly(order, k, times[i]);
            VectorXf get_poly_smooth_2 = -get_poly(order, k, 0);
            for (int ii = 0; ii < order; ii++) {
                A(2*n+2*(half_order-1) + i*2*(half_order-1)+k-1, i*order+ii) = get_poly_smooth_1[ii];
            }
            for (int jj = 0; jj < order; jj++) {
                A(2*n+2*(half_order-1) + i*2*(half_order-1)+k-1, i*order+order+jj) = get_poly_smooth_2[jj];
            }
        }
    }

    // solve for the coefficients
    MatrixXf A_inv = A.inverse();

    // Coefficients
    VectorXf coeff = A_inv * B;

    return coeff;
}

VectorXf pos_waypoint_min(float time, int time_step, VectorXf t_wps, int order,
                          VectorXf cX, VectorXf cY, VectorXf cZ) {
    float desPosX = 0;
    float desPosY = 0;
    float desPosZ = 0;

    float desVelX = 0;
    float desVelY = 0;
    float desVelZ = 0;

    float desAccX = 0;
    float desAccY = 0;
    float desAccZ = 0;

    int t_idx;
    for (int i = 0; i < time_step; i++) {
        if (time >= t_wps[i] and time <t_wps[i+1]) {
            t_idx = i;
            break;
        }
        t_idx = time_step;
    }

    // Scaled time (between 0 and duration of segment)
    float scale = (time - t_wps[t_idx]);

    // Which coefficients to use
    int start = order * t_idx;
    int end = order * (t_idx + 1);

    // Set desired position, velocity and acceleration
    VectorXf get_poly_0 = get_poly(order, 0, scale);
    VectorXf get_poly_1 = get_poly(order, 1, scale);
    VectorXf get_poly_2 = get_poly(order, 2, scale);

    VectorXf Px = cX(seq(start, end-1));
    VectorXf Py = cY(seq(start, end-1));
    VectorXf Pz = cZ(seq(start, end-1));

    for (int i = 0; i < order; i++) {
        desPosX += Px[i] * get_poly_0[i];
        desPosY += Py[i] * get_poly_0[i];
        desPosZ += Pz[i] * get_poly_0[i];
    }

    for (int i = 0; i < order; i++) {
        desVelX += Px[i] * get_poly_1[i];
        desVelY += Py[i] * get_poly_1[i];
        desVelZ += Pz[i] * get_poly_1[i];
    }

    for (int i = 0; i < order; i++) {
        desAccX += Px[i] * get_poly_2[i];
        desAccY += Py[i] * get_poly_2[i];
        desAccZ += Pz[i] * get_poly_2[i];
    }

    VectorXf desPosVelAcc(9);
    desPosVelAcc << desPosX, desPosY, desPosZ, desVelX, desVelY, desVelZ, desAccX, desAccY, desAccZ;

    return desPosVelAcc;
}

int main() {
    int time_step = 5;

    VectorXf t(time_step); // Input
    t << 1.875, 5.46535165 - 1.875, 7.96535165 - 5.46535165, 12.01581459 - 7.96535165, 14.35435046 - 12.01581459;

    // Scaled time (between 0 and duration of segment)
    VectorXf t_wps(time_step + 1);
    t_wps[0] = 0;
    for (int i = 1; i < time_step + 1; i++) {
        t_wps[i] = t_wps[i - 1] + t[i - 1];
    }

    //cout << "t wps: \n" << t_wps << endl;

    MatrixXf wp(time_step+1, 3);
    wp << 0, 0, 0, 2, 2, 1, -2, 3, -3, -2, -1, -3, 3, -2, 1, 0, 0, 0;

    //cout << "way points: \n" << wp << endl;

    int order = 8; // Input

    VectorXf way_points_x = wp(seq(0, last), 0); // x-axis way points
    VectorXf way_points_y = wp(seq(0, last), 1); // y-axis way points
    VectorXf way_points_z = wp(seq(0, last), 2); // z-axis way points

    VectorXf coeff_x = minTraj(way_points_x, t, order, time_step);
    VectorXf coeff_y = minTraj(way_points_y, t, order, time_step);
    VectorXf coeff_z = minTraj(way_points_z, t, order, time_step);

    fstream fs;
    string str_buf;
    fs.open("result.csv", ios::out);

    float time = 0;
    float Ts = 0.005;

    while (time < 14) {
        VectorXf desPosVelAcc = pos_waypoint_min(time, time_step, t_wps, order, coeff_x, coeff_y, coeff_z);
        //cout << "Desired state: \n" << desPosVelAcc << endl;
        fs << desPosVelAcc[0] << "," << desPosVelAcc[1] << "," << desPosVelAcc[2] << "\n";
        time += Ts;
    }
    fs.close();

    return 0;
}
