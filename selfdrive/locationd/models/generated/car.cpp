#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1362453500626495482) {
   out_1362453500626495482[0] = delta_x[0] + nom_x[0];
   out_1362453500626495482[1] = delta_x[1] + nom_x[1];
   out_1362453500626495482[2] = delta_x[2] + nom_x[2];
   out_1362453500626495482[3] = delta_x[3] + nom_x[3];
   out_1362453500626495482[4] = delta_x[4] + nom_x[4];
   out_1362453500626495482[5] = delta_x[5] + nom_x[5];
   out_1362453500626495482[6] = delta_x[6] + nom_x[6];
   out_1362453500626495482[7] = delta_x[7] + nom_x[7];
   out_1362453500626495482[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1829696213295155707) {
   out_1829696213295155707[0] = -nom_x[0] + true_x[0];
   out_1829696213295155707[1] = -nom_x[1] + true_x[1];
   out_1829696213295155707[2] = -nom_x[2] + true_x[2];
   out_1829696213295155707[3] = -nom_x[3] + true_x[3];
   out_1829696213295155707[4] = -nom_x[4] + true_x[4];
   out_1829696213295155707[5] = -nom_x[5] + true_x[5];
   out_1829696213295155707[6] = -nom_x[6] + true_x[6];
   out_1829696213295155707[7] = -nom_x[7] + true_x[7];
   out_1829696213295155707[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4589523013307382378) {
   out_4589523013307382378[0] = 1.0;
   out_4589523013307382378[1] = 0;
   out_4589523013307382378[2] = 0;
   out_4589523013307382378[3] = 0;
   out_4589523013307382378[4] = 0;
   out_4589523013307382378[5] = 0;
   out_4589523013307382378[6] = 0;
   out_4589523013307382378[7] = 0;
   out_4589523013307382378[8] = 0;
   out_4589523013307382378[9] = 0;
   out_4589523013307382378[10] = 1.0;
   out_4589523013307382378[11] = 0;
   out_4589523013307382378[12] = 0;
   out_4589523013307382378[13] = 0;
   out_4589523013307382378[14] = 0;
   out_4589523013307382378[15] = 0;
   out_4589523013307382378[16] = 0;
   out_4589523013307382378[17] = 0;
   out_4589523013307382378[18] = 0;
   out_4589523013307382378[19] = 0;
   out_4589523013307382378[20] = 1.0;
   out_4589523013307382378[21] = 0;
   out_4589523013307382378[22] = 0;
   out_4589523013307382378[23] = 0;
   out_4589523013307382378[24] = 0;
   out_4589523013307382378[25] = 0;
   out_4589523013307382378[26] = 0;
   out_4589523013307382378[27] = 0;
   out_4589523013307382378[28] = 0;
   out_4589523013307382378[29] = 0;
   out_4589523013307382378[30] = 1.0;
   out_4589523013307382378[31] = 0;
   out_4589523013307382378[32] = 0;
   out_4589523013307382378[33] = 0;
   out_4589523013307382378[34] = 0;
   out_4589523013307382378[35] = 0;
   out_4589523013307382378[36] = 0;
   out_4589523013307382378[37] = 0;
   out_4589523013307382378[38] = 0;
   out_4589523013307382378[39] = 0;
   out_4589523013307382378[40] = 1.0;
   out_4589523013307382378[41] = 0;
   out_4589523013307382378[42] = 0;
   out_4589523013307382378[43] = 0;
   out_4589523013307382378[44] = 0;
   out_4589523013307382378[45] = 0;
   out_4589523013307382378[46] = 0;
   out_4589523013307382378[47] = 0;
   out_4589523013307382378[48] = 0;
   out_4589523013307382378[49] = 0;
   out_4589523013307382378[50] = 1.0;
   out_4589523013307382378[51] = 0;
   out_4589523013307382378[52] = 0;
   out_4589523013307382378[53] = 0;
   out_4589523013307382378[54] = 0;
   out_4589523013307382378[55] = 0;
   out_4589523013307382378[56] = 0;
   out_4589523013307382378[57] = 0;
   out_4589523013307382378[58] = 0;
   out_4589523013307382378[59] = 0;
   out_4589523013307382378[60] = 1.0;
   out_4589523013307382378[61] = 0;
   out_4589523013307382378[62] = 0;
   out_4589523013307382378[63] = 0;
   out_4589523013307382378[64] = 0;
   out_4589523013307382378[65] = 0;
   out_4589523013307382378[66] = 0;
   out_4589523013307382378[67] = 0;
   out_4589523013307382378[68] = 0;
   out_4589523013307382378[69] = 0;
   out_4589523013307382378[70] = 1.0;
   out_4589523013307382378[71] = 0;
   out_4589523013307382378[72] = 0;
   out_4589523013307382378[73] = 0;
   out_4589523013307382378[74] = 0;
   out_4589523013307382378[75] = 0;
   out_4589523013307382378[76] = 0;
   out_4589523013307382378[77] = 0;
   out_4589523013307382378[78] = 0;
   out_4589523013307382378[79] = 0;
   out_4589523013307382378[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2951645221556204453) {
   out_2951645221556204453[0] = state[0];
   out_2951645221556204453[1] = state[1];
   out_2951645221556204453[2] = state[2];
   out_2951645221556204453[3] = state[3];
   out_2951645221556204453[4] = state[4];
   out_2951645221556204453[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2951645221556204453[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2951645221556204453[7] = state[7];
   out_2951645221556204453[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8577401023761808986) {
   out_8577401023761808986[0] = 1;
   out_8577401023761808986[1] = 0;
   out_8577401023761808986[2] = 0;
   out_8577401023761808986[3] = 0;
   out_8577401023761808986[4] = 0;
   out_8577401023761808986[5] = 0;
   out_8577401023761808986[6] = 0;
   out_8577401023761808986[7] = 0;
   out_8577401023761808986[8] = 0;
   out_8577401023761808986[9] = 0;
   out_8577401023761808986[10] = 1;
   out_8577401023761808986[11] = 0;
   out_8577401023761808986[12] = 0;
   out_8577401023761808986[13] = 0;
   out_8577401023761808986[14] = 0;
   out_8577401023761808986[15] = 0;
   out_8577401023761808986[16] = 0;
   out_8577401023761808986[17] = 0;
   out_8577401023761808986[18] = 0;
   out_8577401023761808986[19] = 0;
   out_8577401023761808986[20] = 1;
   out_8577401023761808986[21] = 0;
   out_8577401023761808986[22] = 0;
   out_8577401023761808986[23] = 0;
   out_8577401023761808986[24] = 0;
   out_8577401023761808986[25] = 0;
   out_8577401023761808986[26] = 0;
   out_8577401023761808986[27] = 0;
   out_8577401023761808986[28] = 0;
   out_8577401023761808986[29] = 0;
   out_8577401023761808986[30] = 1;
   out_8577401023761808986[31] = 0;
   out_8577401023761808986[32] = 0;
   out_8577401023761808986[33] = 0;
   out_8577401023761808986[34] = 0;
   out_8577401023761808986[35] = 0;
   out_8577401023761808986[36] = 0;
   out_8577401023761808986[37] = 0;
   out_8577401023761808986[38] = 0;
   out_8577401023761808986[39] = 0;
   out_8577401023761808986[40] = 1;
   out_8577401023761808986[41] = 0;
   out_8577401023761808986[42] = 0;
   out_8577401023761808986[43] = 0;
   out_8577401023761808986[44] = 0;
   out_8577401023761808986[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8577401023761808986[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8577401023761808986[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8577401023761808986[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8577401023761808986[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8577401023761808986[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8577401023761808986[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8577401023761808986[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8577401023761808986[53] = -9.8000000000000007*dt;
   out_8577401023761808986[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8577401023761808986[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8577401023761808986[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8577401023761808986[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8577401023761808986[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8577401023761808986[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8577401023761808986[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8577401023761808986[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8577401023761808986[62] = 0;
   out_8577401023761808986[63] = 0;
   out_8577401023761808986[64] = 0;
   out_8577401023761808986[65] = 0;
   out_8577401023761808986[66] = 0;
   out_8577401023761808986[67] = 0;
   out_8577401023761808986[68] = 0;
   out_8577401023761808986[69] = 0;
   out_8577401023761808986[70] = 1;
   out_8577401023761808986[71] = 0;
   out_8577401023761808986[72] = 0;
   out_8577401023761808986[73] = 0;
   out_8577401023761808986[74] = 0;
   out_8577401023761808986[75] = 0;
   out_8577401023761808986[76] = 0;
   out_8577401023761808986[77] = 0;
   out_8577401023761808986[78] = 0;
   out_8577401023761808986[79] = 0;
   out_8577401023761808986[80] = 1;
}
void h_25(double *state, double *unused, double *out_6699396215809902969) {
   out_6699396215809902969[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2083213777335141338) {
   out_2083213777335141338[0] = 0;
   out_2083213777335141338[1] = 0;
   out_2083213777335141338[2] = 0;
   out_2083213777335141338[3] = 0;
   out_2083213777335141338[4] = 0;
   out_2083213777335141338[5] = 0;
   out_2083213777335141338[6] = 1;
   out_2083213777335141338[7] = 0;
   out_2083213777335141338[8] = 0;
}
void h_24(double *state, double *unused, double *out_5802721014004881644) {
   out_5802721014004881644[0] = state[4];
   out_5802721014004881644[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6829723104301650372) {
   out_6829723104301650372[0] = 0;
   out_6829723104301650372[1] = 0;
   out_6829723104301650372[2] = 0;
   out_6829723104301650372[3] = 0;
   out_6829723104301650372[4] = 1;
   out_6829723104301650372[5] = 0;
   out_6829723104301650372[6] = 0;
   out_6829723104301650372[7] = 0;
   out_6829723104301650372[8] = 0;
   out_6829723104301650372[9] = 0;
   out_6829723104301650372[10] = 0;
   out_6829723104301650372[11] = 0;
   out_6829723104301650372[12] = 0;
   out_6829723104301650372[13] = 0;
   out_6829723104301650372[14] = 1;
   out_6829723104301650372[15] = 0;
   out_6829723104301650372[16] = 0;
   out_6829723104301650372[17] = 0;
}
void h_30(double *state, double *unused, double *out_1063869041010769472) {
   out_1063869041010769472[0] = state[4];
}
void H_30(double *state, double *unused, double *out_435119181172107289) {
   out_435119181172107289[0] = 0;
   out_435119181172107289[1] = 0;
   out_435119181172107289[2] = 0;
   out_435119181172107289[3] = 0;
   out_435119181172107289[4] = 1;
   out_435119181172107289[5] = 0;
   out_435119181172107289[6] = 0;
   out_435119181172107289[7] = 0;
   out_435119181172107289[8] = 0;
}
void h_26(double *state, double *unused, double *out_8665038112830089818) {
   out_8665038112830089818[0] = state[7];
}
void H_26(double *state, double *unused, double *out_5824717096209197562) {
   out_5824717096209197562[0] = 0;
   out_5824717096209197562[1] = 0;
   out_5824717096209197562[2] = 0;
   out_5824717096209197562[3] = 0;
   out_5824717096209197562[4] = 0;
   out_5824717096209197562[5] = 0;
   out_5824717096209197562[6] = 0;
   out_5824717096209197562[7] = 1;
   out_5824717096209197562[8] = 0;
}
void h_27(double *state, double *unused, double *out_6745674590238803002) {
   out_6745674590238803002[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2658713252356050506) {
   out_2658713252356050506[0] = 0;
   out_2658713252356050506[1] = 0;
   out_2658713252356050506[2] = 0;
   out_2658713252356050506[3] = 1;
   out_2658713252356050506[4] = 0;
   out_2658713252356050506[5] = 0;
   out_2658713252356050506[6] = 0;
   out_2658713252356050506[7] = 0;
   out_2658713252356050506[8] = 0;
}
void h_29(double *state, double *unused, double *out_3937804226650076173) {
   out_3937804226650076173[0] = state[1];
}
void H_29(double *state, double *unused, double *out_945350525486499473) {
   out_945350525486499473[0] = 0;
   out_945350525486499473[1] = 1;
   out_945350525486499473[2] = 0;
   out_945350525486499473[3] = 0;
   out_945350525486499473[4] = 0;
   out_945350525486499473[5] = 0;
   out_945350525486499473[6] = 0;
   out_945350525486499473[7] = 0;
   out_945350525486499473[8] = 0;
}
void h_28(double *state, double *unused, double *out_8664283432026415929) {
   out_8664283432026415929[0] = state[0];
}
void H_28(double *state, double *unused, double *out_4137048491583031101) {
   out_4137048491583031101[0] = 1;
   out_4137048491583031101[1] = 0;
   out_4137048491583031101[2] = 0;
   out_4137048491583031101[3] = 0;
   out_4137048491583031101[4] = 0;
   out_4137048491583031101[5] = 0;
   out_4137048491583031101[6] = 0;
   out_4137048491583031101[7] = 0;
   out_4137048491583031101[8] = 0;
}
void h_31(double *state, double *unused, double *out_4976512631549297711) {
   out_4976512631549297711[0] = state[8];
}
void H_31(double *state, double *unused, double *out_6450925198442549038) {
   out_6450925198442549038[0] = 0;
   out_6450925198442549038[1] = 0;
   out_6450925198442549038[2] = 0;
   out_6450925198442549038[3] = 0;
   out_6450925198442549038[4] = 0;
   out_6450925198442549038[5] = 0;
   out_6450925198442549038[6] = 0;
   out_6450925198442549038[7] = 0;
   out_6450925198442549038[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_1362453500626495482) {
  err_fun(nom_x, delta_x, out_1362453500626495482);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1829696213295155707) {
  inv_err_fun(nom_x, true_x, out_1829696213295155707);
}
void car_H_mod_fun(double *state, double *out_4589523013307382378) {
  H_mod_fun(state, out_4589523013307382378);
}
void car_f_fun(double *state, double dt, double *out_2951645221556204453) {
  f_fun(state,  dt, out_2951645221556204453);
}
void car_F_fun(double *state, double dt, double *out_8577401023761808986) {
  F_fun(state,  dt, out_8577401023761808986);
}
void car_h_25(double *state, double *unused, double *out_6699396215809902969) {
  h_25(state, unused, out_6699396215809902969);
}
void car_H_25(double *state, double *unused, double *out_2083213777335141338) {
  H_25(state, unused, out_2083213777335141338);
}
void car_h_24(double *state, double *unused, double *out_5802721014004881644) {
  h_24(state, unused, out_5802721014004881644);
}
void car_H_24(double *state, double *unused, double *out_6829723104301650372) {
  H_24(state, unused, out_6829723104301650372);
}
void car_h_30(double *state, double *unused, double *out_1063869041010769472) {
  h_30(state, unused, out_1063869041010769472);
}
void car_H_30(double *state, double *unused, double *out_435119181172107289) {
  H_30(state, unused, out_435119181172107289);
}
void car_h_26(double *state, double *unused, double *out_8665038112830089818) {
  h_26(state, unused, out_8665038112830089818);
}
void car_H_26(double *state, double *unused, double *out_5824717096209197562) {
  H_26(state, unused, out_5824717096209197562);
}
void car_h_27(double *state, double *unused, double *out_6745674590238803002) {
  h_27(state, unused, out_6745674590238803002);
}
void car_H_27(double *state, double *unused, double *out_2658713252356050506) {
  H_27(state, unused, out_2658713252356050506);
}
void car_h_29(double *state, double *unused, double *out_3937804226650076173) {
  h_29(state, unused, out_3937804226650076173);
}
void car_H_29(double *state, double *unused, double *out_945350525486499473) {
  H_29(state, unused, out_945350525486499473);
}
void car_h_28(double *state, double *unused, double *out_8664283432026415929) {
  h_28(state, unused, out_8664283432026415929);
}
void car_H_28(double *state, double *unused, double *out_4137048491583031101) {
  H_28(state, unused, out_4137048491583031101);
}
void car_h_31(double *state, double *unused, double *out_4976512631549297711) {
  h_31(state, unused, out_4976512631549297711);
}
void car_H_31(double *state, double *unused, double *out_6450925198442549038) {
  H_31(state, unused, out_6450925198442549038);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_lib_init(car)
