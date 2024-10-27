#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_1362453500626495482);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1829696213295155707);
void car_H_mod_fun(double *state, double *out_4589523013307382378);
void car_f_fun(double *state, double dt, double *out_2951645221556204453);
void car_F_fun(double *state, double dt, double *out_8577401023761808986);
void car_h_25(double *state, double *unused, double *out_6699396215809902969);
void car_H_25(double *state, double *unused, double *out_2083213777335141338);
void car_h_24(double *state, double *unused, double *out_5802721014004881644);
void car_H_24(double *state, double *unused, double *out_6829723104301650372);
void car_h_30(double *state, double *unused, double *out_1063869041010769472);
void car_H_30(double *state, double *unused, double *out_435119181172107289);
void car_h_26(double *state, double *unused, double *out_8665038112830089818);
void car_H_26(double *state, double *unused, double *out_5824717096209197562);
void car_h_27(double *state, double *unused, double *out_6745674590238803002);
void car_H_27(double *state, double *unused, double *out_2658713252356050506);
void car_h_29(double *state, double *unused, double *out_3937804226650076173);
void car_H_29(double *state, double *unused, double *out_945350525486499473);
void car_h_28(double *state, double *unused, double *out_8664283432026415929);
void car_H_28(double *state, double *unused, double *out_4137048491583031101);
void car_h_31(double *state, double *unused, double *out_4976512631549297711);
void car_H_31(double *state, double *unused, double *out_6450925198442549038);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}