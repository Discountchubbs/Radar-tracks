#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_6095997321851776631);
void live_err_fun(double *nom_x, double *delta_x, double *out_1020174459653040964);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_170284851571879567);
void live_H_mod_fun(double *state, double *out_762102343625566160);
void live_f_fun(double *state, double dt, double *out_5052328704831362084);
void live_F_fun(double *state, double dt, double *out_6576080517812039271);
void live_h_4(double *state, double *unused, double *out_172293272178988194);
void live_H_4(double *state, double *unused, double *out_2162120159528067895);
void live_h_9(double *state, double *unused, double *out_7611091542703511884);
void live_H_9(double *state, double *unused, double *out_1920930512898477250);
void live_h_10(double *state, double *unused, double *out_3864878420190995903);
void live_H_10(double *state, double *unused, double *out_2685759071093232638);
void live_h_12(double *state, double *unused, double *out_1884113930817331366);
void live_H_12(double *state, double *unused, double *out_2857336248503893900);
void live_h_35(double *state, double *unused, double *out_5649298364682089177);
void live_H_35(double *state, double *unused, double *out_5602899280828907609);
void live_h_32(double *state, double *unused, double *out_2637127511937211224);
void live_H_32(double *state, double *unused, double *out_5734064065769648078);
void live_h_13(double *state, double *unused, double *out_5623579500843631738);
void live_H_13(double *state, double *unused, double *out_74919045099890815);
void live_h_14(double *state, double *unused, double *out_7611091542703511884);
void live_H_14(double *state, double *unused, double *out_1920930512898477250);
void live_h_33(double *state, double *unused, double *out_1277502705475652847);
void live_H_33(double *state, double *unused, double *out_8753456285467765213);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}