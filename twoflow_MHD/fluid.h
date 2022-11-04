#pragma once

#include "MHD.h"



void electron_flow();
void electron_flow_v2();
void ion_flow();
void electron_flow_openmp();

struct _U cal_u(int i, int j);
struct _F cal_fr(struct _U uij);
struct _F cal_fz(struct _U uij);
struct _F cal_s(struct _U uij);
struct _U arti_vis(struct  node np, struct  node n, struct  node nn);

double dy_vis(int i, int j);
void cal_tau();
void cal_tau_bar();

void cal_tau_omp();
void cal_tau_bar_omp();

void smooth_ne(int i, int j);
void smooth_pe(int i, int j);