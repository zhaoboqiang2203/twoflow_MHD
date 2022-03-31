#pragma once

#include "MHD.h"
void atom_flow();
struct _AF cal_atom_u(int i, int j);
struct _AF cal_atom_Fqr(struct _AF uij);
struct _AF cal_atom_Fqz(struct _AF uij);
struct _AF cal_atom_s(struct _AF uij);
struct _AF arti_atom_vis(struct _particle np, struct  _particle n, struct  _particle nn);
void atom_boundary_condition();
