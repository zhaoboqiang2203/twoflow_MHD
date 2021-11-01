#pragma once

#include "MHD.h"

void electron_flow();
void ion_flow();

struct _U cal_u(int i, int j);
struct _F cal_fr(struct _U uij);
struct _F cal_fz(struct _U uij);
struct _F cal_s(struct _U uij);