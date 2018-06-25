//
// Created by alberto on 5/2/18.
//
#include "struct.geotop.h"
#include "t_datamanipulation.h"

STATEVAR_3D::STATEVAR_3D(double nan, long nl, long nr, long nc) :
        type{new_shortmatrix(nr, nc)},
        lnum{new_longmatrix(nr, nc)},
        Dzl{new_doubletensor(nl, nr, nc)},
        w_liq{new_doubletensor(nl, nr, nc)},
        w_ice{new_doubletensor(nl, nr, nc)},
        T{new_doubletensor(nl, nr, nc)}
{
  initialize_shortmatrix(type, 2);
  initialize_longmatrix(lnum, 0);
  initialize_doubletensor(Dzl, 0.);
  initialize_doubletensor(T, nan);
  initialize_doubletensor(w_ice, 0.);
  initialize_doubletensor(w_liq, 0.);
}

STATEVAR_3D::~STATEVAR_3D() {
  free_shortmatrix(type);
  free_longmatrix(lnum);
  free_doubletensor(Dzl);
  free_doubletensor(T);
  free_doubletensor(w_ice);
  free_doubletensor(w_liq);
}

SOIL_STATE::SOIL_STATE(const long n, const long nl) :
        T{new_doublematrix(nl, n)},
        P{new_doublematrix0_(nl, n)},
        thi{new_doublematrix(nl, n)}
{
  T = new_doublematrix(nl, n);
  initialize_doublematrix(T, 0.);
  P = new_doublematrix0_(nl, n);
  initialize_doublematrix(P, 0.);
  thi = new_doublematrix(nl, n);
  initialize_doublematrix(thi, 0.);

}

SOIL_STATE::~SOIL_STATE() {

  free_doublematrix(T);
  free_doublematrix(P);
  free_doublematrix(thi);
}
