//
// Created by alberto on 5/2/18.
//
#include "struct.geotop.h"
#include "t_datamanipulation.h"

STATEVAR_3D::STATEVAR_3D(double nan, long nl, long nr, long nc) :
        type{new Matrix<short>{nr, nc}},
        lnum{new Matrix<long>{nr, nc}},
        Dzl{new_doubletensor(nl, nr, nc)},
        w_liq{new_doubletensor(nl, nr, nc)},
        w_ice{new_doubletensor(nl, nr, nc)},
        T{new_doubletensor(nl, nr, nc)}
{
  *type = 2;
  *lnum = 0;
  initialize_doubletensor(Dzl, 0.);
  initialize_doubletensor(T, nan);
  initialize_doubletensor(w_ice, 0.);
  initialize_doubletensor(w_liq, 0.);
}

STATEVAR_3D::~STATEVAR_3D() {
  //free_shortmatrix(type);
 //free_longmatrix(lnum);
  free_doubletensor(Dzl);
  free_doubletensor(T);
  free_doubletensor(w_ice);
  free_doubletensor(w_liq);
}

SOIL_STATE::SOIL_STATE(const long n, const long nl) :
        T{new Matrix<double>{nl,n}},
        P{new Matrix<double>{nl,0,n,1}},
        thi{new Matrix<double>{nl,n}}
{
  T.reset(new Matrix<double>{nl,n});
  P.reset(new Matrix<double>{nl,0,n,1});
  thi.reset(new Matrix<double>{nl,n});
}

