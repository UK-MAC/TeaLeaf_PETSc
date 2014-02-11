/*Crown Copyright 2014 AWE.
*
* This file is part of TeaLeaf.
*
* TeaLeaf is free software: you can redistribute it and/or modify it under
* the terms of the GNU General Public License as published by the
* Free Software Foundation, either version 3 of the License, or (at your option)
* any later version.
*
* TeaLeaf is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public License along with
* TeaLeaf. If not, see http://www.gnu.org/licenses/. */

/**
 *  @brief C heat conduction kernel
 *  @author David Beckingsale, Wayne Gaudin
 *  @details Implicitly calculates the change in temperature using a Jacobi iteration
 */

#include <stdio.h>
#include <stdlib.h>
#include "ftocmacros.h"
#include <math.h>

/* Coefficient constants taken from data.f90 */
#define CONDUCTIVITY 1
#define RECIP_CONDUCTIVITY 2

void tea_leaf_kernel_init_c_(
        int* xmin,
        int* xmax,
        int* ymin,
        int* ymax,
        double* celldx,
        double* celldy,
        double* volume,
        double* density,
        double* energy,
        double* u0,
        double* u1,
        double* un,
        double* heat_capacity,
        double* Kx_tmp,
        double* Ky_tmp,
        double* Kx,
        double* Ky,
        int* coef)
{
    int x_min = *xmin;
    int x_max = *xmax;
    int y_min = *ymin;
    int y_max = *ymax;
    int coefficient = *coef;

    int j,k;

#pragma omp parallel
    {

        if(coefficient == RECIP_CONDUCTIVITY) {
#pragma omp for 
            for(k = y_min-1; k <= y_max+1; k++) {
                for(j = x_min-1; j <= x_max+1; j++) {
                    Kx_tmp[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]=1.0/density[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)];
                    Kx_tmp[FTNREF2D(j+1,k,x_max+5,x_min-2,y_min-2)]=1.0/density[FTNREF2D(j+1,k,x_max+4,x_min-2,y_min-2)];
                    Ky_tmp[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]=1.0/density[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)];
                    Ky_tmp[FTNREF2D(j,k+1,x_max+5,x_min-2,y_min-2)]=1.0/density[FTNREF2D(j,k+1,x_max+4,x_min-2,y_min-2)];
                }
            }
        }
        else if(coefficient == CONDUCTIVITY) {
#pragma omp for
            for(k = y_min-1; k <= y_max+1; k++) {
                for(j = x_min-1; j <= x_max+1; j++) {
                    Kx_tmp[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]=density[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)];
                    Kx_tmp[FTNREF2D(j+1,k,x_max+5,x_min-2,y_min-2)]=density[FTNREF2D(j+1,k,x_max+4,x_min-2,y_min-2)];
                    Ky_tmp[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]=density[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)];
                    Ky_tmp[FTNREF2D(j,k+1,x_max+5,x_min-2,y_min-2)]=density[FTNREF2D(j,k+1,x_max+4,x_min-2,y_min-2)];
                }
            }
        }

#pragma omp for
        for(k = y_min-1; k <= y_max+1; k++) {
            for(j = x_min-1; j <= x_max+1; j++) {
                Kx[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)] = 
                    (Kx_tmp[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]+Kx_tmp[FTNREF2D(j+1,k,x_max+5,x_min-2,y_min-2)])
                    /(2.0*Kx_tmp[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]*Kx_tmp[FTNREF2D(j+1,k,x_max+5,x_min-2,y_min-2)]);

                Ky[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)] = 
                    (Ky_tmp[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]+Ky_tmp[FTNREF2D(j,k+1,x_max+5,x_min-2,y_min-2)])
                    /(2.0*Ky_tmp[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]*Ky_tmp[FTNREF2D(j,k+1,x_max+5,x_min-2,y_min-2)]);
            }
        }

#pragma omp for 
        for(k = y_min-1; k <=  y_max+1; k++) {
            for(j = x_min-1; j <=  x_max+1; j++) {
                u0[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)] =  
                    energy[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)]
                    * density[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)];
            }
        }

#pragma omp for
            for(k = y_min-1; k <=  y_max+1; k++) {
                for(j = x_min-1; j <=  x_max+1; j++) {
                    u1[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)] = u0[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)];
                }
            }
    }
}

void tea_leaf_kernel_solve_c_(
        int *xmin,
        int* xmax,
        int* ymin,
        int* ymax,
        double* rxp,
        double* ryp,
        double* Kx,
        double* Ky,
        double* error,
        double* u0,
        double* u1,
        double* un)
{
    int x_min = *xmin;
    int x_max = *xmax;
    int y_min = *ymin;
    int y_max = *ymax;

    int j,k;

    double rx = *rxp;
    double ry = *ryp;

    *error = 0.0;

#pragma omp parallel
    {
#pragma omp for
        for(k = y_min-1; k <=  y_max+1; k++) {
            for(j = x_min-1; j <=  x_max+1; j++) {
                un[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)] = u1[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)];
            }
        }
    }

    // #pragma omp for reduction(max:error)
    for(k = y_min; k <=  y_max; k++) {
        for(j = x_min; j <=  x_max; j++) {

            u1[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)] = 
                (u0[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)] 
                 + Kx[FTNREF2D(j+1,k,x_max+5,x_min-2,y_min-2)]*rx*un[FTNREF2D(j+1,k,x_max+5,x_min-2,y_min-2)]
                 + Kx[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]*rx*un[FTNREF2D(j-1,k,x_max+5,x_min-2,y_min-2)]
                 + Ky[FTNREF2D(j,k+1,x_max+5,x_min-2,y_min-2)]*ry*un[FTNREF2D(j,k+1,x_max+5,x_min-2,y_min-2)] 
                 + Ky[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]*ry*un[FTNREF2D(j,k-1,x_max+5,x_min-2,y_min-2)])
                / (1+2.0*rx+2.0*ry);

            *error = fmax(*error,
                    fabs(u1[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)] 
                        - u0[FTNREF2D(j,k,x_max+5,x_min-2,y_min-2)]));
        }
    }
    //  }
}

void tea_leaf_kernel_finalise_c_(
        int* xmin,
        int* xmax,
        int* ymin,
        int* ymax,
        double* energy,
        double* density,
        double* u)
{
    int x_min = *xmin;
    int x_max = *xmax;
    int y_min = *ymin;
    int y_max = *ymax;

    int j,k;

#pragma omp parallel for 
    for(k = y_min; k <=  y_max; k++) {
        for(j = x_min; j <=  x_max; j++) {
            energy[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)] = 
                u[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)] / density[FTNREF2D(j,k,x_max+4,x_min-2,y_min-2)];
        }
    }
}
