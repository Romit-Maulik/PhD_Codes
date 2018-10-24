/** \file main.cpp
 *  \brief One dimensional Wavelet Adapted Mesh Refinement Code \n
 *  Systems to solve: \n
 *  Burgers - Baeza, A., Martínez-Gavara, A., & Mulet, P. (2012). Adaptation based on \n
 *  interpolation errors for high order mesh refinement methods applied to conservation laws. \n
 *  Applied Numerical Mathematics, 62(4), 278-296. \n

 *  Shallow Water Equations - Ambrosi, D. "Approximation of shallow water equations by Roe's Riemann solver." \n
 *  International journal for numerical methods in fluids 20.2 (1995): 157-168. \n

 *  Sod Shock Tube - http://www.csun.edu/~jb715473/examples/euler1d.htm \n

 *  Brio-Wu Shock Tube - http://www.csun.edu/~jb715473/examples/mhd1d.htm \n

 *  AMR Routine: \n
 *  Roussel, O., Schneider, K., Tsigulin, A., & Bockhorn, H. (2003). \n
 *  A conservative fully adaptive multiresolution algorithm for parabolic PDEs. \n
 *  Journal of Computational Physics, 188(2), 493-523. \n

 *  Finite Volume Central Flux Computations: \n
 *  Hyman, James M., Robert J. Knapp, and James C. Scovel. \n
 *  "High order finite volume approximations of differential operators on nonuniform grids." \n
 *  Physica D: Nonlinear Phenomena 60.1 (1992): 112-138. \n

 *  Shock Capturing options: \n

 *  MUSCL - MinMod/MC \n
 *  San, Omer, and Kursat Kara. \n
 *  "Numerical assessments of high-order accurate shock capturing schemes: \n
 *  Kelvin–Helmholtz type vortical structures in high-resolutions." \n
 *  Computers & Fluids 89 (2014): 254-276. \n

 *  WENO3/WENO5 \n
 *  San, Omer, and Kursat Kara. "Evaluation of Riemann flux solvers for WENO reconstruction schemes: Kelvin–Helmholtz instability." \n
 *  Computers & Fluids 117 (2015): 24-41. \n

 *  WENO6 - Self Derived \n

 *  Riemann Solver: \n
 *  Rusanov -- Maulik, Romit, and Omer San. \n
 *  "Resolution and Energy Dissipation Characteristics of Implicit LES and Explicit Filtering Models for Compressible Turbulence." \n
 *  Fluids 2.2 (2017): 14. \n

 *  FORCE -- \n
 *  San, Omer, and Kursat Kara. "Evaluation of Riemann flux solvers for WENO reconstruction schemes: Kelvin–Helmholtz instability." \n
 *  Computers & Fluids 117 (2015): 24-41. \n
 *
 *  \author Romit Maulik (romit.maulik@okstate.edu)
 */

#include "Control_Parameters.h"
#include "Global.h"
#include "Function_Def.h"
#include "Init_Domain.h"
#include "Datastruct.h"
#include "Wavelet.h"
#include "Output.h"
#include "Face_Fluxes.h"
#include "Force_flux.h"
#include "Rusanov_Solver.h"
#include "Time_Int.h"
#include "Static_Eval.h"

int main()
{
    unordered_map<int,Cell*>Cellvect;
//    analyse_static(Cellvect);


    create_domain(Cellvect);

    cout<<"Enter number of timestep outputs "<<endl;
    cin>>snap_steps;
    evolution(Cellvect);


    return 0;
}
