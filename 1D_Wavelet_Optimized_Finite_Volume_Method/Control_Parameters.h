#ifndef Control_Pameters_h_inluded
#define Control_Pameters_h_inluded
#include <math.h> //Needed for atan definition
/** \file Control_Parameters.h
 *  \brief Coefficient values,solver flags,choice of problem \n
 *  Switches for dynamic grid refinement etc.
*/

//Datastructure related parameters
const int max_level = 12;/**< preset maximum level of refinement */
const double wc_threshold = 1.0e-3;/**< Wavelet coefficient truncation threshold */
const int n_eq = 3;/**< 1:Burgers, 2:SWE, 3:Euler-1D, 7:MHD-1D*/
const int bc_type = 2;/**< 1 - Periodic, 2 - Open*/
const int refine_switch = 1;/**< 1 - Refine grid dynamically, 0 - Uniform grid*/
const double p_w = 4.0;

//Finite Volume methods
const int mflux = 3;/**<1-WENO3,2-MUSCL,3-WENO5,4-WENO6*/
const int diss_level_thresh=1;/**<Level at which upwinding begins*/
const int selective_weno = 1;/**<0 - Purely WENO, 1-Hybrid*/
const int mlim = 1;/**<Options for MUSCL, 1-MinMod,2-MC*/

//Simulation parameters
const double cfl=0.8;/**<Global CFL*/
double dt;/**<Global dt*/
const double final_time=0.2;/**< Final time*/
int snap_steps;
int temp_min_level;/**< Tracks minimum level of resolution at which a leaf resides*/
int temp_max_level;/**< Tracks maximum level of resolution at which a leaf resides*/


//Riemann Solver
const int riesol = 1;/**< 1- Rusanov; 2-FORCE*/
const double force_dt = 1.0e-4;

//Problem Specific Constants - 1D Euler
//Physical Parameters
double gmm = 7.0/5.0;/**< Ratio of Specific heats (both Euler and MHD)*/

//Problem Specific Constants - 1D MHD
const double bx = 0.75;/**< MHD Constant*/


//Problem Specific Constants - 1D NLE
// const double rho0 = 8.93;
// const double cv =3.9e-4;
// const double t0 = 300.0;
// const double alpha = 1.0;
// const double beta = 3.0;
// const double gamma = 2.0;

// const double cc0 = 4.6;
// const double bb0 = 2.1;

// const double k0 = cc0*cc0 - (4.0/3.0)*bb0*bb0;
// const double b0 = bb0*bb0;

//Problem specific constants - SWE - Classical DamBreak
const double height_1 = 10.0;/**< SWE Constant 1*/
const double height_2 = 1.0;/**< SWE Constant 2*/
const double grav = 1.0;/**< SWE Constant 3*/

#endif // Control_Pameters_h_inluded
