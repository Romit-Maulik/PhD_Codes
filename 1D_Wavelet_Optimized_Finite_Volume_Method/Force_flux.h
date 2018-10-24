#ifndef Force_flux_h_inluded
#define Force_flux_h_inluded


/** \file Force_Flux.h
 *  \brief Header file for FORCE based flux reconstructions from left and right states at a cell face \n
 *  Not used IN IJNMF \n
 *  Requires specification of timestep dt \n
*/


/** \brief This function takes a memory address to a Cell as input and calculates face fluxes \n
 *  Choice of state reconstruction is determined by mflux option (i.e. WENO3,WENO5,MUSCL etc) \n
 *  Note that wavespeeds for the FORCE reconstruction are calculated according to the local stencil. \n
 *  Flux reconstructions are carried out at both left and right faces. \n
 *  Fixed timestep specified by user
 * \param Cell* &cellval
 * \return void
 *
 */
void force_solver(Cell* &cellval)
{
    double f_right[n_eq] = {0.0};
    double f_left[n_eq] = {0.0};
    double q_right[n_eq] = {0.0};
    double q_left[n_eq] = {0.0};

    double q[n_eq] = {0.0};
    double qp1[n_eq] = {0.0};
    double qp2[n_eq] = {0.0};
    double qp3[n_eq] = {0.0};
    double qm1[n_eq] = {0.0};
    double qm2[n_eq] = {0.0};
    double qm3[n_eq] = {0.0};



    for (int j=0;j<=n_eq-1;j++)
    {
        q[j] = cellval->q[j];
        qp1[j] = cellval->right_level->q[j];
        qp2[j] = cellval->right_level->right_level->q[j];
        qp3[j] = cellval->right_level->right_level->right_level->q[j];
        qm1[j] = cellval->left_level->q[j];
        qm2[j] = cellval->left_level->left_level->q[j];
        qm3[j] = cellval->left_level->left_level->left_level->q[j];
    }

    double c=cellval->cell_length/dt;

    if (mflux==1)
    {
        //Calculating right face flux
        weno3(cellval->level,qm1,q,qp1,qp2,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
    }
    else if (mflux==2)
    {
        //Calculating right face flux
        muscl(qm1,q,qp1,qp2,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
    }
    else if (mflux==3)
    {
        //Calculating right face flux
        weno5(cellval->level,qm2,qm1,q,qp1,qp2,qp3,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
    }
    else if (mflux==4)
    {
        //Calculating right face flux
        weno6(cellval->level,qm2,qm1,q,qp1,qp2,qp3,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
    }

    double lf_flux[n_eq] = {0.0};//Only use curly braces for initialization

    for (int j=0;j<=n_eq-1;j++)//Lax flux
    {
        lf_flux[j] = 0.5*(f_right[j]+f_left[j]) + 0.5*c*(q_left[j]-q_right[j]);
    }

    double ri_flux[n_eq] = {0.0};//Only use curly braces for initialization
    double qrich[n_eq] = {0.0};//Only use curly braces for initialization

    for (int j=0;j<=n_eq-1;j++)//Richardson reconstruction
    {
        qrich[j] = 0.5*(q_right[j]+q_left[j]) + 0.5/c*(f_left[j]-f_right[j]);
    }

    calculate_local_fluxes(qrich,ri_flux);//Richardson Flux

    for (int j=0;j<=n_eq-1;j++)
    {
        cellval->f_p[j] = 0.5*(lf_flux[j]+ri_flux[j]);//Averaging Lax and Richardson flux
    }

    if (mflux==1)
    {
        //Calculating left face flux
        weno3(cellval->level,qm2,qm1,q,qp1,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
    }
    else if (mflux==2)
    {
        //Calculating left face flux
        muscl(qm2,qm1,q,qp1,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
    }
    else if (mflux==3)
    {
        //Calculating left face flux
        weno5(cellval->level,qm3,qm2,qm1,q,qp1,qp2,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
    }
    else if (mflux==4)
    {
        //Calculating left face flux
        weno6(cellval->level,qm3,qm2,qm1,q,qp1,qp2,q_right,q_left);
        calculate_local_fluxes(q_right,f_right);
        calculate_local_fluxes(q_left,f_left);
    }

    for (int j=0;j<=n_eq-1;j++)
    {
        lf_flux[j] = 0.5*(f_right[j]+f_left[j]) + 0.5*c*(q_left[j]-q_right[j]);
    }

    for (int j=0;j<=n_eq-1;j++)
    {
        qrich[j] = 0.5*(q_right[j]+q_left[j]) + 0.5/c*(f_left[j]-f_right[j]);
    }

    calculate_local_fluxes(qrich,ri_flux);

    for (int j=0;j<=n_eq-1;j++)
    {
        cellval->f_m[j] = 0.5*(lf_flux[j]+ri_flux[j]);
    }
}







#endif // Force_flux_h_inluded
