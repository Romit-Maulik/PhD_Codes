#ifndef face_fluxes_h_inluded
#define face_fluxes_h_inluded

/** \file Face_Fluxes.h
 *  \brief Header file for datastructure related functions \n
*/

/** \brief This function calculates central fluxes through finite volume interpolation \n
 *  Refer J.M. Hyman et al 1992, Physica D for Finite Volume Flux Constructions \n
 *  These fluxes are nondissipative
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void calculate_central_fluxes(unordered_map<int,Cell*> &leafvector)
{
    for (auto citer_one=leafvector.begin();citer_one!=leafvector.end();citer_one++)
    {
        if ((*citer_one).second->virt == 0)//Only for true leaves
        {

                double f_plus_one[n_eq] = {0.0};
                double f_plus_two[n_eq] = {0.0};
                double f_plus_three[n_eq] = {0.0};
                double f_minus_one[n_eq] = {0.0};
                double f_minus_two[n_eq] = {0.0};
                double f_minus_three[n_eq] = {0.0};

                //Calculating right face flux
                calculate_local_fluxes((*citer_one).second->right_level->q,f_plus_one);
                calculate_local_fluxes((*citer_one).second->right_level->right_level->q,f_plus_two);
                calculate_local_fluxes((*citer_one).second->right_level->right_level->right_level->q,f_plus_three);
                calculate_local_fluxes((*citer_one).second->q,f_minus_one);
                calculate_local_fluxes((*citer_one).second->left_level->q,f_minus_two);
                calculate_local_fluxes((*citer_one).second->left_level->left_level->q,f_minus_three);


                for (int j=0;j<=n_eq-1;j++)
                {
                    //Refer Table 1 Hyman Paper
                    double flux = (f_minus_three[j] - 8.0*f_minus_two[j] + 37.0*f_minus_one[j] + 37.0*f_plus_one[j] - 8.0*f_plus_two[j] + 1.0*f_plus_three[j])/60.0;
                    (*citer_one).second->f_p[j] = flux;
                }

                //Calculating left face flux
                calculate_local_fluxes((*citer_one).second->q,f_plus_one);
                calculate_local_fluxes((*citer_one).second->right_level->q,f_plus_two);
                calculate_local_fluxes((*citer_one).second->right_level->right_level->q,f_plus_three);
                calculate_local_fluxes((*citer_one).second->left_level->q,f_minus_one);
                calculate_local_fluxes((*citer_one).second->left_level->left_level->q,f_minus_two);
                calculate_local_fluxes((*citer_one).second->left_level->left_level->left_level->q,f_minus_three);

                for (int j=0;j<=n_eq-1;j++)
                {
                    //Refer Table 1 Hyman Paper
                    double flux = (f_minus_three[j] - 8.0*f_minus_two[j] + 37.0*f_minus_one[j] + 37.0*f_plus_one[j] - 8.0*f_plus_two[j] + 1.0*f_plus_three[j])/60.0;
                    (*citer_one).second->f_m[j] = flux;
                }

        }
    }
}


/** \brief This function calculates dissipative (or mononotic) fluxes through combination of finite volume interpolation and upwinding \n
 *  Requires a global diss_level_thresh parameter under which only central interpolation happens - refer Control_Parameters.h \n
 *  Above this minimum threshold for upwinding, we may apply WENO5/WENO3 through fully nonlinear weights or scale-selective weights (IJNMF concept) \n
 *  We also have option for MUSCL
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void calculate_monotonic_fluxes(unordered_map<int,Cell*> &leafvector)
{
    for (auto citer_one=leafvector.begin();citer_one!=leafvector.end();citer_one++)
    {
        if ((*citer_one).second->virt == 0)//Only for true leaves
        {
            if ((*citer_one).second->level<diss_level_thresh)
            {
                double f_plus_one[n_eq] = {0.0};
                double f_plus_two[n_eq] = {0.0};
                double f_plus_three[n_eq] = {0.0};
                double f_minus_one[n_eq] = {0.0};
                double f_minus_two[n_eq] = {0.0};
                double f_minus_three[n_eq] = {0.0};

                //Calculating right face flux
                calculate_local_fluxes((*citer_one).second->right_level->q,f_plus_one);
                calculate_local_fluxes((*citer_one).second->right_level->right_level->q,f_plus_two);
                calculate_local_fluxes((*citer_one).second->right_level->right_level->right_level->q,f_plus_three);
                calculate_local_fluxes((*citer_one).second->q,f_minus_one);
                calculate_local_fluxes((*citer_one).second->left_level->q,f_minus_two);
                calculate_local_fluxes((*citer_one).second->left_level->left_level->q,f_minus_three);


                for (int j=0;j<=n_eq-1;j++)
                {
                    //Refer Table 1 Hyman Paper
                    double flux = (f_minus_three[j] - 8.0*f_minus_two[j] + 37.0*f_minus_one[j] + 37.0*f_plus_one[j] - 8.0*f_plus_two[j] + 1.0*f_plus_three[j])/60.0;
                    (*citer_one).second->f_p[j] = flux;
                }

                //Calculating right face flux
                calculate_local_fluxes((*citer_one).second->q,f_plus_one);
                calculate_local_fluxes((*citer_one).second->right_level->q,f_plus_two);
                calculate_local_fluxes((*citer_one).second->right_level->right_level->q,f_plus_three);
                calculate_local_fluxes((*citer_one).second->left_level->q,f_minus_one);
                calculate_local_fluxes((*citer_one).second->left_level->left_level->q,f_minus_two);
                calculate_local_fluxes((*citer_one).second->left_level->left_level->left_level->q,f_minus_three);

                for (int j=0;j<=n_eq-1;j++)
                {
                    //Refer Table 1 Hyman Paper
                    double flux = (f_minus_three[j] - 8.0*f_minus_two[j] + 37.0*f_minus_one[j] + 37.0*f_plus_one[j] - 8.0*f_plus_two[j] + 1.0*f_plus_three[j])/60.0;
                    (*citer_one).second->f_m[j] = flux;
                }

            }
            else//Use smoothing
            {
                if (riesol==1)
                {
                rusanov_solver((*citer_one).second);
                }
                else
                {
                force_solver((*citer_one).second);
                }
            }
        }
    }
}



/** \brief This function calculates local flux quantities given conserved variables \n
 *  Requires an input array of floating point conserved variables and an 'empty' array of floating point fluxes which are populated \n
 *  Flux definitions may be found in IJNMF \n
 * \param double (&q)[n_eq]
 * \param double (&f)[n_eq]
 * \return void
 *
 */
void calculate_local_fluxes(double (&q)[n_eq],double (&f)[n_eq])
{
    if (n_eq==1)
    {
        f[0] = 0.5*q[0]*q[0];
    }
    else if (n_eq==2)
    {
        f[0] = q[1];
        f[1] = q[1]*q[1]/q[0] + 0.5*q[0]*q[0]*grav;
    }
    else if (n_eq==3)
    {
        f[0] = q[1];
        double p = (gmm-1.0)*(q[2]-0.5*q[1]*q[1]/q[0]);
        f[1] = q[1]*q[1]/q[0] + p;
        f[2] = q[1]/q[0]*(q[2]+p);
    }
    else if (n_eq==7)
    {
        f[0] = q[1];

        double p = (gmm-1.0)*(q[6]-0.5*q[1]*q[1]/q[0]-0.5*q[2]*q[2]/q[0]-0.5*q[3]*q[3]/q[0]-0.5*(bx*bx + q[4]*q[4] + q[5]*q[5]));

        if (p<0.0){p=0.0;}

        double pstar = p + 0.5*(bx*bx + q[4]*q[4] + q[5]*q[5]);

        f[1] = q[1]*q[1]/q[0] + pstar - bx*bx;
        f[2] = q[2]*q[1]/q[0] - bx*q[4];
        f[3] = q[3]*q[1]/q[0] - bx*q[5];
        f[4] = q[4]*q[1]/q[0] - bx*q[2]/q[0];
        f[5] = q[5]*q[1]/q[0] - bx*q[3]/q[0];
        f[6] = (q[6]+pstar)*q[1]/q[0] - bx*(bx*q[1]/q[0] + q[4]*q[2]/q[0] + q[5]*q[3]/q[0]);
    }
    else if (n_eq==14)
    {
        //Calculating deformation gradient tensor
        int j=5;
        double def_ten[3][3] = {0.0};//Be careful

        for (int ii=0;ii<=2;ii++)
        {
            for (int jj=0;jj<=2;jj++)
            {
                def_ten[ii][jj] = q[j]/q[0];
                j++;
            }

        }

        double r = q[0];
        double u = q[1]/q[0];
        double v = q[2]/q[0];
        double w = q[3]/q[0];
        double e = q[4]/q[0] - 0.5*(u*u + v*v + w*w);

        double stress[3][3] = {0.0};//Be careful

        stress_calc(r,e,def_ten,stress);


        f[0] = r*u;
        f[1] = r*u*u - stress[0][0];
        f[2] = r*u*v - stress[1][0];
        f[3] = r*u*w - stress[2][0];
        f[4] = u*q[4] - u*stress[0][0]- v*stress[1][0]- w*stress[2][0];

        f[5] = 0.0;
        f[6] = 0.0;
        f[7] = 0.0;

        f[8] = r*(def_ten[1][0]*u - def_ten[0][0]*v);
        f[9] = r*(def_ten[1][1]*u - def_ten[0][1]*v);
        f[10] = r*(def_ten[1][2]*u - def_ten[0][2]*v);

        f[11] = r*(def_ten[2][0]*u - def_ten[0][0]*w);
        f[12] = r*(def_ten[2][1]*u - def_ten[0][1]*w);
        f[13] = r*(def_ten[2][2]*u - def_ten[0][2]*w);

    }
//    else if (n_eq==4)//Eq 12; Cheng-Toro CAF 2015 136-152
//    {
//        f[0] = q[1];
//
//        //Calculation of pressure from Mie-Gruneisen equation of state
//        double eta = q[0]/ep_rho_0;
//        double int_en = q[2]/q[0] - 0.5*q[1]*q[1]/q[0];
//        double mg_f = (eta-1.0)*(eta-ep_gamma_0*(eta-1.0)/2.0)/((eta-ep_s*(eta-1.0))*(eta-ep_s*(eta-1.0)));
//        double p = ep_rho_0*(ep_a_0*ep_a_0)*mg_f + ep_rho_0*ep_gamma_0*int_en;
//
//
//        f[1] = q[1]*q[1]/q[0] + p - q[3];
//        f[2] = q[1]/q[0]*(q[2] + p - q[3]);
//        f[3] = q[3]*q[1]/q[0];
//    }
}


/** \brief This function calculates local left and right state reconstructions using WENO3 \n
 *  Requires input arrays of floating point conserved variables at +- stencil locations \n
 *  Also requires 2 empty arrays for left and right state reconstructions
 *  Full mathematical setup may be found in IJNMF \n
 *  Note that option of scale-selective WENO enforces extra linearized weighting of nonlinear coefficients \n
 *  Formulated in the sense of the right cell face
 * \param double (&qm1)[n_eq]
 * \param double (&q)[n_eq]
 * \param double (&qp1)[n_eq]
 * \param double (&qp2)[n_eq]
 * \param double (&qright)[n_eq]
 * \param double (&qleft)[n_eq]
 * \return void
 *
 */
void weno3(int level, double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq], double (&qright)[n_eq], double (&qleft)[n_eq])
{
    double d1=1.0/3.0;
    double d2=2.0/3.0;

        for (int j=0;j<=n_eq-1;j++)
        {
            double beta1 = (q[j] - qm1[j])*(q[j] - qm1[j])+1.0e-6;
            double beta2 = (qp1[j] - q[j])*(qp1[j] - q[j])+1.0e-6;

            double a1 = d1/((beta1)*(beta1));
            double a2 = d2/((beta2)*(beta2));

            double w1 = a1/(a1+a2);
            double w2 = a2/(a1+a2);

            double alpha_weno = 1.0;

            if (selective_weno==1)
            {
            alpha_weno = (level-temp_min_level)/(temp_max_level-temp_min_level);
            }

            double c1 = (1.0-alpha_weno)*d1 + alpha_weno*w1;
            double c2 = (1.0-alpha_weno)*d2 + alpha_weno*w2;

            qleft[j] = c1/2.0*(-qm1[j]+3.0*q[j]) + c2/2.0*(q[j]+qp1[j]);
        }

        d1 = 2.0/3.0;
        d2 = 1.0/3.0;

        for (int j=0;j<=n_eq-1;j++)
        {
            double beta1 = (qp1[j] - q[j])*(qp1[j] - q[j])+1.0e-6;
            double beta2 = (qp2[j] - qp1[j])*(qp2[j] - qp1[j])+1.0e-6;

            double a1 = d1/((beta1)*(beta1));
            double a2 = d2/((beta2)*(beta2));

            double w1 = a1/(a1+a2);
            double w2 = a2/(a1+a2);

            double alpha_weno = 1.0;

            if (selective_weno==1)
            {
            alpha_weno = (level-temp_min_level)/(temp_max_level-temp_min_level);
            }

            double c1 = (1.0-alpha_weno)*d1 + alpha_weno*w1;
            double c2 = (1.0-alpha_weno)*d2 + alpha_weno*w2;

            qright[j] = c1/2.0*(q[j]+qp1[j]) + c2/2.0*(3.0*qp1[j]-qp2[j]);
        }
}



/** \brief This function calculates local left and right state reconstructions using WENO5 \n
 *  Requires input arrays of floating point conserved variables at +- stencil locations \n
 *  Also requires 2 empty arrays for left and right state reconstructions
 *  Full mathematical setup may be found in IJNMF \n
 *  Note that option of scale-selective WENO enforces extra linearized weighting of nonlinear coefficients \n
 *  Formulated in the sense of the right cell face
 * \param double (&qm2)[n_eq]
 * \param double (&qm1)[n_eq]
 * \param double (&q)[n_eq]
 * \param double (&qp1)[n_eq]
 * \param double (&qp2)[n_eq]
 * \param double (&qp3)[n_eq]
 * \param double (&qright)[n_eq]
 * \param double (&qleft)[n_eq]
 * \return void
 *
 */
void weno5(int level, double (&qm2)[n_eq], double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq], double (&qp3)[n_eq], double (&qright)[n_eq], double (&qleft)[n_eq])
{

        double d1=1.0/10.0;
        double d2=6.0/10.0;
        double d3=3.0/10.0;

        for (int j=0;j<=n_eq-1;j++)
        {
            double beta1 = 13.0/12.0*(qm2[j]-2.0*qm1[j]+q[j])*(qm2[j]-2.0*qm1[j]+q[j]) + 1.0/4.0*(qm2[j]-4.0*qm1[j]+3.0*q[j])*(qm2[j]-4.0*qm1[j]+3.0*q[j])+1.0e-6;
            double beta2 = 13.0/12.0*(qm1[j]-2.0*q[j]+qp1[j])*(qm1[j]-2.0*q[j]+qp1[j]) + 1.0/4.0*(qm1[j]-qp1[j])*(qm1[j]-qp1[j])+1.0e-6;
            double beta3 = 13.0/12.0*(q[j]-2.0*qp1[j]+qp2[j])*(q[j]-2.0*qp1[j]+qp2[j]) + 1.0/4.0*(3.0*q[j]-4.0*qp1[j]+qp2[j])*(3.0*q[j]-4.0*qp1[j]+qp2[j])+1.0e-6;


            double a1 = d1/((beta1)*(beta1));
            double a2 = d2/((beta2)*(beta2));
            double a3 = d3/((beta3)*(beta3));

            double w1 = a1/(a1+a2+a3);
            double w2 = a2/(a1+a2+a3);
            double w3 = a3/(a1+a2+a3);

            double alpha_weno = 1.0;

            if (selective_weno==1 && refine_switch==1)
            {
            alpha_weno = (level-temp_min_level)/(temp_max_level-temp_min_level);
            }

            double c1 = (1.0-alpha_weno)*d1 + alpha_weno*w1;
            double c2 = (1.0-alpha_weno)*d2 + alpha_weno*w2;
            double c3 = (1.0-alpha_weno)*d3 + alpha_weno*w3;

            qleft[j] = c1*(1.0/3.0*qm2[j] - 7.0/6.0*qm1[j] + 11.0/6.0*q[j]) + c2*(-1.0/6.0*qm1[j] + 5.0/6.0*q[j] + 1.0/3.0*qp1[j]) + c3*(1.0/3.0*q[j] + 5.0/6.0*qp1[j] - 1.0/6.0*qp2[j]);
            //qleft - wave coming from left
        }

        d1 = 3.0/10.0;
        d2 = 6.0/10.0;
        d3 = 1.0/10.0;


        for (int j=0;j<=n_eq-1;j++)
        {
            double beta1 = 13.0/12.0*(qm1[j]-2.0*q[j]+qp1[j])*(qm1[j]-2.0*q[j]+qp1[j]) + 1.0/4.0*(qm1[j]-4.0*q[j]+3.0*qp1[j])*(qm1[j]-4.0*q[j]+3.0*qp1[j])+1.0e-6;
            double beta2 = 13.0/12.0*(q[j]-2.0*qp1[j]+qp2[j])*(q[j]-2.0*qp1[j]+qp2[j]) + 1.0/4.0*(q[j]-qp2[j])*(q[j]-qp2[j])+1.0e-6;
            double beta3 = 13.0/12.0*(qp1[j]-2.0*qp2[j]+qp3[j])*(qp1[j]-2.0*qp2[j]+qp3[j]) + 1.0/4.0*(3.0*qp1[j]-4.0*qp2[j]+qp3[j])*(3.0*qp1[j]-4.0*qp2[j]+qp3[j])+1.0e-6;


            double a1 = d1/((beta1)*(beta1));
            double a2 = d2/((beta2)*(beta2));
            double a3 = d3/((beta3)*(beta3));

            double w1 = a1/(a1+a2+a3);
            double w2 = a2/(a1+a2+a3);
            double w3 = a3/(a1+a2+a3);

            double alpha_weno = 1.0;

            if (selective_weno==1 && refine_switch==1)
            {
            alpha_weno = (level-temp_min_level)/(temp_max_level-temp_min_level);
            }

            double c1 = (1.0-alpha_weno)*d1 + alpha_weno*w1;
            double c2 = (1.0-alpha_weno)*d2 + alpha_weno*w2;
            double c3 = (1.0-alpha_weno)*d3 + alpha_weno*w3;

            qright[j] = c1*(1.0/3.0*qp1[j] + 5.0/6.0*q[j] - 1.0/6.0*qm1[j]) + c2*(-1.0/6.0*qp2[j] + 5.0/6.0*qp1[j] + 1.0/3.0*q[j]) + c3*(1.0/3.0*qp3[j] - 7.0/6.0*qp2[j] + 11.0/6.0*qp1[j]);
            //qright - wave coming from right
        }

}


/** \brief This function calculates local left and right state reconstructions using WENO6 \n
 *  Requires input arrays of floating point conserved variables at +- stencil locations \n
 *  Also requires 2 empty arrays for left and right state reconstructions
 *  Full mathematical setup may be found in IJNMF \n
 *  Note that option of scale-selective WENO enforces extra linearized weighting of nonlinear coefficients \n
 *  Formulated in the sense of the right cell face \n
 *  This generally is not effective to damp oscillations
 * \param double (&qm2)[n_eq]
 * \param double (&qm1)[n_eq]
 * \param double (&q)[n_eq]
 * \param double (&qp1)[n_eq]
 * \param double (&qp2)[n_eq]
 * \param double (&qp3)[n_eq]
 * \param double (&qright)[n_eq]
 * \param double (&qleft)[n_eq]
 * \return void
 *
 */
 void weno6(int level, double (&qm2)[n_eq], double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq], double (&qp3)[n_eq], double (&qright)[n_eq], double (&qleft)[n_eq])
{
    double d1=1.0/20.0;
    double d2=9.0/20.0;
    double d3=9.0/20.0;
    double d4=1.0/20.0;

    for (int j=0;j<=n_eq-1;j++)
    {
        double eps = 1.0e-6;
        double h = 13.0/12.0;

        double beta1 = h*(qm2[j]-2.0*qm1[j]+q[j])*(qm2[j]-2.0*qm1[j]+q[j]) + 1.0/4.0*(qm2[j]-4.0*qm1[j]+3.0*q[j])*(qm2[j]-4.0*qm1[j]+3.0*q[j])+eps;
        double beta2 = h*(qm1[j]-2.0*q[j]+qp1[j])*(qm1[j]-2.0*q[j]+qp1[j]) + 1.0/4.0*(qm1[j]-qp1[j])*(qm1[j]-qp1[j])+eps;
        double beta3 = h*(q[j]-2.0*qp1[j]+qp2[j])*(q[j]-2.0*qp1[j]+qp2[j]) + 1.0/4.0*(3.0*q[j]-4.0*qp1[j]+qp2[j])*(3.0*q[j]-4.0*qp1[j]+qp2[j])+eps;
        double beta4 = h*(qp1[j]-2.0*qp2[j]+qp3[j])*(qp1[j]-2.0*qp2[j]+qp3[j]) + 1.0/4.0*(5.0*qp1[j]-8.0*qp2[j]+3.0*qp3[j])*(5.0*qp1[j]-8.0*qp2[j]+3.0*qp3[j])+eps;

        beta2 = max(beta1,beta2);
        beta3 = max(beta2,beta3);
        beta4 = max(beta3,beta4);

        double a1 = d1/((beta1)*(beta1));
        double a2 = d2/((beta2)*(beta2));
        double a3 = d3/((beta3)*(beta3));
        double a4 = d4/((beta4)*(beta4));

        double w1 = a1/(a1+a2+a3+a4);
        double w2 = a2/(a1+a2+a3+a4);
        double w3 = a3/(a1+a2+a3+a4);
        double w4 = a4/(a1+a2+a3+a4);

        double alpha_weno = 1.0;

        if (selective_weno==1 && refine_switch==1)
        {
        alpha_weno = (level-temp_min_level)/(temp_max_level-temp_min_level);
        }

        double c1 = (1.0-alpha_weno)*d1 + alpha_weno*w1;
        double c2 = (1.0-alpha_weno)*d2 + alpha_weno*w2;
        double c3 = (1.0-alpha_weno)*d3 + alpha_weno*w3;
        double c4 = (1.0-alpha_weno)*d4 + alpha_weno*w4;

        double g = 1.0/6.0;

        double q1val = g*(2.0*qm2[j]-7.0*qm1[j]+11.0*q[j]);
        double q2val = g*(-qm1[j]+5.0*q[j]+2.0*qp1[j]);
        double q3val = g*(2.0*q[j]+5.0*qp1[j]-qp2[j]);
        double q4val = g*(11.0*qp1[j]-7.0*qp2[j]+2.0*qp3[j]);

        qleft[j] = c1*q1val + c2*q2val + c3*q3val + c4*q4val;
    }

    for (int j=0;j<=n_eq-1;j++)
    {
        double eps = 1.0e-6;
        double h = 13.0/12.0;

        double beta1 = h*(qp3[j]-2.0*qp2[j]+qp1[j])*(qp3[j]-2.0*qp2[j]+qp1[j]) + 1.0/4.0*(qp3[j]-4.0*qp2[j]+3.0*qp1[j])*(qp3[j]-4.0*qp2[j]+3.0*qp1[j])+eps;
        double beta2 = h*(qp2[j]-2.0*qp1[j]+q[j])*(qp2[j]-2.0*qp1[j]+q[j]) + 1.0/4.0*(qp2[j]-q[j])*(qp2[j]-q[j])+eps;
        double beta3 = h*(qp1[j]-2.0*q[j]+qm1[j])*(qp1[j]-2.0*q[j]+qm1[j]) + 1.0/4.0*(3.0*qp1[j]-4.0*q[j]+qm1[j])*(3.0*qp1[j]-4.0*q[j]+qm1[j])+eps;
        double beta4 = h*(q[j]-2.0*qm1[j]+qm2[j])*(q[j]-2.0*qm1[j]+qm2[j]) + 1.0/4.0*(5.0*q[j]-8.0*qm1[j]+3.0*qm2[j])*(5.0*q[j]-8.0*qm1[j]+3.0*qm2[j])+eps;

        beta2 = max(beta1,beta2);
        beta3 = max(beta2,beta3);
        beta4 = max(beta3,beta4);

        double a1 = d1/((beta1)*(beta1));
        double a2 = d2/((beta2)*(beta2));
        double a3 = d3/((beta3)*(beta3));
        double a4 = d4/((beta4)*(beta4));

        double w1 = a1/(a1+a2+a3+a4);
        double w2 = a2/(a1+a2+a3+a4);
        double w3 = a3/(a1+a2+a3+a4);
        double w4 = a4/(a1+a2+a3+a4);

        double alpha_weno = 1.0;

        if (selective_weno==1 && refine_switch==1)
        {
        alpha_weno = (level-temp_min_level)/(temp_max_level-temp_min_level);
        }

        double c1 = (1.0-alpha_weno)*d1 + alpha_weno*w1;
        double c2 = (1.0-alpha_weno)*d2 + alpha_weno*w2;
        double c3 = (1.0-alpha_weno)*d3 + alpha_weno*w3;
        double c4 = (1.0-alpha_weno)*d4 + alpha_weno*w4;

        double g = 1.0/6.0;

        double q1val = g*(2.0*qp3[j]-7.0*qp2[j]+11.0*qp1[j]);
        double q2val = g*(-qp2[j]+5.0*qp1[j]+2.0*q[j]);
        double q3val = g*(2.0*qp1[j]+5.0*q[j]-qm1[j]);
        double q4val = g*(11.0*q[j]-7.0*qm1[j]+2.0*qm2[j]);

        qright[j] = c1*q1val + c2*q2val + c3*q3val + c4*q4val;
    }

}



/** \brief This function calculates local left and right state reconstructions using MUSCL \n
 *  Requires input arrays of floating point conserved variables at +- stencil locations \n
 *  Also requires 2 empty arrays for left and right state reconstructions
 *  We are calculating face reconstructions using 14,17,18 Dr. San paper CAF \n
 *  Two types of limiters used 1 - MinMod 2 - MC \n
 *  Formulated in the sense of the right cell face
 * \param double (&qm1)[n_eq]
 * \param double (&q)[n_eq]
 * \param double (&qp1)[n_eq]
 * \param double (&qp2)[n_eq]
 * \param double (&qright)[n_eq]
 * \param double (&qleft)[n_eq]
 * \return void
 *
 */
void muscl(double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq], double (&qright)[n_eq], double (&qleft)[n_eq])
{

    if (mlim==1)//Minmod limiter
    {
        double r[n_eq] = {0.0};
        double phi[n_eq] = {0.0};

        for (int j=0;j<=n_eq-1;j++)
        {
            double qrightrightval = qp2[j];
            double qrightval = qp1[j];
            double qval = q[j];
            double qleftval = qm1[j];

             if (fabs(qrightrightval-qrightval)>1.0e-12)
             {
                r[j] = (qrightval-qval)/(qrightrightval-qrightval);
                phi[j]= max(0.0,min(r[j],1.0));
             }
             else
             {
                 phi[j] = 1.0;
             }

             qright[j] = qrightval - 0.5*phi[j]*(qrightrightval-qrightval);


             if (fabs(qrightval-qval)>1.0e-12)
             {
                r[j] = (qval-qleftval)/(qrightval-qval);
                phi[j]= max(0.0,min(r[j],1.0));
             }
             else
             {
                 phi[j] = 1.0;
             }

             qleft[j] = qval + 0.5*phi[j]*(qrightval-qval);
        }
    }
    else if (mlim==2)//MC Limiter
    {
        double r[n_eq] = {0.0};
        double phi[n_eq] = {0.0};

        for (int j=0;j<=n_eq-1;j++)
        {
            double qrightrightval = qp2[j];
            double qrightval = qp1[j];
            double qval = q[j];
            double qleftval = qm1[j];

             if (fabs(qrightrightval-qrightval)>1.0e-12)
             {
                r[j] = (qrightval-qval)/(qrightrightval-qrightval);
                phi[j]= max(0.0,min(2.0*r[j],min(0.5*(1.0+r[j]),2.0)));
             }
             else
             {
                 phi[j] = 2.0;
             }

             qright[j] = qrightval - 0.5*phi[j]*(qrightrightval-qrightval);


             if (fabs(qrightval-qval)>1.0e-12)
             {
                r[j] = (qval-qleftval)/(qrightval-qval);
                phi[j]= max(0.0,min(2.0*r[j],min(0.5*(1.0+r[j]),2.0)));
             }
             else
             {
                 phi[j] = 2.0;
             }

             qleft[j] = qval + 0.5*phi[j]*(qrightval-qval);
        }

    }
}


/** \brief This function corrects fluxes to ensure conservativeness  \n
 * This is to ensure that overlapping cell boundaries have same incoming and outgoing fluxes \n
 * Refer section 2.6 and equation 2.42 for details in Tenaud tutorial
 * \param unordered_map<int,Cell*> &leafvector
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void conservative_correction(unordered_map<int,Cell*> &leafvector,unordered_map<int,Cell*> &Cellvect)
{
    int level = max_level-1;

    while(level>=0)
    {
        for (auto citer_one=leafvector.begin();citer_one!=leafvector.end();citer_one++)
        {
            if ((*citer_one).second->virt==0 && (*citer_one).second->level == level)
            {
                //Check left side
                int i = (*citer_one).second->i;

                if (find_cell_exist(2*i-1,Cellvect) && find_cell(2*i-1,Cellvect)->leaf==1 && find_cell(2*i-1,Cellvect)->virt==0 && find_cell(2*i-1,Cellvect)->level==level+1)
                {
                    for (int j=0;j<=n_eq-1;j++)
                    {
                        (*citer_one).second->f_m[j] = find_cell(2*i-1,Cellvect)->f_p[j];
                    }
                }
                //Check right side
                if (find_cell_exist(2*i+2,Cellvect) && find_cell(2*i+2,Cellvect)->leaf==1 && find_cell(2*i+2,Cellvect)->virt==0 && find_cell(2*i+2,Cellvect)->level==level+1)
                {
                    for (int j=0;j<=n_eq-1;j++)
                    {
                        (*citer_one).second->f_p[j] = find_cell(2*i+2,Cellvect)->f_m[j];
                    }
                }
            }
        }
        level--;
    }
}





#endif // face_fluxes_h_inluded
