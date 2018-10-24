#ifndef Source_terms_h_inluded
#define Source_terms_h_inluded

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
/// \fn    void calculate_sources(double q[],double (&f)[n_eq]);
/// \brief Takes in two arrays by reference, q->conserved variables, s->sources from them
/// \brief Has conditional statements for each system of equations
/// \brief Source term is assumed to be on the right hand side of the system of equations
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
void calculate_sources(unordered_map<int,Cell*> &leafvector)
{
//    if (n_eq==4)
//    {
//        for (auto citer_one=leafvector.begin();citer_one!=leafvector.end();citer_one++)
//        {
//            if ((*citer_one).second->virt == 0)//Only for true leaves
//            {
//                double u_plus_one = (*citer_one).second->right_level->q[1]/(*citer_one).second->right_level->q[0];
//                double u_plus_two = (*citer_one).second->right_level->right_level->q[1]/(*citer_one).second->right_level->right_level->q[0];
//                double u_plus_three = (*citer_one).second->right_level->right_level->right_level->q[1]/(*citer_one).second->right_level->right_level->right_level->q[0];
//                double u_minus_one = (*citer_one).second->q[1]/(*citer_one).second->q[0];
//                double u_minus_two = (*citer_one).second->left_level->q[1]/(*citer_one).second->left_level->q[0];
//                double u_minus_three = (*citer_one).second->left_level->left_level->q[1]/(*citer_one).second->left_level->left_level->q[0];
//
//                double u_rf = 0.0;
//
//                for (int j=0;j<=n_eq-1;j++)
//                {   //Refer Table 1 Hyman Paper
//                    u_rf = (u_minus_three - 8.0*u_minus_two + 37.0*u_minus_one + 37.0*u_plus_one - 8.0*u_plus_two + 1.0*u_plus_three)/60.0;
//                }
//
//                u_plus_one = (*citer_one).second->q[1]/(*citer_one).second->q[0];
//                u_plus_two = (*citer_one).second->right_level->q[1]/(*citer_one).second->right_level->q[0];
//                u_plus_three = (*citer_one).second->right_level->right_level->q[1]/(*citer_one).second->right_level->right_level->q[0];
//                u_minus_one = (*citer_one).second->left_level->q[1]/(*citer_one).second->left_level->q[0];
//                u_minus_two = (*citer_one).second->left_level->left_level->q[1]/(*citer_one).second->left_level->left_level->q[0];
//                u_minus_three = (*citer_one).second->left_level->left_level->left_level->q[1]/(*citer_one).second->left_level->left_level->left_level->q[0];
//
//                double u_lf = 0.0;
//
//                for (int j=0;j<=n_eq-1;j++)
//                {   //Refer Table 1 Hyman Paper
//                    u_lf = (u_minus_three - 8.0*u_minus_two + 37.0*u_minus_one + 37.0*u_plus_one - 8.0*u_plus_two + 1.0*u_plus_three)/60.0;
//                }
//
//                double dudx = (u_rf - u_lf)/(*citer_one).second->cell_length;
//
//                (*citer_one).second->q_source[3] = ((*citer_one).second->q[3]+4.0/3.0*ep_mu*dudx);
//            }
//        }
//    }
}





















#endif
