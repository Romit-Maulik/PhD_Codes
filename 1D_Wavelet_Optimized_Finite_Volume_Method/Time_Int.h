#ifndef Time_Int_h_inluded
#define Time_Int_h_inluded

/** \file Time_Int.h
 *  \brief Header file for Runge-Kutta third order time integration \n
*/


/** \brief This function is the control panel for time integration \n
 *  Also outputs metrics for the size of datastructure with time \n
 *  Also calls functions from Output.h for periodic recording of spatial field. \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void evolution(unordered_map<int,Cell*> &Cellvect)
{
    double t=0.0;
    int t_snap = 1;

    string pt = "Sizes.plt";
    std::ofstream myfile;
    myfile.open(pt);
    myfile << "variables =\"t\",\"length\" \n";

    int file_op=0;
    while(t<final_time)
    {
        file_op++;

        if (refine_switch==1)
        {
            encode(Cellvect);
            hartens_predictive_thresholding(Cellvect);
        }

        unordered_map<int,Cell*>leafvect;
        leafvect = leaf_vector(Cellvect);

        temp_min_level = find_min_level(leafvect);
        temp_max_level = find_max_level(leafvect);

        if (riesol==1)
        {
            double u = calc_max_wave_speed(leafvect);
            double min_cell_size = 1.0/(pow(2,max_level));
            dt = cfl*min_cell_size/u;
        }
        else
        {
            dt = force_dt;
        }

        if (t+dt>final_time)
        {
            dt = fabs(final_time-t);
        }

        runge_kutta_3(leafvect,Cellvect);

        t = t+dt;
        cout<<t<<endl;

        if (t>=t_snap*final_time/snap_steps)
        {
        string filename = NumtoStr(t_snap) + ".plt";
        print_field(filename,Cellvect);
        t_snap++;
        }

        myfile << std::fixed << std::setprecision(10) << t << std::setw(20);
        myfile << std::fixed << std::setprecision(10) << Cellvect.size() << std::endl;


    }

    myfile.close();
}



/** \brief This function converts an integer to string \n
 *  Required for filename modification in periodic output \n
 * \param int a
 * \return string
 *
 */
string NumtoStr (int a)
{
    ostringstream temp;
    temp<<a;
    return temp.str();
}



/** \brief Function for total-variation-diminishing Runge-Kutta third-order ODE integrator \n
 *  Also outputs metrics for the size of datastructure with time \n
 *  Note extra steps related to conservative correction and substage decoding of physical values into virtual leaves
 * \param unordered_map<int,Cell*> &Cellvect
 * \param unordered_map<int,Cell*> &leafvector
 * \return void
 *
 */
void runge_kutta_3(unordered_map<int,Cell*> &leafvector, unordered_map<int,Cell*> &Cellvect)
{

    calculate_monotonic_fluxes(leafvector);//This calculates fluxes  monotonically near shocks
    conservative_correction(leafvector,Cellvect);

    for (auto citer_one = leafvector.begin();citer_one != leafvector.end(); citer_one++)
    {
        if ((*citer_one).second->virt == 0)
        {
            for (int j=0;j<=n_eq-1;j++)
            {
                double q_now = (*citer_one).second->q[j];

                (*citer_one).second->q_temp[j] = q_now;
                (*citer_one).second->q_rhs[j] = - 1.0/((*citer_one).second->cell_length)*((*citer_one).second->f_p[j] - (*citer_one).second->f_m[j]);
                (*citer_one).second->q_rhs[j] = (*citer_one).second->q_rhs[j] + (*citer_one).second->q_source[j];
                (*citer_one).second->q_1[j] = q_now + dt*((*citer_one).second->q_rhs[j]);
                (*citer_one).second->q[j] = (*citer_one).second->q_1[j];
            }
        }
    }

    //substage_decode(Cellvect);//Update virtual cell values
    calculate_monotonic_fluxes(leafvector);//This calculates fluxes  monotonically near shocks
    conservative_correction(leafvector,Cellvect);

    for (auto citer_one = leafvector.begin();citer_one != leafvector.end(); citer_one++)
    {
        if ((*citer_one).second->virt == 0)
        {
            for (int j=0;j<=n_eq-1;j++)
            {
                (*citer_one).second->q_rhs[j] = - 1.0/((*citer_one).second->cell_length)*((*citer_one).second->f_p[j] - (*citer_one).second->f_m[j]);
                (*citer_one).second->q_rhs[j] = (*citer_one).second->q_rhs[j] + (*citer_one).second->q_source[j];
                (*citer_one).second->q_2[j] = 0.75*(*citer_one).second->q_temp[j] + 0.25*(*citer_one).second->q_1[j] + 0.25*dt*((*citer_one).second->q_rhs[j]);
                (*citer_one).second->q[j] = (*citer_one).second->q_2[j];
            }
        }
    }

    //substage_decode(Cellvect);//Update virtual cell values
    calculate_monotonic_fluxes(leafvector);//This calculates fluxes  monotonically near shocks
    conservative_correction(leafvector,Cellvect);

    for (auto citer_one = leafvector.begin();citer_one != leafvector.end(); citer_one++)
    {
        if ((*citer_one).second->virt == 0)
        {
            for (int j=0;j<=n_eq-1;j++)
            {
                (*citer_one).second->q_rhs[j] = - 1.0/((*citer_one).second->cell_length)*((*citer_one).second->f_p[j] - (*citer_one).second->f_m[j]);
                (*citer_one).second->q_rhs[j] = (*citer_one).second->q_rhs[j] + (*citer_one).second->q_source[j];
                (*citer_one).second->q[j] = 1.0/3.0*(*citer_one).second->q_temp[j] + 2.0/3.0*(*citer_one).second->q_2[j] + 2.0/3.0*dt*((*citer_one).second->q_rhs[j]);
                (*citer_one).second->q_rhs[j] = 0.0;
                (*citer_one).second->q_1[j] = 0.0;
                (*citer_one).second->q_2[j] = 0.0;
                (*citer_one).second->q_temp[j] = 0.0;
            }
        }
    }

    substage_decode(Cellvect);
    for (auto citer_one = leafvector.begin();citer_one != leafvector.end(); citer_one++)
    {
        if ((*citer_one).second->virt==1)
        {
            (*citer_one).second->virt=0;
        }
    }

    //check_leaves(Cellvect);
}


/** \brief Function finds the minimum level of resolution among leaves i.e. the least refined level \n
 * \param unordered_map<int,Cell*> &leafvector
 * \return int
 *
 */
int find_min_level(unordered_map<int,Cell*> &leafvector)
{
    int min_level = max_level;

    for (auto citer_one=leafvector.begin();citer_one!=leafvector.end();citer_one++)
    {
        min_level = min((*citer_one).second->level,min_level);
    }

    return min_level;

}

/** \brief Function finds the maximum level of resolution among leaves i.e. the most refined level \n
 * \param unordered_map<int,Cell*> &leafvector
 * \return int
 *
 */
int find_max_level(unordered_map<int,Cell*> &leafvector)
{
    int temp_max_level = 0;

    for (auto citer_one=leafvector.begin();citer_one!=leafvector.end();citer_one++)
    {
        temp_max_level = max((*citer_one).second->level,temp_max_level);
    }

    return temp_max_level;

}

#endif // Time_Int_h_inluded
