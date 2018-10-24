#ifndef Static_Eval_h_inluded
#define Static_Eval_h_inluded

/** \file Static_Eval.h
 *  \brief Header file for analyzing discrete wavelet transform implementation and accuracy for static problem \n
 *  Refer IJNMF for Sine-Wave perturbation demo \n
 *  Can be activated by uncommenting the analyse_static function in main.cpp (and hiding create_domain and evolution functions there)
*/

/** \brief This creates the initial domain for a static error check \n
 *  Also performs thresholding, deletion etc etc to check L1 norm behavior at a snapshot in time \n
 *  Error computed against analytical solution. \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void analyse_static(unordered_map<int,Cell*> &Cellvect)
{
    //Creating root of this particular domain
    Cell *nc = new Cell;
    nc->i = 1;
    nc->level = 0;//Root is kept level 0
    nc->cell_length = 1.0;
    nc->leaf = 0;

    nc->xx = 0.5;

    Cellvect[1]=nc;

    //Adding other instances
    add_init_children(Cellvect);

    //Checking leaves
    check_leaves(Cellvect);

    //Initial physical values to finest cells
    leaf_phys_vals_static(Cellvect);

    //Link Neighbors
    link_level_neighbors(Cellvect);

    //Encode
    encode(Cellvect);

    //Threshold
    hartens_predictive_thresholding(Cellvect);

    //print output
    print_field("Static_Check.plt",Cellvect);

    //Printing L1 norm
    unordered_map<int,Cell*>leaves;
    leaves = leaf_vector(Cellvect);
    double l1 = 0.0;
    double linf = 0.0;

    for (auto citer_one=leaves.begin();citer_one!=leaves.end();citer_one++)
    {
        if ((*citer_one).second->virt==0)
        {
            double x = (*citer_one).second->xx;
            double pi = 4.0*atan(1.0);
            double exact_val = sin(2.0*pi*(x-0.5)) + exp(-10000.0*(x-0.5)*(x-0.5));

            double diff = fabs(exact_val-(*citer_one).second->q[0]);

            l1 = l1 + diff;
            linf = max(linf,diff);
        }
    }

    cout<<"The L1 norm is "<<l1/double(Cellvect.size())<<endl;
    cout<<"The L infinity norm is "<<linf<<endl;

    //Calculating compression ratio
    int compressed_size = Cellvect.size();
    int fine_size = int(pow(2,max_level));

    double cr = 100.0*double(compressed_size)/double(fine_size);

    cout<<"The compression ratio is "<<cr<<endl;

}

/** \brief This specifies initial conditions for the static error check idea  \n
 *  As usual physical values are specified on the leaves only \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void leaf_phys_vals_static(unordered_map<int,Cell*> &Cellvect)
{
    for (auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
    {
        if ((*citer_one).second->leaf == 1)
        {
            double x = (*citer_one).second->xx;
            double pi = 4.0*atan(1.0);
            (*citer_one).second->q[0] = sin(2.0*pi*(x-0.5)) + exp(-10000.0*(x-0.5)*(x-0.5));
        }
    }
}



#endif // Static_Eval_h_inluded
