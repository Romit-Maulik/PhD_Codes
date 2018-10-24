#ifndef Wavelet_h_inluded
#define Wavelet_h_inluded

/** \file Wavelet.h
 *  \brief Header file for wavelet coefficient calculation \n
 * This file also determines projection, prediction and thresholding operations.
*/

/** \brief This is a control panel for calculation of wavelet coefficients \n
 *  Calls two functions for projection and detail calculation respectively \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void encode(unordered_map<int,Cell*> &Cellvect)
{
    projection(Cellvect);

    details_calculation(Cellvect);
}


/** \brief Function enables conservative projection from leaves towards root \n
 *  Calculates average cell value of variable at parent \n
 *  Refer Tenaud tutorial equation 2.28 \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void projection(unordered_map<int,Cell*> &Cellvect)
{
    int level = max_level-1;

    while(level>=0)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level--;
            continue;
        }

        for (auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
            if ((*citer_one).second->leaf==0)
            {
                double current_cell_area = (*citer_one).second->cell_length;
                double child_cell_area = (*citer_one).second->left_child->cell_length;

                for (int j=0;j<=n_eq-1;j++)
                {
                    double left_q = (*citer_one).second->left_child->q[j];
                    double right_q = (*citer_one).second->right_child->q[j];

                    double local_q = child_cell_area/current_cell_area*(left_q+right_q);
                    (*citer_one).second->q[j] = local_q;
                }
            }
        }
        levelvector.clear();
        level--;
    }
}


/** \brief Function uses projection to reinterpolate at finer levels \n
 *  Stores wavelet coefficients as difference of true values and interpolated ones \n
 *  Uses step 6 in Algorithm 1 Tenaud tutorial \n
 *  Note that there are no details calculated for level 0 (or root since it has no parents) \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void details_calculation(unordered_map<int,Cell*> &Cellvect)
{
    int level = max_level;

    while(level>=1)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level--;
            continue;
        }

        for (auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {

                int i = (*citer_one).second->i;
                double xi1 = -22.0/128.0;
                double xi2 = 3.0/128.0;

                for (int j=0;j<=n_eq-1;j++)
                {
                    double parent_val = (*citer_one).second->parent->q[j];

                    double parent_val_right2 = (*citer_one).second->parent->right_level->right_level->q[j];
                    double parent_val_right1 = (*citer_one).second->parent->right_level->q[j];
                    double parent_val_left1 = (*citer_one).second->parent->left_level->q[j];
                    double parent_val_left2 = (*citer_one).second->parent->left_level->left_level->q[j];

                    double interp_val=0.0;

                    if (i%2==0)//left child
                    {
                        interp_val = parent_val+(xi1*(parent_val_right1-parent_val_left1))+(xi2*(parent_val_right2-parent_val_left2));
                    }
                    else if (i%2!=0)//right child
                    {
                        interp_val = parent_val-(xi1*(parent_val_right1-parent_val_left1))-(xi2*(parent_val_right2-parent_val_left2));
                    }

                    (*citer_one).second->det[j] =  fabs((*citer_one).second->q[j] - interp_val);
                }


        }
        level--;
        levelvector.clear();
    }
}

/** \brief Function preserves binary graded structure to ensure appropriate linkages across resolutions \n
 *  This subroutine must be used after marking cells for deletion or retention
 *  Any cell that must be kept has its parent stencil retained \n
 *  This follows steps 13-25: Algorithm 3 Tenaud tutorial \n
 *  Assumes tree is already graded at the end of previous refinement \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void preserve_graded_structure(unordered_map<int,Cell*> &Cellvect)
{
    int level = max_level;

    while(level>1)//Parent of level 1 is root which is always included
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level--;
            continue;
        }

        for(auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
            if (  (*citer_one).second->keep_flag == 1 )
            {
                int i = (*citer_one).second->i;
                int level_ref = (*citer_one).second->level;

                //Keeping parent
                keep_cell(i/2,level_ref-1,Cellvect);

                //Keeping left_level neighbor of parent
                keep_cell(i/2-1,level_ref-1,Cellvect);

                //Keeping right_level neighbor of parent
                keep_cell(i/2+1,level_ref-1,Cellvect);

                //Keeping left_left_level neighbor of parent
                keep_cell(i/2-2,level_ref-1,Cellvect);

                //Keeping right_right_level neighbor of parent
                keep_cell(i/2+2,level_ref-1,Cellvect);
            }
        }
        level--;
        levelvector.clear();//This makes sure Cell* pointers etc are not corrupted
    }
}



/** \brief Function thresholds datastructure in order to bound perturbation error theoretically \n
 *  Used to determine if a cell must be kept in the tree structure or not
 *  In addition, adds virtual cells and decodes them to ensure stencil retention for fluxes \n
 *  This function is key in ensuring error scales linearly with wavelet coefficient \n
 *  Refer IJNMF paper for linear scaling verification
 *  Algorithm 5 in Tenaud tutorial \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void hartens_predictive_thresholding(unordered_map<int,Cell*> &Cellvect)
{
    int level = max_level - 1;

    for(auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)//Step 3 in Tenaud tutorial Algorithm 5
    {
        (*citer_one).second->keep_flag = 0;
        (*citer_one).second->virt = 0;
        (*citer_one).second->leaf = 0;
        (*citer_one).second->new_cell = 0;
    }

    Cellvect[1]->keep_flag = 1;//Root must stay

//    double global_max[n_eq] = {0.0};
//    find_max_vals(Cellvect,global_max);

    while(level>0)
    {

        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level--;
            continue;
        }

        if (n_eq==1)
        {
            for(auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
            {
                //double local_det = (*citer_one).second->det[0]/global_max[0];
                double local_det = (*citer_one).second->det[0];
                double local_threshold = wc_threshold/pow(2,max_level-level);

                if (local_det>=local_threshold)
                {
                    int iref = (*citer_one).second->i;
                    int level_ref = (*citer_one).second->level;

                    //Keeping left_child
                    keep_cell(2*iref,level_ref+1,Cellvect);
                    //Keeping left_child->left_level
                    keep_cell(2*iref-1,level_ref+1,Cellvect);
                   //Keeping right_child
                    keep_cell(2*iref+1,level_ref+1,Cellvect);
                    //Keeping right_child->right_level
                    keep_cell(2*iref+2,level_ref+1,Cellvect);

                    if (local_det>=local_threshold*pow(2,2.0*p_w) && level!=max_level-1)//Refer step 12 of Algo 5
                    {
                            for (int i=4*iref-2;i<=4*iref+5;i++)
                            {
                                keep_cell(i,level_ref+2,Cellvect);
                            }
                    }
                }
            }
        }
        else
        {
            double max_detail = find_max_detail(levelvector);

            for(auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
            {

                double local_det = 0.0;

                for (int j=0;j<=n_eq-1;j++)
                {
                    local_det = local_det+fabs((*citer_one).second->det[j]);
                    //local_det = max(local_det,fabs((*citer_one).second->det[j]/global_max[j]));
                }

                double local_threshold = wc_threshold/pow(2,max_level-level);

                if (local_det/max_detail>local_threshold)
                //if (local_det>=local_threshold)
                {

                    int iref = (*citer_one).second->i;
                    int level_ref = (*citer_one).second->level;

                    //Keeping left_child
                    keep_cell(2*iref,level_ref+1,Cellvect);
                    //Keeping left_child->left_level
                    keep_cell(2*iref-1,level_ref+1,Cellvect);
                   //Keeping right_child
                    keep_cell(2*iref+1,level_ref+1,Cellvect);
                    //Keeping right_child->right_level
                    keep_cell(2*iref+2,level_ref+1,Cellvect);

                    if (local_det/max_detail>=local_threshold*pow(2,2.0*p_w) && level!=max_level-1)//Refer step 12 of Algo 5
                    //if (local_det>=local_threshold*pow(2,2.0*p_w) && level!=max_level-1)//Refer step 12 of Algo 5
                    {
                            for (int i=4*iref-2;i<=4*iref+5;i++)
                            {
                                keep_cell(i,level_ref+2,Cellvect);
                            }
                    }
                }
            }
        }

        level--;
        levelvector.clear();
    }
    preserve_graded_structure(Cellvect);
    conservative_refinement(Cellvect);
    link_level_neighbors(Cellvect);
    decode_new_cells(Cellvect);
    delete_unnecessary_cells(Cellvect);
    check_leaves(Cellvect);
    add_virtual_cells(Cellvect);
    preserve_graded_structure(Cellvect);
    conservative_refinement(Cellvect);
    link_level_neighbors(Cellvect);
    decode_new_cells(Cellvect);
}

/** \brief Function finds the maximum value of the wavelet coefficient for the given unordered map \n
 *  Generally fed a level_vector and not the entire Cellvect unordered_map \n
 *  Needed for Harten's predictive thresholding \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return double
 *
 */
double find_max_detail(unordered_map<int,Cell*> &Cellvect)
{
    double max_detail=0.0;

    if (n_eq!=1)
    {
        for(auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
        {
            double local_det = 0.0;
            for (int j=0;j<=n_eq-1;j++)
            {
                local_det = local_det+fabs((*citer_one).second->det[j]);
            }

            max_detail = max(local_det,max_detail);
        }
    }
    else
    {
        for(auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
        {
            double l1sum = (*citer_one).second->det[0];
            max_detail = max(fabs(l1sum),max_detail);
        }
    }

    return max_detail;

}


/** \brief Function finds the maximum value of each conserved variable globally \n
 *  Can be fed Cellvect or Levelvect unordered_maps \n
 *  Needed for Harten's predictive thresholding \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \param double (&global_max)[n_eq]
 * \return void
 *
 */
void find_max_vals(unordered_map<int,Cell*> &Cellvect,double (&global_max)[n_eq])
{
    for (int j=0;j<=n_eq-1;j++)
    {
        global_max[j] = 0.0;

        for (auto citer_one = Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
        {
            global_max[j] = max(global_max[j],fabs((*citer_one).second->q[j]));
        }
    }
}



/** \brief Function used to ensure that post thresholding the existence of one child ensures existence of all children \n
 *  For example - If any cell has atleast one child - it should have all 2 children. \n
 *  Needed for Harten's predictive thresholding \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void conservative_refinement(unordered_map<int,Cell*> &Cellvect)
{
    int level = 0;

    while(level<=max_level-1)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level++;
            continue;
        }

        for(auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
            int iref = (*citer_one).second->i;
            int level_ref = (*citer_one).second->level;

            if (find_cell_exist(2*iref,Cellvect))//left
            {
                if (find_cell(2*iref,Cellvect)->keep_flag == 1)
                {
                    keep_cell(2*iref,level_ref+1,Cellvect);
                    keep_cell(2*iref+1,level_ref+1,Cellvect);
                }
            }

            if (find_cell_exist(2*iref+1,Cellvect))//right
            {
                if (find_cell(2*iref+1,Cellvect)->keep_flag == 1)
                {
                    keep_cell(2*iref,level_ref+1,Cellvect);
                    keep_cell(2*iref+1,level_ref+1,Cellvect);
                }
            }
        }
    level++;
    levelvector.clear();
    }
}


/** \brief Function adds virtual cells to ensure stability for advective solution and conservative flux calculation \n
 *  Refer step 3(v) in Algorithm 6 of Tenaud tutorial. \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void add_virtual_cells(unordered_map<int,Cell*> &Cellvect)
{
    int level = 0;

    while(level<=max_level)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level++;
            continue;
        }

        for (auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
                if ((*citer_one).second->leaf==1)
                {
                int iref = (*citer_one).second->i;
                int lev = (*citer_one).second->level;
                int ileft = periodize(iref-1,lev);
                int iright = periodize(iref+1,lev);


                if (level<max_level)
                {
                    if (!find_cell_exist(iright,Cellvect))//Right Direction
                    {
                        keep_cell_virt(iref+1,level,Cellvect);
                        keep_cell_virt(iref+2,level,Cellvect);
                        keep_cell_virt(iref+3,level,Cellvect);
                        keep_cell_virt(iref+4,level,Cellvect);
                    }
                    else if (find_cell(iright,Cellvect)->leaf == 0)
                    {
                        keep_cell_virt(2*iref,level+1,Cellvect);
                        keep_cell_virt(2*iref+1,level+1,Cellvect);
                        keep_cell_virt(2*iref-1,level+1,Cellvect);
                        keep_cell_virt(2*iref-2,level+1,Cellvect);
                    }

                    if (!find_cell_exist(ileft,Cellvect))//Left Direction
                    {
                        keep_cell_virt(iref-1,level,Cellvect);
                        keep_cell_virt(iref-2,level,Cellvect);
                        keep_cell_virt(iref-3,level,Cellvect);
                        keep_cell_virt(iref-4,level,Cellvect);
                    }
                    else if (find_cell(ileft,Cellvect)->leaf == 0)
                    {
                        keep_cell_virt(2*iref,level+1,Cellvect);
                        keep_cell_virt(2*iref+1,level+1,Cellvect);
                        keep_cell_virt(2*iref+2,level+1,Cellvect);
                        keep_cell_virt(2*iref+3,level+1,Cellvect);
                    }
                }
                else
                {
                    if (!find_cell_exist(iright,Cellvect))//Right Direction
                    {
                        keep_cell_virt(iref+1,level,Cellvect);
                        keep_cell_virt(iref+2,level,Cellvect);
                        keep_cell_virt(iref+3,level,Cellvect);
                        keep_cell_virt(iref+4,level,Cellvect);
                    }

                    if (!find_cell_exist(ileft,Cellvect))//Left Direction
                    {
                        keep_cell_virt(iref-1,level,Cellvect);
                        keep_cell_virt(iref-2,level,Cellvect);
                        keep_cell_virt(iref-3,level,Cellvect);
                        keep_cell_virt(iref-4,level,Cellvect);
                    }
                }
                }
        }
        level++;
        levelvector.clear();
    }
}



/** \brief Function deletes cell which are deemed to be superfluous by Harten's predictive and conservative thresholding \n
 *  Deletes all cells with keep_flag zero. Simple. \n
 *  Must be used after all thresholding, virtual cells, conservative refinement etc \n
 *  Note Tenaud uses this before adding virtual cells - we do it after. Theoretically no difference (as shown by scaling)
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void delete_unnecessary_cells(unordered_map<int,Cell*> &Cellvect)
{
    int level = max_level;
    while(level>0)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level--;
            continue;
        }

        for(auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
            if ((*citer_one).second->keep_flag == 0)
            {
                int iref = (*citer_one).second->i;
                delete_cell(iref,Cellvect);
            }
        }
        level--;
        levelvector.clear();
    }
}


/** \brief Function interpolates to find solutions at newly created cells \n
 *  Linkages must be complete (and tree must be graded prior to this) \n
 *  Refer algorithm 2 in Tenaud tutorial \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void decode_new_cells(unordered_map<int,Cell*> &Cellvect)
{
    int level = 1;

    while(level<=max_level)
    {
        unordered_map<int,Cell*>levelvector;
        levelvector = level_vector(level,Cellvect);

        if (levelvector.size()==0)
        {
            level++;
			continue;
        }

        for (auto citer_one=levelvector.begin();citer_one!=levelvector.end();citer_one++)
        {
            if ((*citer_one).second->new_cell == 1)
            {
                int i = (*citer_one).second->i;
                double xi1 = -22.0/128.0;
                double xi2 = 3.0/128.0;

                //double parent_cell_length = (*citer_one).second->parent->cell_length;
                //double cell_length = (*citer_one).second->cell_length;

                for (int j=0;j<=n_eq-1;j++)
                {
                    double interp_val=0.0;

                    if (i%2==0)//left child
                    {
                        double parent_val = (*citer_one).second->parent->q[j];
                        double parent_val_right2 = (*citer_one).second->parent->right_level->right_level->q[j];
                        double parent_val_right1 = (*citer_one).second->parent->right_level->q[j];
                        double parent_val_left1 = (*citer_one).second->parent->left_level->q[j];
                        double parent_val_left2 = (*citer_one).second->parent->left_level->left_level->q[j];

                        interp_val = parent_val+(xi1*(parent_val_right1-parent_val_left1))+(xi2*(parent_val_right2-parent_val_left2));
                        interp_val = interp_val + (*citer_one).second->parent->det[j];
                    }
                    else if (i%2!=0)//right child
                    {
                        //interp_val = parent_val-(xi1*(parent_val_right1-parent_val_left1))-(xi2*(parent_val_right2-parent_val_left2));

                        double parent_val = (*citer_one).second->left_level->parent->q[j];
                        double parent_val_right2 = (*citer_one).second->left_level->parent->right_level->right_level->q[j];
                        double parent_val_right1 = (*citer_one).second->left_level->parent->right_level->q[j];
                        double parent_val_left1 = (*citer_one).second->left_level->parent->left_level->q[j];
                        double parent_val_left2 = (*citer_one).second->left_level->parent->left_level->left_level->q[j];

                        interp_val = parent_val+(xi1*(parent_val_right1-parent_val_left1))+(xi2*(parent_val_right2-parent_val_left2));
                        interp_val = interp_val + (*citer_one).second->parent->det[j];

                        interp_val = 2.0*parent_val - interp_val;

                    }

                    (*citer_one).second->q[j] = interp_val;
                }

                (*citer_one).second->new_cell = 0;
            }
        }
        level++;
        levelvector.clear();
    }
}


/** \brief Function interpolates to find solutions at virtual cells at each substage of RK3 \n
 *  This is our addition above what Tenaud did (he does not use multistep time integrators) \n
 *  Virtual leaves (where the solution does not evolve) must reflect update of solution on true leaves \n
 *  Solution from true leaves is encoded towards root and cascaded downwards into virtual leaves \n
 *  Note that this does not affect final solution much - not used in IJNMF
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void substage_decode(unordered_map<int,Cell*> &Cellvect)
{

    encode(Cellvect);

    unordered_map<int,Cell*>leafvector;
    leafvector = leaf_vector(Cellvect);

    for (auto citer_one=leafvector.begin();citer_one!=leafvector.end();citer_one++)
    {
        if ((*citer_one).second->virt == 1)
        {
            int i = (*citer_one).second->i;
                double xi1 = -22.0/128.0;
                double xi2 = 3.0/128.0;

//                double parent_cell_length = (*citer_one).second->parent->cell_length;
//                double cell_length = (*citer_one).second->cell_length;

                for (int j=0;j<=n_eq-1;j++)
                {
                    double interp_val=0.0;

                    if (i%2==0)//left child
                    {
                        double parent_val = (*citer_one).second->parent->q[j];
                        double parent_val_right2 = (*citer_one).second->parent->right_level->right_level->q[j];
                        double parent_val_right1 = (*citer_one).second->parent->right_level->q[j];
                        double parent_val_left1 = (*citer_one).second->parent->left_level->q[j];
                        double parent_val_left2 = (*citer_one).second->parent->left_level->left_level->q[j];

                        interp_val = parent_val+(xi1*(parent_val_right1-parent_val_left1))+(xi2*(parent_val_right2-parent_val_left2));
                        interp_val = interp_val + (*citer_one).second->parent->det[j];
                    }
                    else if (i%2!=0)//right child
                    {
                        //interp_val = parent_val-(xi1*(parent_val_right1-parent_val_left1))-(xi2*(parent_val_right2-parent_val_left2));

                        double parent_val = (*citer_one).second->left_level->parent->q[j];
                        double parent_val_right2 = (*citer_one).second->left_level->parent->right_level->right_level->q[j];
                        double parent_val_right1 = (*citer_one).second->left_level->parent->right_level->q[j];
                        double parent_val_left1 = (*citer_one).second->left_level->parent->left_level->q[j];
                        double parent_val_left2 = (*citer_one).second->left_level->parent->left_level->left_level->q[j];

                        interp_val = parent_val+(xi1*(parent_val_right1-parent_val_left1))+(xi2*(parent_val_right2-parent_val_left2));
                        interp_val = interp_val + (*citer_one).second->parent->det[j];

                        interp_val = 2.0*parent_val - interp_val;

                    }

                    (*citer_one).second->q[j] = interp_val;// + (*citer_one).second->det[j];
                }
        }
    }
}


#endif // Wavelet_h_inluded
