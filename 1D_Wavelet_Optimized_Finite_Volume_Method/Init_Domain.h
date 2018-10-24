#ifndef Init_Domain_h_inluded
#define Init_Domain_h_inluded
/** \file Init_Domain.h
 *  \brief Header file for initializing problem cases and datastructure \n
*/

/** \brief This creates the initial domain\n
* A new cell class instance is created for root and added to unordered_map Cellvect\n
* The root is kept at i=1 and level = 0\n
* The root is at xx=0.5 and has a cell length of 1.0\n
* The root is not a leaf\n
* Children to root are added to unordered_map\n
* Leaves are marked and physical values are embedded in them\n
* Linkages are established between children, parents and neighbors
* \param unordered_map<key,Cell*,hasher_key> &Cellvect
* \return void
*
*/
void create_domain(unordered_map<int,Cell*> &Cellvect)
{
    //Creating root of this particular domain
    Cell *nc = new Cell;
    nc->i = 1;
    nc->level = 0;
    nc->cell_length = 1.0;
    nc->leaf = 0;

    nc->xx = 0.5;

    Cellvect[1]=nc;

    //Adding other instances
    add_init_children(Cellvect);

    //Checking leaves
    check_leaves(Cellvect);

    //Initial physical values to finest cells
    leaf_phys_vals(Cellvect);

    //Link Neighbors
    link_level_neighbors(Cellvect);

}



/** \brief This adds children to roots and cascades downwards successively\n
* Note that the first element of map must be i=1\n
* At the addition of each children a new cell class is instantiated\n
* Children are also linked successively with their parents according to graded logic
* \param unordered_map<key,Cell*,hasher_key> &Cellvect
* \return void
*
*/
void add_init_children(unordered_map<int,Cell*> &Cellvect)
{
    int level = 1;

    int init_max_level=max_level;

    while(level<=init_max_level)
    {
        int upper_limit = int(pow(2,level+1))-1;
        int lower_limit = int(pow(2,level));

        for (int i=lower_limit;i<=upper_limit;i++)
        {
                Cell *nc = new Cell;
                nc->i = i;
                nc->level = level;


                nc->cell_length = 1.0/pow(2,level);


                nc->leaf = 0;
                nc->xx = 0.0;

                nc->parent = Cellvect[i/2];

                if (i%2==0)//left child
                {
                    nc->parent->left_child = nc;
                    nc->xx = nc->parent->xx - (nc->cell_length)/2.0;
                }
                else if (i%2!=0)//right child
                {
                    nc->parent->right_child = nc;
                    nc->xx = nc->parent->xx + (nc->cell_length)/2.0;
                }
                Cellvect[i]=nc;
        }
    level++;
    }

}

/** \brief Initial condition for physical quantities at leaves\n
* Note that only leaves have actual physical quantities which are projected backward\n
* The initial conditions (as well as number of equations) are problem specific.\n
* \param unordered_map<key,Cell*,hasher_key> &Cellvect
* \return void
*
*/
void leaf_phys_vals(unordered_map<int,Cell*> &Cellvect)
{
    for (auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
    {
        if ((*citer_one).second->leaf == 1)
        {
            double x = (*citer_one).second->xx;

            if (n_eq==1)
            {
                //double pi = 4.0*atan(1.0);
                //(*citer_one).second->q[0] = sin(pi*x*2.0);
                //(*citer_one).second->q[0] = exp(-50.0*(x-0.5)*(x-0.5));
                if (x>=0.5 && x<0.6)
                {
                    (*citer_one).second->q[0] = 1.0;
                }
                else
                {
                    (*citer_one).second->q[0] = 0.0;
                }
            }
            else if (n_eq==2)
            {
                if (x<=0.5)
                {
                    (*citer_one).second->q[0] = height_1;
                    (*citer_one).second->q[1] = 0.0;
                }
                else
                {
                    (*citer_one).second->q[0] = height_2;
                    (*citer_one).second->q[1] = 0.0;
                }
            }
            else if (n_eq==3)//SOD Shock tube
            {
                if (x<=0.5)
                {
                    (*citer_one).second->q[0] = 1.0;
                    (*citer_one).second->q[1] = 0.0;
                    (*citer_one).second->q[2] = 1.0/(gmm-1.0);
                }
                else
                {
                    (*citer_one).second->q[0] = 0.125;
                    (*citer_one).second->q[1] = 0.0;
                    (*citer_one).second->q[2] = 0.1/(gmm-1.0);
                }
            }
            else if (n_eq==7)
            {
                if (x<=0.5)
                {
                    (*citer_one).second->q[0] = 1.0;
                    (*citer_one).second->q[1] = 0.0;
                    (*citer_one).second->q[2] = 0.0;
                    (*citer_one).second->q[3] = 0.0;
                    (*citer_one).second->q[4] = 1.0;
                    (*citer_one).second->q[5] = 0.0;
                    (*citer_one).second->q[6] = 1.0/(gmm-1.0) + 0.5*(bx*bx + 1.0);
                }
                else
                {
                    (*citer_one).second->q[0] = 0.125;
                    (*citer_one).second->q[1] = 0.0;
                    (*citer_one).second->q[2] = 0.0;
                    (*citer_one).second->q[3] = 0.0;
                    (*citer_one).second->q[4] = -1.0;
                    (*citer_one).second->q[5] = 0.0;
                    (*citer_one).second->q[6] = 0.1/(gmm-1.0) + 0.5*(bx*bx + 1.0);//Refer http://www.csun.edu/~jb715473/examples/mhd1d.htm
                }
            }
        }
    }
}







#endif
