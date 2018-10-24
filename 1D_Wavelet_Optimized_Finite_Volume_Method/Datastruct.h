#ifndef Datastruct_h_inluded
#define Datastruct_h_inluded
/** \file Datastruct.h
 *  \brief Header file for datastructure related functions \n
*/

/** \brief Marks leaves by sweeping through unordered_map \n
 *  Logic for leaf flag switch is if a cell has children or not.
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void check_leaves(unordered_map<int,Cell*> &Cellvect)
{
    for (auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
    {
        int iref = (*citer_one).second->i;
        (*citer_one).second->leaf = 0;

        if (!find_cell_exist(2*iref+1,Cellvect) && !find_cell_exist(2*iref,Cellvect))//No Children
        {
            (*citer_one).second->leaf = 1;
        }
    }
}



/** \brief Logical function to determine if a particular cell exists \n
 *  Takes in an index i (unique cell index) \n
 *  Logical function to find if cell exists by path traversal from root to cell \n
 *  Returns boolean output \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \param i
 * \return bool
 *
 */
bool find_cell_exist(int i, unordered_map<int,Cell*> &Cellvect)
{
    //Note that first element of Cellvect has to be the root.
    int child_tracker[max_level];
    int level;

    if (i==1)//Root always remains in vector
    {
       return true;
    }


    for (level=0;level<max_level;level++)//Increment after body of loop executes
    {
        if (i%2==0)//Left child
        {
            child_tracker[level]=1;
        }
        else if (i%2!=0)//Right Child
        {
            child_tracker[level]=2;
        }

        i = i/2;

        if (i==1)//We have reached root
        {
            break;
        }
    }

    if (i!=1)
    {
        return false;
    }

    Cell* loc = Cellvect[1];

    do{
        if (child_tracker[level]==1)// left child
        {

            if (loc->left_child == nullptr)//Cell does not exist
            {
                return false;
            }
            else
            {
                loc = loc->left_child;
            }
        }
        else if (child_tracker[level]==2)//Right child
        {
            if (loc->right_child == nullptr)//Cell does not exist
            {
                return false;
            }
            else
            {
                loc = loc->right_child;
            }
        }
        level--;

    }while(level>=0);

    return true;//If it has made it here that means the cell exists

}


/** \brief Function returns pointer to a particular cell with given i value \n
 *  Takes in an index i (unique cell index) \n
 *  Returns Cell* to this i value \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \param i
 * \return Cell*
 *
 */
Cell* find_cell(int i, unordered_map<int,Cell*> &Cellvect)//Use after verifying that cell exists
{
    //Note that first element of Cellvect has to be the root.
    int child_tracker[max_level];
    int level;

    if (i==1)//Root always remains in vector
    {
        Cell* loc = Cellvect[1];
        return loc;
    }


    for (level=0;level<max_level;level++)//Increment after body of loop executes
    {
        if (i%2==0)//left child
        {
            child_tracker[level]=1;
        }
        else if (i%2!=0)// right child
        {
            child_tracker[level]=2;
        }

        i = i/2;

        if (i==1)//We have reached root
        {
            break;
        }
    }

    Cell* loc=Cellvect[1];

    do{

        if (child_tracker[level]==1)//left child
        {
                loc = loc->left_child;
        }
        else if (child_tracker[level]==2)//right child
        {
                loc = loc->right_child;
        }
        level--;
    }while(level>=0);

    return loc;
}



/** \brief Function links to neighboring level if cell exists \n
 *  Required for interpolations \n
 *  Varies at boundaries according to periodic or open boundary conditions
 * \param unordered_map<int,Cell*> &Cellvect
 * \return void
 *
 */
void link_level_neighbors(unordered_map<int,Cell*> &Cellvect)
{

    for (auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
    {
        int iright = (*citer_one).second->i + 1;
        int ileft = (*citer_one).second->i - 1;

        //Linking right cells within domain
        if (find_cell_exist(iright,Cellvect))
        {
            (*citer_one).second->right_level = find_cell(iright,Cellvect);
        }

        //Linking left cells within domain
        if (find_cell_exist(ileft,Cellvect))
        {
            (*citer_one).second->left_level = find_cell(ileft,Cellvect);
        }

        //Linking right direction boundary condition
        if (iright == pow(2,(*citer_one).second->level+1))
        {
            if (bc_type==1)//Periodic
            {
                int per_iright = pow(2,(*citer_one).second->level);
                if (find_cell_exist(per_iright,Cellvect))
                {
                    (*citer_one).second->right_level = find_cell(per_iright,Cellvect);
                }
            }
            else//Open
            {
                (*citer_one).second->right_level = (*citer_one).second;
            }

        }

        //Linking left direction boundary condition
        if (ileft == pow(2,(*citer_one).second->level)-1)
        {
            if (bc_type==1)
            {
                int per_ileft = pow(2,(*citer_one).second->level+1)-1;
                if (find_cell_exist(per_ileft,Cellvect))
                {
                    (*citer_one).second->left_level = find_cell(per_ileft,Cellvect);
                }
            }
            else
            {
                (*citer_one).second->left_level = (*citer_one).second;
            }

        }
    }
}

/** \brief Function either changes existent cell keep_flag to 1 or makes new cell \n
 *  Also sets new_cell flag as 1 (even if cell was already present) \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \param i
 * \return void
 *
 */
void keep_cell(int i, int level, unordered_map<int,Cell*> &Cellvect)
{


    if (bc_type==1)
    {
        //First to ensure that periodicity is respected for i direction
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = pow(2,level) + i - pow(2,level+1);
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = pow(2,level+1) - pow(2,level) + i;
        }
    }
    else
    {
        //Open BC
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = int(pow(2,level+1))-1;
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = int(pow(2,level));
        }
    }


    //Periodicity is taken care of now

    if (find_cell_exist(i,Cellvect))
    {
        Cell* loc = find_cell(i,Cellvect);
        loc->keep_flag = 1;
    }
    else
    {
        int child_tracker[max_level];
        int level;

        for (level=0;level<max_level;level++)//Charting route of cell from root
        {
        if (i%2==0)//left child
        {
            child_tracker[level]=1;
        }
        else if (i%2!=0)//Right child
        {
            child_tracker[level]=2;
        }

        i = i/2;

        if (i==1)//We have reached root
        {
            break;
        }
        }

        Cell* loc = Cellvect[1];
        int insert_level = 1;

        do{

            if (child_tracker[level]==1)//Lower left child
            {

                if (loc->left_child == nullptr)//Cell does not exist
                {
                    Cell* nc = new Cell;

                    nc->i = 2*i;
                    nc->level = insert_level;
                    nc->parent = loc;
                    nc->cell_length = nc->parent->cell_length/2.0;
                    nc->xx = nc->parent->xx - (nc->cell_length)/2.0;;

                    loc->left_child = nc;
                    nc->new_cell = 1;
                    nc->keep_flag = 1;

                    Cellvect[2*i]=nc;//Insertion into unique key

                    loc = nc;

                    i = 2*i;

                }
                else
                {
                    loc = loc->left_child;
                    i = 2*i;
                }


            }
            else if (child_tracker[level]==2)//Right child
            {
                if (loc->right_child == nullptr)//Cell does not exist
                {
                    Cell* nc = new Cell;
                    nc->i = 2*i+1;

                    nc->level = insert_level;
                    nc->parent = loc;
                    nc->cell_length = nc->parent->cell_length/2.0;
                    nc->xx = nc->parent->xx + (nc->cell_length)/2.0;

                    loc->right_child = nc;
                    nc->new_cell = 1;
                    nc->keep_flag = 1;

                    Cellvect[2*i+1]=nc;//Insertion into unique key

                    loc = nc;
                    i = 2*i+1;

                }
                else
                {
                    loc = loc->right_child;
                    i = 2*i+1;
                }
            }
            level--;
            insert_level++;
        }while(level>=0);
    }
}

/** \brief Function returns an unordered_map of memory addresses to Cells of the same level \n
 *  Useful for level operations which are incremented or decremented \n
 *  Note that changes on this unordered_map also affect global binary tree \n
 *  These are just pointers(or memory addresses) \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \param level
 * \return unordered_map<Cell*>
 *
 */
unordered_map<int,Cell*> level_vector(int level, unordered_map<int,Cell*> &Cellvect)
{
    unordered_map<int,Cell*>levelvector;

    if (level==0)
    {
        Cell* loc = Cellvect[1];
        levelvector[1] = loc;
        return levelvector;
    }

    int low_limit = int(pow(2,level));
    int high_limit = int(pow(2,level+1))-1;

    for (int i=low_limit;i<=high_limit;i++)
    {
            if (find_cell_exist(i,Cellvect))
            {
            Cell* loc = find_cell(i,Cellvect);
            levelvector[i] = loc;
            }
    }
    return levelvector;
}

/** \brief Function returns an unordered_map of memory addresses to Cells which are leaves \n
 *  Our governing law essentially evolves on this unordered_map \n
 *  Note that changes on this unordered_map also affect global binary tree \n
 *  These are just pointers(or memory addresses) \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \param level
 * \return unordered_map<Cell*>
 *
 */
unordered_map<int,Cell*> leaf_vector(unordered_map<int,Cell*> &Cellvect)
{
    unordered_map<int,Cell*>leafvector;
    for(auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
    {
        if((*citer_one).second->leaf==1)
        {
            int i = (*citer_one).second->i;
            Cell* loc = find_cell(i,Cellvect);
            leafvector[i] = loc;
        }
    }
    return leafvector;
}


/** \brief Function makes a new virtual cell for conservative flux calculations \n
 *  Also sets new_cell flag as 1 and virt flag as 1 \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \param i
 * \return void
 *
 */
void keep_cell_virt(int i, int level, unordered_map<int,Cell*> &Cellvect)
{

    if (bc_type==1)
    {
        //First to ensure that periodicity is respected for i direction
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = pow(2,level) + i - pow(2,level+1);
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = pow(2,level+1) - pow(2,level) + i;
        }
    }
    else
    {
        //Open BC
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = int(pow(2,level+1))-1;
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = int(pow(2,level));
        }
    }

    //Periodicity is taken care of now

    if (find_cell_exist(i,Cellvect))
    {
        Cell* loc = find_cell(i,Cellvect);
        loc->keep_flag = 1;
    }
    else
    {
        int child_tracker[max_level];
        int level;

        for (level=0;level<max_level;level++)//Charting route of cell from root
        {
        if (i%2==0)//left child
        {
            child_tracker[level]=1;
        }
        else if (i%2!=0)//right child
        {
            child_tracker[level]=2;
        }

        i = i/2;

        if (i==1)//We have reached root
        {
            break;
        }
        }

        Cell* loc = Cellvect[1];
        int insert_level = 1;
        int virt_tracker = 0;

        do{
            if (child_tracker[level]==1)//left child
            {

                if (loc->left_child == nullptr)//Cell does not exist
                {
                    Cell* nc = new Cell;

                    nc->i = 2*i;

                    nc->level = insert_level;
                    nc->parent = loc;
                    nc->cell_length = nc->parent->cell_length/2.0;
                    nc->xx = nc->parent->xx - (nc->cell_length)/2.0;;

                    loc->left_child = nc;
                    nc->new_cell = 1;
                    nc->keep_flag = 1;

                    Cellvect[2*i] = nc;

                    loc = nc;
                    i = 2*i;
                    virt_tracker++;

                }
                else
                {
                    loc = loc->left_child;
                    i = 2*i;
                }
            }
            else if (child_tracker[level]==2)//Upper left child
            {
                if (loc->right_child == nullptr)//Cell does not exist
                {
                    Cell* nc = new Cell;

                    nc->i = 2*i+1;

                    nc->level = insert_level;
                    nc->parent = loc;
                    nc->cell_length = nc->parent->cell_length/2.0;
                    nc->xx = nc->parent->xx + (nc->cell_length)/2.0;;

                    loc->right_child = nc;
                    nc->new_cell = 1;
                    nc->keep_flag = 1;
                    Cellvect[2*i+1] = nc;

                    loc = nc;
                    i = 2*i+1;
                    virt_tracker++;
                }
                else
                {
                    loc = loc->right_child;
                    i = 2*i+1;
                }
            }

            level--;
            insert_level++;

        }while(level>=0);

        if (virt_tracker!=0)
        {
        find_cell(i,Cellvect)->virt = 1;
        find_cell(i,Cellvect)->leaf = 1;
        find_cell(i,Cellvect)->keep_flag = 1;
        }

        find_cell(i,Cellvect)->keep_flag = 1;

    }
}


/** \brief Function deletes an unnecessary cell \n
 *  Also severs pointers for parent and neighbors \n
 *  Does not delete if cell as has children as it must be in graded structure \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \param i
 * \return void
 *
 */
void delete_cell(int i, unordered_map<int,Cell*> &Cellvect)
{
    Cell* loc = find_cell(i,Cellvect);
    //Severing pointer connection with parent

    if (i%2==0)//left child
    {
        loc->parent->left_child = nullptr;
    }
    else if (i%2!=0)//right child
    {
        loc->parent->right_child = nullptr;
    }

    //Delinking neighbor connections
    if (loc->left_level != nullptr)
    {
        loc->left_level->right_level = nullptr;
    }

    if (loc->right_level != nullptr)
    {
        loc->right_level->left_level = nullptr;
    }

    Cellvect.erase(i);

    delete(loc);
}


/** \brief Function updates index values to reflect boundary conditions \n
 *  Adapts i such that boundary is either periodic or open \n
 *  Connects boundary node to itself for Open BC \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \param i
 * \return void
 *
 */
int periodize(int i, int level)
{
    if (bc_type==1)
    {
        //First to ensure that periodicity is respected for i direction
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = pow(2,level) + i - pow(2,level+1);
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = pow(2,level+1) - pow(2,level) + i;
        }
    }
    else
    {
        //Open BC
        if (i>int(pow(2,level+1))-1) //This is beyond the right boundary
        {
            i = int(pow(2,level+1))-1;
        }
        else if (i<int(pow(2,level))) //This is beyond the left boundary
        {
            i = int(pow(2,level));
        }
    }

    return i;
}

#endif // Datastruct_h_inluded
