#ifndef Output_h_inluded
#define Output_h_inluded

/** \file Output.h
 *  \brief Header file for Tecplot output \n
*/

/** \brief This function outputs given datastructure contents in Tecplot format \n
 *  The output is the true leaves and the physical solutions within these cells \n
 * \param unordered_map<int,Cell*> &Cellvect
 * \param string filename
 * \return void
 *
 */
void print_field(string filename, unordered_map<int,Cell*> &Cellvect)
{
    if (n_eq==1)
    {
        std::ofstream myfile;
        myfile.open(filename);
        myfile << "variables =\"x\",\"q1\",\"level\",\"fp\",\"fm\",\"i\",\"det\" \n";

        for (auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
        {
            if ((*citer_one).second->leaf == 1 && (*citer_one).second->virt==0)
            {
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->xx << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[0] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->level << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->f_p[0] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->f_m[0] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->i << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->det[0] << std::endl;
            }
        }

        myfile.close();
    }
    else if (n_eq==2)
    {
        std::ofstream myfile;
        myfile.open(filename);
        myfile << "variables =\"x\",\"q1\",\"q2\",\"level\",\"det\" \n";

        for (auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
        {
            if ((*citer_one).second->leaf == 1 && (*citer_one).second->virt==0)
            {
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->xx << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[0] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[1] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->level << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->det[0] << std::endl;
            }
        }
        myfile.close();
    }
    else if (n_eq==3)
    {
        std::ofstream myfile;
        myfile.open(filename);
        myfile << "variables =\"x\",\"q1\",\"q2\",\"q3\",\"i\",\"level\",\"det\" \n";

        for (auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
        {
            if ((*citer_one).second->leaf == 1 && (*citer_one).second->virt==0)
            {
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->xx << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[0] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[1]<< std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[2]<< std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->i<< std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->level << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->det[0] << std::endl;
            }
        }

        myfile.close();
    }
    else if (n_eq==7)
    {
        std::ofstream myfile;
        myfile.open(filename);
        myfile << "variables =\"x\",\"q1\",\"q2\",\"q3\",\"q4\",\"q5\",\"q6\",\"q7\",\"level\",\"det\" \n";

        for (auto citer_one=Cellvect.begin();citer_one!=Cellvect.end();citer_one++)
        {
            if ((*citer_one).second->leaf == 1 && (*citer_one).second->virt==0)
            {
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->xx << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[0] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[1] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[2] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[3] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[4] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[5] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->q[6] << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->level << std::setw(20);
                myfile << std::fixed << std::setprecision(10) << (*citer_one).second->det[0] << std::endl;
            }
        }

        myfile.close();
    }
}


#endif
