#ifndef Functions_h_inluded
#define Functions_h_inluded

void initialize(int nx, double dx, int nprocs, vector<Cell>&Myvect, double Kappa, double angle)
{
    int pe_res = nx/nprocs;
    int thres = 0;
    int proc = 0;

    //Populating vector
    for (int iter=0;iter<nx;iter++)
    {
        thres++;
        if (thres>pe_res)
        {
            proc++;
            thres = 0;
        }
        Cell object = Cell(iter,proc);
        Myvect.push_back(object);//This is our vector for the heat equation with Cell objects

    }

    //Connecting neighbors
    for (auto it=Myvect.begin()+1;it!=Myvect.end()-1;it++)
    {
        (*it).leftpe = &(*(it-1));
        (*it).rightpe = &(*(it+1));
    }

    //Boundaries
    (*Myvect.begin()).leftpe = &(*(Myvect.end()-1));
    (*Myvect.begin()).rightpe = &(*(Myvect.begin()+1));

    (*(Myvect.end()-1)).rightpe = &(*(Myvect.begin()));
    (*(Myvect.end()-1)).leftpe = &(*(Myvect.end()-2));


    //Initial conditions
    double pi = atan(1.0)*4.0;

    for (int iter=0;iter<nx;iter++)
    {
        Myvect[iter].x = double(iter)*dx;
        Myvect[iter].u_exact = sin(Kappa*double(iter)*dx+pi*angle/180.0);
        Myvect[iter].u_now = sin(Kappa*double(iter)*dx+pi*angle/180.0);
        Myvect[iter].u_prev = sin(Kappa*double(iter)*dx+pi*angle/180.0);
    }

}

void async_time_advance(int mod, int nx, int nprocs, double dx, vector<Cell>&Myvect, double alpha, double c, double r_cb2, double r_alpha, double dt, double final_time, double delay_prob, double Kappa, double angle)
{
    //Time loop
    double t = 0.0;

    while(t<final_time)
    {
        t = t+dt;
        vector<int>proc_del;

        for (int i=0;i<nprocs;i++)
        {
            double rand_val = (double)rand() / (double)RAND_MAX;
            if (rand_val<=delay_prob)
            {
                proc_del.push_back(1);
            }
            else
            {
                proc_del.push_back(0);
            }
        }

        if (mod==1)
        {
            for (int i=0;i<nx;i++)
            {
                if (Myvect[i].leftpe->proc_num != Myvect[i].proc_num)//Left Boundary
                {
                    double async_fac = r_alpha + r_cb2;
                    int kdel = 0;
                    kdel = proc_del[Myvect[i].leftpe->proc_num];
                    double c_new = c/(1.0-kdel*async_fac);
                    double alpha_new = alpha/(1.0-kdel*async_fac);

                    if (kdel==1)
                    {
                        Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_prev,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha_new,c_new,dt,dx);
                    }
                    else
                    {
                        Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha_new,c_new,dt,dx);
                    }


                }
                else if (Myvect[i].rightpe->proc_num != Myvect[i].proc_num)//Right boundary
                {
                    double async_fac = r_alpha - r_cb2;
                    int kdel = 0;
                    kdel = proc_del[Myvect[i].rightpe->proc_num];
                    double c_new = c/(1.0-kdel*async_fac);
                    double alpha_new = alpha/(1.0-kdel*async_fac);

                    if (kdel==1)
                    {
                        Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,Myvect[i].rightpe->u_prev,alpha_new,c_new,dt,dx);
                    }
                    else
                    {
                        Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha_new,c_new,dt,dx);
                    }
                }
                else
                {
                    Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha,c,dt,dx);
                }

            }

            //Updating the exact solution and calculating errors
            for (int i=0;i<nx;i++)
            {
                double x = Myvect[i].x;
                double pi = atan(1.0)*4.0;
                Myvect[i].u_prev = Myvect[i].u_now;
                Myvect[i].u_now = Myvect[i].u_temp;
                Myvect[i].u_exact = exp(-t*alpha*Kappa*Kappa)*sin(Kappa*(x-c*t)+pi*angle/180.0);
                Myvect[i].async_error = fabs(Myvect[i].u_exact-Myvect[i].u_now);
            }
        }
        else if (mod==2)
        {
            for (int i=0;i<nx;i++)
            {
                if (Myvect[i].leftpe->proc_num != Myvect[i].proc_num)//Left Boundary
                {
                    int kdel = 0;
                    kdel = proc_del[Myvect[i].leftpe->proc_num];

                    if (kdel==1)
                    {
                        Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_prev,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha,c,dt,dx);
                    }
                    else
                    {
                        Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha,c,dt,dx);
                    }


                }
                else if (Myvect[i].rightpe->proc_num != Myvect[i].proc_num)//Right boundary
                {
                    int kdel = 0;
                    kdel = proc_del[Myvect[i].rightpe->proc_num];

                    if (kdel==1)
                    {
                        Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,Myvect[i].rightpe->u_prev,alpha,c,dt,dx);
                    }
                    else
                    {
                        Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha,c,dt,dx);
                    }
                }
                else
                {
                    Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha,c,dt,dx);
                }

            }

            //Updating the exact solution and calculating errors
            for (int i=0;i<nx;i++)
            {
                double x = Myvect[i].x;
                double pi = atan(1.0)*4.0;
                Myvect[i].u_prev = Myvect[i].u_now;
                Myvect[i].u_now = Myvect[i].u_temp;
                Myvect[i].u_exact = exp(-t*alpha*Kappa*Kappa)*sin(Kappa*(x-c*t)+pi*angle/180.0);
                Myvect[i].async_error = fabs(Myvect[i].u_exact-Myvect[i].u_now);
            }
        }
        else
        {
            for (int i=0;i<nx;i++)
            {
                if (Myvect[i].leftpe->proc_num != Myvect[i].proc_num)//Left Boundary
                {
                    int kdel = 0;
                    kdel = proc_del[Myvect[i].leftpe->proc_num];

                    if (kdel==1)
                    {
                        double um0 = Myvect[i].u_prev;
                        double um1 = Myvect[i].leftpe->u_prev;
                        double um2 = Myvect[i].leftpe->leftpe->u_prev;
                        double umodleft = um1 + dt*(-c/dx*(um0-um2)+alpha/(dx*dx)*(um0-2.0*um1+um2));

                        Myvect[i].u_temp = update_stencil(umodleft,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha,c,dt,dx);
                    }
                    else
                    {
                        Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha,c,dt,dx);
                    }


                }
                else if (Myvect[i].rightpe->proc_num != Myvect[i].proc_num)//Right boundary
                {

                    int kdel = 0;
                    kdel = proc_del[Myvect[i].rightpe->proc_num];

                    if (kdel==1)
                    {
                        double up0 = Myvect[i].u_prev;
                        double up1 = Myvect[i].rightpe->u_prev;
                        double up2 = Myvect[i].rightpe->rightpe->u_prev;
                        double umodright = up1 + dt*(-c/dx*(up2-up0)+alpha/(dx*dx)*(up2-2.0*up1+up0));

                        Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,umodright,alpha,c,dt,dx);
                    }
                    else
                    {
                        Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha,c,dt,dx);
                    }
                }
                else
                {
                    Myvect[i].u_temp = update_stencil(Myvect[i].leftpe->u_now,Myvect[i].u_now,Myvect[i].rightpe->u_now,alpha,c,dt,dx);
                }

            }

            //Updating the exact solution and calculating errors
            for (int i=0;i<nx;i++)
            {
                double x = Myvect[i].x;
                double pi = atan(1.0)*4.0;
                Myvect[i].u_prev = Myvect[i].u_now;
                Myvect[i].u_now = Myvect[i].u_temp;
                Myvect[i].u_exact = exp(-t*alpha*Kappa*Kappa)*sin(Kappa*(x-c*t)+pi*angle/180.0);
                Myvect[i].async_error = fabs(Myvect[i].u_exact-Myvect[i].u_now);
            }
        }
    }

}

double update_stencil(double uleft, double umid, double uright, double alpha, double c, double dt, double dx)
{
    double unew = umid + dt*((uleft-uright)*c)/(2.0*dx) + dt*(uright-2.0*umid+uleft)*alpha/(dx*dx);
    return unew;
}

void plot_output(vector<Cell>&Myvect, int nx)
{
    ofstream myfile;
    myfile.open("Check.plt");
    myfile << "variables =\"x\",\"u\",\"Exact\" \n";

    for(int i=0; i<nx; i++)
    {
    myfile << fixed << std::setprecision(10) << Myvect[i].x << setw(20);
    myfile << fixed << std::setprecision(10) << Myvect[i].u_now << setw(20);
    myfile << fixed << std::setprecision(10) << Myvect[i].u_exact << endl;
    }
    myfile.close();
}


void avg_error(vector<Cell>&Myvect, int nx)
{
    ofstream myfile;
    myfile.open("Error.plt");
    myfile << "variables =\"x\",\"Error\" \n";

    for(int i=0; i<nx; i++)
    {
    myfile << fixed << std::setprecision(10) << Myvect[i].x << setw(20);
    myfile << fixed << std::setprecision(10) << fabs(Myvect[i].u_now-Myvect[i].u_exact) << endl;
    }
    myfile.close();
}

#endif // Functions_h_inluded
