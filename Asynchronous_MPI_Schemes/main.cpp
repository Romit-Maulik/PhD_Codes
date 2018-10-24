#include "Global.h"
#include "Functions.h"
#include <random>


int main()
{
    vector<Cell>domain;

    int mod = 3; //1-Modified, 2- Unmodified, 3-Our idea

    int nx = 128;
    int nprocs = 4;
    double delay_prob = 1.0;
    double xlength = 2.0*atan(1.0)*4.0;
    double dx = xlength/double(nx);
    double Kappa = 2.0;


    double rand_val = (double)rand() / (double)RAND_MAX;
    double angle = 20.0 * (rand_val - 0.5) * 2.0;
    initialize(nx,dx,nprocs,domain,Kappa,angle);


    double final_time = 1.0;
    double r_alpha = 0.1;
    double alpha = 0.1;
    double dt = r_alpha*dx*dx/alpha;
    double c = 1.0;
    double r_cb2 = c*dt/(2.0*dx);

    async_time_advance(mod,nx,nprocs,dx,domain,alpha,c,r_cb2,r_alpha,dt,final_time,delay_prob,Kappa,angle);
    plot_output(domain,nx);
    avg_error(domain,nx);

    return 0;
}
