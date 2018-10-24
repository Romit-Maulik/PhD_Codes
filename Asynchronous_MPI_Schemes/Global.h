#ifndef Global_h_inluded
#define Global_h_inluded

using namespace std;
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <fstream>      //myfile.open()
#include <iomanip>      // std::setw

class Cell
{
    public:

    double u_now = 0.0;
    double u_prev = 0.0;
    double dudt_prev = 0.0;
    double u_exact = 0.0;
    double x = 0.0;
    double u_temp = 0.0;

    //Error quantification
    double sync_error = 0.0; //Error with no modification
    double async_error = 0.0;//Error after equations are modified

    int proc_num = 0;
    int index = 0;

    Cell* leftpe = nullptr;
    Cell* rightpe = nullptr;

    //Parameterized Constructor
    Cell(int a, int b)
    {
        index = a;
        proc_num = b;
    };
    //Deconstructor
    ~Cell(){};
};

double update_stencil(double uleft, double umid, double uright, double alpha, double c, double dt, double dx);
void plot_output(vector<Cell>&, int );


#endif // Global_h_inluded
