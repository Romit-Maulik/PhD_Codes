#ifndef Global_h_inluded
#define Global_h_inluded
/** \file Global.h
 *  \brief Global header file containing Cell class \n
 *  Also contains include directives for some standard libraries
*/

using namespace std;

//Headers for global use
#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <fstream>      //myfile.open()
#include <iomanip>      // std::setw
#include <memory>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>

/** \class Cell
 *  \brief This is the object for a finite volume cell.\n
 *  Each cell has a unique i,j value to locate it in memory\n
 *  Also initialized with multiple pointers to parents, neighbors and children\n
 *  Refer further documentation related to each instance of storage,pointer etc.
 */
class Cell
{
public:


    double q[n_eq] = {0.0};/**< Storage for conserved variables*/
    int i;/**< Index in hierarchical data representation*/

    //MRA based storage
    double det[n_eq]={0.0};/**< Storage for wavelet coefficients*/

    int level;/**< Current level of resolution*/
    double xx=0.0;/**< Cell-center location in physical space*/
    int keep_flag = 1;/**< Flag to decide if cell is retained*/
    double cell_length=0.0;/**< Cell-width*/
    int leaf = 0;/**< Flag [0]-Not a leaf, [1] Leaf*/
    int new_cell = 0;/**< Flag [0]-Not a newly created cell, [1] Newly created cell*/

    Cell* parent=nullptr;/**< Pointer to parent*/

    Cell* left_child=nullptr;/**< Pointer to left child 2i */
    Cell* right_child=nullptr;/**< Pointer to right child 2i+1 */

    Cell* left_level=nullptr;/**< Pointer to left neighbor i-1 */
    Cell* right_level=nullptr;/**< Pointer to right neighbor i+1 */

    int virt = 0;/**< Flag [0]-Not a virtual leaf, [1] Virtual leaf*/

    //Face fluxes

    double f_p[n_eq] = {0.0};/**< Storage for right face flux*/
    double f_m[n_eq] = {0.0};/**< Storage for left face flux*/

     //Runge Kutta storage
    double q_rhs[n_eq] = {0.0};/**< Runge-Kutta Storage*/
    double q_temp[n_eq] = {0.0};/**< Runge-Kutta Storage*/
    double q_1[n_eq] = {0.0};/**< Runge-Kutta Storage*/
    double q_2[n_eq] = {0.0};/**< Runge-Kutta Storage*/

    //Source Term
    double q_source[n_eq] = {0.0};/**< Source Term*/

    //Perona Malik Storage
    double filt_q[n_eq]={0.0};/**< Perona-Malik Storage*/

    //Constructor
    Cell(){};
    //Deconstructor
    ~Cell(){};

};


#endif // Global_h_inluded
