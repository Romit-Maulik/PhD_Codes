#ifndef Function_Def_h_inluded
#define Function_Def_h_inluded

/** \file Function_Def.h
 *  \brief This header file is essentially used to declare all the functions in this code.
 *
 */


//Init_Domain header
void add_init_children(unordered_map<int,Cell*> &);
void leaf_phys_vals(unordered_map<int,Cell*> &);

//Datastruct header
void check_leaves(unordered_map<int,Cell*> &);
bool find_cell_exist(int , unordered_map<int,Cell*> &);
Cell* find_cell(int , unordered_map<int,Cell*> &);
void link_level_neighbors(unordered_map<int,Cell*> &);
void keep_cell(int , int, unordered_map<int,Cell*> &);
void keep_cell_virt(int, int , unordered_map<int,Cell*> &);
unordered_map<int,Cell*> level_vector(int, unordered_map<int,Cell*> &);
void delete_cell(int, unordered_map<int,Cell*> &);
unordered_map<int,Cell*> leaf_vector(unordered_map<int,Cell*> &);
int periodize(int, int );


//Wavelet header
void encode(unordered_map<int,Cell*> &);
void projection(unordered_map<int,Cell*> &);
void details_calculation(unordered_map<int,Cell*> &);
void preserve_graded_structure(unordered_map<int,Cell*> &);
double find_max_detail(unordered_map<int,Cell*> &);
void conservative_refinement(unordered_map<int,Cell*> &);
void add_virtual_cells(unordered_map<int,Cell*> &);
void delete_unnecessary_cells(unordered_map<int,Cell*> &);
void decode_new_cells(unordered_map<int,Cell*> &);
void substage_decode(unordered_map<int,Cell*> &);
bool thresholding_a(double (&q)[n_eq], double (&q2)[n_eq]);
bool thresholding_b(double (&q)[n_eq], double (&q2)[n_eq]);
void find_max_vals(unordered_map<int,Cell*> &Cellvect,double (&global_max)[n_eq]);
void hartens_predictive_thresholding(unordered_map<int,Cell*> &);
void wamr_stencil(unordered_map<int,Cell*> &,unordered_map<int,Cell*> &);

//Output header
void print_field(string,unordered_map<int,Cell*> &);


//Time Integration header
void advance_in_time(unordered_map<int,Cell*> &);//Only for leaves
void evolution(unordered_map<int,Cell*> &);
string NumtoStr (int );
void runge_kutta_3(unordered_map<int,Cell*> &, unordered_map<int,Cell*> &);
int find_min_level(unordered_map<int,Cell*> &);
int find_max_level(unordered_map<int,Cell*> &);


//Fluxes header
void calculate_central_fluxes(unordered_map<int,Cell*> &);
void calculate_monotonic_fluxes(unordered_map<int,Cell*> &);
void calculate_local_fluxes(double (&q)[n_eq], double (&f)[n_eq]);
void weno3(int , double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq], double (&qright)[n_eq], double (&qleft)[n_eq]);
void weno5(int level, double (&qm2)[n_eq], double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq], double (&qp3)[n_eq], double (&qright)[n_eq], double (&qleft)[n_eq]);
void weno6(int level, double (&qm2)[n_eq], double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq], double (&qp3)[n_eq], double (&qright)[n_eq], double (&qleft)[n_eq]);
void muscl(double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq], double (&qright)[n_eq], double (&qleft)[n_eq]);
void conservative_correction(unordered_map<int,Cell*> &,unordered_map<int,Cell*> &);

//NLE Routines header
void stress_calc(double, double, double (&def_ten)[3][3], double (&stress)[3][3]);
void inverse_calc(double (&a)[3][3], double (&b)[3][3]);//b to be inverse of a
void transpose(double (&a)[3][3], double (&b)[3][3]);//b to be transpose of a
void matrix_mult(double (&a)[3][3], double (&b)[3][3], double (&c)[3][3]);//c = a.b
double invariant_1(double (&a)[3][3]);//Invariant 1 (trace)
double invariant_2(double (&a)[3][3]);//Invariant 2
double invariant_3(double (&a)[3][3]);//Invariant 3
void unit_matrix(double (&a)[3][3]);//Create 3X3 Unit Matrix
double nle_en_calc(double (&def_ten)[3][3], double );


//Rusanov Header
void rusanov_solver(Cell* &);
double swe_wave_speed(double (&q)[n_eq]);
double mhd_wave_speed(double (&q)[n_eq]);
double euler1d_wave_speed(double (&q)[n_eq]);
//double elastoplastic_wave_speed(double (&q)[n_eq]);
double nle_wave_speed(double (&q)[n_eq]);
double calc_wave_speed_weno3(double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq]);
double calc_wave_speed_weno5(double (&qm2)[n_eq], double (&qm1)[n_eq], double (&q)[n_eq], double (&qp1)[n_eq], double (&qp2)[n_eq], double (&qp3)[n_eq]);
double calc_max_wave_speed(unordered_map<int,Cell*> &);


//FORCE Header
void force_solver(Cell* &);

//Filtering header
double euler_pressure_calc(double (&q)[n_eq]);
void relaxation(unordered_map<int,Cell*> &);
void shock_relaxation(unordered_map<int,Cell*> &);



//Source_terms header
void calculate_sources(unordered_map<int,Cell*> &);

//Static Eval header
void analyse_static(unordered_map<int,Cell*> &);
void leaf_phys_vals_static(unordered_map<int,Cell*> &);

#endif // Function_Def_h_inluded
