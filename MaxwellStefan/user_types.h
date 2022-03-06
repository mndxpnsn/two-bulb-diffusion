//
//  user_types.h
//  MaxwellStefan
//
//  Created by mndx on 05/03/2022.
//

#ifndef user_types_h
#define user_types_h

const int NUM_TIME_STEPS = 20;
const double INF = 3e8;
const double RANGE = 1e3;
const double DEC_FAC = 10;
const int NUM_SCALE = 8;
const int NUM_GUESS = 20;

typedef struct grid_properties {
    double L;
    double dz;
} g_props_t;

typedef struct flux_bounds {
    double upper_bound;
    double lower_bound;
} f_bounds_t;

typedef struct mol_frac_results {
    std::vector<std::vector<double>> mol_frac1;
    std::vector<std::vector<double>> mol_frac2;
    std::vector<double> error;
} mol_frac_res_t;

typedef struct physical_properties_and_params {
    double ** D;
    double ct;
} p_params_t;

typedef struct time_params {
    double to;
    double tf;
} t_params_t;

typedef struct bulb_properties {
    double V;
    double d;
} b_props_t;

typedef struct bulb_mol_fractions {
    std::vector<double> mol_frac;
    std::vector<double> mol_frac_E;
} b_fracs_t;

#endif /* user_types_h */
