//
//  lib.hpp
//  MaxwellStefan
//
//  Created by Derek Harrison on 05/03/2022.
//

#ifndef lib_hpp
#define lib_hpp

#include <stdio.h>
#include <vector>

#include "user_types.h"

double ** create_D(int n);

void delete_D(double ** D, int n);

double dx_dz(int component,
             std::vector<double> mol_frac,
             p_params_t p_params,
             std::vector<double> J_vec);

std::vector<double> compute_composition(std::vector<double> mol_frac,
                                        p_params_t p_params,
                                        std::vector<double> J_vec,
                                        g_props_t g_props);

double error(std::vector<double> mol_frac,
             std::vector<double> mol_frac_E);

void compute_fluxes_rec(b_fracs_t b_fracs,
                        p_params_t p_params,
                        g_props_t g_props,
                        int flux_comp,
                        std::vector<double> J_vec_in,
                        f_bounds_t * J_vec_bounds,
                        double & min_dist,
                        double * J_vec);

double * compute_fluxes(b_fracs_t b_fracs,
                        p_params_t p_params,
                        g_props_t g_props);

std::vector<double> convert_to_vec(double * J_vec, int n);

mol_frac_res_t compute_fracs(p_params_t p_params,
                             g_props_t g_props,
                             b_props_t b_props,
                             t_params_t t_params,
                             b_fracs_t b_fracs);

void print_fractions(mol_frac_res_t mol_frac_res, t_params_t t_params, int n);


#endif /* lib_hpp */
