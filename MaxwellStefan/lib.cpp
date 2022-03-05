//
//  lib.cpp
//  MaxwellStefan
//
//  Created by Derek Harrison on 05/03/2022.
//

#include <iostream>
#include <vector>

#include "lib.hpp"
#include "user_types.h"

double ** create_D(int n) {
    double ** D = new double * [n];
    
    for(int i = 0; i < n; ++i) {
        D[i] = new double[n];
        for(int j = 0; j < n; ++j) {
            D[i][j] = 0.0;
        }
    }
    
    return D;
}

void delete_D(double ** D, int n) {
    for(int i = 0; i < n; ++i) {
        delete [] D[i];
    }
    
    delete [] D;
}

double dx_dz(int component,
             std::vector<double> mol_frac,
             p_params_t p_params,
             std::vector<double> J_vec) {
    
    int n = (int) mol_frac.size();
    
    double res = 0.0;

    for(int i = 0; i < n; ++i) {
        if(i != component) {
            res = res + (mol_frac[i] * J_vec[component] - mol_frac[component] * J_vec[i]) /
                  (p_params.ct * p_params.D[component][i]);
        }
    }
    
    return res;
}

std::vector<double> compute_composition(std::vector<double> mol_frac,
                                        p_params_t p_params,
                                        std::vector<double> J_vec,
                                        g_props_t g_props) {
    
    int n = (int) mol_frac.size();
    std::vector<double> mol_frac_E;
    
    int num_steps = g_props.L / g_props.dz;
    
    std::vector<double> mol_frac_in = mol_frac;
    
    for(int c = 0; c < n; ++c) {
        for(int s = 0; s < num_steps; ++s) {
            mol_frac[c] = -dx_dz(c, mol_frac_in, p_params, J_vec) * g_props.dz + mol_frac[c];
        }
    }
    
    for(int c = 0; c < n; ++c) {
        mol_frac_E.push_back(mol_frac[c]);
    }
    
    return mol_frac_E;
}

double error(std::vector<double> mol_frac,
             std::vector<double> mol_frac_E) {
    
    double res = 0.0;
    
    int n = (int) mol_frac.size();
    
    for(int i = 0; i < n; ++i)
        res = res + (mol_frac[i] - mol_frac_E[i]) * (mol_frac[i] - mol_frac_E[i]);
    
    return res;
}

void compute_fluxes_rec(b_fracs_t b_fracs,
                        p_params_t p_params,
                        g_props_t g_props,
                        int flux_comp,
                        std::vector<double> J_vec_in,
                        f_bounds_t * J_vec_bounds,
                        double & min_dist,
                        double * J_vec) {
    
    int n = (int) b_fracs.mol_frac.size();
    int m = (int) J_vec_in.size();
    double min_J = J_vec_bounds[flux_comp].lower_bound;
    double max_J = J_vec_bounds[flux_comp].upper_bound;
    
    // Number of guesses per component
    int ng = NUM_GUESS;
        
    double del_loc = (max_J - min_J) / ng;
    
    // Try all values for flux vector J
    if(m < n - 1) {
        for(int i = 0; i < ng; ++i) {
            double J_elem = min_J + i * del_loc;
            std::vector<double> J_vec_loc = J_vec_in;
            J_vec_loc.push_back(J_elem);
            compute_fluxes_rec(b_fracs, p_params, g_props, flux_comp + 1,
                               J_vec_loc, J_vec_bounds, min_dist, J_vec);
        }
    }
    
    // Compute the last flux component of J
    if(m == n - 1) {
        std::vector<double> J_vec_loc = J_vec_in;
        double j_elem_f = 0.0;
        for(auto j_elem : J_vec_in) {
            j_elem_f = j_elem_f - j_elem;
        }
        J_vec_loc.push_back(j_elem_f);
        compute_fluxes_rec(b_fracs, p_params, g_props, flux_comp + 1,
                           J_vec_loc, J_vec_bounds, min_dist, J_vec);
    }
    
    // Compute the minimum flux vector J
    if(m == n) {
        std::vector<double> mol_frac_E_loc = compute_composition(b_fracs.mol_frac, p_params, J_vec_in, g_props);
        double err = error(b_fracs.mol_frac_E, mol_frac_E_loc);
        if(err < min_dist) {
            for(int k = 0; k < n; ++k) {
                J_vec[k] = J_vec_in[k];
            }
            min_dist = err;
        }
    }
}

double * compute_fluxes(b_fracs_t b_fracs,
                        p_params_t p_params,
                        g_props_t g_props) {
    
    double min_dist = INF;
    int n = (int) b_fracs.mol_frac.size();
    
    double * J_vec = new double[n];
    std::vector<double> J_vec_in;
    
    f_bounds_t * J_vec_bounds = new f_bounds_t[n];
    
    // Set the range, decrease factor and max number of iterations
    double range = RANGE;
    double dec_fac = DEC_FAC;
    int num_iterations = NUM_SCALE;
    
    for(int i = 0; i < n; ++i) {
        J_vec_bounds[i].upper_bound = range;
        J_vec_bounds[i].lower_bound = -range;
    }
    
    // Compute fluxes
    int it = 0;
    while(it < num_iterations) {
        
        // Reset min_dist for calculation
        min_dist = INF;
        
        compute_fluxes_rec(b_fracs, p_params, g_props, 0,
                           J_vec_in, J_vec_bounds, min_dist, J_vec);
        
        for(int i = 0; i < n; ++i) {
            J_vec_bounds[i].upper_bound = J_vec[i] + range;
            J_vec_bounds[i].lower_bound = J_vec[i] - range;
        }
        
        range = range / dec_fac;
        
        ++it;
    }
    
    return J_vec;
}

std::vector<double> convert_to_vec(double * J_vec, int n) {
    std::vector<double> J_vec_inp;
    for(int i = 0; i < n; ++i)
        J_vec_inp.push_back(J_vec[i]);
    
    return J_vec_inp;
}

mol_frac_res_t compute_fracs(p_params_t p_params,
                             g_props_t g_props,
                             b_props_t b_props,
                             t_params_t t_params,
                             b_fracs_t b_fracs) {
    
    // Number of time steps
    int nt = NUM_TIME_STEPS;
    
    double A = 3.14 * b_props.d * b_props.d / 4;
    int n = (int) b_fracs.mol_frac.size();
    double dt = (t_params.tf - t_params.to) / nt;
    double t = t_params.to;
    
    mol_frac_res_t mol_frac_results;
    
    while(t < t_params.tf) {
        double * J_vec = compute_fluxes(b_fracs, p_params, g_props);
        
        for(int i = 0; i < n; ++i) {
            b_fracs.mol_frac[i] = b_fracs.mol_frac[i] - A * J_vec[i] * dt / (p_params.ct * b_props.V);
            b_fracs.mol_frac_E[i] = b_fracs.mol_frac_E[i] + A * J_vec[i] * dt / (p_params.ct * b_props.V);
        }
        
        mol_frac_results.mol_frac1.push_back(b_fracs.mol_frac);
        mol_frac_results.mol_frac2.push_back(b_fracs.mol_frac_E);
        
        t = t + dt;
        
        delete [] J_vec;
    }
    
    return mol_frac_results;
}

void print_fractions(mol_frac_res_t mol_frac_res, t_params_t t_params, int n) {
    
    int nt = (int) mol_frac_res.mol_frac1.size();
    double dt = (t_params.tf - t_params.to) / nt;
    
    for(int i = 0; i < nt; ++i) {
        std::cout << "bulb1 composition at t " << (i + 1) * dt;
        for(int c = 0; c < n; ++c) {
            std::cout << ", " << mol_frac_res.mol_frac1[i][c];
        }
        std::cout << std::endl;
        std::cout << "bulb2 composition at t " << (i + 1) * dt;
        for(int c = 0; c < n; ++c) {
            std::cout << ", " << mol_frac_res.mol_frac2[i][c];
        }
        std::cout << std::endl;
    }
}
