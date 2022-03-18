//
//  lib.cpp
//  MaxwellStefan
//
//  Created by mndx on 05/03/2022.
//

#include <iostream>
#include <stdio.h>
#include <vector>

#include "lib.hpp"
#include "user_types.h"

using namespace std;

double ** create_D(int n) {

    double ** D = new double * [n];
    
    for(int i = 0; i < n; ++i) {
        D[i] = new double[n];
        for(int j = 0; j < n; ++j)
            D[i][j] = 0.0;
    }
    
    return D;
}

void delete_D(double ** D, int n) {

    for(int i = 0; i < n; ++i)
        delete [] D[i];
    
    delete [] D;
}

double dx_dz(int component,
             vec_d_t mol_frac,
             p_params_t p_params,
             vec_d_t J_vec) {
    
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

vector<double> compute_composition(e_params_t e_params,
                                   vec_d_t J_vec) {


    vec_d_t mol_frac = e_params.b_fracs.mol_frac;
    p_params_t p_params = e_params.p_params;
    g_props_t g_props = e_params.g_props;
    
    int n = (int) mol_frac.size();
    vec_d_t mol_frac_E;
    vec_d_t mol_frac_in = mol_frac;
    
    int num_steps = g_props.L / g_props.dz;
    
    for(int c = 0; c < n; ++c) {
        for(int s = 0; s < num_steps; ++s)
            mol_frac[c] = -dx_dz(c, mol_frac_in, p_params, J_vec) * g_props.dz + mol_frac[c];
    }
    
    for(int c = 0; c < n; ++c)
        mol_frac_E.push_back(mol_frac[c]);
    
    return mol_frac_E;
}

double error(vec_d_t mol_frac,
             vec_d_t mol_frac_E) {
    
    double res = 0.0;
    
    int n = (int) mol_frac.size();
    
    for(int i = 0; i < n; ++i)
        res = res + (mol_frac[i] - mol_frac_E[i]) * (mol_frac[i] - mol_frac_E[i]);
    
    return res;
}

void compute_fluxes_rec(e_params_t e_params,
                        vec_d_t J_vec_in,
                        vec_fb_t & J_vec_bounds,
                        double & min_dist,
                        vec_d_t & J_vec) {

    b_fracs_t b_fracs = e_params.b_fracs;

    int n = (int) b_fracs.mol_frac.size();
    int m = (int) J_vec_in.size();
    double min_J = J_vec_bounds[m].lower_bound;
    double max_J = J_vec_bounds[m].upper_bound;
    
    // Number of guesses per component
    int ng = NUM_GUESS;
        
    double del_loc = (max_J - min_J) / ng;
    
    // Try all values for flux vector J
    if(m < n - 1) {
        for(int i = 0; i < ng; ++i) {
            double J_elem = min_J + i * del_loc;
            vec_d_t J_vec_loc = J_vec_in;

            J_vec_loc.push_back(J_elem);

            compute_fluxes_rec(e_params, J_vec_loc, J_vec_bounds, min_dist, J_vec);
        }
    }
    
    // Compute the last flux component of J
    if(m == n - 1) {
        vec_d_t J_vec_loc = J_vec_in;
        double j_elem_f = 0.0;

        for(auto j_elem : J_vec_in)
            j_elem_f = j_elem_f - j_elem;

        J_vec_loc.push_back(j_elem_f);

        compute_fluxes_rec(e_params, J_vec_loc, J_vec_bounds, min_dist, J_vec);
    }
    
    // Compute the minimum flux vector J
    if(m == n) {
        vec_d_t mol_frac_E_loc = compute_composition(e_params, J_vec_in);

        double err = error(b_fracs.mol_frac_E, mol_frac_E_loc);

        if(err < min_dist) {
            for(int k = 0; k < n; ++k)
                J_vec[k] = J_vec_in[k];

            min_dist = err;
        }
    }
}

std::vector<double> compute_fluxes(b_fracs_t b_fracs,
                                   p_params_t p_params,
                                   g_props_t g_props) {

    double min_dist = INF;
    int n = (int) b_fracs.mol_frac.size();
    
    e_params_t e_params;
    e_params.b_fracs = b_fracs;
    e_params.p_params = p_params;
    e_params.g_props = g_props;

    vec_d_t J_vec;
    vec_d_t J_vec_in;
    vec_fb_t J_vec_bounds;
    
    // Set the range, decrease factor and max number of iterations
    double range = RANGE;
    double dec_fac = DEC_FAC;
    int num_iterations = NUM_SCALE;
    
    // Initialize range bounds and flux vector
    for(int i = 0; i < n; ++i) {
        f_bounds_t f_bounds;
        f_bounds.upper_bound = range;
        f_bounds.lower_bound = -range;
        J_vec_bounds.push_back(f_bounds);
        J_vec.push_back(0);
    }
    
    // Compute fluxes
    int it = 0;

    while(it < num_iterations) {
        
        // Reset min_dist for calculation
        min_dist = INF;
        
        compute_fluxes_rec(e_params, J_vec_in, J_vec_bounds, min_dist, J_vec);
        
        range = range / dec_fac;

        for(int i = 0; i < n; ++i) {
            J_vec_bounds[i].upper_bound = J_vec[i] + range;
            J_vec_bounds[i].lower_bound = J_vec[i] - range;
        }
        
        ++it;
    }
    
    return J_vec;
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
        vec_d_t J_vec = compute_fluxes(b_fracs, p_params, g_props);
        
        for(int i = 0; i < n; ++i) {
            b_fracs.mol_frac[i] = b_fracs.mol_frac[i] - A * J_vec[i] * dt / (p_params.ct * b_props.V);
            b_fracs.mol_frac_E[i] = b_fracs.mol_frac_E[i] + A * J_vec[i] * dt / (p_params.ct * b_props.V);
        }
        
        mol_frac_results.mol_frac1.push_back(b_fracs.mol_frac);
        mol_frac_results.mol_frac2.push_back(b_fracs.mol_frac_E);
        
        t = t + dt;
    }
    
    return mol_frac_results;
}

void print_fractions(mol_frac_res_t mol_frac_res,
                     t_params_t t_params,
                     int n) {
    
    int nt = (int) mol_frac_res.mol_frac1.size();
    double dt = (t_params.tf - t_params.to) / nt;
    
    for(int i = 0; i < nt; ++i) {
        cout << "bulb1 composition at t " << (i + 1) * dt;

        for(int c = 0; c < n; ++c)
            cout << ", " << mol_frac_res.mol_frac1[i][c];

        cout << endl;

        cout << "bulb2 composition at t " << (i + 1) * dt;

        for(int c = 0; c < n; ++c)
            cout << ", " << mol_frac_res.mol_frac2[i][c];

        cout << endl;
    }
}

void store_fractions(mol_frac_res_t mol_frac_res, t_params_t t_params, int n) {
    int nt = (int) mol_frac_res.mol_frac1.size();
    double dt = (t_params.tf - t_params.to) / nt;

    FILE * file_ptr1 = fopen("results_bulb1.txt", "w");
    FILE * file_ptr2 = fopen("results_bulb2.txt", "w");

    if(file_ptr1 == nullptr)
        printf("file1 could not be opened\n");

    if(file_ptr2 == nullptr)
        printf("file2 could not be opened\n");

    double t = t_params.to;

    for(int i = 0; i < nt; ++i) {
        
        t = t + dt;
        
        fprintf(file_ptr1, "%f\t", t);
        fprintf(file_ptr2, "%f\t", t);
        
        for(int c = 0; c < n; ++c) {
            fprintf(file_ptr1, "%f\t", mol_frac_res.mol_frac1[i][c]);
            fprintf(file_ptr2, "%f\t", mol_frac_res.mol_frac2[i][c]);
        }
        
        fprintf(file_ptr1, "\n");
        fprintf(file_ptr2, "\n");
    }

    fclose(file_ptr1);
    fclose(file_ptr2);
}
