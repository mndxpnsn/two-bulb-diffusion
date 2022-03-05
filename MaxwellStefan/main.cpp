//
//  main.cpp
//  MaxwellStefan
//
//  Created by mndx on 04/03/2022.
//  One dimensional Maxwell-Stefan diffusion.
//

#include <iostream>
#include <math.h>
#include <vector>

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
    int ng = 4e1;
        
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
    
    double min_dist = 3e8;
    int n = (int) b_fracs.mol_frac.size();
    
    double * J_vec = new double[n];
    std::vector<double> J_vec_in;
    
    f_bounds_t * J_vec_bounds = new f_bounds_t[n];
    
    // Set the range, decrease factor and max number of iterations
    double range = 1e3;
    double dec_fac = 10.0;
    int num_iterations = 6;
    
    for(int i = 0; i < n; ++i) {
        J_vec_bounds[i].upper_bound = range;
        J_vec_bounds[i].lower_bound = -range;
    }
    
    // Compute fluxes
    int it = 0;
    while(it < num_iterations) {
        
        // Reset min_dist for calculation
        min_dist = 3e8;
        
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
    int nt = 40;
    
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

int main(int argc, const char * argv[]) {
    
    // Number of components
    int num_components = 3;
    
    // Bulb1 mol fractions
    double x10 = 0.0; // H2
    double x20 = 0.501; // N2
    double x30 = 1.0 - x10 - x20; // CO2
    
    // Bulb2 mol fractions
    double x1E = 0.501; // H2
    double x2E = 0.499; // N2
    double x3E = 1.0 - x1E - x2E; // CO2
    
    // Diffusivities
    double D12 = 8.33e-5 * 3600; // units are (m2 / h)
    double D13 = 6.8e-5 * 3600; // units are (m2 / h)
    double D23 = 1.68e-5 * 3600; // units are (m2 / h)
    
    // Domain parameters
    g_props_t g_props;
    g_props.L = 1e-2; // units are (m)
    g_props.dz = 1e-2 / 10; // units are (m)
    
    // Bulb parameters
    b_props_t bulb_props;
    bulb_props.V = 5e-4; // Bulb volumes (m3)
    bulb_props.d = 2e-3; // Tube diameter (m)
    
    // Time parameters
    t_params_t time_params;
    time_params.to = 0.0; // Initial time (h)
    time_params.tf = 10.0; // Final time (h)
    
    // Total input concentration
    p_params_t p_params;
    p_params.ct = 1.0;
    
    // Organize data
    b_fracs_t bulb_mol_fracs;
    bulb_mol_fracs.mol_frac.push_back(x10);
    bulb_mol_fracs.mol_frac.push_back(x20);
    bulb_mol_fracs.mol_frac.push_back(x30);
    
    bulb_mol_fracs.mol_frac_E.push_back(x1E);
    bulb_mol_fracs.mol_frac_E.push_back(x2E);
    bulb_mol_fracs.mol_frac_E.push_back(x3E);
    
    p_params.D = create_D(num_components);
    
    p_params.D[0][1] = D12;
    p_params.D[0][2] = D13;
    p_params.D[1][2] = D23;
    
    p_params.D[1][0] = D12;
    p_params.D[2][0] = D13;
    p_params.D[2][1] = D23;
    
    // Perform simulation of two-bulb experiment
    mol_frac_res_t mol_frac_res = compute_fracs(p_params, g_props, bulb_props, time_params, bulb_mol_fracs);

    // Test computation of fluxes
    double * J_vec = compute_fluxes(bulb_mol_fracs, p_params, g_props);
    
    std::vector<double> J_vec_inp = convert_to_vec(J_vec, num_components);

    std::vector<double> exit_frac_E = compute_composition(bulb_mol_fracs.mol_frac, p_params, J_vec_inp, g_props);
    
    // Print results
    for(int i = 0; i < num_components; ++i)
        std::cout << "flux " << i << ": " << J_vec[i] << std::endl;

    for(int i = 0; i < num_components; ++i)
        std::cout << "bulb1 mole frac " << i << ": " << bulb_mol_fracs.mol_frac[i] << std::endl;

    for(int i = 0; i < num_components; ++i)
        std::cout << "bulb2 mole frac " << i << ": " << exit_frac_E[i] << std::endl;

    print_fractions(mol_frac_res, time_params, num_components);
    
    // Free allocated data
    delete_D(p_params.D, num_components);
    delete [] J_vec;
    
    return 0;
}
