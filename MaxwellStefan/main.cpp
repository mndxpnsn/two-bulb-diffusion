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
             double ** D,
             std::vector<double> J_vec,
             double ct) {
    
    int n = (int) mol_frac.size();
    
    double res = 0.0;

    for(int i = 0; i < n; ++i) {
        if(i != component) {
            res = res + (mol_frac[i] * J_vec[component] - mol_frac[component] * J_vec[i]) /
                  (ct * D[component][i]);
        }
    }
    
    return res;
}

std::vector<double> compute_composition(std::vector<double> mol_frac,
                                        double ** D,
                                        std::vector<double> J_vec,
                                        g_props_t g_props,
                                        double ct) {
    
    int n = (int) mol_frac.size();
    std::vector<double> mol_frac_E;
    
    int num_steps = g_props.L / g_props.dz;
    
    std::vector<double> mol_frac_in = mol_frac;
    
    for(int c = 0; c < n; ++c) {
        for(int s = 0; s < num_steps; ++s) {
            mol_frac[c] = -dx_dz(c, mol_frac_in, D, J_vec, ct) * g_props.dz + mol_frac[c];
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

void compute_fluxes_rec(std::vector<double> & mol_frac,
                        std::vector<double> & mol_frac_E,
                        double ** D,
                        g_props_t g_props,
                        double ct,
                        int flux_comp,
                        std::vector<double> J_vec_in,
                        f_bounds_t * J_vec_bounds,
                        double & min_dist,
                        double * J_vec) {
    
    int n = (int) mol_frac.size();
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
            compute_fluxes_rec(mol_frac, mol_frac_E, D, g_props, ct, flux_comp + 1,
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
        compute_fluxes_rec(mol_frac, mol_frac_E, D, g_props, ct, flux_comp + 1,
                           J_vec_loc, J_vec_bounds, min_dist, J_vec);
    }
    
    // Compute the minimum flux vector J
    if(m == n) {
        std::vector<double> mol_frac_E_loc = compute_composition(mol_frac, D, J_vec_in, g_props, ct);
        double err = error(mol_frac_E, mol_frac_E_loc);
        if(err < min_dist) {
            for(int k = 0; k < n; ++k) {
                J_vec[k] = J_vec_in[k];
            }
            min_dist = err;
        }
    }
}

double * compute_fluxes(std::vector<double> & mol_frac,
                        std::vector<double> & mol_frac_E,
                        double ** D,
                        g_props_t g_props,
                        double ct) {
    
    double min_dist = 3e8;
    int n = (int) mol_frac.size();
    
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
        
        compute_fluxes_rec(mol_frac, mol_frac_E, D, g_props, ct, 0,
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

mol_frac_res_t compute_fracs(double ** D,
                             g_props_t g_props,
                             double ct,
                             double V,
                             double d,
                             double to,
                             double tf,
                             std::vector<double> & mol_frac,
                             std::vector<double> & mol_frac_E) {
    
    // Number of time steps
    int nt = 40;
    
    double A = 3.14 * d * d / 4;
    int n = (int) mol_frac.size();
    double dt = (tf - to) / nt;
    double t = to;
    
    mol_frac_res_t mol_frac_results;
    
    while(t < tf) {
        double * J_vec = compute_fluxes(mol_frac, mol_frac_E, D, g_props, ct);
        
        for(int i = 0; i < n; ++i) {
            mol_frac[i] = mol_frac[i] - A * J_vec[i] * dt / (ct * V);
            mol_frac_E[i] = mol_frac_E[i] + A * J_vec[i] * dt / (ct * V);
        }
        
        mol_frac_results.mol_frac1.push_back(mol_frac);
        mol_frac_results.mol_frac2.push_back(mol_frac_E);
        
        t = t + dt;
        
        delete [] J_vec;
    }
    
    return mol_frac_results;
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
    double L = 1e-2; // units are (m)
    double dz = L / 10; // units are (m)
    
    // Bulb parameters
    double V = 5e-4; // Bulb volumes (m3)
    double d = 2e-3; // Tube diameter (m)
    
    // Time parameters
    double to = 0.0; // Initial time (h)
    double tf = 10.0; // Final time (h)
    
    // Total input concentration
    double ct = 1.0;
    
    // Organize data
    std::vector<double> mol_frac;
    mol_frac.push_back(x10);
    mol_frac.push_back(x20);
    mol_frac.push_back(x30);
    
    std::vector<double> mol_frac_E;
    mol_frac_E.push_back(x1E);
    mol_frac_E.push_back(x2E);
    mol_frac_E.push_back(x3E);
    
    double ** D = create_D(num_components);
    
    D[0][1] = D12;
    D[0][2] = D13;
    D[1][2] = D23;
    
    D[1][0] = D12;
    D[2][0] = D13;
    D[2][1] = D23;
    
    g_props_t g_props;
    g_props.L = L;
    g_props.dz = dz;
    
    // Compute fluxes
    double * J_vec = compute_fluxes(mol_frac, mol_frac_E, D, g_props, ct);
    
    // Verify computation of fluxes
    std::vector<double> J_vec_inp = convert_to_vec(J_vec, num_components);

    std::vector<double> exit_frac_E = compute_composition(mol_frac, D, J_vec_inp, g_props, ct);

    // Print results
    for(int i = 0; i < num_components; ++i)
        std::cout << "flux " << i << ": " << J_vec[i] << std::endl;

    for(int i = 0; i < num_components; ++i)
        std::cout << "bulb1 mole frac " << i << ": " << mol_frac[i] << std::endl;

    for(int i = 0; i < num_components; ++i)
        std::cout << "bulb2 mole frac " << i << ": " << exit_frac_E[i] << std::endl;

    // Perform simulation of two-bulb experiment
    mol_frac_res_t mol_frac_res = compute_fracs(D, g_props, ct, V, d, to, tf, mol_frac, mol_frac_E);

    int nt = (int) mol_frac_res.mol_frac1.size();
    double dt = (tf - to) / nt;
    
    for(int i = 0; i < nt; ++i) {
        std::cout << "bulb1 composition at t " << (i + 1) * dt;
        for(int c = 0; c < num_components; ++c) {
            std::cout << ", " << mol_frac_res.mol_frac1[i][c];
        }
        std::cout << std::endl;
        std::cout << "bulb2 composition at t " << (i + 1) * dt;
        for(int c = 0; c < num_components; ++c) {
            std::cout << ", " << mol_frac_res.mol_frac2[i][c];
        }
        std::cout << std::endl;
    }
    
    // Free allocated data
    delete_D(D, num_components);
    delete [] J_vec;
    
    return 0;
}
