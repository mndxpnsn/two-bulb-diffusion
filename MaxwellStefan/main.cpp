//
//  main.cpp
//  MaxwellStefan
//
//  Created by mndx on 04/03/2022.
//  One dimensional Maxwell-Stefan diffusion.
//

#include <iostream>
#include <vector>

#include "lib.hpp"
#include "user_types.h"

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
    
    // Total concentration
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
    mol_frac_res_t mol_frac_res = compute_fracs(p_params,
                                                g_props,
                                                bulb_props,
                                                time_params,
                                                bulb_mol_fracs);
    
    // Print results
    print_fractions(mol_frac_res,
                    time_params,
                    num_components);
    
    // Free allocated data
    delete_D(p_params.D, num_components);
    
    return 0;
}
