//
//  lib.hpp
//  MaxwellStefan
//
//  Created by mndx on 05/03/2022.
//

#ifndef lib_hpp
#define lib_hpp

#include <stdio.h>
#include <vector>

#include "user_types.h"

double ** create_D(int n);

void delete_D(double ** D, int n);

mol_frac_res_t compute_fracs(p_params_t p_params,
                             g_props_t g_props,
                             b_props_t b_props,
                             t_params_t t_params,
                             b_fracs_t b_fracs);

void print_fractions(mol_frac_res_t mol_frac_res, t_params_t t_params, int n);


#endif /* lib_hpp */
