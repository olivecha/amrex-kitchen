
def recipe(field_indexes, box_array, sol_array):
    """
    omega_ratio_H2_O2
    """
    id_O2 = sol_array.species_index('O2')
    id_H2 = sol_array.species_index('H2')
    w_O2 = sol_array.net_production_rates[:, :, :, id_O2]
    w_H2 = sol_array.net_production_rates[:, :, :, id_H2]
    return w_H2/w_O2

