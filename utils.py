import pandas as pd

def load_coeff(DREE_path, atomic_num_path, molar_mass_path, paas_path):
    """
    This function reads the input coeffcients for REEs from spreadsheets.
    """
    REE_list = ['La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd',
                   'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    # DREE (partitioning constants) [mol/mol, no unit]
    D_REE = pd.read_csv(DREE_path)  # (mol/mol)
    D_REE = D_REE[REE_list]
    # REE atomic number [amu]
    atomic_num = pd.read_csv(atomic_num_path)
    # REE molar mass [g/mol]
    molar_mass = pd.read_csv(molar_mass_path)
    molar_mass = molar_mass*10**-3  # mg to g
    # PAAS (ug REE /g ALL(ppm)) [mol/g] from Pourmand et al., 2012
    paas = pd.read_csv(paas_path)  # ug/g
    paas = paas*10**-6  # paas to g/g
    paas_mol = paas/molar_mass.loc[0]  # paas to mol
    paas_mol = paas_mol[REE_list]
    # Variables
    M_Ca = 0.01*10**-3  # mol/g
    X_Ca = 1
    n_Ca = M_Ca/X_Ca
    coeffs = [D_REE, atomic_num, molar_mass, paas_mol, n_Ca]
    return coeffs