# ==============================================================================
# BLOCK 2: STATISTICAL MECHANICS ENGINE & QUASI-RRHO THERMOCHEMISTRY
# ==============================================================================
import inspect
import numpy as np

# Universal Constants (SI & eV)
R_SI = 8.314462618        # Universal gas constant [J / (mol · K)]
KB_SI = 1.380649e-23      # Boltzmann constant [J / K]
H_SI = 6.62607015e-34     # Planck constant [J · s]
HBAR_SI = H_SI / (2.0 * np.pi)  # Reduced Planck constant [J · s]
C_CMS = 2.99792458e10     # Speed of light [cm / s]
NA = 6.02214076e23        # Avogadro constant [mol^-1]
P_REF = 1.0e5             # Standard pressure P° = 1 bar [Pa]

EV_TO_JOULE = 1.602176634e-19
EV_TO_KJ_MOL = EV_TO_JOULE * NA / 1000.0  # ~ 96.4853 kJ/mol
AMU_TO_KG = 1.66053906660e-27
ANG_TO_M = 1.0e-10

def _get_species_database(species_db=None):
    """Retrieves species database from arguments, globals, calling frame, or default snapshot."""
    if species_db is not None:
        return species_db
    if 'species_database' in globals():
        return globals()['species_database']
    frame = inspect.currentframe().f_back
    while frame:
        if 'species_database' in frame.f_globals:
            return frame.f_globals['species_database']
        if 'species_database' in frame.f_locals:
            return frame.f_locals['species_database']
        frame = frame.f_back
    try:
        from kinetics.thermo.species_data import load_default_species_database
        return load_default_species_database()
    except Exception:
        raise NameError("species_database is not defined. Please define it in your notebook or pass species_db explicitly.")

def calculate_gas_thermo(species_name: str, T_K: float, nu_0: float = 100.0, species_db: dict = None, mode: str = 'qRRHO', **kwargs) -> dict:
    """
    Computes H°_gas(T), S°_gas(T), and G°_gas(T) in kJ/mol for a given species
    using canonical statistical mechanics with Grimme Quasi-RRHO damping (or benchmark DFT G).
    """
    db = _get_species_database(species_db)
    data = db[species_name]
    
    # Check if direct benchmark mode is requested or if vibrational data is not populated
    has_freqs = 'frequencies_cm1' in data and len(data['frequencies_cm1']) > 0
    if mode in ('wb97mv', 'b3lyp_benchmark') or not has_freqs:
        G_eV = data.get('G_wb97mv_eV', data.get('G_B3_eV', data.get('E_0K_eV', 0.0)))
        G_kJ_mol = G_eV * EV_TO_KJ_MOL
        return {
            'H_gas_kJ_mol': G_kJ_mol,
            'S_gas_J_mol_K': 0.0,
            'G_gas_kJ_mol': G_kJ_mol,
            'G_gas_eV': G_eV
        }
        
    mass_kg = (data['mass_g_mol'] / 1000.0) / NA
    E_0K_J_mol = data['E_0K_eV'] * EV_TO_KJ_MOL * 1000.0
    sigma_rot = data['sigma_rot']

    # ------------------------------------------------------------------------------
    # 1. TRANSLATIONAL PARTITION FUNCTION (3D PARTICLE IN BOX)
    # For each molecule, the 3 translational degrees of freedom contribute the
    # dominant entropy term through the thermal de Broglie wavelength and the
    # ideal-gas volume. This term reflects the number of accessible microstates
    # in the gas phase and is essential for the Sackur-Tetrode-like contribution.
    lambda_th = np.sqrt(H_SI**2 / (2.0 * np.pi * mass_kg * KB_SI * T_K))
    V_ideal = (KB_SI * T_K) / P_REF
    q_trans = V_ideal / (lambda_th**3)
    
    H_trans = 2.5 * R_SI * T_K
    S_trans = R_SI * (np.log(q_trans) + 2.5)


    # ------------------------------------------------------------------------------
    # 2. ROTATIONAL PARTITION FUNCTION (RIGID ROTOR)
    # A non-linear molecule has 3 rotational degrees of freedom, while a linear
    # molecule has only 2. These modes contribute rotational entropy and thermal
    # energy and depend on the principal moments of inertia, which encode the
    # molecular shape and symmetry.
    I_amu_A2 = data['moments_amu_A2']
    I_kg_m2 = [I * AMU_TO_KG * (ANG_TO_M**2) for I in I_amu_A2]
    
    if any(I == 0.0 for I in I_kg_m2):  # Linear molecule (e.g. CO2, H2, Ethyne)
        I_non_zero = max(I_kg_m2)
        q_rot = (8.0 * np.pi**2 * I_non_zero * KB_SI * T_K) / (sigma_rot * H_SI**2)
        H_rot = 1.0 * R_SI * T_K
        S_rot = R_SI * (np.log(q_rot) + 1.0)
    else:  # Non-linear polyatomic
        I_prod = I_kg_m2[0] * I_kg_m2[1] * I_kg_m2[2]
        q_rot = (np.sqrt(np.pi) / sigma_rot) * ((8.0 * np.pi**2 * KB_SI * T_K / H_SI**2)**1.5) * np.sqrt(I_prod)
        H_rot = 1.5 * R_SI * T_K
        S_rot = R_SI * (np.log(q_rot) + 1.5)


    # ------------------------------------------------------------------------------
    # 3. VIBRATIONAL PARTITION FUNCTION WITH QUASI-RRHO
    # For a non-linear molecule with N atoms, the total number of degrees of
    # freedom is 3N. After subtracting 3 translational and 3 rotational modes,
    # the remaining 3N - 6 vibrational modes describe internal motions. At low
    # frequencies, the harmonic approximation can overestimate entropy, so the
    # quasi-RRHO treatment damps these soft modes toward a free-rotor-like limit.
    freqs = np.array(data['frequencies_cm1'], dtype=float)
    c_cm_s = C_CMS
    
    # ZPE = sum(0.5 * h * c * nu)
    ZPE_J_mol = np.sum(0.5 * H_SI * c_cm_s * freqs) * NA
    
    # Thermal Vibrational Enthalpy (Harmonic Oscillator)
    theta_v = (H_SI * c_cm_s * freqs) / KB_SI
    x = theta_v / T_K
    exp_x = np.exp(np.clip(x, 0, 100))
    U_vib = R_SI * np.sum(theta_v / (exp_x - 1.0))
    H_vib = ZPE_J_mol + U_vib
    
    # Quasi-RRHO Entropy (Grimme 2012)
    S_vib_HO = R_SI * (x / (exp_x - 1.0) - np.log(1.0 - np.exp(-np.clip(x, 0, 100))))
    
    # Free Rotor Entropy for mode k:
    I_eff = HBAR_SI / (4.0 * np.pi * c_cm_s * freqs)
    mu_factor = 8.0 * np.pi**3 * I_eff * KB_SI * T_K / (H_SI**2)
    S_free_rot = R_SI * (0.5 + np.log(np.sqrt(mu_factor)))
    
    # Head-Gordon damping weights: w = 1 / (1 + (nu_0 / nu)^4)
    w = 1.0 / (1.0 + (nu_0 / freqs)**4)
    S_vib_qRRHO = np.sum(w * S_vib_HO + (1.0 - w) * S_free_rot)
    
    # 4. Total Standard Gas-Phase Properties (kJ / mol)
    H_tot_kJ_mol = (E_0K_J_mol + H_trans + H_rot + H_vib) / 1000.0
    S_tot_J_mol_K = S_trans + S_rot + S_vib_qRRHO
    G_tot_kJ_mol = H_tot_kJ_mol - (T_K * S_tot_J_mol_K / 1000.0)
    
    return {
        'H_gas_kJ_mol': H_tot_kJ_mol,
        'S_gas_J_mol_K': S_tot_J_mol_K,
        'G_gas_kJ_mol': G_tot_kJ_mol,
        'G_gas_eV': G_tot_kJ_mol / EV_TO_KJ_MOL
    }
