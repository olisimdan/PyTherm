class Database(object):
    """Tools to access supporting database containing pure compound and mixture properties.

    All pure compound and mixture data could be stored in a SQLite database which is accessible using pugsql.

    Attributes:

    Methods:

    Errors:
    """


class Correlation(object):
    """Correlations for temperature dependent physical properties of pure compounds.

    Attributes:

    Methods:

    Errors:
    """


class Compound(object):
    """A pure chemical compound.

    Attributes:
        name: String representing the compound's name.
        cas_no: String representing the compound's Chemical Abstracts Service Registry Number.
        formula: String representing the compound's chemical formula.
        mw: Float representing the compound's molecular weight [=] (g/mol).
        t_t: Float representing the compound's triple point temperature [=] (K).
        t_c: Float representing the compound's critical point temperature [=] (K).
        p_c: Float representing the compound's critical point pressure [=] (Pa).
        rho_c: Float representing the compound's critical point density [=] (_).
        h_f_ig: Float representing the compound's Ideal Gas Enthalpy of Formation
        g_f_ig: Float representing the compound's Ideal Gas Gibbs Energy of Formation
        s_f_ig: Float representing the compound's Ideal Gas Entropy
        p_sat_data: List of tuples containing vapor pressure data (T, P_sat, P_sat_unc) [=] (K, Pa, Pa)
        rho_l_sat_data: List of tuples containing saturated liquid density data(T, rho, rho_unc) [=] (K, Pa, _)
        rho_v_sat_data: A list of tuples containing saturated vapor density data (T, rho, rho_unc) [=] (K, Pa, _)
        h_vap_data: A list of tuples containing enthalpy of vaporization data (T, h_vap, h_vap_unc) [=] (K, _, _)
        cp_data: A list of tuples containing heat capacity data (T, h_vap, h_vap_unc) [=] (K, _, _)
        cp_l_sat: A list of tuples containing saturated liquid heat capacity data (T, cp, cp_unc) [=] (K, _, _)
        cp_v_sat: A list of tuples containing saturated vapor heat capacity data (T, cp, cp_unc) [=] (K, _, _)


    Methods:
        cp_ig: Returns the ideal gas heat capacity at a given temperature [=] J/mol
        vis_v: Returns the vapor viscosity at a given temperature [=]
        cond_v: Returns the vapor thermal conductivity at a given temperature [=]
    """
    def __init__(self, name, cas_no, formula, mw):
        """Return a Compound object whose name is *name*."""
        self.name = name
        self.cas_no = cas_no
        self.formula = formula
        self.mw = mw


class Phase(object):
    """A mixture of of chemical compounds that form a homogeneous phase.

    A phase manages storage and evaluation of thermodynamic and transport properties.  The
        properties.

        Attributes:

            eos_method -- equation of state for thermodynamic property evaluation
            flash_method -- flash algorithm to use for equilibrum calculations
            transport_method -- approach used for transport property evaluation

            PVT:
            --------------------
            p -- pressure [=] Pa
            t -- temperature [=]K
            rho -- molar density [=] mol/m^3
            rho_mass -- mass density [=] kg/m^3

            Composition:
            --------------------
            x -- mole fraction as numpy array [=] mole fraction
            n -- total moles of phase [=] moles
            w -- mass composition as a numpy array [=] mass fraction
            m -- total mass of phase [=] kg

            Phase Equilibrium:
            --------------------
            f -- fugacity as a numpy array [=] Pa'
            ln_f -- ln(f) as a numpy array [=] Pa
            phi -- fugacity coefficient as a numpy array [=] dimensionless
            ln_phi -- ln(phi) coefficient as a numpy array [=] dimensionless

            Thermal Properties:
            --------------------
            h -- molar enthalpy [=] J/mol
            s -- molar entropy [=] J/(mol.K)
            h_ig -- ideal gas molar enthalpy [=] J/mol
            s_ig -- ideal gas molar entropy [=] J/(mol.K)
            h_r -- residual molar enthalpy [=] J/mol
            s_r -- residual molar entropy [=] J/(mol.K)
            cp -- isobaric heat capacity [=] J/(mol.K)
            cv -- isochoric heat capacity [=] J/(mol.K)
            cp_ig -- ideal gas isobaric heat capacity [=] J/(mol.K)
            cv_ig -- ideal gas isochoric heat capacity [=] J/(mol.K)
            cp_r -- residual isobaric heat capacity [=] J/(mol.K)
            cv_r -- residual isochoric heat capcity [=] J/(mol.K)

            Transport Properties:
            --------------------
            mu -- viscosity [=]
            mu_ig -- ideal gas viscosity [=]
            lambda -- thermal conductivity [=]
            lambda_ig -- ideal gas thermal conductivity [=]

        Methods:
            eos_rho
            eos_fugacity

        Errors:
            Component not available in database.
            Component list not the same size as composition list.

    Attributes:

    Methods:
    """
    def __init__(self, name, compounds):
        self.name = name
        self.compounds = compounds

    def density(self):
        pass

    def residual_enthalpy(self):
        pass

    def residual_entropy(self):
        pass

    def fugacity(self):
        pass


class Mixture(object):
    """A mixture of phases.

    Attributes:

    Methods:

    Errors:
    """
    def __init__(self, name, phases):
        self.name = name
        self.compounds = phases

    def phase_split(self):
        """Two phase solution of the Rachford-Rice equations."""
        pass

    def phase_assignment(self):
        """Assesses if a phase is a gas, liquid, or solid."""
        pass

class Cubic(object):
    """Generalized implementation of the Cubic Plus Association equation of state.

    Attributes:

    Methods:
    """

    def p(self):
        pass


class CubicPlusAssociation(object):
    """Generalized implementation of the Cubic Plus Association equation of state.

    Attributes:

    Methods:
    """

    def p(self):
        pass

class OtherEquationsOfState(object):
    """Generalized implementation of any other equiations of state like PC-SAFT, sPC-SAFT, etc...

    Attributes:

    Methods:
    """

    def p(self):
        pass

class Transport():
    """Expanded Fluid Viscosity (EVF) and Expanded Fluid Thermal Conductivity (EFTC) correlations.

    Attributes:

    Methods:

    """






