import xThermoInterface as xt

import optimizer as opt

# Create a new instance of the thermodynamic property engine (Xiaodong's wrapper for xThermo.dll).
thermo = xt.xThermoInterface()

# Select the CPA equation of state.
thermo.ChooseAModel(1)

# Specify the number of components.
nc = 1
thermo.NoPureComp(nc)
 
# Define the critical properties for Methanol (Tc [=] K, Pc =[=] bar, omega [=] dimensionless).
tc = 512.5
pc = 80.84
omega = 0.565831

# Pass critical properties to the thermodynamic property engine.
thermo.CritProps(1, tc, pc, omega)

# Saturation properties for Methanol (saturation pressure [=] Pa, liquid density [=] kmol/m^3)
type = ['psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat',
        'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden',
        'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'psat', 'psat',
        'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat', 'psat',
        'psat', 'psat', 'psat', 'psat', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden',
        'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden', 'lden']


temp = [205.04,  213.583, 227.822, 242.061, 256.3, 270.539, 284.778, 299.017, 313.256, 327.494, 341.733, 355.972,
        370.211, 384.45, 398.689, 412.928, 427.167, 441.406, 455.644, 461.34, 205.04, 213.583, 227.822, 242.061, 256.3,
        270.539, 284.778, 299.017, 313.256, 327.494, 341.733, 355.972, 370.211, 384.45, 398.689, 412.928, 427.167,
        441.406, 455.644, 461.34, 205.04, 213.583, 227.822, 242.061, 256.3, 270.539, 284.778, 299.017, 313.256, 327.494,
        341.733, 355.972, 370.211, 384.45, 398.689, 412.928, 427.167, 441.406, 455.644, 461.34, 205.04, 213.583,
        227.822, 242.061, 256.3, 270.539, 284.778, 299.017, 313.256, 327.494, 341.733, 355.972, 370.211, 384.45,
        398.689, 412.928, 427.167, 441.406, 455.644, 461.34]

value = [11.22, 29.4, 123.4, 431, 1296, 3434, 8170, 17760, 35690, 67050, 118810, 200100, 322200, 499000, 746700,
         1084000, 1532000, 2116000, 2859000, 3206000, 27.3172710817053, 27.051994257537, 26.6181886274265,
         26.1875039011298, 25.7630609824605, 25.3448598714188, 24.9297796641907, 24.5125148242931, 24.0933774421072,
         23.6658136196242, 23.2288870857, 22.7701142250796, 22.2894950377629, 21.7807877161226, 21.2315086449036,
         20.6322951126646, 19.9737844079645, 19.2341302041071, 18.3914861743961, 18.013856812933, 8.86998014858439,
         24.3985086026671, 108.847413803305, 397.574909611848, 1232.10678345328, 3331.54911912398, 8036.33689533085,
         17605.1058147074, 35537.6249930851, 66889.457291887, 118566.422863296, 199548.040984799, 321078.968760752,
         496809.080175748, 742911.904650219, 1078204.32293734, 1524294.13136298, 2105781.63255288, 2850480.87931292,
         3201048.35320115, 27.1864220324593, 26.9700143796399, 26.6027404329431, 26.2266564112423, 25.8409979218469,
         25.4448868159131, 25.0373063405676, 24.6170689135158, 24.1827736714149, 23.7327817537779, 23.2650109590328,
         22.7770161719119, 22.2657074792961, 21.7271486677825, 21.1561989974983, 20.5459355912848, 19.8866677612897,
         19.1641284158074, 18.3558795402089, 18.0006643827451]

unc_frac = [0.051693, 0.0442176870748299, 0.0364667747163695, 0.0301624129930394, 0.029320987654321, 0.0195107746068725,
            0.0134638922888617, 0.00957207207207207, 0.00700476323900252, 0.00611483967188665, 0.00656510394747917,
            0.00749625187406297, 0.00869025450031037, 0.00961923847695391, 0.0104459622338288, 0.011070110701107,
            0.0117493472584856, 0.0127599243856333, 0.0143406785589367, 0.0149719276356831, 0.00434136867359762,
            0.00369173973234887, 0.00293117598780631, 0.00250268144440472, 0.00218049666868565, 0.00172392562492304,
            0.00125187781672509, 0.0010185503482169, 0.00104922279792746, 0.00131873928524331, 0.0017466075507188,
            0.00219298245614035, 0.00252030243629235, 0.00300902708124373, 0.00352785535793033, 0.00468915443957041,
            0.00640625, 0.0089242252149927, 0.0128966570507382, 0.0148995148995149, 0.397411236335855,
            0.307487527954031, 0.20054568622538, 0.126108651149987, 0.077785323567064, 0.0461275868429737,
            0.0249485131864265, 0.0144844343841984, 0.0105521964417422, 0.00919427402910911, 0.00986788647026127,
            0.0112754802748046, 0.0130808941370731, 0.0144924887392416, 0.0157488390302598, 0.0166944238833719,
            0.0177131168089306, 0.0192327634422858, 0.0215753069758616, 0.0224926311806561, 0.00721954413988909,
            0.0055544420907993, 0.00439931717976281, 0.00374841796801859, 0.00452402842309687, 0.00589668241902396,
            0.00644198750342243, 0.00637082889051615, 0.00554503572598608, 0.00423263493814689, 0.00261584335701284,
            0.00328847689989499, 0.00378449250038222, 0.00452468345645024, 0.00533954474153414, 0.00706329612951163,
            0.00965147036442655, 0.0134352347807032, 0.0193825108727434, 0.022365651812216]

weight = [1.0 / (unc_frac * value ** 2.0) for value, unc_frac in zip(value, unc_frac)]

# Define the association scheme.  Site order is glue, positive, negative (ex. 022 = 4C, 011 = 2B)
assoc_scheme = [11] * len(type)

# Define an args list that we can pass to the least squares optimizer.
cpa_sat_args = list(zip(assoc_scheme, type))


def cpa_sat(t, b, gamma, c1, beta, epsilon, assoc_scheme, prop_type):
    """Calculate the saturation pressure and liquid density.
    Args:
        t (float): Saturation temperature.
        b (float): Pure component CPA parameter.
        gamma (float): Pure component CPA parameter.
        c1 (float): Pure component CPA parameter.
        beta (float): Pure component CPA parameter.
        epsilon (float): Pure component CPA parameter.
        assoc_scheme (int): CPA association scheme.
        prop_type (string): Property to be calculated.
    Returns:
        prop_value: Returns saturation pressure (Pa) or liquid density (kmol/m^3).
    """
    # Pass physical and association parameters to the thermodynamic property engine.
    thermo.CPAParams(1, b, gamma, c1)
    thermo.AssocParams(1, assoc_scheme, beta, epsilon)
    # Fix model inputs passed to the thermodynamic property engine (required before property evaluation).
    thermo.Setup_Thermo()
    # Calculate the saturation pressure.
    psat, ln_k, ierr = thermo.PBubble(t, [1.0])
    # Calculate the liquid compressibility factor at saturated conditions.
    z_l, _ = thermo.FugacityCoeff(t, psat, [1.0], iph=1, job=0)
    # Universal gas constant, J/mol.K
    r = 8.3145
    # Calculate molar liquid density, mol/m^3
    rho_l_molar = (100000.0 * psat) / (z_l * r * t)
    if prop_type == 'psat':
        return psat * 100000.0  # Saturation pressure [=] Pa
    elif prop_type == 'lden':
        return rho_l_molar / 1000.0  # Saturation liquid density [=] kmol/m^3
    else:
        raise ValueError('prop_type not valid')


# Specify initial guess for CPA parameters.
b_0 = 30.978
gamma_0 = 1573.71
c1_0 = 0.431
beta_0 = 16.1
epsilon_0 = 2000.0

# Define bounds for CPA parameters based on the initial guess.
cpa_parms = [b_0, gamma_0, c1_0, beta_0, epsilon_0]
cpa_parms_bounds = [(0.9 * b_0, 1.1 * b_0),
                    (0.1 * gamma_0, 2.0 * gamma_0),
                    (0.1, 1.5),
                    (10, 200),
                    (0.5 * epsilon_0, 2.0 * epsilon_0)]

# Particle swarm optimization identifies a reasonable estimate of the global minimum. Optimization is subject to bounds
# which serve to (1) initialize a search space for the particle swarm algorithm and (2) keep the resulting search
# constrained to realistic CPA parameters.
ps_theta, nm_initial_size = opt.particle_swarm(opt.least_squares_objective_function,
                                               args=(cpa_sat, temp, value, weight, None, cpa_sat_args),
                                               bounds=cpa_parms_bounds,
                                               max_iter=50)

# The ps_theta estimated minimum is passed to the Nelder-Mead simplex algorithm to be refined. The Particle Swarm
# algorithm generates an initial size estimate for the simplex.  This helps ensure that the Nelder-Mead algorithm
# stays confined to the neighborhood around the minimum found by the Particle Swarm algorithm.
nm_theta = opt.nelder_mead(ps_theta,
                           opt.least_squares_objective_function,
                           args=(cpa_sat, temp, value, weight, None, cpa_sat_args),
                           small_tol=10.0**-10,
                           flat_tol=10.0**-9,
                           max_iter=2500,
                           initial_size=nm_initial_size)

# Print CPA parameter regression results.
print("CPA parameter regression (Nelder-Mead final):")
print(nm_theta)

# Define a list of temperatures and types to test regressed CPA parameters.
temp_type_value = list(zip(temp, type, value))

# Test regressed CPA parameters.
assoc_parms = nm_theta
assoc_scheme = 11
abs_dev_psat = []
abs_dev_lden = []
for i in temp_type_value:
    prop_temp = i[0]
    prop_type = i[1]
    prop_value = i[2]
    prop_calc_value = cpa_sat(prop_temp, *assoc_parms, assoc_scheme, prop_type)
    if prop_type == 'psat':
        abs_dev_psat.append(abs((prop_calc_value - prop_value)/prop_value))
    if prop_type == 'lden':
        abs_dev_lden.append(abs((prop_calc_value - prop_value)/prop_value))

print("Liquid Density % AAD:")
print(100.0*sum(abs_dev_lden)/len(abs_dev_lden))
print("Vapor Pressure % AAD:")
print(100.0*sum(abs_dev_psat)/len(abs_dev_psat))
