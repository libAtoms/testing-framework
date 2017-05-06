import numpy as np
from itertools import izip

def mu_range(cur_min_EV, cur_composition, cur_bulk_struct, mcc_compositions, mcc_energies):
    print "mu_range got compositions ", mcc_compositions
    print "mu_range got energies ", mcc_energies
    print "got ref_bulk_model ", cur_min_EV, cur_composition

    Leq = []
    Veq = cur_min_EV * sum( [x[1] for x in cur_composition ] )
    print "equality constraint:",
    cur_elements = []
    i_of_element = {}
    for (i, element) in enumerate(cur_composition):
        if i > 0:
            print "+",
        print "{}*mu_{}".format(element[1], element[0]),
        cur_elements.append(element[0])
        Leq.append(element[1])
        i_of_element[element[0]] = i
    print " = mu_{} = {}".format(cur_bulk_struct, cur_min_EV*len(cur_composition))

    ## print "constraint i_of_element, L, V ", i_of_element, Leq, "=", Veq

    Lne = []
    Vne = []
    print "inequalities:"
    for struct in mcc_energies:
        # skip the current structure, since it is an equality constraint, not a stability limit inequality
        if struct == cur_bulk_struct:
            continue
        # skip if any of the elements in this constraint aren't present in the material we're considering
        skip = False
        for element in mcc_compositions[struct]:
            if not element[0] in cur_elements:
                skip = True
                continue
        if skip:
            continue

        # determine inequality
        Lne_cur = [0] * len(cur_composition)
        # print "   ", struct, multicomponent_constraints_data[struct]["E_per_atom"][model_name], multicomponent_constraints_data[struct]["composition"]
        n_atoms = 0
        for (i, element) in enumerate(mcc_compositions[struct]):
            if i > 0:
                print "+",
            print "   {}*mu_{}".format(element[1], element[0]),
            n_atoms += element[1]
            Lne_cur[i_of_element[element[0]]] = element[1]
        print " <= mu_{} = {}".format(struct, mcc_energies[struct]*n_atoms)
        Lne.append(np.array(Lne_cur))
        Vne.append(mcc_energies[struct]*n_atoms)

    ## for (L, V) in izip(Lne, Vne):
        ## print "inequality i_of_element, L, V ", i_of_element, L, "<=", V

    mu_v_extrema = None
    if len(cur_composition) == 2: # binary
        Leq_3 = np.array([Leq[0], Leq[1], 0.0])
        mu_v_extrema = [ [ -np.finfo(float).max, None ] , [ np.finfo(float).max, None ]  ]
        for (L, V) in izip(Lne, Vne):
            A_lin_sys = np.array( [Leq, L] )
            b_lin_sys = np.array( [Veq, V] )
            x = np.linalg.solve(A_lin_sys, b_lin_sys)
            L_3 = np.array([L[0], L[1], 0.0])
            if np.cross(Leq_3, L_3)[2] < 0.0:
                # new max on component 0
                if x[0] < mu_v_extrema[1][0]:
                    mu_v_extrema[1][0] = x[0]
                    mu_v_extrema[1][1] = x[1]
            else:
                # new min on component 0
                if x[0] > mu_v_extrema[0][0]:
                    mu_v_extrema[0][0] = x[0]
                    mu_v_extrema[0][1] = x[1]

    elif len(cur_composition) == 3: # ternary
        print "****************************************************************************************************"
        print "types ", i_of_element
        print "equality constraint", Leq, Veq
        # do monatomic limits first do define triangle
        for (L1, V1) in izip(Lne, Vne):
            if sum(L1 != 0) == 1:
                i1 = np.where(L1 != 0)[0][0]
                for (L2, V2) in izip(Lne, Vne):
                    i2 = np.where(L2 != 0)[0][0]
                    if sum(L2 != 0) == 1 and i1 < i2:
                        i3 = [0, 1, 2]
                        i3.remove(i1)
                        i3.remove(i2)
                        i3 = i3[0]
                        mu_3 = (Veq - Leq[i1]*V1 - Leq[i2]*V2)/Leq[i3]
                        print "do monatomic inequality pair ", L1, V1, i1, L2, V2, i2, i3
                        print "got mus ", cur_elements[i1], V1, cur_elements[i2], V2, cur_elements[i3], mu_3, "sum",(Leq[i1]*V1+Leq[i2]*V2+Leq[i3]*mu_3)
                print ""
        print "****************************************************************************************************"
        raise Exception("Can't really do mu_range for terunary")
    elif len(cur_composition) > 1: # n >= 4
        raise Exception("Can't do mu_range for %d-ary".format(len(cur_composition)))

    if mu_v_extrema is None:
        mu_extrema = None
    else:
        mu_extrema = []
        for mu_pt in mu_v_extrema:
            mu_pt_Z = {}
            for Z in i_of_element:
                mu_pt_Z[Z] = mu_pt[i_of_element[Z]]
            mu_extrema.append(mu_pt_Z)

    return (mu_extrema)
