import numpy as np
from itertools import izip

debug = False

def intersect_half_plane(polygon, L, V):
    # print "check intersection ", L, V
    # print "polygon", polygon
    last_in_point = None
    first_in_point = None
    # print "number of edges ", len(polygon)
    new_polygon = []
    for e in polygon:
        # print "check edge ", e
        in_0 = sum(L*e[0]) < V
        in_1 = sum(L*e[1]) < V
        # print "in status ", in_0, in_1
        if in_0 and in_1:
            new_polygon.append(e)
        elif not in_0 and not in_1:
            continue
        else: # intersection
            # L . (e0 + x (e1-e0)) = V
            # L.e0 + x L.(e1-e0) = V
            # x = (V-L.e0)/(L.(e1-e0))
            x = (V - np.sum(L*e[0]) ) / np.sum(L*(e[1]-e[0]))
            intersection_p = e[0] + x*(e[1]-e[0])
            if in_0 and not in_1: # in to out
                new_polygon.append( [ e[0], intersection_p ] )
                last_in_point = intersection_p.copy()
            else: # no in_0 and in_1, out to in
                if last_in_point is None: # should only happen on first edge
                    first_in_point = intersection_p.copy()
                else:
                    new_polygon.append( [ last_in_point, intersection_p ] )
                new_polygon.append( [ intersection_p, e[1] ] )
    if first_in_point is not None:
        new_polygon.append( [ last_in_point, first_in_point ] )
    # print "new_polygon ", new_polygon
    return new_polygon

def mu_range(cur_min_EV, cur_composition, cur_bulk_struct, mcc_compositions, mcc_energies):
    if debug:
        print "mu_range got compositions ", mcc_compositions
        print "mu_range got energies ", mcc_energies
        print "got ref_bulk_model ", cur_min_EV, cur_composition

    n_types = len(cur_composition)
    Leq = []
    print "getting Veq from cur_min_EV, cur_composition"
    Veq = cur_min_EV * sum( [x[1] for x in cur_composition ] )
    print "mu_range: equality constraint:",
    cur_elements = []
    i_of_element = {}
    for (i, element) in enumerate(cur_composition):
        if i > 0:
            print "+",
        print "{}*mu_{}".format(element[1], element[0]),
        cur_elements.append(element[0])
        Leq.append(element[1])
        i_of_element[element[0]] = i
    print " = mu_{} = {}".format(cur_bulk_struct, Veq)

    if debug:
        print "constraint i_of_element, L, V ", i_of_element, Leq, "=", Veq

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
        Lne_cur = [0] * n_types
        ## if debug:
            ## print "   ", struct, mcc_energies, mcc_compositions
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

    if debug:
        for (L, V) in izip(Lne, Vne):
            print "inequality i_of_element, L, V ", i_of_element, L, "<=", V

    full_mu_range = None
    stable_mu_range = None
    if n_types == 2: # binary
        if debug:
            print "****************************************************************************************************"
            print "types ", [ (i_of_element[x], x) for x in sorted(i_of_element.keys()) ]
        Leq_3 = np.array([Leq[0], Leq[1], 0.0])
        stable_mu_range = [ [ -np.finfo(float).max, None ] , [ np.finfo(float).max, None ]  ]
        for (L, V) in izip(Lne, Vne):
            if np.all(L == Leq): # skip instance when inequality corresponds to same composition as base structure equilibrium equality
                continue
            if debug:
                print "inequality"
                print "Leq",Leq
                print "L from ne",L
                print "Veq",Veq
                print "V from ne",V
            A_lin_sys = np.array( [Leq, L] )
            b_lin_sys = np.array( [Veq, V] )
            x = np.linalg.solve(A_lin_sys, b_lin_sys)
            L_3 = np.array([L[0], L[1], 0.0])
            if np.cross(Leq_3, L_3)[2] < 0.0:
                # new max on component 0
                if x[0] < stable_mu_range[1][0]:
                    stable_mu_range[1][0] = x[0]
                    stable_mu_range[1][1] = x[1]
            else:
                # new min on component 0
                if x[0] > stable_mu_range[0][0]:
                    stable_mu_range[0][0] = x[0]
                    stable_mu_range[0][1] = x[1]
        full_mu_range = stable_mu_range

    elif n_types == 3: # ternary
        if debug:
            print "****************************************************************************************************"
            print "types ", [ (i_of_element[x], x) for x in sorted(i_of_element.keys()) ]
        # do monatomic limits first
        vertices = []
        for (L1, V1) in izip(Lne, Vne):
            if sum(L1 != 0) == 1: # monatomic
                i1 = np.where(L1 != 0)[0][0]
                for (L2, V2) in izip(Lne, Vne): # monatomic
                    i2 = np.where(L2 != 0)[0][0]
                    if sum(L2 != 0) == 1 and i1 < i2:
                        i3 = [0, 1, 2]
                        i3.remove(i1)
                        i3.remove(i2)
                        i3 = i3[0]
                        mu_3 = (Veq - Leq[i1]*V1 - Leq[i2]*V2)/Leq[i3]
                        pt = [ 0 ] * 3
                        pt[i1] = V1
                        pt[i2] = V2
                        pt[i3] = mu_3
                        vertices.append(np.array(pt))
        polygon = []
        for i in range (3):
            polygon.append([vertices[i], vertices[(i+1)%3]])
        full_mu_range = []
        for e in polygon:
            full_mu_range.append(e[0])
        if debug:
            print "initial polygon ", polygon
        for (L, V) in izip(Lne, Vne):
            if sum(L != 0) > 1:
                polygon = intersect_half_plane (polygon, L, V)
        if debug:
            print "****************************************************************************************************"
        stable_mu_range = []
        for e in polygon:
            stable_mu_range.append(e[0])
    elif n_types > 1: # n >= 4
        raise Exception("Can't do mu_range for %d-ary".format(n_types))

    if stable_mu_range is None:
        stable_mu_range_Z = None
    else:
        stable_mu_range_Z = []
        for mu_pt in stable_mu_range:
            mu_pt_Z = {}
            for Z in i_of_element:
                mu_pt_Z[Z] = mu_pt[i_of_element[Z]]
            stable_mu_range_Z.append(mu_pt_Z)

    if full_mu_range is None:
        full_mu_range_Z = None
    else:
        full_mu_range_Z = []
        for mu_pt in full_mu_range:
            mu_pt_Z = {}
            for Z in i_of_element:
                mu_pt_Z[Z] = mu_pt[i_of_element[Z]]
            full_mu_range_Z.append(mu_pt_Z)

    return (stable_mu_range_Z, full_mu_range_Z)
