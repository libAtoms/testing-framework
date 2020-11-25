import sys
import os
import ase.io
import numpy             as     np
import matplotlib.pyplot as     plt
from   ase.neighborlist  import neighbor_list

#=================================================================
# For persistance of np arrays, saving

def data_save(uid, array_2d, midpoints):
    np.savetxt("midpoints_{}.txt".format(uid), midpoints)
    np.savetxt("array-2d_{}.txt".format(uid), array_2d)

#=================================================================
# For persistance of np arrays, loading

def data_load(uid):
    midpoints = np.loadtxt("midpoints_{}.txt".format(uid))
    array_2d  = np.loadtxt("array-2d_{}.txt".format(uid))
    
    return array_2d, midpoints


#=================================================================
def getCoord(x, y, rho, cutoff, sum=True):
    """
    Returns either the coordination number (default) or the 1d normalised array
     of the integral.


    Arguments:
     x      - np.1darray - Bin midpoints from the RDF
     y      - np.1darray - RDF values for each bin
     rho    - float      - The number density for the system
     cutoff - float      - The value we are intergrating up to for coord calculation
     sum    - bool       - If True return the total integral, if False return the
                            1darray with the value for the area for each bin
    """
    for i in range(len(x)):
        if x[i] > cutoff:
            i -= 1
            break

    # Calculate the bin width
    dr = x[1] - x[0]
    
    # Shorten x & y to the cutoff
    x = x[:i]
    y = y[:i]
    
    # Intergrate
    integral = np.zeros(len(x))
    for i, r in enumerate(x):
        # int(o,r') g(r) * r^2 dr
        integral[i] = (y[i] * dr) * r**2
        
    # Normalise
    norm_integral = 4 * np.pi * rho * integral

    if sum:
        return np.sum(norm_integral)
    else:
        return norm_integral


#=================================================================
def get_rdf(atoms_list, cutoff_size, s1=None, s2=None, bins=200, save=False, comment="", verbose=False):
    """
    Returns a 2d array with each row the normalised histogram for that index
     in the atoms list, and a 1d array for the mid point of each bin (column
     in the 2d array).


    Arguments:
     atoms_list  - list - Containing ase Atoms objects
     cutoff_size - int  - Value for the cutoff of the RDF
     s1          - int  - Species 1 (None for all pairs)
     s2          - int  - Species 2 (None for all pairs)
     bins        - int  - Number of bins to use for the histogram
     save        - bool - If true the rdf data will be saved to text files, and attempted to be loaded
     comment     - str  - Used when saving files as part of the file name
     verbose     - bool - Set to true for increased verbosity
    """

    # Filename 
    if s1 == None:
        species_str    = "all"
    else:
        species_str    = "{}-{}".format(s1, s2)
    uid                = "species-{}_bins-{}_{}".format(species_str, bins, comment)

    # Load Data and exit if it exists 
    if save and os.path.isfile("midpoints_{}.txt".format(uid)) and os.path.isfile("array-2d_{}.txt".format(uid)):
        if verbose:
            print("Found data, loading from files: `midpoints_{uid:}.txt` & `array-2d_{uid:}.txt`".format(uid=uid))
        return data_load(uid)

    if type(atoms_list) is not list:
        raise ValueError("atoms_list must be a list object.")

    if type(s1) is int  and type(s2) is int:
        if verbose:
            print("Setting cutoff to species pair ({},{}), with range {}".format(s1, s2, cutoff_size))
        cutoff = {(s1, s2) : cutoff_size}
    else:
        if verbose:
            print("Setting cutoff to all pairs, with range {}".format(cutoff_size))
        cutoff = cutoff_size

    # 2d array for storing each frames normalised
    count_normd_2d = np.zeros((len(atoms_list), bins))

    # j = count, i = index in atoms file
    for j, atoms in enumerate(atoms_list): 
        list_d = neighbor_list("d", atoms, cutoff)
        
        # n_partices counted in the ase neighbor_list, which we must divide by later
        if type(cutoff) != dict:
            n_particles = len(atoms)
        else:
            n_particles = np.count_nonzero(atoms.get_atomic_numbers() == s1)
            if s1 == s2:
                n_particles = n_particles / 2

        #========================================
        count, bin_edges      = np.histogram(list_d, bins=bins, range=(0, cutoff_size))
       
        # Average number denisty of the system
        rho                 = len(atoms) / atoms.get_volume()
        #                     (( #P   / Volume of the spherical shell) / number density of the system) / number of particles counted
        count_normd_2d[j,:] = ((count / (4.0 * np.pi / 3.0 * (bin_edges[1:]**3 - bin_edges[:-1]**3))) / rho ) / n_particles

        if verbose:
            sys.stdout.write("Index: {:4.0f} | Number of pairs: {:>6} | {:>5.1f}%\r".format(j, len(list_d), ((j+1)/len(atoms_list))*100))
            sys.stdout.flush()
    if verbose:
        # Newline after sys.stdout.write() to move to the line beneith the 
        #   current line of output as we '\r ' rather than '\n' newlining
        print()

    # Calculate the midpoint of each bin
    bins_centre = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    if save:
        data_save(uid, count_normd_2d, bins_centre)

    return count_normd_2d, bins_centre

#=================================================================
def block_average(atoms_list, cutoff_size, n_blocks=10, s1=None, s2=None, bins=200, save=False, comment="", verbose=False):
    """
    Return an n (bins) by m (blocks) array of the average rdf for each block.

    Arguments:
     n_blocks - int - The number of blocks to use for block averaging

     see get_rdf for the further arguments
    """

    array_2d, midpoints = get_rdf(atoms_list, cutoff_size, s1, s2, bins, save, comment, verbose)

    block_avg_2d = np.zeros((n_blocks, bins))

    widths = np.arange(1, n_blocks+1, 1) * (len(atoms_list) / n_blocks)

    if verbose:
        print("Block widths are: {}".format(widths))

    for i, width in enumerate(widths):
        block_avg_2d[i,:] = np.mean(array_2d[:int(width), :], axis=0)

    return block_avg_2d, midpoints

#===========================================
def lammps2atoms(atoms_list, lammps_element={1 : "C", 2 : "Si"}):
    """
    Takes an atoms list, fixes it, and returns the fixed list.
    
    Args:
     atoms_list:     list of ase Atoms objects. (Will accept a single Atoms object)
     lammps_element: dictionary of lammps element number to chemical symbol lookup
    """
    
    new_atoms_list = []
    for j, atoms in enumerate(atoms_list):
        # Convert from lammps numbers to correct symbols
        for i, atom in enumerate(atoms):
            element     = atoms.arrays["element"][i]
            atom.symbol = lammps_element[element]
            atom.mass   = atoms.arrays["mass"][i]

        # Remove unneeded array values
        del atoms.arrays["mass"]
        del atoms.arrays["element"]
        
        # Fix info
        atoms.info["time"] = float(atoms.info["time"])
        
        # Append the new atoms object to the new list
        new_atoms_list.append(atoms)
    
    if len(new_atoms_list) == 1:
        new_atoms_list = new_atoms_list[0]
    
    return new_atoms_list
