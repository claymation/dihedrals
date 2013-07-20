"""
dihedrals.py -- Calculate the dihedral angles between peptide bond planes.

Author: Clay McClure <clay@daemons.net>

Usage:

    python dihedrals.py [--chain=CHAIN] <pdb_file>

Requirements:
    - numpy (for linear algebra)
    - BioPython (for PDB access)

Tested with:
    TIM barrel (1B9B)
    porin (2F1C)
    oxyhemoglobin (1HHO)
    deoxyhemoglobin (2HHB)
"""

import math
import numpy

from Bio.PDB import PDBParser


def calculate_chain_dihedrals(chain):
    """Return a list of (plane, plane, dihedral_angle) tuples for the residues in the given chain"""
    neighbor_residue = None
    neighbor_normal = None
    dihedrals = []

    for residue in chain:
        if not is_amino_acid(residue):
            continue

        if neighbor_residue is not None:
            centroid, normal = calculate_best_plane(neighbor_residue, residue)

            if neighbor_normal is not None:
                dihedral = calculate_dihedral(neighbor_normal, normal)
                dihedrals.append((neighbor_residue, dihedral))

            neighbor_normal = normal

        neighbor_residue = residue

    return dihedrals


def calculate_dihedral(vector1, vector2):
    """Return the dihedral angle between two vectors"""
    return math.degrees(math.acos(numpy.dot(vector1, vector2)))


def compute_covariance_matrix(coords, centroid):
    """Return the 3x3 variance-covariance matrix"""
    covariance = numpy.zeros((3, 3))
    cx, cy, cz = centroid
    for coord in coords:
        x, y, z = coord
        covariance[0][0] += (x-cx)**2
        covariance[0][1] += (x-cx)*(y-cy)
        covariance[0][2] += (x-cx)*(z-cz)
        covariance[1][0] += (y-cy)*(x-cx)
        covariance[1][1] += (y-cy)**2
        covariance[1][2] += (y-cy)*(z-cz)
        covariance[2][0] += (z-cz)*(x-cx)
        covariance[2][1] += (z-cz)*(y-cy)
        covariance[2][2] += (z-cz)**2

    covariance /= len(coords) - 1
    return covariance


def calculate_best_plane(residue1, residue2):
    """Return the centroid and normal vector of the best-fit plane for two neighboring residues"""
    # Best plane algorithm is based on:
    # http://blenderartists.org/forum/showthread.php?152374-Calculating-a-best-fit-plane-to-vertices-SOLVED#post1352722

    planar_atoms = get_planar_atoms(residue1, residue2)
    planar_coords = [atom.coord for atom in planar_atoms]

    centroid = calculate_centroid(planar_coords)
    covariance = compute_covariance_matrix(planar_coords, centroid)
    inverse = numpy.linalg.inv(covariance)

    # Iteratively find the best normal vector
    i = 0
    itermax = 500
    vector = numpy.array([1.0, 1.0, 1.0])
    product = numpy.dot(vector, inverse)
    normal = product / numpy.linalg.norm(product)

    while (vector != normal).any() and i < itermax:
        i += 1
        vector = normal
        product = numpy.dot(vector, inverse)
        normal = product / numpy.linalg.norm(product)

    return centroid, normal


def calculate_centroid(coords):
    """Return the centroid (center of mass) of the the given 3D coordinates"""
    centroid = numpy.zeros(3)
    for coord in coords:
        centroid += coord
    centroid /= len(coords)
    return centroid


def get_planar_atoms(residue, neighbor):
    """Return the (approximately) planar atoms comprising the peptide bond between a residue and its neighbor"""
    # NB: The neighbor N's H atom is also in the plane, but does not appear in crystal structure coordinates.
    return (residue['CA'], residue['N'], neighbor['C'], neighbor['O'], neighbor['CA'])


def is_amino_acid(residue):
    """Return True if the residue looks like an amino acid"""
    for atom in ('N', 'CA', 'C', 'O'):
        if not atom in residue.child_dict:
            return False
    return True


def _resname(residue):
    """Return the residue's 3-letter AA code and sequence number"""
    return '%s%s' % (residue.resname.capitalize(), residue.id[1])


def _print_csv(dihedrals):
    import csv
    import sys

    writer = csv.writer(sys.stdout)
    writer.writerow(['Residue', 'Dihedral angle'])

    for residue, dihedral in dihedrals:
        writer.writerow([_resname(residue), '%.1f' % dihedral])


def _print(dihedrals):
    print
    print 'Residue    Dihedral angle (degrees)'
    print '-------    ------------------------'

    for residue, dihedral in dihedrals:
        print '%-8s   %.1f' % (_resname(residue), dihedral)


if __name__ == '__main__':
    import sys

    from optparse import OptionParser

    parser = OptionParser(usage="usage: %prog [options] pdb_file")
    parser.add_option('-m', '--model', help="The PDB model to analyze (default: 0)", type="int", default=0)
    parser.add_option('-c', '--chain', help="The PDB chain to analyze (default: A)", default="A")
    parser.add_option('-C', '--csv', help="Write output in CSV (comma separated value) format", action="store_true")

    (options, args) = parser.parse_args()

    try:
        pdb_file = args[0]
    except IndexError:
        parser.error("pdb_file is a required argument")

    structure = PDBParser().get_structure('X', pdb_file)

    try:
        model = structure[options.model]
    except IndexError:
        parser.error("model %s not found in structure" % options.model)

    try:
        chain = model[options.chain]
    except KeyError:
        parser.error("chain %s not found in model %s in structure" % (options.chain, options.model))

    dihedrals = calculate_chain_dihedrals(chain)

    if options.csv:
        _print_csv(dihedrals)
    else:
        _print(dihedrals)
