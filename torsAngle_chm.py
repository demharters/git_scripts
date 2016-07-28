#! /usr/bin/env python
import MDAnalysis
from MDAnalysis import *
from MDAnalysis.analysis.distances import *
import numpy as np
import math
import sys

#def stp(vec1, vec2, vec3):
#    r"""Takes the scalar triple product of three vectors.
#
#    Returns the volume *V* of the parallel epiped spanned by the three
#    vectors
#
#    .. math::
#
#        V = \mathbf{v}_3 \cdot (\mathbf{v}_1 \times \mathbf{v}_2)
#
#    .. versionchanged:: 0.11.0
#       Moved into lib.mdamath
#    """
#    return np.dot(vec3, np.cross(vec1, vec2))
#
#def norm(v):
#    r"""Calculate the norm of a vector v.
#
#    .. math:: v = \sqrt{\mathbf{v}\cdot\mathbf{v}}
#
#    This version is faster then numpy.linalg.norm because it only works for a
#    single vector and therefore can skip a lot of the additional fuss
#    linalg.norm does.
#
#    Parameters
#    ----------
#    v: array_like
#        1D array of shape (N) for a vector of length N
#
#    Returns
#    -------
#    float
#        norm of the vector
#
#    """
#    return np.sqrt(np.dot(v, v))
#
#
#def normal(vec1, vec2):
#    r"""Returns the unit vector normal to two vectors.
#
#    .. math::
#
#       \hat{\mathbf{n}} = \frac{\mathbf{v}_1 \times \mathbf{v}_2}{|\mathbf{v}_1 \times \mathbf{v}_2|}
#
#    If the two vectors are collinear, the vector :math:`\mathbf{0}` is returned.
#
#    .. versionchanged:: 0.11.0
#       Moved into lib.mdamath
#    """
#    normal = np.cross(vec1, vec2)
#    n = norm(normal)
#    if n == 0.0:
#        return normal  # returns [0,0,0] instead of [nan,nan,nan]
#    return normal / n  # ... could also use numpy.nan_to_num(normal/norm(normal))
#
#
#def angle(a, b):
#    """Returns the angle between two vectors in radians
#
#    .. versionchanged:: 0.11.0
#       Moved into lib.mdamath
#    """
#    x = np.dot(a, b) / (norm(a) * norm(b))
#    # catch roundoffs that lead to nan otherwise
#    if x > 1.0:
#        return 0.0
#    elif x < -1.0:
#        return -np.pi
#    return np.arccos(x)
#
#def dihedral(ab, bc, cd):
#    """Returns the dihedral angle in radians between vectors connecting A,B,C,D.
#
#    The dihedral measures the rotation around bc::
#
#         ab
#       A---->B
#              \ bc
#              _\'
#                C---->D
#                  cd
#
#    The dihedral angle is restricted to the range -pi <= x <= pi.
#
#    .. versionadded:: 0.8
#    .. versionchanged:: 0.11.0
#       Moved into lib.mdamath
#    """
#    x = angle(normal(ab, bc), normal(bc, cd))
#    return (x if stp(ab, bc, cd) <= 0.0 else -x)


my_traj = sys.argv[1]

u = Universe(my_traj,my_traj)
v = Universe(my_traj)

end = my_traj.find('.pdb')
fout_name = my_traj[0:end] + '_torsAngle.dat'

#A1 = u.selectAtoms("segid A and resid 8 and name HO")
#A2 = u.selectAtoms("segid A and resid 8 and name OH")
#A3 = u.selectAtoms("segid A and resid 8 and name C5M")
#A4 = u.selectAtoms("segid A and resid 8 and name C5")
#A5 = u.selectAtoms("segid A and resid 8 and name C6")
#A6 = u.selectAtoms("segid A and resid 8 and name N1")

A1 = u.selectAtoms("segid A and resid 8 and (name HO or name OH or name C5M or name C5)")
A2 = u.selectAtoms("segid A and resid 8 and (name OH or name C5M or name C5 or name C6)")
A3 = u.selectAtoms("segid A and resid 8 and (name C5M or name C5 or name C6 or name N1)")

f = open(fout_name,'w')

for ts in u.trajectory:
        
    #v1 = A1.coordinates().flatten() - A2.coordinates().flatten() 
    #v2 = A2.coordinates().flatten() - A3.coordinates().flatten()
    #v3 = A3.coordinates().flatten() - A4.coordinates().flatten()
    #v4 = A4.coordinates().flatten() - A5.coordinates().flatten()
    #v5 = A5.coordinates().flatten() - A6.coordinates().flatten()

    #angle1 = dihedral(v1,v2,v3)
    #angle2 = dihedral(v2,v3,v4)
    #angle3 = dihedral(v3,v4,v5)

    angle1 = A1.dihedral()
    angle2 = A2.dihedral()
    angle3 = A3.dihedral()

    f.write('%7.3f\t%7.3f\t%7.3f\n' % (angle1,angle2,angle3))
f.close()

