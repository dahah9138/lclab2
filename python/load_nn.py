# --------------------------------
# Author: Darian Hall            |
# Date: November 19, 2021        |
# --------------------------------
# Loads nn with format (Vx, Vy, Vz, 3)
# and translates into Nemaktis format to
# view the POM
# --------------------------------
import nemaktis as nm
import numpy as np
import scipy.io as sio
from numpy import pi, exp, cos, sin, sign, sqrt
import os.path

#mat_contents = sio.loadmat('low_res_type_II_Q2.mat')
mat_contents = sio.loadmat('../mat/cropped_dimer3.mat')

# Input format: Vx,Vy,Vz,3
nn = mat_contents['nn'];

# Reshape nn into Vz,Vy,Vx,3
nn = np.transpose(nn, (2,1,0,3))

m,n,p,d = nn.shape;
pitch = 4.5 # Micron

dcx,dcy,dcz = 2,2,3

cx = dcx * pitch # Micron
cy = dcy * pitch # Micron
cz = dcz * pitch # Micron


nfield = nm.DirectorField(
    mesh_lengths=(cx, cy, cz),
    mesh_dimensions=(p, n, m))
nfield.vals = nn

    
# We propagate fields through the lc layer
# Currently programmed for 5CB
mat = nm.LCMaterial(
    lc_field=nfield,
    ne=1.77,
    no=1.58)

wavelengths = np.linspace(0.4, 0.8, 11)
sim = nm.LightPropagator(
    material=mat, wavelengths=wavelengths, max_NA_objective=0.85)
output_fields = sim.propagate_fields(method="dtmm")

# Finally, the optical fields are visualized as in a real microscope
viewer = nm.FieldViewer(output_fields)
viewer.plot()