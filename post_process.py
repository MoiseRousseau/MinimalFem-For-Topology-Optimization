import meshio
import numpy as np

### CHANGE THIS BY YOUR <prefix> NAME
prefix = "test/cantilever"


mesh_file = prefix+".mesh"
src = open(mesh_file,'r')

n_vertices = int(src.readline())
vertices = np.zeros((n_vertices,3),dtype='f8')
for i in range(n_vertices):
  vertices[i,0:2] = [float(x) for x in src.readline().split()]

n_elements = int(src.readline())
elems = np.zeros((n_elements,3),dtype='f8')
for i in range(n_elements):
  elems[i] = [int(x) for x in src.readline().split()]

src.close()

data_file = prefix+".displacements" 
disp = np.genfromtxt(data_file) #point

mat_file = prefix+".matprops"
mat = np.genfromtxt(mat_file, skip_header=1) #cell, skip first line which is poisson ratio

stress_file = prefix+".stress"
stress = np.genfromtxt(stress_file) #cell

cells = [("triangle",elems)]
point_data = {"Displacements":disp}
cell_data = {"Young Modulus":[mat], "Von-Misses stress":[stress]}

mesh = meshio.Mesh(vertices, cells, point_data=point_data, cell_data=cell_data)
mesh.write(prefix+'.vtu')
