# MinimalFem for topology optimization

Finite element modelling of elastic problem in two dimensions on triangular grids, and considering heterogeneous distribution on material properties with Jacobian and sensitivity output (for topology optimization)
Derived from the original work at [https://github.com/podgorskiy/MinimalFem](https://github.com/podgorskiy/MinimalFem) and [https://podgorskiy.com/spblog/304/writing-a-fem-solver-in-less-the-180-lines-of-code](https://podgorskiy.com/spblog/304/writing-a-fem-solver-in-less-the-180-lines-of-code).

## Getting started

1. Download or clone this repository
2. Install the Eigen library (`sudo apt install libeigen3-dev` on Ubuntu)
3. Run the `Makefile` (just type `make` in a terminal)
4. Install the Python library `meshio` and `numpy` to allow converting results in `vtu` format: `pip install meshio numpy`

## Use

`MinimalFem` solves the displacements (`<prefix>.displacements` output) and the Von-Misses stresses  (`<prefix>.stress`) of a two-dimensional part according to the static linear elasticity theory and using finite element analysis.
It requires a 2D triangular mesh (`<prefix>.mesh` input), the position of the load and the boundary conditions (`<prefix>.bcs`) and the heterogeneous material property distribution (Poisson ratio and Young modulus in `<prefix>.matprops`).
It also stores the derivative of the linear system solved relative to the displacements  (`<prefix>_jacobian.mtx` file) and to the Yound modulus (`<prefix>_sensitivity.mtx` file) allowing solving adjoint problem and carrying topology optimization.

The solver is called by the command: 
```
./MinimalFem <prefix>
```

After solving, the Python script `post_process.py` can be run (`python post_process.py`) to convert the raw results from the solver in `vtu` format for further visualization in [Paraview](https://www.paraview.org/).



### Inputs

Mesh file (`<prefix>.mesh`) must be organized as follow:
```
n_v   #number of vertices
X1 Y1 #coordinate first vertice
...
Xi Yi #coordinate vertice i
... 
Xn_v Yn_v #coordinate last vertice (n_v)
n_t #number of triangles
u1 v1 w1 #the three vertices of triangle 1
...
uj vj wj #the three vertices of triangle j
...
un_t vn_t wn_t #the three vertices of triangle n_t
```

Loads and boundary conditions (`<prefix>.bcs`) as:
```
n_constraint #number of constraint on the node mouvement
v1 type1 #node number / type = 1: no mouvement in x, type = 2: no mouvement in y, type = 3 no mouvement in x,y
...
vn_constraint typen_constraint
n_loads
v1 dx1 dy1 #vertice number / load in x / load in y
...
vn dxn_loads dyn_loads
```

Material properties (`<prefix>.matprops`) as:
```
poisson_ratio #constant
ym1 #young modulus for triangle 1 (defined by u1 v1 w1 in the mesh file)
...
ymi #young modulus for triangle i
...
ymn_t #young modulus for triangle n_triangles
```

For the two first, also see the original blog [post](https://podgorskiy.com/spblog/304/writing-a-fem-solver-in-less-the-180-lines-of-code) section `input data`.



### Outputs

Displacements are stored in the following format:
```
dx1 dy1 #displacement at first vertice
...
dxi dyi #displacement at vertice i (line i)
...
dxn_v dyn_v #displacement at last vertice
```

Von-Misses stresses as:
```
s1 #stress at triangle 1 (u1 v1 w1)
...
si #stress at triangle i (line i)
..
sn_t #stress at triangle n_t
```

Jacobian (derivative of linear system relative to displacement) and sensitivity (derivative relative to Young modulus) are stored using Eigen `write_binary_sparse` function. 
Both are sparse matrix of size `(2*n_v,n_t)`.
Example Python functions to read the saved matrix are available in file `read_jacobian_sensitivity.py`.

## Example

The `test` folder contains the input of a topology optimized cantilever fixed at the left and loaded at its bottom right. 
The case could be run from the root of this repository after compilation of the solver with:
```
./Minimal test/prefix
```

Post-processing is acheived when running the post-processing script `python post_process.py` by specifying the prefix name as (line 5): 
```
prefix = "test/cantilever"
```

After more processing with Paraview, the results is illustrated in the figure below:
![Post-processing with Paraview](https://github.com/MoiseRousseau/MinimalFem-For-Topology-Optimization/blob/master/test/cantilever.png)
