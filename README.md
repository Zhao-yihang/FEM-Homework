
Hello! Here are several examples of operations based on the finite element method (actually part of a course assignment). Here, you can execute and experience the charm of the finite element method.

The Python package you must have is: FEniCS.


## Assignment 1: Solving with FEniCS
<img width="630" alt="image" src="https://github.com/Zhao-yihang/FEM-Homework/assets/157504045/b879addb-c8d9-4e89-9b21-6b6ac747cdab">

## Assignment 2: Solve Partial Differential Equation
Solve the equation:
\[ -\Delta u = f \]
where \( f(x, y) = 2\pi^2 \sin(\pi x) \sin(\pi y) \). The solution domain is a unit square, and the boundary conditions are homogeneous Dirichlet conditions \( u = 0 \). Experiment with different mesh divisions, and plot the relationship between the error and mesh size \( h \).

## Assignment 3: Cantilever Beam Analysis
### Problem Description
1. **Geometry and Boundary Conditions**:
   - Beam geometry: Length (\(L\)), height (\(H\)), width (\(W\)).
   - Fixed boundary condition: One end fixed.
   - External force: Vertical downward force applied at the free end.
2. **Physical Model**:
   - Elastic body equilibrium equation: \(\nabla \cdot \sigma + \mathbf{f} = \mathbf{0}\), where \(\sigma\) is the stress tensor, \(\mathbf{f}\) is the external force density.
   - Hooke's Law (Linear elastic material): \(\sigma = 2\mu\epsilon(\mathbf{u}) + \lambda \text{tr}(\epsilon(\mathbf{u}))\), where \(\mu\) is the shear modulus, \(\lambda\) is the bulk modulus, \(\epsilon(\mathbf{u})\) is the symmetric gradient of the displacement field (strain tensor), \(\mathbf{I}\) is the identity tensor.
### Solution Approach
1. **Geometry and Mesh Setup**:
   - Use FEniCS to create the geometry of the cantilever beam and generate an appropriate mesh.
2. **Define Function Space**:
   - Use FEniCS to define the function space for the displacement field.
3. **Set Boundary Conditions**:
   - Define fixed boundary conditions for the cantilever beam and apply the external force.
4. **Define Variational Problem**:
   - Construct the variational problem using the equilibrium equation and material constitutive equation.
5. **Solve Variational Problem**:
   - Use the `solve` function in FEniCS to solve the variational problem.
6. **Plot Results**:
   - Use Matplotlib to visualize the bending displacement field of the cantilever beam.
   - Adjust graphical properties such as color, resolution, and axis label positions for a more aesthetic representation.
By following these steps, we can establish a numerical model, solve the bending problem of the cantilever beam, and visualize the displacement field.
