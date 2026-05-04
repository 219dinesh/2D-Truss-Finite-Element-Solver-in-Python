# 2D Truss Finite Element Solver in Python

A lightweight, fully-commented Python script for performing 2D Finite Element Analysis (FEA) on pin-jointed truss structures.

This code bridges the gap between 1D bar elements and complex 2D structural analysis. It demonstrates how to apply geometric coordinate transformations to linear stiffness matrices, allowing you to solve for displacements and reaction forces in planar structures.

## Features

* **Coordinate Transformation:** Mathematically rotates 1D axial stiffness into 2D global $X$ and $Y$ coordinates.
* **Global Assembly:** Uses NumPy's `np.ix_` slicing technique to efficiently assemble the global stiffness matrix.
* **Matrix Partitioning:** Applies boundary conditions by isolating free vs. fixed degrees of freedom (DOFs).
* **Reaction Forces:** Calculates the structural support reactions based on the final deformed shape.
* **Zero Dependencies:** Runs entirely on standard linear algebra using `numpy`.

## The Mathematics of the 2D Truss Element

While a 1D bar only moves along its own axis, a 2D truss member can be oriented at any angle and can move in both the $X$ and $Y$ directions. We must map the local axial stiffness of the bar to the global coordinate system.

### 1. Local Stiffness and Coordinate Transformation
For a single truss element of length $L$, cross-sectional area $A$, and Young's Modulus $E$, the axial stiffness is $k = \frac{AE}{L}$. 

A 2D element has 2 nodes, and each node has 2 DOFs ($X$ and $Y$), resulting in a $4 \times 4$ local stiffness matrix. To transform the local axial stiffness into global coordinates, we use the angle of the element, $\theta$, where $c = \cos(\theta)$ and $s = \sin(\theta)$.

The geometric transformation matrix $\mathbf{T}$ relates global displacements to local axial stretching. Multiplying $\mathbf{T}^T \mathbf{k}_{local} \mathbf{T}$ yields the **Transformed Element Stiffness Matrix ($\mathbf{K}_e$)**:

$$ \mathbf{K}_e = \frac{AE}{L} \begin{bmatrix} c^2 & cs & -c^2 & -cs \\ 
cs & s^2 & -cs & -s^2 \\ 
-c^2 & -cs & c^2 & cs \\ 
-cs & -s^2 & cs & s^2 \end{bmatrix} $$

This matrix is calculated in the script using the differences in nodal coordinates to avoid calculating the actual angle:
&emsp *   $L = \sqrt{(x_2 - x_1)^2 + (y_2 - y_1)^2}$
&emsp *   $c = (x_2 - x_1) / L$
&emsp *   $s = (y_2 - y_1) / L$

### 2. Global Matrix Assembly
The global stiffness matrix $\mathbf{K}$ defines the entire structure. If the structure has $N$ nodes, the global matrix is of size $2N \times 2N$.

For each element, we find its global DOFs (e.g., Node 1 uses DOFs `[2, 3]`, Node 2 uses DOFs `[4, 5]`) and add the $4 \times 4$ element matrix $\mathbf{K}_e$ into the massive global matrix $\mathbf{K}$ at those exact indices.

### 3. Boundary Conditions & Matrix Partitioning
Before solving $\mathbf{K} \mathbf{U} = \mathbf{F}$, we must prevent the structure from flying away by applying boundary conditions (supports). We partition the matrices into **Free** ($f$) DOFs and **Fixed** ($p$) DOFs:

$$ \begin{bmatrix} \mathbf{K}_{ff} & \mathbf{K}_{fp} \\ 
\mathbf{K}_{pf} & \mathbf{K}_{pp} \end{bmatrix} \begin{bmatrix} \mathbf{U}_f \\ 
\mathbf{U}_p \end{bmatrix} = \begin{bmatrix} \mathbf{F}_f \\ 
\mathbf{F}_p \end{bmatrix} $$

Since our wall supports don't move, $\mathbf{U}_p = \mathbf{0}$. We extract the sub-matrix $\mathbf{K}_{ff}$ and the known force vector $\mathbf{F}_f$ to solve for the unknown displacements.

### 4. Solving Displacements and Reactions
We solve the reduced linear system for the displacements at the free nodes:

$$ \mathbf{U}_f = \mathbf{K}_{ff}^{-1} \mathbf{F}_f $$

Finally, we calculate the reaction forces pushing back at the wall supports. Looking at the bottom row of our partitioned matrix, and knowing $\mathbf{U}_p = \mathbf{0}$:

$$ \mathbf{F}_p = \mathbf{K}_{pf} \mathbf{U}_f $$

## Prerequisites
To run this script, you only need Python and the NumPy library installed.
```bash
pip install numpy
```
## Usage & Example Output

Run the script directly from your terminal:
```bash
python 2d_truss_fem.py
```
## Example Output
Based on the default configuration (a 3-node triangular bracket, 10kN downward point load at the tip), you will see the following output:

```Plaintext
--- 2D Truss FEM Solver (2 Elements, 6 DOFs) ---

Nodal Displacements (meters):
Node 0: u_x = 0.00000000, u_y = 0.00000000
Node 1: u_x = 0.00000000, u_y = 0.00000000
Node 2: u_x = 0.00000000, u_y = -0.00001398

Reaction Forces (Newtons):
Node 0: F_x = 10000.00, F_y = 5000.00
Node 1: F_x = -10000.00, F_y = 5000.00
```
