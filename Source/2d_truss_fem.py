import numpy as np

# ==========================================
# 1. DEFINE GEOMETRY & MATERIALS
# ==========================================
# Node coordinates: [X, Y] in meters
nodes = np.array([
    [0.0, 0.0],  # Node 0 (Bottom Wall Support)
    [0.0, 1.0],  # Node 1 (Top Wall Support)
    [1.0, 0.5]   # Node 2 (Free End)
])

# Element connectivity: [Node_i, Node_j]
elements = np.array([
    [0, 2],  # Element 0 connects Node 0 and 2
    [1, 2]   # Element 1 connects Node 1 and 2
])

A = 0.01       # Cross-sectional area (m^2)
E = 200e9      # Young's Modulus (N/m^2)

num_nodes = len(nodes)
num_dofs = 2 * num_nodes  # 2 DOFs per node (X and Y)

print(f"--- 2D Truss FEM Solver ({len(elements)} Elements, {num_dofs} DOFs) ---")

# ==========================================
# 2. INITIALIZE GLOBAL MATRICES
# ==========================================
K_global = np.zeros((num_dofs, num_dofs))
F_global = np.zeros((num_dofs, 1))

# Apply a downward load of 10,000 N at Node 2
# Node 2 DOFs are indices 4 (X) and 5 (Y). Downward is negative Y.
F_global[5, 0] = -10000.0

# ==========================================
# 3. LOCAL STIFFNESS & COORDINATE TRANSFORMATION
# ==========================================
for e in elements:
    n1, n2 = e
    x1, y1 = nodes[n1]
    x2, y2 = nodes[n2]
    
    # Calculate Length (L) and Angle components (cosine and sine)
    L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    c = (x2 - x1) / L  # cos(theta)
    s = (y2 - y1) / L  # sin(theta)
    
    # The 4x4 Transformed Element Stiffness Matrix
    # This matrix distributes the axial stiffness into global X and Y
    k_factor = (A * E) / L
    K_e = k_factor * np.array([
        [ c**2,  c*s, -c**2, -c*s],
        [ c*s,   s**2, -c*s, -s**2],
        [-c**2, -c*s,  c**2,  c*s],
        [-c*s,  -s**2,  c*s,  s**2]
    ])
    
    # Identify the global DOFs for these two nodes
    # Node 1 DOFs: [2*n1, 2*n1+1], Node 2 DOFs: [2*n2, 2*n2+1]
    dofs = [2*n1, 2*n1+1, 2*n2, 2*n2+1]
    
    # Add to global matrix using the ix_ slicer technique
    K_global[np.ix_(dofs, dofs)] += K_e

# ==========================================
# 4. APPLY BOUNDARY CONDITIONS (Partitioning)
# ==========================================
# We want to pin Nodes 0 and 1 to the wall. 
# Therefore DOFs 0, 1, 2, and 3 are fixed.
fixed_dofs = [0, 1, 2, 3]
free_dofs  = [4, 5]  # Only Node 2 can move

# Extract the Free blocks
K_ff = K_global[np.ix_(free_dofs, free_dofs)]
F_f  = F_global[free_dofs]

# ==========================================
# 5. SOLVE THE SYSTEM
# ==========================================
u_global = np.zeros((num_dofs, 1))

# Solve K_ff * d_f = F_f
u_free = np.linalg.solve(K_ff, F_f)
u_global[free_dofs] = u_free

# ==========================================
# 6. POST-PROCESSING (Results)
# ==========================================
print("\nNodal Displacements (meters):")
for i in range(num_nodes):
    print(f"Node {i}: u_x = {u_global[2*i, 0]:.8f}, u_y = {u_global[2*i+1, 0]:.8f}")

# Calculate Reactions: F_p = K_pf * u_f (since u_p = 0)
reactions = np.dot(K_global[np.ix_(fixed_dofs, free_dofs)], u_free)
print("\nReaction Forces (Newtons):")
print(f"Node 0: F_x = {reactions[0,0]:.2f}, F_y = {reactions[1,0]:.2f}")
print(f"Node 1: F_x = {reactions[2,0]:.2f}, F_y = {reactions[3,0]:.2f}")
