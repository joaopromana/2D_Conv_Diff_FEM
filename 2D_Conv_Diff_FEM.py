# FEM for 2D Convection-Diffusion Equation with Homogeneous Dirichlet BC's
# João Pedro Colaço Romana 2022

from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

# Geometric parameters
Lx = Ly = 10

# Time parameters
T = 1

# Parameters for source function
A = 5
alpha = 50

# Source points coordinates
x1 = y1 = 4
x2 = y2 = 4.5
x3 = 7
y3 = 6

# Parameters for w
B = 10
omega = 2 * pi / (4 * T)
D = 0.0025

Nx = int(input('Number of sub-intervals in x-direction: '))
Ny = int(input('Number of sub-intervals in y-direction: '))
num_steps = int(input('Number of time steps: '))
dt = T / num_steps

# Rectangle domain with mesh
mesh = RectangleMesh(Point(0, 0), Point(Lx, Ly), Nx, Ny)

# Finite element with Lagrange basis function
P1 = FiniteElement('Lagrange', triangle, 1)

# Function space V for scalar functions
V = FunctionSpace(mesh, P1)

# Function Space W for vector functions
W = FunctionSpace(mesh, MixedElement([P1, P1]))


# Dirichlet boundary
def boundary(x, on_boundary):
    return on_boundary


# C++ expression for the boundary condition
u_D = Constant(0.0)

# Applying essential Dirichlet BC
bc = DirichletBC(V, u_D, boundary)

# Define initial value
u_k = interpolate(u_D, V)

# Source function
f = Expression('A * (exp(-alpha * pow(x[0] - x1, 2) - alpha * pow(x[1] - y1, 2)) + exp(-alpha * pow(x[0] - x2, 2) - alpha * pow(x[1] - '
               'y2, 2)) + exp(-alpha * pow(x[0] - x3, 2) - alpha * pow(x[1] - y3, 2)))', degree=2, A=A, alpha=alpha, x1=x1, x2=x2, x3=x3, y1=y1, y2=y2, y3=y3)

# Velocity vector field
velocity = Expression(('B * cos(omega * t)', 'B * sin(omega * t)'), degree=2, B=B, omega=omega, t=0.0)
w = Function(W)
w.interpolate(velocity)

# Variational problem
u = TrialFunction(V)
v = TestFunction(V)
a = u * v * dx + dt * D * dot(grad(u), grad(v)) * dx + dt * dot(w, grad(u)) * v * dx
L = (u_k + dt * f) * v * dx

# Computing solution

u = Function(V)
t = 0

plt.show(block=True)

for n in range(num_steps):
    t = t + dt  # updating time
    velocity.t = t  # updating time inside velocity
    w.interpolate(velocity)  # updating w
    solve(a == L, u, bc)

    plt.ion()
    plt.figure(1)
    plt.clf()
    fig, ax = plt.subplots(nrows=1, ncols=1, num=1)
    points = mesh.coordinates()
    u_now = u_k.compute_vertex_values(mesh)
    im = ax.tripcolor(points[:, 0], points[:, 1], u_now)
    ax.set_aspect('equal', adjustable='box')
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal')
    ax.set_title(r'$u(t)$, t=' + str(np.round(t, 2)))

    u_k.assign(u)

    u_now = u_k.compute_vertex_values(mesh)
    im.remove()
    im = ax.tripcolor(points[:, 0], points[:, 1], u_now)
    cbar.update_normal(im)
    ax.set_title(r'$u(t=$' + str(np.round(t, 2)) + r'$)$ for $n=$' + str(num_steps) + r', $N_x=$' + str(Nx) + r' and $N_y=$' + str(Ny))
    
    plt.pause(0.5)
    if n == num_steps // 2:
        plt.savefig('2D_Conv_Diff_FEM_t05_n%02d.png' % num_steps, bbox_inches='tight', dpi=300)
    elif n == num_steps - 1:
        plt.savefig('2D_Conv_Diff_FEM_t1_n%02d.png' % num_steps, bbox_inches='tight', dpi=300)
    
