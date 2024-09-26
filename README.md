# 2D_Conv_Diff_FEM
Solution of 2D Convection-Diffusion equation with Finite Elements Method. The transport of a pollutant $`u(x,y,t)`$ in the atmosphere is described by the convection-diffusion equation. The following boundary-value problem is considered:

```math
\begin{equation}
    \begin{dcases} 
        \frac{\partial u}{\partial t} - D \Delta u + w \cdot \nabla u = f, \ \ \ \ \ (x, y, t) \in \Omega = (0, 10) \times (0, 10) \times I = (0, 1) \\
        u(x, y, t)=0, \ \ \ \ \ (x, y, t) \in \partial \Omega \\
        u(x, y, 0)=0, \ \ \ \ \ (x, y) \in \Omega \\
    \end{dcases}
\end{equation}
```

where $`D=2.5\times10^{-3}`$ is the effective diffusivity caused by the atmospheric turbulence and $`w(x,y,t) = (w_x, w_y)`$ is the vector field of the mean wind velocity with its components given by:

```math
\begin{equation}
  w(x,y,t) = 10 \left( \cos\left(\frac{\pi t}{2} \right), \sin\left(\frac{\pi t}{2} \right)  \right)
\end{equation}
```

The source function describing the emission rate density is given by the following:

```math
\begin{equation}
  f(x,y,t) = 5 \sum_{i=1}^3 e^{-50(x-x_i)^2-50(y-y_i)^2}
\end{equation}
```

where $`x_i=(x_1, x_2, x_3)=(4, 4.5, 7)`$ and $`y_i=(y_1, y_2, y_3)=(4, 4.5, 6)`$. The Backward-Euler method is used as the time integration scheme. Since the FEM is employed, the weak form has to be derived. This formulation is established by calculating the inner product of the strong primal residual with weighting functions $`v \in V_u`$, where the function space $`V_u(\Omega)`$ remains the same for both the primal and weighting functions. Integration by parts is subsequently applied to the resulting equation, leading to the following variational formulation for the time-discretized problem:

```math
\begin{equation}
    \int_\Omega \left( u^{n+1} v + h D \left( \nabla u^{n + 1} \cdot \nabla v \right) + h \left(\nabla w^{n + 1} \cdot \nabla u^{n + 1} \right) v \right) \ d\Omega = \int_{\Omega} \left(u^n v + h f v \right) \ d\Omega
\end{equation}
```

where $`u^n`$ is the solution evaluated at time step $`n`$ and $`h`$ is the time step interval. The legacy FEniCS package is used to solve the boundary-value problem with the FEM.
