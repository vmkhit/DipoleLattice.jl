# DipoleLattice.jl

This package provides set of routines to calculate scattering properties of free standing and supported 2D periodic arrays of magnetoelectric dipolar particles. It is assumed that the scattering properties of these particles can be given by 6x6 dipole polarizability tensor ``\alpha``
```math
\alpha = \begin{bmatrix}
\alpha_{ee} & \alpha_{em}\\
\alpha_{me} & \alpha_{mm}
\end{bmatrix}
```
where, each of ``\alpha_{i, j}``, ``i, j = {e, m}`` are three by three matrices that describe the particle electric and magnetic polarizabilities as well as their cross-coupling, such that the for a given external fields ``\mathbf{E}^{\rm ext}`` and ``\mathbf{H}^{\rm ext}`` the induced electric and magnetic dipole moments are given by

```math
\begin{bmatrix}
\mathbf{p}\\
\mathbf{m}
\end{bmatrix} = \begin{bmatrix}
\alpha_{ee} & \alpha_{em}\\
\alpha_{me} & \alpha_{mm}
\end{bmatrix}
\begin{bmatrix}
\mathbf{E}^{\rm ext} \\
\mathbf{H}^{\rm ext}
\end{bmatrix}
```

The electric and magnetic fields at position ``\mathbf{r}`` produced by such an magnetoelectric scatterer located at some position ``\mathbf{r}_{0}`` can be given using the so called electric and magnetic dyadic green's functions of the electromagnetic field

```math
  \begin{bmatrix}
  \mathbf{E} \\
  \mathbf{H}^{\prime}
  \end{bmatrix} = \begin{bmatrix}
  \mathcal{G}^{\neq}_{E}(\mathbf{r}-\mathbf{r}_0) & \ii k_1 \mathcal{G}_{H}(\mathbf{r}-\mathbf{r}_0)\\
  -\ii k_1 \mathcal{G}^{\neq}_{H}(\mathbf{r}-\mathbf{r}_0) & \mathcal{G}_{E}(\mathbf{r}-\mathbf{r}_0)
  \end{bmatrix}\begin{bmatrix}
  \mathbf{p}^{\prime} \\
  \mathbf{m}^{\prime}
  \end{bmatrix}
```
where for the sake of symmetry the following normalized quantities are introduced ``\mathbf{p}^{\prime} = \mathbf{p}/\varepsilon``, ``\mathbf{m}^{\prime} = Z\mathbf{m}`` and ``\mathbf{H}^{\prime} = Z\mathbf{H}``, with ``\varepsilon`` and ``\mu`` being the dielectric function and permeability of the host medium and ``Z = \sqrt{\mu/\varepsilon}``. The dyadic Green functons are given as

```math
  \mathcal{G}_{E}^{\neq}(\rv-\rv') = \left[\mathds{I} k_1^2 +\nabla\otimes\nabla\right]\dfrac{{\rm e}^{{\rm i} k_1\vert\mathbf{r}-\mathbf{r}'\vert}}{\vert\mathbf{r}-\mathbf{r}'\vert}\\
  \mathcal{G}_{H}^{\neq}(\rv-\rv') = \dfrac{1}{k_{1}^2} \nabla\times \mathcal{G}^{\neq}_{E}(\mathbf{r}-\mathbf{r}') = \mathds{I}  \times \nabla\dfrac{\ee^{{\rm i} k_1\vert\mathbf{r}-\mathbf{r}'\vert}}{\vert\mathbf{r}-\mathbf{r}'\vert}
```
