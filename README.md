# Bifurcation Diagram of a Periodically Forced Delay Differential Equation

This repository contains a Julia script for generating bifurcation diagrams of the periodically forced delay differential equation. 

The script solves the delay differential equation for a sequence of delay values and records the local maxima of the long-time solution to construct a classical bifurcation diagram.

## Model

The equation solved is a periodically forced version of the Mori-Zwangzig ENSO model 

\[
\dot{x}(t)
=
x(t)-x(t)^3
-\alpha x(t-\tau)\left(1-\gamma x(t-\tau)^2\right)
+c\cos(2\pi t),
\]

where

- $\alpha$ is the delayed feedback strength,
- $\tau$ is the delay,
- $\gamma$ controls the nonlinear delayed feedback,
- $c$ is the forcing amplitude.

## Method

For each value of the delay parameter $\tau$:

1. The DDE is integrated using `DifferentialEquations.jl`.
2. An initial transient is discarded.
3. Local maxima are detected directly from the adaptive solver output.
4. The maxima are plotted against $\tau$ to produce the bifurcation diagram.

Unlike many implementations, no uniform interpolation grid is used—the extrema are extracted directly from the adaptive time steps produced by the solver.

## Requirements

- Julia 1.10 or newer (should also work on recent Julia versions)

Packages:

```julia
using DifferentialEquations
using Plots
```

Install them with

```julia
using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("Plots")
```

## Running

Simply execute

```julia
julia bifurcation_diagram.jl
```

The script will

- compute the bifurcation diagram,
- display the figure,
- save it as

```
bifurcation_diagram.png
```

## Parameters

The current script uses

| Parameter | Value |
|-----------|------:|
| α | 0.75 |
| γ | 0.49 |
| c | 1.52 |
| τ | 4.10 → 3.91 |
| simulation time | 5000 |
| transient discarded | 1500 |

These values can easily be modified in the script.

## Notes

The local maxima are detected using neighbouring adaptive solver output points. This approach is computationally efficient and is suitable for producing high-resolution bifurcation diagrams.

## References

The Mori-Zwangzig ENSO model was first porposed in:

Falkena S. K., Quinn C., Sieber J., Frank J., and Dijkstra H. A. “Derivation of delay equation climate models using the Mori-Zwanzig formalism”. In: Proceedings of the Royal Society A 475.2227 (2019), p. 20190075.
