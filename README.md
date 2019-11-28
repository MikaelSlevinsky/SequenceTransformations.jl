# SequenceTransformations.jl

[![Build Status](https://travis-ci.org/MikaelSlevinsky/SequenceTransformations.jl.svg?branch=master)](https://travis-ci.org/MikaelSlevinsky/SequenceTransformations.jl)

Sequence transformations are methods that manipulate sequences and series to accelerate their convergence.

Consider the convergent series for π/4 and Wynn's celebrated ϵ-algorithm:

```julia
using SequenceTransformations
a = Sequence(i->(-1)^(i-1)/(2i-1))
s = cumsum(a)
ϵ(s)(1:20)-π/4
```

Certain sequence transformations are powerful enough to sum divergent series. Consider the Euler series:

```julia
using SpecialFunctions
Gompertz = -.596347362323194074341078499369279376074
a = Sequence(i->(-1)^i*gamma(i))
s = cumsum(a)
Levin(s)(1:20)-Gompertz
Weniger(s)(1:20)-Gompertz
```

# References

E. J. Weniger, Nonlinear sequence transformations for the acceleration of convergence and the summation of divergent series, Comput. Phys. Rep., 10:189--371, 1989.
