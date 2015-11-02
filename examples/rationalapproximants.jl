using SequenceTransformations, Polynomials # use my Polynomials fork Pkg.clone("https://github.com/MikaelSlevinsky/Polynomials.jl.git")

x = RatPoly(Poly([0;big(1)]),one(Poly{BigInt}))

# In this example, we use nonlinear sequence transformations to obtain rational approximants to power series.
# Consider the power series for the inverse tangent function.

a = Sequence(i->x*(-x^2)^(i-1)/(2i-1))

# The ϵ-algorithm corresponds to creating a diagonal sequence in the Padé table.

Epsilon = ϵ(cumsum(a))

Epsilon(1:20)

# Now we can see how the Padé aproximants compare to the function inverse tangent.

myatan(x,n) = convert(typeof(x),Epsilon(n)(x))

myatan(1.0,20)-atan(1.0) # Double precision. Summing 10 million terms can only deliver a little over 8 digits.

myatan(2.0,20)-atan(2.0) # Eight digits beyond region of convergence. Not bad!


# More powerful sequence transformations create useful approximants to factorially divergent power series.
# This series appears when computing the asymptotic expansion of the exponential integral.

a = Sequence(i->(-x)^(i)*factorial(big(i-1)))

D,L,W = Drummond(cumsum(a)),Levin(cumsum(a)),Weniger(cumsum(a))

W(1:20)

mydivergentseries(x,n) = convert(typeof(x),W(n)(x))

Gompertz = -.596347362323194074341078499369279376074

mydivergentseries(1.0,20)-Gompertz # It does quite well!

