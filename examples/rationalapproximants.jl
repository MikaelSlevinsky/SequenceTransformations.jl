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


# The ratio of these sequences is used in Stirling's asymptotic approximation to the Gamma function.
A001163 = big([1,1,1,-139,-571,163879,5246819,-534703531,-4483131259,432261921612371,6232523202521089,-25834629665134204969,-1579029138854919086429,746590869962651602203151,1511513601028097903631961,-8849272268392873147705987190261,-142801712490607530608130701097701])
A001164 = big([1,12,288,51840,2488320,209018880,75246796800,902961561600,86684309913600,514904800886784000,86504006548979712000,13494625021640835072000,9716130015581401251840000,116593560186976815022080000,2798245444487443560529920000,299692087104605205332754432000000,57540880724084199423888850944000000])

# Here, we use the ϵ-algorithm to generate Padé approximants instead of a Poincaré-type asymptotic expansion.

x = RatPoly(Poly([0;big(1)]),one(Poly{BigInt}))

a = Sequence(i->A001163[i]//A001164[i]*x^(i-1))

ϵ(cumsum(a))(1:2:17)

# We can create the Levin d-transformations to the asymptotic expansion for different rational approximants.

Levin(cumsum(a),:d)(1:10)
