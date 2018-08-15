using SequenceTransformations, Test, SpecialFunctions

a = Sequence(i->(-1)^(i-1)/(2i-1))

@test isa(a, Sequence)

s = cumsum(a)

@test isa(s, AbstractSequence)

@test norm(ϵ(s)[19:20].-π/4) < 20eps()

a = Sequence(i->(-1)^(i-1)/(2i-1)^2)
s = cumsum(a)

D,L,W = Drummond(s),Levin(s),Weniger(s)

@test isa(D, Transformation)

@test norm(D[40]-MathConstants.catalan) < 10eps()
@test norm(L[20]-MathConstants.catalan) < 10eps()
@test norm(W[20]-MathConstants.catalan) < 10eps()

const Gompertz = -.596347362323194074341078499369279376074
a = Sequence(i->(-1)^i*gamma(i))
s = cumsum(a)

@test norm(Levin(s)[18]-Gompertz) < 5e-12
@test Weniger(s)[19] ≈ Gompertz

@test Δ(s)[1] == a[2]

E = z-> lgamma(z) - z*log(z)+z-log(2(π/z))/2
z = range(10, stop=100, length=10)
@test norm(stirlingseries(z,Val{15}())-E.(z),Inf) < sqrt(eps())
