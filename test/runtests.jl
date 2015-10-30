using SequenceTransformations, Base.Test

a = Sequence(i->(-1)^(i-1)/(2i-1))

@test isa(a,Sequence)

s = cumsum(a)

@test isa(s,AbstractSequence)

@test norm(ϵ(s)(19:20)-π/4) < 20eps()

a = Sequence(i->(-1)^(i-1)/(2i-1)^2)
s = cumsum(a)

D,L,W = Drummond(s),Levin(s),Weniger(s)

@test isa(D,Transformation)

@test norm(D(40)-catalan) < 10eps()
@test norm(L(20)-catalan) < 10eps()
@test norm(W(20)-catalan) < 10eps()

Gompertz = -.596347362323194074341078499369279376074
a = Sequence(i->(-1)^i*gamma(i))
s = cumsum(a)

@test norm(Levin(s)(18)-Gompertz) < 5e-12
@test_approx_eq Weniger(s)(19) Gompertz
