__precompile__()
module SequenceTransformations
    using Base
    import Base: +,-,.+,.-,*,.*,/,./,muladd

    Base.muladd(x,y,z) = x*y+z

    include("Sequence.jl")
    include("Transformations.jl")
    include("generate.jl")
    include("FactorialSeries.jl")

end # module
