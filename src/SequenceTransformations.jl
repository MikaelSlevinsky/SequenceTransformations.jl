__precompile__()
module SequenceTransformations
    import Base: +, -, *, /

    include("Sequence.jl")
    include("Transformations.jl")
    include("generate.jl")
    include("FactorialSeries.jl")

end # module
