"""
    AutoMinLevel()

Sentinel requesting that a tree constructor choose its default minimum level.
"""
struct AutoMinLevel end

"""
    NoProtrusionCheck()

Disable protrusion-based subdivision checks.
"""
struct NoProtrusionCheck end

"""
    AutoProtrusionCheck()

Sentinel requesting the default protrusion policy for the input type.

Plain point collections resolve to [`NoProtrusionCheck`](@ref). Extensions can
provide input-specific defaults, such as BEAST basis-function protrusions.
"""
struct AutoProtrusionCheck end

"""
    ProtrusionCheck(; max, compute=ComputeProtrusionFunctor())

Subdivision policy that rejects a split when any value protrudes beyond `max`
relative to a child box.

`compute` is called as `compute(center, halfsize, value)`.
"""
struct ProtrusionCheck{M,C}
    max::M
    compute::C
end

function ProtrusionCheck(; max, compute=ComputeProtrusionFunctor())
    return ProtrusionCheck(max, compute)
end

"""
    TwoNTreeBuilder(; minhalfsize=0, minlevel=AutoMinLevel(), root=1, minvalues=70,
        protrusion=AutoProtrusionCheck())

Construction settings for [`TwoNTree`](@ref).

The defaults build an adaptive tree driven mainly by `minvalues`: subdivision
stops once a node holds at most `minvalues` values, unless `minhalfsize`,
`minlevel`, or the protrusion policy impose a stronger constraint.
`AutoProtrusionCheck()` lets the input type choose its default protrusion
policy.
"""
struct TwoNTreeBuilder{M,L,P}
    minhalfsize::M
    minlevel::L
    root::Int
    minvalues::Int
    protrusion::P
end

function TwoNTreeBuilder(;
    minhalfsize=0,
    minlevel=AutoMinLevel(),
    root::Int=1,
    minvalues::Int=70,
    protrusion=AutoProtrusionCheck(),
)
    return TwoNTreeBuilder(minhalfsize, minlevel, root, minvalues, protrusion)
end

"""
    BlockTreeBuilder(; test=TwoNTreeBuilder(), trial=TwoNTreeBuilder())

Construction settings for a [`BlockTree`](@ref).

`test` and `trial` are [`TwoNTreeBuilder`](@ref)s. Their `minhalfsize` and
`root` values must match so both sides use compatible level indexing.
"""
struct BlockTreeBuilder{T,R}
    test::T
    trial::R
end

function BlockTreeBuilder(; test=TwoNTreeBuilder(), trial=TwoNTreeBuilder())
    _validate_block_tree_builders(test, trial)
    return BlockTreeBuilder(test, trial)
end

"""
    defaultprotrusionfunctor(input)

Return the default protrusion evaluator for `input`.

The fallback is [`ComputeProtrusionFunctor`](@ref), which evaluates to zero.
"""
defaultprotrusionfunctor(input) = ComputeProtrusionFunctor()

"""
    defaultprotrusioncheck(input)

Return the default protrusion policy for `input`.

The fallback is [`NoProtrusionCheck`](@ref). Extensions can specialize this for
richer input types.
"""
defaultprotrusioncheck(input) = NoProtrusionCheck()

_resolve_protrusion(protrusion, input) = protrusion
function _resolve_protrusion(::AutoProtrusionCheck, input)
    return defaultprotrusioncheck(input)
end
function _resolve_protrusion(protrusion::ProtrusionCheck, input)
    return _resolve_protrusioncheck(protrusion, protrusion.compute, input)
end
_resolve_protrusioncheck(protrusion, compute, input) = protrusion
function _resolve_protrusioncheck(
    protrusion::ProtrusionCheck, ::ComputeProtrusionFunctor, input
)
    compute = defaultprotrusionfunctor(input)
    compute isa ComputeProtrusionFunctor && return protrusion
    return ProtrusionCheck(; max=protrusion.max, compute=compute)
end

_resolve_builder_protrusion(builder, input) = builder
function _resolve_builder_protrusion(builder::TwoNTreeBuilder, input)
    protrusion = _resolve_protrusion(builder.protrusion, input)
    protrusion === builder.protrusion && return builder
    return TwoNTreeBuilder(;
        minhalfsize=builder.minhalfsize,
        minlevel=builder.minlevel,
        root=builder.root,
        minvalues=builder.minvalues,
        protrusion=protrusion,
    )
end

function _validate_block_tree_builders(test::TwoNTreeBuilder, trial::TwoNTreeBuilder)
    test.minhalfsize == trial.minhalfsize || throw(
        ArgumentError("BlockTreeBuilder requires test.minhalfsize == trial.minhalfsize")
    )
    test.root == trial.root ||
        throw(ArgumentError("BlockTreeBuilder requires test.root == trial.root"))
    return nothing
end

_resolve_minlevel(::AutoMinLevel, default::Int) = default
_resolve_minlevel(minlevel::Int, default::Int) = minlevel
function _resolve_minlevel(minlevel, default::Int)
    return throw(
        ArgumentError("minlevel must be an Int or AutoMinLevel(), got $(typeof(minlevel))")
    )
end

"""
    BoundingBallTreeBuilder(; splitter, numsplits, minvalues=numsplits,
        minlevel=AutoMinLevel(), maxlevel=typemax(Int), root=1,
        balldata=BoundingBallData, updateradii=boundingsphere,
        splitterkwargs=NamedTuple())

Construction settings for [`BoundingBallTree`](@ref). `splitter` is called via
`splitwrapper(points, values, level, numsplits, splitter; splitterkwargs...)`.
For compatibility with existing custom ball splitters, the wrapper also accepts splitters
that ignore `level` and implement `(points, values, numsplits)`. `updateradii` computes the
bounding ball for each node.
"""
struct BoundingBallTreeBuilder{S,L,B,U,K}
    splitter::S
    numsplits::Int
    minvalues::Int
    minlevel::L
    maxlevel::Int
    root::Int
    balldata::B
    updateradii::U
    splitterkwargs::K
end

function BoundingBallTreeBuilder(;
    splitter,
    numsplits::Int,
    minvalues::Int=numsplits,
    minlevel=AutoMinLevel(),
    maxlevel::Int=typemax(Int),
    root::Int=1,
    balldata=BoundingBallData,
    updateradii=boundingsphere,
    splitterkwargs=NamedTuple(),
)
    return BoundingBallTreeBuilder(
        splitter,
        numsplits,
        minvalues,
        minlevel,
        maxlevel,
        root,
        balldata,
        updateradii,
        splitterkwargs,
    )
end

"""
    KMeansTreeBuilder(; numberofclusters=4, minvalues=70,
        minlevel=AutoMinLevel(), root=1, balldata=BoundingBallData,
        updateradii=boundingsphere, splitterkwargs=(; rng=Random.MersenneTwister(1234)))

Construction settings for [`KMeansTree`](@ref).

`splitterkwargs` defaults to a seeded `rng`, so `KMeansTree(points)` is
deterministic unless the caller passes different k-means options.
"""
struct KMeansTreeBuilder{L,B,U,K}
    numberofclusters::Int
    minvalues::Int
    minlevel::L
    root::Int
    balldata::B
    updateradii::U
    splitterkwargs::K
end

function KMeansTreeBuilder(;
    numberofclusters::Int=4,
    minvalues::Int=70,
    minlevel=AutoMinLevel(),
    root::Int=1,
    balldata=BoundingBallData,
    updateradii=boundingsphere,
    splitterkwargs=(; rng=Random.MersenneTwister(1234)),
)
    return KMeansTreeBuilder(
        numberofclusters, minvalues, minlevel, root, balldata, updateradii, splitterkwargs
    )
end

"""
    MetisTreeBuilder(; numdivisions=4, splitunconnectedpartitions=false,
        minvalues=numdivisions, minlevel=AutoMinLevel(),
        maxlevel=typemax(Int), root=1, balldata=BoundingBallData,
        updateradii=boundingsphere)

Construction settings for [`MetisTree`](@ref).

`numdivisions` controls the requested METIS split count at each partitioning
step. `splitunconnectedpartitions=true` post-processes disconnected METIS parts
into connected components.
"""
struct MetisTreeBuilder{L,B,U}
    numdivisions::Int
    splitunconnectedpartitions::Bool
    minvalues::Int
    minlevel::L
    maxlevel::Int
    root::Int
    balldata::B
    updateradii::U
end

function MetisTreeBuilder(;
    numdivisions::Int=4,
    splitunconnectedpartitions::Bool=false,
    minvalues::Int=numdivisions,
    minlevel=AutoMinLevel(),
    maxlevel::Int=typemax(Int),
    root::Int=1,
    balldata=BoundingBallData,
    updateradii=boundingsphere,
)
    return MetisTreeBuilder(
        numdivisions,
        splitunconnectedpartitions,
        minvalues,
        minlevel,
        maxlevel,
        root,
        balldata,
        updateradii,
    )
end

"""
    MetisForestBuilder(; treebuilder=MetisTreeBuilder())

Construction settings for [`MetisForest`](@ref). `treebuilder` is the
[`MetisTreeBuilder`](@ref) used for each connected component.
"""
struct MetisForestBuilder{T}
    treebuilder::T
end

function MetisForestBuilder(; treebuilder=MetisTreeBuilder())
    return MetisForestBuilder(treebuilder)
end

"""
    SimpleHybridTreeBuilder(; hybridhalfsize)

Construction settings for [`SimpleHybridTree`](@ref).

`hybridhalfsize` selects the level at which the tree switches from the original
tree representation to the hybrid wrapper.
"""
struct SimpleHybridTreeBuilder{H}
    hybridhalfsize::H
end

function SimpleHybridTreeBuilder(; hybridhalfsize)
    return SimpleHybridTreeBuilder(hybridhalfsize)
end
