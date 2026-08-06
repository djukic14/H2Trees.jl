# All `Base.show` methods live here: compact one-line for builders, trees, forests, and
# plan/translation results, plus the `text/plain` MIME overrides that switch trees/forests to the
# detailed multi-line summary. The detail itself (`printtree` and its helpers) lives in
# `printing.jl`. Included last so every type it dispatches on is already defined.

# Generic fallback for any `H2ClusterTree` without its own compact `show` below.
function Base.show(io::IO, tree::H2ClusterTree)
    return printtree(io, tree)
end

function Base.show(io::IO, ::MIME"text/plain", tree::H2ClusterTree)
    return printtree(io, tree)
end

_minlevelstr(::AutoMinLevel) = "auto"
_minlevelstr(l::Int) = string(l)

_protrusionstr(::NoProtrusionCheck) = "none"
_protrusionstr(::AutoProtrusionCheck) = "auto"
_protrusionstr(p::ProtrusionCheck) = "max=$(p.max)"

# -- Builders -------------------------------------------------------------------------------
function Base.show(io::IO, b::TwoNTreeBuilder)
    return print(
        io,
        "TwoNTreeBuilder(minhalfsize=",
        b.minhalfsize,
        ", minvalues=",
        b.minvalues,
        ", minlevel=",
        _minlevelstr(b.minlevel),
        ", root=",
        b.root,
        ", protrusion=",
        _protrusionstr(b.protrusion),
        ")",
    )
end

function Base.show(io::IO, b::BlockTreeBuilder)
    return print(io, "BlockTreeBuilder(test=", b.test, ", trial=", b.trial, ")")
end

function Base.show(io::IO, b::KMeansTreeBuilder)
    return print(
        io,
        "KMeansTreeBuilder(numberofclusters=",
        b.numberofclusters,
        ", minvalues=",
        b.minvalues,
        ", minlevel=",
        _minlevelstr(b.minlevel),
        ", root=",
        b.root,
        ")",
    )
end

function Base.show(io::IO, b::MetisTreeBuilder)
    return print(
        io,
        "MetisTreeBuilder(numdivisions=",
        b.numdivisions,
        ", minvalues=",
        b.minvalues,
        ", splitunconnected=",
        b.splitunconnectedpartitions,
        ", minlevel=",
        _minlevelstr(b.minlevel),
        ", root=",
        b.root,
        ")",
    )
end

function Base.show(io::IO, b::MetisForestBuilder)
    return print(io, "MetisForestBuilder(treebuilder=", b.treebuilder, ")")
end

function Base.show(io::IO, b::BoundingBallTreeBuilder)
    return print(
        io,
        "BoundingBallTreeBuilder(numsplits=",
        b.numsplits,
        ", minvalues=",
        b.minvalues,
        ", minlevel=",
        _minlevelstr(b.minlevel),
        ", root=",
        b.root,
        ")",
    )
end

function Base.show(io::IO, b::SimpleHybridTreeBuilder)
    return print(io, "SimpleHybridTreeBuilder(hybridhalfsize=", b.hybridhalfsize, ")")
end

# -- Trees ----------------------------------------------------------------------------------
function _rangestr(r)
    isempty(r) && return "none"
    return string(first(r), ":", last(r))
end

_treelevelstr(tree) = _rangestr(levels(tree))

function _treeheaderstr(tree::Union{TwoNTree{N,D,T},BoundingBallTree{N,D,T}}) where {N,D,T}
    return string(nameof(typeof(tree)), "{", N, ",", T, "}")
end

function Base.show(io::IO, tree::Union{TwoNTree,BoundingBallTree})
    return print(
        io,
        _treeheaderstr(tree),
        "(nodes=",
        numberofnodes(tree),
        ", leaves=",
        length(leaves(tree)),
        ", levels=",
        _treelevelstr(tree),
        ", root=",
        root(tree),
        ")",
    )
end

function Base.show(io::IO, tree::SimpleHybridTree)
    print(io, "SimpleHybridTree(hybridlevel=", hybridlevel(tree), ", tree=")
    show(io, tree.tree)
    return print(io, ")")
end

function Base.show(io::IO, tree::BlockTree)
    print(io, "BlockTree(test=")
    show(io, testtree(tree))
    print(io, ", trial=")
    show(io, trialtree(tree))
    return print(io, ")")
end

# -- Forests --------------------------------------------------------------------------------
function Base.show(io::IO, forest::Forest)
    totalnodes, totalleaves, minlevel, maxlevel = _foreststats(forest)
    return print(
        io,
        "Forest(trees=",
        length(forest),
        ", nodes=",
        totalnodes,
        ", leaves=",
        totalleaves,
        ", levels=",
        _forestlevelstr(minlevel, maxlevel),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", forest::Forest)
    return printforest(io, forest)
end

# -- Plans ----------------------------------------------------------------------------------
function Base.show(io::IO, plans::PlanSet)
    family = ispetrov(plans) ? "Petrov" : "Galerkin"
    print(io, "PlanSet{", family, "}(relevantlevels=", _rangestr(plans.relevantlevels))
    if ispetrov(plans)
        print(io, ", mintranslationlevel=", mintranslationlevel(plans))
    end
    return print(io, ")")
end
