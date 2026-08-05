# Detailed multi-line tree/forest summaries (the `text/plain` MIME `show`). Compact one-line
# `Base.show` methods live in `showmethods.jl`; this file only computes/prints the detail.

"""
    averagenumberofchildren(tree, level)

Return the average number of children per node at `level`.
"""
function averagenumberofchildren(tree, level::Int)
    nodesatlevel = LevelIterator(tree, level)
    numberofchildren = 0
    for node in nodesatlevel
        for _ in children(tree, node)
            numberofchildren += 1
        end
    end

    return numberofchildren / length(nodesatlevel)
end

"""
    averagenumberofpoints(tree, level)

Return the average number of stored values per node at `level`.
"""
function averagenumberofpoints(tree, level::Int)
    nodesatlevel = LevelIterator(tree, level)
    numberofpoints = 0
    for node in nodesatlevel
        numberofpoints += numberofvalues(tree, node)
    end

    return numberofpoints / length(nodesatlevel)
end

_treelevelmetriclabel(::isTwoNTree) = "halfsize"
function _treelevelmetric(tree, level::Int, ::isTwoNTree)
    return halfsize(tree, first(LevelIterator(tree, level)))
end

_treelevelmetriclabel(::isBoundingBallTree) = "radius"
function _treelevelmetric(tree, level::Int, ::isBoundingBallTree)
    return radius(tree, first(LevelIterator(tree, level)))
end

"""
    _levelstats(tree, level, trait)

Return the per-level values printed in the detailed tree summary.
"""
function _levelstats(tree, level::Int, trait)
    nodesatlevel = LevelIterator(tree, level)
    numnodes = length(nodesatlevel)
    numberofchildren = 0
    totalvalues = 0
    for node in nodesatlevel
        totalvalues += numberofvalues(tree, node)
        for _ in children(tree, node)
            numberofchildren += 1
        end
    end
    return (
        level=level,
        numnodes=numnodes,
        avgvalues=totalvalues / numnodes,
        avgchildren=numberofchildren / numnodes,
        metric=_treelevelmetric(tree, level, trait),
    )
end

_statstr(x) = string(round(x; digits=2))
_statstr(x::Integer) = string(x)
_statstr(x::AbstractFloat) = string(round(x; digits=2))

function _printcells(io::IO, cells)
    print(io, rpad(string(first(cells)), 8))
    for cell in Iterators.drop(cells, 1)
        print(io, rpad(string(cell), 14))
    end
    return println(io)
end

function _printtableheader(io::IO, metriclabel)
    return _printcells(io, ("level", "nodes", "avg values", "avg children", metriclabel))
end

function _printlevelrow(io::IO, stats)
    return _printcells(
        io,
        (
            stats.level,
            stats.numnodes,
            _statstr(stats.avgvalues),
            _statstr(stats.avgchildren),
            stats.metric,
        ),
    )
end

"""
    printtree(io, tree)

Print the detailed multi-line tree summary used by `show(io, MIME"text/plain"(), tree)`.
"""
function printtree(io::IO, tree)
    return printtree(io, tree, H2Trees.treetrait(tree))
end

function printtree(io::IO, tree, trait::Union{isTwoNTree,isBoundingBallTree})
    println(io, _treeheaderstr(tree))
    _printtableheader(io, _treelevelmetriclabel(trait))
    for level in H2Trees.levels(tree)
        _printlevelrow(io, _levelstats(tree, level, trait))
    end
end

function _foreststats(forest)
    totalnodes, totalleaves = 0, 0
    minlevel, maxlevel = typemax(Int), typemin(Int)
    for tree in forest
        totalnodes += numberofnodes(tree)
        totalleaves += length(leaves(tree))
        lvls = H2Trees.levels(tree)
        isempty(lvls) && continue
        minlevel = min(minlevel, first(lvls))
        maxlevel = max(maxlevel, last(lvls))
    end
    return totalnodes, totalleaves, minlevel, maxlevel
end

_forestlevelstr(minlevel, maxlevel) = minlevel > maxlevel ? "none" : "$minlevel:$maxlevel"

"""
    printforest(io, forest)

Print the detailed multi-line forest summary used by
`show(io, MIME"text/plain"(), forest)`.
"""
function printforest(io::IO, forest)
    totalnodes, totalleaves, minlevel, maxlevel = _foreststats(forest)
    println(
        io,
        "Forest with ",
        length(forest),
        " tree(s): ",
        totalnodes,
        " total node(s), ",
        totalleaves,
        " total leaves, levels ",
        _forestlevelstr(minlevel, maxlevel),
    )
    for (i, tree) in enumerate(forest)
        print(io, "  [", i, "] ")
        show(io, tree)
        i != length(forest) && println(io)
    end
end
