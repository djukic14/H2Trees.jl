# returns the average number of children per node at a given level
function averagenumberofchildrens(tree, level::Int)
    numberofchildren = 0
    for node in LevelIterator(tree, level)
        numberofchildren += length(collect(children(tree, node)))
    end

    return numberofchildren / length(collect(LevelIterator(tree, level)))
end

# returns the average number of points per node at a given level
function averagenumberofpoints(tree, level::Int)
    numberofpoints = 0
    for node in LevelIterator(tree, level)
        numberofpoints += length(H2Trees.values(tree, node))
    end

    return numberofpoints / length(collect(LevelIterator(tree, level)))
end

function printtree(io::IO, tree)
    return printtree(io, tree, H2Trees.treetrait(tree))
end

function printtree(io::IO, tree, ::isTwoNTree)
    println(io, typeof(tree))
    for level in H2Trees.levels(tree)
        avgnpoints = round(averagenumberofpoints(tree, level); digits=2)
        avgnchildren = round(averagenumberofchildrens(tree, level); digits=2)
        hs = halfsize(tree, H2Trees.LevelIterator(tree, level)[begin])
        numnodes = length(collect(H2Trees.LevelIterator(tree, level)))
        print(io, "-"^(level - H2Trees.minimumlevel(tree)))
        println(
            io,
            " level: $level with $numnodes node(s) with on average $avgnpoints points and $avgnchildren children and halfsize: $hs",
        )
    end
end

function printtree(io::IO, tree, ::isBoundingBallTree)
    println(io, typeof(tree))
    for level in H2Trees.levels(tree)
        avgnpoints = averagenumberofpoints(tree, level)
        avgnchildren = averagenumberofchildrens(tree, level)
        radius = H2Trees.radius(tree, H2Trees.LevelIterator(tree, level)[begin])
        print(io, "-"^(level - H2Trees.minimumlevel(tree)))
        println(
            io,
            " level: $level with on average $avgnpoints points and $avgnchildren childrens and radius: $radius",
        )
    end
end

function Base.show(io::IO, mime::MIME"text/plain", tree::H2ClusterTree)
    return printtree(io, tree)
end

function Base.show(io::IO, tree::H2ClusterTree)
    return printtree(io, tree)
end
