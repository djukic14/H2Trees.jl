
"""
    NearNodeIterator(tree, node::Int; isnear=isnear)

Returns an iterator over the nodes in the tree that are at the same level as the specified
`node` and are near to `node`.
Two nodes are near if the function `isnear(tree, nodea, nodeb)` evaluates to true.

# Arguments

  - `tree`: The tree to search for near nodes.
  - `node::Int`: The node from which to start the search.
  - `isnear`: A function that returns `true` if two nodes are near each other. Defaults to `isnear`.

# Returns

An iterator over the nodes in the tree that are at the same level as the specified `node` and are near to `node`.
"""
function NearNodeIterator(tree, node::Int; isnear=isnear)
    return NodeFilterIterator(tree, node, isnear)
end

"""
    NearNodeIterator(testtree, trialtree, trialnode::Int; isnear=isnear)

Returns an iterator over the nodes in the `testtree` that are at the same level as the specified
`node` in the `trialtree` and are near to `node`.
Two nodes are near if the function `isnear(testtree, trialtree, testnode, trialnode)` evaluates to true.

# Arguments

  - `testtree`: The tree to search for near nodes.
  - `trialtree`: The tree that contains the node to find near nodes for.
  - `trialnode::Int`: The node in the `trialtree` from which to start the search.
  - `isnear`: A function that returns `true` if two nodes are near each other. Defaults to `isnear`.

# Returns

An iterator over the nodes in the `testtree` that are at the same level as the specified `node` in the `trialtree` and are near to `node`.
"""
function NearNodeIterator(testtree, trialtree, trialnode::Int; isnear=isnear)
    return NodeFilterIterator(testtree, trialtree, trialnode, isnear)
end

"""
    See NearNodeIterator.
"""
function FarNodeIterator(tree, node::Int; isfar=isfar)
    return NodeFilterIterator(tree, node, isfar)
end

function FarNodeIterator(testtree, trialtree, trialnode::Int; isfar=isfar)
    return NodeFilterIterator(testtree, trialtree, trialnode, isfar)
end

struct _LeafNearFunctor{IN}
    isnear::IN
end

function (f::_LeafNearFunctor)(tree, nodea, nodeb)
    return isleaf(tree, nodea) && f.isnear(tree, nodea, nodeb)
end

function (f::_LeafNearFunctor)(testtree, trialtree, testnode, trialnode)
    return isleaf(testtree, testnode) && f.isnear(testtree, trialtree, testnode, trialnode)
end

function _getindicesstorage(::Val{:flattened})
    return Int[]
end

function _storeindices!(indices, val, ::Val{:flattened})
    return append!(indices, val)
end

function nearnodevalues(tree, node::Int; isnear=isnear, storevalues=Val{:flattened}())
    indices = _getindicesstorage(storevalues)
    for nearnode in NearNodeIterator(tree, node; isnear=isnear)
        _storeindices!(indices, values(tree, nearnode), storevalues)
    end

    isleafnear = _LeafNearFunctor(isnear)

    for parent in ParentUpwardsIterator(tree, node)
        for nearnode in NearNodeIterator(tree, parent; isnear=isleafnear)
            _storeindices!(indices, values(tree, nearnode), storevalues)
        end
    end
    return indices
end

function nearnodevalues(
    testtree, trialtree, trialnode::Int; isnear=isnear, storevalues=Val{:flattened}()
)
    indices = _getindicesstorage(storevalues)
    for nearnode in NearNodeIterator(testtree, trialtree, trialnode; isnear=isnear)
        _storeindices!(indices, values(testtree, nearnode), storevalues)
    end

    isleafnear = _LeafNearFunctor(isnear)

    for parent in ParentUpwardsIterator(trialtree, trialnode)
        for nearnode in NearNodeIterator(testtree, trialtree, parent; isnear=isleafnear)
            _storeindices!(indices, values(testtree, nearnode), storevalues)
        end
    end
    return indices
end

function farnodevalues(tree, node::Int; isfar=isfar)
    indices = Int[]
    for farnode in FarNodeIterator(tree, node; isfar=isfar)
        append!(indices, values(tree, farnode))
    end

    isleaffar = _LeafNearFunctor(isfar)

    for parent in ParentUpwardsIterator(tree, node)
        for farnode in FarNodeIterator(tree, parent; isfar=isleaffar)
            append!(indices, values(tree, farnode))
        end
    end
    return indices
end

function farnodevalues(testtree, trialtree, trialnode::Int; isfar=isfar)
    indices = Int[]
    for farnode in FarNodeIterator(testtree, trialtree, trialnode; isfar=isfar)
        append!(indices, values(testtree, farnode))
    end

    isleaffar = _LeafNearFunctor(isfar)

    for parent in ParentUpwardsIterator(trialtree, trialnode)
        for farnode in FarNodeIterator(testtree, trialtree, parent; isfar=isleaffar)
            append!(indices, values(testtree, farnode))
        end
    end
    return indices
end
