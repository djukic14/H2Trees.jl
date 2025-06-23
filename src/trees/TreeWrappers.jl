macro treewrapper(treewrapper)
    quote
        function H2Trees.treetrait(tree::$treewrapper)
            return H2Trees.treetrait(tree.tree)
        end

        function H2Trees.nodesatlevel(tree::$treewrapper)
            return H2Trees.nodesatlevel(tree.tree)
        end

        function H2Trees.samelevelnodes(tree::$treewrapper, node::Int)
            return H2Trees.samelevelnodes(tree.tree, node)
        end

        function H2Trees.root(tree::$treewrapper)
            return H2Trees.root(tree.tree)
        end

        function (tree::$treewrapper)(node::Int)
            return tree.tree(node)
        end

        function H2Trees.center(tree::$treewrapper, node::Int)
            return H2Trees.center(tree.tree, node)
        end

        function H2Trees.halfsize(tree::$treewrapper, node::Int)
            return H2Trees.halfsize(tree.tree, node)
        end

        function H2Trees.levels(tree::$treewrapper)
            return H2Trees.levels(tree.tree)
        end

        function H2Trees.leaves(tree::$treewrapper)
            return H2Trees.leaves(tree.tree)
        end

        function H2Trees.numberoflevels(tree::$treewrapper)
            return H2Trees.numberoflevels(tree.tree)
        end

        function H2Trees.values(tree::$treewrapper)
            return H2Trees.values(tree.tree)
        end

        function H2Trees.sector(tree::$treewrapper, node::Int)
            return H2Trees.sector(tree.tree, node)
        end

        function H2Trees.data(tree::$treewrapper, node::Int)
            return H2Trees.data(tree.tree, node)
        end

        function H2Trees.parent(tree::$treewrapper, node::Int)
            return H2Trees.parent(tree.tree, node)
        end

        function H2Trees.nextsibling(tree::$treewrapper, node::Int)
            return H2Trees.nextsibling(tree.tree, node)
        end

        function H2Trees.firstchild(tree::$treewrapper, node::Int)
            return H2Trees.firstchild(tree.tree, node)
        end

        function H2Trees.children(tree::$treewrapper, node::Int)
            return H2Trees.children(tree.tree, node)
        end

        function H2Trees.numberofnodes(tree::$treewrapper)
            return H2Trees.numberofnodes(tree.tree)
        end

        function Base.eltype(tree::$treewrapper)
            return Base.eltype(tree.tree)
        end

        function H2Trees.parentcenterminuschildcenter(tree::$treewrapper, node::Int)
            return H2Trees.parentcenterminuschildcenter(tree.tree, node)
        end

        function H2Trees.oppositesector(tree::$treewrapper, node::Int)
            return H2Trees.oppositesector(tree.tree, node)
        end
    end
end
