using Test
using H2Trees

@testset verbose = true "H2Trees" begin
    @testset "FarMulMode" begin
        @test H2Trees.AggregateMode() == adjoint(H2Trees.AggregateTranslateMode())
        @test H2Trees.AggregateTranslateMode() == adjoint(H2Trees.AggregateMode())
    end

    @testset verbose = true "Hilbert ordering" begin
        include("hilbert/runtests.jl")
    end

    include("sebb/runtests.jl")

    @testset verbose = true "SEBB tree integration" begin
        include("trees/test_boundingballtree_sebb.jl")
    end

    @testset verbose = true "Iterators" begin
        include("trees/test_iterators.jl")
    end

    @testset verbose = true "Near interactions" begin
        include("test_nearinteractions.jl")
    end

    @testset verbose = true "Near lists" begin
        include("test_nearlists.jl")
    end

    include("trees/runtests.jl")
    include("plans/runtests.jl")

    @testset verbose = true "Translations" begin
        include("translations/test_translations.jl")
    end

    @testset verbose = true "Plots" begin
        include("H2PlotlyJSTrees/test_plots.jl")
    end

    @testset verbose = true "H2ParallelKMeansTrees" begin
        include("H2ParallelKMeansTrees/test_kmeanstree.jl")
    end

    @testset verbose = true "Simple Hybrid Tree" begin
        include("simplehybridtree_runtests.jl")
    end

    @testset verbose = true "Performance contracts" begin
        include("performance/runtests.jl")
    end

    @testset verbose = true "H2BEASTTrees" begin
        include("H2BEASTTrees/runtests.jl")
    end

    @testset verbose = true "H2MetisTrees" begin
        include("H2MetisTrees/runtests.jl")
    end

    @testset "Code quality (Aqua.jl)" begin
        using Aqua
        Aqua.test_all(H2Trees)
    end

    @testset "Code formatting (JuliaFormatter.jl)" begin
        using JuliaFormatter
        @test JuliaFormatter.format(pkgdir(H2Trees), overwrite=false)
    end

    @testset "Explicit imports (ExplicitImports.jl)" begin
        using ExplicitImports
        @test ExplicitImports.check_no_stale_explicit_imports(H2Trees) === nothing
        @test ExplicitImports.check_all_explicit_imports_via_owners(H2Trees) === nothing
        # `@treewrapper` (src/trees/treewrappers.jl) generates forwarding methods inside a
        # `quote` block; unqualified names there resolve through macro hygiene, which silently
        # fails to attach several of these as methods of the intended generic (confirmed by
        # `methods(H2Trees.root)` missing a SimpleHybridTree method after de-qualifying), so its
        # `H2Trees.foo(...)` qualifications must stay. `iswellseparated` is qualified in
        # `WellSeparatedIterator.jl` because the unqualified name there resolves to that
        # function's own `iswellseparated` keyword argument, not the global predicate.
        selfqualifiedignore = (
            :treetrait,
            :nodesatlevel,
            :treeindex,
            :samelevelnodes,
            :root,
            :center,
            :halfsize,
            :levels,
            :leaves,
            :numberoflevels,
            :values,
            :sector,
            :data,
            :parent,
            :nextsibling,
            :firstchild,
            :children,
            :numberofnodes,
            :parentcenterminuschildcenter,
            :oppositesector,
            :iswellseparated,
        )
        @test ExplicitImports.check_no_self_qualified_accesses(
            H2Trees; ignore=selfqualifiedignore
        ) === nothing
    end
end
