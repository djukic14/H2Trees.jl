using H2Trees
using CompScienceMeshes, PlotlyJS, ParallelKMeans, BEAST, Metis, Graphs
using Documenter
import H2Trees:
    DepthFirstIterator,
    ParentUpwardsIterator,
    ChildIterator,
    leaves,
    LevelIterator,
    SameLevelIterator,
    NearNodeIterator,
    FarNodeIterator,
    NodeFilterIterator,
    AbstractTranslationTrait,
    AllTranslations,
    DirectionInvariance,
    DirectionInvariancePerLevel,
    translations,
    WellSeparatedIterator,
    TranslatingNodesIterator,
    NearNodeIterator,
    isnear,
    iswellseparated,
    TwoNTree,
    BoundingBallTree,
    tracecube,
    traceball,
    SimpleHybridTree,
    PlanBuilder,
    PlanSet,
    buildplans,
    TwoNTreeBuilder,
    BlockTreeBuilder,
    BoundingBallTreeBuilder,
    KMeansTreeBuilder,
    MetisTreeBuilder,
    MetisForestBuilder,
    SimpleHybridTreeBuilder,
    MetisTree,
    MetisForest

DocMeta.setdocmeta!(H2Trees, :DocTestSetup, :(using H2Trees); recursive=true)

include(joinpath(@__DIR__, "ensureplots.jl"))

makedocs(;
    modules=[
        H2Trees,
        H2Trees.SEBB,
        if isdefined(Base, :get_extension)
            Base.get_extension(H2Trees, :H2MetisTrees)
        else
            H2Trees.H2MetisTrees
        end,
        if isdefined(Base, :get_extension)
            Base.get_extension(H2Trees, :H2BEASTTrees)
        else
            H2Trees.H2BEASTTrees
        end,
        if isdefined(Base, :get_extension)
            Base.get_extension(H2Trees, :H2ParallelKMeansTrees)
        else
            H2Trees.H2ParallelKMeansTrees
        end,
        if isdefined(Base, :get_extension)
            Base.get_extension(H2Trees, :H2PlotlyJSTrees)
        else
            H2Trees.H2PlotlyJSTrees
        end,
    ],
    authors="Danijel Jukić <danijel.jukic14@gmail.com> and contributors",
    sitename="𝓗² Trees.jl",
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://djukic14.github.io/H2Trees.jl",
        edit_link="main",
        assets=String[],
        # The full API reference is one long, legitimately large page (every public and
        # internal docstring across the core package and its four extensions); exempt it
        # from the page-size budget that guards against accidentally-bloated manual pages.
        size_threshold_ignore=["apiref.md"],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Builder Workflow" => "builders.md",
            "Tree Families" => "tree_families.md",
            "TwoNTree" => "twontree.md",
            "BlockTree and Petrov Trees" => "blocktree_petrov.md",
            "Simple Hybrid Tree" => "simplehybridtree.md",
            "Protrusion Policy" => "protrusion.md",
            "Tree Access and Values" => "tree_access.md",
            "Near and Far Predicates" => "near_far.md",
            "Iterators" => "iterators.md",
            "Plans" => [
                "Plans Overview" => "plans/plans.md",
                "Constructing Galerkin Plans" => "plans/galerkinplans.md",
                "Constructing Petrov Plans" => "plans/petrovplans.md",
                "Aggregate Plan" => "plans/aggregateplan.md",
                "Disaggregate Translate Plan" => "plans/disaggregatetranslateplan.md",
                "Aggregate Translate Plan" => "plans/aggregatetranslateplan.md",
                "Disaggregate Plan" => "plans/disaggregateplan.md",
            ],
            "Admissibility Diagnostics" => "admissibility.md",
            "Translations" => "translations.md",
            "Forest" => "forest.md",
            "Printing" => "printing.md",
        ],
        "Hilbert Curve" => "nodeids.md",
        "SEBB" => "sebb.md",
        "Extensions" => [
            "BEAST" => "ext/h2beasttrees.md",
            "Metis" => "ext/h2metistrees.md",
            "ParallelKMeans" => "ext/h2parallelkmeanstrees.md",
            "PlotlyJS" => "ext/h2plotlyjstrees.md",
        ],
        "API Reference" => "apiref.md",
        "Contributing" => "contributing.md",
    ],
)

deploydocs(;
    repo="github.com/djukic14/H2Trees.jl",
    target="build",
    devbranch="main",
    push_preview=true,
    forcepush=true,
    versions=["stable" => "v^", "dev" => "dev"],
)
