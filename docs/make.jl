using H2Trees
using CompScienceMeshes, PlotlyJS, ParallelKMeans, BEAST
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
    QuadPointsTree,
    SimpleHybridTree

DocMeta.setdocmeta!(H2Trees, :DocTestSetup, :(using H2Trees); recursive=true)

makedocs(;
    modules=[
        H2Trees,
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
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "TwoNTree" => "twontree.md",
            "Simple Hybrid Tree" => "simplehybridtree.md",
            "Iterators" => "iterators.md",
            "Plans" => [
                "Aggregate Plan" => "plans/aggregateplan.md",
                "Disaggregate Translate Plan" => "plans/disaggregatetranslateplan.md",
                "Aggregate Translate Plan" => "plans/aggregatetranslateplan.md",
                "Disaggregate Plan" => "plans/disaggregateplan.md",
            ],
            "Translations" => "translations.md",
        ],
        "Extensions" => [
            "BEAST" => "ext/h2beasttrees.md",
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
)
