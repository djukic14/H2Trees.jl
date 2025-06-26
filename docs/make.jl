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
    translations

DocMeta.setdocmeta!(H2Trees, :DocTestSetup, :(using H2Trees); recursive=true)

makedocs(;
    modules=[H2Trees],
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
            "Hybrid Tree" => "hybridtree.md",
            "Quadpoints Tree" => "quadpointtree.md",
            "Iterators" => "iterators.md",
            "Plans" => [
                "Aggregate Plan" => "plans/aggregateplan.md",
                "Disaggregate Translate Plan" => "plans/disaggregatetranslateplan.md",
                "Aggregate Translate Plan" => "plans/aggregatetranslateplan.md",
                "Disaggregate Plan" => "plans/disaggregateplan.md",
            ],
            "Translations" => "translations.md",
        ],
        "Visualization" => [
            "TwoNTree" => "visualization/twontree_plot.md",
            "BoundingBallTree" => "visualization/boundingballtree_plot.md",
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
