"""
    splitplan(tree, plan)

Split a translating plan into two subplans at `hybridlevel(tree)`.

The function returns `(upperplan, lowerplan)`, both with the same concrete plan
type as `plan`, where level sets are partitioned into levels at or above the
hybrid level and levels below it.
"""
function splitplan(tree, plan::AggregateTranslatePlan)
    startlevel, endlevel = levels(plan)[end], levels(plan)[begin]
    hybridlevel_ = hybridlevel(tree)

    upperrange = hybridlevel_:-1:startlevel
    lowerrange = endlevel:-1:(hybridlevel_ + 1)
    upperids = leveltolevelid.(Ref(plan), upperrange)
    lowerids = leveltolevelid.(Ref(plan), lowerrange)

    return AggregateTranslatePlan(
        plan.receivingnodes[upperids],
        plan.receivingnodes_by_level[upperids],
        plan.nodes[upperids],
        upperrange,
        plan.istranslatingnode,
        plan.rootoffset,
        tree,
    ),
    AggregateTranslatePlan(
        plan.receivingnodes[lowerids],
        plan.receivingnodes_by_level[lowerids],
        plan.nodes[lowerids],
        lowerrange,
        plan.istranslatingnode,
        plan.rootoffset,
        tree,
    )
end

function splitplan(tree, plan::DisaggregateTranslatePlan)
    startlevel, endlevel = levels(plan)[begin], levels(plan)[end]
    hybridlevel_ = hybridlevel(tree)

    upperrange = startlevel:hybridlevel_
    lowerrange = (hybridlevel_ + 1):endlevel
    upperids = leveltolevelid.(Ref(plan), upperrange)
    lowerids = leveltolevelid.(Ref(plan), lowerrange)

    return DisaggregateTranslatePlan(
        plan.translatingnodes[upperids],
        plan.receivingnodes_by_level[upperids],
        plan.nodes[upperids],
        upperrange,
        plan.isdisaggregationnode,
        plan.rootoffset,
        plan.tree,
    ),
    DisaggregateTranslatePlan(
        plan.translatingnodes[lowerids],
        plan.receivingnodes_by_level[lowerids],
        plan.nodes[lowerids],
        lowerrange,
        plan.isdisaggregationnode,
        plan.rootoffset,
        plan.tree,
    )
end
