function splitplan(tree, plan::AggregateTranslatePlan)
    startlevel, endlevel = levels(plan)[end], levels(plan)[begin]
    hybridlevel = H2Trees.hybridlevel(tree)

    upperrange = hybridlevel:-1:startlevel
    lowerrange = endlevel:-1:(hybridlevel + 1)
    return AggregateTranslatePlan(
        plan.receivingnodes[leveltolevelid.(Ref(plan), upperrange)],
        plan.nodes[leveltolevelid.(Ref(plan), upperrange)],
        upperrange,
        plan.rootoffset,
        tree,
    ),
    AggregateTranslatePlan(
        plan.receivingnodes[leveltolevelid.(Ref(plan), lowerrange)],
        plan.nodes[leveltolevelid.(Ref(plan), lowerrange)],
        lowerrange,
        plan.rootoffset,
        tree,
    )
end

function splitplan(tree, plan::DisaggregateTranslatePlan)
    startlevel, endlevel = levels(plan)[begin], levels(plan)[end]
    hybridlevel = H2Trees.hybridlevel(tree)

    upperrange = startlevel:hybridlevel
    lowerrange = (hybridlevel + 1):endlevel

    return DisaggregateTranslatePlan(
        plan.translatingnodes[leveltolevelid.(Ref(plan), upperrange)],
        plan.nodes[leveltolevelid.(Ref(plan), upperrange)],
        upperrange,
        plan.isdisaggregationnode,
        plan.rootoffset,
        plan.tree,
    ),
    DisaggregateTranslatePlan(
        plan.translatingnodes[leveltolevelid.(Ref(plan), lowerrange)],
        plan.nodes[leveltolevelid.(Ref(plan), lowerrange)],
        lowerrange,
        plan.isdisaggregationnode,
        plan.rootoffset,
        plan.tree,
    )
end
