# --- hide-from-docs ---
# What one box has to translate FROM, and what a symmetry-reduced collection actually stores.
#
# ONE BOX AT CHILD SLOT (0,0), and that choice is load-bearing rather than cosmetic. The canonical
# representative of an orbit is picked globally, so for a general child slot some representatives
# are NOT themselves offsets of the box being drawn, so the "stored" panel would then contain arrows
# the left panel never showed. At slot (0,0) all seven representatives lie inside this box's own
# 27 offsets, so the right panel is a genuine SUBSET of the left one. (Measured: the other slots
# give 3, 6 and 6 of 7 inside.)
#
# `FullLatticeSymmetry` for a second reason of the same kind: one box's list meets every orbit the
# whole LEVEL has, so 7 is simultaneously this box's count and the level-wide stored count. Under
# `OppositeSymmetry` those differ (19 against 20), and the figure would need an asterisk.
using H2Trees
using PlotlyJS

const NEAR2D = Set(Iterators.product(-1:1, -1:1))

# The interaction list of a box at child slot `slot`: the children of the parent's neighbours,
# minus the box's own neighbours. In units of the box halfsize.
function boxoffsets2d(slot)
    offsets = NTuple{2,Int}[]
    for parent in Iterators.product(-1:1, -1:1), child in Iterators.product(0:1, 0:1)
        q = ntuple(i -> 2parent[i] + child[i] - slot[i], 2)
        q in NEAR2D || push!(offsets, q)
    end
    return sort!(offsets)
end

group2d = H2Trees.symmetrygroup(H2Trees.FullLatticeSymmetry(), Val(2))
offsets2d = boxoffsets2d((0, 0))
canonical2d = [first(H2Trees.canonicalizetranslation(q, group2d)) for q in offsets2d]
representatives2d = unique(canonical2d)

# One colour per orbit, shared by both panels: an arrow on the left carries the colour of the
# stored arrow on the right that serves it.
palette = [
    "#1f77b4",
    "#d62728",
    "#2ca02c",
    "#9467bd",
    "#ff7f0e",
    "#17becf",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
]
orbitcolor2d = Dict(
    r => palette[mod1(i, length(palette))] for (i, r) in enumerate(representatives2d)
)

# A box outline, centered on `q`, in halfsize units (so a box spans one unit each way).
function boxtrace2d(q; color, width=1, dash="solid", name="", showlegend=false)
    x = q[1] .+ [-1, 1, 1, -1, -1]
    y = q[2] .+ [-1, -1, 1, 1, -1]
    return scatter(;
        x=x,
        y=y,
        mode="lines",
        line=attr(; color=color, width=width, dash=dash),
        name=name,
        showlegend=showlegend,
        hoverinfo="skip",
    )
end

function arrowtrace2d(q; color, width, opacity, col)
    # Centre to centre: the arrow runs from the receiving box's midpoint to the midpoint of the box
    # it translates from, which is what the offset `q` actually measures. Boxes are drawn at `2q`
    # with half-width 1, so `2q` is the target's centre.
    return scatter(;
        x=[0, 2q[1]],
        y=[0, 2q[2]],
        mode="lines+markers",
        line=attr(; color=color, width=width),
        marker=attr(; size=[0, 6], color=color, symbol="arrow-up", angleref="previous"),
        opacity=opacity,
        showlegend=false,
        hovertemplate="offset ($(q[1]), $(q[2]))<extra></extra>",
        xaxis=col == 1 ? "x" : "x2",
        yaxis=col == 1 ? "y" : "y2",
    )
end

traces = GenericTrace[]

# LEFT: every translation this box needs, which is what `DirectionInvariancePerLevel` stores for it.
for (q, c) in zip(offsets2d, canonical2d)
    t = boxtrace2d(2 .* q; color=orbitcolor2d[c], width=1)
    t[:opacity] = 0.35
    t[:xaxis], t[:yaxis] = "x", "y"
    push!(traces, t)
    push!(traces, arrowtrace2d(q; color=orbitcolor2d[c], width=1, opacity=0.5, col=1))
end

# RIGHT: only the stored representatives, one per orbit, each serving every arrow of its colour
# on the left.
for r in representatives2d
    t = boxtrace2d(2 .* r; color=orbitcolor2d[r], width=2)
    t[:xaxis], t[:yaxis] = "x2", "y2"
    push!(traces, t)
    push!(traces, arrowtrace2d(r; color=orbitcolor2d[r], width=3, opacity=1.0, col=2))
end

# THE SAME 27 BOXES, FAINT, ON THE RIGHT TOO. Without them the reduced panel is seven arrows into
# empty space and the eye has nothing to compare against; with them the stored set is visibly a
# SUBSET of the interaction list rather than a different picture.
for q in offsets2d
    t = boxtrace2d(2 .* q; color="#dddddd", width=1)
    t[:xaxis], t[:yaxis] = "x2", "y2"
    pushfirst!(traces, t)
end

# The receiving box itself, and its near-field ring, on both panels for orientation.
for (xa, ya) in (("x", "y"), ("x2", "y2"))
    for n in NEAR2D
        n == (0, 0) && continue
        t = boxtrace2d(2 .* n; color="#cccccc", width=1, dash="dot")
        t[:xaxis], t[:yaxis] = xa, ya
        push!(traces, t)
    end
    t = boxtrace2d((0, 0); color="#111111", width=3)
    t[:xaxis], t[:yaxis] = xa, ya
    push!(traces, t)
end

# ONE EXPLICIT RANGE FOR BOTH PANELS. Left autoscales over all 27 offsets, right over only the
# seven representatives, and canonicalization maps every orbit into ONE octant, so the reduced
# set is not centred on the box. Left to themselves the two panels come out at different scales and
# different origins, which is precisely the comparison this figure is supposed to make.
# `q` spans -2:3 per axis and a box is drawn at `2q` with half-width 1, so the data occupies
# [-5, 7]; one box of margin either side.
const EXTENT2D = [-6, 8]
axis2d = attr(; showgrid=false, zeroline=false, showticklabels=false, range=EXTENT2D)
p = PlotlyJS.plot(
    traces,
    Layout(;
        title=attr(;
            text="one box's translations (left, 27) and what a lattice symmetry stores (right, 7)",
            x=0.5,
        ),
        xaxis=attr(; axis2d..., domain=[0.0, 0.47], scaleanchor="y", constrain="domain"),
        yaxis=attr(; axis2d...),
        xaxis2=attr(; axis2d..., domain=[0.53, 1.0], scaleanchor="y2", constrain="domain"),
        yaxis2=attr(; axis2d..., anchor="x2"),
        margin=attr(; l=10, r=10, t=50, b=10),
        plot_bgcolor="white",
        paper_bgcolor="white",
    ),
)
# --- end-hide-from-docs ---
p
