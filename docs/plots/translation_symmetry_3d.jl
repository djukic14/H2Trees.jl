# --- hide-from-docs ---
# The three-dimensional twin of `translation_symmetry_2d.jl`. Read that file's header first: the
# choice of child slot (0,0,0) and of `FullLatticeSymmetry` is load-bearing for the same two
# reasons, and holds here too. All 16 representatives lie inside this box's own 189 offsets
# (other slots give 6, 12 and 15 of 16), and one box's list meets every orbit the level has, so 16
# is both this box's count and the level-wide stored count.
#
# INTERACTIVE RATHER THAN STATIC, and that is why this is a PlotlyJS figure and not an image: 189
# offsets form a hollow shell (the 3x3x3 near field is cut out of a 7x7x7 block), so in any fixed
# projection the front and back faces overlap into a solid mass. Rotating it is what makes the
# shell, and the sparseness of the 16 stored directions inside it, legible at all.
#
# Boxes are drawn as MARKERS, not wireframe cubes: 189 cubes is 189 x 12 line segments, which is
# both slow to render and visually solid. The near-field ring the 2-D figure draws is omitted for
# the same reason.
using H2Trees
using PlotlyJS

const NEAR3D = Set(Iterators.product(-1:1, -1:1, -1:1))

function boxoffsets3d(slot)
    offsets = NTuple{3,Int}[]
    for parent in Iterators.product(-1:1, -1:1, -1:1),
        child in Iterators.product(0:1, 0:1, 0:1)

        q = ntuple(i -> 2parent[i] + child[i] - slot[i], 3)
        q in NEAR3D || push!(offsets, q)
    end
    return sort!(offsets)
end

group3d = H2Trees.symmetrygroup(H2Trees.FullLatticeSymmetry(), Val(3))
offsets3d = boxoffsets3d((0, 0, 0))
canonical3d = [first(H2Trees.canonicalizetranslation(q, group3d)) for q in offsets3d]
representatives3d = unique(canonical3d)

palette3d = [
    "#1f77b4",
    "#d62728",
    "#2ca02c",
    "#9467bd",
    "#ff7f0e",
    "#17becf",
    "#8c564b",
    "#e377c2",
    "#bcbd22",
    "#393b79",
    "#637939",
    "#8c6d31",
    "#843c39",
    "#7b4173",
    "#3182bd",
    "#31a354",
]
orbitcolor3d = Dict(
    r => palette3d[mod1(i, length(palette3d))] for (i, r) in enumerate(representatives3d)
)

# One `scatter3d` per orbit rather than per offset: 189 traces would make the figure sluggish in a
# browser, and grouping by orbit is also what the colours are trying to say.
function orbittraces(qs, canonicals, colorof; scene, arrows::Bool, opacity, width)
    traces = GenericTrace[]
    for r in unique(canonicals)
        members = [q for (q, c) in zip(qs, canonicals) if c == r]
        push!(
            traces,
            scatter3d(;
                x=[2q[1] for q in members],
                y=[2q[2] for q in members],
                z=[2q[3] for q in members],
                mode="markers",
                marker=attr(; size=arrows ? 5 : 3, color=colorof[r]),
                opacity=opacity,
                showlegend=false,
                scene=scene,
                hovertemplate="offset (%{x:.0f}, %{y:.0f}, %{z:.0f}) / 2<extra></extra>",
            ),
        )
        arrows || continue
        # Arrows only on the reduced panel: 189 of them would be the hairball this figure exists
        # to avoid. `nothing` separators keep one trace per orbit rather than one per segment.
        xs, ys, zs = Float64[], Float64[], Float64[]
        for q in members
            append!(xs, [0.0, 2q[1], NaN])
            append!(ys, [0.0, 2q[2], NaN])
            append!(zs, [0.0, 2q[3], NaN])
        end
        push!(
            traces,
            scatter3d(;
                x=xs,
                y=ys,
                z=zs,
                mode="lines",
                line=attr(; color=colorof[r], width=width),
                showlegend=false,
                scene=scene,
                hoverinfo="skip",
            ),
        )
    end
    return traces
end

traces = GenericTrace[]
append!(
    traces,
    orbittraces(
        offsets3d,
        canonical3d,
        orbitcolor3d;
        scene="scene",
        arrows=false,
        opacity=0.55,
        width=1,
    ),
)
append!(
    traces,
    orbittraces(
        representatives3d,
        representatives3d,
        orbitcolor3d;
        scene="scene2",
        arrows=true,
        opacity=1.0,
        width=4,
    ),
)

# The receiving box at the origin, on both panels.
for sc in ("scene", "scene2")
    push!(
        traces,
        scatter3d(;
            x=[0],
            y=[0],
            z=[0],
            mode="markers",
            marker=attr(; size=9, color="#111111", symbol="square"),
            showlegend=false,
            scene=sc,
            hovertemplate="receiving box<extra></extra>",
        ),
    )
end

# THE FULL SHELL, FAINT, BEHIND THE STORED ONES. Sixteen arrows into an empty scene say nothing
# about what they replaced; drawn inside the 189 they came from, the sparseness is the message.
push!(
    traces,
    scatter3d(;
        x=[2q[1] for q in offsets3d],
        y=[2q[2] for q in offsets3d],
        z=[2q[3] for q in offsets3d],
        mode="markers",
        marker=attr(; size=2, color="#cccccc"),
        opacity=0.5,
        showlegend=false,
        scene="scene2",
        hoverinfo="skip",
    ),
)

# ONE EXPLICIT RANGE FOR BOTH SCENES, for the reason the 2-D figure spells out: canonicalization
# maps every orbit into one octant, so the reduced scene would otherwise autoscale to that octant
# and put the receiving box somewhere else entirely. `q` spans -2:3 per axis, drawn at `2q`.
const EXTENT3D = [-5, 7]
blank = attr(;
    showgrid=false,
    zeroline=false,
    showticklabels=false,
    title=attr(; text=""),
    range=EXTENT3D,
)
function sceneattr(domain)
    return attr(;
        domain=domain,
        xaxis=blank,
        yaxis=blank,
        zaxis=blank,
        aspectmode="cube",
        camera=attr(; eye=attr(; x=1.5, y=1.5, z=1.1)),
    )
end

p = PlotlyJS.plot(
    traces,
    Layout(;
        title=attr(;
            text="one box's translations (left, 189) and what a lattice symmetry stores (right, 16)",
            x=0.5,
        ),
        scene=sceneattr(attr(; x=[0.0, 0.48], y=[0.0, 1.0])),
        scene2=sceneattr(attr(; x=[0.52, 1.0], y=[0.0, 1.0])),
        margin=attr(; l=0, r=0, t=50, b=0),
        plot_bgcolor="white",
        paper_bgcolor="white",
    ),
)
# --- end-hide-from-docs ---
p
