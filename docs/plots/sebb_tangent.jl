using H2Trees
using LinearAlgebra
using StaticArrays
using PlotlyJS

const SEBB = H2Trees.SEBB

# A configuration where EVERY input ball touches the enclosing ball from the inside, in each
# supported dimension. Pick the answer first, then build inputs that must produce it.
#
# Put the target ball at the origin with radius 1 and choose contact directions `u`. Placing an
# input ball of radius `r` at `c = (1 - r) * u` makes it internally tangent, since
# `‖c‖ + r = (1 - r) + r = 1`. The enclosing ball is then exactly the target one provided the
# origin lies in the convex hull of the contact directions. The directions below satisfy that in
# 1D, 2D, and 3D.
function tangentballs(directions, radii)
    return [SEBB.Ball((1 - r) * SVector(u...), r) for (u, r) in zip(directions, radii)]
end

configurations = (
    (label="1D", balls=tangentballs(((1.0,), (-1.0,)), (0.55, 0.40))),
    (
        label="2D",
        balls=tangentballs(
            map(θ -> (cos(θ), sin(θ)), (π / 2, π / 2 + 2π / 3, π / 2 + 4π / 3)),
            (0.45, 0.30, 0.55),
        ),
    ),
    (
        label="3D",
        balls=tangentballs(
            map(
                v -> v ./ sqrt(3),
                ((1.0, 1.0, 1.0), (1.0, -1.0, -1.0), (-1.0, 1.0, -1.0), (-1.0, -1.0, 1.0)),
            ),
            (0.40, 0.50, 0.30, 0.45),
        ),
    ),
)

# The figure claims tangency, so the figure checks it: solve each case and confirm the answer is
# the unit ball at the origin and that every input really does touch its boundary. A drifting
# solver turns this into a build failure rather than a quietly wrong picture.
for configuration in configurations
    enclosing = SEBB.smallest_enclosing_ball(configuration.balls)
    @assert isapprox(SEBB.radius(enclosing), 1; atol=1e-12)
    @assert all(isapprox.(SEBB.center(enclosing), 0; atol=1e-12))
    for ball in configuration.balls
        gap = 1 - (norm(SEBB.center(ball)) + SEBB.radius(ball))
        @assert isapprox(gap, 0; atol=1e-12) "input ball is not on the boundary (gap $gap)"
    end
end

# --- hide-from-docs ---

const ENCLOSINGCOLOR = "rgb(70,110,190)"
const BALLCOLORS = [
    "rgb(214,96,77)", "rgb(244,165,130)", "rgb(146,197,222)", "rgb(33,102,172)"
]

# Keep all dimensions in one 3D scene, so the selector only toggles visibility and camera.
function spheretrace(center, radius, color, opacity, name; showscale=false)
    us = range(0, 2π; length=48)
    vs = range(0, π; length=24)
    x = [center[1] + radius * cos(u) * sin(v) for v in vs, u in us]
    y = [center[2] + radius * sin(u) * sin(v) for v in vs, u in us]
    z = [center[3] + radius * cos(v) for v in vs, u in us]
    return surface(;
        x=x,
        y=y,
        z=z,
        opacity=opacity,
        showscale=showscale,
        colorscale=[[0, color], [1, color]],
        name=name,
        hoverinfo="name",
    )
end

function circletrace(center, radius, color, width, name; dash="solid")
    θ = range(0, 2π; length=181)
    return scatter3d(;
        x=center[1] .+ radius .* cos.(θ),
        y=center[2] .+ radius .* sin.(θ),
        z=fill(0.0, length(θ)),
        mode="lines",
        line=attr(; color=color, width=width, dash=dash),
        name=name,
        hoverinfo="name",
    )
end

function segmenttrace(lo, hi, yoffset, color, width, name)
    return scatter3d(;
        x=[lo, hi],
        y=[yoffset, yoffset],
        z=[0.0, 0.0],
        mode="lines+markers",
        line=attr(; color=color, width=width),
        marker=attr(; size=4, color=color),
        name=name,
        hoverinfo="name",
    )
end

traces = GenericTrace[]
groups = Vector{Int}[]   # trace indices belonging to each configuration

# 1D: balls are intervals; y offsets only separate overlaps.
start1d = length(traces)
push!(traces, segmenttrace(-1.0, 1.0, 0.0, ENCLOSINGCOLOR, 10, "enclosing ball"))
for (i, ball) in enumerate(configurations[1].balls)
    c, r = SEBB.center(ball)[1], SEBB.radius(ball)
    push!(
        traces,
        segmenttrace(
            c - r, c + r, 0.22 * i, BALLCOLORS[i], 7, "ball $i (r = $(round(r; digits=2)))"
        ),
    )
end
push!(groups, collect((start1d + 1):length(traces)))

# 2D: circles in the z = 0 plane.
start2d = length(traces)
push!(traces, circletrace((0.0, 0.0), 1.0, ENCLOSINGCOLOR, 6, "enclosing ball"))
for (i, ball) in enumerate(configurations[2].balls)
    c, r = SEBB.center(ball), SEBB.radius(ball)
    push!(
        traces,
        circletrace(
            (c[1], c[2]), r, BALLCOLORS[i], 4, "ball $i (r = $(round(r; digits=2)))"
        ),
    )
end
push!(groups, collect((start2d + 1):length(traces)))

# 3D: the enclosing sphere is translucent so tangent balls stay visible.
start3d = length(traces)
push!(traces, spheretrace((0.0, 0.0, 0.0), 1.0, ENCLOSINGCOLOR, 0.18, "enclosing ball"))
for (i, ball) in enumerate(configurations[3].balls)
    c, r = SEBB.center(ball), SEBB.radius(ball)
    push!(
        traces,
        spheretrace(
            (c[1], c[2], c[3]),
            r,
            BALLCOLORS[i],
            0.95,
            "ball $i (r = $(round(r; digits=2)))",
        ),
    )
end
push!(groups, collect((start3d + 1):length(traces)))

visibility(active) = [j in groups[active] for j in 1:length(traces)]

# Straight-on cameras for the degenerate cases: down the y-axis for 1D so the intervals lie on a
# line, from directly above for 2D so the circles are undistorted.
const CAMERAS = (
    attr(; eye=attr(; x=0.0, y=-2.2, z=0.01), up=attr(; x=0, y=0, z=1)),
    attr(; eye=attr(; x=0.0, y=0.0, z=2.2), up=attr(; x=0, y=1, z=0)),
    attr(; eye=attr(; x=1.6, y=1.5, z=1.0), up=attr(; x=0, y=0, z=1)),
)

buttons = [
    attr(;
        label=configurations[k].label,
        method="update",
        args=[attr(; visible=visibility(k)), attr(; scene=attr(; camera=CAMERAS[k]))],
    ) for k in (3, 2, 1)
]

# 3D is the default view, so only that group starts visible.
for (j, trace) in enumerate(traces)
    trace[:visible] = j in groups[3]
end

p = PlotlyJS.plot(
    traces,
    Layout(;
        scene=attr(;
            aspectmode="cube",
            camera=CAMERAS[3],
            xaxis=attr(; range=[-1.25, 1.25], title=""),
            yaxis=attr(; range=[-1.25, 1.25], title=""),
            zaxis=attr(; range=[-1.25, 1.25], title=""),
        ),
        margin=attr(; l=0, r=0, t=40, b=0),
        showlegend=false,
        updatemenus=[
            attr(;
                type="buttons",
                direction="right",
                showactive=true,
                x=0.0,
                xanchor="left",
                y=1.12,
                yanchor="top",
                buttons=buttons,
            ),
        ],
    ),
)

p
