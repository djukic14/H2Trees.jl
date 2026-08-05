# Plans: General Introduction

In 𝓗² methods, these plans describe how information moves through the tree hierarchy during matrix-vector products.
It defines where values are aggregated, where translations are applied, and how results are disaggregated back to the required space.

[`buildplans`](@ref) with a [`PlanBuilder`](@ref) is the public entry point: it builds a
Galerkin [`PlanSet`](@ref) for a single tree or a Petrov [`PlanSet`](@ref) for a `BlockTree`,
dispatching automatically on which one `tree` is.

## Core Ingredients

Independent of Galerkin or Petrov formulation, plan construction uses:

- A tree representation for test and trial spaces.
- A translating-nodes iterator that separates near and far interactions.
- A node-selection rule (`aggregatenode`) that marks translating nodes.
- An aggregation mode that chooses the forward aggregation/disaggregation pattern.

## Plan Families

Each constructed [`PlanSet`](@ref) contains these public entries:

- `trialaggregationplan`
- `testdisaggregationplan`
- `testaggregationplan`
- `trialdisaggregationplan`
- `relevantlevels`

For Petrov planning, `mintranslationlevel` is additionally provided.

The entries can be accessed by property, by symbol key, or through the accessor
functions:

```julia
plans.trialaggregationplan
plans[:trialaggregationplan]
H2Trees.trialaggregationplan(plans)
```

`PlanSet` is a concrete return type so dispatch can distinguish Galerkin and
Petrov plan metadata without sentinel symbols. Code that needs raw named-tuple
interoperability can use `NamedTuple(plans)`, keyword splatting with
`(; plans...)`, or `merge` with named tuples in either argument order.

```julia
nt = NamedTuple(plans)
kwargs = (; plans...)
extended = merge(plans, (label=:forward,))
```

!!! note
    Older code that only accessed plan fields by name should continue to work.
    Code that dispatches on `NamedTuple` should convert explicitly with
    `NamedTuple(plans)`.

## Aggregation Modes

The selected mode determines the forward plans.
The transpose plans are then built automatically so transpose operations remain consistent.

### AggregateMode

```@raw html
<table>
    <tr>
        <th style="padding: 8px 18px; text-align: left;">AggregatePlan</th>
        <th style="padding: 8px 18px; text-align: left;">DisaggregateTranslatePlan</th>
    </tr>
    <tr>
        <td style="padding: 12px 18px; vertical-align: top;">
            <img src="../../assets/aggregationtree.svg" alt="AggregatePlan" style="height: 300%;" />
        </td>
        <td style="padding: 12px 18px; vertical-align: top;">
        <img src="../../assets/disaggregationtranslationtree.svg" alt="DisaggregateTranslatePlan" style="height: 300%;" />
        </td>
        </tr>
</table>
```

### AggregateTranslateMode

```@raw html
<table>
    <tr>
        <th style="padding: 8px 18px; text-align: left;">AggregateTranslatePlan</th>
        <th style="padding: 8px 18px; text-align: left;">DisaggregatePlan</th>
    </tr>
    <tr>
        <td style="padding: 12px 18px; vertical-align: top;">
            <img src="../../assets/aggregationtranslationtree.svg" alt="AggregateTranslatePlan" style="height: 300%;" />
        </td>
        <td style="padding: 12px 18px; vertical-align: top;">
            <img src="../../assets/disaggregationtree.svg" alt="DisaggregatePlan" style="height: 300%;" />
        </td>
    </tr>
</table>
```

The two modes produce different forward and adjoint plan assignments:

Forward plans:

| Mode | Trial-side forward plan (`trialaggregationplan`) | Test-side forward plan (`testdisaggregationplan`) |
| --- | --- | --- |
| `AggregateMode()` | [`AggregatePlan`](aggregateplan.md) | [`DisaggregateTranslatePlan`](disaggregatetranslateplan.md) |
| `AggregateTranslateMode()` | [`AggregateTranslatePlan`](aggregatetranslateplan.md) | [`DisaggregatePlan`](disaggregateplan.md) |

Adjoint (transpose) plans:

| Mode | Test-side adjoint plan (`testaggregationplan`) | Trial-side adjoint plan (`trialdisaggregationplan`) |
| --- | --- | --- |
| `AggregateMode()` | [`AggregateTranslatePlan`](aggregatetranslateplan.md) | [`DisaggregatePlan`](disaggregateplan.md) |
| `AggregateTranslateMode()` | [`AggregatePlan`](aggregateplan.md) | [`DisaggregateTranslatePlan`](disaggregatetranslateplan.md) |

## Galerkin and Petrov Context

The main distinction is how test and trial trees are related:

- Galerkin: test tree and trial tree are the same.
- Petrov: test tree and trial tree are generally different.

For full construction workflows and runnable examples, see:

- [Constructing Galerkin Plans](galerkinplans.md)
- [Constructing Petrov Plans](petrovplans.md)
