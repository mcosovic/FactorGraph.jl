# [LDPC-Style Decoding](@id ldpc-decoding-example)

```@meta
CurrentModule = FactorGraph
```

This example builds a small low-density parity-check style decoding problem. The
transmitted codeword is observed through a binary symmetric channel. Unary factors encode
the channel likelihoods, and parity-check factors enforce even parity constraints.

The graph contains cycles, so iterative sum-product belief propagation is used to estimate the posterior
probability of each bit.

---

## Variable Nodes

Each code bit is a binary discrete variable with state order `[0, 1]`:

```@example ldpc_decoding
using FactorGraph

states = [0, 1]
bitIds = [:x1, :x2, :x3, :x4, :x5, :x6]

variables = [
    DiscreteVariable(bitId, 2; label = string(bitId), states = states)
    for bitId in bitIds
]

nothing # hide
```

---

## Channel Factors

Assume a binary symmetric channel with bit-flip probability `p`. The example transmits a
valid nonzero codeword and flips one received bit. If the received bit is `0`, the
likelihood over states `[0, 1]` is `[1 - p, p]`; if it is `1`, the likelihood is
`[p, 1 - p]`.

```@example ldpc_decoding
p = 0.08
transmitted = [1, 0, 1, 1, 1, 0]
received = [1, 0, 1, 0, 1, 0]

channelLikelihood(y, p) = y == 0 ? [1.0 - p, p] : [p, 1.0 - p]

channelFactors = [
    DiscreteFactor(
        bitIds[index],
        channelLikelihood(received[index], p);
        label = "channel_$(bitIds[index])",
        initialize = true
    )
    for index in eachindex(bitIds)
]

nothing # hide
```

---

## Parity-Check Factors

The parity-check matrix defines which bits participate in each even-parity constraint:

```@example ldpc_decoding
H = [
    1 1 0 1 0 0
    0 1 1 0 1 0
    1 0 1 0 0 1
]

function evenParityTable(degree)
    table = zeros(Float64, ntuple(_ -> 2, degree))

    for assignment in Iterators.product(ntuple(_ -> 0:1, degree)...)
        table[map(bit -> bit + 1, assignment)...] = iseven(sum(assignment)) ? 1.0 : 0.0
    end

    return table
end

parityFactors = [
    begin
        checkVariables = bitIds[findall(H[checkIndex, :] .== 1)]

        DiscreteFactor(
            checkVariables...,
            evenParityTable(length(checkVariables));
            label = "check_$checkIndex"
        )
    end
    for checkIndex in axes(H, 1)
]

nothing # hide
```

For each parity factor, table axes follow the order of `checkVariables`.

---

## Factor Graph Construction

Build the factor graph:

```@example ldpc_decoding
graph = factorGraph(variables, vcat(channelFactors, parityFactors))

nothing # hide
```

The graph can be rendered as an SVG factor graph figure:

```@example ldpc_decoding
saveGraphFigure("../ldpcd.svg", graph)

nothing # hide
```

```@raw html
<div class="graph-figure" style="text-align: center;">
  <object
    data="ldpcd.svg"
    type="image/svg+xml"
    aria-label="LDPC decoding factor graph"
    style="width: 55%; height: auto;">
    <a href="ldpcd.svg">LDPC decoding factor graph</a>
  </object>
</div>
```

---

## Running Belief Propagation

Run loopy discrete belief propagation on the graph:

```@example ldpc_decoding
inference = sumproduct(graph)

gbp!(graph, inference; iterations = 25, tolerance = 1e-8, schedule = :flooding)

nothing # hide
```

The posterior marginal of each bit gives the decoder's soft information:

```@example ldpc_decoding
posterior = [
    marginal(graph, inference, bitId)
    for bitId in bitIds
]
```

A hard decision chooses the most likely state for each bit:

```@example ldpc_decoding
decoded = [
    states[argmax(marginal(graph, inference, bitId))]
    for bitId in bitIds
]
```

The decoded word satisfies the parity checks:

```@example ldpc_decoding
mod.(H * decoded, 2)
```
