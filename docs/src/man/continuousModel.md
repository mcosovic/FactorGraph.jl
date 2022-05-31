# [Graphical Model](@id graphicalModelContinuous)

The FactorGraph supports the composite type `ContinuousModel` related with the [synchronous message passing] (@ref synchronousSchedule), with three fields:
- `ContinuousGraph`;
- `ContinuousInference`;
- `ContinuousSystem`.

The subtype `ContinuousGraph` describes the factor graph obtained based on the input data. The GBP inference and marginal values are kept in the subtype `ContinuousInference`. The system of the linear equations being solved is preserved in the subtype `ContinuousSystem`. Note that the function `continuousModel()` returns the main FactorGraph composite type `ContinuousModel` with all subtypes.

In addition, we also provide several functions for factor graph manipulation.

---

#### Build graphical model

Input arguments of the function `continuousModel()` describe the graphical model, while the function returns `ContinuousModel` type.

Loads the system data passing arguments:
```julia-repl
gbp = continuousModel(jacobian, observation, variances)
```

---

#### Virtual factor nodes

The function `continuousModel()` receives arguments by keyword to set the mean and variance of the virtual factor nodes. We advise the reader to read the Section [initialisation procedure] (@ref initialisationProcedure) which provides a detailed description of the virtual factor nodes.

```julia-repl
gbp = continuousModel(DATA; mean = value, variance = value)
```
Default setting of the mean value is `mean = 0.0`, while the default variance is equal to `variance = 1e10`.

---

#### Randomized damping parametars

The function `continuousModel()` receives arguments by keyword to set damping parametars. We advise the reader to read the section [the GBP with randomized damping] (@ref dampGBP) which provides a detailed description of the input parameters.
```julia-repl
gbp = continuousModel(DATA; prob = value, alpha = value)
```
The keyword `prob` represents the probability of the Bernoulli random variable, independently sampled for each mean value message from a factor node to a variable node, applied for randomised damping iteration scheme with `value` between 0 and 1. Default setting is set to `prob = 0.6`. The damped message is evaluated as a linear combination of the message from the previous and the current iteration, with weights `alpha = value` and `1 - alpha`, applied for randomised damping iteration scheme where `alpha` is between 0 and 1. Default setting is set to `alpha = 0.4`.

Using the function `continuousModel()`, the set of damp messages are fixed through GBP iterations. However, we provide the function that changes damp parameters `prob` and `alpha` on the fly:
```julia-repl
damping!(gbp; prob = value, alpha = value)
```

---

#### Freeze factor node, variable node or edge
The functions freeze target factor or variable node, whereby all messages sent by the factor or variable node retain latest obtained values at the time of freezing.
```julia-repl
freezeFactor!(gbp; factor = index)
```
```julia-repl
freezeVariable!(gbp; variable = index)
```

We provide functions that freeze the target edge. More precisely, the function freezes the message from variable node to factor node, or the message from factor node to variable node. Hence, the frozen message keeps the last value obtained at the time of freezing.
```julia-repl
freezeVariableFactor!(gbp; variable = index, factor = index)
```
```julia-repl
freezeFactorVariable!(gbp; factor = index, variable = index)
```
The functions accept following parameters: composite type `ContinuousModel`; the factor node index corresponding to the row index of the Jacobian matrix; and the variable node index corresponding to the column index of the Jacobian matrix. Note that the singly connected factor nodes can not be frozen because they always send the same message.

---

#### Defreeze factor node, variable node or edge
The functions unfreeze the target frozen factor node or frozen variable node, allowing the factor or variable node to calculate outgoing messages.
```julia-repl
defreezeFactor!(gbp; factor = index)
```
```julia-repl
defreezeVariable!(gbp; variable = index)
```

Also, we provide functions that unfreeze the target edge, starting the calculation of messages either from variable node to factor node or from factor node to variable node calculates.
```julia-repl
defreezeVariableFactor!(gbp; variable = index, factor = index)
```
```julia-repl
defreezeFactorVariable!(gbp; factor = index, variable = index)
```

The functions accept following parameters: composite type `ContinuousModel`; the factor node index corresponding to the row index of the Jacobian matrix; and the variable node index corresponding to the column index of the Jacobian matrix. Since singly connected factors cannot be frozen, they cannot be unfreezed.

---

#### Hide factor node
Utilising a hiding mechanism, the function softly deletes factor node. Hence, the function obliterates the target factor node from the graph during the calculation. Soft delete actually removes the node, while preserving node numbering and keeping the same dimensions of the internal variables.
```julia-repl
hideFactor!(gbp; factor = index)
```
If the function targets the singly connected factor node, the function obliterates the target factor only if there are two or more singly connected factor nodes at the same variable node. If there is only one singly connected factor node at the variable node, the function transforms the target factor node to the virtual factor node. Note that to maintain consistency, the function also affects `ContinuousSystem.observation`, `ContinuousSystem.jacobian` and `ContinuousSystem.jacobianTranspose` fields by setting non-zero elements to zero.

---

#### Add factor nodes
The function adds new factor nodes to the existing factor graph.
```julia-repl
addFactors!(gbp; mean = vector, variance = vector, jacobian = matrix)
```
The function supports addition of the multiple factor nodes to initial (existing) formation of the factor graph using the same input data format. The function accepts the following parameters: composite type `ContinuousModel`; the `mean` and `variance` vectors representing new measurement values and variances, respectively. The keyword `jacobian` with corresponding coefficients defines the set of equations describing new factor nodes. Also, function initializes messages from variable nodes to a new factor node using results from the last GBP iteration. Note that the function also affects `ContinuousSystem.observation`, `ContinuousSystem.variance`, `ContinuousSystem.jacobian` and `ContinuousSystem.jacobianTranspose` fields.

