# [Graphical Model](@id graphicalModel)

The GaussBP supports the composite type `GraphicalModel` with three fields:
- `FactorGraph`;
- `Inference`;
- `SystemModel`.

The subtype `FactorGraph` describes the factor graph obtained based on the input data. The GBP inference and marginal values are kept in the subtype `Inference`. The system of the linear equations being solved is preserved in the subtype `SystemModel`. Note that the function `graphicalModel()` returns the main GaussBP composite type `GraphicalModel` with all subtypes.

In addition, we also provide several functions for factor graph manipulation.

---

#### Build graphical model

Input arguments DATA of the function `graphicalModel()` describe the graphical model, while the function returns `GraphicalModel` type.

Loads the system data using h5-file from the package:
```julia-repl
gbp = graphicalModel("data33_14.h5")
```

Loads the system data using xlsx-file from the package:
```julia-repl
gbp = graphicalModel("data33_14.xlsx")
```

Loads the system data from a custom path:
```julia-repl
gbp = graphicalModel("C:/name.h5")
```

Loads the system data passing arguments directly:
```julia-repl
gbp = graphicalModel(jacobian, observation, variances)
```

---

#### Virtual factor nodes

The GBP function `graphicalModel()` receives arguments by keyword to set the mean and variance of the virtual factor nodes. We advise the reader to read the section [message passing schedule] (@ref schedule) which provides a detailed description of the virtual factor nodes.

```julia-repl
gbp = graphicalModel(DATA; mean = value, variance = value)
```
Default setting of the mean value is `mean = 0.0`, while the default variance is equal to `variance = 1e10`.

---

#### Randomized damping parametars

The GBP function `graphicalModel()` receives arguments by keyword to set damping parametars. We advise the reader to read the section [the GBP with randomized damping] (@ref dampGBP) which provides a detailed description of the input parameters.
```julia-repl
gbp = graphicalModel(DATA; prob = value, alpha = value)
```
The keyword `prob` represents the probability of the Bernoulli random variable, independently sampled for each mean value message from a factor node to a variable node, applied for randomised damping iteration scheme with `value` between 0 and 1. Default setting is set to `prob = 0.6`. The damped message is evaluated as a linear combination of the message from the previous and the current iteration, with weights `alpha = value` and `1 - alpha`, applied for randomised damping iteration scheme where `alpha` is between 0 and 1. Default setting is set to `alpha = 0.4`.

Using the function `graphicalModel()`, the set of damp messages are fixed through GBP iterations. However, we provide the function that changes damp parameters `prob` and `alpha` on the fly:
```julia-repl
damping!(gbp; prob = value, alpha = value)
```

---

#### Freeze factor node, variable node or edge
The functions freeze the target factor or variable node, whereby all messages sent by the factor or variable node retain the values that were at the time of freezing.
```julia-repl
freezeFactor!(gbp; factor = value)
```
```julia-repl
freezeVariable!(gbp; variable = value)
```

Additionally, we provide functions that freeze the target edge. More precisely, the function freezes the message from variable node to factor node or the message from factor node to variable node. Hence, the frozen message retains the value that was at the time of freezing.
```julia-repl
freezeVariableFactor!(gbp; variable = value, factor = value)
```
```julia-repl
freezeFactorVariable!(gbp; factor = value, variable = value)
```

The functions accept the composite type `GraphicalModel` and the factor node index corresponding to the row number of the jacobian matrix, while the variable node index corresponding to the column number of the jacobian matrix.


---

#### Defreeze factor node, variable node or edge
The functions refreeze the target frozen factor node or frozen variable node, whereby the factor or variable node begins to calculate outgoing messages.
```julia-repl
defreezeFactor!(gbp; factor = value)
```
```julia-repl
defreezeVariable!(gbp; variable = value)
```

Also, we provide functions that refreeze the target edge, whereby the message from variable node to factor node or the message from factor node to variable node begins to calculate.
```julia-repl
defreezeVariableFactor!(gbp; variable = value, factor = value)
```
```julia-repl
defreezeFactorVariable!(gbp; factor = value, variable = value)
```

The functions accept the composite type `GraphicalModel` and the factor node index corresponding to the row number of the jacobian matrix, while the variable node index corresponding to the column number of the jacobian matrix.