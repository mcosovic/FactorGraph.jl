# [Graphical Model](@id graphicalModel)

The GaussBP supports the composite type `GraphicalModel` with three fields:
- `FactorGraph`;
- `Inference`;
- `SystemModel`.

The subtype `FactorGraph` describes the factor graph obtained based on the input data. The GBP inference and marginal values are kept in the subtype `Inference`. The system of the linear equations being solved is preserved in the subtype `SystemModel`. Note that the function `graphicalModel()` returns the main GaussBP composite type `GraphicalModel` with all subtypes.

---

Input arguments of the function `graphicalModel()` describe the graphical model, while the function returns `GraphicalModel` type.  The order of inputs and their appearance is arbitrary, with only DATA input required. Still, for the methodological reasons, the syntax examples follow a certain order.

#### Syntax
```julia-repl
graphicalModel(DATA)
graphicalModel(DATA; VIRTUAL)
graphicalModel(DATA; VIRTUAL, DAMP)

```
```@raw html
&nbsp;
```
#### Description
```julia-repl
graphicalModel(DATA) forms the graphical model using DATA input
graphicalModel(DATA; VIRTUAL) sets virtual factor nodes
graphicalModel(DATA; VIRTUAL, DAMP) sets damping parameters
```
```@raw html
&nbsp;
```
#### Output
```julia-repl
gbp = graphicalModel() returns the graphical model
```

---

#### Variable arguments
The function `graphicalModel()` receives the variable argument DATA.


| DATA            | input system data                                                           |
|:----------------|:----------------------------------------------------------------------------|
|                 |                                                                             |
| **Example**     | `"data33_14.h5"`                                                            |
| **Description** | loads the system data using h5-file from the package                        |
|                 |                                                                             |
| **Example**     | `"data33_14.xlsx"`                                                          |
| **Description** |  loads the system data using xlsx-file from the package                     |
|                 |                                                                             |
| **Example**     | `"C:/name.h5"`                                                              |
| **Description** | loads the system data using h5-file from a custom path                      |
|                 |                                                                             |
| **Example**     | `"C:/name.xlsx"`                                                            |
| **Description** | loads the system data using xlsx-file from a custom path                    |
|                 |                                                                             |
| **Example**     | `jacobian, observation, variances`                                          |
| **Description** | loads the system data passing arguments directly                            |


---

#### Keyword arguments
The GBP function `graphicalModel()` receives a group of arguments by keyword: VIRTUAL and DAMP.

| VIRTUAL         | sets virtual factor nodes                                                           |
|:----------------|:------------------------------------------------------------------------------------|
|                 |                                                                                     |
| **Command**     | `mean = value`                                                                      |
| **Description** |  the mean value of the virtual factor nodes, `default setting: mean = 0.0`          |
|                 |                                                                                     |
| **Command**     | `variance = value`                                                                  |
| **Description** |  the variance value of the virtual factor nodes, `default setting: variance = 1e10` |

We advise the reader to read the section [message passing schedule] (@ref schedule) which provides a detailed description of the virtual factor nodes.

```@raw html
&nbsp;
```

| DAMP            | sets damping parameters                                                     |
|:----------------|:----------------------------------------------------------------------------|
| **Command**     | `prob = value`                                                              |
| **Description** |  a Bernoulli random variable with probability `prob = value` independently sampled for each mean value message from a factor node to a variable node, applied for randomised damping iteration scheme with `value` between 0 and 1, `default setting: prob = 0.6`     |
|                 |                                                                             |
| **Command**     | `alpha = value`                                                             |
| **Description** |  the damped message is evaluated as a linear combination of the message from the previous and the current iteration, with weights `alpha = value` and `1 - alpha`, applied for randomised damping iteration scheme where `alpha` is between 0 and 1, `default setting: alpha = 0.4` |

We advise the reader to read the section [the GBP with randomized damping] (@ref dampGBP) which provides a detailed description of the input parameters.











