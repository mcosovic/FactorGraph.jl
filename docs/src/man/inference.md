# [Inference](@id vanilla)

We advise the reader to read the [theoretical bacground] (@ref theoretical) which provides a detailed description of the inference algorithms. To exchange information over the factor graph, the GaussBP provides three inference approaches:
- [vanilla GBP algorithm] (@ref vanillaGBP),
- [computation-efficient GBP algorithm] (@ref efficientGBP),
- [computation-efficient kahan–babuška GBP algorithm] (@ref kahanGBP).

Each of the inference functions accepts only the composite type `GraphicalModel`, i.e., an output variable of the function `gbp = graphicalModel()`.

---

#### Message inference
The set of functions that can be used to preform message inference:
```julia-repl
messageFactorVariableVanilla(gbp); messageVariableFactorVanilla(gbp)
messageFactorVariableEfficient(gbp); messageVariableFactorEfficient(gbp)
messageFactorVariableKahan(gbp); messageVariableFactorKahan(gbp)
```
---

#### Mean inference
The set of functions that can be used to preform only mean inference:
```julia-repl
meanFactorVariableVanilla(gbp); meanVariableFactorVanilla(gbp)
meanFactorVariableEfficient(gbp); meanVariableFactorEfficient(gbp)
meanFactorVariableKahan(gbp); meanVariableFactorKahan(gbp)
```
---

#### Variance inference
The set of functions that can be used to preform only variance inference:
```julia-repl
varianceFactorVariableVanilla(gbp); varianceVariableFactorVanilla(gbp)
varianceFactorVariableEfficient(gbp); varianceVariableFactorEfficient(gbp)
varianceFactorVariableKahan(gbp); varianceVariableFactorKahan(gbp)
```
---


#### Randomised damping inference
Additionaly, we provide the set of functions to preform damping inference:
```julia-repl
messageDampFactorVariableVanilla(gbp); meanDampFactorVariableVanilla(gbp)
messageDampFactorVariableEfficient(gbp); meanDampFactorVariableEfficient(gbp)
messageDampFactorVariableKahan(gbp); meanDampFactorVariableKahan(gbp)
```
---

#### Marginal inference
To compute marginals the GaussBP provides the function:
```julia-repl
marginal(gbp)
```
Same as before, the function accepts only the composite type `GraphicalModel`.

---

#### Dynamic inference
This framework is an extension to the real-time model that operates continuously and accepts asynchronous measurement mean and variance values. More precisely, in each GBP iteration user can change the mean and variance values of the corresponding factor nodes and continue the GBP iteration process. We advise the reader to read the section [dynamic GBP algorithm] (@ref dynamicGBP) which provides a detailed description of the input parameters.
```julia-repl
dynamicFactor!(gbp; factor = index, mean = value, variance = value)
```
The function accepts the composite type `GraphicalModel` and keywords `factor`, `mean` and `variance`, which defines the dynamic update scheme of the factor nodes. The factor node index corresponding to the row index of the jacobian matrix. Note that during each function call, `SystemModel.observation` and `SystemModel.variance` fields also change values according to the scheme.

---

#### Dynamic inference with variance ageing
The ageing framework represents an extension of the dynamic model and establishes a model for measurement arrival processes and for the process of measurement deterioration or ageing over time (or GBP iterations). We integrate these measurements regularly into the running instances of the GBP algorithm. We advise the reader to read the section [ageing GBP algorithm] (@ref ageingGBP) which provides a detailed description of the input parameters.
```julia-repl
ageingVariance!(gbp; factor = index, initial = value, limit = value,
                model = value, a = value, b = value, tau = value)
```
This function should be integrated into the iteration loop to ensure variance ageing over iterations. The function accepts the composite type `GraphicalModel` and the keywords `factor`, `initial`, `limit`, `model`, `a`, `b` and `tau`. The variance growth model can be linear `model = 1`, logarithmic `model = 2` and exponential `model = 3`, where parameters `a` and `b` control the rate of the growth. The `initial` defines the initial value of the variance, while the variance upper limit value is defined according to `limit`. The ageing model increases the value of variance over iterations, thus the current iteration step should be forwarded using `tau` keyword. Also, during each function call, `SystemModel.variance` field changes values according to the ageing model.