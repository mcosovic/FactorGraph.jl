# [Inference](@id inferenceContinuous)

We advise the reader to read the Section [continuous Gaussian random variables] (@ref continuousVariables) which provides a detailed description of the inference algorithms. To exchange information over the factor graph, the FactorGraph provides three inference approaches:
 - [vanilla GBP algorithm] (@ref vanillaGBP);
 - [broadcast GBP algorithm] (@ref broadcastGBP);
 - [broadcast GBP with Kahan–Babuška algorithm] (@ref kahanGBP).

Each of the inference functions accepts only the composite type `ContinuousModel`, i.e., an output variable of the function `gbp = continuousModel()` and applies the [synchronous message passing schedule] (@ref synchronousSchedule).

---

#### Message inference
The set of functions that can be used to preform message inference:
```julia-repl
messageFactorVariable(gbp); messageVariableFactor(gbp)
messageFactorVariableBroadcast(gbp); messageVariableFactorBroadcast(gbp)
messageFactorVariableKahan(gbp); messageVariableFactorKahan(gbp)
```
---

#### Mean inference
The set of functions that can be used to preform only mean inference:
```julia-repl
meanFactorVariable(gbp); meanVariableFactor(gbp)
meanFactorVariableBroadcast(gbp); meanVariableFactorBroadcast(gbp)
meanFactorVariableKahan(gbp); meanVariableFactorKahan(gbp)
```
---

#### Variance inference
The set of functions that can be used to preform only variance inference:
```julia-repl
varianceFactorVariable(gbp); varianceVariableFactor(gbp)
varianceFactorVariableBroadcast(gbp); varianceVariableFactorBroadcast(gbp)
varianceFactorVariableKahan(gbp); varianceVariableFactorKahan(gbp)
```
---

#### Randomised damping inference
Additionaly, we provide the set of functions to preform damping inference:
```julia-repl
messageDampFactorVariable(gbp); meanDampFactorVariable(gbp)
messageDampFactorVariableBroadcast(gbp); meanDampFactorVariableBroadcast(gbp)
messageDampFactorVariableKahan(gbp); meanDampFactorVariableKahan(gbp)
```
---

#### Marginal inference
To compute marginals the FactorGraph provides the function:
```julia-repl
marginal(gbp)
```
Same as before, the function accepts only the composite type `ContinuousModel`.

---

#### Dynamic inference
This framework is an extension to the real-time model that operates continuously and accepts asynchronous observation and variance values. More precisely, in each GBP iteration user can change the observation and variance values of the corresponding factor nodes and continue the GBP iteration process. We advise the reader to read the section [dynamic GBP algorithm] (@ref dynamicGBP) which provides a detailed description of the input parameters.
```julia-repl
dynamicFactor!(gbp; factor = index, observation = value, variance = value)
```
The function accepts the composite type `ContinuousModel` and keywords `factor`, `observation` and `variance`, which defines the dynamic update scheme of the factor nodes. The factor node index corresponding to the row index of the coefficient matrix. Note that during each function call, `ContinuousSystem.observation` and `ContinuousSystem.variance` fields also change values according to the scheme.

---

#### Dynamic inference with variance ageing
The ageing framework represents an extension of the dynamic model and establishes a model for data arrival processes and for the process of data deterioration or ageing over time (or GBP iterations). We integrate these data regularly into the running instances of the GBP algorithm. We advise the reader to read the section [ageing GBP algorithm] (@ref ageingGBP) which provides a detailed description of the input parameters.
```julia-repl
ageingVariance!(gbp; factor = index, initial = value, limit = value,
                model = value, a = value, b = value, tau = value)
```
This function should be integrated into the iteration loop to ensure variance ageing over iterations. The function accepts the composite type `ContinuousModel` and the keywords `factor`, `initial`, `limit`, `model`, `a`, `b` and `tau`. The variance growth model can be linear `model = 1`, logarithmic `model = 2` and exponential `model = 3`, where parameters `a` and `b` control the rate of the growth. The `initial` defines the initial value of the variance, while the variance upper limit value is defined according to `limit`. The ageing model increases the value of variance over iterations, thus the current iteration step should be forwarded using `tau` keyword. Also, during each function call, `ContinuousSystem.variance` field changes values according to the ageing model.