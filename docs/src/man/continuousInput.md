# [Input Data](@id inputdataContinuous)

The FactorGraph package requires knowledge about the joint probability density function ``g(\mathcal{X})`` of the set of random variables ``\mathcal{X} = \{x_1,\dots,x_n\}`` that can be factorised proportionally  (``\propto``) to a product of local functions:
```math
    g(\mathcal{X}) \propto \prod_{i=1}^m \psi_i(\mathcal{X}_i).
```
The FactorGraph package supports continuous random variables, where local function ``\psi_i(\mathcal{X}_i)`` is defined as a continuous Gaussian distribution:
```math
  \mathcal{N}(z_i|\mathcal{X}_i,v_i) \propto \exp\Bigg\{-\cfrac{[z_i-h_i(\mathcal{X}_i)]^2}{2v_i}\Bigg\}.
```
Hence, the local function is associated with observation ``z_i``, variance  ``v_i``, and linear equation ``h_i(\mathcal{X}_i)``. To describe the joint probability density function ``g(\mathcal{X})``, it is enough to define the coefficient matrix containing coefficients of the equations, and vectors of observation and variance values. Then, the ``i``-th row of the coefficient matrix, with consistent entries of observation and variance, defines the local function ``\psi_i(\mathcal{X}_i)``.

Thus, the input data structure includes the `coefficient` variable which describes coefficients of the equations, while variables `observation` and `variance` represent observation and variance vectors, respectively. The functions `continuousModel()` and `continuousTreeModel()` accept variables `coefficient`, `observation` and `variance`.

---


#### Build the graphical model
Let us observe the following joint probability density function:
```math
    g(\mathcal{X}) \propto  \exp\Bigg\{-\cfrac{[2.5 - 0.2x_1]^2}{2\cdot 1.1}\Bigg\}\exp\Bigg\{-\cfrac{[0.6 - (2.1x_1 + 3.4x_2)]^2}{2\cdot 3.5}\Bigg\}
```

We can describe the joint probability density function using variables `coefficient`, `observation` and `variance`. Passing data directly via command-line support the following format, where the coefficient matrix can be defined as a full or sparse matrix:
```julia-repl
coefficient = zeros(2, 2)
coefficient[1, 1] = 0.2; coefficient[2, 1] = 2.1; coefficient[2, 2] = 3.4
observation = [2.5; 0.6]
variance = [1.1; 3.5]

gbp = continuousModel(coefficient, observation, variance)
```
Here, the variable `gbp` holds the main composite type related to the continuous model. In the case of a tree factor graph, when you want to use a forward-backward GBP algorithm, then the following command can be used:
```julia-repl
gbp = continuousTreeModel(coefficient, observation, variance)
```