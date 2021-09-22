# [Input Data](@id inputdataDiscrete)

The FactorGraph package requires knowledge about the joint probability density function ``g(\mathcal{X})`` of the set of random variables ``\mathcal{X} = \{x_1,\dots,x_n\}`` that can be factorised proportionally  (``\propto``) to a product of local functions:
```math
    g(\mathcal{X}) \propto \prod_{i=1}^m \psi_i(\mathcal{X}_i).
```
The FactorGraph package supports disrecte random variables, where each random variable having ``K`` possible states, while local function ``\psi_i(\mathcal{X}_i)`` is defined as the conditional probability distribution:
```math
  p_i(x_i|\mathcal{X}_i \setminus x_i),
```
Hence, the local function is associated with the conditional distribution ``p_i(x_i|\mathcal{X}_i \setminus x_i)`` and the set of random variables ``\mathcal{X}_i``. The conditional distribution takes probability values over all possible state combinations of the random variables from the set ``\mathcal{X}_i`` forming the conditional probability table. To describe the joint probability density function ``g(\mathcal{X})``, it is enough to define the container that saves indices of discrete variables and conditional probability tables.


The input data used to describe the joint probability density function can be given using HDF5 and XLSX input files or by passing data directly via command-line arguments. The examples of the input files are located in the `src/example`. Note that, with large-scale systems, we strongly recommend using the HDF5 file data format.

Thus, the parameters that describe the factor graph structure are represented by variable `probability` containing indices of random variables, while variable `table` saves conditional probability tables. The function `discreteTreeModel()` accepts HDF5 and XLSX input files with variables `probability` and `table`. Also, this function accepts variables `probability` and `table` directly via function arguments.

---


#### Build the graphical model by passing arguments
Let us observe the following joint probability density function:
```math
    g(\mathcal{X})  =  p_1(x_1)p_2(x_1|x_2)p_3(x_1|x_2,x_3),
```
where all variables have two states ``1`` and ``2``. The conditional probability tables can be written in the compact form:

|             |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
| :---------: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: |
| ``x_1``     | 1   | 2   | 1   | 2   | 1   | 2   | 1   | 2   | 1   | 2   | 1   | 2   | 1   | 2   |
| ``x_2``     |     |     | 1   | 1   | 2   | 2   | 1   | 1   | 2   | 2   | 1   | 1   | 2   | 2   |
| ``x_3``     |     |     |     |     |     |     | 1   | 1   | 1   | 1   | 2   | 2   | 2   | 2   |
| probability | 0.6 | 0.4 | 0.8 | 0.5 | 0.2 | 0.5 | 1.0 | 0.1 | 0.5 | 0.6 | 0.0 | 0.9 | 0.5 | 0.4 |
|             |     |     |     |     |     |     |     |     |     |     |     |     |     |     |

Let us define the variable `probability` that contains discrete variable indices:
```julia-repl
probability1 = [1]
probability2 = [1; 2]
probability3 = [1; 2; 3]

probability = [probability1, probability2, probability3]
```

Further, we define the variable `table` that holds conditional probability tables:
```julia-repl
table1 = zeros(2)
table1[1] = 0.6; table1[1] = 0.4

table2 = zeros(2, 2)
table2[1, 1] = 0.8; table2[2, 1] = 0.5; table2[1, 2] = 0.2; table2[2, 2] = 0.5

table3 = zeros(2, 2, 2)
table3[1, 1, 1] = 1.0; table3[2, 1, 1] = 0.1; table3[1, 2, 1] = 0.5; table3[2, 2, 1] = 0.6
table3[1, 1, 2] = 0.0; table3[2, 1, 2] = 0.9; table3[1, 2, 2] = 0.5; table3[2, 2, 2] = 0.4

table = [table1, table2, table3]
```
Passing data directly via command-line arguments supports the above format, where the function `discreteTreeModel()` accepts variables `probability` and `table`:
```julia-repl
bp = discreteTreeModel(probability, table)
```
Here, the variable `bp` holds the main composite type related to the discrete model. We can also pack the input data as a dictionary and pass to the function `discreteTreeModel()`:
```julia-repl
inputData = Dict("probability" => (probability1, probability2, probability3),
                 "table"  => (table1, table2, table3))

bp = discreteTreeModel(inputData)
```

----

#### Build the graphical model using HDF5 file
The HDF5 input data stored above format in the fields `probability` and `table`. Inside these fields, each subfield holds data about specific conditional probability distribution with a unique subfield name. Subfield names through fields `probability` and `table` must match.

---

#### Build the graphical model using XLSX file
The input data stored in the XLSX file requires conditional probability tables in the native form with discrete variable indices. Further, the number of possible states should be given for each random variable.