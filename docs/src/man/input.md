# [Input Data](@id inputdata)

The package GaussBP supports HDF5 and CSV input files or accept arguments passed directly. The basic input data structure used to describe a linear system of equations includes the matrix `H` that contains coefficients of the equations, while vectors `z` and `v` represent observation or measurement values and observation variances, respectively. Note that, in the case of large-scale systems, we strongly recommend to use the HDF5 files for the input data.  

---

#### HDF5
The HDF5 input file must contain the following elements:
- `model_name.h5/H`, `model_name.h5/z`, `model_name.h5/v` 
  - coefficient data `H::Array{Float64, 2} = [row column coefficient]`
  - observation values `z::Array{Float64, 1}`
  - observation variances `v::Array{Float64, 1}`
---

#### CSV
The structure of the CSV input file must contain:
- `model_name.csv` 
  - data with columns `row | column | coefficient | observation | variance`
---

#### Passing arguments
The structure of the arguments should be:
- `gbp(H, z, v)`
  - coefficient data `H::Union{Array{Float64, 2}, SparseMatrixCSC{Float64, Int64}}`
  - observation values `z::Array{Float64, 1}`
  - observation variances `v::Array{Float64, 1}` 
---

#### Use Cases
The pre-defined data are located in the `src/example` as the `.h5` or `.csv` files.

| Case                | Variables     | Equations | 
|:--------------------|--------------:|----------:|
| data33_14.csv       | 14            | 33        |
| data33_14.h5        | 33            | 14        | 
| data897_300.h5      | 300           | 897       | 
| data3119_1354.h5    | 1354          | 3119      | 
| data5997_2000.csv   | 2000          | 5997      | 
| data5997_2000.h5    | 2000          | 5997      | 
| data7149_2000.h5    | 2000          | 7149      |
| data29997_10000.h5  | 10000         | 29997     |
| data283803_70000.h5 | 70000         | 283803    | 
