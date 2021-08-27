# [Input Data](@id inputdata)

The GaussBP package supports HDF5 and XLSX input files or passing data directly via command-line arguments. The basic input data structure describing a linear system of equations includes the `jacobian` matrix containing coefficients of the equations, while vectors `observation` and `variance` represent measurement values and measurement variances, respectively. The function `graphicalModel()` accepts `jacobian`, `observation` and `variance` variables to form appropriate probabilistic graphical model. Note that, with large-scale systems, we strongly recommend using the HDF5 file data format.



---

#### HDF5
The HDF5 input file must contain the following elements:
- coefficient data `model.h5/jacobian`;
- measurement values `model.h5/observation`;
- measurement variances `model.h5/variance`.

The type and structure of the input data must be:
- `jacobian::Array{Float64, 2} = [row column coefficient]`;
- `observation::Array{Float64, 1}`;
- `variance::Array{Float64, 1}`.

---

#### XLSX
The XLSX input file must contain the following sheets:
- coefficient data `sheet: jacobian`;
- measurement data `sheet: measurement`.

The type and structure of the input data must be:
- `jacobian - row | column | coefficient`;
- `measurement -  observation | variance`.

---

#### Passing arguments
The structure of the arguments must be:
- `gbp = graphicalModel(jacobian, observation, variance)`.

The type and structure of the arguments must be:
  - `jacobian::Union{Array{Float64, 2}, SparseMatrixCSC{Float64, Int64}}`;
  - `observation::Array{Float64, 1}`;
  - `variance::Array{Float64, 1}`;
---

#### Data structure
The GaussBP package uses `jacobian` input data format for all analyses. Jacobian input data contains coefficients of the linear system of the equations. The structure of the `jacobian` variable being loaded from HDF5 or XLSX input files is given in the table below, while passing data directly via command-line arguments allows the use of a sparse or full matrix to describe the `jacobian` variable.

| Column   | Description                                               |
|:--------:|:----------------------------------------------------------|
| 1        | row indices of the corresponding jacobian matrix          |
| 2        | column indices of the corresponding jacobian matrix       |
| 3        | coefficient values of the corresponding jacobian matrix   |

```@raw html
&nbsp;
```
The `observation` input data is used for all analyses available in the GaussBP package and contains measurement mean values.

| Column   | Description                                               |
|:--------:|:----------------------------------------------------------|
| 1        | measurement mean value                                    |

```@raw html
&nbsp;
```
The `variance` input data is used for all analyses available in the GaussBP package and contains measurement variance values.

| Column   | Description                                               |
|:--------:|:----------------------------------------------------------|
| 1        | measurement variance value                                |

---


#### Use cases
The pre-defined data are located in the `src/example` as the `.h5` or `.xlsx` files.

| Case                        | Variables     | Equations |
|:----------------------------|--------------:|----------:|
| data33_14.xlsx              | 14            | 33        |
| data33_14.h5                | 33            | 14        |
| data897_300.h5              | 300           | 897       |
| data3119_1354.h5            | 1354          | 3119      |
| data5997_2000.h5            | 2000          | 5997      |
| data7149_2000.h5            | 2000          | 7149      |
| data29997_10000.h5          | 10000         | 29997     |
| data283803_70000.h5         | 70000         | 283803    |