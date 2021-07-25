# [Input Data](@id inputdata)

The package GaussBP supports HDF5 and XLSX input files, or passing data directly via command line arguments. The basic input data structure describing a linear system of equations includes the `jacobian` matrix that contains coefficients of the equations, while vectors `observation` and `variance` represent measurement values and measurement variances, respectively. Note that, with large-scale systems, we strongly recommend using the HDF5 data format.  

Running the GBP algorithm in a dynamic framework requires a variable `dynamic`, which defines the dynamic update scheme of the factor nodes.

---

#### HDF5
The HDF5 input file must contain the following elements:
- coefficient data `model.h5/jacobian`; 
- measurement values `model.h5/observation`;
- measurement variances `model.h5/variance`. 


The input data for the dynamic GBP algorithm must contain the additional variable:
- dynamic data `model.h5/dynamic`.


The type and structure of the input data must be:
- `jacobian::Array{Float64, 2} = [row column coefficient]`;
- `observation::Array{Float64, 1}`;
- `variance::Array{Float64, 1}`;
- `dynamic::Array{Float64, 2} = [iteration factor observation variance]`.

---

#### XLSX
The XLSX input file must contain the following sheets:
- coefficient data `sheet: jacobian`; 
- measurement data `sheet: measurement`.

The input data for the dynamic GBP algorithm must contain the additional sheet:
- dynamic data `sheet: dynamic`.

The type and structure of the input data must be:
- `jacobian - row | column | coefficient`;
- `measurement -  observation | variance`;
- `dynamic - iteration | factor | observation | variance`.
---

#### Passing arguments
The structure of the arguments should be:
- `gbp(jacobian, observation, variance)`

The structure for the dynamic GBP algorithm must contain the additional argument:
- `gbp(jacobian, observation, variance, dynamic)`

The type and structure of the arguments must be:
  - `jacobian::Union{Array{Float64, 2}, SparseMatrixCSC{Float64, Int64}}`
  - `observation::Array{Float64, 1}`
  - `variance::Array{Float64, 1}` 
  - `dynamic::Array{Float64, 2}`
---

#### Data Structure
The `jacobian` input data is used for all analysis available in the GaussBP package and contains coefficients of the linear system of the equations.   

| Column   | Description                                               |
|:--------:|:----------------------------------------------------------|   
| 1        | row indices of the corresponding jacobian matrix          |
| 2        | column indices of the corresponding jacobian matrix       |
| 3        | coefficient values of the corresponding jacobian matrix   | 

```@raw html
&nbsp;
```

The `dynamic` input data is used only for the dynamic GBP algorithm and contains factor nodes update scheme. 

| Column   | Description                                                                            |
|:--------:|:---------------------------------------------------------------------------------------|    
| 1        | the iteration number when the factor node receives new observation and variance value  |
| 2        | factor node index corresponding to the row number of the jacobian matrix               |
| 3        | new observation value                                                                  |  
| 4        | new variance value                                                                     | 

---


#### Use Cases
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
| dataDynamic33_14.xlsx       | 14            | 33        |
| dataDynamic33_14.h5         | 14            | 33        