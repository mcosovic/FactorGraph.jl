# Gaussian Belief Propagation Linear Systems Solver
The solver provides the solution of the linear system of equations with/without Gaussian noise using belief propagation (BP) algorithm applied over the factor graph.

## The System Model
We observe a noisy linear system of equations with real coefficients and variables:

![equation](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bb%7D%20%3D%20%5Ctextbf%7Bf%7D%28%5Ctextbf%7Bx%7D%29%20&plus;%20%5Ctextbf%7Bu%7D)

where **x** is the vector of the state variables (i.e., unknowns), **f**(**x**) is the vector of linear functions, **b** is the vector of observation values and **u** is the vector of uncorrelated observation errors. Note that the linear system of equations represents an overdetermined system.

The solution can be obtained by solving linear weighted least-squares (WLS) problem:

![wls](https://latex.codecogs.com/gif.latex?%28%5Ctextbf%7BA%7D%5ET%5Ctextbf%7BW%7D%5Ctextbf%7BA%7D%29%5Ctextbf%7Bx%7D%3D%5Ctextbf%7BA%7D%5ET%5Ctextbf%7BW%7D%5Ctextbf%7Bb%7D)

where **A** is the Jacobian matrix of linear functions or the coefficient  matrix for our system, and **W** is a diagonal matrix containing inverses of observation variances.

Further, the solution to the problem can be found via maximization of the likelihood function which is defined via likelihoods of independent observations, and that can be efficiently solved utilizing factor graphs and the Gaussian belief propagation (BP) algorithm.

## Installation
To install GaussianBP, you can run the following:
```
(v1.2) pkg> add https://github.com/mcosovic/GaussianBP
```

## Input Data
MODEL.h5 file located in `\src\data\` with variables:
- `MODEL.h5/H` - coefficient list of type Array{Float64,2} in the form [row column coefficient];
- `MODEL.h5/b` - observation values of type Array{Float64,1};
- `MODEL.h5/v` - observation variances of type Array{Float64,1};


 ## User Options
1. Design of Iteration Scheme:
   - `DAMP` - applied randomized damping at the BP iteration;
   - `MAXI` - the upper limit on BP iterations;

2. Convergence Parameters:
   - `PROB` - a Bernoulli random variable with probability "PROB" independently sampled for each mean value message from indirect factor node to a variable node, with values between 0 and 1;
   - `ALPH` - the damped message is evaluated as a linear combination of the message from the previous and the current iteration, with weights "ALPH" and 1 - "ALPH", where "ALPH" is between 0 and 1;

Note: We use an improved BP algorithm that applies synchronous scheduling  with randomized damping. The randomized damping parameter pairs lead to a trade-off between the number of non-converging simulations and the rate of convergence. In general, for the selection of "prob" and "alph" for which only a small fraction of messages are combined with their values in a previous iteration, and that is a case for "prob" close to 0 or "alph" close to 1, we observe a large number of non-converging simulations.

3. Virtual Factor Nodes
   - `MEAN` - the mean value of virtual factor nodes;
   - `VARI` - the variance value of the virtual factor nodes;

Note: The virtual factor node is a singly-connected factor node used if the variable node x is not directly observed. In a usual scenario, without prior knowledge, the variance of virtual factor nodes tend to infinity.

4. Post-Processing Options:
  - `TIME = "on"` - shows belief propagation time
  - `error = "on"` - shows belief propagation evaluation versus weighted least-squares  

## Algorithms
1. Belief propagation with simply summing messages
```
julia> using GaussianBP
julia> xbp = bp(MODEL, MAXI, DAMP, PROB, ALPH, MEAN, VARI; TIME, ERROR)
```

2. Improved Kahan-Babuska algorithm for summing variance messages with executive function
```
julia> using GaussianBP
julia> xbp = bpn(MODEL, MAXI, DAMP, PROB, ALPH, MEAN, VARI; TIME, ERROR)
```

## Quick Start
```
julia> using GaussianBP
julia> bp("data33_14", 100; TIME = "on")
julia> bp("data33_14", 100; ERROR = "on")
julia> bpn("data33_14", 100, 10)
julia> bpn("data33_14", 100, 10; TIME = "on")
```

## More information:
- M. Cosovic and D. Vukobratovic, "Distributed Gauss-Newton Method for State Estimation Using Belief Propagation," in IEEE Transactions on  Power Systems, vol. 34, no. 1, pp. 648-658, Jan. 2019. [arxiv.org](https://arxiv.org/pdf/1702.05781.pdf)
- M. Cosovic, "Design and Analysis of Distributed State Estimation Algorithms Based on Belief Propagation and Applications in Smart Grids." arXiv preprint arXiv:1811.08355 (2018). [arxiv.org](https://arxiv.org/pdf/1811.08355.pdf)
