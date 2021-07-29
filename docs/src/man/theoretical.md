# [Theoretical Background](@id theoretical)

As an input, we observe a noisy linear system of equations with real coefficients and variables:
```math
        \mathbf{z}=\mathbf{h}(\mathbf{x})+\mathbf{u},
```
where ``\mathbf {x}=[x_1,\dots,x_{n}]^{{T}}`` is the vector of the state variables, ``\mathbf{h}(\mathbf{x})= [h_1(\mathbf{x})``, ``\dots``, ``h_k(\mathbf{x})]^{{T}}`` is the vector of observation or measurement functions,  ``\mathbf{z} = [z_1,\dots,z_m]^{{T}}`` is the vector of measurement values, and ``\mathbf{u} = [u_1,\dots,u_k]^{{T}}`` is the vector of uncorrelated measurement errors. The linear system of equations is an overdetermined ``m>n`` arising in many technical fields, such as statistics, signal processing, and control theory. 

Each observation is associated with measured value ``z_i``, measurement error  ``u_i``, and measurement function ``h_i(\mathbf{x})``. Under the assumption that measurement errors ``u_i`` follow a zero-mean Gaussian distribution, the probability density function associated with the ``i``-th measurement is proportional to:
```math
    \mathcal{N}(z_i|\mathbf{x},v_i) \propto \exp\Bigg\{-\cfrac{[z_i-h_i(\mathbf{x})]^2}{2v_i}\Bigg\},
```
where ``v_i`` is the measurement variance defined by the measurement error ``u_i``, and the measurement function ``h_i(\mathbf{x})`` connects the vector of state variables ``\mathbf{x}`` to the value of the ``i``-th measurement.

The goal is to determine state variables ``\mathbf{x}`` according to the noisy observed data ``\mathbf{z}`` and a prior knowledge: 
```math
    p(\mathbf{x}|\mathbf{z})= \cfrac{p(\mathbf{z}|\mathbf{x})p(\mathbf{x})}{p(\mathbf{z})}.
```
Assuming that the prior probability distribution ``p(\mathbf{x})`` is uniform, and given that ``p(\mathbf{z})`` does not depend on ``\mathbf{x}``, the maximum a posteriori solution reduces to the maximum likelihood solution, as given below:
```math
    \hat{\mathbf{x}}= \mathrm{arg}\max_{\mathbf{x}}p(\mathbf{x}|\mathbf{z})= \mathrm{arg}\max_{\mathbf{x}}p(\mathbf{z}|\mathbf{x})=
    \mathrm{arg}\max_{\mathbf{x}}\mathcal{L}(\mathbf{z}|\mathbf{x}).
```

One can find the solution via maximization of the likelihood function ``\mathcal{L}(\mathbf{z}|\mathbf{x})``, which is defined via likelihoods of ``m`` independent measurements:  
```math
    \hat{\mathbf x}= \mathrm{arg} \max_{\mathbf{x}}\mathcal{L}(\mathbf{z}|\mathbf{x})= 
    \mathrm{arg} \max_{\mathbf{x}} \prod_{i=1}^m \mathcal{N}(z_i|\mathbf{x},v_i).
```

It can be shown that the maximum a posteriori solution can be obtained by solving the following optimization problem, known as the weighted least-squares (WLS) problem:
```math
    \hat{\mathbf x} = \mathrm{arg}\min_{\mathbf{x}} \sum_{i=1}^m  \cfrac{[z_i-h_i(\mathbf x)]^2}{v_i}.
```
The state estimate ``\hat{\mathbf x}`` representing the solution of the optimization problem is known as the WLS estimator. The maximum likelihood and WLS estimator are equivalent to the maximum a posteriori solution.

---

## [Linear GBP Algorithm] (@id vanillaGBP)

In the standard setup, the goal of the belief propagation (BP) algorithm is to efficiently evaluate the marginals of a set of random variables ``\mathcal{X} = \{x_1,\dots,x_n\}`` described via the joint probability density function ``g(\mathcal{X})``. Assuming the function ``g(\mathcal{X})`` can be factorised proportionally (``\propto``) to a product of local functions:
```math
    g(\mathcal{X}) \propto \prod_{i=1}^m \psi(\mathcal{X}_i),
```
where ``\mathcal{X}_i \subseteq \mathcal{X}``. The first step is forming a factor graph, which is a bipartite graph that describes the structure of the factorisation. Factor graph allows a graph-based representation of probability density functions using variable and factor nodes connected by edges. In contrast to directed and undirected graphical models, factor graphs provide the details of the factorisation more explicitly. The factor graph structure comprises the set of factor nodes ``\mathcal{F}=\{f_1,\dots,f_m\}``, where each factor node  ``f_i`` represents local function ``\psi(\mathcal{X}_i)``, and the set of variable nodes ``\mathcal{X}``. The factor node ``f_i`` connects to the variable node ``x_s`` if and only if ``x_s \in \mathcal{X}_i``. 

The BP algorithm on factor graphs proceeds by passing two types of messages along the edges of the factor graph: 
- a variable node ``x_s \in \mathcal{X}`` to a factor node ``f_i \in \mathcal{F}`` message ``\mu_{x_s \to f_i}(x_s)``, and  
- a factor node ``f_i \in \mathcal{F}`` to a variable node ``x_s \in \mathcal{X}`` message ``\mu_{f_i \to x_s}(x_s)``.
Both variable and factor nodes in a factor graph process the incoming messages and calculate outgoing messages, where an output message on any edge depends on incoming messages from all other edges. The BP messages represent ``beliefs" about variable nodes, thus a message that arrives or departs a certain variable node is a function (distribution) of the random variable corresponding to the variable node.

The Gaussian belief propagation (GBP) represents a class of the BP, where local function ``\psi(\mathcal{X}_i)`` is defined as a continuous Gaussian distribution:
```math
    \mathcal{N}(z_i|\mathcal{X}_i,v_i) \propto \exp\Bigg\{-\cfrac{[z_i-h(\mathcal{X}_i)]^2}{2v_i}\Bigg\},
```
where ``v_i`` is the variance, and the function ``h(\mathcal{X}_i)`` connects the set of state variables ``\mathcal{X}_i`` to the known ``z_i`` value. The \emph{linear}-GBP model implies the linear function ``h(\mathcal{X}_i)``. If the linear-GBP algorithm converges, it will converge to a fixed point representing a true means \cite{bickson}, regardless of the structure of the factor graph. Unlike means, the variances of the linear-GBP algorithm may not converge to correct values for graphical models with loops, while for models without loops (i.e., tree factor graph) variances will have exact values.

Under the **native GBP algorithm** , we imply the algorithm in which messages are calculated as described below.

#### Message from a variable node to a factor node
Consider a part of a factor graph with a group of factor nodes ``\mathcal{F}_s=\{f_i,f_w,...,f_W\}`` ``\subseteq`` ``\mathcal{F}`` that are neighbours of the variable node ``x_s \in \mathcal{X}``. The message ``\mu_{x_s \to f_i}(x_s)`` from the variable node ``x_s`` to the factor node ``f_i`` is equal to the product of all incoming factor node to variable node messages arriving at all the other incident edges: 
```math
    \mu_{x_s \to f_i}(x_s) =\prod_{f_a \in \mathcal{F}_s \setminus f_i} \mu_{f_a \to x_s}(x_s),
```
where ``\mathcal{F}_s \setminus f_i`` represents the set of factor nodes incident to the variable node ``x_s``, excluding the factor node ``f_i``. Note that each message is a function of the variable ``x_s``.

Let us assume that the incoming messages ``\mu_{f_w \to x_s}(x_s)``, ``\dots``, ``\mu_{f_W \to x_s}(x_s)`` into the variable node ``x_s`` are Gaussian and represented by their mean-variance pairs ``(z_{f_w \to x_s},v_{f_w \to x_s})``, ``\dots``, ``(z_{f_W \to x_s},v_{f_W \to x_s})``. Note that these messages carry beliefs about the variable node ``x_s`` provided by its neighbouring factor nodes ``\mathcal{F}_s\setminus f_i``. It can be shown that the message ``\mu_{x_s \to f_i}(x_s)`` from the variable node ``x_s`` to the factor node ``f_i`` is proportional to:  
```math
    \mu_{x_s \to f_i}(x_s) \propto \mathcal{N}(x_s|z_{x_s \to f_i}, v_{x_s \to f_i}),       
```
with mean ``z_{x_s \to f_i}`` and variance ``v_{x_s \to f_i}`` obtained as: 
```math
    z_{x_s \to f_i} = \Bigg( \sum_{f_a \in \mathcal{F}_s\setminus f_i} \cfrac{z_{f_a \to x_s}}{v_{f_a \to x_s}}\Bigg) v_{x_s \to f_i} \\
    \cfrac{1}{v_{x_s \to f_i}} = \sum_{f_a \in \mathcal{F}_s\setminus f_i} \cfrac{1}{v_{f_a \to x_s}}.
```

After the variable node ``x_s`` receives the messages from all of the neighbouring factor nodes from the set ``\mathcal{F}_s\setminus f_i``, it evaluates the message ``\mu_{x_s \to f_i}(x_s)``, and sends it to the factor node ``f_i``. 

#### Message from a factor node to a variable node
Consider a part of a factor graph that consists of a group of variable nodes ``\mathcal{X}_i = \{x_s, x_l,...,x_L\}`` ``\subseteq`` ``\mathcal X`` that are neighbours of the factor node ``f_i`` ``\in`` ``\mathcal{F}``. The message ``\mu_{f_i \to x_s}(x_s)`` from the factor node ``f_i`` to the variable node ``x_s`` is defined as a product of all incoming variable node to factor node messages arriving at other incident edges, multiplied by the function ``\psi_i(\mathcal{X}_i)`` associated to the factor node ``f_i``, and marginalised over all of the variables associated with the incoming messages:
```math
    \mu_{f_i \to x_s}(x_s)= \int\limits_{x_l}\dots\int\limits_{x_L} \psi_i(\mathcal{X}_i)
    \prod_{x_b \in \mathcal{X}_i\setminus x_s} \big[\mu_{x_b \to f_i}(x_b) \cdot \mathrm{d}x_b\big], 
``` 
where ``\mathcal{X}_i\setminus x_s`` is the set of variable nodes incident to the factor node ``f_i``, excluding the variable node ``x_s``.

Due to linearity of measurement functions ``h_i(\mathcal{X}_i)``, closed form expressions for these messages is easy to obtain and follow a Gaussian form:
```math
    \mu_{f_i \to x_s}(x_s) \propto \mathcal{N}(x_s|z_{f_i \to x_s},v_{f_i \to x_s}).
```
The message ``\mu_{f_i \to x_s}(x_s)`` can be computed only when all other incoming messages (variable to factor node messages) are known. Let us assume that the messages into factor nodes are Gaussian, denoted by: 
```math
        \mu_{x_l \to f_i}(x_l) \propto \mathcal{N}(x_l|z_{x_l \to f_i}, v_{x_l \to f_i})\\
        \vdots\\
        \mu_{x_L \to f_i}(x_L) \propto \mathcal{N}(x_L|z_{x_L \to f_i}, v_{x_L \to f_i}).
``` 
The Gaussian function associated with the factor node ``f_i`` is:
```math
    \mathcal{N}(z_i|\mathcal{X}_i, v_i) \propto \exp\Bigg\{-\cfrac{[z_i-h_i(\mathcal{X}_i)]^2} {2v_i}\Bigg\}.
```
The model contains only linear functions which we represent in a general form as:
```math
    h_i(\mathcal{X}_i) = C_{x_s} x_s + \sum_{x_b \in \mathcal{X}_i\setminus x_s} C_{x_b} x_b,
```
where ``\mathcal{X}_i\setminus x_s`` is the set of variable nodes incident to the factor node ``f_i``, excluding the variable node ``x_s``. 

It can be shown that the message ``\mu_{f_i \to x_s}(x_s)`` from the factor node ``f_i`` to the variable node ``x_s`` is represented by the Gaussian function \eqref{BP_Gauss_fv}, with mean ``z_{f_i \to x_s}`` and variance ``v_{f_i \to x_s}`` obtained as: 
```math
        z_{f_i \to x_s} = \cfrac{1}{C_{x_s}} \Bigg(z_i - \sum_{x_b \in \mathcal{X}_i \setminus x_s} 
        C_{x_b} z_{x_b \to f_i} \Bigg)\\
        v_{f_i \to x_s} = \cfrac{1}{C_{x_s}^2} \Bigg( v_i + \sum_{x_b \in \mathcal{X}_i \setminus x_s} C_{x_b}^2 v_{x_b \to f_i}  \Bigg).
```

To summarise, after the factor node ``f_i`` receives the messages from all of the neighbouring variable nodes from the set ``\mathcal{X}_i\setminus x_s``, it evaluates the message ``\mu_{f_i \to x_s}(x_s)``, and sends it to the variable node ``x_s``. 

#### Marginal inference
The marginal of the variable node ``x_s`` is obtained as the product of all incoming messages into the variable node ``x_s``:
```math
    p(x_s) =\prod_{f_c \in \mathcal{F}_s} \mu_{f_c \to x_s}(x_s),
```
where ``\mathcal{F}_s`` is the set of factor nodes incident to the variable node ``x_s``. It can be shown that the marginal of the state variable ``x_s`` is represented by: 
```math
    p(x_s) \propto \mathcal{N}(x_s|\hat x_s,v_{x_s}),
```
with the mean value ``\hat x_s`` and variance ``v_{x_s}``:      
```math
    \hat x_s = \Bigg( \sum_{f_c \in \mathcal{F}_s} \cfrac{z_{f_c \to x_s}}{v_{f_c \to x_s}}\Bigg) v_{x_s} \\
    \cfrac{1}{v_{x_s}} = \sum_{f_c \in \mathcal{F}_s} \cfrac{1}{v_{f_c \to x_s}}.
```

Finally, the mean-value ``\hat x_s`` is adopted as the estimated value of the state variable ``x_s``. 

---

## [Computation-efficient GBP Algorithm]  (@id efficientGBP)
We can make a substantial improvement to the vanilla GBP algorithm's complexity by reducing the number of calculations per outgoing messages. We achieve this reduction by summarisation of all incoming messages for each variable and factor node instead of summarising all incoming messages per each outgoing message. This simple trick, allow a single variable or factor node to share these summations across all outgoing messages, hence calculating these summations only once. As a result, each outgoing message involves a constant number of operations improving the worst-case running complexity to ``\mathcal{O}(nm)``. In this framework, we calculate the message from the variable node to the factor node as:
```math
        z_{x_s \to f_i} = \Bigg(\alpha_{x_s} - \cfrac{z_{f_i \to x_s}}{v_{f_i \to x_s}}\Bigg) v_{x_s \to f_i} \\
        \cfrac{1}{v_{x_s \to f_i}} = \beta_{x_s} - \cfrac{1}{v_{f_i \to x_s}},
```
where:
```math
    \alpha_{x_s} = \sum_{f_a \in \mathcal{F}_s} \cfrac{z_{f_a \to x_s}}{v_{f_a \to x_s}};  \quad
    \beta_{x_s} = \sum_{f_a \in \mathcal{F}_s} \cfrac{1}{v_{f_a \to x_s}}. 
```
Likewise, the message from the factor node to the variable node is:
```math
    z_{f_i \to x_s} = \cfrac{1}{C_{x_s}} \left(z_i - \alpha_{f_i} \right) + z_{x_s \to f_i} \\
    v_{f_i \to x_s} = \cfrac{1}{C_{x_s}^2} \left( v_i +  \beta_{f_i}  \right) - v_{x_s \to f_i},
```
where:
```math
    \alpha_{f_i} = \sum_{x_b \in \mathcal{X}_i} C_{x_b} z_{x_b \to f_i};  \quad
    \beta_{f_i} = \sum_{x_b \in \mathcal{X}_i} C_{x_b}^2 v_{x_b \to f_i}. 
```

---

## [The GBP and Kahan–Babuška Algorithm]  (@id kahanGBP)
The major drawback of the computation-efficient GBP algorithm is sensitivity to numerical errors because of the summation of floating-point numbers, due to possible significant differences in the values of incoming means and variances. However, this limitation can be alleviated with a compensated summation algorithm, such as the Kahan summation or the improved Kahan–Babuška algorithm. These algorithms increase the complexity of the operations by a constant factor, which means the time complexity of the worst-case remains unaffected. More precisely, we do summation that exists in the messages as:  
```julia-repl
function kahan(summands, total, epsilon)
    t = total + summands
    if abs(total) >= abs(summands)
        epsilon += (total - t) + summands
    else
        epsilon += (summands - t) + total
    end
    total = t
    
    return total, epsilon
end
```

---

## [The Dynamic GBP Algorithm]  (@id dynamicGBP)
To recall, each factor node is associated with the measurement value ``z_i`` and the measurement variance  ``v_i``. The dynamic framework allows the update of these values in any GBP iteration ``\tau``. This framework is an extension to the real-time model that operates continuously and accepts asynchronous measurements from different measurement subsystems. Such measurements are continuously integrated into the running instances of the GBP algorithm. Hence, the GBP algorithm can update the state estimate vector in a time-continuous process.

Additionally, this framework allows for the artificial addition and removal of factor nodes. Then, the initial factor graph, described with the Jacobian matrix, should include all possible measurements. Measurements that are not active are then taken into account via extremely large values of variances (e.g., ``10^{60}``). Consequently, estimates will have a unique solution according to measurement variances whose values are much smaller than ``10^{60}``.

---

## [The Ageing GBP Algorithm]  (@id ageingGBP)
The ageing framework represents an extension of the dynamic model and establishes a model for measurement arrival processes and for the process of measurement deterioration or ageing over time (or GBP iterations). We integrate these measurements regularly into the running instances of the GBP algorithm. 

Let ``\alpha`` denote iteration number when the factor node ``f_i`` receives the new measurement value ``z_i`` with the predefined variance ``v_{i}``. At first, depending on the network dynamics observed by the group of measurements, we may mark predefined variance ``v_i`` for a certain period of iterations ``\alpha \leq \tau \leq \rho`` as a constant value, which is especially advantageous in networks whose dynamics change slowly. Then, after iteration instant ``\rho``, the ageing model increases variance value over iterations ``v_i(\tau)``. More precisely, we associate the Gaussian distribution ``\mathcal{N}(z_i|\mathcal{X}_i, v_i(\tau))`` to the corresponding factor node ``f_i``, where the variance ``v_i(\tau)`` increases its value starting from the predefined variance ``v_i(\tau) = v_i``. Finally, in practice ageing model requires defining a limit from above ``\bar {v}_i`` of a function ``v_i(\tau)``, instead of allowing variance to take on extremely large values. 

Depending on the measurements arriving dynamic, an adaptive mechanism for increasing the variance over the time ``v_i(t)`` can be derived. The logarithmic growth model represents a promising solution for systems with a high sampling rate of the measurements, where a rapid increase in variance is required: 
```math
    v_i(\tau) = \begin{cases} 
      v_i, & \alpha \leq \tau \leq \rho   \\
      a \, \text{log} \left(\cfrac{\tau - \rho + 1 + b}{1 + b} \right ) + v_i, & \rho \leq \tau \leq \theta \\
      \bar {v}_i, & \tau \geq \theta,
  \end{cases}
```
where ``a`` and ``b`` control the rate of growth.  In contrast, the exponential growth model corresponds to systems with a low sampling rate of the measurements:
```math
    v_i(\tau) = \begin{cases} 
      v_i, & \alpha \leq \tau \leq \rho   \\
      v_i(1+b)^{a(\tau - \rho)}, & \rho \leq \tau \leq \theta \\
      \bar {v}_i, & \tau \geq \theta.
  \end{cases}
```
Finally, the linear growth model can be observed as a compromise between logarithmic and exponential growth models:
```math
    v_i(\tau) = \begin{cases} 
      v_i, & \alpha \leq \tau \leq \rho   \\
      a(\tau-\rho) + v_i, & \rho \leq \tau \leq \theta \\
      \bar {v}_i, & \tau \geq \theta.
  \end{cases}
```
GaussBP uses the above data according to the Table.

| Column   |  Label            | Description                                                                            |
|:--------:|:-----------------:|:---------------------------------------------------------------------------------------|  
| 1        | ``\alpha``        | the iteration number when the factor node receives new observation and variance value  |
| 2        | ``i``             | factor node index corresponding to the row number of the jacobian matrix               |
| 3        | ``z_i``           |  new observation value                                                                 | 
| 4        | ``v_i``           | new variance value                                                                     |    
| 5        |                   | the growth model, where linear = 1, logarithmic = 2, exponential = 3                   |
| 6        | ``\rho``          | the iteration number to which the variance maintains a constant value                  |
| 7        | ``a``             | the parameter that controls the rate of the growth                                     |
| 8        | ``b``             | the parameter that controls the rate of the growth                                     |
| 9        | ``\bar {v}_i``    | the variance upper limit value                                                         |