# [Continuous Gaussian random variables] (@id continuousVariables)
The Gaussian belief propagation (GBP) represents a class of the BP, where local function ``\psi(\mathcal{X}_i)`` is defined as a continuous Gaussian distribution:
```math
    \mathcal{N}(z_i|\mathcal{X}_i,v_i) \propto \exp\Bigg\{-\cfrac{[z_i-h(\mathcal{X}_i)]^2}{2v_i}\Bigg\},
```
where ``v_i`` is the variance, and the function ``h(\mathcal{X}_i)`` connects the set of state variables ``\mathcal{X}_i`` to the known ``z_i`` value. The linear-GBP model implies the linear function ``h(\mathcal{X}_i)``. If the linear-GBP algorithm converges, it will converge to a fixed point representing a true means [1], regardless of the structure of the factor graph. Unlike means, the variances of the linear-GBP algorithm may not converge to correct values for graphical models with loops, while for models without loops (i.e., tree factor graph) variances will have exact values.

Thus, as an input, we observe a noisy linear system of equations with real coefficients and variables:
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

### [Vanilla GBP algorithm] (@id vanillaGBP)
Under the vanilla GBP algorithm, we imply the algorithm in which messages are calculated as described below.

#### Message from a factor node to a variable node
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
It can be shown that the marginal of the state variable ``x_s`` is represented by:
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

### [Broadcast GBP algorithm]  (@id broadcastGBP)
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

### [Broadcast GBP and Kahan–Babuška algorithm]  (@id kahanGBP)
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

### [The GBP with randomized damping]  (@id dampGBP)
We propose a randomized damping approach, where each mean value message from factor node to a variable node is damped independently with probability ``p``, otherwise, the message is calculated as in the standard the GBP algorithm. The damped message is evaluated as a linear combination of the message from the previous and the current iteration, with weights ``\alpha`` and ``1 - \alpha``, respectively. More, precisly, the proposed randomized damping scheduling updates of selected factor to variable
node means in every iteration by combining them with their values from the previous iteration using convergence parameters ``p`` and ``\alpha``:
```math
    z_{f_{i} \rightarrow x_{s}}^{(\tau)}=\left(1-q_{i s}\right) \cdot z_{f_{i} \rightarrow x_{s}}^{(\tau)}+q_{i s} \cdot\left(\alpha \cdot z_{f_{x} \rightarrow x_{s}}^{(\tau-1)}+(1-\alpha) \cdot z_{f_{i} \rightarrow x_{a}}^{(\tau)}\right),
```
where ``q_{i s} \sim \operatorname{Ber}(p) \in\{0,1\}`` is independently sampled with probability ``p`` for the mean from factor node ``f_i`` to the variable node ``x_s``.

The randomised damping parameter pairs lead to a trade-off between the number of non-converging simulations and the rate of convergence. In general, we observe a large number of non-converging simulations for the selection of `prob` and `alpha` for which only a small fraction of messages are combined with their values in a previous iteration, and that is a case for `prob` close to 0 or `alpha` close to 1.

---

### [Dynamic GBP algorithm]  (@id dynamicGBP)
To recall, each factor node is associated with the measurement value ``z_i`` and the measurement variance  ``v_i``. The dynamic framework allows the update of these values in any GBP iteration ``\tau``. This framework is an extension to the real-time model that operates continuously and accepts asynchronous measurements from different measurement subsystems. Such measurements are continuously integrated into the running instances of the GBP algorithm. Hence, the GBP algorithm can update the state estimate vector in a time-continuous process.

Additionally, this framework allows for the artificial addition and removal of factor nodes. Then, the initial factor graph, described with the Jacobian matrix, should include all possible measurements. Measurements that are not active are then taken into account via extremely large values of variances (e.g., ``10^{60}``). Consequently, estimates will have a unique solution according to measurement variances whose values are much smaller than ``10^{60}``.

---

### [Ageing GBP algorithm]  (@id ageingGBP)
The ageing framework represents an extension of the dynamic model and establishes a model for measurement arrival processes and for the process of measurement deterioration or ageing over time (or GBP iterations). We integrate these measurements regularly into the running instances of the GBP algorithm.

Let us assume that factor node ``f_i`` receives the new variance ``v_{i}``. After that moment, the ageing model increases variance value over iterations ``v_i(\tau)``. More precisely, we associate the Gaussian distribution ``\mathcal{N}(z_i|\mathcal{X}_i, v_i(\tau))`` to the corresponding factor node ``f_i``, where the variance ``v_i(\tau)`` increases its value starting from the predefined variance ``v_i(\tau) = v_i``. Finally, in practice ageing model requires defining a limit from above ``\bar {v}_i`` of a function ``v_i(\tau)``, instead of allowing variance to take on extremely large values.

Depending on the measurements arriving dynamic, an adaptive mechanism for increasing the variance over iterations ``v_i(\tau)`` can be derived. The logarithmic growth model represents a promising solution for systems with a high sampling rate of the measurements, where a rapid increase in variance is required:
```math
    v_i(\tau) = \begin{cases}
      a \, \text{log} \left(\cfrac{\tau + 1 + b}{1 + b} \right ) + v_i, & 1 \leq \tau \leq \theta \\
      \bar {v}_i, & \tau \geq \theta,
  \end{cases}
```
where ``a`` and ``b`` control the rate of growth.  In contrast, the exponential growth model corresponds to systems with a low sampling rate of the measurements:
```math
    v_i(\tau) = \begin{cases}
      v_i(1+b)^{a\tau}, & 1 \leq \tau \leq \theta \\
      \bar {v}_i, & \tau \geq \theta.
  \end{cases}
```
Finally, the linear growth model can be observed as a compromise between logarithmic and exponential growth models:
```math
    v_i(\tau) = \begin{cases}
      a\tau + v_i, & 1 \leq \tau \leq \theta \\
      \bar {v}_i, & \tau \geq \theta.
  \end{cases}
```

---

### [References](@id refsBelief)
[1] D. Bickson, *Gaussian Belief Propagation: Theory and Aplication,* ArXiv e-prints, Nov. 2008.


