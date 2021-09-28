# [Inference in Factor Graphs] (@id inferenceFactorGraphs)

In the standard setup, we observe the set of random variables ``\mathcal{X} = \{x_1,\dots,x_n\}`` described via the joint probability density function ``g(\mathcal{X})``. Assuming the function ``g(\mathcal{X})`` can be factorised proportionally (``\propto``) to a product of local functions:
```math
    g(\mathcal{X}) \propto \prod_{i=1}^m \psi(\mathcal{X}_i),
```
where ``\mathcal{X}_i \subseteq \mathcal{X}``. The first step is forming a factor graph, which is a bipartite graph that describes the structure of the factorisation. Factor graph allows a graph-based representation of probability density functions using variable and factor nodes connected by edges. In contrast to directed and undirected graphical models, factor graphs provide the details of the factorisation more explicitly. The factor graph structure comprises the set of factor nodes ``\mathcal{F}=\{f_1,\dots,f_m\}``, where each factor node  ``f_i`` represents local function ``\psi(\mathcal{X}_i)``, and the set of variable nodes ``\mathcal{X}``. The factor node ``f_i`` connects to the variable node ``x_s`` if and only if ``x_s \in \mathcal{X}_i``.


The message passing algorithm on factor graphs proceeds by passing two types of messages along the edges of the factor graph:
- a variable node ``x_s \in \mathcal{X}`` to a factor node ``f_i \in \mathcal{F}`` message ``\mu_{x_s \to f_i}(x_s)``, and
- a factor node ``f_i \in \mathcal{F}`` to a variable node ``x_s \in \mathcal{X}`` message ``\mu_{f_i \to x_s}(x_s)``.
Both variable and factor nodes in a factor graph process the incoming messages and calculate outgoing messages, where an output message on any edge depends on incoming messages from all other edges. The messages represent "beliefs" about variable nodes, thus a message that arrives or departs a certain variable node is a function (distribution) of the random variable corresponding to the variable node.

Here we shall focus on the problem of evaluating local marginals over variable nodes ``\mathcal{X} = \{x_1,\dots,x_n\}`` described via the joint probability density function ``g(\mathcal{X})``, which will lead us to the belief propagation (BP) algorithm, also known as the sum-product algorithm [1].

---

### [Message from a variable node to a factor node] (@id variableFactorMessage)
Consider a part of a factor graph with a group of factor nodes ``\mathcal{F}_s=\{f_i,f_w,...,f_W\}`` ``\subseteq`` ``\mathcal{F}`` that are neighbours of the variable node ``x_s \in \mathcal{X}``. The message ``\mu_{x_s \to f_i}(x_s)`` from the variable node ``x_s`` to the factor node ``f_i`` is equal to the product of all incoming factor node to variable node messages arriving at all the other incident edges:
```math
    \mu_{x_s \to f_i}(x_s) =\prod_{f_a \in \mathcal{F}_s \setminus f_i} \mu_{f_a \to x_s}(x_s),
```
where ``\mathcal{F}_s \setminus f_i`` represents the set of factor nodes incident to the variable node ``x_s``, excluding the factor node ``f_i``. Note that each message is a function of the variable ``x_s``.

---

### [Message from a factor node to a variable node] (@id factorVariableMessage)
Consider a part of a factor graph that consists of a group of variable nodes ``\mathcal{X}_i = \{x_s, x_l,...,x_L\}`` ``\subseteq`` ``\mathcal X`` that are neighbours of the factor node ``f_i`` ``\in`` ``\mathcal{F}``. The message ``\mu_{f_i \to x_s}(x_s)`` from the factor node ``f_i`` to the variable node ``x_s`` is defined as a product of all incoming variable node to factor node messages arriving at other incident edges, multiplied by the function ``\psi_i(\mathcal{X}_i)`` associated to the factor node ``f_i``, and marginalised over all of the variables associated with the incoming messages:
```math
    \mu_{f_i \to x_s}(x_s)= \sum\limits_{x_l}\dots\sum\limits_{x_L} \psi_i(\mathcal{X}_i)
    \prod_{x_b \in \mathcal{X}_i\setminus x_s} \mu_{x_b \to f_i}(x_b),
```
where ``\mathcal{X}_i\setminus x_s`` is the set of variable nodes incident to the factor node ``f_i``, excluding the variable node ``x_s``. For continuous variables the summations are simply replaced by integrations:
```math
    \mu_{f_i \to x_s}(x_s)= \int\limits_{x_l}\dots\int\limits_{x_L} \psi_i(\mathcal{X}_i)
    \prod_{x_b \in \mathcal{X}_i\setminus x_s} \big[\mu_{x_b \to f_i}(x_b) \cdot \mathrm{d}x_b\big].
```

---

### [Marginal inference] (@id marginal)
The marginal of the variable node ``x_s`` is obtained as the product of all incoming messages into the variable node ``x_s``:
```math
    p(x_s) =\prod_{f_c \in \mathcal{F}_s} \mu_{f_c \to x_s}(x_s),
```
where ``\mathcal{F}_s`` is the set of factor nodes incident to the variable node ``x_s``.

---

### [Message passing schedule]  (@id MessagePassingSchedule)
The message passing algorithm is an iterative algorithm, and requires a message-passing schedule. Typically, the message updating schedule can be implemented using:
 - synchronous, or
 - forward-backward schedule.

#### [Synchronous schedule]  (@id synchronousSchedule)
The scheduling where messages from variable to factor nodes, and messages from factor nodes to variable nodes, are updated in parallel in respective half-iterations, is known as synchronous scheduling. Synchronous scheduling updates all messages in a given iteration using the output of the previous iteration as an input. The synchronous scheduling allows inference for an arbitrary factor graph structure.

#### [Forward-backward schedule]  (@id treeSchedule)
The forward–backward schedule allows exact inference in tree factor graph. We start by viewing an arbitrary variable node as the root of the factor graph and initiating messages at the leaves of the tree factor graph using. The message passing steps from variable nodes to factor nodes and from factor nodes to variable nodes are then applied recursively until messages have been propagated along every link, and the root node has received messages from all of its neighbours. Each node can send a message towards the root once it has received messages from all of its other neighbours. This step is known as the forward recursion.

The backward recursion starts when the root node received messages from all of its neighbours. It can therefore send out messages to all of its neighbours. These in turn will then have received messages from all of their neighbours and so can send out messages along the links going away from the root, and so on. In this way, messages are passed outwards from the root all the way to the leaves.

By now, a message will have passed in both directions across every link in the graph, and every node will have received
a message from all of its neighbours. Every variable node will have received messages from all of its neighbours, we can readily calculate the marginal distribution for every variable in the graph. The number of messages that have to be computed is given by twice the number of links in the graph and so involves only twice the computation involved in finding a single marginal [1].

#### [Initialisation procedure]  (@id initialisationProcedure)
The Initialisation step starts with messages from singly connected factor nodes to variable nodes. Then, variable nodes forward the incoming messages received from singly connected factor nodes along remaining edges. To ensure this, we are using virtual factor nodes. Hence, the virtual factor node is a singly connected factor node used if the variable node is not directly observed. In a typical scenario, without prior knowledge, the variance of virtual factor nodes tend to infinity for continuous variables. Then, we also improve convergence performance using virtual factor nodes. For discrete variables messages from virtual factor nodes are set to unity.

---

### [References](@id refs)
[1] C. M. Bishop, *Pattern Recognition and Machine Learning* (Information Science and Statistics). Berlin, Heidelberg: Springer-Verlag, 2006.


