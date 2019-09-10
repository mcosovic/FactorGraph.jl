#############################################################
#  Initialize the algorithm and produce damping parameters  #
#############################################################


#------------------------------------------------------------------------
# Pass messages from singly-connected factor nodes to all indirect links
# Input Data:
#   - Hi: vector of coefficient of indirect factor nodes
#   - Ji: indices of indirect factors according to variable nodes (columns)
#   - Nli: number of links between indirect factor and variable nodes
#   - md: like means from singly-connected factors
#   - vid: inverse variance from singly-connected factors
#   - m_vf: mean messages from variable node to factor node
#   - v_vf: variance messages from variable node to factor node
# Output Data:
#   - m_vf: mean messages from variable node to factor node
#   - v_vf: variance messages from variable node to factor node
#------------------------------------------------------------------------
function forward_directs(Hi, Ji, Nli, md, vid, v_vf, m_vf)
    @inbounds for i = 1:Nli
        v_vf[i] = 1 / vid[Ji[i]]
        m_vf[i] = md[Ji[i]] / vid[Ji[i]]
    end

    return m_vf, v_vf
end
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# Set damping parameters
# Input Data:
#   - Nli: number of links between indirect factor and variable node
#   - ALPH: damping weights
#   - PROB: a Bernoulli random variable with probability
# Output Data:
#   - ah1: current iteration message weights
#   - ah2: previous iteration message weights
#------------------------------------------------------------------------
function damping(Nli, ALPH, PROB)
    wow = randsubseq(collect(1:Nli), PROB)

    ah1 = fill(1.0, Nli)
    ah1[wow] .= 1.0 - ALPH

    ah2 = fill(0.0, Nli)
    ah2[wow] .= ALPH

    return ah1, ah2
end
#------------------------------------------------------------------------
