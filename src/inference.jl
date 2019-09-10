################################################################################
# Compute messages and marginals
################################################################################


################################################################################
# Compute messages from factor nodes to variable nodes using state-of-the-art
# equations, or damping equations, or using sum compensation without or with
# damping
# Input Data:
#   - m_vf: mean messages from variable node to factor node
#   - v_vf: variance messages from variable node to factor node
#   - m_fv: mean messages from factor node to variable node
#   - vi_fv: inverse variance messages from factor node to variable node
#   - msr: sum vector of row mean messages
#   - vsr: sum vector of row variance messages
#   - Hi: vector of coefficient of indirect factor nodes
#   - bi: vector of measurement means of indecies factor nodes
#   - vi: vector of measurement variances of indirect factor nodes
#   - Ii: indices of indirect factors according to factor nodes (row indices)
#   - Nli: number of links between indirect factor and variable nodes
#   - ah1: current iteration message weights
#   - ah2: previous iteration message weights
#   - evr: error vector of variance summation
# Output Data:
#   - m_fv: mean messages from factor node to variable node
#   - vi_fv: inverse variance messages from factor node to variable node
#-------------------------------------------------------------------------------
 function factor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, Hi, bi, vi, Ii, Nli)
    @inbounds for i = 1:Nli
        m_fv[i] = (bi[Ii[i]] - msr[Ii[i]]) / Hi[i] + m_vf[i]
        vi_fv[i] =  (Hi[i]^2) / (vi[Ii[i]] + vsr[Ii[i]] - Hi[i]^2 * v_vf[i])
    end

    return m_fv, vi_fv
 end

 function dfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, Hi, bi, vi, Ii, Nli, ah1, ah2)
    @inbounds for i = 1:Nli
        m_fv[i] = ah1[i] * ((bi[Ii[i]] - msr[Ii[i]]) / Hi[i] + m_vf[i]) + ah2[i] * m_fv[i]
        vi_fv[i] =  (Hi[i]^2) / (vi[Ii[i]] + vsr[Ii[i]] - Hi[i]^2 * v_vf[i])
    end

    return m_fv, vi_fv
 end

 function nfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, evr, Hi, bi, vi, Ii, Nli)
    @inbounds for i = 1:Nli
        m_fv[i] = (bi[Ii[i]] - msr[Ii[i]]) / Hi[i] + m_vf[i]
        vi_fv[i] =  (Hi[i]^2) / (vi[Ii[i]] + (vsr[Ii[i]] - Hi[i]^2 * v_vf[i]) + evr[Ii[i]])
    end

    return m_fv, vi_fv
 end

 function ndfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, evr, Hi, bi, vi, Ii, Nli, ah1, ah2)
     @inbounds for i = 1:Nli
         m_fv[i] = ah1[i] * ((bi[Ii[i]] - msr[Ii[i]]) / Hi[i] + m_vf[i]) + ah2[i] * m_fv[i]
         vi_fv[i] =  (Hi[i]^2) / (vi[Ii[i]] + (vsr[Ii[i]] - Hi[i]^2 * v_vf[i]) + evr[Ii[i]])
     end

     return m_fv, vi_fv
 end
################################################################################


################################################################################
# Compute messages from variable nodes to factor nodes using state-of-the-art
# equations, or sum compensation
# Input Data:
#   - m_vf: mean messages from variable node to factor node
#   - v_vf: variance messages from variable node to factor node
#   - m_fv: mean messages from factor node to variable node
#   - vi_fv: inverse variance messages from factor node to variable node
#   - md: total means from singly-connected factors to variable nodes
#   - vid: total inverse variance from singly-connected factors to variable nodes
#   - msc: sum vector of column mean messages
#   - vsc: sum vector of column variance messages
#   - Ji: indices of indirect factors according to variable nodes (column indices)
#   - Nli: number of links between indirect factor and variable nodes
#   - evc: error vector of variance summation
# Output Data:
#   - m_vf: mean messages from variable node to factor node
#   - v_vf: variance messages from variable node to factor node
#-------------------------------------------------------------------------------
function variable_to_factor(m_vf, v_vf, m_fv, vi_fv, md, vid, msc, vsc, Ji, Nli)
    @inbounds for i = 1:Nli
        v_vf[i] = 1 / (vsc[Ji[i]] + vid[Ji[i]] - vi_fv[i])
        m_vf[i] = (msc[Ji[i]] - m_fv[i] * vi_fv[i] + md[Ji[i]]) * v_vf[i]
    end

    return m_vf, v_vf
end

function nvariable_to_factor(m_vf, v_vf, m_fv, vi_fv, md, vid, msc, vsc, evc, Ji, Nli)
    @inbounds for i = 1:Nli
        v_vf[i] = 1 / ((vsc[Ji[i]] - vi_fv[i]) + evc[Ji[i]] + vid[Ji[i]])
        m_vf[i] = (msc[Ji[i]] - m_fv[i] * vi_fv[i] + md[Ji[i]]) * v_vf[i]
    end

    return m_vf, v_vf
end
################################################################################


################################################################################
# Compute marginals using state-of-the-art equations, or sum compensation
# Input Data:
#   - md: total means from singly-connected factors to variable nodes
#   - vid: total inverse variance from singly-connected factors to variable nodes
#   - msc: sum vector of column mean messages
#   - vsc: sum vector of column variance messages
#   - Ji: indices of indirect factors according to variable nodes (column indices)
#   - Nv: number of variable nodes
#   - evc: error vector of variance summation
# Output Data:
#   - xbp: estimate values
#-------------------------------------------------------------------------------
function marginal(md, vid, msc, vsc, Ji, Nv)
    xbp = fill(0.0, Nv)

    @inbounds for i = 1:Nv
        xbp[i] = (msc[i] + md[i]) / (vsc[i] + vid[i])
    end

    return xbp
end

function nmarginal(md, vid, msc, vsc, evc, Ji, Nv)
    xbp = fill(0.0, Nv)

    @inbounds for i = 1:Nv
        xbp[i] = (msc[i] + md[i]) / (vsc[i] + vid[i] + evc[i])
    end

    return xbp
end
################################################################################
