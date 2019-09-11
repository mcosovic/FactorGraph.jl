##########################
#  Summarizing messages  #
##########################


#------------------------------------------------------------------------
# Summarizing messages from variable nodes to factor nodes using simply
# summation, or using Kahan-Babuska algorithm
# Input Data:
#   - Hi: vector of coefficient of indirect factor nodes
#   - Ii: indices of indirect factors according to factor nodes (rows)
#   - Nli: number of links between indirect factor and variable nodes
#   - m_vf: mean messages from variable node to factor node
#   - v_vf: variance messages from variable node to factor node
#   - msr: sum vector of row mean messages
#   - vsr: sum vector of row variance messages
#   - evr: error vector of variance summation
#   - iter: current iteration
#   - BUMP: canceled variances
# Output Data:
#   - msr: sum vector of row mean messages
#   - vsr: sum vector of row variance messages
#------------------------------------------------------------------------
function sum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr, iter, BUMP)
    @inbounds for i = 1:Nli
        msr[Ii[i]] += (Hi[i] * m_vf[i])

        if iter < BUMP
            vsr[Ii[i]] += (Hi[i]^2 * v_vf[i])
        end
    end

    return msr, vsr
end

function nsum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr, evr, iter, BUMP)
    @inbounds for i = 1:Nli
        msr[Ii[i]] += (Hi[i] * m_vf[i])

        if iter < BUMP
            x = Hi[i]^2 * v_vf[i]
            t = vsr[Ii[i]] + x
            if abs(vsr[Ii[i]]) >= abs(x)
                evr[Ii[i]] += (vsr[Ii[i]] - t) + x
            else
                evr[Ii[i]] += (x - t) + vsr[Ii[i]]
            end
            vsr[Ii[i]] = t
        end
    end

   return msr, vsr, evr
end
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# Summarizing messages from factor nodes to variable nodes using simply
# summation, or using Kahan-Babuska algorithm
# Input Data:
#   - Ji: indices of indirect factors according to variable nodes (columns)
#   - Nli: number of links between indirect factor and variable nodes
#   - m_fv: mean messages from factor node to variable node
#   - vi_fv: inverse variance messages from factor node to variable node
#   - msc: sum vector of column mean messages
#   - vsc: sum vector of column variance messages
#   - evc: error vector of variance summation
#   - iter: current iteration
#   - BUMP: canceled variances
# Output Data:
#   - msc: sum vector of column mean messages
#   - vsc: sum vector of column variance messages
#------------------------------------------------------------------------
function sum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc, iter, BUMP)
    @inbounds for i = 1:Nli
        msc[Ji[i]] += m_fv[i] * vi_fv[i]

        if iter < BUMP
            vsc[Ji[i]] += vi_fv[i]
        end
    end

    return msc, vsc
end

function nsum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc, evc, iter, BUMP)
    @inbounds for i = 1:Nli
        msc[Ji[i]] += m_fv[i] * vi_fv[i]

        if iter < BUMP
            t = vsc[Ji[i]] + vi_fv[i]
            if abs(vsc[Ji[i]]) >= abs(vi_fv[i])
                evc[Ji[i]] += (vsc[Ji[i]] - t) + vi_fv[i]
            else
                evc[Ji[i]] += (vi_fv[i] - t) + vsc[Ji[i]]
            end
            vsc[Ji[i]] = t
        end
    end

   return msc, vsc, evc
end
#------------------------------------------------------------------------
