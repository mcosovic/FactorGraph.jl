#############################################################
#  The belief propagation using simply summation procedure  #
#############################################################


function bps(H, b, v, MAXI, DAMP, PROB, ALPH, MEAN, VARI, TIME)

    fgraph = @elapsed begin
        Nf, Nv, T = graph(H)
        Nld, Nli, dir = links(Nf, T)
        vir = virtuals(Nv, dir)
        Ii, Ji, Ni, bi, vi, Hi, md, vid = factors(Nf, Nv, Nld, Nli, T, b, v, vir, MEAN, VARI)
    end

    init = @elapsed begin
        m_fv, vi_fv, m_vf, v_vf = load_messages(Hi)
        ah1, ah2 = damping(Nli, ALPH, PROB)
        m_vf, v_vf = forward_directs(Hi, Ji, Nli, md, vid, v_vf, m_vf)

        msr, vsr, msc, vsc = load_sum(Nv, Ni)
        msr, vsr = sum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr)
        m_fv, vi_fv = factor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, Hi, bi, vi, Ii, Nli)
        msr, vsr, msc, vsc = clear_sum(msr, vsr, msc, vsc)
    end

    infe = @elapsed begin
        for i = 1:MAXI
            msc, vsc = sum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc)
            m_vf, v_vf = variable_to_factor(m_vf, v_vf, m_fv, vi_fv, md, vid, msc, vsc, Ji, Nli)
            msr, vsr = sum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr)
            if i < DAMP
                m_fv, vi_fv = factor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, Hi, bi, vi, Ii, Nli)
            else
                m_fv, vi_fv = dfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, Hi, bi, vi, Ii, Nli, ah1, ah2)
            end
            msr, vsr, msc, vsc = clear_sum(msr, vsr, msc, vsc)
        end
    end

    solu = @elapsed begin
        msc, vsc = sum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc)
        xbp = marginal(md, vid, msc, vsc, Ji, Nv)
    end

    bp_time(fgraph, init, infe, solu, TIME)

    return xbp
 end
