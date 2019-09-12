####################################################################
#  The belief propagation using Kahan-Babuska summation algorithm  #
####################################################################


function bpn(H, b, v, MAXI, DAMP, BUMP, PROB, ALPH, MEAN, VARI, TIME)

    fgraph = @elapsed begin
        Nf, Nv, T = graph(H)
        Nld, Nli, dir = links(Nf, T)
        vir = virtuals(Nv, dir)
        Ii, Ji, Ni, bi, vi, Hi, Hr, md, vid = factors(Nf, Nv, Nld, Nli, T, b, v, vir, MEAN, VARI)
    end

    init = @elapsed begin
        m_fv, vi_fv, m_vf, v_vf = load_messages(Hi)
        ah1, ah2 = damping(Nli, ALPH, PROB)
        m_vf, v_vf = forward_directs(Hi, Ji, Nli, md, vid, v_vf, m_vf)

        msr, vsr, evr, msc, vsc, evc = nload_sum(Nv, Ni)
        msr, vsr, evr = nsum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr, evr, -1, BUMP)
        m_fv, vi_fv, msr, vsr, evr = nfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, evr, Hi, Hr, bi, vi, Ii, Nli, -1, BUMP)
    end

    infe = @elapsed begin
        for i = 1:MAXI
            msc, vsc, evc = nsum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc, evc, i, BUMP)
            m_vf, v_vf, msc, vsc, evc = nvariable_to_factor(m_vf, v_vf, m_fv, vi_fv, md, vid, msc, vsc, evc, Ji, Nli, i, BUMP)
            msr, vsr, evr = nsum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr, evr, i, BUMP)
            if i < DAMP
                m_fv, vi_fv, msr, vsr, evr = nfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, evr, Hi, Hr, bi, vi, Ii, Nli, i, BUMP)
            else
                m_fv, vi_fv, msr, vsr, evr = ndfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, evr, Hi, Hr, bi, vi, Ii, Nli, ah1, ah2, i, BUMP)
            end
        end
    end

    solu = @elapsed begin
        msc, vsc, evc = nsum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc, evc, -1, BUMP)
        xbp = nmarginal(md, vid, msc, vsc, evc, Ji, Nv)
    end

    graph_statistic(Nf, Nv, Nld, Nli, vir, v)
    bp_time(fgraph, init, infe, solu, TIME)

    return xbp
end
