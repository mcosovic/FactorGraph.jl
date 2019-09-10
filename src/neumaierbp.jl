####################################################################
#  The belief propagation using Kahan-Babuska summation algorithm  #
####################################################################


function bpn(H, b, v, MAXI, DAMP, PROB, ALPH, MEAN, VARI, TIME)

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

        msr, vsr, evr, msc, vsc, evc = nload_sum(Nv, Ni)
        msr, vsr, evr = nsum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr, evr)
        m_fv, vi_fv = nfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, evr, Hi, bi, vi, Ii, Nli)                      
        msr, vsr, evr, msc, vsc, evc = nclear_sum(msr, vsr, evr, msc, vsc, evc)
    end

    infe = @elapsed begin
        for i = 1:MAXI
            msc, vsc, evc = nsum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc, evc)
            m_vf, v_vf = nvariable_to_factor(m_vf, v_vf, m_fv, vi_fv, md, vid, msc, vsc, evc, Ji, Nli)                      
            msr, vsr, evr = nsum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr, evr)
            if i < DAMP
                m_fv, vi_fv = nfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, evr, Hi, bi, vi, Ii, Nli)              
            else
                m_fv, vi_fv = ndfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, evr, Hi, bi, vi, Ii, Nli, ah1, ah2)   
            end
            msr, vsr, evr, msc, vsc, evc = nclear_sum(msr, vsr, evr, msc, vsc, evc)
        end
    end

    solu = @elapsed begin
        msc, vsc, evc = nsum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc, evc)
        xbp = nmarginal(md, vid, msc, vsc, evc, Ji, Nv)                                                                     
    end

    bp_time(fgraph, init, infe, solu, TIME)

    return xbp
end
