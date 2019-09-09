function bpn(data::String="data_33_14", MAXI::Int64=20, DAMP::Int64=10, PROB::Float64=0.6, ALPH::Float64=0.4, MEAN::Float64=0.0, VARI::Float64=1e3; TIME::String = "off", ERROR::String = "off")
    H, b, v = model(data)

    fgraph = @elapsed begin
        Nf, Nv, T = graph(H)                                                                                                # FactorGraph
        Nld, Nli, dir = links(Nf, T)                                                                                        # FactorGraph
        vir = virtuals(Nv, dir)                                                                                             # FactorGraph
        Ii, Ji, Ni, bi, vi, Hi, md, vid = factors(Nf, Nv, Nld, Nli, T, b, v, vir, MEAN, VARI)                               # FactorGraph
    end

    init = @elapsed begin
        m_fv, vi_fv, m_vf, v_vf = load_messages(Hi)                                                                         # AuxiliaryFunction
        ah1, ah2 = damping(Nli, ALPH, PROB)                                                                                 # InitializeMessages
        m_vf, v_vf = forward_directs(Hi, Ji, Nli, md, vid, v_vf, m_vf)                                                      # InitializeMessages

        msr, vsr, evr, msc, vsc, evc = nload_sum(Nv, Ni)                                                                    # AuxiliaryFunction
        msr, vsr, evr = nsum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr, evr)                                                   # SummationBelief
        m_fv, vi_fv = nfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, evr, Hi, bi, vi, Ii, Nli)                      # BeliefPropagation
        msr, vsr, evr, msc, vsc, evc = nclear_sum(msr, vsr, evr, msc, vsc, evc)                                             # AuxiliaryFunction
    end

    infe = @elapsed begin
        for i = 1:MAXI
            msc, vsc, evc = nsum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc, evc)                                                  # SummationBelief
            m_vf, v_vf = nvariable_to_factor(m_vf, v_vf, m_fv, vi_fv, md, vid, msc, vsc, evc, Ji, Nli)                      # BeliefPropagation
            msr, vsr, evr = nsum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr, evr)                                               # SummationBelief
            if i < DAMP
                m_fv, vi_fv = nfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, evr, Hi, bi, vi, Ii, Nli)              # BeliefPropagation
            else
                m_fv, vi_fv = ndfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, evr, Hi, bi, vi, Ii, Nli, ah1, ah2)   # BeliefPropagation
            end
            msr, vsr, evr, msc, vsc, evc = nclear_sum(msr, vsr, evr, msc, vsc, evc)                                         # AuxiliaryFunction
        end
    end

    solu = @elapsed begin
        msc, vsc, evc = nsum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc, evc)                                                      # SummationBelief
        xbp = nmarginal(md, vid, msc, vsc, evc, Ji, Nv)                                                                     # BeliefPropagation
    end

    if TIME == "on"
       bp_time(fgraph, init, infe, solu)
    end

    if ERROR == "on"
        wls = @elapsed begin
            xwls = wlsMldivide(H, b, v)
        end

        wrss_bp, wrss_wls = wrss(H, b, v, xbp, xwls)
        rmse_bp, rmse_wls = rmse(H, b, v, xbp, xwls)

        wls_vs_bp(wls, wrss_wls, wrss_bp, rmse_wls, rmse_bp, xbp, xwls)
    end


    return xbp
end
