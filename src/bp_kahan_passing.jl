####################################################################
#  The belief propagation using Kahan-Babuska summation algorithm  #
####################################################################


function bp_kahan_passing(
    jacobian, observation, noise,
    MAXI, DAMP, BUMP,
    PROB, ALPH,
    MEAN, VARI,
    TIME, ERROR, STATISTIC)

factorgraph = @elapsed begin
    Nfac, Nvar, jacobianT = graph(jacobian)
    Ndir, Nlink, dir_position = links(Nfac, jacobianT)
    virtual = virtuals(Nvar, dir_position)

    row, row_colptr, col, col_colptr, Nind, Mind, Vind, coeff, coeffInv, Mdir, Wdir  =
        factors(Nfac, Nvar, Ndir, Nlink, jacobianT, observation, noise, virtual, MEAN, VARI)
end

initialize = @elapsed begin
    Mfac_var, Wfac_var, Mvar_fac, Vvar_fac = load_messages(coeff)
    alpha1, alpha2 = damping(Nlink, ALPH, PROB)
    to_fac, to_var = keep_order(Nlink, row, col)
    Mvar_fac, Vvar_fac = forward_directs_to_links(Mvar_fac, Vvar_fac, Mdir, Wdir, Nvar, col_colptr, to_fac)
end

inference = @elapsed begin
    for i = 1:BUMP
        if i < DAMP
            Mfac_var, Wfac_var = factor_to_variable_kahan(Mvar_fac, Vvar_fac, Mfac_var, Wfac_var, Mind, Vind,
                                                          coeff, coeffInv, row, Nind, row_colptr, to_var)
        else
            Mfac_var, Wfac_var = factor_to_variable_damp_kahan(Mvar_fac, Vvar_fac, Mfac_var, Wfac_var, Mind, Vind, coeff,
                                                               coeffInv, row, Nind, row_colptr, to_var, alpha1, alpha2)
        end
        Mvar_fac, Vvar_fac = variable_to_factor_kahan(Mvar_fac, Vvar_fac, Mfac_var, Wfac_var,
                                                      Mdir, Wdir, col, Nvar, col_colptr, to_fac)
    end

    for i = (BUMP + 1):MAXI
        if i < DAMP
            Mfac_var = factor_to_variable_mean(Mvar_fac, Mfac_var, Mind, coeff,
                                               coeffInv, row, Nind, row_colptr, to_var)
        else
            Mfac_var = factor_to_variable_damp_mean(Mvar_fac, Mfac_var, Mind, coeff, coeffInv,
                                                    row, Nind, row_colptr, to_var, alpha1, alpha2)
        end
        Mvar_fac = variable_to_factor_mean(Mvar_fac, Vvar_fac, Mfac_var, Wfac_var, Mdir,
                                           col, Nvar, col_colptr, to_fac)
    end
end

solution = @elapsed begin
    Xbp = marginal(Mfac_var, Wfac_var, Mdir, Wdir, col, Nvar, col_colptr)
end

    if STATISTIC == "on"
        graph_statistic(Nfac, Nvar, Ndir, Nlink, virtual, noise)
    end
    if TIME == "on"
        bp_time(factorgraph, initialize, inference, solution)
    end
    if ERROR == "on"
        errors(jacobian, observation, noise, Xbp, TIME)
    end

    return Xbp
end
