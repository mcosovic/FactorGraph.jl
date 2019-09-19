####################################################################
#  The belief propagation using Kahan-Babuska summation algorithm  #
####################################################################


function bp_simple_recursion(
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
    Mfac_var, Wfac_var = factor_to_variable(Mvar_fac, Vvar_fac, Mfac_var, Wfac_var, Mind, Vind,
                                            coeff, coeffInv, row, Nind, row_colptr, to_var)
    Mcol, Wcol = load_sum_col_recursion(Nvar)
    end

inference = @elapsed begin
    for i = 1:BUMP
        if i < DAMP
            Mfac_var, Wfac_var, Mcol, Wcol = factor_recursion(
                    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var, Mdir, Wdir, Mind, Vind,
                    coeff, coeffInv, row, col, Mcol, Wcol,
                    Nind, Nlink, row_colptr, to_var, Nvar, col_colptr)
        else
            Mfac_var, Wfac_var, Mcol, Wcol = factor_recursion_damp(
                    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var, Mdir, Wdir, Mind, Vind,
                    coeff, coeffInv, row, col, Mcol, Wcol,
                    Nind, Nlink, row_colptr, to_var, alpha1, alpha2, Nvar, col_colptr)
        end
    end

    for i = (BUMP + 1):MAXI
        if i < DAMP
            Mfac_var, Mcol = factor_recursion_mean(
                    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var, Mdir, Mind,
                    coeff, coeffInv, row, col, Mcol, Wcol,
                    Nind, Nlink, row_colptr, to_var, Nvar, col_colptr)
        else
            Mfac_var, Mcol = factor_recursion_damp_mean(
                    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var, Mdir, Mind,
                    coeff, coeffInv, row, col, Mcol, Wcol,
                    Nind, Nlink, row_colptr, to_var, alpha1, alpha2, Nvar, col_colptr)

        end
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
