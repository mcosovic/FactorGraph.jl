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
    Nfactor, Nvariable, jacobianT = graph(jacobian)
    Ndir, Nlink, dir_position = links(Nfactor, jacobianT)
    virtual = virtuals(Nvariable, dir_position)

    row, col, Nind, Mind, Vind, coeff, coeffInv, Mdir, VdirInv,
    variable_colptr, factor_colptr =
        factors(Nfactor, Nvariable, Ndir, Nlink, jacobianT, observation,
                noise, virtual, MEAN, VARI)
end

initialize = @elapsed begin
    Mfac, VfacInv, Mvar, Vvar = load_messages(coeff)
    alpha1, alpha2 = damping(Nlink, ALPH, PROB)
    fv, vf = keep_order(Nlink, row, col)
    Mvar, Vvar = forward_directs_to_links(Mvar, Vvar, Mdir, VdirInv, Nvariable, variable_colptr, fv)
    Mfac, VfacInv = factor_to_variable(Mvar, Vvar, Mfac, VfacInv, Mind, Vind, coeff,
                                       coeffInv, row, Nind, factor_colptr, vf)
    Mcol, VcolInv = load_sum_col_recursion(Nvariable)
    end

inference = @elapsed begin
    for i = 1:BUMP
        if i < DAMP
            Mfac, VfacInv, Mcol, VcolInv = factor_recursion(
                    Mvar, Vvar, Mfac, VfacInv, Mdir, VdirInv, Mind, Vind,
                    coeff, coeffInv, row, col, Mcol, VcolInv,
                    Nind, Nlink, factor_colptr, vf)
        else
            Mfac, VfacInv, Mcol, VcolInv = factor_recursion_damp(
                    Mvar, Vvar, Mfac, VfacInv, Mdir, VdirInv, Mind, Vind,
                    coeff, coeffInv, row, col, Mcol, VcolInv,
                    Nind, Nlink, factor_colptr, vf, alpha1, alpha2)
        end
    end

    for i = (BUMP + 1):MAXI
        if i < DAMP
            Mfac, Mcol = factor_recursion_mean(
                    Mvar, Vvar, Mfac, VfacInv, Mdir, Mind,
                    coeff, coeffInv, row, col, Mcol, VcolInv,
                    Nind, Nlink, factor_colptr, vf)
        else
            Mfac, Mcol = factor_recursion_damp_mean(
                    Mvar, Vvar, Mfac, VfacInv, Mdir, Mind,
                    coeff, coeffInv, row, col, Mcol, VcolInv,
                    Nind, Nlink, factor_colptr, vf, alpha1, alpha2)

        end
    end
end

solution = @elapsed begin
    Xbp = marginal(Mfac, VfacInv, Mdir, VdirInv, col, Nvariable, variable_colptr)
end

    if STATISTIC == "on"
        graph_statistic(Nfactor, Nvariable, Ndir, Nlink, virtual, noise)
    end
    if TIME == "on"
        bp_time(factorgraph, initialize, inference, solution)
    end
    if ERROR == "on"
        errors(jacobian, observation, noise, Xbp, TIME)
    end

    return Xbp
end
