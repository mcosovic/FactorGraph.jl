#############################################################
#  The belief propagation using simply summation procedure  #
#############################################################


function bp_simple_passing(
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
end

inference = @elapsed begin
    for i = 1:BUMP
        if i < DAMP
            Mfac, VfacInv = factor_to_variable(Mvar, Vvar, Mfac, VfacInv, Mind, Vind,
                            coeff, coeffInv, row, Nind, factor_colptr, vf)
        else
            Mfac, VfacInv = factor_to_variable_damp(Mvar, Vvar, Mfac, VfacInv, Mind, Vind,
                            coeff, coeffInv, row, Nind, factor_colptr, vf, alpha1, alpha2)
        end
        Mvar, Vvar = variable_to_factor(Mvar, Vvar, Mfac, VfacInv, Mdir, VdirInv,
                     col, Nvariable, variable_colptr, fv)
    end

    for i = (BUMP + 1):MAXI
        if i < DAMP
            Mfac = factor_to_variable_mean(Mvar, Mfac, Mind, coeff, coeffInv,
                                           row, Nind, factor_colptr, vf)
        else
            Mfac = factor_to_variable_damp_mean(Mvar, Mfac, Mind, coeff, coeffInv,
                                                row, Nind, factor_colptr, vf, alpha1, alpha2)
        end
        Mvar = variable_to_factor_mean(Mvar, Vvar, Mfac, VfacInv, Mdir,
                                       col, Nvariable, variable_colptr, fv)
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
