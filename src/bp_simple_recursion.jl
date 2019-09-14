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

    row, col, Nind, Mind, Vind, coeff, coeffInv, Mdir, VdirInv =
        factors(Nfactor, Nvariable, Ndir, Nlink, jacobianT, observation,
                noise, virtual, MEAN, VARI)
end

initialize = @elapsed begin
    Mfac, VfacInv, Mvar, Vvar = load_messages(coeff)
    alpha1, alpha2 = damping(Nlink, ALPH, PROB)
    Mvar, Vvar = forward_directs_to_links(Mvar, Vvar, Mind, Vind, col, Nlink)

    Mrow, Vrow, Mcol, VcolInv, Vaux, Maux =
        load_sum_recursion(Nvariable, Nind, Nlink)

    Mrow, Vrow = sum_rows(Mvar, Vvar, Mrow, Vrow, coeff, row, Nlink)

    Mfac, VfacInv, Mrow, Vrow =
        factor_to_variable(Mvar, Vvar, Mfac, VfacInv, Mrow, Vrow, Mind, Vind,
                           coeff, coeffInv, row, Nlink)
end

inference = @elapsed begin
    for i = 1:BUMP
        Mcol, VcolInv =
            sum_cols(Mfac, VfacInv, Mcol, VcolInv, col, Nlink)

        Mrow, Vrow, Vaux, Maux =
            sum_rows_recursion(Mfac, VfacInv, Mrow, Vrow, Mcol, VcolInv, Maux,
                               Vaux, Mdir, VdirInv, coeff, row, col, Nlink)

        if i < DAMP
            Mfac, VfacInv, Mrow, Vrow, Mcol, VcolInv =
                factor_recursion(Mvar, Vvar, Mfac, VfacInv, Mrow, Vrow, Mcol,
                                 VcolInv, Vaux, Maux, Mind, Vind, coeff,
                                 coeffInv, row, Nlink)
        else
            Mfac, VfacInv, Mrow, Vrow, Mcol, VcolInv =
                factor_recursion_damp(Mvar, Vvar, Mfac, VfacInv, Mrow, Vrow,
                                      Mcol, VcolInv, Vaux, Maux, Mind, Vind,
                                      coeff, coeffInv, row, Nlink, alpha1,
                                      alpha2)
        end
    end
    
    for i = (BUMP + 1):MAXI
        Mcol = sum_cols_mean(Mfac, VfacInv, Mcol, col, Nlink)

        Mrow, Maux =
            sum_rows_mean_recursion(Mfac, VfacInv, Mrow, Mcol, Maux, Vaux, Mdir,
                                    coeff, row, col, Nlink)
        if i < DAMP
            Mfac, Mrow, Mcol =
                factor_mean_recursion(Mfac, Mrow, Mcol, Maux, Mind,coeffInv,
                                      row, Nlink)
        else
            Mfac, Mrow, Mcol =
                factor_recursion_damp(Mfac, Mrow, Mcol, Maux, Mind, coeffInv,
                                      row, Nlink, alpha1, alpha2)
        end
    end
end

solution = @elapsed begin
    fill!(VcolInv, 0.0)
    Mcol, VcolInv = sum_cols(Mfac, VfacInv, Mcol, VcolInv, col, Nlink)
    Xbp = marginal(Mcol, VcolInv, Mdir, VdirInv, col, Nvariable)
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
