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

    Mrow, Vrow, error_row, Mcol, VcolInv, error_col =
        load_sum_kahan(Nvariable, Nind)

    Mrow, Vrow, error_row =
        sum_rows_kahan(Mvar, Vvar, Mrow, Vrow, error_row,
                       coeff, row, Nlink)

    Mfac, VfacInv, Mrow, Vrow, error_row =
        factor_to_variable_kahan(Mvar, Vvar, Mfac, VfacInv, Mrow, Vrow,
                                 error_row, Mind, Vind, coeff, coeffInv, row,
                                 Nlink)

end

inference = @elapsed begin
    for i = 1:BUMP
        Mcol, VcolInv, error_col =
            sum_cols_kahan(Mfac, VfacInv, Mcol, VcolInv, error_col, col, Nlink)

        Mvar, Vvar, Mcol, VcolInv, error_col =
            variable_to_factor_kahan(Mvar, Vvar, Mfac, VfacInv, Mcol, VcolInv,
                                     error_col, Mdir, VdirInv, col, Nlink)

        Mrow, Vrow, error_row =
            sum_rows_kahan(Mvar, Vvar, Mrow, Vrow, error_row,
                           coeff, row, Nlink)

            if i < DAMP
                Mfac, VfacInv, Mrow, Vrow, error_row =
                    factor_to_variable_kahan(Mvar, Vvar, Mfac, VfacInv, Mrow,
                                             Vrow, error_row, Mind, Vind, coeff,
                                             coeffInv, row, Nlink)
            else
                Mfac, VfacInv, Mrow, Vrow, error_row =
                factor_to_variable_kahan_damp(Mvar, Vvar, Mfac, VfacInv, Mrow,
                                              Vrow, error_row, Mind, Vind, coeff,
                                              coeffInv, row, Nlink, alpha1,
                                              alpha2)
            end
    end

    for i = 1:(MAXI - BUMP)
        Mcol = sum_cols_mean(Mfac, VfacInv, Mcol, col, Nlink)

        Mvar, Mcol =
            variable_to_factor_mean(Mvar, Vvar, Mfac, VfacInv, Mcol, Mdir,
                                    col, Nlink)

        Mrow = sum_rows_mean(Mvar, Mrow, coeff, row, Nlink)

        if i < DAMP
            Mfac, Mrow =
                factor_to_variable_mean(Mvar, Mfac, Mrow, Mind, coeffInv,
                                        row, Nlink)
        else
            Mfac, Mrow =
                factor_to_variable_mean_damp(Mvar, Mfac, Mrow, Mind, coeffInv,
                                            row, Nlink, alpha1, alpha2)
        end
    end
end

solution = @elapsed begin
    Mcol, VcolInv, error_col =
        sum_cols_kahan(Mfac, VfacInv, Mcol, VcolInv, error_col, col, Nlink)

    Xbp = marginal_kahan(Mcol, VcolInv, error_col, Mdir, VdirInv, col, Nvariable)
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
