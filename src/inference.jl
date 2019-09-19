####################################
#  Compute messages and marginals  #
####################################


#-------------------------------------------------------------------------------
# Compute messages using state-of-the-art equations with simply summation
#-------------------------------------------------------------------------------
function factor_to_variable(
    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var,
    Mind, Vind, coeff, coeffInv,
    row, Nind, row_colptr, to_var)

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0

        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mrow += (coeff[j] * Mvar_fac[j])
            Vrow += (coeff[j]^2 * Vvar_fac[j])
        end
        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mfac_var[to_var[j]] = (Mind[row[j]] - Mrow) * coeffInv[j] + Mvar_fac[j]
            Wfac_var[to_var[j]] = (coeff[j]^2) / (Vind[row[j]] + Vrow - coeff[j]^2 * Vvar_fac[j])
        end
    end

    return Mfac_var, Wfac_var
end

function factor_to_variable_damp(
    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var,
    Mind, Vind, coeff, coeffInv,
    row, Nind, row_colptr, to_var,
    alpha1, alpha2)

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0

        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mrow += (coeff[j] * Mvar_fac[j])
            Vrow += (coeff[j]^2 * Vvar_fac[j])
        end
        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mfac_var[to_var[j]] = alpha1[j] * ((Mind[row[j]] - Mrow) * coeffInv[j] + Mvar_fac[j]) + alpha2[j] * Mfac_var[to_var[j]]
            Wfac_var[to_var[j]] = (coeff[j]^2) / (Vind[row[j]] + Vrow - coeff[j]^2 * Vvar_fac[j])
        end
    end

    return Mfac_var, Wfac_var
end

function variable_to_factor(
    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var,
    Mdir, VdirInv,
    col, Nvariable, col_colptr, to_fac)

    @inbounds for i = 1:Nvariable
        Mcol = 0.0
        Wcol = 0.0

        for j in (col_colptr[i]):col_colptr[i+1]-1
            Mcol += Mfac_var[j] * Wfac_var[j]
            Wcol += Wfac_var[j]
        end
        for j in (col_colptr[i]):col_colptr[i+1]-1
            Vvar_fac[to_fac[j]] = 1 / (Wcol + VdirInv[i] - Wfac_var[j])
            Mvar_fac[to_fac[j]] = (Mcol - Mfac_var[j] * Wfac_var[j] + Mdir[i]) * Vvar_fac[to_fac[j]]
        end
    end

    return Mvar_fac, Vvar_fac
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Compute messages using sum compensation, Kahan-Babuska algorithm
#------------------------------------------------------------------------------
function factor_to_variable_kahan(
    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var,
    Mind, Vind, coeff, coeffInv,
    row, Nind, row_colptr, to_var)

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0
        error = 0.0

        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mrow += (coeff[j] * Mvar_fac[j])

            x = coeff[j]^2 * Vvar_fac[j]
            t = Vrow + x
            if abs(Vrow) >= abs(x)
                error += (Vrow - t) + x
            else
                error += (x - t) + Vrow
            end
            Vrow = t
        end
        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mfac_var[to_var[j]] = (Mind[row[j]] - Mrow) * coeffInv[j] + Mvar_fac[j]
            Wfac_var[to_var[j]] = (coeff[j]^2) / (Vind[row[j]] + (Vrow - coeff[j]^2 * Vvar_fac[j]) + error)
        end
    end

    return Mfac_var, Wfac_var
end

function factor_to_variable_damp_kahan(
    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var,
    Mind, Vind, coeff, coeffInv,
    row, Nind, row_colptr, to_var,
    alpha1, alpha2)

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0
        error = 0.0

        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mrow += (coeff[j] * Mvar_fac[j])

            x = coeff[j]^2 * Vvar_fac[j]
            t = Vrow + x
            if abs(Vrow) >= abs(x)
                error += (Vrow - t) + x
            else
                error += (x - t) + Vrow
            end
            Vrow = t
        end
        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mfac_var[to_var[j]] = alpha1[j] * ((Mind[row[j]] - Mrow) * coeffInv[j] + Mvar_fac[j]) + alpha2[j] * Mfac_var[to_var[j]]
            Wfac_var[to_var[j]] = (coeff[j]^2) / (Vind[row[j]] + (Vrow - coeff[j]^2 * Vvar_fac[j]) + error)
        end
    end

    return Mfac_var, Wfac_var
end

function variable_to_factor_kahan(
    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var,
    Mdir, VdirInv,
    col, Nvariable, col_colptr, to_fac)

    @inbounds for i = 1:Nvariable
        Mcol = 0.0
        Wcol = 0.0
        error = 0.0

        for j in (col_colptr[i]):col_colptr[i+1]-1
            Mcol += Mfac_var[j] * Wfac_var[j]

            t = Wcol + Wfac_var[j]
            if abs(Wcol) >= abs(Wfac_var[j])
                error += (Wcol - t) + Wfac_var[j]
            else
                error += (Wfac_var[j] - t) + Wcol
            end
            Wcol = t
        end
        for j in (col_colptr[i]):col_colptr[i+1]-1
            Vvar_fac[to_fac[j]] = 1 / ((Wcol - Wfac_var[j]) + error + VdirInv[i])
            Mvar_fac[to_fac[j]] = (Mcol - Mfac_var[j] * Wfac_var[j] + Mdir[i]) * Vvar_fac[to_fac[j]]
        end
    end

    return Mvar_fac, Vvar_fac
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Compute messages using recursion with simply summation
#-------------------------------------------------------------------------------
function factor_recursion(
    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var,
    Mdir, VdirInv, Mind, Vind,
    coeff, coeffInv, row, col, Mcol, Wcol,
    Nind, Nlink, row_colptr, to_var)

    @inbounds for i = 1:Nlink
        Mcol[col[i]] += Mfac_var[to_var[i]] * Wfac_var[to_var[i]]
        Wcol[col[i]] += Wfac_var[to_var[i]]
    end

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0

        for j in (row_colptr[i]):row_colptr[i+1]-1
            Vvar_fac[j] = 1 / (Wcol[col[j]] + VdirInv[col[j]] - Wfac_var[to_var[j]])
            Mvar_fac[j] = (Mcol[col[j]] - Mfac_var[to_var[j]] * Wfac_var[to_var[j]] + Mdir[col[j]]) * Vvar_fac[j]

            Mrow += (coeff[j] * Mvar_fac[j])
            Vrow += (coeff[j]^2 * Vvar_fac[j])
        end
        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mfac_var[to_var[j]] = (Mind[row[j]] - Mrow) * coeffInv[j] + Mvar_fac[j]
            Wfac_var[to_var[j]] = (coeff[j]^2) / (Vind[row[j]] + Vrow - coeff[j]^2 * Vvar_fac[j])
        end
    end

    fill!(Mcol, 0.0)
    fill!(Wcol, 0.0)

    return Mfac_var, Wfac_var, Mcol, Wcol
end

function factor_recursion_damp(
    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var,
    Mdir, VdirInv, Mind, Vind,
    coeff, coeffInv, row, col, Mcol, Wcol,
    Nind, Nlink, row_colptr, to_var, alpha1, alpha2)

    @inbounds for i = 1:Nlink
        Mcol[col[i]] += Mfac_var[to_var[i]] * Wfac_var[to_var[i]]
        Wcol[col[i]] += Wfac_var[to_var[i]]
    end

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0

        for j in (row_colptr[i]):row_colptr[i+1]-1
            Vvar_fac[j] = 1 / (Wcol[col[j]] + VdirInv[col[j]] - Wfac_var[to_var[j]])
            Mvar_fac[j] = (Mcol[col[j]] - Mfac_var[to_var[j]] * Wfac_var[to_var[j]] + Mdir[col[j]]) * Vvar_fac[j]

            Mrow += (coeff[j] * Mvar_fac[j])
            Vrow += (coeff[j]^2 * Vvar_fac[j])
        end
        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mfac_var[to_var[j]] = alpha1[j] * ((Mind[row[j]] - Mrow) * coeffInv[j] + Mvar_fac[j]) + alpha2[j] * Mfac_var[to_var[j]]
            Wfac_var[to_var[j]] = (coeff[j]^2) / (Vind[row[j]] + Vrow - coeff[j]^2 * Vvar_fac[j])
        end
    end

    fill!(Mcol, 0.0)
    fill!(Wcol, 0.0)

    return Mfac_var, Wfac_var, Mcol, Wcol
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Compute marginals
#-------------------------------------------------------------------------------
function marginal(Mfac_var, Wfac_var, Mdir, VdirInv, col, Nvariable, col_colptr)
    Xbp = fill(0.0, Nvariable)

    @inbounds for i = 1:Nvariable
        Mcol = 0.0
        Wcol = 0.0

        for j in (col_colptr[i]):col_colptr[i+1]-1
            Mcol += Mfac_var[j] * Wfac_var[j]
            Wcol += Wfac_var[j]
        end
        Xbp[i] = (Mcol + Mdir[i]) / (Wcol + VdirInv[i])
    end

    return Xbp
end
#-------------------------------------------------------------------------------
