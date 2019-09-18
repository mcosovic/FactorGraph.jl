####################################
#  Compute messages and marginals  #
####################################


#-------------------------------------------------------------------------------
# Compute messages using state-of-the-art equations with simply summation
#-------------------------------------------------------------------------------
function factor_to_variable(
    Mvar, Vvar, Mfac, VfacInv,
    Mind, Vind, coeff, coeffInv,
    row, Nind, factor_colptr, vf)

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0

        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mrow += (coeff[j] * Mvar[j])
            Vrow += (coeff[j]^2 * Vvar[j])
        end
        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mfac[vf[j]] = (Mind[row[j]] - Mrow) * coeffInv[j] + Mvar[j]
            VfacInv[vf[j]] = (coeff[j]^2) / (Vind[row[j]] + Vrow - coeff[j]^2 * Vvar[j])
        end
    end

    return Mfac, VfacInv
end

function factor_to_variable_damp(
    Mvar, Vvar, Mfac, VfacInv,
    Mind, Vind, coeff, coeffInv,
    row, Nind, factor_colptr, vf,
    alpha1, alpha2)

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0

        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mrow += (coeff[j] * Mvar[j])
            Vrow += (coeff[j]^2 * Vvar[j])
        end
        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mfac[vf[j]] = alpha1[vf[j]] * ((Mind[row[j]] - Mrow) * coeffInv[j] + Mvar[j]) + alpha2[vf[j]] * Mfac[vf[j]]
            VfacInv[vf[j]] = (coeff[j]^2) / (Vind[row[j]] + Vrow - coeff[j]^2 * Vvar[j])
        end
    end

    return Mfac, VfacInv
end

function variable_to_factor(
    Mvar, Vvar, Mfac, VfacInv,
    Mdir, VdirInv,
    col, Nvariable, variable_colptr, fv)

    @inbounds for i = 1:Nvariable
        Mcol = 0.0
        VcolInv = 0.0

        for j in (variable_colptr[i]):variable_colptr[i+1]-1
            Mcol += Mfac[j] * VfacInv[j]
            VcolInv += VfacInv[j]
        end
        for j in (variable_colptr[i]):variable_colptr[i+1]-1
            Vvar[fv[j]] = 1 / (VcolInv + VdirInv[i] - VfacInv[j])
            Mvar[fv[j]] = (Mcol - Mfac[j] * VfacInv[j] + Mdir[i]) * Vvar[fv[j]]
        end
    end

    return Mvar, Vvar
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Compute messages using sum compensation, Kahan-Babuska algorithm
#------------------------------------------------------------------------------
function factor_to_variable_kahan(
    Mvar, Vvar, Mfac, VfacInv,
    Mind, Vind, coeff, coeffInv,
    row, Nind, factor_colptr, vf)

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0
        error = 0.0

        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mrow += (coeff[j] * Mvar[j])

            x = coeff[j]^2 * Vvar[j]
            t = Vrow + x
            if abs(Vrow) >= abs(x)
                error += (Vrow - t) + x
            else
                error += (x - t) + Vrow
            end
            Vrow = t
        end
        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mfac[vf[j]] = (Mind[row[j]] - Mrow) * coeffInv[j] + Mvar[j]
            VfacInv[vf[j]] = (coeff[j]^2) / (Vind[row[j]] + (Vrow - coeff[j]^2 * Vvar[j]) + error)
        end
    end

    return Mfac, VfacInv
end

function factor_to_variable_damp_kahan(
    Mvar, Vvar, Mfac, VfacInv,
    Mind, Vind, coeff, coeffInv,
    row, Nind, factor_colptr, vf,
    alpha1, alpha2)

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0
        error = 0.0

        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mrow += (coeff[j] * Mvar[j])

            x = coeff[j]^2 * Vvar[j]
            t = Vrow + x
            if abs(Vrow) >= abs(x)
                error += (Vrow - t) + x
            else
                error += (x - t) + Vrow
            end
            Vrow = t
        end
        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mfac[vf[j]] = alpha1[vf[j]] * ((Mind[row[j]] - Mrow) * coeffInv[j] + Mvar[j]) + alpha2[vf[j]] * Mfac[vf[j]]
            VfacInv[vf[j]] = (coeff[j]^2) / (Vind[row[j]] + (Vrow - coeff[j]^2 * Vvar[j]) + error)
        end
    end

    return Mfac, VfacInv
end

function variable_to_factor_kahan(
    Mvar, Vvar, Mfac, VfacInv,
    Mdir, VdirInv,
    col, Nvariable, variable_colptr, fv)

    @inbounds for i = 1:Nvariable
        Mcol = 0.0
        VcolInv = 0.0
        error = 0.0

        for j in (variable_colptr[i]):variable_colptr[i+1]-1
            Mcol += Mfac[j] * VfacInv[j]

            t = VcolInv + VfacInv[j]
            if abs(VcolInv) >= abs(VfacInv[j])
                error += (VcolInv - t) + VfacInv[j]
            else
                error += (VfacInv[j] - t) + VcolInv
            end
            VcolInv = t
        end
        for j in (variable_colptr[i]):variable_colptr[i+1]-1
            Vvar[fv[j]] = 1 / ((VcolInv - VfacInv[j]) + error + VdirInv[i])
            Mvar[fv[j]] = (Mcol - Mfac[j] * VfacInv[j] + Mdir[i]) * Vvar[fv[j]]
        end
    end

    return Mvar, Vvar
end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Compute messages using recursion with simply summation
#-------------------------------------------------------------------------------
function factor_recursion(
    Mvar, Vvar, Mfac, VfacInv,
    Mdir, VdirInv, Mind, Vind,
    coeff, coeffInv, row, col, Mcol, VcolInv,
    Nind, Nlink, factor_colptr, vf)

    @inbounds for i = 1:Nlink
        Mcol[col[i]] += Mfac[vf[i]] * VfacInv[vf[i]]
        VcolInv[col[i]] += VfacInv[vf[i]]
    end

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0

        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Vvar[j] = 1 / (VcolInv[col[j]] + VdirInv[col[j]] - VfacInv[vf[j]])
            Mvar[j] = (Mcol[col[j]] - Mfac[vf[j]] * VfacInv[vf[j]] + Mdir[col[j]]) * Vvar[j]

            Mrow += (coeff[j] * Mvar[j])
            Vrow += (coeff[j]^2 * Vvar[j])
        end
        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mfac[vf[j]] = (Mind[row[j]] - Mrow) * coeffInv[j] + Mvar[j]
            VfacInv[vf[j]] = (coeff[j]^2) / (Vind[row[j]] + Vrow - coeff[j]^2 * Vvar[j])
        end
    end

    fill!(Mcol, 0.0)
    fill!(VcolInv, 0.0)

    return Mfac, VfacInv, Mcol, VcolInv
end

function factor_recursion_damp(
    Mvar, Vvar, Mfac, VfacInv,
    Mdir, VdirInv, Mind, Vind,
    coeff, coeffInv, row, col, Mcol, VcolInv,
    Nind, Nlink, factor_colptr, vf, alpha1, alpha2)

    @inbounds for i = 1:Nlink
        Mcol[col[i]] += Mfac[vf[i]] * VfacInv[vf[i]]
        VcolInv[col[i]] += VfacInv[vf[i]]
    end

    @inbounds for i = 1:Nind
        Mrow = 0.0
        Vrow = 0.0

        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Vvar[j] = 1 / (VcolInv[col[j]] + VdirInv[col[j]] - VfacInv[vf[j]])
            Mvar[j] = (Mcol[col[j]] - Mfac[vf[j]] * VfacInv[vf[j]] + Mdir[col[j]]) * Vvar[j]

            Mrow += (coeff[j] * Mvar[j])
            Vrow += (coeff[j]^2 * Vvar[j])
        end
        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mfac[vf[j]] = alpha1[vf[j]] * ((Mind[row[j]] - Mrow) * coeffInv[j] + Mvar[j]) + alpha2[vf[j]] * Mfac[vf[j]]
            VfacInv[vf[j]] = (coeff[j]^2) / (Vind[row[j]] + Vrow - coeff[j]^2 * Vvar[j])
        end
    end

    fill!(Mcol, 0.0)
    fill!(VcolInv, 0.0)

    return Mfac, VfacInv, Mcol, VcolInv
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Compute marginals
#-------------------------------------------------------------------------------
function marginal(Mfac, VfacInv, Mdir, VdirInv, col, Nvariable, variable_colptr)
    Xbp = fill(0.0, Nvariable)

    @inbounds for i = 1:Nvariable
        Mcol = 0.0
        VcolInv = 0.0

        for j in (variable_colptr[i]):variable_colptr[i+1]-1
            Mcol += Mfac[j] * VfacInv[j]
            VcolInv += VfacInv[j]
        end
        Xbp[i] = (Mcol + Mdir[i]) / (VcolInv + VdirInv[i])
    end

    return Xbp
end
#-------------------------------------------------------------------------------
