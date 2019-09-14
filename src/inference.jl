####################################
#  Compute messages and marginals  #
####################################


#-------------------------------------------------------------------------------
# Compute messages from factor nodes to variable nodes using
# state-of-the-art equations, or damping equations, or using sum
# compensation without or with damping
#-------------------------------------------------------------------------------
 function factor_to_variable(
    Mvar, Vvar, Mfac, VfacInv,
    Mrow, Vrow,
    Mind, Vind,
    coeff, coeffInv,
    row, Nlink)

    @inbounds for i = 1:Nlink
        Mfac[i] = (Mind[row[i]] - Mrow[row[i]]) * coeffInv[i] + Mvar[i]

        VfacInv[i] = (coeff[i]^2) / (Vind[row[i]] + Vrow[row[i]] -
                      coeff[i]^2 * Vvar[i])
    end

    fill!(Mrow, 0.0)
    fill!(Vrow, 0.0)

    return Mfac, VfacInv, Mrow, Vrow
 end

 function factor_to_variable_damp(
    Mvar, Vvar, Mfac, VfacInv,
    Mrow, Vrow,
    Mind, Vind,
    coeff, coeffInv,
    row, Nlink, alpha1, alpha2)

    @inbounds for i = 1:Nlink
        Mfac[i] = alpha1[i] * ((Mind[row[i]] - Mrow[row[i]]) * coeffInv[i] +
                  Mvar[i]) + alpha2[i] * Mfac[i]

        VfacInv[i] = (coeff[i]^2) / (Vind[row[i]] + Vrow[row[i]] -
                     coeff[i]^2 * Vvar[i])
    end

    fill!(Mrow, 0.0)
    fill!(Vrow, 0.0)

    return Mfac, VfacInv, Mrow, Vrow
 end


 function factor_to_variable_kahan(
    Mvar, Vvar, Mfac, VfacInv,
    Mrow, Vrow, error_row,
    Mind, Vind, coeff, coeffInv, row, Nlink)

    @inbounds for i = 1:Nlink
        Mfac[i] = (Mind[row[i]] - Mrow[row[i]]) * coeffInv[i] + Mvar[i]

        VfacInv[i] = (coeff[i]^2) / (Vind[row[i]] + (Vrow[row[i]] -
                     coeff[i]^2 * Vvar[i]) + error_row[row[i]])
    end

    fill!(Mrow, 0.0)
    fill!(Vrow, 0.0)
    fill!(error_row, 0.0)

    return Mfac, VfacInv, Mrow, Vrow, error_row
 end

 function factor_to_variable_kahan_damp(
    Mvar, Vvar, Mfac, VfacInv,
    Mrow, Vrow, error_row,
    Mind, Vind, coeff, coeffInv,
    row, Nlink, alpha1, alpha2)

    @inbounds for i = 1:Nlink
        Mfac[i] = alpha1[i] * ((Mind[row[i]] - Mrow[row[i]]) * coeffInv[i] +
                  Mvar[i]) + alpha2[i] * Mfac[i]

        VfacInv[i] = (coeff[i]^2) / (Vind[row[i]] + (Vrow[row[i]] -
                     coeff[i]^2 * Vvar[i]) + error_row[row[i]])
    end

    fill!(Mrow, 0.0)
    fill!(Vrow, 0.0)
    fill!(error_row, 0.0)

    return Mfac, VfacInv, Mrow, Vrow, error_row
 end

 function factor_recursion(
    Mvar, Vvar, Mfac, VfacInv,
    Mrow, Vrow, Mcol, VcolInv, Vaux, Maux,
    Mind, Vind,
    coeff, coeffInv,
    row, Nlink)

    @inbounds for i = 1:Nlink
        Mfac[i] = (Mind[row[i]] - Mrow[row[i]]) * coeffInv[i] + Maux[i]

        VfacInv[i] = (coeff[i]^2) / (Vind[row[i]] + Vrow[row[i]] -
                     coeff[i]^2 * Vaux[i])
    end

    fill!(Mrow, 0.0)
    fill!(Mcol, 0.0)
    fill!(Vrow, 0.0)
    fill!(VcolInv, 0.0)

   return Mfac, VfacInv, Mrow, Vrow, Mcol, VcolInv
end

function factor_recursion_damp(
   Mvar, Vvar, Mfac, VfacInv,
   Mrow, Vrow, Mcol, VcolInv, Vaux, Maux,
   Mind, Vind,
   coeff, coeffInv,
   row, Nlink, alpha1, alpha2)

    @inbounds for i = 1:Nlink
        Mfac[i] = alpha1[i] * ((Mind[row[i]] - Mrow[row[i]]) * coeffInv[i] +
                  Maux[i]) + alpha2[i] * Mfac[i]

        VfacInv[i] = (coeff[i]^2) / (Vind[row[i]] + Vrow[row[i]] -
                     coeff[i]^2 * Vaux[i])
    end

    fill!(Mrow, 0.0)
    fill!(Mcol, 0.0)
    fill!(Vrow, 0.0)
    fill!(VcolInv, 0.0)

  return Mfac, VfacInv, Mrow, Vrow, Mcol, VcolInv
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Compute messages from variable nodes to factor nodes using
# state-of-the-art equations, or sum compensation
#-------------------------------------------------------------------------------
function variable_to_factor(
    Mvar, Vvar, Mfac, VfacInv,
    Mcol, VcolInv,
    Mdir, VdirInv,
    col, Nlink)

    @inbounds for i = 1:Nlink
        Vvar[i] = 1 / (VcolInv[col[i]] + VdirInv[col[i]] - VfacInv[i])

        Mvar[i] = (Mcol[col[i]] - Mfac[i] * VfacInv[i] + Mdir[col[i]]) * Vvar[i]
    end

    fill!(Mcol, 0.0)
    fill!(VcolInv, 0.0)

    return Mvar, Vvar, Mcol, VcolInv
end

function variable_to_factor_kahan(
    Mvar, Vvar, Mfac, VfacInv,
    Mcol, VcolInv, error_col,
    Mdir, VdirInv,
    col, Nlink)

    @inbounds for i = 1:Nlink
        Vvar[i] = 1 / ((VcolInv[col[i]] - VfacInv[i]) +
                  error_col[col[i]] + VdirInv[col[i]])

        Mvar[i] = (Mcol[col[i]] - Mfac[i] * VfacInv[i] + Mdir[col[i]]) * Vvar[i]
    end

    fill!(Mcol, 0.0)
    fill!(VcolInv, 0.0)
    fill!(error_col, 0.0)

    return Mvar, Vvar, Mcol, VcolInv, error_col
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Compute marginals using state-of-the-art equations, or sum compensation
#-------------------------------------------------------------------------------
function marginal(Mcol, VcolInv, Mdir, VdirInv, col, Nvariable)
    Xbp = fill(0.0, Nvariable)

    @inbounds for i = 1:Nvariable
        Xbp[i] = (Mcol[i] + Mdir[i]) / (VcolInv[i] + VdirInv[i])
    end

    return Xbp
end

function marginal_kahan(Mcol, VcolInv, error_col, Mdir, VdirInv, col, Nvariable)
    Xbp = fill(0.0, Nvariable)

    @inbounds for i = 1:Nvariable
        Xbp[i] = (Mcol[i] + Mdir[i]) / (VcolInv[i] + VdirInv[i] + error_col[i])
    end

    return Xbp
end
#-------------------------------------------------------------------------------
