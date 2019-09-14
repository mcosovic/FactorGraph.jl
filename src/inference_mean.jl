####################################
#  Compute messages and marginals  #
####################################


#-------------------------------------------------------------------------------
# Compute only means from factor nodes to variable nodes using
# state-of-the-art equations, or damping equations, or using sum
# compensation without or with damping
#-------------------------------------------------------------------------------
function factor_to_variable_mean(Mvar, Mfac, Mrow, Mind, coeffInv, row, Nlink)
    @inbounds for i = 1:Nlink
        Mfac[i] = (Mind[row[i]] - Mrow[row[i]]) * coeffInv[i] + Mvar[i]
    end

    fill!(Mrow, 0.0)

    return Mfac, Mrow
end

function factor_to_variable_mean_damp(
    Mvar, Mfac, Mrow, Mind, coeffInv, row, Nlink,
    alpha1, alpha2)

    @inbounds for i = 1:Nlink
        Mfac[i] = alpha1[i] * ((Mind[row[i]] - Mrow[row[i]]) * coeffInv[i] +
                  Mvar[i]) + alpha2[i] * Mfac[i]
    end

    fill!(Mrow, 0.0)

    return Mfac, Mrow
end

function factor_mean_recursion(
    Mfac, Mrow, Mcol, Maux, Mind,
    coeffInv, row, Nlink)

   @inbounds for i = 1:Nlink
       Mfac[i] = (Mind[row[i]] - Mrow[row[i]]) * coeffInv[i] + Maux[i]
   end

   fill!(Mrow, 0.0)
   fill!(Mcol, 0.0)

   return Mfac, Mrow, Mcol
end

function factor_recursion_damp(
    Mfac, Mrow, Mcol, Maux, Mind,
    coeffInv, row, Nlink, alpha1, alpha2)

    @inbounds for i = 1:Nlink
        Mfac[i] = alpha1[i] * ((Mind[row[i]] - Mrow[row[i]]) * coeffInv[i] +
                  Maux[i]) + alpha2[i] * Mfac[i]
    end

    fill!(Mrow, 0.0)
    fill!(Mcol, 0.0)

    return Mfac, Mrow, Mcol
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Compute only means from variable nodes to factor nodes using
# state-of-the-art equations, or sum compensation
#-------------------------------------------------------------------------------
function variable_to_factor_mean(Mvar, Vvar, Mfac, VfacInv, Mcol, Mdir, col, Nlink)

    @inbounds for i = 1:Nlink
        Mvar[i] = (Mcol[col[i]] - Mfac[i] * VfacInv[i] + Mdir[col[i]]) * Vvar[i]
    end

    fill!(Mcol, 0.0)

    return Mvar, Mcol
end
#-------------------------------------------------------------------------------
