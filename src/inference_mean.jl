##########################################
#  Compute mean values of messages only  #
##########################################


#-------------------------------------------------------------------------------
# Compute means using state-of-the-art equations
#-------------------------------------------------------------------------------
function factor_to_variable_mean(
   Mvar, Mfac, Mind, coeff, coeffInv,
   row, Nind, factor_colptr, vf)

   @inbounds for i = 1:Nind
       Mrow = 0.0
       Vrow = 0.0

       for j in (factor_colptr[i]):factor_colptr[i+1]-1
           Mrow += (coeff[j] * Mvar[j])
       end
       for j in (factor_colptr[i]):factor_colptr[i+1]-1
           Mfac[vf[j]] = (Mind[row[j]] - Mrow) * coeffInv[j] + Mvar[j]
       end
   end

   return Mfac
end

function factor_to_variable_damp_mean(
    Mvar, Mfac, Mind, coeff, coeffInv,
    row, Nind, factor_colptr, vf,
    alpha1, alpha2)

   @inbounds for i = 1:Nind
       Mrow = 0.0

       for j in (factor_colptr[i]):factor_colptr[i+1]-1
           Mrow += (coeff[j] * Mvar[j])
       end
       for j in (factor_colptr[i]):factor_colptr[i+1]-1
           Mfac[vf[j]] = alpha1[vf[j]] * ((Mind[row[j]] - Mrow) * coeffInv[j] + Mvar[j]) + alpha2[vf[j]] * Mfac[vf[j]]
       end
   end

   return Mfac
end

function variable_to_factor_mean(
    Mvar, Vvar, Mfac, VfacInv, Mdir,
    col, Nvariable, variable_colptr, fv)

    @inbounds for i = 1:Nvariable
        Mcol = 0.0

        for j in (variable_colptr[i]):variable_colptr[i+1]-1
            Mcol += Mfac[j] * VfacInv[j]
        end
        for j in (variable_colptr[i]):variable_colptr[i+1]-1
            Mvar[fv[j]] = (Mcol - Mfac[j] * VfacInv[j] + Mdir[i]) * Vvar[fv[j]]
        end
    end
    return Mvar
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Compute means using recursion equations
#-------------------------------------------------------------------------------
function factor_recursion_mean(
    Mvar, Vvar, Mfac, VfacInv,
    Mdir, Mind,
    coeff, coeffInv, row, col, Mcol, VcolInv,
    Nind, Nlink, factor_colptr, vf)

    @inbounds for i = 1:Nlink
        Mcol[col[i]] += Mfac[vf[i]] * VfacInv[vf[i]]
    end

    @inbounds for i = 1:Nind
        Mrow = 0.0

        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mvar[j] = (Mcol[col[j]] - Mfac[vf[j]] * VfacInv[vf[j]] + Mdir[col[j]]) * Vvar[j]

            Mrow += (coeff[j] * Mvar[j])
        end
        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mfac[vf[j]] = (Mind[row[j]] - Mrow) * coeffInv[j] + Mvar[j]
        end
    end

    fill!(Mcol, 0.0)

    return Mfac, Mcol
end

function factor_recursion_damp_mean(
    Mvar, Vvar, Mfac, VfacInv,
    Mdir, Mind,
    coeff, coeffInv, row, col, Mcol, VcolInv,
    Nind, Nlink, factor_colptr, vf, alpha1, alpha2)

    @inbounds for i = 1:Nlink
        Mcol[col[i]] += Mfac[vf[i]] * VfacInv[vf[i]]
    end

    @inbounds for i = 1:Nind
        Mrow = 0.0

        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mvar[j] = (Mcol[col[j]] - Mfac[vf[j]] * VfacInv[vf[j]] + Mdir[col[j]]) * Vvar[j]

            Mrow += (coeff[j] * Mvar[j])
        end
        for j in (factor_colptr[i]):factor_colptr[i+1]-1
            Mfac[vf[j]] = alpha1[vf[j]] * ((Mind[row[j]] - Mrow) * coeffInv[j] + Mvar[j]) + alpha2[vf[j]] * Mfac[vf[j]]
        end
    end

    fill!(Mcol, 0.0)

    return Mfac, Mcol
end
#-------------------------------------------------------------------------------
