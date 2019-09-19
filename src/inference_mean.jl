##########################################
#  Compute mean values of messages only  #
##########################################


#-------------------------------------------------------------------------------
# Compute means using state-of-the-art equations
#-------------------------------------------------------------------------------
function factor_to_variable_mean(
   Mvar_fac, Mfac_var, Mind, coeff, coeffInv,
   row, Nind, row_colptr, to_var)

   @inbounds for i = 1:Nind
       Mrow = 0.0
       Vrow = 0.0

       for j in (row_colptr[i]):row_colptr[i+1]-1
           Mrow += (coeff[j] * Mvar_fac[j])
       end
       for j in (row_colptr[i]):row_colptr[i+1]-1
           Mfac_var[to_var[j]] = (Mind[row[j]] - Mrow) * coeffInv[j] + Mvar_fac[j]
       end
   end

   return Mfac_var
end

function factor_to_variable_damp_mean(
    Mvar_fac, Mfac_var, Mind, coeff, coeffInv,
    row, Nind, row_colptr, to_var,
    alpha1, alpha2)

   @inbounds for i = 1:Nind
       Mrow = 0.0

       for j in (row_colptr[i]):row_colptr[i+1]-1
           Mrow += (coeff[j] * Mvar_fac[j])
       end
       for j in (row_colptr[i]):row_colptr[i+1]-1
           Mfac_var[to_var[j]] = alpha1[j] * ((Mind[row[j]] - Mrow) * coeffInv[j] + Mvar_fac[j]) + alpha2[j] * Mfac_var[to_var[j]]
       end
   end

   return Mfac_var
end

function variable_to_factor_mean(
    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var, Mdir,
    col, Nvar, col_colptr, to_fac)

    @inbounds for i = 1:Nvar
        Mcol = 0.0

        for j in (col_colptr[i]):col_colptr[i+1]-1
            Mcol += Mfac_var[j] * Wfac_var[j]
        end
        for j in (col_colptr[i]):col_colptr[i+1]-1
            Mvar_fac[to_fac[j]] = (Mcol - Mfac_var[j] * Wfac_var[j] + Mdir[i]) * Vvar_fac[to_fac[j]]
        end
    end
    return Mvar_fac
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Compute means using recursion equations
#-------------------------------------------------------------------------------
function factor_recursion_mean(
    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var,
    Mdir, Mind,
    coeff, coeffInv, row, col, Mcol, Wcol,
    Nind, Nlink, row_colptr, to_var)

    @inbounds for i = 1:Nlink
        Mcol[col[i]] += Mfac_var[to_var[i]] * Wfac_var[to_var[i]]
    end

    @inbounds for i = 1:Nind
        Mrow = 0.0

        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mvar_fac[j] = (Mcol[col[j]] - Mfac_var[to_var[j]] * Wfac_var[to_var[j]] + Mdir[col[j]]) * Vvar_fac[j]

            Mrow += (coeff[j] * Mvar_fac[j])
        end
        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mfac_var[to_var[j]] = (Mind[row[j]] - Mrow) * coeffInv[j] + Mvar_fac[j]
        end
    end

    fill!(Mcol, 0.0)

    return Mfac_var, Mcol
end

function factor_recursion_damp_mean(
    Mvar_fac, Vvar_fac, Mfac_var, Wfac_var,
    Mdir, Mind,
    coeff, coeffInv, row, col, Mcol, Wcol,
    Nind, Nlink, row_colptr, to_var, alpha1, alpha2)

    @inbounds for i = 1:Nlink
        Mcol[col[i]] += Mfac_var[to_var[i]] * Wfac_var[to_var[i]]
    end

    @inbounds for i = 1:Nind
        Mrow = 0.0

        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mvar_fac[j] = (Mcol[col[j]] - Mfac_var[to_var[j]] * Wfac_var[to_var[j]] + Mdir[col[j]]) * Vvar_fac[j]

            Mrow += (coeff[j] * Mvar_fac[j])
        end
        for j in (row_colptr[i]):row_colptr[i+1]-1
            Mfac_var[to_var[j]] = alpha1[j] * ((Mind[row[j]] - Mrow) * coeffInv[j] + Mvar_fac[j]) + alpha2[j] * Mfac_var[to_var[j]]
        end
    end

    fill!(Mcol, 0.0)

    return Mfac_var, Mcol
end
#-------------------------------------------------------------------------------
