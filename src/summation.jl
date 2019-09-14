##########################
#  Summarizing messages  #
##########################


#-------------------------------------------------------------------------------
# Summarizing messages from variable nodes to factor nodes using simply
# summation, or using Kahan-Babuska algorithm
#-------------------------------------------------------------------------------
function sum_rows(
    Mvar, Vvar, Mrow, Vrow,
    coeff, row, Nlink)

    @inbounds for i = 1:Nlink
        Mrow[row[i]] += (coeff[i] * Mvar[i])
        Vrow[row[i]] += (coeff[i]^2 * Vvar[i])
    end

    return Mrow, Vrow
end

function sum_rows_kahan(
    Mvar, Vvar, Mrow, Vrow, error_row,
    coeff, row, Nlink)

    @inbounds for i = 1:Nlink
        Mrow[row[i]] += (coeff[i] * Mvar[i])

        x = coeff[i]^2 * Vvar[i]
        t = Vrow[row[i]] + x
        if abs(Vrow[row[i]]) >= abs(x)
            error_row[row[i]] += (Vrow[row[i]] - t) + x
        else
            error_row[row[i]] += (x - t) + Vrow[row[i]]
        end
        Vrow[row[i]] = t
    end

   return Mrow, Vrow, error_row
end

function sum_rows_recursion(
    Mfac, VfacInv, Mrow, Vrow, Mcol, VcolInv, Maux, Vaux,
    Mdir, VdirInv,
    coeff, row, col, Nlink)

    @inbounds for i = 1:Nlink
        Vaux[i] = 1 / (VcolInv[col[i]] + VdirInv[col[i]] - VfacInv[i])
        Vrow[row[i]] += (coeff[i]^2 * Vaux[i])

        Maux[i] = (Mcol[col[i]] - Mfac[i] * VfacInv[i] + Mdir[col[i]]) * Vaux[i]
        Mrow[row[i]] += (coeff[i] * Maux[i])
    end

    return Mrow, Vrow, Vaux, Maux
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Summarizing messages from factor nodes to variable nodes using simply
# summation, or using Kahan-Babuska algorithm
#-------------------------------------------------------------------------------
function sum_cols(
    Mfac, VfacInv, Mcol, VcolInv,
    col, Nlink)

    @inbounds for i = 1:Nlink
        Mcol[col[i]] += Mfac[i] * VfacInv[i]
        VcolInv[col[i]] += VfacInv[i]
    end

    return Mcol, VcolInv
end

function sum_cols_kahan(
    Mfac, VfacInv, Mcol, VcolInv, error_col,
    col, Nlink)

    @inbounds for i = 1:Nlink
        Mcol[col[i]] += Mfac[i] * VfacInv[i]

        t = VcolInv[col[i]] + VfacInv[i]
        if abs(VcolInv[col[i]]) >= abs(VfacInv[i])
            error_col[col[i]] += (VcolInv[col[i]] - t) + VfacInv[i]
        else
            error_col[col[i]] += (VfacInv[i] - t) + VcolInv[col[i]]
        end
        VcolInv[col[i]] = t
    end

   return Mcol, VcolInv, error_col
end
#-------------------------------------------------------------------------------
