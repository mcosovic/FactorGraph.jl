##########################
#  Summarizing messages  #
##########################


#-------------------------------------------------------------------------------
# Summarizing messages from variable nodes to factor nodes using simply
# summation, or using Kahan-Babuska algorithm
#-------------------------------------------------------------------------------
function sum_rows_mean(Mvar, Mrow, coeff, row, Nlink)

    @inbounds for i = 1:Nlink
        Mrow[row[i]] += (coeff[i] * Mvar[i])
    end

    return Mrow
end

function sum_rows_mean_recursion(
    Mfac, VfacInv, Mrow, Mcol, Maux, Vaux, Mdir,
    coeff, row, col, Nlink)

    @inbounds for i = 1:Nlink
        Maux[i] = (Mcol[col[i]] - Mfac[i] * VfacInv[i] + Mdir[col[i]]) * Vaux[i]
        Mrow[row[i]] += (coeff[i] * Maux[i])
    end

    return Mrow, Maux
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Summarizing messages from factor nodes to variable nodes using simply
# summation, or using Kahan-Babuska algorithm
#-------------------------------------------------------------------------------
function sum_cols_mean(Mfac, VfacInv, Mcol, col, Nlink)

    @inbounds for i = 1:Nlink
        Mcol[col[i]] += Mfac[i] * VfacInv[i]
    end

    return Mcol
end
#-------------------------------------------------------------------------------
