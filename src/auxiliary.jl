#########################################################
#  Auxiliary functions for arrays loading and clearing  #
#########################################################


#-------------------------------------------------------------------------------
# Load arrays for messages
#-------------------------------------------------------------------------------
function load_messages(coeff)
    Mfac = similar(coeff)
    VfacInv = similar(coeff)
    Mvar = similar(coeff)
    Vvar = similar(coeff)

    return Mfac, VfacInv, Mvar, Vvar
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Load arrays for message summation in METHOD = "passing" and "recursion" and
# ALGORITHM = "sum" and "kahan"
#-------------------------------------------------------------------------------
function load_sum(Nvariable, Nind)
    Mrow = fill(0.0, Nind)
    Vrow = fill(0.0, Nind)
    Mcol = fill(0.0, Nvariable)
    VcolInv = fill(0.0, Nvariable)

    return Mrow, Vrow, Mcol, VcolInv
end

function load_sum_recursion(Nvariable, Nind, Nlink)
    Mrow = fill(0.0, Nind)
    Vrow = fill(0.0, Nind)
    Mcol = fill(0.0, Nvariable)
    VcolInv = fill(0.0, Nvariable)
    Vaux = fill(0.0, Nlink)
    Maux = similar(Vaux)

    return Mrow, Vrow, Mcol, VcolInv, Vaux, Maux
end

function load_sum_kahan(Nvariable, Nind)
    Mrow = fill(0.0, Nind)
    Vrow = fill(0.0, Nind)
    error_row = fill(0.0, Nind)
    Mcol = fill(0.0, Nvariable)
    VcolInv = fill(0.0, Nvariable)
    error_col = fill(0.0, Nvariable)

    return Mrow, Vrow, error_row, Mcol, VcolInv, error_col
end
#-------------------------------------------------------------------------------
