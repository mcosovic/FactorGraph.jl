########################
#  Auxiliary functions #
########################


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
# Load arrays for column summation, METHOD = "recursion"
#-------------------------------------------------------------------------------
function load_sum_col_recursion(Nvariable)
    Mcol = fill(0.0, Nvariable)
    VcolInv = fill(0.0, Nvariable)

    return Mcol, VcolInv
end
#-------------------------------------------------------------------------------
