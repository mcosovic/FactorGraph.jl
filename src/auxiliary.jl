########################
#  Auxiliary functions #
########################


#-------------------------------------------------------------------------------
# Load arrays for messages
#-------------------------------------------------------------------------------
function load_messages(coeff)
    Mfac_var = similar(coeff)
    Wfac_var = similar(coeff)
    Mvar_fac = similar(coeff)
    Vvar_fac = similar(coeff)

    return Mfac_var, Wfac_var, Mvar_fac, Vvar_fac
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Load arrays for column summation, METHOD = "recursion"
#-------------------------------------------------------------------------------
function load_sum_col_recursion(Nvar)
    Mcol = fill(0.0, Nvar)
    Wcol = fill(0.0, Nvar)

    return Mcol, Wcol
end
#-------------------------------------------------------------------------------
