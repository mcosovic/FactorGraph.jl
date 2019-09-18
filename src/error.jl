##################
#  Check Inputs  #
##################


#-------------------------------------------------------------------------------
# Check MAXI, DAMP, BUMP
#-------------------------------------------------------------------------------
function check_iteration_scheme(MAXI, DAMP, BUMP)
    if MAXI < 0
        error("The upper limit on BP iterations MAXI has invalid value.")
    end
    if DAMP > MAXI || DAMP < 0
        error("Applied randomized damping parameter DAMP has invalid value.")
    end
    if BUMP > MAXI || BUMP <= 0
        error("Applied variance computation parameter BUMP has invalid value.")
    end
end
#-------------------------------------------------------------------------------
