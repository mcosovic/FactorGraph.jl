################################################################################
# Compare the accuracy of the belief propagation algorithm to that of the
# weighted-least squares method and print results in the REPL
################################################################################


################################################################################
# Compute the weighted-least squares solution
#-------------------------------------------------------------------------------
function wlsMldivide(H, b, v)
    W = spdiagm(0 =>  @. 1.0 / sqrt(v))
    Hti = W * H
    rti = W * b

    xml = (Hti' * Hti) \ (Hti' * rti)

    return xml
end
################################################################################



################################################################################
# Compute weighted residual sum of squares (WRSS) and  root mean square error
# (RMSE) metrics
#-------------------------------------------------------------------------------
function wrss(H, b, v, xbp, xwls)
    fbp = H * xbp
    fwls = H * xwls

    wrss_bp = sum(((b - fbp).^2) ./ v)
    wrss_wls = sum(((b - fwls).^2) ./ v)

    return wrss_bp, wrss_wls
end

function rmse(H, b, v, xbp, xwls)
    fbp = H * xbp
    fwls = H * xwls
    Nf, ~ = size(H)

    rmse_bp = ((sum(b - fbp).^2) / Nf)^(1/2)
    rmse_wls = ((sum(b - fwls).^2) / Nf)^(1/2)

    return rmse_bp, rmse_wls
end
################################################################################


################################################################################
# The belief propagation algorithm is evaluated for the option ERROR == "on"
# and results are shown in the REPL
#-------------------------------------------------------------------------------
function errors(H, b, v, xbp, ERROR)
    if ERROR == "on"
        wls = @elapsed begin
            xwls = wlsMldivide(H, b, v)
        end

        wrss_bp, wrss_wls = wrss(H, b, v, xbp, xwls)
        rmse_bp, rmse_wls = rmse(H, b, v, xbp, xwls)

        print("WLS:          $(@sprintf("%.6f", wls * 1000)) ms \n")
        println(" ")
        print("WRSS WLS:     $(@sprintf("%.6f", wrss_wls)) \n")
        print("WRSS BP:      $(@sprintf("%.6f", wrss_bp)) \n")
        println(" ")
        print("RMSE WLS:     $(@sprintf("%.6f", rmse_wls)) \n")
        print("RMSE BP:      $(@sprintf("%.6f", rmse_bp)) \n")
        println(" ")
        println("MAX DIFF: ", sort(abs.(xbp-xwls))[end])
    end
end
################################################################################


################################################################################
# The time evolution of the belief propagation algorithm for the option
# TIME = "on"
#-------------------------------------------------------------------------------
function bp_time(fg, it, ic, so)
    print("Preprocesing: $(@sprintf("%.6f", fg * 1000)) ms \n")
    print("Initialize:   $(@sprintf("%.6f", it * 1000)) ms \n")
    print("Inference:    $(@sprintf("%.6f", ic * 1000)) ms \n")
    print("Marginal:     $(@sprintf("%.6f", so * 1000)) ms \n")
    println(" ")
    print("BP:           $(@sprintf("%.6f", (fg + it + ic + so) * 1000)) ms \n")
end
################################################################################
