function wlsMldivide(H, b, v)
    W = spdiagm(0 =>  @. 1.0 / sqrt(v))
    Hti = W * H
    rti = W * b

    xml = (Hti' * Hti) \ (Hti' * rti)

    return xml
end

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

function errors(H, b, v, xbp, ERROR)
    if ERROR == "on"
        wls = @elapsed begin
            xwls = wlsMldivide(H, b, v)
        end

        wrss_bp, wrss_wls = wrss(H, b, v, xbp, xwls)
        rmse_bp, rmse_wls = rmse(H, b, v, xbp, xwls)

        wls_vs_bp(wls, wrss_wls, wrss_bp, rmse_wls, rmse_bp, xbp, xwls)
    end
end

function bp_time(factorgraph, initialize, inference, solution)
    print("Preprocesing: $(@sprintf("%.6f", factorgraph*1000)) ms \n")
    print("Initialize:   $(@sprintf("%.6f", initialize*1000)) ms \n")
    print("Inference:    $(@sprintf("%.6f", inference*1000)) ms \n")
    print("Marginal:     $(@sprintf("%.6f", solution*1000)) ms \n")
    println(" ")
    print("BP:           $(@sprintf("%.6f", (factorgraph+initialize+inference+solution)*1000)) ms \n")
end

function wls_vs_bp(wls, wrss_wls, wrss_bp, rmse_wls, rmse_bp, xbp, xwls)
    print("WLS:          $(@sprintf("%.6f", wls*1000)) ms \n")
    println(" ")
    print("WRSS WLS:     $(@sprintf("%.6f", wrss_wls)) \n")
    print("WRSS BP:      $(@sprintf("%.6f", wrss_bp)) \n")
    println(" ")
    print("RMSE WLS:     $(@sprintf("%.6f", rmse_wls)) \n")
    print("RMSE BP:      $(@sprintf("%.6f", rmse_bp)) \n")
    println(" ")
    println("MAX DIFF: ", sort(abs.(xbp-xwls))[end])
end
