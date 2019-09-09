

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
    println("MAX DIFF: ", sort(abs.(xbp-xwls))[end-5:end])
end
