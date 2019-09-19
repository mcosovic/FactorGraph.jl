########################################################################
#  Compare the accuracy of the belief propagation algorithm to that    #
#  of the weighted-least squares method and print results in the REPL  #
########################################################################


#-------------------------------------------------------------------------------
# Compute the weighted-least squares solution
#-------------------------------------------------------------------------------
function wls_mldivide(jacobian, observation, noise)
    W = spdiagm(0 =>  @. 1.0 / sqrt(noise))
    Hti = W * jacobian
    rti = W * observation

    Xml = (Hti' * Hti) \ (Hti' * rti)

    return Xml
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Compute weighted residual sum of squares (WRSS) and root mean square
# error (RMSE) metrics
#-------------------------------------------------------------------------------
function wrss(jacobian, observation, noise, Xbp, Xwls)
    observation_estimate_bp = jacobian * Xbp
    observation_estimate_wls = jacobian * Xwls

    wrss_bp = sum(((observation - observation_estimate_bp).^2) ./ noise)
    wrss_wls = sum(((observation - observation_estimate_wls).^2) ./ noise)

    return wrss_bp, wrss_wls
end

function rmse(jacobian, observation, noise, Xbp, Xwls)
    observation_estimate_bp = jacobian * Xbp
    observation_estimate_wls = jacobian * Xwls
    N, ~ = size(jacobian)

    rmse_bp = ((sum(observation - observation_estimate_bp).^2) / N)^(1/2)
    rmse_wls = ((sum(observation - observation_estimate_wls).^2) / N)^(1/2)

    return rmse_bp, rmse_wls
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# The belief propagation algorithm is evaluated for the option
# ERROR == "on" and results are shown in the REPL
#-------------------------------------------------------------------------------
function errors(jacobian, observation, noise, Xbp, TIME)
    wls = @elapsed begin
        Xwls = wls_mldivide(jacobian, observation, noise)
    end

    wrss_bp, wrss_wls = wrss(jacobian, observation, noise, Xbp, Xwls)
    rmse_bp, rmse_wls = rmse(jacobian, observation, noise, Xbp, Xwls)

    if TIME == "on"
        pretty_table(
            ["WLS" 1000 * wls],
            ["" "Time (ms)"];
            screen_size = (-1,-1),
            alignment=[:l,:r],
            formatter = ft_printf("%3.6f", [2]))
    end

    col0 = ["WRSS", "RMSE"]
    col1 = [wrss_bp, rmse_bp]
    col2 = [wrss_wls, rmse_wls]
    data = [col1 col1 col2 abs.(col1 ./ col2)]

    high = Highlighter((data,i,j) -> i == 1 && data[1,4] > 1.1,
           Crayon(bold = true, background = :red))

    pretty_table(
        [col0 data[:,2:4]],
        ["Error" "BP" "WLS" "Ratio (BP/WLS)"];
        screen_size = (-1,-1),
        alignment=[:r,:r,:r, :r],
        formatter = ft_printf(["%3.6f","%3.6f","%3.6e"], [2,3,4]),
        highlighters = high)

    A = [collect(1:length(Xbp)) Xbp Xwls Xbp ./ Xwls]
    pretty_table(
        A[reverse(sortperm(A[:, 4])),  :],
        ["State Variable" "BP Estimate" "WLS Estimate" "Max to Min Ratio"],
        alignment=[:r,:r,:r, :r],
        formatter = ft_printf(["%3.0f", "%3.6f","%3.6f","%3.6e"], [1,2,3,4]))
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# The time evolution of the belief propagation algorithm for the option
# TIME = "on"
#-------------------------------------------------------------------------------
function bp_time(factorgraph, initialize, inference, solution)
    data = ["Preprocesing" 1000 * factorgraph;
            "Initialize"  1000 * initialize;
            "Inference"  1000 * inference;
            "Marginal"  1000 * solution;
            "Total" 1000 * (factorgraph + initialize + inference + solution)]

    pretty_table(
        data, ["BP Phase" "Time (ms)"];
        screen_size = (-1,-1),
        alignment=[:l,:r],
        formatter = ft_printf("%3.6f", [2]),
        hlines = [4])
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Graph Data
#-------------------------------------------------------------------------------
function graph_statistic(Nfac, Nvar, Ndir, Nlink, virtual, noise)
    data = ["Number of Factor Nodes" Nfac;
            "Number of Virtual Nodes" length(findall(!iszero, virtual));
            "Number of Variable Nodes" Nvar;
            "Number of Direct Links" Ndir;
            "Number of Indirect Links" Nlink;
            "Minimum Variance Value" minimum(noise);
            "Maximum Variance Value" maximum(noise)]

    pretty_table(
        data, ["Graph Data" ""];
        screen_size = (-1,-1),
        alignment=[:l,:r],
        hlines = [3,5])
end
#-------------------------------------------------------------------------------
