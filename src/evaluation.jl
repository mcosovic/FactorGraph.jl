########################################################################
#  Compare the accuracy of the belief propagation algorithm to that    #
#  of the weighted-least squares method and print results in the REPL  #
########################################################################


#------------------------------------------------------------------------
# Compute the weighted-least squares solution
#------------------------------------------------------------------------
function wlsMldivide(H, b, v)
    W = spdiagm(0 =>  @. 1.0 / sqrt(v))
    Hti = W * H
    rti = W * b

    xml = (Hti' * Hti) \ (Hti' * rti)

    return xml
end
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# Compute weighted residual sum of squares (WRSS) and root mean square
# error (RMSE) metrics
#------------------------------------------------------------------------
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
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# The belief propagation algorithm is evaluated for the option
# ERROR == "on" and results are shown in the REPL
#------------------------------------------------------------------------
function errors(H, b, v, xbp, ERROR, TIME)
    if ERROR == "on"
        wls = @elapsed begin
            xwls = wlsMldivide(H, b, v)
        end

        wrss_bp, wrss_wls = wrss(H, b, v, xbp, xwls)
        rmse_bp, rmse_wls = rmse(H, b, v, xbp, xwls)

        if TIME == "on"
            pretty_table(["WLS" 1000 * wls],
                         ["" "Time (ms)"];
                         screen_size = (-1,-1),
                         alignment=[:l,:r],
                         formatter = ft_printf("%3.6f", [2]))
        end

        col0 = ["WRSS", "RMSE"]
        col1 = [wrss_bp, rmse_bp]
        col2 = [wrss_wls, rmse_wls]
        data = [col1 col1 col2 abs.(col1 - col2)]

        if data[1,4] > 1
            high = Highlighter((data,i,j) -> j in (4) && data[i,j] > 1e-2,
                   Crayon(bold = true, background = :red))
        else
            high = ()
        end
        pretty_table([col0 col1 col2 abs.(col1 - col2)],
                     ["Error" "BP" "WLS" "Distance"];
                     screen_size = (-1,-1),
                     alignment=[:r,:r,:r, :r],
                     formatter = ft_printf(["%3.6f","%3.6f","%3.6e"], [2,3,4]),
                     highlighters = high)

        A = [collect(1:length(xbp)) xbp xwls abs.(xbp - xwls)]
        pretty_table(A[reverse(sortperm(A[:, 4])),  :],
                     ["State Variable" "BP Estimate" "WLS Estimate" "Max to Min Distance"],
                     alignment=[:r,:r,:r, :r],
                     formatter = ft_printf(["%3.0f", "%3.6f","%3.6f","%3.6e"], [1,2,3,4]))
    end
end
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# The time evolution of the belief propagation algorithm for the option
# TIME = "on"
#------------------------------------------------------------------------
function bp_time(fg, it, ic, so, TIME)
    if TIME == "on"
        col1 = ["Preprocesing", "Initialize", "Inference", "Marginal", "Total"]
        col2 = 1000 .* [fg, it, ic, so, fg + it + ic + so]

        pretty_table([col1 col2], ["BP Phase" "Time (ms)"];
                     screen_size = (-1,200),
                     alignment=[:l,:r],
                     formatter = ft_printf("%3.6f", [2]),
                     hlines = [4])
    end
end
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# Graph Data
#------------------------------------------------------------------------
function graph_statistic(Nf, Nv, Ndl, Nli, vir, v)
    col1 = ["Number of Factor Nodes",
            "Number of Virtual Nodes",
            "Number of Variable Nodes",
            "Number of Direct Links",
            "Number of Indirect Links",
            "Minimum Variance Value",
            "Maximum Variance Value"]

    col2 = [Nf, length(findall(!iszero, vir)), Nv, Ndl, Nli, minimum(v), maximum(v)]

    pretty_table([col1 col2], ["Graph Data" ""];
                 screen_size = (-1,200),
                 alignment=[:l,:r],
                 hlines = [3,5])
end
#------------------------------------------------------------------------
