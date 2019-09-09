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

function errors(H, b, v, xbp)
    if ERROR == "on"
        wls = @elapsed begin
            xwls = wlsMldivide(H, b, v)
        end

        wrss_bp, wrss_wls = wrss(H, b, v, xbp, xwls)
        rmse_bp, rmse_wls = rmse(H, b, v, xbp, xwls)

        wls_vs_bp(wls, wrss_wls, wrss_bp, rmse_wls, rmse_bp, xbp, xwls)
    end
end
