### Efficient GBP with Kahan-Babuska algorithm: Factor to variable messages
@inline function kahan_factor_to_variable(graph, bp)
    @inbounds for i = 1:graph.Nind
        Mrow = 0.0; Vrow = 0.0; errorV = 0.0; errorM = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            summands = graph.coeff[j]^2 * bp.Vvar_fac[j]
            Vrow, errorV = kahan(summands, Vrow, errorV)

            summands = graph.coeff[j] * bp.Mvar_fac[j]
            Mrow, errorM = kahan(summands, Mrow, errorM)
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = (graph.Mind[j] - (Mrow - graph.coeff[j] * bp.Mvar_fac[j]) - errorM) * graph.coeffInv[j] 
            bp.Wfac_var[bp.to_var[j]] = (graph.coeff[j]^2) / (graph.Vind[j] + (Vrow - graph.coeff[j]^2 * bp.Vvar_fac[j]) + errorV)
        end
    end
end

### Efficient GBP with Kahan-Babuska algorithm: Factor to variable messages with damping
@inline function kahan_factor_to_variable_damp(graph, bp)
    @inbounds for i = 1:graph.Nind
        Mrow = 0.0; Vrow = 0.0; errorV = 0.0; errorM = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            summands = graph.coeff[j]^2 * bp.Vvar_fac[j]
            Vrow, errorV = kahan(summands, Vrow, errorV)

            summands = graph.coeff[j] * bp.Mvar_fac[j]
            Mrow, errorM = kahan(summands, Mrow, errorM)
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = bp.alpha1[j] * (graph.Mind[j] - (Mrow - graph.coeff[j] * bp.Mvar_fac[j]) - errorM) * graph.coeffInv[j] + bp.alpha2[j] * bp.Mfac_var[bp.to_var[j]]
            bp.Wfac_var[bp.to_var[j]] = (graph.coeff[j]^2) / (graph.Vind[j] + (Vrow - graph.coeff[j]^2 * bp.Vvar_fac[j]) + errorV)
        end
    end
end

### Efficient GBP with Kahan-Babuska algorithm: Factor to variable messages means only
@inline function kahan_factor_to_variable_mean(graph, bp)
    @inbounds for i = 1:graph.Nind
        Mrow = 0.0; errorM = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            summands = graph.coeff[j] * bp.Mvar_fac[j]
            Mrow, errorM = kahan(summands, Mrow, errorM)
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = (graph.Mind[j] - (Mrow - graph.coeff[j] * bp.Mvar_fac[j]) - errorM) * graph.coeffInv[j] 
        end
    end
end

### Efficient GBP with Kahan-Babuska algorithm: Factor to variable messages means only with damping
@inline function kahan_factor_to_variable_mean_damp(graph, bp)
    @inbounds for i = 1:graph.Nind
        Mrow = 0.0; errorM = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            summands = graph.coeff[j] * bp.Mvar_fac[j]
            Mrow, errorM = kahan(summands, Mrow, errorM)
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = bp.alpha1[j] * (graph.Mind[j] - (Mrow - graph.coeff[j] * bp.Mvar_fac[j]) - errorM) * graph.coeffInv[j] + bp.alpha2[j] * bp.Mfac_var[bp.to_var[j]]
        end
    end
end

###  Efficient GBP with Kahan-Babuska algorithm: Variable to factor messages
@inline function kahan_variable_to_factor(graph, bp)
    @inbounds for i = 1:graph.Nvar
        Mcol = 0.0; Wcol = 0.0; errorV = 0.0; errorM = 0.0

        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            Wcol, errorV = kahan(bp.Wfac_var[j], Wcol, errorV)

            summands = bp.Mfac_var[j] * bp.Wfac_var[j]
            Mcol, errorM = kahan(summands, Mcol, errorM)
        end
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            bp.Vvar_fac[bp.to_fac[j]] = 1 / ((Wcol - bp.Wfac_var[j]) + errorV + graph.Wdir[i])
            bp.Mvar_fac[bp.to_fac[j]] = ((Mcol - bp.Mfac_var[j] * bp.Wfac_var[j]) + errorM + graph.Mdir[i]) * bp.Vvar_fac[bp.to_fac[j]]
        end
    end
end

###  Efficient GBP with Kahan-Babuska algorithm: Variable to factor messages means only
@inline function kahan_variable_to_factor_mean(graph, bp)
    @inbounds for i = 1:graph.Nvar
        Mcol = 0.0; errorM = 0.0

        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            summands = bp.Mfac_var[j] * bp.Wfac_var[j]
            Mcol, errorM = kahan(summands, Mcol, errorM)
        end
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            bp.Mvar_fac[bp.to_fac[j]] = ((Mcol - bp.Mfac_var[j] * bp.Wfac_var[j]) + errorM + graph.Mdir[i]) * bp.Vvar_fac[bp.to_fac[j]]
        end
    end
end

###  Kahan-Babuska algorithm
@inline function kahan(summands, total, epsilon)
    t = total + summands
    if abs(total) >= abs(summands)
        epsilon += (total - t) + summands
    else
        epsilon += (summands - t) + total
    end
    total = t
    
    return total, epsilon
end
