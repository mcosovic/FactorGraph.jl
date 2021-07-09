### Vanilla GBP: Factor to variable messages
@inline function vanilla_factor_to_variable(graph, bp)
    @inbounds for i = 1:graph.Nind
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            Mrow = graph.Mind[j]; Vrow = graph.Vind[j]
            for k in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
                if j != k
                    Mrow -= (graph.coeff[k] * bp.Mvar_fac[k])
                    Vrow += (graph.coeff[k]^2 * bp.Vvar_fac[k])
                end
            end
            bp.Mfac_var[bp.to_var[j]] = Mrow * graph.coeffInv[j] 
            bp.Wfac_var[bp.to_var[j]] = graph.coeff[j]^2 / Vrow
        end
    end
end

### Vanilla GBP: Factor to variable messages with damping
@inline function vanilla_factor_to_variable_damp(graph, bp)
    @inbounds for i = 1:graph.Nind
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            Mrow = graph.Mind[j]; Vrow = graph.Vind[j]
            for k in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
                if j != k
                    Mrow -= (graph.coeff[k] * bp.Mvar_fac[k])
                    Vrow += (graph.coeff[k]^2 * bp.Vvar_fac[k])
                end
            end
            bp.Mfac_var[bp.to_var[j]] = bp.alpha1[j] * Mrow * graph.coeffInv[j] + bp.alpha2[j] * bp.Mfac_var[bp.to_var[j]]  
            bp.Wfac_var[bp.to_var[j]] = graph.coeff[j]^2 / Vrow
        end
    end
end

### Vanilla GBP: Factor to variable means only
@inline function vanilla_factor_to_variable_mean(graph, bp)
    @inbounds for i = 1:graph.Nind
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            Mrow = graph.Mind[j]
            for k in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
                if j != k
                    Mrow -= (graph.coeff[k] * bp.Mvar_fac[k])
                end
            end
            bp.Mfac_var[bp.to_var[j]] = Mrow * graph.coeffInv[j] 
        end
    end
end

### Vanilla GBP: Factor to variable means only with damping
@inline function vanilla_factor_to_variable_mean_damp(graph, bp)
    @inbounds for i = 1:graph.Nind
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            Mrow = graph.Mind[j]
            for k in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
                if j != k
                    Mrow -= (graph.coeff[k] * bp.Mvar_fac[k])
                end
            end
            bp.Mfac_var[bp.to_var[j]] = bp.alpha1[j] * Mrow * graph.coeffInv[j] + bp.alpha2[j] * bp.Mfac_var[bp.to_var[j]]  
        end
    end
end

### Vanilla GBP: Variable to factor messages
@inline function vanilla_variable_to_factor(graph, bp)
    @inbounds for i = 1:graph.Nvar
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            Mcol = graph.Mdir[i]; Wcol = graph.Wdir[i]
            for k in graph.colptr[i]:(graph.colptr[i + 1] - 1)
                if j != k
                    Mcol += bp.Mfac_var[k] * bp.Wfac_var[k]
                    Wcol += bp.Wfac_var[k]
                end
            end
            bp.Vvar_fac[bp.to_fac[j]] = 1 / Wcol
            bp.Mvar_fac[bp.to_fac[j]] = Mcol * bp.Vvar_fac[bp.to_fac[j]]
        end
    end
end

### Vanilla GBP: Variable to factor means only
@inline function vanilla_variable_to_factor_mean(graph, bp)
    @inbounds for i = 1:graph.Nvar
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            Mcol = graph.Mdir[i]
            for k in graph.colptr[i]:(graph.colptr[i + 1] - 1)
                if j != k
                    Mcol += bp.Mfac_var[k] * bp.Wfac_var[k]
                end
            end
            bp.Mvar_fac[bp.to_fac[j]] = Mcol * bp.Vvar_fac[bp.to_fac[j]]
        end
    end
end

### Compute marginals
function marginal(graph, bp, results)
    bp.iterCount[1] += 1
    @inbounds for i = 1:graph.Nvar
        Mcol = 0.0; Wcol = 0.0

        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            Mcol += bp.Mfac_var[j] * bp.Wfac_var[j]
            Wcol += bp.Wfac_var[j]
        end
        results.variances[i, bp.iterCount[1]] = 1 / (Wcol + graph.Wdir[i])
        results.means[i, bp.iterCount[1]] = (Mcol + graph.Mdir[i]) * results.variances[i, bp.iterCount[1]]
    end
end