### Efficient GBP: Factor to variable messages
@inline function efficient_factor_to_variable(graph, bp)
    @inbounds for i = 1:graph.Nind
        Mrow = 0.0; Vrow = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            Mrow += (graph.coeff[j] * bp.Mvar_fac[j])
            Vrow += (graph.coeff[j]^2 * bp.Vvar_fac[j])
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = (graph.Mind[j] - Mrow) * graph.coeffInv[j] + bp.Mvar_fac[j]
            bp.Wfac_var[bp.to_var[j]] = (graph.coeff[j]^2) / (graph.Vind[j] + Vrow - graph.coeff[j]^2 * bp.Vvar_fac[j])
        end
    end
end

### Efficient GBP: Factor to variable messages with damping
@inline function efficient_factor_to_variable_damp(graph, bp)
    @inbounds for i = 1:graph.Nind
        Mrow = 0.0; Vrow = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            Mrow += (graph.coeff[j] * bp.Mvar_fac[j])
            Vrow += (graph.coeff[j]^2 * bp.Vvar_fac[j])
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = bp.alpha1[j] * ((graph.Mind[j] - Mrow) * graph.coeffInv[j] + bp.Mvar_fac[j]) + bp.alpha2[j] * bp.Mfac_var[bp.to_var[j]]
            bp.Wfac_var[bp.to_var[j]] = (graph.coeff[j]^2) / (graph.Vind[j] + Vrow - graph.coeff[j]^2 * bp.Vvar_fac[j])
        end
    end
end

### Efficient GBP: Factor to variable means only
@inline function efficient_factor_to_variable_mean(graph, bp)
    @inbounds for i = 1:graph.Nind
        Mrow = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            Mrow += (graph.coeff[j] * bp.Mvar_fac[j])
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = (graph.Mind[j] - Mrow) * graph.coeffInv[j] + bp.Mvar_fac[j]
        end
    end
end

### Efficient GBP: Factor to variable means only with damping
@inline function efficient_factor_to_variable_mean_damp(graph, bp)
    @inbounds for i = 1:graph.Nind
        Mrow = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            Mrow += (graph.coeff[j] * bp.Mvar_fac[j])
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = bp.alpha1[j] * ((graph.Mind[j] - Mrow) * graph.coeffInv[j] + bp.Mvar_fac[j]) + bp.alpha2[j] * bp.Mfac_var[bp.to_var[j]]
        end
    end
end

### Efficient GBP: Variable to factor messages
@inline function efficient_variable_to_factor(graph, bp)
    @inbounds for i = 1:graph.Nvar
        Mcol = 0.0; Wcol = 0.0

        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            Mcol += bp.Mfac_var[j] * bp.Wfac_var[j]
            Wcol += bp.Wfac_var[j]
        end
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            bp.Vvar_fac[bp.to_fac[j]] = 1 / (Wcol + graph.Wdir[i] - bp.Wfac_var[j])
            bp.Mvar_fac[bp.to_fac[j]] = (Mcol - bp.Mfac_var[j] * bp.Wfac_var[j] + graph.Mdir[i]) * bp.Vvar_fac[bp.to_fac[j]]
        end
    end
end

### Efficient GBP: Variable to factor means only
@inline function efficient_variable_to_factor_mean(graph, bp)
    @inbounds for i = 1:graph.Nvar
        Mcol = 0.0

        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            Mcol += bp.Mfac_var[j] * bp.Wfac_var[j]
        end
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            bp.Mvar_fac[bp.to_fac[j]] = (Mcol - bp.Mfac_var[j] * bp.Wfac_var[j] + graph.Mdir[i]) * bp.Vvar_fac[bp.to_fac[j]]
        end
    end
end