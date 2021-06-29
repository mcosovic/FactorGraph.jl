### Factor to variable messages
@inline function factor_to_variable(graph, bp)
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

### Factor to variable messages with damping
@inline function factor_to_variable_damp(graph, bp)
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

### Factor to variable means only
@inline function factor_to_variable_mean(graph, bp)
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

### Factor to variable means only with damping
@inline function factor_to_variable_mean_damp(graph, bp)
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


### Variable to factor messages
@inline function variable_to_factor(graph, bp)
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

### Variable to factor means only
@inline function variable_to_factor_mean(graph, bp)
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

### Factor to variable messages with Kahan-Babuska summation algorithm
@inline function factor_to_variable_kahan(graph, bp)
    @inbounds for i = 1:graph.Nind
        Mrow = 0.0; Vrow = 0.0; error = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            Mrow += (graph.coeff[j] * bp.Mvar_fac[j])

            x = graph.coeff[j]^2 * bp.Vvar_fac[j]
            t = Vrow + x
            if abs(Vrow) >= abs(x)
                error += (Vrow - t) + x
            else
                error += (x - t) + Vrow
            end
            Vrow = t
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = (graph.Mind[j] - Mrow) * graph.coeffInv[j] + bp.Mvar_fac[j]
            bp.Wfac_var[bp.to_var[j]] = (graph.coeff[j]^2) / (graph.Vind[j] + (Vrow - graph.coeff[j]^2 * bp.Vvar_fac[j]) + error)
        end
    end
end

### Factor to variable messages with damping and Kahan-Babuska summation algorithm
@inline function factor_to_variable_kahan_damp(graph, bp)
    @inbounds for i = 1:graph.Nind
        Mrow = 0.0; Vrow = 0.0; error = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            Mrow += (graph.coeff[j] * bp.Mvar_fac[j])

            x = graph.coeff[j]^2 * bp.Vvar_fac[j]
            t = Vrow + x
            if abs(Vrow) >= abs(x)
                error += (Vrow - t) + x
            else
                error += (x - t) + Vrow
            end
            Vrow = t
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = bp.alpha1[j] * ((graph.Mind[j] - Mrow) * graph.coeffInv[j] + bp.Mvar_fac[j]) + bp.alpha2[j] * bp.Mfac_var[bp.to_var[j]]
            bp.Wfac_var[bp.to_var[j]] = (graph.coeff[j]^2) / (graph.Vind[j] + (Vrow - graph.coeff[j]^2 * bp.Vvar_fac[j]) + error)
        end
    end
end

### Variable to factor messages using Kahan-Babuska summation algorithm
@inline function variable_to_factor_kahan(graph, bp)
    @inbounds for i = 1:graph.Nvar
        Mcol = 0.0; Wcol = 0.0; error = 0.0

        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            Mcol += bp.Mfac_var[j] * bp.Wfac_var[j]

            t = Wcol + bp.Wfac_var[j]
            if abs(Wcol) >= abs(bp.Wfac_var[j])
                error += (Wcol - t) + bp.Wfac_var[j]
            else
                error += (bp.Wfac_var[j] - t) + Wcol
            end
            Wcol = t
        end
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            bp.Vvar_fac[bp.to_fac[j]] = 1 / ((Wcol - bp.Wfac_var[j]) + error + graph.Wdir[i])
            bp.Mvar_fac[bp.to_fac[j]] = (Mcol -bp.Mfac_var[j] * bp.Wfac_var[j] + graph.Mdir[i]) * bp.Vvar_fac[bp.to_fac[j]]
        end
    end
end

### Factor to variable messages in recursion mode
function factor_recursion(graph, bp, rec)
    @inbounds for i = 1:graph.Nvar
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            rec.Mcol[i] += bp.Mfac_var[j] * bp.Wfac_var[j]
            rec.Wcol[i] += bp.Wfac_var[j]
        end
    end

    @inbounds for i = 1:graph.Nind
        Mrow = 0.0; Vrow = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Vvar_fac[j] = 1 / (rec.Wcol[graph.col[j]] + graph.Wdir[graph.col[j]] - bp.Wfac_var[bp.to_var[j]])
            bp.Mvar_fac[j] = (rec.Mcol[graph.col[j]] - bp.Mfac_var[bp.to_var[j]] * bp.Wfac_var[bp.to_var[j]] + graph.Mdir[graph.col[j]]) * bp.Vvar_fac[j]

            Mrow += (graph.coeff[j] * bp.Mvar_fac[j])
            Vrow += (graph.coeff[j]^2 * bp.Vvar_fac[j])
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = (graph.Mind[j] - Mrow) * graph.coeffInv[j] + bp.Mvar_fac[j]
            bp.Wfac_var[bp.to_var[j]] = (graph.coeff[j]^2) / (graph.Vind[j] + Vrow - graph.coeff[j]^2 * bp.Vvar_fac[j])
        end
    end

    @inbounds for i = 1:graph.Nvar
        rec.Mcol[i] = 0.0
        rec.Wcol[i] = 0.0
    end
end


### Factor to variable messages in recursion mode
function factor_recursion_damp(graph, bp, rec)
    @inbounds for i = 1:graph.Nvar
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            rec.Mcol[i] += bp.Mfac_var[j] * bp.Wfac_var[j]
            rec.Wcol[i] += bp.Wfac_var[j]
        end
    end

    @inbounds for i = 1:graph.Nind
        Mrow = 0.0; Vrow = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Vvar_fac[j] = 1 / (rec.Wcol[graph.col[j]] + graph.Wdir[graph.col[j]] - bp.Wfac_var[bp.to_var[j]])
            bp.Mvar_fac[j] = (rec.Mcol[graph.col[j]] - bp.Mfac_var[bp.to_var[j]] * bp.Wfac_var[bp.to_var[j]] + graph.Mdir[graph.col[j]]) * bp.Vvar_fac[j]

            Mrow += (graph.coeff[j] * bp.Mvar_fac[j])
            Vrow += (graph.coeff[j]^2 * bp.Vvar_fac[j])
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = bp.alpha1[j] * ((graph.Mind[j] - Mrow) * graph.coeffInv[j] + bp.Mvar_fac[j]) + bp.alpha2[j] * bp.Mfac_var[bp.to_var[j]]
            bp.Wfac_var[bp.to_var[j]] = (graph.coeff[j]^2) / (graph.Vind[j] + Vrow - graph.coeff[j]^2 * bp.Vvar_fac[j])
        end
    end

    @inbounds for i = 1:graph.Nvar
        rec.Mcol[i] = 0.0
        rec.Wcol[i] = 0.0
    end
end

### Factor to variable means only in recursion mode
function factor_recursion_mean(graph, bp, rec)
    @inbounds for i = 1:graph.Nvar
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            rec.Mcol[i] += bp.Mfac_var[j] * bp.Wfac_var[j]
        end
    end

    @inbounds for i = 1:graph.Nind
        Mrow = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mvar_fac[j] = (rec.Mcol[graph.col[j]] - bp.Mfac_var[bp.to_var[j]] * bp.Wfac_var[bp.to_var[j]] + graph.Mdir[graph.col[j]]) * bp.Vvar_fac[j]

            Mrow += (graph.coeff[j] * bp.Mvar_fac[j])
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = (graph.Mind[j] - Mrow) * graph.coeffInv[j] + bp.Mvar_fac[j]
        end
    end

    @inbounds for i = 1:graph.Nvar
        rec.Mcol[i] = 0.0
    end
end

### Factor to variable means only in recursion mode with damping
function factor_recursion_mean_damp(graph, bp, rec)
    @inbounds for i = 1:graph.Nvar
        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            rec.Mcol[i] += bp.Mfac_var[j] * bp.Wfac_var[j]
        end
    end

    @inbounds for i = 1:graph.Nind
        Mrow = 0.0

        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mvar_fac[j] = (rec.Mcol[graph.col[j]] - bp.Mfac_var[bp.to_var[j]] * bp.Wfac_var[bp.to_var[j]] + graph.Mdir[graph.col[j]]) * bp.Vvar_fac[j]

            Mrow += (graph.coeff[j] * bp.Mvar_fac[j])
        end
        for j in graph.rowptr[i]:(graph.rowptr[i + 1] - 1)
            bp.Mfac_var[bp.to_var[j]] = bp.alpha1[j] * ((graph.Mind[j] - Mrow) * graph.coeffInv[j] + bp.Mvar_fac[j]) + bp.alpha2[j] * bp.Mfac_var[bp.to_var[j]]
        end
    end

    @inbounds for i = 1:graph.Nvar
        rec.Mcol[i] = 0.0
    end
end



### Compute marginals
function marginal(graph, bp)
    Xbp = fill(0.0, graph.Nvar)
    @inbounds for i = 1:graph.Nvar
        Mcol = 0.0; Wcol = 0.0

        for j in graph.colptr[i]:(graph.colptr[i + 1] - 1)
            Mcol += bp.Mfac_var[j] * bp.Wfac_var[j]
            Wcol += bp.Wfac_var[j]
        end
        Xbp[i] = (Mcol + graph.Mdir[i]) / (Wcol + graph.Wdir[i])
    end

    return Xbp
end
