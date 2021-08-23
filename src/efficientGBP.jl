########## Efficient GBP messages: Factor to variable ##########
function messageFactorVariableEfficient(gbp::GraphicalModel)
    @inbounds for i = 1:gbp.graph.Nindirect
        Mrow = 0.0; Vrow = 0.0
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            Mrow += (gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j])
            Vrow += (gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j])
        end
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.meanIndirect[j] - Mrow) / gbp.graph.coefficient[j] + gbp.inference.meanVariableFactor[j]
            gbp.inference.varianceFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.varianceIndirect[j] + Vrow - gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]) / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Efficient GBP means: Factor to variable ##########
function meanFactorVariableEfficient(gbp::GraphicalModel)
    @inbounds for i = 1:gbp.graph.Nindirect
        Mrow = 0.0
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            Mrow += (gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j])
        end
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.meanIndirect[j] - Mrow) / gbp.graph.coefficient[j] + gbp.inference.meanVariableFactor[j]
        end
    end
end

########## Efficient GBP variances: Factor to variable ##########
function varianceFactorVariableEfficient(gbp::GraphicalModel)
    @inbounds for i = 1:gbp.graph.Nindirect
        Vrow = 0.0
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            Vrow += (gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j])
        end
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            gbp.inference.varianceFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.varianceIndirect[j] + Vrow - gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]) / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Efficient GBP damp messages: Factor to variable ##########
function messageDampFactorVariableEfficient(gbp::GraphicalModel)
    @inbounds for i = 1:gbp.graph.Nindirect
        Mrow = 0.0; Vrow = 0.0
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            Mrow += (gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j])
            Vrow += (gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j])
        end
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = gbp.graph.alphaNew[j] * ((gbp.graph.meanIndirect[j] - Mrow) / gbp.graph.coefficient[j] + gbp.inference.meanVariableFactor[j]) + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]]
            gbp.inference.varianceFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.varianceIndirect[j] + Vrow - gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]) / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Efficient GBP damp means: Factor to variable ##########
function meanDampFactorVariableEfficient(gbp::GraphicalModel)
    @inbounds for i = 1:gbp.graph.Nindirect
        Mrow = 0.0
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            Mrow += (gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j])
        end
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = gbp.graph.alphaNew[j] * ((gbp.graph.meanIndirect[j] - Mrow) / gbp.graph.coefficient[j] + gbp.inference.meanVariableFactor[j]) + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]]
        end
    end
end

########## Efficient GBP messages: Variable to factor ##########
function messageVariableFactorEfficient(gbp::GraphicalModel)
    @inbounds for i = 1:gbp.graph.Nvariable
        Mcol = gbp.graph.meanDirect[i]; Wcol = gbp.graph.weightDirect[i]
        for j in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
            Mcol += gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
            Wcol += 1 / gbp.inference.varianceFactorVariable[j]
        end
        for j in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
            gbp.inference.varianceVariableFactor[gbp.graph.sendToFactor[j]] = 1 / (Wcol - 1 / gbp.inference.varianceFactorVariable[j])
            gbp.inference.meanVariableFactor[gbp.graph.sendToFactor[j]] = (Mcol - gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]) * gbp.inference.varianceVariableFactor[gbp.graph.sendToFactor[j]]
        end
    end
end

########## Efficient GBP means: Variable to factor ##########
function meanVariableFactorEfficient(gbp::GraphicalModel)
    @inbounds for i = 1:gbp.graph.Nvariable
        Mcol = gbp.graph.meanDirect[i]
        for j in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
            Mcol += gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
        end
        for j in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
            gbp.inference.meanVariableFactor[gbp.graph.sendToFactor[j]] = (Mcol - gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]) * gbp.inference.varianceVariableFactor[gbp.graph.sendToFactor[j]]
        end
    end
end

########## Efficient GBP variances: Variable to factor ##########
function varianceVariableFactorEfficient(gbp::GraphicalModel)
    @inbounds for i = 1:gbp.graph.Nvariable
        Wcol = gbp.graph.weightDirect[i]
        for j in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
            Wcol += 1 / gbp.inference.varianceFactorVariable[j]
        end
        for j in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
            gbp.inference.varianceVariableFactor[gbp.graph.sendToFactor[j]] = 1 / (Wcol - 1 / gbp.inference.varianceFactorVariable[j])
        end
    end
end