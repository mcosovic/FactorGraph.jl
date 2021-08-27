function messageFactorVariableVanilla(gbp::GraphicalModel)
    @inbounds for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            Mrow = gbp.graph.meanIndirect[j]; Vrow = gbp.graph.varianceIndirect[j]
            for k in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
                if j != k
                    Mrow -= (gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k])
                    Vrow += (gbp.graph.coefficient[k]^2 * gbp.inference.varianceVariableFactor[k])
                end
            end
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = Mrow / gbp.graph.coefficient[j]
            gbp.inference.varianceFactorVariable[gbp.graph.sendToVariable[j]] = Vrow / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Vanilla GBP means: Factor to variable ##########
function meanFactorVariableVanilla(gbp::GraphicalModel)
    @inbounds for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            Mrow = gbp.graph.meanIndirect[j]
            for k in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
                if j != k
                    Mrow -= (gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k])
                end
            end
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = Mrow / gbp.graph.coefficient[j]
        end
    end
end

########## Vanilla GBP variances: Factor to variable ##########
function varianceFactorVariableVanilla(gbp::GraphicalModel)
    @inbounds for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            Vrow = gbp.graph.varianceIndirect[j]
            for k in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
                if j != k
                    Vrow += (gbp.graph.coefficient[k]^2 * gbp.inference.varianceVariableFactor[k])
                end
            end
            gbp.inference.varianceFactorVariable[gbp.graph.sendToVariable[j]] = Vrow / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Vanilla GBP damp messages: Factor to variable ##########
function messageDampFactorVariableVanilla(gbp::GraphicalModel)
    @inbounds for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            Mrow = gbp.graph.meanIndirect[j]; Vrow = gbp.graph.varianceIndirect[j]
            for k in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
                if j != k
                    Mrow -= (gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k])
                    Vrow += (gbp.graph.coefficient[k]^2 * gbp.inference.varianceVariableFactor[k])
                end
            end
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = gbp.graph.alphaNew[j] * Mrow / gbp.graph.coefficient[j] + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]]
            gbp.inference.varianceFactorVariable[gbp.graph.sendToVariable[j]] = Vrow / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Vanilla GBP damp means: Factor to variable ##########
function meanDampFactorVariableVanilla(gbp::GraphicalModel)
    @inbounds for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
            Mrow = gbp.graph.meanIndirect[j]
            for k in gbp.graph.rowptr[i]:(gbp.graph.rowptr[i + 1] - 1)
                if j != k
                    Mrow -= (gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k])
                end
            end
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = gbp.graph.alphaNew[j] * Mrow / gbp.graph.coefficient[j] + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]]
        end
    end
end

############# Vanilla GBP messages: Variable to factor ##########
function messageVariableFactorVanilla(gbp::GraphicalModel)
    @inbounds for i in gbp.graph.iterateVariable
        for j in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
            Mcol = gbp.graph.meanDirect[i]; Wcol = gbp.graph.weightDirect[i]
            for k in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
                if j != k
                    Mcol += gbp.inference.meanFactorVariable[k] / gbp.inference.varianceFactorVariable[k]
                    Wcol += 1 / gbp.inference.varianceFactorVariable[k]
                end
            end
            gbp.inference.varianceVariableFactor[gbp.graph.sendToFactor[j]] = 1 / Wcol
            gbp.inference.meanVariableFactor[gbp.graph.sendToFactor[j]] = Mcol / Wcol
        end
    end
end

############# Vanilla GBP means: Variable to factor ##########
function meanVariableFactorVanilla(gbp::GraphicalModel)
    @inbounds for i = 1:gbp.graph.Nvariable
        for j in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
            Mcol = gbp.graph.meanDirect[i]
            for k in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
                if j != k
                    Mcol += gbp.inference.meanFactorVariable[k] / gbp.inference.varianceFactorVariable[k]
                end
            end
            gbp.inference.meanVariableFactor[gbp.graph.sendToFactor[j]] = Mcol * gbp.inference.varianceVariableFactor[gbp.graph.sendToFactor[j]]
        end
    end
end

############# Vanilla GBP variances: Variable to factor ##########
function varianceVariableFactorVanilla(gbp::GraphicalModel)
    @inbounds for i = 1:gbp.graph.Nvariable
        for j in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
            Wcol = gbp.graph.weightDirect[i]
            for k in gbp.graph.colptr[i]:(gbp.graph.colptr[i + 1] - 1)
                if j != k
                    Wcol += 1 / gbp.inference.varianceFactorVariable[k]
                end
            end
            gbp.inference.varianceVariableFactor[gbp.graph.sendToFactor[j]] = 1 / Wcol
        end
    end
end