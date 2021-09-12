########## Efficient GBP messages: Factor to variable ##########
function messageFactorVariableEfficient(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        Mrow = 0.0; Vrow = 0.0
        for j in gbp.graph.rowptr[i]
            Mrow += (gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j])
            Vrow += (gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j])
        end
        for j in gbp.graph.rowptr[i]
            q = gbp.graph.toVariable.nzval[j]
            gbp.inference.meanFactorVariable[q] = (gbp.graph.meanIndirect[j] - Mrow) / gbp.graph.coefficient[j] + gbp.inference.meanVariableFactor[j]
            gbp.inference.varianceFactorVariable[q] = (gbp.graph.varianceIndirect[j] + Vrow - gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]) / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Efficient GBP means: Factor to variable ##########
function meanFactorVariableEfficient(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        Mrow = 0.0
        for j in gbp.graph.rowptr[i]
            Mrow += (gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j])
        end
        for j in gbp.graph.rowptr[i]
            gbp.inference.meanFactorVariable[gbp.graph.toVariable.nzval[j]] = (gbp.graph.meanIndirect[j] - Mrow) / gbp.graph.coefficient[j] + gbp.inference.meanVariableFactor[j]
        end
    end
end

########## Efficient GBP variances: Factor to variable ##########
function varianceFactorVariableEfficient(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        Vrow = 0.0
        for j in gbp.graph.rowptr[i]
            Vrow += (gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j])
        end
        for j in gbp.graph.rowptr[i]
            gbp.inference.varianceFactorVariable[gbp.graph.toVariable.nzval[j]] = (gbp.graph.varianceIndirect[j] + Vrow - gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]) / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Efficient GBP damp messages: Factor to variable ##########
function messageDampFactorVariableEfficient(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        Mrow = 0.0; Vrow = 0.0
        for j in gbp.graph.rowptr[i]
            Mrow += (gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j])
            Vrow += (gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j])
        end
        for j in gbp.graph.rowptr[i]
            q = gbp.graph.toVariable.nzval[j]
            gbp.inference.meanFactorVariable[q] = gbp.graph.alphaNew[j] * ((gbp.graph.meanIndirect[j] - Mrow) / gbp.graph.coefficient[j] + gbp.inference.meanVariableFactor[j]) + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[q]
            gbp.inference.varianceFactorVariable[q] = (gbp.graph.varianceIndirect[j] + Vrow - gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]) / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Efficient GBP damp means: Factor to variable ##########
function meanDampFactorVariableEfficient(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        Mrow = 0.0
        for j in gbp.graph.rowptr[i]
            Mrow += (gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j])
        end
        for j in gbp.graph.rowptr[i]
            q = gbp.graph.toVariable.nzval[j]
            gbp.inference.meanFactorVariable[q] = gbp.graph.alphaNew[j] * ((gbp.graph.meanIndirect[j] - Mrow) / gbp.graph.coefficient[j] + gbp.inference.meanVariableFactor[j]) + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[q]
        end
    end
end

########## Efficient GBP messages: Variable to factor ##########
function messageVariableFactorEfficient(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateVariable
        Mcol = gbp.graph.meanDirect[i]; Wcol = gbp.graph.weightDirect[i]
        for j in gbp.graph.colptr[i]
            Mcol += gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
            Wcol += 1 / gbp.inference.varianceFactorVariable[j]
        end
        for j in gbp.graph.colptr[i]
            q = gbp.graph.toFactor.nzval[j]
            gbp.inference.varianceVariableFactor[q] = 1 / (Wcol - 1 / gbp.inference.varianceFactorVariable[j])
            gbp.inference.meanVariableFactor[q] = (Mcol - gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]) * gbp.inference.varianceVariableFactor[q]
        end
    end
end

########## Efficient GBP means: Variable to factor ##########
function meanVariableFactorEfficient(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateVariable
        Mcol = gbp.graph.meanDirect[i]
        for j in gbp.graph.colptr[i]
            Mcol += gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
        end
        for j in gbp.graph.colptr[i]
            q = gbp.graph.toFactor.nzval[j]
            gbp.inference.meanVariableFactor[q] = (Mcol - gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]) * gbp.inference.varianceVariableFactor[q]
        end
    end
end

########## Efficient GBP variances: Variable to factor ##########
function varianceVariableFactorEfficient(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateVariable
        Wcol = gbp.graph.weightDirect[i]
        for j in gbp.graph.colptr[i]
            Wcol += 1 / gbp.inference.varianceFactorVariable[j]
        end
        for j in gbp.graph.colptr[i]
            gbp.inference.varianceVariableFactor[gbp.graph.toFactor.nzval[j]] = 1 / (Wcol - 1 / gbp.inference.varianceFactorVariable[j])
        end
    end
end