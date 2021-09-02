function messageFactorVariableVanilla(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]
            Mrow = gbp.graph.meanIndirect[j]; Vrow = gbp.graph.varianceIndirect[j]
            for k in gbp.graph.rowptr[i]
                if j != k
                    Mrow -= (gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k])
                    Vrow += (gbp.graph.coefficient[k]^2 * gbp.inference.varianceVariableFactor[k])
                end
            end
            gbp.inference.meanFactorVariable[gbp.graph.toVariable.nzval[j]] = Mrow / gbp.graph.coefficient[j]
            gbp.inference.varianceFactorVariable[gbp.graph.toVariable.nzval[j]] = Vrow / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Vanilla GBP means: Factor to variable ##########
function meanFactorVariableVanilla(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]
            Mrow = gbp.graph.meanIndirect[j]
            for k in gbp.graph.rowptr[i]
                if j != k
                    Mrow -= (gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k])
                end
            end
            gbp.inference.meanFactorVariable[gbp.graph.toVariable.nzval[j]] = Mrow / gbp.graph.coefficient[j]
        end
    end
end

########## Vanilla GBP variances: Factor to variable ##########
function varianceFactorVariableVanilla(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]
            Vrow = gbp.graph.varianceIndirect[j]
            for k in gbp.graph.rowptr[i]
                if j != k
                    Vrow += (gbp.graph.coefficient[k]^2 * gbp.inference.varianceVariableFactor[k])
                end
            end
            gbp.inference.varianceFactorVariable[gbp.graph.toVariable.nzval[j]] = Vrow / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Vanilla GBP damp messages: Factor to variable ##########
function messageDampFactorVariableVanilla(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]
            Mrow = gbp.graph.meanIndirect[j]; Vrow = gbp.graph.varianceIndirect[j]
            for k in gbp.graph.rowptr[i]
                if j != k
                    Mrow -= (gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k])
                    Vrow += (gbp.graph.coefficient[k]^2 * gbp.inference.varianceVariableFactor[k])
                end
            end
            gbp.inference.meanFactorVariable[gbp.graph.toVariable.nzval[j]] = gbp.graph.alphaNew[j] * Mrow / gbp.graph.coefficient[j] + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[gbp.graph.toVariable.nzval[j]]
            gbp.inference.varianceFactorVariable[gbp.graph.toVariable.nzval[j]] = Vrow / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Vanilla GBP damp means: Factor to variable ##########
function meanDampFactorVariableVanilla(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]
            Mrow = gbp.graph.meanIndirect[j]
            for k in gbp.graph.rowptr[i]
                if j != k
                    Mrow -= (gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k])
                end
            end
            gbp.inference.meanFactorVariable[gbp.graph.toVariable.nzval[j]] = gbp.graph.alphaNew[j] * Mrow / gbp.graph.coefficient[j] + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[gbp.graph.toVariable.nzval[j]]
        end
    end
end

############# Vanilla GBP messages: Variable to factor ##########
function messageVariableFactorVanilla(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateVariable
        for j in gbp.graph.colptr[i]
            Mcol = gbp.graph.meanDirect[i]; Wcol = gbp.graph.weightDirect[i]
            for k in gbp.graph.colptr[i]
                if j != k
                    Mcol += gbp.inference.meanFactorVariable[k] / gbp.inference.varianceFactorVariable[k]
                    Wcol += 1 / gbp.inference.varianceFactorVariable[k]
                end
            end
            gbp.inference.varianceVariableFactor[gbp.graph.toFactor.nzval[j]] = 1 / Wcol
            gbp.inference.meanVariableFactor[gbp.graph.toFactor.nzval[j]] = Mcol / Wcol
        end
    end
end

############# Vanilla GBP means: Variable to factor ##########
function meanVariableFactorVanilla(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateVariable
        for j in gbp.graph.colptr[i]
            Mcol = gbp.graph.meanDirect[i]
            for k in gbp.graph.colptr[i]
                if j != k
                    Mcol += gbp.inference.meanFactorVariable[k] / gbp.inference.varianceFactorVariable[k]
                end
            end
            gbp.inference.meanVariableFactor[gbp.graph.toFactor.nzval[j]] = Mcol * gbp.inference.varianceVariableFactor[gbp.graph.toFactor.nzval[j]]
        end
    end
end

############# Vanilla GBP variances: Variable to factor ##########
function varianceVariableFactorVanilla(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateVariable
        for j in gbp.graph.colptr[i]
            Wcol = gbp.graph.weightDirect[i]
            for k in gbp.graph.colptr[i]
                if j != k
                    Wcol += 1 / gbp.inference.varianceFactorVariable[k]
                end
            end
            gbp.inference.varianceVariableFactor[gbp.graph.toFactor.nzval[j]] = 1 / Wcol
        end
    end
end