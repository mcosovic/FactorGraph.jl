########## Vanilla GBP messages: Factor to variable ##########
function messageFactorVariable(gbp::ContinuousModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]
            Mrow = gbp.graph.meanIndirect[j]; Vrow = gbp.graph.varianceIndirect[j]; q = gbp.graph.toVariable.nzval[j]
            for k in gbp.graph.rowptr[i]
                if j != k
                    Mrow -= gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k]
                    Vrow += gbp.graph.coefficient[k]^2 * gbp.inference.varianceVariableFactor[k]
                end
            end
            gbp.inference.meanFactorVariable[q] = Mrow / gbp.graph.coefficient[j]
            gbp.inference.varianceFactorVariable[q] = Vrow / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Vanilla GBP means: Factor to variable ##########
function meanFactorVariable(gbp::ContinuousModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]
            Mrow = gbp.graph.meanIndirect[j]
            for k in gbp.graph.rowptr[i]
                if j != k
                    Mrow -= gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k]
                end
            end
            gbp.inference.meanFactorVariable[gbp.graph.toVariable.nzval[j]] = Mrow / gbp.graph.coefficient[j]
        end
    end
end

########## Vanilla GBP variances: Factor to variable ##########
function varianceFactorVariable(gbp::ContinuousModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]
            Vrow = gbp.graph.varianceIndirect[j]
            for k in gbp.graph.rowptr[i]
                if j != k
                    Vrow += gbp.graph.coefficient[k]^2 * gbp.inference.varianceVariableFactor[k]
                end
            end
            gbp.inference.varianceFactorVariable[gbp.graph.toVariable.nzval[j]] = Vrow / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Vanilla GBP damp messages: Factor to variable ##########
function messageDampFactorVariable(gbp::ContinuousModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]
            Mrow = gbp.graph.meanIndirect[j]; Vrow = gbp.graph.varianceIndirect[j]; q = gbp.graph.toVariable.nzval[j]
            for k in gbp.graph.rowptr[i]
                if j != k
                    Mrow -= gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k]
                    Vrow += gbp.graph.coefficient[k]^2 * gbp.inference.varianceVariableFactor[k]
                end
            end
            gbp.inference.meanFactorVariable[q] = gbp.graph.alphaNew[j] * Mrow / gbp.graph.coefficient[j] + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[q]
            gbp.inference.varianceFactorVariable[q] = Vrow / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Vanilla GBP damp means: Factor to variable ##########
function meanDampFactorVariable(gbp::ContinuousModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        for j in gbp.graph.rowptr[i]
            Mrow = gbp.graph.meanIndirect[j]; q = gbp.graph.toVariable.nzval[j]
            for k in gbp.graph.rowptr[i]
                if j != k
                    Mrow -= (gbp.graph.coefficient[k] * gbp.inference.meanVariableFactor[k])
                end
            end
            gbp.inference.meanFactorVariable[q] = gbp.graph.alphaNew[j] * Mrow / gbp.graph.coefficient[j] + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[q]
        end
    end
end

############# Vanilla GBP messages: Variable to factor ##########
function messageVariableFactor(gbp::ContinuousModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateVariable
        for j in gbp.graph.colptr[i]
            Mcol = gbp.graph.meanDirect[i]; Wcol = gbp.graph.weightDirect[i]; q = gbp.graph.toFactor.nzval[j]
            for k in gbp.graph.colptr[i]
                if j != k
                    Mcol += gbp.inference.meanFactorVariable[k] / gbp.inference.varianceFactorVariable[k]
                    Wcol += 1 / gbp.inference.varianceFactorVariable[k]
                end
            end
            gbp.inference.varianceVariableFactor[q] = 1 / Wcol
            gbp.inference.meanVariableFactor[q] = Mcol / Wcol
        end
    end
end

############# Vanilla GBP means: Variable to factor ##########
function meanVariableFactor(gbp::ContinuousModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateVariable
        for j in gbp.graph.colptr[i]
            Mcol = gbp.graph.meanDirect[i]; q = gbp.graph.toFactor.nzval[j]
            for k in gbp.graph.colptr[i]
                if j != k
                    Mcol += gbp.inference.meanFactorVariable[k] / gbp.inference.varianceFactorVariable[k]
                end
            end
            gbp.inference.meanVariableFactor[q] = Mcol * gbp.inference.varianceVariableFactor[q]
        end
    end
end

############# Vanilla GBP variances: Variable to factor ##########
function varianceVariableFactor(gbp::ContinuousModel)
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