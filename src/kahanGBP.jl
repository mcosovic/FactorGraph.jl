########## Efficient Kahan-Babuska GBP messages: Factor to variable ##########
function messageFactorVariableKahan(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        Mrow = 0.0; Vrow = 0.0; errorV = 0.0; errorM = 0.0

        for j in gbp.graph.rowptr[i]
            summands = gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]
            Vrow, errorV = kahanbabuska(summands, Vrow, errorV)

            summands = gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j]
            Mrow, errorM = kahanbabuska(summands, Mrow, errorM)
        end
        for j in gbp.graph.rowptr[i]
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.meanIndirect[j] - (Mrow - gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j]) - errorM) / gbp.graph.coefficient[j]
            gbp.inference.varianceFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.varianceIndirect[j] + (Vrow - gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]) + errorV) / (gbp.graph.coefficient[j]^2)
        end
    end
end

# ########## Efficient Kahan-Babuska GBP means: Factor to variable ##########
function meanFactorVariableKahan(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        Mrow = 0.0; errorM = 0.0

        for j in gbp.graph.rowptr[i]
            summands = gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j]
            Mrow, errorM = kahanbabuska(summands, Mrow, errorM)
        end
        for j in gbp.graph.rowptr[i]
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.meanIndirect[j] - (Mrow - gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j]) - errorM) / gbp.graph.coefficient[j]
        end
    end
end

########## Efficient Kahan-Babuska GBP variances: Factor to variable ##########
function varianceFactorVariableKahan(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        Vrow = 0.0; errorV = 0.0

        for j in gbp.graph.rowptr[i]
            summands = gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]
            Vrow, errorV = kahanbabuska(summands, Vrow, errorV)
        end
        for j in gbp.graph.rowptr[i]
            gbp.inference.varianceFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.varianceIndirect[j] + (Vrow - gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]) + errorV) / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Efficient Kahan-Babuska GBP damp messages: Factor to variable ##########
function messageDampFactorVariableKahan(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        Mrow = 0.0; Vrow = 0.0; errorV = 0.0; errorM = 0.0

        for j in gbp.graph.rowptr[i]
            summands = gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]
            Vrow, errorV = kahanbabuska(summands, Vrow, errorV)

            summands = gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j]
            Mrow, errorM = kahanbabuska(summands, Mrow, errorM)
        end
        for j in gbp.graph.rowptr[i]
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.alphaNew[j] * (gbp.graph.meanIndirect[j] - (Mrow - gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j]) - errorM) / gbp.graph.coefficient[j]) + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]]
            gbp.inference.varianceFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.varianceIndirect[j] + (Vrow - gbp.graph.coefficient[j]^2 * gbp.inference.varianceVariableFactor[j]) + errorV) / (gbp.graph.coefficient[j]^2)
        end
    end
end

########## Efficient Kahan-Babuska GBP damp means: Factor to variable ##########
function meanDampFactorVariableKahan(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateFactor
        Mrow = 0.0; errorM = 0.0

        for j in gbp.graph.rowptr[i]
            summands = gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j]
            Mrow, errorM = kahanbabuska(summands, Mrow, errorM)
        end
        for j in gbp.graph.rowptr[i]
            gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]] = (gbp.graph.alphaNew[j] * (gbp.graph.meanIndirect[j] - (Mrow - gbp.graph.coefficient[j] * gbp.inference.meanVariableFactor[j]) - errorM) / gbp.graph.coefficient[j]) + gbp.graph.alphaOld[j] * gbp.inference.meanFactorVariable[gbp.graph.sendToVariable[j]]
        end
    end
end

########## Efficient Kahan-Babuska GBP messages: Variable to factor ##########
function messageVariableFactorKahan(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateVariable
        Mcol = 0.0; Wcol = 0.0; errorV = 0.0; errorM = 0.0

        for j in gbp.graph.colptr[i]
            Wcol, errorV = kahanbabuska(1 / gbp.inference.varianceFactorVariable[j], Wcol, errorV)

            summands = gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
            Mcol, errorM = kahanbabuska(summands, Mcol, errorM)
        end
        for j in gbp.graph.colptr[i]
            gbp.inference.varianceVariableFactor[gbp.graph.sendToFactor[j]] = 1 / ((Wcol - 1 / gbp.inference.varianceFactorVariable[j]) + errorV + gbp.graph.weightDirect[i])
            gbp.inference.meanVariableFactor[gbp.graph.sendToFactor[j]] = ((Mcol - gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]) + errorM + gbp.graph.meanDirect[i]) * gbp.inference.varianceVariableFactor[gbp.graph.sendToFactor[j]]
        end
    end
end

########## Efficient Kahan-Babuska GBP means: Variable to factor ##########
function meanVariableFactorKahan(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateVariable
        Mcol = 0.0; errorM = 0.0

        for j in gbp.graph.colptr[i]
            summands = gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]
            Mcol, errorM = kahanbabuska(summands, Mcol, errorM)
        end
        for j in gbp.graph.colptr[i]
            gbp.inference.meanVariableFactor[gbp.graph.sendToFactor[j]] = ((Mcol - gbp.inference.meanFactorVariable[j] / gbp.inference.varianceFactorVariable[j]) + errorM + gbp.graph.meanDirect[i]) * gbp.inference.varianceVariableFactor[gbp.graph.sendToFactor[j]]
        end
    end
end

########## Efficient Kahan-Babuska GBP variances: Variable to factor ##########
function varianceVariableFactorKahan(gbp::GraphicalModel)
    @inbounds Threads.@threads for i in gbp.graph.iterateVariable
        Wcol = 0.0; errorV = 0.0

        for j in gbp.graph.colptr[i]
            Wcol, errorV = kahanbabuska(1 / gbp.inference.varianceFactorVariable[j], Wcol, errorV)
        end
        for j in gbp.graph.colptr[i]
            gbp.inference.varianceVariableFactor[gbp.graph.sendToFactor[j]] = 1 / ((Wcol - 1 / gbp.inference.varianceFactorVariable[j]) + errorV + gbp.graph.weightDirect[i])
        end
    end
end

########## Kahan-Babuska algorithm ##########
@inline function kahanbabuska(summands, total, epsilon)
    t = total + summands
    if abs(total) >= abs(summands)
        epsilon += (total - t) + summands
    else
        epsilon += (summands - t) + total
    end
    total = t

    return total, epsilon
end