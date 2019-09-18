###########################################################
#  Produce the factor graph and define singly-connected,  #
#  virtual and indirect factor nodes                      #
###########################################################


#-------------------------------------------------------------------------------
# Define the node numbers and the transpose system matrix
#-------------------------------------------------------------------------------
function graph(jacobian)
    Nfactor, Nvariable = size(jacobian)
    jacobianT = SparseMatrixCSC(jacobian')

    return Nfactor, Nvariable, jacobianT
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Find number of links
#-------------------------------------------------------------------------------
function links(Nfactor, jacobianT)
    Ndir = 0
    Nlink = 0
    variables_in_column = 0
    dir_position = fill(0, Nfactor)

    @inbounds for i = 1:Nfactor
        variables_in_column = jacobianT.colptr[i + 1] - jacobianT.colptr[i]

        if variables_in_column == 1
            Ndir += variables_in_column
            dir_position[i] = jacobianT.rowval[jacobianT.colptr[i]]
        else
            Nlink += variables_in_column
        end
    end

    return Ndir, Nlink, dir_position
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Define position of virtual factor nodes
#-------------------------------------------------------------------------------
function virtuals(Nvariable, dir_position)
    idx = findall(!iszero, dir_position)
    virtual = fill(1, Nvariable)

    @inbounds for i in idx
        virtual[dir_position[i]] = 0
    end

    return virtual
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Define singly-connected, virtual and indirect factor nodes arrays
#-------------------------------------------------------------------------------
function factors(
    Nfactor, Nvariable, Ndir, Nlink,
    jacobianT, observation, noise, virtual,
    MEAN, VARI)

    Nind = Nfactor - Ndir

    Mdir = fill(0.0, Nvariable)
    VdirInv = fill(0.0, Nvariable)

    row = fill(0, Nlink)
    col = similar(row)
    coeff = fill(0.0, Nlink)
    coeffInv = similar(coeff)
    Mind = fill(0.0, Nind)
    Vind = similar(Mind)

    factor_colptr = fill(0, Nind)
    variable_colptr = fill(0, Nvariable)

    idxi = 1
    idxr = 1

    idxT = findall(!iszero, jacobianT)
    @inbounds for i in idxT
        if (jacobianT.colptr[i[2] + 1] - jacobianT.colptr[i[2]]) == 1
            Mdir[i[1]] += observation[i[2]] * jacobianT[i] / noise[i[2]]
            VdirInv[i[1]] += (jacobianT[i]^2) / noise[i[2]]
        else
            row[idxi] = idxr
            col[idxi] = i[1]

            coeff[idxi] = jacobianT[i]
            coeffInv[idxi] = 1 / jacobianT[i]

            idxi += 1

            variable_colptr[i[1]] += 1
            factor_colptr[idxr] = idxi

            if idxT[jacobianT.colptr[i[2] + 1] - 1] == i
                Mind[idxr] = observation[i[2]]
                Vind[idxr] = noise[i[2]]
                idxr += 1
            end
        end
        if virtual[i[1]] !== 0
            Mdir[i[1]] = MEAN / VARI
            VdirInv[i[1]] = 1 / VARI
        end
    end

    pushfirst!(variable_colptr, 1)
    variable_colptr = cumsum(variable_colptr)
    pushfirst!(factor_colptr, 1)

    return row, col, Nind, Mind, Vind, coeff, coeffInv, Mdir, VdirInv,
           variable_colptr, factor_colptr
end
#-------------------------------------------------------------------------------
