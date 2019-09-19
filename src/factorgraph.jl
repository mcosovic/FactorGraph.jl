###########################################################
#  Produce the factor graph and define singly-connected,  #
#  virtual and indirect factor nodes                      #
###########################################################


#-------------------------------------------------------------------------------
# Define the node numbers and the transpose system matrix
#-------------------------------------------------------------------------------
function graph(jacobian)
    Nfac, Nvar = size(jacobian)
    jacobianT = SparseMatrixCSC(jacobian')

    return Nfac, Nvar, jacobianT
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Find number of links
#-------------------------------------------------------------------------------
function links(Nfac, jacobianT)
    Ndir = 0
    Nlink = 0
    var_in_column = 0
    dir_position = fill(0, Nfac)

    @inbounds for i = 1:Nfac
        var_in_column = jacobianT.colptr[i + 1] - jacobianT.colptr[i]

        if var_in_column == 1
            Ndir += var_in_column
            dir_position[i] = jacobianT.rowval[jacobianT.colptr[i]]
        else
            Nlink += var_in_column
        end
    end

    return Ndir, Nlink, dir_position
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Define position of virtual factor nodes
#-------------------------------------------------------------------------------
function virtuals(Nvar, dir_position)
    idx = findall(!iszero, dir_position)
    virtual = fill(1, Nvar)

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
    Nfac, Nvar, Ndir, Nlink,
    jacobianT, observation, noise, virtual,
    MEAN, VARI)

    Nind = Nfac - Ndir

    Mdir = fill(0.0, Nvar)
    Wdir = fill(0.0, Nvar)

    row = fill(0, Nlink)
    col = similar(row)
    coeff = fill(0.0, Nlink)
    coeffInv = similar(coeff)
    Mind = fill(0.0, Nind)
    Vind = similar(Mind)

    row_colptr = fill(0, Nind)
    col_colptr = fill(0, Nvar)

    idxi = 1
    idxr = 1

    idxT = findall(!iszero, jacobianT)
    @inbounds for i in idxT
        if (jacobianT.colptr[i[2] + 1] - jacobianT.colptr[i[2]]) == 1
            Mdir[i[1]] += observation[i[2]] * jacobianT[i] / noise[i[2]]
            Wdir[i[1]] += (jacobianT[i]^2) / noise[i[2]]
        else
            row[idxi] = idxr
            col[idxi] = i[1]

            coeff[idxi] = jacobianT[i]
            coeffInv[idxi] = 1 / jacobianT[i]

            idxi += 1

            col_colptr[i[1]] += 1
            row_colptr[idxr] = idxi

            if idxT[jacobianT.colptr[i[2] + 1] - 1] == i
                Mind[idxr] = observation[i[2]]
                Vind[idxr] = noise[i[2]]
                idxr += 1
            end
        end
        if virtual[i[1]] !== 0
            Mdir[i[1]] = MEAN / VARI
            Wdir[i[1]] = 1 / VARI
        end
    end

    pushfirst!(col_colptr, 1)
    pushfirst!(row_colptr, 1)
    col_colptr = cumsum(col_colptr)

    return row, row_colptr, col, col_colptr,
           Nind, Mind, Vind, coeff, coeffInv, Mdir, Wdir
end
#-------------------------------------------------------------------------------
