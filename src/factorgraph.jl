


#-------------------------------------------------------------------------------
# Input Data:
#   - H: system model sparse matrix
# Output Data:
#   - Nf: number of factor nodes
#   - Nv: number of variable nodes
#   - T: system model transpose sparse matrix
#-------------------------------------------------------------------------------
function graph(H)
    Nf, Nv = size(H)
    T = SparseMatrixCSC(H')

    return Nf, Nv, T
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Input Data:
#   - Nf: number of factor nodes
#   - T: system model transpose sparse matrix
# Output Data:
#   - Nld: number of links between singly-connected factor and variable nodes
#   - Nli: number of links between indirect factor and variable nodes
#   - dir: position of each singly-connected factor according to variable nodes
#-------------------------------------------------------------------------------
function links(Nf, T)
    Nld = 0
    Nli = 0
    num = 0
    dir = fill(0, Nf)

    @inbounds for i = 1:Nf
        num = T.colptr[i + 1] - T.colptr[i]

        if num == 1
            Nld += num
            dir[i] = T.rowval[T.colptr[i]]
        else
            Nli += num
        end
    end

    return Nld, Nli, dir
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Input Data:
#   - Nf: number of factor nodes
#   - dir: position of each singly-connected factor according to variable nodes
# Output Data:
#   - vir: virtual factors according to variable nodes
#-------------------------------------------------------------------------------
function virtuals(Nv, dir)
    idx = findall(!iszero, dir)
    vir = fill(1, Nv)

    @inbounds for i in idx
        vir[dir[i]] = 0
    end

    return vir
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Input Data:
#   - Nf: number of factor nodes
#   - Nv: number of variable nodes
#   - Nld: number of links between singly-connected factor and variable nodes
#   - Nli: number of links between indirect factor and variable nodes
#   - T: system model transpose sparse matrix
#   - b: mean vector of measurements
#   - v: variance vector of measurements
#   - vir: virtual factors according to variable nodes
#   - MEAN: virtual factor node mean value
#   - VARI: virtual factor node variance value
# Output Data:
#   - Hi: vector of coefficient of indirect factor nodes
#   - bi: vector of measurement means of indecies factor nodes
#   - vi: vector of measurement variances of indirect factor nodes
#   - Ii: indices of indirect factors according to factor nodes (row indices)
#   - Ji: indices of indirect factors according to variable nodes (column indices)
#   - Ni: number of indirect factor nodes
#   - md: total means from singly-connected factors to variable nodes
#   - vid: total inverse variance from singly-connected factors to variable nodes
#-------------------------------------------------------------------------------
function factors(Nf, Nv, Nld, Nli, T, b, v, vir, MEAN, VARI)
    Ni = Nf - Nld

    md = fill(0.0, Nv)
    vid = fill(0.0, Nv)

    Ii = fill(0, Nli)
    Ji = similar(Ii)
    Hi = fill(0.0, Nli)
    bi = fill(0.0, Ni)
    vi = similar(bi)

    idxi = 1
    idxr = 1

    idxT = findall(!iszero, T)
    @inbounds for i in idxT
        if (T.colptr[i[2] + 1] - T.colptr[i[2]]) == 1
            md[i[1]] += b[i[2]] * T[i] / v[i[2]]
            vid[i[1]] += (T[i]^2) / v[i[2]]
        else
            Ii[idxi] = idxr
            Ji[idxi] = i[1]
            Hi[idxi] = T[i]
            idxi += 1
            if idxT[T.colptr[i[2] + 1] - 1] == i
                bi[idxr] = b[i[2]]
                vi[idxr] = v[i[2]]
                idxr += 1
            end
        end
        if vir[i[1]] !== 0
            md[i[1]] = MEAN
            vid[i[1]] = 1 / VARI
        end
    end

    return Ii, Ji, Ni, bi, vi, Hi, md, vid
end
#-------------------------------------------------------------------------------
