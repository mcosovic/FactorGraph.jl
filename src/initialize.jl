#############################################################
#  Initialize the algorithm and produce damping parameters  #
#############################################################


#-------------------------------------------------------------------------------
# Pass messages from singly-connected factor nodes to all indirect links
#-------------------------------------------------------------------------------
function forward_directs_to_links(Mvar, Vvar, Mdir, VdirInv, Nvariable, variable_colptr, fv)
    @inbounds for i = 1:Nvariable
        for j in (variable_colptr[i]):variable_colptr[i+1]-1
            Vvar[fv[j]] = 1 / VdirInv[i]
            Mvar[fv[j]] = Mdir[i] * Vvar[fv[j]]
        end

    end

    return Mvar, Vvar
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Keep order of messages
#-------------------------------------------------------------------------------
function keep_order(Nlink, row, col)
    temp = collect(1:Nlink)
    sort_col_idx = sortperm(col)
    set_around_factor = temp[sort_col_idx]

    new_row = row[sort_col_idx]
    sort_row_idx = sortperm(new_row)
    set_around_variable = temp[sort_row_idx]

    return set_around_factor, set_around_variable
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Set damping parameters
#-------------------------------------------------------------------------------
function damping(Nlink, ALPH, PROB)
    bernoulli_sample = randsubseq(collect(1:Nlink), PROB)

    alpha1 = fill(1.0, Nlink)
    alpha1[bernoulli_sample] .= 1.0 - ALPH

    alpha2 = fill(0.0, Nlink)
    alpha2[bernoulli_sample] .= ALPH

    return alpha1, alpha2
end
#-------------------------------------------------------------------------------
