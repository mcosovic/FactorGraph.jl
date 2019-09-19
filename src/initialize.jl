#############################################################
#  Initialize the algorithm and produce damping parameters  #
#############################################################


#-------------------------------------------------------------------------------
# Pass messages from singly-connected factor nodes to all indirect links
#-------------------------------------------------------------------------------
function forward_directs_to_links(Mvar_fac, Vvar_fac, Mdir, Wdir, Nvar, col_colptr, to_fac)
    @inbounds for i = 1:Nvar
        for j in col_colptr[i]:(col_colptr[i + 1] - 1)
            Vvar_fac[to_fac[j]] = 1 / Wdir[i]
            Mvar_fac[to_fac[j]] = Mdir[i] * Vvar_fac[to_fac[j]]
        end

    end

    return Mvar_fac, Vvar_fac
end
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Keep order of messages
#-------------------------------------------------------------------------------
function keep_order(Nlink, row, col)
    temp = collect(1:Nlink)
    sort_col_idx = sortperm(col)
    to_fac = temp[sort_col_idx]

    new_row = row[sort_col_idx]
    sort_row_idx = sortperm(new_row)
    to_var = temp[sort_row_idx]

    return to_fac, to_var
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
