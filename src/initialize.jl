#############################################################
#  Initialize the algorithm and produce damping parameters  #
#############################################################


#-------------------------------------------------------------------------------
# Pass messages from singly-connected factor nodes to all indirect links
#-------------------------------------------------------------------------------
function forward_directs_to_links(Mvar, Vvar, Mind, Vind, col, Nlink)
    @inbounds for i = 1:Nlink
        temp = 1 / Vind[col[i]]
        Mvar[i] = Mind[col[i]] * temp
        Vvar[i] = temp
    end

    return Mvar, Vvar
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
