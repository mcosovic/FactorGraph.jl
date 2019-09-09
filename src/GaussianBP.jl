module GaussianBP

export runbp

using SparseArrays
using HDF5
using Random


######################
# Default Parameters #
######################

const data = "data_33_14"

const MAXI = 15
const DAMP = 10

const PROB = 0.6
const ALPH = 0.4

const MEAN = 0.0
const VARI = 1e3


############
# Includes #
############

include("input.jl")
include("factorgraph.jl")
include("auxiliary.jl")
include("initialize.jl")
include("inference.jl")
include("summation.jl")


######################
# Simply Gaussian BP #
######################

function runbp(data, MAXI, DAMP, PROB, ALPH, MEAN, VARI)
    H, b, v = model(data)                                                                                      # input.jl

    Nf, Nv, T = graph(H)                                                                                       # factorgraph.jl
    Nld, Nli, dir = links(Nf, T)                                                                               # factorgraph.jl
    vir = virtuals(Nv, dir)                                                                                    # factorgraph.jl
    Ii, Ji, Ni, bi, vi, Hi, md, vid = factors(Nf, Nv, Nld, Nli, T, b, v, vir, MEAN, VARI)                      # factorgraph.jl

    m_fv, vi_fv, m_vf, v_vf = load_messages(Hi)                                                                # auxiliary.jl
    ah1, ah2 = damping(Nli, ALPH, PROB)                                                                        # initialize.jl
    m_vf, v_vf = forward_directs(Hi, Ji, Nli, md, vid, v_vf, m_vf)                                             # initialize.jl

    msr, vsr, msc, vsc = load_sum(Nv, Ni)                                                                      # auxiliary.jl
    msr, vsr = sum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr)                                                     # summation.jl
    m_fv, vi_fv = factor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, Hi, bi, vi, Ii, Nli)                   # inference.jl
    msr, vsr, msc, vsc = clear_sum(msr, vsr, msc, vsc)                                                         # auxiliary.jl

for i = 1:MAXI
    msc, vsc = sum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc)                                                         # summation.jl
    m_vf, v_vf = variable_to_factor(m_vf, v_vf, m_fv, vi_fv, md, vid, msc, vsc, Ji, Nli)                        # inference.jl
    msr, vsr = sum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr)                                                      # summation.jl
    if i < DAMP
        m_fv, vi_fv = factor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, Hi, bi, vi, Ii, Nli)                # inference.jl
    else
        m_fv, vi_fv = dfactor_to_variable(m_vf, v_vf, m_fv, vi_fv, msr, vsr, Hi, bi, vi, Ii, Nli, ah1, ah2)     # inference.jl
    end
    msr, vsr, msc, vsc = clear_sum(msr, vsr, msc, vsc)                                                          # auxiliary.jl
end

    msc, vsc = sum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc)                                                         # summation.jl
    xbp = marginal(md, vid, msc, vsc, Ji, Nv)                                                                   # inference.jl

    return xbp
 end

end # SimplyGBP
