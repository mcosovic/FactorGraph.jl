module Instability

using SparseArrays
using HDF5
using Test
using Random

################################################################################
# Read data from SimpleTest.h5
#-------------------------------------------------------------------------------
function model(DATA, PATH)
    system = string(PATH, DATA, ".h5")

    Hlist = h5read(system, "/H")
    H = sparse(Hlist[:,1], Hlist[:,2], Hlist[:,3])

    b = h5read(system, "/b")
    v = h5read(system, "/v")

    return H, b, v
end
################################################################################

H, b, v = model("SimpleTest", "test/")

include("../src/factorgraph.jl")
Nf, Nv, T = @inferred graph(H)
Nld, Nli, dir =  @inferred links(Nf, T)
vir = @inferred virtuals(Nv, dir)                                                                                         # factorgraph.jl
Ii, Ji, Ni, bi, vi, Hi, md, vid = @inferred factors(Nf, Nv, Nld, Nli, T, b, v, vir, 0.0, 1e3)                           # factorgraph.jl

include("../src/auxiliary.jl")
m_fv, vi_fv, m_vf, v_vf = @inferred load_messages(Hi)
msr, vsr, msc, vsc = @inferred load_sum(Nv, Ni)
msr, vsr, msc, vsc = @inferred clear_sum(msr, vsr, msc, vsc)
msr, vsr, evr, msc, vsc, evc = @inferred nload_sum(Nv, Ni)
msr, vsr, evr, msc, vsc, evc = @inferred nclear_sum(msr, vsr, evr, msc, vsc, evc)                                             # AuxiliaryFunction

include("../src/initialize.jl")
ah1, ah2 =  @inferred damping(Nli, 0.5, 0.5)                                                                                 # InitializeMessages
m_vf, v_vf =  @inferred forward_directs(Hi, Ji, Nli, md, vid, v_vf, m_vf)


end
