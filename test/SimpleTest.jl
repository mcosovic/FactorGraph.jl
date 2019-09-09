module SimpleTest

using GaussianBP
using SparseArrays
using HDF5
using Test

function model(DATA, PATH)
    system = string(PATH, DATA, ".h5")

    Hlist = h5read(system, "/H")
    H = sparse(Hlist[:,1], Hlist[:,2], Hlist[:,3])

    b = h5read(system, "/b")
    v = h5read(system, "/v")

    return H, b, v
end

function wlsMldivide(H, b, v)
    W = spdiagm(0 =>  @. 1.0 / sqrt(v))
    Hti = W * H
    rti = W * b

    xml = (Hti' * Hti) \ (Hti' * rti)

    return xml
end

H, b, v = model("SimpleTest", "test/")


bp("SimpleTest", 1000, 10, 0.6, 0.4, 0.0, 1e6, ALGORITHM = "sum", PATH = "test/") ≈ wlsMldivide(H, b, v)
# println(@test xbp ≈ wlsMldivide(H, b, v))
#
# xbp, H, b, v = bp("SimpleTest", 1000, 10, 0.6, 0.4, 0.0, 1e60, ALGORITHM = "kahan", PATH = "test/")
# println(@test xbp ≈ wlsMldivide(H, b, v))

end # SimpleTest
