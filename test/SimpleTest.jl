module SimpleTest

using GaussianBP
using SparseArrays
using HDF5
using Test

function wlsMldivide(H, b, v)
    W = spdiagm(0 =>  @. 1.0 / sqrt(v))
    Hti = W * H
    rti = W * b

    xml = (Hti' * Hti) \ (Hti' * rti)

    return xml
end

xbp, H, b, v = bp("SimpleTest", 1000, 10, 0.6, 0.4, 0.0, 1e6, ALGORITHM = "sum", PATH = "test/")
println(@test xbp ≈ wlsMldivide(H, b, v))

xbp, H, b, v = bp("SimpleTest", 1000, 10, 0.6, 0.4, 0.0, 1e60, ALGORITHM = "kahan", PATH = "test/")
println(@test xbp ≈ wlsMldivide(H, b, v))

end # SimpleTest
