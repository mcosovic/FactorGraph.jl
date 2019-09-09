module SimpleTest

using GaussianBP
using SparseArrays
using HDF5
using Test

bp("SimpleTest", 20, 10, 0.6, 0.4, 0.0, 1e3, ALGORITHM = "sum", ERROR = "on", PATH = "test/")


# function wlsMldivide(H, b, v)
#     W = spdiagm(0 =>  @. 1.0 / sqrt(v))
#     Hti = W * H
#     rti = W * b
#
#     xml = (Hti' * Hti) \ (Hti' * rti)
#
#     return xml
# end



# xml = wlsMldivide(H, b, v)
#
# a = @test bp("SimpleTest", 1000, 10, 0.6, 0.4, 0.0, 1e6)  ≈ xml
# b = @test bpn("SimpleTest", 1000, 10, 0.6, 0.4, 0.0, 1e60)  ≈ xml
# println(a)
# println(b)



end # SimpleTest
