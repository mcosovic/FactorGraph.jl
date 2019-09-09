module SimpleTest

using GaussianBP
using SparseArrays
using HDF5
using Test


################################################################################
# The belief propagation accuracy test, where the algorithm converged
################################################################################


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


################################################################################
# Compute the weighted least-squares solution
#-------------------------------------------------------------------------------
function wlsMldivide(H, b, v)
    W = spdiagm(0 =>  @. 1.0 / sqrt(v))
    Hti = W * H
    rti = W * b

    xml = (Hti' * Hti) \ (Hti' * rti)

    return xml
end
################################################################################


################################################################################
# Test the belief propagation accuracy
#-------------------------------------------------------------------------------
H, b, v = model("SimpleTest", "test/")

@testset "SimplyBP" begin
    @test bp("SimpleTest", 1000, 10, 0.6, 0.4, 0.0, 1e6, ALGORITHM = "sum", PATH = "test/") ≈ wlsMldivide(H, b, v)
    @test bp("SimpleTest", 1000, 50, 0.0, 0.4, 0.0, 1e6, ALGORITHM = "sum", PATH = "test/") ≈ wlsMldivide(H, b, v)
    @test bp("SimpleTest", 1000, 10, 0.6, 0.4, 10.0, 1e8, ALGORITHM = "sum", PATH = "test/") ≈ wlsMldivide(H, b, v)
end

@testset "KahanBP" begin
    @test bp("SimpleTest", 1000, 10, 0.6, 0.4, 0.0, 1e30, ALGORITHM = "kahan", PATH = "test/") ≈ wlsMldivide(H, b, v)
    @test bp("SimpleTest", 1000, 50, 0.0, 0.0, 0.0, 1e60, ALGORITHM = "kahan", PATH = "test/") ≈ wlsMldivide(H, b, v)
    @test bp("SimpleTest", 500, 10, 0.6, 0.4, 10.0, 1e80, ALGORITHM = "kahan", PATH = "test/") ≈ wlsMldivide(H, b, v)
end
################################################################################


end # SimpleTest
