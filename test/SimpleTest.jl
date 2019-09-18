module SimpleTest

using GaussianBP
using SparseArrays
using HDF5
using Test


#########################################################################
#  The belief propagation accuracy test, where the algorithm converged  #
#########################################################################


##################################
#  Read data from SimpleTest.h5  #
##################################
function model(DATA, PATH)
    system = string(PATH, DATA, ".h5")

    list = h5read(system, "/H")
    jacobian = sparse(list[:,1], list[:,2], list[:,3])

    observation = h5read(system, "/b")
    noise = h5read(system, "/v")

    return jacobian, observation, noise
end


#################################################
#  Compute the weighted least-squares solution  #
#################################################
function wls_mldivide(jacobian, observation, noise)
    W = spdiagm(0 =>  @. 1.0 / sqrt(noise))
    Hti = W * jacobian
    rti = W * observation

    Xml = (Hti' * Hti) \ (Hti' * rti)

    return Xml
end


##############################################
#  Test for the belief propagation accuracy  #
##############################################
jacobian, observation, noise = model("SimpleTest", "test/")

@testset "SimplyBP" begin
    @test bp("SimpleTest", 1000, 10, 50, 0.6, 0.4, 0.0, 1e6, ALGORITHM = "sum", PATH = "test/") ≈ wls_mldivide(jacobian, observation, noise)
    @test bp("SimpleTest", 1000, 50, 80, 0.0, 0.4, 0.0, 1e6, ALGORITHM = "sum", PATH = "test/") ≈ wls_mldivide(jacobian, observation, noise)
    @test bp("SimpleTest", 1000, 10, 80, 0.2, 0.3, 10.0, 1e8, ALGORITHM = "sum", PATH = "test/") ≈ wls_mldivide(jacobian, observation, noise)
end

@testset "KahanBP" begin
    @test bp("SimpleTest", 1000, 10, 50, 0.6, 0.4, 0.0, 1e30, ALGORITHM = "kahan", PATH = "test/") ≈ wls_mldivide(jacobian, observation, noise)
    @test bp("SimpleTest", 1000, 50, 50, 0.0, 0.0, 0.0, 1e60, ALGORITHM = "kahan", PATH = "test/") ≈ wls_mldivide(jacobian, observation, noise)
    @test bp("SimpleTest", 500, 10, 50, 0.6, 0.4, 10.0, 1e80, ALGORITHM = "kahan", PATH = "test/") ≈ wls_mldivide(jacobian, observation, noise)
end

end # SimpleTest
