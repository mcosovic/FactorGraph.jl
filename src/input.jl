


#--------------------Load System Model from HDF5 File---------------------------
function model(data::String)
    system = string("src/data/", data, ".h5")

    Hlist = h5read(system, "/H")
    H = sparse(Hlist[:,1], Hlist[:,2], Hlist[:,3])

    b = h5read(system, "/b")
    v = h5read(system, "/v")

    return H, b, v
end
#-------------------------------------------------------------------------------
