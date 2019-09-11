#########################################################
#  Auxiliary functions for arrays loading and clearing  #
#########################################################


#------------------------------------------------------------------------
# Load arrays for messages
#------------------------------------------------------------------------
function load_messages(Hi)
    m_fv = similar(Hi)
    vi_fv = similar(Hi)
    m_vf = similar(Hi)
    v_vf = similar(Hi)

    return m_fv, vi_fv, m_vf, v_vf
end
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# Load arrays for message summation in ALGORITHM = "sum" and "kahan"
#------------------------------------------------------------------------
function load_sum(Nv, Ni)
    msr = fill(0.0, Ni)
    vsr = fill(0.0, Ni)
    msc = fill(0.0, Nv)
    vsc = fill(0.0, Nv)

    return msr, vsr, msc, vsc
end

function nload_sum(Nv, Ni)
    msr = fill(0.0, Ni)
    vsr = fill(0.0, Ni)
    evr = fill(0.0, Ni)
    msc = fill(0.0, Nv)
    vsc = fill(0.0, Nv)
    evc = fill(0.0, Nv)

    return msr, vsr, evr, msc, vsc, evc
end
#------------------------------------------------------------------------
