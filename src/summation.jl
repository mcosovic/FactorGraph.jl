##########################
#  Summarizing messages  #
##########################


#--------------------------------------------------------------------------
# Summarizing messages from variable nodes to factor nodes using simply
# summation, or using Kahan-Babuska algorithm
# Input Data:
#	- Hi: vector of coefficient of indirect factor nodes
#	- Ii: indices of indirect factors according to factor nodes (rows)
#	- Nli: number of links between indirect factor and variable nodes
#	- m_vf: mean messages from variable node to factor node
#	- v_vf: variance messages from variable node to factor node
#	- msr: sum vector of row mean messages
#	- vsr: sum vector of row variance messages
#	- evr: error vector of variance summation
# Output Data:
#	- msr: sum vector of row mean messages
#	- vsr: sum vector of row variance messages
#--------------------------------------------------------------------------
function sum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr)
	@inbounds for i = 1:Nli
		msr[Ii[i]] += (Hi[i] * m_vf[i])
		vsr[Ii[i]] += (Hi[i]^2 * v_vf[i])
	end

	return msr, vsr
end

function nsum_rows(Hi, Ii, Nli, m_vf, v_vf, msr, vsr, evr)
	@inbounds for i = 1:Nli
		msr[Ii[i]] += (Hi[i] * m_vf[i])

		x = Hi[i]^2 * v_vf[i]
		t = vsr[Ii[i]] + x
		if abs(vsr[Ii[i]]) >= abs(x)
			evr[Ii[i]] += (vsr[Ii[i]] - t) + x
		else
			evr[Ii[i]] += (x - t) + vsr[Ii[i]]
		end
		vsr[Ii[i]] = t
	end

	return msr, vsr, evr
end
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
# Summarizing messages from factor nodes to variable nodes using simply
# summation, or using Kahan-Babuska algorithm
# Input Data:
#	- Ji: indices of indirect factors according to variable nodes (columns)
#	- Nli: number of links between indirect factor and variable nodes
#	- m_fv: mean messages from factor node to variable node
#	- vi_fv: inverse variance messages from factor node to variable node
#	- msc: sum vector of column mean messages
#	- vsc: sum vector of column variance messages
#	- evc: error vector of variance summation
# Output Data:
#	- msc: sum vector of column mean messages
#	- vsc: sum vector of column variance messages
#--------------------------------------------------------------------------
function sum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc)
	@inbounds for i = 1:Nli
		msc[Ji[i]] += m_fv[i] * vi_fv[i]
		vsc[Ji[i]] += vi_fv[i]
	end

	return msc, vsc
end

function nsum_cols(Ji, Nli, m_fv, vi_fv, msc, vsc, evc)
	@inbounds for i = 1:Nli
		msc[Ji[i]] += m_fv[i] * vi_fv[i]

		t = vsc[Ji[i]] + vi_fv[i]
		if abs(vsc[Ji[i]]) >= abs(vi_fv[i])
			evc[Ji[i]] += (vsc[Ji[i]] - t) + vi_fv[i]
		else
			evc[Ji[i]] += (vi_fv[i] - t) + vsc[Ji[i]]
		end
		vsc[Ji[i]] = t
	end

	return msc, vsc, evc
end
#--------------------------------------------------------------------------
