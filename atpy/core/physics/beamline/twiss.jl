
_symplectic_conjugate_2by2(M)=return RealType[M[2,2] -M[1,2];-M[2,1] M[1,1]]

_matrixTransform_2by2(M)=begin
    m11,m21,m12,m22=M
    [m11*m11 -2m11*m12 m12*m12
    -m11*m21 1.0+2m12*m21 -m12*m22
    m21*m21 -2m21*m22 m22*m22]
end

"""
propagate the twiss function from begin to the end
period twiss function and local transfer matrixes from begin (findorbit) need be calculated before calling function
"""
function propagate_twiss!(lat::Lattice,tin::AbstractTwiss,beam::AbstractBeam=_beam[])
    inv_pi=1/2pi
	for i in 2:length(lat.Line)
        lat.Line[i].twiss.s=lat.Line[i-1].twiss.s+lat.Line[i].elem.L
        mux,muy=propagate_twiss!(lat.Line[i].twiss,tin,lat.Line[i].T )
        nux,nuy=inv_pi*mux,inv_pi*muy
        int_nux,frac_nux=fldmod(lat.Line[i-1].twiss.nux,1.0)
        lat.Line[i].twiss.nux = int_nux+ (frac_nux>nux+DoubleEPS ? 1.0+nux : nux)
        int_nuy,frac_nuy=fldmod(lat.Line[i-1].twiss.nuy,1.0)
        lat.Line[i].twiss.nuy = int_nuy+ (frac_nuy>nuy+DoubleEPS ? 1.0+nuy : nuy)
	end
end



function propagate_twiss!(tout::EdwardsTengTwiss,tin::EdwardsTengTwiss,M::Matrix{RealType})
	A=@view M[1:2,1:2]
	B=@view M[1:2,3:4]
	C=@view M[3:4,1:2]
	D=@view M[3:4,3:4]

	R1=reshape((@view tin.data[_k_R1:_k_R4]), (2,2) )

	_R1=_symplectic_conjugate_2by2(R1)
	if tin.mode == RealType(1)
		X=A-B*R1
		begin t=det(X)
			if t>RealType(0.1)
				R=(D*R1-C)*_symplectic_conjugate_2by2(X)
				R/=t
				X/=sqrt(t)
				Y=D+C*_R1
				Y/=sqrt(det(Y))
				mode=RealType(1)
			else
				X=C-D*R1
				X/=sqrt(det(X))
				Y=B+A*_R1
				t=det(Y)
				R=-(D+C*_R1)*_symplectic_conjugate_2by2(Y)
				R/=t
				Y/=sqrt(t)
				mode=RealType(2)
			end
		end
	elseif tin.mode == RealType(2) 
		X=B+A*_R1
		begin t=det(X)
			if t>RealType(0.1)
				R=-(D+C*_R1)*_symplectic_conjugate_2by2(X)
				R/=t
				X/=sqrt(t)
				Y=C-D*R1
				Y/=sqrt(det(Y))
				mode=RealType(1)
			else
				X=D+C*_R1
				X/=sqrt(det(X))
				Y=A-B*R1
				t=det(Y)
				R=(D*R1-C)*_symplectic_conjugate_2by2(Y)
				R/=t
				Y/=sqrt(t)
				mode=RealType(2)
			end
		end
	else
		#throw(AssertionError("Mode should be integer 1 or 2."))
		println(stderr,"Invalid mode.")
		return EdwardsTengTwiss(;betax=RealType(1),betay=RealType(1),mode=RealType(0))
	end
	tout.R1=R[1,1]
	tout.R2=R[1,2]
	tout.R3=R[2,1]
	tout.R4=R[2,2]
    return _decoupled_twiss_propagate!(tout,tin,X,Y,M)
end


@inline function _decoupled_twiss_propagate!(tout::AbstractTwiss,tin::AbstractTwiss,X::AbstractArray{RealType,2},Y::AbstractArray{RealType,2},M::AbstractArray{RealType,2})
    mux=_propagate_decoupled_twiss_2x2!( @view(tout.data[_k_betax:_k_gammax]), @view(tin.data[_k_betax:_k_gammax]), X )
    muy=_propagate_decoupled_twiss_2x2!( @view(tout.data[_k_betay:_k_gammay]), @view(tin.data[_k_betay:_k_gammay]), Y )
	(@view tout.data[_k_etax:_k_etapy]) .= (@view M[1:4,1:4])*(@view tin.data[_k_etax:_k_etapy]) .+ (@view M[1:4,6])
    return mux,muy
end

@inline _propagate_decoupled_twiss_2x2!(tout::AbstractArray{RealType,1},tin::AbstractArray{RealType,1},X::AbstractArray{RealType,2})=begin
    m11,m21,m12,m22=X
    tout[1] = m11*m11*tin[1] -      2m11*m12*tin[2] +  m12*m12*tin[3]
    tout[2] =-m11*m21*tin[1] +(1.0+2m12*m21)*tin[2] -  m12*m22*tin[3]
    tout[3] = m21*m21*tin[1] -      2m21*m22*tin[2] +  m22*m22*tin[3]
	sin_dmu=m12/sqrt(tout[1]*tin[1])
	cos_dmu=m11*sqrt(tin[1]/tout[1])-tin[2]*sin_dmu
    mu = mod(atan(sin_dmu,cos_dmu),2π )
    return mu
end


function periodic_EdwardsTengTwiss(M::Matrix{RealType};output_warning=true)
	A=@view M[1:2,1:2]
	B=@view M[1:2,3:4]
	C=@view M[3:4,1:2]
	D=@view M[3:4,3:4]
	invalid_ret=EdwardsTengTwiss(;betax=RealType(1),betay=RealType(1),mode=RealType(0))

	Bbar_and_C=_symplectic_conjugate_2by2(B)+C
	t1=0.5*(tr(A)-tr(D))
	Δ=t1*t1+det(Bbar_and_C)
	Δ<RealType(0) && begin
		if output_warning
			println(stderr,"Failed to decouple periodic transfer matrix. The linear matrix is unstable.")
		end
		return invalid_ret
	end

	_sign= t1>RealType(0) ? RealType(-1) : RealType(1)

	t2=abs(t1)+sqrt(Δ)
	if t2==RealType(0)
		R=RealType[0 0;0 0]
	else
		R=Bbar_and_C*(_sign/t2)
	end

	X=A-B*R
	Y=D+C*_symplectic_conjugate_2by2(R)

	# It should be equal to 1
	(det(X)<RealType(0.9) || det(Y)<RealType(0.9))  && begin
		if output_warning
			println(stderr,"Failed to decouple the periodic transfer matrix with mode 1.")
		end
		return invalid_ret
	end

	cmux=RealType(0.5)*(X[1,1]+X[2,2])
	cmuy=RealType(0.5)*(Y[1,1]+Y[2,2])
	(RealType(-1)<cmux<RealType(1) && RealType(-1)<cmuy<RealType(1)) || begin
		if output_warning
			println(stderr,"Failed to get beta functions. The linear matrix is unstable.")
		end
		return invalid_ret
	end

	smux=sqrt(RealType(1)-cmux*cmux)*sign(X[1,2])
	smuy=sqrt(RealType(1)-cmuy*cmuy)*sign(Y[1,2])
	betax=X[1,2]/smux
	gammax=-X[2,1]/smux
	betay=Y[1,2]/smuy
	gammay=-Y[2,1]/smuy

	alphax=RealType(0.5)*(X[1,1]-X[2,2])/smux
	alphay=RealType(0.5)*(Y[1,1]-Y[2,2])/smuy

	println("nux: ",atan(smux, cmux)/(2π), ", nuy: ", atan(smuy, cmuy)/(2π) )

	eta=inv(Matrix{RealType}(I,(4,4))-(@view M[1:4,1:4]))*(@view M[1:4,6])

	return EdwardsTengTwiss(s=0,betax=betax,alphax=alphax,nux=0,betay=betay,alphay=alphay,nuy=0,
                    etax=eta[1],etapx=eta[2],etay=eta[3],etapy=eta[4],R1=R[1,1],R2=R[1,2],R3=R[2,1],R4=R[2,2],mode=RealType(1))
end


function normalMatrix(tin::EdwardsTengTwiss)
	(tin.mode==RealType(1) || tin.mode==RealType(2)) || begin
		println(stderr,"Warning: return identity matrix for unknown mode $(tin.mode) as the normal matrix (transformation matrix from normal space to physical space).")
		return RealType(1)*Matrix{RealType}(I,6,6)
	end
	D=RealType[1 0 0 0 0 tin.etax
			   0 1 0 0 0 tin.etapx
			   0 0 1 0 0 tin.etay
			   0 0 0 1 0 tin.etapy
			   -tin.etapx tin.etax -tin.etapy tin.etay 1 0
			   0 0 0 0 0 1]
	sbx,sby=sqrt.((tin.betax,tin.betay) )
	B=RealType[sbx 0 0 0 0 0
			   -tin.alphax/sbx 1/sbx 0 0 0 0
			   0 0 sby 0 0 0
			   0 0 -tin.alphay/sby 1/sby 0 0
			   0 0 0 0 1 0
			   0 0 0 0 0 1]
    R=reshape((@view tin.data[_k_R1:_k_R4]), (2,2) )
	λ=RealType(1)/sqrt(abs(RealType(1)+det(R)))
	R=λ*R
	_R=_symplectic_conjugate_2by2(R)
	O=RealType[0 0;0 0]
	U=RealType[λ 0;0 λ]
	if tin.mode==RealType(1)
		V=[U _R O;-R U O;O O I]
	else
		V=[_R U O;U -R O;O O I]
	end
	return D*V*B
end
