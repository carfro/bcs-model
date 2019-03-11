FUNCTION WD(U,V,DMAT,N) result(WW)
	COMPLEX(dp) 	:: WW(2*N,2*N),U(N,N),V(N,N)
	integer 	:: N 

	WTW(1:N,1:N) = matmul(transpose(VM),UM)
	WTW(1:N,(N+1):2*N) = matmul(transpose(VM),matmul(DMAT,conjg(VD)))
	WTW((N+1):(2*N),1:N) = -1*matmul(transpose(conjg(VD)),matmul(transpose(DMAT),VM))
	WTW((N+1):(2*N),(N+1):(2*N)) = matmul(transpose(conjg(UD)),conjg(VD))
END FUNCTION WD

program main
	


end program main
