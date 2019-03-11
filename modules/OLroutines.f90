module OLroutines
	use MO_module, only: dp, qp, ten_quad
	implicit none
	!real(kind=qp), parameter :: ten_quad = 10._qp
	real(8), parameter :: PI=4.D0*DATAN(1.D0)
	private 
	public sum_check,prod_calc,WTW,ten_quad,PI,M_create,M_flip
contains

! 	Calculates the sum of all v_i² ,i=1,...,N - should be the same as <N>
	FUNCTION sum_check(V,N) result(summ)
		COMPLEX(dp) 		:: V(N,N) 
		integer 		:: N,i
		!COMPLEX(kind=qp) 	:: prod
		real(kind=qp) 	:: summ
		
		summ=0
		do i=1,N
			if(i+1<=N) then
				!prod=prod*V(i,i+1)*V(i,i+1)
			!	write(*,*) V(i,i+1)
				summ=summ+V(i,i+1)**2
			end if
		end do
	END FUNCTION sum_check

! 	Calculates the prod of v_i² , i=1,..,N/2, since we're using paired states
	FUNCTION prod_calc(V,N) result(prod)
		COMPLEX(dp) 		:: V(N,N) 
		integer 		:: N,i
		!COMPLEX(kind=qp) 	:: prod
		real(kind=qp) 	:: prod
		
		prod=1
		do i=1,N,2
			if(i+1<=N) then
				!prod=prod*V(i,i+1)*V(i,i+1)
				prod=prod*V(i,i+1)**2
			end if
		end do
	END FUNCTION prod_calc

! 	Creates the matrix of <phi
	FUNCTION WTW(U,V,N) result(WW)
		COMPLEX(dp) 	:: WW(2*N,2*N),U(N,N),V(N,N)
		integer 	:: N 

		WW(1:N,1:N) = matmul(transpose(V),U)
		WW(N+1:2*N,1:N) = -matmul(transpose(V),V)
		WW(1:N,N+1:2*N) = matmul(transpose(V),V)
		WW(N+1:2*N,N+1:2*N) = matmul(transpose(U),V)
	END FUNCTION WTW

	function M_create(UD,UM,VD,VM,DMAT,N) 
		complex(kind=dp),dimension(N,N), intent(in) 	:: DMAT, UD, VD, UM, VM
		complex(dp) 					:: M_create(2*N,2*N)
		integer 	:: N 

		M_create(1:N,1:N) = matmul(transpose(VM),UM)
		M_create(1:N,(N+1):2*N) = matmul(transpose(VM),matmul(DMAT,conjg(VD)))
		M_create((N+1):(2*N),1:N) = -1*matmul(transpose(conjg(VD)),matmul(transpose(DMAT),VM))
		M_create((N+1):(2*N),(N+1):(2*N)) = matmul(transpose(conjg(UD)),conjg(VD))
	end function M_create

	function M_flip(UD,UM,VD,VM,DMAT,N) 
		complex(kind=dp),dimension(N,N), intent(in) 	:: DMAT, UD, VD, UM, VM
		complex(dp) 					:: M_flip(2*N,2*N)
		integer 	:: N 

		M_flip(1:N,1:N) = matmul(transpose(VM),UM)
		M_flip(1:N,(N+1):2*N) = matmul(transpose(VM),matmul(transpose(conjg(DMAT)),conjg(VD)))
		M_flip((N+1):(2*N),1:N) = -1*matmul(transpose(conjg(VD)),matmul(conjg(DMAT),VM))
		M_flip((N+1):(2*N),(N+1):(2*N)) = matmul(transpose(conjg(UD)),conjg(VD))
	end function M_flip
! ------------------------------------------------ Superflous routines -----------------------------------------------

	FUNCTION DMAT_CREATOR(phi,N) result(DMAT)
		! computes the overlap of two wave functions \bra{UM,VM} DMAT \ket{UD,VD}   with ZSKPFA from PFAPACK,
		implicit none

		integer, intent(in):: N
		real(dp), intent(in) :: phi
		
		COMPLEX(dp), dimension(N,N) :: DMAT

		COMPLEX(dp) :: z_no
		integer :: i

		z_no = cmplx(0,phi,16)
		DMAT=0

		do i=1,N
			DMAT(i,i) = exp(z_no)
		end do 
	END FUNCTION DMAT_CREATOR

	FUNCTION norm_fac(U,V,N) result(norm)
		COMPLEX(dp) 	:: WW(2*N,2*N),U(N,N),V(N,N)
		integer 	:: N 
		COMPLEX(kind=dp)	:: 	Pf2P(2)
		REAL(kind=qp) 		:: 	norm

		WW=WTW(U,V,N)
		call ZSKPF10_F95(WW,Pf2P)
		norm=(-1)**(N*(N-1)/2)*real(PF2P(1))*ten_quad**real(PF2P(2))

	END FUNCTION norm_fac
	
	FUNCTION prod_sqrt(V,N) result(prod)
		COMPLEX(dp) 		:: V(N,N) 
		integer 		:: N,i
		COMPLEX(kind=qp) 	:: prod
		
		prod=1
		do i=1,N,2
			if(i+1<=N) then
				!prod=prod*V(i,i+1)*V(i,i+1)
				prod=prod*V(i,i+1)
			end if
		end do
	END FUNCTION prod_sqrt

end module

