module part_no_module
	use mroutines
	use MO_module
	use analytic_module
	use pfaffian_module
	use F95_PFAPACK
	implicit none
	!integer, parameter :: N_tot=330
	real(8), parameter :: PI=4.D0*DATAN(1.D0)
	real(kind=qp), parameter :: ten_quad = 10._qp

contains

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
	
	FUNCTION WTW(U,V,N) result(WW)
		COMPLEX(dp) 	:: WW(2*N,2*N),U(N,N),V(N,N)
		integer 	:: N 

		WW(1:N,1:N) = matmul(transpose(V),U)
		WW(N+1:2*N,1:N) = -matmul(transpose(V),V)
		WW(1:N,N+1:2*N) = matmul(transpose(V),V)
		WW(N+1:2*N,N+1:2*N) = matmul(transpose(U),V)
	END FUNCTION WTW

	FUNCTION MAT_NORM2(N,WW) result(NORM)
		complex(dp) 	:: WW(N,N)
		real(dp) 	:: NORM
		integer 	:: N,i,j

		NORM=0
		WW=matmul(conjg(WW),WW)
		do i=1,N
			NORM=NORM+WW(i,i)	
		end do

		NORM=sqrt(NORM)
	END FUNCTION MAT_NORM2

	function M_create(UD,UM,VD,VM,DMAT,N) 
		complex(kind=dp),dimension(N,N), intent(in) 	:: DMAT, UD, VD, UM, VM
		complex(dp) 					:: M_create(2*N,2*N)
		integer 	:: N 

		M_create(1:N,1:N) = matmul(transpose(VM),UM)
		M_create(1:N,(N+1):2*N) = matmul(transpose(VM),matmul(DMAT,conjg(VD)))
		M_create((N+1):(2*N),1:N) = -1*matmul(transpose(conjg(VD)),matmul(transpose(DMAT),VM))
		M_create((N+1):(2*N),(N+1):(2*N)) = matmul(transpose(conjg(UD)),conjg(VD))
	end function M_create

	complex(qp) function spart(N,UD,UM,VD,VM,DMAT,OP) result(trace)
		integer, intent(in) 				:: N
		complex(kind=dp),dimension(N,N), intent(in) 	:: DMAT, UD, VD, UM, VM, OP
		complex(kind=dp),dimension(N,N) 	 	:: X
		integer 					:: i
		
		X = matmul(transpose(UD),matmul(DMAT,UM))+&
			matmul(transpose(VD),matmul(conjg(DMAT),VM))

		X = transpose(inv(N,X))
		!X = inv(N,X)

		X = matmul(matmul(conjg(VM),X),matmul(transpose(VD),transpose(conjg(DMAT))))
		X = matmul(OP,X)
		
		trace = 0	

		do i=1,N
			trace=trace+X(i,i)
		end do

	end function

	function ext_pfaffian(N,UD,UM,VD,VM,DMAT,get_norm) result(pf)
	    ! computes the overlap of two wave functions with pfaffian formula
		integer, intent(in) 				:: N
		COMPLEX(kind=dp),dimension(N,N), intent(in) 	:: DMAT, UD, VD, UM, VM
		REAL(kind=qp), intent(in), optional		:: get_norm

		COMPLEX(kind=dp) 	:: 	WTW(2*N,2*N)
		COMPLEX(kind=qp)	:: 	pf
		REAL(kind=qp) 		:: 	norm
		integer :: IPIV(2*N,2)
		integer 		:: 	i,j,S_n

		S_n=(-1)**(N*(N-1)/2)

		if( present(get_norm)) then
			norm=get_norm
		else
			norm=prod_calc(VD,N)
		end if

		WTW(1:N,1:N) = matmul(transpose(VM),UM)
		WTW(1:N,(N+1):2*N) = matmul(transpose(VM),matmul(DMAT,conjg(VD)))
		WTW((N+1):(2*N),1:N) = -1*matmul(transpose(conjg(VD)),matmul(transpose(DMAT),VM))
		WTW((N+1):(2*N),(N+1):(2*N)) = matmul(transpose(conjg(UD)),conjg(VD))

		call Zpfaffian_ext(WTW,2*N,2*N,IPIV,pf)
		
		pf=S_n*pf/(norm)

	    return
	end function ext_pfaffian

	SUBROUTINE ZSKPF10_OLP(Pf2P_out,N,UD,UM,VD,VM,DMAT,get_norm) 
		! computes the overlap of two wave functions \bra{UM,VM} DMAT \ket{UD,VD}   with ZSKPFA from PFAPACK,
		implicit none

		integer, intent(in) 				:: N
		COMPLEX(kind=dp),dimension(N,N), intent(in) 	:: DMAT, UD, VD, UM, VM
		REAL(kind=qp), intent(in), optional		:: get_norm

		COMPLEX(kind=dp) 	:: 	WTW(2*N,2*N)
		COMPLEX(kind=dp)	:: 	Pf2P(2),Pf2P_out
		REAL(kind=qp) 		:: 	norm
		integer 		:: 	i,j,S_n

		S_n=(-1)**(N*(N-1)/2)

		if( present(get_norm)) then
			norm=get_norm
		else
			norm=prod_calc(VD,N)
		end if

		WTW(1:N,1:N) = matmul(transpose(VM),UM)
		WTW(1:N,(N+1):2*N) = matmul(transpose(VM),matmul(DMAT,conjg(VD)))
		WTW((N+1):(2*N),1:N) = -1*matmul(transpose(conjg(VD)),matmul(transpose(DMAT),VM))
		WTW((N+1):(2*N),(N+1):(2*N)) = matmul(transpose(conjg(UD)),conjg(VD))

		!IF(.NOT. (MAT_NORM2(2*N,WTW+transpose(WTW))<ten_quad**(-4))) write(*,*) MAT_NORM2(N,WTW+transpose(WTW))

		call ZSKPF10_F95(WTW,Pf2P) 

		Pf2P_out=S_n*cmplx(real(PF2P(1)),aimag(PF2P(1)),kind=qp)*ten_quad**(real(PF2P(2)))/(norm)

	END SUBROUTINE ZSKPF10_OLP

	SUBROUTINE ZSKPF10_MOLP(Pf2P_out,N,WDW,norm) 
		! computes the overlap of two wave functions \bra{UM,VM} DMAT \ket{UD,VD}   with ZSKPFA from PFAPACK,
		implicit none

		integer, intent(in) 					:: N
		COMPLEX(kind=dp),dimension(2*N,2*N), intent(in) 	:: WDW
		REAL(kind=qp), intent(in) 				:: norm

		COMPLEX(kind=dp)	:: 	Pf2P(2),Pf2P_out
		COMPLEX(kind=dp) 	:: 	WTW(2*N,2*N)
		integer 		:: 	i,j,S_n

		WTW=WDW
		S_n=(-1)**(N*(N-1)/2)

		!IF(.NOT. (MAT_NORM2(2*N,WTW+transpose(WTW))<ten_quad**(-4))) write(*,*) MAT_NORM2(N,WTW+transpose(WTW))

		call ZSKPF10_F95(WTW,Pf2P) 

		Pf2P_out=S_n*cmplx(real(PF2P(1)),aimag(PF2P(1)),kind=qp)*ten_quad**(real(PF2P(2)))/(norm)

	END SUBROUTINE ZSKPF10_MOLP

	SUBROUTINE ONISHI_OVERLAP(overlap_out,N,UD,UM,VD,VM,DMAT,phi) 
		! computes the overlap of two wave functions \bra{UM,VM} DMAT \ket{UD,VD}   with ZSKPFA from PFAPACK,
		implicit none

		integer, intent(in) :: N
		COMPLEX(kind=dp),dimension(N,N), intent(in) :: DMAT,UD, VD, UM, VM
		real(kind=dp), intent(in) :: phi

		COMPLEX(kind=qp), intent(out) :: overlap_out

		COMPLEX(kind=qp),dimension(N,N)	:: 	X
		COMPLEX(kind=qp) 	:: 	det_X,det_test,expan_factor
		integer 		:: 	i,j,S_n

		X = matmul(transpose(conjg(UD)),matmul(transpose(conjg(DMAT)),UM))+&
			matmul(transpose(conjg(VD)),matmul(transpose(DMAT),VM))

		!expan_factor=minval(real(X))
		!expan_factor=1._qp
		expan_factor=ten_quad**(2)
		X=X*expan_factor

		!det_X = det(N,X)/(expan_factor**(N))
		!write(*,*) 'LU-dec det 	', det_X
		
		det_test = det_eig(N,X)/(expan_factor**(N))	
		!write(*,*) 
		!write(*,*) 'Eig-dec det	 ', det_test

	!	do i=1,N
	!		write(7,'(100F41.30)') (real(X(i,j)),j=1,N)
	!	end do 	
	!	write(8,*) 

		!write(3,'(3F12.7)') phi,abs(sqrt(det_X)),abs(sqrt(det_test))

		overlap_out=sqrt(det_test)*exp(cmplx(0,-0.5*N*phi,kind=qp))

	END SUBROUTINE ONISHI_OVERLAP

	
end module part_no_module
! Constructs the BCS-model and computes the system for \lambda=[1..-1] and plots the result versus number of particles
program main
	use part_no_module
	implicit none
! 	Nucleus, filled using MO_module
	type(Nucleon), dimension(:,:), allocatable :: nucleus

! 	Variables used for the particle projection
	COMPLEX(kind=dp) :: 	Pf2P(2),Pf2P_ol,Pf2P_ol2
	COMPLEX(kind=qp) ::	on_ol,ext_ol,prod_N,prod_Z,&
		 		summ,summ2,summ3,TrcRhoR,&
				Nexp,Nexp2

	REAL(kind=qp) 	:: 	norm
	REAL(kind=dp) 	:: 	dPhi,N_exp,N_op,tol,lam_sN,lam_sZ,scaleFactor(7),deltaR,phi_arg,phi_arg2

	COMPLEX(kind=dp), dimension(:,:), allocatable :: 	U_N,V_N,U_Z,V_Z,DMAT,Mmat,EYE
	COMPLEX(kind=dp), dimension(:), allocatable   :: 	ROL,ROL2,ROL3,NROL,NROL2
	REAL(kind=dp), dimension(:,:), allocatable    :: 	EV_N,EV_Z


	! Nbr of neutrons/protons and loop integer(s),
	integer :: i,j,N_loop,k,l
	integer :: N,Z,step

	integer, parameter :: cut = 100

	N=24	! Number of NEUTRONS to find
	Z=24 	! Number of PROTONS to find
	deltaR=12.d0/sqrt(real(N+Z))
	
	scaleFactor = (/0.5*deltaR,1*deltaR,1.5*deltaR,3*deltaR,6*deltaR,12*deltaR,24*deltaR/)

!------------------------- Test-loop for the particle nbr op exp value ----------------------------------

	allocate(nucleus(cut,2),U_N(cut,cut),V_N(cut,cut),U_Z(cut,cut),&
		V_Z(cut,cut),DMAT(cut,cut),Mmat(2*cut,2*cut),EYE(cut,cut))

	call nucleus_creator(N,Z,cut,nucleus)
	call qpart_creator(nucleus,N,Z,cut,scaleFactor(1),U_N,V_N,prod_N,U_Z,V_Z,prod_Z)

! 	overlap testing files
	open(unit=2,file='data/part_no/overlap_proj_test.dat',status='replace')
	open(unit=3,file='data/part_no/Nexp.dat',status='replace')
	open(unit=9,file='data/part_no/overlap_real.dat',status='replace')
	open(unit=10,file='data/part_no/overlap_imag.dat',status='replace')
	open(unit=11,file='data/part_no/sum.dat',status='replace')

	norm=prod_calc(V_N,cut)
	write(*,*) 'norm 	', norm

	N_loop=10!*20 	!5*30*30
	allocate(ROL(cut),ROL2(cut),ROL3(cut),NROL(cut),NROL2(cut))
	
	N_exp=2*sum_check(V_N,cut)
	write(*,*) 'Nbr of part 	', N_exp

	dPhi=2*PI/(2*N_loop+1)	

	Mmat=0
	DMAT=0
	EYE=0

	Pf2P_ol=0
	ext_ol=0
	on_ol=0

	forall(i=1:cut) EYE(i,i)=cmplx(1,0,16)

	do j=ceiling(N_exp)-6,ceiling(N_exp)+6
		write(*,*) j
		ROL =0
		ROL2=0
		ROL3=0

		!!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(DMAT,EYE,U_N,V_N,norm,N_loop,phi_arg)
		do i=-1*N_loop,N_loop
			!write(*,*) i!, i*dPhi

			DMAT=EYE*exp(cmplx(0,i*dPhi,16))

			Mmat=M_create(U_N,U_N,V_N,V_N,DMAT,cut)

			call ZSKPF10_MOLP(Pf2P_ol,cut,Mmat,norm)

			ext_ol=ext_pfaffian(cut,U_N,U_N,V_N,V_N,DMAT,norm)
			
			call ONISHI_OVERLAP(on_ol,cut,U_N,U_N,V_N,V_N,DMAT,i*dPhi)
			
			TrcRhoR 		= 	spart(cut,U_N,U_N,V_N,V_N,DMAT,EYE)

			ROL(i+N_loop+1) 	= 	dPhi*exp(cmplx(0,-1*i*dPhi*j,16))*Pf2P_ol
			ROL2(i+N_loop+1) 	= 	dPhi*exp(cmplx(0,-1*i*dPhi*j*0.5,16))*on_ol
			!ROL3(i+N_loop+1) 	= 	dPhi*exp(cmplx(0,-1*i*dPhi*j*0.5,16))*on_ol
			NROL(i+N_loop+1) 	= 	ROL(i+N_loop+1)*TrcRhoR
			NROL2(i+N_loop+1) 	= 	ROL2(i+N_loop+1)*TrcRhoR
			
		end do
		!!$OMP END PARALLEL DO

		summ 	= 1/(2*PI)*sum(ROL)	
		summ2 	= 1/(2*PI)*sum(ROL2)
		!summ3 	= 1/(2*N_loop+1)*sum(ROL3)
 		NROL 	= NROL/(sum(ROL))
 		NROL2 	= NROL2/(sum(ROL2))
		Nexp 	= sum(NROL)
		Nexp2 	= sum(NROL2)
 		
		write(11,'(I6,5E52.35E6)') j, abs(summ),abs(summ2)!,abs(summ3)
		write(3,'(I6,5E52.35E6)') j, abs(Nexp), abs(Nexp2)
		write(*,*) 'done'
	end do

 	deallocate(nucleus,U_N,V_N,U_Z,V_Z,DMAT,Mmat,EYE)
	close(unit=2)
	close(unit=3)
end program main
