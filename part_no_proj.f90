module part_no_module
	use mroutines
	use OLroutines
	use MO_module
	use analytic_module
	use pfaffian_module
	use F95_PFAPACK
	implicit none
contains

	complex(qp) function spart(N,UD,UM,VD,VM,DMAT,OP) result(trace)
		integer, intent(in) 				:: N
		complex(kind=dp),dimension(N,N), intent(in) 	:: DMAT, UD, VD, UM, VM, OP
		complex(kind=dp),dimension(N,N) 	 	:: X,Xflip
		real(kind=dp) 					:: y
		integer 					:: i
		
		X = transpose(matmul(transpose(conjg(UD)),matmul(transpose(conjg(DMAT)),UM))+&
			matmul(transpose(conjg(VD)),matmul(transpose(DMAT),VM)))

		!Xflip = transpose(matmul(transpose(conjg(UD)),matmul(transpose(DMAT),UM))+&
		!	matmul(transpose(conjg(VD)),matmul(conjg(DMAT),VM)))

		X = matmul(matmul(DMAT,conjg(VD)),matmul(diag_inv(N,X),transpose(VM)))
		!Xflip = matmul(conjg(VD),matmul(diag_inv(N,Xflip),matmul(transpose(VM),transpose(conjg(DMAT)))))

		!y = NORM2(transfer(X-transpose(conjg(Xflip)),[real(kind(X-transpose(conjg(Xflip))))::]))
		!if(y>ten_quad**-4) write(*,*) y
		!write(*,*) y

		trace = 0	

		do i=1,N
			trace=trace+X(i,i)
		end do

	end function

	function ext_pfaffian(N,UD,UM,VD,VM,DMAT,get_norm) result(pf)
	    ! computes the overlap of two wave functions with pfaffian formula
		integer, intent(in) 				:: N
		complex(kind=dp),dimension(N,N), intent(in) 	:: DMAT, UD, VD, UM, VM
		REAL(kind=qp), intent(in), optional		:: get_norm

		complex(kind=dp) 	:: 	WTW(2*N,2*N)
		complex(kind=qp)	:: 	pf
		real(kind=qp) 		:: 	norm
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
	complex(kind=dp) :: 	Pf2P(2),Pf2P_ol,Pf2P_ol2
	complex(kind=qp) ::	on_ol,ext_ol,prod_N,prod_Z,&
		 		summ,summ2,summ3,TrcRhoR,&
				Nexp,Nexp2,Nexp3

	real(kind=qp) 	:: 	norm
	real(kind=dp) 	:: 	dPhi,N_exp,N_op,tol,lam_sN,lam_sZ,scaleFactor(7),deltaR,phi_arg,phi_arg2

	complex(kind=dp), dimension(:,:), allocatable :: 	U_N,V_N,U_Z,V_Z,DMAT,Mmat,EYE
	complex(kind=dp), dimension(:), allocatable   :: 	ROL,ROL2,ROL3,NROL,NROL2,NROL3
	real(kind=dp), dimension(:,:), allocatable    :: 	EV_N,EV_Z

	! Nbr of neutrons/protons and loop integer(s),
	integer, parameter :: cut = 50
	character(len=1024) :: outfmt
	character(len=3) :: x1
	
	integer :: i,j,N_loop,k,l,w
	integer :: N,Z,step

	N=20 	! Number of NEUTRONS to find
	Z=20 	! Number of PROTONS to find
	deltaR=12.d0/sqrt(real(N+Z))
	
	scaleFactor = [0.5*deltaR,1*deltaR,1.5*deltaR,3*deltaR,6*deltaR,12*deltaR,24*deltaR]

!------------------------- Test-loop for the particle nbr op exp value ----------------------------------
	allocate(nucleus(cut,2),U_N(cut,cut),V_N(cut,cut),U_Z(cut,cut),&
		V_Z(cut,cut),DMAT(cut,cut),Mmat(2*cut,2*cut),EYE(cut,cut))

	call nucleus_creator(N,Z,cut,nucleus)
	call qpart_creator(nucleus,N,Z,cut,scaleFactor(2),U_N,V_N,prod_N,U_Z,V_Z,prod_Z)

! 	overlap testing files
	open(unit=2,file='data/part_no/overlap_proj_test.dat',status='replace')
	open(unit=3,file='data/part_no/Nexp.dat',status='replace')
	open(unit=9,file='data/part_no/overlap_real.dat',status='replace')
	open(unit=10,file='data/part_no/overlap_imag.dat',status='replace')
	open(unit=11,file='data/part_no/sum.dat',status='replace')
! 	U,V matrices
	open(unit=4,file='matlab/U_mat.dat',status='replace')
	open(unit=5,file='matlab/V_mat.dat',status='replace')
	open(unit=6,file='matlab/A_mat.dat',status='replace')

	write(x1,'(I3)') cut
	outfmt='"('//trim(adjustl(x1))//'F41.30)"'

	do l=1,cut
		write(4,"(100F41.30)") (real(U_N(l,k)),k=1,cut)
	end do 	
	do l=1,cut
		write(5,'(100F41.30)') (real(V_N(l,k)),k=1,cut)
	end do 	

	norm=prod_calc(V_N,cut)
	write(*,*) 'norm 	', norm

	N_loop=1*10!*20 	!5*30*30
	w=6
	allocate(ROL(2*N_loop+1),ROL2(2*N_loop+1),ROL3(2*N_loop+1),NROL(2*N_loop+1),NROL2(2*N_loop+1),NROL3(2*N_loop+1))
	
	N_exp=2*sum_check(V_N,cut)
	write(*,*) 'Nbr of part 	', N_exp

	dPhi=2*PI/(2*N_loop+1)	

	Mmat=0
	DMAT=0
	EYE=0

	Pf2P_ol=0
	ext_ol=0
	on_ol=0

	do i=1,cut
		EYE(i,i)=cmplx(1,0,16)
	end do

	!write(*,*) 'Trace(\rho_) ', spart(cut,U_N,U_N,V_N,V_N,EYE,EYE)

	do j=ceiling(N_exp)-w,ceiling(N_exp)+w
		write(*,*) j
		ROL =0
		ROL2=0
		ROL3=0
		summ=0
		summ2=0
		summ3=0

		!!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(DMAT,EYE,U_N,V_N,norm,N_loop,phi_arg)
		do i=0,2*N_loop
			DMAT=EYE*exp(cmplx(0,i*dPhi,16))

			Mmat=M_create(U_N,U_N,V_N,V_N,DMAT,cut)
			call ZSKPF10_MOLP(Pf2P_ol,cut,Mmat,norm)

			ext_ol=ext_pfaffian(cut,U_N,U_N,V_N,V_N,DMAT,norm)
			
			call ONISHI_OVERLAP(on_ol,cut,U_N,U_N,V_N,V_N,DMAT,i*dPhi)
			
			!TrcRhoR 		= 	spart(cut,U_N,U_N,V_N,V_N,EYE,EYE)
			!if(abs(TrcRhoR - N_exp) .gt. ten_quad**-10) write(*,*) 'Trace(\rho_) changed; ', TrcRhoR
			
			TrcRhoR 		= 	spart(cut,U_N,U_N,V_N,V_N,DMAT,EYE)

			ROL(i+1) 	= 	dPhi*exp(cmplx(0,-1*i*dPhi*j,16))*Pf2P_ol
			ROL2(i+1) 	= 	dPhi*exp(cmplx(0,-1*i*dPhi*j*0.5,16))*on_ol
			ROL3(i+1) 	= 	dPhi*exp(cmplx(0,-1*i*dPhi*j,16))*ext_ol
			NROL(i+1) 	= 	ROL(i+1)*TrcRhoR
			NROL2(i+1) 	= 	ROL2(i+1)*TrcRhoR
			NROL3(i+1) 	= 	ROL3(i+1)*TrcRhoR
			
			!write(*,*) abs(conjg(Pf2P_ol)-Pf2P_ol2)
			!if(abs(Pf2P_ol-Pf2P_ol2)>ten_quad**-4) write(*,*) &
			!	char(9),char(9), 'ol - ol*  ',char(9), abs(Pf2P_ol-Pf2P_ol2)

		end do
		!!$OMP END PARALLEL DO

		summ 	= 1/(2*PI)*sum(ROL)	
		summ2 	= 1/(2*PI)*sum(ROL2)
		summ3 	= 1/(2*PI)*sum(ROL3)
 		!NROL 	= NROL/(sum(ROL))
 		!NROL2 	= NROL2/(sum(ROL2))
		Nexp 	= sum(NROL)/(sum(ROL))
		Nexp2 	= sum(NROL2)/(sum(ROL2))
		Nexp3 	= sum(NROL3)/(sum(ROL3))

 		
		write(11,'(I6,5E52.35E6)') j, abs(summ),abs(summ2),abs(summ3)
		write(3,'(I6,5E52.35E6)') j, abs(Nexp), abs(Nexp2),abs(Nexp3)
		write(*,*) 'done'
	end do

 	deallocate(nucleus,U_N,V_N,U_Z,V_Z,DMAT,Mmat,EYE)
	close(unit=2)
	close(unit=3)
end program main
