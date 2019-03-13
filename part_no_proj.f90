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

		complex(kind=dp) 	:: 	Mmat(2*N,2*N)
		complex(kind=qp)	:: 	pf
		real(kind=qp) 		:: 	norm
		integer 		:: 	S_n,IPIV(2*N,2)

		S_n=(-1)**(N*(N-1)/2)

		if( present(get_norm)) then
			norm=get_norm
		else
			norm=prod_calc(VD,N)
		end if

		Mmat=M_create(UD,UM,VD,VM,DMAT,N)

		call Zpfaffian_ext(Mmat,2*N,2*N,IPIV,pf)
		
		pf=S_n*pf/(norm)

	end function ext_pfaffian

	function ZSKPF10_OLP(N,UD,UM,VD,VM,DMAT,get_norm) result(Pf2P_out)
		! computes the overlap of two wave functions \bra{UM,VM} DMAT \ket{UD,VD}   with ZSKPFA from PFAPACK,
		implicit none

		integer, intent(in) 				:: N
		COMPLEX(kind=dp),dimension(N,N), intent(in) 	:: DMAT, UD, VD, UM, VM
		REAL(kind=qp), intent(in), optional		:: get_norm

		COMPLEX(kind=dp) 	:: 	Mmat(2*N,2*N)
		COMPLEX(kind=dp)	:: 	Pf2P(2),Pf2P_out
		REAL(kind=qp) 		:: 	norm
		integer 		:: 	S_n

		S_n=(-1)**(N*(N-1)/2)

		if( present(get_norm)) then
			norm=get_norm
		else
			norm=prod_calc(VD,N)
		end if

		Mmat=M_create(UD,UM,VD,VM,DMAT,N)

		!IF(.NOT. (MAT_NORM2(2*N,WTW+transpose(WTW))<ten_quad**(-4))) write(*,*) MAT_NORM2(N,WTW+transpose(WTW))

		call ZSKPF10_F95(Mmat,Pf2P) 

		Pf2P_out=S_n*cmplx(real(PF2P(1)),aimag(PF2P(1)),kind=qp)*ten_quad**(real(PF2P(2)))/(norm)

	end function ZSKPF10_OLP

	function ZSKPF10_MOLP(N,WDW,norm) result(Pf2P_out)
		! computes the overlap of two wave functions \bra{UM,VM} DMAT \ket{UD,VD}   with ZSKPFA from PFAPACK,
		implicit none

		integer, intent(in) 					:: N
		COMPLEX(kind=dp),dimension(2*N,2*N), intent(in) 	:: WDW
		REAL(kind=qp), intent(in) 				:: norm

		COMPLEX(kind=dp)	:: 	Pf2P(2),Pf2P_out
		COMPLEX(kind=dp) 	:: 	Mmat(2*N,2*N)
		integer 		:: 	i,j,S_n

		Mmat=WDW
		S_n=(-1)**(N*(N-1)/2)

		!IF(.NOT. (MAT_NORM2(2*N,WTW+transpose(WTW))<ten_quad**(-4))) write(*,*) MAT_NORM2(N,WTW+transpose(WTW))

		call ZSKPF10_F95(Mmat,Pf2P) 

		Pf2P_out=S_n*cmplx(real(PF2P(1)),aimag(PF2P(1)),kind=qp)*ten_quad**(real(PF2P(2)))/(norm)

	end function ZSKPF10_MOLP

	function ONISHI_OVERLAP(N,UD,UM,VD,VM,DMAT,phi) result(overlap_out) 
		! computes the overlap of two wave functions \bra{UM,VM} DMAT \ket{UD,VD}   with ZSKPFA from PFAPACK,
		implicit none

		integer, intent(in) :: N
		COMPLEX(kind=dp),dimension(N,N), intent(in) :: DMAT,UD, VD, UM, VM
		real(kind=dp), intent(in) :: phi

		COMPLEX(kind=qp) :: overlap_out

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

	end function ONISHI_OVERLAP

	
end module part_no_module
! Constructs the BCS-model and computes the system for \lambda=[1..-1] and plots the result versus number of particles
program main
	use part_no_module
	implicit none
! 	Nucleus, filled using MO_module
	type(Nucleon), dimension(:,:), allocatable :: nucleus

! 	Variables used for the particle projection
	!complex(kind=dp) :: 	Pf2P(2),Pf2P_ol,Pf2P_ol2
	complex(kind=qp) ::	prod_N,prod_Z,summ,summ2,&
				summ3,TrcRhoR,Nexp,Nexp2,Nexp3

	real(kind=qp) 	:: 	norm
	real(kind=dp) 	:: 	dPhi,N_exp,N_op,tol,lam_sN,lam_sZ,scaleFactor(7),deltaR,halvr

	complex(kind=dp), dimension(:,:), allocatable :: 	U_N,V_N,U_Z,V_Z,DMAT,Mmat,EYE
	complex(kind=qp), dimension(:), allocatable   :: 	ROL,ROL2,ROL3,NROL,NROL2,NROL3,&
								ON_OL,PF_OL,PFE_OL
	real(kind=dp), dimension(:,:), allocatable    :: 	EV_N,EV_Z,OL_VECTOR

	! Nbr of neutrons/protons and loop integer(s),
	integer, parameter :: cut = 100
	character(len=1024) :: outfmt
	character(len=4) :: x1
	
	integer :: i,j
	integer :: N_loop,k,L,w !N
	integer :: N,Z,step,N_P

	N=24 	! Number of NEUTRONS to find
	Z=24 	! Number of PROTONS to find
	deltaR=12.d0/sqrt(real(N+Z))
	
	scaleFactor = [0.25**32*deltaR,1*deltaR,2*deltaR,4*deltaR,8*deltaR,12*deltaR,24*deltaR]

!------------------------- Test-loop for the particle nbr op exp value ----------------------------------
	allocate(nucleus(cut,2),U_N(cut,cut),V_N(cut,cut),U_Z(cut,cut),&
		V_Z(cut,cut),DMAT(cut,cut),Mmat(2*cut,2*cut),EYE(cut,cut))

	call nucleus_creator(N,Z,cut,nucleus)
! 	overlap testing files
	!open(unit=1,file='data/part_no/overlap_proj_test.dat',status='replace')
	open(unit=2,file='data/part_no/Nexp.dat',status='replace')
	open(unit=3,file='data/part_no/rot_ol_sum.dat',status='replace')
	open(unit=6,file='data/part_no/pairing_test.dat',status='replace')

	N_loop=1*10!*20 	!5*30*30
	L=2*N_loop+1
	dPhi=2*PI/(L)	
	w=6

	allocate(ROL(L),ROL2(L),ROL3(L),NROL(L),NROL2(L),NROL3(L),ON_OL(L),PF_OL(L),PFE_OL(L),OL_VECTOR(0:2*w,0:10))

	EYE=0
	do i=1,cut
		EYE(i,i)=cmplx(1,0,16)
	end do

	do k=0,10
		write(*,*) 'k 	', k
		halvr=2.d0**(-1*k)
		call qpart_creator(nucleus,N,Z,cut,halvr*deltaR,U_N,V_N,prod_N,U_Z,V_Z,prod_Z)

	! 	U,V matrices
	!	open(unit=4,file='matlab/U_mat_00625d.dat',status='replace')
	!	open(unit=5,file='matlab/V_mat_00625d.dat',status='replace')
	!
	!! 	U,V matrices printing - for Matlab
	!	write(x1,"(I3)") cut
	!	outfmt='('//trim(adjustl(x1))//'F41.30)'
	!	do l=1,cut
	!		write(4,outfmt) (real(U_N(l,k)),k=1,cut)
	!		write(5,outfmt) (real(V_N(l,k)),k=1,cut)
	!	end do 	

		norm=prod_calc(V_N,cut)
		!write(*,*) 'norm 	', norm

		N_exp=2*sum_check(V_N,cut)
		!write(*,*) 'Nbr of part 	', N_exp
		N_exp = N!ceiling(N_exp)

		DMAT=0
		PF_OL=0
		ON_OL=0
		PFE_OL=0

		do j=0,2*W
			N_P=N_exp-W+j
			!write(*,*) N_P

			ROL   = 0
			ROL2  = 0
			ROL3  = 0
			summ  = 0
			summ2 = 0
			summ3 = 0

			!!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(DMAT,EYE,U_N,V_N,norm,N_loop,phi_arg)
			do i=0,2*N_loop
				DMAT=EYE*exp(cmplx(0,i*dPhi,16))

				if(j==0) PF_OL(i+1)=ZSKPF10_OLP(cut,U_N,U_N,V_N,V_N,DMAT,norm)

				!DMAT=EYE*exp(cmplx(0,i*dPhi*0.5,16))
			!	if(j==0) PFE_OL(i+1)=ext_pfaffian(cut,U_N,U_N,V_N,V_N,DMAT,norm)
			!	
			!	if(j==0) ON_OL(i+1)=ONISHI_OVERLAP(cut,U_N,U_N,V_N,V_N,DMAT,i*dPhi)

				ROL(i+1) 	= 	PF_OL(i+1) *dPhi*exp(cmplx(0,-1*i*dPhi*N_P,16))
			!	ROL2(i+1) 	= 	ON_OL(i+1) *dPhi*exp(cmplx(0,-1*i*dPhi*N_P*0.5,16))
			!	ROL3(i+1) 	= 	PFE_OL(i+1)*dPhi*exp(cmplx(0,-1*i*dPhi*N_P,16))

			!	NROL(i+1) 	= 	ROL(i+1)*spart(cut,U_N,U_N,V_N,V_N,DMAT,EYE)
			!	NROL2(i+1) 	= 	ROL2(i+1)*spart(cut,U_N,U_N,V_N,V_N,DMAT,EYE)
			!	NROL3(i+1) 	= 	ROL3(i+1)*spart(cut,U_N,U_N,V_N,V_N,DMAT,EYE)
			end do
			!!$OMP END PARALLEL DO

			summ 	= 1/(2*PI)*sum(ROL)	
			!summ2 	= 1/(2*PI)*sum(ROL2)
			!summ3 	= 1/(2*PI)*sum(ROL3)
			!Nexp 	= sum(NROL)/(sum(ROL))
			!Nexp2 	= sum(NROL2)/(sum(ROL2))
			!Nexp3 	= sum(NROL3)/(sum(ROL3))
			
			OL_VECTOR(j,k) = abs(summ)
			!write(3,'(I6,5E52.35E6)') j, abs(summ),abs(summ2),abs(summ3)
			!write(2,'(I6,5E52.35E6)')  j, abs(Nexp), abs(Nexp2),abs(Nexp3)
			!write(*,*) 'done'
		end do
		write(*,*) 
	end do

	do j=0,2*W
		N_P=N_exp-W+j
		write(*,*) N_P, (OL_VECTOR(j,k),k=0,10)
		write(6,'(I6,11E52.35E6)') N_P, (OL_VECTOR(j,k),k=0,10)
	end do 

 	deallocate(nucleus,U_N,V_N,U_Z,V_Z,DMAT,Mmat,EYE)
	close(unit=2)
	close(unit=3)
	close(unit=4)
	close(unit=5)
end program main
