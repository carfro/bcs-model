module part_no_module
	use MO_module
	use analytic_module
	use pfaffian_module
	USE F95_PFAPACK
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
		COMPLEX(dp) 	:: WW(N,N)
		REAL(dp) 	:: NORM
		integer 	:: N,i,j

		NORM=0
		WW=matmul(conjg(WW),WW)
		do i=1,N
			NORM=NORM+WW(i,i)	
		end do

		NORM=sqrt(NORM)
	END FUNCTION MAT_NORM2

	FUNCTION WD(UD,UM,VD,VM,DMAT,N) 
		COMPLEX(kind=dp),dimension(N,N), intent(in) 	:: DMAT, UD, VD, UM, VM
		COMPLEX(dp) 	:: WD(2*N,2*N)
		integer 	:: N 

		WD(1:N,1:N) = matmul(transpose(VM),UM)
		WD(1:N,(N+1):2*N) = matmul(transpose(VM),matmul(DMAT,conjg(VD)))
		WD((N+1):(2*N),1:N) = -1*matmul(transpose(conjg(VD)),matmul(transpose(DMAT),VM))
		WD((N+1):(2*N),(N+1):(2*N)) = matmul(transpose(conjg(UD)),conjg(VD))
	END FUNCTION WD

	FUNCTION ext_pfaffian(N,UD,UM,VD,VM,DMAT,get_norm) result(pf)
	    ! computes the overlap of two wave functions with pfaffian formula
		implicit none

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
	END FUNCTION ext_pfaffian

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


	complex(kind=qp) function det_eig(N, mat_in)
		! Computes the determinant of a matrix mat_NxN using eig.val. decomposition (ZGEEV) from lapack
		implicit none

		integer, intent(in) :: N 
		complex(kind=qp), intent(inout), dimension(:,:) :: mat_in
		integer 	   ::          lda, ldvl, ldvr
            	integer, parameter ::         lwmax=3000

		integer(kind=dp) :: i,j, info,lwork
		integer, allocatable :: ipiv(:)
		real(kind=qp) :: sgn

		real(kind=dp), allocatable :: rwork(:)
		complex(kind=dp), allocatable :: vl(:,:), vr(:,:), w(:), work(:),mat_used(:,:)
		
      		lda = N; ldvl = N; ldvr = N
		allocate(ipiv(N),mat_used(N,N),rwork(2*N),vl(ldvl,N),vr(ldvr,N),w(N),work(lwmax))
		ipiv = 0
		
 		mat_used=mat_in
		    
		info=0
		lwork=-1
		call zgeev('vectors','vectors',N,mat_used,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
		lwork=min(lwmax,int(work(1)))
		!write(*,*) 	
		!write(*,*) 'lwork' , lwork

		call zgeev('vectors','vectors',N,mat_used,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)

		if(info.gt.0) then
			write(*,*) 'the algorithm failed to compute eigenvalues.'
			write(*,*) 'info 	',info
			stop
		end if
		
		!call print_matrix('eigenvalues',1,N,w,1)
		!write(*,*) 'info 	',info
		!write(*,*) 

		det_eig = 1._qp

		do i = 1, N
			det_eig = det_eig*w(i)
			!write(*,*) 'det_eig 	', det_eig
		end do

	end function det_eig

	complex(kind=qp) function det(N, mat)
		! Computes the determinant of a matrix mat_NxN using LU-decomposition from lapack
		implicit none

		INTEGER, intent(in) :: N 
		complex(kind=qp), intent(in), dimension(:,:) :: mat
		complex(kind=qp),allocatable, dimension(:,:) :: mat_used
		integer(kind=dp) :: i, info
		integer, allocatable :: ipiv(:)
		real(kind=qp) :: sgn

		allocate(ipiv(N),mat_used(N,N))
		mat_used=mat
		ipiv = 0
		    
		call zgetrf(N, N, mat_used, N, ipiv, info)
		    
		det = 1._qp

		do i = 1, N
			det = det*mat_used(i, i)
			if(abs(mat_used(i,i))<10._qp**(-4) ) write(*,*) real(mat_used(i,i))
		end do

		sgn = 1._qp

		do i = 1, N
			if(ipiv(i) /= i) sgn = -sgn
		end do

		det = sgn*det   

	end function det


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
	
end module part_no_module
! Constructs the BCS-model and computes the system for \lambda=[1..-1] and plots the result versus number of particles
program main
	use part_no_module
	implicit none
! 	Nucleus, filled using MO_module
	type(Nucleon), dimension(:,:), allocatable :: nucleus
! 	Matrices used for analytic solution, and particle projection
	!COMPLEX(kind=dp) :: 	U_N(N_tot,N_tot),V_N(N_tot,N_tot),&
	!			U_Z(N_tot,N_tot),V_Z(N_tot,N_tot),&
	!			WW_N(2*N_tot,2*N_tot)!, WW_D(2*N_tot,2*N_tot)!,&
				!DMAT(N_tot,N_tot)

! 	Variables used for the particle projection
	COMPLEX(kind=dp) :: 	Pf2P(2),Pf2P_ol,Pf2P_ol2
	COMPLEX(kind=qp) ::	on_ol,ext_ol,prod_N,prod_Z,&
				summ,summ2,summ3

	REAL(kind=qp) 	:: 	norm
	REAL(kind=dp) 	:: 	dPhi,N_exp,N_op,tol,lam_sN,lam_sZ,scaleFactor(7),deltaR,phi_arg,phi_arg2

	COMPLEX(kind=dp), dimension(:,:), allocatable :: 	U_N,V_N,U_Z,V_Z,DMAT,WDW,EYE
	REAL(kind=dp), dimension(:,:), allocatable    :: 	EV_N,EV_Z


	! Nbr of neutrons/protons and loop integer(s),
	integer :: i,j,N_loop,k,l
	integer :: N,Z,step

	integer, parameter :: cut = 100
	!CHARACTER(LEN=3), DIMENSION(2*cut,2*cut) :: imag_unit
	!!imag_unit = '+i*'
	!imag_unit = '+'

	N=24	! Number of NEUTRONS to find
	Z=24 	! Number of PROTONS to find
	deltaR=12.d0/sqrt(real(N+Z))
	
	scaleFactor = (/0.5*deltaR,1*deltaR,1.5*deltaR,3*deltaR,6*deltaR,12*deltaR,24*deltaR/)

!------------------------- Test-loop for the particle nbr op exp value ----------------------------------

	allocate(nucleus(cut,2),U_N(cut,cut),V_N(cut,cut),U_Z(cut,cut),&
		V_Z(cut,cut),DMAT(cut,cut),WDW(2*cut,2*cut),EYE(cut,cut))

	call nucleus_creator(N,Z,cut,nucleus)
	call qpart_creator(nucleus,N,Z,cut,5d0,U_N,V_N,prod_N,U_Z,V_Z,prod_Z)
! 	overlap testing files
	open(unit=1,file='data/part_no/part_no_test.dat',status='replace')
	open(unit=2,file='data/part_no/overlap_proj_test.dat',status='replace')
	open(unit=9,file='data/part_no/overlap_real.dat',status='replace')
	open(unit=10,file='data/part_no/overlap_imag.dat',status='replace')
! 	det test
	open(unit=3,file='data/part_no/det_test.dat',status='replace')
! 	U,V matrices
	open(unit=4,file='matlab/U_mat.dat',status='replace')
	open(unit=5,file='matlab/V_mat.dat',status='replace')
! 	W matrices

	do l=1,cut
		write(4,'(100F41.30)') (real(U_N(l,k)),k=1,cut)
	end do 	
	do l=1,cut
		write(5,'(100F41.30)') (real(V_N(l,k)),k=1,cut)
	end do 	

	norm=prod_calc(V_N,cut)
	write(*,*) 'norm 	', norm

	!N_loop=5*30*30
	N_loop=10
	
	N_exp=2*sum_check(V_N,cut)
	write(*,*) 'sum_check 	', N_exp

	dPhi=4*4*PI/((cut+N_exp)*(N_loop))	
	phi_arg=dPhi!*N_exp

	WDW=0
	DMAT=0
	EYE=0

	Pf2P_ol=0
	ext_ol=0
	on_ol=0

	forall(i=1:cut) EYE(i,i)=cmplx(1,0,16)

	do j=N_exp,N_exp!N_exp-2,N_exp+2
		write(*,*) j
		!summ=0
		!summ2=0
		!!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(DMAT,EYE,U_N,V_N,norm,N_loop,phi_arg)
		!do i=0.25*N_loop,0.25*N_loop
		do i=0,N_loop
			write(*,*) i!, i*dPhi

			!DMAT=EYE*exp(cmplx(0,i*phi_arg,16))
			if(i>0) DMAT=DMAT*exp(cmplx(0,phi_arg,16))

			WDW=WD(U_N,U_N,V_N,V_N,DMAT,cut)

			call ZSKPF10_MOLP(Pf2P_ol,cut,WDW,norm)
			!call ZSKPF10_OLP(Pf2P_ol,N,U_N,U_N,V_N,V_N,DMAT,norm)

			ext_ol=ext_pfaffian(cut,U_N,U_N,V_N,V_N,DMAT,norm)
			
			call ONISHI_OVERLAP(on_ol,cut,U_N,U_N,V_N,V_N,DMAT,i*phi_arg)

			write(2,'(I6,F12.7,3E52.35E6)') i, i*dPhi, abs(Pf2P_ol), abs(on_ol), abs(ext_ol)
			write(9,'(I6,F12.7,3E52.35E6)') i, i*dPhi, real(Pf2P_ol), real(on_ol), real(ext_ol)
			write(10,'(I6,F12.7,3E52.35E6)') i, i*dPhi, aimag(Pf2P_ol), aimag(on_ol), aimag(ext_ol)

		!	summ= summ + dPhi*exp(cmplx(0,-1*i*dPhi*j,16))*Pf2P_ol
		!	summ2= summ2 + dPhi*exp(cmplx(0,-1*i*dPhi*j,16))*exp(cmplx(0,1*dPhi*i*N_exp,16))
		end do
		!summ=1/(2*PI)*summ
		!write(1,'(I2,3ES48.38)') j,real(summ),real(summ2),real(summ3)
	end do
	!!$OMP END PARALLEL DO

 	deallocate(nucleus,U_N,V_N,U_Z,V_Z,DMAT,WDW,EYE)
	close(unit=2)
	close(unit=3)
end program main
