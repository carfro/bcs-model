module BCS
	use MO_module
	use analytic_module
	implicit none
	integer, parameter :: qp = selected_real_kind(33, 4931) ! 128-bit real w/ 33 sig. fig, exponent range 4931 
contains
	FUNCTION overlap_pfaffian(lambda,NDIM,dimD,dimM,UD,UM,VD,VM) result(pf)
	    ! computes the overlap of two wave functions with pfaffian formula

	    implicit none


	    integer, intent(in) :: NDIM, dimD, dimM
	    complex*16 :: W((dimM+dimD),(dimM+dimD))
	    complex*16, intent(in) :: UD(NDIM,dimD), VD(NDIM,dimD), UM(NDIM,dimM), VM(NDIM,dimM)
	    integer :: IPIV((dimM+dimD),2)
	    complex(kind=qp) :: pf
	    real(kind=qp),intent(in) :: lambda

	    W(1:dimM,1:dimM) = matmul(transpose(VM),UM)
	    W(1:dimM,(dimM+1):(dimM+dimD)) = matmul(transpose(VM),conjg(VD))
	    W((dimM+1):(dimM+dimD),1:dimM) = -matmul(transpose(conjg(VD)),VM)
	    W((dimM+1):(dimM+dimD),(dimM+1):(dimM+dimD)) = matmul(transpose(conjg(UD)),conjg(VD))

	    W=lambda*W

	    call Zpfaffian_ext(W,(dimM+dimD),(dimM+dimD),IPIV,pf)
	    return
	END FUNCTION overlap_pfaffian
end module BCS

! Constructs the BCS-model and computes the system for \lambda=[1..-1] and plots the result versus number of particles
program main
	use BCS
	implicit none
	integer, parameter :: N_tot=330

! 	Nucleus, filled using MO_module
	type(Nucleon), dimension(:,:), allocatable :: nucleus

! 	Variables used for analytic solution
	REAL(8), dimension(:,:), allocatable :: EV_N,EV_Z
	COMPLEX(8) ::  U(N_tot,N_tot) , V(N_tot,N_tot) , WW(2*N_tot,2*N_tot), &
		A(N_tot,N_tot),B(N_tot,N_tot),C(N_tot,N_tot), &
		WWD(2*N_tot,2*N_tot),WW2(2*24,2*24),Pf2,WW22(2*24,2*24)
		 

	complex(kind=qp) :: Pf22,Pf
	real(kind=qp) 	:: prod,factor,factor2

	Integer 	:: Ipiv(2*N_tot,2),Ipiv2(2*24,2),Ipiv22(2*24,2),Ipiv3(4,2)
				

	REAL(8) :: 	V_N_tot(N_tot,7), Theta_N(N_tot,7), & 
			lam_sN,lam_sZ,tol,scaleFactor(7),&
		 	prod2

	! Nbr of neutrons/protons and loop integer(s),
	integer :: i,j
	integer :: N,Z,siz,step

!-------Analytical solution using BCS-equations

	step=1000 	! Nbr of points in \lambda vector
	N=24		! Number of NEUTRONS to find
	Z=24 		! Number of PROTONS to find
	tol=0.001	! tolerance to find root

	allocate(nucleus(N_tot,2),EV_N(N_tot,2),EV_Z(N_tot,2))

	call nucleus_creator(N,Z,nucleus)

	scaleFactor = (/0.5d0,0.75d0,1.d0,1.75d0,2.5d0,3d0,4d0/)

	do i=1,7
		call analytic_solve_sweep(nucleus,N,Z,step,tol,&
			EV_N,lam_sN,EV_Z,lam_sZ,scaleFactor(i))
		V_N_tot(:,i) = EV_N(:,2) 	! Contains all v_i^2 
	end do

	U=0;V=0;prod=1;prod2=1;
	do i=1,N_tot,2
		Theta_N(i,:)=(/  (0.5*dacos(-1.d0+2*V_N_tot(i,j)), j=1,7) /)
		Theta_N(i+1,:)=(/  (0.5*dacos(-1.d0+2*V_N_tot(i+1,j)), j=1,7) /)

		if(i+1<=24) then
			prod2=prod2*V_N_tot(i,3)
		end if

		U(i,i)=cmplx(dsin(theta_N(i,3)),0,8)
		U(i+1,i+1)=cmplx(dsin(theta_N(i,3)),0,8)
		if(i+1<=N_tot) then
			prod=prod*V_N_tot(i,3)
			V(i,i+1)=cmplx(dcos(theta_N(i,3)),0,8)
			V(i+1,i)=cmplx(-1*dcos(theta_N(i,3)),0,8)
		end if
		!write(*,*) i, real(U(i,i)), real(V(i,i+1))
	end do
	!write(*,*) real(V(N_tot,N_tot-1)), real(V(N_tot-1,N_tot))
	!write(*,*) Theta_N(330,3), dcos(Theta_N(330,3)) 
	write(*,*) 'Prod (=Pf)', prod
	write(*,*) 'Prod2 (=Pf2)', prod2

	A=matmul(transpose(V),U)
	B=matmul(transpose(U),V)
	C=matmul(transpose(V),V)

	WW(1:N_tot,1:N_tot) = A
	WW(N_tot+1:2*N_tot,1:N_tot) = -C
	WW(1:N_tot,N_tot+1:2*N_tot) = C
	WW(N_tot+1:2*N_tot,N_tot+1:2*N_tot) = B

		
	WW2(1:24,1:24) = A(1:24,1:24)
	WW2(24+1:2*24,1:24) = -C(1:24,1:24)
	WW2(1:24,24+1:2*24) = C(1:24,1:24)
	WW2(24+1:2*24,24+1:2*24) = B(1:24,1:24)
	
	WW22=WW2

	open(unit=1,file='data/WW.dat',status='replace')
	open(unit=2,file='data/WWD.dat',status='replace')
	open(unit=3,file='data/A.dat',status='replace')
	open(unit=4,file='data/B.dat',status='replace')
	open(unit=5,file='data/U.dat',status='replace')
	open(unit=6,file='data/V.dat',status='replace')

	WWD=WW+transpose(WW)

	do i=1,2*N_tot
		do j=1,2*N_tot
			if(.not.(abs(WW(i,j))<1.d-5)) then 
				write(1,"(2I5,F18.12)") j,2*N_tot+1-1*i,real(WW(i,j),8) 
			end if

			if(i<=2*24 .and. j<=2*24) then
				if(.not.(abs(WW2(i,j))<1.d-3)) then 
					write(2,"(2I5,F18.12)") j,2*24+1-1*i,real(WW2(i,j),8)
				end if
			end if
			if (i<=N_tot .and. j<=N_tot) then
				if(.not.(abs(A(i,j))<1.d-3)) then 
					!write(3,"(2I5)") j,N_tot+1-1*i
					write(3,"(2I5,F17.12)") j,N_tot+1-1*i, real(A(i,j),8) 
				end if
				if(.not.(abs(B(i,j))<1.d-3)) then 
					write(4,"(2I5)") j,N_tot+1-1*i 
				end if
				if(.not.(abs(U(i,j))<1.d-3)) then 
					write(5,"(2I5,2I5,F18.12)") j,N_tot+1-1*i,i,j,real(U(i,j),8) 
				end if
				if(.not.(abs(V(i,j))<1.d-3)) then 
					write(6,"(2I5,2I5,F18.12)") j,N_tot+1-1*i,i,j,real(V(i,j),8) 
				end if 
				!if(i==j) then
				!	write(*,*) i, j, real(C(i,j),8)
				!end if
			end if 
		end do
	end do
	factor=1
	Pf=overlap_pfaffian(factor,N_tot,N_tot,N_tot,U,U,V,V)
	!call ZPfaffian_EXT(WW,2*N_tot,2*N_tot,Ipiv,Pf)
	call ZPfaffianF(WW2,2*24,2*24,Ipiv2,Pf2)

	call ZPfaffian_EXT(WW22,2*24,2*24,Ipiv22,Pf22)
	
	write(*,*)
	factor2=factor**(N_tot/2)
	write(*,*) 'Pf*factor: ', real(Pf)
	write(*,*) 'factor2: ', real(factor2)
	write(*,*) 'Pf/factor2: ', real(Pf/factor2)
	write(*,*)
	write(*,*) 'Pf2: ', real(Pf2)
	write(*,*) 'Pf22: ', real(Pf22)

	close(unit=1)
	close(unit=2)
	close(unit=3)
	close(unit=4)
	close(unit=5)
	close(unit=6)
	deallocate(nucleus,EV_N)
end program main
