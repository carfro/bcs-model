module BCS
	use MO_module
	use analytic_module
	implicit none
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
	COMPLEX(8) ::  U(N_tot,N_tot) , V(N_tot,N_tot) , SK(2*N_tot,2*N_tot), Pf, &
		A(N_tot,N_tot),B(N_tot,N_tot),C(N_tot,N_tot), W(2*N_tot,2*N_tot),&
		SK3(4,4),Pf3,Ipiv(2*N_tot,2),Ipiv3(4,2),SKD(2*N_tot,2*N_tot),&
		 SK2(2*24,2*24),Pf2,Ipiv2(2*24,2)

	REAL(8) :: V_N_tot(N_tot,7), Theta_N(N_tot,7) 
	REAL(8) :: lam_sN,lam_sZ,tol,factor(7)
	REAL(8) :: prod,prod2

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

	factor = (/0.5d0,0.75d0,1.d0,1.75d0,2.5d0,3d0,4d0/)

	do i=1,7
		call analytic_solve_sweep(nucleus,N,Z,step,tol,&
			EV_N,lam_sN,EV_Z,lam_sZ,factor(i))
		V_N_tot(:,i) = EV_N(:,2)
	end do

	U=0;V=0;prod2=1
	do i=1,N_tot,2
		Theta_N(i,:)=(/  (0.5*dacos(-1.d0+2*V_N_tot(i,j)), j=1,7) /)
		Theta_N(i+1,:)=(/  (0.5*dacos(-1.d0+2*V_N_tot(i+1,j)), j=1,7) /)

		if(i+1<=24) then
			prod2=prod2*V_N_tot(i,3)
			write(*,*) i, prod,prod2
		end if

		U(i,i)=cmplx(dsin(theta_N(i,3)),0,8)
		U(i+1,i+1)=cmplx(dsin(theta_N(i,3)),0,8)
		if(i+1<=N_tot) then
			!write(*,*) i
			V(i,i+1)=cmplx(dcos(theta_N(i,3)),0,8)
			V(i+1,i)=cmplx(-1*dcos(theta_N(i,3)),0,8)
		end if
		!write(*,*) i, real(U(i,i)), real(V(i,i+1))
	end do
	!write(*,*) real(V(N_tot,N_tot-1)), real(V(N_tot-1,N_tot))
	!write(*,*) Theta_N(330,3), dcos(Theta_N(330,3)) 

	A=matmul(transpose(V),U)
	B=matmul(transpose(U),V)
	C=matmul(transpose(V),V)

	SK(1:N_tot,1:N_tot) = A
	SK(N_tot+1:2*N_tot,1:N_tot) = -C
	SK(1:N_tot,N_tot+1:2*N_tot) = C
	SK(N_tot+1:2*N_tot,N_tot+1:2*N_tot) = B

	W=reshape((/  U ,V ,&
			V , U /), shape(W))

	SK3=reshape((/ 	(0,0),(-1,0),(0,0),(0,0),&
			(1,0),(0,0),(-2,0),(0,0),&
			(0,0),(2,0),(0,0),(-3,0),&
			(0,0),(0,0),(3,0),(0,0)/), shape(SK3))	
		
	SK2(1:24,1:24) = A(1:24,1:24)
	SK2(24+1:2*24,1:24) = -C(1:24,1:24)
	SK2(1:24,24+1:2*24) = C(1:24,1:24)
	SK2(24+1:2*24,24+1:2*24) = B(1:24,1:24)

	open(unit=1,file='data/SK.dat',status='replace')
	open(unit=2,file='data/SKD.dat',status='replace')
	open(unit=3,file='data/A.dat',status='replace')
	open(unit=4,file='data/B.dat',status='replace')
	open(unit=5,file='data/U.dat',status='replace')
	open(unit=6,file='data/V.dat',status='replace')

	SKD=SK+transpose(SK)

	do i=1,2*N_tot
		do j=1,2*N_tot
			if(.not.(abs(SK(i,j))<1.d-5)) then 
				write(1,"(2I5,F18.12)") j,2*N_tot+1-1*i,real(SK(i,j),8) 
			end if

			if(i<=2*24 .and. j<=2*24) then
				if(.not.(abs(SK2(i,j))<1.d-3)) then 
					write(2,"(2I5,F18.12)") j,2*24+1-1*i,real(SK2(i,j),8)
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
	call ZPfaffianF(SK,2*N_tot,2*N_tot,Ipiv,Pf)
	call ZPfaffianF(SK2,2*24,2*24,Ipiv2,Pf2)
	call ZPfaffianF(SK3,4,4,Ipiv3,Pf3)

	write(*,*) Pf
	write(*,*) Pf2
	write(*,*) Pf3
	close(unit=1)
	close(unit=2)
	close(unit=3)
	close(unit=4)
	close(unit=5)
	close(unit=6)
	deallocate(nucleus,EV_N)
end program main
