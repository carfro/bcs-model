module BCS
	use MO_module
	use analytic_module
	use pfaffian_module
	USE F95_PFAPACK
	implicit none
	!integer, parameter :: N_tot=330
	real(8), parameter :: PI=4.D0*DATAN(1.D0)
	!real(kind=qp), parameter :: ten_quad = 10._qp

contains
	
end module BCS
! Constructs the BCS-model and computes the system for \lambda=[1..-1] and plots the result versus number of particles
program main
	use BCS
	implicit none
! 	Nucleus, filled using MO_module
	type(Nucleon), dimension(:,:), allocatable :: nucleus
! 	Matrices used for analytic solution, and particle projection
	COMPLEX(kind=dp) :: 	U_N(N_tot,N_tot),V_N(N_tot,N_tot),&
				U_Z(N_tot,N_tot),V_Z(N_tot,N_tot),&
				DMAT(N_tot,N_tot),&
				WW_N(2*N_tot,2*N_tot), WW_D(2*N_tot,2*N_tot)

! 	Variables used for the particle projection
	COMPLEX(kind=qp) :: 	Pf2P(2),Pf2P_ol,Pf2P_ol2,on_ol,&
				prod_N,prod_Z,&
				summ,summ2,summ3

	REAL(kind=qp) 	:: 	norm,dPhi,N_exp,N_op
	REAL(kind=dp) 	:: 	tol,lam_sN,lam_sZ,scaleFactor(7),deltaR


	COMPLEX(kind=dp), dimension(:,:), allocatable :: 	U_test,V_test
	REAL(kind=dp), dimension(:,:), allocatable    :: 	EV_N,EV_Z


	! Nbr of neutrons/protons and loop integer(s),
	integer :: i,j,N_loop
	integer :: N,Z,cut,step

	N=24		! Number of NEUTRONS to find
	Z=24 		! Number of PROTONS to find
	deltaR=12.d0/sqrt(real(N+Z))
	
	scaleFactor = (/0.5*deltaR,1*deltaR,1.5*deltaR,3*deltaR,6*deltaR,12*deltaR,24*deltaR/)

!-------Analytical solution using BCS-equations

	step=10000 ! Nbr of points in \lambda vector
	tol=0.001 ! tolerance to find root

	!allocate(nucleus(N_tot,2),EV_N(N_tot,2),EV_Z(N_tot,2))

	call nucleus_creator(N,Z,N_tot,nucleus)
	open(unit=3,file='data/occupation.dat',status='replace')

	do i=1,7
		call analytic_solve_sweep(nucleus,N,Z,step,tol,EV_N,lam_sN,EV_Z,lam_sZ,scaleFactor(i))
				
		write(*,*) 'Delta_pairing  ', scaleFactor(i), ' MeV'
		write(*,*) 'lam_sN, sum(EV_N(:,2),1), nucleus(N,1)%E'
		write(*,*) lam_sN, sum(EV_N(:,2),1), nucleus(N,1)%E
		write(*,*) 'lam_sZ, sum(EV_Z(:,2),1), nucleus(Z,1)%E'
		write(*,*) lam_sZ, sum(EV_Z(:,2),1), nucleus(Z,2)%E
		write(*,*) 
		
		write(3,*) '\Delta_{pairing} ', scaleFactor(i), 'MeV'
		do j=1,N_tot
			write(3,'(2(F20.10))') nucleus(j,1)%E,EV_N(j,2)  
			!write(3,'(I3,3(F20.10))') j,EV_N(j,1),nucleus(j,1)%E,EV_N(j,2)  
		end do 
		write(3,'(/)') 
	end do 

	deallocate(nucleus,EV_N,EV_Z)
	close(unit=3)
end program main
