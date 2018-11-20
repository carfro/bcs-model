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
	REAL(8) :: V_N_tot(N_tot,7), theta_N(N_tot,7), U(N_tot,N_tot) , V(N_tot,N_tot)
	REAL(8) :: lam_sN,lam_sZ,tol
	REAL(8) :: factor(7)

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

	!call print_n(nucleus)
	!call analytic_solve_sweep(nucleus,N,Z,step,tol,EV_N,lam_sN,EV_Z,lam_sZ)
	!write(*,*) lam_sN, sum(EV_N(:,2),1), nucleus(N,1)%E
	!write(*,*) lam_sZ, sum(EV_Z(:,2),1), nucleus(Z,2)%E

	open(unit=1,file='data/ev_neutrons.dat')
	open(unit=2,file='data/eu_neutrons.dat')
	factor = (/0.5d0,0.75d0,1.d0,1.25d0,1.5d0,1.75d0,2d0/)

	do i=1,7
		call analytic_solve_sweep(nucleus,N,Z,step,tol,EV_N,lam_sN,EV_Z,lam_sZ,factor(i))
		V_N_tot(:,i) = EV_N(:,2)
	end do

	do i=1,N_tot
		Theta_N(i,:)=(/  (0.5*dacos(-1.d0+2*V_N_tot(i,j)), j=1,7) /)
		write(1,'(15F20.10)') nucleus(i,1)%E , ( V_N_tot(i,j) , j=1,7)
		write(2,'(15F20.10)') nucleus(i,1)%E , ( dsin(theta_N(i,j))**2 , j=1,7)
	end do

	U=0;V=0;
	do i=1,N_tot
		U(i,i)=dsin(theta_N(i,3))
		V(i,i)=dcos(theta_N(i,3))
	end do

	write(*,*) shape(U), shape(V)
	deallocate(nucleus,EV_N)
end program main
