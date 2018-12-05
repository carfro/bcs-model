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
	COMPLEX(8) ::  U(N_tot,N_tot) , V(N_tot,N_tot) , SK(2*N_tot,2*N_tot)
	REAL(8) :: V_N_tot(N_tot,7), theta_N(N_tot,7) , Pf
	REAL(8) :: lam_sN,lam_sZ,tol
	REAL(8) :: factor(7)
	Integer,dimension(N_tot,2)::Ipiv


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

	factor = (/0.5d0,0.75d0,1.d0,1.25d0,1.5d0,1.75d0,2d0/)

	do i=1,7
		call analytic_solve_sweep(nucleus,N,Z,step,tol,EV_N,lam_sN,EV_Z,lam_sZ,factor(i))
		V_N_tot(:,i) = EV_N(:,2)
	end do

	U=0;V=0;
	do i=1,N_tot
		U(i,i)=cmplx(dsin(theta_N(i,3)),0,8)
		V(i,i)=cmplx(dcos(theta_N(i,3)),0,8)
	end do

	SK=reshape((/ matmul(transpose(V),U),-matmul(transpose(V),V),&
		matmul(transpose(V),V),matmul(transpose(U),V)/), shape(SK))

	call ZPfaffianF(SK,2*N_tot,2*N_tot,Ipiv,Pf)

	write(*,*) Pf
	deallocate(nucleus,EV_N)
end program main
