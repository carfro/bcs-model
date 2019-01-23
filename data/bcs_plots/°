! Module containing all the subroutines and functions for the modified oscillator calculations
include 'MO_module.f03'
! Module containing print routines aswell as extraction of nucleons from .dat files
include 'print_module.f03'
! Module -||- to solve BCS-model by numerical diagonalization
include 'num_matrix_module.f03'

! Module contanining all the subroutines and functions for the BCS-calculations, using a MO-basis
module BCS 
	use modified_oscillator
	use print_module
	implicit none
contains

	function analytic_EV(nucleons,delta,lam) result(EV)
		type(Nucleon),dimension(:), intent(in) :: nucleons
		REAL(16), intent(in) :: delta,lam
		REAL(16), dimension(size(nucleons,1),2) :: EV
		integer :: i

		do i=1,size(nucleons,1)
			EV(i,1)=sqrt((nucleons(i)%E-lam)**2 + delta**2)
			EV(i,2)=0.5*(1-1*(nucleons(i)%E-lam)/EV(i,1))
		end do
	end function analytic_EV
		
	subroutine analytic_solve_sweep(nucleus,N,Z,step,tol,EV_N,lam_sN,&
			EV_Z,lam_sZ,factor)
		type(Nucleon),dimension(:,:), intent(in) :: nucleus
		integer, intent(in) 			:: N,Z,step
		REAL(16), intent(in) 			:: tol
		REAL(16), intent(in), optional 		:: factor
		REAL(16), dimension(:,:), intent(out) :: EV_N,EV_Z
		REAL(16), intent(out) 		:: lam_sN,lam_sZ
		
		! 	Variables used for solution	
		REAL(16) 				:: deltaR,Vsum, scal
		REAL(16), dimension(:),allocatable 	:: lam_vector
		integer 				:: siz,i

		if(present(factor)) then
			scal=factor
		else
			scal=12.d0
		end if

		deltaR=scal/sqrt(real(N+Z))

		allocate(lam_vector(step))
		lam_vector=linspace( lower_energy(nucleus(:,1),N),&
			higher_energy(nucleus(:,1),N),step)

		i=0
		do
			if (i>=step) exit
			siz = size(lam_vector,1)
			lam_sN=lam_vector(siz/2)

			EV_N = analytic_EV(nucleus(:,1),deltaR,lam_sN)
			Vsum=sum(EV_N(:,2),1)

			if (abs(Vsum-N)<tol) then
			       exit
			else if (Vsum>N) then
				lam_vector=lam_vector(1:siz/2)	       
			else
				lam_vector=lam_vector(siz/2:siz)
			end if
			i=i+1
		end do
		
		lam_vector=linspace( lower_energy(nucleus(:,2),Z),&
			higher_energy(nucleus(:,2),Z),step)
		i=0
		do
			if (i>=step) exit
			siz = size(lam_vector,1)
			lam_sZ=lam_vector(siz/2)

			EV_Z = analytic_EV(nucleus(:,1),deltaR,lam_sZ)
			Vsum=sum(EV_Z(:,2),1)

			if (abs(Vsum-N)<tol) then
			       exit
			else if (Vsum>N) then
				lam_vector=lam_vector(1:siz/2)	       
			else
				lam_vector=lam_vector(siz/2:siz)
			end if
			i=i+1
		end do
		deallocate(lam_vector)
	end subroutine analytic_solve_sweep

	function linspace(st,en,points) result(vector)
		real(16), intent(in) 		:: st,en
		integer, intent(in) 		:: points
		real(16), dimension(points) 	:: vector
		real(16) 			:: step
		integer 			:: i

		step = (en-st)/( points )
		!vector(1) = st
		vector =  (/((st+i*step),i=1,points)/)
		!vector(points) = en
	end function linspace

	function lower_energy(nucleons,N) result(E_L)
		type(Nucleon),dimension(:), intent(in) :: nucleons
		integer, intent(in) :: N
		integer :: R
		real(16) :: E_L,E0

		E0=nucleons(N)%E
		R=N
		do 
			R=R-1
			if(nucleons(R)%E<E0) exit
		end do
		E_L=nucleons(R)%E
	end function lower_energy

	function higher_energy(nucleons,N) result(E_H)
		type(Nucleon),dimension(:), intent(in) :: nucleons
		integer, intent(in) :: N
		integer :: R
		real(16) :: E_H,E0

		E0=nucleons(N)%E
		R=N
		do 
			R=R+1
			if(nucleons(R)%E>E0) exit
		end do
		E_H=nucleons(R)%E
	end function higher_energy      	

end module BCS

! Constructs the BCS-model and diagonalizes the system for \lambda=[1..-1] and plots the result versus number of particles 
program main
	use BCS
	implicit none
	integer, parameter :: N_tot=330

	type(Nucleon), dimension(:,:), allocatable :: nucleus

! 	Variables used for numerical solution	
	COMPLEX(16), dimension(:,:), allocatable :: sHN,deltaMN,hamiltonianN
	REAL(16), dimension(:,:), allocatable :: En
	COMPLEX(16), dimension(:,:,:), allocatable :: UVn
	COMPLEX(16), dimension(1000) :: Ln
	COMPLEX(16), dimension(:,:), allocatable :: sHZ,deltaMZ
	COMPLEX(16) :: delta
	
! 	Variables used for analytic solution	
	REAL(16), dimension(:,:), allocatable :: EV_N,EV_Z
	REAL(16) :: lam_sN,lam_sZ,tol
	real(16) :: factor(6) = (\9,10,11,12,13,14,15\)

	! Nbr of neutrons/protons and loop integer(s)
	integer :: i
	integer :: N,Z,siz,step 

!-------Analytical solution using BCS-equations	

	step=10000 	! Nbr of points in \lambda vector
	N=24		! Number of NEUTRONS to find
	Z=24 		! Number of PROTONS to find 
	tol=0.001	! tolerance to find root
	
	allocate(nucleus(N_tot,2),EV_N(N_tot,2),EV_Z(N_tot,2))

	call nucleus_creator(N,Z,nucleus)
	call print_n(nucleus)

	call analytic_solve_sweep(nucleus,N,Z,step,tol,EV_N,lam_sN,EV_Z,lam_sZ)
	
	write(*,*) lam_sN, sum(EV_N(:,2),1), nucleus(N,1)%E
	write(*,*) lam_sZ, sum(EV_Z(:,2),1), nucleus(Z,2)%E

	

	open(unit=1,file='ev_neutrons.dat')
	write(1,'(1F20.10)') EV_N(:,2)  

	deallocate(nucleus,EV_N)	

!-------Numerical solution using diagonalization

	!X=330 ! Neutrons
	!Y=330 ! Protons

	!allocate(nucleus(MAX0(X,Y),2),sHN(X,X),deltaMN(X,X),hamiltonianN(2*X,2*X),&
	!	UVn(2*X,2*X,size(Ln,1)),En(2*X,size(Ln,1)))

	!nucleus=nucleus_extractor(X,Y)

	!delta=12/sqrt(CMPLX(X+Y,0.0))
	!sHN=single_particle_hamiltonian(nucleus(:,1))
	!deltaMN=pairing_matrix(nucleus(:,1),delta)
	!hamiltonianN=hamiltonian_assembler(deltaMN,sHN,CMPLX(1.0,0.0,16))
	!call PRINT_MATRIX('Hamiltonian', 2*X,2*X,hamiltonianN,2*X)

	!call diagonalize(nucleus(:,1),delta,2*X,UVn,En,Ln)
	!call PRINT_RMATRIX('Eigenvalues', 1, 2*X, En(:,1),1)

	!do i=1,1000
	!	call PRINT_RMATRIX('Eigenvalues', 1, 2*X, En(:,i),1)
	!	call PRINT_MATRIX('Eigenvectors (stored columnwise)', 2*X, 2*X, &
	!		UVn(:,:,i),2*X)
	!end do

	!deallocate (nucleus,sHN,deltaMN,hamiltonianN,UVn,En) 
end program main

