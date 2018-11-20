module analytic_module
	use MO_module
	implicit none
contains
	function analytic_EV(nucleons,delta,lam) result(EV)
		type(Nucleon),dimension(:), intent(in) :: nucleons
		real(8), intent(in) :: delta,lam
		real(8), dimension(size(nucleons,1),2) :: EV
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
		real(8), intent(in) 			:: tol
		real(8), intent(in), optional 		:: factor
		real(8), dimension(:,:), intent(out) :: EV_N,EV_Z
		real(8), intent(out) 		:: lam_sN,lam_sZ
		
		! 	Variables used for solution	
		real(8) 				:: deltaR,Vsum, scal
		real(8), dimension(:),allocatable 	:: lam_vector
		integer 				:: siz,i

		if(present(factor)) then
			scal=factor*12.d0
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
		vector =  (/((st+i*step),i=1,points)/)
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
end module analytic_module

