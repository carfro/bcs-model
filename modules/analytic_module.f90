module analytic_module
	use MO_module
	implicit none
	integer, parameter :: dp=selected_real_kind(2*kind(1.0))
contains
	function analytic_EV(nucleons,delta,lam,NN) result(EV)
		type(Nucleon),dimension(:), intent(in) :: nucleons
		real(dp), intent(in) :: delta,lam
		integer,intent(in) :: NN
		real(dp)  :: summ
		real(dp), dimension(size(nucleons,1),2) :: EV
		integer :: i

		if(delta==0._dp) then
			summ=0
			i=1
			do while(summ<NN)
				EV(i,1)=sqrt((nucleons(i)%E-lam)**2 + delta**2)
				EV(i,2)=0.5*(1-1*(nucleons(i)%E-lam)/EV(i,1))

				summ=summ+EV(i,2)
				i=i+1
			end do
			EV(i:size(nucleons,1),1)=0
			EV(i:size(nucleons,1),2)=0
		else
			do i=1,size(nucleons,1)
				EV(i,1)=sqrt((nucleons(i)%E-lam)**2 + delta**2)
				EV(i,2)=0.5*(1-1*(nucleons(i)%E-lam)/EV(i,1))
			end do
		end if

	end function analytic_EV
		
	subroutine analytic_solve_sweep(nucleus,N,Z,step,tol,EV_N,lam_sN,&
			EV_Z,lam_sZ,factor)
		type(Nucleon),dimension(:,:), intent(in) :: nucleus
		integer, intent(in) 			:: N,Z,step
		real(dp), intent(in) 			:: tol
		real(dp), intent(in), optional 		:: factor
		real(dp), dimension(:,:), intent(out) :: EV_N,EV_Z
		real(dp), intent(out) 		:: lam_sN,lam_sZ
		
		! 	Variables used for solution	
		real(dp) 				:: deltaR,Vsum, scal
		real(dp), dimension(:),allocatable 	:: lam_vector,lam_vector_old
		integer 				:: siz,i,midP

		if(present(factor)) then
			scal=factor*12.d0
		else
			scal=12.d0
		end if
		
		!open(unit=1,file='data/sweep_test_N.dat',status='replace')

		deltaR=scal/sqrt(real(N+Z))

		allocate(lam_vector(step))
		lam_vector=linspace( lower_energy(nucleus(:,1),N),&
			higher_energy(nucleus(:,1),N),step)

		DO i=1,step
			write(1,*) i

			siz = size(lam_vector,1)
			midP=ceiling(siz/2d0)
			allocate(lam_vector_old(siz))
			lam_vector_old=lam_vector(1:siz)

			!write(1,*) 'size of lam_vector 	' , siz
			!write(1,*) 'midp 	', midP

			lam_sN=lam_vector(midP)

			!write(1,*) 'lam_sN', lam_sN
			!write(1,*) 'E_N 	', nucleus(N,1)%E

			EV_N = analytic_EV(nucleus(:,1),deltaR,lam_sN,N)
			Vsum=sum(EV_N(:,2),1)

			write(1,*) 'Vsum' , Vsum

			if (abs(Vsum-N)<tol) then
				!write(1,*) 'GOLD'
			       	exit 
			else if (Vsum>N) then
				if(midP==1) exit
				deallocate(lam_vector)
				allocate(lam_vector(midP-1))
				lam_vector=lam_vector_old(1:midP-1)	       
				deallocate(lam_vector_old)
				!write(1,*) 'Vsum>N -- first case '
				!write(1,*) 'size of reduced lam_vec  ',size(lam_vector,1)
			else
				if(midP==siz) exit
				deallocate(lam_vector)
				allocate(lam_vector(siz-midP))
				lam_vector=lam_vector_old(midP+1:siz)
				deallocate(lam_vector_old)
				!write(1,*) 'Vsum<N -- second case '
				!write(1,*) 'size of reduced lam_vec  ',size(lam_vector,1)
			end if
			!write(1,*) 
			!write(1,*) 'lam_vector' 
			!write(1,*) lam_vector
		END DO 
		
		if(allocated(lam_vector)) deallocate(lam_vector)
		if(allocated(lam_vector_old)) deallocate(lam_vector_old)

		allocate(lam_vector(step))
		lam_vector=linspace( lower_energy(nucleus(:,2),Z),&
			higher_energy(nucleus(:,2),Z),step)

		do i=1,step
			siz = size(lam_vector,1)
			midP=ceiling(siz/2d0)
			allocate(lam_vector_old(siz))
			lam_vector_old=lam_vector(1:siz)
			lam_sZ=lam_vector(midP)

			EV_Z = analytic_EV(nucleus(:,2),deltaR,lam_sZ,Z)
			Vsum=sum(EV_Z(:,2),1)

			if (abs(Vsum-Z)<tol) then
				!write(1,*) 'GOLD'
			       	exit 
			else if (Vsum>Z) then
				if(midP==1) exit
				deallocate(lam_vector)
				allocate(lam_vector(midP-1))
				lam_vector=lam_vector_old(1:midP-1)	       
				deallocate(lam_vector_old)
			else
				if(midP==siz) exit
				deallocate(lam_vector)
				allocate(lam_vector(siz-midP))
				lam_vector=lam_vector_old(midP+1:siz)
				deallocate(lam_vector_old)
			end if
		end do
		!close(unit=1)
		if(allocated(lam_vector)) deallocate(lam_vector)
		if(allocated(lam_vector_old)) deallocate(lam_vector_old)
	end subroutine analytic_solve_sweep

	function linspace(st,en,points) result(vector)
		real(dp), intent(in) 		:: st,en
		integer, intent(in) 		:: points
		real(dp), dimension(points) 	:: vector
		real(dp) 			:: step
		integer 			:: i

		step = (en-st)/( points )
		vector =  (/((st+i*step),i=1,points)/)
	end function linspace

	function lower_energy(nucleons,N) result(E_L)
		type(Nucleon),dimension(:), intent(in) :: nucleons
		integer, intent(in) :: N
		integer :: R
		real(dp) :: E_L,E0

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
		real(dp) :: E_H,E0

		E0=nucleons(N)%E
		R=N

		do 
			R=R+1
			if(nucleons(R)%E>E0) exit
		end do
		E_H=nucleons(R)%E
	end function higher_energy      	
end module analytic_module

