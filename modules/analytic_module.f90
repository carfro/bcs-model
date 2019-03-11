module analytic_module
	use MO_module
	implicit none
	integer, parameter :: N_tot=330
contains
! Creates the quasiparticle U,V matrices for neutrons (_N) and protons (_Z) using analytical
! solutions from the analytic_module, also calculates the \prod_1^N/2 v_i^2
	Subroutine qpart_creator(nucleus,N,Z,lvls,del_pairing,U_N,V_N,prod_N,U_Z,V_Z,prod_Z)

		type(Nucleon), dimension(:,:),intent(in) 	:: nucleus
		real(dp),intent(in)				:: del_pairing
		integer,intent(in) 				:: N,Z,lvls

		COMPLEX(kind=qp),intent(out)			:: prod_N,prod_Z
		COMPLEX(dp),intent(out) 			:: U_N(lvls,lvls),V_N(lvls,lvls),&
								   U_Z(lvls,lvls),V_Z(lvls,lvls)

		REAL(dp) 	:: 	Theta_N(lvls),Theta_Z(lvls),&
				 	EV_N(lvls,2),EV_Z(lvls,2),&
					lam_sN,lam_sZ,tol

		integer 	:: i,j,step

		step=100000 	! Nbr of points in \lambda vector
		tol=0.0000001	! tolerance to find root

		call analytic_solve_sweep(nucleus,N,Z,step,tol,&
					  EV_N,lam_sN,EV_Z,lam_sZ,del_pairing)

		U_N=0 	; V_N=0 	; prod_N=1;
		U_Z=0 	; V_Z=0 	; prod_Z=1;
		do i=1,lvls,2
			Theta_N(i) 	= 0.5*dacos(-1.d0+2*EV_N(i,2))
			Theta_N(i+1) 	= 0.5*dacos(-1.d0+2*EV_N(i+1,2))
			Theta_Z(i) 	= 0.5*dacos(-1.d0+2*EV_Z(i,2))
			Theta_Z(i+1) 	= 0.5*dacos(-1.d0+2*EV_Z(i+1,2))
			
			U_N(i,i) 	= cmplx(dsin(theta_N(i)),0,8)
			U_N(i+1,i+1) 	= cmplx(dsin(theta_N(i)),0,8)
			U_Z(i,i) 	= cmplx(dsin(theta_Z(i)),0,8)
			U_Z(i+1,i+1) 	= cmplx(dsin(theta_Z(i)),0,8)

			if(i+1<=lvls) then
				prod_N 		= prod_N*EV_N(i,2)
				V_N(i,i+1) 	= cmplx(dcos(theta_N(i)),0,8)
				V_N(i+1,i) 	= cmplx(-1*dcos(theta_N(i)),0,8)

				prod_Z 		= prod_Z*EV_Z(i,2)
				V_Z(i,i+1) 	= cmplx(dcos(theta_Z(i)),0,8)
				V_Z(i+1,i) 	= cmplx(-1*dcos(theta_Z(i)),0,8)
			end if
		end do
	END Subroutine qpart_creator

		
! Solves the BCS-model for both protons EV_Z and neutrons EV_N using lambda sweep until the energy gives the correct value of
! <N>² = sum(EV_x(:,2))
	subroutine analytic_solve_sweep(nucleus,N,Z,step,tol,EV_N,lam_sN,&
			EV_Z,lam_sZ,del_pairing)
		type(Nucleon),dimension(:,:), intent(in) :: nucleus
		integer, intent(in) 			 :: N,Z,step
		real(dp), intent(in) 			 :: tol
		real(dp), intent(in), optional 		 :: del_pairing
		real(dp), dimension(:,:), intent(out) :: EV_N,EV_Z
		real(dp), intent(out) 			 :: lam_sN,lam_sZ
		
		! 	Variables used for solution	
		real(dp) 				 :: deltaR,Vsum, scal
		real(dp), dimension(:),allocatable 	 :: lam_vector,lam_vector_old
		integer 				 :: siz,i,midP

		if(present(del_pairing)) then
			deltaR=del_pairing
		else
			write(*,*) 'Empirical pairing activated'
			deltaR=12.d0/sqrt(real(N+Z))
		end if
		
		open(unit=1,file='data/sweep_test_N.dat',status='replace')

		allocate(lam_vector(step),lam_vector_old(step))
		lam_vector=linspace( lower_energy(nucleus(:,1),N),&
			higher_energy(nucleus(:,1),N),step)

		DO i=1,step
			write(1,*) i

			siz = size(lam_vector,1)
			midP=ceiling(siz/2d0)

			lam_sN=lam_vector(midP)

			EV_N = analytic_EV(nucleus(:,1),deltaR,lam_sN,N)
			Vsum=sum(EV_N(:,2),1)

			if (abs(Vsum-N)<tol) then
			       	exit 
			else if (Vsum>N) then
				lam_vector=linspace(lam_vector(1),&
				lam_vector(midP),step)
				if(lam_vector(1)==lam_vector(siz)) then
					lam_vector=linspace(lam_vector(1)-1d0,&
						lam_vector(midP),step)
				end if
			else
				lam_vector=linspace(lam_vector(midP),&
				lam_vector(size(lam_vector,1)),step)
				if(lam_vector(1)==lam_vector(siz)) then
					lam_vector=linspace(lam_vector(1),&
						lam_vector(midP)+1d0,step)
				end if
			end if
		END DO 
		
		if(allocated(lam_vector)) deallocate(lam_vector)
		if(allocated(lam_vector_old)) deallocate(lam_vector_old)

		allocate(lam_vector(step))
		lam_vector=linspace( lower_energy(nucleus(:,2),z),&
			higher_energy(nucleus(:,2),z),step)

		do i=1,step
			siz = size(lam_vector,1)
			midp=ceiling(siz/2d0)

			lam_sz=lam_vector(midp)

			ev_z = analytic_ev(nucleus(:,2),deltar,lam_sz,z)
			vsum=sum(ev_z(:,2),1)

			if (abs(Vsum-N)<tol) then
			       	exit 
			else if (Vsum>N) then
				lam_vector=linspace(lam_vector(1),&
				lam_vector(midP),step)
				if(lam_vector(1)==lam_vector(siz)) then
					lam_vector=linspace(lam_vector(1)-1d0,&
						lam_vector(midP),step)
				end if
			else
				lam_vector=linspace(lam_vector(midP),&
				lam_vector(size(lam_vector,1)),step)
				if(lam_vector(1)==lam_vector(siz)) then
					lam_vector=linspace(lam_vector(1),&
						lam_vector(midP)+1d0,step)
				end if
			end if
		end do
		!close(unit=1)
		if(allocated(lam_vector)) deallocate(lam_vector)
		if(allocated(lam_vector_old)) deallocate(lam_vector_old)
	end subroutine analytic_solve_sweep

!Computes the quasi-particle energy E(i,1) and occupation prob E(i,2)=v_i²
	function analytic_EV(nucleons,delta,lam,NN) result(EV)
		type(Nucleon),dimension(:), intent(in) :: nucleons
		real(dp), intent(in) :: delta,lam
		integer,intent(in) :: NN
		real(dp)  :: summ
		real(dp), dimension(size(nucleons,1),2) :: EV
		integer :: i

		if(delta==0._dp) then
			!write(*,*) 'Delta 0'
			summ=0
			i=1
			do while(summ<NN)
				write(*,*) 'analytic_EV FAIL'
				EV(i,1)=sqrt((nucleons(i)%E-lam)**2 + delta**2)
				EV(i,2)=0.5*(1-1*(nucleons(i)%E-lam)/EV(i,1))

				summ=summ+EV(i,2)
				i=i+1
			end do
			EV(i:size(nucleons,1),1)=0
			EV(i:size(nucleons,1),2)=0
		else
			!write(*,*) 'delta'
			!write(*,*) delta
			do i=1,size(nucleons,1)
				EV(i,1)=sqrt((nucleons(i)%E-lam)**2 + delta**2)
				EV(i,2)=0.5*(1-1*(nucleons(i)%E-lam)/EV(i,1))
				!write(*,*) EV(i,1)
			end do
		end if

	end function analytic_EV

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

		do while(R>2)
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

		do while(R+1<=size(nucleons,1)) 
			R=R+1
			if(nucleons(R)%E>E0) exit
		end do
		E_H=nucleons(R)%E
	end function higher_energy      	
end module analytic_module

