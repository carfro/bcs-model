module BCS
	use MO_module
	use analytic_module
	use pfaffian_module
	implicit none
	integer, parameter :: N_tot=330
	real(kind=qp), parameter :: ten_quad = 10.
contains
! Creates the quasiparticle U,V matrices for neutrons (_N) and protons (_Z) using analytical
! solutions from the analytic_module, also calculates the \prod_1^N/2 v_i^2
	Subroutine qpart_creator(nucleus,N,Z,scalefactor,U_N,V_N,prod_N,U_Z,V_Z,prod_Z)

		type(Nucleon), dimension(:,:),intent(in) 	:: nucleus
		real(dp),intent(in)				:: scaleFactor
		integer,intent(in) 				:: N,Z

		complex(kind=qp),intent(out)			:: prod_N,prod_Z
		complex(dp),intent(out) 				:: U_N(N_tot,N_tot),V_N(N_tot,N_tot),&
										U_Z(N_tot,N_tot),V_Z(N_tot,N_tot)

		real(dp) 	:: Theta_N(N_tot),Theta_Z(N_tot),&
					lam_sN,lam_sZ,tol,EV_N(N_tot,2),EV_Z(N_tot,2)

		integer 	:: i,j,step

		step=1000 	! Nbr of points in \lambda vector
		tol=0.001	! tolerance to find root

		call analytic_solve_sweep(nucleus,N,Z,step,tol,&
			EV_N,lam_sN,EV_Z,lam_sZ,scaleFactor)

		U_N=0;V_N=0;prod_N=1;
		U_Z=0;V_Z=0;prod_Z=1;
		do i=1,N_tot,2
			Theta_N(i)= 0.5*dacos(-1.d0+2*EV_N(i,2))
			Theta_N(i+1)= 0.5*dacos(-1.d0+2*EV_N(i+1,2))
			Theta_Z(i)= 0.5*dacos(-1.d0+2*EV_Z(i,2))
			Theta_Z(i+1)= 0.5*dacos(-1.d0+2*EV_Z(i+1,2))
			
			U_N(i,i)=cmplx(dsin(theta_N(i)),0,8)
			U_N(i+1,i+1)=cmplx(dsin(theta_N(i)),0,8)
			U_Z(i,i)=cmplx(dsin(theta_Z(i)),0,8)
			U_Z(i+1,i+1)=cmplx(dsin(theta_Z(i)),0,8)
			if(i+1<=N_tot) then
				prod_N=prod_N*EV_N(i,2)
				V_N(i,i+1)=cmplx(dcos(theta_N(i)),0,8)
				V_N(i+1,i)=cmplx(-1*dcos(theta_N(i)),0,8)

				prod_Z=prod_Z*EV_Z(i,2)
				V_Z(i,i+1)=cmplx(dcos(theta_Z(i)),0,8)
				V_Z(i+1,i)=cmplx(-1*dcos(theta_Z(i)),0,8)
			end if
		end do
	END Subroutine qpart_creator

	FUNCTION prod_calc(V,N) result(prod)
		complex(dp)	:: V(N,N) 
		integer 		:: N,i
		!complex(kind=qp) 	:: prod
		real(kind=qp) 	:: prod
		
		prod=1
		do i=1,N,2
			if(i+1<=N) then
				!prod=prod*V(i,i+1)*V(i,i+1)
				prod=prod*V(i,i+1)**2
			end if
		end do
	END FUNCTION prod_calc

	FUNCTION WTW(U,V,N) result(WW)
		complex(dp) 	:: WW(2*N,2*N),U(N,N),V(N,N)
		integer 	:: N 

		WW(1:N,1:N) = matmul(transpose(V),U)
		WW(N+1:2*N,1:N) = -matmul(transpose(V),V)
		WW(1:N,N+1:2*N) = matmul(transpose(V),V)
		WW(N+1:2*N,N+1:2*N) = matmul(transpose(U),V)
	END FUNCTION WTW

end module BCS
! Constructs the BCS-model and computes the system for \lambda=[1..-1] and plots the result versus number of particles
program main
	use BCS
	USE F95_PFAPACK
	implicit none
! 	Nucleus, filled using MO_module
	type(Nucleon), dimension(:,:), allocatable :: nucleus
! 	Variables used for analytic solution
	COMPLEX(dp) :: 	Pf2H(2),Pf2P(2),&
			U_N(N_tot,N_tot),V_N(N_tot,N_tot),&
			U_Z(N_tot,N_tot),V_Z(N_tot,N_tot)

	COMPLEX(dp), dimension(:,:), allocatable :: 	U_test,V_test,WW_N

	complex(kind=qp) :: Pf,prod_N,prod_Z
	real(kind=qp) 	:: factor,pf2p_real,pf2h_real,prod_test

	REAL(dp) 	:: scaleFactor(7),t1,t2,t3,t4,t5,st,step,mult

	character(45) :: prod_str,pf_str
	character(50) :: pf2h_str
	character(45) :: pf2p_str
	character(39) :: char1
	character(4) :: char2

	! Nbr of neutrons/protons and loop integer(s),
	integer :: i
	integer :: N,Z,N_mult

	N=24		! Number of NEUTRONS to find
	Z=24 		! Number of PROTONS to find
	
	scaleFactor = (/0.5d0,0.75d0,1.d0,1.75d0,2.5d0,3d0,4d0/)

!-------Analytical solution using BCS-equations
	allocate(nucleus(N_tot,2))

	call nucleus_creator(N,Z,nucleus)

	call qpart_creator(nucleus,N,Z,scalefactor(3),U_N,V_N,prod_N,U_Z,V_Z,prod_Z)
	
!-------Test-loop for the pfaffians below	

	open(unit=1,file='data/perf_plots/performance_err.dat',status='replace')
	open(unit=2,file='data/perf_plots/performance_time.dat',status='replace')
	
	write(1,'(A4,4A48)') "N" , "Prod" ,"Pf$_{\\text{G-C,Robledo,Bertsch}}$",&
	       	"Pf$_{\\text{Wimmer;Parlett-Reid}}$", "Pf$_{\\text{Wimmer;Housholder}}$"

	write(2,'(A4,4A48)') "N" , "Prod" ,"Pf$_{\\text{G-C,Robledo,Bertsch}}$",&
	       	"Pf$_{\\text{Wimmer;Parlett-Reid}}$", "Pf$_{\\text{Wimmer;Housholder}}$"

	! set start (st) and step-size for the testing loop
	st=0.1
	step=0.05
	do i=1,19
		mult=st + (i-1)*step
		N_mult=floor(N_tot*mult) + mod(floor(N_tot*mult),2)

		allocate(U_test(N_mult,N_mult),V_test(N_mult,N_mult),WW_N(2*N_mult,2*N_mult))

		U_test=U_N(1:N_mult,1:N_mult)
		V_test=V_N(1:N_mult,1:N_mult)
		write(*,*) 'N_mult: ', N_mult, ',   i: ', i

		WW_N = WTW(U_test,V_test,N_mult)

		factor=1
		call cpu_time(t1)
		Pf=overlap_pfaffian(factor,N_mult,N_mult,&
			N_mult,U_test,U_test,V_test,V_test)
		call cpu_time(t2)
		call ZSKPF10_F95(WW_N,Pf2P) 
		call cpu_time(t3)
		WW_N = WTW(U_test,V_test,N_mult)
		call cpu_time(t3)
		call ZSKPF10_F95(WW_N,Pf2H,MTHD='H') 
		call cpu_time(t4)
		prod_test=prod_calc(V_test,N_mult)
		call cpu_time(t5)

		write(*,*) 'prod_test'
		write(*,*) prod_test
		write(prod_str,'(ES43.36)') prod_test

		pf2p_real=real(Pf2P(1))*ten_quad**real(Pf2P(2))
		write(*,*) 'Pf2P'
		write(*,*) pf2p_real
		pf2p_real=abs(abs(Pf2p_real)-prod_test)/prod_test
		write(pf2p_str,'(ES43.36)') pf2p_real

		pf2h_real=real(Pf2h(1))*ten_quad**real(Pf2h(2))
		write(*,*) 'Pf2H'
		write(*,*) pf2h_real 
		pf2h_real=abs(abs(Pf2h_real)-prod_test)/prod_test
	 	write(pf2h_str,*) pf2h_real
		
		write(*,*) 'Pf'
		write(*,*) real(Pf)
		Pf=abs(abs(Pf)-prod_test)/prod_test
		write(pf_str,'(ES43.36)') real(Pf)

! Writes product calculation, and pfaffian using: robledo, pfa_parlett-reid, pfa_householder
		write(1,'(I4,4ES48.38)') N_mult , prod_test,real(Pf), pf2p_real, pf2h_real
		write(*,*)

! Writes time of product calculation, robledo, pfa_parlett-reid, pfa_householder
		write(2,'(I4,4E24.16)') N_mult, t5-t4, t2-t1, t3-t2, t4-t3 
		deallocate(U_test,V_test,WW_N)
	end do

	close(unit=1)
	close(unit=2)
	deallocate(nucleus)
end program main
