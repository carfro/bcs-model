module BCS
	use MO_module
	use analytic_module
	use pfaffian_module
	USE F95_PFAPACK
	implicit none
	integer, parameter :: N_tot=330
	real(8), parameter :: PI=4.D0*DATAN(1.D0)
	real(kind=qp), parameter :: ten_quad = 10.

contains
! Creates the quasiparticle U,V matrices for neutrons (_N) and protons (_Z) using analytical
! solutions from the analytic_module, also calculates the \prod_1^N/2 v_i^2
	Subroutine qpart_creator(nucleus,N,Z,scalefactor,U_N,V_N,prod_N,U_Z,V_Z,prod_Z)

		type(Nucleon), dimension(:,:),intent(in) 	:: nucleus
		real(dp),intent(in)				:: scaleFactor
		integer,intent(in) 				:: N,Z

		COMPLEX(kind=qp),intent(out)			:: prod_N,prod_Z
		COMPLEX(dp),intent(out) 				:: U_N(N_tot,N_tot),V_N(N_tot,N_tot),&
								U_Z(N_tot,N_tot),V_Z(N_tot,N_tot)

		real(dp) 	:: Theta_N(N_tot),Theta_Z(N_tot),&
				lam_sN,lam_sZ,tol,EV_N(N_tot,2),EV_Z(N_tot,2)
		integer 	:: i,j,step

		step=100000 	! Nbr of points in \lambda vector
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

	FUNCTION prod_sqrt(V,N) result(prod)
		COMPLEX(dp) 		:: V(N,N) 
		integer 		:: N,i
		COMPLEX(kind=qp) 	:: prod
		
		prod=1
		do i=1,N,2
			if(i+1<=N) then
				!prod=prod*V(i,i+1)*V(i,i+1)
				prod=prod*V(i,i+1)
			end if
		end do
	END FUNCTION prod_sqrt

	FUNCTION sum_check(V,N) result(summ)
		COMPLEX(dp) 		:: V(N,N) 
		integer 		:: N,i
		!COMPLEX(kind=qp) 	:: prod
		real(kind=qp) 	:: summ
		
		summ=0
		do i=1,N
			if(i+1<=N) then
				!prod=prod*V(i,i+1)*V(i,i+1)
			!	write(*,*) V(i,i+1)
				summ=summ+V(i,i+1)**2
			end if
		end do
	END FUNCTION sum_check

	FUNCTION prod_calc(V,N) result(prod)
		COMPLEX(dp) 		:: V(N,N) 
		integer 		:: N,i
		!COMPLEX(kind=qp) 	:: prod
		real(kind=qp) 	:: prod
		
		prod=1
		do i=1,N,2
			if(i+1<=N) then
				!prod=prod*V(i,i+1)*V(i,i+1)
				prod=prod*V(i,i+1)**2
			end if
		end do
	END FUNCTION prod_calc

	FUNCTION norm_fac(U,V,N) result(norm)
		COMPLEX(dp) 	:: WW(2*N,2*N),U(N,N),V(N,N)
		integer 	:: N 
		COMPLEX(kind=dp)	:: 	Pf2P(2)
		REAL(kind=qp) 		:: 	norm

		WW=WTW(U,V,N)
		call ZSKPF10_F95(WW,Pf2P)
		norm=(-1)**(N*(N-1)/2)*real(PF2P(1))*ten_quad**real(PF2P(2))

	END FUNCTION norm_fac
	
	FUNCTION WTW(U,V,N) result(WW)
		COMPLEX(dp) 	:: WW(2*N,2*N),U(N,N),V(N,N)
		integer 	:: N 

		WW(1:N,1:N) = matmul(transpose(V),U)
		WW(N+1:2*N,1:N) = -matmul(transpose(V),V)
		WW(1:N,N+1:2*N) = matmul(transpose(V),V)
		WW(N+1:2*N,N+1:2*N) = matmul(transpose(U),V)
	END FUNCTION WTW

	SUBROUTINE ZSKPF10_OVERLAPP(Pf2P_out,N,UD,UM,VD,VM,DMAT,get_norm) 
		! computes the overlap of two wave functions \bra{UM,VM} DMAT \ket{UD,VD}   with ZSKPFA from PFAPACK,
		implicit none

		integer, intent(in) :: N
		COMPLEX(kind=dp),dimension(N,N), intent(in) :: DMAT, UD, VD, UM, VM
		REAL(kind=qp), intent(in), optional	:: 	get_norm

		COMPLEX(kind=dp) 	:: 	WTW(2*N,2*N)
		COMPLEX(kind=dp)	:: 	Pf2P(2),Pf2P_out
		REAL(kind=qp) 		:: 	norm
		integer 		:: 	i,j,S_n

		S_n=(-1)**(N*(N-1)/2)

		if( present(get_norm)) then
			norm=get_norm
		else
			norm=norm_fac(UD,VD,N)
		end if

		WTW(1:N,1:N) = matmul(transpose(VM),UM)
		WTW((N+1):(2*N),1:N) = -matmul(transpose(conjg(VD)),matmul(transpose(DMAT),VM))
		WTW(1:N,(N+1):2*N) = matmul(transpose(VM),matmul(DMAT,conjg(VD)))
		WTW((N+1):(2*N),(N+1):(2*N)) = matmul(transpose(conjg(UD)),conjg(VD))

		call ZSKPF10_F95(WTW,Pf2P) 
		Pf2P_out=S_n*real(PF2P(1))*ten_quad**real(PF2P(2))/(norm)

	END SUBROUTINE ZSKPF10_OVERLAPP

	
	FUNCTION DMAT_CREATOR(phi,N) result(DMAT)
		! computes the overlap of two wave functions \bra{UM,VM} DMAT \ket{UD,VD}   with ZSKPFA from PFAPACK,
		implicit none

		integer, intent(in):: N
		real(dp), intent(in) :: phi
		
		COMPLEX(dp), dimension(N,N) :: DMAT

		COMPLEX(dp) :: z_no
		integer :: i

		z_no = cmplx(0,phi,16)
		DMAT=0

		do i=1,N
			DMAT(i,i) = exp(z_no)
		end do 

	END FUNCTION DMAT_CREATOR
	
end module BCS
! Constructs the BCS-model and computes the system for \lambda=[1..-1] and plots the result versus number of particles
program main
	use BCS
	implicit none
! 	Nucleus, filled using MO_module
	type(Nucleon), dimension(:,:), allocatable :: nucleus
! 	Variables used for analytic solution
	COMPLEX(dp) :: 	Pf2P(2),Pf2P_ol,&
			U_N(N_tot,N_tot),V_N(N_tot,N_tot),&
			U_Z(N_tot,N_tot),V_Z(N_tot,N_tot),&
			DMAT(N_tot,N_tot),&
		       	WW_N(2*N_tot,2*N_tot), WW_D(2*N_tot,2*N_tot)

	COMPLEX(kind=qp) :: prod_N,prod_Z,summ,summ2
	REAL(kind=qp) 	:: norm


	REAL(kind=dp) 	:: scaleFactor(7),dPhi,S_n,N_exp,N_op


	! Nbr of neutrons/protons and loop integer(s),
	integer :: i,j,N_loop
	integer :: N,Z

	N=24		! Number of NEUTRONS to find
	Z=24 		! Number of PROTONS to find
	
	scaleFactor = (/0.5d0,0.75d0,1.d0,1.75d0,2.5d0,3d0,4d0/)

!-------Analytical solution using BCS-equations
	allocate(nucleus(N_tot,2))

	call nucleus_creator(N,Z,nucleus)

!-------Test-loop for the overlapp pfaffian, should be =1 for U=U', V='V, DMAT=1

	!DMAT=DMAT_CREATOR(0d0,N_tot)
	!do i=1,7
	!	call qpart_creator(nucleus,N,Z,scaleFactor(i),U_N,V_N,prod_N,U_Z,V_Z,prod_Z)
	!	call ZSKPF10_OVERLAPP(Pf2P_ol,N_tot,U_N,U_N,V_N,V_N,DMAT)
	!	write(*,*) 'Pairing strength scaled by : ',scaleFactor(i)
	!	write(*,*) 'Pfaffian_olverlap : ', real(Pf2P_ol)
	!	write(*,*) 
        !end do

!-------Test-loop for the particle nbr op exp value

	call qpart_creator(nucleus,N,Z,1d0,U_N,V_N,prod_N,U_Z,V_Z,prod_Z)
	open(unit=1,file='data/part_no_test.dat',status='replace')
	open(unit=2,file='data/phase_proj_test.dat',status='replace')
	
	norm=norm_fac(U_N,V_N,N_tot)

	N_loop=8*24
	dPhi=2*PI/(2*N_loop+1)	

	N_op=24
	N_exp=24

	do j=N_exp,N_exp!N_exp-2,N_exp+2
		write(*,*) j
		summ=0
		summ2=0
		do i=-1*N_loop,N_loop
			DMAT=DMAT_CREATOR(1*i*dPhi,N_tot)
			call ZSKPF10_OVERLAPP(Pf2P_ol,N_tot,U_N,U_N,V_N,V_N,DMAT,norm)
			write(2,'(I4,2ES48.38)') i,i*dPhi, real(Pf2P_ol)**2
			summ= summ + dPhi*exp(cmplx(0,-1*i*dPhi*j,16))*Pf2P_ol
			summ2= summ2 + dPhi*exp(cmplx(0,-1*i*dPhi*j,16))*exp(cmplx(0,1*dPhi*i*N_exp,16))
		end do
		summ=1/(2*PI)*summ
		summ2=1/(2*PI)*summ2
		write(1,'(I2,2ES48.38)') j,real(summ),real(summ2)
		!if(j .NE. N_exp+2) write(2,'(2/)')
		!write(2,*)
	end do

	deallocate(nucleus)
	close(unit=1)
	close(unit=2)
end program main
