module BCS
	use MO_module
	use analytic_module
	use pfaffian_module
	USE F95_PFAPACK
	implicit none
	integer, parameter :: N_tot=330
<<<<<<< 39081309af0c5314a057e32da0b2ba81764cc1d9
<<<<<<< ed64cb91f3bfcf01ea739c71ab7d2bf93469f696
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
<<<<<<< 4a8e749bb2e5597a330efac1e9fb799123643e5d
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
	
=======
=======
	real(8), parameter :: PI=4.D0*DATAN(1.D0)
	real(kind=qp), parameter :: ten_quad = 10.

>>>>>>> The pfaffian performance test is now finished. Currently in the process of creating particle nbr projector
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

	FUNCTION WTW(U,V,N) result(WW)
		COMPLEX(dp) 	:: WW(2*N,2*N),U(N,N),V(N,N)
		integer 	:: N 

		WW(1:N,1:N) = matmul(transpose(V),U)
		WW(N+1:2*N,1:N) = -matmul(transpose(V),V)
		WW(1:N,N+1:2*N) = matmul(transpose(V),V)
		WW(N+1:2*N,N+1:2*N) = matmul(transpose(U),V)
	END FUNCTION WTW

<<<<<<< 39081309af0c5314a057e32da0b2ba81764cc1d9
>>>>>>> Last commit was uncomplete. Added libs and scrapbook containing old cod that might be useful
=======
	FUNCTION ZSKPF10_OVERLAPP(N,UD,UM,VD,VM,DMAT) result(Pf2P_out)
=======
>>>>>>> Trying to figure out why the part. no. op. is not working. A check using the onishi formula might be usefull. Printing results to visualize, will remove later.
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
	
>>>>>>> The pfaffian performance test is now finished. Currently in the process of creating particle nbr projector
end module BCS
! Constructs the BCS-model and computes the system for \lambda=[1..-1] and plots the result versus number of particles
program main
	use BCS
	implicit none
! 	Nucleus, filled using MO_module
	type(Nucleon), dimension(:,:), allocatable :: nucleus
! 	Variables used for analytic solution
<<<<<<< 39081309af0c5314a057e32da0b2ba81764cc1d9
<<<<<<< 56f1d0abafcb13c8470f19f2d610ccb98fef0fc6
<<<<<<< f8ecfa0a313073c6540362884910bffeee3ffa04
<<<<<<< ed64cb91f3bfcf01ea739c71ab7d2bf93469f696
	COMPLEX(dp) :: 	Pf2P(2),Pf2P_ol,&
			U_N(N_tot,N_tot),V_N(N_tot,N_tot),&
			U_Z(N_tot,N_tot),V_Z(N_tot,N_tot),&
			DMAT(N_tot,N_tot),&
		       	WW_N(2*N_tot,2*N_tot), WW_D(2*N_tot,2*N_tot)

	COMPLEX(kind=qp) :: prod_N,prod_Z,summ,summ2
	REAL(kind=qp) 	:: norm
<<<<<<< 4a8e749bb2e5597a330efac1e9fb799123643e5d


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
=======
	COMPLEX(8) :: 	Pf2,WW_N(2*N_tot,2*N_tot),&
=======
	COMPLEX(8) :: 	Pf2H(2),Pf2P(2),WW_N(2*N_tot,2*N_tot),&
>>>>>>> Finally managed to get the linking of the pfapack-library to workgit status!
=======
	COMPLEX(8) :: 	Pf2H(2),Pf2P(2),&
>>>>>>> The test of pfaffian routines aswell as plots of their abs. error and runtimes have now been made. The error measure is not a very good one, will look into making it a relative one - right now the number doesn't say much. It would also be nice to run precision measurements on more realistic wavefunctions W.
=======
	COMPLEX(dp) :: 	Pf2P(2),Pf2P_ol,&
>>>>>>> The pfaffian performance test is now finished. Currently in the process of creating particle nbr projector
			U_N(N_tot,N_tot),V_N(N_tot,N_tot),&
			U_Z(N_tot,N_tot),V_Z(N_tot,N_tot),&
			DMAT(N_tot,N_tot),&
		       	WW_N(2*N_tot,2*N_tot), WW_D(2*N_tot,2*N_tot)

	COMPLEX(kind=qp) :: prod_N,prod_Z,summ
	!REAL(dp) 	:: Pf2P_out
=======
>>>>>>> Trying to figure out why the part. no. op. is not working. A check using the onishi formula might be usefull. Printing results to visualize, will remove later.


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

<<<<<<< 39081309af0c5314a057e32da0b2ba81764cc1d9
	prod_test=prod_calc(V_N,N_tot)
	
<<<<<<< 56f1d0abafcb13c8470f19f2d610ccb98fef0fc6
<<<<<<< f8ecfa0a313073c6540362884910bffeee3ffa04
	write(*,*) 'Product: ', real(prod_N)
	write(*,*) 'Pf_Extended: ', real(Pf)
	write(*,*) 'Pf_SKPF10: ', real(Pf2)
>>>>>>> Last commit was uncomplete. Added libs and scrapbook containing old cod that might be useful
=======
	write(*,*) 'Product: ' 			, real(prod_N)
	write(*,*)
	write(*,*) 'Pf_Extended: ' 		, real(Pf)
<<<<<<< 16dc980778b554c97900b8315ec9ef82fcbfb481
	write(*,*) 'Pf_SKPF10, Parlett-Reid : ' , real(Pf2P)
	write(*,*) 'Pf_SKPF10, Householder : ' 	, real(Pf2H)
>>>>>>> Finally managed to get the linking of the pfapack-library to workgit status!
=======
	write(*,*) 'elapsed time: ', t2-t1
	write(*,*)
	write(*,*) 'Pf_SKPF10, Parlett-Reid : ' , real(Pf2P) 	
	write(*,*) 'elapsed time: ', t3-t2
	write(*,*)
	write(*,*) 'Pf_SKPF10, Householder : ' 	, real(Pf2H) 	
	write(*,*) 'elapsed time: ', t4-t3
>>>>>>> Too tiered to get anything done, simply added some lines for timing the different pfaffian routines.
=======
!-------Test-loop for the pfaffians below	

	open(unit=1,file='data/performance_err.dat',status='replace')
	open(unit=2,file='data/performance_time.dat',status='replace')

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
		call ZSKPF10_F95(WW_N,Pf2H,MTHD='H') 
		call cpu_time(t4)
		prod_test=prod_calc(V_test,N_mult)
		call cpu_time(t5)

		write(prod_str,'(ES43.36)') prod_test
		!prod_str=prod_str(1:18)//'E'//prod_str(index(prod_str(3:len(prod_str)),'-')+2:len_trim(prod_str))
		!read(prod_str,*)  prod_test

		write(char1,'(F39.36)') real(Pf2P(1)); write(char2,'(I4)') int(real(Pf2P(2)))
		pf2p_str=char1//'E'//trim(adjustl(char2))
		read(pf2p_str,*) pf2p_real
		pf2p_real=abs(abs(Pf2p_real)-prod_test)/prod_test
		write(pf2p_str,'(ES43.36)') pf2p_real

		!write(char1,'(F39.36)') real(Pf2h(1)); write(char2,'(I4)') int(real(Pf2h(2)))
		!pf2h_str=char1//'E'//trim(adjustl(char2))
		!read(pf2h_str,*) pf2h_real
		pf2h_real=real(Pf2h(1))*10**real(Pf2h(2))
		pf2h_real=abs(abs(Pf2h_real)-prod_test)/prod_test
		write(*,*) pf2h_real
	 	write(pf2h_str,*) pf2h_real
		write(*,*) pf2h_str
		!write(*,*) pf2h_str(1:index(pf2h_str,'E')-1)
		!pf2h_str=pf2h_str(1:39)//'E'//pf2h_str(index(pf2h_str,'E')+1:len_trim(pf2h_str))
		!write(*,*) pf2h_str
		!pf2h_str=char1//'E'//trim(adjustl(char2))
		!write(pf2h_str,'(ES43.36)') pf2h_real
	!	
		Pf=abs(abs(Pf)-prod_test)/prod_test
		write(pf_str,'(ES43.36)') real(Pf)

! Writes product calculation, and pfaffian using: robledo, pfa_parlett-reid, pfa_householder
		!write(*,'(I4,ES51.36)') N_mult , prod_test
		write(*,'(I4,3A45,A50)') N_mult , prod_str,pf_str, pf2p_str, pf2h_str
		write(1,'(I4,3A45,A50)') N_mult , prod_str,pf_str, pf2p_str, pf2h_str
		!write(1,'(I4,5A,3(ES48.38,A5))') N_mult ,char(9), prod_test,char(9), abs(abs(Pf)-prod_test)&
		!	, char(9),abs(abs(pf2h_real)-prod_test), char(9)!,abs(abs(Pf2p_real)-prod_test)

		write(*,*)
! Writes time of product calculation, robledo, pfa_parlett-reid, pfa_householder
		write(2,'(I4,4E24.16)') N_mult, t5-t4, t2-t1, t3-t2, t4-t3 
		deallocate(U_test,V_test,WW_N)
=======
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
<<<<<<< 4a8e749bb2e5597a330efac1e9fb799123643e5d
		write(*,*) summ
>>>>>>> The pfaffian performance test is now finished. Currently in the process of creating particle nbr projector
=======
		summ=1/(2*PI)*summ
		summ2=1/(2*PI)*summ2
		write(1,'(I2,2ES48.38)') j,real(summ),real(summ2)
		!if(j .NE. N_exp+2) write(2,'(2/)')
		!write(2,*)
>>>>>>> Trying to figure out why the part. no. op. is not working. A check using the onishi formula might be usefull. Printing results to visualize, will remove later.
	end do
>>>>>>> The test of pfaffian routines aswell as plots of their abs. error and runtimes have now been made. The error measure is not a very good one, will look into making it a relative one - right now the number doesn't say much. It would also be nice to run precision measurements on more realistic wavefunctions W.

<<<<<<< 39081309af0c5314a057e32da0b2ba81764cc1d9
	deallocate(nucleus)
	close(unit=1)
<<<<<<< 56f1d0abafcb13c8470f19f2d610ccb98fef0fc6
<<<<<<< ed64cb91f3bfcf01ea739c71ab7d2bf93469f696
	close(unit=2)
=======
=======
	close(unit=2)
>>>>>>> The test of pfaffian routines aswell as plots of their abs. error and runtimes have now been made. The error measure is not a very good one, will look into making it a relative one - right now the number doesn't say much. It would also be nice to run precision measurements on more realistic wavefunctions W.
=======
>>>>>>> The pfaffian performance test is now finished. Currently in the process of creating particle nbr projector
	deallocate(nucleus)
<<<<<<< 4a8e749bb2e5597a330efac1e9fb799123643e5d
>>>>>>> Last commit was uncomplete. Added libs and scrapbook containing old cod that might be useful
=======
	close(unit=1)
	close(unit=2)
>>>>>>> Trying to figure out why the part. no. op. is not working. A check using the onishi formula might be usefull. Printing results to visualize, will remove later.
end program main
