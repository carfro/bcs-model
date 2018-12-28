module BCS
	use MO_module
	use analytic_module
	use pfaffian_module
	implicit none
	integer, parameter :: N_tot=330
contains

	Subroutine qpart_creator(nucleus,N,Z,scalefactor,U_N,V_N,prod_N,U_Z,V_Z,prod_Z)

		type(Nucleon), dimension(:,:),intent(in) 	:: nucleus
		real(8),intent(in)				:: scaleFactor
		integer,intent(in) 				:: N,Z

		complex(kind=qp),intent(out)			:: prod_N,prod_Z
		complex(8),intent(out) 				:: U_N(N_tot,N_tot),V_N(N_tot,N_tot),&
								U_Z(N_tot,N_tot),V_Z(N_tot,N_tot)

		real(8) 	:: Theta_N(N_tot),Theta_Z(N_tot),&
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

	FUNCTION WTW(U,V) result(WW)
		complex(8) 	:: WW(2*N_tot,2*N_tot),U(N_tot,N_tot),V(N_tot,N_tot)

		WW(1:N_tot,1:N_tot) = matmul(transpose(V),U)
		WW(N_tot+1:2*N_tot,1:N_tot) = -matmul(transpose(V),V)
		WW(1:N_tot,N_tot+1:2*N_tot) = matmul(transpose(V),U)
		WW(N_tot+1:2*N_tot,N_tot+1:2*N_tot) = matmul(transpose(U),V)
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
	COMPLEX(8) :: 	Pf2H(2),Pf2P(2),WW_N(2*N_tot,2*N_tot),&
			U_N(N_tot,N_tot),V_N(N_tot,N_tot),&
			U_Z(N_tot,N_tot),V_Z(N_tot,N_tot)

	complex(kind=qp) :: Pf22,Pf,prod_N,prod_Z
	real(kind=qp) 	:: factor

	Integer 	:: Ipiv(2*N_tot,2),Ipiv2(2*24,2),Ipiv22(2*24,2),Ipiv3(4,2)

	REAL(8) 	:: scaleFactor(7)

	! Nbr of neutrons/protons and loop integer(s),
	!integer :: i,j
	integer :: N,Z,siz,step

	N=24		! Number of NEUTRONS to find
	Z=24 		! Number of PROTONS to find
	
!-------Analytical solution using BCS-equations
	allocate(nucleus(N_tot,2))

	call nucleus_creator(N,Z,nucleus)

	scaleFactor = (/0.5d0,0.75d0,1.d0,1.75d0,2.5d0,3d0,4d0/)
	
	call qpart_creator(nucleus,N,Z,scalefactor(3),U_N,V_N,prod_N,U_Z,V_Z,prod_Z)
	WW_N = WTW(U_N,V_N)

	open(unit=1,file='data/performance.dat',status='replace')

	factor=1
	Pf=overlap_pfaffian(factor,N_tot,N_tot,N_tot,U_N,U_N,V_N,V_N)

	!call ZPfaffian_EXT(WW,2*N_tot,2*N_tot,Ipiv,Pf)

	call ZSKPF10_F95(WW_N,Pf2P) 
	call ZSKPF10_F95(WW_N,Pf2H,MTHD='H') 
	
	write(*,*) 'Product: ' 			, real(prod_N)
	write(*,*) 'Pf_Extended: ' 		, real(Pf)
	write(*,*) 'Pf_SKPF10, Parlett-Reid : ' , real(Pf2P)
	write(*,*) 'Pf_SKPF10, Householder : ' 	, real(Pf2H)

	close(unit=1)
	deallocate(nucleus)
end program main
