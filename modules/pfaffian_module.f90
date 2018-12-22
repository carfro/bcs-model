! Pfaffian routines used: 
module pfaffian_module 
	implicit none
	
	! 128-bit real w/ 33 sig. fig, exponent range 4931 
	integer, parameter 	:: qp = selected_real_kind(33, 4931)
	
contains 
	FUNCTION overlap_pfaffian(lambda,NDIM,dimD,dimM,UD,UM,VD,VM) result(pf)
	    ! computes the overlap of two wave functions with pfaffian formula
	    implicit none

	    integer, intent(in) :: NDIM, dimD, dimM
	    complex*16 :: W((dimM+dimD),(dimM+dimD))
	    complex*16, intent(in) :: UD(NDIM,dimD), VD(NDIM,dimD), UM(NDIM,dimM), VM(NDIM,dimM)
	    integer :: IPIV((dimM+dimD),2)
	    complex(kind=qp) :: pf
	    real(kind=qp),intent(in) :: lambda

	    W(1:dimM,1:dimM) = matmul(transpose(VM),UM)
	    W(1:dimM,(dimM+1):(dimM+dimD)) = matmul(transpose(VM),conjg(VD))
	    W((dimM+1):(dimM+dimD),1:dimM) = -matmul(transpose(conjg(VD)),VM)
	    W((dimM+1):(dimM+dimD),(dimM+1):(dimM+dimD)) = matmul(transpose(conjg(UD)),conjg(VD))

	    W=lambda*W

	    call Zpfaffian_ext(W,(dimM+dimD),(dimM+dimD),IPIV,pf)
	    return
	END FUNCTION overlap_pfaffian

	Subroutine ZPfaffian_EXT(SK,LDS,N,Ipiv,Pf)
	    ! CALCULATES THE PFAFFIAN OF A SKEW-SYMMETRIC MATRIX SK USING EXTENDED PRECISION
	    Implicit none

	    Integer NB,IB,i,j,k,ip,NR,NC,I1,I2
	    Integer, intent(in)::N,LDS
	    Integer, intent(inout), dimension(N,2):: Ipiv
	    real*16 :: epsln,phas,big
	    Complex*16 :: SS,SW !,one,zero
	    complex(kind=qp), intent(inout)::Pf
	    Complex*16, intent(inout),dimension(LDS,N):: SK
	    complex(kind=qp) :: one,zero

	    one=(1.0_qp,0.0_qp)
	    zero=(0.0_qp,0.0_qp)

	    epsln = 1.0d-13   ! smallest number such that 1+epsln=1

	    if(mod(N,2) == 1) then ! N odd

	       Pf = zero

	    else

	       NB= N / 2   ! Number of 2x2 blocks

	       if (NB==1) then

		  ! The pfaffian of a 2x2 matrix is Pf=Sk(1,2)

		  Pf = SK(1,2)

	       else   ! NB>1

		  Pf = one

		  do IB=1,NB-1

		     NR = IB*2 - 1 ! row numb of the 1,2 element of the 2x2 block
		     NC = NR + 1

		     big = 0.0d+00
		     !big=0.0_

		     I1 = NR
		     I2 = NC

		     do i=NR,N-1 ! all rows
			do j=i+1,N

			   if(abs(SK( i , j ))> big) then
			      big = abs(SK( i , j ))
			      I1 = i
			      I2 = j
			   end if

			end do
		     end do

		     Ipiv( IB , 1) = I1 ! to initialize
		     Ipiv( IB , 2) = I2 ! to initialize

		     phas = 1.0d+00
		     if(I1==NR) phas = -phas
		     if(I2==NC) phas = -phas

		     ! Pivoting of element NR,NC with ip1,ip2

		     !
		     !   Pivoting for a skew-symmetric matrix (Upper)
		     !
		     if(I1/=NR) Call Zexch_ext(SK,LDS,N,NR,I1)
		     if(I2/=NC) Call Zexch_ext(SK,LDS,N,NC,I2)

		     ss = Sk( NR , NC )

		     !       write(6,*) ' IB I2, I2 ', IB,I1,I2,big,ss

		     if(abs(ss)>epsln) then
			!
			!      Updating the Schur complement  matrix
			!

			do i = 2*IB+1 , N-1

			   do j = i+1 , N

			      Sk(i,j) = Sk(i,j) + (Sk(NC,i)*Sk(NR,j)-Sk(NR,i)*Sk(NC,j))/SS
			   end do
			end do
			!
			! Storing X and Y vectors in the lower part of the matrix
			!
			do i = 2*IB+1 , N

			   Sk( i , NR ) =  -Sk( NC , i )/SS
			   Sk( i , NC ) =   Sk( NR , i )/SS

			end do

			if (IB>=2) then  ! swap

			   !
			   !   Swapping
			   !
			   do j= 1 , 2*IB

			      SW           = Sk( NR , j )
			      Sk( NR , j ) = Sk( I1 , j )
			      Sk( I1 , j ) = SW

			      SW           = Sk( NC , j )
			      Sk( NC , j ) = Sk( I2 , j )
			      Sk( I2 , j ) = SW

			   end do

			end if

		     else

			Pf = zero

			return

		     end if ! dabs(ss)>epsln

		     Pf = Pf * SS * phas

		  end do ! IB

		  Pf = Pf * Sk(N-1,N)

	       end if !    NB==1

	    end if !    mod(N,2)==1

	    return
	end Subroutine ZPfaffian_EXT

	Subroutine ZExch_ext(SK,LDS,N,I1,I2)
	    !  Exchange of rows i1 i2 and columns i1 i2
	    !  SK is assumed upper triangular
	    !  i1 < i2

	    Implicit none

	    Integer,intent(in)::LDS,N,I1,I2
	    Integer :: I
	    Complex*16,intent(inout),dimension(LDS,N):: SK
	    Complex*16 :: SS

	    SK(i1,i2) = -SK(i1,i2)

	    if(I1/=1) then

	       do i=1,I1-1
		  SS           = SK( i , I1 )
		  SK( i , I1 ) = SK( i , I2 )
		  SK( i , I2 ) = SS
	       end do
	    end if

	    if(I2/=N) then

	       do i=I2+1,N
		  SS           = SK( I1 , i )
		  SK( I1 , i ) = SK( I2 , i )
		  SK( I2 , i ) = SS
	       end do
	    end if

	    if(I2>=I1+2) then

	       do i=I1+1,I2-1
		  SS           =  SK( I1 , i )
		  SK( I1 , i ) = -SK( i , I2 )
		  SK( i , I2 ) = -SS
	       end do
	    end if

	    return
	end Subroutine ZExch_ext
!+---------------------------------------------------------------------+
!|   Computes the Pfaffian of a skew-symmetric matrix SK using         |
!|   the ideas of Aitken's block diagonalization formula.              |
!|   Full pivoting is implemented to make the whole procedure stable   |
!+---------------------------------------------------------------------+
!|   Only the upper triangle of SK is used                             |
!|   In output the lower triangle contains the transformation matrix   |
!+---------------------------------------------------------------------+
!|   SK .....Skew symmetric input matrix                               |
!|   LDS ....Leading dimension of matrix SK                            |
!|   N  .....Dimension of SK                                           |
!|   Ipiv ...Integer vector to hold the pivoted columns                |
!+---------------------------------------------------------------------+
!| This program may be freely used providing such use cites the source,|
!|      Numeric and symbolic evaluation of the pfaffian of general     |
!|      skew-symmetric matrices                                        |
!|                                                                     |
!|      C. Gonzalez-Ballestero, L.M.Robledo and G.F. Bertsch           |
!|      Computer Physics Communications (2011)                         |
!+---------------------------------------------------------------------+
!
!|   C o m p l e x    V e r s i o n                                   |
!|                                                                    |
!+--------------------------------------------------------------------+
	Subroutine ZPfaffianF(SK,LDS,N,Ipiv,Pf)

		Implicit none

		Integer NB,IB,i,j,k,ip,NR,NC,I1,I2
		Integer, intent(in)::N,LDS
		Integer, intent(inout), dimension(N,2):: Ipiv
		Real*8 epsln,phas,big
		Double Complex SS,SW,one,zero
		Double Complex, intent(inout)::Pf
		Double Complex, intent(inout),dimension(LDS,N):: SK

		one = dcmplx(1.0d+00,0.0d+00)
		zero= dcmplx(0.0d+00,0.0d+00)

		epsln = 1.0d-13   ! smallest number such that 1+epsln=1

		if(mod(N,2) == 1) then ! N odd

		Pf = zero

		else

		NB= N / 2   ! Number of 2x2 blocks

		if (NB==1) then

		! The pfaffian of a 2x2 matrix is Pf=Sk(1,2)

		Pf = SK(1,2)

		else   ! NB>1

		Pf = one

		do IB=1,NB-1

		   NR = IB*2 - 1 ! row numb of the 1,2 element of the 2x2 block
		   NC = NR + 1

		   big = 0.0d+00

		   I1 = NR
		   I2 = NC

		   do i=NR,N-1 ! all rows
		      do j=i+1,N

			 if(abs(SK( i , j ))> big) then
			    big = abs(SK( i , j ))
			    I1 = i
			    I2 = j
			 end if

		      end do
		   end do

		   Ipiv( IB , 1) = I1 ! to initialize
		   Ipiv( IB , 2) = I2 ! to initialize

		   phas = 1.0d+00
		   if(I1==NR) phas = -phas
		   if(I2==NC) phas = -phas

		! Pivoting of element NR,NC with ip1,ip2

		!
		!   Pivoting for a skew-symmetric matrix (Upper)
		!
		   if(I1/=NR) Call Zexch(SK,LDS,N,NR,I1)
		   if(I2/=NC) Call Zexch(SK,LDS,N,NC,I2)

		   ss = Sk( NR , NC )

		!       write(6,*) ' IB I2, I2 ', IB,I1,I2,big,ss

		   if(abs(ss)>epsln) then
		!
		!      Updating the Schur complement  matrix
		!

		   do i = 2*IB+1 , N-1

		      do j = i+1 , N

			 Sk(i,j) = Sk(i,j) + (Sk(NC,i)*Sk(NR,j)-Sk(NR,i)*Sk(NC,j))/SS
		      end do
		   end do
		!
		! Storing X and Y vectors in the lower part of the matrix
		!
		   do i = 2*IB+1 , N

		      Sk( i , NR ) =  -Sk( NC , i )/SS
		      Sk( i , NC ) =   Sk( NR , i )/SS

		   end do

		   if (IB>=2) then  ! swap

		!
		!   Swapping
		!
		   do j= 1 , 2*IB

		      SW           = Sk( NR , j )
		      Sk( NR , j ) = Sk( I1 , j )
		      Sk( I1 , j ) = SW

		      SW           = Sk( NC , j )
		      Sk( NC , j ) = Sk( I2 , j )
		      Sk( I2 , j ) = SW

		   end do

		   end if

		   else

		   Pf = zero

		   return

		   end if ! dabs(ss)>epsln

		   Pf = Pf * SS * phas

		   end do ! IB

		Pf = Pf * Sk(N-1,N)

		end if !    NB==1

		end if !    mod(N,2)==1

		return
	end Subroutine ZPfaffianF
	!
	!  Exchange of rows i1 i2 and columns i1 i2
	!  SK is assumed upper triangular
	!  i1 < i2
	!
	Subroutine ZExch(SK,LDS,N,I1,I2)
		Implicit none

		Integer,intent(in)::LDS,N,I1,I2
		Integer :: I
		Double Complex,intent(inout),dimension(LDS,N):: SK
		Double Complex :: SS

		SK(i1,i2) = -SK(i1,i2)

		if(I1/=1) then

		   do i=1,I1-1
		      SS           = SK( i , I1 )
		      SK( i , I1 ) = SK( i , I2 )
		      SK( i , I2 ) = SS
		   end do
		end if

		if(I2/=N) then

		   do i=I2+1,N
		      SS           = SK( I1 , i )
		      SK( I1 , i ) = SK( I2 , i )
		      SK( I2 , i ) = SS
		   end do
		end if

		if(I2>=I1+2) then

		   do i=I1+1,I2-1
		      SS           =  SK( I1 , i )
		      SK( I1 , i ) = -SK( i , I2 )
		      SK( i , I2 ) = -SS
		   end do
		end if

		return
	end Subroutine ZExch
!+---------------------------------------------------------------------+
!|   Computes the Pfaffian of a skew-symmetric matrix SK using         |
!|   Householder transformations to bring the matrix  to tridiagonal   |
!|   form. The pfaffian of a tridiagonal skew-symmetric matrix is      |
!|   straitforward to compute.                                         |
!+---------------------------------------------------------------------+
!|   The full matrix SK is used                                        |
!+---------------------------------------------------------------------+
!+---------------------------------------------------------------------+
!|   Complex version                                                   |
!+---------------------------------------------------------------------+
!|   SK ....Skew symmetric input matrix                                |
!|   LDS ...Leading dimension of matrix SK                             |
!|   N  ....Dimension of SK                                            |
!|   Pf  ...Pfaffian                                                   |
!|   WORK ..Working matrix of dimension LDS,N                          |
!|   IH.....Switch to use full (1) or half (2) Householder             |
!+---------------------------------------------------------------------+
!|   BLAS subroutines are used                                         |
!|                                                                     |
!|   zcopy .... to copy vectors and matrices                           |
!|   zdot  .... scalar product of two vectors                          |
!|   zgemv .... Matrix vector multiplication                           |
!|   zger  .... Matrix plus diadic product of two vectors              |
!+---------------------------------------------------------------------+

!	Subroutine ZPfaffianH (SK,LDS,N,Pf,WORK,IH)
!		Implicit none
!
!		Integer IS,M,i,j,k,IDEN,IH
!		Integer,intent(in)::LDS,N
!		Real*8 epsln,xmod,u2
!		Double Complex ki,one,zero,al
!		Double Complex,intent(inout)::Pf
!		Double Complex,intent(inout),dimension(LDS,N):: SK
!		Double Complex,intent(inout),dimension(LDS,N):: WORK
!		Double Complex ZDOTC
!
!		one = dcmplx(1.0d+00,0.0d+00)
!		zero= dcmplx(0.0d+00,0.0d+00)
!
!		epsln = 1.0d-13   ! smallest number such that 1+epsln=1
!
!		if(IH.ne.1.and.IH.ne.2) then
!		   write(6,'(" Value of IH undefined in pfaffianH",I2)') IH
!		end if
!
!		if(mod(N,2)==1) then 
!		!                                                   N odd   Pf=0              
!		Pf = zero 
!
!		else
!		!                                                   N even                               
!		if (N==2) then 
!
!		!                          The pfaffian of a 2x2 matrix is Pf=Sk(1,2)
!		       
!		Pf = Sk(1,2)
!
!		else   ! N.gt.2
!
!		!
!		!                                      N even and > 2       
!		Pf = one
!		IDen=0
!
!		do IS=1,N-2,IH
!
!		M = N - IS
!
!		! Column N-IS+1 of Sk is copied into column IS of working vector
!		! This is vector x with dimension M of the Householder transformation   
!		    
!		call zcopy(M,Sk(1,N-IS+1),1,WORK(1,IS),1) 
!
!		! Norm of x
!
!		xmod = dsqrt( abs( ZDOTC(M,WORK(1,IS),1,WORK(1,IS),1) ) ) ! |x|
!
!		! if Norm of x=0 the pfaffian is zero
!		 
!		if (xmod>epsln) then
!
!		if(abs(Sk(M,M+1))>epsln) then
!
!		!   u = x -+ sign(Sk(N-IS,N-IS+1))*xmod e (N-IS)
!		!   e (M) is the unit cartesian vector with 1 in position M and
!		!   zero elsewhere
!
!		ki    = xmod*Sk(M,M+1)/abs(Sk(M,M+1))
!
!		else
!
!		ki    = dcmplx(xmod,0.0d+00)
!
!		end if
!
!		!
!		!      u = x -+ e**[i*arg(Sk(N-IS,N-IS+1))]*xmod e (N-IS)
!		!       
!		work(M,IS) = work(M,IS) - ki
!
!		!   Norm of u
!			      
!		u2 = abs(ZDOTC(M,WORK(1,IS),1,WORK(1,IS),1))  ! |u|**2
!
!		!  If the norm of u is too small the negative sign in the definition
!		!  of u is taken
!		     
!		if (u2 < epsln**2) then 
!		       
!		!  ki changes sign
!		 
!		ki    = -ki
!
!		work(M,IS) = work(M,IS) - 2.0d+00*ki
!
!		u2 = abs(ZDOTC(M,WORK(1,IS),1,WORK(1,IS),1))  ! |u|**2
!
!		end if
!
!		!     *
!		!    u  --> work(1,N-1)
!
!		do i=1,M
!		   work( i , N-1 ) = dconjg(work( i , IS ))
!		end do
!		       
!		!                *
!		!    v = A(N-1) u --> work(1,N)
!
!			 
!		call Zgemv('N',M,M,one,Sk,LDS,work(1,N-1),1,zero,work(1,N),1)  
!
!
!		!      Update of A(N-IS)
!		       
!		al = dcmplx(2.0d+00/u2,0.0d+00)
!
!		!     Sk -> Sk + al*u v'        
!
!		call zgeru(M,M, al,work(1,IS),1,work(1,N),1,Sk,LDS)
!
!		!     Sk -> Sk - al*v u'         
!
!		call zgeru(M,M,-al,work(1,N),1,work(1,IS),1,Sk,LDS)
!
!		else ! xmod.gt.epsln
!		       
!		ki = 0.0d+00
!		Iden = Iden + 1
!		       
!		end if   ! xmod.gt.epsln  
!			    
!		if(mod(IS,2)==1) Pf = Pf * ki
!
!		end do ! IS
!
!		Pf = Pf * Sk(1,2) * dfloat(1-2*mod(Iden,2))
!
!		if (IH==2) Pf = Pf * float(1-2*mod(N/2-1,2))
!
!		end if !    N==2
!		   
!		end if !    mod(N,2)==1
!
!
!		return
!	end
	
end module pfaffian_module		
