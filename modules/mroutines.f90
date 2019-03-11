module mroutines 
	use MO_module, only: dp, qp, ten_quad
	implicit none
	private 
	public det, det_eig, inv,inverse, mat_norm2, print_matrix,diag_inv
contains
	subroutine inverse(a,c,n)
		!============================================================
		! Inverse matrix
		! Method: Based on Doolittle LU factorization for Ax=b
		! Alex G. December 2009
		!-----------------------------------------------------------
		! input ...
		! a(n,n) - array of coefficients for matrix A
		! n      - dimension
		! output ...
		! c(n,n) - inverse matrix of A
		! comments ...
		! the original matrix a(n,n) will be destroyed 
		! during the calculation
		!===========================================================
		implicit none 
		integer 		:: 	n
		complex(kind=dp) 	:: 	a(n,n), c(n,n)
		complex(kind=dp) 	:: 	L(n,n), U(n,n), b(n), d(n), x(n)
		complex(kind=dp) 	:: 	coeff
		integer 		:: 	i, j, k

		! step 0: initialization for matrices L and U and b
		! Fortran 90/95 allows such operations on matrices
		L=0.0
		U=0.0
		b=0.0

		! step 1: forward elimination
		do k=1, n-1
		   do i=k+1,n
		      coeff=a(i,k)/a(k,k)
		      L(i,k) = coeff
		      do j=k+1,n
			 a(i,j) = a(i,j)-coeff*a(k,j)
		      end do
		   end do
		end do

		! Step 2: prepare L and U matrices 
		! L matrix is a matrix of the elimination coefficient
		! + the diagonal elements are 1.0
		do i=1,n
		  L(i,i) = 1.0
		end do
		! U matrix is the upper triangular part of A
		do j=1,n
		  do i=1,j
		    U(i,j) = a(i,j)
		  end do
		end do

		! Step 3: compute columns of the inverse matrix C
		do k=1,n
		  b(k)=1.0
		  d(1) = b(1)
		! Step 3a: Solve Ld=b using the forward substitution
		  do i=2,n
		    d(i)=b(i)
		    do j=1,i-1
		      d(i) = d(i) - L(i,j)*d(j)
		    end do
		  end do
		! Step 3b: Solve Ux=d using the back substitution
		  x(n)=d(n)/U(n,n)
		  do i = n-1,1,-1
		    x(i) = d(i)
		    do j=n,i+1,-1
		      x(i)=x(i)-U(i,j)*x(j)
		    end do
		    x(i) = x(i)/u(i,i)
		  end do
		! Step 3c: fill the solutions x(n) into column k of C
		  do i=1,n
		    c(i,k) = x(i)
		  end do
		  b(k)=0.0
		end do
	end subroutine inverse
	
	function diag_inv(N,mat) result(matinv)
		integer, intent(in)  		:: N 
		complex(kind=dp), intent(in) 	:: mat(N,N)
		complex(kind=dp) 		:: matinv(N,N)
		integer(kind=dp) 		:: i 

		if (abs(mat(1,2))>ten_quad**-4) stop 'Matrix not diagonal!'

		matinv 	= 0
		    
		do i=1,N
			matinv(i,i) = mat(i,i)**-1
		end do 

	end function

	function inv(N,mat) result(matinv)
		integer, intent(in)  		:: N 
		complex(kind=dp), intent(in) 	:: mat(N,N)
		complex(kind=dp) 		:: matinv(N,N), work(N)
		integer(kind=dp) 		:: i, info, ipiv(N)
		!real(kind=qp) 			:: sgn

		matinv 	= mat
		ipiv 	= 0
		info 	= 0
		    
		call zgetrf(N, N, matinv, N, ipiv, info)
		if (info.ne.0) stop 'Matrix is numerically singular!' 

		call zgetri(N, matinv, N, ipiv, work, N, info)
		if (info.ne.0) stop 'Matrix inversion failed!'

	end function

	complex(kind=qp) function det_eig(N, mat_in)
		! Computes the determinant of a matrix mat_NxN using eig.val. decomposition (ZGEEV) from lapack
		integer, intent(in) :: N 
		complex(kind=qp), intent(inout), dimension(:,:) :: mat_in
		integer 	   ::          lda, ldvl, ldvr
            	integer, parameter ::         lwmax=3000

		integer(kind=dp) :: i,j, info,lwork
		integer, allocatable :: ipiv(:)
		real(kind=qp) :: sgn

		real(kind=dp), allocatable :: rwork(:)
		complex(kind=dp), allocatable :: vl(:,:), vr(:,:), w(:), work(:),mat_used(:,:)
		
      		lda = N; ldvl = N; ldvr = N
		allocate(ipiv(N),mat_used(N,N),rwork(2*N),vl(ldvl,N),vr(ldvr,N),w(N),work(lwmax))
		ipiv = 0
		
 		mat_used=mat_in
		    
		info=0
		lwork=-1
		call zgeev('vectors','vectors',N,mat_used,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
		lwork=min(lwmax,int(work(1)))

		call zgeev('vectors','vectors',N,mat_used,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)

		if(info.gt.0) then
			write(*,*) 'the algorithm failed to compute eigenvalues.'
			write(*,*) 'info 	',info
			stop
		end if
		
		det_eig = 1._qp

		do i = 1, N
			det_eig = det_eig*w(i)
		end do

	end function det_eig

	complex(kind=qp) function det(N, mat)
		! Computes the determinant of a matrix mat_NxN using LU-decomposition from lapack
		integer, intent(in) :: N 
		complex(kind=qp), intent(in), dimension(:,:) :: mat
		complex(kind=qp),allocatable, dimension(:,:) :: mat_used
		integer(kind=dp) :: i, info
		integer, allocatable :: ipiv(:)
		real(kind=qp) :: sgn

		allocate(ipiv(N),mat_used(N,N))
		mat_used=mat
		ipiv = 0
		    
		call zgetrf(N, N, mat_used, N, ipiv, info)
		    
		det = 1._qp

		do i = 1, N
			det = det*mat_used(i, i)
			if(abs(mat_used(i,i))<10._qp**(-4) ) write(*,*) real(mat_used(i,i))
		end do

		sgn = 1._qp

		do i = 1, N
			if(ipiv(i) /= i) sgn = -sgn
		end do

		det = sgn*det   

	end function det

!     Auxiliary routine: printing a matrix.
!
	subroutine PRINT_MATRIX( DESC, M, N, A, LDA )
	      CHARACTER*(*)    DESC
	      INTEGER          M, N, LDA
	      COMPLEX*16          A( LDA, * )
	
	      INTEGER          I, J
	
	      WRITE(*,*)
	      WRITE(*,*) DESC
	      DO I = 1, M
		 WRITE(*,9998) ( A( I, J ), J = 1, N )
	      END DO
	
	      !write(FMT,*) 
	 9998 FORMAT( 4(:,1X,'(',F10.6,',',F10.6,')') )
	      RETURN
	end subroutine PRINT_MATRIX

!
!     Auxiliary routine: printing a matrix.
!
!	subroutine PRINT_ZMATRIX( DESC, M, N, A, LDA )
!	      CHARACTER*(*)    DESC
!	      INTEGER          M, N, LDA
!	      COMPLEX(kind=qp)          A( LDA, * )
!	
!	      INTEGER          I, J
!	
!	      WRITE(*,*)
!	      WRITE(*,*) DESC
!	      DO I = 1, M
!		 WRITE(*,9998) ( A( I, J ), J = 1, N )
!	      END DO
!	
!	 9998 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
!	      RETURN
!	end subroutine PRINT_ZMATRIX

!
!     Auxiliary routine: printing a real matrix.
!
	subroutine PRINT_REAL_MATRIX( DESC, M, N, A, LDA )
	      CHARACTER*(*)    DESC
	      INTEGER          M, N, LDA
	      REAL(8)             A( LDA, * )
	
	      INTEGER          I, J
	
	      WRITE(*,*)
	      WRITE(*,*) DESC
	      DO I = 1, M
		 WRITE(*,9998) ( A( I, J ), J = 1, N )
	      END DO
	
	 9998 FORMAT( 11(:,1X,F6.2) )
	      RETURN
	end subroutine PRINT_REAL_MATRIX  

	FUNCTION MAT_NORM2(N,WW) result(NORM)
		complex(dp) 	:: WW(N,N)
		real(dp) 	:: NORM
		integer 	:: N,i,j

		NORM=0
		WW=matmul(transpose(conjg(WW)),WW)
		do i=1,N
			NORM=NORM+WW(i,i)	
		end do

		NORM=sqrt(NORM)
	END FUNCTION MAT_NORM2


end module
