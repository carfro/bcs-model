module mroutines 
	use MO_module, only: dp, qp
	implicit none
	private 
	public det, det_eig, inv
contains

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

end module
