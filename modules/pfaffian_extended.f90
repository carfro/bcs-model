! FUNCTION overlap_pfaffian(lambda,NDIM,dimD,dimM,UD,UM,VD,VM)
!     ! computes the overlap of two wave functions with pfaffian formula
!
!     implicit none
!
!     integer, parameter :: qp = selected_real_kind(33, 4931)
!
!     integer, intent(in) :: NDIM, dimD, dimM, lambda
!     complex*16 :: W((dimM+dimD),(dimM+dimD))
!     complex*16, intent(in) :: UD(NDIM,dimD), VD(NDIM,dimD), UM(NDIM,dimM), VM(NDIM,dimM)
!     integer :: IPIV((dimM+dimD),2)
!     complex(kind=qp),intent(out) :: overlap_pfaffian
!
!     W(1:dimM,1:dimM) = matmul(transpose(VM),UM)
!     W(1:dimM,(dimM+1):(dimM+dimD)) = matmul(transpose(VM),conjg(VD))
!     W((dimM+1):(dimM+dimD),1:dimM) = -matmul(transpose(conjg(VD)),VM)
!     W((dimM+1):(dimM+dimD),(dimM+1):(dimM+dimD)) = matmul(transpose(conjg(UD)),conjg(VD))
!
!     W=lambda*W
!
!     call Zpfaffian_ext(W,(dimM+dimD),(dimM+dimD),IPIV,overlap_pfaffian)
!     return
! END FUNCTION overlap_pfaffian

! FUNCTION overlap_pfaffian(lambda,NDIM,dimD,dimM,UD,UM,VD,VM,DMAT)
!     ! computes the overlap of two wave functions with pfaffian formula
!
!     implicit none
!
!     integer, parameter :: qp = selected_real_kind(33, 4931)
!
!     integer, intent(in) :: NDIM, dimD, dimM, lambda
!     complex*16 :: W((dimM+dimD),(dimM+dimD))
!     complex*16, intent(in) :: DMAT(NDIM,NDIM), UD(NDIM,dimD), VD(NDIM,dimD), UM(NDIM,dimM), VM(NDIM,dimM)
!     integer :: IPIV((dimM+dimD),2)
!     complex(kind=qp) :: overlap_pfaffian
!
!     W(1:dimM,1:dimM) = matmul(transpose(VM),UM)
!     W(1:dimM,(dimM+1):(dimM+dimD)) = matmul(transpose(VM),matmul(DMAT,conjg(VD)))
!     W((dimM+1):(dimM+dimD),1:dimM) = -matmul(transpose(conjg(VD)),matmul(transpose(DMAT),VM))
!     W((dimM+1):(dimM+dimD),(dimM+1):(dimM+dimD)) = matmul(transpose(conjg(UD)),conjg(VD))
!
!     W=lambda*W
!
!     call Zpfaffian_ext(W,(dimM+dimD),(dimM+dimD),IPIV,overlap_pfaffian)
!
! END FUNCTION overlap_pfaffian

Subroutine ZPfaffian_EXT(SK,LDS,N,Ipiv,Pf)
    ! CALCULATES THE PFAFFIAN OF A SKEW-SYMMETRIC MATRIX SK USING EXTENDED PRECISION

    Implicit none

    integer, parameter :: qp = selected_real_kind(33, 4931)

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
