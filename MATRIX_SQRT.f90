
!***********************************************************
!***********************************************************
!   Sqrt of a matrix : P -> sqrt(P)
SUBROUTINE MATRIX_SQRT(MatSize,MAT,MATSQRT)

  implicit none 

  integer MatSize
  double complex MAT(1:MatSize,1:MatSize),MATSQRT(1:MatSize,1:MatSize)

  integer i,j,k
  double complex eig_sqrt(1:MatSize)
  !lapack
  character jobvl,jobvr
  integer lda,ldvl,ldvr,info,lwork
  double complex w(1:MatSize),VL(1:MatSize),VR(1:MatSize,1:MATSIZE),&
       work(1:2*MatSize),rwork(1:2*MatSize)

  jobvl='N'
  jobvr='V'
  lda=MatSize
  ldvl=1
  ldvr=MatSize
  lwork=2*MatSize

  !MAT2=MAT

  call ZGEEV(jobvl,jobvr,MatSize,MAT,lda,w,vl,ldvl,vr,ldvr,work,lwork,rwork,info)
  ! in this case w(i) is real. 
  do i=1,MatSize
      eig_sqrt(i)=dcmplx(dsqrt(dble(w(i))))
      !write(*,*)eig_sqrt(i)
  end do
  MATSQRT=(0d0,0d0)
  do i=1,MatSize
     do j=1,MatSize
        do k=1,MatSize
           MATSQRT(i,j)=MATSQRT(i,j)+VR(i,k)*eig_sqrt(k)*dconjg(VR(j,k))
       end do
     end do
  end do

  return

END SUBROUTINE MATRIX_SQRT
