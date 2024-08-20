!***********************************************************
!***********************************************************
!   Log of a matrix : P -> Log(P)
SUBROUTINE MATRIX_LOG(MatSize,MAT,MATLOG)

  implicit none 

  integer MatSize
  double complex MAT(1:MatSize,1:MatSize),MATLOG(1:MatSize,1:MatSize)

  integer i,j,k
  double complex eig_log(1:MatSize)
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
     !write(*,*)w(i)
     eig_log(i)=dcmplx(dlog(dble(w(i))))
     !write(*,*)eig_log(i) 
 end do
  MATLOG=(0d0,0d0)
  do i=1,MatSize
     do j=1,MatSize
        do k=1,MatSize
           MATLOG(i,j)=MATLOG(i,j)+VR(i,k)*eig_log(k)*dconjg(VR(j,k))
       end do
     end do
  end do

  return

END SUBROUTINE MATRIX_LOG

