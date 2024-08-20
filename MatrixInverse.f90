!***********************************************************
!***********************************************************
! output -> (MAT)^{-1}
SUBROUTINE MatrixInverse(NMAT,MAT)

  implicit none

  integer NMAT,DIM,lwork,ipiv(1:NMAT),info 
  doublecomplex MAT(1:NMAT,1:NMAT)
  !work(lwork);lwork must be equal to or larger than 2*(NMAT).    
  doublecomplex WORK(16*(NMAT)) 
  lwork=16*(NMAT)

  call zgetrf(NMAT,NMAT,MAT,NMAT,ipiv,info)
!  write(*,*)"info=",info
  call zgetri(NMAT,MAT,NMAT,ipiv,work,lwork,info)
!  write(*,*)"info=",info
  return

END SUBROUTINE MatrixInverse

