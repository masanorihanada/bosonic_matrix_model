!***********************************************************
!***********************************************************
!    Polaakov loop 
SUBROUTINE Calc_Polyakov(nmat,nsite,umat,Pol,unitarity_breaking)
  
  implicit none
  
  integer nmat,nsite
  double complex umat(1:nmat,1:nmat,1:nsite)
  double complex trace,MAT1(1:nmat,1:nmat),MAT2(1:nmat,1:nmat)
  double precision Pol,unitarity_breaking
  integer imat,jmat,kmat
  integer isite

  trace=(0d0,0d0)
 
  do imat=1,nmat
     do jmat=1,nmat
        MAT1(imat,jmat)=umat(imat,jmat,1)
     end do
  end do
  do isite=2,nsite
     MAT2=(0d0,0d0)
     do imat=1,nmat
        do jmat=1,nmat
           do kmat=1,nmat
              MAT2(imat,jmat)=MAT2(imat,jmat)&
                   +MAT1(imat,kmat)*umat(kmat,jmat,isite)
           end do
        end do
     end do
     MAT1=MAT2
  end do
  trace=(0d0,0d0)
  do imat=1,nmat
     trace=trace+MAT2(imat,imat)
  end do
  
 
  Pol=abs(trace)/dble(NMAT)
  
  
unitarity_breaking=0d0
do imat=1,nmat
do jmat=1,nmat
unitarity_breaking=unitarity_breaking+dble(mat2(imat,jmat)*dconjg(mat2(imat,jmat)))
end do
end do
unitarity_breaking=unitarity_breaking-dble(nmat)

  return
  
END SUBROUTINE Calc_Polyakov
