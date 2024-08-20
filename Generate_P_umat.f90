! Generate P_umat with Gaussian weight.
SUBROUTINE Generate_P_umat(nmat,nsite,P_umat)
  
  implicit none
  integer nmat,nsite
  integer imat,jmat,idim,isite
  double precision r1,r2

  double complex P_umat(1:nmat,1:nmat,1:nsite)
  
  do isite=1,nsite
     do imat=1,nmat-1
        do jmat=imat+1,nmat
           call BoxMuller(r1,r2)
           P_umat(imat,jmat,isite)=&
                dcmplx(r1/dsqrt(2d0))+dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
           P_umat(jmat,imat,isite)=&
                dcmplx(r1/dsqrt(2d0))-dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
        end do
     end do
     do imat=1,nmat
        call BoxMuller(r1,r2)
        P_umat(imat,imat,isite)=dcmplx(r1)
     end do
  end do
  
  return
  
END SUBROUTINE Generate_P_umat
