subroutine Calc_X2(nmat,ndim,nsite,nbc,xmat,trx2)

  implicit none

  integer nmat,ndim,nsite,nbc
  double complex xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double precision trx2


  integer isite
  integer idim,jdim
  integer imat,jmat

  trx2=0d0
  do isite=1,nsite
     do idim=1,ndim
        do imat=1,nmat
           do jmat=1,nmat
              trx2=trx2+dble(xmat(imat,jmat,idim,isite)&
                   &*dconjg(xmat(imat,jmat,idim,isite)))
           end do
        end do
     end do
  end do
  trx2=trx2/dble(nmat*nsite)
  
  return

END subroutine Calc_X2
