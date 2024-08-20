subroutine Calc_energy(nmat,ndim,nsite,temperature,xmat,umat,energy,mass,tHooft)

  implicit none

  integer nmat,ndim,nsite,nbc
  double complex xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex umat(1:nmat,1:nmat,1:nsite)
  double precision energy
  double precision temperature,mass,tHooft
  double precision   trx2,trf2
  double complex commutator(1:nmat,1:nmat),ux(1:nmat,1:nmat),&
       uxumx(1:nmat,1:nmat)
  integer isite,isite_p1
  integer idim,jdim
  integer imat,jmat,kmat
  !potential term
  trf2=0d0
  do isite=1,nsite
     do idim=1,ndim-1
        do jdim=idim+1,ndim
           commutator=(0d0,0d0)
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    commutator(imat,jmat)=commutator(imat,jmat)&
                         +xmat(imat,kmat,idim,isite)*xmat(kmat,jmat,jdim,isite)&
                         -xmat(imat,kmat,jdim,isite)*xmat(kmat,jmat,idim,isite)
                 end do
              end do
           end do
           do imat=1,nmat
              do jmat=1,nmat
                 trf2=trf2&
                      +dble(commutator(imat,jmat)*dconjg(commutator(imat,jmat)))
              end do
           end do

        end do
     end do
  end do
  energy=trf2*0.75d0/dble(nmat*nsite)*2d0*tHooft
  !******************************  
  !*** Plane wave deformation ***
  !******************************
  !*****************
  !*** mass term ***
  !*****************
  trx2=0d0
  do isite=1,nsite
     do imat=1,nmat
        do jmat=1,nmat
           do idim=1,ndim
              trx2=trx2&
                   +dble(xmat(imat,jmat,idim,isite)&
                   *dconjg(xmat(imat,jmat,idim,isite)))
           end do
        end do
     end do
  end do
  energy=energy+mass*mass*(trx2)/dble(nmat*nsite)

  return

END subroutine Calc_Energy
