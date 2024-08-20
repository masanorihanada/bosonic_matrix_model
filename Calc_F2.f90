subroutine Calc_F2(nmat,ndim,nsite,nbc,xmat,trf2)

  implicit none

  integer nmat,ndim,nsite,nbc
  double complex xmat(1:nmat,1:nmat,1:ndim,1:nsite),commutator(1:nmat,1:nmat)
  double precision potential,trf2

  integer isite
  integer idim,jdim
  integer imat,jmat,kmat

  potential=0d0
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
                 potential=potential&
                      +dble(commutator(imat,jmat)*dconjg(commutator(imat,jmat)))
              end do
           end do
           
        end do
     end do
  end do
  !potential=potential*0.5d0*dble(nmat)*lattice_spacing
  trf2=potential/dble(nmat*nsite)*2d0


  return

END subroutine Calc_F2
