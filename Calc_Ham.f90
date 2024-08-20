subroutine Calc_Ham(nmat,ndim,nsite,nbc,temperature,xmat,umat,P_xmat,P_umat,ham,mass,tHooft)

  implicit none

  integer nmat,ndim,nsite,nbc
  double complex xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex umat(1:nmat,1:nmat,1:nsite)
  double complex P_umat(1:nmat,1:nmat,1:nsite)
  double precision ham
  double precision temperature,mass,tHooft

  integer isite
  integer idim
  integer imat,jmat

  call Calc_action(nmat,ndim,nsite,nbc,temperature,xmat,umat,ham,mass,tHooft)

  do isite=1,nsite
     do idim=1,ndim
        do imat=1,nmat
           do jmat=1,nmat
              ham=ham&
                   +dble(P_xmat(imat,jmat,idim,isite)&
                   *P_xmat(jmat,imat,idim,isite))*0.5d0
           end do
        end do
     end do
     do imat=1,nmat
        do jmat=1,nmat
           ham=ham+dble(P_umat(imat,jmat,isite)*P_umat(jmat,imat,isite))*0.5d0
        end do
     end do
  end do
  
  return
  
END subroutine Calc_Ham
