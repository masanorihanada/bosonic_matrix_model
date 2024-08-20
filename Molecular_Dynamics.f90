! umat -> exp(i*P_umat*dtau_umat)*umat
! xmat -> xmat + P_xmat*dtau_xmat
! P_umat -> P_umat - delh_umat*dtau_umat
! P_xmat -> P_xmat - delh_xmat*dtau_xmat
! delh_xmat(imat,jmat)=dS/dxmat(jmat,imat)

subroutine Molecular_Dynamics(nmat,ndim,nsite,nbc,temperature,&
     ntau,dtau_xmat,dtau_umat,xmat,umat,P_xmat,P_umat,mass,tHooft)

  implicit none

  integer nbc,nmat,nsite,ndim

  double complex xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex umat(1:nmat,1:nmat,1:nsite)
  double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex P_umat(1:nmat,1:nmat,1:nsite)
  double complex delh_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex delh_umat(1:nmat,1:nmat,1:nsite)

  integer ntau
  double precision dtau_xmat,dtau_umat

  doubleprecision temperature,mass,tHooft
  double complex MAT(1:NMAT,1:NMAT),EXP_iP(1:NMAT,1:NMAT)

  !matrix indices
  integer imat,jmat,kmat,lmat
  integer isite,jsite,ksite,lsite
  integer idim,jdim,kdim,ldim
  integer ispin,jspin,kspin,lspin
  integer step

  ! first step of leap frog
  xmat=xmat+P_xmat*dcmplx(0.5d0*dtau_xmat)
  do isite=1,nsite  
     do imat=1,nmat
        do jmat=1,nmat
           MAT(imat,jmat)=P_umat(imat,jmat,isite)&
                *dcmplx(0.5d0*dtau_umat)
        end do
     end do
     call MATRIX_iEXP(NMAT,MAT,EXP_iP)
     do imat=1,nmat
        do jmat=1,nmat
           MAT(imat,jmat)= umat(imat,jmat,isite)
           umat(imat,jmat,isite)=(0d0,0d0)
        end do
     end do
     do imat=1,nmat
        do jmat=1,nmat
           do kmat=1,nmat
              umat(imat,jmat,isite)=umat(imat,jmat,isite)&
                   +EXP_iP(imat,kmat)*MAT(kmat,jmat) 
           end do
        end do
     end do
  end do
  
  ! second,...,Ntau-th   
  do step=2,ntau
     call Calc_DELH(nmat,ndim,nsite,nbc,temperature,xmat,umat,delh_xmat,delh_umat,mass,tHooft)
     P_umat=P_umat-delh_umat*dcmplx(dtau_umat)
     P_xmat=P_xmat-delh_xmat*dcmplx(dtau_xmat)
     
     xmat=xmat+P_xmat*dcmplx(dtau_xmat)
     do isite=1,nsite  
        do imat=1,nmat
           do jmat=1,nmat
              MAT(imat,jmat)=P_umat(imat,jmat,isite)&
                   *dcmplx(dtau_umat)
           end do
        end do
        call MATRIX_iEXP(NMAT,MAT,EXP_iP)
        do imat=1,nmat
           do jmat=1,nmat
              MAT(imat,jmat)= umat(imat,jmat,isite)
              umat(imat,jmat,isite)=(0d0,0d0)
           end do
        end do
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 umat(imat,jmat,isite)=umat(imat,jmat,isite)&
                      +EXP_iP(imat,kmat)*MAT(kmat,jmat) 
              end do
           end do
        end do
     end do
  end do
  ! last step
  call Calc_DELH(nmat,ndim,nsite,nbc,temperature,xmat,umat,delh_xmat,delh_umat,mass,tHooft)
  P_umat=P_umat-delh_umat*dcmplx(dtau_umat)
  P_xmat=P_xmat-delh_xmat*dcmplx(dtau_xmat)
  
  xmat=xmat+P_xmat*dcmplx(0.5d0*dtau_xmat)
  do isite=1,nsite  
     do imat=1,nmat
        do jmat=1,nmat
           MAT(imat,jmat)=P_umat(imat,jmat,isite)&
                *dcmplx(0.5d0*dtau_umat)
        end do
     end do
     call MATRIX_iEXP(NMAT,MAT,EXP_iP)
     do imat=1,nmat
        do jmat=1,nmat
           MAT(imat,jmat)= umat(imat,jmat,isite)
           umat(imat,jmat,isite)=(0d0,0d0)
        end do
     end do
     do imat=1,nmat
        do jmat=1,nmat
           do kmat=1,nmat
              umat(imat,jmat,isite)=umat(imat,jmat,isite)&
                   +EXP_iP(imat,kmat)*MAT(kmat,jmat) 
           end do
        end do
     end do
  end do

  return

END subroutine Molecular_Dynamics
