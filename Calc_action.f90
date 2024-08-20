subroutine Calc_action(nmat,ndim,nsite,nbc,temperature,xmat,umat,action,mass,tHooft)

  implicit none

  integer nmat,ndim,nsite,nbc
  double complex xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex umat(1:nmat,1:nmat,1:nsite)
  double precision action,kinetic,potential,trx2
  double precision temperature,lattice_spacing,mass,tHooft
  double complex commutator(1:nmat,1:nmat),ux(1:nmat,1:nmat),&
       uxumx(1:nmat,1:nmat),udx(1:nmat,1:nmat)
  integer isite,isite_p1,isite_m1
  integer idim,jdim
  integer imat,jmat,kmat

  lattice_spacing=1d0/temperature/dble(nsite)

  !potential term
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
  potential=potential*0.5d0*dble(nmat)*lattice_spacing*tHooft
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
  potential=potential+mass*mass*0.5d0*trx2&
       *lattice_spacing*dble(nmat)





  !**********************
  !**** kinetic term ****
  !**********************
  kinetic=0d0
  do isite=1,nsite
     
     if(isite.LT.nsite)then
        isite_p1=isite+1
     else
        isite_p1=1
     end if

     if(isite.GT.1)then
        isite_m1=isite-1
     else
        isite_m1=nsite
     end if
     
     do idim=1,ndim
        !u(t)*x(t+a)
        ux=(0d0,0d0)
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 ux(imat,jmat)=ux(imat,jmat)&
                      +umat(imat,kmat,isite)*xmat(kmat,jmat,idim,isite_p1)
              end do
           end do
        end do
        !u^dag(t-a)*x(t-a)
        udx=(0d0,0d0)
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 udx(imat,jmat)=udx(imat,jmat)&
                      +dconjg(umat(kmat,imat,isite_m1))*xmat(kmat,jmat,idim,isite_m1)
              end do
           end do
        end do
        !-0.5*u(t)*x(t+a)*u^dagger(t) + 2*x(t) - 1.5*u^dagger(t-a)*x(t-a)*u(t-a)
        do imat=1,nmat
           do jmat=1,nmat
              uxumx(imat,jmat)=(2d0,0d0)*xmat(imat,jmat,idim,isite)
           end do
        end do
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 uxumx(imat,jmat)=uxumx(imat,jmat)&
                      -(0.5d0,0d0)*ux(imat,kmat)*dconjg(umat(jmat,kmat,isite))&
                      -(1.5d0,0d0)*udx(imat,kmat)*umat(kmat,jmat,isite_m1)
              end do
           end do
        end do
        
        do imat=1,nmat
           do jmat=1,nmat
              kinetic=kinetic+dble(uxumx(imat,jmat)*uxumx(jmat,imat))
           end do
        end do
        
     end do
  end do
  
  
  
  kinetic=kinetic*0.5d0*dble(nmat)/lattice_spacing
  action=kinetic+potential


  return

END subroutine Calc_action
