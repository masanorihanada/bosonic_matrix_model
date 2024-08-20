! umat -> exp(i*P_umat*dtau_umat)*umat
! xmat -> xmat + P_xmat*dtau_xmat
! P_umat -> P_umat - delh_umat*dtau_umat
! P_xmat -> P_xmat - delh_xmat*dtau_xmat
! delh_xmat(imat,jmat)=dS/dxmat(jmat,imat)

SUBROUTINE Calc_DELH(nmat,ndim,nsite,nbc,temperature,xmat,umat,&
delh_xmat,delh_umat,mass,tHooft)

  implicit none

  integer nmat,ndim,nsite,nbc
  double complex xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex umat(1:nmat,1:nmat,1:nsite)
  double complex delh_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex delh_umat(1:nmat,1:nmat,1:nsite),delh_umat_temp(1:nmat,1:nmat,1:nsite)

  double precision temperature,lattice_spacing,mass,tHooft
  integer imat,jmat,kmat
  integer idim,jdim
  integer isite,isite_p1,isite_m1
  double complex commutator(1:nmat,1:nmat,1:ndim,1:ndim),ux(1:nmat,1:nmat),&
       udx(1:nmat,1:nmat),uxudx(1:nmat,1:nmat),xxx,&
       temp1(1:nmat,1:nmat),temp2(1:nmat,1:nmat),xududx(1:nmat,1:nmat),&
       uu(1:nmat,1:nmat),udud(1:nmat,1:nmat)
  
  double complex com12(1:nmat,1:nmat),com23(1:nmat,1:nmat),&
       &com31(1:nmat,1:nmat),temp,uxumx(1:nmat,1:nmat)

  lattice_spacing=1d0/temperature/dble(nsite)
  !***************************
  !*** calculate delh_umat ***
  !***************************
  delh_umat=(0d0,0d0)
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
        !u(t)^¥dagger*x(t)
        udx=(0d0,0d0)
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 ux(imat,jmat)=ux(imat,jmat)&
                      +umat(imat,kmat,isite)*xmat(kmat,jmat,idim,isite_p1)
                 udx(imat,jmat)=udx(imat,jmat)&
                      +dconjg(umat(kmat,imat,isite))*xmat(kmat,jmat,idim,isite)
              end do
           end do
        end do
        !u(t)*x(t+a)*u^dagger(t)*x(t)
        uxudx=(0d0,0d0)
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 uxudx(imat,jmat)=uxudx(imat,jmat)+ux(imat,kmat)*udx(kmat,jmat)
              end do
           end do
        end do
        
        do imat=1,nmat
           do jmat=1,nmat
              delh_umat(imat,jmat,isite)=&
                   delh_umat(imat,jmat,isite)&
                   -(4d0,0d0)*uxudx(imat,jmat)!+dconjg(uxudx(jmat,imat)))
           end do
        end do

        !temp1=x(t+a)*u^dag(t)
        !temp2=u^dag(t-a)*x(t-a)
        !xududx=x(t+a)*u^dag(t)*u^dag(t-a)*x(t-a)
        temp1=(0d0,0d0)
        temp2=(0d0,0d0)
        xududx=(0d0,0d0)
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 temp1(imat,jmat)=temp1(imat,jmat)&
                      +xmat(imat,kmat,idim,isite_p1)*dconjg(umat(jmat,kmat,isite))
                 temp2(imat,jmat)=temp2(imat,jmat)&
                      +dconjg(umat(kmat,imat,isite_m1))*xmat(kmat,jmat,idim,isite_m1)             
              end do
           end do
        end do
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 xududx(imat,jmat)=xududx(imat,jmat)+temp1(imat,kmat)*temp2(kmat,jmat)          
              end do
           end do
        end do
        temp1=(0d0,0d0)   
        temp2=(0d0,0d0)
        !temp1=u(t)*xududx
        !temp2=xududx*u(t-a)
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 temp1(imat,jmat)=temp1(imat,jmat)&
                      +umat(imat,kmat,isite)*xududx(kmat,jmat)
                 temp2(imat,jmat)=temp2(imat,jmat)&
                      +xududx(imat,kmat)*umat(kmat,jmat,isite_m1)
              end do
           end do
        end do
        
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 delh_umat(imat,jmat,isite)=&
                      delh_umat(imat,jmat,isite)&
                      +(0.75d0,0d0)*umat(imat,kmat,isite)*temp2(kmat,jmat)
                 delh_umat(imat,jmat,isite_m1)=&
                      delh_umat(imat,jmat,isite_m1)&
                      +(0.75d0,0d0)*umat(imat,kmat,isite_m1)*temp1(kmat,jmat)
           end do
        end do
     end do
        
     end do
  end do

  do isite=1,nsite
     do imat=1,nmat
        do jmat=1,nmat
           delh_umat_temp(imat,jmat,isite)=-(delh_umat(imat,jmat,isite)-dconjg(delh_umat(jmat,imat,isite)))
        end do
     end do
  end do
  delh_umat=delh_umat_temp*(0d0,-1d0)*dcmplx(NMAT)/dcmplx(lattice_spacing)


  !***************************
  !*** calculate delh_xmat ***
  !***************************

  delh_xmat=(0d0,0d0)
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
        !u(t-a)^¥dagger*x(t-a)
        udx=(0d0,0d0)
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 ux(imat,jmat)=ux(imat,jmat)&
                      +umat(imat,kmat,isite)*xmat(kmat,jmat,idim,isite_p1)
                 udx(imat,jmat)=udx(imat,jmat)&
                      +dconjg(umat(kmat,imat,isite_m1))*xmat(kmat,jmat,idim,isite_m1)
              end do
           end do
        end do
        
        do imat=1,nmat
           do jmat=1,nmat
              delh_xmat(imat,jmat,idim,isite)=&
                   delh_xmat(imat,jmat,idim,isite)+(6.5d0,0d0)*xmat(imat,jmat,idim,isite)
              do kmat=1,nmat
                 delh_xmat(imat,jmat,idim,isite)=&
                      delh_xmat(imat,jmat,idim,isite)&
                      -(4d0,0d0)*ux(imat,kmat)*dconjg(umat(jmat,kmat,isite))&
                      -(4d0,0d0)*udx(imat,kmat)*umat(kmat,jmat,isite_m1)
              end do
           end do
        end do

        !uu=u(t-a)*u(t), udud = (uu)^dagger
        uu=(0d0,0d0)
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 uu(imat,jmat)=uu(imat,jmat)+umat(imat,kmat,isite_m1)*umat(kmat,jmat,isite)
              end do
           end do
        end do
        do imat=1,nmat
           do jmat=1,nmat
              udud(imat,jmat)=dconjg(uu(jmat,imat))
           end do
        end do
        !temp1=uu*x(t+a)
        !temp2=udud*x(t-a)
        temp1=(0d0,0d0)
        temp2=(0d0,0d0)
        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 temp1(imat,jmat)=temp1(imat,jmat)+uu(imat,kmat)*xmat(kmat,jmat,idim,isite_p1)
                 temp2(imat,jmat)=temp2(imat,jmat)+udud(imat,kmat)*xmat(kmat,jmat,idim,isite_m1)
              end do
           end do
        end do

        do imat=1,nmat
           do jmat=1,nmat
              do kmat=1,nmat
                 delh_xmat(imat,jmat,idim,isite_m1)=&
                      delh_xmat(imat,jmat,idim,isite_m1)&
                      +(0.75d0,0d0)*temp1(imat,kmat)*udud(kmat,jmat)
                 delh_xmat(imat,jmat,idim,isite_p1)=&
                      delh_xmat(imat,jmat,idim,isite_p1)&
                      +(0.75d0,0d0)*temp2(imat,kmat)*uu(kmat,jmat)
              end do
           end do
        end do
        
     end do
  end do
  delh_xmat=delh_xmat*dcmplx(nmat)/dcmplx(lattice_spacing)


  !commutator term
  do isite=1,nsite 
     commutator=(0d0,0d0)
     do idim=1,ndim-1
        do jdim=idim+1,ndim
           do imat=1,nmat
              do jmat=1,nmat
                 do kmat=1,nmat
                    commutator(imat,jmat,idim,jdim)=&
                         commutator(imat,jmat,idim,jdim)&
                         +xmat(imat,kmat,idim,isite)*xmat(kmat,jmat,jdim,isite)&
                         -xmat(imat,kmat,jdim,isite)*xmat(kmat,jmat,idim,isite)
                 end do
              end do
           end do
        end do
     end do
     
     do idim=1,ndim
        do jdim=1,ndim
           do imat=1,nmat
              do jmat=1,nmat
                 xxx=(0d0,0d0)
                 if(jdim.GT.idim)then
                    do kmat=1,nmat
                       xxx=xxx&
                            +xmat(imat,kmat,jdim,isite)&
                            *commutator(kmat,jmat,idim,jdim)&
                            -commutator(imat,kmat,idim,jdim)&
                            *xmat(kmat,jmat,jdim,isite)
                    end do
                 else if(jdim.LT.idim)then
                    do kmat=1,nmat
                       xxx=xxx&
                            -xmat(imat,kmat,jdim,isite)&
                            *commutator(kmat,jmat,jdim,idim)&
                            +commutator(imat,kmat,jdim,idim)&
                            *xmat(kmat,jmat,jdim,isite)
                    end do
                 end if
                 
                 delh_xmat(imat,jmat,idim,isite)=&
                      delh_xmat(imat,jmat,idim,isite)&
                      -dcmplx(nmat)*dcmplx(lattice_spacing)*xxx*dcmplx(tHooft)
              end do
           end do
           
        end do
     end do
     
  end do

  !***************************
  !*** deriv. of mass term ***
  !***************************
  do isite=1,nsite
     do imat=1,nmat
        do jmat=1,nmat
           temp=dcmplx(mass*mass*lattice_spacing*dble(nmat))
           do idim=1,ndim
              delh_xmat(imat,jmat,idim,isite)=&
                   delh_xmat(imat,jmat,idim,isite)&
                   +xmat(imat,jmat,idim,isite)*temp
           end do
        end do
     end do
  end do

  
  return

END SUBROUTINE Calc_DELH
