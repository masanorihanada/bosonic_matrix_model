!##############################################################################
!######              BFSS matrix model on lattice                     #########
!######                                                               #########
!######                 written by Masanori Hanada                    #########
!######                                                               #########
!######            ver.1  bosonic, HMC, no parallelization.           #########
!##############################################################################
!Mersenne twister.
include 'mt19937.f90'
program BFSS

  use mtmod !Mersenne twistor

  implicit none

  include 'size.h'
  include 'include.h'

  double precision unitarity_breaking,norm
  double complex inner_prod
  !---------------------------------




  integer CheckHam,info



  open(unit=10,status='OLD',file='input_v1.dat',action='READ')
  read(10,*) input_config
  read(10,*) output_config
  read(10,*) data_output
  read(10,*) init
  read(10,*) temperature
  read(10,*) mass
  read(10,*) tHooft
  read(10,*) ntraj
  read(10,*) nskip
  read(10,*) ntau
  read(10,*) dtau_xmat
  read(10,*) dtau_umat
  read(10,*) checkham
  read(10,*) mersenne_seed

  close(10)


  !call random_seed() this function does not exist. Enrico.

  !*************************************
  !*** Set the initial configuration ***
  !*************************************

  if(init.EQ.0)then
     !continue from old config
     open(unit=9,status='OLD',file=input_config,action='READ')
     call mtgetu(9)
     read(9,*) itraj
     read(9,*) xmat
     read(9,*) umat
     close(9)

  else if(init.EQ.1)then
     !initialize random number generator
     call sgrnd(mersenne_seed)
     !new config, cold start
     itraj=1
     xmat=(0d0,0d0)
     umat=(0d0,0d0)
     do isite=1,nsite
        do idim=1,ndim
           do imat=1,nmat
              xmat(imat,imat,idim,isite)=(0d0,0d0)
           end do
        end do
     end do
     do isite=1,nsite
        do imat=1,nmat
           umat(imat,imat,isite)=(1d0,0d0)
        end do
     end do
  end if

  !**************************************************
  !**************************************************
  nacceptance=0 !number of acceptance
  ntrial=0 !number of trial

  !************************************
  !************************************
  !     Make the output file
  !************************************
  !************************************

  open(unit=10,status='REPLACE',file=data_output,action='WRITE')
  write(10,*) "#size of the gauge group: nmat=",nmat
  write(10,*) "#temperature, mass, 't Hooft coupling:",temperature, mass, tHooft
  write(10,*) "#ntau=",ntau
  write(10,*) "#dtau for xmat=",Dtau_xmat
  write(10,*) "#dtau for umat=",Dtau_umat
  write(10,*) "# traj, ham_fin-ham_init, energy, |Pol.|, trx2, trf2, unitarity_breaking, acceptance"
  write(10,*) "#------------------------------------------------"



  nacceptance=0
  ntrial=0
  do while (itraj.LE.ntraj)

     backup_xmat=xmat
     backup_umat=umat

     call Generate_P_xmat(nmat,ndim,nsite,P_xmat)
     call Generate_P_umat(nmat,nsite,P_umat)


     call Calc_Ham(nmat,ndim,nsite,nbc,temperature,xmat,umat,&
          P_xmat,P_umat,ham_init,mass,tHooft)


     call Molecular_Dynamics(nmat,ndim,nsite,nbc,temperature,&
          ntau,dtau_xmat,dtau_umat,xmat,umat,P_xmat,P_umat,mass,tHooft)


     call Calc_Ham(nmat,ndim,nsite,nbc,temperature,xmat,umat,&
          P_xmat,P_umat,ham_fin,mass,tHooft)


     !write(*,*)ham_init-ham_fin,ham_init,ham_fin

     !call random_number(metropolis)
     metropolis=grnd()
     ntrial=ntrial+1
     If(dexp(ham_init-ham_fin) > metropolis)THEN
        !accept
        nacceptance=nacceptance+1
     else
        !reject
        xmat=backup_xmat
        umat=backup_umat

     end If



     ! measurements
     if(MOD(itraj,nskip).EQ.0)then

        !projection to traceless, Hermitian & unitary
        do isite=1,nsite
           do idim=1,ndim
              trace=0d0
              do imat=1,nmat
                 trace=trace+dble(xmat(imat,imat,idim,isite))
              end do
              trace=trace/dble(nmat)
              do imat=1,nmat
                 xmat(imat,imat,idim,isite)=dcmplx(dble(xmat(imat,imat,idim,isite)))-dcmplx(trace)
              end do
           end do
        end do

        do isite=1,nsite
           do idim=1,ndim
              do imat=1,nmat-1
                 do jmat=imat+1,nmat
                    xmat(jmat,imat,idim,isite)=dconjg(xmat(imat,jmat,idim,isite))
                 end do
              end do
           end do
        end do

        do isite=1,nsite
           norm=0d0
           do imat=1,nmat
              norm=norm+dble(umat(imat,1,isite)*dconjg(umat(imat,1,isite)))
           end do
           norm=dsqrt(norm)
           do imat=1,nmat
             umat(imat,1,isite)=umat(imat,1,isite)/dcmplx(norm)
           end do
           do jmat=2,nmat
              do kmat=1,jmat-1
                 inner_prod=(0d0,0d0)
                 do imat=1,nmat
                    inner_prod=inner_prod+umat(imat,jmat,isite)*dconjg(umat(imat,kmat,isite))
                 end do
                 do imat=1,nmat
                    umat(imat,jmat,isite)=umat(imat,jmat,isite)-inner_prod*umat(imat,kmat,isite)
                 end do
              end do
              norm=0d0
              do imat=1,nmat
                 norm=norm+dble(umat(imat,jmat,isite)*dconjg(umat(imat,jmat,isite)))
              end do
              norm=dsqrt(norm)
              do imat=1,nmat
                 umat(imat,jmat,isite)=umat(imat,jmat,isite)/dcmplx(norm)
              end do
           end do
        end do
        
        call Calc_Polyakov(nmat,nsite,umat,Pol,unitarity_breaking)
        !call Calc_action(nmat,ndim,nsite,nbc,temperature,xmat,umat,action,flux)
        call Calc_F2(nmat,ndim,nsite,nbc,xmat,trf2)
        call Calc_X2(nmat,ndim,nsite,nbc,xmat,trx2)
        call Calc_energy(nmat,ndim,nsite,temperature,xmat,umat,energy,mass,tHooft)
       

        !output
        write(10,*)itraj,-ham_init+ham_fin,energy,abs(Pol),trx2,trf2,unitarity_breaking,dble(nacceptance)/dble(ntrial)
        write(*,*)itraj,-ham_init+ham_fin,energy,abs(Pol),trx2,trf2,unitarity_breaking,dble(nacceptance)/dble(ntrial)

     end if


     itraj=itraj+1
  end do
  !**************************************************
  !**************************************************
  !   End of iteration
  !**************************************************
  !**************************************************

  close(10)




  open(UNIT = 22, File = output_config, STATUS = "REPLACE", ACTION = "WRITE")
  call mtsaveu(22)
  write(22,*) itraj
  write(22,*) xmat
  write(22,*) umat
  close(22)






end program BFSS

include 'BoxMuller.f90'
include 'MATRIX_iEXP.f90'
include 'MatrixInverse.f90'
include 'MATRIX_LOG.f90'
include 'MATRIX_SQRT.f90'
include 'Calc_Polyakov.f90'
include 'Molecular_Dynamics.f90'
include 'Calc_action.f90'
include 'Calc_Ham.f90'
include 'Generate_P_xmat.f90'
include 'Generate_P_umat.f90'
include 'Calc_DELH.f90'
include 'Calc_F2.f90'
include 'Calc_X2.f90'
include 'Calc_energy.f90'
