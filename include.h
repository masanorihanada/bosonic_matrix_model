
  integer nbc !boundary condition for fermions; 0 -> pbc, 1 -> apbc
  integer init !initial condition; 0 -> continue, 1 -> new

!matrices
  double complex xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex umat(1:nmat,1:nmat,1:nsite)
  double complex backup_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex backup_umat(1:nmat,1:nmat,1:nsite)
  double complex P_xmat(1:nmat,1:nmat,1:ndim,1:nsite)
  double complex P_umat(1:nmat,1:nmat,1:nsite)
  !matrix indices
  integer imat,jmat,kmat,lmat
  integer isite,jsite,ksite,lsite
  integer idim,jdim,kdim,ldim
  integer ispin,jspin,kspin,lspin
  
  double precision temperature,mass,tHooft

  !parameters for molecular evolution
  integer ntau
  double precision dtau_xmat,dtau_umat

  !number of trajectories
  integer ntraj     !total number of trajectories at the end of the run
  integer itraj

  double precision ham_init,ham_fin,metropolis

  !measurements
  integer nskip !measurement is performed every nskiptrajectories
  integer nacceptance !number of acceptance
  integer ntrial !number of trial

  !absolute value of Polyakov loop
  double precision Pol

  !energy, action, etc
    double precision energy,action, trx2, trf2,trace
character(150) input_config,data_output,output_config

!For Mersenne Twister
integer mersenne_seed
