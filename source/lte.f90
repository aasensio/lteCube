program lteCube
use globalModule, only : config, lineList, atmosphere3D, atmosphere, PH, identity, packageAtmosphere, packageSizeAtmosphere, packageStokes, &
	packageSizeStokes, myrank, atomicE1, atomicAbundance, atomicAffinity, atomicName
use ioModule, only : readConfigFile, initTransition, initAtmospheres, read3DModels, check
use synthModule, only : generateAtomicZeemanComponents, synthAllRegionsAndPixels
use mpiModule, only : mpiBroadcastGeneral, Master2Slave_SendAtmosphere, SlaveFromMaster_GetAtmosphere, Slave2Master_SendStokes, MasterFromSlave_GetStokes, killSlave
use chemicalModule, only : initChemical
use netcdf
implicit none

	integer :: mpi_status, nprocs, ierr, i, kill_flag, slave, packagesize, j
	character(len=8) :: date
	character(len=10) :: time
	character(len=5) :: zone
	integer :: values(8), xPos, yPos, loop, stopFlag
	include 'mpif.h'

	call MPI_INIT(mpi_status)
	call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpi_status)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpi_status)

	write(*,FMT='(A,I3,A,I3)') 'Node ', myrank, '/', nprocs

! Number of processes
	config%nprocs = nprocs
	
	config%activeSlaves = config%nprocs-1
	
! Read configuration file. This needs to be done only by the master
	if (myrank == 0) then
		call readConfigFile(config, lineList)
	
! Read 3D models and set the common logTau500 axis
		call read3DModels(config, atmosphere3D)
		atmosphere%nDepths = config%nDepths
		allocate(atmosphere%lTau500(atmosphere%nDepths))
		atmosphere%lTau500 = atmosphere3D%lTau500		

		call initChemical(atomicE1, atomicAbundance, atomicAffinity, atomicName)
	endif
					
! Broadcast information to all nodes	
	call mpiBroadcastGeneral(myrank, config, lineList, atmosphere)
		
! Allocate memory for atmospheres in all nodes
	call initAtmospheres(config, atmosphere)

! lTau500 has been already allocated and set in the master. Allocate memory for the slaves
	if (myrank /= 0) then
		if (associated(atmosphere%lTau500)) deallocate(atmosphere%lTau500)
		allocate(atmosphere%lTau500(atmosphere%nDepths))
	endif

! Broadcast common logTau500 axis
	call MPI_Bcast(atmosphere%lTau500,atmosphere%nDepths,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		
! Precompute Zeeman components and allocate the memory needed for the synthesis
	do i = 1, config%nRegions
		do j = 1, lineList(i)%nLines
			call generateAtomicZeemanComponents(lineList(i)%transition(j))
			call initTransition(lineList(i)%transition(j), atmosphere%nDepths)
		enddo
		
		allocate(lineList(i)%opacity(7,atmosphere%nDepths,lineList(i)%nLambdaTotal))
		allocate(lineList(i)%source(atmosphere%nDepths,lineList(i)%nLambdaTotal))
		allocate(lineList(i)%boundary(4,lineList(i)%nLambdaTotal))
		allocate(lineList(i)%opacityContinuum(atmosphere%nDepths,lineList(i)%nLambdaTotal))
		allocate(lineList(i)%stokesOut(config%nPixelsChunk,4,lineList(i)%nLambdaTotal))
	enddo
		
! Generate identity matrix
	identity = 0.d0
	do i = 1, 4
		identity(i,i) = 1.d0
	enddo
	
! Sizes of the packages
	packageSizeAtmosphere = atmosphere%nDepths * config%nPixelsChunk * 6 * sizeof(PH) + (config%nPixelsChunk * 2 + 2) * sizeof(mpi_status)
	allocate(packageAtmosphere(packageSizeAtmosphere))
	
	packageSizeStokes = sizeof(mpi_status) + config%nPixelsChunk * 2 * sizeof(mpi_status) + sizeof(mpi_status)
	do i = 1, config%nRegions
		packageSizeStokes = packageSizeStokes + 4 * lineList(i)%nLambdaTotal * config%nPixelsChunk * sizeof(PH)
	enddo
	allocate(packageStokes(packageSizeStokes))
	
! Barrier to start computations
	call MPI_Barrier(MPI_COMM_WORLD, ierr)

! We start from a different pixel
	if (config%initialPixel == 0) then
		xPos = 1
		yPos = 1
	else
		yPos = config%initialPixel / config%nx
		xPos = config%initialPixel - yPos * config%nx + 1
	endif


! FOR PROFILING
	! if (myrank == 0) then
	! 	do i = 1, 288
	! 		atmosphere%TChunk(:,i) = atmosphere3D%T(1,:,1)
	! 		atmosphere%vmacChunk(:,i) = atmosphere3D%vmac(1,:,1)
	! 		atmosphere%thetaBChunk(:,i) = atmosphere3D%thetaB(1,:,1)
	! 		atmosphere%chiBChunk(:,i) = atmosphere3D%chiB(1,:,1)
	! 		atmosphere%BChunk(:,i) = atmosphere3D%B(1,:,1)
	! 		atmosphere%PeChunk(:,i) = atmosphere3D%Pe(1,:,1)
	! 	enddo

	! 	call synthAllRegionsAndPixels(config, atmosphere, lineList)
	! 	stop
	! endif

! Initial distribution of works
	if (nprocs >= 2) then
		if (myrank == 0) then
			do slave = 1, min(config%nprocs-1, config%remainingColumns)
				
				call Master2Slave_SendAtmosphere(config, atmosphere, xPos, yPos, slave)
								
			enddo
		endif
	else
		write(*,*) 'Master-slave strategy cannot be implemented in serial computations'
		stop
	endif
	
! Main loop. Now the MPI master-slave logic comes
	
	stopFlag = 0
	do while (stopFlag == 0)
! Master		
		if (myrank == 0) then
			call MasterFromSlave_GetStokes(slave, config, lineList, atmosphere)
			
			if (config%remainingColumns /= 0) then
				call Master2Slave_SendAtmosphere(config, atmosphere, xPos, yPos, slave)
			else
				call killSlave(slave)
				config%activeSlaves = config%activeSlaves - 1
				write(*,*) 'Slave ', slave, 'killed. Remaining: ', config%activeSlaves
			endif
			
			if (config%activeSlaves == 0) then
				stopFlag = 1
			endif
		else
			call SlaveFromMaster_GetAtmosphere(config, atmosphere, stopFlag)

			if (stopFlag == 0) then
! Do work		
 				call synthAllRegionsAndPixels(config, atmosphere, lineList)
				call Slave2Master_SendStokes(myrank, config, lineList, atmosphere)
			endif
		endif

	enddo
	
	if (myrank == 0) then
		call check( nf90_close(config%outputID) )
		write(*,*) 'Done'
	endif
	
! Finalize MPI
	call MPI_FINALIZE(mpi_status)
	
   
end program lteCube