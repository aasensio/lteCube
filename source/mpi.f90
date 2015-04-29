module mpiModule
use globalModule, only : configType, packageAtmosphere, packageSizeAtmosphere, lineListType, atmosphereType, atmosphere3D, packageStokes, packageSizeStokes,&
	atomicE1, atomicE2, atomicWeight, atomicAbundance, atomicAffinity, atomicName
use ioModule, only : check
use variables
use netcdf
implicit none

contains

!------------------------------------------------------------
! Broadcast some general information
!------------------------------------------------------------
	subroutine mpiBroadcastGeneral(myrank, config, lineList, atmosphere)
	type(configType) :: config
	type(lineListType), pointer :: lineList(:)
	type(atmosphereType) :: atmosphere
	integer :: ierr, myrank, i, j
	include 'mpif.h'
		
		call MPI_Bcast(config%nRegions,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_Bcast(config%zeemanSynthesis,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		
		if (myrank /= 0) then
			allocate(lineList(config%nRegions))
		endif
		
		do j = 1, config%nRegions
			call MPI_Bcast(lineList(j)%nLines,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
			
			if (myrank /= 0) then
				allocate(lineList(j)%transition(lineList(j)%nLines))
			endif
			
			call MPI_Bcast(lineList(j)%nLambdaTotal,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
			
			if (myrank /= 0) then
				allocate(lineList(j)%lambda(lineList(j)%nLambdaTotal))
				allocate(lineList(j)%frequency(lineList(j)%nLambdaTotal))
				allocate(lineList(j)%lambdaForInterpolation(lineList(j)%nLines))
				allocate(lineList(j)%contOpacityForInterpolation(lineList(j)%nLines))
			endif
			
			call MPI_Bcast(lineList(j)%lambda,lineList(j)%nLambdaTotal,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			call MPI_Bcast(lineList(j)%frequency,lineList(j)%nLambdaTotal,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
						
			do i = 1, lineList(j)%nLines
				
				call MPI_Bcast(lineList(j)%transition(i)%lambda0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%frequency0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%Elow,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%gf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%element,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%ionization,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%mass,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%sigmaABO,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%alphaABO,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%Jl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%Ju,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%gl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%gu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
				call MPI_Bcast(lineList(j)%transition(i)%molecule,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)			
				
			enddo
		enddo
		
		call MPI_Bcast(atmosphere%nDepths,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)									
				
		call MPI_Bcast(config%nPixelsChunk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Broadcast information for chemical equilibrium
		call MPI_Bcast(equil,273,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

		do i = 1, 21
			call MPI_Bcast(elements(i),2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
		enddo

		call MPI_Bcast(abund_atom,21,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_Bcast(pot_ion,21,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_Bcast(afinidad,21,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

		call MPI_Bcast(estequio,273*21,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_Bcast(charge,273,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

		do i = 1, 273
			call MPI_Bcast(nombre_mol(i),16,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
		enddo

		call MPI_Bcast(composicion,273*4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
		call MPI_Bcast(n_atoms_mol,273,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

		call MPI_Bcast(n_included,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

		if (myrank /= 0) then
			allocate(which_included(n_included))
		endif

		call MPI_Bcast(which_included,n_included,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

		call MPI_Bcast(atomicE1,92,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_Bcast(atomicE2,92,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_Bcast(atomicWeight,92,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_Bcast(atomicAbundance,92,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_Bcast(atomicAffinity,92,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		do i = 1, 92
			call MPI_Bcast(atomicName(i),2,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
		enddo
		
 	end subroutine mpiBroadcastGeneral
 	
!------------------------------------------------------------
! MASTER: send a calculation to a slave
!------------------------------------------------------------
	subroutine Master2Slave_SendAtmosphere(config, atmosphere, xPos, yPos, slave)
	type(configType) :: config
	type(atmosphereType) :: atmosphere
	integer :: i, xPos, yPos, slave, nPacked, stopFlag
	integer :: ierr, pos
	character(len=8) :: date
	character(len=10) :: time
	character(len=5) :: zone
	integer :: values(8)
	include 'mpif.h'
		
		nPacked = 0		
		do i = 1, min(config%nPixelsChunk, config%remainingColumns)
			if (xPos > config%nx) then 
				xPos = 1
				yPos = yPos + 1
			endif
			
			atmosphere%xPosChunk(i) = xPos
			atmosphere%yPosChunk(i) = yPos
			
			atmosphere%TChunk(:,i) = atmosphere3D%T(xPos,:,yPos)
			atmosphere%vmacChunk(:,i) = atmosphere3D%vmac(xPos,:,yPos)
			atmosphere%thetaBChunk(:,i) = atmosphere3D%thetaB(xPos,:,yPos)
			atmosphere%chiBChunk(:,i) = atmosphere3D%chiB(xPos,:,yPos)
			atmosphere%BChunk(:,i) = atmosphere3D%B(xPos,:,yPos)
			atmosphere%PeChunk(:,i) = atmosphere3D%Pe(xPos,:,yPos)
						
			xPos = xPos + 1
			nPacked = nPacked + 1						
		enddo
		
		config%remainingColumns = config%remainingColumns - nPacked
						
		pos = 0
		call MPI_Pack(slave, 1, MPI_INTEGER, packageAtmosphere, packageSizeAtmosphere, pos, MPI_COMM_WORLD,ierr)
		call MPI_Pack(nPacked, 1, MPI_INTEGER, packageAtmosphere, packageSizeAtmosphere, pos, MPI_COMM_WORLD,ierr)
		call MPI_Pack(atmosphere%xPosChunk, config%nPixelsChunk, MPI_INTEGER, packageAtmosphere, packageSizeAtmosphere, pos, MPI_COMM_WORLD,ierr)
		call MPI_Pack(atmosphere%yPosChunk, config%nPixelsChunk, MPI_INTEGER, packageAtmosphere, packageSizeAtmosphere, pos, MPI_COMM_WORLD,ierr)
		call MPI_Pack(atmosphere%TChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, packageAtmosphere, packageSizeAtmosphere, pos, MPI_COMM_WORLD,ierr)
		call MPI_Pack(atmosphere%vmacChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, packageAtmosphere, packageSizeAtmosphere, pos, MPI_COMM_WORLD,ierr)
		call MPI_Pack(atmosphere%thetaBChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, packageAtmosphere, packageSizeAtmosphere, pos, MPI_COMM_WORLD,ierr)
		call MPI_Pack(atmosphere%chiBChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, packageAtmosphere, packageSizeAtmosphere, pos, MPI_COMM_WORLD,ierr)
		call MPI_Pack(atmosphere%BChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, packageAtmosphere, packageSizeAtmosphere, pos, MPI_COMM_WORLD,ierr)
		call MPI_Pack(atmosphere%PeChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, packageAtmosphere, packageSizeAtmosphere, pos, MPI_COMM_WORLD,ierr)
		
		call MPI_Send(packageAtmosphere, packageSizeAtmosphere, MPI_PACKED, slave, 10, MPI_COMM_WORLD, ierr)
						
		call date_and_time(date, time, zone, values)		
		write(*,FMT='(A,A,A,I5,A,I7,A,I5,A,F5.1,A)') '(',time,') * Master -> Slave ', slave, ' - sent ', nPacked,' columns - remaining ', config%remainingColumns, ' - ',&
			100.d0 * config%remainingColumns / (config%nx * config%ny), '%'
			
	end subroutine Master2Slave_SendAtmosphere
	
!------------------------------------------------------------
! SLAVE: get a calculation from the master
!------------------------------------------------------------
	subroutine SlaveFromMaster_GetAtmosphere(config, atmosphere, stopFlag)
	type(configType) :: config
	type(atmosphereType) :: atmosphere
	integer :: i, slave, nPacked, stopFlag
	integer :: ierr, pos
	character(len=8) :: date
	character(len=10) :: time
	character(len=5) :: zone
	integer :: values(8)
	include 'mpif.h'
	integer :: status(MPI_STATUS_SIZE)
		
			
		call MPI_Recv(packageAtmosphere, packageSizeAtmosphere, MPI_PACKED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
		
		if (status(MPI_TAG) == 10) then
			pos = 0
			call MPI_Unpack(packageAtmosphere, packageSizeAtmosphere, pos, slave, 1, MPI_INTEGER, MPI_COMM_WORLD,ierr)
			call MPI_Unpack(packageAtmosphere, packageSizeAtmosphere, pos, atmosphere%modelsInChunk, 1, MPI_INTEGER, MPI_COMM_WORLD,ierr)
			call MPI_Unpack(packageAtmosphere, packageSizeAtmosphere, pos, atmosphere%xPosChunk, config%nPixelsChunk, MPI_INTEGER, MPI_COMM_WORLD,ierr)
			call MPI_Unpack(packageAtmosphere, packageSizeAtmosphere, pos, atmosphere%yPosChunk, config%nPixelsChunk, MPI_INTEGER, MPI_COMM_WORLD,ierr)
			call MPI_Unpack(packageAtmosphere, packageSizeAtmosphere, pos, atmosphere%TChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,ierr)
			call MPI_Unpack(packageAtmosphere, packageSizeAtmosphere, pos, atmosphere%vmacChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,ierr)
			call MPI_Unpack(packageAtmosphere, packageSizeAtmosphere, pos, atmosphere%thetaBChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,ierr)
			call MPI_Unpack(packageAtmosphere, packageSizeAtmosphere, pos, atmosphere%chiBChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,ierr)
			call MPI_Unpack(packageAtmosphere, packageSizeAtmosphere, pos, atmosphere%BChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,ierr)
			call MPI_Unpack(packageAtmosphere, packageSizeAtmosphere, pos, atmosphere%PeChunk, config%nPixelsChunk*atmosphere%nDepths, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,ierr)
		endif
				
		if (status(MPI_TAG) == 999) then
			stopFlag = 1
		endif
									
	end subroutine SlaveFromMaster_GetAtmosphere
	
!------------------------------------------------------------
! SLAVE: send Stokes profiles to master
!------------------------------------------------------------
	subroutine Slave2Master_SendStokes(slave, config, lineList, atmosphere)
	type(configType) :: config
	type(lineListType), pointer :: lineList(:)
	type(atmosphereType) :: atmosphere
	integer :: slave, ierr, pos, i
	include 'mpif.h'
		
		pos = 0
		call MPI_Pack(slave, 1, MPI_INTEGER, packageStokes, packageSizeStokes, pos, MPI_COMM_WORLD,ierr)
		
		call MPI_Pack(atmosphere%modelsInChunk, 1, MPI_INTEGER, packageStokes, packageSizeStokes, pos, MPI_COMM_WORLD,ierr)
		call MPI_Pack(atmosphere%xPosChunk, config%nPixelsChunk, MPI_INTEGER, packageStokes, packageSizeStokes, pos, MPI_COMM_WORLD,ierr)
		call MPI_Pack(atmosphere%yPosChunk, config%nPixelsChunk, MPI_INTEGER, packageStokes, packageSizeStokes, pos, MPI_COMM_WORLD,ierr)		
		
		do i = 1, config%nRegions
			call MPI_Pack(lineList(i)%stokesOut, 4*lineList(i)%nLambdaTotal*config%nPixelsChunk, MPI_DOUBLE_PRECISION, packageStokes, packageSizeStokes, pos, MPI_COMM_WORLD,ierr)
		enddo
		
		call MPI_Send(packageStokes, packageSizeStokes, MPI_PACKED, 0, 3, MPI_COMM_WORLD, ierr)
							
	end subroutine Slave2Master_SendStokes
	
!------------------------------------------------------------
! MASTER: get Stokes profiles from slave
!------------------------------------------------------------
	subroutine MasterFromSlave_GetStokes(slave, config, lineList, atmosphere)
	type(configType) :: config
	type(lineListType), pointer :: lineList(:)
	type(atmosphereType) :: atmosphere
	integer :: slave, ierr, pos, i, j
	include 'mpif.h'
	integer :: status(MPI_STATUS_SIZE)
	
		call MPI_Recv(packageStokes, packageSizeStokes, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)		
		
		pos = 0
		call MPI_Unpack(packageStokes, packageSizeStokes, pos, slave, 1, MPI_INTEGER, MPI_COMM_WORLD,ierr)
	
		call MPI_Unpack(packageStokes, packageSizeStokes, pos, atmosphere%modelsInChunk, 1, MPI_INTEGER, MPI_COMM_WORLD,ierr)
		call MPI_Unpack(packageStokes, packageSizeStokes, pos, atmosphere%xPosChunk, config%nPixelsChunk, MPI_INTEGER, MPI_COMM_WORLD,ierr)
		call MPI_Unpack(packageStokes, packageSizeStokes, pos, atmosphere%yPosChunk, config%nPixelsChunk, MPI_INTEGER, MPI_COMM_WORLD,ierr)
		
		do i = 1, config%nRegions
			call MPI_Unpack(packageStokes, packageSizeStokes, pos, lineList(i)%stokesOut, 4*lineList(i)%nLambdaTotal*config%nPixelsChunk, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,ierr)			
		enddo
		
! Write data
		do i = 1, config%nRegions
			do j = 1, atmosphere%modelsInChunk
				lineList(i)%valuesNetCDF = lineList(i)%stokesOut(j,:,:)
				lineList(i)%startNetCDF = (/ atmosphere%xPosChunk(j), atmosphere%yPosChunk(j), 1, 1 /)
				lineList(i)%countNetCDF = (/ 1, 1, 4, lineList(i)%nLambdaTotal /)
 				call check( nf90_put_var(config%outputID, config%regionStokes_id(i), lineList(i)%valuesNetCDF, start=lineList(i)%startNetCDF, count=lineList(i)%countNetCDF) )
			enddo
		enddo
	
	end subroutine MasterFromSlave_GetStokes
	
!------------------------------------------------------------
! Send the kill signal to a slave
!------------------------------------------------------------
	subroutine killSlave(slave)
	integer :: slave
	integer :: ierr, i
	include 'mpif.h'

! Send an empty message with the tag 999
		call MPI_Send(0, 0, MPI_INTEGER, slave, 999, MPI_COMM_WORLD, ierr)

 	end subroutine killSlave

end module mpiModule