module ioModule
use globalModule, only : PC, PH, PK, PI, lineListType, configType, atmosphereType, atmosphere3DType, transitionType
use netcdf
implicit none
contains

!------------------------------------------------------------
! Check if the NetCDF action ended correctly
!------------------------------------------------------------
	subroutine check(status)
	integer, intent ( in) :: status

		if (status /= nf90_noerr) then
			print *, trim(nf90_strerror(status))
			stop
		endif
	end subroutine check

!---------------------------------------------------------
! Read the input file and set up the variables
!---------------------------------------------------------
	subroutine readConfigFile(config, lineList)
	type(configType) :: config
	type(lineListType), pointer :: lineList(:)
	integer :: nargs, j, i, k
	character(len=100) :: config_file
	real(kind=8) :: lambdaInitial, lambdaStep
	integer :: nLambda
	character(len=1) :: str
	
! Verify if a configuration file is passed as an argument
! If not, use the default 'config' file
		nargs = iargc()

		if (nargs == 0) then
			write(*,*) 'Config file not give. Trying config.input...'			
			config_file = 'conf.input'
		endif

		if (nargs == 1) then	      
			call getarg(1,config_file)
			write(*,FMT='(A,A)') 'Using config file : ', config_file
		endif

! Read the configuration file
		open(unit=12,file=config_file,action='read',status='old')
		
		read(12,*) config%verbose
		read(12,*) config%zeemanSynthesis
		read(12,*) config%lTau500File
		read(12,*) config%PeFile
		read(12,*) config%TFile
		read(12,*) config%BFile
		read(12,*) config%thetaBFile
		read(12,*) config%chiBFile
		read(12,*) config%vFile
		read(12,*) config%nx, config%ny, config%nDepths		
		read(12,*) config%out_file
		read(12,*) config%initialPixel
		
		read(12,*) config%nPixelsChunk
		read(12,*) config%nRegions
		
		config%remainingColumns = config%nx * config%ny - config%initialPixel
				
		allocate(lineList(config%nRegions))
		
! Read line information
		do j = 1, config%nRegions
			read(12,*) lineList(j)%nLines
						
			allocate(lineList(j)%transition(lineList(j)%nLines))
			read(12,*) nLambda
			read(12,*) lambdaInitial
			read(12,*) lambdaStep
			lineList(j)%nLambdaTotal = nLambda
			
			allocate(lineList(j)%lambda(nLambda))
			allocate(lineList(j)%frequency(nLambda))
			
			allocate(lineList(j)%valuesNetCDF(4,nLambda))
			
			do k = 1, nLambda
				lineList(j)%lambda(k) = lambdaInitial + lambdaStep * (k-1.d0)
			enddo
			
			lineList(j)%frequency = PC / (lineList(j)%lambda * 1.d-8)
				
			do i = 1, lineList(j)%nLines
				lineList(j)%transition(i)%nLambda = nLambda								
				
				read(12,*) lineList(j)%transition(i)%lambda0, lineList(j)%transition(i)%Elow, lineList(j)%transition(i)%gf, lineList(j)%transition(i)%element, &
					lineList(j)%transition(i)%ionization, lineList(j)%transition(i)%mass, lineList(j)%transition(i)%sigmaABO, lineList(j)%transition(i)%alphaABO, &
					lineList(j)%transition(i)%Jl, lineList(j)%transition(i)%gl, lineList(j)%transition(i)%Ju, lineList(j)%transition(i)%gu
				
				if (config%verbose == 1) then
					write(*,*) 'Region ', j, ' - Line ', i, ' - l0=', lineList(j)%transition(i)%lambda0,' A'
				endif
				
				lineList(j)%transition(i)%gf = 10.d0**lineList(j)%transition(i)%gf
				lineList(j)%transition(i)%Elow = lineList(j)%transition(i)%Elow * PH * PC         ! Transform to erg
				lineList(j)%transition(i)%frequency0 = PC / (lineList(j)%transition(i)%lambda0 * 1.d-8)								
			enddo
		enddo
				
		close(12)

		if (config%initialPixel == 0) then
		
! Open new output file
			call check( nf90_create(config%out_file, NF90_CLOBBER, config%outputID) )
			call check( nf90_def_dim(config%outputID, "nx", config%nx, config%nx_id) )
			call check( nf90_def_dim(config%outputID, "ny", config%ny, config%ny_id) )
			call check( nf90_def_dim(config%outputID, "nstokes", 4, config%nstokes_id) )
			call check( nf90_def_dim(config%outputID, "nregions", config%nRegions, config%nRegions_id) )
			
			allocate(config%regionSize_id(config%nRegions))
			allocate(config%regionLambda_id(config%nRegions))
			allocate(config%regionStokes_id(config%nRegions))
			
			do i = 1, config%nRegions
				write(str,FMT='(I1)') i
				call check( nf90_def_dim(config%outputID, "nLambdaRegion"//str, lineList(i)%nLambdaTotal, config%regionSize_id(i)) )		
	  			call check( nf90_def_var(config%outputID, "lambda"//str, NF90_DOUBLE, (/ config%regionSize_id(i) /), config%regionLambda_id(i)) )
	 			call check( nf90_def_var(config%outputID, "stokes"//str, NF90_FLOAT, (/ config%nx_id, config%ny_id, config%nstokes_id, config%regionSize_id(i) /), config%regionStokes_id(i)) ) 			
			enddo
					
	 		call check( nf90_enddef(config%outputID) )
	 		
! Save wavelength axis
	 		do i = 1, config%nRegions
				call check( nf90_put_var(config%outputID, config%regionLambda_id(i), lineList(i)%lambda) )
			enddo
		else
! Open existing output file
			call check( nf90_open(config%out_file, NF90_WRITE, config%outputID) )
			call check( nf90_inq_dimid(config%outputID, "nx", config%nx_id) )
			call check( nf90_inq_dimid(config%outputID, "ny", config%ny_id) )
			call check( nf90_inq_dimid(config%outputID, "nstokes", config%nstokes_id) )
			call check( nf90_inq_dimid(config%outputID, "nregions", config%nRegions_id) )

			allocate(config%regionSize_id(config%nRegions))
			allocate(config%regionLambda_id(config%nRegions))
			allocate(config%regionStokes_id(config%nRegions))

			do i = 1, config%nRegions
				write(str,FMT='(I1)') i
				call check( nf90_inq_dimid(config%outputID, "nLambdaRegion"//str, config%regionSize_id(i)) )
	  			call check( nf90_inq_varid(config%outputID, "lambda"//str, config%regionLambda_id(i)) )
	 			call check( nf90_inq_varid(config%outputID, "stokes"//str, config%regionStokes_id(i)) ) 			
			enddo

		endif
						
	end subroutine readConfigFile

!---------------------------------------------------------
! Initialize arrays for atmosphere
!---------------------------------------------------------
	subroutine initAtmospheres(config, atmosphere)
	type(configType) :: config
	type(atmosphereType) :: atmosphere
	integer :: i
   							
		if (associated(atmosphere%height)) deallocate(atmosphere%height)
		allocate(atmosphere%height(atmosphere%nDepths))
		
		if (associated(atmosphere%T)) deallocate(atmosphere%T)
		allocate(atmosphere%T(atmosphere%nDepths))
						
		if (associated(atmosphere%vmac)) deallocate(atmosphere%vmac)
		allocate(atmosphere%vmac(atmosphere%nDepths))
		
		if (associated(atmosphere%B)) deallocate(atmosphere%B)
		allocate(atmosphere%B(atmosphere%nDepths))
		
		if (associated(atmosphere%thetaB)) deallocate(atmosphere%thetaB)
		allocate(atmosphere%thetaB(atmosphere%nDepths))
		
		if (associated(atmosphere%chiB)) deallocate(atmosphere%chiB)
		allocate(atmosphere%chiB(atmosphere%nDepths))
				
		if (associated(atmosphere%niovern)) deallocate(atmosphere%niovern)
		allocate(atmosphere%niovern(atmosphere%nDepths))
		
		if (associated(atmosphere%ui)) deallocate(atmosphere%ui)
		allocate(atmosphere%ui(atmosphere%nDepths))
		
		if (associated(atmosphere%u1)) deallocate(atmosphere%u1)
		allocate(atmosphere%u1(atmosphere%nDepths))
		
		if (associated(atmosphere%u2)) deallocate(atmosphere%u2)
		allocate(atmosphere%u2(atmosphere%nDepths))
		
		if (associated(atmosphere%u3)) deallocate(atmosphere%u3)
		allocate(atmosphere%u3(atmosphere%nDepths))
		
		if (associated(atmosphere%n1overn0)) deallocate(atmosphere%n1overn0)
		allocate(atmosphere%n1overn0(atmosphere%nDepths))
		
		if (associated(atmosphere%n2overn1)) deallocate(atmosphere%n2overn1)
		allocate(atmosphere%n2overn1(atmosphere%nDepths))
		
		if (associated(atmosphere%Pe)) deallocate(atmosphere%Pe)
		allocate(atmosphere%Pe(atmosphere%nDepths))
		
		if (associated(atmosphere%Pg)) deallocate(atmosphere%Pg)
		allocate(atmosphere%Pg(atmosphere%nDepths))
		
		if (associated(atmosphere%PH)) deallocate(atmosphere%PH)
		allocate(atmosphere%PH(atmosphere%nDepths))
		
		if (associated(atmosphere%PHminus)) deallocate(atmosphere%PHminus)
		allocate(atmosphere%PHminus(atmosphere%nDepths))
		
		if (associated(atmosphere%PHplus)) deallocate(atmosphere%PHplus)
		allocate(atmosphere%PHplus(atmosphere%nDepths))
		
		if (associated(atmosphere%PH2)) deallocate(atmosphere%PH2)
		allocate(atmosphere%PH2(atmosphere%nDepths))
		
		if (associated(atmosphere%PH2plus)) deallocate(atmosphere%PH2plus)
		allocate(atmosphere%PH2plus(atmosphere%nDepths))
		
		if (associated(atmosphere%PTotal)) deallocate(atmosphere%PTotal)
		allocate(atmosphere%PTotal(atmosphere%nDepths))
		
		if (associated(atmosphere%nHtot)) deallocate(atmosphere%nHtot)
		allocate(atmosphere%nHtot(atmosphere%nDepths))
		
		if (associated(atmosphere%order)) deallocate(atmosphere%order)
		allocate(atmosphere%order(atmosphere%nDepths))		
		
		if (associated(atmosphere%opacity500)) deallocate(atmosphere%opacity500)
		allocate(atmosphere%opacity500(atmosphere%nDepths))
		
		if (associated(atmosphere%splitting)) deallocate(atmosphere%splitting)
		allocate(atmosphere%splitting(atmosphere%nDepths))
		
		if (associated(atmosphere%profile)) deallocate(atmosphere%profile)
		allocate(atmosphere%profile(2,atmosphere%nDepths))
		
		if (associated(atmosphere%zeeman_voigt)) deallocate(atmosphere%zeeman_voigt)
		allocate(atmosphere%zeeman_voigt(3,atmosphere%nDepths))
		
		if (associated(atmosphere%zeeman_faraday)) deallocate(atmosphere%zeeman_faraday)
		allocate(atmosphere%zeeman_faraday(3,atmosphere%nDepths))
		
		if (associated(atmosphere%coefficients)) deallocate(atmosphere%coefficients)
		allocate(atmosphere%coefficients(7,atmosphere%nDepths))
				
		if (associated(atmosphere%lTau500Chunk)) deallocate(atmosphere%lTau500Chunk)
		allocate(atmosphere%lTau500Chunk(atmosphere%nDepths,config%nPixelsChunk))
		
		if (associated(atmosphere%TChunk)) deallocate(atmosphere%TChunk)
		allocate(atmosphere%TChunk(atmosphere%nDepths,config%nPixelsChunk))
				
		if (associated(atmosphere%vmacChunk)) deallocate(atmosphere%vmacChunk)
		allocate(atmosphere%vmacChunk(atmosphere%nDepths,config%nPixelsChunk))
		
		if (associated(atmosphere%BChunk)) deallocate(atmosphere%BChunk)
		allocate(atmosphere%BChunk(atmosphere%nDepths,config%nPixelsChunk))
		
		if (associated(atmosphere%thetaBChunk)) deallocate(atmosphere%thetaBChunk)
		allocate(atmosphere%thetaBChunk(atmosphere%nDepths,config%nPixelsChunk))
		
		if (associated(atmosphere%chiBChunk)) deallocate(atmosphere%chiBChunk)
		allocate(atmosphere%chiBChunk(atmosphere%nDepths,config%nPixelsChunk))
						
		if (associated(atmosphere%PeChunk)) deallocate(atmosphere%PeChunk)
		allocate(atmosphere%PeChunk(atmosphere%nDepths,config%nPixelsChunk))
		
		if (associated(atmosphere%xPosChunk)) deallocate(atmosphere%xPosChunk)
		allocate(atmosphere%xPosChunk(config%nPixelsChunk))
		
		if (associated(atmosphere%yPosChunk)) deallocate(atmosphere%yPosChunk)
		allocate(atmosphere%yPosChunk(config%nPixelsChunk))

				
	end subroutine initAtmospheres
	
!---------------------------------------------------------
! Allocate some variables for the lines that will be used in the synthesis
!---------------------------------------------------------
	subroutine initTransition(transition, nDepths)
	type(transitionType) :: transition
	integer :: nDepths
	integer :: i
			
 		allocate(transition%lineOpacity(nDepths))
		allocate(transition%dopplerWidth(nDepths))
		allocate(transition%deltaNu(nDepths))
		allocate(transition%damping(nDepths))						
				
	end subroutine initTransition
	

!---------------------------------------------------------
! Allocate some variables for the lines that will be used in the synthesis
!---------------------------------------------------------
	subroutine read3DModels(config, atmosphere3D)
	type(configType) :: config
	type(atmosphere3DType) :: atmosphere3D
	integer :: nDepths
	integer :: i
	real(kind=4), allocatable :: dummy(:,:,:)
		
		if (config%verbose == 1) write(*,*) 'Reading logTau500 axis...'
		allocate(atmosphere3D%lTau500(config%nDepths))
		open(unit=12,file=config%lTau500File,action='read',access='stream',status='old')
		read(12) atmosphere3D%lTau500	
		close(12)
				
		allocate(dummy(config%nx, config%nDepths, config%ny))
	
		if (config%verbose == 1) write(*,*) 'Reading T cube...'
		allocate(atmosphere3D%T(config%nx, config%nDepths, config%ny))
		open(unit=12,file=config%TFile,action='read',access='stream',status='old')
		read(12) dummy
		atmosphere3D%T = dummy
		close(12)
		
		if (config%verbose == 1) write(*,*) 'Reading v cube...'
		allocate(atmosphere3D%vmac(config%nx, config%nDepths, config%ny))		
		open(unit=12,file=config%vFile,action='read',access='stream',status='old')
		read(12) dummy
		atmosphere3D%vmac = dummy
		close(12)
				
		if (config%verbose == 1) write(*,*) 'Reading B cube...'
		allocate(atmosphere3D%B(config%nx, config%nDepths, config%ny))		
		open(unit=12,file=config%BFile,action='read',access='stream',status='old')
		read(12) dummy
		atmosphere3D%B = dummy
		close(12)
		
		if (config%verbose == 1) write(*,*) 'Reading thetaB cube...'
		allocate(atmosphere3D%thetaB(config%nx, config%nDepths, config%ny))		
		open(unit=12,file=config%thetaBFile,action='read',access='stream',status='old')
		read(12) dummy
		atmosphere3D%thetaB = dummy
		close(12)
		atmosphere3D%thetaB = atmosphere3D%thetaB * PI / 180.0
		
		if (config%verbose == 1) write(*,*) 'Reading chiB cube...'
		allocate(atmosphere3D%chiB(config%nx, config%nDepths, config%ny))		
		open(unit=12,file=config%chiBFile,action='read',access='stream',status='old')
		read(12) dummy
		atmosphere3D%chiB = dummy
		close(12)
		atmosphere3D%chiB = atmosphere3D%chiB * PI / 180.0
		
		if (config%verbose == 1) write(*,*) 'Reading Ne cube...'
		allocate(atmosphere3D%Pe(config%nx, config%nDepths, config%ny))		
		open(unit=12,file=config%PeFile,action='read',access='stream',status='old')
		read(12) dummy
		atmosphere3D%Pe = dummy
		close(12)
		
! The file contains the electron density. Transform to pressure
		atmosphere3D%Pe = atmosphere3D%Pe * PK * atmosphere3D%T
		
		deallocate(dummy)
		
	end subroutine read3DModels

	
end module ioModule