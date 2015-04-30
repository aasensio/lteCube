module globalModule
implicit none

	real(kind=8), parameter :: PK = 1.3806503d-16, UMA = 1.66053873d-24, PC = 2.99792458d10
	real(kind=8), parameter :: PH = 6.62606876d-27, PHK = PH / PK, PHC = 2.d0 * PH / PC**2
	real(kind=8), parameter :: PME = 9.10938188d-28, PE = 4.8032d-10, PI = 3.14159265359d0
	real(kind=8), parameter :: PHC2 = 2.d0 * PH * PC**2
	real(kind=8), parameter :: SQRTPI = 1.77245385091d0, EARTH_SUN_DISTANCE = 1.495979d13
	real(kind=8), parameter :: LARMOR = PE/(4.d0*PI*PME*PC), PMP = 1.67262158d-24
	real(kind=8), parameter :: NUCLEAR_MAGNETON = PE*PH/(2.d0*PI*PMP)
	real(kind=8), parameter :: BOHR_MAGNETON = PE*PH/(2.d0*PI*PME)
	real(kind=8), parameter :: PH2C3 = 2.d0 * PH**2 * PC**3, PERTURBATION = 0.005d0
	real(kind=8), parameter :: OPA = PI * PE**2 / (PME * PC), EV_ERG = 1.60217646d-12, GRAV = 2.7398d4
	real(kind=8), parameter :: BOHR_R = (PH/(2.d0*PI))**2/(PME*PE**2), CEV_ERG = 1.60217646d-12


	type atmosphereType		
		character(len=256) :: file
		integer :: nDepths, nNodes, modelsInChunk
		real(kind=8) :: abundance
		real(kind=8), dimension(:), pointer :: lTau500, height, T, nodes, niovern, ui, u1, u2, u3, n2overn1, n1overn0, microturbulence, vmac, B, thetaB, chiB
		real(kind=8), dimension(:), pointer :: opacity500, splitting, molecularAbundance, molecularPartition, molecularEquilibrium
 		real(kind=8), dimension(:), pointer :: Pe, Pg, PH, PHminus, PHplus, PH2, PH2plus, PTotal, nHtot, TOrdered, lTau500Ordered
 		real(kind=8), dimension(:,:), pointer :: lTau500Chunk, TChunk, vmacChunk, BChunk, thetaBChunk, chiBChunk, PeChunk
 		real(kind=8), dimension(:,:), pointer :: profile, zeeman_voigt, zeeman_faraday, coefficients
 		integer, dimension(:), pointer :: order, xPosChunk, yPosChunk
	end type atmosphereType
	
	type atmosphere3DType		
		character(len=256) :: file
		integer :: nDepths
		real(kind=8), dimension(:), pointer :: lTau500
		real(kind=8), dimension(:,:,:), pointer :: T, vmac, B, thetaB, chiB, Pe		
	end type atmosphere3DType

	type transitionType
		integer :: nLambda, element, ionization, nComponents(3)
		logical :: active, molecule
		real(kind=8) :: lambda0, gf, Elow, gu, Ju, gl, Jl, sigmaABO, alphaABO, lambdaLeft, lambdaRight, mass
		real(kind=8) :: frequency0, boundary
		real(kind=8), dimension(:), pointer :: lineOpacity, backOpacity, dopplerWidth, deltaNu, damping
		real(kind=8), dimension(:,:), pointer :: source, splitting, strength, opacityContinuum
		real(kind=8), dimension(:,:,:), pointer :: opacity
	end type transitionType
	
	type lineListType		
		integer :: nLines, nLambdaTotal
		type(transitionType), dimension(:), pointer :: transition
		real(kind=8), dimension(:,:), pointer :: source, opacityContinuum, boundary
		real(kind=8), dimension(:,:,:), pointer :: opacity, stokesOut
		real(kind=8), dimension(:), pointer :: lambda, frequency, contOpacityForInterpolation, lambdaForInterpolation
		real(kind=8), pointer :: valuesNetCDF(:,:)
		integer :: startNetCDF(4), countNetCDF(4)
	end type lineListType
	
	type configType
		integer :: nPixelsChunk, nx, ny, nDepths, nRegions, nprocs, verbose, remainingColumns, activeSlaves, zeemanSynthesis, chemicalEquilibriumOption
		integer :: outputID, nx_id, ny_id, nstokes_id, nRegions_id, initialPixel
		integer, pointer :: regionSize_id(:), regionLambda_id(:), regionStokes_id(:)
		character(len=100) :: PeFile, TFile, BFile, thetaBFile, chiBFile, vFile, out_file, lTau500File
	end type configType

	real(kind=8), parameter :: CHPartitionT(15) = (/1000.d0,1500.d0,2000.d0,2500.d0,3000.d0,3500.d0,4000.d0,4500.d0,5000.d0,5500.d0,6000.d0,6500.d0,7000.d0,7500.d0,8000.d0/)
	real(kind=8), parameter :: CHPartitionU(15) = (/135.80691,218.34335,313.03912,422.82003,549.45839,694.31549,858.62851,1043.5906,1250.3427,1479.9248,1733.2079,&
		2010.8152,2313.0385,2639.7540,2990.3400/)
	real(kind=8), parameter :: CHEquilibrium(9) = (/7.04641,-8.49562,-9.63745,-8.30710,-1.72875,-4.09372,0.94498,-1.01589,0.04693/)

	real(kind=8), dimension(92) :: atomicE1, atomicE2, atomicWeight, atomicAbundance, atomicAffinity
	character(len=2), dimension(92) :: atomicName

			
	type(atmosphereType) :: atmosphere	
	type(configType) :: config
	type(atmosphere3DType) :: atmosphere3D
	type(lineListType), pointer :: lineList(:)
	
	real(kind=8) :: identity(4,4)
	
	integer :: packageSizeAtmosphere, packageSizeStokes, myrank
	character, allocatable :: packageAtmosphere(:), packageStokes(:)		
	
end module globalModule