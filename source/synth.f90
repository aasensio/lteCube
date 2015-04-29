module synthModule
use globalModule, only : OPA, PH, PC, PK, PHK, UMA, EV_ERG, SQRTPI, LARMOR, atmosphereType, lineListType, transitionType, &
	configType, myrank, CHEquilibrium, CHPartitionT, CHPartitionU
use atomicPartitionModule, only : partitionAtomic
use mathsModule, only : saha, vecfvoigt2, calculateDamping, planckFrequency, shortCharacteristics, computeHeight, gasPressure, vecfvoigt2, &
	vecfvoigt_zeeman, vecfvoigt_zeeman2, strength_zeeman, formal_sol_polarized, linInterpol
use backgroundOpacityModule, only : backgroundOpacity
use chemicalModule, only : getAbundance
implicit none
contains

!-----------------------------------------------------------------
! Fill the strength and splitting variables of the atomic_transition structure
!-----------------------------------------------------------------
	subroutine generateAtomicZeemanComponents(transition)
	type(transitionType) :: transition
	integer :: i, n, nlow, nup, iup, ilow, i_pi, i_red, i_blue, cual
	real(kind=8) :: Mup, Mlow, strength

		nup = 2*transition%Ju+1
		nlow = 2*transition%Jl+1

		i_pi = 0
		i_blue = 0
		i_red = 0

! First count the number of components of each type (
		transition%nComponents = 0
		do iup = 1, nup
			Mup = transition%Ju + 1 - iup  ! Mup=J...-J
			do ilow = 1, 3
				Mlow = Mup-2+ilow			! Mlow=Mup-1,Mup,Mup+1 (the allowed transitions)
				if (abs(Mlow) <= transition%Jl) then
					transition%nComponents(ilow) = transition%nComponents(ilow)+1
				endif
			enddo
		enddo

		allocate(transition%splitting(3,maxval(transition%nComponents)))
		allocate(transition%strength(3,maxval(transition%nComponents)))

! Now generate all the data
		do iup = 1, nup
			Mup = transition%Ju + 1 - iup  ! Mup=J...-J
			do ilow = 1, 3
				Mlow = Mup-2+ilow			! Mlow=Mup-1,Mup,Mup+1 (the allowed transitions)
				if (abs(Mlow) <= transition%Jl) then
					if (ilow == 1) then
						i_blue = i_blue + 1
						cual = i_blue
					endif
					if (ilow == 2) then
						i_pi = i_pi + 1
						cual = i_pi
					endif
					if (ilow == 3) then
						i_red = i_red + 1
						cual = i_red
					endif

					transition%strength(ilow,cual) = strength_zeeman(transition%Ju,transition%Jl,Mup,Mlow)
					transition%splitting(ilow,cual) = (transition%gu*Mup - transition%gl*Mlow)

				endif
			enddo
		enddo

	end subroutine generateAtomicZeemanComponents

!-----------------------------------------------------------------
! Return the seven independent elements of the absorption matrix
! Remember that, zeeman_voigt(q,:) and zeeman_faraday(q,:) have
!  q=1  Mlow=Mup-1  (sigma blue)
!  q=2  Mlow=Mup    (sigma pi)
!  q=3  Mlow=Mup+1  (sigma red)
!-----------------------------------------------------------------
	subroutine zeeman_opacity(transition, atmos, frequency)	
	real(kind=8) :: frequency
	type(atmosphereType) :: atmos
	type(transitionType) :: transition	
	integer :: n, nlow, nup, iup, ilow, i_pi, i_blue, i_red, cual
	real(kind=8) :: Mup, Mlow, strength
	
		nup = 2*transition%Ju+1
		nlow = 2*transition%Jl+1

		atmos%zeeman_voigt = 0.d0
		atmos%zeeman_faraday = 0.d0

		i_red = 0
		i_pi = 0
		i_blue = 0

		do iup = 1, nup
			Mup = transition%Ju + 1 - iup  ! Mup=J...-J
			do ilow = 1, 3
				Mlow = Mup-2+ilow			! Mlow=Mup-1,Mup,Mup+1 (the allowed transitions)
				if (abs(Mlow) <= transition%Jl) then

					if (ilow == 1) then
						i_blue = i_blue + 1
						cual = i_blue
					endif
					if (ilow == 2) then
						i_pi = i_pi + 1
						cual = i_pi
					endif
					if (ilow == 3) then
						i_red = i_red + 1
						cual = i_red
					endif

					atmos%splitting = LARMOR * atmos%B * transition%splitting(ilow,cual)

 					atmos%profile = vecfvoigt_zeeman2(transition%damping, (transition%frequency0 - frequency - transition%frequency0 * atmos%vmac / PC + atmos%splitting) / transition%deltaNu)
 					
					atmos%zeeman_voigt(ilow,:) = atmos%zeeman_voigt(ilow,:) + transition%strength(ilow,cual) * atmos%profile(1,:)
					atmos%zeeman_faraday(ilow,:) = atmos%zeeman_faraday(ilow,:) + transition%strength(ilow,cual) * atmos%profile(2,:)
															
				endif
			enddo
		enddo
		
! Classical absorption coefficients
		atmos%coefficients(1,:) = 0.5d0 * (atmos%zeeman_voigt(2,:)*sin(atmos%thetaB)**2 + &
			0.5d0*(atmos%zeeman_voigt(1,:)+atmos%zeeman_voigt(3,:))*(1.d0+cos(atmos%thetaB)**2))  ! eta_I
 		atmos%coefficients(2,:) = 0.5d0 * (atmos%zeeman_voigt(2,:) - &
 			0.5d0*(atmos%zeeman_voigt(1,:)+atmos%zeeman_voigt(3,:))) * sin(atmos%thetaB)**2*cos(2.d0*atmos%chiB)  ! eta_Q
 		atmos%coefficients(3,:) = 0.5d0 * (atmos%zeeman_voigt(2,:) - &
 			0.5d0*(atmos%zeeman_voigt(1,:)+atmos%zeeman_voigt(3,:))) * sin(atmos%thetaB)**2*sin(2.d0*atmos%chiB)  ! eta_U
 		atmos%coefficients(4,:) = 0.5d0 * (atmos%zeeman_voigt(3,:)-atmos%zeeman_voigt(1,:)) * cos(atmos%thetaB)         ! eta_V
! 
! ! Magneto-optical coefficients
 		atmos%coefficients(5,:) = 0.5d0 * (atmos%zeeman_faraday(2,:) - &
 			0.5d0*(atmos%zeeman_faraday(1,:)+atmos%zeeman_faraday(3,:))) * sin(atmos%thetaB)**2*cos(2.d0*atmos%chiB)  ! rho_Q
 		atmos%coefficients(6,:) = 0.5d0 * (atmos%zeeman_faraday(2,:) - &
 			0.5d0*(atmos%zeeman_faraday(1,:)+atmos%zeeman_faraday(3,:))) * sin(atmos%thetaB)**2*sin(2.d0*atmos%chiB)  ! rho_U
 		atmos%coefficients(7,:) = 0.5d0 * (atmos%zeeman_faraday(3,:)-atmos%zeeman_faraday(1,:)) * cos(atmos%thetaB)  ! rho_V 		 		 		

	end subroutine zeeman_opacity
	
!------------------------------------------------
! Synthesize all lines in a region for Stokes I, Q, U and V
!------------------------------------------------
	subroutine synthLinesZeeman(atmosphere, lineList, chunkIndex)
	type(atmosphereType) :: atmosphere
	type(lineListType) :: lineList
	real(kind=8) :: n1overn0, n2overn1, ei1, ei2, weight, abundance
	real(kind=8), allocatable :: voigtProfile(:)
	integer :: i, j, loop, fromLambda, toLambda, k, chunkIndex
					 		 		 				
		fromLambda = 1
		loop = 1
		
		lineList%opacity = 0.d0
		lineList%opacityContinuum = 0.d0
		lineList%boundary = 0.d0
		lineList%source = 0.d0
		
! Evaluate continuum opacity at the central wavelengths of the lines and then use linear interpolation
		do j = 1, atmosphere%nDepths
			do k = 1, lineList%nLines
				lineList%lambdaForInterpolation(k) = lineList%transition(k)%lambda0
				lineList%contOpacityForInterpolation(k) = backgroundOpacity(atmosphere%T(j), atmosphere%Pe(j), atmosphere%PH(j), atmosphere%PHminus(j), &
					atmosphere%PHplus(j), atmosphere%PH2(j), atmosphere%PH2plus(j), lineList%transition(k)%lambda0)
			enddo
			call linInterpol(lineList%lambdaForInterpolation,lineList%contOpacityForInterpolation,lineList%lambda,lineList%opacityContinuum(j,:))
		enddo
		
		do i = 1, lineList%nLines
			call partitionAtomic(atmosphere%T, lineList%transition(i)%element, atmosphere%u1, atmosphere%u2, &
				atmosphere%u3, ei1, ei2, weight, atmosphere%abundance)
			
			atmosphere%n1overn0 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u1, atmosphere%u2, ei1)
			atmosphere%n2overn1 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u2, atmosphere%u3, ei2)
			atmosphere%niovern = 1.d0 / (1.d0 + atmosphere%n1overn0 + atmosphere%n2overn1 * atmosphere%n1overn0)						
						
			select case (lineList%transition(i)%ionization)
				case(1)
					atmosphere%ui = atmosphere%u1
				case(2)
					atmosphere%ui = atmosphere%u2
				case(3)
					atmosphere%ui = atmosphere%u3
			end select
		
! Compute line opacity
			lineList%transition(i)%lineOpacity = OPA * lineList%transition(i)%gf / atmosphere%ui * dexp(-lineList%transition(i)%Elow / (PK * atmosphere%T)) *&
				(1.d0 - dexp(-PHK * lineList%transition(i)%frequency0 / atmosphere%T)) * atmosphere%niovern * &
				(atmosphere%nhtot * atmosphere%abundance)
																							
! Doppler width			
			lineList%transition(i)%dopplerWidth = dsqrt(2.d0 * PK * atmosphere%T / (lineList%transition(i)%mass * UMA))
			lineList%transition(i)%deltaNu = lineList%transition(i)%dopplerWidth * lineList%transition(i)%frequency0 / PC
									
! Compute the damping
			lineList%transition(i)%damping = calculateDamping(atmosphere%T, atmosphere%nHtot, atmosphere%Pe, atmosphere%PTotal, &
				lineList%transition(i)%dopplerWidth, ei1, lineList%transition(i)%Elow / EV_ERG, 0.d0, 0.d0, 0.d0, lineList%transition(i)%lambda0, &
				lineList%transition(i)%alphaABO, lineList%transition(i)%sigmaABO, lineList%transition(i)%mass)
				
! Add to total opacity in the region
			do j = 1, lineList%nLambdaTotal
				call zeeman_opacity(lineList%transition(i), atmosphere, lineList%frequency(j))
				do k = 1, 7
 					lineList%opacity(k,:,j) = lineList%opacity(k,:,j) + lineList%transition(i)%lineOpacity * atmosphere%coefficients(k,:) / (lineList%transition(i)%deltaNu * SQRTPI)
				enddo
				lineList%source(:,j) = planckFrequency(lineList%frequency(j), atmosphere%T)				
				lineList%boundary(1,j) = lineList%source(1,j)
			enddo
			
		enddo
		
		
! Add continuum opacity
		lineList%opacity(1,:,:) = lineList%opacity(1,:,:) + lineList%opacityContinuum		
										
! Solve RT equation
		do i = 1, lineList%nLambdaTotal			
			lineList%stokesOut(chunkIndex,:,i) = formal_sol_polarized(atmosphere%height, lineList%opacity(:,:,i), lineList%source(:,i), 1.d0, lineList%boundary(:,i), 1)
		enddo
							
	end subroutine synthLinesZeeman

!------------------------------------------------
! Compute molecular abundance
!------------------------------------------------
! 	subroutine molecularAbundance(molCode, atmosphere)
! 	integer :: molCode
! 	type(atmosphereType) :: atmosphere
! 	real(kind=8) ::  ei1, ei2, weight, abundance
		
! ! CH
! 		if (molCode == 101) then
			
! ! Hydrogen
! 			call partitionAtomic(atmosphere%T, 1, atmosphere%u1, atmosphere%u2, atmosphere%u3, ei1, ei2, weight, abundance)			
! 			atmosphere%n1overn0 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u1, atmosphere%u2, ei1)
! 			atmosphere%n2overn1 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u2, atmosphere%u3, ei2)
! 			atmosphere%niovern = atmosphere%nhtot * abundance / (1.d0 + atmosphere%n1overn0 + atmosphere%n2overn1 * atmosphere%n1overn0)

! ! Carbon
! 			call partitionAtomic(atmosphere%T, 6, atmosphere%u1, atmosphere%u2, atmosphere%u3, ei1, ei2, weight, abundance)
! 			atmosphere%n1overn0 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u1, atmosphere%u2, ei1)
! 			atmosphere%n2overn1 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u2, atmosphere%u3, ei2)
! 			atmosphere%niovern = atmosphere%niovern * atmosphere%nhtot * abundance / (1.d0 + atmosphere%n1overn0 + atmosphere%n2overn1 * atmosphere%n1overn0)	
		
! 			atmosphere%niovern = atmosphere%niovern * PK * atmosphere%T / atmosphere%equilibriumConstant

! 		endif

! 		print *, atmosphere%niovern / atmosphere%nhtot
! 		stop

! 	end subroutine molecularAbundance

!------------------------------------------------
! Synthesize all lines in a region for Stokes I
!------------------------------------------------
	subroutine synthLines(atmosphere, lineList, chunkIndex)
	type(atmosphereType) :: atmosphere
	type(lineListType) :: lineList
	real(kind=8) :: n1overn0, n2overn1, ei1, ei2, weight, abundance
	real(kind=8), allocatable :: voigtProfile(:)
	integer :: i, j, loop, fromLambda, toLambda, k, chunkIndex
					 		 		 				
		fromLambda = 1
		loop = 1
		
		lineList%opacity = 0.d0
		lineList%opacityContinuum = 0.d0
		lineList%boundary = 0.d0
		lineList%source = 0.d0
		
! Evaluate continuum opacity at the central wavelengths of the lines and then use linear interpolation
		do j = 1, atmosphere%nDepths
			do k = 1, lineList%nLines
				lineList%lambdaForInterpolation(k) = lineList%transition(k)%lambda0
				lineList%contOpacityForInterpolation(k) = backgroundOpacity(atmosphere%T(j), atmosphere%Pe(j), atmosphere%PH(j), atmosphere%PHminus(j), &
					atmosphere%PHplus(j), atmosphere%PH2(j), atmosphere%PH2plus(j), lineList%transition(k)%lambda0)
			enddo
			call linInterpol(lineList%lambdaForInterpolation,lineList%contOpacityForInterpolation,lineList%lambda,lineList%opacityContinuum(j,:))
		enddo
		
		do i = 1, lineList%nLines

! Molecular line
			if (lineList%transition(i)%molecule) then

! Compute line opacity
				lineList%transition(i)%lineOpacity = OPA * lineList%transition(i)%gf / atmosphere%molecularPartition * dexp(-lineList%transition(i)%Elow / (PK * atmosphere%T)) *&
					(1.d0 - dexp(-PHK * lineList%transition(i)%frequency0 / atmosphere%T)) * atmosphere%molecularAbundance


			else
! Atomic line
				call partitionAtomic(atmosphere%T, lineList%transition(i)%element, atmosphere%u1, atmosphere%u2, &
					atmosphere%u3, ei1, ei2, weight, atmosphere%abundance)
				
				atmosphere%n1overn0 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u1, atmosphere%u2, ei1)
				atmosphere%n2overn1 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u2, atmosphere%u3, ei2)
				atmosphere%niovern = 1.d0 / (1.d0 + atmosphere%n1overn0 + atmosphere%n2overn1 * atmosphere%n1overn0)						
							
				select case (lineList%transition(i)%ionization)
					case(1)
						atmosphere%ui = atmosphere%u1
					case(2)
						atmosphere%ui = atmosphere%u2
					case(3)
						atmosphere%ui = atmosphere%u3
				end select
! Compute line opacity
				lineList%transition(i)%lineOpacity = OPA * lineList%transition(i)%gf / atmosphere%ui * dexp(-lineList%transition(i)%Elow / (PK * atmosphere%T)) *&
					(1.d0 - dexp(-PHK * lineList%transition(i)%frequency0 / atmosphere%T)) * atmosphere%niovern * &
					(atmosphere%nhtot * atmosphere%abundance)

			endif		
																							
! Doppler width			
			lineList%transition(i)%dopplerWidth = dsqrt(2.d0 * PK * atmosphere%T / (lineList%transition(i)%mass * UMA))
			lineList%transition(i)%deltaNu = lineList%transition(i)%dopplerWidth * lineList%transition(i)%frequency0 / PC
									
! Compute the damping
			lineList%transition(i)%damping = calculateDamping(atmosphere%T, atmosphere%nHtot, atmosphere%Pe, atmosphere%PTotal, &
				lineList%transition(i)%dopplerWidth, ei1, lineList%transition(i)%Elow / EV_ERG, 0.d0, 0.d0, 0.d0, lineList%transition(i)%lambda0, &
				lineList%transition(i)%alphaABO, lineList%transition(i)%sigmaABO, lineList%transition(i)%mass)
				
! Add to total opacity in the region
			do j = 1, lineList%nLambdaTotal
				atmosphere%coefficients(1,:) = vecfvoigt2(lineList%transition(i)%damping, (lineList%transition(i)%frequency0 - lineList%frequency(j) - &
					lineList%transition(i)%frequency0 * atmosphere%vmac / PC) / lineList%transition(i)%deltaNu)
				lineList%opacity(1,:,j) = lineList%opacity(1,:,j) + lineList%transition(i)%lineOpacity * atmosphere%coefficients(1,:) / (lineList%transition(i)%deltaNu * SQRTPI)
				lineList%source(:,j) = planckFrequency(lineList%frequency(j), atmosphere%T)				
				lineList%boundary(1,j) = lineList%source(1,j)
			enddo
			
		enddo
			
! Add continuum opacity
		lineList%opacity(1,:,:) = lineList%opacity(1,:,:) + lineList%opacityContinuum		
										
! Solve RT equation
		lineList%stokesOut(chunkIndex,1,:) = shortCharacteristics(atmosphere%height, lineList%opacity(1,:,:), lineList%source, 1.d0, lineList%boundary(1,:), 1)		
							
	end subroutine synthLines
	
!------------------------------------------------
! Synthesize all regions
!------------------------------------------------
	subroutine synthAllRegionsAndPixels(config, atmosphere, lineList)
	type(configType) :: config
	type(atmosphereType) :: atmosphere
	type(lineListType) :: lineList(:)
	integer :: i, j, k
		
		! atmosphere%modelsInChunk = 10
		do i = 1, atmosphere%modelsInChunk
			atmosphere%T = atmosphere%TChunk(:,i)
			atmosphere%vmac = atmosphere%vmacChunk(:,i)
			atmosphere%B = atmosphere%BChunk(:,i)
			atmosphere%thetaB = atmosphere%thetaBChunk(:,i)
			atmosphere%chiB = atmosphere%chiBChunk(:,i)
			atmosphere%Pe = atmosphere%PeChunk(:,i)
			
! Compute the chemical equilibrium
			call getAbundance(9, atmosphere%nDepths, atmosphere%T, atmosphere%Pe, atmosphere%PH, atmosphere%PHminus, atmosphere%PHplus, &
					atmosphere%PH2, atmosphere%PH2plus, atmosphere%PTotal, atmosphere%molecularAbundance)

			! call gasPressure(atmosphere%Pe, atmosphere%T, atmosphere%PH, atmosphere%PHminus, atmosphere%PHplus, atmosphere%PH2, atmosphere%PH2plus, atmosphere%PTotal)

! Total hydrogen density
			atmosphere%nhtot = (atmosphere%PH + atmosphere%PHminus + atmosphere%PHplus + atmosphere%PH2 + atmosphere%PH2plus) / (PK * atmosphere%T)

! Compute partition function
			 call linInterpol(CHPartitionT, CHPartitionU, atmosphere%T, atmosphere%molecularPartition)
				
! Generate height axis
			do k = 1, atmosphere%nDepths
				atmosphere%opacity500(k) = backgroundOpacity(atmosphere%T(k), atmosphere%Pe(k), atmosphere%PH(k), atmosphere%PHminus(k), &
					atmosphere%PHplus(k), atmosphere%PH2(k), atmosphere%PH2plus(k), 5000.d0)		
			enddo
	
! Compute the height depth scale			
			call computeHeight(atmosphere%lTau500, atmosphere%opacity500, atmosphere%height)

! CH equilibrium constant
			! atmosphere%equilibriumConstant = 0.d0
			! do k = 0, 8
				! atmosphere%equilibriumConstant = atmosphere%equilibriumConstant + CHEquilibrium(k+1) * (dlog10(5040.d0 / atmosphere%T))**i
			! enddo

! Transform the units from SI to cgs multipliying by 10 the necessary times depending on the
! units of the equilibrium constant
			! atmosphere%equilibriumConstant = 10.d0**atmosphere%equilibriumConstant * 10.d0
	
 			do j = 1, config%nRegions
 				if (config%zeemanSynthesis == 1) then
 					call synthLinesZeeman(atmosphere, lineList(j), i)
 				else
 					call synthLines(atmosphere, lineList(j), i)
 				endif
 			enddo
		enddo   	
		
	end subroutine synthAllRegionsAndPixels
	
end module synthModule