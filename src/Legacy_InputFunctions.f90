!==========================================================================================
    subroutine ReadInput_Umbrella(fileUnit)
       use MiscelaniousVars
       use SimpleDistPair, only: nDistPair, pairArrayIndx, CalcDistPairs_New
       use SimParameters, only: NMAX, NMIN, NPART, NPart_New, nMolTypes, maxMol, echoInput, NTotal, NTotal_New
       use Q6Functions, only: q6ArrayIndx, CalcQ6_Disp, CalcQ6_SwapIn, CalcQ6_SwapOut
       use ParallelVar, only: nout
       use WHAM_Module, only: refBin, refSizeNumbers
       implicit none
       integer, intent(in) :: fileUnit
       integer :: iUmbrella, AllocateStatus
       integer :: indxVar, stat
       integer :: iDisp, iSwapIn, iSwapOut
       real(dp) :: binSize, valMax, valMin
       real(dp), allocatable :: refVals(:)
       character(len=100), allocatable :: inputLines(:)
       character(len=100) :: refLine
       character(len=30) :: labelField
       character(len=30) :: umbrellaName

       read (fileUnit, *) labelField, nBiasVariables
       if (nBiasVariables .lt. 0) then
          write (*, *) "ERROR! The user has specified an invalid number of Umbrella Sampling Variables"
          write (*, *) labelField, nBiasVariables
          stop
       end if
       if (nBiasVariables .eq. 0) then
          useUmbrella = .false.
          write (nout, *) "Umbrella Sampling Used? = ", useUmbrella
          return
       else
          useUmbrella = .true.
          write (nout, *) "Umbrella Sampling Used? = ", useUmbrella
          write (nout, *) "Number of Umbrella Variables:", nBiasVariables
       end if

       read (fileUnit, "(A)") refLine
       allocate (inputLines(1:nBiasVariables))
       do iUmbrella = 1, nBiasVariables
          read (fileUnit, "(A)") inputLines(iUmbrella)
          if (echoinput) then
             write (35, "(A)") inputLines(iUmbrella)
          end if
       end do

       allocate (refVals(1:nBiasVariables))

       nDispFunc = 0
       nSwapInFunc = 0
       nSwapOutFunc = 0
       do iUmbrella = 1, nBiasVariables
          read (inputLines(iUmbrella), *) umbrellaName
          select case (trim(adjustl(umbrellaName)))
          case ("clustersize")
!          nSwapInFunc = nSwapInFunc + 1
!          nSwapOutFunc = nSwapOutFunc + 1
          case ("totalclustersize")
!          nSwapInFunc = nSwapInFunc + 1
!          nSwapOutFunc = nSwapOutFunc + 1
          case ("pairdist")
             nDispFunc = nDispFunc + 1
             nSwapInFunc = nSwapInFunc + 1
             nSwapOutFunc = nSwapOutFunc + 1
          case ("q6")
             nDispFunc = nDispFunc + 1
             nSwapInFunc = nSwapInFunc + 1
             nSwapOutFunc = nSwapOutFunc + 1
          case default
             write (*, *) "ERROR! Invalid variable type specified in input file"
             write (*, *) umbrellaName
             stop
          end select
       end do

       call AllocateUmbrellaVariables

       iDisp = 0
       iSwapIn = 0
       iSwapOut = 0
       do iUmbrella = 1, nBiasVariables
          read (inputLines(iUmbrella), *) umbrellaName
          select case (trim(adjustl(umbrellaName)))
          case ("clustersize")
             read (inputLines(iUmbrella), *) labelField, indxVar
             if (indxVar .le. 0) then
                write (*, *) "Error! An invalid molecule type has been chosen!"
                write (*, *) "Defined Mol Types:", nMolTypes
                write (*, *) "Chosen Mol Type:", indxVar
                stop
             end if
             if (indxVar .gt. nMolTypes) then
                write (*, *) "Error! An invalid molecule type has been chosen!"
                write (*, *) "Defined Mol Types:", nMolTypes
                write (*, *) "Chosen Mol Type:", indxVar
                stop
             end if
             biasvar(iUmbrella)%varType = 1
             biasvar(iUmbrella)%intVar => NPart(indxVar)
             biasvarnew(iUmbrella)%varType = 1
             biasvarnew(iUmbrella)%intVar => NPart_New(indxVar)
             binMax(iUmbrella) = NMAX(indxVar)
             binMin(iUmbrella) = NMin(indxVar)
             UBinSize(iUmbrella) = 1E0
             outputFormat(iUmbrella) = "2x,F5.1,"

!          iSwapIn = iSwapIn + 1
!          iSwapOut = iSwapOut + 1
          case ("totalclustersize")
             biasvar(iUmbrella)%varType = 1
             biasvar(iUmbrella)%intVar => NTotal
             biasvarnew(iUmbrella)%varType = 1
             biasvarnew(iUmbrella)%intVar => NTotal_New
             binMax(iUmbrella) = maxMol
             binMin(iUmbrella) = 1
             UBinSize(iUmbrella) = 1E0
             outputFormat(iUmbrella) = "2x,F5.1,"

!          iSwapIn = iSwapIn + 1
!          iSwapOut = iSwapOut + 1
          case ("pairdist")
             indxVar = 0
             read (inputLines(iUmbrella), *) labelField, indxVar, valMin, valMax, binSize
             if (nDistPair .le. 0) then
                write (*, *) "Error! An invalid distance variable has been chosen!"
                write (*, *) "Defined Distance Pairs:", nDistPair
                write (*, *) "Chosen Distance Pair:", indxVar
                stop
             end if
             if (indxVar .gt. nDistPair) then
                write (*, *) "Error! An invalid distance variable has been chosen!"
                write (*, *) "Defined Distance Pairs:", nDistPair
                write (*, *) "Chosen Distance Pair:", indxVar
                stop
             end if
             biasvar(iUmbrella)%varType = 2
             biasvar(iUmbrella)%realVar => miscCoord(pairArrayIndx(indxVar))
             biasvarnew(iUmbrella)%varType = 2
             biasvarnew(iUmbrella)%realVar => miscCoord_New(pairArrayIndx(indxVar))
             UBinSize(iUmbrella) = binSize
             binMin(iUmbrella) = nint(valMin/binSize)
             binMax(iUmbrella) = nint(valMax/binSize)
             outputFormat(iUmbrella) = "2x,F12.6,"

             iDisp = iDisp + 1
             DispUmbrella(iDisp)%func => CalcDistPairs_New
             iSwapIn = iSwapIn + 1
             SwapInUmbrella(iSwapIn)%func => CalcDistPairs_SwapIn
             iSwapOut = iSwapOut + 1
             SwapOutUmbrella(iSwapOut)%func => CalcDistPairs_SwapOut

          case ("q6")
             read (inputLines(iUmbrella), *) labelField, valMin, valMax, binSize
             biasvar(iUmbrella)%varType = 2
             biasvar(iUmbrella)%realVar => miscCoord(q6ArrayIndx)
             biasvarnew(iUmbrella)%varType = 2
             biasvarnew(iUmbrella)%realVar => miscCoord_New(q6ArrayIndx)
             UBinSize(iUmbrella) = binSize
             binMin(iUmbrella) = nint(valMin/binSize)
             binMax(iUmbrella) = nint(valMax/binSize)
             outputFormat(iUmbrella) = "2x,F12.6,"

             iDisp = iDisp + 1
             DispUmbrella(iDisp)%func => CalcQ6_Disp
             iSwapIn = iSwapIn + 1
             SwapInUmbrella(iSwapIn)%func => CalcQ6_SwapIn
             iSwapOut = iSwapOut + 1
             SwapOutUmbrella(iSwapOut)%func => CalcQ6_SwapOut
!          write(*,*) labelField, valMin, valMax, binSize
          case default
             write (*, *) "ERROR! Invalid variable type specified in input file"
             write (*, *) umbrellaName
             stop
          end select
       end do
       write (screenFormat, *) "(", (trim(adjustl(outputFormat(iUmbrella))), iUmbrella=1, nBiasVariables), "1x)"

       allocate (varValues(1:nBiasVariables))
       allocate (UArray(1:nBiasVariables))
       deallocate (inputLines)

       call AllocateUmbrellaArray
       call ReadInitialBias

!      Use the input to specify the reference bin
       allocate (refSizeNumbers(1:nBiasVariables), STAT=AllocateStatus)
       read (refLine, *) labelField, (refVals(iUmbrella), iUmbrella=1, nBiasVariables)
       call getUIndexArray(refVals, refBin, stat)
       do iUmbrella = 1, nBiasVariables
          refSizeNumbers(iUmbrella) = refVals(iUmbrella)
       end do
       if (stat .eq. 1) then
          write (*, *) "ERROR! An invalid ref bin was specified!"
          write (*, *) refLine
          stop
       end if

       deallocate (refVals)

    end subroutine

!======================================================
    subroutine ReadAnalysisInput(fileUnit)
       use MiscelaniousVars
       use SimpleDistPair, only: nDistPair, pairArrayIndx
       use SimParameters, only: NMAX, NMIN, NPART, NPart_New, nMolTypes
       use Q6Functions, only: CalcQ6, Initialize_q6, useQ6, q6Dist, q6DistSq
       implicit none
       integer, intent(in) :: fileUnit
       integer :: iAnalysis, AllocateStatus
       integer :: indxVar, nBins
       integer :: iRadial, iDistPair, iRadDens, iPostMove, iOutput
       integer :: type1, type2, mol1, mol2, atom1, atom2
       real(dp) :: binSize, realVar
       character(len=200), allocatable :: inputLines(:)
       character(len=30) :: labelField
       character(len=30) :: analysisName, fileName

       read (fileUnit, *) labelField, nAnalysisVar
       if (nAnalysisVar .lt. 0) then
          write (*, *) "ERROR! The user has specified an invalid number of Analysis Variables"
          write (*, *) labelField, nAnalysisVar
          stop
       end if
       if (nAnalysisVar .eq. 0) then
          useAnalysis = .false.
          return
       else
          useAnalysis = .true.
       end if

       allocate (inputLines(1:nAnalysisVar), stat=AllocateStatus)
!      Due to the nature how the analysis arrays are structured, it is nessiary
!      to first gather the lines
       do iAnalysis = 1, nAnalysisVar
          read (fileUnit, "(A)") inputLines(iAnalysis)
!        write(*,*) inputLines(iAnalysis)
       end do

!      Begin by counting how many entries are required for each function array.
!      These will be used in the next step to allocate the arrays.
       nPostMove = 0
       nOutput = 0
       nDistPair = 0
       nRadialDist = 0
       nRadialDens = 0
       do iAnalysis = 1, nAnalysisVar
          read (inputLines(iAnalysis), *) analysisName

          select case (trim(adjustl(analysisName)))
          case ("radialdistribution")
             nRadialDist = nRadialDist + 1
             if (nRadialDist .eq. 1) then
                nPostMove = nPostMove + 1
                nOutput = nOutput + 1
             end if
          case ("pairdist")
             nDistPair = nDistPair + 1
             if (nDistPair .eq. 1) then
                nPostMove = nPostMove + 1
             end if
          case ("q6")
             if (useQ6 .eqv. .false.) then
                useQ6 = .true.
                nPostMove = nPostMove + 1
             else
                stop "ERROR! The Q6 analysis function has been defined more than once in the input script."
             end if
          case ("radialdensity")
             nRadialDens = nRadialDens + 1
             if (nRadialDens .eq. 1) then
                nPostMove = nPostMove + 1
                nOutput = nOutput + 1
             end if
          case default
             write (*, *) "ERROR! Invalid variable type specified in input file"
             write (*, *) analysisName
             stop
          end select
       end do

!      Now that we know how much memory is required, allocate all relevent arrays
       if (nPostMove .ne. 0) then
          allocate (postMoveArray(1:nPostMove))
       end if
       if (nOutput .ne. 0) then
          allocate (outputArray(1:nOutput))
       end if
       allocate (loadUmbArray(1:nAnalysisVar))
       call Initialize_RadialDist
       call Initialize_DistPair
       call Initialize_RadialDens
       call Initialize_Q6
       call AllocateMiscArrays

       iOutPut = 0
       iPostMove = 0
       iRadial = 0
       iDistPair = 0
       iRadDens = 0
       do iAnalysis = 1, nAnalysisVar
          read (inputLines(iAnalysis), *) analysisName

          select case (trim(adjustl(analysisName)))
          case ("radialdistribution")
             iRadial = iRadial + 1
             read (inputLines(iAnalysis), *) analysisName, type1, atom1, type2, atom2, binSize, nBins, fileName

             call SetRadialParameters(iRadial, type1, type2, atom1, atom2)
             call SetRadialHist(iRadial, binSize, nBins, fileName)

             if (iRadial .eq. 1) then
                iPostMove = iPostMove + 1
                postMoveArray(iPostMove)%func => Calc_RadialDist
                iOutPut = iOutPut + 1
                outputArray(iOutPut)%func => Output_RadialDist
             end if
          case ("pairdist")
             iDistPair = iDistPair + 1
             read (inputLines(iAnalysis), *) analysisName, type1, mol1, atom1, type2, mol2, atom2
             call SetPairVariables(iDistPair, Type1, Mol1, Atom1, Type2, Mol2, Atom2)
             if (iDistPair .eq. 1) then
                iPostMove = iPostMove + 1
                postMoveArray(iPostMove)%func => CalcDistPairs
             end if
          case ("q6")
             read (inputLines(iAnalysis), *) analysisName, realVar
             iPostMove = iPostMove + 1
             q6Dist = realVar
             q6DistSq = q6Dist*q6Dist
             postMoveArray(iPostMove)%func => CalcQ6
             loadUmbArray(iAnalysis)%func => UmbrellaVar_Q6
          case ("radialdensity")
             read (inputLines(iAnalysis), *) analysisName, type1, binSize, nBins, fileName
             iRadDens = iRadDens + 1
             call SetDensityParameters(iRadDens, type1)
             call SetDensityHist(iRadDens, binSize, nBins, fileName)
             if (iRadDens .eq. 1) then
                iPostMove = iPostMove + 1
                iOutPut = iOutPut + 1
                postMoveArray(iPostMove)%func => Calc_RadialDensity
                outputArray(iOutPut)%func => Output_RadialDensity
             end if
          case default
             write (*, *) "ERROR! Invalid variable type specified in input file"
             write (*, *) analysisName
             stop
          end select
       end do

       call AllocateHistBins

       deallocate (inputLines)

    end subroutine
!======================================================
    subroutine ReadInput_MCMove(fileUnit)
       implicit none
       integer, intent(in) :: fileUnit
       integer :: iMoves, AllocateStatus
       real(dp) :: norm
       character(len=30) :: labelField
       character(len=30) :: moveName_temp

       norm = 0d0
       read (fileUnit, *) labelField, nMoveTypes
       if (nMoveTypes .le. 0) then
          write (*, *) "ERROR! The user has specified an invalid number of Monte Carlo moves"
          write (*, *) "Please specify at least one valid Monte Carlo move to continue"
          write (*, *) labelField, nMoveTypes
          stop
       end if

       allocate (mcMoveArray(1:nMoveTypes), STAT=AllocateStatus)
       allocate (moveProbability(1:nMoveTypes), STAT=AllocateStatus)
       allocate (movesAccepted(1:nMoveTypes), STAT=AllocateStatus)
       allocate (movesAttempt(1:nMoveTypes), STAT=AllocateStatus)
       allocate (accptRate(1:nMoveTypes), STAT=AllocateStatus)
       allocate (moveName(1:nMoveTypes), STAT=AllocateStatus)
       norm = 0d0
       avbmcUsed = .false.
       cbmcUsed = .false.
       do iMoves = 1, nMoveTypes
          read (fileUnit, *) moveName_temp, moveProbability(iMoves)
          norm = norm + moveProbability(iMoves)
          select case (trim(adjustl(moveName_temp)))
          case ("translation")
             mcMoveArray(iMoves)%moveFunction => Translation
             moveName(iMoves) = "Translation"
          case ("rotation")
             mcMoveArray(iMoves)%moveFunction => Rotation
             moveName(iMoves) = "Rotation"
          case ("avbmc")
             mcMoveArray(iMoves)%moveFunction => AVBMC
             moveName(iMoves) = "AVBMC"
             avbmcUsed = .true.
          case ("cbmc")
             mcMoveArray(iMoves)%moveFunction => CBMC
             moveName(iMoves) = "CBMC"
             cbmcUsed = .true.
          case ("exchange")
             mcMoveArray(iMoves)%moveFunction => Exchange
             moveName(iMoves) = "Exchange"
          case ("singleatom_translation")
             mcMoveArray(iMoves)%moveFunction => SingleAtom_Translation
             moveName(iMoves) = "Single Atom Translation"
          case default
             write (*, *) "ERROR! Invalid move type specified in input file"
             write (*, *) moveName, moveProbability(iMoves)
             stop
          end select
!        moveName(i) = moveName_temp
       end do

       do iMoves = 1, nMoveTypes
          moveProbability(iMoves) = moveProbability(iMoves)/norm
       end do
       if (nMoveTypes .gt. 1) then
          do iMoves = 2, nMoveTypes
             moveProbability(iMoves) = moveProbability(iMoves) + moveProbability(iMoves - 1)
          end do
       end if

    end subroutine
!===========================================================================
