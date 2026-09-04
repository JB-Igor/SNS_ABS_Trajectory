#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma DefaultTab={3,20,4}

//==============================================================================
// SNS_SpatialMaps.ipf
//
// Spatial LDOS workflows built on SNS ray tracing and DOS assembly.
//
// Responsibilities:
//   - 1D line maps
//   - 2D area maps
//   - fixed-B LDOS maps
//   - point/angle sweeps
//   - spatial sampling helpers
//
// Dependencies:
//   SNS_Core.ipf
//   SNS_Logging.ipf
//   SNS_GeometryFromMask.ipf
//   SNS_DOS.ipf
//   SNS_Broadening.ipf
//
// Notes:
//   This file orchestrates repeated spatial calculations.
//   It should not contain low-level ABS solver code, interface transparency,
//   broadening kernels, or display-only helpers.
//==============================================================================

//==============================================================================
// 1D LINE WORKFLOWS
//
// Functions in this section compute LDOS/DOS along a line through the mask.
// They repeatedly build local channel ensembles at line positions and assemble
// spatially resolved spectra.
//==============================================================================

//==============================================================================
// SNS_LineLDOS_BSweep_FromMask
//
// Purpose:
//   Compute LDOS(E, r, B) along a straight 1D line through a 2D N-region mask.
//
//   For each spatial point r along the line:
//     1. build SNS ray channels from Nmask,
//     2. solve the ABS/DOS for each B value,
//     3. apply experimental T + lock-in modulation broadening,
//     4. store the result in a 3D LDOS wave.
//
//   Output dimension order:
//       x : energy E [eV]
//       y : line position r [nm]
//       z : magnetic field B [T]
//
// Inputs:
//   Nmask      : 2D binary N-region mask.
//                Inside N: > 0.5, outside: <= 0.5.
//                Wave axes must be scaled in nm.
//
//   phiB       : in-plane magnetic-field angle [rad]
//
//   xStart,
//   yStart     : line start point [nm]
//
//   xEnd,
//   yEnd       : line end point [nm]
//
//   Bmin       : minimum magnetic field [T]
//   Bmax       : maximum magnetic field [T]
//
//   NBcenters  : number of B points between Bmin and Bmax
//
//   stepFac    : ray/channel sampling refinement factor.
//                Used by SNS_BuildChannelsFromMask2D.
//
//   Zbarrier   : legacy input, currently unused.
//                Interface transparency is now read from SNS_Settings via
//                SNS_ChannelTransmissionFromCos().
//
// Outputs:
//   LDOS_B...  : LDOS(E,r,B) wave
//   E_axis_... : energy axis [eV]
//   R_axis_... : line-position axis [nm]
//   B_axis_... : magnetic-field axis [T]
//
// Returns:
//   Nr         : number of line positions.
//
// Dependencies:
//   SNS_Core.ipf
//     - SNS_Params
//     - SNS_LoadParams
//
//   SNS_GeometryFromMask.ipf
//     - SNS_BuildLinePositions
//     - SNS_BuildChannelsFromMask2D
//
//   SNS_Broadening.ipf
//     - SNS_ApplyDOS_Broadening_TplusMod
//
//   SNS_Maps.ipf / legacy residual
//     - Compute_DOS_SNS_Map_FromChannels_Coh
//
// Notes:
//   This is a spatial workflow function and belongs in SNS_SpatialMaps.ipf.
//   Long term, replace Compute_DOS_SNS_Map_FromChannels_Coh with
//   SNS_ComputeDOS_FromChannels / SNS_ComputeDOS_FromSettings.
//==============================================================================
Function SNS_LineLDOS_BSweep_FromMask(Nmask, phiB, xStart, yStart, xEnd, yEnd, Bmin, Bmax, NBcenters, stepFac, Zbarrier)
    Wave   Nmask
    Variable phiB
    Variable xStart, yStart, xEnd, yEnd      // nm
    Variable Bmin, Bmax                      // T
    Variable NBcenters
    Variable stepFac, Zbarrier

    //---------- 0) Get current folder reference ----------
    DFREF dfrCaller = GetDataFolderDFR()    
    
    // --- Load standard SNS settings (Igor Pro 9) ---
    STRUCT SNS_Params params
    SNS_LoadParams(params)

    Variable lambdaF    = params.LambdaF      // [m]
    Variable Delta      = params.Delta        // [eV]
    Variable vF         = params.vF           // [m/s]
    Variable lambdaL    = params.lambdaL      // [m]
    Variable Broadening = params.Broadening   // [eV]
    Variable NE         = params.NE

    // If you need Nphi later, you can re-enable this:
    // Variable Nphi = round(SNS_EstimateNphi_FromMask(Nmask, lambdaF*1e9))

    // ---------- 1) Auto-generated names ----------
    String tag
    tag = "B" + num2str(round(Bmin*1e3)) + "to" + num2str(round(Bmax*1e3)) + "mT_" \
          + num2str(round(phiB*180/pi)) + "deg"

    String nameLDOS   = "LDOS_"   + tag
    String nameEaxis  = "E_axis_" + tag
    String nameRaxis  = "R_axis_" + tag
    String nameBaxis  = "B_axis_" + tag

    // ---------- 2) Build line positions r across N-region ----------
    String nameX = "SNS_Line_X"
    String nameY = "SNS_Line_Y"

    Variable Nr = SNS_BuildLinePositions(xStart, yStart, xEnd, yEnd, lambdaF, \
                                         nameX, nameY, nameRaxis)
    Wave Xline = $nameX
    Wave Yline = $nameY
    Wave Rline = $nameRaxis

    if (Nr > 1)
        Variable drnm = (Rline[Nr-1] - Rline[0]) / (Nr - 1)
        SetScale/P x, Rline[0], drnm, "nm", Rline
    else
        SetScale/P x, 0, 1, "nm", Rline
    endif

    // ---------- 3) B centers axis ----------
    Make/O/D/N=(NBcenters) $nameBaxis
    Wave B_centers = $nameBaxis

    Variable dBcent
    if (NBcenters > 1)
        dBcent = (Bmax - Bmin) / (NBcenters - 1)
    else
        dBcent = 0
    endif
    B_centers = Bmin + p*dBcent
    SetScale/P x, Bmin, (NBcenters>1 ? dBcent : 1), "T", B_centers

    // ---------- 4) Allocate LDOS(E,r,B) ----------
    Make/O/D/N=(NE, Nr, NBcenters) $nameLDOS
    Wave LDOS = $nameLDOS

    // ---------- 5) Temp folder & single-B array ----------
    String chanFolder = "root:SNS_LineTmp"
    NewDataFolder/O $chanFolder
    SetDataFolder $chanFolder

    // Single B value per solve
    Make/O/D/N=1 Bwin

    String nameDOS_EB    = "DOS_local_EB"
    String nameEaxisLoc  = "E_axis_local"
    String nameDOS_EBbr  = "DOS_local_EB_broad"

    Variable ir, k, iE
    Variable Nch

    // ---------- 6) Loop over positions along the line ----------
    for (ir = 0; ir < Nr; ir += 1)

        Variable r0x = Xline[ir]
        Variable r0y = Yline[ir]

        // Build channels for this position
        SNS_BuildChannelsFromMask2D(Nmask, r0x, r0y, phiB, \
                                         stepFac, chanFolder)

        Wave L_N_List_nm
        Wave W_eff_List_nm
        Wave wChan
        Wave T_eff_List

        Nch = DimSize(L_N_List_nm, 0)
        if ( (Nch <= 0) || (Nch != DimSize(W_eff_List_nm,0)) || \
             (Nch != DimSize(wChan,0)) || (Nch != DimSize(T_eff_List,0)) )

            // No valid channels → mark whole column in B as NaN
            SetDataFolder root:
            LDOS[][ir][] = NaN
            SetDataFolder $chanFolder
            continue
        endif

        // nm → m for solver
        Duplicate/O L_N_List_nm,   L_N_List_m
        Duplicate/O W_eff_List_nm, W_eff_List_m
        L_N_List_m   *= 1e-9
        W_eff_List_m *= 1e-9

        // ---------- 7) Loop over B centers (no averaging) ----------
        for (k = 0; k < NBcenters; k += 1)

            Variable B0 = B_centers[k]

            // Single-field evaluation
            Bwin[0] = B0

            // Raw DOS(E,B) from coherent solver; B dimension has size 1
            Compute_DOS_SNS_Map_FromChannels_Coh(Bwin, L_N_List_m, W_eff_List_m, \
                                                 wChan, T_eff_List, \
                                                 Delta, vF, lambdaL, Broadening, NE, \
                                                 nameDOS_EB, nameEaxisLoc)

            Wave DOS_EB_raw   = $nameDOS_EB      // [E, B=0]
            Wave E_axis_local = $nameEaxisLoc    // [E]

            // Experimental T+mod broadening
            SNS_ApplyDOS_Broadening_TplusMod(DOS_EB_raw, nameDOS_EBbr)

            Wave DOS_EB_broad = $nameDOS_EBbr    // [E, B=0] broadened

            // Copy energy axis once to root
            if ( (ir == 0) && (k == 0) )
                SetDataFolder root:
                Duplicate/O E_axis_local, $nameEaxis
                SetDataFolder $chanFolder
            endif

            // Store LDOS(E, r_ir, B_k)
            SetDataFolder root:
            for (iE = 0; iE < NE; iE += 1)
                LDOS[iE][ir][k] = DOS_EB_broad[iE][0]
            endfor
            SetDataFolder $chanFolder

        endfor // k over B centers

    endfor // ir over r

    // ---------- 8) Attach axes to LDOS ----------
    SetDataFolder root:

    Wave E_axis = $nameEaxis
    // x: energy
    SetScale/P x, DimOffset(E_axis, 0), DimDelta(E_axis, 0), WaveUnits(E_axis, 0), LDOS
    // y: r
    SetScale/P y, DimOffset(Rline, 0),  DimDelta(Rline, 0),  WaveUnits(Rline, 0),  LDOS
    // z: B
    SetScale/P z, DimOffset(B_centers, 0), DimDelta(B_centers, 0), WaveUnits(B_centers, 0), LDOS

    SetDataFolder dfrCaller
    return Nr
End

//==============================================================================
// LINE-CACHE HELPERS
//
// Small utilities used by line-resolved workflows to preserve per-position
// trajectory state for later UI-based trajectory lookup.
//==============================================================================

//==============================================================================
// SNS_LineDOS_DupWaveIfExists
//
// Purpose:
//   Duplicate a wave only if the source wave exists.
//
// Inputs:
//   srcPath : full or relative Igor path to source wave.
//   dstPath : full or relative Igor path/name for destination wave.
//
// Outputs:
//   dstPath : duplicated wave, only if srcPath exists.
//
// Returns:
//   1 if the wave was duplicated.
//   0 if srcPath does not exist.
//
// Notes:
//   Generic helper used by spatial DOS/LDOS workflows.
//   Does not change the current data folder.
//==============================================================================
Function SNS_LineDOS_DupWaveIfExists(srcPath, dstPath)
    String srcPath, dstPath

    Wave/Z w = $srcPath
    if (WaveExists(w))
        Duplicate/O w, $dstPath
        return 1
    endif

    return 0
End

//==============================================================================
// SNS_LineDOS_SaveTrajectoryState
//
// Purpose:
//   Save per-line-position channel and branch state needed by the interactive
//   trajectory picker.
//
//   For line position index ir, this creates/replaces:
//
//       lineCacheFolder:r####:
//
//   and copies the channel waves and branch-energy waves from chanFolder.
//
// Inputs:
//   chanFolder       : source folder containing current-position channel waves.
//   lineCacheFolder  : destination parent folder for cached line-position states.
//   ir               : line-position index.
//
// Copied channel/state waves, if present:
//   Hit1x_List, Hit1y_List
//   Hit2x_List, Hit2y_List
//   wChan
//   T_eff_List
//   L_N_List, W_eff_List
//   L_N_List_m, W_eff_List_m
//   betaExtra_List
//
// Copied branch waves:
//   E_allBranches:ch###
//
// Outputs:
//   lineCacheFolder:r####: channel/state waves
//   lineCacheFolder:r####:E_allBranches:ch###
//
// Returns:
//    0  success
//   -1  missing source wChan; branch copy skipped
//
// Dependencies:
//   SNS_LineDOS_TrailingColon
//   SNS_LineDOS_DupWaveIfExists
//
// Notes:
//   This is used by SNS_UI for LDOS(E,r) click-to-trajectory lookup.
//   It does not solve spectra or compute DOS.
//==============================================================================
Function SNS_LineDOS_SaveTrajectoryState(chanFolder, lineCacheFolder, ir)
    String chanFolder, lineCacheFolder
    Variable ir

    String srcDF = SNS_LineDOS_TrailingColon(chanFolder)
    String baseDF = SNS_LineDOS_TrailingColon(lineCacheFolder)

    String dstDF
    sprintf dstDF, "%sr%04d", baseDF, ir

    KillDataFolder/Z $dstDF
    NewDataFolder/O $dstDF
    dstDF = SNS_LineDOS_TrailingColon(dstDF)

    // Channel geometry/state
    SNS_LineDOS_DupWaveIfExists(srcDF + "Hit1x_List_nm", dstDF + "Hit1x_List_nm")
    SNS_LineDOS_DupWaveIfExists(srcDF + "Hit1y_List_nm", dstDF + "Hit1y_List_nm")
    SNS_LineDOS_DupWaveIfExists(srcDF + "Hit2x_List_nm", dstDF + "Hit2x_List_nm")
    SNS_LineDOS_DupWaveIfExists(srcDF + "Hit2y_List_nm", dstDF + "Hit2y_List_nm")
    SNS_LineDOS_DupWaveIfExists(srcDF + "wChan",         dstDF + "wChan")
    SNS_LineDOS_DupWaveIfExists(srcDF + "T_eff_List",    dstDF + "T_eff_List")
    SNS_LineDOS_DupWaveIfExists(srcDF + "L_N_List_nm",   dstDF + "L_N_List_nm")
    SNS_LineDOS_DupWaveIfExists(srcDF + "W_eff_List_nm", dstDF + "W_eff_List_nm")
    SNS_LineDOS_DupWaveIfExists(srcDF + "betaExtra_List", dstDF + "betaExtra_List")

    // Branch energies
    NewDataFolder/O $(dstDF + "E_allBranches")

    Wave/Z wChan = $(srcDF + "wChan")
    if (!WaveExists(wChan))
        return -1
    endif

    Variable ch, Nch = numpnts(wChan)
    String srcW, dstW

    for (ch = 0; ch < Nch; ch += 1)
        sprintf srcW, "%sE_allBranches:ch%03d", srcDF, ch
        sprintf dstW, "%sE_allBranches:ch%03d", dstDF, ch
        SNS_LineDOS_DupWaveIfExists(srcW, dstW)
    endfor

    return 0
End

//==============================================================================
// FIXED-FIELD LINE LDOS
//
// 2D and 3D fixed-B line workflows. These produce DOS(E,r) and broadened
// DOS_Conv(E,r), plus trajectory diagnostics such as W_max, L_max, and selected
// ABS energies.
//==============================================================================

//==============================================================================
// SNS_LineDOS_BFixed_FromMask
//
// Purpose:
//   Compute a 1D line-resolved LDOS/DOS profile at fixed magnetic field B0:
//
//       DOS(E, r) and experimentally broadened DOS_Conv(E, r)
//
//   For each position r along the line, the function:
//     1. builds SNS channels from the 2D N-region mask,
//     2. optionally computes an extra phase betaExtra_List from a vortex /
//        phase-field / Q-field model,
//     3. computes DOS(E,B0) from the channel ensemble,
//     4. applies thermal + lock-in modulation broadening,
//     5. stores line-position diagnostics used for trajectory inspection.
//
// Special B=0 behavior:
//   If B0 == 0 and ScreeningModel is omitted, no phase-field / vortex correction
//   is used. The function skips phase-field folder/wave checks and uses simple
//   output names:
//
//       LDOS_B0mT
//       LDOS_Conv_B0mT
//       E_axis_B0mT
//       R_axis_B0mT
//
//   If ScreeningModel is explicitly supplied, the selected model is used even at
//   B0 == 0. In that case maskFolder is required and the output tag keeps the
//   angle/model suffix.
//
// Inputs:
//   Nmask        : 2D binary N-region mask.
//                  Inside N: > 0.5, outside: <= 0.5.
//                  Wave axes must be scaled in nm.
//
//   phiB         : in-plane magnetic-field angle [rad]
//
//   xStart,
//   yStart       : line start point [nm]
//
//   xEnd,
//   yEnd         : line end point [nm]
//
//   B0           : fixed magnetic field [T]
//
//   stepFac      : ray/channel sampling refinement factor passed to
//                  SNS_BuildChannelsFromMask2D.
//
// Optional Inputs:
//   BTK_barrier  : legacy argument, currently only stored in metadata.
//                  Interface transparency is now controlled globally through
//                  SNS_Settings and SNS_ChannelTransmissionFromCos().
//
//   maskFolder   : optional parent folder containing phase-field model
//                  subfolders. Required for ScreeningModel 1--4. For model 0,
//                  supplied PhaseReFree/PhaseImFree waves are reused; when
//                  absent, the analytical free-vortex phase is generated in a
//                  temporary calculation folder.
//
//   ScreeningModel:
//                  phase / Q-field selector:
//
//                    0 = free vortex phase
//                    1 = stiffness-corrected vortex phase
//                    2 = corrected phase-gradient Q
//                    3 = screened corrected Q
//                    4 = local-London corrected Q
//
//                  If omitted and B0 == 0, no phase model is used.
//                  If supplied, the selected model is used even at B0 == 0.
//
//   NLinePts     : optional number of line positions.
//                  If supplied, overrides the default LambdaF-based line step.
//                  Must be finite integer >= 2.
//
// Outputs:
//   LDOS_*                 : raw DOS(E,r)
//   LDOS_Conv_*            : broadened DOS(E,r)
//   E_axis_*               : energy axis [eV]
//   R_axis_*               : line-position axis [nm]
//   W_max_*                : maximum W_eff per line point
//   L_max_*                : maximum L_N per line point
//   betaExtra_min_*        : smallest |betaExtra| per line point
//   Energy_betaExtra_min_* : ABS energy for channel with smallest |betaExtra|
//   Energy_wmax_*          : lowest |E| branch for largest-W_eff channel
//   Energy_Lmax_*          : lowest |E| branch for longest-L_N channel
//
// Returns:
//   Nr                     : number of line positions.
//
// Dependencies:
//   SNS_Core.ipf
//     - SNS_Params
//     - SNS_LoadParams
//
//   SNS_GeometryFromMask.ipf
//     - SNS_BuildLinePositions
//     - SNS_BuildChannelsFromMask2D
//     - SNS_MaskAreaPerim_FromParticles
//
//   SNS_DOS.ipf
//     - SNS_ComputeDOS_FromChannels
//
//   SNS_Broadening.ipf
//     - SNS_ApplyDOS_Broadening_TplusMod
//
//   SNS_PhaseFields / residual SNS_Maps.ipf
//     - SNS_ComputeBetaExtraFromExteriorPhase2D
//     - SNS_ComputeBetaExtraFromQField2D
//
//   SNS_UI.ipf
//     - consumes trajectory cache written by SNS_LineDOS_SaveTrajectoryState
//
// Notes:
//   This is a 1D spatial fixed-field workflow and belongs in
//   SNS_SpatialMaps.ipf.
//==============================================================================
Function SNS_LineDOS_BFixed_FromMask(Nmask, phiB, xStart, yStart, xEnd, yEnd, B0, stepFac, [maskFolder, ScreeningModel, NLinePts])

    Wave   Nmask
    Variable phiB
    Variable xStart, yStart, xEnd, yEnd
    Variable B0
    Variable stepFac
    String maskFolder
    Variable ScreeningModel
    Variable NLinePts

    DFREF dfrCaller = GetDataFolderDFR()

    Variable isZeroField = (abs(B0) <= 1e-15)
    Variable hasScreeningModel = !ParamIsDefault(ScreeningModel)

    Variable smodel = 0
    if (hasScreeningModel)
        smodel = ScreeningModel
    endif
    smodel = round(smodel)

    if ((smodel < 0) || (smodel > 4))
        Abort "SNS_LineDOS_BFixed_FromMask: ScreeningModel must be 0, 1, 2, 3, or 4."
    endif

    if (hasScreeningModel && smodel != 0 && ParamIsDefault(maskFolder))
        SetDataFolder dfrCaller
        Abort "SNS_LineDOS_BFixed_FromMask: maskFolder is required for ScreeningModel 1, 2, 3, or 4."
    endif

    // Use phase-field logic if a maskFolder is supplied and either:
    //   1) B0 != 0, or
    //   2) ScreeningModel was explicitly selected.
    // Only the model-free B0 case skips phase-field checks.
    Variable usePhaseField = (hasScreeningModel && smodel == 0) || ((!ParamIsDefault(maskFolder)) && ((!isZeroField) || hasScreeningModel))

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    Variable lambdaF_nm = params.LambdaF
    Variable NE         = params.NE

    // ---------- geometry folder ----------
    DFREF dfrGeom = GetWavesDataFolderDFR(Nmask)
    SetDataFolder dfrGeom

    Wave/Z Vortex_ptx = Vortex_ptx
    Wave/Z Vortex_pty = Vortex_pty

    Variable xV_nm = 0
    Variable yV_nm = 0
    if (WaveExists(Vortex_ptx) && WaveExists(Vortex_pty))
        xV_nm = Vortex_ptx[0]
        yV_nm = Vortex_pty[0]
    endif

    Wave/Z w_area_nm2 = SNS_MaskAreaPerim_FromParticles(Nmask)

    Variable N_cont_statesPer_eV = 0
    if (WaveExists(w_area_nm2))
        Variable area_nm2 = w_area_nm2[0]
        Variable area_m2  = area_nm2 * 1e-18
        N_cont_statesPer_eV = params.DOS2D_eV_Area * area_m2
    endif

    SetDataFolder dfrCaller

    // ---------- output tag ----------
    String bTag = "B" + num2str(round(B0*1e3)) + "mT"
    String angleTag = num2str(round(phiB*180/pi)) + "deg"
    String tag

    if (isZeroField && !hasScreeningModel)
        tag = "B0mT"
    elseif (smodel == 1)
        tag = bTag + "_" + angleTag + "_PhaseCorr"
    elseif (smodel == 2)
        tag = bTag + "_" + angleTag + "_QCorr"
    elseif (smodel == 3)
        tag = bTag + "_" + angleTag + "_ScreenedCorr"
    elseif (smodel == 4)
        tag = bTag + "_" + angleTag + "_LondonCorr"
    else
        tag = bTag + "_" + angleTag + "_PhaseFree"
    endif

    String nameDOSLine           = "LDOS_"                 + tag
    String nameDOSLineConv       = "LDOS_Conv_"            + tag
    String nameEaxis             = "E_axis_"               + tag
    String nameRaxis             = "R_axis_"               + tag
    String nameWmax              = "W_max_"                + tag
    String nameLmax              = "L_max_"                + tag
    String nameBetaExtraMin      = "betaExtra_min_"        + tag
    String nameEnergyBetaMin     = "Energy_betaExtra_min_" + tag
    String nameEnergyWmax        = "Energy_wmax_"          + tag
    String nameEnergyLmax        = "Energy_Lmax_"          + tag

    // ---------- model subfolders ----------
    String maskBase = ""
    String stiffnessFolder = ""
    String screenedFolder = ""
    String londonFolder = ""

    if (!ParamIsDefault(maskFolder))
        maskBase = SNS_LineDOS_TrailingColon(maskFolder)
        stiffnessFolder = maskBase + "StiffnessPhase"
        screenedFolder  = maskBase + "ScreenedStiffnessPhase"
        londonFolder    = maskBase + "LocalLondonStiffnessPhase"
    endif

    // ---------- line positions ----------
    String nameX = "SNS_Line_X"
    String nameY = "SNS_Line_Y"

    Variable lineStep_nm = lambdaF_nm

    if (!ParamIsDefault(NLinePts))
        NLinePts = round(NLinePts)

        if ((numtype(NLinePts) != 0) || (NLinePts < 2))
            SetDataFolder dfrCaller
            Abort "SNS_LineDOS_BFixed_FromMask: NLinePts must be a finite integer >= 2."
        endif

        Variable lineLength_nm = sqrt((xEnd - xStart)^2 + (yEnd - yStart)^2)

        if ((numtype(lineLength_nm) != 0) || (lineLength_nm <= 0))
            SetDataFolder dfrCaller
            Abort "SNS_LineDOS_BFixed_FromMask: cannot use NLinePts for a zero-length or invalid line."
        endif

        lineStep_nm = lineLength_nm / NLinePts
    endif

    Variable lineStep_m = lineStep_nm * 1e-9

    Variable Nr = SNS_BuildLinePositions(xStart, yStart, xEnd, yEnd, lineStep_m, nameX, nameY, nameRaxis)

    Wave Xline = $nameX
    Wave Yline = $nameY
    Wave Rline = $nameRaxis

    if (Nr > 1)
        Variable drnm = (Rline[Nr-1] - Rline[0]) / (Nr - 1)
        SetScale/P x, Rline[0], drnm, "nm", Rline
    else
        SetScale/P x, 0, 1, "nm", Rline
    endif

    // ---------- allocate outputs ----------
    Make/O/D/N=(NE, Nr) $nameDOSLine
    Wave DOS_Line = $nameDOSLine
    DOS_Line = NaN

    Make/O/D/N=(NE, Nr) $nameDOSLineConv
    Wave DOS_Line_Conv = $nameDOSLineConv
    DOS_Line_Conv = NaN

    Make/O/D/N=(Nr) $nameWmax
    Make/O/D/N=(Nr) $nameLmax
    Make/O/D/N=(Nr) $nameBetaExtraMin
    Make/O/D/N=(Nr) $nameEnergyBetaMin
    Make/O/D/N=(Nr) $nameEnergyWmax
    Make/O/D/N=(Nr) $nameEnergyLmax

    Wave W_max_Line            = $nameWmax
    Wave L_max_Line            = $nameLmax
    Wave betaExtra_min_Line    = $nameBetaExtraMin
    Wave Energy_beta_min_Line  = $nameEnergyBetaMin
    Wave Energy_wmax_Line      = $nameEnergyWmax
    Wave Energy_lmax_Line      = $nameEnergyLmax

    W_max_Line           = NaN
    L_max_Line           = NaN
    betaExtra_min_Line   = NaN
    Energy_beta_min_Line = NaN
    Energy_wmax_Line     = NaN
    Energy_lmax_Line     = NaN

    // ---------- temp/cache folders ----------
    String chanFolder = "root:SNS_LineTmp"
    NewDataFolder/O $chanFolder

    String freePhaseFolder = ""
    if (usePhaseField && smodel == 0)
        freePhaseFolder = SNS_EnsureFreeVortexPhase2D(Nmask, xV_nm, yV_nm, params.SNS_nFlux, stiffnessFolder, chanFolder + ":FreeVortexPhase")
    endif

    String cacheName = CleanupName("LineDOSCache_" + tag, 0)
    String lineCacheFolder = chanFolder + ":" + cacheName

    KillDataFolder/Z $lineCacheFolder
    NewDataFolder/O $lineCacheFolder

    SetDataFolder $chanFolder

    Make/O/D/N=1 Bwin
    Bwin[0] = B0

    String nameDOS_EB    = "DOS_local_EB"
    String nameEaxisLoc  = "E_axis_local"
    String nameDOS_EBbr  = "DOS_local_EB_broad"

    Variable ir, iE
    Variable Nch
    Variable haveEaxis = 0

    // ---------- main loop ----------
    for (ir = 0; ir < Nr; ir += 1)

        Variable r0x = Xline[ir]
        Variable r0y = Yline[ir]

        SetDataFolder $chanFolder

        KillWaves/Z L_N_List_nm, W_eff_List_nm, wChan, T_eff_List
        KillWaves/Z Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm
        KillWaves/Z betaExtra_List
        KillWaves/Z L_N_List_m, W_eff_List_m

        SNS_BuildChannelsFromMask2D(Nmask, r0x, r0y, phiB, stepFac, chanFolder)

        Wave/Z L_N_List_nm   = L_N_List_nm
        Wave/Z W_eff_List_nm = W_eff_List_nm
        Wave/Z wChan      = wChan
        Wave/Z T_eff_List = T_eff_List

        Wave/Z Hit1x_List_nm = Hit1x_List_nm
        Wave/Z Hit1y_List_nm = Hit1y_List_nm
        Wave/Z Hit2x_List_nm = Hit2x_List_nm
        Wave/Z Hit2y_List_nm = Hit2y_List_nm

        if (!WaveExists(L_N_List_nm) || !WaveExists(W_eff_List_nm) || !WaveExists(wChan) || !WaveExists(T_eff_List) || !WaveExists(Hit1x_List_nm) || !WaveExists(Hit1y_List_nm) || !WaveExists(Hit2x_List_nm) || !WaveExists(Hit2y_List_nm))

            SetDataFolder dfrCaller
            DOS_Line[][ir]             = NaN
            DOS_Line_Conv[][ir]        = NaN
            W_max_Line[ir]             = NaN
            L_max_Line[ir]             = NaN
            betaExtra_min_Line[ir]     = NaN
            Energy_beta_min_Line[ir]   = NaN
            Energy_wmax_Line[ir]       = NaN
            Energy_lmax_Line[ir]       = NaN
            SetDataFolder $chanFolder
            continue
        endif

        Nch = DimSize(L_N_List_nm, 0)

        if ((Nch <= 0) || (Nch != DimSize(W_eff_List_nm,0)) || (Nch != DimSize(wChan,0)) || (Nch != DimSize(T_eff_List,0)) || (Nch != DimSize(Hit1x_List_nm,0)) || (Nch != DimSize(Hit1y_List_nm,0)) || (Nch != DimSize(Hit2x_List_nm,0)) || (Nch != DimSize(Hit2y_List_nm,0)))

            SetDataFolder dfrCaller
            DOS_Line[][ir]             = NaN
            DOS_Line_Conv[][ir]        = NaN
            W_max_Line[ir]             = NaN
            L_max_Line[ir]             = NaN
            betaExtra_min_Line[ir]     = NaN
            Energy_beta_min_Line[ir]   = NaN
            Energy_wmax_Line[ir]       = NaN
            Energy_lmax_Line[ir]       = NaN
            SetDataFolder $chanFolder
            continue
        endif

        Variable chWmax = -1
        Variable chLmax = -1
        Variable tmpMaxW = -Inf
        Variable tmpMaxL = -Inf
        Variable jScan

        for (jScan = 0; jScan < Nch; jScan += 1)
            if (!numtype(W_eff_List_nm[jScan]) && (W_eff_List_nm[jScan] > tmpMaxW))
                tmpMaxW = W_eff_List_nm[jScan]
                chWmax = jScan
            endif

            if (!numtype(L_N_List_nm[jScan]) && (L_N_List_nm[jScan] > tmpMaxL))
                tmpMaxL = L_N_List_nm[jScan]
                chLmax = jScan
            endif
        endfor

        SetDataFolder dfrCaller
        W_max_Line[ir] = tmpMaxW
        L_max_Line[ir] = tmpMaxL
        SetDataFolder $chanFolder

        if (usePhaseField)

            if (smodel == 0)
                Wave/Z PhaseReFree_line0 = $(freePhaseFolder + "PhaseReFree")
                Wave/Z PhaseImFree_line0 = $(freePhaseFolder + "PhaseImFree")
                if (!WaveExists(PhaseReFree_line0) || !WaveExists(PhaseImFree_line0))
                    Abort "SNS_LineDOS_BFixed_FromMask: free-vortex phase generation did not produce PhaseReFree/PhaseImFree."
                endif
                SNS_ComputeBetaExtraFromExteriorPhase2D(Nmask, PhaseReFree_line0, PhaseImFree_line0, Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, chanFolder)
            elseif (smodel == 1)
                Wave/Z PhaseReCorr_line1 = $(stiffnessFolder + ":PhaseReCorr")
                Wave/Z PhaseImCorr_line1 = $(stiffnessFolder + ":PhaseImCorr")
                if (!WaveExists(PhaseReCorr_line1) || !WaveExists(PhaseImCorr_line1))
                    Abort "SNS_LineDOS_BFixed_FromMask: ScreeningModel=1 requires PhaseReCorr/PhaseImCorr in StiffnessPhase."
                endif
                SNS_ComputeBetaExtraFromExteriorPhase2D(Nmask, PhaseReCorr_line1, PhaseImCorr_line1, Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, chanFolder)
            elseif (smodel == 2)
                Wave/Z QxCorr_line2 = $(stiffnessFolder + ":QxCorr")
                Wave/Z QyCorr_line2 = $(stiffnessFolder + ":QyCorr")
                Wave/Z PhaseReCorr_line2 = $(stiffnessFolder + ":PhaseReCorr")
                Wave/Z PhaseImCorr_line2 = $(stiffnessFolder + ":PhaseImCorr")
                if (!WaveExists(QxCorr_line2) || !WaveExists(QyCorr_line2) || !WaveExists(PhaseReCorr_line2) || !WaveExists(PhaseImCorr_line2))
                    Abort "SNS_LineDOS_BFixed_FromMask: ScreeningModel=2 requires QxCorr/QyCorr and PhaseReCorr/PhaseImCorr in StiffnessPhase."
                endif
                SNS_ComputeBetaExtraFromQField2D(QxCorr_line2, QyCorr_line2, Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, chanFolder, coreHandling=4, xV_nm=xV_nm, yV_nm=yV_nm, rCore_nm=15, PhaseReCore=PhaseReCorr_line2, PhaseImCore=PhaseImCorr_line2)
            elseif (smodel == 3)
                Wave/Z QxScreen_line3 = $(screenedFolder + ":QxScreen")
                Wave/Z QyScreen_line3 = $(screenedFolder + ":QyScreen")
                Wave/Z PhaseReCorr_line3 = $(stiffnessFolder + ":PhaseReCorr")
                Wave/Z PhaseImCorr_line3 = $(stiffnessFolder + ":PhaseImCorr")
                if (!WaveExists(QxScreen_line3) || !WaveExists(QyScreen_line3) || !WaveExists(PhaseReCorr_line3) || !WaveExists(PhaseImCorr_line3))
                    Abort "SNS_LineDOS_BFixed_FromMask: ScreeningModel=3 requires QxScreen/QyScreen and PhaseReCorr/PhaseImCorr."
                endif
                SNS_ComputeBetaExtraFromQField2D(QxScreen_line3, QyScreen_line3, Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, chanFolder, coreHandling=4, xV_nm=xV_nm, yV_nm=yV_nm, rCore_nm=15, PhaseReCore=PhaseReCorr_line3, PhaseImCore=PhaseImCorr_line3)
            elseif (smodel == 4)
                Wave/Z QxLondon_line4 = $(londonFolder + ":QxLondonCorr")
                Wave/Z QyLondon_line4 = $(londonFolder + ":QyLondonCorr")
                if (!WaveExists(QxLondon_line4) || !WaveExists(QyLondon_line4))
                    Abort "SNS_LineDOS_BFixed_FromMask: ScreeningModel=4 requires QxLondonCorr/QyLondonCorr in LocalLondonStiffnessPhase."
                endif
                SNS_ComputeBetaExtraFromQField2D(QxLondon_line4, QyLondon_line4, Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, chanFolder, coreHandling=0, xV_nm=xV_nm, yV_nm=yV_nm, rCore_nm=15)
            endif

            Wave/Z betaExtra_List = betaExtra_List

            if (!WaveExists(betaExtra_List))
                Abort "SNS_LineDOS_BFixed_FromMask: betaExtra_List was not created."
            endif

            if (DimSize(betaExtra_List,0) != Nch)
                Abort "SNS_LineDOS_BFixed_FromMask: betaExtra_List length mismatch."
            endif

            SNS_ComputeDOS_FromChannels(Bwin, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, xV_nm, yV_nm, N_cont_statesPer_eV, nameDOS_EB, nameEaxisLoc, betaExtra_List=betaExtra_List)

            Variable jCh, kBr
            Variable minAbsBeta = Inf
            Variable absBeta
            Variable bestBetaCh = -1

            for (jCh = 0; jCh < Nch; jCh += 1)
                absBeta = abs(betaExtra_List[jCh])
                if (!numtype(absBeta) && (absBeta < minAbsBeta))
                    minAbsBeta = absBeta
                    bestBetaCh = jCh
                endif
            endfor

            SetDataFolder dfrCaller
            if (bestBetaCh >= 0)
                betaExtra_min_Line[ir] = minAbsBeta
            else
                betaExtra_min_Line[ir] = NaN
            endif
            Energy_beta_min_Line[ir] = NaN
            SetDataFolder $chanFolder

            if (bestBetaCh >= 0)
                String chWavePathBeta
                sprintf chWavePathBeta, "%s:E_allBranches:ch%03d", chanFolder, bestBetaCh
                Wave/Z E_beta_ch = $chWavePathBeta

                if (!WaveExists(E_beta_ch))
                    sprintf chWavePathBeta, "%s:E_allBranches_ch%03d", chanFolder, bestBetaCh
                    Wave/Z E_beta_ch = $chWavePathBeta
                endif

                if (WaveExists(E_beta_ch))
                    Variable nBrBeta = DimSize(E_beta_ch, 1)
                    Variable minAbsE_beta = Inf
                    Variable absE_beta

                    for (kBr = 0; kBr < nBrBeta; kBr += 1)
                        absE_beta = abs(E_beta_ch[0][kBr])
                        if (!numtype(absE_beta) && (absE_beta < minAbsE_beta))
                            minAbsE_beta = absE_beta
                        endif
                    endfor

                    SetDataFolder dfrCaller
                    if (minAbsE_beta < Inf)
                        Energy_beta_min_Line[ir] = minAbsE_beta
                    else
                        Energy_beta_min_Line[ir] = NaN
                    endif
                    SetDataFolder $chanFolder
                endif
            endif

            Variable kBr2
            Variable minAbsE_W = Inf
            Variable minAbsE_L = Inf
            Variable absE2

            if (chWmax >= 0)
                String chWavePathW
                sprintf chWavePathW, "%s:E_allBranches:ch%03d", chanFolder, chWmax
                Wave/Z E_w_ch = $chWavePathW

                if (!WaveExists(E_w_ch))
                    sprintf chWavePathW, "%s:E_allBranches_ch%03d", chanFolder, chWmax
                    Wave/Z E_w_ch = $chWavePathW
                endif

                if (WaveExists(E_w_ch))
                    for (kBr2 = 0; kBr2 < DimSize(E_w_ch, 1); kBr2 += 1)
                        absE2 = abs(E_w_ch[0][kBr2])
                        if (!numtype(absE2) && (absE2 < minAbsE_W))
                            minAbsE_W = absE2
                        endif
                    endfor
                endif
            endif

            if (chLmax >= 0)
                String chWavePathL
                sprintf chWavePathL, "%s:E_allBranches:ch%03d", chanFolder, chLmax
                Wave/Z E_l_ch = $chWavePathL

                if (!WaveExists(E_l_ch))
                    sprintf chWavePathL, "%s:E_allBranches_ch%03d", chanFolder, chLmax
                    Wave/Z E_l_ch = $chWavePathL
                endif

                if (WaveExists(E_l_ch))
                    for (kBr2 = 0; kBr2 < DimSize(E_l_ch, 1); kBr2 += 1)
                        absE2 = abs(E_l_ch[0][kBr2])
                        if (!numtype(absE2) && (absE2 < minAbsE_L))
                            minAbsE_L = absE2
                        endif
                    endfor
                endif
            endif

            SetDataFolder dfrCaller
            Energy_wmax_Line[ir] = (minAbsE_W < Inf) ? minAbsE_W : NaN
            Energy_lmax_Line[ir] = (minAbsE_L < Inf) ? minAbsE_L : NaN
            SetDataFolder $chanFolder

        else

            SNS_ComputeDOS_FromChannels(Bwin, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, xV_nm, yV_nm, N_cont_statesPer_eV, nameDOS_EB, nameEaxisLoc)

            SetDataFolder dfrCaller
            betaExtra_min_Line[ir]   = NaN
            Energy_beta_min_Line[ir] = NaN
            SetDataFolder $chanFolder

            Variable kBrNoMask
            Variable minAbsE_W_NoMask = Inf
            Variable minAbsE_L_NoMask = Inf
            Variable absE_NoMask

            if (chWmax >= 0)
                String chWavePathW_NoMask
                sprintf chWavePathW_NoMask, "%s:E_allBranches:ch%03d", chanFolder, chWmax
                Wave/Z E_w_ch_NoMask = $chWavePathW_NoMask

                if (!WaveExists(E_w_ch_NoMask))
                    sprintf chWavePathW_NoMask, "%s:E_allBranches_ch%03d", chanFolder, chWmax
                    Wave/Z E_w_ch_NoMask = $chWavePathW_NoMask
                endif

                if (WaveExists(E_w_ch_NoMask))
                    for (kBrNoMask = 0; kBrNoMask < DimSize(E_w_ch_NoMask, 1); kBrNoMask += 1)
                        absE_NoMask = abs(E_w_ch_NoMask[0][kBrNoMask])
                        if (!numtype(absE_NoMask) && (absE_NoMask < minAbsE_W_NoMask))
                            minAbsE_W_NoMask = absE_NoMask
                        endif
                    endfor
                else
                    SNS_Log("E_allBranches wave for chWmax not found.", level="WARN")
                endif
            endif

            if (chLmax >= 0)
                String chWavePathL_NoMask
                sprintf chWavePathL_NoMask, "%s:E_allBranches:ch%03d", chanFolder, chLmax
                Wave/Z E_l_ch_NoMask = $chWavePathL_NoMask

                if (!WaveExists(E_l_ch_NoMask))
                    sprintf chWavePathL_NoMask, "%s:E_allBranches_ch%03d", chanFolder, chLmax
                    Wave/Z E_l_ch_NoMask = $chWavePathL_NoMask
                endif

                if (WaveExists(E_l_ch_NoMask))
                    for (kBrNoMask = 0; kBrNoMask < DimSize(E_l_ch_NoMask, 1); kBrNoMask += 1)
                        absE_NoMask = abs(E_l_ch_NoMask[0][kBrNoMask])
                        if (!numtype(absE_NoMask) && (absE_NoMask < minAbsE_L_NoMask))
                            minAbsE_L_NoMask = absE_NoMask
                        endif
                    endfor
                else
                    SNS_Log("E_allBranches wave for chLmax not found.", level="WARN")
                endif
            endif

            SetDataFolder dfrCaller
            Energy_wmax_Line[ir] = (minAbsE_W_NoMask < Inf) ? minAbsE_W_NoMask : NaN
            Energy_lmax_Line[ir] = (minAbsE_L_NoMask < Inf) ? minAbsE_L_NoMask : NaN
            SetDataFolder $chanFolder

        endif

        SNS_LineDOS_SaveTrajectoryState(chanFolder, lineCacheFolder, ir)

        Wave/Z DOS_EB_raw   = $nameDOS_EB
        Wave/Z E_axis_local = $nameEaxisLoc

        if (!WaveExists(DOS_EB_raw) || !WaveExists(E_axis_local))
            SetDataFolder dfrCaller
            DOS_Line[][ir]      = NaN
            DOS_Line_Conv[][ir] = NaN
            SetDataFolder $chanFolder
            continue
        endif

        SNS_ApplyDOS_Broadening_TplusMod(DOS_EB_raw, nameDOS_EBbr)
        Wave/Z DOS_EB_broad = $nameDOS_EBbr

        if (!WaveExists(DOS_EB_broad))
            SetDataFolder dfrCaller
            DOS_Line[][ir]      = NaN
            DOS_Line_Conv[][ir] = NaN
            SetDataFolder $chanFolder
            continue
        endif

        if (!haveEaxis)
            SetDataFolder dfrCaller
            Duplicate/O E_axis_local, $nameEaxis
            haveEaxis = 1
            SetDataFolder $chanFolder
        endif

        SetDataFolder dfrCaller
        for (iE = 0; iE < NE; iE += 1)
            DOS_Line[iE][ir]      = DOS_EB_raw[iE][0]
            DOS_Line_Conv[iE][ir] = DOS_EB_broad[iE][0]
        endfor
        SetDataFolder $chanFolder

    endfor

    // ---------- axes and metadata ----------
    SetDataFolder dfrCaller

    Wave/Z E_axis = $nameEaxis
    Wave/Z Rline2 = $nameRaxis

    if (!WaveExists(E_axis))
        Abort "SNS_LineDOS_BFixed_FromMask: no energy axis was created; all line positions failed."
    endif

    if (!WaveExists(Rline2))
        Abort "SNS_LineDOS_BFixed_FromMask: R axis missing."
    endif

    SetScale/P x, DimOffset(E_axis, 0), DimDelta(E_axis, 0), WaveUnits(E_axis, 0), DOS_Line
    SetScale/P x, DimOffset(E_axis, 0), DimDelta(E_axis, 0), WaveUnits(E_axis, 0), DOS_Line_Conv

    SetScale/P y, DimOffset(Rline2, 0), DimDelta(Rline2, 0), WaveUnits(Rline2, 0), DOS_Line
    SetScale/P y, DimOffset(Rline2, 0), DimDelta(Rline2, 0), WaveUnits(Rline2, 0), DOS_Line_Conv

    if (Nr > 1)
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", W_max_Line
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", L_max_Line
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", betaExtra_min_Line
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", Energy_beta_min_Line
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", Energy_wmax_Line
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", Energy_lmax_Line
    else
        SetScale/P x, 0, 1, "nm", W_max_Line
        SetScale/P x, 0, 1, "nm", L_max_Line
        SetScale/P x, 0, 1, "nm", betaExtra_min_Line
        SetScale/P x, 0, 1, "nm", Energy_beta_min_Line
        SetScale/P x, 0, 1, "nm", Energy_wmax_Line
        SetScale/P x, 0, 1, "nm", Energy_lmax_Line
    endif

    String meta
    meta  = "SNS_LineCacheFolder=" + SNS_LineDOS_TrailingColon(lineCacheFolder) + ";"
    meta += "SNS_B0_T=" + num2str(B0) + ";"
    meta += "SNS_phiB_rad=" + num2str(phiB) + ";"
    meta += "SNS_HasScreeningModel=" + num2str(hasScreeningModel) + ";"
    meta += "SNS_ScreeningModel=" + num2str(smodel) + ";"
    meta += "SNS_UsePhaseField=" + num2str(usePhaseField) + ";"
    meta += "SNS_LineStep_nm=" + num2str(lineStep_nm) + ";"
    meta += "SNS_LineStep_m=" + num2str(lineStep_m) + ";"

    if (!ParamIsDefault(NLinePts))
        meta += "SNS_NLinePts_requested=" + num2str(NLinePts) + ";"
    else
        meta += "SNS_NLinePts_requested=;"
    endif

    meta += "SNS_Raxis=" + GetWavesDataFolder(Rline2, 2) + ";"
    meta += "SNS_betaExtra_min=" + GetWavesDataFolder(betaExtra_min_Line, 2) + ";"
    meta += "SNS_Energy_betaExtra_min=" + GetWavesDataFolder(Energy_beta_min_Line, 2) + ";"
    meta += "SNS_Energy_wmax=" + GetWavesDataFolder(Energy_wmax_Line, 2) + ";"
    meta += "SNS_Energy_Lmax=" + GetWavesDataFolder(Energy_lmax_Line, 2) + ";"

    if (usePhaseField)
        meta += "SNS_MaskFolder=" + maskBase + ";"
        meta += "SNS_StiffnessPhaseFolder=" + stiffnessFolder + ";"
        meta += "SNS_ScreenedStiffnessPhaseFolder=" + screenedFolder + ";"
        meta += "SNS_LocalLondonStiffnessPhaseFolder=" + londonFolder + ";"
    else
        meta += "SNS_MaskFolder=;"
        meta += "SNS_StiffnessPhaseFolder=;"
        meta += "SNS_ScreenedStiffnessPhaseFolder=;"
        meta += "SNS_LocalLondonStiffnessPhaseFolder=;"
    endif
    meta += "SNS_FreePhaseFolder=" + freePhaseFolder + ";"

    Note/K DOS_Line
    Note DOS_Line, meta

    Note/K DOS_Line_Conv
    Note DOS_Line_Conv, meta

    Note/K betaExtra_min_Line
    Note betaExtra_min_Line, meta

    Note/K Energy_beta_min_Line
    Note Energy_beta_min_Line, meta

    Note/K Energy_wmax_Line
    Note Energy_wmax_Line, meta

    Note/K Energy_lmax_Line
    Note Energy_lmax_Line, meta

    SetDataFolder dfrCaller
    return Nr
End

//==============================================================================
// SNS_LineDOS_BFixed_FromMask3D
//
// Purpose:
//   Compute a 1D line-resolved 3D LDOS/DOS profile at fixed magnetic field B0:
//
//       DOS(E, r) and experimentally broadened DOS_Conv(E, r)
//
//   For each position r along the line, the function:
//     1. builds 3D SNS channels using the Cu(111)-faceted tracer,
//     2. optionally adds betaExtra_List from a vortex / phase-field / Q-field,
//     3. computes DOS(E,B0) from the channel ensemble,
//     4. applies thermal + lock-in modulation broadening,
//     5. stores line-position diagnostics for trajectory inspection.
//
// 3D-specific behavior:
//   - Uses SNS_BuildChannelsFromMask3D(...) instead of 2D ray tracing.
//   - Requires H_nm and maxPath_nm as positional inputs.
//   - Uses DOS3D_eV_Vol * area * H_nm for the continuum DOS.
//   - L_N_List and W_eff_List are returned in meters and passed directly to
//     SNS_ComputeDOS_FromChannels(...).
//   - W_max_* and L_max_* outputs are converted to nm.
//   - For ScreeningModel = 2/3/4, betaExtra_List is generated inside the
//     3D channel builder from QxPhase/QyPhase line integration.
//   - For endpoint phase modes 0/1, betaExtra_List is computed after channel
//     construction using fullMask_loc, matching the existing 3D convention.
//
// Interface transparency:
//   This function no longer accepts BTK_barrier / Zbarrier.
//   The 3D channel builder/tracer reads the barrier from SNS_Settings via
//   SNS_Params and applies the 3D bottom-interface incidence formula.
//
// Special B=0 behavior:
//   If B0 == 0 and ScreeningModel is omitted, no phase-field / vortex correction
//   is used. The function skips phase-field folder/wave checks and uses simple
//   output names:
//
//       LDOS_3D_B0mT
//       LDOS_Conv_3D_B0mT
//
//   If ScreeningModel is explicitly supplied, the selected model is used even at
//   B0 == 0. In that case maskFolder is required and the output tag keeps the
//   angle/model suffix.
//
// Inputs:
//   Nmask      : 2D mask defining the N footprint; axes scaled in nm.
//   phiB       : in-plane magnetic-field angle [rad]
//   xStart,
//   yStart     : line start point [nm]
//   xEnd,
//   yEnd       : line end point [nm]
//   B0         : fixed magnetic field [T]
//   stepFac    : angular/ray sampling step factor
//   H_nm       : 3D normal-region height [nm]
//   maxPath_nm : maximum path / phase-coherence cutoff [nm]
//
// Optional Inputs:
//   maskFolder     : optional parent folder containing phase/Q-field
//                    subfolders. Required for ScreeningModel 1--4. For model 0,
//                    compatible PhaseReFree/PhaseImFree waves are reused; when
//                    absent, the analytical free-vortex phase is generated in a
//                    temporary calculation folder.
//
//   ScreeningModel : phase / Q-field selector:
//
//                      0 = free vortex phase
//                      1 = stiffness-corrected vortex phase
//                      2 = corrected phase-gradient Q
//                      3 = screened corrected Q
//                      4 = local-London corrected Q
//
//                    If omitted and B0 == 0, no phase model is used.
//                    If supplied, the selected model is used even at B0 == 0.
//
//   NLinePts       : optional number of line positions.
//                    If supplied, overrides the default LambdaF-based line step.
//                    Must be finite integer >= 2.
//
// Outputs:
//   LDOS_*                  : raw DOS(E,r) at B = B0 [states/eV]
//   LDOS_Conv_*             : broadened DOS(E,r) [states/eV]
//   E_axis_*                : energy axis [eV]
//   R_axis_*                : distance along line [nm]
//   W_max_*                 : maximum effective channel width vs r [nm]
//   L_max_*                 : maximum trajectory length vs r [nm]
//   betaExtra_min_*         : min(abs(betaExtra_List)) vs r [rad]
//   Energy_betaExtra_min_*  : min(abs(E_ABS)) for betaExtra_min channel [eV]
//   Energy_maxW_*           : min(abs(E_ABS)) for W_max channel [eV]
//   Energy_Lmax_*           : min(abs(E_ABS)) for L_max channel [eV]
//
// Returns:
//   Nr : number of positions along the line
//
// Notes:
//   This is a 1D spatial fixed-field workflow and belongs in
//   SNS_SpatialMaps.ipf.
//==============================================================================
Function SNS_LineDOS_BFixed_FromMask3D(Nmask, phiB, xStart, yStart, xEnd, yEnd, B0, stepFac, H_nm, maxPath_nm, [maskFolder, ScreeningModel, NLinePts])

    Wave Nmask
    Variable phiB
    Variable xStart, yStart, xEnd, yEnd
    Variable B0
    Variable stepFac
    Variable H_nm, maxPath_nm
    String maskFolder
    Variable ScreeningModel
    Variable NLinePts

    DFREF dfrCaller = GetDataFolderDFR()

    Variable isZeroField = (abs(B0) <= 1e-15)
    Variable hasScreeningModel = !ParamIsDefault(ScreeningModel)

    Variable smodel = 0
    if (hasScreeningModel)
        smodel = ScreeningModel
    endif
    smodel = round(smodel)

    if ((smodel < 0) || (smodel > 4))
        Abort "SNS_LineDOS_BFixed_FromMask3D: ScreeningModel must be 0, 1, 2, 3, or 4."
    endif

    if (hasScreeningModel && smodel != 0 && ParamIsDefault(maskFolder))
        SetDataFolder dfrCaller
        Abort "SNS_LineDOS_BFixed_FromMask3D: maskFolder is required for ScreeningModel 1, 2, 3, or 4."
    endif

    Variable usePhaseField = (hasScreeningModel && smodel == 0) || ((!ParamIsDefault(maskFolder)) && ((!isZeroField) || hasScreeningModel))

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    if (numtype(H_nm) != 0 || H_nm <= 0)
        SetDataFolder dfrCaller
        Abort "SNS_LineDOS_BFixed_FromMask3D: H_nm must be finite and positive."
    endif

    if (numtype(maxPath_nm) != 0 || maxPath_nm <= 0)
        SetDataFolder dfrCaller
        Abort "SNS_LineDOS_BFixed_FromMask3D: maxPath_nm must be finite and positive."
    endif

    Variable lambdaF_nm = params.LambdaF
    Variable NE         = params.NE

    // ---------- geometry folder ----------
    DFREF dfrGeom = GetWavesDataFolderDFR(Nmask)
    SetDataFolder dfrGeom

    Wave/Z Vortex_ptx = Vortex_ptx
    Wave/Z Vortex_pty = Vortex_pty

    Variable xV_nm = 0
    Variable yV_nm = 0
    if (WaveExists(Vortex_ptx) && WaveExists(Vortex_pty))
        xV_nm = Vortex_ptx[0]
        yV_nm = Vortex_pty[0]
    endif

    Wave/Z w_area_nm2 = SNS_MaskAreaPerim_FromParticles(Nmask)

    Variable N_cont_statesPer_eV = 0
    if (WaveExists(w_area_nm2))
        Variable area_nm2 = w_area_nm2[0]
        Variable area_m2  = area_nm2 * 1e-18
        N_cont_statesPer_eV = params.DOS3D_eV_Vol * area_m2 * H_nm * 1e-9
    endif

    SetDataFolder dfrCaller

    // ---------- model subfolders ----------
    String maskBase = ""
    String stiffnessFolder = ""
    String screenedFolder = ""
    String londonFolder = ""

    if (!ParamIsDefault(maskFolder))
        maskBase = SNS_LineDOS_TrailingColon(maskFolder)
        stiffnessFolder = maskBase + "StiffnessPhase"
        screenedFolder  = maskBase + "ScreenedStiffnessPhase"
        londonFolder    = maskBase + "LocalLondonStiffnessPhase"
    endif

    // ---------- output names ----------
    String bTag = "B" + num2str(round(B0*1e3)) + "mT"
    String angleTag = num2str(round(phiB*180/pi)) + "deg"
    String tag

    if (isZeroField && !hasScreeningModel)
        tag = "3D_B0mT"
    elseif (smodel == 1)
        tag = "3D_" + bTag + "_" + angleTag + "_PhaseCorr"
    elseif (smodel == 2)
        tag = "3D_" + bTag + "_" + angleTag + "_QCorr"
    elseif (smodel == 3)
        tag = "3D_" + bTag + "_" + angleTag + "_ScreenedCorr"
    elseif (smodel == 4)
        tag = "3D_" + bTag + "_" + angleTag + "_LondonCorr"
    else
        tag = "3D_" + bTag + "_" + angleTag + "_PhaseFree"
    endif

    String nameDOSLine       = "LDOS_"                 + tag
    String nameDOSLineConv   = "LDOS_Conv_"            + tag
    String nameEaxis         = "E_axis_"               + tag
    String nameRaxis         = "R_axis_"               + tag
    String nameWmax          = "W_max_"                + tag
    String nameLmax          = "L_max_"                + tag
    String nameBetaExtraMin  = "betaExtra_min_"        + tag
    String nameEnergyBetaMin = "Energy_betaExtra_min_" + tag
    String nameEnergyMaxW    = "Energy_maxW_"          + tag
    String nameEnergyLmax    = "Energy_Lmax_"          + tag

    // ---------- line positions ----------
    String nameX = "SNS_Line_X"
    String nameY = "SNS_Line_Y"

    Variable lineStep_nm = lambdaF_nm

    if (!ParamIsDefault(NLinePts))
        NLinePts = round(NLinePts)

        if ((numtype(NLinePts) != 0) || (NLinePts < 2))
            SetDataFolder dfrCaller
            Abort "SNS_LineDOS_BFixed_FromMask3D: NLinePts must be a finite integer >= 2."
        endif

        Variable lineLength_nm = sqrt((xEnd - xStart)^2 + (yEnd - yStart)^2)

        if ((numtype(lineLength_nm) != 0) || (lineLength_nm <= 0))
            SetDataFolder dfrCaller
            Abort "SNS_LineDOS_BFixed_FromMask3D: cannot use NLinePts for a zero-length or invalid line."
        endif

        lineStep_nm = lineLength_nm / NLinePts
    endif

    Variable lineStep_m = lineStep_nm * 1e-9

    Variable Nr = SNS_BuildLinePositions(xStart, yStart, xEnd, yEnd, lineStep_m, nameX, nameY, nameRaxis)

    Wave Xline = $nameX
    Wave Yline = $nameY
    Wave Rline = $nameRaxis

    if (Nr > 1)
        Variable drnm = (Rline[Nr-1] - Rline[0]) / (Nr - 1)
        SetScale/P x, Rline[0], drnm, "nm", Rline
    else
        SetScale/P x, 0, 1, "nm", Rline
    endif

    // ---------- allocate outputs ----------
    Make/O/D/N=(NE, Nr) $nameDOSLine
    Make/O/D/N=(NE, Nr) $nameDOSLineConv
    Make/O/D/N=(Nr)     $nameWmax
    Make/O/D/N=(Nr)     $nameLmax
    Make/O/D/N=(Nr)     $nameBetaExtraMin
    Make/O/D/N=(Nr)     $nameEnergyBetaMin
    Make/O/D/N=(Nr)     $nameEnergyMaxW
    Make/O/D/N=(Nr)     $nameEnergyLmax

    Wave DOS_Line             = $nameDOSLine
    Wave DOS_Line_Conv        = $nameDOSLineConv
    Wave W_max_Line           = $nameWmax
    Wave L_max_Line           = $nameLmax
    Wave betaExtra_min_Line   = $nameBetaExtraMin
    Wave Energy_beta_min_Line = $nameEnergyBetaMin
    Wave Energy_maxW_Line     = $nameEnergyMaxW
    Wave Energy_lmax_Line     = $nameEnergyLmax

    DOS_Line             = NaN
    DOS_Line_Conv        = NaN
    W_max_Line           = NaN
    L_max_Line           = NaN
    betaExtra_min_Line   = NaN
    Energy_beta_min_Line = NaN
    Energy_maxW_Line     = NaN
    Energy_lmax_Line     = NaN

    // ---------- temp/cache folders ----------
    String chanFolder = "root:SNS_LineTmp3D"
    NewDataFolder/O $chanFolder

    String freePhaseFolder = ""
    if (usePhaseField && smodel == 0)
        freePhaseFolder = SNS_EnsureFreeVortexPhase2D(Nmask, xV_nm, yV_nm, params.SNS_nFlux, stiffnessFolder, chanFolder + ":FreeVortexPhase")
    endif

    String cacheName = CleanupName("LineDOSCache_" + tag, 0)
    String lineCacheFolder = chanFolder + ":" + cacheName

    KillDataFolder/Z $lineCacheFolder
    NewDataFolder/O $lineCacheFolder

    SetDataFolder $chanFolder

    Make/O/D/N=1 Bwin
    Bwin[0] = B0

    String nameDOS_EB   = "DOS_local_EB"
    String nameEaxisLoc = "E_axis_local"
    String nameDOS_EBbr = "DOS_local_EB_broad"

    Variable ir, iE, Nch
    Variable haveEaxis = 0

    // ---------- main loop ----------
    for (ir = 0; ir < Nr; ir += 1)

        Variable r0x = Xline[ir]
        Variable r0y = Yline[ir]

        SetDataFolder $chanFolder

        KillWaves/Z L_N_List_nm, W_eff_List_nm, wChan, T_eff_List
        KillWaves/Z Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm
        KillWaves/Z betaExtra_List

        // ---------- build channels ----------
        if (usePhaseField && smodel == 2)

            Wave/Z QxCorr_loc = $(stiffnessFolder + ":QxCorr")
            Wave/Z QyCorr_loc = $(stiffnessFolder + ":QyCorr")
            Wave/Z PhaseReCorr_loc = $(stiffnessFolder + ":PhaseReCorr")
            Wave/Z PhaseImCorr_loc = $(stiffnessFolder + ":PhaseImCorr")

            if (!WaveExists(QxCorr_loc) || !WaveExists(QyCorr_loc) || !WaveExists(PhaseReCorr_loc) || !WaveExists(PhaseImCorr_loc))
                Abort "SNS_LineDOS_BFixed_FromMask3D: StiffnessPhase folder missing QxCorr/QyCorr or PhaseReCorr/PhaseImCorr."
            endif

            SNS_BuildChannelsFromMask3D(Nmask, r0x, r0y, phiB, stepFac, H_nm, maxPath_nm, chanFolder, QxPhase=QxCorr_loc, QyPhase=QyCorr_loc, qNstep=1, xV_nm=xV_nm, yV_nm=yV_nm, rCore_nm=15, PhaseReCore=PhaseReCorr_loc, PhaseImCore=PhaseImCorr_loc)

        elseif (usePhaseField && smodel == 3)

            Wave/Z QxScreen_loc = $(screenedFolder + ":QxScreen")
            Wave/Z QyScreen_loc = $(screenedFolder + ":QyScreen")
            Wave/Z PhaseReCorr_loc = $(stiffnessFolder + ":PhaseReCorr")
            Wave/Z PhaseImCorr_loc = $(stiffnessFolder + ":PhaseImCorr")

            if (!WaveExists(QxScreen_loc) || !WaveExists(QyScreen_loc))
                Abort "SNS_LineDOS_BFixed_FromMask3D: ScreenedStiffnessPhase folder missing QxScreen/QyScreen."
            endif
            if (!WaveExists(PhaseReCorr_loc) || !WaveExists(PhaseImCorr_loc))
                Abort "SNS_LineDOS_BFixed_FromMask3D: StiffnessPhase folder missing PhaseReCorr/PhaseImCorr required for model 3 core handling."
            endif

            SNS_BuildChannelsFromMask3D(Nmask, r0x, r0y, phiB, stepFac, H_nm, maxPath_nm, chanFolder, QxPhase=QxScreen_loc, QyPhase=QyScreen_loc, qNstep=1, xV_nm=xV_nm, yV_nm=yV_nm, rCore_nm=15, PhaseReCore=PhaseReCorr_loc, PhaseImCore=PhaseImCorr_loc)

        elseif (usePhaseField && smodel == 4)

            Wave/Z QxLondon_loc = $(londonFolder + ":QxLondonCorr")
            Wave/Z QyLondon_loc = $(londonFolder + ":QyLondonCorr")

            if (!WaveExists(QxLondon_loc) || !WaveExists(QyLondon_loc))
                Abort "SNS_LineDOS_BFixed_FromMask3D: LocalLondonStiffnessPhase folder missing QxLondonCorr/QyLondonCorr."
            endif

            SNS_BuildChannelsFromMask3D(Nmask, r0x, r0y, phiB, stepFac, H_nm, maxPath_nm, chanFolder, QxPhase=QxLondon_loc, QyPhase=QyLondon_loc, qNstep=1)

        else

            SNS_BuildChannelsFromMask3D(Nmask, r0x, r0y, phiB, stepFac, H_nm, maxPath_nm, chanFolder)

        endif

        Wave/Z L_N_List_nm   = L_N_List_nm
        Wave/Z W_eff_List_nm = W_eff_List_nm
        Wave/Z wChan         = wChan
        Wave/Z T_eff_List    = T_eff_List

        Wave/Z Hit1x_List_nm = Hit1x_List_nm
        Wave/Z Hit1y_List_nm = Hit1y_List_nm
        Wave/Z Hit2x_List_nm = Hit2x_List_nm
        Wave/Z Hit2y_List_nm = Hit2y_List_nm

        if (!WaveExists(L_N_List_nm) || !WaveExists(W_eff_List_nm) || !WaveExists(wChan) || !WaveExists(T_eff_List) || !WaveExists(Hit1x_List_nm) || !WaveExists(Hit1y_List_nm) || !WaveExists(Hit2x_List_nm) || !WaveExists(Hit2y_List_nm))

            SetDataFolder dfrCaller
            DOS_Line[][ir]             = NaN
            DOS_Line_Conv[][ir]        = NaN
            W_max_Line[ir]             = NaN
            L_max_Line[ir]             = NaN
            betaExtra_min_Line[ir]     = NaN
            Energy_beta_min_Line[ir]   = NaN
            Energy_maxW_Line[ir]       = NaN
            Energy_lmax_Line[ir]       = NaN
            SetDataFolder $chanFolder
            continue
        endif

        Nch = numpnts(L_N_List_nm)

        if ((Nch <= 0) || (Nch != numpnts(W_eff_List_nm)) || (Nch != numpnts(wChan)) || (Nch != numpnts(T_eff_List)) || (Nch != numpnts(Hit1x_List_nm)) || (Nch != numpnts(Hit1y_List_nm)) || (Nch != numpnts(Hit2x_List_nm)) || (Nch != numpnts(Hit2y_List_nm)))

            SetDataFolder dfrCaller
            DOS_Line[][ir]             = NaN
            DOS_Line_Conv[][ir]        = NaN
            W_max_Line[ir]             = NaN
            L_max_Line[ir]             = NaN
            betaExtra_min_Line[ir]     = NaN
            Energy_beta_min_Line[ir]   = NaN
            Energy_maxW_Line[ir]       = NaN
            Energy_lmax_Line[ir]       = NaN
            SetDataFolder $chanFolder
            continue
        endif

        // ---------- identify max-W and max-L channels ----------
        Variable chWmax = -1
        Variable chLmax = -1
        Variable tmpMaxW = -Inf
        Variable tmpMaxL = -Inf
        Variable jScan

        for (jScan = 0; jScan < Nch; jScan += 1)
            if (!numtype(W_eff_List_nm[jScan]) && (W_eff_List_nm[jScan] > tmpMaxW))
                tmpMaxW = W_eff_List_nm[jScan]
                chWmax = jScan
            endif

            if (!numtype(L_N_List_nm[jScan]) && (L_N_List_nm[jScan] > tmpMaxL))
                tmpMaxL = L_N_List_nm[jScan]
                chLmax = jScan
            endif
        endfor

        SetDataFolder dfrCaller
        W_max_Line[ir] = tmpMaxW
        L_max_Line[ir] = tmpMaxL
        SetDataFolder $chanFolder

        // ---------- betaExtra handling ----------
        if (usePhaseField)

            if (smodel == 2 || smodel == 3 || smodel == 4)

                Wave/Z betaExtra_List = betaExtra_List

                if (!WaveExists(betaExtra_List))
                    Abort "SNS_LineDOS_BFixed_FromMask3D: Q-mode requested but betaExtra_List was not created by SNS_BuildChannelsFromMask3D."
                endif

                if (numpnts(betaExtra_List) != Nch)
                    Abort "SNS_LineDOS_BFixed_FromMask3D: betaExtra_List length mismatch in Q-mode."
                endif

            else

                KillWaves/Z betaExtra_List

                Duplicate/O Nmask fullMask_loc
                fullMask_loc = 0

                if (smodel == 1)
                    Wave/Z PhaseReCorr_line3D = $(stiffnessFolder + ":PhaseReCorr")
                    Wave/Z PhaseImCorr_line3D = $(stiffnessFolder + ":PhaseImCorr")
                    if (!WaveExists(PhaseReCorr_line3D) || !WaveExists(PhaseImCorr_line3D))
                        Abort "SNS_LineDOS_BFixed_FromMask3D: ScreeningModel=1 requires PhaseReCorr/PhaseImCorr in StiffnessPhase."
                    endif
                    SNS_ComputeBetaExtraFromExteriorPhase2D(fullMask_loc, PhaseReCorr_line3D, PhaseImCorr_line3D, Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, chanFolder)
                else
                    Wave/Z PhaseReFree_line3D = $(freePhaseFolder + "PhaseReFree")
                    Wave/Z PhaseImFree_line3D = $(freePhaseFolder + "PhaseImFree")
                    if (!WaveExists(PhaseReFree_line3D) || !WaveExists(PhaseImFree_line3D))
                        Abort "SNS_LineDOS_BFixed_FromMask3D: free-vortex phase generation did not produce PhaseReFree/PhaseImFree."
                    endif
                    SNS_ComputeBetaExtraFromExteriorPhase2D(fullMask_loc, PhaseReFree_line3D, PhaseImFree_line3D, Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, chanFolder)
                endif

                Wave/Z betaExtra_List = betaExtra_List

                if (!WaveExists(betaExtra_List))
                    Abort "SNS_LineDOS_BFixed_FromMask3D: endpoint betaExtra_List was not created."
                endif

                if (numpnts(betaExtra_List) != Nch)
                    Abort "SNS_LineDOS_BFixed_FromMask3D: endpoint betaExtra_List length mismatch."
                endif
            endif

            SNS_ComputeDOS_FromChannels(Bwin, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, xV_nm, yV_nm, N_cont_statesPer_eV, nameDOS_EB, nameEaxisLoc, betaExtra_List=betaExtra_List)

            // ---------- betaExtra_min and corresponding ABS energy ----------
            Variable jCh, kBr
            Variable minAbsBeta = Inf
            Variable absBeta
            Variable bestBetaCh = -1

            for (jCh = 0; jCh < Nch; jCh += 1)
                absBeta = abs(betaExtra_List[jCh])
                if (!numtype(absBeta) && (absBeta < minAbsBeta))
                    minAbsBeta = absBeta
                    bestBetaCh = jCh
                endif
            endfor

            SetDataFolder dfrCaller
            if (bestBetaCh >= 0)
                betaExtra_min_Line[ir] = minAbsBeta
            else
                betaExtra_min_Line[ir] = NaN
            endif
            Energy_beta_min_Line[ir] = NaN
            SetDataFolder $chanFolder

            if (bestBetaCh >= 0)
                String chWavePathBeta
                sprintf chWavePathBeta, "%s:E_allBranches:ch%03d", chanFolder, bestBetaCh
                Wave/Z E_beta_ch = $chWavePathBeta

                if (!WaveExists(E_beta_ch))
                    sprintf chWavePathBeta, "%s:E_allBranches_ch%03d", chanFolder, bestBetaCh
                    Wave/Z E_beta_ch = $chWavePathBeta
                endif

                if (WaveExists(E_beta_ch))
                    Variable minAbsE_beta = Inf
                    Variable absE_beta

                    for (kBr = 0; kBr < DimSize(E_beta_ch, 1); kBr += 1)
                        absE_beta = abs(E_beta_ch[0][kBr])
                        if (!numtype(absE_beta) && (absE_beta < minAbsE_beta))
                            minAbsE_beta = absE_beta
                        endif
                    endfor

                    SetDataFolder dfrCaller
                    if (minAbsE_beta < Inf)
                        Energy_beta_min_Line[ir] = minAbsE_beta
                    else
                        Energy_beta_min_Line[ir] = NaN
                    endif
                    SetDataFolder $chanFolder
                endif
            endif

            // ---------- Energy corresponding to W_max and L_max channels ----------
            Variable kBr2
            Variable minAbsE_W = Inf
            Variable minAbsE_L = Inf
            Variable absE2

            if (chWmax >= 0)
                String chWavePathW
                sprintf chWavePathW, "%s:E_allBranches:ch%03d", chanFolder, chWmax
                Wave/Z E_w_ch = $chWavePathW

                if (!WaveExists(E_w_ch))
                    sprintf chWavePathW, "%s:E_allBranches_ch%03d", chanFolder, chWmax
                    Wave/Z E_w_ch = $chWavePathW
                endif

                if (WaveExists(E_w_ch))
                    for (kBr2 = 0; kBr2 < DimSize(E_w_ch, 1); kBr2 += 1)
                        absE2 = abs(E_w_ch[0][kBr2])
                        if (!numtype(absE2) && (absE2 < minAbsE_W))
                            minAbsE_W = absE2
                        endif
                    endfor
                endif
            endif

            if (chLmax >= 0)
                String chWavePathL
                sprintf chWavePathL, "%s:E_allBranches:ch%03d", chanFolder, chLmax
                Wave/Z E_l_ch = $chWavePathL

                if (!WaveExists(E_l_ch))
                    sprintf chWavePathL, "%s:E_allBranches_ch%03d", chanFolder, chLmax
                    Wave/Z E_l_ch = $chWavePathL
                endif

                if (WaveExists(E_l_ch))
                    for (kBr2 = 0; kBr2 < DimSize(E_l_ch, 1); kBr2 += 1)
                        absE2 = abs(E_l_ch[0][kBr2])
                        if (!numtype(absE2) && (absE2 < minAbsE_L))
                            minAbsE_L = absE2
                        endif
                    endfor
                endif
            endif

            SetDataFolder dfrCaller
            if (minAbsE_W < Inf)
                Energy_maxW_Line[ir] = minAbsE_W
            else
                Energy_maxW_Line[ir] = NaN
            endif

            if (minAbsE_L < Inf)
                Energy_lmax_Line[ir] = minAbsE_L
            else
                Energy_lmax_Line[ir] = NaN
            endif
            SetDataFolder $chanFolder

        else

            SNS_ComputeDOS_FromChannels(Bwin, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, xV_nm, yV_nm, N_cont_statesPer_eV, nameDOS_EB, nameEaxisLoc)

            SetDataFolder dfrCaller
            betaExtra_min_Line[ir]   = NaN
            Energy_beta_min_Line[ir] = NaN
            Energy_maxW_Line[ir]     = NaN
            Energy_lmax_Line[ir]     = NaN
            SetDataFolder $chanFolder

        endif

        // ---------- save trajectory state ----------
        SNS_LineDOS_SaveTrajectoryState(chanFolder, lineCacheFolder, ir)

        Wave/Z DOS_EB_raw   = $nameDOS_EB
        Wave/Z E_axis_local = $nameEaxisLoc

        if (!WaveExists(DOS_EB_raw) || !WaveExists(E_axis_local))
            SetDataFolder dfrCaller
            DOS_Line[][ir]      = NaN
            DOS_Line_Conv[][ir] = NaN
            SetDataFolder $chanFolder
            continue
        endif

        SNS_ApplyDOS_Broadening_TplusMod(DOS_EB_raw, nameDOS_EBbr)

        Wave/Z DOS_EB_broad = $nameDOS_EBbr
        if (!WaveExists(DOS_EB_broad))
            SetDataFolder dfrCaller
            DOS_Line[][ir]      = NaN
            DOS_Line_Conv[][ir] = NaN
            SetDataFolder $chanFolder
            continue
        endif

        if (!haveEaxis)
            SetDataFolder dfrCaller
            Duplicate/O E_axis_local, $nameEaxis
            haveEaxis = 1
            SetDataFolder $chanFolder
        endif

        SetDataFolder dfrCaller
        for (iE = 0; iE < NE; iE += 1)
            DOS_Line[iE][ir]      = DOS_EB_raw[iE][0]
            DOS_Line_Conv[iE][ir] = DOS_EB_broad[iE][0]
        endfor
        SetDataFolder $chanFolder

    endfor

    // ---------- axes and metadata ----------
    SetDataFolder dfrCaller

    Wave/Z E_axis = $nameEaxis
    Wave/Z Rline2 = $nameRaxis

    if (!WaveExists(E_axis))
        Abort "SNS_LineDOS_BFixed_FromMask3D: no energy axis was created; all line positions failed."
    endif
    if (!WaveExists(Rline2))
        Abort "SNS_LineDOS_BFixed_FromMask3D: R axis missing."
    endif

    String eUnits = WaveUnits(E_axis, 0)
    String rUnits = WaveUnits(Rline2, 0)

    SetScale/P x, DimOffset(E_axis, 0), DimDelta(E_axis, 0), eUnits, DOS_Line
    SetScale/P y, DimOffset(Rline2, 0), DimDelta(Rline2, 0), rUnits, DOS_Line

    SetScale/P x, DimOffset(E_axis, 0), DimDelta(E_axis, 0), eUnits, DOS_Line_Conv
    SetScale/P y, DimOffset(Rline2, 0), DimDelta(Rline2, 0), rUnits, DOS_Line_Conv

    if (Nr > 1)
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", W_max_Line
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", L_max_Line
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", betaExtra_min_Line
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", Energy_beta_min_Line
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", Energy_maxW_Line
        SetScale/P x, Rline2[0], Rline2[1]-Rline2[0], "nm", Energy_lmax_Line
    else
        SetScale/P x, 0, 1, "nm", W_max_Line
        SetScale/P x, 0, 1, "nm", L_max_Line
        SetScale/P x, 0, 1, "nm", betaExtra_min_Line
        SetScale/P x, 0, 1, "nm", Energy_beta_min_Line
        SetScale/P x, 0, 1, "nm", Energy_maxW_Line
        SetScale/P x, 0, 1, "nm", Energy_lmax_Line
    endif

    String meta
    meta  = "SNS_LineCacheFolder=" + SNS_LineDOS_TrailingColon(lineCacheFolder) + ";"
    meta += "SNS_ChannelBuilder=3D;"
    meta += "SNS_B0_T=" + num2str(B0) + ";"
    meta += "SNS_phiB_rad=" + num2str(phiB) + ";"
    meta += "SNS_H_nm=" + num2str(H_nm) + ";"
    meta += "SNS_maxPath_nm=" + num2str(maxPath_nm) + ";"
    meta += "SNS_HasScreeningModel=" + num2str(hasScreeningModel) + ";"
    meta += "SNS_ScreeningModel=" + num2str(smodel) + ";"
    meta += "SNS_UsePhaseField=" + num2str(usePhaseField) + ";"
    meta += "SNS_LineStep_nm=" + num2str(lineStep_nm) + ";"
    meta += "SNS_LineStep_m=" + num2str(lineStep_m) + ";"

    if (!ParamIsDefault(NLinePts))
        meta += "SNS_NLinePts_requested=" + num2str(NLinePts) + ";"
    else
        meta += "SNS_NLinePts_requested=;"
    endif

    meta += "SNS_Raxis=" + GetWavesDataFolder(Rline2, 2) + ";"
    meta += "SNS_betaExtra_min=" + GetWavesDataFolder(betaExtra_min_Line, 2) + ";"
    meta += "SNS_Energy_betaExtra_min=" + GetWavesDataFolder(Energy_beta_min_Line, 2) + ";"
    meta += "SNS_Energy_maxW=" + GetWavesDataFolder(Energy_maxW_Line, 2) + ";"
    meta += "SNS_Energy_Lmax=" + GetWavesDataFolder(Energy_lmax_Line, 2) + ";"

    if (usePhaseField)
        meta += "SNS_MaskFolder=" + maskBase + ";"
        meta += "SNS_StiffnessPhaseFolder=" + stiffnessFolder + ";"
        meta += "SNS_ScreenedStiffnessPhaseFolder=" + screenedFolder + ";"
        meta += "SNS_LocalLondonStiffnessPhaseFolder=" + londonFolder + ";"
    else
        meta += "SNS_MaskFolder=;"
        meta += "SNS_StiffnessPhaseFolder=;"
        meta += "SNS_ScreenedStiffnessPhaseFolder=;"
        meta += "SNS_LocalLondonStiffnessPhaseFolder=;"
    endif
    meta += "SNS_FreePhaseFolder=" + freePhaseFolder + ";"

    Note/K DOS_Line
    Note DOS_Line, meta

    Note/K DOS_Line_Conv
    Note DOS_Line_Conv, meta

    Note/K betaExtra_min_Line
    Note betaExtra_min_Line, meta

    Note/K Energy_beta_min_Line
    Note Energy_beta_min_Line, meta

    Note/K Energy_maxW_Line
    Note Energy_maxW_Line, meta

    Note/K Energy_lmax_Line
    Note Energy_lmax_Line, meta

    SetDataFolder dfrCaller
    return Nr
End

//==============================================================================
// SPATIAL SAMPLING HELPERS
//
// Generic helpers for testing or sampling spatial coordinates against scaled
// mask waves.
//==============================================================================

//==============================================================================
// SNS_SampleMaskNearest
//
// Purpose:
//   Sample a 2D scaled mask at physical coordinates (x_nm, y_nm) using nearest-
//   pixel lookup.
//
// Inputs:
//   Nmask : 2D mask wave with x/y scaling in nm.
//   x_nm  : x coordinate [nm].
//   y_nm  : y coordinate [nm].
//
// Returns:
//   Mask value at the nearest pixel.
//   Returns 0 if the coordinate lies outside the wave bounds.
//
// Notes:
//   Generic spatial-map helper.
//   Does not change the current data folder.
//==============================================================================
Function SNS_SampleMaskNearest(Nmask, x_nm, y_nm)

    Wave Nmask
    Variable x_nm, y_nm

    Variable ix = round((x_nm - DimOffset(Nmask, 0)) / DimDelta(Nmask, 0))
    Variable iy = round((y_nm - DimOffset(Nmask, 1)) / DimDelta(Nmask, 1))

    if ((ix < 0) || (ix >= DimSize(Nmask, 0)) || (iy < 0) || (iy >= DimSize(Nmask, 1)))
        return 0
    endif

    return Nmask[ix][iy]
End

//==============================================================================
// FIXED-FIELD POINT-LIST LDOS MAPS
//
// 2D and 3D fixed-B map workflows. These evaluate local DOS at explicit point
// lists and store scalar energy-window-integrated LDOS values per point.
//==============================================================================

//==============================================================================
// SNS_LDOSmap_BFixed_FromMask
//
// Purpose:
//   Compute a 2D point-list LDOS map at fixed magnetic field B0.
//
//   For each point in Pos_nm, the function:
//     1. checks whether the point lies inside Nmask,
//     2. builds a local 2D SNS channel ensemble,
//     3. optionally adds betaExtra_List from a vortex / phase-field / Q-field,
//     4. computes DOS(E,B0),
//     5. applies thermal + lock-in modulation broadening,
//     6. integrates DOS over a symmetric energy window around zero.
//
//   In addition to scalar integrated LDOS values, this function stores full
//   local spectra as LDOS(E, point) waves. Use these together with XY_* to
//   build XYS triplets at a chosen energy:
//       x = XY_*[][0]
//       y = XY_*[][1]
//       s = LDOS_E_raw_*[iE][] or LDOS_E_conv_*[iE][]
//   It does not persist the full trajectory cache.
//
// Inputs:
//   Nmask   : 2D binary N-region mask.
//             Inside N: > 0.5, outside: <= 0.5.
//             Wave axes must be scaled in nm.
//
//   Pos_nm  : 2D point list with columns:
//               Pos_nm[][0] = x [nm]
//               Pos_nm[][1] = y [nm]
//
//   phiB    : in-plane magnetic-field angle [rad]
//   B0      : fixed magnetic field [T]
//   stepFac : ray/channel sampling refinement factor
//
// Optional Inputs:
//   Eint_meV       : half-width of symmetric integration window [meV].
//                   Default: 0.1 meV.
//                   The integrated window is [-Eint_meV, +Eint_meV].
//
//   maskFolder     : optional parent folder containing phase/Q-field
//                   subfolders. Required for ScreeningModel 1--4. For model 0,
//                   supplied PhaseReFree/PhaseImFree waves are reused; when
//                   absent, the analytical free-vortex phase is generated in a
//                   temporary calculation folder.
//
//   ScreeningModel : phase / Q-field selector:
//
//                      0 = free vortex phase
//                      1 = stiffness-corrected vortex phase
//                      2 = corrected phase-gradient Q
//                      3 = screened corrected Q
//                      4 = local-London corrected Q
//
// Interface transparency:
//   This function no longer accepts BTK_barrier / Zbarrier.
//   Interface transparency is controlled globally through SNS_Settings and
//   applied inside SNS_BuildChannelsFromMask2D.
//
// Optional model-0 reuse location:
//   maskFolder:StiffnessPhase:
//       PhaseReFree, PhaseImFree
//
// Required phase-folder layout for ScreeningModel 1--4:
//   maskFolder:StiffnessPhase:
//       PhaseReCorr, PhaseImCorr
//       QxCorr, QyCorr
//
//   maskFolder:ScreenedStiffnessPhase:
//       QxScreen, QyScreen
//
//   maskFolder:LocalLondonStiffnessPhase:
//       QxLondonCorr, QyLondonCorr
//
// Outputs:
//   LDOS_Int_*       : raw integrated LDOS per input point
//   LDOS_Int_Conv_*  : broadened integrated LDOS per input point
//   W_max_*          : maximum W_eff per input point
//   L_max_*          : maximum L_N per input point
//   ValidPt_*        : 1 for valid computed point, 0 otherwise
//   XY_*             : copy of input point list
//   Ewin_*           : integration window [eV]
//   LDOS_E_raw_*     : raw local spectra, [energy][point]
//   LDOS_E_conv_*    : broadened local spectra, [energy][point]
//
// Returns:
//   Nvalid           : number of valid computed points.
//
// Notes:
//   Model 4 uses QxLondonCorr/QyLondonCorr directly with coreHandling=0.
//   This is a spatial fixed-field map workflow and belongs in SNS_SpatialMaps.ipf.
//==============================================================================
Function SNS_LDOSmap_BFixed_FromMask(Nmask, Pos_nm, phiB, B0, stepFac, [Eint_meV, maskFolder, ScreeningModel])

    Wave   Nmask
    Wave   Pos_nm
    Variable phiB
    Variable B0
    Variable stepFac
    Variable Eint_meV
    String maskFolder
    Variable ScreeningModel

    DFREF dfrCaller = GetDataFolderDFR()
    Variable runtimeTimerRef = StartMSTimer
    Variable runtimeStartDateTime = DateTime

    Variable hasScreeningModel = !ParamIsDefault(ScreeningModel)
    Variable smodel = 0
    if (hasScreeningModel)
        smodel = round(ScreeningModel)
    endif
	if (hasScreeningModel && smodel != 0 && ParamIsDefault(maskFolder))
	    SetDataFolder dfrCaller
	    Abort "SNS_LDOSmap_BFixed_FromMask: maskFolder is required for ScreeningModel 1, 2, 3, or 4."
	endif

    if ((smodel < 0) || (smodel > 4))
        Abort "SNS_LDOSmap_BFixed_FromMask: ScreeningModel must be 0, 1, 2, 3, or 4."
    endif

    if (ParamIsDefault(Eint_meV))
        Eint_meV = 0.1
    endif
    Variable Eint_eV = Eint_meV * 1e-3

    if (DimSize(Pos_nm, 1) < 2)
        Abort "SNS_LDOSmap_BFixed_FromMask: Pos_nm must be a 2D wave with columns [x_nm, y_nm]."
    endif

    Variable Np = DimSize(Pos_nm, 0)
    if (Np <= 0)
        Abort "SNS_LDOSmap_BFixed_FromMask: Pos_nm contains no points."
    endif
    SNS_Log("SNS_LDOSmap_BFixed_FromMask: start; input points=" + num2istr(Np), level="INFO")

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    Variable NE = params.NE

    DFREF dfrGeom = GetWavesDataFolderDFR(Nmask)
    SetDataFolder dfrGeom

    Wave/Z Vortex_ptx = Vortex_ptx
    Wave/Z Vortex_pty = Vortex_pty

    Variable xV_nm, yV_nm
    if (WaveExists(Vortex_ptx) && WaveExists(Vortex_pty))
        xV_nm = Vortex_ptx[0]
        yV_nm = Vortex_pty[0]
    else
        xV_nm = 0
        yV_nm = 0
    endif

    Wave/Z w_area_nm2 = SNS_MaskAreaPerim_FromParticles(Nmask)

    Variable N_cont_statesPer_eV
    if (WaveExists(w_area_nm2))
        Variable area_nm2 = w_area_nm2[0]
        Variable area_m2  = area_nm2 * 1e-18
        N_cont_statesPer_eV = params.DOS2D_eV_Area * area_m2
    else
        N_cont_statesPer_eV = 0
    endif

    SetDataFolder dfrCaller

    String maskBase = ""
    String stiffnessFolder = ""
    String screenedFolder = ""
    String londonFolder = ""

    if (!ParamIsDefault(maskFolder))
        maskBase = SNS_LineDOS_TrailingColon(maskFolder)
        stiffnessFolder = maskBase + "StiffnessPhase"
        screenedFolder  = maskBase + "ScreenedStiffnessPhase"
        londonFolder    = maskBase + "LocalLondonStiffnessPhase"
    endif

    String tag
    if (smodel == 1)
        tag = "B" + num2str(round(B0*1e3)) + "mT_" + num2str(round(phiB*180/pi)) + "deg_PhaseCorr"
    elseif (smodel == 2)
        tag = "B" + num2str(round(B0*1e3)) + "mT_" + num2str(round(phiB*180/pi)) + "deg_QCorr"
    elseif (smodel == 3)
        tag = "B" + num2str(round(B0*1e3)) + "mT_" + num2str(round(phiB*180/pi)) + "deg_ScreenedCorr"
    elseif (smodel == 4)
        tag = "B" + num2str(round(B0*1e3)) + "mT_" + num2str(round(phiB*180/pi)) + "deg_LondonCorr"
    else
        tag = "B" + num2str(round(B0*1e3)) + "mT_" + num2str(round(phiB*180/pi)) + "deg_PhaseFree"
    endif

    String nameLDOSIntRaw  = "LDOS_Int_"      + tag
    String nameLDOSIntConv = "LDOS_Int_Conv_" + tag
    String nameWmax        = "W_max_"         + tag
    String nameLmax        = "L_max_"         + tag
    String nameValid       = "ValidPt_"       + tag
    String nameXY          = "XY_"            + tag
    String nameEwin        = "Ewin_"          + tag
    String nameLDOSEraw    = "LDOS_E_raw_"    + tag
    String nameLDOSEconv   = "LDOS_E_conv_"   + tag

    Make/O/D/N=(Np) $nameLDOSIntRaw
    Make/O/D/N=(Np) $nameLDOSIntConv
    Make/O/D/N=(Np) $nameWmax
    Make/O/D/N=(Np) $nameLmax
    Make/O/B/N=(Np) $nameValid
    Duplicate/O Pos_nm, $nameXY
    Make/O/D/N=2 $nameEwin

    Wave LDOS_Int      = $nameLDOSIntRaw
    Wave LDOS_Int_Conv = $nameLDOSIntConv
    Wave W_max_List    = $nameWmax
    Wave L_max_List    = $nameLmax
    Wave ValidPt       = $nameValid
    Wave XY_out        = $nameXY
    Wave Ewin          = $nameEwin

    LDOS_Int      = NaN
    LDOS_Int_Conv = NaN
    W_max_List    = NaN
    L_max_List    = NaN
    ValidPt       = 0

    Ewin[0] = -Eint_eV
    Ewin[1] = +Eint_eV

    SetScale/P x, 0, 1, "", LDOS_Int
    SetScale/P x, 0, 1, "", LDOS_Int_Conv
    SetScale/P x, 0, 1, "", W_max_List
    SetScale/P x, 0, 1, "", L_max_List
    SetScale/P x, 0, 1, "", ValidPt

    String chanFolder = "root:SNS_MapTmp"
    NewDataFolder/O $chanFolder
    Variable useFreePhase = (smodel == 0) && (hasScreeningModel || !ParamIsDefault(maskFolder))
    String freePhaseFolder = ""
    if (useFreePhase)
        freePhaseFolder = SNS_EnsureFreeVortexPhase2D(Nmask, xV_nm, yV_nm, params.SNS_nFlux, stiffnessFolder, chanFolder + ":FreeVortexPhase")
    endif
    SetDataFolder $chanFolder

    Make/O/D/N=1 Bwin
    Bwin[0] = B0

    String nameDOS_EB   = "DOS_local_EB"
    String nameEaxisLoc = "E_axis_local"
    String nameDOS_EBbr = "DOS_local_EB_broad"

    Variable ip, iE
    Variable x0_nm, y0_nm
    Variable Nch
    Variable iLo, iHi, dE_eV
    Variable sumRaw, sumConv
    Variable Nvalid = 0
    Variable ldosInitialized = 0

    for (ip = 0; ip < Np; ip += 1)

        x0_nm = Pos_nm[ip][0]
        y0_nm = Pos_nm[ip][1]

        if (SNS_SampleMaskNearest(Nmask, x0_nm, y0_nm) <= 0.5)
            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        SetDataFolder $chanFolder

        KillWaves/Z L_N_List_nm, W_eff_List_nm, wChan, T_eff_List
        KillWaves/Z Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm
        KillWaves/Z betaExtra_List
        SNS_BuildChannelsFromMask2D(Nmask, x0_nm, y0_nm, phiB, stepFac, chanFolder)

        Wave/Z L_N_List_nm   = L_N_List_nm
        Wave/Z W_eff_List_nm = W_eff_List_nm
        Wave/Z wChan      = wChan
        Wave/Z T_eff_List = T_eff_List

        Wave/Z Hit1x_List_nm = Hit1x_List_nm
        Wave/Z Hit1y_List_nm = Hit1y_List_nm
        Wave/Z Hit2x_List_nm = Hit2x_List_nm
        Wave/Z Hit2y_List_nm = Hit2y_List_nm

        if (!WaveExists(L_N_List_nm) || !WaveExists(W_eff_List_nm) || !WaveExists(wChan) || !WaveExists(T_eff_List) || \
            !WaveExists(Hit1x_List_nm) || !WaveExists(Hit1y_List_nm) || !WaveExists(Hit2x_List_nm) || !WaveExists(Hit2y_List_nm))

            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        Nch = DimSize(L_N_List_nm, 0)

        if ((Nch <= 0) || \
            (Nch != DimSize(W_eff_List_nm,0)) || \
            (Nch != DimSize(wChan,0)) || \
            (Nch != DimSize(T_eff_List,0)) || \
            (Nch != DimSize(Hit1x_List_nm,0)) || \
            (Nch != DimSize(Hit1y_List_nm,0)) || \
            (Nch != DimSize(Hit2x_List_nm,0)) || \
            (Nch != DimSize(Hit2y_List_nm,0)))

            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        SetDataFolder dfrCaller
        W_max_List[ip] = WaveMax(W_eff_List_nm)
        L_max_List[ip] = WaveMax(L_N_List_nm)
        SetDataFolder $chanFolder

        if (useFreePhase || !ParamIsDefault(maskFolder))

            if (smodel == 0)
                Wave/Z PhaseReFree_loc = $(freePhaseFolder + "PhaseReFree")
                Wave/Z PhaseImFree_loc = $(freePhaseFolder + "PhaseImFree")
                if (!WaveExists(PhaseReFree_loc) || !WaveExists(PhaseImFree_loc))
                    Abort "SNS_LDOSmap_BFixed_FromMask: free-vortex phase generation did not produce PhaseReFree/PhaseImFree."
                endif
                SNS_ComputeBetaExtraFromExteriorPhase2D( \
                    Nmask, PhaseReFree_loc, PhaseImFree_loc, \
                    Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
                    chanFolder)
            elseif (smodel == 1)
                Wave/Z PhaseReCorr_model1 = $(stiffnessFolder + ":PhaseReCorr")
                Wave/Z PhaseImCorr_model1 = $(stiffnessFolder + ":PhaseImCorr")
                if (!WaveExists(PhaseReCorr_model1) || !WaveExists(PhaseImCorr_model1))
                    Abort "SNS_LDOSmap_BFixed_FromMask: ScreeningModel=1 requires PhaseReCorr/PhaseImCorr in StiffnessPhase."
                endif
                SNS_ComputeBetaExtraFromExteriorPhase2D( \
                    Nmask, PhaseReCorr_model1, PhaseImCorr_model1, \
                    Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
                    chanFolder)
            elseif (smodel == 2)
                Wave/Z QxCorr_model2 = $(stiffnessFolder + ":QxCorr")
                Wave/Z QyCorr_model2 = $(stiffnessFolder + ":QyCorr")
                Wave/Z PhaseReCorr_model2 = $(stiffnessFolder + ":PhaseReCorr")
                Wave/Z PhaseImCorr_model2 = $(stiffnessFolder + ":PhaseImCorr")
                if (!WaveExists(QxCorr_model2) || !WaveExists(QyCorr_model2) || !WaveExists(PhaseReCorr_model2) || !WaveExists(PhaseImCorr_model2))
                    Abort "SNS_LDOSmap_BFixed_FromMask: ScreeningModel=2 requires QxCorr/QyCorr and PhaseReCorr/PhaseImCorr in StiffnessPhase."
                endif
                SNS_ComputeBetaExtraFromQField2D( \
                    QxCorr_model2, QyCorr_model2, \
                    Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
                    chanFolder, coreHandling=4, xV_nm=xV_nm, yV_nm=yV_nm, rCore_nm=15, PhaseReCore=PhaseReCorr_model2, PhaseImCore=PhaseImCorr_model2, useInterp=1)
            elseif (smodel == 3)
                Wave/Z QxScreen_model3 = $(screenedFolder + ":QxScreen")
                Wave/Z QyScreen_model3 = $(screenedFolder + ":QyScreen")
                Wave/Z PhaseReCorr_model3 = $(stiffnessFolder + ":PhaseReCorr")
                Wave/Z PhaseImCorr_model3 = $(stiffnessFolder + ":PhaseImCorr")
                if (!WaveExists(QxScreen_model3) || !WaveExists(QyScreen_model3) || !WaveExists(PhaseReCorr_model3) || !WaveExists(PhaseImCorr_model3))
                    Abort "SNS_LDOSmap_BFixed_FromMask: ScreeningModel=3 requires QxScreen/QyScreen and PhaseReCorr/PhaseImCorr."
                endif
                SNS_ComputeBetaExtraFromQField2D( \
                    QxScreen_model3, QyScreen_model3, \
                    Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
                    chanFolder, coreHandling=4, xV_nm=xV_nm, yV_nm=yV_nm, rCore_nm=15, PhaseReCore=PhaseReCorr_model3, PhaseImCore=PhaseImCorr_model3, useInterp=1)
            elseif (smodel == 4)
                Wave/Z QxLondon_model4 = $(londonFolder + ":QxLondonCorr")
                Wave/Z QyLondon_model4 = $(londonFolder + ":QyLondonCorr")
                if (!WaveExists(QxLondon_model4) || !WaveExists(QyLondon_model4))
                    Abort "SNS_LDOSmap_BFixed_FromMask: ScreeningModel=4 requires QxLondonCorr/QyLondonCorr in LocalLondonStiffnessPhase."
                endif
                SNS_ComputeBetaExtraFromQField2D( \
                    QxLondon_model4, QyLondon_model4, \
                    Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
                    chanFolder, coreHandling=0, useInterp=1)
            endif

            Wave/Z betaExtra_List = betaExtra_List

            if (!WaveExists(betaExtra_List))
                Abort "SNS_LDOSmap_BFixed_FromMask: betaExtra_List was not created."
            endif

            if (DimSize(betaExtra_List,0) != Nch)
                Abort "SNS_LDOSmap_BFixed_FromMask: betaExtra_List length mismatch."
            endif

            SNS_ComputeDOS_FromChannels( \
                Bwin, \
                L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
                Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
                xV_nm, yV_nm, \
                N_cont_statesPer_eV, \
                nameDOS_EB, nameEaxisLoc, \
                betaExtra_List=betaExtra_List)

        else

            SNS_ComputeDOS_FromChannels( \
                Bwin, \
                L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
                Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
                xV_nm, yV_nm, \
                N_cont_statesPer_eV, \
                nameDOS_EB, nameEaxisLoc)

        endif

        Wave/Z DOS_EB_raw   = $nameDOS_EB
        Wave/Z E_axis_local = $nameEaxisLoc

        if (!WaveExists(DOS_EB_raw) || !WaveExists(E_axis_local))
            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        SNS_ApplyDOS_Broadening_TplusMod(DOS_EB_raw, nameDOS_EBbr)
        Wave/Z DOS_EB_broad = $nameDOS_EBbr

        if (!WaveExists(DOS_EB_broad))
            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        if (!ldosInitialized)
            SetDataFolder dfrCaller
            Make/O/D/N=(NE, Np) $nameLDOSEraw
            Make/O/D/N=(NE, Np) $nameLDOSEconv
            Wave LDOS_E_Raw_Init = $nameLDOSEraw
            Wave LDOS_E_Conv_Init = $nameLDOSEconv
            LDOS_E_Raw_Init = NaN
            LDOS_E_Conv_Init = NaN

            SetScale/P x, E_axis_local[0], E_axis_local[1]-E_axis_local[0], "eV", LDOS_E_Raw_Init
            SetScale/P x, E_axis_local[0], E_axis_local[1]-E_axis_local[0], "eV", LDOS_E_Conv_Init
            SetScale/P y, 0, 1, "point", LDOS_E_Raw_Init
            SetScale/P y, 0, 1, "point", LDOS_E_Conv_Init

            ldosInitialized = 1
            SetDataFolder $chanFolder
        endif

        SetDataFolder dfrCaller
        Wave LDOS_E_Raw_Fill = $nameLDOSEraw
        Wave LDOS_E_Conv_Fill = $nameLDOSEconv
        for (iE = 0; iE < NE; iE += 1)
            LDOS_E_Raw_Fill[iE][ip] = DOS_EB_raw[iE][0]
            LDOS_E_Conv_Fill[iE][ip] = DOS_EB_broad[iE][0]
        endfor
        SetDataFolder $chanFolder

        iLo = floor(x2pnt(E_axis_local, -Eint_eV))
        iHi = ceil( x2pnt(E_axis_local, +Eint_eV))

        iLo = max(0, iLo)
        iHi = min(NE-1, iHi)

        if (iHi < iLo)
            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        dE_eV = abs(DimDelta(E_axis_local, 0))

        sumRaw  = 0
        sumConv = 0
        for (iE = iLo; iE <= iHi; iE += 1)
            sumRaw  += DOS_EB_raw[iE][0]
            sumConv += DOS_EB_broad[iE][0]
        endfor

        SetDataFolder dfrCaller
        LDOS_Int[ip]      = sumRaw  * dE_eV
        LDOS_Int_Conv[ip] = sumConv * dE_eV
        ValidPt[ip]       = 1
        SetDataFolder $chanFolder

        Nvalid += 1

    endfor

    SetDataFolder dfrCaller

    Variable runtime_s = StopMSTimer(runtimeTimerRef)/1e6
    Variable runtimeEndDateTime = DateTime
    Variable/G v_SNS_LDOSmap_Runtime_s = runtime_s
    Variable/G v_SNS_LDOSmap_SecondsPerInput = runtime_s/Np
    Variable/G v_SNS_LDOSmap_NInput = Np
    Variable/G v_SNS_LDOSmap_NValid = Nvalid
    Variable/G v_SNS_LDOSmap_StartDateTime = runtimeStartDateTime
    Variable/G v_SNS_LDOSmap_EndDateTime = runtimeEndDateTime
    SNS_Log("SNS_LDOSmap_BFixed_FromMask: finished; input points=" + num2istr(Np) + "; valid points=" + num2istr(Nvalid) + "; runtime_s=" + num2str(runtime_s) + "; seconds_per_input=" + num2str(runtime_s/Np), level="INFO")

    String meta
    meta  = "SNS_ChannelBuilder=2D;"
    meta += "SNS_B0_T=" + num2str(B0) + ";"
    meta += "SNS_phiB_rad=" + num2str(phiB) + ";"
    meta += "SNS_ScreeningModel=" + num2str(smodel) + ";"
    meta += "SNS_Eint_meV=" + num2str(Eint_meV) + ";"
    meta += "SNS_Runtime_s=" + num2str(runtime_s) + ";"
    meta += "SNS_SecondsPerInput=" + num2str(runtime_s/Np) + ";"
    meta += "SNS_NInput=" + num2istr(Np) + ";"
    meta += "SNS_NValid=" + num2istr(Nvalid) + ";"
    meta += "SNS_LDOSDims=energy,point;"
    meta += "SNS_LDOSPositionWave=" + nameXY + ";"

    if (!ParamIsDefault(maskFolder))
        meta += "SNS_MaskFolder=" + maskBase + ";"
        meta += "SNS_StiffnessPhaseFolder=" + stiffnessFolder + ";"
        meta += "SNS_ScreenedStiffnessPhaseFolder=" + screenedFolder + ";"
        meta += "SNS_LocalLondonStiffnessPhaseFolder=" + londonFolder + ";"
    else
        meta += "SNS_MaskFolder=;"
        meta += "SNS_StiffnessPhaseFolder=;"
        meta += "SNS_ScreenedStiffnessPhaseFolder=;"
        meta += "SNS_LocalLondonStiffnessPhaseFolder=;"
    endif
    meta += "SNS_FreePhaseFolder=" + freePhaseFolder + ";"

    Note/K LDOS_Int
    Note LDOS_Int, meta

    Note/K LDOS_Int_Conv
    Note LDOS_Int_Conv, meta

    Wave/Z LDOS_E_Raw_Note = $nameLDOSEraw
    if (WaveExists(LDOS_E_Raw_Note))
        Note/K LDOS_E_Raw_Note
        Note LDOS_E_Raw_Note, meta
    endif

    Wave/Z LDOS_E_Conv_Note = $nameLDOSEconv
    if (WaveExists(LDOS_E_Conv_Note))
        Note/K LDOS_E_Conv_Note
        Note LDOS_E_Conv_Note, meta
    endif

    return Nvalid
End


//==============================================================================
// SNS_LDOSmap_BFixed_FromMask3D
//
// Purpose:
//   Compute a 3D point-list LDOS map at fixed magnetic field B0.
//
//   For each point in Pos_nm, the function:
//     1. checks whether the point lies inside Nmask,
//     2. builds a local 3D SNS channel ensemble using the Cu(111)-faceted tracer,
//     3. optionally adds betaExtra_List from a vortex / phase-field / Q-field,
//     4. computes DOS(E,B0),
//     5. applies thermal + lock-in modulation broadening,
//     6. integrates DOS over a symmetric energy window around zero.
//
//   Unlike the line-DOS workflow, this function stores only scalar integrated
//   LDOS values per input point. It does not persist the full trajectory cache.
//
// 3D-specific behavior:
//   - Uses SNS_BuildChannelsFromMask3D(...) instead of 2D ray tracing.
//   - Requires H_nm and maxPath_nm as positional inputs.
//   - Uses DOS3D_eV_Vol * area * H_nm for the continuum DOS.
//   - L_N_List and W_eff_List are returned in meters and passed directly to
//     SNS_ComputeDOS_FromChannels(...).
//   - W_max_* and L_max_* outputs are converted to nm.
//   - For ScreeningModel = 2/3/4, betaExtra_List is generated inside the
//     3D channel builder from QxPhase/QyPhase line integration.
//   - For endpoint phase modes 0/1, betaExtra_List is computed after channel
//     construction using fullMask_loc, matching the existing 3D convention.
//
// Interface transparency:
//   This function no longer accepts BTK_barrier / Zbarrier.
//   The 3D channel builder/tracer reads the barrier from SNS_Settings via
//   SNS_Params and applies the 3D bottom-interface incidence formula.
//
// Inputs:
//   Nmask      : 2D mask defining the N footprint; axes scaled in nm.
//   Pos_nm     : 2D point list with columns:
//                  Pos_nm[][0] = x [nm]
//                  Pos_nm[][1] = y [nm]
//   phiB       : in-plane magnetic-field angle [rad]
//   B0         : fixed magnetic field [T]
//   stepFac    : ray/channel sampling refinement factor
//   H_nm       : 3D normal-region height [nm]
//   maxPath_nm : maximum path / phase-coherence cutoff [nm]
//
// Optional Inputs:
//   Eint_meV       : half-width of symmetric integration window [meV].
//                   Default: 0.1 meV.
//                   The integrated window is [-Eint_meV, +Eint_meV].
//
//   maskFolder     : parent folder containing phase/Q-field subfolders.
//                   Required for ScreeningModel 1--4. For model 0, compatible
//                   free-phase waves are reused when present and otherwise
//                   generated analytically in a temporary calculation folder.
//
//   ScreeningModel : phase / Q-field selector:
//
//                      0 = free vortex phase
//                      1 = stiffness-corrected vortex phase
//                      2 = corrected phase-gradient Q
//                      3 = screened corrected Q
//                      4 = local-London corrected Q
//
// Optional model-0 reuse location:
//   maskFolder:StiffnessPhase:
//       PhaseReFree, PhaseImFree
//
// Required phase-folder layout for ScreeningModel 1--4:
//   maskFolder:StiffnessPhase:
//       PhaseReCorr, PhaseImCorr
//       QxCorr, QyCorr
//
//   maskFolder:ScreenedStiffnessPhase:
//       QxScreen, QyScreen
//
//   maskFolder:LocalLondonStiffnessPhase:
//       QxLondonCorr, QyLondonCorr
//
// Outputs:
//   LDOS_Int_*       : raw integrated LDOS per input point
//   LDOS_Int_Conv_*  : broadened integrated LDOS per input point
//   W_max_*          : maximum W_eff per input point [nm]
//   L_max_*          : maximum L_N per input point [nm]
//   ValidPt_*        : 1 for valid computed point, 0 otherwise
//   XY_*             : copy of input point list
//   Ewin_*           : integration window [eV]
//
// Returns:
//   Nvalid           : number of valid computed points.
//
// Notes:
//   Model 4 passes QxLondonCorr/QyLondonCorr directly to the 3D builder,
//   without phase-core fill arguments. This corresponds to using the London-Q
//   field as computed, analogous to coreHandling=0.
//   This is a spatial fixed-field map workflow and belongs in SNS_SpatialMaps.ipf.
//==============================================================================
Function SNS_LDOSmap_BFixed_FromMask3D(Nmask, Pos_nm, phiB, B0, stepFac, H_nm, maxPath_nm, [Eint_meV, maskFolder, ScreeningModel])

    Wave   Nmask
    Wave   Pos_nm
    Variable phiB
    Variable B0
    Variable stepFac
    Variable H_nm, maxPath_nm
    Variable Eint_meV
    String maskFolder
    Variable ScreeningModel

    DFREF dfrCaller = GetDataFolderDFR()

    Variable hasScreeningModel = !ParamIsDefault(ScreeningModel)
    Variable smodel = 0
    if (hasScreeningModel)
        smodel = ScreeningModel
    endif
    smodel = round(smodel)

    if (hasScreeningModel && smodel != 0 && ParamIsDefault(maskFolder))
    	SetDataFolder dfrCaller
		Abort "SNS_LDOSmap_BFixed_FromMask3D: maskFolder is required for ScreeningModel 1, 2, 3, or 4."
	endif

    if ((smodel < 0) || (smodel > 4))
        Abort "SNS_LDOSmap_BFixed_FromMask3D: ScreeningModel must be 0, 1, 2, 3, or 4."
    endif

    if (ParamIsDefault(Eint_meV))
        Eint_meV = 0.1
    endif
    Variable Eint_eV = Eint_meV * 1e-3

    if (DimSize(Pos_nm, 1) < 2)
        Abort "SNS_LDOSmap_BFixed_FromMask3D: Pos_nm must be a 2D wave with columns [x_nm, y_nm]."
    endif

    Variable Np = DimSize(Pos_nm, 0)
    if (Np <= 0)
        Abort "SNS_LDOSmap_BFixed_FromMask3D: Pos_nm contains no points."
    endif

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    Variable NE = params.NE

    DFREF dfrGeom = GetWavesDataFolderDFR(Nmask)
    SetDataFolder dfrGeom

    Wave/Z Vortex_ptx = Vortex_ptx
    Wave/Z Vortex_pty = Vortex_pty

    Variable xV_nm, yV_nm
    if (WaveExists(Vortex_ptx) && WaveExists(Vortex_pty))
        xV_nm = Vortex_ptx[0]
        yV_nm = Vortex_pty[0]
    else
        xV_nm = 0
        yV_nm = 0
    endif

    Wave/Z w_area_nm2 = SNS_MaskAreaPerim_FromParticles(Nmask)

    Variable N_cont_statesPer_eV
    if (WaveExists(w_area_nm2))
        Variable area_nm2 = w_area_nm2[0]
        Variable area_m2  = area_nm2 * 1e-18
        N_cont_statesPer_eV = params.DOS3D_eV_Vol * area_m2 * H_nm * 1e-9
    else
        N_cont_statesPer_eV = 0
    endif

    SetDataFolder dfrCaller

    String maskBase = ""
    String stiffnessFolder = ""
    String screenedFolder = ""
    String londonFolder = ""

    if (!ParamIsDefault(maskFolder))
        maskBase = SNS_LineDOS_TrailingColon(maskFolder)
        stiffnessFolder = maskBase + "StiffnessPhase"
        screenedFolder  = maskBase + "ScreenedStiffnessPhase"
        londonFolder    = maskBase + "LocalLondonStiffnessPhase"
    endif

    String tag
    if (smodel == 1)
        tag = "3D_B" + num2str(round(B0*1e3)) + "mT_" + num2str(round(phiB*180/pi)) + "deg_PhaseCorr"
    elseif (smodel == 2)
        tag = "3D_B" + num2str(round(B0*1e3)) + "mT_" + num2str(round(phiB*180/pi)) + "deg_QCorr"
    elseif (smodel == 3)
        tag = "3D_B" + num2str(round(B0*1e3)) + "mT_" + num2str(round(phiB*180/pi)) + "deg_ScreenedCorr"
    elseif (smodel == 4)
        tag = "3D_B" + num2str(round(B0*1e3)) + "mT_" + num2str(round(phiB*180/pi)) + "deg_LondonCorr"
    else
        tag = "3D_B" + num2str(round(B0*1e3)) + "mT_" + num2str(round(phiB*180/pi)) + "deg_PhaseFree"
    endif

    String nameLDOSIntRaw  = "LDOS_Int_"      + tag
    String nameLDOSIntConv = "LDOS_Int_Conv_" + tag
    String nameWmax        = "W_max_"         + tag
    String nameLmax        = "L_max_"         + tag
    String nameValid       = "ValidPt_"       + tag
    String nameXY          = "XY_"            + tag
    String nameEwin        = "Ewin_"          + tag

    Make/O/D/N=(Np) $nameLDOSIntRaw
    Make/O/D/N=(Np) $nameLDOSIntConv
    Make/O/D/N=(Np) $nameWmax
    Make/O/D/N=(Np) $nameLmax
    Make/O/B/N=(Np) $nameValid
    Duplicate/O Pos_nm, $nameXY
    Make/O/D/N=2 $nameEwin

    Wave LDOS_Int      = $nameLDOSIntRaw
    Wave LDOS_Int_Conv = $nameLDOSIntConv
    Wave W_max_List    = $nameWmax
    Wave L_max_List    = $nameLmax
    Wave ValidPt       = $nameValid
    Wave XY_out        = $nameXY
    Wave Ewin          = $nameEwin

    LDOS_Int      = NaN
    LDOS_Int_Conv = NaN
    W_max_List    = NaN
    L_max_List    = NaN
    ValidPt       = 0

    Ewin[0] = -Eint_eV
    Ewin[1] = +Eint_eV

    SetScale/P x, 0, 1, "", LDOS_Int
    SetScale/P x, 0, 1, "", LDOS_Int_Conv
    SetScale/P x, 0, 1, "", W_max_List
    SetScale/P x, 0, 1, "", L_max_List
    SetScale/P x, 0, 1, "", ValidPt

    String chanFolder = "root:SNS_MapTmp3D"
    NewDataFolder/O $chanFolder
    Variable useFreePhase = (smodel == 0) && (hasScreeningModel || !ParamIsDefault(maskFolder))
    String freePhaseFolder = ""
    if (useFreePhase)
        freePhaseFolder = SNS_EnsureFreeVortexPhase2D(Nmask, xV_nm, yV_nm, params.SNS_nFlux, stiffnessFolder, chanFolder + ":FreeVortexPhase")
    endif
    SetDataFolder $chanFolder

    Make/O/D/N=1 Bwin
    Bwin[0] = B0

    String nameDOS_EB   = "DOS_local_EB"
    String nameEaxisLoc = "E_axis_local"
    String nameDOS_EBbr = "DOS_local_EB_broad"

    Variable ip, iE
    Variable x0_nm, y0_nm
    Variable Nch
    Variable iLo, iHi, dE_eV
    Variable sumRaw, sumConv
    Variable Nvalid = 0

    for (ip = 0; ip < Np; ip += 1)

        x0_nm = Pos_nm[ip][0]
        y0_nm = Pos_nm[ip][1]

        if (SNS_SampleMaskNearest(Nmask, x0_nm, y0_nm) <= 0.5)
            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        SetDataFolder $chanFolder

        KillWaves/Z L_N_List, W_eff_List, wChan, T_eff_List
        KillWaves/Z Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
        KillWaves/Z betaExtra_List

        if (smodel == 2)

            if (ParamIsDefault(maskFolder))
                Abort "SNS_LDOSmap_BFixed_FromMask3D: ScreeningModel=2 requires maskFolder."
            endif

            Wave/Z QxCorr_loc = $(stiffnessFolder + ":QxCorr")
            Wave/Z QyCorr_loc = $(stiffnessFolder + ":QyCorr")
            Wave/Z PhaseReCorr_loc = $(stiffnessFolder + ":PhaseReCorr")
            Wave/Z PhaseImCorr_loc = $(stiffnessFolder + ":PhaseImCorr")

            if (!WaveExists(QxCorr_loc) || !WaveExists(QyCorr_loc) || \
                !WaveExists(PhaseReCorr_loc) || !WaveExists(PhaseImCorr_loc))
                Abort "SNS_LDOSmap_BFixed_FromMask3D: StiffnessPhase folder missing QxCorr/QyCorr or PhaseReCorr/PhaseImCorr."
            endif

            SNS_BuildChannelsFromMask3D(Nmask, x0_nm, y0_nm, phiB, stepFac, H_nm, maxPath_nm, chanFolder, \
                QxPhase=QxCorr_loc, QyPhase=QyCorr_loc, qNstep=1, \
                xV_nm=xV_nm, yV_nm=yV_nm, rCore_nm=15, PhaseReCore=PhaseReCorr_loc, PhaseImCore=PhaseImCorr_loc)

        elseif (smodel == 3)

            if (ParamIsDefault(maskFolder))
                Abort "SNS_LDOSmap_BFixed_FromMask3D: ScreeningModel=3 requires maskFolder."
            endif

            Wave/Z QxScreen_loc = $(screenedFolder + ":QxScreen")
            Wave/Z QyScreen_loc = $(screenedFolder + ":QyScreen")
            Wave/Z PhaseReCorr_loc = $(stiffnessFolder + ":PhaseReCorr")
            Wave/Z PhaseImCorr_loc = $(stiffnessFolder + ":PhaseImCorr")

            if (!WaveExists(QxScreen_loc) || !WaveExists(QyScreen_loc))
                Abort "SNS_LDOSmap_BFixed_FromMask3D: ScreenedStiffnessPhase folder missing QxScreen/QyScreen."
            endif
            if (!WaveExists(PhaseReCorr_loc) || !WaveExists(PhaseImCorr_loc))
                Abort "SNS_LDOSmap_BFixed_FromMask3D: StiffnessPhase folder missing PhaseReCorr/PhaseImCorr required for model 3."
            endif

            SNS_BuildChannelsFromMask3D(Nmask, x0_nm, y0_nm, phiB, stepFac, H_nm, maxPath_nm, chanFolder, \
                QxPhase=QxScreen_loc, QyPhase=QyScreen_loc, qNstep=1, \
                xV_nm=xV_nm, yV_nm=yV_nm, rCore_nm=15, PhaseReCore=PhaseReCorr_loc, PhaseImCore=PhaseImCorr_loc)

        elseif (smodel == 4)

            if (ParamIsDefault(maskFolder))
                Abort "SNS_LDOSmap_BFixed_FromMask3D: ScreeningModel=4 requires maskFolder."
            endif

            Wave/Z QxLondon_loc = $(londonFolder + ":QxLondonCorr")
            Wave/Z QyLondon_loc = $(londonFolder + ":QyLondonCorr")

            if (!WaveExists(QxLondon_loc) || !WaveExists(QyLondon_loc))
                Abort "SNS_LDOSmap_BFixed_FromMask3D: LocalLondonStiffnessPhase folder missing QxLondonCorr/QyLondonCorr."
            endif

            SNS_BuildChannelsFromMask3D(Nmask, x0_nm, y0_nm, phiB, stepFac, H_nm, maxPath_nm, chanFolder, \
                QxPhase=QxLondon_loc, QyPhase=QyLondon_loc, qNstep=1)

        else

            SNS_BuildChannelsFromMask3D(Nmask, x0_nm, y0_nm, phiB, stepFac, H_nm, maxPath_nm, chanFolder)

        endif

        Wave/Z L_N_List_nm   = L_N_List_nm
        Wave/Z W_eff_List_nm = W_eff_List_nm
        Wave/Z wChan         = wChan
        Wave/Z T_eff_List    = T_eff_List

        Wave/Z Hit1x_List_nm = Hit1x_List_nm
        Wave/Z Hit1y_List_nm = Hit1y_List_nm
        Wave/Z Hit2x_List_nm = Hit2x_List_nm
        Wave/Z Hit2y_List_nm = Hit2y_List_nm

        if (!WaveExists(L_N_List_nm) || !WaveExists(W_eff_List_nm) || !WaveExists(wChan) || !WaveExists(T_eff_List) || \
            !WaveExists(Hit1x_List_nm) || !WaveExists(Hit1y_List_nm) || !WaveExists(Hit2x_List_nm) || !WaveExists(Hit2y_List_nm))

            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        Nch = DimSize(L_N_List_nm, 0)

        if ((Nch <= 0) || \
            (Nch != DimSize(W_eff_List_nm,0)) || \
            (Nch != DimSize(wChan,0)) || \
            (Nch != DimSize(T_eff_List,0)) || \
            (Nch != DimSize(Hit1x_List_nm,0)) || \
            (Nch != DimSize(Hit1y_List_nm,0)) || \
            (Nch != DimSize(Hit2x_List_nm,0)) || \
            (Nch != DimSize(Hit2y_List_nm,0)))

            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        SetDataFolder dfrCaller
        W_max_List[ip] = WaveMax(W_eff_List_nm)
        L_max_List[ip] = WaveMax(L_N_List_nm)
        SetDataFolder $chanFolder

        if (useFreePhase || !ParamIsDefault(maskFolder))

            if (smodel == 2 || smodel == 3 || smodel == 4)

                Wave/Z betaExtra_List = betaExtra_List

                if (!WaveExists(betaExtra_List))
                    Abort "SNS_LDOSmap_BFixed_FromMask3D: Q-mode requested but betaExtra_List was not created by SNS_BuildChannelsFromMask3D."
                endif

                if (DimSize(betaExtra_List,0) != Nch)
                    Abort "SNS_LDOSmap_BFixed_FromMask3D: betaExtra_List length mismatch in Q-mode."
                endif

            else

                KillWaves/Z betaExtra_List

                Duplicate/O Nmask fullMask_loc
                fullMask_loc = 0

                if (smodel == 1)
                    Wave/Z PhaseReCorr_map3D = $(stiffnessFolder + ":PhaseReCorr")
                    Wave/Z PhaseImCorr_map3D = $(stiffnessFolder + ":PhaseImCorr")
                    if (!WaveExists(PhaseReCorr_map3D) || !WaveExists(PhaseImCorr_map3D))
                        Abort "SNS_LDOSmap_BFixed_FromMask3D: ScreeningModel=1 requires PhaseReCorr/PhaseImCorr in StiffnessPhase."
                    endif
                    SNS_ComputeBetaExtraFromExteriorPhase2D( \
                        fullMask_loc, PhaseReCorr_map3D, PhaseImCorr_map3D, \
                        Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
                        chanFolder)
                else
                    Wave/Z PhaseReFree_map3D = $(freePhaseFolder + "PhaseReFree")
                    Wave/Z PhaseImFree_map3D = $(freePhaseFolder + "PhaseImFree")
                    if (!WaveExists(PhaseReFree_map3D) || !WaveExists(PhaseImFree_map3D))
                        Abort "SNS_LDOSmap_BFixed_FromMask3D: free-vortex phase generation did not produce PhaseReFree/PhaseImFree."
                    endif
                    SNS_ComputeBetaExtraFromExteriorPhase2D( \
                        fullMask_loc, PhaseReFree_map3D, PhaseImFree_map3D, \
                        Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
                        chanFolder)
                endif

                Wave/Z betaExtra_List = betaExtra_List

                if (!WaveExists(betaExtra_List))
                    Abort "SNS_LDOSmap_BFixed_FromMask3D: endpoint betaExtra_List was not created."
                endif

                if (DimSize(betaExtra_List,0) != Nch)
                    Abort "SNS_LDOSmap_BFixed_FromMask3D: endpoint betaExtra_List length mismatch."
                endif
            endif

            SNS_ComputeDOS_FromChannels( \
                Bwin, \
                L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
                Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
                xV_nm, yV_nm, \
                N_cont_statesPer_eV, \
                nameDOS_EB, nameEaxisLoc, \
                betaExtra_List=betaExtra_List)

        else

            SNS_ComputeDOS_FromChannels( \
                Bwin, \
                L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
                Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
                xV_nm, yV_nm, \
                N_cont_statesPer_eV, \
                nameDOS_EB, nameEaxisLoc)

        endif

        Wave/Z DOS_EB_raw   = $nameDOS_EB
        Wave/Z E_axis_local = $nameEaxisLoc

        if (!WaveExists(DOS_EB_raw) || !WaveExists(E_axis_local))
            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        SNS_ApplyDOS_Broadening_TplusMod(DOS_EB_raw, nameDOS_EBbr)
        Wave/Z DOS_EB_broad = $nameDOS_EBbr

        if (!WaveExists(DOS_EB_broad))
            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        iLo = floor(x2pnt(E_axis_local, -Eint_eV))
        iHi = ceil( x2pnt(E_axis_local, +Eint_eV))

        iLo = max(0, iLo)
        iHi = min(NE-1, iHi)

        if (iHi < iLo)
            SetDataFolder dfrCaller
            LDOS_Int[ip]      = NaN
            LDOS_Int_Conv[ip] = NaN
            W_max_List[ip]    = NaN
            L_max_List[ip]    = NaN
            ValidPt[ip]       = 0
            SetDataFolder $chanFolder
            continue
        endif

        dE_eV = abs(DimDelta(E_axis_local, 0))

        sumRaw  = 0
        sumConv = 0
        for (iE = iLo; iE <= iHi; iE += 1)
            sumRaw  += DOS_EB_raw[iE][0]
            sumConv += DOS_EB_broad[iE][0]
        endfor

        SetDataFolder dfrCaller
        LDOS_Int[ip]      = sumRaw  * dE_eV
        LDOS_Int_Conv[ip] = sumConv * dE_eV
        ValidPt[ip]       = 1
        SetDataFolder $chanFolder

        Nvalid += 1

    endfor

    SetDataFolder dfrCaller

    String meta
    meta  = "SNS_ChannelBuilder=3D;"
    meta += "SNS_B0_T=" + num2str(B0) + ";"
    meta += "SNS_phiB_rad=" + num2str(phiB) + ";"
    meta += "SNS_H_nm=" + num2str(H_nm) + ";"
    meta += "SNS_maxPath_nm=" + num2str(maxPath_nm) + ";"
    meta += "SNS_ScreeningModel=" + num2str(smodel) + ";"
    meta += "SNS_Eint_meV=" + num2str(Eint_meV) + ";"

    if (!ParamIsDefault(maskFolder))
        meta += "SNS_MaskFolder=" + maskBase + ";"
        meta += "SNS_StiffnessPhaseFolder=" + stiffnessFolder + ";"
        meta += "SNS_ScreenedStiffnessPhaseFolder=" + screenedFolder + ";"
        meta += "SNS_LocalLondonStiffnessPhaseFolder=" + londonFolder + ";"
    else
        meta += "SNS_MaskFolder=;"
        meta += "SNS_StiffnessPhaseFolder=;"
        meta += "SNS_ScreenedStiffnessPhaseFolder=;"
        meta += "SNS_LocalLondonStiffnessPhaseFolder=;"
    endif
    meta += "SNS_FreePhaseFolder=" + freePhaseFolder + ";"

    Note/K LDOS_Int
    Note LDOS_Int, meta

    Note/K LDOS_Int_Conv
    Note LDOS_Int_Conv, meta

    return Nvalid
End

//==============================================================================
// POINT ANGLE-SWEEP WORKFLOWS
//
// Fixed-position workflows that rebuild the channel ensemble while rotating the
// in-plane magnetic-field direction.
//==============================================================================

//==============================================================================
// SNS_PointLDOS_PhiSweep_FromMask
//
// Purpose:
//   Compute LDOS(E,phiB) at one fixed point (x0,y0) inside the N region by
//   rotating the in-plane magnetic-field direction at fixed field magnitude B0.
//
//   For each angle phiB, the function:
//     1. rebuilds the local 2D SNS channel ensemble,
//     2. computes DOS(E,B0) from the channel ensemble,
//     3. applies thermal + lock-in modulation broadening,
//     4. stores LDOS(E,phiB) and zero-bias LDOS vs angle.
//
// Inputs:
//   Nmask      : 2D binary N-region mask.
//                Inside N: > 0.5, outside: <= 0.5.
//                Wave axes must be scaled in nm.
//
//   x0, y0     : spectroscopy position [nm]
//   B0         : fixed magnetic-field magnitude [T]
//   phiStepDeg : angular step size [deg]
//   stepFac    : ray/channel sampling refinement factor
//
// Optional Inputs:
//   phiStartDeg : start angle [deg]. Default: 0.
//   phiEndDeg   : end angle [deg]. Default: 360 - phiStepDeg.
//
// Interface transparency:
//   This function no longer accepts BTK_barrier / Zbarrier.
//   Interface transparency is controlled globally through SNS_Settings and
//   applied inside SNS_BuildChannelsFromMask2D.
//
// Outputs:
//   LDOS_phi<tag>        : broadened DOS(E,phi)
//   LDOS_Gamma_phi<tag>  : broadened gamma-diagnostic DOS(E,phi), if available
//   E_axis_phi<tag>      : energy axis [eV]
//   Phi_axis<tag>        : angle axis [deg]
//   ZBC_phi<tag>         : zero-bias LDOS vs angle
//
// Returns:
//   Nphi                : number of angle points.
//
// Notes:
//   Channels are rebuilt for every phiB because W_eff and trajectory geometry
//   depend on field direction.
//   This is a point/angle spatial workflow and belongs in SNS_SpatialMaps.ipf.
//==============================================================================
Function SNS_PointLDOS_PhiSweep_FromMask(Nmask, x0, y0, B0, phiStepDeg, stepFac, [phiStartDeg, phiEndDeg])
    Wave Nmask
    Variable x0, y0                  // [nm]
    Variable B0                      // [T]
    Variable phiStepDeg              // [deg]
    Variable stepFac
    Variable phiStartDeg, phiEndDeg

    DFREF dfrCaller = GetDataFolderDFR()

    if (ParamIsDefault(phiStartDeg))
        phiStartDeg = 0
    endif
    if (ParamIsDefault(phiEndDeg))
        phiEndDeg = 360 - phiStepDeg
    endif

    if (phiStepDeg <= 0)
        Abort "SNS_PointLDOS_PhiSweep_FromMask: phiStepDeg must be > 0."
    endif
    if (phiEndDeg < phiStartDeg)
        Abort "SNS_PointLDOS_PhiSweep_FromMask: phiEndDeg must be >= phiStartDeg."
    endif

    // ---------- load standard settings ----------
    STRUCT SNS_Params params
    SNS_LoadParams(params)

    Variable NE = params.NE

    // geometry / vortex folder = folder of Nmask
    DFREF dfrGeom = GetWavesDataFolderDFR(Nmask)
    SetDataFolder dfrGeom

	 Wave/Z Vortex_ptx = Vortex_ptx
	 Wave/Z Vortex_pty = Vortex_pty
	 Variable xV_nm, yV_nm
	 if (WaveExists(Vortex_ptx) && WaveExists(Vortex_pty))
	    xV_nm = Vortex_ptx[0]
	    yV_nm = Vortex_pty[0]
	 else
	    xV_nm = NaN
	    yV_nm = NaN
	 endif

    Wave/Z w_area_nm2 = w_area_nm2
    Variable N_cont_statesPer_eV
    if (WaveExists(w_area_nm2))
        Variable area_nm2 = w_area_nm2[0]
        Variable area_m2  = area_nm2 * 1e-18
        N_cont_statesPer_eV = params.DOS2D_eV_Area * area_m2
    else
        N_cont_statesPer_eV = 0
    endif

    // ---------- angle axis ----------
    SetDataFolder dfrCaller

    Variable Nphi = floor((phiEndDeg - phiStartDeg)/phiStepDeg + 0.5) + 1
    String tag = "ptx_" + num2str(round(x0)) + "nm_pty" + num2str(round(y0)) \
               + "nm_B" + num2str(round(B0*1e3)) + "mT"

    String nameLDOS = "LDOS_phi" + tag
    String nameLDOS_Gamma = "LDOS_Gamma_phi" + tag
    String nameE    = "E_axis_phi" + tag
    String namePhi  = "Phi_axis" + tag
    String nameZBC  = "ZBC_phi" + tag

    Make/O/D/N=(Nphi) $namePhi
    Wave Phi_axis = $namePhi
    Phi_axis = phiStartDeg + p*phiStepDeg
    SetScale/P x, phiStartDeg, phiStepDeg, "deg", Phi_axis

    Make/O/D/N=(NE, Nphi) $nameLDOS
    Wave LDOS_phi = $nameLDOS
    LDOS_phi = NaN
    
    Make/O/D/N=(NE, Nphi) $nameLDOS_Gamma
    Wave LDOS_Gamma_phi = $nameLDOS_Gamma
    LDOS_Gamma_phi = NaN

    // ---------- temp folder ----------
    String chanFolder = "root:SNS_PhiTmp"
    NewDataFolder/O $chanFolder
    SetDataFolder $chanFolder

    Make/O/D/N=1 Bwin
    Bwin[0] = B0

    String nameDOS_EB        = "DOS_local_EB"
    String nameEaxisLoc      = "E_axis_local"
    String nameDOS_EBbr      = "DOS_local_EB_broad"
    String nameDOS_EBGammabr = "DOS_local_EB_Gamma_broad"

    Variable iphi, iE, Nch
    for (iphi = 0; iphi < Nphi; iphi += 1)

        Variable phiDeg = Phi_axis[iphi]
        Variable phiRad = phiDeg * pi / 180

        // build channels for this point and field direction
        SNS_BuildChannelsFromMask2D(Nmask, x0, y0, phiRad, stepFac, chanFolder)

        Wave L_N_List_nm
        Wave W_eff_List_nm
        Wave wChan
        Wave T_eff_List
        Wave Hit1x_List_nm
        Wave Hit1y_List_nm
        Wave Hit2x_List_nm
        Wave Hit2y_List_nm

        Nch = DimSize(L_N_List_nm, 0)
        if ( (Nch <= 0) || (Nch != DimSize(W_eff_List_nm,0)) || \
             (Nch != DimSize(wChan,0)) || (Nch != DimSize(T_eff_List,0)) || \
             (Nch != DimSize(Hit1x_List_nm,0)) || (Nch != DimSize(Hit2x_List_nm,0)) )

				SetDataFolder dfrCaller
				LDOS_phi[][iphi] = NaN
				LDOS_Gamma_phi[][iphi] = NaN
				SetDataFolder $chanFolder
				continue
        endif

        // raw DOS
        SNS_ComputeDOS_FromChannels( \
            Bwin, \
            L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
            Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
            xV_nm, yV_nm, \
            N_cont_statesPer_eV, \
            nameDOS_EB, nameEaxisLoc)
			
        Wave DOS_EB_raw   = $nameDOS_EB
        Wave E_axis_local = $nameEaxisLoc
        Wave/Z DOS_EB_Gamma = $(nameDOS_EB + "_Gamma")



        // broaden
        SNS_ApplyDOS_Broadening_TplusMod(DOS_EB_raw, nameDOS_EBbr)
        Wave DOS_EB_broad = $nameDOS_EBbr
		 if (WaveExists(DOS_EB_Gamma))
		    SNS_ApplyDOS_Broadening_TplusMod(DOS_EB_Gamma, nameDOS_EBGammabr)
		    Wave DOS_EB_Gamma_broad = $nameDOS_EBGammabr
		 endif


        // copy energy axis once
        if (iphi == 0)
            SetDataFolder dfrCaller
            Duplicate/O E_axis_local, $nameE
            SetDataFolder $chanFolder
        endif

		SetDataFolder dfrCaller
		for (iE = 0; iE < NE; iE += 1)
		    LDOS_phi[iE][iphi] = DOS_EB_broad[iE][0]
		    if (WaveExists(DOS_EB_Gamma))
		        LDOS_Gamma_phi[iE][iphi] = DOS_EB_Gamma_broad[iE][0]
		    else
		        LDOS_Gamma_phi[iE][iphi] = NaN
		    endif
		endfor
		SetDataFolder $chanFolder
    endfor

    // ---------- attach axes ----------
    SetDataFolder dfrCaller
    Wave E_axis = $nameE

    // x = energy, y = phi
    SetScale/P x, DimOffset(E_axis,0), DimDelta(E_axis,0), WaveUnits(E_axis,0), LDOS_phi
    SetScale/P y, DimOffset(Phi_axis,0), DimDelta(Phi_axis,0), WaveUnits(Phi_axis,0), LDOS_phi
    SetScale/P x, DimOffset(E_axis,0), DimDelta(E_axis,0), WaveUnits(E_axis,0), LDOS_Gamma_phi
    SetScale/P y, DimOffset(Phi_axis,0), DimDelta(Phi_axis,0), WaveUnits(Phi_axis,0), LDOS_Gamma_phi

    // optional ZBC(phi)
    Make/O/D/N=(Nphi) $nameZBC
    Wave ZBC_phi = $nameZBC

    Variable i0 = BinarySearchInterp(E_axis, 0)
    if (numtype(i0) != 0)
        ZBC_phi = NaN
    else
        Variable iLo = floor(i0)
        Variable iHi = ceil(i0)
        if (iLo == iHi)
            ZBC_phi = LDOS_phi[iLo][p]
        else
            Variable wHi = i0 - iLo
            Variable wLo = 1 - wHi
            ZBC_phi = wLo*LDOS_phi[iLo][p] + wHi*LDOS_phi[iHi][p]
        endif
    endif
    SetScale/P x, DimOffset(Phi_axis,0), DimDelta(Phi_axis,0), "deg", ZBC_phi

    return Nphi
End

//==============================================================================
// SNS_PointLDOS_PhiSweep_360
//
// Purpose:
//   Convenience wrapper for SNS_PointLDOS_PhiSweep_FromMask using a full
//   0–360 degree rotation of the in-plane magnetic field.
//
// Inputs:
//   Nmask      : 2D binary N-region mask; axes scaled in nm.
//   x0, y0     : spectroscopy position [nm]
//   B0_mT      : fixed magnetic-field magnitude [mT]
//   phiStepDeg : angular step [deg]
//   stepFac    : ray/channel sampling refinement factor
//
// Interface transparency:
//   This function no longer accepts BTK_barrier / Zbarrier.
//   Interface transparency is controlled globally through SNS_Settings and
//   applied inside SNS_BuildChannelsFromMask2D.
//
// Outputs:
//   Same as SNS_PointLDOS_PhiSweep_FromMask.
//
// Returns:
//   Nphi : number of angle points.
//==============================================================================
Function SNS_PointLDOS_PhiSweep_360(Nmask, x0, y0, B0_mT, phiStepDeg, stepFac)

    Wave Nmask
    Variable x0, y0
    Variable B0_mT
    Variable phiStepDeg
    Variable stepFac

    Variable B0 = B0_mT * 1e-3

    return SNS_PointLDOS_PhiSweep_FromMask( \
        Nmask, \
        x0, y0, \
        B0, \
        phiStepDeg, \
        stepFac, \
        phiStartDeg = 0, \
        phiEndDeg = 360 - phiStepDeg)

End

//==============================================================================
// SPATIAL RAY-DIAGNOSTIC WORKFLOWS
//
// Diagnostic line workflows that repeatedly call ray-trace histogram helpers and
// collect supported extrema of L_N, W_eff, and related quantities.
//==============================================================================


//==============================================================================
// SNS_LineExtremeStats_FromMask
//
// Purpose:
//   Evaluate supported ray-tracing extrema along a line through a 2D N-region
//   mask.
//
//   For each line position, this function runs:
//
//       SNS_RayTrace_Hist_LN_Weff(..., doDisplay=0)
//
//   and reads the supported extrema from:
//
//       <maskFolder>:RayTraceHist
//
//   The stored extrema are determined by the upper-tail discriminator used by
//   SNS_RayTrace_Hist_LN_Weff.
//
// Inputs:
//   Nmask   : 2D binary N-region mask; axes scaled in nm.
//   phiB    : in-plane magnetic-field angle [rad]
//   xStart,
//   yStart  : line start point [nm]
//   xEnd,
//   yEnd    : line end point [nm]
//   stepFac : ray/channel sampling refinement factor
//
// Optional Inputs:
//   NLinePts  : optional number of line positions.
//               If supplied, overrides the default LambdaF-based line step.
//               Must be finite integer >= 2.
//
//   nBins     : number of histogram bins. Default: 50.
//   binL_nm   : optional fixed L_N histogram bin width [nm].
//   binW_nm   : optional fixed W_eff histogram bin width [nm].
//   absW      : pass-through flag for |W_eff| handling. Default: 1.
//   normalize : pass-through histogram normalization flag. Default: 1.
//   tagSuffix : optional suffix for output wave names.
//
// Interface transparency:
//   This function no longer accepts BTK_barrier / Zbarrier.
//   Interface transparency is controlled globally through SNS_Settings and
//   should be handled inside SNS_RayTrace_Hist_LN_Weff / channel construction.
//
// Outputs:
//   R_axis_ExtremeStats*
//   L_extreme_ExtremeStats*
//   W_extreme_ExtremeStats*
//   L_sigma_ExtremeStats*
//   W_sigma_ExtremeStats*
//   L_trace_ExtremeStats*
//   W_trace_ExtremeStats*
//   L_tailN_ExtremeStats*
//   W_tailN_ExtremeStats*
//   L_clusterStart_ExtremeStats*
//   W_clusterStart_ExtremeStats*
//   L_clusterSize_ExtremeStats*
//   W_clusterSize_ExtremeStats*
//   L_minClusterSize_ExtremeStats*
//   W_minClusterSize_ExtremeStats*
//
// Returns:
//   Nr : number of line positions.
//
// Notes:
//   This is a spatial diagnostic workflow and belongs in SNS_SpatialMaps.ipf.
//   Optional binL_nm / binW_nm are passed through to SNS_RayTrace_Hist_LN_Weff
//   so supported extrema can be evaluated with fixed histogram bin widths.
//==============================================================================
Function SNS_LineExtremeStats_FromMask(Nmask, phiB, xStart, yStart, xEnd, yEnd, stepFac, [BTK_barrier, NLinePts, nBins, binL_nm, binW_nm, absW, normalize, tagSuffix])
    Wave Nmask
    Variable phiB
    Variable xStart, yStart, xEnd, yEnd
    Variable stepFac
    Variable BTK_barrier
    Variable NLinePts
    Variable nBins, binL_nm, binW_nm, absW, normalize
    String tagSuffix

    DFREF dfrCaller = GetDataFolderDFR()

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    // ---------------- Defaults ----------------
    if (ParamIsDefault(BTK_barrier))
        BTK_barrier = params.BTK_barrier
    endif
    if (ParamIsDefault(nBins))
        nBins = 50
    endif
    if (ParamIsDefault(absW))
        absW = 1
    endif
    if (ParamIsDefault(normalize))
        normalize = 1
    endif
    if (ParamIsDefault(tagSuffix))
        tagSuffix = ""
    endif

    if (numtype(BTK_barrier) != 0 || BTK_barrier < 0)
        SetDataFolder dfrCaller
        Abort "SNS_LineExtremeStats_FromMask: BTK_barrier must be finite and non-negative."
    endif

    if (!ParamIsDefault(binL_nm))
        if (numtype(binL_nm) != 0 || binL_nm <= 0)
            SetDataFolder dfrCaller
            Abort "SNS_LineExtremeStats_FromMask: binL_nm must be finite and positive."
        endif
    endif

    if (!ParamIsDefault(binW_nm))
        if (numtype(binW_nm) != 0 || binW_nm <= 0)
            SetDataFolder dfrCaller
            Abort "SNS_LineExtremeStats_FromMask: binW_nm must be finite and positive."
        endif
    endif

    nBins = max(2, round(nBins))

    Variable lambdaF_nm = params.LambdaF

    // ---------------- Geometry / mask folder ----------------
    DFREF dfrGeom = GetWavesDataFolderDFR(Nmask)
    SetDataFolder dfrGeom

    String maskFolder = GetDataFolder(1)

    // SNS_RayTrace_Hist_LN_Weff expects a wave named w_mask in dfPath.
    // If the input mask has a different name, create/update a local duplicate.
    if (cmpstr(NameOfWave(Nmask), "w_mask") != 0)
        Duplicate/O Nmask, w_mask
    endif

    SetDataFolder dfrCaller

    // ---------------- Output tag ----------------
    String suffix = ""
    if (strlen(tagSuffix) > 0)
        suffix = "_" + CleanupName(tagSuffix, 0)
    endif

    String nameRaxis       = "R_axis_ExtremeStats" + suffix
    String nameLext        = "L_extreme_ExtremeStats" + suffix
    String nameWext        = "W_extreme_ExtremeStats" + suffix
    String nameLsig        = "L_sigma_ExtremeStats" + suffix
    String nameWsig        = "W_sigma_ExtremeStats" + suffix
    String nameLtrace      = "L_trace_ExtremeStats" + suffix
    String nameWtrace      = "W_trace_ExtremeStats" + suffix
    String nameLNtail      = "L_tailN_ExtremeStats" + suffix
    String nameWNtail      = "W_tailN_ExtremeStats" + suffix
    String nameLclStart    = "L_clusterStart_ExtremeStats" + suffix
    String nameWclStart    = "W_clusterStart_ExtremeStats" + suffix
    String nameLclSize     = "L_clusterSize_ExtremeStats" + suffix
    String nameWclSize     = "W_clusterSize_ExtremeStats" + suffix
    String nameLminClSize  = "L_minClusterSize_ExtremeStats" + suffix
    String nameWminClSize  = "W_minClusterSize_ExtremeStats" + suffix

    // ---------------- Line positions ----------------
    String nameX = "SNS_ExtremeStats_Line_X" + suffix
    String nameY = "SNS_ExtremeStats_Line_Y" + suffix

    Variable lineStep_nm = lambdaF_nm

    if (!ParamIsDefault(NLinePts))
        NLinePts = round(NLinePts)

        if ((numtype(NLinePts) != 0) || (NLinePts < 2))
            SetDataFolder dfrCaller
            Abort "SNS_LineExtremeStats_FromMask: NLinePts must be a finite integer >= 2."
        endif

        Variable lineLength_nm = sqrt((xEnd - xStart)^2 + (yEnd - yStart)^2)

        if ((numtype(lineLength_nm) != 0) || (lineLength_nm <= 0))
            SetDataFolder dfrCaller
            Abort "SNS_LineExtremeStats_FromMask: cannot use NLinePts for a zero-length or invalid line."
        endif

        lineStep_nm = lineLength_nm / (NLinePts - 1)
    endif

    Variable lineStep_m = lineStep_nm * 1e-9

    Variable Nr = SNS_BuildLinePositions(xStart, yStart, xEnd, yEnd, lineStep_m, nameX, nameY, nameRaxis)

    Wave Xline = $nameX
    Wave Yline = $nameY
    Wave Rline = $nameRaxis

    if (Nr > 1)
        Variable drnm = (Rline[Nr-1] - Rline[0]) / (Nr - 1)
        SetScale/P x, Rline[0], drnm, "nm", Rline
    else
        SetScale/P x, 0, 1, "nm", Rline
    endif

    // ---------------- Allocate outputs ----------------
    Make/O/D/N=(Nr) $nameLext, $nameWext
    Make/O/D/N=(Nr) $nameLsig, $nameWsig
    Make/O/D/N=(Nr) $nameLtrace, $nameWtrace
    Make/O/D/N=(Nr) $nameLNtail, $nameWNtail
    Make/O/D/N=(Nr) $nameLclStart, $nameWclStart
    Make/O/D/N=(Nr) $nameLclSize, $nameWclSize
    Make/O/D/N=(Nr) $nameLminClSize, $nameWminClSize

    Wave L_extreme_Line        = $nameLext
    Wave W_extreme_Line        = $nameWext
    Wave L_sigma_Line          = $nameLsig
    Wave W_sigma_Line          = $nameWsig
    Wave L_trace_Line          = $nameLtrace
    Wave W_trace_Line          = $nameWtrace
    Wave L_tailN_Line          = $nameLNtail
    Wave W_tailN_Line          = $nameWNtail
    Wave L_clusterStart_Line   = $nameLclStart
    Wave W_clusterStart_Line   = $nameWclStart
    Wave L_clusterSize_Line    = $nameLclSize
    Wave W_clusterSize_Line    = $nameWclSize
    Wave L_minClusterSize_Line = $nameLminClSize
    Wave W_minClusterSize_Line = $nameWminClSize

    L_extreme_Line        = NaN
    W_extreme_Line        = NaN
    L_sigma_Line          = NaN
    W_sigma_Line          = NaN
    L_trace_Line          = NaN
    W_trace_Line          = NaN
    L_tailN_Line          = NaN
    W_tailN_Line          = NaN
    L_clusterStart_Line   = NaN
    W_clusterStart_Line   = NaN
    L_clusterSize_Line    = NaN
    W_clusterSize_Line    = NaN
    L_minClusterSize_Line = NaN
    W_minClusterSize_Line = NaN

    // ---------------- Main loop ----------------
    String histFolder
    Variable ir

    for (ir = 0; ir < Nr; ir += 1)

        Variable r0x = Xline[ir]
        Variable r0y = Yline[ir]

        // Skip line positions outside the N mask.
        // Nmask axes are in nm.
        Variable px = round((r0x - DimOffset(Nmask, 0)) / DimDelta(Nmask, 0))
        Variable py = round((r0y - DimOffset(Nmask, 1)) / DimDelta(Nmask, 1))

        if (px < 0 || px >= DimSize(Nmask, 0) || py < 0 || py >= DimSize(Nmask, 1))
            SNS_Log("SNS_LineExtremeStats_FromMask: line point outside mask bounds at index " + num2str(ir) + ". Leaving NaNs.", level="WARN")
            continue
        endif

        if (Nmask[px][py] <= 0)
            SNS_Log("SNS_LineExtremeStats_FromMask: line point outside N region at index " + num2str(ir) + ". Leaving NaNs.", level="WARN")
            continue
        endif

        // Run geometry/statistics helper only. No histogram display.
        if (!ParamIsDefault(binL_nm) && !ParamIsDefault(binW_nm))
            SNS_RayTrace_Hist_LN_Weff(maskFolder, Bangle_deg=phiB*180/pi, STSx=r0x, STSy=r0y, nBins=nBins, binL_nm=binL_nm, binW_nm=binW_nm, absW=absW, normalize=normalize, showFitInfo=0, doDisplay=0)
        elseif (!ParamIsDefault(binL_nm))
            SNS_RayTrace_Hist_LN_Weff(maskFolder, Bangle_deg=phiB*180/pi, STSx=r0x, STSy=r0y, nBins=nBins, binL_nm=binL_nm, absW=absW, normalize=normalize, showFitInfo=0, doDisplay=0)
        elseif (!ParamIsDefault(binW_nm))
            SNS_RayTrace_Hist_LN_Weff(maskFolder, Bangle_deg=phiB*180/pi, STSx=r0x, STSy=r0y, nBins=nBins, binW_nm=binW_nm, absW=absW, normalize=normalize, showFitInfo=0, doDisplay=0)
        else
            SNS_RayTrace_Hist_LN_Weff(maskFolder, Bangle_deg=phiB*180/pi, STSx=r0x, STSy=r0y, nBins=nBins, absW=absW, normalize=normalize, showFitInfo=0, doDisplay=0)
        endif

        histFolder = maskFolder + "RayTraceHist"

        if (!DataFolderExists(histFolder))
            SNS_Log("SNS_LineExtremeStats_FromMask: missing RayTraceHist folder after SNS_RayTrace_Hist_LN_Weff at line index " + num2str(ir) + ". Leaving NaNs.", level="WARN")
            SetDataFolder dfrCaller
            continue
        endif

        SetDataFolder dfrCaller

        Variable diagVal

        diagVal = SNS_RayDiagValue(histFolder + ":", "L", "tailCluster", "max")
        if (numtype(diagVal) == 0)
            L_extreme_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "W", "tailCluster", "max")
        if (numtype(diagVal) == 0)
            W_extreme_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "L", "tailCluster", "sigma")
        if (numtype(diagVal) == 0)
            L_sigma_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "W", "tailCluster", "sigma")
        if (numtype(diagVal) == 0)
            W_sigma_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "L", "tailCluster", "max_trace")
        if (numtype(diagVal) == 0)
            L_trace_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "W", "tailCluster", "max_trace")
        if (numtype(diagVal) == 0)
            W_trace_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "L", "tailCluster", "n")
        if (numtype(diagVal) == 0)
            L_tailN_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "W", "tailCluster", "n")
        if (numtype(diagVal) == 0)
            W_tailN_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "L", "tailCluster", "start")
        if (numtype(diagVal) == 0)
            L_clusterStart_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "W", "tailCluster", "start")
        if (numtype(diagVal) == 0)
            W_clusterStart_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "L", "tailCluster", "size")
        if (numtype(diagVal) == 0)
            L_clusterSize_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "W", "tailCluster", "size")
        if (numtype(diagVal) == 0)
            W_clusterSize_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "L", "tailCluster", "minSize")
        if (numtype(diagVal) == 0)
            L_minClusterSize_Line[ir] = diagVal
        endif
        diagVal = SNS_RayDiagValue(histFolder + ":", "W", "tailCluster", "minSize")
        if (numtype(diagVal) == 0)
            W_minClusterSize_Line[ir] = diagVal
        endif
    endfor

    // ---------------- Axes ----------------
    if (Nr > 1)
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", L_extreme_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", W_extreme_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", L_sigma_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", W_sigma_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", L_trace_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", W_trace_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", L_tailN_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", W_tailN_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", L_clusterStart_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", W_clusterStart_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", L_clusterSize_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", W_clusterSize_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", L_minClusterSize_Line
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", W_minClusterSize_Line
    else
        SetScale/P x, 0, 1, "nm", L_extreme_Line
        SetScale/P x, 0, 1, "nm", W_extreme_Line
        SetScale/P x, 0, 1, "nm", L_sigma_Line
        SetScale/P x, 0, 1, "nm", W_sigma_Line
        SetScale/P x, 0, 1, "nm", L_trace_Line
        SetScale/P x, 0, 1, "nm", W_trace_Line
        SetScale/P x, 0, 1, "nm", L_tailN_Line
        SetScale/P x, 0, 1, "nm", W_tailN_Line
        SetScale/P x, 0, 1, "nm", L_clusterStart_Line
        SetScale/P x, 0, 1, "nm", W_clusterStart_Line
        SetScale/P x, 0, 1, "nm", L_clusterSize_Line
        SetScale/P x, 0, 1, "nm", W_clusterSize_Line
        SetScale/P x, 0, 1, "nm", L_minClusterSize_Line
        SetScale/P x, 0, 1, "nm", W_minClusterSize_Line
    endif

    // ---------------- Metadata ----------------
    String meta
    meta  = "SNS_SourceFunction=SNS_LineExtremeStats_FromMask;"
    meta += "SNS_ExtremeDefinition=supported cluster maximum;"
    meta += "SNS_HistOutputFolder=RayTraceHist;"
    meta += "SNS_SupportCriterion=q>=0.9*qmax;contiguous_cluster_size>=max(5,round(0.01*nRay));"
    meta += "SNS_SelectedCluster=largest_cluster_size_tiebreak_larger_median;"
    meta += "SNS_MaskFolder=" + maskFolder + ";"
    meta += "SNS_phiB_rad=" + num2str(phiB) + ";"
    meta += "SNS_Bangle_deg=" + num2str(phiB*180/pi) + ";"
    meta += "SNS_BTK_barrier=" + num2str(BTK_barrier) + ";"
    meta += "SNS_LineStep_nm=" + num2str(lineStep_nm) + ";"
    meta += "SNS_LineStep_m=" + num2str(lineStep_m) + ";"
    meta += "SNS_absW=" + num2str(absW) + ";"
    meta += "SNS_nBins=" + num2str(nBins) + ";"
    meta += "SNS_normalize=" + num2str(normalize) + ";"
    if (!ParamIsDefault(binL_nm))
        meta += "SNS_binL_nm=" + num2str(binL_nm) + ";"
    endif
    if (!ParamIsDefault(binW_nm))
        meta += "SNS_binW_nm=" + num2str(binW_nm) + ";"
    endif
    meta += "SNS_Raxis=" + GetWavesDataFolder(Rline, 2) + ";"

    Note/K L_extreme_Line
    Note L_extreme_Line, meta

    Note/K W_extreme_Line
    Note W_extreme_Line, meta

    Note/K L_sigma_Line
    Note L_sigma_Line, meta

    Note/K W_sigma_Line
    Note W_sigma_Line, meta

    Note/K L_trace_Line
    Note L_trace_Line, meta

    Note/K W_trace_Line
    Note W_trace_Line, meta

    Note/K L_tailN_Line
    Note L_tailN_Line, meta

    Note/K W_tailN_Line
    Note W_tailN_Line, meta

    Note/K L_clusterStart_Line
    Note L_clusterStart_Line, meta

    Note/K W_clusterStart_Line
    Note W_clusterStart_Line, meta

    Note/K L_clusterSize_Line
    Note L_clusterSize_Line, meta

    Note/K W_clusterSize_Line
    Note W_clusterSize_Line, meta

    Note/K L_minClusterSize_Line
    Note L_minClusterSize_Line, meta

    Note/K W_minClusterSize_Line
    Note W_minClusterSize_Line, meta

    SetDataFolder dfrCaller

    return Nr
End
