#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma DefaultTab={3,20,4}

//==============================================================================
// SNS_DOS.ipf
//
// DOS assembly from precomputed SNS channel waves.
//
// Responsibilities:
//   - field/channel-dependent Lorentzian broadening Γ(B,L,W)
//   - DOS(E,B) construction from solved ABS branches
//   - continuum-background normalization
//   - high-level wrapper SNS_ComputeDOS_FromSettings
//
// Dependencies:
//   SNS_Core.ipf
//     - Structure SNS_Params
//     - SNS_LoadParams
//     - constants HBAR_eVs, etc.
//
//   SNS_Logging.ipf
//     - SNS_Log
//
//   SNS_Solver.ipf
//     - SNS_beta
//     - Solve_AllBranches_SNS_dGSJ_betaextra
//
//   SNS_RayTrace2D.ipf or future SNS_ChannelSchema/SNS_GeometryUtils.ipf
//     - SNS_MaskAreaPerim_FromParticles
//
//   SNS_Utilities.ipf or SNS_DOS_Display.ipf
//     - MoveBranchWavesToSubfolder
//
// Notes:
//   This file should not perform ray tracing, build masks, draw plots, or
//   implement interface transparency.
//==============================================================================


//==============================================================================
// SNS_Gamma_of_B
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   Bval : input
//   B0 : input
//   Gamma0 : input
//   alpha1 : input
//   alpha2 : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_Gamma_of_B(Bval, B0, Gamma0, alpha1, alpha2)
    Variable Bval, B0, Gamma0, alpha1, alpha2

    Variable bDim
    if (B0 > 0)
        bDim = abs(Bval) / B0
    else
        bDim = 0
    endif

    return Gamma0 * (1 + alpha1*bDim + alpha2*bDim*bDim)
End

//==============================================================================
// SNS_ClampGamma
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   GammaIn : input
//   GammaMin : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_ClampGamma(GammaIn, GammaMin)
    Variable GammaIn, GammaMin
    if (GammaIn < GammaMin)
        return GammaMin
    endif
    return GammaIn
End

//============================================================
// SNS_ComputeGammaTot
//
// Purpose:
//   Compute total Lorentzian width Γ_tot(B, Lch, Wch) for one
//   channel and one field value.
//
//   The zero-field width is fixed by the LDOS fit:
//
//       Γ_base = SNS_GammaBase_eV, if valid
//              = SNS_p.Broadening otherwise.
//
//   Field-dependent depairing is represented by a quadratic
//   additive term:
//
//       Γ_pair(B) = SNS_GammaPairScale * B^2
//
//   with SNS_GammaPairScale in eV/T^2.
//
//   The other optional mechanisms remain quadrature terms:
//
//       Γ_quad^2 = Γ_base^2
//                + Γ_Dopp^2
//                + Γ_Zeeman^2
//                + Γ_User^2
//                + Γ_Coh^2
//
//       Γ_tot = sqrt(Γ_quad^2) + Γ_pair(B)
//
// Inputs:
//   Bval    : magnetic field [T]
//   Lch     : channel chord length [m]
//   Wch     : effective magnetic width [m]
//   SNS_p   : SNS_Params structure
//
// Reads from root:SNS_Settings:
//   SNS_GammaBase_eV
//   SNS_useGammaDopp,   SNS_GammaDoppAlpha
//   SNS_useGammaZeeman, SNS_GammaZeemanScale
//   SNS_useGammaPair,   SNS_GammaPairScale
//   SNS_useGammaUser
//   SNS_useGammaCoh
//
// Returns:
//   Γ_tot [eV]
//============================================================
Function SNS_ComputeGammaTot(Bval, Lch, Wch, SNS_p)
    Variable Bval, Lch, Wch
    STRUCT SNS_Params &SNS_p

    String dfOld = GetDataFolder(1)
    SetDataFolder root:SNS_Settings

    NVAR SNS_GammaBase_eV
    NVAR SNS_useGammaDopp,   SNS_GammaDoppAlpha
    NVAR SNS_useGammaZeeman, SNS_GammaZeemanScale
    NVAR SNS_useGammaPair,   SNS_GammaPairScale
    NVAR SNS_useGammaUser,   SNS_useGammaCoh

    SetDataFolder $dfOld

    // --- base width: fixed by zero-field LDOS fit ---
    Variable Gamma_base
    if (numtype(SNS_GammaBase_eV) == 0 && SNS_GammaBase_eV > 0)
        Gamma_base = SNS_GammaBase_eV
    else
        Gamma_base = SNS_p.Broadening
    endif

    // --- Doppler-related diagnostic broadening ---
    Variable Gamma_Dopp = 0
    if (SNS_useGammaDopp)
        Variable betaOrb = SNS_beta(Bval, Lch, Wch, SNS_p.lambdaL)
        Gamma_Dopp = SNS_GammaDoppAlpha * SNS_p.Delta * abs(betaOrb)
    endif

    // --- Zeeman broadening ---
    Variable Gamma_Zeeman = 0
    if (SNS_useGammaZeeman)
        Gamma_Zeeman = SNS_GammaZeemanScale * abs(Bval)
    endif

    // --- Quadratic depairing broadening ---
    // SNS_GammaPairScale has units eV/T^2.
    // This is additive on top of the zero-field broadening.
    Variable Gamma_Pair = 0
    if (SNS_useGammaPair)
        Gamma_Pair = SNS_GammaPairScale * Bval * Bval
    endif

    // --- User-defined broadening placeholder ---
    Variable Gamma_User = 0
    if (SNS_useGammaUser)
        Gamma_User = 0    // extend later
    endif

    // --- Coherence-length broadening diagnostic ---
    Variable Gamma_Coh = 0
    if (SNS_useGammaCoh && SNS_p.Delta > 0 && SNS_p.vF > 0)
        Variable xi = 5 * HBAR_eVs * SNS_p.vF / (pi * SNS_p.Delta)
        Gamma_Coh = (SNS_p.Delta/pi) * (1 - exp(-Lch/xi))
    endif

    // --- combine ---
    Variable GammaQuad
    GammaQuad = sqrt( Gamma_base*Gamma_base \
                    + Gamma_Dopp*Gamma_Dopp \
                    + Gamma_Zeeman*Gamma_Zeeman \
                    + Gamma_User*Gamma_User \
                    + Gamma_Coh*Gamma_Coh )

    Variable GammaTot
    GammaTot = GammaQuad + Gamma_Pair

    return GammaTot
End


//==============================================================================
// SNS_DynesDOS
//
// Scalar Dynes quasiparticle density of states. All three energy arguments
// must use the same units.
//==============================================================================
Function SNS_DynesDOS(energy, delta, gamma)
    Variable energy, delta, gamma

    Complex energyGamma = cmplx(abs(energy), gamma)
    Complex denominator = energyGamma*energyGamma - delta*delta
    return real(energyGamma/sqrt(denominator))
End


//==============================================================================
// SNS_ComputeDOS_FromChannels
//
// Purpose:
//   Compute DOS(E,B) in absolute units [states/eV] from a set of
//   precomputed SNS channels.
//
//   The function:
//     • Uses global SNS_Settings via SNS_LoadParams()
//       (Delta, vF, lambdaL, Broadening, NE).
//     • Uses SNS_Settings for vortex flags (SNS_useVortex, SNS_nFlux).
//     • Takes the magnetic-field grid explicitly (B_T).
//     • Adds a constant continuum DOS above |E| > Delta using
//       N_cont_statesPer_eV.
//
//   Legacy mode:
//     each ABS branch contributes a Lorentzian with the full channel weight.
//   branchNormalize mode:
//     all ABS branches generated from one channel share that channel weight,
//     so sum(branch weights for channel j) = normalized wChan[j].
//   Channel weights are absolute (not normalized):
//       weight_j = 2 * wChan[j]
//   → In 2D:     wChan[j] = 1
//   → In 3D:     wChan[j] = sin(theta_j) (or appropriate angular weight)
//
// Inputs:
//   B_T                   : magnetic field axis [T]
//   L_N_List_nm           : channel chord lengths [nm]
//   W_eff_List_nm         : effective magnetic widths [nm]
//   wChan                 : geometric channel weight (dimensionless)
//   T_eff_List            : normal-state transparency per channel (0…1)
//   rS1x_nm, rS1y_nm      : first S-contact hit point per channel [nm]
//   rS2x_nm, rS2y_nm      : second S-contact hit point per channel [nm]
//   xV_nm, yV_nm          : vortex center position [nm]
//   N_cont_statesPer_eV   : continuum DOS for this region [states/eV]
//   nameDOS               : output DOS wave name (2D: E × B)
//   nameEaxis             : output energy axis wave name [eV]
//
// Optional Inputs:
//   betaExtra_List        : precomputed per-channel extra phase [rad]
//                           Accepted forms:
//                             • 1D wave betaExtra_List[ch]
//                               -> field-independent extra phase
//                             • 2D wave betaExtra_List[ch][iB]
//                               -> field-dependent extra phase
//                           If supplied, it overrides the legacy point-vortex
//                           endpoint formula.
//                           If omitted, betaExtra is computed from xV_nm, yV_nm
//                           when SNS_useVortex != 0.
//
//   is3D                  : flag used only for cleanup/archive naming.
//                           If nonzero, branch waves are moved into
//                           E_allBranches_3D, m_allBranches_3D,
//                           and s_allBranches_3D.
//
//   Delta_fit_eV          : optional effective gap used for generated fit
//                           spectra [eV].
//                           If supplied, this value replaces SNS_Settings:Delta
//                           inside this call only. This is intended for fitting,
//                           where Delta_eff is varied by the optimizer while
//                           preserving the established DOS construction and
//                           normalization of SNS_ComputeDOS_FromChannels.
//                           If omitted, SNS_Settings:Delta is used.
//
//   branchNormalize       : optional diagnostic/physical sum-rule flag.
//                           0 or omitted = legacy per-branch channel weight.
//                           nonzero = divide channel weight by nBr so all ABS
//                           branches corresponding to a channel/wavefunction
//                           carry one channel's total spectral weight.
//
// Outputs:
//   nameDOS               : DOS(E,B) in [states/eV]
//   nameDOS + "_Gamma"    : DOS(E,B) with field/coherence broadening
//   nameEaxis             : energy axis [eV]
//
// Returns:
//   Nch                   : number of channels processed
//
// Notes:
//   • Subgap ABS energies are obtained from
//       Solve_AllBranches_SNS_dGSJ_betaextra().
//   • Continuum DOS is modeled as constant above |E| > Delta.
//   • No A(x,y) or Delta(x,y) maps are used (phase-winding vortex model).
//==============================================================================
Function SNS_ComputeDOS_FromChannels(B_T, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, rS1x_nm, rS1y_nm, rS2x_nm, rS2y_nm, xV_nm, yV_nm, N_cont_statesPer_eV, nameDOS, nameEaxis, [betaExtra_List, is3D, Delta_fit_eV, branchNormalize])
    Wave    B_T, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List
    Wave    rS1x_nm, rS1y_nm, rS2x_nm, rS2y_nm
    Variable xV_nm, yV_nm
    Variable N_cont_statesPer_eV
    String  nameDOS, nameEaxis
    Wave betaExtra_List
    Variable is3D
    Variable Delta_fit_eV
    Variable branchNormalize

    String currentDF = GetDataFolder(1)

    // --- log: start ---
    SNS_Log("SNS_ComputeDOS_FromChannels: start; DOS=" + nameDOS + \
            ", Eaxis=" + nameEaxis, level="INFO")

    // --- load global SNS parameters ---
    STRUCT SNS_Params SNS_p
    SNS_LoadParams(SNS_p)

    SetDataFolder root:SNS_Settings
    NVAR SNS_useVortex = SNS_useVortex
    NVAR SNS_nFlux     = SNS_nFlux
    SetDataFolder $currentDF

    Variable useVortex = SNS_useVortex
    Variable nFlux     = SNS_nFlux

    Variable Delta = SNS_p.Delta
    Variable useBranchNormalize = 0

    if (!ParamIsDefault(branchNormalize))
        useBranchNormalize = (branchNormalize != 0)
    endif

	if (!ParamIsDefault(Delta_fit_eV))
	    if (numtype(Delta_fit_eV) != 0 || Delta_fit_eV <= 0)
	        Abort "SNS_ComputeDOS_FromChannels: invalid Delta_fit_eV."
	    endif
	    Delta = Delta_fit_eV
	endif
    Variable vF         = SNS_p.vF
    Variable lambdaL    = SNS_p.lambdaL
    Variable Broadening = SNS_p.Broadening
    Variable NE         = SNS_p.NE

    // ---- Delta(B) parameters ----
    Variable SNS_DeltaBmax_T   = 0.5
    Variable SNS_DeltaFracDrop = 0.001

    Variable nB  = numpnts(B_T)
    Variable Nch = numpnts(L_N_List_nm)

    if ((Nch <= 0) || (Nch != numpnts(W_eff_List_nm)) || (Nch != numpnts(wChan)) \
        || (Nch != numpnts(T_eff_List)) \
        || (Nch != numpnts(rS1x_nm)) || (Nch != numpnts(rS1y_nm)) \
        || (Nch != numpnts(rS2x_nm)) || (Nch != numpnts(rS2y_nm)))

        String msg1 = "SNS_ComputeDOS_FromChannels: channel waves inconsistent (Nch=" \
                      + num2str(Nch) + ", nW=" + num2str(numpnts(W_eff_List_nm)) + ")"
        SNS_Log(msg1, level="ERR")
        Abort msg1
    endif

    if (nB <= 0)
        String msg2 = "SNS_ComputeDOS_FromChannels: B_T has no points."
        SNS_Log(msg2, level="ERR")
        Abort msg2
    endif

    // ---- basic sanity on continuum reference ----
    if (N_cont_statesPer_eV < 0)
        String msg3 = "SNS_ComputeDOS_FromChannels: N_cont_statesPer_eV < 0 (" + \
                      num2str(N_cont_statesPer_eV) + ")"
        SNS_Log(msg3, level="ERR")
        Abort msg3
    endif
    if (N_cont_statesPer_eV == 0)
        SNS_Log("SNS_ComputeDOS_FromChannels: N_cont_statesPer_eV = 0 → continuum normalization disabled; ABS only.",\
                level="WARN")
    endif

    // ---- normalize geometric weights only ----
    Variable sumW = sum(wChan)
    if (sumW <= 0)
        String msg4 = "SNS_ComputeDOS_FromChannels: all channel weights are zero."
        SNS_Log(msg4, level="ERR")
        Abort msg4
    endif

    // ---------------------------
    // 2. Precompute Delta(B) scaling only
    // ---------------------------
    Make/FREE/D/N=(nB) deltaScale_wv
    Variable iB, Bval, bDimDelta, deltaScale

    for (iB = 0; iB < nB; iB += 1)
        Bval = B_T[iB]

        // Delta(B) (global scale factor)
        if (SNS_DeltaBmax_T > 0)
            bDimDelta = abs(Bval)/SNS_DeltaBmax_T
        else
            bDimDelta = 0
        endif
        deltaScale = 1 - SNS_DeltaFracDrop*bDimDelta*bDimDelta
        if (deltaScale < 0)
            deltaScale = 0
        endif
        deltaScale_wv[iB] = deltaScale
    endfor

    //================================================
    // 3. Energy axis and DOS wave
    //================================================
    Variable dEnergy = 2.5e-3

    Make/O/D/N=(NE) $nameEaxis
    Wave E_axis = $nameEaxis
    // Legacy: limited to subgap energies
    // SetScale/I x, -Delta, Delta, "eV", E_axis
    SetScale/I x, -dEnergy, dEnergy, "eV", E_axis
    E_axis = x

    Make/O/D/N=(NE, nB) $nameDOS
    Wave DOS_EB = $nameDOS
    DOS_EB = 0

    String nameDOS_broadening = nameDOS + "_Gamma"
    Make/O/D/N=(NE, nB) $nameDOS_broadening
    Wave DOS_EB_broadening = $nameDOS_broadening
    DOS_EB_broadening = 0

    SetScale/P x, DimOffset(E_axis,0), DimDelta(E_axis,0), "eV", DOS_EB
    SetScale/P y, DimOffset(B_T,0),    DimDelta(B_T,0),    "T",  DOS_EB
    SetScale/P x, DimOffset(E_axis,0), DimDelta(E_axis,0), "eV", DOS_EB_broadening
    SetScale/P y, DimOffset(B_T,0),    DimDelta(B_T,0),    "T",  DOS_EB_broadening

    //================================================
    // 4. Channel loop + ABS solver + DOS
    //================================================
    Variable k, iE, j
    Variable Lch_nm, Wch_nm, Lch, Wch, Teff, weight, branchWeight
    Variable nBr, E0, E0_scaled, dE
    Variable GammaB_loc, GammaC_loc, GammaTot, deltaScale_loc
    Variable betaExtra, betaChord
    Variable DeltaCh, DeltaUsed

    String nameE2D, nameM, nameS

    SNS_Log("SNS_ComputeDOS_FromChannels: starting channel loop; Nch=" + num2str(Nch) + \
            ", nB=" + num2str(nB) + ", NE=" + num2str(NE), level="INFO")

	for (j = 0; j < Nch; j += 1)
	
	    Lch_nm = L_N_List_nm[j]
	    Wch_nm = W_eff_List_nm[j]
	    Lch    = Lch_nm * 1e-9
	    Wch    = Wch_nm * 1e-9
	    Teff = T_eff_List[j]
	    weight = wChan[j] / sumW
	
	    if ((weight <= 0) || numtype(weight) || numtype(Lch_nm) || numtype(Wch_nm) || numtype(Teff))
	        SNS_Log("SNS_ComputeDOS_FromChannels: skipping channel " + num2str(j) + \
	                " due to NaN/Inf or zero weight.", level="DBG")
	        continue
	    endif
	
	    if (Lch_nm <= 0 || Wch_nm < 0 || Teff < 0 || Teff > 1)
	        SNS_Log("SNS_ComputeDOS_FromChannels: skipping invalid channel " + num2str(j) + \
	                " L_nm=" + num2str(Lch_nm) + \
	                " W_nm=" + num2str(Wch_nm) + \
	                " T=" + num2str(Teff), level="WARN")
	        continue
	    endif

        // --- per-channel effective Delta (no Δ-map) ---
        DeltaCh   = Delta
        DeltaUsed = DeltaCh

        // --- per-channel extra phase from vortex phase winding at endpoints ---
        Variable useBetaExtraWave = 0
        Variable betaExtraIs2D = 0
        Variable betaExtraLegacy = 0
        Variable betaExtraUse

        if (!ParamIsDefault(betaExtra_List))

            useBetaExtraWave = 1

            if (WaveDims(betaExtra_List) == 1)

                if (numpnts(betaExtra_List) != Nch)
                    Abort "SNS_ComputeDOS_FromChannels: 1D betaExtra_List has wrong length."
                endif
                betaExtraIs2D = 0

            elseif (WaveDims(betaExtra_List) == 2)

                if ((DimSize(betaExtra_List, 0) != Nch) || (DimSize(betaExtra_List, 1) != nB))
                    Abort "SNS_ComputeDOS_FromChannels: 2D betaExtra_List has wrong dimensions."
                endif
                betaExtraIs2D = 1

            else
                Abort "SNS_ComputeDOS_FromChannels: betaExtra_List must be 1D or 2D."
            endif

        elseif (useVortex)

            Variable th1 = atan2(rS1y_nm[j] - yV_nm, rS1x_nm[j] - xV_nm)
            Variable th2 = atan2(rS2y_nm[j] - yV_nm, rS2x_nm[j] - xV_nm)
            Variable dth = th2 - th1

            // wrap to [-pi, pi]
            if (dth > pi)
                dth -= 2*pi
            elseif (dth < -pi)
                dth += 2*pi
            endif

            betaExtraLegacy = nFlux * dth
        endif

        sprintf nameE2D, "E_allBranches_ch%03d", j
        sprintf nameM,   "m_allBranches_ch%03d", j
        sprintf nameS,   "s_allBranches_ch%03d", j

        if (!useBetaExtraWave || !betaExtraIs2D)

            // ------------------------------------------
            // Case 1:
            //   no betaExtra_List      -> legacy betaExtraLegacy
            //   1D betaExtra_List[ch]  -> field-independent betaExtra
            //   one solve for full B_T
            // ------------------------------------------
            if (useBetaExtraWave)
                betaExtraUse = betaExtra_List[j]
                if (numtype(betaExtraUse) != 0)
                    continue
                endif
            else
                betaExtraUse = betaExtraLegacy
            endif

            nBr = Solve_AllBranches_SNS_dGSJ_betaextra( \
                        B_T, Lch, Wch, DeltaUsed, vF, lambdaL, Teff, \
                        betaExtraUse, nameE2D, nameM, nameS)

            if (nBr <= 0)
                SNS_Log("SNS_ComputeDOS_FromChannels: channel " + num2str(j) + \
                        " returned no branches (nBr=" + num2str(nBr) + ").", level="DBG")
                continue
            endif
            branchWeight = weight
            if (useBranchNormalize)
                branchWeight = weight / nBr
            endif

            Wave E_all = $nameE2D

            for (iB = 0; iB < nB; iB += 1)

                Bval           = B_T[iB]
                deltaScale_loc = deltaScale_wv[iB]

                // --- total broadening from modular Gamma block ---
                GammaTot = SNS_ComputeGammaTot(Bval, Lch, Wch, SNS_p)

                for (k = 0; k < nBr; k += 1)

                    E0 = E_all[iB][k]
                    if (numtype(E0))
                        continue
                    endif

                    E0_scaled = E0 * deltaScale_loc

                    for (iE = 0; iE < NE; iE += 1)
                        dE = E_axis[iE] - E0_scaled

                        // Reference DOS: fixed base Broadening only
                        DOS_EB[iE][iB] += branchWeight * (Broadening/pi) / (dE*dE + Broadening*Broadening)

                        // DOS with all field/lifetime broadening mechanisms
                        DOS_EB_broadening[iE][iB] += branchWeight * (GammaTot/pi) / (dE*dE + GammaTot*GammaTot)
                    endfor

                endfor
            endfor

        else

            // ------------------------------------------
            // Case 2:
            //   2D betaExtra_List[ch][iB]
            //   solve one B-point at a time
            // ------------------------------------------
            Make/FREE/D/N=1 B_one

            for (iB = 0; iB < nB; iB += 1)

                betaExtraUse = betaExtra_List[j][iB]
                if (numtype(betaExtraUse) != 0)
                    continue
                endif

                Bval = B_T[iB]
                B_one[0] = Bval
                deltaScale_loc = deltaScale_wv[iB]

                // unique temp names per channel/B point
                sprintf nameE2D, "E_allBranches_ch%03d_B%03d", j, iB
                sprintf nameM,   "m_allBranches_ch%03d_B%03d", j, iB
                sprintf nameS,   "s_allBranches_ch%03d_B%03d", j, iB

                nBr = Solve_AllBranches_SNS_dGSJ_betaextra( \
                            B_one, Lch, Wch, DeltaUsed, vF, lambdaL, Teff, \
                            betaExtraUse, nameE2D, nameM, nameS)

                if (nBr <= 0)
                    continue
                endif
                branchWeight = weight
                if (useBranchNormalize)
                    branchWeight = weight / nBr
                endif

                Wave E_all_1B = $nameE2D

                // --- total broadening from modular Gamma block ---
                GammaTot = SNS_ComputeGammaTot(Bval, Lch, Wch, SNS_p)

                for (k = 0; k < nBr; k += 1)

                    E0 = E_all_1B[0][k]
                    if (numtype(E0))
                        continue
                    endif

                    E0_scaled = E0 * deltaScale_loc

                    for (iE = 0; iE < NE; iE += 1)
                        dE = E_axis[iE] - E0_scaled

                        // Reference DOS: fixed base Broadening only
                        DOS_EB[iE][iB] += branchWeight * (Broadening/pi) / (dE*dE + Broadening*Broadening)

                        // DOS with all field/lifetime broadening mechanisms
                        DOS_EB_broadening[iE][iB] += branchWeight * (GammaTot/pi) / (dE*dE + GammaTot*GammaTot)
                    endfor

                endfor
            endfor

        endif
    endfor
    //================================================
    // 5. Normalize ABS DOS so that the *total* ABS
    //    weight in the full window [Emin,Emax]
    //    equals the normal-state continuum weight
    //    in the subgap region (|E|<Delta, clipped
    //    by the energy window).
    //
    //    Then add a continuum background with a tanh
    //    onset at |E| ~ Delta and amplitude
    //    N_cont_statesPer_eV.
    //================================================

    Variable NE_local = DimSize(E_axis, 0)
    Variable dEgrid
    Variable E_sm = 3 * Broadening   // [eV] broadening of gap edge

    if (NE_local < 2)
        String msg6 = "SNS_ComputeDOS_FromChannels: NE_local < 2; invalid energy axis."
        SNS_Log(msg6, level="ERR")
        Abort msg6
    endif

    dEgrid = E_axis[1] - E_axis[0]
    if (dEgrid <= 0 || numtype(dEgrid))
        String msg7 = "SNS_ComputeDOS_FromChannels: invalid energy step dE=" + num2str(dEgrid)
        SNS_Log(msg7, level="ERR")
        Abort msg7
    endif

    Variable Emin = E_axis[0]
    Variable Emax = E_axis[NE_local - 1]
    Variable Ewin = Emax - Emin

    if (N_cont_statesPer_eV <= 0)
        SNS_Log("SNS_ComputeDOS_FromChannels: N_cont_statesPer_eV <= 0; " + \
                "skipping ABS normalization and continuum background.", \
                level="WARN")
    else

        // ---------- 5a. Define subgap interval and target gap weight ----------

        // subgap window inside the current E window: [Egap_low, Egap_high]
        Variable Egap_low  = max(Emin, -Delta)
        Variable Egap_high = min(Emax,  Delta)
        Variable EgapSpan  = Egap_high - Egap_low

        if (EgapSpan <= 0)
            String msg8 = "SNS_ComputeDOS_FromChannels: energy window does not overlap |E|<Delta; " + \
                          "ABS normalization to subgap continuum skipped."
            SNS_Log(msg8, level="WARN")
        else
            // normal-state continuum spectral weight that would lie in |E|<Delta
            Variable N_gap_ref = N_cont_statesPer_eV * EgapSpan

            // ---------- 5b. ABS weight over full window & scaling per B ----------

            Make/FREE/D/N=(nB) N_abs_B, alpha_B
            Variable acc, alpha
            Variable nAlphaOver = 0
            Variable nAlphaUnder = 0
            Variable firstAlphaOver = NaN, lastAlphaOver = NaN
            Variable firstAlphaUnder = NaN, lastAlphaUnder = NaN
            Variable alphaOverMin = Inf, alphaOverMax = -Inf
            Variable alphaUnderMin = Inf, alphaUnderMax = -Inf

            for (iB = 0; iB < nB; iB += 1)

                // integrate ABS over the *entire* [Emin,Emax] window
                acc = 0
                for (iE = 0; iE < NE_local; iE += 1)
                    acc += DOS_EB[iE][iB]
                endfor
                N_abs_B[iB] = acc * dEgrid

                if (N_abs_B[iB] <= 0 || numtype(N_abs_B[iB]))
                    // no ABS weight in window: nothing to normalize
                    alpha_B[iB] = 0
                    SNS_Log("SNS_ComputeDOS_FromChannels: no ABS weight in window at B index " + \
                            num2str(iB) + "; skipping ABS normalization for this slice.", \
                            level="DBG")
                else
                    // scale total ABS weight in [Emin,Emax] to match
                    // continuum gap weight N_gap_ref
                    alpha_B[iB] = N_gap_ref / N_abs_B[iB]

                    // log significant over-/under-counting
                    if (alpha_B[iB] < 0.5)
                        nAlphaOver += 1
                        if (numtype(firstAlphaOver) != 0)
                            firstAlphaOver = iB
                        endif
                        lastAlphaOver = iB
                        alphaOverMin = min(alphaOverMin, alpha_B[iB])
                        alphaOverMax = max(alphaOverMax, alpha_B[iB])
                    elseif (alpha_B[iB] > 2)
                        nAlphaUnder += 1
                        if (numtype(firstAlphaUnder) != 0)
                            firstAlphaUnder = iB
                        endif
                        lastAlphaUnder = iB
                        alphaUnderMin = min(alphaUnderMin, alpha_B[iB])
                        alphaUnderMax = max(alphaUnderMax, alpha_B[iB])
                    else
                        SNS_Log("SNS_ComputeDOS_FromChannels: ABS normalized at B index " + \
                                num2str(iB) + " with alpha=" + num2str(alpha_B[iB]) + \
                                " to match subgap continuum weight.", \
                                level="DBG")
                    endif
                endif
            endfor

            if (nAlphaOver > 0)
                SNS_Log("SNS_ComputeDOS_FromChannels: ABS overcount normalization in " + \
                        num2str(nAlphaOver) + "/" + num2str(nB) + \
                        " B slices; alpha range=[" + num2str(alphaOverMin) + "," + \
                        num2str(alphaOverMax) + "], first/last B index=" + \
                        num2str(firstAlphaOver) + "/" + num2str(lastAlphaOver) + ".", \
                        level="WARN")
            endif
            if (nAlphaUnder > 0)
                SNS_Log("SNS_ComputeDOS_FromChannels: ABS undercount normalization in " + \
                        num2str(nAlphaUnder) + "/" + num2str(nB) + \
                        " B slices; alpha range=[" + num2str(alphaUnderMin) + "," + \
                        num2str(alphaUnderMax) + "], first/last B index=" + \
                        num2str(firstAlphaUnder) + "/" + num2str(lastAlphaUnder) + ".", \
                        level="WARN")
            endif

            // Keep alpha_B as an internal/free wave for now. If detailed
            // per-B normalization diagnostics become important, export a named
            // alpha_B wave next to nameDOS rather than logging every slice.

            // apply ABS scaling (whole slice)
            for (iB = 0; iB < nB; iB += 1)
                alpha = alpha_B[iB]
                if (alpha <= 0 || numtype(alpha))
                    continue        // leave ABS as computed
                endif

                for (iE = 0; iE < NE_local; iE += 1)
                    DOS_EB[iE][iB]            *= alpha
                    DOS_EB_broadening[iE][iB] *= alpha
                endfor
            endfor

            SNS_Log("SNS_ComputeDOS_FromChannels: ABS DOS normalized so that total ABS weight in " + \
                    "[" + num2str(Emin) + "," + num2str(Emax) + "] equals subgap continuum weight " + \
                    "N_gap_ref=" + num2str(N_gap_ref) + " states per B-slice.", level="INFO")
        endif

        // ---------- 5c. Add continuum background via tanh profile ----------

        // continuum shape f(E); tanh onset at |E| ~ Delta
        Make/FREE/D/N=(NE_local) contShape
        Variable p, E
        for (p = 0; p < NE_local; p += 1)
            E = E_axis[p]
            contShape[p] = 0.5 * (1 + tanh( (abs(E) - Delta)/E_sm ))
        endfor

        for (iB = 0; iB < nB; iB += 1)
            for (iE = 0; iE < NE_local; iE += 1)
                Variable fE = contShape[iE]
                DOS_EB[iE][iB]            += N_cont_statesPer_eV * fE
                DOS_EB_broadening[iE][iB] += N_cont_statesPer_eV * fE
            endfor
        endfor

        SNS_Log("SNS_ComputeDOS_FromChannels: added continuum background with N_cont_statesPer_eV=" + \
                num2str(N_cont_statesPer_eV) + " using tanh(|E|-Delta)/E_sm profile.", \
                level="INFO")
    endif

    //================================================
    // 6. Cleanup: move branch waves into subfolders
    //================================================
    String brSuffix = ""
    if (!ParamIsDefault(is3D) && is3D)
        brSuffix = "_3D"
    endif

    String folderE = "E_allBranches" + brSuffix
    String folderM = "m_allBranches" + brSuffix
    String folderS = "s_allBranches" + brSuffix

    NewDataFolder/O $folderE
    NewDataFolder/O $folderM
    NewDataFolder/O $folderS

    MoveBranchWavesToSubfolder("E_allBranches_", folderE)
    MoveBranchWavesToSubfolder("m_allBranches_", folderM)
    MoveBranchWavesToSubfolder("s_allBranches_", folderS)

    SNS_Log("SNS_ComputeDOS_FromChannels: finished successfully; Nch=" + num2str(Nch), level="INFO")

    return Nch
End

//==============================================================================
// Green-function 1D SNS helper block
//
// This is a separate finite-chain SNS solver.  It is intentionally independent
// of the branch/Lorentzian ABS construction above: the LDOS is obtained from the
// retarded Green function of a short 1D envelope chain with BCS reservoir
// self-energies attached to both ends.
//==============================================================================
Function SNS_Green1D_AddComplexElement(M, rowC, colC, a)
    Wave M
    Variable rowC, colC
    Complex a

    Variable rowR = 2*rowC
    Variable colR = 2*colC

    M[rowR][colR]     += real(a)
    M[rowR][colR+1]   += -imag(a)
    M[rowR+1][colR]   += imag(a)
    M[rowR+1][colR+1] += real(a)
End

Function/C SNS_Green1D_BCSDenom(E_eV, eta_eV, Delta_eV)
    Variable E_eV, eta_eV, Delta_eV

    Complex zE = cmplx(E_eV, eta_eV)
    Complex denom = sqrt(cmplx(Delta_eV*Delta_eV, 0) - zE*zE)

    if ((abs(E_eV) < Delta_eV) && (real(denom) < 0))
        denom = -denom
    endif

    if ((abs(E_eV) > Delta_eV) && (E_eV > 0) && (imag(denom) > 0))
        denom = -denom
    endif
    if ((abs(E_eV) > Delta_eV) && (E_eV < 0) && (imag(denom) < 0))
        denom = -denom
    endif

    return denom
End

Function SNS_Green1D_AddBCSSelfEnergy(A, site, E_eV, etaS_eV, Delta_eV, Gamma_eV, phi_rad)
    Wave A
    Variable site
    Variable E_eV, etaS_eV, Delta_eV, Gamma_eV, phi_rad

    Complex zE = cmplx(E_eV, etaS_eV)
    Complex denom = SNS_Green1D_BCSDenom(E_eV, etaS_eV, Delta_eV)
    Complex phaseP = cmplx(cos(phi_rad), sin(phi_rad))
    Complex phaseM = cmplx(cos(phi_rad), -sin(phi_rad))

    Complex sigmaDiag = -Gamma_eV * zE / denom
    Complex sigmaEH   = -Gamma_eV * Delta_eV * phaseP / denom
    Complex sigmaHE   = -Gamma_eV * Delta_eV * phaseM / denom

    Variable eIdx = 2*site
    Variable hIdx = eIdx + 1

    // A = zI - H - Sigma.
    SNS_Green1D_AddComplexElement(A, eIdx, eIdx, -sigmaDiag)
    SNS_Green1D_AddComplexElement(A, hIdx, hIdx, -sigmaDiag)
    SNS_Green1D_AddComplexElement(A, eIdx, hIdx, -sigmaEH)
    SNS_Green1D_AddComplexElement(A, hIdx, eIdx, -sigmaHE)
End

Function SNS_Green1D_FillMatrix(A, Nsite, t_eV, E_eV, etaN_eV, etaS_eV, Delta_eV, Gamma_eV, beta_rad)
    Wave A
    Variable Nsite, t_eV
    Variable E_eV, etaN_eV, etaS_eV, Delta_eV, Gamma_eV, beta_rad

    A = 0

    Complex zN = cmplx(E_eV, etaN_eV)
    Variable i
    Variable eIdx, hIdx

    for (i = 0; i < Nsite; i += 1)
        eIdx = 2*i
        hIdx = eIdx + 1

        SNS_Green1D_AddComplexElement(A, eIdx, eIdx, zN)
        SNS_Green1D_AddComplexElement(A, hIdx, hIdx, zN)

        if (i < Nsite-1)
            // Envelope tight-binding at band center:
            // electron H hopping = -t, hole H hopping = +t.
            // Therefore A = zI - H contributes +t and -t respectively.
            SNS_Green1D_AddComplexElement(A, eIdx,   eIdx+2, cmplx(t_eV, 0))
            SNS_Green1D_AddComplexElement(A, eIdx+2, eIdx,   cmplx(t_eV, 0))
            SNS_Green1D_AddComplexElement(A, hIdx,   hIdx+2, cmplx(-t_eV, 0))
            SNS_Green1D_AddComplexElement(A, hIdx+2, hIdx,   cmplx(-t_eV, 0))
        endif
    endfor

    SNS_Green1D_AddBCSSelfEnergy(A, 0,       E_eV, etaS_eV, Delta_eV, Gamma_eV, -0.5*beta_rad)
    SNS_Green1D_AddBCSSelfEnergy(A, Nsite-1, E_eV, etaS_eV, Delta_eV, Gamma_eV,  0.5*beta_rad)
End

Function SNS_Green1D_ChannelLDOS(E_eV, beta_rad, Lch_nm, Teff, weight, Delta_eV, vF_mps, Gamma0_eV, etaN_eV, etaS_eV, NsiteMax, avgSites)
    Variable E_eV, beta_rad, Lch_nm, Teff, weight
    Variable Delta_eV, vF_mps, Gamma0_eV, etaN_eV, etaS_eV
    Variable NsiteMax, avgSites

    if ((Lch_nm <= 0) || (weight <= 0) || (Teff <= 0))
        return 0
    endif

    Variable Nsite = round(NsiteMax)
    if (Nsite < 3)
        Nsite = 3
    endif

    Variable Lch_m = Lch_nm * 1e-9
    Variable dx_m = Lch_m / (Nsite - 1)
    Variable t_eV = HBAR_eVs * vF_mps / (2*dx_m)
    Variable Gamma_eV = Gamma0_eV * Teff

    Variable dimC = 2*Nsite
    Variable dimR = 2*dimC

    Make/FREE/D/N=(dimR, dimR) M_Green1D_A
    Make/FREE/D/N=(dimR) w_Green1D_b

    SNS_Green1D_FillMatrix(M_Green1D_A, Nsite, t_eV, E_eV, etaN_eV, etaS_eV, Delta_eV, Gamma_eV, beta_rad)

    Variable nProbe
    if (avgSites <= 0)
        nProbe = Nsite
    else
        nProbe = round(avgSites)
        if (nProbe < 1)
            nProbe = 1
        endif
        if (nProbe > Nsite)
            nProbe = Nsite
        endif
    endif

    Variable ip, site, eIdx, rowR
    Variable ldos = 0
    Variable frac
    Complex Gdiag

    for (ip = 0; ip < nProbe; ip += 1)
        if (nProbe == 1)
            site = round(0.5*(Nsite-1))
        elseif (avgSites <= 0)
            site = ip
        else
            frac = ip / (nProbe - 1)
            site = round(frac * (Nsite - 1))
        endif

        eIdx = 2*site
        rowR = 2*eIdx

        w_Green1D_b = 0
        w_Green1D_b[rowR] = 1

        Duplicate/FREE M_Green1D_A, M_Green1D_A_work
        MatrixLinearSolve/O/Z M_Green1D_A_work, w_Green1D_b
        if (V_flag != 0)
            return 0
        endif

        Gdiag = cmplx(w_Green1D_b[rowR], w_Green1D_b[rowR+1])
        ldos += -imag(Gdiag) / pi
    endfor

    return weight * ldos / nProbe
End

//==============================================================================
// SNS_ComputeDOS_FromChannels_Green1D
//
// Purpose:
//   Compute LDOS(E,B) from a retarded 1D SNS Green-function model per channel.
//   This avoids explicit ABS root finding and Lorentzian placement.
//
// Notes:
//   The normal segment is a band-center envelope tight-binding chain with group
//   velocity vF.  Superconductors enter only through fixed BCS self-energies at
//   the two endpoints.  avgSites>0 samples that many sites for a cheap
//   length-average; avgSites<=0 uses all chain sites.
//==============================================================================
Function SNS_ComputeDOS_FromChannels_Green1D(B_T, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, nameDOS, nameEaxis, [betaExtra_List, Delta_fit_eV, Gamma0_eV, etaN_eV, etaS_eV, NsiteMax, avgSites])
    Wave B_T, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List
    String nameDOS, nameEaxis
    Wave betaExtra_List
    Variable Delta_fit_eV, Gamma0_eV, etaN_eV, etaS_eV, NsiteMax, avgSites

    STRUCT SNS_Params SNS_p
    SNS_LoadParams(SNS_p)

    Variable Delta = SNS_p.Delta
    if (!ParamIsDefault(Delta_fit_eV))
        Delta = Delta_fit_eV
    endif

    Variable Gamma0 = Delta
    if (!ParamIsDefault(Gamma0_eV))
        Gamma0 = Gamma0_eV
    endif

    Variable etaN = SNS_p.Broadening
    if (!ParamIsDefault(etaN_eV))
        etaN = etaN_eV
    endif

    Variable etaS = max(SNS_p.Broadening, 1e-9)
    if (!ParamIsDefault(etaS_eV))
        etaS = etaS_eV
    endif

    Variable NsiteUse = 17
    if (!ParamIsDefault(NsiteMax))
        NsiteUse = NsiteMax
    endif

    Variable avgSitesUse = 7
    if (!ParamIsDefault(avgSites))
        avgSitesUse = avgSites
    endif

    Variable nB = numpnts(B_T)
    Variable Nch = numpnts(L_N_List_nm)
    Variable NE = SNS_p.NE

    if ((Nch <= 0) || (Nch != numpnts(W_eff_List_nm)) || (Nch != numpnts(wChan)) || (Nch != numpnts(T_eff_List)))
        Abort "SNS_ComputeDOS_FromChannels_Green1D: channel waves inconsistent."
    endif
    if (nB <= 0)
        Abort "SNS_ComputeDOS_FromChannels_Green1D: B_T has no points."
    endif

    Variable betaExtraIs2D = 0
    Variable useBetaExtraWave = 0
    if (!ParamIsDefault(betaExtra_List))
        useBetaExtraWave = 1
        if (WaveDims(betaExtra_List) == 1)
            if (numpnts(betaExtra_List) != Nch)
                Abort "SNS_ComputeDOS_FromChannels_Green1D: 1D betaExtra_List has wrong length."
            endif
        elseif (WaveDims(betaExtra_List) == 2)
            if ((DimSize(betaExtra_List, 0) != Nch) || (DimSize(betaExtra_List, 1) != nB))
                Abort "SNS_ComputeDOS_FromChannels_Green1D: 2D betaExtra_List has wrong dimensions."
            endif
            betaExtraIs2D = 1
        else
            Abort "SNS_ComputeDOS_FromChannels_Green1D: betaExtra_List must be 1D or 2D."
        endif
    endif

    Variable sumW = sum(wChan)
    if (sumW <= 0)
        Abort "SNS_ComputeDOS_FromChannels_Green1D: all channel weights are zero."
    endif

    Variable dEnergy = 2.5e-3

    Make/O/D/N=(NE) $nameEaxis
    Wave E_axis = $nameEaxis
    SetScale/I x, -dEnergy, dEnergy, "eV", E_axis
    E_axis = x

    Make/O/D/N=(NE, nB) $nameDOS
    Wave DOS_EB = $nameDOS
    DOS_EB = 0
    SetScale/P x, DimOffset(E_axis,0), DimDelta(E_axis,0), "eV", DOS_EB
    SetScale/P y, DimOffset(B_T,0),    DimDelta(B_T,0),    "T",  DOS_EB

    Variable j, iB, iE
    Variable Lch_nm, Wch_nm, Teff, weight
    Variable betaExtraUse, betaChord, betaTot

    SNS_Log("SNS_ComputeDOS_FromChannels_Green1D: start; channels=" + num2str(Nch) + ", nB=" + num2str(nB) + ", NE=" + num2str(NE), level="INFO")

    for (j = 0; j < Nch; j += 1)
        Lch_nm = L_N_List_nm[j]
        Wch_nm = W_eff_List_nm[j]
        Teff = T_eff_List[j]
        weight = wChan[j] / sumW

        if ((weight <= 0) || (Lch_nm <= 0) || (Teff <= 0) || numtype(Lch_nm) || numtype(Wch_nm) || numtype(Teff) || numtype(weight))
            continue
        endif

        for (iB = 0; iB < nB; iB += 1)
            betaChord = SNS_beta(B_T[iB], Lch_nm*1e-9, Wch_nm*1e-9, SNS_p.lambdaL)

            betaExtraUse = 0
            if (useBetaExtraWave)
                if (betaExtraIs2D)
                    betaExtraUse = betaExtra_List[j][iB]
                else
                    betaExtraUse = betaExtra_List[j]
                endif
                if (numtype(betaExtraUse) != 0)
                    continue
                endif
            endif

            betaTot = betaChord + betaExtraUse

            for (iE = 0; iE < NE; iE += 1)
                DOS_EB[iE][iB] += SNS_Green1D_ChannelLDOS(E_axis[iE], betaTot, Lch_nm, Teff, weight, Delta, SNS_p.vF, Gamma0, etaN, etaS, NsiteUse, avgSitesUse)
            endfor
        endfor
    endfor

    String meta
    meta  = "SNS_Solver=Green1D;"
    meta += "SNS_Green1D_Gamma0_eV=" + num2str(Gamma0) + ";"
    meta += "SNS_Green1D_etaN_eV=" + num2str(etaN) + ";"
    meta += "SNS_Green1D_etaS_eV=" + num2str(etaS) + ";"
    meta += "SNS_Green1D_NsiteMax=" + num2str(NsiteUse) + ";"
    meta += "SNS_Green1D_avgSites=" + num2str(avgSitesUse) + ";"
    meta += "SNS_Green1D_note=band-center envelope chain with endpoint BCS self-energies;"

    Note/K DOS_EB
    Note DOS_EB, meta

    return Nch
End

//==============================================================================
// SNS_ComputeDOS_FromChannels_TwoPass
//
// Purpose:
//   Diagnostic two-pass DOS builder for phase-space-proximity broadening tests.
//
//   Pass 1 solves all channel branches and stores a flattened branch table:
//       Branch_E_eV[record][B]
//       Branch_ch, Branch_br, Branch_m, Branch_s
//       Branch_L_nm, Branch_W_nm, Branch_T, Branch_weight
//       Branch_phi_rad and branch endpoint coordinates
//
//   Pass 2 computes a local branch-neighborhood density
//
//       rho_i(B) = sum_j w_j K_E(E_i-E_j) K_W(W_i-W_j)
//
//   and accumulates three DOS waves:
//       nameDOS                 : fixed base Broadening
//       nameDOS + "_Gamma"      : current SNS_ComputeGammaTot broadening
//       nameDOS + "_Mix"        : optional Gamma_mix diagnostic broadening
//
// Inputs:
//   Same core channel inputs as SNS_ComputeDOS_FromChannels(...).
//
// Optional Inputs:
//   betaExtra_List       : accepted only as 1D betaExtra_List[ch].
//                          2D betaExtra_List[ch][B] is intentionally not
//                          supported in this first two-pass diagnostic because
//                          it breaks simple branch identity across B.
//   sigmaE_eV            : energy kernel width for rho. Default: Broadening.
//   sigmaW_nm            : W_eff kernel width for rho. Default: 25 nm.
//   GammaMixMax_eV       : maximum candidate mixing broadening.
//                          Default: 0, so nameDOS+"_Mix" equals Gamma path.
//   rho0                 : saturation density scale. Default: 1.
//   useMixBroadening     : if nonzero, add Gamma_mix in quadrature.
//                          Default: 0, diagnostic only.
//   outPrefix            : prefix for diagnostic branch/rho waves.
//                          Default: nameDOS + "_TwoPass".
//   phiList_rad          : optional compacted channel-angle list [rad].
//                          When supplied, Branch_phi_rad[rec] =
//                          phiList_rad[Branch_ch[rec]].
//
// Outputs:
//   nameDOS, nameDOS+"_Gamma", nameDOS+"_Mix", nameEaxis.
//   Diagnostic waves named from outPrefix:
//       <outPrefix>_Branch_E_eV
//       <outPrefix>_Branch_ch
//       <outPrefix>_Branch_br
//       <outPrefix>_Branch_m
//       <outPrefix>_Branch_s
//       <outPrefix>_Branch_L_nm
//       <outPrefix>_Branch_W_nm
//       <outPrefix>_Branch_T
//       <outPrefix>_Branch_weight
//       <outPrefix>_Branch_phi_rad
//       <outPrefix>_rho_branch
//       <outPrefix>_Gamma_mix_branch
//
// Returns:
//   Number of flattened branch records.
//
// Notes:
//   This function is deliberately separate from SNS_ComputeDOS_FromChannels.
//   It is a diagnostic/model-development path and should reproduce the current
//   DOS when useMixBroadening=0 and GammaMixMax_eV=0, up to ordinary branch
//   bookkeeping differences.
//==============================================================================
Function SNS_ComputeDOS_FromChannels_TwoPass(B_T, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, rS1x_nm, rS1y_nm, rS2x_nm, rS2y_nm, xV_nm, yV_nm, N_cont_statesPer_eV, nameDOS, nameEaxis, [betaExtra_List, is3D, Delta_fit_eV, sigmaE_eV, sigmaW_nm, GammaMixMax_eV, rho0, useMixBroadening, outPrefix, phiList_rad])
    Wave    B_T, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List
    Wave    rS1x_nm, rS1y_nm, rS2x_nm, rS2y_nm
    Variable xV_nm, yV_nm
    Variable N_cont_statesPer_eV
    String  nameDOS, nameEaxis
    Wave betaExtra_List
    Variable is3D, Delta_fit_eV
    Variable sigmaE_eV, sigmaW_nm, GammaMixMax_eV, rho0, useMixBroadening
    String outPrefix
    Wave phiList_rad

    String currentDF = GetDataFolder(1)

    STRUCT SNS_Params SNS_p
    SNS_LoadParams(SNS_p)

    SetDataFolder root:SNS_Settings
    NVAR SNS_useVortex = SNS_useVortex
    NVAR SNS_nFlux     = SNS_nFlux
    SetDataFolder $currentDF

    Variable Delta = SNS_p.Delta
    if (!ParamIsDefault(Delta_fit_eV))
        if (numtype(Delta_fit_eV) != 0 || Delta_fit_eV <= 0)
            Abort "SNS_ComputeDOS_FromChannels_TwoPass: invalid Delta_fit_eV."
        endif
        Delta = Delta_fit_eV
    endif

    Variable vF         = SNS_p.vF
    Variable lambdaL    = SNS_p.lambdaL
    Variable Broadening = SNS_p.Broadening
    Variable NE         = SNS_p.NE

    if (ParamIsDefault(sigmaE_eV) || sigmaE_eV <= 0 || numtype(sigmaE_eV))
        sigmaE_eV = Broadening
    endif
    if (sigmaE_eV <= 0 || numtype(sigmaE_eV))
        sigmaE_eV = 1e-9
    endif
    if (ParamIsDefault(sigmaW_nm) || sigmaW_nm <= 0 || numtype(sigmaW_nm))
        sigmaW_nm = 25
    endif
    if (ParamIsDefault(GammaMixMax_eV) || GammaMixMax_eV < 0 || numtype(GammaMixMax_eV))
        GammaMixMax_eV = 0
    endif
    if (ParamIsDefault(rho0) || rho0 <= 0 || numtype(rho0))
        rho0 = 1
    endif
    if (ParamIsDefault(useMixBroadening))
        useMixBroadening = 0
    endif
    if (ParamIsDefault(outPrefix) || strlen(outPrefix) <= 0)
        outPrefix = nameDOS + "_TwoPass"
    endif

    Variable nB  = numpnts(B_T)
    Variable Nch = numpnts(L_N_List_nm)

    Variable nW = numpnts(W_eff_List_nm)
    Variable nWeight = numpnts(wChan)
    Variable nT = numpnts(T_eff_List)
    Variable nX1 = numpnts(rS1x_nm)
    Variable nY1 = numpnts(rS1y_nm)
    Variable nX2 = numpnts(rS2x_nm)
    Variable nY2 = numpnts(rS2y_nm)

    if ((Nch <= 0) || (Nch != nW) || (Nch != nWeight) \
        || (Nch != nT) \
        || (Nch != nX1) || (Nch != nY1) \
        || (Nch != nX2) || (Nch != nY2))

        String msgInconsistent
        msgInconsistent = "SNS_ComputeDOS_FromChannels_TwoPass: channel waves inconsistent. " + \
            "nL=" + num2str(Nch) + \
            ", nW=" + num2str(nW) + \
            ", nwChan=" + num2str(nWeight) + \
            ", nT=" + num2str(nT) + \
            ", nHit1x=" + num2str(nX1) + \
            ", nHit1y=" + num2str(nY1) + \
            ", nHit2x=" + num2str(nX2) + \
            ", nHit2y=" + num2str(nY2) + "."
        SNS_Log(msgInconsistent, level="ERR")
        Abort msgInconsistent
    endif
    if (nB <= 0)
        Abort "SNS_ComputeDOS_FromChannels_TwoPass: B_T has no points."
    endif
    if (N_cont_statesPer_eV < 0)
        Abort "SNS_ComputeDOS_FromChannels_TwoPass: N_cont_statesPer_eV < 0."
    endif
    if (!ParamIsDefault(phiList_rad) && (numpnts(phiList_rad) != Nch))
        Abort "SNS_ComputeDOS_FromChannels_TwoPass: phiList_rad length mismatch."
    endif

    Variable useBetaExtraWave = 0
    if (!ParamIsDefault(betaExtra_List))
        useBetaExtraWave = 1
        if (WaveDims(betaExtra_List) != 1)
            Abort "SNS_ComputeDOS_FromChannels_TwoPass: betaExtra_List must be 1D in this diagnostic implementation."
        endif
        if (numpnts(betaExtra_List) != Nch)
            Abort "SNS_ComputeDOS_FromChannels_TwoPass: betaExtra_List length mismatch."
        endif
    endif

    Variable sumW = sum(wChan)
    if (sumW <= 0)
        Abort "SNS_ComputeDOS_FromChannels_TwoPass: all channel weights are zero."
    endif

    // ---- Delta(B) scaling ----
    Variable SNS_DeltaBmax_T   = 0.5
    Variable SNS_DeltaFracDrop = 0.001
    Make/FREE/D/N=(nB) deltaScale_wv
    Variable iB, Bval, bDimDelta, deltaScale
    for (iB = 0; iB < nB; iB += 1)
        Bval = B_T[iB]
        if (SNS_DeltaBmax_T > 0)
            bDimDelta = abs(Bval) / SNS_DeltaBmax_T
        else
            bDimDelta = 0
        endif
        deltaScale = 1 - SNS_DeltaFracDrop*bDimDelta*bDimDelta
        if (deltaScale < 0)
            deltaScale = 0
        endif
        deltaScale_wv[iB] = deltaScale
    endfor

    // ---- Energy axis and output DOS waves ----
    Variable dEnergy = 2.5e-3

    Make/O/D/N=(NE) $nameEaxis
    Wave E_axis = $nameEaxis
    SetScale/I x, -dEnergy, dEnergy, "eV", E_axis
    E_axis = x

    Make/O/D/N=(NE, nB) $nameDOS
    Wave DOS_EB = $nameDOS
    DOS_EB = 0

    String nameDOSGamma = nameDOS + "_Gamma"
    Make/O/D/N=(NE, nB) $nameDOSGamma
    Wave DOS_EB_Gamma = $nameDOSGamma
    DOS_EB_Gamma = 0

    String nameDOSMix = nameDOS + "_Mix"
    Make/O/D/N=(NE, nB) $nameDOSMix
    Wave DOS_EB_Mix = $nameDOSMix
    DOS_EB_Mix = 0

    SetScale/P x, DimOffset(E_axis,0), DimDelta(E_axis,0), "eV", DOS_EB
    SetScale/P y, DimOffset(B_T,0),    DimDelta(B_T,0),    "T",  DOS_EB
    SetScale/P x, DimOffset(E_axis,0), DimDelta(E_axis,0), "eV", DOS_EB_Gamma
    SetScale/P y, DimOffset(B_T,0),    DimDelta(B_T,0),    "T",  DOS_EB_Gamma
    SetScale/P x, DimOffset(E_axis,0), DimDelta(E_axis,0), "eV", DOS_EB_Mix
    SetScale/P y, DimOffset(B_T,0),    DimDelta(B_T,0),    "T",  DOS_EB_Mix

    // ---- Branch table waves ----
    String nameBranchE = outPrefix + "_Branch_E_eV"
    String nameBranchCh = outPrefix + "_Branch_ch"
    String nameBranchBr = outPrefix + "_Branch_br"
    String nameBranchM = outPrefix + "_Branch_m"
    String nameBranchS = outPrefix + "_Branch_s"
    String nameBranchL = outPrefix + "_Branch_L_nm"
    String nameBranchW = outPrefix + "_Branch_W_nm"
    String nameBranchT = outPrefix + "_Branch_T"
    String nameBranchWeight = outPrefix + "_Branch_weight"
    String nameBranchX1 = outPrefix + "_Branch_x1_nm"
    String nameBranchY1 = outPrefix + "_Branch_y1_nm"
    String nameBranchX2 = outPrefix + "_Branch_x2_nm"
    String nameBranchY2 = outPrefix + "_Branch_y2_nm"
    String nameBranchPhi = outPrefix + "_Branch_phi_rad"

    Make/O/D/N=(0, nB) $nameBranchE
    Make/O/D/N=0 $nameBranchCh, $nameBranchBr, $nameBranchM, $nameBranchS
    Make/O/D/N=0 $nameBranchL, $nameBranchW, $nameBranchT, $nameBranchWeight
    Make/O/D/N=0 $nameBranchX1, $nameBranchY1, $nameBranchX2, $nameBranchY2
    Make/O/D/N=0 $nameBranchPhi

    Wave Branch_E = $nameBranchE
    Wave Branch_ch = $nameBranchCh
    Wave Branch_br = $nameBranchBr
    Wave Branch_m = $nameBranchM
    Wave Branch_s = $nameBranchS
    Wave Branch_L = $nameBranchL
    Wave Branch_W = $nameBranchW
    Wave Branch_T = $nameBranchT
    Wave Branch_weight = $nameBranchWeight
    Wave Branch_x1 = $nameBranchX1
    Wave Branch_y1 = $nameBranchY1
    Wave Branch_x2 = $nameBranchX2
    Wave Branch_y2 = $nameBranchY2
    Wave Branch_phi = $nameBranchPhi

    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", Branch_E

    //==========================================================================
    // Pass 1: solve/store all full-B branches
    //==========================================================================
    Variable j, k, rec0, rec, nRec = 0
    Variable Lch_nm, Wch_nm, Lch, Wch, Teff, weight
    Variable betaExtraUse, betaExtraLegacy
    Variable nBr, E0
    String nameE2D, nameM, nameS

    for (j = 0; j < Nch; j += 1)

        Lch_nm = L_N_List_nm[j]
        Wch_nm = W_eff_List_nm[j]
        Lch = Lch_nm * 1e-9
        Wch = Wch_nm * 1e-9
        Teff = T_eff_List[j]
        weight = wChan[j] / sumW

        if ((weight <= 0) || numtype(weight) || numtype(Lch_nm) || numtype(Wch_nm) || numtype(Teff))
            continue
        endif
        if (Lch_nm <= 0 || Wch_nm < 0 || Teff < 0 || Teff > 1)
            continue
        endif

        if (useBetaExtraWave)
            betaExtraUse = betaExtra_List[j]
            if (numtype(betaExtraUse) != 0)
                continue
            endif
        elseif (SNS_useVortex)
            Variable th1 = atan2(rS1y_nm[j] - yV_nm, rS1x_nm[j] - xV_nm)
            Variable th2 = atan2(rS2y_nm[j] - yV_nm, rS2x_nm[j] - xV_nm)
            Variable dth = th2 - th1
            if (dth > pi)
                dth -= 2*pi
            elseif (dth < -pi)
                dth += 2*pi
            endif
            betaExtraUse = SNS_nFlux * dth
        else
            betaExtraUse = 0
        endif

        sprintf nameE2D, "%s_tmp_E_ch%03d", outPrefix, j
        sprintf nameM,   "%s_tmp_m_ch%03d", outPrefix, j
        sprintf nameS,   "%s_tmp_s_ch%03d", outPrefix, j

        nBr = Solve_AllBranches_SNS_dGSJ_betaextra( \
                    B_T, Lch, Wch, Delta, vF, lambdaL, Teff, \
                    betaExtraUse, nameE2D, nameM, nameS)

        if (nBr <= 0)
            continue
        endif

        Wave E_all = $nameE2D
        Wave m_all = $nameM
        Wave s_all = $nameS

        rec0 = nRec
        nRec += nBr
        Redimension/N=(nRec, nB) Branch_E
        Redimension/N=(nRec) Branch_ch, Branch_br, Branch_m, Branch_s
        Redimension/N=(nRec) Branch_L, Branch_W, Branch_T, Branch_weight
        Redimension/N=(nRec) Branch_x1, Branch_y1, Branch_x2, Branch_y2
        Redimension/N=(nRec) Branch_phi

        for (k = 0; k < nBr; k += 1)
            rec = rec0 + k
            Branch_ch[rec] = j
            Branch_br[rec] = k
            Branch_m[rec] = m_all[k]
            Branch_s[rec] = s_all[k]
            Branch_L[rec] = Lch_nm
            Branch_W[rec] = Wch_nm
            Branch_T[rec] = Teff
            Branch_weight[rec] = weight
            Branch_x1[rec] = rS1x_nm[j]
            Branch_y1[rec] = rS1y_nm[j]
            Branch_x2[rec] = rS2x_nm[j]
            Branch_y2[rec] = rS2y_nm[j]
            if (!ParamIsDefault(phiList_rad))
                Branch_phi[rec] = phiList_rad[j]
            else
                Branch_phi[rec] = SNS__ChordPhiFromEndpoints_rad(rS1x_nm[j], rS1y_nm[j], rS2x_nm[j], rS2y_nm[j])
            endif

            for (iB = 0; iB < nB; iB += 1)
                E0 = E_all[iB][k]
                if (numtype(E0))
                    Branch_E[rec][iB] = NaN
                else
                    Branch_E[rec][iB] = E0 * deltaScale_wv[iB]
                endif
            endfor
        endfor

        KillWaves/Z $nameE2D, $nameM, $nameS
    endfor

    if (nRec <= 0)
        Abort "SNS_ComputeDOS_FromChannels_TwoPass: no branch records generated."
    endif

    //==========================================================================
    // Pass 2a: compute branch-neighborhood rho and candidate Gamma_mix
    //==========================================================================
    String nameRho = outPrefix + "_rho_branch"
    String nameGammaMix = outPrefix + "_Gamma_mix_branch"
    Make/O/D/N=(nRec, nB) $nameRho
    Make/O/D/N=(nRec, nB) $nameGammaMix
    Wave rho_branch = $nameRho
    Wave Gamma_mix_branch = $nameGammaMix
    rho_branch = 0
    Gamma_mix_branch = 0
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", rho_branch
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", Gamma_mix_branch

    Variable rec2, dEij, dWij, kE, kW, rho
    for (iB = 0; iB < nB; iB += 1)
        for (rec = 0; rec < nRec; rec += 1)
            E0 = Branch_E[rec][iB]
            if (numtype(E0))
                rho_branch[rec][iB] = NaN
                Gamma_mix_branch[rec][iB] = NaN
                continue
            endif

            rho = 0
            for (rec2 = 0; rec2 < nRec; rec2 += 1)
                if (rec2 == rec)
                    continue
                endif
                if (numtype(Branch_E[rec2][iB]))
                    continue
                endif

                dEij = (E0 - Branch_E[rec2][iB]) / sigmaE_eV
                dWij = (Branch_W[rec] - Branch_W[rec2]) / sigmaW_nm
                kE = exp(-0.5*dEij*dEij)
                kW = exp(-0.5*dWij*dWij)
                rho += Branch_weight[rec2] * kE * kW
            endfor

            rho_branch[rec][iB] = rho
            if (GammaMixMax_eV > 0)
                Gamma_mix_branch[rec][iB] = GammaMixMax_eV * rho / (rho + rho0)
            else
                Gamma_mix_branch[rec][iB] = 0
            endif
        endfor
    endfor

    //==========================================================================
    // Pass 2b: accumulate DOS from branch table
    //==========================================================================
    Variable iE, dE, GammaTot, GammaMix, GammaUsed
    for (rec = 0; rec < nRec; rec += 1)
        Lch = Branch_L[rec] * 1e-9
        Wch = Branch_W[rec] * 1e-9
        weight = Branch_weight[rec]

        for (iB = 0; iB < nB; iB += 1)
            E0 = Branch_E[rec][iB]
            if (numtype(E0))
                continue
            endif

            Bval = B_T[iB]
            GammaTot = SNS_ComputeGammaTot(Bval, Lch, Wch, SNS_p)
            GammaMix = Gamma_mix_branch[rec][iB]
            if (useMixBroadening && GammaMix > 0)
                GammaUsed = sqrt(GammaTot*GammaTot + GammaMix*GammaMix)
            else
                GammaUsed = GammaTot
            endif

            for (iE = 0; iE < NE; iE += 1)
                dE = E_axis[iE] - E0
                DOS_EB[iE][iB] += weight * (Broadening/pi) / (dE*dE + Broadening*Broadening)
                DOS_EB_Gamma[iE][iB] += weight * (GammaTot/pi) / (dE*dE + GammaTot*GammaTot)
                DOS_EB_Mix[iE][iB] += weight * (GammaUsed/pi) / (dE*dE + GammaUsed*GammaUsed)
            endfor
        endfor
    endfor

    // ---- Normalize ABS weight and add continuum background, matching current path. ----
    Variable NE_local = DimSize(E_axis, 0)
    Variable dEgrid = E_axis[1] - E_axis[0]
    Variable Emin = E_axis[0]
    Variable Emax = E_axis[NE_local - 1]
    Variable E_sm = 3 * Broadening
    if (E_sm <= 0 || numtype(E_sm))
        E_sm = 3e-9
    endif

    if (N_cont_statesPer_eV > 0 && NE_local >= 2 && dEgrid > 0)
        Variable Egap_low  = max(Emin, -Delta)
        Variable Egap_high = min(Emax,  Delta)
        Variable EgapSpan  = Egap_high - Egap_low

        if (EgapSpan > 0)
            Variable N_gap_ref = N_cont_statesPer_eV * EgapSpan
            Make/FREE/D/N=(nB) alpha_B
            Variable acc, alpha

            for (iB = 0; iB < nB; iB += 1)
                acc = 0
                for (iE = 0; iE < NE_local; iE += 1)
                    acc += DOS_EB[iE][iB]
                endfor
                acc *= dEgrid
                if (acc > 0 && numtype(acc) == 0)
                    alpha_B[iB] = N_gap_ref / acc
                else
                    alpha_B[iB] = 0
                endif
            endfor

            for (iB = 0; iB < nB; iB += 1)
                alpha = alpha_B[iB]
                if (alpha <= 0 || numtype(alpha))
                    continue
                endif
                for (iE = 0; iE < NE_local; iE += 1)
                    DOS_EB[iE][iB] *= alpha
                    DOS_EB_Gamma[iE][iB] *= alpha
                    DOS_EB_Mix[iE][iB] *= alpha
                endfor
            endfor
        endif

        Make/FREE/D/N=(NE_local) contShape
        Variable p, E
        for (p = 0; p < NE_local; p += 1)
            E = E_axis[p]
            contShape[p] = 0.5 * (1 + tanh((abs(E) - Delta) / E_sm))
        endfor

        for (iB = 0; iB < nB; iB += 1)
            for (iE = 0; iE < NE_local; iE += 1)
                DOS_EB[iE][iB] += N_cont_statesPer_eV * contShape[iE]
                DOS_EB_Gamma[iE][iB] += N_cont_statesPer_eV * contShape[iE]
                DOS_EB_Mix[iE][iB] += N_cont_statesPer_eV * contShape[iE]
            endfor
        endfor
    endif

    String meta
    meta  = "SNS_DOSBuilder=TwoPass;"
    meta += "SNS_sigmaE_eV=" + num2str(sigmaE_eV) + ";"
    meta += "SNS_sigmaW_nm=" + num2str(sigmaW_nm) + ";"
    meta += "SNS_GammaMixMax_eV=" + num2str(GammaMixMax_eV) + ";"
    meta += "SNS_rho0=" + num2str(rho0) + ";"
    meta += "SNS_useMixBroadening=" + num2str(useMixBroadening) + ";"
    meta += "SNS_Delta_eV=" + num2str(Delta) + ";"
    meta += "SNS_N_cont_statesPer_eV=" + num2str(N_cont_statesPer_eV) + ";"
    meta += "SNS_BranchTablePrefix=" + outPrefix + ";"

    Note/K DOS_EB
    Note DOS_EB, meta
    Note/K DOS_EB_Gamma
    Note DOS_EB_Gamma, meta
    Note/K DOS_EB_Mix
    Note DOS_EB_Mix, meta
    Note/K Branch_E
    Note Branch_E, meta
    Note/K Branch_phi
    Note Branch_phi, meta
    Note/K rho_branch
    Note rho_branch, meta
    Note/K Gamma_mix_branch
    Note Gamma_mix_branch, meta

    SNS_Log("SNS_ComputeDOS_FromChannels_TwoPass: finished; branch records=" + num2str(nRec), level="INFO")
    return nRec
End

//==============================================================================
// SNS__ChordPhiFromEndpoints_rad
//
// Purpose:
//   Return the chord direction angle modulo pi [rad]. This is only a fallback
//   when a legacy branch table lacks the sampled channel angle phiList_rad.
//==============================================================================
Function SNS__ChordPhiFromEndpoints_rad(x1, y1, x2, y2)
    Variable x1, y1, x2, y2

    Variable phi = atan2(y2 - y1, x2 - x1)
    if (phi < 0)
        phi += pi
    endif
    if (phi >= pi)
        phi -= pi
    endif

    return phi
End

//==============================================================================
// SNS__EndpointDistance2D_nm
//
// Purpose:
//   Endpoint-pair distance between two 2D SNS chords, independent of endpoint
//   ordering. Used as a geometry kernel for local branch hybridization.
//==============================================================================
Function SNS__EndpointDistance2D_nm(x1a, y1a, x2a, y2a, x1b, y1b, x2b, y2b)
    Variable x1a, y1a, x2a, y2a
    Variable x1b, y1b, x2b, y2b

    Variable dSame = 0.5*((x1a-x1b)^2 + (y1a-y1b)^2 + (x2a-x2b)^2 + (y2a-y2b)^2)
    Variable dSwap = 0.5*((x1a-x2b)^2 + (y1a-y2b)^2 + (x2a-x1b)^2 + (y2a-y1b)^2)

    return sqrt(min(dSame, dSwap))
End

//==============================================================================
// SNS__HybridCoupling_ij
//
// Purpose:
//   Phenomenological channel-coupling matrix element V_ij [eV].
//
// Notes:
//   V0_eV is the fit/control parameter. The remaining kernels only suppress
//   coupling between geometrically unlike branches.
//==============================================================================
Function SNS__HybridCoupling_ij(recA, recB, Branch_L, Branch_W, Branch_T, Branch_x1, Branch_y1, Branch_x2, Branch_y2, Branch_phi, V0_eV, sigmaW_nm, sigmaL_nm, sigmaEnd_nm, sigmaPhi_rad)
    Variable recA, recB
    Wave Branch_L, Branch_W, Branch_T
    Wave Branch_x1, Branch_y1, Branch_x2, Branch_y2
    Wave Branch_phi
    Variable V0_eV, sigmaW_nm, sigmaL_nm, sigmaEnd_nm, sigmaPhi_rad

    Variable arg = 0
    Variable dW, dL, dEnd, dPhi

    if (sigmaW_nm > 0)
        dW = (Branch_W[recA] - Branch_W[recB]) / sigmaW_nm
        arg += dW*dW
    endif

    if (sigmaL_nm > 0)
        dL = (Branch_L[recA] - Branch_L[recB]) / sigmaL_nm
        arg += dL*dL
    endif

    if (sigmaEnd_nm > 0)
        dEnd = SNS__EndpointDistance2D_nm( \
            Branch_x1[recA], Branch_y1[recA], Branch_x2[recA], Branch_y2[recA], \
            Branch_x1[recB], Branch_y1[recB], Branch_x2[recB], Branch_y2[recB]) / sigmaEnd_nm
        arg += dEnd*dEnd
    endif

    if (sigmaPhi_rad > 0)
        dPhi = abs(Branch_phi[recA] - Branch_phi[recB])
        dPhi = min(dPhi, pi - dPhi)
        dPhi /= sigmaPhi_rad
        arg += dPhi*dPhi
    endif

    Variable Tfac = sqrt(max(0, Branch_T[recA]) * max(0, Branch_T[recB]))
    return abs(V0_eV) * exp(-0.5*arg) * Tfac
End

//==============================================================================
// SNS__JacobiEigenSymmetric
//
// Purpose:
//   Diagonalize a small real symmetric matrix A using Jacobi rotations.
//
// Inputs/Outputs:
//   A      : on input, symmetric matrix; overwritten during diagonalization.
//   eval   : eigenvalues [n]
//   evec   : eigenvectors, column-major evec[row][eigenIndex]
//   n      : active matrix size
//
// Returns:
//   number of sweeps/rotations attempted.
//==============================================================================
Function SNS__JacobiEigenSymmetric(A, eval, evec, n)
    Wave A, eval, evec
    Variable n

    Variable i, j, k
    evec = 0
    for (i = 0; i < n; i += 1)
        evec[i][i] = 1
    endfor

    Variable maxIter = max(50, 50*n*n)
    Variable tol = 1e-18
    Variable iter, p, q, maxOff, val
    Variable app, aqq, apq, tau, t, c, s
    Variable akp, akq, vip, viq

    for (iter = 0; iter < maxIter; iter += 1)
        maxOff = 0
        p = 0
        q = 1

        for (i = 0; i < n-1; i += 1)
            for (j = i+1; j < n; j += 1)
                val = abs(A[i][j])
                if (val > maxOff)
                    maxOff = val
                    p = i
                    q = j
                endif
            endfor
        endfor

        if (maxOff < tol)
            break
        endif

        app = A[p][p]
        aqq = A[q][q]
        apq = A[p][q]
        if (apq == 0)
            continue
        endif

        tau = (aqq - app) / (2*apq)
        if (tau >= 0)
            t = 1 / (tau + sqrt(1 + tau*tau))
        else
            t = -1 / (-tau + sqrt(1 + tau*tau))
        endif
        c = 1 / sqrt(1 + t*t)
        s = t*c

        A[p][p] = app - t*apq
        A[q][q] = aqq + t*apq
        A[p][q] = 0
        A[q][p] = 0

        for (k = 0; k < n; k += 1)
            if ((k == p) || (k == q))
                continue
            endif

            akp = A[k][p]
            akq = A[k][q]
            A[k][p] = c*akp - s*akq
            A[p][k] = A[k][p]
            A[k][q] = s*akp + c*akq
            A[q][k] = A[k][q]
        endfor

        for (k = 0; k < n; k += 1)
            vip = evec[k][p]
            viq = evec[k][q]
            evec[k][p] = c*vip - s*viq
            evec[k][q] = s*vip + c*viq
        endfor
    endfor

    for (i = 0; i < n; i += 1)
        eval[i] = A[i][i]
    endfor

    return iter
End

//==============================================================================
// SNS_ComputeDOS_FromBranchTable_Hybridized
//
// Purpose:
//   Build a DOS(E,B) by local hybridization of an existing two-pass branch table.
//
//   At each field, mutually nearest nearby branches are paired and represented
//   by a 2x2 real symmetric Hamiltonian
//
//       H_ij = E_i(B) delta_ij + V_ij
//
//   with V_ij set by V0_eV and geometry kernels. The paired eigenvalues are
//   deposited with the ordinary narrow DOS broadening. This models level
//   repulsion / avoided crossings rather than large Lorentzian tails.
//
// Inputs:
//   branchPrefix : prefix used by SNS_ComputeDOS_FromChannels_TwoPass, e.g.
//                  "TwoPassDiag".
//   nameDOS      : output hybridized DOS wave.
//   nameEaxis    : output energy-axis wave.
//
// Optional Inputs:
//   V0_eV              : channel-coupling energy scale [eV]. Default 0.
//   clusterEnergy_eV   : seed energy window for cluster candidates [eV].
//                        Default max(6*V0_eV, 4*Broadening).
//   sigmaW_nm          : W_eff coupling kernel width [nm]. Default 25.
//   sigmaL_nm          : L_N coupling kernel width [nm]. Default 0 (disabled).
//   sigmaEnd_nm        : endpoint-pair coupling kernel width [nm]. Default 0.
//   maxCluster         : retained for call compatibility; ignored by this
//                        pairwise implementation.
//   GammaHybrid_eV     : Lorentzian HWHM used to deposit eigenvalues.
//                        Default root:SNS_Settings:Broadening.
//   N_cont_statesPer_eV: optional continuum background scale. Default 0.
//   Delta_fit_eV       : optional Delta for continuum shape. Default settings.
//   outPrefix          : diagnostic prefix. Default nameDOS + "_Hybrid".
//   slopeMin_eVperT    : minimum relative slope |dEi/dB-dEj/dB| for an
//                        avoided-crossing candidate. Default 0 disables.
//   crossingWindow_T   : require the linearized crossing field to be within
//                        this field distance of the current B. Default:
//                        disabled when <= 0.
//   sigmaPhi_rad       : channel-angle coupling kernel width [rad].
//                        Default 0 disables angular filtering.
//   sourceEaxis        : optional DOS energy axis to reuse. If omitted, this
//                        function uses E_axis_2pass when present, then falls
//                        back to the standard +/-2.5 meV axis.
//   minVijActive_eV    : optional diagnostic threshold for active coupling.
//                        Events with peak or instantaneous V_ij below this
//                        value are ignored. Default 0 disables this threshold.
//   eventPruneFactor   : controls global pruning of overlapping crossing
//                        events that share a branch. Default 1 reproduces the
//                        previous behavior. Set <= 0 to keep all events and
//                        let the local mutual-best selection choose the active
//                        protected pair at each field.
//
// Outputs:
//   nameDOS
//   nameEaxis
//   <outPrefix>_clusterSize
//   <outPrefix>_hybridShift_eV
//   <outPrefix>_partner
//   <outPrefix>_Vij_eV
//   <outPrefix>_Bcross_T
//   <outPrefix>_eventIndex
//
// Returns:
//   number of clusters diagonalized.
//==============================================================================
Function SNS_ComputeDOS_FromBranchTable_Hybridized(branchPrefix, nameDOS, nameEaxis, [V0_eV, clusterEnergy_eV, sigmaW_nm, sigmaL_nm, sigmaEnd_nm, maxCluster, GammaHybrid_eV, N_cont_statesPer_eV, Delta_fit_eV, outPrefix, slopeMin_eVperT, crossingWindow_T, sigmaPhi_rad, sourceEaxis, minVijActive_eV, eventPruneFactor])
    String branchPrefix, nameDOS, nameEaxis
    Variable V0_eV, clusterEnergy_eV
    Variable sigmaW_nm, sigmaL_nm, sigmaEnd_nm, sigmaPhi_rad
    Variable maxCluster, GammaHybrid_eV, N_cont_statesPer_eV, Delta_fit_eV
    String outPrefix
    Variable slopeMin_eVperT, crossingWindow_T, minVijActive_eV, eventPruneFactor
    Wave sourceEaxis

    STRUCT SNS_Params SNS_p
    SNS_LoadParams(SNS_p)

    DFREF dfrSet = root:SNS_Settings
    Wave B_T = dfrSet:B_T

    Wave Branch_E = $(branchPrefix + "_Branch_E_eV")
    Wave Branch_L = $(branchPrefix + "_Branch_L_nm")
    Wave Branch_W = $(branchPrefix + "_Branch_W_nm")
    Wave Branch_T = $(branchPrefix + "_Branch_T")
    Wave Branch_weight = $(branchPrefix + "_Branch_weight")
    Wave Branch_x1 = $(branchPrefix + "_Branch_x1_nm")
    Wave Branch_y1 = $(branchPrefix + "_Branch_y1_nm")
    Wave Branch_x2 = $(branchPrefix + "_Branch_x2_nm")
    Wave Branch_y2 = $(branchPrefix + "_Branch_y2_nm")
    Wave/Z Branch_phi_existing = $(branchPrefix + "_Branch_phi_rad")
    Wave/Z defaultSourceEaxis = E_axis_2pass
    String branchMeta = note(Branch_E)

    Variable nRec = DimSize(Branch_E, 0)
    Variable nB = DimSize(Branch_E, 1)
    Variable nBset = numpnts(B_T)
    if (nRec <= 0 || nB <= 0 || nB != nBset)
        Abort "SNS_ComputeDOS_FromBranchTable_Hybridized: invalid branch table or B_T mismatch."
    endif

    Variable NE = SNS_p.NE
    Variable Delta = SNS_p.Delta
    if (!ParamIsDefault(Delta_fit_eV))
        Delta = Delta_fit_eV
    else
        String deltaMetaStr = StringByKey("SNS_Delta_eV", branchMeta, "=", ";")
        if (strlen(deltaMetaStr) > 0)
            Variable deltaMeta = str2num(deltaMetaStr)
            if (numtype(deltaMeta) == 0 && deltaMeta > 0)
                Delta = deltaMeta
            endif
        endif
    endif

    if (ParamIsDefault(V0_eV) || V0_eV < 0 || numtype(V0_eV))
        V0_eV = 0
    endif
    if (ParamIsDefault(GammaHybrid_eV) || GammaHybrid_eV <= 0 || numtype(GammaHybrid_eV))
        GammaHybrid_eV = SNS_p.Broadening
    endif
    if (GammaHybrid_eV <= 0 || numtype(GammaHybrid_eV))
        GammaHybrid_eV = 1e-9
    endif
    if (ParamIsDefault(clusterEnergy_eV) || clusterEnergy_eV <= 0 || numtype(clusterEnergy_eV))
        clusterEnergy_eV = max(6*V0_eV, 4*GammaHybrid_eV)
    endif
    if (clusterEnergy_eV <= 0 || numtype(clusterEnergy_eV))
        clusterEnergy_eV = 4e-6
    endif
    if (ParamIsDefault(sigmaW_nm) || sigmaW_nm < 0 || numtype(sigmaW_nm))
        sigmaW_nm = 25
    endif
    if (ParamIsDefault(sigmaL_nm) || sigmaL_nm < 0 || numtype(sigmaL_nm))
        sigmaL_nm = 0
    endif
    if (ParamIsDefault(sigmaEnd_nm) || sigmaEnd_nm < 0 || numtype(sigmaEnd_nm))
        sigmaEnd_nm = 0
    endif
    if (ParamIsDefault(sigmaPhi_rad) || sigmaPhi_rad < 0 || numtype(sigmaPhi_rad))
        sigmaPhi_rad = 0
    endif
    if (ParamIsDefault(maxCluster) || maxCluster < 1 || numtype(maxCluster))
        maxCluster = 12
    endif
    maxCluster = round(maxCluster)
    if (ParamIsDefault(N_cont_statesPer_eV))
        String nContMetaStr = StringByKey("SNS_N_cont_statesPer_eV", branchMeta, "=", ";")
        if (strlen(nContMetaStr) > 0)
            N_cont_statesPer_eV = str2num(nContMetaStr)
        else
            N_cont_statesPer_eV = 0
        endif
    endif
    if (N_cont_statesPer_eV < 0 || numtype(N_cont_statesPer_eV))
        N_cont_statesPer_eV = 0
    endif
    if (ParamIsDefault(outPrefix) || strlen(outPrefix) <= 0)
        outPrefix = nameDOS + "_Hybrid"
    endif
    if (ParamIsDefault(slopeMin_eVperT) || slopeMin_eVperT < 0 || numtype(slopeMin_eVperT))
        slopeMin_eVperT = 0
    endif
    if (ParamIsDefault(crossingWindow_T) || crossingWindow_T < 0 || numtype(crossingWindow_T))
        crossingWindow_T = 0
    endif
    if (ParamIsDefault(minVijActive_eV) || minVijActive_eV < 0 || numtype(minVijActive_eV))
        minVijActive_eV = 0
    endif
    if (ParamIsDefault(eventPruneFactor) || numtype(eventPruneFactor))
        eventPruneFactor = 1
    endif

    Make/O/D/N=(NE) $nameEaxis
    Wave E_axis = $nameEaxis
    if (!ParamIsDefault(sourceEaxis) && WaveExists(sourceEaxis) && (numpnts(sourceEaxis) == NE))
        SetScale/P x, DimOffset(sourceEaxis,0), DimDelta(sourceEaxis,0), WaveUnits(sourceEaxis,0), E_axis
    elseif (WaveExists(defaultSourceEaxis) && (numpnts(defaultSourceEaxis) == NE))
        SetScale/P x, DimOffset(defaultSourceEaxis,0), DimDelta(defaultSourceEaxis,0), WaveUnits(defaultSourceEaxis,0), E_axis
    else
        SetScale/I x, -2.5e-3, 2.5e-3, "eV", E_axis
    endif
    E_axis = x

    Make/O/D/N=(NE, nB) $nameDOS
    Wave DOS_Hybrid = $nameDOS
    DOS_Hybrid = 0
    SetScale/P x, DimOffset(E_axis,0), DimDelta(E_axis,0), "eV", DOS_Hybrid
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", DOS_Hybrid

    String nameCl = outPrefix + "_clusterSize"
    String nameShift = outPrefix + "_hybridShift_eV"
    String namePartner = outPrefix + "_partner"
    String nameVij = outPrefix + "_Vij_eV"
    String nameBcross = outPrefix + "_Bcross_T"
    String nameEventIndex = outPrefix + "_eventIndex"
    Make/O/D/N=(nRec, nB) $nameCl
    Make/O/D/N=(nRec, nB) $nameShift
    Make/O/D/N=(nRec, nB) $namePartner
    Make/O/D/N=(nRec, nB) $nameVij
    Make/O/D/N=(nRec, nB) $nameBcross
    Make/O/D/N=(nRec, nB) $nameEventIndex
    Wave clusterSize = $nameCl
    Wave hybridShift = $nameShift
    Wave hybridPartner = $namePartner
    Wave hybridVij = $nameVij
    Wave hybridBcross = $nameBcross
    Wave hybridEventIndex = $nameEventIndex
    clusterSize = NaN
    hybridShift = NaN
    hybridPartner = NaN
    hybridVij = NaN
    hybridBcross = NaN
    hybridEventIndex = NaN
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", clusterSize
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", hybridShift
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", hybridPartner
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", hybridVij
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", hybridBcross
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", hybridEventIndex

    Make/FREE/D/N=(nRec) Branch_phi
    if (WaveExists(Branch_phi_existing) && (numpnts(Branch_phi_existing) == nRec))
        Branch_phi = Branch_phi_existing[p]
    else
        Branch_phi = SNS__ChordPhiFromEndpoints_rad(Branch_x1[p], Branch_y1[p], Branch_x2[p], Branch_y2[p])
    endif

    Make/FREE/B/N=(nRec) assigned
    Make/FREE/D/N=(nRec) bestPartner, bestScore, bestVij, bestBcross
    Make/FREE/D/N=(nRec) bestEvent
    Make/FREE/D/N=0 eventA, eventB, eventBcross, eventVij, eventSigmaB

    Variable iB, rec, rec2, iE
    Variable Ei, Ej, score, Vij, dE
    Variable nClusters = 0
    Variable alpha
    Variable partner
    Variable Eavg, det, split, Eplus, Eminus, cos2, sin2, wPlus, wMinus
    Variable kBseg, Eia, Eib, Eja, Ejb, D0, D1, B0seg, B1seg, dBseg
    Variable dSlope, Bcross, sigmaB, Bnow, Veff, nEvents = 0
    Variable ev1, ev2, shareBranch, eventSep_T

    // Build field-continuous crossing events from bare branch crossings.
    for (rec = 0; rec < nRec-1; rec += 1)
        for (rec2 = rec + 1; rec2 < nRec; rec2 += 1)

            Vij = SNS__HybridCoupling_ij(rec, rec2, Branch_L, Branch_W, Branch_T, Branch_x1, Branch_y1, Branch_x2, Branch_y2, Branch_phi, V0_eV, sigmaW_nm, sigmaL_nm, sigmaEnd_nm, sigmaPhi_rad)
            if (Vij <= 0 || numtype(Vij))
                continue
            endif
            if ((minVijActive_eV > 0) && (Vij < minVijActive_eV))
                continue
            endif

            for (kBseg = 0; kBseg < nB-1; kBseg += 1)
                Eia = Branch_E[rec][kBseg]
                Eib = Branch_E[rec][kBseg+1]
                Eja = Branch_E[rec2][kBseg]
                Ejb = Branch_E[rec2][kBseg+1]
                if (numtype(Eia) || numtype(Eib) || numtype(Eja) || numtype(Ejb))
                    continue
                endif

                D0 = Eia - Eja
                D1 = Eib - Ejb
                if ((abs(D0) > clusterEnergy_eV) && (abs(D1) > clusterEnergy_eV))
                    continue
                endif
                if (D0*D1 > 0)
                    continue
                endif

                B0seg = B_T[kBseg]
                B1seg = B_T[kBseg+1]
                dBseg = B1seg - B0seg
                if (dBseg == 0)
                    continue
                endif

                dSlope = (D1 - D0) / dBseg
                if (abs(dSlope) < slopeMin_eVperT)
                    continue
                endif
                if (dSlope == 0)
                    continue
                endif

                Bcross = B0seg - D0/dSlope
                if ((Bcross < min(B0seg, B1seg)) || (Bcross > max(B0seg, B1seg)))
                    continue
                endif

                if (crossingWindow_T > 0)
                    sigmaB = crossingWindow_T
                else
                    sigmaB = max(abs(dBseg), 3*abs(Vij/dSlope))
                endif
                if (sigmaB <= 0 || numtype(sigmaB))
                    continue
                endif

                Redimension/N=(nEvents + 1) eventA, eventB, eventBcross, eventVij, eventSigmaB
                eventA[nEvents] = rec
                eventB[nEvents] = rec2
                eventBcross[nEvents] = Bcross
                eventVij[nEvents] = Vij
                eventSigmaB[nEvents] = sigmaB
                nEvents += 1
            endfor
        endfor
    endfor

    Make/FREE/B/N=(nEvents) eventKeep
    eventKeep = 1

    if (eventPruneFactor > 0)
        // Remove overlapping crossing events that compete for the same branch.
        // Setting eventPruneFactor<=0 keeps all events; the field-local
        // mutual-best selection below then lets a branch participate in
        // multiple protected pair encounters along the sweep.
        for (ev1 = 0; ev1 < nEvents-1; ev1 += 1)
            if (!eventKeep[ev1])
                continue
            endif
            for (ev2 = ev1 + 1; ev2 < nEvents; ev2 += 1)
                if (!eventKeep[ev2])
                    continue
                endif

                shareBranch = (eventA[ev1] == eventA[ev2]) || (eventA[ev1] == eventB[ev2]) || \
                              (eventB[ev1] == eventA[ev2]) || (eventB[ev1] == eventB[ev2])
                if (!shareBranch)
                    continue
                endif

                eventSep_T = eventPruneFactor * max(eventSigmaB[ev1], eventSigmaB[ev2])
                if (abs(eventBcross[ev1] - eventBcross[ev2]) > eventSep_T)
                    continue
                endif

                if (eventVij[ev1] >= eventVij[ev2])
                    eventKeep[ev2] = 0
                else
                    eventKeep[ev1] = 0
                    break
                endif
            endfor
        endfor
    endif

    for (iB = 0; iB < nB; iB += 1)
        assigned = 0
        bestPartner = -1
        bestScore = Inf
        bestVij = 0
        bestBcross = NaN
        bestEvent = NaN
        Bnow = B_T[iB]

        // Activate stored crossing events smoothly around their crossing field.
        for (rec2 = 0; rec2 < nEvents; rec2 += 1)
            if (!eventKeep[rec2])
                continue
            endif

            rec = eventA[rec2]
            partner = eventB[rec2]

            Ei = Branch_E[rec][iB]
            Ej = Branch_E[partner][iB]
            if (numtype(Ei) || numtype(Ej))
                continue
            endif

            dBseg = (Bnow - eventBcross[rec2]) / eventSigmaB[rec2]
            if (abs(dBseg) > 4)
                continue
            endif

            Veff = eventVij[rec2] * exp(-0.5*dBseg*dBseg)
            if (Veff <= 0 || numtype(Veff))
                continue
            endif
            if ((minVijActive_eV > 0) && (Veff < minVijActive_eV))
                continue
            endif

            score = abs(Ei - Ej) / Veff + 0.25*abs(dBseg)

            if (score < bestScore[rec])
                bestScore[rec] = score
                bestPartner[rec] = partner
                bestVij[rec] = Veff
                bestBcross[rec] = eventBcross[rec2]
                bestEvent[rec] = rec2
            endif
            if (score < bestScore[partner])
                bestScore[partner] = score
                bestPartner[partner] = rec
                bestVij[partner] = Veff
                bestBcross[partner] = eventBcross[rec2]
                bestEvent[partner] = rec2
            endif
        endfor

        for (rec = 0; rec < nRec; rec += 1)
            if (assigned[rec])
                continue
            endif
            Ei = Branch_E[rec][iB]
            if (numtype(Ei))
                assigned[rec] = 1
                continue
            endif

            partner = bestPartner[rec]
            if ((partner >= 0) && (partner < nRec) && !assigned[partner] && (bestPartner[partner] == rec))
                Ej = Branch_E[partner][iB]
                Vij = bestVij[rec]

                if (!numtype(Ej) && Vij > 0)
                    Eavg = 0.5*(Ei + Ej)
                    det = 0.5*(Ei - Ej)
                    split = sqrt(det*det + Vij*Vij)
                    Eminus = Eavg - split
                    Eplus = Eavg + split

                    // Eigenvector weights for the lower/upper eigenvalues.
                    if (split > 0)
                        cos2 = 0.5*(1 + det/split)
                        sin2 = 0.5*(1 - det/split)
                    else
                        cos2 = 0.5
                        sin2 = 0.5
                    endif

                    wMinus = Branch_weight[rec]*sin2 + Branch_weight[partner]*cos2
                    wPlus = Branch_weight[rec]*cos2 + Branch_weight[partner]*sin2

                    for (iE = 0; iE < NE; iE += 1)
                        dE = E_axis[iE] - Eminus
                        DOS_Hybrid[iE][iB] += wMinus * (GammaHybrid_eV/pi) / (dE*dE + GammaHybrid_eV*GammaHybrid_eV)
                        dE = E_axis[iE] - Eplus
                        DOS_Hybrid[iE][iB] += wPlus * (GammaHybrid_eV/pi) / (dE*dE + GammaHybrid_eV*GammaHybrid_eV)
                    endfor

                    clusterSize[rec][iB] = 2
                    clusterSize[partner][iB] = 2
                    hybridShift[rec][iB] = min(abs(Eminus - Ei), abs(Eplus - Ei))
                    hybridShift[partner][iB] = min(abs(Eminus - Ej), abs(Eplus - Ej))
                    hybridPartner[rec][iB] = partner
                    hybridPartner[partner][iB] = rec
                    hybridVij[rec][iB] = Vij
                    hybridVij[partner][iB] = Vij
                    hybridBcross[rec][iB] = bestBcross[rec]
                    hybridBcross[partner][iB] = bestBcross[rec]
                    hybridEventIndex[rec][iB] = bestEvent[rec]
                    hybridEventIndex[partner][iB] = bestEvent[rec]

                    assigned[rec] = 1
                    assigned[partner] = 1
                    nClusters += 1
                    continue
                endif
            endif

            // Unpaired branch: deposit unchanged.
            for (iE = 0; iE < NE; iE += 1)
                dE = E_axis[iE] - Ei
                DOS_Hybrid[iE][iB] += Branch_weight[rec] * (GammaHybrid_eV/pi) / (dE*dE + GammaHybrid_eV*GammaHybrid_eV)
            endfor
            clusterSize[rec][iB] = 1
            hybridShift[rec][iB] = 0
            assigned[rec] = 1
        endfor
    endfor

    // Match the ABS normalization convention used by the branch-table builder.
    Variable dEgrid = E_axis[1] - E_axis[0]
    Variable Emin = E_axis[0]
    Variable Emax = E_axis[NE - 1]
    Variable Egap_low = max(Emin, -Delta)
    Variable Egap_high = min(Emax, Delta)
    Variable EgapSpan = Egap_high - Egap_low
    Variable acc, E, contShapeVal

    if (N_cont_statesPer_eV > 0 && dEgrid > 0 && EgapSpan > 0)
        Variable N_gap_ref = N_cont_statesPer_eV * EgapSpan
        for (iB = 0; iB < nB; iB += 1)
            acc = 0
            for (iE = 0; iE < NE; iE += 1)
                acc += DOS_Hybrid[iE][iB]
            endfor
            acc *= dEgrid
            if (acc > 0 && numtype(acc) == 0)
                alpha = N_gap_ref / acc
                for (iE = 0; iE < NE; iE += 1)
                    DOS_Hybrid[iE][iB] *= alpha
                endfor
            endif
        endfor

        Variable E_sm = 3*GammaHybrid_eV
        if (E_sm <= 0 || numtype(E_sm))
            E_sm = 3e-9
        endif
        for (iE = 0; iE < NE; iE += 1)
            E = E_axis[iE]
            contShapeVal = 0.5 * (1 + tanh((abs(E) - Delta) / E_sm))
            for (iB = 0; iB < nB; iB += 1)
                DOS_Hybrid[iE][iB] += N_cont_statesPer_eV * contShapeVal
            endfor
        endfor
    endif

    String meta
    meta  = "SNS_DOSBuilder=BranchTableHybridized;"
    meta += "SNS_BranchTablePrefix=" + branchPrefix + ";"
    meta += "SNS_V0_eV=" + num2str(V0_eV) + ";"
    meta += "SNS_clusterEnergy_eV=" + num2str(clusterEnergy_eV) + ";"
    meta += "SNS_sigmaW_nm=" + num2str(sigmaW_nm) + ";"
    meta += "SNS_sigmaL_nm=" + num2str(sigmaL_nm) + ";"
    meta += "SNS_sigmaEnd_nm=" + num2str(sigmaEnd_nm) + ";"
    meta += "SNS_sigmaPhi_rad=" + num2str(sigmaPhi_rad) + ";"
    meta += "SNS_maxCluster=" + num2str(maxCluster) + ";"
    meta += "SNS_GammaHybrid_eV=" + num2str(GammaHybrid_eV) + ";"
    meta += "SNS_Delta_eV=" + num2str(Delta) + ";"
    meta += "SNS_N_cont_statesPer_eV=" + num2str(N_cont_statesPer_eV) + ";"
    meta += "SNS_slopeMin_eVperT=" + num2str(slopeMin_eVperT) + ";"
    meta += "SNS_crossingWindow_T=" + num2str(crossingWindow_T) + ";"
    meta += "SNS_minVijActive_eV=" + num2str(minVijActive_eV) + ";"
    meta += "SNS_eventPruneFactor=" + num2str(eventPruneFactor) + ";"

    Note/K DOS_Hybrid
    Note DOS_Hybrid, meta
    Note/K clusterSize
    Note clusterSize, meta
    Note/K hybridShift
    Note hybridShift, meta
    Note/K hybridPartner
    Note hybridPartner, meta
    Note/K hybridVij
    Note hybridVij, meta
    Note/K hybridBcross
    Note hybridBcross, meta
    Note/K hybridEventIndex
    Note hybridEventIndex, meta

    SNS_Log("SNS_ComputeDOS_FromBranchTable_Hybridized: finished; clusters=" + num2str(nClusters), level="INFO")
    return nClusters
End

//==============================================================================
// SNS_ComputeDOS_FromBranchTable_ManifoldHybridized
//
// Purpose:
//   Compute DOS(E,B) from a two-pass branch table using local coherent
//   hybridization inside nearly-parallel phase-space manifolds.
//
//   Unlike SNS_ComputeDOS_FromBranchTable_Hybridized(...), this model does not
//   require a branch crossing. At each field point it builds small local
//   Hamiltonian blocks from branches that are close in energy, geometry, angle,
//   and field slope dE/dB. The blocks are diagonalized and their eigenvectors
//   redistribute branch spectral weight among the local eigenvalues.
//
// Inputs:
//   branchPrefix : prefix used by SNS_ComputeDOS_FromChannels_TwoPass, e.g.
//                  "TwoPassDiag".
//   nameDOS      : output DOS wave name.
//   nameEaxis    : output energy-axis wave name.
//
// Optional Inputs:
//   V0_eV              : maximum off-diagonal manifold coupling [eV].
//                        Default: 0.
//   clusterEnergy_eV   : local energy window / kernel scale [eV].
//                        Default: max(6*V0_eV, 4*GammaHybrid_eV).
//   sigmaSlope_eVperT  : slope kernel width for dE/dB similarity [eV/T].
//                        Default: 0, which disables slope filtering.
//   sigmaW_nm          : W_eff kernel width [nm]. Default: 25.
//   sigmaL_nm          : L_N kernel width [nm]. Default: 0, disabled.
//   sigmaEnd_nm        : endpoint-pair kernel width [nm]. Default: 0, disabled.
//   sigmaPhi_rad       : chord-angle kernel width [rad]. Default: 0, disabled.
//   maxCluster         : maximum local Hamiltonian size. Default: 6.
//   GammaHybrid_eV     : Lorentzian deposition width. Default: Broadening.
//   N_cont_statesPer_eV: continuum DOS scale. If omitted, inherited from the
//                        branch table note when available.
//   Delta_fit_eV       : gap used for continuum placement. If omitted,
//                        inherited from branch table note when available.
//   outPrefix          : diagnostic output prefix. Default: nameDOS+"_Manifold".
//   sourceEaxis        : optional energy axis to reuse.
//   minVijActive_eV    : optional diagnostic threshold for active coupling.
//                        Default: 0, disabled.
//
// Outputs:
//   nameDOS, nameEaxis.
//   Diagnostic waves:
//       <outPrefix>_clusterSize
//       <outPrefix>_hybridShift_eV
//       <outPrefix>_partner
//       <outPrefix>_Vij_eV
//       <outPrefix>_slope_eVperT
//
// Returns:
//   number of local manifold blocks diagonalized.
//==============================================================================
Function SNS_ComputeDOS_FromBranchTable_ManifoldHybridized(branchPrefix, nameDOS, nameEaxis, [V0_eV, clusterEnergy_eV, sigmaSlope_eVperT, sigmaW_nm, sigmaL_nm, sigmaEnd_nm, sigmaPhi_rad, maxCluster, GammaHybrid_eV, N_cont_statesPer_eV, Delta_fit_eV, outPrefix, sourceEaxis, minVijActive_eV])
    String branchPrefix, nameDOS, nameEaxis
    Variable V0_eV, clusterEnergy_eV, sigmaSlope_eVperT
    Variable sigmaW_nm, sigmaL_nm, sigmaEnd_nm, sigmaPhi_rad
    Variable maxCluster, GammaHybrid_eV, N_cont_statesPer_eV, Delta_fit_eV
    String outPrefix
    Wave sourceEaxis
    Variable minVijActive_eV

    STRUCT SNS_Params SNS_p
    SNS_LoadParams(SNS_p)

    DFREF dfrSet = root:SNS_Settings
    Wave B_T = dfrSet:B_T

    Wave Branch_E = $(branchPrefix + "_Branch_E_eV")
    Wave Branch_L = $(branchPrefix + "_Branch_L_nm")
    Wave Branch_W = $(branchPrefix + "_Branch_W_nm")
    Wave Branch_T = $(branchPrefix + "_Branch_T")
    Wave Branch_weight = $(branchPrefix + "_Branch_weight")
    Wave Branch_x1 = $(branchPrefix + "_Branch_x1_nm")
    Wave Branch_y1 = $(branchPrefix + "_Branch_y1_nm")
    Wave Branch_x2 = $(branchPrefix + "_Branch_x2_nm")
    Wave Branch_y2 = $(branchPrefix + "_Branch_y2_nm")
    Wave/Z Branch_phi_existing = $(branchPrefix + "_Branch_phi_rad")
    Wave/Z defaultSourceEaxis = E_axis_2pass
    String branchMeta = note(Branch_E)

    Variable nRec = DimSize(Branch_E, 0)
    Variable nB = DimSize(Branch_E, 1)
    Variable nBset = numpnts(B_T)
    if (nRec <= 0 || nB <= 0 || nB != nBset)
        Abort "SNS_ComputeDOS_FromBranchTable_ManifoldHybridized: invalid branch table or B_T mismatch."
    endif

    Variable NE = SNS_p.NE
    Variable Delta = SNS_p.Delta
    if (!ParamIsDefault(Delta_fit_eV))
        Delta = Delta_fit_eV
    else
        String deltaMetaStr = StringByKey("SNS_Delta_eV", branchMeta, "=", ";")
        if (strlen(deltaMetaStr) > 0)
            Variable deltaMeta = str2num(deltaMetaStr)
            if (numtype(deltaMeta) == 0 && deltaMeta > 0)
                Delta = deltaMeta
            endif
        endif
    endif

    if (ParamIsDefault(V0_eV) || V0_eV < 0 || numtype(V0_eV))
        V0_eV = 0
    endif
    if (ParamIsDefault(GammaHybrid_eV) || GammaHybrid_eV <= 0 || numtype(GammaHybrid_eV))
        GammaHybrid_eV = SNS_p.Broadening
    endif
    if (GammaHybrid_eV <= 0 || numtype(GammaHybrid_eV))
        GammaHybrid_eV = 1e-9
    endif
    if (ParamIsDefault(clusterEnergy_eV) || clusterEnergy_eV <= 0 || numtype(clusterEnergy_eV))
        clusterEnergy_eV = max(6*V0_eV, 4*GammaHybrid_eV)
    endif
    if (clusterEnergy_eV <= 0 || numtype(clusterEnergy_eV))
        clusterEnergy_eV = 4e-6
    endif
    if (ParamIsDefault(sigmaSlope_eVperT) || sigmaSlope_eVperT < 0 || numtype(sigmaSlope_eVperT))
        sigmaSlope_eVperT = 0
    endif
    if (ParamIsDefault(sigmaW_nm) || sigmaW_nm < 0 || numtype(sigmaW_nm))
        sigmaW_nm = 25
    endif
    if (ParamIsDefault(sigmaL_nm) || sigmaL_nm < 0 || numtype(sigmaL_nm))
        sigmaL_nm = 0
    endif
    if (ParamIsDefault(sigmaEnd_nm) || sigmaEnd_nm < 0 || numtype(sigmaEnd_nm))
        sigmaEnd_nm = 0
    endif
    if (ParamIsDefault(sigmaPhi_rad) || sigmaPhi_rad < 0 || numtype(sigmaPhi_rad))
        sigmaPhi_rad = 0
    endif
    if (ParamIsDefault(maxCluster) || maxCluster < 2 || numtype(maxCluster))
        maxCluster = 6
    endif
    maxCluster = round(maxCluster)
    if (maxCluster < 2)
        maxCluster = 2
    endif
    if (ParamIsDefault(N_cont_statesPer_eV))
        String nContMetaStr = StringByKey("SNS_N_cont_statesPer_eV", branchMeta, "=", ";")
        if (strlen(nContMetaStr) > 0)
            N_cont_statesPer_eV = str2num(nContMetaStr)
        else
            N_cont_statesPer_eV = 0
        endif
    endif
    if (N_cont_statesPer_eV < 0 || numtype(N_cont_statesPer_eV))
        N_cont_statesPer_eV = 0
    endif
    if (ParamIsDefault(outPrefix) || strlen(outPrefix) <= 0)
        outPrefix = nameDOS + "_Manifold"
    endif
    if (ParamIsDefault(minVijActive_eV) || minVijActive_eV < 0 || numtype(minVijActive_eV))
        minVijActive_eV = 0
    endif

    Make/O/D/N=(NE) $nameEaxis
    Wave E_axis = $nameEaxis
    if (!ParamIsDefault(sourceEaxis) && WaveExists(sourceEaxis) && (numpnts(sourceEaxis) == NE))
        SetScale/P x, DimOffset(sourceEaxis,0), DimDelta(sourceEaxis,0), WaveUnits(sourceEaxis,0), E_axis
    elseif (WaveExists(defaultSourceEaxis) && (numpnts(defaultSourceEaxis) == NE))
        SetScale/P x, DimOffset(defaultSourceEaxis,0), DimDelta(defaultSourceEaxis,0), WaveUnits(defaultSourceEaxis,0), E_axis
    else
        SetScale/I x, -2.5e-3, 2.5e-3, "eV", E_axis
    endif
    E_axis = x

    Make/O/D/N=(NE, nB) $nameDOS
    Wave DOS_Manifold = $nameDOS
    DOS_Manifold = 0
    SetScale/P x, DimOffset(E_axis,0), DimDelta(E_axis,0), "eV", DOS_Manifold
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", DOS_Manifold

    String nameCl = outPrefix + "_clusterSize"
    String nameShift = outPrefix + "_hybridShift_eV"
    String namePartner = outPrefix + "_partner"
    String nameVij = outPrefix + "_Vij_eV"
    String nameSlope = outPrefix + "_slope_eVperT"
    Make/O/D/N=(nRec, nB) $nameCl
    Make/O/D/N=(nRec, nB) $nameShift
    Make/O/D/N=(nRec, nB) $namePartner
    Make/O/D/N=(nRec, nB) $nameVij
    Make/O/D/N=(nRec, nB) $nameSlope
    Wave clusterSize = $nameCl
    Wave hybridShift = $nameShift
    Wave hybridPartner = $namePartner
    Wave hybridVij = $nameVij
    Wave BranchSlope = $nameSlope
    clusterSize = NaN
    hybridShift = NaN
    hybridPartner = NaN
    hybridVij = NaN
    BranchSlope = NaN
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", clusterSize
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", hybridShift
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", hybridPartner
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", hybridVij
    SetScale/P y, DimOffset(B_T,0), DimDelta(B_T,0), "T", BranchSlope

    Make/FREE/D/N=(nRec) Branch_phi
    if (WaveExists(Branch_phi_existing) && (numpnts(Branch_phi_existing) == nRec))
        Branch_phi = Branch_phi_existing[p]
    else
        Branch_phi = SNS__ChordPhiFromEndpoints_rad(Branch_x1[p], Branch_y1[p], Branch_x2[p], Branch_y2[p])
    endif

    Variable rec, rec2, iB, iE
    Variable B0seg, B1seg, dBseg, Eprev, Enext, Ei, Ej
    for (rec = 0; rec < nRec; rec += 1)
        for (iB = 0; iB < nB; iB += 1)
            if (iB == 0)
                B0seg = B_T[0]
                B1seg = B_T[1]
                Eprev = Branch_E[rec][0]
                Enext = Branch_E[rec][1]
            elseif (iB == nB-1)
                B0seg = B_T[nB-2]
                B1seg = B_T[nB-1]
                Eprev = Branch_E[rec][nB-2]
                Enext = Branch_E[rec][nB-1]
            else
                B0seg = B_T[iB-1]
                B1seg = B_T[iB+1]
                Eprev = Branch_E[rec][iB-1]
                Enext = Branch_E[rec][iB+1]
            endif
            dBseg = B1seg - B0seg
            if (dBseg != 0 && !numtype(Eprev) && !numtype(Enext))
                BranchSlope[rec][iB] = (Enext - Eprev) / dBseg
            endif
        endfor
    endfor

    Make/FREE/B/N=(nRec) assigned
    Make/FREE/D/N=(maxCluster) candRec, candScore
    Make/FREE/D/N=(maxCluster, maxCluster) H, Hwork, evec
    Make/FREE/D/N=(maxCluster) eval, eigWeight

    Variable nClusters = 0
    Variable nCand, pos, insertAt, cIdx, cJdx, cand
    Variable dE, dSlope, dSlopeNorm, Vij, Veff, score, worstScore, worstIdx
    Variable maxVforRec, bestPartner, k, a, b
    Variable eig, wEig, shiftMin, dEgrid, alpha
    Variable slopeOK

    for (iB = 0; iB < nB; iB += 1)
        assigned = 0

        for (rec = 0; rec < nRec; rec += 1)
            if (assigned[rec])
                continue
            endif

            Ei = Branch_E[rec][iB]
            if (numtype(Ei))
                assigned[rec] = 1
                continue
            endif

            nCand = 1
            candRec[0] = rec
            candScore[0] = -Inf

            for (rec2 = 0; rec2 < nRec; rec2 += 1)
                if ((rec2 == rec) || assigned[rec2])
                    continue
                endif

                Ej = Branch_E[rec2][iB]
                if (numtype(Ej))
                    continue
                endif

                dE = abs(Ei - Ej)
                if (dE > clusterEnergy_eV)
                    continue
                endif

                Vij = SNS__HybridCoupling_ij(rec, rec2, Branch_L, Branch_W, Branch_T, Branch_x1, Branch_y1, Branch_x2, Branch_y2, Branch_phi, V0_eV, sigmaW_nm, sigmaL_nm, sigmaEnd_nm, sigmaPhi_rad)
                if (Vij <= 0 || numtype(Vij))
                    continue
                endif

                slopeOK = 1
                dSlopeNorm = 0
                if (sigmaSlope_eVperT > 0)
                    if (numtype(BranchSlope[rec][iB]) || numtype(BranchSlope[rec2][iB]))
                        slopeOK = 0
                    else
                        dSlope = BranchSlope[rec][iB] - BranchSlope[rec2][iB]
                        dSlopeNorm = dSlope / sigmaSlope_eVperT
                    endif
                endif
                if (!slopeOK)
                    continue
                endif

                Veff = Vij * exp(-0.5*(dE/clusterEnergy_eV)^2) * exp(-0.5*dSlopeNorm*dSlopeNorm)
                if (Veff <= 0 || numtype(Veff))
                    continue
                endif
                if ((minVijActive_eV > 0) && (Veff < minVijActive_eV))
                    continue
                endif

                score = Veff / (1 + dE/clusterEnergy_eV + abs(dSlopeNorm))

                if (nCand < maxCluster)
                    candRec[nCand] = rec2
                    candScore[nCand] = score
                    nCand += 1
                else
                    worstIdx = 1
                    worstScore = candScore[1]
                    for (pos = 2; pos < nCand; pos += 1)
                        if (candScore[pos] < worstScore)
                            worstScore = candScore[pos]
                            worstIdx = pos
                        endif
                    endfor
                    if (score > worstScore)
                        candRec[worstIdx] = rec2
                        candScore[worstIdx] = score
                    endif
                endif
            endfor

            if (nCand < 2)
                for (iE = 0; iE < NE; iE += 1)
                    dE = E_axis[iE] - Ei
                    DOS_Manifold[iE][iB] += Branch_weight[rec] * (GammaHybrid_eV/pi) / (dE*dE + GammaHybrid_eV*GammaHybrid_eV)
                endfor
                clusterSize[rec][iB] = 1
                hybridShift[rec][iB] = 0
                assigned[rec] = 1
                continue
            endif

            H = 0
            Hwork = 0
            for (cIdx = 0; cIdx < nCand; cIdx += 1)
                cand = candRec[cIdx]
                H[cIdx][cIdx] = Branch_E[cand][iB]
            endfor

            for (cIdx = 0; cIdx < nCand-1; cIdx += 1)
                a = candRec[cIdx]
                for (cJdx = cIdx + 1; cJdx < nCand; cJdx += 1)
                    b = candRec[cJdx]
                    dE = abs(Branch_E[a][iB] - Branch_E[b][iB])
                    dSlopeNorm = 0
                    if (sigmaSlope_eVperT > 0)
                        dSlopeNorm = (BranchSlope[a][iB] - BranchSlope[b][iB]) / sigmaSlope_eVperT
                    endif
                    Vij = SNS__HybridCoupling_ij(a, b, Branch_L, Branch_W, Branch_T, Branch_x1, Branch_y1, Branch_x2, Branch_y2, Branch_phi, V0_eV, sigmaW_nm, sigmaL_nm, sigmaEnd_nm, sigmaPhi_rad)
                    Veff = Vij * exp(-0.5*(dE/clusterEnergy_eV)^2) * exp(-0.5*dSlopeNorm*dSlopeNorm)
                    if ((minVijActive_eV > 0) && (Veff < minVijActive_eV))
                        Veff = 0
                    endif
                    H[cIdx][cJdx] = Veff
                    H[cJdx][cIdx] = Veff
                endfor
            endfor

            Hwork = H
            SNS__JacobiEigenSymmetric(Hwork, eval, evec, nCand)

            for (eig = 0; eig < nCand; eig += 1)
                wEig = 0
                for (cIdx = 0; cIdx < nCand; cIdx += 1)
                    cand = candRec[cIdx]
                    wEig += Branch_weight[cand] * evec[cIdx][eig] * evec[cIdx][eig]
                endfor

                for (iE = 0; iE < NE; iE += 1)
                    dE = E_axis[iE] - eval[eig]
                    DOS_Manifold[iE][iB] += wEig * (GammaHybrid_eV/pi) / (dE*dE + GammaHybrid_eV*GammaHybrid_eV)
                endfor
            endfor

            for (cIdx = 0; cIdx < nCand; cIdx += 1)
                cand = candRec[cIdx]
                shiftMin = Inf
                for (eig = 0; eig < nCand; eig += 1)
                    dE = abs(eval[eig] - Branch_E[cand][iB])
                    if (dE < shiftMin)
                        shiftMin = dE
                    endif
                endfor

                maxVforRec = 0
                bestPartner = -1
                for (cJdx = 0; cJdx < nCand; cJdx += 1)
                    if (cJdx == cIdx)
                        continue
                    endif
                    if (abs(H[cIdx][cJdx]) > maxVforRec)
                        maxVforRec = abs(H[cIdx][cJdx])
                        bestPartner = candRec[cJdx]
                    endif
                endfor

                clusterSize[cand][iB] = nCand
                hybridShift[cand][iB] = shiftMin
                hybridPartner[cand][iB] = bestPartner
                hybridVij[cand][iB] = maxVforRec
                assigned[cand] = 1
            endfor

            nClusters += 1
        endfor
    endfor

    dEgrid = E_axis[1] - E_axis[0]
    Variable Emin = E_axis[0]
    Variable Emax = E_axis[NE - 1]
    Variable Egap_low = max(Emin, -Delta)
    Variable Egap_high = min(Emax, Delta)
    Variable EgapSpan = Egap_high - Egap_low
    Variable acc, E, contShapeVal

    if (N_cont_statesPer_eV > 0 && dEgrid > 0 && EgapSpan > 0)
        Variable N_gap_ref = N_cont_statesPer_eV * EgapSpan
        for (iB = 0; iB < nB; iB += 1)
            acc = 0
            for (iE = 0; iE < NE; iE += 1)
                acc += DOS_Manifold[iE][iB]
            endfor
            acc *= dEgrid
            if (acc > 0 && numtype(acc) == 0)
                alpha = N_gap_ref / acc
                for (iE = 0; iE < NE; iE += 1)
                    DOS_Manifold[iE][iB] *= alpha
                endfor
            endif
        endfor

        Variable E_sm = 3*GammaHybrid_eV
        if (E_sm <= 0 || numtype(E_sm))
            E_sm = 3e-9
        endif
        for (iE = 0; iE < NE; iE += 1)
            E = E_axis[iE]
            contShapeVal = 0.5 * (1 + tanh((abs(E) - Delta) / E_sm))
            for (iB = 0; iB < nB; iB += 1)
                DOS_Manifold[iE][iB] += N_cont_statesPer_eV * contShapeVal
            endfor
        endfor
    endif

    String meta
    meta  = "SNS_DOSBuilder=BranchTableManifoldHybridized;"
    meta += "SNS_BranchTablePrefix=" + branchPrefix + ";"
    meta += "SNS_V0_eV=" + num2str(V0_eV) + ";"
    meta += "SNS_clusterEnergy_eV=" + num2str(clusterEnergy_eV) + ";"
    meta += "SNS_sigmaSlope_eVperT=" + num2str(sigmaSlope_eVperT) + ";"
    meta += "SNS_sigmaW_nm=" + num2str(sigmaW_nm) + ";"
    meta += "SNS_sigmaL_nm=" + num2str(sigmaL_nm) + ";"
    meta += "SNS_sigmaEnd_nm=" + num2str(sigmaEnd_nm) + ";"
    meta += "SNS_sigmaPhi_rad=" + num2str(sigmaPhi_rad) + ";"
    meta += "SNS_maxCluster=" + num2str(maxCluster) + ";"
    meta += "SNS_GammaHybrid_eV=" + num2str(GammaHybrid_eV) + ";"
    meta += "SNS_Delta_eV=" + num2str(Delta) + ";"
    meta += "SNS_N_cont_statesPer_eV=" + num2str(N_cont_statesPer_eV) + ";"
    meta += "SNS_minVijActive_eV=" + num2str(minVijActive_eV) + ";"

    Note/K DOS_Manifold
    Note DOS_Manifold, meta
    Note/K clusterSize
    Note clusterSize, meta
    Note/K hybridShift
    Note hybridShift, meta
    Note/K hybridPartner
    Note hybridPartner, meta
    Note/K hybridVij
    Note hybridVij, meta
    Note/K BranchSlope
    Note BranchSlope, meta

    SNS_Log("SNS_ComputeDOS_FromBranchTable_ManifoldHybridized: finished; clusters=" + num2str(nClusters), level="INFO")
    return nClusters
End

//==============================================================================
// SNS_ComputeDOS_FromSettings_TwoPass
//
// Purpose:
//   High-level wrapper for SNS_ComputeDOS_FromChannels_TwoPass(...).
//   Reads canonical channel waves from dataFolder and creates the two-pass DOS
//   and branch-density diagnostics in the caller's current data folder.
//
// Inputs:
//   dataFolder          : folder containing canonical channel geometry waves and
//                         w_mask.
//   nameDOS             : output DOS base name.
//   nameEaxis           : output energy-axis wave name.
//
// Optional Inputs:
//   h_eff_3D_nm         : if supplied, use 3D continuum DOS volume estimate.
//   betaExtra_List      : optional 1D betaExtra_List[ch].
//   sigmaE_eV           : energy kernel width for rho.
//   sigmaW_nm           : W_eff kernel width for rho.
//   GammaMixMax_eV      : maximum candidate mixing broadening.
//   rho0                : saturation density scale.
//   useMixBroadening    : nonzero to apply Gamma_mix in quadrature.
//   outPrefix           : branch/rho diagnostic output prefix.
//
// Returns:
//   Number of flattened branch records.
//==============================================================================
Function SNS_ComputeDOS_FromSettings_TwoPass(dataFolder, nameDOS, nameEaxis, [h_eff_3D_nm, betaExtra_List, sigmaE_eV, sigmaW_nm, GammaMixMax_eV, rho0, useMixBroadening, outPrefix])
    String dataFolder
    String nameDOS, nameEaxis
    Variable h_eff_3D_nm
    Wave betaExtra_List
    Variable sigmaE_eV, sigmaW_nm, GammaMixMax_eV, rho0, useMixBroadening
    String outPrefix

    String savedDF = GetDataFolder(1)

    STRUCT SNS_Params SNS_p
    SNS_LoadParams(SNS_p)

    DFREF dfrSet  = root:SNS_Settings
    DFREF dfrGeom = $dataFolder

    Wave B_T = dfrSet:B_T

    Wave L_N_List_nm   = dfrGeom:L_N_List_nm
    Wave W_eff_List_nm = dfrGeom:W_eff_List_nm
    Wave/Z phiList_rad = dfrGeom:phiList_rad
    Wave wChan         = dfrGeom:wChan
    Wave T_eff_List    = dfrGeom:T_eff_List

    Wave Hit1x_List_nm = dfrGeom:Hit1x_List_nm
    Wave Hit1y_List_nm = dfrGeom:Hit1y_List_nm
    Wave Hit2x_List_nm = dfrGeom:Hit2x_List_nm
    Wave Hit2y_List_nm = dfrGeom:Hit2y_List_nm

    Wave/Z w_mask = dfrGeom:w_mask
    if (!WaveExists(w_mask))
        SetDataFolder $savedDF
        Abort "SNS_ComputeDOS_FromSettings_TwoPass: wave w_mask is missing."
    endif
    if (!WaveExists(phiList_rad))
        SetDataFolder $savedDF
        Abort "SNS_ComputeDOS_FromSettings_TwoPass: wave phiList_rad is missing."
    endif

    SetDataFolder dfrGeom
    Wave/Z w_area_nm2 = SNS_MaskAreaPerim_FromParticles(w_mask)
    SetDataFolder $savedDF

    if (!WaveExists(w_area_nm2))
        Abort "SNS_ComputeDOS_FromSettings_TwoPass: could not create/read w_area_nm2."
    endif

    Variable area_nm2 = w_area_nm2[0]
    Variable area_m2  = area_nm2 * 1e-18

    Variable N_cont_statesPer_eV
    if (ParamIsDefault(h_eff_3D_nm))
        N_cont_statesPer_eV = SNS_p.DOS2D_eV_Area * area_m2
    else
        N_cont_statesPer_eV = SNS_p.DOS3D_eV_Vol * area_m2 * h_eff_3D_nm * 1e-9
    endif

    Wave/Z Vortex_ptx = dfrGeom:Vortex_ptx
    Wave/Z Vortex_pty = dfrGeom:Vortex_pty

    Variable xV_nm = 0
    Variable yV_nm = 0
    if (WaveExists(Vortex_ptx) && WaveExists(Vortex_pty))
        xV_nm = Vortex_ptx[0]
        yV_nm = Vortex_pty[0]
    endif

    NVAR SNS_useVortex = dfrSet:SNS_useVortex
    NVAR SNS_nFlux     = dfrSet:SNS_nFlux

    Variable nChIn = numpnts(L_N_List_nm)
    Variable nRec
    Variable is3D = !ParamIsDefault(h_eff_3D_nm)

    SetDataFolder $savedDF

    if (!ParamIsDefault(betaExtra_List))

        nRec = SNS_ComputeDOS_FromChannels_TwoPass( \
            B_T, \
            L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
            Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
            xV_nm, yV_nm, \
            N_cont_statesPer_eV, \
            nameDOS, nameEaxis, \
            betaExtra_List=betaExtra_List, \
            is3D=is3D, \
            sigmaE_eV=sigmaE_eV, \
            sigmaW_nm=sigmaW_nm, \
            GammaMixMax_eV=GammaMixMax_eV, \
            rho0=rho0, \
            useMixBroadening=useMixBroadening, \
            outPrefix=outPrefix, \
            phiList_rad=phiList_rad)

    elseif (SNS_useVortex)

        Make/FREE/D/N=(nChIn) betaExtra_local

        Variable j, th1, th2, dth
        for (j = 0; j < nChIn; j += 1)
            th1 = atan2(Hit1y_List_nm[j] - yV_nm, Hit1x_List_nm[j] - xV_nm)
            th2 = atan2(Hit2y_List_nm[j] - yV_nm, Hit2x_List_nm[j] - xV_nm)
            dth = th2 - th1
            if (dth > pi)
                dth -= 2*pi
            elseif (dth < -pi)
                dth += 2*pi
            endif
            betaExtra_local[j] = SNS_nFlux * dth
        endfor

        nRec = SNS_ComputeDOS_FromChannels_TwoPass( \
            B_T, \
            L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
            Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
            xV_nm, yV_nm, \
            N_cont_statesPer_eV, \
            nameDOS, nameEaxis, \
            betaExtra_List=betaExtra_local, \
            is3D=is3D, \
            sigmaE_eV=sigmaE_eV, \
            sigmaW_nm=sigmaW_nm, \
            GammaMixMax_eV=GammaMixMax_eV, \
            rho0=rho0, \
            useMixBroadening=useMixBroadening, \
            outPrefix=outPrefix, \
            phiList_rad=phiList_rad)

    else

        nRec = SNS_ComputeDOS_FromChannels_TwoPass( \
            B_T, \
            L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
            Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
            xV_nm, yV_nm, \
            N_cont_statesPer_eV, \
            nameDOS, nameEaxis, \
            is3D=is3D, \
            sigmaE_eV=sigmaE_eV, \
            sigmaW_nm=sigmaW_nm, \
            GammaMixMax_eV=GammaMixMax_eV, \
            rho0=rho0, \
            useMixBroadening=useMixBroadening, \
            outPrefix=outPrefix, \
            phiList_rad=phiList_rad)

    endif

    SetDataFolder $savedDF
    return nRec
End

//==============================================================================
// SNS_ComputeDOS_FromSettings_Green1D
//
// Purpose:
//   Folder-based wrapper for SNS_ComputeDOS_FromChannels_Green1D().
//   It uses the same canonical channel waves as SNS_ComputeDOS_FromSettings()
//   and constructs the same vortex betaExtra list when vortex mode is enabled.
//==============================================================================
Function SNS_ComputeDOS_FromSettings_Green1D(dataFolder, nameDOS, nameEaxis, [betaExtra_List, Delta_fit_eV, Gamma0_eV, etaN_eV, etaS_eV, NsiteMax, avgSites])
    String dataFolder
    String nameDOS, nameEaxis
    Wave betaExtra_List
    Variable Delta_fit_eV, Gamma0_eV, etaN_eV, etaS_eV, NsiteMax, avgSites

    String savedDF = GetDataFolder(1)

    STRUCT SNS_Params SNS_p
    SNS_LoadParams(SNS_p)

    Variable DeltaUse = SNS_p.Delta
    if (!ParamIsDefault(Delta_fit_eV))
        DeltaUse = Delta_fit_eV
    endif

    Variable Gamma0Use = DeltaUse
    if (!ParamIsDefault(Gamma0_eV))
        Gamma0Use = Gamma0_eV
    endif

    Variable etaNUse = SNS_p.Broadening
    if (!ParamIsDefault(etaN_eV))
        etaNUse = etaN_eV
    endif

    Variable etaSUse = max(SNS_p.Broadening, 1e-9)
    if (!ParamIsDefault(etaS_eV))
        etaSUse = etaS_eV
    endif

    Variable NsiteUse = 17
    if (!ParamIsDefault(NsiteMax))
        NsiteUse = NsiteMax
    endif

    Variable avgSitesUse = 7
    if (!ParamIsDefault(avgSites))
        avgSitesUse = avgSites
    endif

    DFREF dfrSet  = root:SNS_Settings
    DFREF dfrGeom = $dataFolder

    Wave B_T = dfrSet:B_T

    Wave L_N_List_nm   = dfrGeom:L_N_List_nm
    Wave W_eff_List_nm = dfrGeom:W_eff_List_nm
    Wave wChan         = dfrGeom:wChan
    Wave T_eff_List    = dfrGeom:T_eff_List

    Wave/Z Hit1x_List_nm = dfrGeom:Hit1x_List_nm
    Wave/Z Hit1y_List_nm = dfrGeom:Hit1y_List_nm
    Wave/Z Hit2x_List_nm = dfrGeom:Hit2x_List_nm
    Wave/Z Hit2y_List_nm = dfrGeom:Hit2y_List_nm

    Wave/Z Vortex_ptx = dfrGeom:Vortex_ptx
    Wave/Z Vortex_pty = dfrGeom:Vortex_pty

    Variable xV_nm = 0
    Variable yV_nm = 0
    if (WaveExists(Vortex_ptx) && WaveExists(Vortex_pty))
        xV_nm = Vortex_ptx[0]
        yV_nm = Vortex_pty[0]
    endif

    NVAR SNS_useVortex = dfrSet:SNS_useVortex
    NVAR SNS_nFlux     = dfrSet:SNS_nFlux

    Variable nChIn = numpnts(L_N_List_nm)
    Variable nCh

    SetDataFolder $savedDF

    if (!ParamIsDefault(betaExtra_List))

        nCh = SNS_ComputeDOS_FromChannels_Green1D( \
            B_T, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
            nameDOS, nameEaxis, \
            betaExtra_List=betaExtra_List, \
            Delta_fit_eV=DeltaUse, \
            Gamma0_eV=Gamma0Use, \
            etaN_eV=etaNUse, \
            etaS_eV=etaSUse, \
            NsiteMax=NsiteUse, \
            avgSites=avgSitesUse)

    elseif (SNS_useVortex)

        if (!WaveExists(Hit1x_List_nm) || !WaveExists(Hit1y_List_nm) || !WaveExists(Hit2x_List_nm) || !WaveExists(Hit2y_List_nm))
            Abort "SNS_ComputeDOS_FromSettings_Green1D: endpoint waves required for vortex betaExtra."
        endif

        Make/FREE/D/N=(nChIn) betaExtra_local

        Variable j
        Variable th1, th2, dth

        for (j = 0; j < nChIn; j += 1)
            th1 = atan2(Hit1y_List_nm[j] - yV_nm, Hit1x_List_nm[j] - xV_nm)
            th2 = atan2(Hit2y_List_nm[j] - yV_nm, Hit2x_List_nm[j] - xV_nm)
            dth = th2 - th1

            if (dth > pi)
                dth -= 2*pi
            elseif (dth < -pi)
                dth += 2*pi
            endif

            betaExtra_local[j] = SNS_nFlux * dth
        endfor

        nCh = SNS_ComputeDOS_FromChannels_Green1D( \
            B_T, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
            nameDOS, nameEaxis, \
            betaExtra_List=betaExtra_local, \
            Delta_fit_eV=DeltaUse, \
            Gamma0_eV=Gamma0Use, \
            etaN_eV=etaNUse, \
            etaS_eV=etaSUse, \
            NsiteMax=NsiteUse, \
            avgSites=avgSitesUse)

    else

        nCh = SNS_ComputeDOS_FromChannels_Green1D( \
            B_T, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
            nameDOS, nameEaxis, \
            Delta_fit_eV=DeltaUse, \
            Gamma0_eV=Gamma0Use, \
            etaN_eV=etaNUse, \
            etaS_eV=etaSUse, \
            NsiteMax=NsiteUse, \
            avgSites=avgSitesUse)

    endif

    SetDataFolder $savedDF
    return nCh
End

//==============================================================================
// SNS_PrintChannelWaveCounts
//
// Purpose:
//   Print canonical channel-wave lengths in a geometry data folder.
//   Useful when SNS_ComputeDOS_FromChannels_TwoPass reports inconsistent
//   channel waves.
//
// Inputs:
//   dataFolder : folder containing channel waves.
//
// Returns:
//   0
//==============================================================================
Function SNS_PrintChannelWaveCounts(dataFolder)
    String dataFolder

    DFREF dfrGeom = $dataFolder

    Wave/Z L_N_List_nm   = dfrGeom:L_N_List_nm
    Wave/Z W_eff_List_nm = dfrGeom:W_eff_List_nm
    Wave/Z phiList_rad   = dfrGeom:phiList_rad
    Wave/Z wChan         = dfrGeom:wChan
    Wave/Z T_eff_List    = dfrGeom:T_eff_List
    Wave/Z Hit1x_List_nm = dfrGeom:Hit1x_List_nm
    Wave/Z Hit1y_List_nm = dfrGeom:Hit1y_List_nm
    Wave/Z Hit2x_List_nm = dfrGeom:Hit2x_List_nm
    Wave/Z Hit2y_List_nm = dfrGeom:Hit2y_List_nm

    Print "SNS channel wave counts for ", dataFolder
    if (WaveExists(L_N_List_nm))
        Print "  L_N_List_nm   = ", numpnts(L_N_List_nm)
    else
        Print "  L_N_List_nm   = MISSING"
    endif
    if (WaveExists(W_eff_List_nm))
        Print "  W_eff_List_nm = ", numpnts(W_eff_List_nm)
    else
        Print "  W_eff_List_nm = MISSING"
    endif
    if (WaveExists(phiList_rad))
        Print "  phiList_rad   = ", numpnts(phiList_rad)
    else
        Print "  phiList_rad   = MISSING"
    endif
    if (WaveExists(wChan))
        Print "  wChan         = ", numpnts(wChan)
    else
        Print "  wChan         = MISSING"
    endif
    if (WaveExists(T_eff_List))
        Print "  T_eff_List    = ", numpnts(T_eff_List)
    else
        Print "  T_eff_List    = MISSING"
    endif
    if (WaveExists(Hit1x_List_nm))
        Print "  Hit1x_List_nm = ", numpnts(Hit1x_List_nm)
    else
        Print "  Hit1x_List_nm = MISSING"
    endif
    if (WaveExists(Hit1y_List_nm))
        Print "  Hit1y_List_nm = ", numpnts(Hit1y_List_nm)
    else
        Print "  Hit1y_List_nm = MISSING"
    endif
    if (WaveExists(Hit2x_List_nm))
        Print "  Hit2x_List_nm = ", numpnts(Hit2x_List_nm)
    else
        Print "  Hit2x_List_nm = MISSING"
    endif
    if (WaveExists(Hit2y_List_nm))
        Print "  Hit2y_List_nm = ", numpnts(Hit2y_List_nm)
    else
        Print "  Hit2y_List_nm = MISSING"
    endif

    return 0
End

//==============================================================================
// SNS_ComputeDOS_FromSettings
//
// Purpose:
//   High-level wrapper to compute DOS(E,B) in absolute units [states/eV]
//   for a given SNS channel geometry stored in a data folder.
//
//   This function:
//     • Loads global SNS parameters from SNS_Settings via SNS_LoadParams()
//       (Delta, vF, lambdaL, Broadening, NE, DOS2D_eV_Area,
//        DOS3D_eV_Vol, …).
//     • Uses the global magnetic-field axis B_T stored in root:SNS_Settings.
//     • Reads channel geometry and weights from `dataFolder`:
//         L_N_List_nm, W_eff_List_nm, wChan, T_eff_List,
//         Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm.
//     • Reads the N-region mask w_mask and obtains the N-region area using
//         SNS_MaskAreaPerim_FromParticles(w_mask).
//     • Constructs the scalar continuum DOS contribution
//         N_cont_statesPer_eV
//       from either the 2D normal-state DOS or, when h_eff_3D_nm is supplied,
//       the 3D normal-state DOS times the effective volume.
//     • Optionally reads a vortex center from
//         Vortex_ptx / Vortex_pty.
//     • Reads vortex controls from root:SNS_Settings:
//         SNS_useVortex, SNS_nFlux, SNS_nIntSteps, SNS_useDeltamap.
//       In the present implementation only SNS_useVortex and SNS_nFlux are
//       used directly in this wrapper.
//     • If a precomputed betaExtra_List is supplied, forwards it directly to
//         SNS_ComputeDOS_FromChannels(...).
//       This takes precedence over internally generated vortex phases.
//     • If no betaExtra_List is supplied but SNS_useVortex is true, constructs
//       a controlled endpoint vortex phase
//         betaExtra_local[ch]
//       from the difference of endpoint polar angles around the vortex center,
//       wraps that difference once into [-pi, pi], multiplies by SNS_nFlux,
//       and passes the resulting phase list explicitly to
//         SNS_ComputeDOS_FromChannels(...).
//       This avoids triggering lower-level legacy vortex handling implicitly.
//     • If neither betaExtra_List nor SNS_useVortex is active, calls
//         SNS_ComputeDOS_FromChannels(...)
//       without an extra phase wave.
//     • Creates output waves in the caller's current data folder, not inside
//       `dataFolder`.
//
// Inputs:
//   dataFolder          : path/name of data folder containing the channel
//                         geometry and mask waves.
//
//                         Required waves:
//                           L_N_List_nm   [ch]    N-region path length [nm]
//                           W_eff_List_nm [ch]    effective orbital lever arm [nm]
//                           wChan         [ch]    channel weight
//                           T_eff_List    [ch]    effective transparency
//                           Hit1x_List_nm [ch]    first endpoint x coordinate [nm]
//                           Hit1y_List_nm [ch]    first endpoint y coordinate [nm]
//                           Hit2x_List_nm [ch]    second endpoint x coordinate [nm]
//                           Hit2y_List_nm [ch]    second endpoint y coordinate [nm]
//                           w_mask                binary/particle mask of N region
//
//                         Optional waves:
//                           Vortex_ptx            vortex x coordinate [nm]
//                           Vortex_pty            vortex y coordinate [nm]
//
//   nameDOS             : output DOS wave name to create in the caller folder.
//                         This function does not modify the supplied name.
//
//   nameEaxis           : output energy-axis wave name [eV] to create in the
//                         caller folder. This function does not modify the
//                         supplied name.
//
// Optional Inputs:
//   h_eff_3D_nm         : effective height of the N region contributing to
//                         surface LDOS for a 3D S-N-I / 3D-carrier
//                         implementation [nm].
//                         If omitted:
//                           N_cont_statesPer_eV = DOS2D_eV_Area * area.
//                         If supplied:
//                           N_cont_statesPer_eV = DOS3D_eV_Vol * area
//                                                * h_eff_3D_nm.
//
//   betaExtra_List      : precomputed per-channel extra phase [rad].
//
//                         Accepted forms:
//                           • 1D wave betaExtra_List[ch]
//                             -> field-independent extra phase
//
//                           • 2D wave betaExtra_List[ch][iB]
//                             -> field-dependent extra phase
//
//                         If supplied, this wave is passed directly to
//                         SNS_ComputeDOS_FromChannels(...), overriding the
//                         internally generated endpoint vortex phase.
//
// Outputs:
//   nameDOS             : DOS(E,B) in [states/eV], created in the caller's
//                         current data folder by SNS_ComputeDOS_FromChannels.
//
//   nameDOS + "_Gamma"  : additionally broadened DOS(E,B), if generated by
//                         SNS_ComputeDOS_FromChannels.
//
//   nameEaxis           : energy axis [eV], created in the caller's current
//                         data folder.
//
// Returns:
//   nCh                 : number of input channels passed to
//                         SNS_ComputeDOS_FromChannels(...).
//
// Notes:
//   • This wrapper intentionally does not alter output wave names. Any naming
//     convention distinguishing vortex/non-vortex runs must be handled by the
//     caller through nameDOS/nameEaxis.
//   • The vortex phase generated here is endpoint-based:
//         betaExtra_local[ch] = SNS_nFlux *
//             wrapToPi( atan2(Hit2 - Vortex) - atan2(Hit1 - Vortex) )
//     with a single bounded wrap step, avoiding while-loop phase unwrapping.
//   • If Vortex_ptx/Vortex_pty are absent and SNS_useVortex is true, the vortex
//     center defaults to (0,0). For strict input checking, validate dataFolder
//     before calling this function.
//   • The original data folder is restored before normal return and before the
//     explicit aborts in this wrapper.
//==============================================================================
Function SNS_ComputeDOS_FromSettings(dataFolder, nameDOS, nameEaxis, [h_eff_3D_nm, betaExtra_List, branchNormalize])

    String dataFolder
    String nameDOS, nameEaxis
    Variable h_eff_3D_nm
    Wave betaExtra_List
    Variable branchNormalize

    String savedDF = GetDataFolder(1)

    // ---- load global SNS settings ----
    STRUCT SNS_Params SNS_p
    SNS_LoadParams(SNS_p)

    // ---- use data-folder references, avoid leaking current DF state ----
    DFREF dfrSet  = root:SNS_Settings
    DFREF dfrGeom = $dataFolder

    Wave B_T = dfrSet:B_T

    Wave L_N_List_nm   = dfrGeom:L_N_List_nm
    Wave W_eff_List_nm = dfrGeom:W_eff_List_nm
    Wave wChan      = dfrGeom:wChan
    Wave T_eff_List = dfrGeom:T_eff_List

    Wave Hit1x_List_nm = dfrGeom:Hit1x_List_nm
    Wave Hit1y_List_nm = dfrGeom:Hit1y_List_nm
    Wave Hit2x_List_nm = dfrGeom:Hit2x_List_nm
    Wave Hit2y_List_nm = dfrGeom:Hit2y_List_nm

    Wave/Z w_mask = dfrGeom:w_mask
    if (!WaveExists(w_mask))
        String msg1 = "SNS_ComputeDOS_FromSettings: wave w_mask is missing."
        SNS_Log(msg1, level="ERR")
        SetDataFolder $savedDF
        Abort msg1
    endif

    // area / perimeter helper waves
    SetDataFolder dfrGeom
    Wave/Z w_area_nm2 = SNS_MaskAreaPerim_FromParticles(w_mask)
    SetDataFolder $savedDF

    if (!WaveExists(w_area_nm2))
        String msg2 = "SNS_ComputeDOS_FromSettings: could not create/read w_area_nm2."
        SNS_Log(msg2, level="ERR")
        Abort msg2
    endif

    // ---- continuum states/eV for this patch ----
    Variable area_nm2 = w_area_nm2[0]
    Variable area_m2  = area_nm2 * 1e-18

    Variable N_cont_statesPer_eV
    if (ParamIsDefault(h_eff_3D_nm))
        N_cont_statesPer_eV = SNS_p.DOS2D_eV_Area * area_m2
    else
        N_cont_statesPer_eV = SNS_p.DOS3D_eV_Vol * area_m2 * h_eff_3D_nm * 1e-9
    endif

    // ---- vortex center, if present ----
    Wave/Z Vortex_ptx = dfrGeom:Vortex_ptx
    Wave/Z Vortex_pty = dfrGeom:Vortex_pty

    Variable xV_nm = 0
    Variable yV_nm = 0

    if (WaveExists(Vortex_ptx) && WaveExists(Vortex_pty))
        xV_nm = Vortex_ptx[0]
        yV_nm = Vortex_pty[0]
    endif

    // ---- vortex controls ----
    NVAR SNS_useVortex   = dfrSet:SNS_useVortex
    NVAR SNS_nFlux       = dfrSet:SNS_nFlux
    NVAR SNS_nIntSteps   = dfrSet:SNS_nIntSteps
    NVAR SNS_useDeltamap = dfrSet:SNS_useDeltamap

    Variable useVortex = SNS_useVortex
    Variable nFlux     = SNS_nFlux

    Variable nChIn = numpnts(L_N_List_nm)
    Variable nCh

    Variable is3D = !ParamIsDefault(h_eff_3D_nm)
    Variable useBranchNormalize = 0
    if (!ParamIsDefault(branchNormalize))
        useBranchNormalize = (branchNormalize != 0)
    endif

    // ---- output should be created in caller folder ----
    SetDataFolder $savedDF

    //==========================================================================
    // Critical part:
    // If caller supplied betaExtra_List, use it.
    // Else, if vortex is enabled, create a controlled endpoint phase list here.
    // Else, call with no betaExtra_List.
    //==========================================================================
    if (!ParamIsDefault(betaExtra_List))

        nCh = SNS_ComputeDOS_FromChannels( \
            B_T, \
            L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
            Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
            xV_nm, yV_nm, \
            N_cont_statesPer_eV, \
            nameDOS, nameEaxis, \
            betaExtra_List=betaExtra_List, \
            is3D=is3D, \
            branchNormalize=useBranchNormalize)

    elseif (useVortex)

        Make/FREE/D/N=(nChIn) betaExtra_local

        Variable j
        Variable th1, th2, dth

        for (j = 0; j < nChIn; j += 1)

            th1 = atan2(Hit1y_List_nm[j] - yV_nm, Hit1x_List_nm[j] - xV_nm)
            th2 = atan2(Hit2y_List_nm[j] - yV_nm, Hit2x_List_nm[j] - xV_nm)

            dth = th2 - th1

            // wrap once into [-pi, pi]; no while-loop
            if (dth > pi)
                dth -= 2*pi
            elseif (dth < -pi)
                dth += 2*pi
            endif

            betaExtra_local[j] = nFlux * dth
        endfor

        nCh = SNS_ComputeDOS_FromChannels( \
            B_T, \
            L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
            Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
            xV_nm, yV_nm, \
            N_cont_statesPer_eV, \
            nameDOS, nameEaxis, \
            betaExtra_List=betaExtra_local, \
            is3D=is3D, \
            branchNormalize=useBranchNormalize)

    else

        nCh = SNS_ComputeDOS_FromChannels( \
            B_T, \
            L_N_List_nm, W_eff_List_nm, wChan, T_eff_List, \
            Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, \
            xV_nm, yV_nm, \
            N_cont_statesPer_eV, \
            nameDOS, nameEaxis, \
            is3D=is3D, \
            branchNormalize=useBranchNormalize)

    endif

    SetDataFolder $savedDF
    return nCh
End
