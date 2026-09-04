#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma ModuleName=SNS_Maps
//==============================================================================
// SNS_Maps.ipf
//
// Assemble DOS(E,B) / LDOS(E,r,B) from channel ensembles and dGSJ eigenbranches,
// including thermal + lock-in modulation broadening.
//
// Depends on:
//   - SNS_Core.ipf
//   - SNS_Solver.ipf
//   - SNS_GeometryFromMask.ipf (for mask-derived channel inputs and vortex/Delta maps)
//==============================================================================






// Global parameters for Gamma(B)
Variable/G gSNS_B0_T        = 0.15     // Tesla
Variable/G gSNS_GammaAlpha1 = 0
Variable/G gSNS_GammaAlpha2 = 1
Variable/G gSNS_GammaMin    = 1e-12

//============================================================
// 7. DOS / MAP builder using all branches and angle sampling
//============================================================

//==============================================================================
// Compute_DOS_SNS_Map_AllBranches
//
// Purpose:
//   Compute derived quantities and assemble output maps.
//
// Inputs:
//   B_T : [T]
//   L : input
//   Wmax : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   T : input
//   Broadening : input
//   lambdaF : input
//   NE : input
//   nameDOS : input
//   nameEaxis : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Compute_DOS_SNS_Map_AllBranches(B_T, L, Wmax, Delta, vF, lambdaL, T, Broadening, lambdaF, NE, nameDOS, nameEaxis)
    Wave B_T
    Variable L, Wmax, Delta, vF, lambdaL, T
    Variable Broadening      // base Lorentzian broadening [eV]
    Variable lambdaF         // Fermi wavelength [m]
    Variable NE              // number of energy points
    String nameDOS, nameEaxis

    // ---- field–dependent Gamma(B) parameters (LOCAL) ----
    // Change these numbers here if you want to tune later.
    Variable SNS_B0_T        = 0.1      // Tesla, scale field for Gamma(B)
    Variable SNS_GammaAlpha1 = 1         // linear term in b = |B|/SNS_B0_T
    Variable SNS_GammaAlpha2 = 0.01         // quadratic term in b
    Variable SNS_GammaMin    = 75e-6     // lower bound for Gamma

    // ---- field–dependent Delta(B) parameters (LOCAL) ----
    // We want a 20% drop at B = 0.5 T: Δ(B) = Δ0 * [1 - 0.2 * (|B|/0.5)^2]
    Variable SNS_DeltaBmax_T   = 0.5     // Tesla: field where drop is 20%
    Variable SNS_DeltaFracDrop = 0.2     // 20% reduction at Bmax

    Variable nB = numpnts(B_T)
    if (nB <= 0)
        Abort "Compute_DOS_SNS_Map_AllBranches: B_T has no points."
    endif

    // ---------------------------
    // 1. Energy axis
    // ---------------------------
    Make/O/D/N=(NE) $nameEaxis
    Wave E_axis = $nameEaxis
    SetScale/I x, -Delta, Delta, E_axis
    E_axis = x

    // DOS(E,B)
    Make/O/D/N=(NE, nB) $nameDOS
    Wave DOS_EB = $nameDOS
    DOS_EB = 0

    SetScale/P x, DimOffset(E_axis, 0), DimDelta(E_axis, 0), "eV", DOS_EB
    SetScale/P y, DimOffset(B_T, 0),    DimDelta(B_T, 0),    "T",  DOS_EB

    // ---------------------------
    // 2. Channel sampling
    // ---------------------------
    Variable Nch = trunc(2*Wmax / lambdaF)
    if (Nch < 1)
        Nch = 1
    endif

    Variable dw = Wmax / Nch
    Make/O/D/N=(Nch) thetaList, wList

    Variable j
    for (j = 0; j < Nch; j += 1)
        wList[j]     = (j + 0.5) * dw
        thetaList[j] = atan2(wList[j], L)  // [rad]
    endfor

    // ---------------------------
    // 3. Solve ABS + build DOS
    // ---------------------------
    Variable iB, k, iE
    Variable theta, W_theta, nBr, E0, E0_scaled, dE
    Variable weight = 1.0 / Nch

    String nameE2D, nameM, nameS

    for (j = 0; j < Nch; j += 1)

        theta   = thetaList[j]
        W_theta = L * tan(theta)

        sprintf nameE2D, "E_allBranches_th%03d", j
        sprintf nameM,   "m_allBranches_th%03d", j
        sprintf nameS,   "s_allBranches_th%03d", j

        nBr = Solve_AllBranches_SNS_dGSJ(B_T, L, W_theta, Delta, vF, lambdaL, T, \
                                         nameE2D, nameM, nameS)
        if (nBr <= 0)
            continue
        endif

        Wave E_all = $nameE2D

        for (iB = 0; iB < nB; iB += 1)

            Variable Bval   = B_T[iB]
            Variable GammaB, bDimGamma, gammaScale
            Variable bDimDelta, deltaScale

            // -------------------------------
            // Inline field-dependent Gamma(B)
            // -------------------------------
            if ((SNS_GammaAlpha1 == 0) && (SNS_GammaAlpha2 == 0) || (SNS_B0_T <= 0))
                GammaB = Broadening
            else
                bDimGamma = abs(Bval)/SNS_B0_T
                gammaScale = 1 + SNS_GammaAlpha1*bDimGamma + SNS_GammaAlpha2*bDimGamma*bDimGamma
                GammaB = Broadening * gammaScale
                if (GammaB < SNS_GammaMin)
                    GammaB = SNS_GammaMin
                endif
            endif

            // -------------------------------
            // Field-dependent Delta(B) scaling
            // -------------------------------
            // deltaScale(B) = Δ(B)/Δ0 = 1 - SNS_DeltaFracDrop * (|B|/SNS_DeltaBmax_T)^2
            // clamp so it never goes negative
            if (SNS_DeltaBmax_T > 0)
                bDimDelta = abs(Bval)/SNS_DeltaBmax_T
            else
                bDimDelta = 0
            endif
            deltaScale = 1 - SNS_DeltaFracDrop * bDimDelta*bDimDelta
            if (deltaScale < 0)
                deltaScale = 0
            endif

            for (k = 0; k < nBr; k += 1)

                E0 = E_all[iB][k]
                if (numtype(E0) != 0)
                    continue
                endif

                // apply Delta(B) scaling to the ABS energy
                E0_scaled = E0 * deltaScale

                for (iE = 0; iE < NE; iE += 1)
                    dE = E_axis[iE] - E0_scaled
                    DOS_EB[iE][iB] += weight * (GammaB/pi) / (dE*dE + GammaB*GammaB)
                endfor

            endfor
        endfor
    endfor

    return Nch
End


//============================================================
// Robust full-B solver: scan all branches at every field
//============================================================
//
// Signature:
//   Solve_AllBranches_SNS_dGSJ_FullScan(B_T, L, W, Delta, vF, lambdaL, T, mMax, \
//                                       nameE2D, nameM, nameS [, tol, maxIters])
//
// - B_T      : 1D wave of magnetic field values [T]
// - L, W     : geometry [m]
// - Delta    : gap [eV]
// - vF       : Fermi velocity [m/s]
// - lambdaL  : London depth [m]
// - T        : transparency (0..1)
// - mMax     : max |m| to scan (integer)
// - nameE2D  : name of output 2D wave E2D[iB][jBranch]
// - nameM    : name of 1D wave storing m per branch index j
// - nameS    : name of 1D wave storing sSign (=±1) per branch index j
//
// Optional:
//   tol      : root-finder tolerance; default 1e-9
//   maxIters : max iterations; default 200
//
// Behaviour:
//   - For EVERY (B, m, sSign) it calls SolveBranchAtB_dGSJ.
//   - If there is a valid subgap root, it writes E2D[iB][j] = E;
//     otherwise E2D[iB][j] = NaN.
//   - j runs over all (m, sSign) in a fixed, deterministic order.
//============================================================

//==============================================================================
// Compute_DOS_SNS_Ensemble
//
// Purpose:
//   Compute derived quantities and assemble output maps.
//
// Inputs:
//   B_T : [T]
//   L : input
//   Wmax : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   T : input
//   Broadening : input
//   lambdaF : input
//   NE : input
//   nameDOS_ensemble : input
//   nameEaxis : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Compute_DOS_SNS_Ensemble(B_T, L, Wmax, Delta, vF, lambdaL, T, Broadening, lambdaF, NE, nameDOS_ensemble, nameEaxis)
    Wave B_T
    Variable L, Wmax, Delta, vF, lambdaL, T
    Variable Broadening, lambdaF, NE
    String nameDOS_ensemble, nameEaxis

    // Define a few geometries and weights
    // Example: three Wmax values: nominal, -10%, +10%
    Variable Ngeom = 11
    Make/O/D/N=(Ngeom) L_geom, W_geom, p_geom

    L_geom[0] = L
    W_geom[0] = Wmax * 0.9
    p_geom[0] = 0.3

    L_geom[1] = L
    W_geom[1] = Wmax
    p_geom[1] = 0.4

    L_geom[2] = L
    W_geom[2] = Wmax * 1.1
    p_geom[2] = 0.3

    // First geometry: compute DOS and axis
    String nameDOSg, nameEaxisLocal
    sprintf nameDOSg, "DOS_SNS_geom%02d", 0
    nameEaxisLocal = nameEaxis

    Variable Nch0 = Compute_DOS_SNS_Map_AllBranches(B_T, L_geom[0], W_geom[0], Delta, vF, lambdaL, T, \
                                                    Broadening, lambdaF, NE, nameDOSg, nameEaxisLocal)
    Wave DOSg = $nameDOSg
    Wave E_axis = $nameEaxisLocal

    // Allocate ensemble DOS and initialize with first geometry
    Make/O/D/N=(DimSize(DOSg,0), DimSize(DOSg,1)) $nameDOS_ensemble
    Wave DOS_ensemble = $nameDOS_ensemble
    DOS_ensemble = p_geom[0] * DOSg

    Variable ig
    for (ig = 1; ig < Ngeom; ig += 1)

        sprintf nameDOSg, "DOS_SNS_geom%02d", ig
        Compute_DOS_SNS_Map_AllBranches(B_T, L_geom[ig], W_geom[ig], Delta, vF, lambdaL, T, \
                                        Broadening, lambdaF, NE, nameDOSg, nameEaxisLocal)
        Wave DOSg_i = $nameDOSg

        DOS_ensemble += p_geom[ig] * DOSg_i
    endfor

    // Set axis scales on DOS_ensemble to match E_axis and B_T
    SetScale/P x, DimOffset(E_axis, 0), DimDelta(E_axis, 0), "eV", DOS_ensemble
    SetScale/P y, DimOffset(B_T,    0), DimDelta(B_T,    0), "T",  DOS_ensemble

    return Ngeom
End

//============================================================
// Channel-based DOS(E,B) with per-channel transparency
//
// B_T        : 1D wave of B fields [T]
// L_N_List   : chord lengths per channel [m]
// W_eff_List : effective magnetic widths [m]
// wChan      : geometric weights (any scale; normalized inside)
// T_eff_List : per-channel transparencies (0..1), only used as T in solver
//
// Delta      : gap [eV]
// vF         : Fermi velocity [m/s]
// lambdaL    : London penetration depth [m]
// Broadening : base Lorentzian broadening [eV]
// NE         : number of energy points
// nameDOS    : output DOS(E,B)
// nameEaxis  : output energy axis
//============================================================

//==============================================================================
// Compute_DOS_SNS_Map_FromChannels
//
// Purpose:
//   Compute derived quantities and assemble output maps.
//
// Inputs:
//   B_T : [T]
//   L_N_List : input
//   W_eff_List : input
//   wChan : input
//   T_eff_List : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   Broadening : input
//   NE : input
//   nameDOS : input
//   nameEaxis : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Compute_DOS_SNS_Map_FromChannels(B_T, L_N_List, W_eff_List, wChan, T_eff_List,Delta, vF, lambdaL, Broadening, NE, nameDOS, nameEaxis)

    Wave    B_T, L_N_List, W_eff_List, wChan, T_eff_List
    Variable Delta, vF, lambdaL
    Variable Broadening, NE
    String  nameDOS, nameEaxis

    // ---- Gamma(B) parameters ----
    Variable SNS_B0_T        = 0.1
    Variable SNS_GammaAlpha1 = 1
    Variable SNS_GammaAlpha2 = 0.01
    Variable SNS_GammaMin    = Broadening

    // ---- Delta(B) parameters ----
    Variable SNS_DeltaBmax_T   = 0.5
    Variable SNS_DeltaFracDrop = 0.2

    Variable nB  = numpnts(B_T)
    Variable Nch = numpnts(L_N_List)

    if ((Nch <= 0) || (Nch != numpnts(W_eff_List)) || (Nch != numpnts(wChan)) || (Nch != numpnts(T_eff_List)))
        Abort "Compute_DOS_SNS_Map_FromChannels: channel waves inconsistent."
    endif

    // ---- normalize geometric weights ONLY ----
    Variable sumW = sum(wChan)
    if (sumW <= 0)
        Abort "Compute_DOS_SNS_Map_FromChannels: all channel weights are zero."
    endif

    // ---------------------------
    // Energy axis
    // ---------------------------
    Make/O/D/N=(NE) $nameEaxis
    Wave E_axis = $nameEaxis
    SetScale/I x, -Delta, Delta, E_axis
    E_axis = x

    // DOS(E,B)
    Make/O/D/N=(NE, nB) $nameDOS
    Wave DOS_EB = $nameDOS
    DOS_EB = 0

    SetScale/P x, DimOffset(E_axis,0), DimDelta(E_axis,0), "eV", DOS_EB
    SetScale/P y, DimOffset(B_T,0),    DimDelta(B_T,0),    "T",  DOS_EB

    // ---------------------------
    // Precompute Gamma(B), Delta(B) scaling
    // ---------------------------
    Make/FREE/D/N=(nB) GammaB_wv, deltaScale_wv
    Variable iB, Bval, bDimGamma, gammaScale, bDimDelta, deltaScale

    for (iB = 0; iB < nB; iB += 1)
        Bval = B_T[iB]

        // Gamma(B)
        if (SNS_B0_T > 0)
            bDimGamma  = abs(Bval)/SNS_B0_T
            gammaScale = 1 + SNS_GammaAlpha1*bDimGamma + SNS_GammaAlpha2*bDimGamma*bDimGamma
            GammaB_wv[iB] = max(Broadening * gammaScale, SNS_GammaMin)
        else
            GammaB_wv[iB] = Broadening
        endif

        // Delta(B)
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

    // ---------------------------
    // Channel loop
    // ---------------------------
    Variable j, k, iE
    Variable Lch, Wch, Teff, weight
    Variable nBr, E0, E0_scaled, dE, GammaB_loc, deltaScale_loc

    String nameE2D, nameM, nameS

    for (j = 0; j < Nch; j += 1)

        Lch  = L_N_List[j]
        Wch  = W_eff_List[j]
        Teff = T_eff_List[j]
        weight = wChan[j]/sumW    // geometric only

        if ((weight <= 0) || numtype(Lch) || numtype(Wch) || numtype(Teff))
            continue
        endif

        sprintf nameE2D, "E_allBranches_ch%03d", j
        sprintf nameM,   "m_allBranches_ch%03d", j
        sprintf nameS,   "s_allBranches_ch%03d", j

        // per-channel transparency ONLY here
        nBr = Solve_AllBranches_SNS_dGSJ(B_T, Lch, Wch, Delta, vF, lambdaL, Teff, \
                                         nameE2D, nameM, nameS)
        if (nBr <= 0)
            continue
        endif

        Wave E_all = $nameE2D

        for (iB = 0; iB < nB; iB += 1)

            GammaB_loc    = GammaB_wv[iB]
            deltaScale_loc = deltaScale_wv[iB]

            for (k = 0; k < nBr; k += 1)

                E0 = E_all[iB][k]
                if (numtype(E0))
                    continue
                endif

                E0_scaled = E0 * deltaScale_loc * Delta

                for (iE = 0; iE < NE; iE += 1)
                    dE = E_axis[iE] - E0_scaled
                    DOS_EB[iE][iB] += weight * (GammaB_loc/pi) / (dE*dE + GammaB_loc*GammaB_loc)
                endfor

            endfor
        endfor
    endfor

    return Nch
End


//============================================================
// DOS(E,B) from channels with per-channel transparency +
// coherence-length–dependent broadening.
//
// B_T        : 1D wave of B fields [T]
// L_N_List   : chord lengths per channel [m]
// W_eff_List : effective magnetic widths [m]
// wChan      : geometric weights (any scale; normalized inside)
// T_eff_List : per-channel transparencies (0..1)
//
// Delta      : gap [eV]
// vF         : Fermi velocity [m/s]
// lambdaL    : London penetration depth [m]
// Broadening : base Lorentzian broadening at B=0 [eV]
// NE         : number of energy points
// nameDOS    : output DOS(E,B)
// nameEaxis  : output energy axis
//
// Broadening model:
//   xi = ħ vF / (π Δ)
//   Gamma_coh(L) = (Δ/π) * (1 - exp(-L/xi))
//   Gamma_tot(B,L) = sqrt( Gamma_B(B)^2 + Gamma_coh(L)^2 )
//============================================================

//==============================================================================
// Compute_DOS_SNS_Map_FromChannels_Coh
//
// Purpose:
//   Compute derived quantities and assemble output maps.
//
// Inputs:
//   B_T : [T]
//   L_N_List : input
//   W_eff_List : input
//   wChan : input
//   T_eff_List : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   Broadening : input
//   NE : input
//   nameDOS : input
//   nameEaxis : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Compute_DOS_SNS_Map_FromChannels_Coh(B_T, L_N_List, W_eff_List, wChan, T_eff_List,Delta, vF, lambdaL, Broadening, NE, nameDOS, nameEaxis)

    Wave    B_T, L_N_List, W_eff_List, wChan, T_eff_List
    Variable Delta, vF, lambdaL
    Variable Broadening, NE
    String  nameDOS, nameEaxis

    // ---- Gamma(B) parameters ----
    Variable SNS_B0_T        = 0
    Variable SNS_GammaAlpha1 = 0
    Variable SNS_GammaAlpha2 = 0
    Variable SNS_GammaMin    = Broadening

    // ---- Delta(B) parameters ----
    Variable SNS_DeltaBmax_T   = 0.5
    Variable SNS_DeltaFracDrop = 0.001

    Variable nB  = numpnts(B_T)
    Variable Nch = numpnts(L_N_List)

    if ((Nch <= 0) || (Nch != numpnts(W_eff_List)) || (Nch != numpnts(wChan)) || (Nch != numpnts(T_eff_List)))
        Abort "Compute_DOS_SNS_Map_FromChannels_Coh: channel waves inconsistent."
    endif

    if (nB <= 0)
        Abort "Compute_DOS_SNS_Map_FromChannels_Coh: B_T has no points."
    endif

    // ---- normalize geometric weights only ----
    Variable sumW = sum(wChan)
    if (sumW <= 0)
        Abort "Compute_DOS_SNS_Map_FromChannels_Coh: all channel weights are zero."
    endif

    //================================================
    // 1. Coherence length and Gamma_coh(L_j)
    //================================================
    // xi = ħ vF / (π Δ)
    Variable hbar_eVs = 6.582119569e-16    // eV·s
    Variable xi = hbar_eVs * vF / (Pi * Delta)   // [m]

    if (xi <= 0)
        Abort "Compute_DOS_SNS_Map_FromChannels_Coh: invalid coherence length."
    endif

	Make/FREE/D/N=(Nch) GammaCohCh
	Variable j, Lj, xLoverXi

	for (j = 0; j < Nch; j += 1)
    	Lj = L_N_List[j]          // [m]
    	if (numtype(Lj))
        	GammaCohCh[j] = 0
    	else
        	xLoverXi = Lj/(xi)
        	if (xLoverXi < 0)
            	xLoverXi = 0
        	endif
        	// Gamma_coh(L) = (Delta/pi) * (1 - exp(-L/xi))
        	GammaCohCh[j] = (Delta/Pi) * (1 - exp(-xLoverXi))
    	endif
	endfor


        // ---------------------------
    // Precompute Gamma(B), Delta(B) scaling
    // ---------------------------
    Make/FREE/D/N=(nB) GammaB_wv, deltaScale_wv
    Variable iB, Bval, bDimGamma, gammaScale, bDimDelta, deltaScale

    for (iB = 0; iB < nB; iB += 1)
        Bval = B_T[iB]

        // Gamma(B)
        if (SNS_B0_T > 0)
            bDimGamma  = abs(Bval)/SNS_B0_T
            gammaScale = 1 + SNS_GammaAlpha1*bDimGamma + SNS_GammaAlpha2*bDimGamma*bDimGamma
            GammaB_wv[iB] = max(Broadening * gammaScale, SNS_GammaMin)
        else
            GammaB_wv[iB] = Broadening
        endif

        // Delta(B)
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
    Make/O/D/N=(NE) $nameEaxis
    Wave E_axis = $nameEaxis
    SetScale/I x, -Delta, Delta, "eV", E_axis
    E_axis = x

    Make/O/D/N=(NE, nB) $nameDOS
    Wave DOS_EB = $nameDOS
    DOS_EB = 0
    
    string nameDOS_broadening=nameDOS
    nameDOS_broadening+= "_Gamma"
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
    Variable k, iE
    Variable Lch, Wch, Teff, weight
    Variable nBr, E0, E0_scaled, dE
    Variable GammaB_loc, GammaC_loc, GammaTot, deltaScale_loc

    String nameE2D, nameM, nameS

    for (j = 0; j < Nch; j += 1)

        Lch  = L_N_List[j]
        Wch  = W_eff_List[j]
        Teff = T_eff_List[j]
        weight = wChan[j]/sumW

        if ((weight <= 0) || numtype(Lch) || numtype(Wch) || numtype(Teff))
            continue
        endif

        sprintf nameE2D, "E_allBranches_ch%03d", j
        sprintf nameM,   "m_allBranches_ch%03d", j
        sprintf nameS,   "s_allBranches_ch%03d", j

        // per-channel transparency in the solver
        nBr = Solve_AllBranches_SNS_dGSJ(B_T, Lch, Wch, Delta, vF, lambdaL, Teff, \
                                         nameE2D, nameM, nameS)
        if (nBr <= 0)
            continue
        endif

        Wave E_all = $nameE2D
        GammaC_loc = GammaCohCh[j]       // fixed per channel

	for (iB = 0; iB < nB; iB += 1)

    		Bval = B_T[iB]
    		GammaB_loc     = GammaB_wv[iB]
    		deltaScale_loc = deltaScale_wv[iB]

    		// --- field-dependent factor for coherence broadening ---
    		// use the same scale as Delta(B) suppression:
    		// SNS_DeltaBmax_T is your "pair-breaking" field scale
    		Variable bCo = 0
    		if (SNS_DeltaBmax_T > 0)
        		bCo = abs(Bval)/SNS_DeltaBmax_T
        		if (bCo > 1)
            		bCo = 1    // saturate for large B
        		endif
    		endif
    		// coherence broadening vanishes at B=0, reaches full GammaCohCh at B ~ SNS_DeltaBmax_T
    		GammaC_loc = GammaCohCh[j] * (bCo*bCo)    // ∝ B^2, as for orbital pair-breaking

    		// total broadening
    		GammaTot = sqrt(GammaB_loc*GammaB_loc + GammaC_loc*GammaC_loc)

            
            for (k = 0; k < nBr; k += 1)

                E0 = E_all[iB][k]
                if (numtype(E0))
                    continue
                endif

                E0_scaled = E0 * deltaScale_loc

                for (iE = 0; iE < NE; iE += 1)
                    dE = E_axis[iE] - E0_scaled
                    DOS_EB[iE][iB] += weight * (Broadening/Pi) / (dE*dE + Broadening*Broadening)
                    DOS_EB_broadening[iE][iB] += weight * (GammaTot/Pi) / (dE*dE + GammaTot*GammaTot)
                endfor

            endfor
        endfor
    endfor
    
    //================================================
    // 5. Cleanup: move branch waves into subfolders
    //================================================
    NewDataFolder/O E_allBranches
    NewDataFolder/O m_allBranches
    NewDataFolder/O s_allBranches

    // E-branches: "E_allBranches_chXXX" -> :E_allBranches:chXXX
    MoveBranchWavesToSubfolder("E_allBranches_", "E_allBranches")

    // m-branches: "m_allBranches_chXXX" -> :m_allBranches:chXXX
    MoveBranchWavesToSubfolder("m_allBranches_", "m_allBranches")

    // s-branches: "s_allBranches_chXXX" -> :s_allBranches:chXXX
    MoveBranchWavesToSubfolder("s_allBranches_", "s_allBranches")   

    return Nch
End



// Wrapper using SNS_Params instead of long parameter list

//==============================================================================
// Compute_DOS_SNS_Map_FromChannels_Coh_Params
//
// Purpose:
//   Compute derived quantities and assemble output maps.
//
// Inputs:
//   B_T : [T]
//   L_N_List : input
//   W_eff_List : input
//   wChan : input
//   T_eff_List : input
//   p : input
//   dosName : input
//   eAxisName : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Compute_DOS_SNS_Map_FromChannels_Coh_Params(B_T, L_N_List, W_eff_List, wChan, T_eff_List, p, dosName, eAxisName)
    Wave   B_T, L_N_List, W_eff_List, wChan, T_eff_List
    STRUCT SNS_Params &p
    String dosName, eAxisName

    return Compute_DOS_SNS_Map_FromChannels_Coh(B_T, L_N_List, W_eff_List, wChan, T_eff_List,p.Delta, p.vF, p.lambdaL, p.Broadening, p.NE,dosName, eAxisName)
End



//==============================================================================
// SNS_ComputeDOS_FromSettings_TI
//
// Minimal TI wrapper around the existing DOS machinery.
//
// What it does:
//   1. Calls SNS_InitDefaultSettings with TI Dirac mode enabled.
//   2. Sets hbar*vD, mu relative to Dirac point, and Dirac barrier ZD.
//   3. Rebuilds ray channels so T_eff_List uses Dirac/Klein transparency.
//   4. Calls the existing SNS_ComputeDOS_FromSettings(...).
//
// Notes:
//   - ZD is stored in root:SNS_Settings:BTK_barrier for compatibility.
//     In TI mode this is interpreted as the Dirac delta-barrier phase
//     Z_D = U0*d/(hbar*vD).
//   - hbarvD_eVA is hbar*vD in eV Angstrom.
//   - muDirac_eV is chemical potential relative to the Dirac point [eV].
//   - This wrapper leaves SNS_bandModel = 1 after execution. Normal runs should
//     call SNS_InitDefaultSettings(..., bandModel_in=0), or simply
//     SNS_InitDefaultSettings(...) if default bandModel_in behavior is 0.
//==============================================================================
Function SNS_ComputeDOS_FromSettings_TI(dataFolder, nameDOS, nameEaxis, STSx, STSy, Bangle_deg, [hbarvD_eVA, muDirac_eV, ZD, stepFac])
    String dataFolder
    String nameDOS, nameEaxis
    Variable STSx, STSy
    Variable Bangle_deg
    Variable hbarvD_eVA, muDirac_eV, ZD, stepFac

    String savedDF = GetDataFolder(1)

    if (ParamIsDefault(hbarvD_eVA))
        hbarvD_eVA = 3.0
    endif
    if (ParamIsDefault(muDirac_eV))
        muDirac_eV = 0.2
    endif
    if (ParamIsDefault(ZD))
        ZD = 0.3
    endif
    if (ParamIsDefault(stepFac))
        stepFac = 0.5
    endif

    // ------------------------------
    // Initialize TI settings
    // ------------------------------
    SNS_InitDefaultSettings( \
        bandModel_in=1, \
        hbarvD_eVA_in=hbarvD_eVA, \
        muDirac_eV_in=muDirac_eV, \
        BTKbarrier_in=ZD)

    // Optional diagnostic bindings
    NVAR/Z LambdaF = root:SNS_Settings:LambdaF
    NVAR/Z vF      = root:SNS_Settings:vF

    // ------------------------------
    // Rebuild channels using TI transparency
    // ------------------------------
    DFREF dfrGeom = $dataFolder

    Wave/Z w_mask = dfrGeom:w_mask
    if (!WaveExists(w_mask))
        SetDataFolder $savedDF
        Abort "SNS_ComputeDOS_FromSettings_TI: missing w_mask in dataFolder."
    endif

    SetDataFolder dfrGeom

    // phiB is expected in radians by SNS_BuildChannelsFromMask2D.
    Variable phiB = Bangle_deg*pi/180

    SNS_BuildChannelsFromMask2D(w_mask, STSx, STSy, phiB, stepFac, "")

    // Ensure mask area/perimeter helper waves exist for DOS normalization.
    Wave/Z w_area_nm2 = SNS_MaskAreaPerim_FromParticles(w_mask)

    SetDataFolder $savedDF

    // ------------------------------
    // Existing DOS machinery
    // ------------------------------
    Variable nCh
    nCh = SNS_ComputeDOS_FromSettings(dataFolder, nameDOS, nameEaxis)

    Print "SNS TI DOS: folder = ", dataFolder
    Print "SNS TI DOS: nCh = ", nCh
    Print "SNS TI DOS: hbarvD_eVA = ", hbarvD_eVA
    Print "SNS TI DOS: muDirac_eV = ", muDirac_eV
    Print "SNS TI DOS: ZD = ", ZD

    if (NVAR_Exists(LambdaF))
        Print "SNS TI DOS: lambdaF [nm] = ", LambdaF*1e9
    endif
    if (NVAR_Exists(vF))
        Print "SNS TI DOS: vD [m/s] = ", vF
    endif

    SetDataFolder $savedDF
    return nCh
End

//============================================================
// SNS_VortexPhaseOnChord
//
// Line-integral of A(x,y) along a straight SNS chord, converted
// to the total orbital phase β_vortex seen by the electron+hole
// loop for that channel.
//
// Ax,Ay : 2D waves from SNS_MakeVortexA (same scaling as Nmask)
// x1,y1 : S-contact 1 position [nm]
// x2,y2 : S-contact 2 position [nm]
// nStep : number of integration steps along chord
// betaOut : returned β_vortex [radians]
//
// NOTES:
// - Uses wave scaling (DimOffset/DimDelta) to map (x,y)→(ix,iy).
// - Assumes Ax/Ay as constructed in SNS_MakeVortexA (Φ_v/(2πr) with r in nm).
//   With that convention, the sum(A·dl_nm) is numerically equal to the
//   physical flux in Wb, so β = (2e/ħ) * ∫A·dl is consistent.
//============================================================
//============================================================
// SNS_VortexPhaseOnChord  (REPLACEMENT)
//
// Line-integral of A(x,y) along a straight SNS chord, converted
// to β_vortex = (2e/ħ) ∫ A·dl.
//
// IMPORTANT FIXES vs old version:
//  1) Consistent axis mapping: dim0 = x, dim1 = y (matches SNS_MakeVortexA: Ax[ix][iy]).
//  2) Proper SI units: SNS_MakeVortexA used r in nm → Ax/Ay are Wb/nm.
//     Convert Ax/Ay to Wb/m (×1e9) and dl to meters (×1e-9) so ∫A·dl is Wb.
//============================================================

//==============================================================================
// SNS_VortexPhaseOnChord
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   Ax : input
//   Ay : input
//   x1 : input
//   y1 : input
//   x2 : input
//   y2 : input
//   nStep : input
//   betaOut : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_VortexPhaseOnChord(Ax, Ay, x1, y1, x2, y2, nStep, betaOut)
    Wave Ax, Ay
    Variable x1, y1, x2, y2
    Variable nStep
    Variable &betaOut

    if (nStep < 1)
        betaOut = 0
        return -1
    endif

    // segment (nm)
    Variable dxSeg_nm = (x2 - x1)/nStep
    Variable dySeg_nm = (y2 - y1)/nStep

    // segment (m)
    Variable dxSeg_m  = dxSeg_nm * 1e-9
    Variable dySeg_m  = dySeg_nm * 1e-9

    // scaling of A-waves (nm): dim0 = x, dim1 = y
    Variable x0_nm = DimOffset(Ax, 0)
    Variable y0_nm = DimOffset(Ax, 1)
    Variable dx_nm = DimDelta(Ax, 0)
    Variable dy_nm = DimDelta(Ax, 1)

    Variable nX = DimSize(Ax, 0)   // x index range
    Variable nY = DimSize(Ax, 1)   // y index range

    Variable k, t, x_nm, y_nm
    Variable ix, iy
    Variable phaseFlux_Wb = 0      // ∫A·dl in Wb

    for (k = 0; k < nStep; k += 1)
        t    = (k + 0.5)/nStep
        x_nm = x1 + t*(x2 - x1)
        y_nm = y1 + t*(y2 - y1)

        ix = round((x_nm - x0_nm)/dx_nm)
        iy = round((y_nm - y0_nm)/dy_nm)

        if (ix < 0 || ix >= nX || iy < 0 || iy >= nY)
            continue
        endif

        // Ax/Ay are Wb/nm (because r was in nm). Convert to Wb/m via ×1e9.
        phaseFlux_Wb += (Ax[ix][iy]*1e9)*dxSeg_m + (Ay[ix][iy]*1e9)*dySeg_m
    endfor

    // Convert ∫A·dl → β_vortex = (2e/ħ) * ∫A·dl
    Variable q_e  = 1.602176634e-19   // C
    Variable hbar = 1.054571817e-34   // J·s

    betaOut = (2*q_e/hbar) * phaseFlux_Wb

    return 0
End


//============================================================
// SNS_EffectiveDeltaOnChord
//
// Compute average Delta along a straight chord from (x1,y1)
// to (x2,y2) in nm on the DeltaMap grid.
//
// DeltaMap : 2D wave [eV], same scaling as Nmask
// x1,y1,x2,y2 : endpoints [nm]
// nStep   : number of sample points along chord
//============================================================

//==============================================================================
// SNS_EffectiveDeltaOnChord
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   DeltaMap : [eV]
//   x1 : input
//   y1 : input
//   x2 : input
//   y2 : input
//   nStep : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_EffectiveDeltaOnChord(DeltaMap, x1, y1, x2, y2, nStep)
    Wave DeltaMap
    Variable x1, y1, x2, y2
    Variable nStep

    if (nStep < 1)
        return NaN
    endif

    if (WaveDims(DeltaMap) != 2)
        Abort "SNS_EffectiveDeltaOnChord: DeltaMap must be 2D."
    endif

    Variable nX = DimSize(DeltaMap, 0)
    Variable nY = DimSize(DeltaMap, 1)

    Variable x0 = DimOffset(DeltaMap,0)
    Variable dx = DimDelta(DeltaMap,0)
    Variable y0 = DimOffset(DeltaMap,1)
    Variable dy = DimDelta(DeltaMap,1)

    if (dx == 0 || dy == 0)
        Abort "SNS_EffectiveDeltaOnChord: invalid scaling."
    endif

    Variable k, t, x, y
    Variable ix, iy
    Variable sumDelta = 0
    Variable nValid   = 0
    Variable DeltaLoc

    for (k = 0; k < nStep; k += 1)
        t = (k + 0.5)/nStep
        x = x1 + t*(x2 - x1)
        y = y1 + t*(y2 - y1)

        ix = round((x - x0)/dx)
        iy = round((y - y0)/dy)

        if (ix < 0 || ix >= nX || iy < 0 || iy >= nY)
            continue
        endif

        DeltaLoc = DeltaMap[ix][iy]
        if (numtype(DeltaLoc) != 0)
            continue
        endif

        sumDelta += DeltaLoc
        nValid   += 1
    endfor

    if (nValid <= 0)
        return NaN
    endif

    return sumDelta/nValid
End



//============================================================
// Build lock-in modulation kernel (sinusoidal modulation).
// E_axis  : 1D energy axis [eV], symmetric about 0
// Eac_eV  : modulation amplitude in energy [eV] (peak)
//           e.g. 50 µeV -> 50e-6
// nameKM  : name of output kernel wave
//============================================================



// Return integer half-width (in points) of a symmetric kernel,
// defined as max distance from center where |K| > eps.





//============================================================
// Apply thermal + modulation broadening to DOS(E,B) using
// 1D convolutions along the energy axis, with mirror padding
// to avoid edge artifacts.
//
// Uses SNS settings:
//   root:SNS_Settings:T_K    (temperature in K)
//   root:SNS_Settings:V_mod  (modulation amplitude in eV)
// Energy axis is taken from x-scaling of DOS_EB_in.
//============================================================


//==============================================================================
// SNS_PointLDOS_ZeroField_ResolutionMatrix_FromSettings
//
// Zero-field local LDOS resolution sweep for one tip position.
//
// Channel-source priority:
//   1) dataFolder contains full solver-ready channel waves.
//   2) dataFolder:RayTraceHist contains full solver-ready channel waves.
//   3) dataFolder:RayTraceHist_3D contains full tagged solver-ready channel waves.
//   4) otherwise run SNS_ExtractModesForFolder(...).
//
// The selected channels are copied into a temporary internal solver folder:
//   dataFolder:_ZeroFieldRes_WorkChannels
//
// This avoids modifying archived RayTraceHist folders.
//
// Output folder:
//   dataFolder:outFolderName
//
// Output waves:
//   <outPrefix>_LDOS_ERes
//   <outPrefix>_E_eV
//   <outPrefix>_T_K
//   <outPrefix>_Lphi_um
//   <outPrefix>_Gamma_eV
//   <outPrefix>_Gamma_ueV
//   <outPrefix>_label
//
// Defaults:
//   T_delta = 0.01 K
//   phase_coherence_delta = 1 um
//   outFolderName = "ZeroField_ResolutionMatrix"
//   outPrefix = "ZeroFieldRes"
//==============================================================================

Function SNS_PointLDOS_ZeroField_ResolutionMatrix_FromSettings(dataFolder, [T_delta, phase_coherence_delta, outFolderName, outPrefix, Tmin_K, Tmax_K, LphiMin_um, LphiMax_um, Bangle_deg, BTK_barrier, STSx, STSy, forceRayTrace])
	String dataFolder
	Variable T_delta, phase_coherence_delta
	String outFolderName, outPrefix
	Variable Tmin_K, Tmax_K, LphiMin_um, LphiMax_um
	Variable Bangle_deg, BTK_barrier, STSx, STSy, forceRayTrace

	if (ParamIsDefault(T_delta))
		T_delta = 0.01
	endif
	if (ParamIsDefault(phase_coherence_delta))
		phase_coherence_delta = 1
	endif
	if (ParamIsDefault(outFolderName))
		outFolderName = "ZeroField_ResolutionMatrix"
	endif
	if (ParamIsDefault(outPrefix))
		outPrefix = "ZeroFieldRes"
	endif
	if (ParamIsDefault(Tmin_K))
		Tmin_K = 0.01
	endif
	if (ParamIsDefault(Tmax_K))
		Tmax_K = 4.2
	endif
	if (ParamIsDefault(LphiMin_um))
		LphiMin_um = 1
	endif
	if (ParamIsDefault(LphiMax_um))
		LphiMax_um = 100
	endif
	if (ParamIsDefault(Bangle_deg))
		Bangle_deg = 225
	endif
	if (ParamIsDefault(forceRayTrace))
		forceRayTrace = 0
	endif

	if (T_delta <= 0)
		Abort "SNS_PointLDOS_ZeroField_ResolutionMatrix_FromSettings: T_delta must be > 0."
	endif
	if (phase_coherence_delta <= 0)
		Abort "SNS_PointLDOS_ZeroField_ResolutionMatrix_FromSettings: phase_coherence_delta must be > 0."
	endif
	if (Tmax_K < Tmin_K)
		Abort "SNS_PointLDOS_ZeroField_ResolutionMatrix_FromSettings: Tmax_K must be >= Tmin_K."
	endif
	if (LphiMax_um < LphiMin_um)
		Abort "SNS_PointLDOS_ZeroField_ResolutionMatrix_FromSettings: LphiMax_um must be >= LphiMin_um."
	endif

	String savedDF = GetDataFolder(1)

	// -------------------------------------------------------------------------
	// Ensure settings exist.
	// -------------------------------------------------------------------------
	if (!DataFolderExists("root:SNS_Settings:"))
		SNS_InitDefaultSettings()
	endif
	if (Exists("root:SNS_Settings:vF") != 2)
		SNS_InitDefaultSettings()
	endif
	if (Exists("root:SNS_Settings:T_K") != 2)
		SNS_InitDefaultSettings()
	endif
	if (Exists("root:SNS_Settings:V_mod") != 2)
		SNS_InitDefaultSettings()
	endif
	if (Exists("root:SNS_Settings:Broadening") != 2)
		SNS_InitDefaultSettings()
	endif
	if (Exists("root:SNS_Settings:SNS_GammaBase_eV") != 2)
		SNS_InitDefaultSettings()
	endif
	if (Exists("root:SNS_Settings:B_T") != 1)
		SNS_InitDefaultSettings()
	endif

	// Compatibility with older gamma code.
	Variable/G root:SNS_Settings:SNS_GammaDoppScale = 0

	// -------------------------------------------------------------------------
	// Check mask.
	// -------------------------------------------------------------------------
	if (Exists(dataFolder + ":w_mask") != 1)
		SetDataFolder $savedDF
		Abort "SNS_PointLDOS_ZeroField_ResolutionMatrix_FromSettings: missing w_mask in dataFolder."
	endif

	// -------------------------------------------------------------------------
	// Select channel source.
	// -------------------------------------------------------------------------
	String sourceDF = ""
	Variable sourceIs3DTagged = 0

	String hist2D = dataFolder + ":RayTraceHist"
	String hist3D = dataFolder + ":RayTraceHist_3D"

	if (!forceRayTrace && SNS__ZeroFieldRes_HasGenericSolverChannels(dataFolder))
		sourceDF = dataFolder
		sourceIs3DTagged = 0
	elseif (!forceRayTrace && DataFolderExists(hist2D + ":") && SNS__ZeroFieldRes_HasGenericSolverChannels(hist2D))
		sourceDF = hist2D
		sourceIs3DTagged = 0
	elseif (!forceRayTrace && DataFolderExists(hist3D + ":") && SNS__ZeroFieldRes_HasTagged3DSolverChannels(hist3D))
		sourceDF = hist3D
		sourceIs3DTagged = 1
	else
		if (ParamIsDefault(BTK_barrier) && ParamIsDefault(STSx) && ParamIsDefault(STSy))
			SNS_ExtractModesForFolder(dataFolder, Bangle_deg=Bangle_deg, doDisplay=0)
		elseif (ParamIsDefault(BTK_barrier) && !ParamIsDefault(STSx) && !ParamIsDefault(STSy))
			SNS_ExtractModesForFolder(dataFolder, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, doDisplay=0)
		elseif (!ParamIsDefault(BTK_barrier) && ParamIsDefault(STSx) && ParamIsDefault(STSy))
			SNS_ExtractModesForFolder(dataFolder, Bangle_deg=Bangle_deg, doDisplay=0)
		elseif (!ParamIsDefault(BTK_barrier) && !ParamIsDefault(STSx) && !ParamIsDefault(STSy))
			SNS_ExtractModesForFolder(dataFolder, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, doDisplay=0)
		else
			SetDataFolder $savedDF
			Abort "SNS_PointLDOS_ZeroField_ResolutionMatrix_FromSettings: provide both STSx and STSy, or neither."
		endif

		if (!SNS__ZeroFieldRes_HasGenericSolverChannels(dataFolder))
			SetDataFolder $savedDF
			Abort "SNS_PointLDOS_ZeroField_ResolutionMatrix_FromSettings: ray tracing did not produce a full solver-ready channel set."
		endif

		sourceDF = dataFolder
		sourceIs3DTagged = 0
	endif

	// -------------------------------------------------------------------------
	// Copy selected source into temporary solver folder.
	// SNS_ComputeDOS_FromSettings expects generic names and w_mask in same DF.
	// -------------------------------------------------------------------------
	String workDF = dataFolder + ":ZeroFieldResWorkChannels"

	KillDataFolder/Z $workDF

	SetDataFolder $dataFolder
	NewDataFolder/O/S ZeroFieldResWorkChannels
	workDF = GetDataFolder(1)

	// Copy mask into the actual current work folder.
	WAVE srcMask = $(dataFolder + ":w_mask")
	Duplicate/O srcMask, w_mask

	if (sourceIs3DTagged)

		WAVE srcL3D    = $(sourceDF + ":L_N_List_3D")
		WAVE srcW3D    = $(sourceDF + ":W_eff_List_3D")
		WAVE srcC3D    = $(sourceDF + ":wChan_3D")
		WAVE srcT3D    = $(sourceDF + ":T_eff_List_3D")
		WAVE srcH1x3D  = $(sourceDF + ":Hit1x_List_3D")
		WAVE srcH1y3D  = $(sourceDF + ":Hit1y_List_3D")
		WAVE srcH2x3D  = $(sourceDF + ":Hit2x_List_3D")
		WAVE srcH2y3D  = $(sourceDF + ":Hit2y_List_3D")

		Duplicate/O srcL3D,   L_N_List
		Duplicate/O srcW3D,   W_eff_List
		Duplicate/O srcC3D,   wChan
		Duplicate/O srcT3D,   T_eff_List
		Duplicate/O srcH1x3D, Hit1x_List
		Duplicate/O srcH1y3D, Hit1y_List
		Duplicate/O srcH2x3D, Hit2x_List
		Duplicate/O srcH2y3D, Hit2y_List

	else

		WAVE srcL    = $(sourceDF + ":L_N_List")
		WAVE srcW    = $(sourceDF + ":W_eff_List")
		WAVE srcC    = $(sourceDF + ":wChan")
		WAVE srcT    = $(sourceDF + ":T_eff_List")
		WAVE srcH1x  = $(sourceDF + ":Hit1x_List")
		WAVE srcH1y  = $(sourceDF + ":Hit1y_List")
		WAVE srcH2x  = $(sourceDF + ":Hit2x_List")
		WAVE srcH2y  = $(sourceDF + ":Hit2y_List")

		Duplicate/O srcL,   L_N_List
		Duplicate/O srcW,   W_eff_List
		Duplicate/O srcC,   wChan
		Duplicate/O srcT,   T_eff_List
		Duplicate/O srcH1x, Hit1x_List
		Duplicate/O srcH1y, Hit1y_List
		Duplicate/O srcH2x, Hit2x_List
		Duplicate/O srcH2y, Hit2y_List

	endif

	if (Exists(sourceDF + ":Vortex_ptx") == 1)
		WAVE srcVx = $(sourceDF + ":Vortex_ptx")
		Duplicate/O srcVx, Vortex_ptx
	endif
	if (Exists(sourceDF + ":Vortex_pty") == 1)
		WAVE srcVy = $(sourceDF + ":Vortex_pty")
		Duplicate/O srcVy, Vortex_pty
	endif

	// Sanity check: this is exactly the folder passed to SNS_ComputeDOS_FromSettings.
	if (Exists(workDF + ":w_mask") != 1)
		SetDataFolder $savedDF
		Abort "SNS_PointLDOS_ZeroField_ResolutionMatrix_FromSettings: internal work folder was created, but w_mask was not copied."
	endif

	SetDataFolder $savedDF

	// -------------------------------------------------------------------------
	// Bind settings.
	// -------------------------------------------------------------------------
	NVAR T_K        = root:SNS_Settings:T_K
	NVAR V_mod      = root:SNS_Settings:V_mod
	NVAR Broadening = root:SNS_Settings:Broadening
	NVAR GammaBase  = root:SNS_Settings:SNS_GammaBase_eV
	NVAR vF         = root:SNS_Settings:vF

	NVAR useDopp    = root:SNS_Settings:SNS_useGammaDopp
	NVAR useZeeman  = root:SNS_Settings:SNS_useGammaZeeman
	NVAR usePair    = root:SNS_Settings:SNS_useGammaPair
	NVAR useUser    = root:SNS_Settings:SNS_useGammaUser
	NVAR useCoh     = root:SNS_Settings:SNS_useGammaCoh
	NVAR useVortex  = root:SNS_Settings:SNS_useVortex
	NVAR nFlux      = root:SNS_Settings:SNS_nFlux

	Variable save_T_K        = T_K
	Variable save_V_mod      = V_mod
	Variable save_Broadening = Broadening
	Variable save_GammaBase  = GammaBase
	Variable save_useDopp    = useDopp
	Variable save_useZeeman  = useZeeman
	Variable save_usePair    = usePair
	Variable save_useUser    = useUser
	Variable save_useCoh     = useCoh
	Variable save_useVortex  = useVortex
	Variable save_nFlux      = nFlux

	Duplicate/O root:SNS_Settings:B_T, root:SNS_Settings:B_T__before_ZeroFieldResSweep

	Make/O/D/N=1 root:SNS_Settings:B_T
	WAVE BfieldWaveZeroFieldRes = root:SNS_Settings:B_T
	BfieldWaveZeroFieldRes[0] = 0
	SetScale/P x, 0, 1, "T", BfieldWaveZeroFieldRes

	V_mod     = 0
	useDopp   = 0
	useZeeman = 0
	usePair   = 0
	useUser   = 0
	useCoh    = 0
	useVortex = 0
	nFlux     = 0

	// -------------------------------------------------------------------------
	// Determine number of columns.
	// -------------------------------------------------------------------------
	Variable nT = floor((Tmax_K - Tmin_K)/T_delta + 1e-9) + 1
	if ((Tmax_K - (nT - 1)*T_delta) > (Tmin_K + 1e-9))
		nT += 1
	endif

	Variable nLphi = floor((LphiMax_um - LphiMin_um)/phase_coherence_delta + 1e-9) + 1
	if ((LphiMin_um + (nLphi - 1)*phase_coherence_delta) < (LphiMax_um - 1e-9))
		nLphi += 1
	endif

	Variable nColsAlloc = nT + nLphi - 1

	// -------------------------------------------------------------------------
	// Output folder and output waves.
	// -------------------------------------------------------------------------
	SetDataFolder $dataFolder
	NewDataFolder/O/S $outFolderName

	Make/O/D/N=(1, nColsAlloc) $(outPrefix + "_LDOS_ERes")
	Make/O/D/N=1               $(outPrefix + "_E_eV")
	Make/O/D/N=(nColsAlloc)    $(outPrefix + "_T_K")
	Make/O/D/N=(nColsAlloc)    $(outPrefix + "_Lphi_um")
	Make/O/D/N=(nColsAlloc)    $(outPrefix + "_Gamma_eV")
	Make/O/D/N=(nColsAlloc)    $(outPrefix + "_Gamma_ueV")
	Make/O/T/N=(nColsAlloc)    $(outPrefix + "_label")

	WAVE LDOSmat    = $(outPrefix + "_LDOS_ERes")
	WAVE Eout       = $(outPrefix + "_E_eV")
	WAVE Tcol       = $(outPrefix + "_T_K")
	WAVE Lphicol    = $(outPrefix + "_Lphi_um")
	WAVE Gcol       = $(outPrefix + "_Gamma_eV")
	WAVE GueVcol    = $(outPrefix + "_Gamma_ueV")
	WAVE/T labelcol = $(outPrefix + "_label")

	Variable col = 0
	Variable i, Tnow, LphiNow, GammaNow, nE = 0
	String dosName, eName, gammaName, broadName

	//==========================================================================
	// 1) Temperature sweep: Tmax_K -> Tmin_K at fixed Lphi = LphiMin_um
	//==========================================================================
	for (i = 0; i < nT; i += 1)

		Tnow = Tmax_K - i*T_delta
		if (Tnow < Tmin_K)
			Tnow = Tmin_K
		endif

		LphiNow  = LphiMin_um
		GammaNow = SNS_GammaHWHM_FromLphi_um(vF, LphiNow)

		T_K        = Tnow
		Broadening = GammaNow
		GammaBase  = GammaNow

		sprintf dosName,   "%s_raw_col%04d", outPrefix, col
		sprintf eName,     "%s_Eraw_col%04d", outPrefix, col
		sprintf broadName, "%s_broad_col%04d", outPrefix, col

		SNS_ComputeDOS_FromSettings(workDF, dosName, eName)

		gammaName = dosName + "_Gamma"
		WAVE DOSgamma = $gammaName

		SNS_ApplyDOS_Broadening_TplusMod(DOSgamma, broadName)

		WAVE DOSbroad = $broadName
		WAVE Eraw     = $eName

		if (col == 0)
			nE = DimSize(DOSbroad, 0)

			Redimension/N=(nE, nColsAlloc) LDOSmat
			Redimension/N=(nE) Eout

			Eout[] = Eraw[p]

			SetScale/P x, DimOffset(DOSbroad,0), DimDelta(DOSbroad,0), WaveUnits(DOSbroad,0), LDOSmat
			SetScale/P y, 0, 1, "resolution index", LDOSmat
			SetScale/P x, DimOffset(Eraw,0), DimDelta(Eraw,0), WaveUnits(Eraw,0), Eout
		endif

		LDOSmat[][col] = DOSbroad[p][0]
		Tcol[col]      = Tnow
		Lphicol[col]   = LphiNow
		Gcol[col]      = GammaNow
		GueVcol[col]   = GammaNow * 1e6
		labelcol[col]  = "T sweep; Lphi=" + num2str(LphiNow) + " um; T=" + num2str(Tnow) + " K"

		KillWaves/Z $dosName, $gammaName, $eName, $broadName
		KillDataFolder/Z E_allBranches
		KillDataFolder/Z m_allBranches
		KillDataFolder/Z s_allBranches

		col += 1

		if (Tnow <= Tmin_K + 1e-12)
			break
		endif
	endfor

	//==========================================================================
	// 2) Lphi sweep at Tmin_K.
	//==========================================================================
	for (i = 1; i < nLphi; i += 1)

		LphiNow = LphiMin_um + i*phase_coherence_delta
		if (LphiNow > LphiMax_um)
			LphiNow = LphiMax_um
		endif

		Tnow     = Tmin_K
		GammaNow = SNS_GammaHWHM_FromLphi_um(vF, LphiNow)

		T_K        = Tnow
		Broadening = GammaNow
		GammaBase  = GammaNow

		sprintf dosName,   "%s_raw_col%04d", outPrefix, col
		sprintf eName,     "%s_Eraw_col%04d", outPrefix, col
		sprintf broadName, "%s_broad_col%04d", outPrefix, col

		SNS_ComputeDOS_FromSettings(workDF, dosName, eName)

		gammaName = dosName + "_Gamma"
		WAVE DOSgamma2 = $gammaName

		SNS_ApplyDOS_Broadening_TplusMod(DOSgamma2, broadName)

		WAVE DOSbroad2 = $broadName

		LDOSmat[][col] = DOSbroad2[p][0]
		Tcol[col]      = Tnow
		Lphicol[col]   = LphiNow
		Gcol[col]      = GammaNow
		GueVcol[col]   = GammaNow * 1e6
		labelcol[col]  = "Lphi sweep; T=" + num2str(Tnow) + " K; Lphi=" + num2str(LphiNow) + " um"

		KillWaves/Z $dosName, $gammaName, $eName, $broadName
		KillDataFolder/Z E_allBranches
		KillDataFolder/Z m_allBranches
		KillDataFolder/Z s_allBranches

		col += 1

		if (LphiNow >= LphiMax_um - 1e-12)
			break
		endif
	endfor

	// -------------------------------------------------------------------------
	// Trim output waves.
	// -------------------------------------------------------------------------
	Redimension/N=(nE, col) LDOSmat
	Redimension/N=(col) Tcol, Lphicol, Gcol, GueVcol, labelcol

	// -------------------------------------------------------------------------
	// Restore original settings and B grid.
	// -------------------------------------------------------------------------
	Duplicate/O root:SNS_Settings:B_T__before_ZeroFieldResSweep, root:SNS_Settings:B_T
	KillWaves/Z root:SNS_Settings:B_T__before_ZeroFieldResSweep

	T_K        = save_T_K
	V_mod      = save_V_mod
	Broadening = save_Broadening
	GammaBase  = save_GammaBase
	useDopp    = save_useDopp
	useZeeman  = save_useZeeman
	usePair    = save_usePair
	useUser    = save_useUser
	useCoh     = save_useCoh
	useVortex  = save_useVortex
	nFlux      = save_nFlux

	KillDataFolder/Z $workDF

	SetDataFolder $savedDF

//	Print "SNS_PointLDOS_ZeroField_ResolutionMatrix_FromSettings:"
//	Print "  source folder:  ", sourceDF
//	Print "  output folder:  ", dataFolder + ":" + outFolderName
//	Print "  matrix:         ", outPrefix + "_LDOS_ERes"
//	Print "  columns:        ", col
//	Print "  vF [m/s]:       ", vF
//	Print "  Gamma(1 um):    ", SNS_GammaHWHM_FromLphi_um(vF, 1)*1e6, " ueV"
//	Print "  Gamma(10 um):   ", SNS_GammaHWHM_FromLphi_um(vF, 10)*1e6, " ueV"
//	Print "  Gamma(100 um):  ", SNS_GammaHWHM_FromLphi_um(vF, 100)*1e6, " ueV"

	return col
End


//==============================================================================
// SNS__ZeroFieldRes_HasGenericSolverChannels
//==============================================================================

Function SNS__ZeroFieldRes_HasGenericSolverChannels(dfPath)
	String dfPath

	if (!DataFolderExists(dfPath + ":"))
		return 0
	endif

	if (Exists(dfPath + ":L_N_List") != 1)
		return 0
	endif
	if (Exists(dfPath + ":W_eff_List") != 1)
		return 0
	endif
	if (Exists(dfPath + ":wChan") != 1)
		return 0
	endif
	if (Exists(dfPath + ":T_eff_List") != 1)
		return 0
	endif
	if (Exists(dfPath + ":Hit1x_List") != 1)
		return 0
	endif
	if (Exists(dfPath + ":Hit1y_List") != 1)
		return 0
	endif
	if (Exists(dfPath + ":Hit2x_List") != 1)
		return 0
	endif
	if (Exists(dfPath + ":Hit2y_List") != 1)
		return 0
	endif

	return 1
End


//==============================================================================
// SNS__ZeroFieldRes_HasTagged3DSolverChannels
//==============================================================================

Function SNS__ZeroFieldRes_HasTagged3DSolverChannels(dfPath)
	String dfPath

	if (!DataFolderExists(dfPath + ":"))
		return 0
	endif

	if (Exists(dfPath + ":L_N_List_3D") != 1)
		return 0
	endif
	if (Exists(dfPath + ":W_eff_List_3D") != 1)
		return 0
	endif
	if (Exists(dfPath + ":wChan_3D") != 1)
		return 0
	endif
	if (Exists(dfPath + ":T_eff_List_3D") != 1)
		return 0
	endif
	if (Exists(dfPath + ":Hit1x_List_3D") != 1)
		return 0
	endif
	if (Exists(dfPath + ":Hit1y_List_3D") != 1)
		return 0
	endif
	if (Exists(dfPath + ":Hit2x_List_3D") != 1)
		return 0
	endif
	if (Exists(dfPath + ":Hit2y_List_3D") != 1)
		return 0
	endif

	return 1
End


//==============================================================================
// SNS_GammaHWHM_FromLphi_um
//
// Gamma = hbar vF / (2 Lphi)
//
// Inputs:
//   vF_mps   : Fermi velocity [m/s]
//   Lphi_um  : phase-coherence length [um]
//
// Output:
//   Gamma [eV]
//==============================================================================

Function SNS_GammaHWHM_FromLphi_um(vF_mps, Lphi_um)
	Variable vF_mps, Lphi_um

	if (Lphi_um <= 0)
		return NaN
	endif
	if (numtype(Lphi_um) != 0)
		return NaN
	endif
	if (numtype(vF_mps) != 0)
		return NaN
	endif

	return HBAR_eVs * vF_mps / (2 * Lphi_um * 1e-6)
End

//==============================================================================
// SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings
//
// Cheap zero-energy resolution matrix.
// Uses existing ray-tracer channel waves directly.
//
// Does NOT calculate full LDOS(E,B).
// Instead computes:
//
//   ZBC(B) = sum_ch weight_ch * Lorentzian[ E_ch(B), Gamma_eff ]
//
// with
//
//   E_ch(B) ≈ Sphi_ch * ( |beta_ch(B)| - pi )
//   beta_ch(B) = 2*pi*B*hEff*abs(W_ch)/Phi0
//   Sphi_ch = 1 / ( 2*(1/Delta + L_ch/(hbar*vF)) )
//
// Field range:
//   Bfirst_T = Phi0/(2*hEff*wMax)
//   Bmax_T   = 1.1 * Bfirst_T
//   calculated range = -Bmax_T ... +Bmax_T
//
// Resolution logic is identical to zero-field:
//   1) Temperature sweep from Tmax_K to Tmin_K at fixed LphiMin_um.
//   2) Lphi sweep from LphiMin_um to LphiMax_um at fixed Tmin_K,
//      skipping duplicate LphiMin_um column.
//
// Reuses existing:
//   SNS__ZeroFieldRes_HasGenericSolverChannels
//   SNS__ZeroFieldRes_HasTagged3DSolverChannels
//   SNS_GammaHWHM_FromLphi_um
//
// Requires reduced helpers:
//   SNS_ZeroEnergyRes_MakeTrace
//   SNS_ZeroEnergyRes_GammaEff_eV
//   SNS_ZeroEnergyRes_dWrep_nm
//   SNS_ZeroEnergyRes_MaxAbsIndex
//   SNS_ZeroEnergyRes_NormalizeColumns
//==============================================================================

Function SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings(dataFolder, [T_delta, phase_coherence_delta, outFolderName, outPrefix, Tmin_K, Tmax_K, LphiMin_um, LphiMax_um, hEff_nm, dB_T, Bangle_deg, BTK_barrier, STSx, STSy, forceRayTrace, useWeights, DeltaEff_eV])
	String dataFolder
	Variable T_delta, phase_coherence_delta
	String outFolderName, outPrefix
	Variable Tmin_K, Tmax_K, LphiMin_um, LphiMax_um
	Variable hEff_nm, dB_T
	Variable Bangle_deg, BTK_barrier, STSx, STSy
	Variable forceRayTrace, useWeights
	Variable DeltaEff_eV

	if (ParamIsDefault(T_delta))
		T_delta = 0.01
	endif
	if (ParamIsDefault(phase_coherence_delta))
		phase_coherence_delta = 1
	endif
	if (ParamIsDefault(outFolderName))
		outFolderName = "ZeroEnergy_ResolutionMatrix"
	endif
	if (ParamIsDefault(outPrefix))
		outPrefix = "ZeroEnergyRes"
	endif
	if (ParamIsDefault(Tmin_K))
		Tmin_K = 0.01
	endif
	if (ParamIsDefault(Tmax_K))
		Tmax_K = 4.2
	endif
	if (ParamIsDefault(LphiMin_um))
		LphiMin_um = 1
	endif
	if (ParamIsDefault(LphiMax_um))
		LphiMax_um = 100
	endif
	if (ParamIsDefault(dB_T))
		dB_T = 0.0005
	endif
	if (ParamIsDefault(Bangle_deg))
		Bangle_deg = 225
	endif
	if (ParamIsDefault(forceRayTrace))
		forceRayTrace = 0
	endif
	if (ParamIsDefault(useWeights))
		useWeights = 0
	endif

	if (T_delta <= 0)
		Abort "SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings: T_delta must be > 0."
	endif
	if (phase_coherence_delta <= 0)
		Abort "SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings: phase_coherence_delta must be > 0."
	endif
	if (Tmax_K < Tmin_K)
		Abort "SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings: Tmax_K must be >= Tmin_K."
	endif
	if (LphiMax_um < LphiMin_um)
		Abort "SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings: LphiMax_um must be >= LphiMin_um."
	endif
	if (dB_T <= 0)
		Abort "SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings: dB_T must be > 0."
	endif

	String savedDF = GetDataFolder(1)

	if (!DataFolderExists("root:SNS_Settings:"))
		SNS_InitDefaultSettings()
	endif

	STRUCT SNS_Params SNS_p
	SNS_LoadParams(SNS_p)

	NVAR T_K        = root:SNS_Settings:T_K
	NVAR V_mod      = root:SNS_Settings:V_mod
	NVAR Broadening = root:SNS_Settings:Broadening
	NVAR GammaBase  = root:SNS_Settings:SNS_GammaBase_eV
	NVAR vF         = root:SNS_Settings:vF

	if (ParamIsDefault(DeltaEff_eV))
		DeltaEff_eV = SNS_p.Delta
	endif

	if (ParamIsDefault(hEff_nm))
		if (Exists("root:v_cfg_h_eff_nm") == 2)
			NVAR hEffCfg = root:v_cfg_h_eff_nm
			hEff_nm = hEffCfg
		else
			hEff_nm = SNS_p.lambdaL * 1e9
		endif
	endif

	if (hEff_nm <= 0)
		SetDataFolder $savedDF
		Abort "SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings: invalid hEff_nm."
	endif

	// -------------------------------------------------------------------------
	// Select existing ray-tracer output.
	// -------------------------------------------------------------------------
	String sourceDF = ""
	Variable sourceIs3DTagged = 0

	String hist2D = dataFolder + ":RayTraceHist"
	String hist3D = dataFolder + ":RayTraceHist_3D"

	if (!forceRayTrace && SNS__ZeroFieldRes_HasGenericSolverChannels(dataFolder))
		sourceDF = dataFolder
		sourceIs3DTagged = 0
	elseif (!forceRayTrace && DataFolderExists(hist2D + ":") && SNS__ZeroFieldRes_HasGenericSolverChannels(hist2D))
		sourceDF = hist2D
		sourceIs3DTagged = 0
	elseif (!forceRayTrace && DataFolderExists(hist3D + ":") && SNS__ZeroFieldRes_HasTagged3DSolverChannels(hist3D))
		sourceDF = hist3D
		sourceIs3DTagged = 1
	else
		if (ParamIsDefault(BTK_barrier) && ParamIsDefault(STSx) && ParamIsDefault(STSy))
			SNS_ExtractModesForFolder(dataFolder, Bangle_deg=Bangle_deg, doDisplay=0)
		elseif (ParamIsDefault(BTK_barrier) && !ParamIsDefault(STSx) && !ParamIsDefault(STSy))
			SNS_ExtractModesForFolder(dataFolder, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, doDisplay=0)
		elseif (!ParamIsDefault(BTK_barrier) && ParamIsDefault(STSx) && ParamIsDefault(STSy))
			SNS_ExtractModesForFolder(dataFolder, Bangle_deg=Bangle_deg, doDisplay=0)
		elseif (!ParamIsDefault(BTK_barrier) && !ParamIsDefault(STSx) && !ParamIsDefault(STSy))
			SNS_ExtractModesForFolder(dataFolder, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, doDisplay=0)
		else
			SetDataFolder $savedDF
			Abort "SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings: provide both STSx and STSy, or neither."
		endif

		if (!SNS__ZeroFieldRes_HasGenericSolverChannels(dataFolder))
			SetDataFolder $savedDF
			Abort "SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings: ray tracing did not produce solver-ready channel waves."
		endif

		sourceDF = dataFolder
		sourceIs3DTagged = 0
	endif

	String Lpath, Wpath, Cpath, Tpath

	if (sourceIs3DTagged)
		Lpath = sourceDF + ":L_N_List_3D"
		Wpath = sourceDF + ":W_eff_List_3D"
		Cpath = sourceDF + ":wChan_3D"
		Tpath = sourceDF + ":T_eff_List_3D"
	else
		Lpath = sourceDF + ":L_N_List"
		Wpath = sourceDF + ":W_eff_List"
		Cpath = sourceDF + ":wChan"
		Tpath = sourceDF + ":T_eff_List"
	endif

	WAVE L_N_List   = $Lpath
	WAVE W_eff_List = $Wpath
	WAVE wChan      = $Cpath
	WAVE T_eff_List = $Tpath

	Variable nCh = numpnts(W_eff_List)
	if (nCh <= 0)
		SetDataFolder $savedDF
		Abort "SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings: no channels."
	endif

	// -------------------------------------------------------------------------
	// Symmetric B range from first zero-bias resonance of max-|w| channel.
	// -------------------------------------------------------------------------
	Variable idxW = SNS_ZeroEnergyRes_MaxAbsIndex(W_eff_List)
	Variable Wrep_m = abs(W_eff_List[idxW])
	Variable Lrep_m = L_N_List[idxW]
	Variable hEff_m = hEff_nm * 1e-9

	Variable Phi0 = 2.067833848e-15
	Variable Bfirst_T = Phi0 / (2*hEff_m*Wrep_m)
	Variable Bmax_T   = 1.5 * Bfirst_T
	Variable Bmin_T   = -Bmax_T

	if (numtype(Bmax_T) != 0 || Bmax_T <= 0)
		SetDataFolder $savedDF
		Abort "SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings: invalid Bmax_T."
	endif

	// Force a symmetric grid with B=0 exactly included.
	Variable nBhalf = ceil(Bmax_T/dB_T)
	if (nBhalf < 1)
		nBhalf = 1
	endif

	Variable nB = 2*nBhalf + 1
	Variable dBactual_T = Bmax_T / nBhalf

	Bmin_T = -Bmax_T

	// -------------------------------------------------------------------------
	// Resolution columns.
	// Same as zero-field:
	//   T sweep: Tmax_K -> Tmin_K at fixed LphiMin_um.
	//   Lphi sweep: Tmin_K at LphiMin_um+step ... LphiMax_um.
	// -------------------------------------------------------------------------
	Variable nT = floor((Tmax_K - Tmin_K)/T_delta + 1e-9) + 1
	if ((Tmax_K - (nT - 1)*T_delta) > (Tmin_K + 1e-9))
		nT += 1
	endif

	Variable nLphi = floor((LphiMax_um - LphiMin_um)/phase_coherence_delta + 1e-9) + 1
	if ((LphiMin_um + (nLphi - 1)*phase_coherence_delta) < (LphiMax_um - 1e-9))
		nLphi += 1
	endif

	Variable nColsAlloc = nT + nLphi - 1

	SetDataFolder $dataFolder
	NewDataFolder/O/S $outFolderName

	Make/O/D/N=(nB, nColsAlloc) $(outPrefix + "_ZBC_BRes")
	Make/O/D/N=(nB, nColsAlloc) $(outPrefix + "_ZBC_BRes_n")
	Make/O/D/N=(nB)             $(outPrefix + "_B_T")

	Make/O/D/N=(nColsAlloc) $(outPrefix + "_T_K")
	Make/O/D/N=(nColsAlloc) $(outPrefix + "_Lphi_um")
	Make/O/D/N=(nColsAlloc) $(outPrefix + "_GammaLife_eV")
	Make/O/D/N=(nColsAlloc) $(outPrefix + "_GammaEff_eV")
	Make/O/D/N=(nColsAlloc) $(outPrefix + "_GammaEff_ueV")
	Make/O/D/N=(nColsAlloc) $(outPrefix + "_dWrep_nm")
	Make/O/T/N=(nColsAlloc) $(outPrefix + "_label")

	WAVE ZBCmat     = $(outPrefix + "_ZBC_BRes")
	WAVE ZBCmatN    = $(outPrefix + "_ZBC_BRes_n")
	WAVE Bwave      = $(outPrefix + "_B_T")
	WAVE Tcol       = $(outPrefix + "_T_K")
	WAVE Lphicol    = $(outPrefix + "_Lphi_um")
	WAVE GlifeCol   = $(outPrefix + "_GammaLife_eV")
	WAVE GeffCol    = $(outPrefix + "_GammaEff_eV")
	WAVE GueVcol    = $(outPrefix + "_GammaEff_ueV")
	WAVE dWrepCol   = $(outPrefix + "_dWrep_nm")
	WAVE/T labelcol = $(outPrefix + "_label")

	Bwave[] = Bmin_T + p*dBactual_T

	SetScale/P x, Bmin_T, dBactual_T, "T", ZBCmat
	SetScale/P y, 0, 1, "resolution index", ZBCmat
	SetScale/P x, Bmin_T, dBactual_T, "T", ZBCmatN
	SetScale/P y, 0, 1, "resolution index", ZBCmatN
	SetScale/P x, Bmin_T, dBactual_T, "T", Bwave

	Variable/G $(outPrefix + "_Bfirst_T") = Bfirst_T
	Variable/G $(outPrefix + "_Bmax_T")   = Bmax_T
	Variable/G $(outPrefix + "_Bmin_T")   = Bmin_T
	Variable/G $(outPrefix + "_dB_T")     = dBactual_T
	Variable/G $(outPrefix + "_Wrep_nm")  = Wrep_m * 1e9
	Variable/G $(outPrefix + "_Lrep_nm")  = Lrep_m * 1e9
	Variable/G $(outPrefix + "_hEff_nm")  = hEff_nm
	Variable/G $(outPrefix + "_DeltaEff_eV") = DeltaEff_eV

	// Save settings affected during sweep.
	Variable save_T_K        = T_K
	Variable save_Broadening = Broadening
	Variable save_GammaBase  = GammaBase

	Variable col = 0
	Variable i, Tnow, LphiNow, GammaLife, GammaEff, dWrep_nm
	String traceName

	//==========================================================================
	// 1) Temperature sweep: Tmax_K -> Tmin_K at fixed Lphi = LphiMin_um
	//==========================================================================
	for (i = 0; i < nT; i += 1)

		Tnow = Tmax_K - i*T_delta
		if (Tnow < Tmin_K)
			Tnow = Tmin_K
		endif

		LphiNow = LphiMin_um

		GammaLife = SNS_GammaHWHM_FromLphi_um(vF, LphiNow)
		GammaEff  = SNS_ZeroEnergyRes_GammaEff_eV(GammaLife, Tnow, V_mod)
		dWrep_nm  = SNS_ZeroEnergyRes_dWrep_nm(GammaEff, Lrep_m, Wrep_m, vF, DeltaEff_eV)

		T_K        = Tnow
		Broadening = GammaLife
		GammaBase  = GammaLife

		sprintf traceName, "%s_trace_col%04d", outPrefix, col

		SNS_ZeroEnergyRes_MakeTrace(L_N_List, W_eff_List, wChan, T_eff_List, traceName, Bmin_T, Bmax_T, dBactual_T, hEff_nm, DeltaEff_eV, vF, GammaEff, useWeights)

		WAVE tr = $traceName
		ZBCmat[][col] = tr[p]

		Tcol[col]     = Tnow
		Lphicol[col]  = LphiNow
		GlifeCol[col] = GammaLife
		GeffCol[col]  = GammaEff
		GueVcol[col]  = GammaEff * 1e6
		dWrepCol[col] = dWrep_nm
		labelcol[col] = "T sweep; Lphi=" + num2str(LphiNow) + " um; T=" + num2str(Tnow) + " K"

		KillWaves/Z $traceName

		col += 1

		if (Tnow <= Tmin_K + 1e-12)
			break
		endif
	endfor

	//==========================================================================
	// 2) Lphi sweep at Tmin_K.
	//==========================================================================
	for (i = 1; i < nLphi; i += 1)

		LphiNow = LphiMin_um + i*phase_coherence_delta
		if (LphiNow > LphiMax_um)
			LphiNow = LphiMax_um
		endif

		Tnow = Tmin_K

		GammaLife = SNS_GammaHWHM_FromLphi_um(vF, LphiNow)
		GammaEff  = SNS_ZeroEnergyRes_GammaEff_eV(GammaLife, Tnow, V_mod)
		dWrep_nm  = SNS_ZeroEnergyRes_dWrep_nm(GammaEff, Lrep_m, Wrep_m, vF, DeltaEff_eV)

		T_K        = Tnow
		Broadening = GammaLife
		GammaBase  = GammaLife

		sprintf traceName, "%s_trace_col%04d", outPrefix, col

		SNS_ZeroEnergyRes_MakeTrace(L_N_List, W_eff_List, wChan, T_eff_List, traceName, Bmin_T, Bmax_T, dBactual_T, hEff_nm, DeltaEff_eV, vF, GammaEff, useWeights)

		WAVE tr2 = $traceName
		ZBCmat[][col] = tr2[p]

		Tcol[col]     = Tnow
		Lphicol[col]  = LphiNow
		GlifeCol[col] = GammaLife
		GeffCol[col]  = GammaEff
		GueVcol[col]  = GammaEff * 1e6
		dWrepCol[col] = dWrep_nm
		labelcol[col] = "Lphi sweep; T=" + num2str(Tnow) + " K; Lphi=" + num2str(LphiNow) + " um"

		KillWaves/Z $traceName

		col += 1

		if (LphiNow >= LphiMax_um - 1e-12)
			break
		endif
	endfor

	Redimension/N=(nB, col) ZBCmat, ZBCmatN
	Redimension/N=(col) Tcol, Lphicol, GlifeCol, GeffCol, GueVcol, dWrepCol, labelcol

	SNS_ZeroEnergyRes_NormalizeColumns(ZBCmat, ZBCmatN)

	T_K        = save_T_K
	Broadening = save_Broadening
	GammaBase  = save_GammaBase

	SetDataFolder $savedDF

	Print "SNS_PointLDOS_ZeroEnergy_ResolutionMatrix_FromSettings:"
	Print "  source folder: ", sourceDF
	Print "  output folder: ", dataFolder + ":" + outFolderName
	Print "  Bfirst [T]:    ", Bfirst_T
	Print "  B range [T]:   ", Bmin_T, " to ", Bmax_T
	Print "  dB actual [T]: ", dBactual_T
	Print "  Wrep [nm]:     ", Wrep_m * 1e9
	Print "  columns:       ", col

	return col
End
//==============================================================================
// SNS_ZeroEnergyRes_MakeTrace
//
// Reduced ZBC(B) estimator from channel rays.
// Symmetric in B: first zero-bias resonance occurs at ±B0.
//
// Uses:
//   betaAbs = | 2*pi*B*hEff*W / Phi0 |
//   E_ch(B) ≈ Sphi_ch * ( betaAbs - pi )
//
// Therefore each channel contributes near both +B0 and -B0.
//==============================================================================

Function SNS_ZeroEnergyRes_MakeTrace(L_N_List, W_eff_List, wChan, T_eff_List, outName, Bmin_T, Bmax_T, dB_T, hEff_nm, DeltaEff_eV, vF_mps, GammaEff_eV, useWeights)
	WAVE L_N_List, W_eff_List, wChan, T_eff_List
	String outName
	Variable Bmin_T, Bmax_T, dB_T
	Variable hEff_nm, DeltaEff_eV, vF_mps, GammaEff_eV
	Variable useWeights

	Variable nB = floor((Bmax_T - Bmin_T)/dB_T + 1e-9) + 1
	if ((Bmin_T + (nB - 1)*dB_T) < (Bmax_T - 1e-12))
		nB += 1
	endif

	Make/O/D/N=(nB) $outName
	WAVE rho = $outName
	rho = 0
	SetScale/P x, Bmin_T, dB_T, "T", rho

	Variable Phi0 = 2.067833848e-15
	Variable hEff_m = hEff_nm * 1e-9

	Variable nCh = numpnts(W_eff_List)
	Variable sumW = sum(wChan)
	if (sumW <= 0)
		sumW = nCh
	endif

	if (GammaEff_eV <= 0)
		GammaEff_eV = 1e-12
	endif

	Variable iB, ch
	Variable Bnow, Lch, WchAbs, betaAbs
	Variable Sphi_eV, Ech_eV, weight, val

	for (iB = 0; iB < nB; iB += 1)

		Bnow = Bmin_T + iB*dB_T
		val = 0

		for (ch = 0; ch < nCh; ch += 1)

			Lch    = L_N_List[ch]
			WchAbs = abs(W_eff_List[ch])

			if (numtype(Lch) != 0 || numtype(WchAbs) != 0 || Lch <= 0 || WchAbs <= 0)
				continue
			endif

			if (sumW > 0)
				weight = wChan[ch] / sumW
			else
				weight = 1 / nCh
			endif

			if (useWeights)
				weight *= T_eff_List[ch]
			endif

			if (weight <= 0 || numtype(weight) != 0)
				continue
			endif

			Sphi_eV = 1 / (2 * (1/DeltaEff_eV + Lch/(HBAR_eVs*vF_mps)))

			// Critical correction:
			// zero-bias resonance at both +B0 and -B0.
			betaAbs = abs(2*pi*Bnow*hEff_m*WchAbs/Phi0)

			Ech_eV = Sphi_eV * (betaAbs - pi)

			val += weight * (GammaEff_eV/pi) / (Ech_eV*Ech_eV + GammaEff_eV*GammaEff_eV)
		endfor

		rho[iB] = val
	endfor

	return nB
End

//==============================================================================
// SNS_ZeroEnergyRes_GammaEff_eV
//
// New reduced-width helper for the cheap ZBC estimator.
// This is only for the reduced trace, not a replacement for SNS_ComputeGammaTot.
//==============================================================================

Function SNS_ZeroEnergyRes_GammaEff_eV(GammaLife_eV, T_K, Vmod_eV)
	Variable GammaLife_eV, T_K, Vmod_eV

	Variable kB_eVperK = 8.617333262e-5

	// HWHM estimate of thermal derivative: FWHM ≈ 3.53 kBT.
	Variable GammaTherm_eV = 1.765 * kB_eVperK * T_K

	// Effective HWHM estimate for lock-in modulation.
	Variable GammaMod_eV = 0.85 * abs(Vmod_eV)

	if (GammaLife_eV < 0 || numtype(GammaLife_eV) != 0)
		GammaLife_eV = 0
	endif

	return sqrt(GammaLife_eV^2 + GammaTherm_eV^2 + GammaMod_eV^2)
End

//==============================================================================
// SNS_ZeroEnergyRes_dWrep_nm
//
// Representative w-resolution at the first zero-bias resonance.
//==============================================================================

Function SNS_ZeroEnergyRes_dWrep_nm(Gamma_eV, L_m, W_m, vF_mps, Delta_eV)
	Variable Gamma_eV, L_m, W_m, vF_mps, Delta_eV

	if (Gamma_eV <= 0 || L_m <= 0 || W_m <= 0 || vF_mps <= 0 || Delta_eV <= 0)
		return NaN
	endif

	Variable dBeta = 2 * Gamma_eV * (1/Delta_eV + L_m/(HBAR_eVs*vF_mps))

	return W_m * dBeta / pi * 1e9
End

//==============================================================================
// SNS_ZeroEnergyRes_MaxAbsIndex
//==============================================================================

Function SNS_ZeroEnergyRes_MaxAbsIndex(w)
	WAVE w

	Variable n = numpnts(w)
	Variable i
	Variable best = 0
	Variable bestAbs = -Inf
	Variable a

	for (i = 0; i < n; i += 1)
		a = abs(w[i])
		if (numtype(a) == 0 && a > bestAbs)
			bestAbs = a
			best = i
		endif
	endfor

	return best
End

//==============================================================================
// SNS_ZeroEnergyRes_NormalizeColumns
//==============================================================================

Function SNS_ZeroEnergyRes_NormalizeColumns(M, Mout)
	WAVE M, Mout

	Variable nx = DimSize(M, 0)
	Variable ny = DimSize(M, 1)
	Variable j, vmax

	Redimension/N=(nx, ny) Mout

	for (j = 0; j < ny; j += 1)
		WaveStats/Q/RMD=[][j,j] M
		vmax = V_max

		if (numtype(vmax) == 0 && vmax > 0)
			Mout[][j] = M[p][j] / vmax
		else
			Mout[][j] = 0
		endif
	endfor

	SetScale/P x, DimOffset(M,0), DimDelta(M,0), WaveUnits(M,0), Mout
	SetScale/P y, DimOffset(M,1), DimDelta(M,1), WaveUnits(M,1), Mout
End








//------------------------------------------------------------------------------
// Normalize data folder path to trailing-colon form.
//------------------------------------------------------------------------------
Function/S SNS_LineDOS_TrailingColon(dfPath)
    String dfPath

    if (strlen(dfPath) == 0)
        return ""
    endif

    if (cmpstr(dfPath[strlen(dfPath)-1, strlen(dfPath)-1], ":") != 0)
        dfPath += ":"
    endif

    return dfPath
End














//======================================================================
// SNS_QPI_UnitCircle_FromSolver
// ---------------------------------------------------------------------
// Build a toy QPI pattern in k-space from the *actual* ABS spectrum
// for each real trajectory (channel).
//
// Inputs:
//   B_T        : 1D wave of B values [T], same as for DOS solver
//   phiList    : 1D wave of channel angles (direction of k) [rad]
//   L_N_List   : 1D wave of N-path lengths per channel [m]
//   W_eff_List : 1D wave of effective magnetic widths per channel [m]
//   wChan      : geometric weights per channel (any scale; normalized inside)
//   T_eff_List : 1D wave of transparencies per channel (0..1)
//   Delta      : gap [eV]
//   vF         : Fermi velocity [m/s]
//   lambdaL    : London penetration depth [m]
//   lambdaF    : Fermi wavelength [m] → sets kF = 2π / λF
//   Ewin       : half-width of energy window around zero [eV]
//                (e.g. 100e-6 for ±100 µeV)
//   Nk         : number of pixels along kx,ky (square grid, −kF..+kF)
//   nameQPI    : name of output 3D wave QPI[kx][ky][B]
//   nameKaxis  : name of k-axis wave (shared by x and y of QPI)
//
// Output:
//   QPI[kx][ky][iB] : accumulated weight of channels with |E|<Ewin
//                     at that B and that direction (±k)
//   kAxis[k]        : k-axis (1/m), scaled −kF..+kF
//
// Uses: Solve_AllBranches_SNS_dGSJ(...) from ABS_in_SNS_SSpS_withBfield_v9.ipf
//======================================================================

//==============================================================================
// SNS_QPI_UnitCircle_FromSolver
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B_T : [T]
//   phiList : input
//   L_N_List : input
//   W_eff_List : input
//   wChan : input
//   T_eff_List : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   lambdaF : input
//   Ewin : input
//   Nk : input
//   nameQPI : input
//   nameKaxis : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_QPI_UnitCircle_FromSolver(B_T, phiList, L_N_List, W_eff_List, wChan, T_eff_List, Delta, vF, lambdaL, lambdaF, Ewin, Nk, nameQPI, nameKaxis)
    Wave    B_T, phiList, L_N_List, W_eff_List, wChan, T_eff_List
    Variable Delta, vF, lambdaL, lambdaF
    Variable Ewin
    Variable Nk
    String  nameQPI, nameKaxis

    Variable nB  = numpnts(B_T)
    Variable Nch = numpnts(phiList)

    if (nB <= 0)
        Abort "SNS_QPI_UnitCircle_FromSolver: B_T has no points."
    endif
    if ((Nch <= 0) || (Nch != numpnts(L_N_List)) || (Nch != numpnts(W_eff_List)) || \
        (Nch != numpnts(wChan)) || (Nch != numpnts(T_eff_List)))
        Abort "SNS_QPI_UnitCircle_FromSolver: channel waves inconsistent."
    endif
    if (Nk < 2)
        Abort "SNS_QPI_UnitCircle_FromSolver: Nk must be >= 2."
    endif
    if (lambdaF <= 0)
        Abort "SNS_QPI_UnitCircle_FromSolver: lambdaF must be > 0."
    endif

    // ---- sanity check: at least one non-zero channel weight ----
    Variable sumW = sum(wChan)
    if (sumW <= 0)
        Abort "SNS_QPI_UnitCircle_FromSolver: all channel weights are zero."
    endif

    // ---- k-axis: unit circle in k-space ----
    Variable kF = 2*pi/lambdaF

    Make/O/D/N=(Nk) $nameKaxis
    Wave kAxis = $nameKaxis
    SetScale/I x, -kF, +kF, kAxis
    kAxis = x

    // 3D QPI cube: kx, ky, B
    Make/O/D/N=(Nk, Nk, nB) $nameQPI
    Wave QPI = $nameQPI
    QPI = 0;

    SetScale/P x, DimOffset(kAxis,0), DimDelta(kAxis,0), "1/m", QPI
    SetScale/P y, DimOffset(kAxis,0), DimDelta(kAxis,0), "1/m", QPI
    SetScale/P z, DimOffset(B_T,0),   DimDelta(B_T,0),   "T",   QPI

    // ---- main channel loop ----
    Variable jCh, iB, k
    Variable Lch, Wch, Teff, weight
    Variable phi, kx, ky
    Variable ix, iy, ixm, iym

    String nameE2D, nameM, nameS

    for (jCh = 0; jCh < Nch; jCh += 1)

        Lch   = L_N_List[jCh]
        Wch   = W_eff_List[jCh]
        Teff  = T_eff_List[jCh]
        weight = wChan[jCh]/sumW

        if ((weight <= 0) || numtype(Lch) || numtype(Wch) || numtype(Teff))
            continue
        endif

        // actual direction of this trajectory
        phi = phiList[jCh]
        kx  = kF*cos(phi)
        ky  = kF*sin(phi)

        // bin indices for +k and -k
        ix  = round((kx - DimOffset(kAxis,0))/DimDelta(kAxis,0))
        iy  = round((ky - DimOffset(kAxis,0))/DimDelta(kAxis,0))
        ixm = round((-kx - DimOffset(kAxis,0))/DimDelta(kAxis,0))
        iym = round((-ky - DimOffset(kAxis,0))/DimDelta(kAxis,0))

        // skip if both ±k are out of range (should not happen for kF in axis)
        if ((ix  < 0 || ix  >= Nk || iy  < 0 || iy  >= Nk) && \
            (ixm < 0 || ixm >= Nk || iym < 0 || iym >= Nk))
            continue
        endif

        // ---- full ABS spectrum for THIS channel from the SNS solver ----
        sprintf nameE2D, "E_allBranches_QPI_ch%03d", jCh
        sprintf nameM,   "m_allBranches_QPI_ch%03d", jCh
        sprintf nameS,   "s_allBranches_QPI_ch%03d", jCh

        Variable nBr = Solve_AllBranches_SNS_dGSJ(B_T, Lch, Wch, Delta, vF, lambdaL, Teff, \
                                                 nameE2D, nameM, nameS)
        if (nBr <= 0)
            continue
        endif

        Wave E_all = $nameE2D   // [iB][branch]

        // ---- for each B: check if ANY branch sits inside |E| <= Ewin ----
        for (iB = 0; iB < nB; iB += 1)

            Variable active = 0
            for (k = 0; k < nBr; k += 1)
                Variable E0 = E_all[iB][k]
                if (numtype(E0) == 0 && abs(E0) <= Ewin)
                    active = 1
                    break
                endif
            endfor

            if (!active)
                continue
            endif

            // add geometric weight at ±k for this B
            if (ix >= 0 && ix < Nk && iy >= 0 && iy < Nk)
                QPI[ix][iy][iB] += abs(weight)
            endif
            if (ixm >= 0 && ixm < Nk && iym >= 0 && iym < Nk)
                QPI[ixm][iym][iB] += abs(weight)
            endif

        endfor // iB

    endfor // jCh
    

    return Nch
End

//==============================================================================
// SNS_QPI_UnitCircle_FromChannelDOS
//
// Build a toy k-space QPI pattern on a unit circle using the *experimental*
// DOS from Compute_DOS_SNS_Map_FromChannels_Coh, channel by channel,
// including T+mod broadening.
//
// For each channel ch:
//   1) Build a 1-channel set: {L_N[ch], W_eff[ch], wChan[ch], T_eff[ch]}.
//   2) Compute DOS_ch(E,B) with Compute_DOS_SNS_Map_FromChannels_Coh.
//   3) Apply SNS_ApplyDOS_Broadening_TplusMod → DOS_ch_broad(E,B).
//   4) Zero DOS outside E ∈ [0, Ewin).
//   5) SumDimension over energy → W_ch(B) = ∫_0^{Ewin} DOS_ch_broad(E,B) dE.
//   6) Deposit W_ch(B) * (wChan[ch]/Σ wChan) at ±(cos φ_ch, sin φ_ch).
//
// OUTPUT:
//   QPI_unit[qx,qy,B] : 3D wave (unit-circle q-space vs B) in root:
//   B_axis            : copy of B_T
//
// Notes:
//   - L_N_List_m and W_eff_List_m must be in *meters*.
//   - Ewin is the positive-bias energy window upper bound [0,Ewin) in eV.
//   - q-space is dimensionless (unit circle [-1,1]) as a directional diagnostic.
//==============================================================================

//==============================================================================
// SNS_QPI_UnitCircle_FromChannelDOS
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B_T : [T]
//   phiList : input
//   L_N_List_m : input
//   W_eff_List_m : input
//   wChan : input
//   T_eff_List : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   Broadening : input
//   NE : input
//   T_K : input
//   Eac_eV : input
//   Ewin : input
//   Nq : input
//   nameQPI : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_QPI_UnitCircle_FromChannelDOS(B_T, phiList, L_N_List_m, W_eff_List_m, wChan, T_eff_List, Delta, vF, lambdaL, Broadening, NE, T_K, Eac_eV, Ewin, Nq,nameQPI)

    Wave    B_T           // [nB] Tesla
    Wave    phiList       // [Nch] radians
    Wave    L_N_List_m    // [Nch] meters
    Wave    W_eff_List_m  // [Nch] meters
    Wave    wChan         // [Nch]
    Wave    T_eff_List    // [Nch]

    Variable Delta        // [eV]
    Variable vF           // [m/s]
    Variable lambdaL      // [m]
    Variable Broadening   // intrinsic Lorentzian [eV]
    Variable NE           // # energy points
    Variable T_K          // temperature [K]
    Variable Eac_eV       // modulation amplitude [eV]
    Variable Ewin         // upper bound of energy window [0,Ewin) in eV
    Variable Nq           // q-grid size (e.g. 256)
    String  nameQPI       // output: QPI_unit

    DFREF dfrCaller = GetDataFolderDFR()
    SetDataFolder root:

    Variable nB  = numpnts(B_T)
    Variable Nch = numpnts(phiList)

    if (nB <= 0)
        Abort "SNS_QPI_UnitCircle_FromChannelDOS: B_T has no points."
    endif
    if ( (Nch <= 0) || (Nch != numpnts(L_N_List_m)) || \
         (Nch != numpnts(W_eff_List_m)) || (Nch != numpnts(wChan)) || \
         (Nch != numpnts(T_eff_List)) )
        Abort "SNS_QPI_UnitCircle_FromChannelDOS: channel waves inconsistent."
    endif
    if (Nq < 2)
        Abort "SNS_QPI_UnitCircle_FromChannelDOS: Nq must be >= 2."
    endif
    if (Delta <= 0 || Ewin <= 0)
        Abort "SNS_QPI_UnitCircle_FromChannelDOS: Delta and Ewin must be > 0."
    endif

    // --- normalize channel weights ---
    Variable sumW = sum(wChan)
    if (sumW <= 0)
        Abort "SNS_QPI_UnitCircle_FromChannelDOS: all channel weights are zero."
    endif

    // --- copy B axis ---
    Duplicate/O B_T, B_axis
    SetScale/P x, DimOffset(B_T,0), DimDelta(B_T,0), "T", B_axis

    // --- allocate QPI cube: qx,qy in [-1,1], z = B index ---
    Make/O/D/N=(Nq, Nq, nB) $nameQPI
    Wave QPI_unit = $nameQPI
    QPI_unit = 0

    SetScale/I x, -1, 1, QPI_unit
    SetScale/I y, -1, 1, QPI_unit
    SetScale/P z, DimOffset(B_axis,0), DimDelta(B_axis,0), "T", QPI_unit

    // --- temp folder for single-channel DOS ---
    String chanFolder = "root:SNS_QPI_ChanTmp"
    NewDataFolder/O $chanFolder

    Variable ch, iB
    Variable phi, nx, ny
    Variable ixPos, iyPos, ixNeg, iyNeg
    Variable weightNorm
    Variable dE

    String nameDOS, nameEaxis, nameDOSbr

    for (ch = 0; ch < Nch; ch += 1)

        if (wChan[ch] <= 0 || numtype(wChan[ch]) != 0)
            continue
        endif

        weightNorm = wChan[ch]/sumW

        phi = phiList[ch]
        nx  = cos(phi)
        ny  = sin(phi)

        // map ±(nx,ny) from [-1,1] to [0..Nq-1]
        ixPos = round( (nx + 1) * (Nq - 1) / 2 )
        iyPos = round( (ny + 1) * (Nq - 1) / 2 )
        ixNeg = round( (-nx + 1) * (Nq - 1) / 2 )
        iyNeg = round( (-ny + 1) * (Nq - 1) / 2 )

        if ( (ixPos < 0 || ixPos >= Nq || iyPos < 0 || iyPos >= Nq) && \
             (ixNeg < 0 || ixNeg >= Nq || iyNeg < 0 || iyNeg >= Nq) )
            continue
        endif

        // --- build single-channel waves in temp folder ---
        SetDataFolder $chanFolder

        Make/O/D/N=1 L_single_m, W_single_m, w_single, T_single
        Wave L_single_m, W_single_m, w_single, T_single

        L_single_m[0] = L_N_List_m[ch]
        W_single_m[0] = W_eff_List_m[ch]
        w_single[0]   = 1           // use weightNorm outside
        T_single[0]   = T_eff_List[ch]

        sprintf nameDOS,   "DOS_ch_%03d", ch
        sprintf nameEaxis, "E_axis_ch_%03d", ch
        sprintf nameDOSbr, "DOS_ch_%03d_broad", ch

        // --- compute coherent DOS(E,B) for this channel ---
        Compute_DOS_SNS_Map_FromChannels_Coh(B_T, L_single_m, W_single_m, \
                                             w_single, T_single, \
                                             Delta, vF, lambdaL, Broadening, NE, \
                                             nameDOS, nameEaxis)

        Wave DOS_ch    = $nameDOS      // [E,B]
        Wave E_axis_ch = $nameEaxis    // [E]

        // --- apply experimental T+mod broadening ---
        SNS_ApplyDOS_Broadening_TplusMod(DOS_ch, nameDOSbr)
        Wave DOS_ch_broad = $nameDOSbr   // [E,B], broadened as in experiment

        // energy step
        dE = abs(DimDelta(E_axis_ch, 0))

        // --- build masked broadened DOS over E ∈ [0, Ewin) ---
        Duplicate/FREE DOS_ch_broad, DOS_win
        DOS_win = ((E_axis_ch[p] >= 0) && (E_axis_ch[p] < Ewin)) ? DOS_ch_broad : 0

        // --- sum along energy dimension (dim 0) for all B at once ---
        Make/FREE/D/N=(nB) EInt
        SumDimension/D=0/DEST=EInt DOS_win   // sum over E → EInt[B]

        // approximate integral over [0,Ewin) by sum * dE
        EInt *= dE

        // --- deposit into QPI_unit for all B ---
        SetDataFolder root:

        for (iB = 0; iB < nB; iB += 1)
            Variable wB = EInt[iB]
            if (wB <= 0 || numtype(wB) != 0)
                continue
            endif

            if (ixPos >= 0 && ixPos < Nq && iyPos >= 0 && iyPos < Nq)
                QPI_unit[ixPos][iyPos][iB] += weightNorm * wB
            endif
            if (ixNeg >= 0 && ixNeg < Nq && iyNeg >= 0 && iyNeg < Nq)
                QPI_unit[ixNeg][iyNeg][iB] += weightNorm * wB
            endif
        endfor

        // back to temp folder, clean up this channel
        SetDataFolder $chanFolder
        KillWaves/Z DOS_ch, DOS_ch_broad, E_axis_ch, L_single_m, W_single_m, w_single, T_single

    endfor

    KillDataFolder/Z $chanFolder
    SetDataFolder dfrCaller
    return Nch
End



////==============================================================================
//// SNS_MapDOS_BFixed_FromMask
////
//// Purpose:
////   Build an energy-window LDOS map at fixed B using the SNS channel solver.
////   Vortex effects enter via the global SNS_useVortex/SNS_nFlux settings and
////   the vortex center (Vortex_ptx/pty); if SNS_useVortex == 0, no vortex phase
////   is applied.
////
////   For each pixel (x,y) inside the N region, the function:
////     • Builds local SNS channels from the mask via
////         SNS_BuildChannelsFromMask2D(...).
////     • Calls SNS_ComputeDOS_FromChannels(...) with Bwin = {B0} to obtain
////         DOS(E, B0) in [states/eV].
////     • Applies experimental T+mod broadening via
////         SNS_ApplyDOS_Broadening_TplusMod(...).
////     • Integrates DOS(E, B0) over an energy window [E0 - dE, E0 + dE] to
////         produce a scalar LDOS value at that pixel.
////
//// Inputs:
////   Nmask    : 2D mask wave defining N-region geometry (1 inside N, 0 outside);
////              axes in nm.
////   phiB     : in-plane field direction [rad].
////   x0, y0   : lower-left corner of map (nm).
////   x1, y1   : upper-right corner of map (nm).
////   dx_nm    : pixel step in x [nm].
////   dy_nm    : pixel step in y [nm].
////   B0       : magnetic field [T].
////   E0       : center energy [eV].
////   dE       : half-width of energy window [eV].
////   stepFac  : step factor passed to SNS_BuildChannelsFromMask2D.
////   Zbarrier : BTK-like interface barrier parameter passed to channel builder.
////
//// Optional:
////   nameOut  : output 2D wave name; if omitted, an auto-generated name
////              "LDOSmap_B...mT_...deg_E...uV_dE...uV" is used.
////
//// Outputs (in caller’s current data folder):
////   nameOut      : 2D LDOS map over (x,y), each value
////                  ≈ ∫_{E0-dE}^{E0+dE} DOS_broadened(E; x,y,B0) dE  [states].
////   E_axis_map   : 1D energy axis [eV] used for the integration.
////
//// Returns:
////   Number of pixels for which channels were successfully built (nDone).
////==============================================================================
//
//Function SNS_MapDOS_BFixed_FromMask(Nmask, phiB, x0, y0, x1, y1, dx_nm, dy_nm, B0, E0, dE, stepFac, Zbarrier, [nameOut])
//    Wave     Nmask
//    Variable phiB
//    Variable x0, y0, x1, y1          // nm
//    Variable dx_nm, dy_nm            // nm
//    Variable B0                      // T
//    Variable E0, dE                  // eV
//    Variable stepFac, Zbarrier
//    String   nameOut
//
//    DFREF dfrCaller = GetDataFolderDFR()
//
//    // ---------------- Logging start ----------------
//    Variable tStart = DateTime
//    String runTag = "MapDOS_Bfixed"
//    runTag += " B=" + num2str(B0) + "T"
//    runTag += " phi=" + num2str(round(phiB*180/pi)) + "deg"
//    runTag += " E0=" + num2str(E0) + " dE=" + num2str(dE)
//    SNS_LogStart(runTag)
//    SNS_Log("Input: x0=" + num2str(x0) + " y0=" + num2str(y0) + \
//            " x1=" + num2str(x1) + " y1=" + num2str(y1) + \
//            " dx=" + num2str(dx_nm) + " dy=" + num2str(dy_nm), level="INFO")
//    SNS_Log("Params: stepFac=" + num2str(stepFac) + " Zbarrier=" + num2str(Zbarrier), level="INFO")
//    // ------------------------------------------------
//
//    // --- Load standard SNS settings ---
//    STRUCT SNS_Params params
//    SNS_LoadParams(params)
//
//    Variable lambdaF = params.LambdaF      // [m]
//    Variable NE      = params.NE
//
//    // geometry / vortex folder = folder of Nmask
//    DFREF dfrGeom = GetWavesDataFolderDFR(Nmask)
//    SetDataFolder dfrGeom
//
//    // vortex center (nm) MUST be available here:
//    Wave Vortex_ptx = Vortex_ptx
//    Wave Vortex_pty = Vortex_pty
//    Variable xV_nm = Vortex_ptx[0]
//    Variable yV_nm = Vortex_pty[0]
//
//    // continuum DOS for this mask (2D DOS × area)
//    Wave/Z w_area_nm2 = w_area_nm2
//    Variable N_cont_statesPer_eV
//    if (WaveExists(w_area_nm2))
//        Variable area_nm2 = w_area_nm2[0]
//        Variable area_m2  = area_nm2 * 1e-18
//        N_cont_statesPer_eV = params.DOS2D_eV_Area * area_m2   // [states/eV]
//    else
//        N_cont_statesPer_eV = 0
//    endif
//
//    // estimate Nphi from mask + lambdaF (in nm) – not used, kept for diagnostics
//    Variable Nphi = round(SNS_EstimateNphi_FromMask(Nmask, lambdaF*1e9))
//    SNS_Log("Nphi estimate from mask: " + num2istr(Nphi), level="INFO")
//
//    // --- output name ---
//    SetDataFolder dfrCaller
//    if (ParamIsDefault(nameOut))
//        String tag
//        tag = "B" + num2str(round(B0*1e3)) + "mT_" \
//            + num2str(round(phiB*180/pi)) + "deg_" \
//            + "E" + num2str(round(E0*1e6)) + "uV_" \
//            + "dE" + num2str(round(dE*1e6)) + "uV"
//        nameOut = "LDOSmap_" + tag
//    endif
//    SNS_Log("Output wave: " + nameOut, level="INFO")
//
//    // --- grid size (at least 1x1) ---
//    Variable Nx = floor( abs(x1 - x0)/dx_nm + 0.5 ) + 1
//    Variable Ny = floor( abs(y1 - y0)/dy_nm + 0.5 ) + 1
//    if (Nx < 1)
//        Nx = 1
//    endif
//    if (Ny < 1)
//        Ny = 1
//    endif
//
//    // Make sure dx,dy follow the direction x0->x1 and y0->y1
//    Variable dxEff = (x1 >= x0) ? abs(dx_nm) : -abs(dx_nm)
//    Variable dyEff = (y1 >= y0) ? abs(dy_nm) : -abs(dy_nm)
//
//    SNS_Log("Grid: Nx=" + num2istr(Nx) + " Ny=" + num2istr(Ny) + \
//            " dxEff=" + num2str(dxEff) + " dyEff=" + num2str(dyEff), level="INFO")
//
//    // --- allocate output map in caller folder ---
//    Make/O/D/N=(Nx, Ny) $nameOut
//    Wave LDOSmap = $nameOut
//    LDOSmap = NaN
//
//    SetScale/P x, x0, dxEff, "nm", LDOSmap
//    SetScale/P y, y0, dyEff, "nm", LDOSmap
//
//    // --- temp folder (like line routine) ---
//    String tmpFolder = "root:SNS_MapTmp"
//    NewDataFolder/O $tmpFolder
//    SetDataFolder $tmpFolder
//
//    // Single B value
//    Make/O/D/N=1 Bwin
//    Bwin[0] = B0
//
//    String nameDOS_EB    = "DOS_local_EB"
//    String nameEaxisLoc  = "E_axis_local"
//    String nameDOS_EBbr  = "DOS_local_EB_broad"
//
//    Variable haveEaxis = 0
//    Variable iE1 = 0, iE2 = 0
//    Variable dEgrid = 1
//    Variable tmp
//
//    Variable ix, iy, iE
//    Variable Nch
//    Variable nDone = 0
//    Variable nSkip = 0
//    Variable nSkipOutside = 0
//
//    Variable logEveryNy = 10
//
//    // We will use mask scaling (nm) to test membership before channel-building
//    Variable px, py
//
//    // Rebind Nmask in its own folder for ScaleToIndex calls
//    SetDataFolder dfrGeom
//    Wave Nmask_loc = Nmask
//
//    for (iy = 0; iy < Ny; iy += 1)
//        Variable r0y = y0 + iy*dyEff
//
//        if (mod(iy, logEveryNy) == 0)
//            SNS_Log("Progress: row " + num2istr(iy) + "/" + num2istr(Ny-1) + \
//                    " (done=" + num2istr(nDone) + ", skip=" + num2istr(nSkip) + \
//                    ", outside=" + num2istr(nSkipOutside) + ")", level="INFO")
//        endif
//
//        for (ix = 0; ix < Nx; ix += 1)
//            Variable r0x = x0 + ix*dxEff
//
//            //============================================================
//            // Early skip: point outside N-region mask -> keep NaN
//            //============================================================
//            px = ScaleToIndex(Nmask_loc, r0x, 0)    // dim0 (x)
//            py = ScaleToIndex(Nmask_loc, r0y, 1)    // dim1 (y)
//
//            // out of bounds => outside
//            if (px < 0 || px >= DimSize(Nmask_loc,0) || py < 0 || py >= DimSize(Nmask_loc,1))
//                nSkipOutside += 1
//                continue
//            endif
//
//            // mask convention: inside N-region = 1, outside = 0
//            if (Nmask_loc[px][py] <= 0)
//                nSkipOutside += 1
//                continue
//            endif
//            //============================================================
//
//            // Build channels for this position (stored in tmpFolder)
//            SetDataFolder $tmpFolder
//            SNS_BuildChannelsFromMask2D(Nmask, r0x, r0y, phiB, stepFac, Zbarrier, tmpFolder)
//
//            // Access channel waves created by builder
//            Wave L_N_List
//            Wave W_eff_List
//            Wave wChan
//            Wave T_eff_List
//            Wave Hit1x_List
//            Wave Hit1y_List
//            Wave Hit2x_List
//            Wave Hit2y_List
//
//            Nch = DimSize(L_N_List, 0)
//
//            // Guard against empty / inconsistent channels
//            if ( (Nch <= 0) || (Nch != DimSize(W_eff_List,0)) || \
//                 (Nch != DimSize(wChan,0)) || (Nch != DimSize(T_eff_List,0)) || \
//                 (Nch != DimSize(Hit1x_List,0)) || (Nch != DimSize(Hit2x_List,0)) )
//                nSkip += 1
//                continue
//            endif
//
//            // Convert nm -> m for solver
//            Duplicate/O L_N_List,   L_N_List_m
//            Duplicate/O W_eff_List, W_eff_List_m
//            L_N_List_m   *= 1e-9
//            W_eff_List_m *= 1e-9
//
//            // Compute raw DOS(E,B) (B dim = 1) with continuum normalization
//            SNS_ComputeDOS_FromChannels( \
//                Bwin, \
//                L_N_List_m, W_eff_List_m, wChan, T_eff_List, \
//                Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, \
//                xV_nm, yV_nm, \
//                N_cont_statesPer_eV, \
//                nameDOS_EB, nameEaxisLoc)
//
//            Wave DOS_EB_raw   = $nameDOS_EB
//            Wave E_axis_local = $nameEaxisLoc
//
//            // Apply experimental T+mod broadening
//            SNS_ApplyDOS_Broadening_TplusMod(DOS_EB_raw, nameDOS_EBbr)
//            Wave DOS_EB_broad = $nameDOS_EBbr
//
//            // Copy energy axis and precompute integration indices once
//            if (!haveEaxis)
//                SetDataFolder dfrCaller
//                Duplicate/O E_axis_local, $"E_axis_map"
//                Wave E_axis_map = $"E_axis_map"
//
//                Variable E1 = E0 - dE
//                Variable E2 = E0 + dE
//                if (E2 < E1)
//                    tmp = E1
//                    E1 = E2
//                    E2 = tmp
//                endif
//
//                iE1 = x2pnt(E_axis_map, E1)
//                iE2 = x2pnt(E_axis_map, E2)
//
//                if (iE1 < 0)
//                    iE1 = 0
//                endif
//                if (iE2 < 0)
//                    iE2 = 0
//                endif
//                if (iE1 > NE-1)
//                    iE1 = NE-1
//                endif
//                if (iE2 > NE-1)
//                    iE2 = NE-1
//                endif
//                if (iE2 < iE1)
//                    tmp = iE1
//                    iE1 = iE2
//                    iE2 = tmp
//                endif
//
//                dEgrid = abs(DimDelta(DOS_EB_broad, 0))    // eV per point
//                if (dEgrid <= 0)
//                    dEgrid = 1
//                endif
//
//                SetDataFolder $tmpFolder
//                haveEaxis = 1
//
//                SNS_Log("Energy window indices: iE1=" + num2istr(iE1) + \
//                        " iE2=" + num2istr(iE2) + " dEgrid=" + num2str(dEgrid), level="INFO")
//            endif
//
//            // Integrate over E-window (discrete sum * dEgrid)
//            Variable acc = 0
//            for (iE = iE1; iE <= iE2; iE += 1)
//                acc += DOS_EB_broad[iE][0]
//            endfor
//            acc *= dEgrid
//
//            // Store into map (caller DF)
//            SetDataFolder dfrCaller
//            LDOSmap[ix][iy] = acc
//            SetDataFolder $tmpFolder
//
//            nDone += 1
//        endfor
//    endfor
//
//    // --- Logging end ---
//    Variable tEnd = DateTime
//    SNS_Log("Finished: done=" + num2istr(nDone) + " skip=" + num2istr(nSkip) + \
//            " outside=" + num2istr(nSkipOutside) + \
//            " runtime_s=" + num2str(tEnd - tStart), level="INFO")
//    SNS_LogEnd()
//    // ------------------
//
//    SetDataFolder dfrCaller
//    return nDone
//End




////==============================================================================
//// SNS_GetMapBounds_FromCoordsAndSTS
////
//// Determine map bounds and step size (x0..x1, y0..y1, dx, dy) from:
////   - coords : 2×(nx*ny) wave of actual spectrum positions in nm (e.g. :pos:STSpos)
////   - STS    : 2D or 3D STS matrix; only dims 0,1 are used for nx,ny and XY grid
////
//// Output:
////   boundsW must be a 1D wave with at least 6 points.
////   boundsW[0]=x0_nm, [1]=y0_nm, [2]=x1_nm, [3]=y1_nm, [4]=dx_nm, [5]=dy_nm
////
//// Method:
////   Uses STMtools DetermineActualCoordinatesOfGrid(img2D, coords) to create img2D_xy,
////   then reads DimOffset/DimDelta from that regular XY grid. :contentReference[oaicite:0]{index=0}
////
//// Optional:
////   killXY : delete the generated *_xy wave after reading scaling (default 0).
////==============================================================================
//Function SNS_GetMapBounds_FromCoordsAndSTS(STS, coords, [killXY])
//    Wave STS
//    Wave coords
//    Variable killXY
//
//    if (ParamIsDefault(killXY))
//        killXY = 0
//    endif
//
//    // dims 0/1 define the spatial grid
//    Variable nx = DimSize(STS, 0)
//    Variable ny = DimSize(STS, 1)
//
//    // coords must be 2×(nx*ny)
//    if (DimSize(coords,0) != 2)
//        Abort "SNS_GetMapBounds_FromCoordsAndSTS: coords must be 2×N."
//    endif
//    if (DimSize(coords,1) != nx*ny)
//        Abort "SNS_GetMapBounds_FromCoordsAndSTS: coords columns must equal nx*ny from STS dims."
//    endif
//
//    // DetermineActualCoordinatesOfGrid expects a 2D image wave.
//    // If STS is 3D, use layer 0 as a dummy image (content doesn't matter for dx/dy).
//    String baseName = NameOfWave(STS)
//    String img2DName = baseName + "_tmp2D_forGrid"
//
//    if (WaveDims(STS) == 2)
//        // Use directly
//        DetermineActualCoordinatesOfGrid(STS, coords)
//        // output wave name will be baseName+"_xy"
//    else
//        // Make a 2D slice [nx,ny] from layer 0
//        Make/O/D/N=(nx,ny) $img2DName
//        Wave img2D = $img2DName
//        img2D[][] = STS[p][q][0]
//
//        DetermineActualCoordinatesOfGrid(img2D, coords)
//        // output wave name will be img2DName+"_xy"
//        baseName = img2DName
//    endif
//
//    String xyName = baseName + "_xy"
//    Wave/Z STSxy = $xyName
//    if (!WaveExists(STSxy))
//        Abort "SNS_GetMapBounds_FromCoordsAndSTS: failed to create " + xyName
//    endif
//
//    // Read regular XY grid scaling (nm). :contentReference[oaicite:1]{index=1}
//    Variable x0 = DimOffset(STSxy, 0)
//    Variable dx = DimDelta (STSxy, 0)
//    Variable y0 = DimOffset(STSxy, 1)
//    Variable dy = DimDelta (STSxy, 1)
//
//    Variable x1 = x0 + dx*(nx-1)
//    Variable y1 = y0 + dy*(ny-1)
//
//	Make/O/D/N=6 wBounds
//    wBounds[0] = x0
//    wBounds[1] = y0
//    wBounds[2] = x1
//    wBounds[3] = y1
//    wBounds[4] = dx
//    wBounds[5] = dy
//
//    if (killXY)
//        KillWaves/Z STSxy
//    endif
//    // clean temp 2D if we created one
//    KillWaves/Z $img2DName
//
//    return 0
//End





//==============================================================================
// SNS_FindPointVortexMinEnergy2D
//
// Scan a trial grid of point-like vortex positions (xV,yV) and find the minimum
// of the zero-temperature ABS free-energy proxy
//
//   C(xV,yV) = - sum_{ch,br} wChan[ch] * |E_ch,br(xV,yV)|
//
// using the existing endpoint phase winding and
// Solve_AllBranches_SNS_dGSJ_betaExtra(...).
//
// Inputs
//   Bval_T         : magnetic field [T] at which to evaluate the cost
//   L_N_List       : channel N lengths [m]
//   W_eff_List     : channel effective widths [m]
//   wChan          : channel weights / multiplicities
//   T_eff_List     : per-channel effective transparency
//   rS1x_nm/y_nm   : channel endpoint 1 [nm]
//   rS2x_nm/y_nm   : channel endpoint 2 [nm]
//   xCand_nm       : 1D wave of trial x positions [nm]
//   yCand_nm       : 1D wave of trial y positions [nm]
//   folder         : output data folder
//
// Optional
//   nFlux          : winding number (default: root:SNS_Settings:SNS_nFlux)
//
// Outputs in 'folder'
//   CostMap                : 2D cost wave, dims = (nX, nY)
//   xCand_nm, yCand_nm     : copies of candidate axes
//   xV_best_nm             : best x [nm]
//   yV_best_nm             : best y [nm]
//   cost_best              : minimum cost
//
// Notes
//   - This is the simplest point-vortex implementation.
//   - Lower (more negative) cost is better.
//   - Start with a coarse x/y candidate grid, then refine near the minimum.
//==============================================================================

Function SNS_FindPointVortexMinEnergy2D(Bval_T, L_N_List, W_eff_List, wChan, T_eff_List, rS1x_nm, rS1y_nm, rS2x_nm, rS2y_nm, xCand_nm, yCand_nm, folder, [nFlux])

    Variable Bval_T
    Wave L_N_List, W_eff_List, wChan, T_eff_List
    Wave rS1x_nm, rS1y_nm, rS2x_nm, rS2y_nm
    Wave xCand_nm, yCand_nm
    String folder
    Variable nFlux

    // ---------- optional input ----------
    if (ParamIsDefault(nFlux))
        NVAR SNS_nFlux = root:SNS_Settings:SNS_nFlux
        nFlux = SNS_nFlux
    endif

    // ---------- sanity ----------
    Variable Nch = numpnts(L_N_List)
    if ((Nch <= 0) || (Nch != numpnts(W_eff_List)) || (Nch != numpnts(wChan)) \
        || (Nch != numpnts(T_eff_List)) \
        || (Nch != numpnts(rS1x_nm)) || (Nch != numpnts(rS1y_nm)) \
        || (Nch != numpnts(rS2x_nm)) || (Nch != numpnts(rS2y_nm)))
        Abort "SNS_FindPointVortexMinEnergy2D: inconsistent channel waves."
    endif

    Variable nX = numpnts(xCand_nm)
    Variable nY = numpnts(yCand_nm)
    if ((nX <= 0) || (nY <= 0))
        Abort "SNS_FindPointVortexMinEnergy2D: empty candidate grid."
    endif

    Variable sumW = sum(wChan)
    if (sumW <= 0)
        Abort "SNS_FindPointVortexMinEnergy2D: sum(wChan) <= 0."
    endif

    // ---------- load standard SNS params ----------
    STRUCT SNS_Params SNS_p
    SNS_LoadParams(SNS_p)

    Variable Delta   = SNS_p.Delta
    Variable vF      = SNS_p.vF
    Variable lambdaL = SNS_p.lambdaL

    // ---------- output folder ----------
    String oldDF = GetDataFolder(1)
    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    Make/O/D/N=(nX, nY) CostMap
    CostMap = NaN

    // ---------- one-point B wave for existing solver ----------
    Make/FREE/D/N=1 B_one
    B_one[0] = Bval_T

    // ---------- reusable temp folder for solver outputs ----------
    String tmpFolder = "root:Packages:SNS:tmpVortexMin2D"
    NewDataFolder/O root:Packages
    NewDataFolder/O root:Packages:SNS
    NewDataFolder/O $tmpFolder

    String nameE2D = tmpFolder + ":E2D_tmp"
    String nameM   = tmpFolder + ":M_tmp"
    String nameS   = tmpFolder + ":S_tmp"

    // ---------- scan ----------
    Variable ix, iy, j, k
    Variable xV_nm, yV_nm
    Variable Lch, Wch, Teff, weight
    Variable betaExtra, th1, th2, dth
    Variable nBr, Eval, costHere

    Variable bestX_local = NaN
    Variable bestY_local = NaN
    Variable bestCost_local = NaN

    for (ix = 0; ix < nX; ix += 1)
        xV_nm = xCand_nm[ix]

        for (iy = 0; iy < nY; iy += 1)
            yV_nm = yCand_nm[iy]

            costHere = 0

            for (j = 0; j < Nch; j += 1)

                weight = wChan[j] / sumW
                if (weight <= 0)
                    continue
                endif

                Lch  = L_N_List[j]
                Wch  = W_eff_List[j]
                Teff = T_eff_List[j]

                if (numtype(Lch) || numtype(Wch) || numtype(Teff))
                    continue
                endif

                // --- same point-vortex endpoint phase as in SNS_ComputeDOS_FromChannels ---
                th1 = atan2(rS1y_nm[j] - yV_nm, rS1x_nm[j] - xV_nm)
                th2 = atan2(rS2y_nm[j] - yV_nm, rS2x_nm[j] - xV_nm)
                dth = th2 - th1

                // wrap to [-pi, pi]
                if (dth > pi)
                    dth -= 2*pi
                elseif (dth < -pi)
                    dth += 2*pi
                endif

                betaExtra = nFlux * dth

                nBr = Solve_AllBranches_SNS_dGSJ_betaExtra( \
                            B_one, Lch, Wch, Delta, vF, lambdaL, Teff, \
                            betaExtra, nameE2D, nameM, nameS)

                if (nBr <= 0)
                    continue
                endif

                Wave E2D_tmp = $nameE2D

                for (k = 0; k < nBr; k += 1)
                    Eval = E2D_tmp[0][k]
                    if (numtype(Eval) == 0)
                        costHere -= weight * abs(Eval)
                    endif
                endfor

            endfor

            CostMap[ix][iy] = costHere

            if ((numtype(bestCost_local) != 0) || (costHere < bestCost_local))
                bestCost_local = costHere
                bestX_local    = xV_nm
                bestY_local    = yV_nm
            endif

        endfor
    endfor

    // ---------- save best result ----------
    Variable/G xV_best_nm = bestX_local
    Variable/G yV_best_nm = bestY_local
    Variable/G cost_best  = bestCost_local

    // optional regular scaling if candidate axes are uniform
    if ((nX > 1) && (nY > 1))
        Variable dxX = xCand_nm[1] - xCand_nm[0]
        Variable dyY = yCand_nm[1] - yCand_nm[0]
        SetScale/P x, xCand_nm[0], dxX, "nm", CostMap
        SetScale/P y, yCand_nm[0], dyY, "nm", CostMap
    endif

    KillDataFolder/Z $tmpFolder
    SetDataFolder $oldDF

    return 0
End


//==============================================================================
// SNS_Usadel2D_SinglePos_FromSettings
//
// Purpose:
//   Diffusive Usadel LDOS(E,B) for one experimental STS position.
//   Reads standard SNS package settings from root:SNS_Settings.
//
// Uses:
//   root:SNS_Settings:B_T
//   root:SNS_Settings:NE
//   root:SNS_Settings:Delta
//   root:SNS_Settings:Broadening
//   root:SNS_Settings:lambdaL   // [m]
//
// Inputs:
//   mask         : cropped/coarse binary mask, scaled in nm
//   STSx, STSy   : real-space STS position, same units as mask scaling
//   outDOS       : output DOS(E,B) wave name
//   outEaxis     : output energy-axis wave name
//   D_nm2_s      : diffusion constant [nm^2/s]
//   GammaEdge_eV : phenomenological edge coupling [eV]
//
// Notes:
//   Energy range is +/- 5*Delta.
//   For this v0 model LDOS is particle-hole symmetric, so only E>=0 is solved.
//==============================================================================

Function SNS_Usadel2D_SinglePos_FromSettings(mask, STSx, STSy, outDOS, outEaxis, D_nm2_s, GammaEdge_eV)
	Wave mask
	String outDOS, outEaxis
	Variable STSx, STSy
	Variable D_nm2_s, GammaEdge_eV

	WAVE Bwave_T = root:SNS_Settings:B_T
	NVAR NE_cfg = root:SNS_Settings:NE
	NVAR Delta = root:SNS_Settings:Delta
	NVAR Broadening = root:SNS_Settings:Broadening
	NVAR lambdaL_m = root:SNS_Settings:lambdaL

	Variable NE = round(NE_cfg)

	// Force odd NE so E=0 is exactly represented.
	if (mod(NE, 2) == 0)
		NE += 1
	endif

	Variable Emin_eV = -2 * Delta
	Variable Emax_eV =  2 * Delta
	Variable eta_eV = Broadening

	Variable nx = DimSize(mask,0)
	Variable ny = DimSize(mask,1)
	Variable nB = DimSize(Bwave_T,0)

	Variable dx_nm = abs(DimDelta(mask,0))
	Variable dy_nm = abs(DimDelta(mask,1))
	Variable x0_nm = DimOffset(mask,0)
	Variable y0_nm = DimOffset(mask,1)

	if (dx_nm <= 0 || dy_nm <= 0)
		Abort "SNS_Usadel2D_SinglePos_FromSettings: mask has invalid wave scaling."
	endif

	if (abs(dx_nm - dy_nm) > 1e-6*dx_nm)
		Abort "SNS_Usadel2D_SinglePos_FromSettings: mask x/y scaling differs; anisotropic Laplacian not implemented."
	endif

	Variable ixSTS = round((STSx - x0_nm) / DimDelta(mask,0))
	Variable iySTS = round((STSy - y0_nm) / DimDelta(mask,1))

	if (ixSTS < 0 || ixSTS >= nx || iySTS < 0 || iySTS >= ny)
		Abort "SNS_Usadel2D_SinglePos_FromSettings: STS position outside mask wave bounds."
	endif

	if (mask[ixSTS][iySTS] < 0.5)
		Abort "SNS_Usadel2D_SinglePos_FromSettings: STS position is outside island mask."
	endif

	Make/O/D/N=(NE) $outEaxis
	Wave E_axis_local = $outEaxis
	SetScale/I x, Emin_eV, Emax_eV, "eV", E_axis_local
	E_axis_local = x

	Make/O/D/N=(NE,nB) $outDOS
	Wave DOS_EB_local = $outDOS
	DOS_EB_local = NaN
	SetScale/I x, Emin_eV, Emax_eV, "eV", DOS_EB_local

	if (nB > 1)
		SetScale/I y, Bwave_T[0], Bwave_T[nB-1], "T", DOS_EB_local
	endif

	Duplicate/O mask, w_usadel_mask
	Duplicate/O mask, w_edgeWeight_Usadel
	SNS_Usadel2D_MakeEdgeWeight(w_usadel_mask, w_edgeWeight_Usadel)

	Make/C/O/N=(nx,ny) theta
	Make/O/D/N=(NE,nB) w_iterCount_Usadel = NaN
	Make/O/D/N=(NE,nB) w_converged_Usadel = 0

	Variable coeff_eVnm2 = 0.5 * HBAR_eVs * D_nm2_s

	// Fast-test defaults. Tighten later if needed.
	Variable maxIter = 50
	Variable tol = 1e-4
	Variable mix = 0.45
	Variable maxStepAllowed = 0.75

	Variable iB, ie, ieMirror
	Variable q_nmInv, q2_nmInv2
	Variable iterUsed, convFlag

	Variable iZero = floor(NE / 2)

	for (iB=0; iB<nB; iB+=1)

		// Uniform in-plane Meissner current:
		// Q = (2e/hbar) B lambdaL.
		// lambdaL is stored in meters. Convert Q from 1/m to 1/nm.
		q_nmInv = (2 * q_e / HBAR_SI) * Bwave_T[iB] * lambdaL_m * 1e-9
		q2_nmInv2 = q_nmInv * q_nmInv

		// Critical fix:
		// At E=0 the proximitized solution should start near theta=pi/2,
		// not theta=0. theta=0 is the normal-state branch and gives
		// unphysical high zero-bias DOS.
		theta = cmplx(pi/2, 0)

		for (ie=iZero; ie<NE; ie+=1)

			SNS_Usadel2D_SolveOneEnergy_UniformQ(theta, w_usadel_mask, w_edgeWeight_Usadel, dx_nm, coeff_eVnm2, q2_nmInv2, E_axis_local[ie], eta_eV, Delta, GammaEdge_eV, maxIter, tol, mix, maxStepAllowed, iterUsed, convFlag)

			DOS_EB_local[ie][iB] = real(cos(theta[ixSTS][iySTS]))
			w_iterCount_Usadel[ie][iB] = iterUsed
			w_converged_Usadel[ie][iB] = convFlag

			// Particle-hole mirror: N(-E,B)=N(+E,B)
			ieMirror = NE - 1 - ie

			if (ieMirror >= 0 && ieMirror < iZero)
				DOS_EB_local[ieMirror][iB] = DOS_EB_local[ie][iB]
				w_iterCount_Usadel[ieMirror][iB] = iterUsed
				w_converged_Usadel[ieMirror][iB] = convFlag
			endif
		endfor
	endfor
End


//==============================================================================
// SNS_Usadel2D_SolveOneEnergy_UniformQ
//
// Purpose:
//   Nonlinear point-relaxed Newton solve for one energy and one uniform Q^2.
//
// Fixes:
//   Uses retarded-safe BCS reservoir Green functions:
//      cosS = z / sqrt(z^2 - Delta^2)
//      sinS = Delta / sqrt(z^2 - Delta^2)
//   with branch chosen so Re(cosS) >= 0 above the gap.
//
// Equation:
//   coeff*(lap - q2*sin(theta)*cos(theta))
//   + i*zE*sin(theta)
//   + Gamma_edge*(sinS*cos(theta) - cosS*sin(theta)) = 0
//==============================================================================

Function SNS_Usadel2D_SolveOneEnergy_UniformQ(theta, mask, edgeW, dx_nm, coeff, q2, E_eV, eta_eV, Delta_eV, GammaEdge_eV, maxIter, tol, mix, maxStepAllowed, iterUsed, convFlag)
	Wave/C theta
	Wave mask, edgeW
	Variable dx_nm, coeff, q2, E_eV, eta_eV, Delta_eV, GammaEdge_eV
	Variable maxIter, tol, mix, maxStepAllowed
	Variable &iterUsed, &convFlag

	Variable nx = DimSize(mask,0)
	Variable ny = DimSize(mask,1)
	Variable dx2 = dx_nm * dx_nm

	Variable iter, ix, iy, nNbr
	Variable maxStep, stepAbs, gammaLoc

	Complex zE = cmplx(E_eV, eta_eV)

	// Retarded-safe BCS reservoir.
	Complex rootBCS = sqrt(zE*zE - cmplx(Delta_eV^2,0))
	Complex cosS = zE / rootBCS
	Complex sinS = Delta_eV / rootBCS

	// Choose physical retarded branch above gap:
	// Re(cosS) should be positive for E>Delta.
	if (E_eV > Delta_eV && real(cosS) < 0)
		rootBCS = -rootBCS
		cosS = zE / rootBCS
		sinS = Delta_eV / rootBCS
	endif

	// For subgap positive energies, keep Im(sin/cos) continuous with same branch.
	if (E_eV >= 0 && imag(rootBCS) < 0)
		rootBCS = -rootBCS
		cosS = zE / rootBCS
		sinS = Delta_eV / rootBCS
	endif

	Complex th, lap, s, c, R, dR, dth

	convFlag = 0
	iterUsed = maxIter

	for (iter=0; iter<maxIter; iter+=1)

		maxStep = 0

		for (ix=0; ix<nx; ix+=1)
			for (iy=0; iy<ny; iy+=1)

				if (mask[ix][iy] < 0.5)
					continue
				endif

				th = theta[ix][iy]
				lap = SNS_Usadel2D_Lap(theta, mask, ix, iy, dx2, nNbr)

				s = sin(th)
				c = cos(th)

				gammaLoc = GammaEdge_eV * edgeW[ix][iy]

				R = coeff * (lap - q2*s*c) + cmplx(0,1)*zE*s + gammaLoc * (sinS*c - cosS*s)

				dR = coeff * (-nNbr/dx2 - q2*cos(2*th)) + cmplx(0,1)*zE*c + gammaLoc * (-sinS*s - cosS*c)

				if (cabs(dR) == 0)
					continue
				endif

				dth = -R / dR

				stepAbs = cabs(dth)
				if (stepAbs > maxStepAllowed)
					dth *= maxStepAllowed / stepAbs
					stepAbs = maxStepAllowed
				endif

				theta[ix][iy] = th + mix*dth
				maxStep = max(maxStep, mix*stepAbs)
			endfor
		endfor

		if (maxStep < tol)
			iterUsed = iter + 1
			convFlag = 1
			break
		endif
	endfor
End

//==============================================================================
// Masked Laplacian
//==============================================================================

Function/C SNS_Usadel2D_Lap(theta, mask, ix, iy, dx2, nNbr)
	Wave/C theta
	Wave mask
	Variable ix, iy, dx2
	Variable &nNbr

	Variable nx = DimSize(mask,0)
	Variable ny = DimSize(mask,1)

	Complex th0 = theta[ix][iy]
	Complex sum = cmplx(0,0)
	nNbr = 0

	if (ix > 0 && mask[ix-1][iy] > 0.5)
		sum += theta[ix-1][iy] - th0
		nNbr += 1
	endif
	if (ix < nx-1 && mask[ix+1][iy] > 0.5)
		sum += theta[ix+1][iy] - th0
		nNbr += 1
	endif
	if (iy > 0 && mask[ix][iy-1] > 0.5)
		sum += theta[ix][iy-1] - th0
		nNbr += 1
	endif
	if (iy < ny-1 && mask[ix][iy+1] > 0.5)
		sum += theta[ix][iy+1] - th0
		nNbr += 1
	endif

	return sum / dx2
End


//==============================================================================
// Edge coupling weight
//==============================================================================

Function SNS_Usadel2D_MakeEdgeWeight(mask, edgeW)
	Wave mask, edgeW

	Variable nx = DimSize(mask,0)
	Variable ny = DimSize(mask,1)
	Variable ix, iy, miss

	edgeW = 0

	for (ix=0; ix<nx; ix+=1)
		for (iy=0; iy<ny; iy+=1)

			if (mask[ix][iy] < 0.5)
				continue
			endif

			miss = 0

			if (ix <= 0 || mask[ix-1][iy] < 0.5)
				miss += 1
			endif
			if (ix >= nx-1 || mask[ix+1][iy] < 0.5)
				miss += 1
			endif
			if (iy <= 0 || mask[ix][iy-1] < 0.5)
				miss += 1
			endif
			if (iy >= ny-1 || mask[ix][iy+1] < 0.5)
				miss += 1
			endif

			edgeW[ix][iy] = miss / 4
		endfor
	endfor
End

//==============================================================================
// SNS_Usadel2D_MakeCroppedCoarseMask
//
// Purpose:
//   Build a smaller mask for the Usadel solver:
//     1. Crop around the connected mask component containing STSx, STSy.
//     2. Add a margin.
//     3. Optionally coarse-grain by integer binFactor.
//
// Inputs:
//   mask        : 2D binary mask, scaled in real units, usually nm
//   STSx, STSy  : position in same real units as mask scaling
//   outMaskName : output wave name, e.g. "w_usadel_mask_crop"
//   margin_nm   : real-space crop margin around connected component
//   binFactor   : integer coarse-graining factor; 1 = no coarse-graining
//
// Output:
//   Creates outMaskName in current data folder.
//   Output keeps real-space scaling.
//==============================================================================

Function SNS_Usadel2D_MakeCroppedCoarseMask(mask, STSx, STSy, outMaskName, margin_nm, binFactor)
	Wave mask
	Variable STSx, STSy
	String outMaskName
	Variable margin_nm, binFactor

	Variable nx = DimSize(mask,0)
	Variable ny = DimSize(mask,1)

	Variable dx = DimDelta(mask,0)
	Variable dy = DimDelta(mask,1)
	Variable x0 = DimOffset(mask,0)
	Variable y0 = DimOffset(mask,1)

	if (dx == 0 || dy == 0)
		Abort "SNS_Usadel2D_MakeCroppedCoarseMask: mask scaling is invalid."
	endif

	Variable adx = abs(dx)
	Variable ady = abs(dy)

	if (abs(adx - ady) > 1e-6*adx)
		Abort "SNS_Usadel2D_MakeCroppedCoarseMask: x/y scaling differs; only square pixels supported."
	endif

	Variable ix0 = round((STSx - x0) / dx)
	Variable iy0 = round((STSy - y0) / dy)

	if (ix0 < 0 || ix0 >= nx || iy0 < 0 || iy0 >= ny)
		Abort "SNS_Usadel2D_MakeCroppedCoarseMask: STS position outside mask bounds."
	endif

	if (mask[ix0][iy0] < 0.5)
		Abort "SNS_Usadel2D_MakeCroppedCoarseMask: STS position is outside mask."
	endif

	binFactor = max(1, round(binFactor))
	Variable marginPix = max(0, ceil(margin_nm / adx))

	// Flood-fill connected component containing STS point.
	Make/O/B/U/N=(nx,ny) w_usadel_component_tmp = 0
	Make/O/D/N=(nx*ny) w_usadel_queueX_tmp, w_usadel_queueY_tmp

	Variable head = 0
	Variable tail = 0
	Variable ix, iy, jx, jy

	w_usadel_queueX_tmp[tail] = ix0
	w_usadel_queueY_tmp[tail] = iy0
	tail += 1
	w_usadel_component_tmp[ix0][iy0] = 1

	Variable minX = ix0
	Variable maxX = ix0
	Variable minY = iy0
	Variable maxY = iy0

	do
		ix = w_usadel_queueX_tmp[head]
		iy = w_usadel_queueY_tmp[head]
		head += 1

		minX = min(minX, ix)
		maxX = max(maxX, ix)
		minY = min(minY, iy)
		maxY = max(maxY, iy)

		// left
		jx = ix - 1
		jy = iy
		if (jx >= 0)
			if (mask[jx][jy] > 0.5 && w_usadel_component_tmp[jx][jy] == 0)
				w_usadel_component_tmp[jx][jy] = 1
				w_usadel_queueX_tmp[tail] = jx
				w_usadel_queueY_tmp[tail] = jy
				tail += 1
			endif
		endif

		// right
		jx = ix + 1
		jy = iy
		if (jx < nx)
			if (mask[jx][jy] > 0.5 && w_usadel_component_tmp[jx][jy] == 0)
				w_usadel_component_tmp[jx][jy] = 1
				w_usadel_queueX_tmp[tail] = jx
				w_usadel_queueY_tmp[tail] = jy
				tail += 1
			endif
		endif

		// down
		jx = ix
		jy = iy - 1
		if (jy >= 0)
			if (mask[jx][jy] > 0.5 && w_usadel_component_tmp[jx][jy] == 0)
				w_usadel_component_tmp[jx][jy] = 1
				w_usadel_queueX_tmp[tail] = jx
				w_usadel_queueY_tmp[tail] = jy
				tail += 1
			endif
		endif

		// up
		jx = ix
		jy = iy + 1
		if (jy < ny)
			if (mask[jx][jy] > 0.5 && w_usadel_component_tmp[jx][jy] == 0)
				w_usadel_component_tmp[jx][jy] = 1
				w_usadel_queueX_tmp[tail] = jx
				w_usadel_queueY_tmp[tail] = jy
				tail += 1
			endif
		endif

	while (head < tail)

	minX = max(0, minX - marginPix)
	maxX = min(nx-1, maxX + marginPix)
	minY = max(0, minY - marginPix)
	maxY = min(ny-1, maxY + marginPix)

	Variable cropNx = maxX - minX + 1
	Variable cropNy = maxY - minY + 1

	if (cropNx <= 1 || cropNy <= 1)
		Abort "SNS_Usadel2D_MakeCroppedCoarseMask: crop is too small."
	endif

	// First create exact crop of the connected component only.
	Make/O/D/N=(cropNx,cropNy) w_usadel_crop_tmp = 0

	for (ix=0; ix<cropNx; ix+=1)
		for (iy=0; iy<cropNy; iy+=1)
			w_usadel_crop_tmp[ix][iy] = w_usadel_component_tmp[minX+ix][minY+iy]
		endfor
	endfor

	Variable cropX0 = x0 + minX*dx
	Variable cropY0 = y0 + minY*dy

	SetScale/P x, cropX0, dx, "", w_usadel_crop_tmp
	SetScale/P y, cropY0, dy, "", w_usadel_crop_tmp

	if (binFactor <= 1)
		Duplicate/O w_usadel_crop_tmp, $outMaskName
		Wave outMask = $outMaskName
		SetScale/P x, cropX0, dx, "", outMask
		SetScale/P y, cropY0, dy, "", outMask
	else
		Variable coarseNx = floor(cropNx / binFactor)
		Variable coarseNy = floor(cropNy / binFactor)

		if (coarseNx < 2 || coarseNy < 2)
			Abort "SNS_Usadel2D_MakeCroppedCoarseMask: binFactor too large for crop."
		endif

		Make/O/D/N=(coarseNx,coarseNy) $outMaskName
		Wave outMask = $outMaskName
		outMask = 0

		Variable cx, cy, bx, by
		Variable sumVal, nVal

		for (cx=0; cx<coarseNx; cx+=1)
			for (cy=0; cy<coarseNy; cy+=1)

				sumVal = 0
				nVal = 0

				for (bx=0; bx<binFactor; bx+=1)
					for (by=0; by<binFactor; by+=1)
						sumVal += w_usadel_crop_tmp[cx*binFactor + bx][cy*binFactor + by]
						nVal += 1
					endfor
				endfor

				outMask[cx][cy] = (sumVal/nVal > 0.5) ? 1 : 0
			endfor
		endfor

		SetScale/P x, cropX0, dx*binFactor, "", outMask
		SetScale/P y, cropY0, dy*binFactor, "", outMask
	endif

	KillWaves/Z w_usadel_component_tmp, w_usadel_queueX_tmp, w_usadel_queueY_tmp, w_usadel_crop_tmp
End

//==============================================================================
// SNS_Usadel2D_PhaseBias_SinglePos_FromSettings
//
// Purpose:
//   Diffusive phase-biased proximity benchmark for one STS position.
//   Maps in-plane Meissner current to superconducting edge phase:
//
//       chiS(r,B) = (2e/hbar) * B * lambdaL * (r · e_perp)
//
//   Solves positive energies only and mirrors DOS to negative energy.
//
// Fix relative to previous version:
//   Computes LDOS from two anomalous amplitudes Fplus and Fminus:
//
//       g = sqrt(1 - Fplus*Fminus)
//       DOS = Re(g)
//
//   NOT sqrt(1-|F|^2), which incorrectly turns coherence peaks into dips.
//==============================================================================

Function SNS_Usadel2D_PhaseBias_SinglePos_FromSettings(mask, STSx, STSy, outDOS, outEaxis, D_nm2_s, GammaEdge_eV, Bangle_deg)
	Wave mask
	String outDOS, outEaxis
	Variable STSx, STSy
	Variable D_nm2_s, GammaEdge_eV, Bangle_deg

	WAVE Bwave_T = root:SNS_Settings:B_T
	NVAR NE_cfg = root:SNS_Settings:NE
	NVAR Delta = root:SNS_Settings:Delta
	NVAR Broadening = root:SNS_Settings:Broadening
	NVAR lambdaL_m = root:SNS_Settings:lambdaL

	Variable NE = round(NE_cfg)

	if (mod(NE, 2) == 0)
		NE += 1
	endif

	Variable Emin_eV = -2 * Delta
	Variable Emax_eV =  2 * Delta
	Variable eta_eV = Broadening

	Variable nx = DimSize(mask,0)
	Variable ny = DimSize(mask,1)
	Variable nB = DimSize(Bwave_T,0)

	Variable dx_nm = abs(DimDelta(mask,0))
	Variable dy_nm = abs(DimDelta(mask,1))
	Variable x0_nm = DimOffset(mask,0)
	Variable y0_nm = DimOffset(mask,1)

	if (dx_nm <= 0 || dy_nm <= 0)
		Abort "SNS_Usadel2D_PhaseBias_SinglePos_FromSettings: mask has invalid wave scaling."
	endif

	if (abs(dx_nm - dy_nm) > 1e-6*dx_nm)
		Abort "SNS_Usadel2D_PhaseBias_SinglePos_FromSettings: mask x/y scaling differs; anisotropic Laplacian not implemented."
	endif

	Variable ixSTS = round((STSx - x0_nm) / DimDelta(mask,0))
	Variable iySTS = round((STSy - y0_nm) / DimDelta(mask,1))

	if (ixSTS < 0 || ixSTS >= nx || iySTS < 0 || iySTS >= ny)
		Abort "SNS_Usadel2D_PhaseBias_SinglePos_FromSettings: STS position outside mask wave bounds."
	endif

	if (mask[ixSTS][iySTS] < 0.5)
		Abort "SNS_Usadel2D_PhaseBias_SinglePos_FromSettings: STS position is outside island mask."
	endif

	Make/O/D/N=(NE) $outEaxis
	Wave E_axis_local = $outEaxis
	SetScale/I x, Emin_eV, Emax_eV, "eV", E_axis_local
	E_axis_local = x

	Make/O/D/N=(NE,nB) $outDOS
	Wave DOS_EB_local = $outDOS
	DOS_EB_local = NaN
	SetScale/I x, Emin_eV, Emax_eV, "eV", DOS_EB_local

	if (nB > 1)
		SetScale/I y, Bwave_T[0], Bwave_T[nB-1], "T", DOS_EB_local
	endif

	Duplicate/O mask, w_usadel_mask
	Duplicate/O mask, w_edgeWeight_Usadel
	SNS_Usadel2D_MakeEdgeWeight(w_usadel_mask, w_edgeWeight_Usadel)

	Make/C/O/N=(nx,ny) Fplus_UsadelPhase
	Make/C/O/N=(nx,ny) Fminus_UsadelPhase
	Make/C/O/N=(nx,ny) FplusZeroSeed_UsadelPhase
	Make/C/O/N=(nx,ny) FminusZeroSeed_UsadelPhase
	Make/C/O/N=(nx,ny) ePhaseEdge_UsadelPhase
	Make/O/D/N=(nx,ny) w_phaseEdge_UsadelPhase

	Make/O/D/N=(NE,nB) w_iterCount_UsadelPhase = NaN
	Make/O/D/N=(NE,nB) w_converged_UsadelPhase = 0

	Variable coeff_eVnm2 = 0.5 * HBAR_eVs * D_nm2_s

	Variable maxIter = 200
	Variable tol = 1e-5
	Variable mix = 0.55

	Variable iZero = floor(NE / 2)
	Variable iB, ie, ieMirror
	Variable iterUsedP, convFlagP
	Variable iterUsedM, convFlagM
	Variable Nloc
	Complex gLoc, prodFM

	Variable phiB = Bangle_deg * pi / 180
	Variable ePerpX = -sin(phiB)
	Variable ePerpY =  cos(phiB)

	Variable uRef_nm = STSx * ePerpX + STSy * ePerpY

	FplusZeroSeed_UsadelPhase = cmplx(0,0)
	FminusZeroSeed_UsadelPhase = cmplx(0,0)

	for (iB=0; iB<nB; iB+=1)

		SNS_Usadel2D_MakeMeissnerEdgePhase(mask, w_edgeWeight_Usadel, ePhaseEdge_UsadelPhase, w_phaseEdge_UsadelPhase, Bwave_T[iB], lambdaL_m, ePerpX, ePerpY, uRef_nm)

		// B-continuation at E=0.
		Fplus_UsadelPhase = FplusZeroSeed_UsadelPhase
		Fminus_UsadelPhase = FminusZeroSeed_UsadelPhase

		for (ie=iZero; ie<NE; ie+=1)

			SNS_Usadel2D_SolveOneEnergy_PhaseF_Sign(Fplus_UsadelPhase, w_usadel_mask, w_edgeWeight_Usadel, ePhaseEdge_UsadelPhase, +1, dx_nm, coeff_eVnm2, E_axis_local[ie], eta_eV, Delta, GammaEdge_eV, maxIter, tol, mix, iterUsedP, convFlagP)

			SNS_Usadel2D_SolveOneEnergy_PhaseF_Sign(Fminus_UsadelPhase, w_usadel_mask, w_edgeWeight_Usadel, ePhaseEdge_UsadelPhase, -1, dx_nm, coeff_eVnm2, E_axis_local[ie], eta_eV, Delta, GammaEdge_eV, maxIter, tol, mix, iterUsedM, convFlagM)

			if (ie == iZero)
				FplusZeroSeed_UsadelPhase = Fplus_UsadelPhase
				FminusZeroSeed_UsadelPhase = Fminus_UsadelPhase
			endif

			prodFM = Fplus_UsadelPhase[ixSTS][iySTS] * Fminus_UsadelPhase[ixSTS][iySTS]
			gLoc = sqrt(1 - prodFM)

			// Physical retarded branch: Re(g) >= 0.
			if (real(gLoc) < 0)
				gLoc = -gLoc
			endif

			Nloc = real(gLoc)

			// Small negative values are numerical branch/noise artifacts.
			if (Nloc < 0)
				Nloc = 0
			endif

			DOS_EB_local[ie][iB] = Nloc
			w_iterCount_UsadelPhase[ie][iB] = max(iterUsedP, iterUsedM)
			w_converged_UsadelPhase[ie][iB] = convFlagP * convFlagM

			ieMirror = NE - 1 - ie

			if (ieMirror >= 0 && ieMirror < iZero)
				DOS_EB_local[ieMirror][iB] = DOS_EB_local[ie][iB]
				w_iterCount_UsadelPhase[ieMirror][iB] = w_iterCount_UsadelPhase[ie][iB]
				w_converged_UsadelPhase[ieMirror][iB] = w_converged_UsadelPhase[ie][iB]
			endif
		endfor
	endfor
End


//==============================================================================
// SNS_Usadel2D_SolveOneEnergy_PhaseF_Sign
//
// Purpose:
//   Linearized complex anomalous Usadel solve for either anomalous sector.
//
// phaseSign = +1:
//      source = fS * exp(+i chiS)
//
// phaseSign = -1:
//      source = fS * exp(-i chiS)
//
// This pair is needed because LDOS uses Fplus*Fminus, not |F|^2.
//==============================================================================

Function SNS_Usadel2D_SolveOneEnergy_PhaseF_Sign(F, mask, edgeW, ePhase, phaseSign, dx_nm, coeff, E_eV, eta_eV, Delta_eV, GammaEdge_eV, maxIter, tol, mix, iterUsed, convFlag)
	Wave/C F
	Wave mask, edgeW
	Wave/C ePhase
	Variable phaseSign
	Variable dx_nm, coeff, E_eV, eta_eV, Delta_eV, GammaEdge_eV
	Variable maxIter, tol, mix
	Variable &iterUsed, &convFlag

	Variable nx = DimSize(mask,0)
	Variable ny = DimSize(mask,1)
	Variable dx2 = dx_nm * dx_nm

	Complex fS = SNS_Usadel2D_BCS_fS_PositiveE(E_eV, eta_eV, Delta_eV)
	Complex zE = cmplx(E_eV, eta_eV)

	Variable iter, ix, iy, nNbr
	Variable maxStep, stepAbs, gammaLoc
	Complex sumNbr, rhs, denom, Fnew, dF, phaseFactor

	convFlag = 0
	iterUsed = maxIter

	for (iter=0; iter<maxIter; iter+=1)

		maxStep = 0

		for (ix=0; ix<nx; ix+=1)
			for (iy=0; iy<ny; iy+=1)

				if (mask[ix][iy] < 0.5)
					continue
				endif

				sumNbr = SNS_Usadel2D_NeighborSum(F, mask, ix, iy, nNbr)
				gammaLoc = GammaEdge_eV * edgeW[ix][iy]

				if (phaseSign >= 0)
					phaseFactor = ePhase[ix][iy]
				else
					phaseFactor = conj(ePhase[ix][iy])
				endif

				rhs = coeff * sumNbr / dx2 + gammaLoc * fS * phaseFactor
				denom = coeff * nNbr / dx2 - cmplx(0,1) * zE + gammaLoc

				if (cabs(denom) == 0)
					continue
				endif

				Fnew = rhs / denom
				dF = Fnew - F[ix][iy]

				F[ix][iy] += mix * dF

				stepAbs = cabs(dF)
				maxStep = max(maxStep, mix * stepAbs)
			endfor
		endfor

		if (maxStep < tol)
			iterUsed = iter + 1
			convFlag = 1
			break
		endif
	endfor
End

//==============================================================================
// SNS_Usadel2D_MakeMeissnerEdgePhase
//
// Purpose:
//   Build exp(i chiS) on edge pixels for the current B field.
//   chiS = (2e/hbar) B lambdaL (r dot e_perp - reference)
//
// Units:
//   B          : tesla
//   lambdaL_m : meters
//   coordinates from mask scaling are assumed nm
//==============================================================================

Function SNS_Usadel2D_MakeMeissnerEdgePhase(mask, edgeW, ePhase, phaseW, B_T, lambdaL_m, ePerpX, ePerpY, uRef_nm)
	Wave mask, edgeW, phaseW
	Wave/C ePhase
	Variable B_T, lambdaL_m, ePerpX, ePerpY, uRef_nm

	Variable nx = DimSize(mask,0)
	Variable ny = DimSize(mask,1)

	Variable ix, iy
	Variable x_nm, y_nm, u_nm
	Variable q_nmInv, phaseVal

	// Qphase = (2e/hbar) B lambdaL, converted from 1/m to 1/nm.
	q_nmInv = (2 * q_e / HBAR_SI) * B_T * lambdaL_m * 1e-9

	for (ix=0; ix<nx; ix+=1)
		for (iy=0; iy<ny; iy+=1)

			if (mask[ix][iy] > 0.5 && edgeW[ix][iy] > 0)
				x_nm = DimOffset(mask,0) + ix * DimDelta(mask,0)
				y_nm = DimOffset(mask,1) + iy * DimDelta(mask,1)

				u_nm = x_nm * ePerpX + y_nm * ePerpY
				phaseVal = q_nmInv * (u_nm - uRef_nm)

				phaseW[ix][iy] = phaseVal
				ePhase[ix][iy] = cmplx(cos(phaseVal), sin(phaseVal))
			else
				phaseW[ix][iy] = NaN
				ePhase[ix][iy] = cmplx(0,0)
			endif
		endfor
	endfor
End


//==============================================================================
// SNS_Usadel2D_SolveOneEnergy_PhaseF
//
// Purpose:
//   Linearized complex anomalous Usadel solve for one positive energy.
//   Solves:
//
//      coeff * lap(F) + i(E+iη) F
//      - Gamma_edge F
//      + Gamma_edge fS exp(i chiS) = 0
//
//   The phase-biased reservoirs enter via exp(i chiS).
//   This avoids theta-branch jumping and is suitable as a first robust
//   le-Sueur-style phase-bias benchmark.
//
// Branch choice:
//   fS is chosen continuous on positive energy.
//==============================================================================

Function SNS_Usadel2D_SolveOneEnergy_PhaseF(F, mask, edgeW, ePhase, dx_nm, coeff, E_eV, eta_eV, Delta_eV, GammaEdge_eV, maxIter, tol, mix, iterUsed, convFlag)
	Wave/C F
	Wave mask, edgeW
	Wave/C ePhase
	Variable dx_nm, coeff, E_eV, eta_eV, Delta_eV, GammaEdge_eV
	Variable maxIter, tol, mix
	Variable &iterUsed, &convFlag

	Variable nx = DimSize(mask,0)
	Variable ny = DimSize(mask,1)
	Variable dx2 = dx_nm * dx_nm

	Complex fS = SNS_Usadel2D_BCS_fS_PositiveE(E_eV, eta_eV, Delta_eV)
	Complex zE = cmplx(E_eV, eta_eV)

	Variable iter, ix, iy, nNbr
	Variable maxStep, stepAbs, gammaLoc
	Complex sumNbr, rhs, denom, Fnew, dF

	convFlag = 0
	iterUsed = maxIter

	for (iter=0; iter<maxIter; iter+=1)

		maxStep = 0

		for (ix=0; ix<nx; ix+=1)
			for (iy=0; iy<ny; iy+=1)

				if (mask[ix][iy] < 0.5)
					continue
				endif

				sumNbr = SNS_Usadel2D_NeighborSum(F, mask, ix, iy, nNbr)

				gammaLoc = GammaEdge_eV * edgeW[ix][iy]

				rhs = coeff * sumNbr / dx2 + gammaLoc * fS * ePhase[ix][iy]
				denom = coeff * nNbr / dx2 - cmplx(0,1) * zE + gammaLoc

				if (cabs(denom) == 0)
					continue
				endif

				Fnew = rhs / denom
				dF = Fnew - F[ix][iy]

				F[ix][iy] += mix * dF

				stepAbs = cabs(dF)
				maxStep = max(maxStep, mix * stepAbs)
			endfor
		endfor

		if (maxStep < tol)
			iterUsed = iter + 1
			convFlag = 1
			break
		endif
	endfor
End


//==============================================================================
// SNS_Usadel2D_BCS_fS_PositiveE
//
// Purpose:
//   Positive-energy retarded BCS anomalous amplitude for boundary source.
//   Uses fS = Delta / sqrt(Delta^2 - (E+iη)^2)
//   with branch chosen so fS is positive real near E=0.
//==============================================================================

Function/C SNS_Usadel2D_BCS_fS_PositiveE(E_eV, eta_eV, Delta_eV)
	Variable E_eV, eta_eV, Delta_eV

	Complex zE = cmplx(E_eV, eta_eV)
	Complex denom = sqrt(cmplx(Delta_eV^2,0) - zE*zE)
	Complex fS = Delta_eV / denom

	if (E_eV < Delta_eV && real(fS) < 0)
		fS = -fS
	endif

	return fS
End


//==============================================================================
// SNS_Usadel2D_NeighborSum
//
// Purpose:
//   Sum nearest-neighbour F values inside mask.
//==============================================================================

Function/C SNS_Usadel2D_NeighborSum(F, mask, ix, iy, nNbr)
	Wave/C F
	Wave mask
	Variable ix, iy
	Variable &nNbr

	Variable nx = DimSize(mask,0)
	Variable ny = DimSize(mask,1)

	Complex sumNbr = cmplx(0,0)
	nNbr = 0

	if (ix > 0 && mask[ix-1][iy] > 0.5)
		sumNbr += F[ix-1][iy]
		nNbr += 1
	endif

	if (ix < nx-1 && mask[ix+1][iy] > 0.5)
		sumNbr += F[ix+1][iy]
		nNbr += 1
	endif

	if (iy > 0 && mask[ix][iy-1] > 0.5)
		sumNbr += F[ix][iy-1]
		nNbr += 1
	endif

	if (iy < ny-1 && mask[ix][iy+1] > 0.5)
		sumNbr += F[ix][iy+1]
		nNbr += 1
	endif

	return sumNbr
End


//==============================================================================
// SNS_EstimateDiffusivePhaseBias_FromMask
//
// Purpose:
//   Fast mask-only estimate of the diffusive phase-bias field scale.
//   Uses the connected mask component containing STSx, STSy.
//   Evaluates edge phase spread along e_perp, the Meissner-current direction.
//
// Physics:
//   chiS(r,B) = (2e/hbar) B lambdaL (r · e_perp)
//   First diffusive minigap closing estimate:
//      Delta chi ~ pi
//      Bpi = pi*hbar / (2e lambdaL Wperp)
//
// Inputs:
//   mask       : 2D binary mask, scaled in nm
//   STSx, STSy : STS position in same physical units as mask scaling
//   Bangle_deg : in-plane B direction in degrees
//   outPrefix  : prefix for scalar outputs, e.g. "usadelPhaseEst"
//
// Uses:
//   root:SNS_Settings:lambdaL   // [m]
//
// Outputs in current data folder:
//   <outPrefix>_Wrange_nm
//   <outPrefix>_Wrms_nm
//   <outPrefix>_Bpi_range_T
//   <outPrefix>_Bpi_rms_T
//   <outPrefix>_uMean_nm
//   <outPrefix>_uMin_nm
//   <outPrefix>_uMax_nm
//   <outPrefix>_nEdge
//
// Optional diagnostic waves:
//   w_diffPhase_component
//   w_diffPhase_edgeWeight
//==============================================================================

Function SNS_EstimateDiffusivePhaseBias_FromMask(mask, STSx, STSy, Bangle_deg, outPrefix)
	Wave mask
	Variable STSx, STSy, Bangle_deg
	String outPrefix

	NVAR lambdaL_m = root:SNS_Settings:lambdaL

	Variable nx = DimSize(mask,0)
	Variable ny = DimSize(mask,1)

	Variable dx = DimDelta(mask,0)
	Variable dy = DimDelta(mask,1)
	Variable x0 = DimOffset(mask,0)
	Variable y0 = DimOffset(mask,1)

	if (dx == 0 || dy == 0)
		Abort "SNS_EstimateDiffusivePhaseBias_FromMask: mask scaling is invalid."
	endif

	Variable adx = abs(dx)
	Variable ady = abs(dy)

	if (abs(adx - ady) > 1e-6*adx)
		Abort "SNS_EstimateDiffusivePhaseBias_FromMask: x/y scaling differs; only square pixels supported."
	endif

	Variable ix0 = round((STSx - x0) / dx)
	Variable iy0 = round((STSy - y0) / dy)

	if (ix0 < 0 || ix0 >= nx || iy0 < 0 || iy0 >= ny)
		Abort "SNS_EstimateDiffusivePhaseBias_FromMask: STS position outside mask bounds."
	endif

	if (mask[ix0][iy0] < 0.5)
		Abort "SNS_EstimateDiffusivePhaseBias_FromMask: STS position is outside mask."
	endif

	// -------------------------------------------------------------------------
	// Flood-fill connected component containing STS point.
	// -------------------------------------------------------------------------
	Make/O/B/U/N=(nx,ny) w_diffPhase_component = 0
	Make/O/D/N=(nx*ny) w_diffPhase_queueX_tmp, w_diffPhase_queueY_tmp

	Variable head = 0
	Variable tail = 0
	Variable ix, iy, jx, jy

	w_diffPhase_queueX_tmp[tail] = ix0
	w_diffPhase_queueY_tmp[tail] = iy0
	tail += 1
	w_diffPhase_component[ix0][iy0] = 1

	do
		ix = w_diffPhase_queueX_tmp[head]
		iy = w_diffPhase_queueY_tmp[head]
		head += 1

		// left
		jx = ix - 1
		jy = iy
		if (jx >= 0)
			if (mask[jx][jy] > 0.5 && w_diffPhase_component[jx][jy] == 0)
				w_diffPhase_component[jx][jy] = 1
				w_diffPhase_queueX_tmp[tail] = jx
				w_diffPhase_queueY_tmp[tail] = jy
				tail += 1
			endif
		endif

		// right
		jx = ix + 1
		jy = iy
		if (jx < nx)
			if (mask[jx][jy] > 0.5 && w_diffPhase_component[jx][jy] == 0)
				w_diffPhase_component[jx][jy] = 1
				w_diffPhase_queueX_tmp[tail] = jx
				w_diffPhase_queueY_tmp[tail] = jy
				tail += 1
			endif
		endif

		// down
		jx = ix
		jy = iy - 1
		if (jy >= 0)
			if (mask[jx][jy] > 0.5 && w_diffPhase_component[jx][jy] == 0)
				w_diffPhase_component[jx][jy] = 1
				w_diffPhase_queueX_tmp[tail] = jx
				w_diffPhase_queueY_tmp[tail] = jy
				tail += 1
			endif
		endif

		// up
		jx = ix
		jy = iy + 1
		if (jy < ny)
			if (mask[jx][jy] > 0.5 && w_diffPhase_component[jx][jy] == 0)
				w_diffPhase_component[jx][jy] = 1
				w_diffPhase_queueX_tmp[tail] = jx
				w_diffPhase_queueY_tmp[tail] = jy
				tail += 1
			endif
		endif

	while (head < tail)

	KillWaves/Z w_diffPhase_queueX_tmp, w_diffPhase_queueY_tmp

	// -------------------------------------------------------------------------
	// Build edge weights for connected component only.
	// Missing-neighbour count / 4, same philosophy as Usadel edge coupling.
	// -------------------------------------------------------------------------
	Make/O/D/N=(nx,ny) w_diffPhase_edgeWeight = 0

	Variable miss
	for (ix=0; ix<nx; ix+=1)
		for (iy=0; iy<ny; iy+=1)

			if (w_diffPhase_component[ix][iy] < 0.5)
				continue
			endif

			miss = 0

			if (ix <= 0 || w_diffPhase_component[ix-1][iy] < 0.5)
				miss += 1
			endif
			if (ix >= nx-1 || w_diffPhase_component[ix+1][iy] < 0.5)
				miss += 1
			endif
			if (iy <= 0 || w_diffPhase_component[ix][iy-1] < 0.5)
				miss += 1
			endif
			if (iy >= ny-1 || w_diffPhase_component[ix][iy+1] < 0.5)
				miss += 1
			endif

			w_diffPhase_edgeWeight[ix][iy] = miss / 4
		endfor
	endfor

	// -------------------------------------------------------------------------
	// Evaluate edge u-distribution.
	// -------------------------------------------------------------------------
	Variable phiB = Bangle_deg * pi / 180
	Variable ePerpX = -sin(phiB)
	Variable ePerpY =  cos(phiB)

	Variable x_nm, y_nm, u_nm
	Variable w, wSum = 0
	Variable uSum = 0
	Variable u2Sum = 0
	Variable uMin = Inf
	Variable uMax = -Inf
	Variable nEdge = 0

	for (ix=0; ix<nx; ix+=1)
		for (iy=0; iy<ny; iy+=1)

			w = w_diffPhase_edgeWeight[ix][iy]

			if (w <= 0)
				continue
			endif

			x_nm = x0 + ix*dx
			y_nm = y0 + iy*dy
			u_nm = x_nm*ePerpX + y_nm*ePerpY

			wSum += w
			uSum += w*u_nm
			u2Sum += w*u_nm*u_nm
			uMin = min(uMin, u_nm)
			uMax = max(uMax, u_nm)
			nEdge += 1
		endfor
	endfor

	if (wSum <= 0)
		Abort "SNS_EstimateDiffusivePhaseBias_FromMask: no edge pixels found."
	endif

	Variable uMean = uSum / wSum
	Variable uVar = u2Sum / wSum - uMean*uMean
	uVar = max(uVar, 0)

	Variable Wrange_nm = uMax - uMin
	Variable Wrms_nm = 2 * sqrt(uVar)

	Variable Bpi_range_T = NaN
	Variable Bpi_rms_T = NaN

	if (Wrange_nm > 0)
		Bpi_range_T = pi * HBAR_SI / (2 * q_e * lambdaL_m * Wrange_nm * 1e-9)
	endif

	if (Wrms_nm > 0)
		Bpi_rms_T = pi * HBAR_SI / (2 * q_e * lambdaL_m * Wrms_nm * 1e-9)
	endif

	// -------------------------------------------------------------------------
	// Store scalar outputs as globals.
	// -------------------------------------------------------------------------
	Variable/G $(outPrefix + "_Wrange_nm") = Wrange_nm
	Variable/G $(outPrefix + "_Wrms_nm") = Wrms_nm
	Variable/G $(outPrefix + "_Bpi_range_T") = Bpi_range_T
	Variable/G $(outPrefix + "_Bpi_rms_T") = Bpi_rms_T
	Variable/G $(outPrefix + "_uMean_nm") = uMean
	Variable/G $(outPrefix + "_uMin_nm") = uMin
	Variable/G $(outPrefix + "_uMax_nm") = uMax
	Variable/G $(outPrefix + "_nEdge") = nEdge

	Print "SNS diffusive phase-bias estimate:"
	Print "  Wrange_nm   = ", Wrange_nm
	Print "  Wrms_nm     = ", Wrms_nm
	Print "  Bpi_range_T = ", Bpi_range_T
	Print "  Bpi_rms_T   = ", Bpi_rms_T
	Print "  nEdge       = ", nEdge
End

//==============================================================================
// SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings
//
// Purpose:
//   1D phase-biased Usadel SNS benchmark over root:SNS_Settings:B_T,
//   or for one optional single B value.
//
//   Uses phase folding:
//      phiFold = abs(atan2(sin(phiRaw), cos(phiRaw)))
//
//   Therefore the solver only follows the physical branch from phi=0 to pi.
//   This avoids branch jumps after gap closing and fills the full B map by
//   the symmetry N(E,phi)=N(E,-phi)=N(E,2pi-phi).
//
// Important length separation:
//   Ldiff_nm  : Usadel diffusion length.
//               Sets dx, diffusion operator, Thouless scale, zero-field LDOS.
//
//   Lphase_nm : magnetic phase lever arm.
//               Sets phi(B) only.
//               If omitted, Lphase_nm = Ldiff_nm for backward compatibility.
//
// Output waves left in caller folder:
//   outDOS
//   outDOS+"_sym"
//   outEaxis
//
// Diagnostics moved into:
//   :Usadel1D:
//
// Optional:
//   etaOverride_eV : overrides root:SNS_Settings:Broadening
//   energySign     : +1 default, -1 tests opposite convention
//   singleB_T      : if supplied, calculate only this single B value
//   Lphase_nm      : phase lever arm for phi(B); default = Ldiff_nm
//
// Requires existing core helpers:
//   SNS_Usadel1D_Phi_SolveOneEnergyNewton
//   SNS_Usadel1D_Phi_SetLinearChi
//   SNS_Usadel1D_Phi_BCS_Theta_PositiveE
//==============================================================================

Function SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings(Ldiff_nm, xProbe_nm, outDOS, outEaxis, D_nm2_s, Ngrid, [etaOverride_eV, energySign, singleB_T, Lphase_nm])
    Variable Ldiff_nm, xProbe_nm, D_nm2_s, Ngrid
    Variable etaOverride_eV, energySign, singleB_T, Lphase_nm
    String outDOS, outEaxis

    Wave/Z Bwave_T = root:SNS_Settings:B_T
    NVAR NE_cfg      = root:SNS_Settings:NE
    NVAR Delta       = root:SNS_Settings:Delta
    NVAR Broadening  = root:SNS_Settings:Broadening
    NVAR lambdaL_m   = root:SNS_Settings:lambdaL

    Variable useSingleB = !ParamIsDefault(singleB_T)

    if (!useSingleB && !WaveExists(Bwave_T))
        Abort "SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings: root:SNS_Settings:B_T does not exist."
    endif

    Variable signE = 1
    if (!ParamIsDefault(energySign))
        if (energySign < 0)
            signE = -1
        endif
    endif

    Variable NE = round(NE_cfg)
    if (mod(NE, 2) == 0)
        NE += 1
    endif

    Ngrid = max(21, round(Ngrid))
    if (mod(Ngrid, 2) == 0)
        Ngrid += 1
    endif

    if (Ldiff_nm <= 0)
        Abort "SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings: Ldiff_nm must be positive."
    endif

    Variable LphaseUse_nm = Ldiff_nm
    if (!ParamIsDefault(Lphase_nm))
        LphaseUse_nm = Lphase_nm
    endif

    if (LphaseUse_nm <= 0)
        Abort "SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings: Lphase_nm must be positive."
    endif

    Variable eta_eV = Broadening
    if (!ParamIsDefault(etaOverride_eV))
        eta_eV = etaOverride_eV
    endif

    Variable nB
    if (useSingleB)
        nB = 1
    else
        nB = DimSize(Bwave_T,0)
    endif

    if (nB <= 0)
        Abort "SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings: no B values."
    endif

    Variable Emin_eV = -2.5e-3
    Variable Emax_eV =  2.5e-3
    Variable dx_nm = Ldiff_nm / (Ngrid - 1)

    if (xProbe_nm < 0 || xProbe_nm > Ldiff_nm)
        Abort "SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings: xProbe_nm outside [0,Ldiff_nm]."
    endif

    Variable iProbe = round(xProbe_nm / dx_nm)
    iProbe = max(0, min(Ngrid-1, iProbe))

    Make/O/D/N=(NE) $outEaxis
    Wave E_axis_local = $outEaxis
    SetScale/I x, Emin_eV, Emax_eV, "eV", E_axis_local
    E_axis_local = x

    Make/O/D/N=(NE,nB) $outDOS
    Wave DOS_EB_local = $outDOS
    DOS_EB_local = NaN
    SetScale/I x, Emin_eV, Emax_eV, "eV", DOS_EB_local

    if (useSingleB)
        SetScale/P y, singleB_T, 1, "T", DOS_EB_local
    else
        if (nB > 1)
            SetScale/I y, Bwave_T[0], Bwave_T[nB-1], "T", DOS_EB_local
        else
            SetScale/P y, Bwave_T[0], 1, "T", DOS_EB_local
        endif
    endif

    Make/C/O/N=(Ngrid) theta_Phi
    Make/C/O/N=(Ngrid) chi_Phi

    Make/C/O/N=(NE,Ngrid) thetaSeed_Usadel1D_Phi
    Make/C/O/N=(NE,Ngrid) chiSeed_Usadel1D_Phi

    Make/O/D/N=(NE,nB) w_iterCount_Usadel1D_Phi = NaN
    Make/O/D/N=(NE,nB) w_converged_Usadel1D_Phi = 0
    Make/O/D/N=(NE,nB) w_resMax_Usadel1D_Phi = NaN
    Make/O/D/N=(NE,nB) w_stepMax_Usadel1D_Phi = NaN

    Make/O/D/N=(nB) w_B_Usadel1D_Phi = NaN
    Make/O/D/N=(nB) w_phiRaw_Usadel1D_Phi = NaN
    Make/O/D/N=(nB) w_phiWrapped_Usadel1D_Phi = NaN
    Make/O/D/N=(nB) w_phiFold_Usadel1D_Phi = NaN
    Make/O/D/N=(nB) w_calcOrder_Usadel1D_Phi = NaN
    Make/O/D/N=(nB) w_solved_Usadel1D_Phi = 0

    Variable/G v_eta_Usadel1D_Phi = eta_eV
    Variable/G v_energySign_Usadel1D_Phi = signE

    Variable/G v_L_Usadel1D_Phi_nm = Ldiff_nm
    Variable/G v_Ldiff_Usadel1D_Phi_nm = Ldiff_nm
    Variable/G v_Lphase_Usadel1D_Phi_nm = LphaseUse_nm

    Variable/G v_xProbe_Usadel1D_Phi_nm = xProbe_nm
    Variable/G v_D_Usadel1D_Phi_nm2_s = D_nm2_s
    Variable/G v_singleBMode_Usadel1D_Phi = useSingleB
    Variable/G v_symmetryFolded_Usadel1D_Phi = 1

    Variable coeff_eVnm2 = 0.5 * HBAR_eVs * D_nm2_s

//    Variable maxNewton = 50
//    Variable resTol_eV = 2e-7
//    Variable stepTol = 1e-6
//    Variable damp = 1.0
//testing resolution
 	 Variable maxNewton = 30
	 Variable resTol_eV = 1e-6
	 Variable stepTol = 7e-5
	 Variable damp = 1.0

    Variable iZero = floor(NE/2)
    Variable iB, iStart, iNext, nDone
    Variable Bval, phiRaw, phiWrapped, phiFold
    Variable bestPhi, bestDiff

    for (iB=0; iB<nB; iB+=1)

        if (useSingleB)
            Bval = singleB_T
        else
            Bval = Bwave_T[iB]
        endif

        w_B_Usadel1D_Phi[iB] = Bval

        phiRaw = (2*q_e/HBAR_SI) * Bval * lambdaL_m * (LphaseUse_nm*1e-9)
        phiWrapped = atan2(sin(phiRaw), cos(phiRaw))
        phiFold = abs(phiWrapped)

        w_phiRaw_Usadel1D_Phi[iB] = phiRaw
        w_phiWrapped_Usadel1D_Phi[iB] = phiWrapped
        w_phiFold_Usadel1D_Phi[iB] = phiFold
    endfor

    bestPhi = Inf
    iStart = 0

    for (iB=0; iB<nB; iB+=1)
        if (w_phiFold_Usadel1D_Phi[iB] < bestPhi)
            bestPhi = w_phiFold_Usadel1D_Phi[iB]
            iStart = iB
        endif
    endfor

    SNS_Usadel1D_Phi_SolveOneB_GenericSeed(DOS_EB_local, E_axis_local, iStart, iZero, iProbe, theta_Phi, chi_Phi, thetaSeed_Usadel1D_Phi, chiSeed_Usadel1D_Phi, dx_nm, coeff_eVnm2, w_phiFold_Usadel1D_Phi[iStart], eta_eV, Delta, signE, maxNewton, resTol_eV, stepTol, damp)

    w_solved_Usadel1D_Phi[iStart] = 1
    w_calcOrder_Usadel1D_Phi[iStart] = 0
    nDone = 1

    do
        if (nDone >= nB)
            break
        endif

        bestDiff = Inf
        iNext = -1

        for (iB=0; iB<nB; iB+=1)
            if (w_solved_Usadel1D_Phi[iB] < 0.5)
                if (w_phiFold_Usadel1D_Phi[iB] < bestDiff)
                    bestDiff = w_phiFold_Usadel1D_Phi[iB]
                    iNext = iB
                endif
            endif
        endfor

        if (iNext < 0)
            break
        endif

        SNS_Usadel1D_Phi_SolveOneB_FromStoredSeed(DOS_EB_local, E_axis_local, iNext, iZero, iProbe, theta_Phi, chi_Phi, thetaSeed_Usadel1D_Phi, chiSeed_Usadel1D_Phi, dx_nm, coeff_eVnm2, w_phiFold_Usadel1D_Phi[iNext], eta_eV, Delta, signE, maxNewton, resTol_eV, stepTol, damp)

        w_solved_Usadel1D_Phi[iNext] = 1
        w_calcOrder_Usadel1D_Phi[iNext] = nDone
        nDone += 1

    while (1)

    String symName = outDOS + "_sym"
    Duplicate/O DOS_EB_local, $symName

    // -------------------------------------------------------------------------
    // Move diagnostics into :Usadel1D:, leaving only result waves in caller.
    // Result waves left in caller:
    //   outDOS
    //   outDOS+"_sym"
    //   outEaxis
    // -------------------------------------------------------------------------

    NewDataFolder/O :Usadel1D

    Duplicate/O theta_Phi, :Usadel1D:theta_Phi
    Duplicate/O chi_Phi, :Usadel1D:chi_Phi

    Duplicate/O thetaSeed_Usadel1D_Phi, :Usadel1D:thetaSeed_Usadel1D_Phi
    Duplicate/O chiSeed_Usadel1D_Phi, :Usadel1D:chiSeed_Usadel1D_Phi

    Duplicate/O w_iterCount_Usadel1D_Phi, :Usadel1D:w_iterCount_Usadel1D_Phi
    Duplicate/O w_converged_Usadel1D_Phi, :Usadel1D:w_converged_Usadel1D_Phi
    Duplicate/O w_resMax_Usadel1D_Phi, :Usadel1D:w_resMax_Usadel1D_Phi
    Duplicate/O w_stepMax_Usadel1D_Phi, :Usadel1D:w_stepMax_Usadel1D_Phi

    Duplicate/O w_B_Usadel1D_Phi, :Usadel1D:w_B_Usadel1D_Phi
    Duplicate/O w_phiRaw_Usadel1D_Phi, :Usadel1D:w_phiRaw_Usadel1D_Phi
    Duplicate/O w_phiWrapped_Usadel1D_Phi, :Usadel1D:w_phiWrapped_Usadel1D_Phi
    Duplicate/O w_phiFold_Usadel1D_Phi, :Usadel1D:w_phiFold_Usadel1D_Phi
    Duplicate/O w_calcOrder_Usadel1D_Phi, :Usadel1D:w_calcOrder_Usadel1D_Phi
    Duplicate/O w_solved_Usadel1D_Phi, :Usadel1D:w_solved_Usadel1D_Phi

    Variable/G :Usadel1D:v_eta_Usadel1D_Phi = v_eta_Usadel1D_Phi
    Variable/G :Usadel1D:v_energySign_Usadel1D_Phi = v_energySign_Usadel1D_Phi

    Variable/G :Usadel1D:v_L_Usadel1D_Phi_nm = v_L_Usadel1D_Phi_nm
    Variable/G :Usadel1D:v_Ldiff_Usadel1D_Phi_nm = v_Ldiff_Usadel1D_Phi_nm
    Variable/G :Usadel1D:v_Lphase_Usadel1D_Phi_nm = v_Lphase_Usadel1D_Phi_nm

    Variable/G :Usadel1D:v_xProbe_Usadel1D_Phi_nm = v_xProbe_Usadel1D_Phi_nm
    Variable/G :Usadel1D:v_D_Usadel1D_Phi_nm2_s = v_D_Usadel1D_Phi_nm2_s
    Variable/G :Usadel1D:v_singleBMode_Usadel1D_Phi = v_singleBMode_Usadel1D_Phi
    Variable/G :Usadel1D:v_symmetryFolded_Usadel1D_Phi = v_symmetryFolded_Usadel1D_Phi

    KillWaves/Z theta_Phi, chi_Phi
    KillWaves/Z thetaSeed_Usadel1D_Phi, chiSeed_Usadel1D_Phi
    KillWaves/Z w_iterCount_Usadel1D_Phi, w_converged_Usadel1D_Phi, w_resMax_Usadel1D_Phi, w_stepMax_Usadel1D_Phi
    KillWaves/Z w_B_Usadel1D_Phi, w_phiRaw_Usadel1D_Phi, w_phiWrapped_Usadel1D_Phi, w_phiFold_Usadel1D_Phi, w_calcOrder_Usadel1D_Phi, w_solved_Usadel1D_Phi

    KillVariables/Z v_eta_Usadel1D_Phi, v_energySign_Usadel1D_Phi
    KillVariables/Z v_L_Usadel1D_Phi_nm, v_Ldiff_Usadel1D_Phi_nm, v_Lphase_Usadel1D_Phi_nm
    KillVariables/Z v_xProbe_Usadel1D_Phi_nm, v_D_Usadel1D_Phi_nm2_s, v_singleBMode_Usadel1D_Phi, v_symmetryFolded_Usadel1D_Phi

    return 0
End


//==============================================================================
// SNS_Usadel1D_Phi_SolveOneB_GenericSeed
//
// Purpose:
//   Solve one field from E=0 upward using generic seed.
//   Stores theta/chi solution at each positive energy for continuation.
//==============================================================================

Function SNS_Usadel1D_Phi_SolveOneB_GenericSeed(DOS_EB, E_axis, iB, iZero, iProbe, theta, chi, thetaSeed, chiSeed, dx_nm, coeff, phiUse, eta_eV, Delta, signE, maxNewton, resTol_eV, stepTol, damp)
	Wave DOS_EB, E_axis
	Variable iB, iZero, iProbe
	Wave/C theta, chi, thetaSeed, chiSeed
	Variable dx_nm, coeff, phiUse, eta_eV, Delta, signE
	Variable maxNewton, resTol_eV, stepTol, damp

	Wave w_iterCount_Usadel1D_Phi
	Wave w_converged_Usadel1D_Phi
	Wave w_resMax_Usadel1D_Phi
	Wave w_stepMax_Usadel1D_Phi

	Variable NE = DimSize(E_axis,0)
	Variable Ngrid = DimSize(theta,0)
	Variable ie, ieMirror, ix
	Variable iterUsed, convFlag, resMax, stepMax, Nloc
	Complex thetaS0

	thetaS0 = SNS_Usadel1D_Phi_BCS_Theta_PositiveE(0, eta_eV, Delta)
	theta = thetaS0
	SNS_Usadel1D_Phi_SetLinearChi(chi, phiUse)

	for (ie=iZero; ie<NE; ie+=1)

		SNS_Usadel1D_Phi_SolveOneEnergyNewton(theta, chi, dx_nm, coeff, phiUse, E_axis[ie], eta_eV, Delta, signE, maxNewton, resTol_eV, stepTol, damp, iterUsed, convFlag, resMax, stepMax)

		Nloc = real(cos(theta[iProbe]))
		if (Nloc < 0)
			Nloc = 0
		endif

		DOS_EB[ie][iB] = Nloc
		w_iterCount_Usadel1D_Phi[ie][iB] = iterUsed
		w_converged_Usadel1D_Phi[ie][iB] = convFlag
		w_resMax_Usadel1D_Phi[ie][iB] = resMax
		w_stepMax_Usadel1D_Phi[ie][iB] = stepMax

		for (ix=0; ix<Ngrid; ix+=1)
			thetaSeed[ie][ix] = theta[ix]
			chiSeed[ie][ix] = chi[ix]
		endfor

		ieMirror = NE - 1 - ie
		if (ieMirror >= 0 && ieMirror < iZero)
			DOS_EB[ieMirror][iB] = DOS_EB[ie][iB]
			w_iterCount_Usadel1D_Phi[ieMirror][iB] = iterUsed
			w_converged_Usadel1D_Phi[ieMirror][iB] = convFlag
			w_resMax_Usadel1D_Phi[ieMirror][iB] = resMax
			w_stepMax_Usadel1D_Phi[ieMirror][iB] = stepMax
		endif
	endfor
End


//==============================================================================
// SNS_Usadel1D_Phi_SolveOneB_FromStoredSeed
//
// Purpose:
//   Solve one field using previous folded-phase solution as seed.
//   Updates seed waves to this new folded phase.
//==============================================================================

Function SNS_Usadel1D_Phi_SolveOneB_FromStoredSeed(DOS_EB, E_axis, iB, iZero, iProbe, theta, chi, thetaSeed, chiSeed, dx_nm, coeff, phiUse, eta_eV, Delta, signE, maxNewton, resTol_eV, stepTol, damp)
	Wave DOS_EB, E_axis
	Variable iB, iZero, iProbe
	Wave/C theta, chi, thetaSeed, chiSeed
	Variable dx_nm, coeff, phiUse, eta_eV, Delta, signE
	Variable maxNewton, resTol_eV, stepTol, damp

	Wave w_iterCount_Usadel1D_Phi
	Wave w_converged_Usadel1D_Phi
	Wave w_resMax_Usadel1D_Phi
	Wave w_stepMax_Usadel1D_Phi

	Variable NE = DimSize(E_axis,0)
	Variable Ngrid = DimSize(theta,0)
	Variable ie, ieMirror, ix
	Variable iterUsed, convFlag, resMax, stepMax, Nloc

	for (ie=iZero; ie<NE; ie+=1)

		for (ix=0; ix<Ngrid; ix+=1)
			theta[ix] = thetaSeed[ie][ix]
			chi[ix] = chiSeed[ie][ix]
		endfor

		chi[0] = cmplx(-0.5*phiUse, 0)
		chi[Ngrid-1] = cmplx(0.5*phiUse, 0)

		SNS_Usadel1D_Phi_SolveOneEnergyNewton(theta, chi, dx_nm, coeff, phiUse, E_axis[ie], eta_eV, Delta, signE, maxNewton, resTol_eV, stepTol, damp, iterUsed, convFlag, resMax, stepMax)

		Nloc = real(cos(theta[iProbe]))
		if (Nloc < 0)
			Nloc = 0
		endif

		DOS_EB[ie][iB] = Nloc
		w_iterCount_Usadel1D_Phi[ie][iB] = iterUsed
		w_converged_Usadel1D_Phi[ie][iB] = convFlag
		w_resMax_Usadel1D_Phi[ie][iB] = resMax
		w_stepMax_Usadel1D_Phi[ie][iB] = stepMax

		for (ix=0; ix<Ngrid; ix+=1)
			thetaSeed[ie][ix] = theta[ix]
			chiSeed[ie][ix] = chi[ix]
		endfor

		ieMirror = NE - 1 - ie
		if (ieMirror >= 0 && ieMirror < iZero)
			DOS_EB[ieMirror][iB] = DOS_EB[ie][iB]
			w_iterCount_Usadel1D_Phi[ieMirror][iB] = iterUsed
			w_converged_Usadel1D_Phi[ieMirror][iB] = convFlag
			w_resMax_Usadel1D_Phi[ieMirror][iB] = resMax
			w_stepMax_Usadel1D_Phi[ieMirror][iB] = stepMax
		endif
	endfor
End

//==============================================================================
// SNS_Usadel1D_Phi_SolveOneEnergyNewton
//
// Unknowns per interior point:
//   Re(theta), Im(theta), Re(chi), Im(chi)
//
// Residual per interior point:
//   Rtheta complex
//   Rchi complex
//==============================================================================

Function SNS_Usadel1D_Phi_SolveOneEnergyNewton(theta, chi, dx_nm, coeff, phi, E_eV, eta_eV, Delta_eV, signE, maxNewton, resTol_eV, stepTol, damp, iterUsed, convFlag, resMax, stepMax)
	Wave/C theta, chi
	Variable dx_nm, coeff, phi, E_eV, eta_eV, Delta_eV, signE
	Variable maxNewton, resTol_eV, stepTol, damp
	Variable &iterUsed, &convFlag, &resMax, &stepMax

	Variable n = DimSize(theta,0)
	Variable nInt = n - 2
	Variable nUnknown = 4 * nInt

	Variable iter, k, i, idx
	Variable scale

	Complex thetaS = SNS_Usadel1D_Phi_BCS_Theta_PositiveE(E_eV, eta_eV, Delta_eV)

	Make/O/D/N=(nUnknown) w_Phi_R, w_Phi_b, w_Phi_delta
	Make/O/D/N=(nUnknown,nUnknown) M_Phi_J

	theta[0] = thetaS
	theta[n-1] = thetaS
	chi[0] = cmplx(-0.5*phi, 0)
	chi[n-1] = cmplx(0.5*phi, 0)

	convFlag = 0
	iterUsed = maxNewton
	resMax = NaN
	stepMax = NaN

	for (iter=0; iter<maxNewton; iter+=1)

		theta[0] = thetaS
		theta[n-1] = thetaS
		chi[0] = cmplx(-0.5*phi, 0)
		chi[n-1] = cmplx(0.5*phi, 0)

		M_Phi_J = 0
		w_Phi_R = 0

		SNS_Usadel1D_Phi_BuildResidualAndJac(theta, chi, w_Phi_R, M_Phi_J, dx_nm, coeff, E_eV, eta_eV, Delta_eV, signE)

		resMax = SNS_Usadel1D_Phi_MaxAbsRealWave(w_Phi_R)

		if (resMax < resTol_eV)
			iterUsed = iter + 1
			convFlag = 1
			stepMax = 0
			break
		endif

		w_Phi_b = -w_Phi_R
		MatrixLinearSolve/O/Z M_Phi_J, w_Phi_b

		if (V_flag != 0)
			convFlag = 0
			iterUsed = iter + 1
			break
		endif

		w_Phi_delta = w_Phi_b
		stepMax = SNS_Usadel1D_Phi_MaxAbsRealWave(w_Phi_delta)

		scale = damp
		if (stepMax > 0.5)
			scale = damp * 0.5 / stepMax
		endif

		for (k=0; k<nInt; k+=1)
			i = k + 1
			idx = 4*k

			theta[i] += scale * cmplx(w_Phi_delta[idx],   w_Phi_delta[idx+1])
			chi[i]   += scale * cmplx(w_Phi_delta[idx+2], w_Phi_delta[idx+3])
		endfor

		if (stepMax*scale < stepTol)
			theta[0] = thetaS
			theta[n-1] = thetaS
			chi[0] = cmplx(-0.5*phi, 0)
			chi[n-1] = cmplx(0.5*phi, 0)

			SNS_Usadel1D_Phi_BuildResidual(theta, chi, w_Phi_R, dx_nm, coeff, E_eV, eta_eV, Delta_eV, signE)
			resMax = SNS_Usadel1D_Phi_MaxAbsRealWave(w_Phi_R)

			if (resMax < resTol_eV)
				iterUsed = iter + 1
				convFlag = 1
				break
			endif
		endif
	endfor

	theta[0] = thetaS
	theta[n-1] = thetaS
	chi[0] = cmplx(-0.5*phi, 0)
	chi[n-1] = cmplx(0.5*phi, 0)

	SNS_Usadel1D_Phi_BuildResidual(theta, chi, w_Phi_R, dx_nm, coeff, E_eV, eta_eV, Delta_eV, signE)
	resMax = SNS_Usadel1D_Phi_MaxAbsRealWave(w_Phi_R)

	if (resMax < resTol_eV)
		convFlag = 1
	endif
End


//==============================================================================
// SNS_Usadel1D_Phi_BuildResidualAndJac
//==============================================================================

Function SNS_Usadel1D_Phi_BuildResidualAndJac(theta, chi, Rvec, Jmat, dx_nm, coeff, E_eV, eta_eV, Delta_eV, signE)
	Wave/C theta, chi
	Wave Rvec, Jmat
	Variable dx_nm, coeff, E_eV, eta_eV, Delta_eV, signE

	Variable n = DimSize(theta,0)
	Variable nInt = n - 2
	Variable dx2 = dx_nm * dx_nm
	Variable k, i, idx, col

	Complex zE = cmplx(E_eV, eta_eV)

	Complex th, thM, thP
	Complex chM, ch0, chP
	Complex s, c, sM, cM, sP, cP
	Complex q, lap
	Complex Rtheta, Rchi
	Complex Aplus, Aminus, dAi, dAp, dAm
	Complex dTheta_i, dTheta_m, dTheta_p
	Complex dChi_m, dChi_i, dChi_p
	Complex dp, dm

	for (k=0; k<nInt; k+=1)

		i = k + 1
		idx = 4*k

		th  = theta[i]
		thM = theta[i-1]
		thP = theta[i+1]

		ch0 = chi[i]
		chM = chi[i-1]
		chP = chi[i+1]

		s  = sin(th)
		c  = cos(th)
		sM = sin(thM)
		cM = cos(thM)
		sP = sin(thP)
		cP = cos(thP)

		q = (chP - chM) / (2*dx_nm)
		lap = (thM - 2*th + thP) / dx2

		Rtheta = coeff * (lap - q*q*s*c) + signE * cmplx(0,1)*zE*s

		dp = chP - ch0
		dm = ch0 - chM

		Aplus  = 0.5 * (s*s + sP*sP)
		Aminus = 0.5 * (sM*sM + s*s)

		Rchi = coeff * (Aplus*dp - Aminus*dm) / dx2

		Rvec[idx]   = real(Rtheta)
		Rvec[idx+1] = imag(Rtheta)
		Rvec[idx+2] = real(Rchi)
		Rvec[idx+3] = imag(Rchi)

		// Rtheta Jacobian
		dTheta_i = coeff * (-2/dx2 - q*q*cos(2*th)) + signE * cmplx(0,1)*zE*c
		dTheta_m = coeff / dx2
		dTheta_p = coeff / dx2

		dChi_m =  coeff * q*s*c / dx_nm
		dChi_p = -coeff * q*s*c / dx_nm

		SNS_Usadel1D_Phi_AddComplexBlock(Jmat, idx, idx, dTheta_i)

		if (i > 1)
			col = 4*(k-1)
			SNS_Usadel1D_Phi_AddComplexBlock(Jmat, idx, col, dTheta_m)
		endif

		if (i < n-2)
			col = 4*(k+1)
			SNS_Usadel1D_Phi_AddComplexBlock(Jmat, idx, col, dTheta_p)
		endif

		if (i > 1)
			col = 4*(k-1) + 2
			SNS_Usadel1D_Phi_AddComplexBlock(Jmat, idx, col, dChi_m)
		endif

		if (i < n-2)
			col = 4*(k+1) + 2
			SNS_Usadel1D_Phi_AddComplexBlock(Jmat, idx, col, dChi_p)
		endif

		// Rchi Jacobian
		dAi = s*c
		dAp = sP*cP
		dAm = sM*cM

		dTheta_i = coeff * dAi * (dp - dm) / dx2
		dTheta_p = coeff * dAp * dp / dx2
		dTheta_m = -coeff * dAm * dm / dx2

		dChi_i = -coeff * (Aplus + Aminus) / dx2
		dChi_p =  coeff * Aplus / dx2
		dChi_m =  coeff * Aminus / dx2

		SNS_Usadel1D_Phi_AddComplexBlock(Jmat, idx+2, idx, dTheta_i)

		if (i > 1)
			col = 4*(k-1)
			SNS_Usadel1D_Phi_AddComplexBlock(Jmat, idx+2, col, dTheta_m)
		endif

		if (i < n-2)
			col = 4*(k+1)
			SNS_Usadel1D_Phi_AddComplexBlock(Jmat, idx+2, col, dTheta_p)
		endif

		SNS_Usadel1D_Phi_AddComplexBlock(Jmat, idx+2, idx+2, dChi_i)

		if (i > 1)
			col = 4*(k-1) + 2
			SNS_Usadel1D_Phi_AddComplexBlock(Jmat, idx+2, col, dChi_m)
		endif

		if (i < n-2)
			col = 4*(k+1) + 2
			SNS_Usadel1D_Phi_AddComplexBlock(Jmat, idx+2, col, dChi_p)
		endif
	endfor
End


//==============================================================================
// SNS_Usadel1D_Phi_BuildResidual
//==============================================================================

Function SNS_Usadel1D_Phi_BuildResidual(theta, chi, Rvec, dx_nm, coeff, E_eV, eta_eV, Delta_eV, signE)
	Wave/C theta, chi
	Wave Rvec
	Variable dx_nm, coeff, E_eV, eta_eV, Delta_eV, signE

	Variable n = DimSize(theta,0)
	Variable nInt = n - 2
	Variable dx2 = dx_nm * dx_nm
	Variable k, i, idx

	Complex zE = cmplx(E_eV, eta_eV)

	Complex th, thM, thP
	Complex chM, ch0, chP
	Complex s, c, sM, sP
	Complex q, lap
	Complex Rtheta, Rchi
	Complex Aplus, Aminus
	Complex dp, dm

	for (k=0; k<nInt; k+=1)

		i = k + 1
		idx = 4*k

		th  = theta[i]
		thM = theta[i-1]
		thP = theta[i+1]

		ch0 = chi[i]
		chM = chi[i-1]
		chP = chi[i+1]

		s  = sin(th)
		c  = cos(th)
		sM = sin(thM)
		sP = sin(thP)

		q = (chP - chM) / (2*dx_nm)
		lap = (thM - 2*th + thP) / dx2

		Rtheta = coeff * (lap - q*q*s*c) + signE * cmplx(0,1)*zE*s

		dp = chP - ch0
		dm = ch0 - chM

		Aplus  = 0.5 * (s*s + sP*sP)
		Aminus = 0.5 * (sM*sM + s*s)

		Rchi = coeff * (Aplus*dp - Aminus*dm) / dx2

		Rvec[idx]   = real(Rtheta)
		Rvec[idx+1] = imag(Rtheta)
		Rvec[idx+2] = real(Rchi)
		Rvec[idx+3] = imag(Rchi)
	endfor
End


//==============================================================================
// SNS_Usadel1D_Phi_AddComplexBlock
//==============================================================================

Function SNS_Usadel1D_Phi_AddComplexBlock(M, row, col, a)
	Wave M
	Variable row, col
	Complex a

	M[row][col]     += real(a)
	M[row][col+1]   += -imag(a)
	M[row+1][col]   += imag(a)
	M[row+1][col+1] += real(a)
End


//==============================================================================
// SNS_Usadel1D_Phi_SetLinearChi
//==============================================================================

Function SNS_Usadel1D_Phi_SetLinearChi(chi, phi)
	Wave/C chi
	Variable phi

	Variable n = DimSize(chi,0)
	Variable i

	for (i=0; i<n; i+=1)
		chi[i] = cmplx(-0.5*phi + phi*i/(n-1), 0)
	endfor
End


//==============================================================================
// SNS_Usadel1D_Phi_BCS_Theta_PositiveE
//==============================================================================

Function/C SNS_Usadel1D_Phi_BCS_Theta_PositiveE(E_eV, eta_eV, Delta_eV)
	Variable E_eV, eta_eV, Delta_eV

	Complex zE = cmplx(E_eV, eta_eV)
	Complex rootBCS = sqrt(zE*zE - cmplx(Delta_eV^2,0))
	Complex gS = zE / rootBCS

	if (E_eV > Delta_eV && real(gS) < 0)
		rootBCS = -rootBCS
		gS = zE / rootBCS
	endif

	if (E_eV >= 0 && real(gS) < 0)
		gS = -gS
	endif

	return acos(gS)
End


//==============================================================================
// SNS_Usadel1D_Phi_MaxAbsRealWave
//==============================================================================

Function SNS_Usadel1D_Phi_MaxAbsRealWave(w)
	Wave w

	Variable n = DimSize(w,0)
	Variable i
	Variable m = 0

	for (i=0; i<n; i+=1)
		m = max(m, abs(w[i]))
	endfor

	return m
End


//==============================================================================
// SNS_Usadel1D_PhaseBiasNewton_XRange_OneB_FromSettings
//
// Purpose:
//   Distance-dependence wrapper for the 1D phase-biased Usadel solver.
//   Varies xProbe_nm for one fixed B value.
//
// Calls:
//   SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings(..., singleB_T=B_T_single)
//
// Input:
//   L_nm          : 1D wire length [nm]
//   xProbeWave_nm : 1D wave of probe positions [nm]
//   B_T_single    : fixed field [T]
//   outDOS        : output DOS(E,xProbe)
//   outEaxis      : output energy axis
//   D_nm2_s       : diffusion constant [nm^2/s]
//   Ngrid         : spatial grid for the Usadel solver
//
// Optional:
//   etaOverride_eV : overrides root:SNS_Settings:Broadening
//   energySign     : +1 default, -1 alternate convention
//
// Outputs:
//   outDOS[energy][xProbe]
//   outEaxis[energy]
//   w_xProbe_Usadel1D_X_nm
//   w_resMax_Usadel1D_X
//   w_converged_Usadel1D_X
//   w_iterCount_Usadel1D_X
//==============================================================================

Function SNS_Usadel1D_PhaseBiasNewton_XRange_OneB_FromSettings(L_nm, xProbeWave_nm, B_T_single, outDOS, outEaxis, D_nm2_s, Ngrid, [etaOverride_eV, energySign])
	Variable L_nm, B_T_single, D_nm2_s, Ngrid
	Wave xProbeWave_nm
	String outDOS, outEaxis
	Variable etaOverride_eV, energySign

	NVAR NE_cfg = root:SNS_Settings:NE
	NVAR Delta  = root:SNS_Settings:Delta

	Variable signE = 1
	if (!ParamIsDefault(energySign))
		if (energySign < 0)
			signE = -1
		endif
	endif

	Variable NE = round(NE_cfg)
	if (mod(NE, 2) == 0)
		NE += 1
	endif

	Variable nX = DimSize(xProbeWave_nm, 0)
	if (nX <= 0)
		Abort "SNS_Usadel1D_PhaseBiasNewton_XRange_OneB_FromSettings: xProbeWave_nm is empty."
	endif

	Variable Emin_eV = -5 * Delta
	Variable Emax_eV =  5 * Delta

	Make/O/D/N=(NE) $outEaxis
	Wave E_axis_out = $outEaxis
	SetScale/I x, Emin_eV, Emax_eV, "eV", E_axis_out
	E_axis_out = x

	Make/O/D/N=(NE,nX) $outDOS
	Wave DOS_EX = $outDOS
	DOS_EX = NaN
	SetScale/I x, Emin_eV, Emax_eV, "eV", DOS_EX

	if (nX > 1)
		SetScale/I y, xProbeWave_nm[0], xProbeWave_nm[nX-1], "nm", DOS_EX
	else
		SetScale/P y, xProbeWave_nm[0], 1, "nm", DOS_EX
	endif

	Duplicate/O xProbeWave_nm, w_xProbe_Usadel1D_X_nm

	Make/O/D/N=(NE,nX) w_resMax_Usadel1D_X = NaN
	Make/O/D/N=(NE,nX) w_converged_Usadel1D_X = 0
	Make/O/D/N=(NE,nX) w_iterCount_Usadel1D_X = NaN
	Make/O/D/N=(NE,nX) w_stepMax_Usadel1D_X = NaN
	
	Wave w_resMax_Usadel1D_Phi
	Wave w_converged_Usadel1D_Phi
	Wave w_iterCount_Usadel1D_Phi
	Wave w_stepMax_Usadel1D_Phi

	String tmpDOS = "__tmp_DOS_Usadel1D_X"
	String tmpE   = "__tmp_E_Usadel1D_X"

	Variable ixp, ie
	Variable xProbe_nm

	for (ixp=0; ixp<nX; ixp+=1)

		xProbe_nm = xProbeWave_nm[ixp]

		if (xProbe_nm < 0 || xProbe_nm > L_nm)
			Abort "SNS_Usadel1D_PhaseBiasNewton_XRange_OneB_FromSettings: xProbe value outside [0,L_nm]."
		endif

		if (!ParamIsDefault(etaOverride_eV))
			SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings(L_nm, xProbe_nm, tmpDOS, tmpE, D_nm2_s, Ngrid, etaOverride_eV=etaOverride_eV, energySign=signE, singleB_T=B_T_single)
		else
			SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings(L_nm, xProbe_nm, tmpDOS, tmpE, D_nm2_s, Ngrid, energySign=signE, singleB_T=B_T_single)
		endif

		Wave tmpDOSW = $tmpDOS
		Wave tmpEW   = $tmpE

		for (ie=0; ie<NE; ie+=1)
			DOS_EX[ie][ixp] = tmpDOSW[ie][0]
			E_axis_out[ie] = tmpEW[ie]

			w_resMax_Usadel1D_X[ie][ixp] = w_resMax_Usadel1D_Phi[ie][0]
			w_converged_Usadel1D_X[ie][ixp] = w_converged_Usadel1D_Phi[ie][0]
			w_iterCount_Usadel1D_X[ie][ixp] = w_iterCount_Usadel1D_Phi[ie][0]
			w_stepMax_Usadel1D_X[ie][ixp] = w_stepMax_Usadel1D_Phi[ie][0]
		endfor
	endfor

	KillWaves/Z $tmpDOS, $tmpE
End


//==============================================================================
// SNS_Usadel1D_FromMaskFolder
//
// Purpose:
//   One-call wrapper for the effective 1D Usadel calculation from a 2D mask
//   folder and one STS position.
//
//   Required workflow performed internally:
//     1) Build ballistic ray ensemble at STSx/STSy for the specified B angle.
//     2) Compress ray geometry to Usadel 1D parameters:
//          Ldiff_nm, xProbe_nm
//     3) Compress ray magnetic lever arm:
//          Lphi_use_nm
//     4) Estimate effective mean free path from bulk + edge scattering:
//
//          1/l_eff = 1/l_bulk + 1/l_edge
//
//        with
//
//          l_edge = edgePrefactor * Area / Perimeter
//
//        using SNS_MaskAreaPerim_FromParticles(w_mask).
//
//     5) Compute
//
//          D = vF * l_eff / 2
//
//        in nm^2/s and run SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings.
//
// Required inputs:
//   maskFolder : folder containing w_mask and where ray outputs are stored.
//                Can be supplied with or without trailing colon.
//   STSx       : STS x position [nm].
//   STSy       : STS y position [nm].
//   Bangle_deg : in-plane B direction [deg]. Required because W_eff_List and
//                Lphi depend on the magnetic-field direction.
//   outDOS     : output DOS wave name created in caller folder.
//   outEaxis   : output energy-axis wave name created in caller folder.
//
// SNS settings loaded through SNS_Params:
//   params.vF                  : Fermi velocity [m/s].
//   params.lmfp_bulk_nm        : bulk / impurity mean free path [nm].
//   params.lmfp_edge_prefactor : edge mfp prefactor.
//   params.BTK_barrier         : default BTK barrier if optional input omitted.
//
// Optional inputs:
//   BTK_barrier : optional override. Default: params.BTK_barrier.
//   Ngrid       : default 101.
//   useWeights  : default 0. Passed to geometry and lever-arm helpers.
//                 0 = equal channel weights.
//                 1 = weight by max(T_eff_List,0), if available.
//   useMax      : default 0. Passed to lever-arm helper.
//                 0 = use Lphi_rms.
//                 1 = use Lphi_max.
//   doDisplay   : default 0. Passed to SNS_ExtractModesForFolder.
//   doPrint     : default 0. Passed to Usadel geometry/lever-arm helpers.
//   energySign  : default 1. Passed to Usadel solver.
//
// Outputs in caller folder:
//   outDOS
//   outDOS+"_sym"
//   outEaxis
//
// Diagnostics stored in maskFolder:
//   v_Usadel_area_nm2
//   v_Usadel_perim_nm
//   v_Usadel_lmfp_bulk_nm
//   v_Usadel_lmfp_edge_nm
//   v_Usadel_lmfp_eff_nm
//   v_Usadel_lmfp_edge_prefactor
//   v_Usadel_D_nm2_s
//   v_Usadel_D_m2_s
//   v_Usadel_D_cm2_s
//   v_Usadel_Bangle_deg
//   v_Usadel_BTK_barrier
//
//==============================================================================

Function SNS_Usadel1D_FromMaskFolder(maskFolder, STSx, STSy, Bangle_deg, outDOS, outEaxis, [BTK_barrier, Ngrid, useWeights, useMax, doDisplay, doPrint, energySign])
    String maskFolder
    Variable STSx, STSy, Bangle_deg
    String outDOS, outEaxis

    Variable BTK_barrier, Ngrid, useWeights, useMax, doDisplay, doPrint, energySign

    String oldDF = GetDataFolder(1)

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    // ------------------------------
    // Defaults
    // ------------------------------
    if (ParamIsDefault(Ngrid))
        Ngrid = 101
    endif

    if (ParamIsDefault(useWeights))
        useWeights = 0
    endif

    if (ParamIsDefault(useMax))
        useMax = 0
    endif

    if (ParamIsDefault(doDisplay))
        doDisplay = 0
    endif

    if (ParamIsDefault(doPrint))
        doPrint = 0
    endif

    if (ParamIsDefault(energySign))
        energySign = 1
    endif

    if (ParamIsDefault(BTK_barrier))
        BTK_barrier = params.BTK_barrier
    endif

    if (numtype(BTK_barrier) != 0 || BTK_barrier < 0)
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_FromMaskFolder: BTK_barrier must be finite and non-negative."
    endif

    if (numtype(params.vF) != 0 || params.vF <= 0)
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_FromMaskFolder: invalid params.vF. Run SNS_InitDefaultSettings first."
    endif

    if (numtype(params.lmfp_bulk_nm) != 0 || params.lmfp_bulk_nm <= 0)
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_FromMaskFolder: invalid params.lmfp_bulk_nm."
    endif

    if (numtype(params.lmfp_edge_prefactor) != 0 || params.lmfp_edge_prefactor <= 0)
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_FromMaskFolder: invalid params.lmfp_edge_prefactor."
    endif

    // ------------------------------
    // Folder normalization
    // ------------------------------
    String maskFolderColon = maskFolder
    if (cmpstr(maskFolderColon[strlen(maskFolderColon)-1, strlen(maskFolderColon)-1], ":") != 0)
        maskFolderColon += ":"
    endif

    String maskFolderNoColon = maskFolderColon[0, strlen(maskFolderColon)-2]

    if (!DataFolderExists(maskFolderColon))
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_FromMaskFolder: maskFolder does not exist."
    endif

    Wave/Z Nmask = $(maskFolderColon + "w_mask")
    if (!WaveExists(Nmask))
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_FromMaskFolder: missing w_mask in maskFolder."
    endif

    // ------------------------------
    // Build ray ensemble and Usadel geometry
    // ------------------------------
    SNS_ExtractModesForFolder(maskFolderNoColon, STSx=STSx, STSy=STSy, Bangle_deg=Bangle_deg, doDisplay=doDisplay)

    SNS_Usadel1D_GeometryFromRayFolder(maskFolderColon, STSx=STSx, STSy=STSy, useWeights=useWeights, doPrint=doPrint)

    SNS_Usadel1D_LeverArmFromRayFolder(maskFolderColon, useWeights=useWeights, useMax=useMax, doPrint=doPrint)

    // ------------------------------
    // Compute edge + bulk mean free path
    // ------------------------------
    SetDataFolder $maskFolderColon

    Wave w_area_perim = SNS_MaskAreaPerim_FromParticles(Nmask)

    Variable area_nm2 = w_area_perim[0]
    Variable perim_nm = w_area_perim[1]

    if (numtype(area_nm2) != 0 || numtype(perim_nm) != 0 || area_nm2 <= 0 || perim_nm <= 0)
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_FromMaskFolder: invalid mask area/perimeter."
    endif

    Variable lmfpEdge_nm = params.lmfp_edge_prefactor * area_nm2 / perim_nm
    Variable lmfpEff_nm  = 1 / (1/params.lmfp_bulk_nm + 1/lmfpEdge_nm)

    Variable D_nm2_s = 0.5 * params.vF * lmfpEff_nm * 1e9
    Variable D_m2_s  = D_nm2_s * 1e-18
    Variable D_cm2_s = D_m2_s * 1e4

    Variable/G v_Usadel_area_nm2 = area_nm2
    Variable/G v_Usadel_perim_nm = perim_nm
    Variable/G v_Usadel_lmfp_bulk_nm = params.lmfp_bulk_nm
    Variable/G v_Usadel_lmfp_edge_nm = lmfpEdge_nm
    Variable/G v_Usadel_lmfp_eff_nm = lmfpEff_nm
    Variable/G v_Usadel_lmfp_edge_prefactor = params.lmfp_edge_prefactor
    Variable/G v_Usadel_D_nm2_s = D_nm2_s
    Variable/G v_Usadel_D_m2_s = D_m2_s
    Variable/G v_Usadel_D_cm2_s = D_cm2_s
    Variable/G v_Usadel_Bangle_deg = Bangle_deg
    Variable/G v_Usadel_BTK_barrier = BTK_barrier

    NVAR v_Usadel_Ldiff_nm = v_Usadel_Ldiff_nm
    NVAR v_Usadel_xProbe_nm = v_Usadel_xProbe_nm
    NVAR v_Usadel_Lphi_use_nm = v_Usadel_Lphi_use_nm

    // ------------------------------
    // Run Usadel solver in caller folder
    // ------------------------------
    SetDataFolder $oldDF

    SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings(v_Usadel_Ldiff_nm, v_Usadel_xProbe_nm, outDOS, outEaxis, D_nm2_s, Ngrid, energySign=energySign, Lphase_nm=v_Usadel_Lphi_use_nm)

    return 0
End



//==============================================================================
// SNS_FiniteMedianFromWave
//
// Purpose:
//   Return the median of finite entries in a 1D wave.
//   NaN/Inf entries are ignored.
//
// Returns:
//   median value, or NaN if no finite entries exist.
//==============================================================================
Function SNS_FiniteMedianFromWave(wIn)
    Wave wIn

    Variable n = numpnts(wIn)
    if (n <= 0)
        return NaN
    endif

    Make/FREE/D/N=(n) wTmp

    Variable i, nValid = 0
    for (i = 0; i < n; i += 1)
        if (numtype(wIn[i]) == 0)
            wTmp[nValid] = wIn[i]
            nValid += 1
        endif
    endfor

    if (nValid <= 0)
        return NaN
    endif

    Redimension/N=(nValid) wTmp
    Sort wTmp, wTmp

    if (mod(nValid, 2) == 1)
        return wTmp[floor(nValid/2)]
    else
        return 0.5 * (wTmp[nValid/2 - 1] + wTmp[nValid/2])
    endif
End


//==============================================================================
// SNS_Usadel1D_LineDOS_GlobalWireFromMask
//
// Purpose:
//   Map a ballistic STM LineSTS geometry onto one global effective 1D Usadel
//   wire and compute a line-resolved Usadel LDOS(E,r) at fixed B0.
//
//   Workflow:
//     1) Build STM line positions from xStart/yStart to xEnd/yEnd.
//     2) For each STM position:
//          - ray trace using SNS_ExtractModesForFolder(...)
//          - compute local Usadel geometry using
//              SNS_Usadel1D_GeometryFromRayFolder(...)
//          - compute local magnetic lever arm using
//              SNS_Usadel1D_LeverArmFromRayFolder(...)
//          - store local Ldiff, fProbe, xProbe, Lphi.
//     3) Define one global 1D Usadel wire:
//          Ldiff_global = median_s[Ldiff_local(s)]
//          Lphi_global  = median_s[Lphi_local(s)]
//     4) Map each STM position onto this global wire:
//          xProbe_global(s) = fProbe_local(s) * Ldiff_global
//     5) Compute Usadel LDOS at each xProbe_global(s), using the same
//        Ldiff_global, Lphi_global, D, and B0.
//
// Required inputs:
//   maskFolder : folder containing w_mask and image data for ray tracing.
//                Can be supplied with or without trailing colon.
//   xStart     : STM line start x [nm].
//   yStart     : STM line start y [nm].
//   xEnd       : STM line end x [nm].
//   yEnd       : STM line end y [nm].
//   B0         : fixed magnetic field [T].
//   Bangle_deg : in-plane B-field angle [deg].
//
// Optional inputs:
//   outDOS      : raw output DOS name.
//                 Default at B != 0:
//                     LDOS_Usadel1D_Global_B###mT_###deg
//                 Default at B == 0:
//                     LDOS_Usadel1D_Global_B0mT
//   outEaxis    : output energy-axis name.
//                 Default at B != 0:
//                     E_axis_Usadel1D_Global_B###mT_###deg
//                 Default at B == 0:
//                     E_axis_Usadel1D_Global_B0mT
//   BTK_barrier : default params.BTK_barrier.
//   Ngrid       : Usadel spatial grid. Default 101.
//   useWeights  : default 0. Passed to geometry/lever-arm helpers.
//                 0 = equal channel weights.
//                 1 = weight by max(T_eff_List,0), if available.
//   useMax      : default 0. Passed to lever-arm helper.
//                 0 = local Lphi uses RMS.
//                 1 = local Lphi uses max.
//   doDisplay   : default 0. Passed to SNS_ExtractModesForFolder.
//   doPrint     : default 0. Prints helper summaries.
//   energySign  : default 1. Passed to Usadel solver.
//   doBroadening: default 1. Also create LDOS_Conv_Usadel1D_Global_*.
//   NLinePts    : optional number of positions along the STM line.
//                 If supplied, this overrides the default LambdaF-based
//                 line discretization by using
//                     lineStep_nm = lineLength_nm / (NLinePts - 1)
//                 The value passed into SNS_BuildLinePositions(...) is
//                     lineStep_m = lineStep_nm * 1e-9
//                 Must be a finite integer >= 2.
//
// Outputs in caller folder:
//   LDOS_Usadel1D_Global_*
//   LDOS_Conv_Usadel1D_Global_*
//   E_axis_Usadel1D_Global_*
//   R_axis_Usadel1D_Global_*
//
// Diagnostics in caller folder:
//   Usadel1D_LineGlobal:
//      SNS_UsadelLine_X_nm
//      SNS_UsadelLine_Y_nm
//      R_axis_Usadel1D_Global_*
//      w_Ldiff_local_nm
//      w_Lphi_local_nm
//      w_fProbe_local
//      w_xProbe_local_nm
//      w_xProbe_global_nm
//      scalar v_UsadelLine_* diagnostics
//
// Notes:
//   Local compression remains RMS through the existing helper definitions.
//   The global wire uses median across STM positions to avoid outlier line
//   positions defining the whole effective 1D wire.
//
//==============================================================================
Function SNS_Usadel1D_LineDOS_GlobalWireFromMask(maskFolder, xStart, yStart, xEnd, yEnd, B0, Bangle_deg, [outDOS, outEaxis, BTK_barrier, Ngrid, useWeights, useMax, doDisplay, doPrint, energySign, doBroadening, NLinePts])
    String maskFolder
    Variable xStart, yStart, xEnd, yEnd
    Variable B0, Bangle_deg

    String outDOS, outEaxis
    Variable BTK_barrier, Ngrid, useWeights, useMax, doDisplay, doPrint, energySign, doBroadening
    Variable NLinePts

    DFREF dfrCaller = GetDataFolderDFR()
    String callerPath = GetDataFolder(1)

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    // ------------------------------
    // Defaults
    // ------------------------------
    if (ParamIsDefault(BTK_barrier))
        BTK_barrier = params.BTK_barrier
    endif
    if (ParamIsDefault(Ngrid))
        Ngrid = 101
    endif
    if (ParamIsDefault(useWeights))
        useWeights = 0
    endif
    if (ParamIsDefault(useMax))
        useMax = 0
    endif
    if (ParamIsDefault(doDisplay))
        doDisplay = 0
    endif
    if (ParamIsDefault(doPrint))
        doPrint = 0
    endif
    if (ParamIsDefault(energySign))
        energySign = 1
    endif
    if (ParamIsDefault(doBroadening))
        doBroadening = 1
    endif

    Variable isZeroField = (abs(B0) <= 1e-15)

    // ------------------------------
    // Default output names
    // ------------------------------
    String tag
    if (isZeroField)
        tag = "B0mT"
    else
        tag = "B" + num2str(round(B0*1e3)) + "mT_" + num2str(round(Bangle_deg)) + "deg"
    endif

    String nameDOSDefault     = "LDOS_Usadel1D_Global_"      + tag
    String nameDOSConvDefault = "LDOS_Conv_Usadel1D_Global_" + tag
    String nameEaxisDefault   = "E_axis_Usadel1D_Global_"    + tag
    String nameRaxisDefault   = "R_axis_Usadel1D_Global_"    + tag

    if (ParamIsDefault(outDOS))
        outDOS = nameDOSDefault
    endif

    if (ParamIsDefault(outEaxis))
        outEaxis = nameEaxisDefault
    endif

    String outDOSConv = nameDOSConvDefault
    String nameRaxis = nameRaxisDefault

    if (numtype(BTK_barrier) != 0 || BTK_barrier < 0)
        Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: BTK_barrier must be finite and non-negative."
    endif
    if (numtype(params.vF) != 0 || params.vF <= 0)
        Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: invalid params.vF. Run SNS_InitDefaultSettings first."
    endif
    if (numtype(params.lmfp_bulk_nm) != 0 || params.lmfp_bulk_nm <= 0)
        Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: invalid params.lmfp_bulk_nm."
    endif
    if (numtype(params.lmfp_edge_prefactor) != 0 || params.lmfp_edge_prefactor <= 0)
        Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: invalid params.lmfp_edge_prefactor."
    endif

    // ------------------------------
    // Folder normalization
    // ------------------------------
    String maskFolderColon = maskFolder
    if (cmpstr(maskFolderColon[strlen(maskFolderColon)-1, strlen(maskFolderColon)-1], ":") != 0)
        maskFolderColon += ":"
    endif

    String maskFolderNoColon = maskFolderColon[0, strlen(maskFolderColon)-2]

    if (!DataFolderExists(maskFolderColon))
        Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: maskFolder does not exist."
    endif

    Wave/Z Nmask = $(maskFolderColon + "w_mask")
    if (!WaveExists(Nmask))
        Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: missing w_mask in maskFolder."
    endif

    // ------------------------------
    // Diagnostics folder and line positions
    // ------------------------------
    String diagFolder = callerPath + "Usadel1D_LineGlobal"
    NewDataFolder/O $diagFolder

    SetDataFolder $diagFolder

    String nameX = "SNS_UsadelLine_X_nm"
    String nameY = "SNS_UsadelLine_Y_nm"

    Variable lineStep_nm = params.LambdaF

    if (!ParamIsDefault(NLinePts))
        NLinePts = round(NLinePts)

        if ((numtype(NLinePts) != 0) || (NLinePts < 2))
            SetDataFolder dfrCaller
            Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: NLinePts must be a finite integer >= 2."
        endif

        Variable lineLength_nm = sqrt((xEnd - xStart)^2 + (yEnd - yStart)^2)

        if ((numtype(lineLength_nm) != 0) || (lineLength_nm <= 0))
            SetDataFolder dfrCaller
            Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: cannot use NLinePts for a zero-length or invalid line."
        endif

        lineStep_nm = lineLength_nm / NLinePts
    endif

    Variable lineStep_m = lineStep_nm * 1e-9

    Variable Nr = SNS_BuildLinePositions(xStart, yStart, xEnd, yEnd, lineStep_m, nameX, nameY, nameRaxis)

    Wave Xline = $nameX
    Wave Yline = $nameY
    Wave Rline = $nameRaxis

    if (Nr <= 0)
        SetDataFolder dfrCaller
        Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: no line positions."
    endif

    if (Nr > 1)
        SetScale/P x, Rline[0], Rline[1]-Rline[0], "nm", Rline
    else
        SetScale/P x, 0, 1, "nm", Rline
    endif

    Make/O/D/N=(Nr) w_Ldiff_local_nm = NaN
    Make/O/D/N=(Nr) w_Lphi_local_nm  = NaN
    Make/O/D/N=(Nr) w_fProbe_local   = NaN
    Make/O/D/N=(Nr) w_xProbe_local_nm = NaN
    Make/O/D/N=(Nr) w_xProbe_global_nm = NaN

    Wave LdiffLocal = w_Ldiff_local_nm
    Wave LphiLocal  = w_Lphi_local_nm
    Wave fProbeLocal = w_fProbe_local
    Wave xProbeLocal = w_xProbe_local_nm
    Wave xProbeGlobal = w_xProbe_global_nm

    // ------------------------------
    // Local ray tracing and local Usadel compression
    // ------------------------------
    Variable ir
    for (ir = 0; ir < Nr; ir += 1)

        SNS_ExtractModesForFolder(maskFolderNoColon, STSx=Xline[ir], STSy=Yline[ir], Bangle_deg=Bangle_deg, doDisplay=doDisplay)

        SNS_Usadel1D_GeometryFromRayFolder(maskFolderColon, STSx=Xline[ir], STSy=Yline[ir], useWeights=useWeights, doPrint=doPrint)

        SNS_Usadel1D_LeverArmFromRayFolder(maskFolderColon, useWeights=useWeights, useMax=useMax, doPrint=doPrint)

        NVAR/Z vLdiff = $(maskFolderColon + "v_Usadel_Ldiff_nm")
        NVAR/Z vLphi  = $(maskFolderColon + "v_Usadel_Lphi_use_nm")
        NVAR/Z vfProbe = $(maskFolderColon + "v_Usadel_fProbe")
        NVAR/Z vxProbe = $(maskFolderColon + "v_Usadel_xProbe_nm")

        if (NVAR_Exists(vLdiff))
            LdiffLocal[ir] = vLdiff
        endif
        if (NVAR_Exists(vLphi))
            LphiLocal[ir] = vLphi
        endif
        if (NVAR_Exists(vfProbe))
            fProbeLocal[ir] = vfProbe
        endif
        if (NVAR_Exists(vxProbe))
            xProbeLocal[ir] = vxProbe
        endif

        SetDataFolder $diagFolder
    endfor

    // ------------------------------
    // Global wire from median over STM positions
    // ------------------------------
    Variable LdiffGlobal_nm = SNS_FiniteMedianFromWave(LdiffLocal)
    Variable LphiGlobal_nm  = SNS_FiniteMedianFromWave(LphiLocal)

    if (numtype(LdiffGlobal_nm) != 0 || LdiffGlobal_nm <= 0)
        SetDataFolder dfrCaller
        Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: invalid global Ldiff."
    endif

    if (numtype(LphiGlobal_nm) != 0 || LphiGlobal_nm <= 0)
        SetDataFolder dfrCaller
        Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: invalid global Lphi."
    endif

    xProbeGlobal = fProbeLocal[p] * LdiffGlobal_nm
    xProbeGlobal = xProbeGlobal[p] < 0 ? 0 : xProbeGlobal[p]
    xProbeGlobal = xProbeGlobal[p] > LdiffGlobal_nm ? LdiffGlobal_nm : xProbeGlobal[p]

    // ------------------------------
    // Effective mfp and diffusion constant
    // ------------------------------
    Wave w_area_perim = SNS_MaskAreaPerim_FromParticles(Nmask)

    Variable area_nm2 = w_area_perim[0]
    Variable perim_nm = w_area_perim[1]

    if (numtype(area_nm2) != 0 || numtype(perim_nm) != 0 || area_nm2 <= 0 || perim_nm <= 0)
        SetDataFolder dfrCaller
        Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: invalid mask area/perimeter."
    endif

    Variable lmfpEdge_nm = params.lmfp_edge_prefactor * area_nm2 / perim_nm
    Variable lmfpEff_nm  = 1 / (1/params.lmfp_bulk_nm + 1/lmfpEdge_nm)

    Variable D_nm2_s = 0.5 * params.vF * lmfpEff_nm * 1e9
    Variable D_m2_s  = D_nm2_s * 1e-18
    Variable D_cm2_s = D_m2_s * 1e4

    Variable/G v_UsadelLine_Ldiff_global_nm = LdiffGlobal_nm
    Variable/G v_UsadelLine_Lphi_global_nm = LphiGlobal_nm
    Variable/G v_UsadelLine_B0_T = B0
    Variable/G v_UsadelLine_Bangle_deg = Bangle_deg
    Variable/G v_UsadelLine_BTK_barrier = BTK_barrier
    Variable/G v_UsadelLine_Ngrid = Ngrid
    Variable/G v_UsadelLine_useWeights = useWeights
    Variable/G v_UsadelLine_useMax = useMax
    Variable/G v_UsadelLine_LineStep_nm = lineStep_nm
    Variable/G v_UsadelLine_LineStep_m = lineStep_m

    if (!ParamIsDefault(NLinePts))
        Variable/G v_UsadelLine_NLinePts_requested = NLinePts
    else
        Variable/G v_UsadelLine_NLinePts_requested = NaN
    endif

    Variable/G v_UsadelLine_area_nm2 = area_nm2
    Variable/G v_UsadelLine_perim_nm = perim_nm
    Variable/G v_UsadelLine_lmfp_bulk_nm = params.lmfp_bulk_nm
    Variable/G v_UsadelLine_lmfp_edge_nm = lmfpEdge_nm
    Variable/G v_UsadelLine_lmfp_eff_nm = lmfpEff_nm
    Variable/G v_UsadelLine_lmfp_edge_prefactor = params.lmfp_edge_prefactor
    Variable/G v_UsadelLine_D_nm2_s = D_nm2_s
    Variable/G v_UsadelLine_D_m2_s = D_m2_s
    Variable/G v_UsadelLine_D_cm2_s = D_cm2_s

    // ------------------------------
    // Allocate output waves in caller folder
    // ------------------------------
    Variable NE = round(params.NE)
    if (mod(NE, 2) == 0)
        NE += 1
    endif

    SetDataFolder dfrCaller

    Make/O/D/N=(NE, Nr) $outDOS
    Wave DOS_Line = $outDOS
    DOS_Line = NaN

    if (doBroadening)
        Make/O/D/N=(NE, Nr) $outDOSConv
        Wave DOS_Line_Conv = $outDOSConv
        DOS_Line_Conv = NaN
    endif

    Duplicate/O $(diagFolder + ":" + nameRaxis), $nameRaxis
    Wave RaxisOut = $nameRaxis

    // ------------------------------
    // Solve Usadel at each mapped probe position
    // ------------------------------
    String tmpFolder = "root:SNS_UsadelLineTmp"
    NewDataFolder/O $tmpFolder

    Variable haveEaxis = 0
    Variable iE

    for (ir = 0; ir < Nr; ir += 1)

        SetDataFolder $tmpFolder

        KillWaves/Z DOS_tmp, DOS_tmp_sym, DOS_tmp_conv, E_tmp
        KillDataFolder/Z :Usadel1D

        SNS_Usadel1D_PhaseBiasNewton_BRange_FromSettings(LdiffGlobal_nm, xProbeGlobal[ir], "DOS_tmp", "E_tmp", D_nm2_s, Ngrid, energySign=energySign, singleB_T=B0, Lphase_nm=LphiGlobal_nm)

        Wave/Z DOS_tmp = DOS_tmp
        Wave/Z E_tmp = E_tmp

        if (!WaveExists(DOS_tmp) || !WaveExists(E_tmp))
            SetDataFolder dfrCaller
            DOS_Line[][ir] = NaN
            if (doBroadening)
                DOS_Line_Conv[][ir] = NaN
            endif
            SetDataFolder $tmpFolder
            continue
        endif

        if (doBroadening)
            SNS_ApplyDOS_Broadening_TplusMod(DOS_tmp, "DOS_tmp_conv")
            Wave/Z DOS_tmp_conv = DOS_tmp_conv
        endif

        if (!haveEaxis)
            SetDataFolder dfrCaller
            Duplicate/O $(tmpFolder + ":E_tmp"), $outEaxis
            Wave EaxisOut = $outEaxis

            SetScale/P x, DimOffset(EaxisOut, 0), DimDelta(EaxisOut, 0), WaveUnits(EaxisOut, 0), DOS_Line
            SetScale/P y, DimOffset(RaxisOut, 0), DimDelta(RaxisOut, 0), WaveUnits(RaxisOut, 0), DOS_Line

            if (doBroadening)
                SetScale/P x, DimOffset(EaxisOut, 0), DimDelta(EaxisOut, 0), WaveUnits(EaxisOut, 0), DOS_Line_Conv
                SetScale/P y, DimOffset(RaxisOut, 0), DimDelta(RaxisOut, 0), WaveUnits(RaxisOut, 0), DOS_Line_Conv
            endif

            haveEaxis = 1
            SetDataFolder $tmpFolder
        endif

        if (DimSize(DOS_tmp, 0) != NE)
            SetDataFolder dfrCaller
            Abort "SNS_Usadel1D_LineDOS_GlobalWireFromMask: DOS_tmp NE mismatch."
        endif

        SetDataFolder dfrCaller
        for (iE = 0; iE < NE; iE += 1)
            DOS_Line[iE][ir] = DOS_tmp[iE][0]
            if (doBroadening && WaveExists(DOS_tmp_conv))
                DOS_Line_Conv[iE][ir] = DOS_tmp_conv[iE][0]
            endif
        endfor

    endfor

    // ------------------------------
    // Metadata
    // ------------------------------
    SetDataFolder dfrCaller

    String meta
    meta  = "SNS_UsadelLineMode=GlobalSingleWire;"
    meta += "SNS_DiagnosticsFolder=" + diagFolder + ":;"
    meta += "SNS_MaskFolder=" + maskFolderColon + ";"
    meta += "SNS_B0_T=" + num2str(B0) + ";"
    meta += "SNS_IsZeroField=" + num2str(isZeroField) + ";"
    meta += "SNS_Bangle_deg=" + num2str(Bangle_deg) + ";"
    meta += "SNS_BTK_barrier=" + num2str(BTK_barrier) + ";"
    meta += "SNS_Ldiff_global_nm=" + num2str(LdiffGlobal_nm) + ";"
    meta += "SNS_Lphi_global_nm=" + num2str(LphiGlobal_nm) + ";"
    meta += "SNS_D_nm2_s=" + num2str(D_nm2_s) + ";"
    meta += "SNS_lmfp_bulk_nm=" + num2str(params.lmfp_bulk_nm) + ";"
    meta += "SNS_lmfp_edge_nm=" + num2str(lmfpEdge_nm) + ";"
    meta += "SNS_lmfp_eff_nm=" + num2str(lmfpEff_nm) + ";"
    meta += "SNS_LineStep_nm=" + num2str(lineStep_nm) + ";"
    meta += "SNS_LineStep_m=" + num2str(lineStep_m) + ";"

    if (!ParamIsDefault(NLinePts))
        meta += "SNS_NLinePts_requested=" + num2str(NLinePts) + ";"
    else
        meta += "SNS_NLinePts_requested=;"
    endif

    meta += "SNS_Raxis=" + GetWavesDataFolder(RaxisOut, 2) + ";"
    meta += "SNS_xProbe_global=" + diagFolder + ":w_xProbe_global_nm;"

    Note/K DOS_Line
    Note DOS_Line, meta

    if (doBroadening)
        Note/K DOS_Line_Conv
        Note DOS_Line_Conv, meta
    endif

    Wave/Z EaxisOutFinal = $outEaxis
    if (WaveExists(EaxisOutFinal))
        Note/K EaxisOutFinal
        Note EaxisOutFinal, meta
    endif

    SetDataFolder dfrCaller
    return Nr
End

