#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma ModuleName=SNS_Legacy_SSpS

//==========================================================================================
// LEGACY: S–S'–S / SSpS solver and DOS utilities
//
// This optional module is not required for the current S–N–S
// mask-derived geometry workflow.
//
// Keep this file only if you still need to reproduce old SSpS results.
// Otherwise it can be removed from the Igor Procedures path.
//==========================================================================================

//==============================================================================
// SSpS_delta_eV
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B : [T]
//   vF : input
//   lambdaL : [m]
//   theta : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SSpS_delta_eV(B, vF, lambdaL, theta)
    Variable B, vF, lambdaL, theta
    return vF * B * lambdaL * sin(theta)
    // return vF * abs(B) * lambdaL * sin(theta)   // use this if you insist on even-in-B physics
End


// Momentum splitting in S' when |δ_theta| > Δp (Δp = DeltaP = gap of S')
// Δk = sgn(delta) * sqrt(delta^2 - Δp^2) / (ħ vF)
//
// Returns DeltaK [1/m] or NaN if gapped (no propagating PF modes).
// S' longitudinal wavevector shift at given energy E (eV)
// Returns real DeltaK if |E - delta_theta| > DeltaP; NaN otherwise (evanescent)
// (You can later extend the evanescent case instead of NaN.)
//==============================================================================
// SSpS_DeltaK_E
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B : [T]
//   vF : input
//   lambdaL : [m]
//   DeltaP : [eV]
//   theta : input
//   E : [eV]
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SSpS_DeltaK_E(B, vF, lambdaL, DeltaP, theta, E)
    Variable B, vF, lambdaL, DeltaP, theta, E

    Variable delta = SSpS_delta_eV(B, vF, lambdaL, theta)
    Variable x = E - delta

    if (abs(x) <= DeltaP)
        return NaN
    endif

    return sign(x) * sqrt(x*x - DeltaP*DeltaP) / (HBAR_eVs * vF)
End



//============================================================
// SS'S quantization root function
//============================================================
//
// Parameter wave layout (MUST match every FindRoots call):
//   pw = {DeltaS, DeltaP, tau, Leff, vF, lambdaL, B, m, sSign, theta}
//
// Interpretation (current model):
//   - Outer S leads are treated as "stiff": Andreev phase uses DeltaS only.
//   - Field dependence enters via propagation phase in S' through DeltaK(B,theta).
//   - sSign = ±1 labels the two propagation directions / phase branches.
//   - m is an integer winding index (allowed to jump during tracking).
//
// Root equation solved:
//   f(E) = 2*acos(E/DeltaS) - ( sSign*(2*DeltaK*Leff) + 2*pi*m ) = 0
//
// Domain:
//   - Solve only for |E| < DeltaS (subgap of outer S leads).
//============================================================
//==============================================================================
// SSpS_BranchFunc
//
// Purpose:
//   Internal solver helper.
//
// Inputs:
//   pw : input
//   E : [eV]
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SSpS_BranchFunc(pw, E)
    Wave pw
    Variable E

    Variable DeltaS  = pw[0]
    Variable DeltaP  = pw[1]
    Variable Leff    = pw[3]
    Variable vF      = pw[4]
    Variable lambdaL = pw[5]
    Variable B       = pw[6]
    Variable m       = pw[7]
    Variable sSign   = pw[8]
    Variable theta   = pw[9]

    if (abs(E) >= DeltaS)
        return 1e6
    endif

	Variable DeltaK = SSpS_DeltaK_E(B, vF, lambdaL, DeltaP, theta, E)
	if (numtype(DeltaK) != 0)
    	return 1e6
	endif


    // Andreev phase from the *lead* S only (stiff-S approximation)
    Variable phiA = 2*AcosSafe(E/DeltaS)

    // Propagation phase in S' (PF)
    Variable phiP = 2*DeltaK*Leff

    return phiA - (sSign*phiP + 2*pi*m)
End


//============================================================
// Robust single-branch solver using FindRoots + bracket
// (compatible with E-dependent DeltaK inside SSpS_BranchFunc)
//============================================================
//
// Purpose:
//   - Solve one (m,sSign) branch at a given B.
//   - Does NOT use PF-only DeltaK(B) checks (invalid once DeltaK depends on E).
//   - Uses a conservative bracketing strategy and lets FindRoots do the work.
//   - Still guards against edge-hit artifacts at |E| ~ DeltaS.
//
// Returns E (eV) or NaN
//============================================================
//==============================================================================
// SolveBranchAtB_SSpS
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B : [T]
//   L : input
//   W : input
//   DeltaS : [eV]
//   DeltaP : [eV]
//   vF : input
//   lambdaL : [m]
//   m : input
//   sSign : input
//   tolLocal : input
//   maxI : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SolveBranchAtB_SSpS(B, L, W, DeltaS, DeltaP, vF, lambdaL, m, sSign, tolLocal, maxI)
    Variable B, L, W, DeltaS, DeltaP, vF, lambdaL, m, sSign, tolLocal, maxI

    Variable theta = SNS_theta(L, W)
    Variable Leff  = L/cos(theta)

    // We cannot compute a meaningful phiP from DeltaK(B) anymore because DeltaK depends on E.
    // Use a neutral initial guess inside the gap and bracket around it.
    Variable edgePad = 1e-9
    Variable Eguess  = 0

    // If you expect solutions mainly above DeltaP, bias the guess slightly outward:
    // (helps FindRoots for the |E| > DeltaP regime without forcing anything)
    if (DeltaP > 0 && DeltaP < 0.9*DeltaS)
        Eguess = 0.5*(DeltaP + DeltaS)
        // keep sign convention neutral; FindRoots will converge to the root if bracket contains it
        Eguess = sign(sSign) * Eguess
    endif

    // Start with a modest bracket; widen aggressively if needed.
    Variable dEBracket = 2e-3*DeltaS
    Variable lo = max(-DeltaS + edgePad, Eguess - dEBracket)
    Variable hi = min( DeltaS - edgePad, Eguess + dEBracket)

    Variable tries
    for (tries = 0; tries < 14; tries += 1)

        Make/FREE/D/N=10 pw
        pw[0] = DeltaS
        pw[1] = DeltaP
        pw[2] = 0
        pw[3] = Leff
        pw[4] = vF
        pw[5] = lambdaL
        pw[6] = B
        pw[7] = m
        pw[8] = sSign
        pw[9] = theta

        FindRoots /Q /B=1 /T=(tolLocal) /I=(maxI) /L=(lo) /H=(hi) SSpS_BranchFunc, pw

        if (V_flag == 0 && numtype(V_Root) == 0)
            if (abs(V_Root) < DeltaS - 2e-9)     // edge-hit guard
                return V_Root
            endif
        endif

        // widen bracket and retry
        dEBracket *= 2.5
        lo = max(-DeltaS + edgePad, Eguess - dEBracket)
        hi = min( DeltaS - edgePad, Eguess + dEBracket)

        // if we essentially cover the full gap and still fail, stop
        if (lo <= -DeltaS + 5e-9 && hi >= DeltaS - 5e-9)
            break
        endif
    endfor

    return NaN
End


//============================================================
// Tracked single-branch solver using FindRoots (continuity)
// (compatible with E-dependent DeltaK inside SSpS_BranchFunc)
//============================================================
//
// Purpose:
//   - Solve one branch at B while enforcing continuity from Eprev.
//   - Does NOT use PF-only DeltaK(B) checks (invalid once DeltaK depends on E).
//   - Brackets around Eprev and expands until convergence or full-gap.
//
// Returns E (eV) or NaN
//============================================================
//==============================================================================
// SolveBranchAtB_SSpS_Tracked
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B : [T]
//   L : input
//   W : input
//   DeltaS : [eV]
//   DeltaP : [eV]
//   vF : input
//   lambdaL : [m]
//   mBase : input
//   sSign : input
//   Eprev : input
//   tolLocal : input
//   maxI : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SolveBranchAtB_SSpS_Tracked(B, L, W, DeltaS, DeltaP, vF, lambdaL, mBase, sSign, Eprev, tolLocal, maxI)
    Variable B, L, W, DeltaS, DeltaP, vF, lambdaL, mBase, sSign, Eprev, tolLocal, maxI

    Variable theta = SNS_theta(L, W)
    Variable Leff  = L/cos(theta)

    Variable edgePad = 1e-9

    // If previous is NaN, fall back to center-ish guess
    if (numtype(Eprev) != 0)
        Eprev = 0
    endif

    // Bracket around previous E (continuity)
    Variable dEBracket = 5e-5*DeltaS
    Variable lo = max(-DeltaS + edgePad, Eprev - dEBracket)
    Variable hi = min( DeltaS - edgePad, Eprev + dEBracket)

    if (hi <= lo)
        lo = -DeltaS + edgePad
        hi =  DeltaS - edgePad
    endif

    Variable tries
    for (tries = 0; tries < 14; tries += 1)

        Make/FREE/D/N=10 pw
        pw[0] = DeltaS
        pw[1] = DeltaP
        pw[2] = 0
        pw[3] = Leff
        pw[4] = vF
        pw[5] = lambdaL
        pw[6] = B
        pw[7] = mBase
        pw[8] = sSign
        pw[9] = theta

        FindRoots /Q /B=1 /T=(tolLocal) /I=(maxI) /L=(lo) /H=(hi) SSpS_BranchFunc, pw

        if (V_flag == 0 && numtype(V_Root) == 0)
            if (abs(V_Root) < DeltaS - 2e-9)     // edge-hit guard
                return V_Root
            endif
        endif

        // widen bracket around Eprev
        dEBracket *= 2
        lo = max(-DeltaS + edgePad, Eprev - dEBracket)
        hi = min( DeltaS - edgePad, Eprev + dEBracket)

        if (lo <= -DeltaS + 5e-9 && hi >= DeltaS - 5e-9)
            break
        endif
    endfor

    return NaN
End


//============================================================
// SS'S ROOT-FINDING AND BRANCH TRACKING (FindRoots-based)
// Patched for E-dependent propagation in S' (DeltaK depends on E inside SSpS_BranchFunc)
//============================================================


//------------------------------------------------------------
// Low-level: try FindRoots in [Elo,Ehi] for given (m,s) at B (and given theta)
// NOTE: Do NOT PF-kill based on DeltaK(B) here; validity is decided inside SSpS_BranchFunc.
//------------------------------------------------------------
//==============================================================================
// SSpS_TryRootBracket
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B : [T]
//   Leff : input
//   DeltaS : [eV]
//   DeltaP : [eV]
//   vF : input
//   lambdaL : [m]
//   theta : input
//   mTry : input
//   sTry : input
//   Elo : input
//   Ehi : input
//   tolLocal : input
//   maxI : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SSpS_TryRootBracket(B, Leff, DeltaS, DeltaP, vF, lambdaL, theta, mTry, sTry, Elo, Ehi, tolLocal, maxI)
    Variable B, Leff, DeltaS, DeltaP, vF, lambdaL, theta, mTry, sTry, Elo, Ehi, tolLocal, maxI

    Variable edgePad = 1e-9
    if (Elo < -DeltaS + edgePad)
        Elo = -DeltaS + edgePad
    endif
    if (Ehi >  DeltaS - edgePad)
        Ehi =  DeltaS - edgePad
    endif
    if (Ehi <= Elo)
        return NaN
    endif

    Make/FREE/D/N=10 pw
    pw[0] = DeltaS
    pw[1] = DeltaP
    pw[2] = 0
    pw[3] = Leff
    pw[4] = vF
    pw[5] = lambdaL
    pw[6] = B
    pw[7] = mTry
    pw[8] = sTry
    pw[9] = theta

    FindRoots /Q /B=1 /T=(tolLocal) /I=(maxI) /L=(Elo) /H=(Ehi) SSpS_BranchFunc, pw

    if (V_flag == 0 && numtype(V_Root) == 0)
        if (abs(V_Root) < DeltaS - 2e-9)
            // Soft residual check
            Variable fAtRoot = SSpS_BranchFunc(pw, V_Root)
            if (numtype(fAtRoot) == 0 && abs(fAtRoot) < 1e-4)
                return V_Root
            endif
        endif
    endif

    return NaN
End


//------------------------------------------------------------
// Continuation step for one branch at fixed theta
// Given previous (Eprev, mPrev, sPrev), try to find a root at this B near Eprev.
// Allows m jumps and (optionally) sSign flips.
// NOTE: No PF-kill based on DeltaK(B); SSpS_BranchFunc handles validity.
//------------------------------------------------------------
//==============================================================================
// SSpS_ContinueOne
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B : [T]
//   Leff : input
//   DeltaS : [eV]
//   DeltaP : [eV]
//   vF : input
//   lambdaL : [m]
//   theta : input
//   Eprev : input
//   mPrev : input
//   sPrev : input
//   mJump : input
//   allowSFlip : input
//   tolLocal : input
//   maxI : input
//   outMS : input
//
// Outputs:
//   return : numeric
//   outMS : output wave(s)
//
//==============================================================================
Function SSpS_ContinueOne(B, Leff, DeltaS, DeltaP, vF, lambdaL, theta, Eprev, mPrev, sPrev, mJump, allowSFlip, tolLocal, maxI, outMS)
    Variable B, Leff, DeltaS, DeltaP, vF, lambdaL, theta, Eprev, mPrev, sPrev, mJump, allowSFlip, tolLocal, maxI
    Wave outMS

    outMS[0] = NaN
    outMS[1] = NaN

    if (numtype(Eprev) != 0)
        return NaN
    endif

    Variable dE = 5e-6*DeltaS
    Variable tries

    Variable bestE = NaN
    Variable bestErr = 1e99
    Variable bestM = NaN
    Variable bestS = NaN

    Variable mTry, sTry, Etest, err, Elo, Ehi
    Variable sLoopStart = sPrev
    Variable sLoopStop  = sPrev
    Variable sLoopStep  = 2

    if (allowSFlip != 0)
        sLoopStart = -1
        sLoopStop  = 1
        sLoopStep  = 2
    endif

    for (tries = 0; tries < 10; tries += 1)

        Elo = Eprev - dE
        Ehi = Eprev + dE

        for (sTry = sLoopStart; sTry <= sLoopStop; sTry += sLoopStep)
            for (mTry = mPrev - mJump; mTry <= mPrev + mJump; mTry += 1)

                Etest = SSpS_TryRootBracket(B, Leff, DeltaS, DeltaP, vF, lambdaL, theta, mTry, sTry, Elo, Ehi, tolLocal, maxI)
                if (numtype(Etest) == 0)
                    err = abs(Etest - Eprev)
                    if (err < bestErr)
                        bestErr = err
                        bestE = Etest
                        bestM = mTry
                        bestS = sTry
                    endif
                endif

            endfor
        endfor

        if (numtype(bestE) == 0)
            outMS[0] = bestM
            outMS[1] = bestS
            return bestE
        endif

        dE *= 3
    endfor

    return NaN
End


//------------------------------------------------------------
// Discovery at B for fixed theta:
// With E-dependent propagation, phiP(B) is not single-valued (depends on E),
// so "mCenter from phiP" is invalid.
// Instead: scan a fixed m-range for both sSign values and let FindRoots decide.
//------------------------------------------------------------
//==============================================================================
// SSpS_DiscoverAtB
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B : [T]
//   Leff : input
//   DeltaS : [eV]
//   DeltaP : [eV]
//   vF : input
//   lambdaL : [m]
//   theta : input
//   mWindow : input
//   tolLocal : input
//   maxI : input
//   Elist : input
//   Mlist : input
//   Slist : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SSpS_DiscoverAtB(B, Leff, DeltaS, DeltaP, vF, lambdaL, theta, mWindow, tolLocal, maxI, Elist, Mlist, Slist)
    Variable B, Leff, DeltaS, DeltaP, vF, lambdaL, theta, mWindow, tolLocal, maxI
    Wave Elist, Mlist, Slist

    Redimension/N=0 Elist, Mlist, Slist

    Variable sTry, mTry, Etest

    for (sTry = -1; sTry <= 1; sTry += 2)

        // Scan m in a fixed window around 0. This is robust when m bookkeeping is not important.
        for (mTry = -mWindow; mTry <= mWindow; mTry += 1)

            Etest = SSpS_TryRootBracket(B, Leff, DeltaS, DeltaP, vF, lambdaL, theta, mTry, sTry, -DeltaS, DeltaS, tolLocal, maxI)

            if (numtype(Etest) == 0)
                Redimension/N=(numpnts(Elist)+1) Elist, Mlist, Slist
                Elist[numpnts(Elist)-1] = Etest
                Mlist[numpnts(Mlist)-1] = mTry
                Slist[numpnts(Slist)-1] = sTry
            endif

        endfor
    endfor

    return numpnts(Elist)
End


//------------------------------------------------------------
// Utility: check if E is already in E2D row within eps
//------------------------------------------------------------
//==============================================================================
// SSpS_IsEnergyAlreadyPresent
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   E2D : input
//   row : input
//   Etest : input
//   epsE : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SSpS_IsEnergyAlreadyPresent(E2D, row, Etest, epsE)
    Wave E2D
    Variable row, Etest, epsE

    Variable nCol = DimSize(E2D, 1)
    Variable j
    for (j = 0; j < nCol; j += 1)
        if (numtype(E2D[row][j]) == 0)
            if (abs(E2D[row][j] - Etest) < epsE)
                return 1
            endif
        endif
    endfor
    return 0
End


//============================================================
// Solve_AllBranches_SSpS (continuation + discovery)
// Patched seeding: choose seed B where discovery finds roots (not where DeltaK(B) is real).
//============================================================
//==============================================================================
// Solve_AllBranches_SSpS
//
// Purpose:
//   Solve subgap bound-state energies and/or track branches.
//
// Inputs:
//   B_T : [T]
//   L : input
//   W : input
//   DeltaS : [eV]
//   DeltaP : [eV]
//   vF : input
//   lambdaL : [m]
//   nameE2D : input
//   nameM : input
//   nameS : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Solve_AllBranches_SSpS(B_T, L, W, DeltaS, DeltaP, vF, lambdaL, nameE2D, nameM, nameS)
    Wave B_T
    Variable L, W, DeltaS, DeltaP, vF, lambdaL
    String nameE2D, nameM, nameS

    Variable nB = numpnts(B_T)
    if (nB <= 0)
        Abort "Solve_AllBranches_SSpS: B_T has no points."
    endif

    Variable theta = SNS_theta(L, W)
    Variable Leff  = L/cos(theta)

    Variable tolLocal = 1e-10
    Variable maxI     = 400

    Variable mJump      = 3
    Variable allowSFlip = 1
    Variable epsMergeE  = 2e-7

    Variable mWindowSeed = 14
    Variable mWindowStep = 8

    // ------------------------------------------------------------
    // Seed at smallest |B| where discovery finds at least one root
    // ------------------------------------------------------------
    Variable i, iSeed = -1
    Variable bestAbsB = 1e99

    for (i = 0; i < nB; i += 1)

        Make/FREE/D/N=0 Etmp, Mtmp, Stmp
        Variable nTmp = SSpS_DiscoverAtB(B_T[i], Leff, DeltaS, DeltaP, vF, lambdaL, theta, mWindowSeed, tolLocal, maxI, Etmp, Mtmp, Stmp)

        if (nTmp > 0)
            if (abs(B_T[i]) < bestAbsB)
                bestAbsB = abs(B_T[i])
                iSeed = i
            endif
        endif

    endfor

    if (iSeed < 0)
        Print "Solve_AllBranches_SSpS: no roots found anywhere in B_T."
        Make/O/D/N=(nB, 0) $nameE2D
        Make/O/D/N=0 $nameM, $nameS
        return 0
    endif

    Variable Bseed = B_T[iSeed]

    // Discover initial roots at seed
    Make/FREE/D/N=0 Eseed, Mseed, Sseed
    Variable nSeed = SSpS_DiscoverAtB(Bseed, Leff, DeltaS, DeltaP, vF, lambdaL, theta, mWindowSeed, tolLocal, maxI, Eseed, Mseed, Sseed)

    if (nSeed <= 0)
        Print "Solve_AllBranches_SSpS: no roots found at seed B=", Bseed, " T (iSeed=", iSeed, ")."
        Make/O/D/N=(nB, 0) $nameE2D
        Make/O/D/N=0 $nameM, $nameS
        return 0
    endif

    // Create outputs
    Make/O/D/N=(nB, nSeed) $nameE2D
    Wave E2D = $nameE2D
    E2D = NaN

    Make/O/D/N=(nSeed) $nameM, $nameS
    Wave mBranch = $nameM
    Wave sBranch = $nameS

    Variable j
    for (j = 0; j < nSeed; j += 1)
        E2D[iSeed][j] = Eseed[j]
        mBranch[j] = Mseed[j]
        sBranch[j] = Sseed[j]
    endfor

    // Save seed labels for backward sweep
    Make/FREE/D/N=(nSeed) mBranchSeed, sBranchSeed
    mBranchSeed = mBranch
    sBranchSeed = sBranch

    Make/FREE/D/N=2 outMS

    // -------- Forward sweep --------
    Variable iStep, Eprev, Enew
    for (iStep = iSeed + 1; iStep < nB; iStep += 1)

        for (j = 0; j < DimSize(E2D, 1); j += 1)
            Eprev = E2D[iStep - 1][j]
            if (numtype(Eprev) != 0)
                continue
            endif

            Enew = SSpS_ContinueOne(B_T[iStep], Leff, DeltaS, DeltaP, vF, lambdaL, theta, Eprev, mBranch[j], sBranch[j], mJump, allowSFlip, tolLocal, maxI, outMS)

            E2D[iStep][j] = Enew
            if (numtype(Enew) == 0)
                mBranch[j] = outMS[0]
                sBranch[j] = outMS[1]
            endif
        endfor

        Make/FREE/D/N=0 Edisc, Mdisc, Sdisc
        Variable nDisc = SSpS_DiscoverAtB(B_T[iStep], Leff, DeltaS, DeltaP, vF, lambdaL, theta, mWindowStep, tolLocal, maxI, Edisc, Mdisc, Sdisc)

        for (j = 0; j < nDisc; j += 1)
            if (SSpS_IsEnergyAlreadyPresent(E2D, iStep, Edisc[j], epsMergeE) == 0)
                Variable nColNew = DimSize(E2D, 1) + 1
                Redimension/N=(nB, nColNew) E2D
                E2D[][nColNew-1] = NaN
                E2D[iStep][nColNew-1] = Edisc[j]

                Redimension/N=(nColNew) mBranch, sBranch
                mBranch[nColNew-1] = Mdisc[j]
                sBranch[nColNew-1] = Sdisc[j]
            endif
        endfor

    endfor

    // Restore seed labels before backward sweep (only original seed branches)
    for (j = 0; j < nSeed; j += 1)
        mBranch[j] = mBranchSeed[j]
        sBranch[j] = sBranchSeed[j]
    endfor

    // -------- Backward sweep --------
    for (iStep = iSeed - 1; iStep >= 0; iStep -= 1)

        for (j = 0; j < DimSize(E2D, 1); j += 1)
            Eprev = E2D[iStep + 1][j]
            if (numtype(Eprev) != 0)
                continue
            endif

            Enew = SSpS_ContinueOne(B_T[iStep], Leff, DeltaS, DeltaP, vF, lambdaL, theta, Eprev, mBranch[j], sBranch[j], mJump, allowSFlip, tolLocal, maxI, outMS)

            E2D[iStep][j] = Enew
            if (numtype(Enew) == 0)
                mBranch[j] = outMS[0]
                sBranch[j] = outMS[1]
            endif
        endfor

        Make/FREE/D/N=0 Edisc2, Mdisc2, Sdisc2
        Variable nDisc2 = SSpS_DiscoverAtB(B_T[iStep], Leff, DeltaS, DeltaP, vF, lambdaL, theta, mWindowStep, tolLocal, maxI, Edisc2, Mdisc2, Sdisc2)

        for (j = 0; j < nDisc2; j += 1)
            if (SSpS_IsEnergyAlreadyPresent(E2D, iStep, Edisc2[j], epsMergeE) == 0)
                Variable nColNew2 = DimSize(E2D, 1) + 1
                Redimension/N=(nB, nColNew2) E2D
                E2D[][nColNew2-1] = NaN
                E2D[iStep][nColNew2-1] = Edisc2[j]

                Redimension/N=(nColNew2) mBranch, sBranch
                mBranch[nColNew2-1] = Mdisc2[j]
                sBranch[nColNew2-1] = Sdisc2[j]
            endif
        endfor

    endfor

    return DimSize(E2D, 1)
End


//============================================================
// Compute_DOS_SSpS_Map_AllBranches
//============================================================
//
// Builds a DOS(E,B) map by:
//   1) Sampling "channels" via a list of trajectory angles thetaList[j]
//      derived from a width coordinate wList[j] and SNS_theta(L, w).
//   2) For each channel, solving SS'S ABS branches vs B using Solve_AllBranches_SSpS.
//   3) Accumulating a Lorentzian DOS contribution at ±E0 for particle-hole symmetry.
//
// IMPORTANT (geometry / bookkeeping):
//   - Solve_AllBranches_SSpS internally computes theta = SNS_theta(L, W_arg).
//   - In this DOS builder we want theta to be the channel angle thetaList[j].
//   - We therefore pass W_arg = W_theta = L*tan(thetaList[j]) so that
//       SNS_theta(L, W_arg) returns thetaList[j] again.
//   - W_theta here is NOT the physical junction width; it's a proxy that
//     reproduces the desired channel angle inside the solver.
//============================================================
//==============================================================================
// Compute_DOS_SSpS_Map_AllBranches
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
Function Compute_DOS_SSpS_Map_AllBranches(B_T, L, Wmax, Delta, vF, lambdaL, T, Broadening, lambdaF, NE, nameDOS, nameEaxis)
    Wave B_T
    Variable L, Wmax, Delta, vF, lambdaL, T
    Variable Broadening      // Lorentzian broadening [eV]
    Variable lambdaF         // Fermi wavelength [m]
    Variable NE              // number of energy points
    String nameDOS, nameEaxis

    Variable nB = numpnts(B_T)
    if (nB <= 0)
        Abort "Compute_DOS_SSpS_Map_AllBranches: B_T has no points."
    endif

    // In this SS'S model we use DeltaS = Delta and DeltaP = 0.5*Delta by default
    Variable DeltaSP = 0.5*Delta

    // ---------------------------
    // 1) Build energy axis
    // ---------------------------
    Make/O/D/N=(NE) $nameEaxis
    Wave E_axis = $nameEaxis
    SetScale/I x, -Delta, Delta, E_axis
    E_axis = x

    // DOS(E,B) 2D map: [iE][iB]
    Make/O/D/N=(NE, nB) $nameDOS
    Wave DOS_EB = $nameDOS
    DOS_EB = 0

    // Axis scales
    SetScale/P x, DimOffset(E_axis, 0), DimDelta(E_axis, 0), "eV", DOS_EB
    SetScale/P y, DimOffset(B_T,    0), DimDelta(B_T,    0), "T",  DOS_EB

    // ---------------------------
    // 2) Build channel sampling from Fermi wavelength
    // ---------------------------
    // Nch ≈ 2 * Wmax / lambdaF (simple transverse-mode count proxy)
    Variable Ntheta = trunc(2*Wmax / lambdaF)
    if (Ntheta < 1)
        Ntheta = 1
    endif

    Variable Nch = Ntheta
    Variable dw  = Wmax / Nch

    Make/O/D/N=(Nch) thetaList, wList
    Variable j
    for (j = 0; j < Nch; j += 1)
        wList[j]     = (j + 0.5) * dw
        thetaList[j] = SNS_theta(L, wList[j])     // typically atan(w/L)
    endfor

    // ---------------------------
    // 3) Loop over channels, solve branches, accumulate DOS
    // ---------------------------
    Variable iB, k, iE
    Variable theta, W_theta, nBr, E0, dE
    Variable weight = 1.0 / Nch

    String nameE2D, nameM, nameS

    for (j = 0; j < Nch; j += 1)

        theta   = thetaList[j]

        // Proxy width argument that reproduces this theta inside Solve_AllBranches_SSpS
        // because Solve_AllBranches_SSpS uses theta = SNS_theta(L, W_arg).
        W_theta = L * tan(theta)

        sprintf nameE2D, "E_allBranches_th%03d", j
        sprintf nameM,   "m_allBranches_th%03d", j
        sprintf nameS,   "s_allBranches_th%03d", j

        nBr = Solve_AllBranches_SSpS(B_T, L, W_theta, Delta, DeltaSP, vF, lambdaL, nameE2D, nameM, nameS)
        if (nBr <= 0)
            Print "Compute_DOS_SSpS_Map_AllBranches: no branches at channel index ", j
            continue
        endif

        Wave E_all = $nameE2D

        for (iB = 0; iB < nB; iB += 1)
            for (k = 0; k < nBr; k += 1)

                E0 = E_all[iB][k]
                if (numtype(E0) != 0)
                    continue
                endif

                // Add Lorentzian peaks at +E0 and -E0 (enforce particle-hole symmetry in DOS)
                for (iE = 0; iE < NE; iE += 1)

                    dE = E_axis[iE] - E0
                    DOS_EB[iE][iB] += weight * (Broadening/pi) / (dE*dE + Broadening*Broadening)

                    dE = E_axis[iE] + E0
                    DOS_EB[iE][iB] += weight * (Broadening/pi) / (dE*dE + Broadening*Broadening)

                endfor

            endfor
        endfor

    endfor

    return Nch
End



//==============================================================================
// SNS_theta
//
// Purpose:
//   Return the channel opening angle θ from chord length and width.
//
// Inputs:
//   L : chord length [m]
//   W : transverse width [m]
//
// Outputs:
//   return : θ = atan2(W, L) [rad]
//
// Notes:
//   Used only in legacy SSpS geometry calculations.
//
//==============================================================================
Static Function SNS_theta(L, W)
    Variable L, W
    return atan2(W, L)
End

