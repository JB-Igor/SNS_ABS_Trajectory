#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma ModuleName=SNS_Solver
#include "SNS_Core"

//==============================================================================
// SNS_Solver.ipf
//
// Model:
//   Subgap ABS quantization in a ballistic SNS channel using the
//   de Gennes–Saint-James (dGSJ) phase-accumulation condition, generalized
//   to finite normal-state transparency T.
//
// Reference:
//   P.G. de Gennes and D. Saint-James, Phys. Lett. 4, 151 (1963).
//
// Quantization (T = 1, Bohr–Sommerfeld form):
//   2*acos(E/Δ) - 2*E*τ  ∓  β  =  2π*n
//
//   Variables : input
//     E   : quasiparticle energy [eV]
//     Δ   : superconducting gap magnitude [eV]
//     τ   : flight-time parameter [1/eV] (e.g. L/(ħ v_F) in a ballistic model)
//     β   : total magnetic phase [rad]
//           (includes orbital/superflow contribution and optional vortex phase)
//     n   : integer branch index (…,-1,0,1,…)
//
// Finite transparency (0 ≤ T ≤ 1):
//   The condition is implemented in cosine form:
//
//     cos( Φ(E) ) = (1 - T) + T*cos(β)
//
//     with Φ(E) = 2*acos(E/Δ) - 2*E*τ.
//
//   Equivalently (phase-accumulation form):
//
//     Φ(E) = 2π*n  ±  δ(T,β)
//
//     δ(T,β) = acos( (1 - T) + T*cos(β) ).
//
//   Here T is the channel transmission probability (dimensionless).
//
// Conventions used in code:
//   - acos() arguments are clamped to [-1,1] via AcosSafe().
//   - β_total may be composed as β(B) + β_extra, where β_extra can encode
//     vortex-induced phase winding or other phase offsets.
//
//==============================================================================


//==============================================================================
// SNS_tau_eVs
//
// Purpose:
//   Convert a ballistic flight time L/vF into τ in units of 1/eV.
//
// Inputs:
//   L  : chord length [m]
//   vF : Fermi velocity [m/s]
//
// Outputs:
//   return : τ = (L/vF)/ħ  [1/eV], with ħ in eV*s
//
//==============================================================================
Function SNS_tau_eVs(L, vF)
    Variable L, vF
    return (L / vF) / HBAR_eVs // + ((200e-9 / vF) / HBAR_eVs) testing constant offset penetration into "S"
End


// Orbital phase β(B) with the *correct* prefactor:
// β = - (2 e B λ_L W) / ħ  (your v4 geometry, not the old 4e version)

//==============================================================================
// SNS_beta
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B : [T]
//   L : input
//   W : input
//   lambdaL : [m]
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_beta(B, L, W, lambdaL)
    Variable B, L, W, lambdaL

    // B [T], λ_L [m], W [m] → β dimensionless
    return - (2.0 * q_e * B * lambdaL * W) / HBAR_SI
End

//============================================================
// 3. ABS QUANTIZATION EQUATION (ONE EQUATION FOR ALL T)
//============================================================

// Finite-transparency dGSJ condition (T ∈ (0,1], T=1 → transparent limit)
// F(E) = cos(2 arccos(E/Δ) - 2Eτ) - [(1-T) + T cosβ] = 0

//==============================================================================
// PhaseEq_dGSJ
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
// Notes:
//   Quantization is evaluated via cos(Φ(E)) = (1-T) + T*cos(β), with Φ(E)=2*acos(E/Δ)-2*E*τ.
//   All variables and conventions are defined in the module header.
//
//==============================================================================
Function PhaseEq_dGSJ(pw, E)
    WAVE pw
    Variable E

    Variable Delta = pw[0]
    Variable tau   = pw[1]
    Variable beta  = pw[2]
    Variable T     = pw[3]

    // Guard against nonsense
    if (numtype(Delta + tau + beta + T) != 0)
        return NaN
    endif

    // Subgap guard (optional; solver should enforce anyway)
    if (abs(E) >= Delta)
        return NaN
    endif

    Variable arg = 2.0*AcosSafe((E)/Delta) - 2.0*E*tau
    
    
    Variable lhs = cos(arg)
    Variable rhs = (1.0 - T) + T*cos(beta)

    return lhs - rhs
End

// Helper to fill the parameter wave pw = [Δ, τ, β, T]

//==============================================================================
// FillParams_dGSJ
//
// Purpose:
//   Internal solver helper.
//
// Inputs:
//   pw : input
//   Delta : [eV]
//   tau : input
//   beta : [rad]
//   T : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function FillParams_dGSJ(pw, Delta, tau, beta, T)
    WAVE pw
    Variable Delta, tau, beta, T

    if (numpnts(pw) < 4)
        Redimension/N=4 pw
    endif

    pw[0] = Delta
    pw[1] = tau
    pw[2] = beta
    pw[3] = T
End

//============================================================
// 4. ROOT-FINDING AND BRANCH TRACKING  (FindRoots-based)
//============================================================

// dGSJ branch function for FindRoots
// We solve: Φ(E) = sSign*θ0 + 2π m
// with Φ(E)   = 2*acos(E/Δ) - 2*E*τ
//     C      = (1-T) + T*cos(β)
//     θ0     = acos(C)
// so f(E) = Φ(E) - (sSign*θ0 + 2*pi*m) = 0
//
// pw[0] = Delta
// pw[1] = tau
// pw[2] = beta
// pw[3] = T
// pw[4] = m          (integer index, stored as double)
// pw[5] = sSign      (+1 or -1)

//==============================================================================
// dGSJ_BranchFunc
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
Function dGSJ_BranchFunc(pw, E)
    Wave pw
    Variable E

    Variable Delta = pw[0]
    Variable tau   = pw[1]
    Variable beta  = pw[2]
    Variable T     = pw[3]
    Variable m     = pw[4]
    Variable sSign = pw[5]

    // Basic sanity
    if (numtype(Delta + tau + beta + T + m + sSign) != 0)
        return 1e6
    endif

    // Enforce subgap domain
    if (abs(E) >= Delta)
        return 1e6
    endif

    // RHS of original cos equation
    Variable C = (1 - T) + T*cos(beta)
    if (abs(C) > 1)
        // no subgap solution for these parameters
        return 1e6
    endif

    Variable theta0 = AcosSafe(C)                       // θ0 ∈ [0,π]
    Variable phi    = 2*AcosSafe(E/Delta) - 2*E*tau     // Φ(E)
    Variable target = sSign*theta0 + 2*pi*m

    // f(E) = 0 ⇒ solution
    return phi - target
End

// Helper: solve *one* branch (m, sSign) at a given B.
// Returns the root E (or NaN if no valid subgap solution).

//==============================================================================
// SolveBranchAtB_dGSJ
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B : [T]
//   L : input
//   W : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   T : input
//   m : input
//   sSign : input
//   tol : input
//   maxIters : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SolveBranchAtB_dGSJ(B, L, W, Delta, vF, lambdaL, T, m, sSign, tol, maxIters)
    Variable B, L, W, Delta, vF, lambdaL, T
    Variable m, sSign
    Variable tol, maxIters

    // Geometry / physics
    Variable tau  = SNS_tau_eVs(L, vF)
    Variable beta = SNS_beta(B, L, W, lambdaL)

    // Parameter wave for dGSJ_BranchFunc
    Make/FREE/D/N=6 pw
    pw[0] = Delta
    pw[1] = tau
    pw[2] = beta
    pw[3] = T
    pw[4] = m
    pw[5] = sSign

    // Subgap bracket
    Variable eps = 1e-9
    Variable lo  = -Delta + eps
    Variable hi  =  Delta - eps

    // Use Igor's root finder
    // /B=1: let Igor try to bracket if needed
    // /L, /H: initial bracket in the full subgap window
    FindRoots /Q /B=1 /T=(tol) /I=(maxIters) /L=(lo) /H=(hi) dGSJ_BranchFunc, pw

    if (V_flag != 0)
        // root search failed
        return NaN
    endif

    // V_Root is the solution candidate; check it's subgap
    if (abs(V_Root) >= Delta)
        return NaN
    endif

    // Optionally check residual; if it's huge, reject
    if (numtype(V_YatRoot) != 0 || abs(V_YatRoot) > 1e-4)
        return NaN
    endif

    return V_Root
End

//============================================================
// Low-level: solve one branch (m,sSign) at field B
// with an extra phase betaExtra added to SNS_beta.
//
// Returns E(B,m,s) in [-Delta,Delta]. NaN if no root.
//============================================================

//==============================================================================
// SolveBranchAtB_dGSJ_betaExtra
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   B : [T]
//   L : input
//   W : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   T : input
//   betaExtra : input
//   m : input
//   sSign : input
//   tolLocal : input
//   maxI : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SolveBranchAtB_dGSJ_betaExtra(B, L, W, Delta, vF, lambdaL, T, betaExtra, m, sSign, tolLocal, maxI)
    Variable B, L, W, Delta, vF, lambdaL, T
    Variable betaExtra, m, sSign, tolLocal, maxI

    Variable tau  = SNS_tau_eVs(L, vF)
    Variable beta = SNS_beta(B, L, W, lambdaL) + betaExtra

    Make/FREE/D/N=6 pw
    pw[0] = Delta
    pw[1] = tau
    pw[2] = beta
    pw[3] = T
    pw[4] = m
    pw[5] = sSign

    Variable eps = 1e-9
    Variable lo  = -Delta + eps
    Variable hi  =  Delta - eps

    FindRoots /Q /B=1 /T=(tolLocal) /I=(maxI) /L=(lo) /H=(hi) dGSJ_BranchFunc, pw

    if (V_flag == 0 && numtype(V_Root) == 0 && abs(V_Root) < Delta && abs(V_YatRoot) < 1e-4)
        return V_Root
    endif

    // failed: return NaN so caller can handle it
    return NaN
End


//============================================================
// Solve all subgap branches and return an E(B, branch) matrix,
// including an extra phase offset betaExtra for this channel.
//
// NEW compared to your original:
//   - extra argument betaExtra
//   - uses SolveBranchAtB_dGSJ_betaExtra instead of SolveBranchAtB_dGSJ
//
// Everything else is identical to your running code.
//============================================================

//==============================================================================
// Solve_AllBranches_SNS_dGSJ_betaExtra
//
// Purpose:
//   Solve subgap bound-state energies and/or track branches.
//
// Inputs:
//   B_T : [T]
//   L : input
//   W : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   T : input
//   betaExtra : input
//   nameE2D : input
//   nameM : input
//   nameS : input
//
// Outputs:
//   return : numeric
//
// Notes:
//   Core SNS ABS solver. Uses β_total(B)=β_orbital(B)+betaExtra.
//
//==============================================================================
Function Solve_AllBranches_SNS_dGSJ_betaExtra(B_T, L, W, Delta, vF, lambdaL, T, betaExtra, nameE2D, nameM, nameS)
    Wave B_T
    Variable L, W, Delta, vF, lambdaL, T
    Variable betaExtra
    String nameE2D, nameM, nameS
    
    Variable nB = numpnts(B_T)
    if (nB <= 0)
        Abort "Solve_AllBranches_SNS_dGSJ_betaExtra: B_T has no points."
    endif

    // ---------- 1. Scan all branches at starting field B0 ----------

    Variable B0    = B_T[0]
    Variable tau0  = SNS_tau_eVs(L, vF)
    Variable beta0 = SNS_beta(B0, L, W, lambdaL) + betaExtra

    // Parameter scan range in m. Increase if needed.
    Variable mMax = 10

    Make/FREE/D/N=0 Elist, Mlist, Slist

    Variable m, sSign
    Variable tolLocal = 1e-9
    Variable maxI     = 200

    for (m = -mMax; m <= mMax; m += 1)
        for (sSign = -1; sSign <= 1; sSign += 2)

            Make/FREE/D/N=6 pw
            pw[0] = Delta
            pw[1] = tau0
            pw[2] = beta0
            pw[3] = T
            pw[4] = m
            pw[5] = sSign

            Variable eps = 1e-9
            Variable lo  = -Delta + eps
            Variable hi  =  Delta - eps

            FindRoots /Q /B=1 /T=(tolLocal) /I=(maxI) /L=(lo) /H=(hi) dGSJ_BranchFunc, pw

            if (V_flag == 0 && numtype(V_Root) == 0)
                if (abs(V_Root) < Delta && abs(V_YatRoot) < 1e-4)
                    Redimension/N=(numpnts(Elist)+1) Elist, Mlist, Slist
                    Elist[numpnts(Elist)-1] = V_Root
                    Mlist[numpnts(Mlist)-1] = m
                    Slist[numpnts(Slist)-1] = sSign
                endif
            endif

        endfor
    endfor

    if (numpnts(Elist) == 0)
        SNS_Log("Solve_AllBranches_SNS_dGSJ_betaExtra: no subgap roots found at B[0].", level="ERR")
        return 0
    endif

    // ---------- 2. Sort branches by |E(B0)| and prepare outputs ----------

    Variable nBranches = numpnts(Elist)

    Make/FREE/D/N=(nBranches) absE
    absE = abs(Elist)
    Sort absE, Elist, Mlist, Slist

    // E(B, branch) 2D wave
    Make/O/D/N=(nB, nBranches) $nameE2D
    Wave E2D = $nameE2D
    E2D = NaN

    // Store (m,s) labels for each branch
    Make/O/D/N=(nBranches) $nameM, $nameS
    Wave mBranch = $nameM
    Wave sBranch = $nameS
    mBranch = Mlist
    sBranch = Slist

    // Set row 0 = B0 energies
    Variable j
    for (j = 0; j < nBranches; j += 1)
        E2D[0][j] = Elist[j]
    endfor

    // ---------- 3. Track each branch across all B ----------

    Variable i, B, Enew

    for (j = 0; j < nBranches; j += 1)
        m     = mBranch[j]
        sSign = sBranch[j]

        for (i = 1; i < nB; i += 1)
            B    = B_T[i]
            Enew = SolveBranchAtB_dGSJ_betaExtra(B, L, W, Delta, vF, lambdaL, T, \
                                                 betaExtra, m, sSign, tolLocal, maxI)
            E2D[i][j] = Enew
        endfor
    endfor

    return nBranches
End



// Track a single ABS branch over B, specified by "branchIndex":
// - At B[0], we scan over a small range of (m, sSign) combinations,
//   collect all subgap roots, sort them by |E|, and select the
//   "branchIndex"-th one.
// - For all subsequent B, we keep that (m, sSign) fixed and just
//   rerun FindRoots for that phase branch.
//
// branchIndex = 0 → smallest |E| at B[0] ("ground" ABS)
// branchIndex = 1 → next one, etc.

//==============================================================================
// TrackBranchOverB_dGSJ
//
// Purpose:
//   Solve subgap bound-state energies and/or track branches.
//
// Inputs:
//   B_T : [T]
//   L : input
//   W : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   T : input
//   branchIndex : input
//   tol : input
//   maxIters : input
//   nameE : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function TrackBranchOverB_dGSJ(B_T, L, W, Delta, vF, lambdaL, T, branchIndex, tol, maxIters, nameE)
    Wave B_T
    Variable L, W, Delta, vF, lambdaL, T
    Variable branchIndex
    Variable tol, maxIters
    String nameE

    Variable nB = numpnts(B_T)
    if (nB <= 0)
        Abort "TrackBranchOverB_dGSJ: B_T has no points."
    endif

    // Output wave
    Make/O/D/N=(nB) $nameE
    Wave Ebranch = $nameE
    Ebranch = NaN

    // ---------- 1. Branch identification at B[0] ----------

    Variable B0    = B_T[0]
    Variable tau0  = SNS_tau_eVs(L, vF)
    Variable beta0 = SNS_beta(B0, L, W, lambdaL)

    // Search window over (m, sSign)
    Variable mMax = 5

    Make/FREE/D/N=0 Elist, Mlist, Slist

    Variable m, sSign
    Variable tolLocal = (tol > 0) ? tol : 1e-9
    Variable maxI     = (maxIters > 0) ? maxIters : 200

    for (m = -mMax; m <= mMax; m += 1)
        for (sSign = -1; sSign <= 1; sSign += 2)

            Make/FREE/D/N=6 pw
            pw[0] = Delta
            pw[1] = tau0
            pw[2] = beta0
            pw[3] = T
            pw[4] = m
            pw[5] = sSign

            Variable eps = 1e-9
            Variable lo  = -Delta + eps
            Variable hi  =  Delta - eps

            FindRoots /Q /B=1 /T=(tolLocal) /I=(maxI) /L=(lo) /H=(hi) dGSJ_BranchFunc, pw

            if (V_flag == 0 && numtype(V_Root) == 0)
                if (abs(V_Root) < Delta && abs(V_YatRoot) < 1e-4)
                    Redimension/N=(numpnts(Elist)+1) Elist, Mlist, Slist
                    Elist[numpnts(Elist)-1] = V_Root
                    Mlist[numpnts(Mlist)-1] = m
                    Slist[numpnts(Slist)-1] = sSign
                endif
            endif

        endfor
    endfor

    if (numpnts(Elist) == 0)
        SNS_Log("Solve_AllBranches_SNS_dGSJ: no subgap roots found at B[0].", level="ERR")
        return 0
    endif

    // Sort by |E|
    Make/FREE/D/N=(numpnts(Elist)) absE
    absE = abs(Elist)
    Sort absE, Elist, Mlist, Slist

    if (branchIndex < 0 || branchIndex >= numpnts(Elist))
        SNS_Log("TrackBranchOverB_dGSJ: branchIndex " + num2istr(branchIndex) + \
                " out of range; only " + num2istr(numpnts(Elist)) + " roots.", level="ERR")
        return 0
    endif

    Variable mSel = Mlist[branchIndex]
    Variable sSel = Slist[branchIndex]
    Variable E0   = Elist[branchIndex]
    Ebranch[0]    = E0

    // ---------- 2. Sweep over B for that (mSel, sSel) ----------

    Variable i, B, Enew
    for (i = 1; i < nB; i += 1)
        B    = B_T[i]
        Enew = SolveBranchAtB_dGSJ(B, L, W, Delta, vF, lambdaL, T,   \
                                   mSel, sSel, tolLocal, maxI)
        Ebranch[i] = Enew
    endfor

    return 0
End


//============================================================
// 5. HIGH-LEVEL SOLVER API (WHAT YOU CALL FROM NOTEBOOK)
//============================================================

// Solve a single ABS branch for given transparency T
// branchIndex = 0 → lowest |E| at B[0], 1 → next, etc.

//==============================================================================
// Solve_Ebranch_SNS_dGSJ
//
// Purpose:
//   Solve subgap bound-state energies and/or track branches.
//
// Inputs:
//   B_T : [T]
//   L : input
//   W : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   T : input
//   branchIndex : input
//   nameE : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Solve_Ebranch_SNS_dGSJ(B_T, L, W, Delta, vF, lambdaL, T, branchIndex, nameE)
    WAVE B_T
    Variable L, W, Delta, vF, lambdaL, T
    Variable branchIndex
    String nameE

    Variable tolX     = 1e-9
    Variable maxIters = 100

    TrackBranchOverB_dGSJ(B_T, L, W, Delta, vF, lambdaL, T, branchIndex, \
                          tolX, maxIters, nameE)
End


// Solve all subgap branches and return an E(B, branch) matrix.
// Rows  : B index (same length as B_T)
// Cols  : branch index, sorted by |E(B0)| at the starting field B_T[0]
//
// Outputs:
//   nameE2D : 2D wave [nB x nBranches] with energies
//   nameM   : 1D wave [nBranches] with m index for each branch
//   nameS   : 1D wave [nBranches] with sSign (+1/-1) for each branch
//
// Return value = number of branches found.

//==============================================================================
// Solve_AllBranches_SNS_dGSJ
//
// Purpose:
//   Solve subgap bound-state energies and/or track branches.
//
// Inputs:
//   B_T : [T]
//   L : input
//   W : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   T : input
//   nameE2D : input
//   nameM : input
//   nameS : input
//
// Outputs:
//   return : numeric
//
// Notes:
//   Convenience wrapper calling Solve_AllBranches_SNS_dGSJ_betaExtra with betaExtra=0.
//
//==============================================================================
Function Solve_AllBranches_SNS_dGSJ(B_T, L, W, Delta, vF, lambdaL, T, nameE2D, nameM, nameS)
    Wave B_T
    Variable L, W, Delta, vF, lambdaL, T
    String nameE2D, nameM, nameS
    
    Variable nB = numpnts(B_T)
    if (nB <= 0)
        Abort "Solve_AllBranches_SNS_dGSJ: B_T has no points."
    endif

    // ---------- 1. Scan all branches at starting field B0 ----------

    Variable B0    = B_T[0]
    Variable tau0  = SNS_tau_eVs(L, vF)
    Variable beta0 = SNS_beta(B0, L, W, lambdaL)

    // Parameter scan range in m. Increase if needed.
    Variable mMax = 10

    Make/FREE/D/N=0 Elist, Mlist, Slist

    Variable m, sSign
    Variable tolLocal = 1e-9
    Variable maxI     = 200

    for (m = -mMax; m <= mMax; m += 1)
        for (sSign = -1; sSign <= 1; sSign += 2)

            Make/FREE/D/N=6 pw
            pw[0] = Delta
            pw[1] = tau0
            pw[2] = beta0
            pw[3] = T
            pw[4] = m
            pw[5] = sSign

            Variable eps = 1e-9
            Variable lo  = -Delta + eps
            Variable hi  =  Delta - eps

            FindRoots /Q /B=1 /T=(tolLocal) /I=(maxI) /L=(lo) /H=(hi) dGSJ_BranchFunc, pw

            if (V_flag == 0 && numtype(V_Root) == 0)
                if (abs(V_Root) < Delta && abs(V_YatRoot) < 1e-4)
                    Redimension/N=(numpnts(Elist)+1) Elist, Mlist, Slist
                    Elist[numpnts(Elist)-1] = V_Root
                    Mlist[numpnts(Mlist)-1] = m
                    Slist[numpnts(Slist)-1] = sSign
                endif
            endif

        endfor
    endfor

    if (numpnts(Elist) == 0)
        Print "Solve_AllBranches_SNS_dGSJ: no subgap roots found at B[0]."
        return 0
    endif

    // ---------- 2. Sort branches by |E(B0)| and prepare outputs ----------

    Variable nBranches = numpnts(Elist)

    Make/FREE/D/N=(nBranches) absE
    absE = abs(Elist)
    Sort absE, Elist, Mlist, Slist

    // E(B, branch) 2D wave
    Make/O/D/N=(nB, nBranches) $nameE2D
    Wave E2D = $nameE2D
    E2D = NaN

    // Store (m,s) labels for each branch
    Make/O/D/N=(nBranches) $nameM, $nameS
    Wave mBranch = $nameM
    Wave sBranch = $nameS
    mBranch = Mlist
    sBranch = Slist

    // Set row 0 = B0 energies
    Variable j
    for (j = 0; j < nBranches; j += 1)
        E2D[0][j] = Elist[j]
    endfor

    // ---------- 3. Track each branch across all B ----------

    Variable i, B, Enew

    for (j = 0; j < nBranches; j += 1)
        m = mBranch[j]
        sSign = sBranch[j]

        for (i = 1; i < nB; i += 1)
            B    = B_T[i]
            Enew = SolveBranchAtB_dGSJ(B, L, W, Delta, vF, lambdaL, T, \
                                       m, sSign, tolLocal, maxI)
            E2D[i][j] = Enew
        endfor
    endfor

    return nBranches
End

// Solve the lowest-energy branch and its particle-hole partner
// If you want strictly symmetric ±E, you can just mirror by hand;
// here we track only the positive-energy (or lowest-|E|) branch and
// then define E_minus = -E_plus by construction.

//==============================================================================
// Solve_Epair_SNS_dGSJ
//
// Purpose:
//   Solve subgap bound-state energies and/or track branches.
//
// Inputs:
//   B_T : [T]
//   L : input
//   W : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   T : input
//   namePlus : input
//   nameMinus : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Solve_Epair_SNS_dGSJ(B_T, L, W, Delta, vF, lambdaL, T, namePlus, nameMinus)
    WAVE B_T
    Variable L, W, Delta, vF, lambdaL, T
    String namePlus, nameMinus

    // Lowest-|E| branch
    Solve_Ebranch_SNS_dGSJ(B_T, L, W, Delta, vF, lambdaL, T, 0, namePlus)
    WAVE Eplus = $namePlus

    Make/O/D/N=(numpnts(Eplus)) $nameMinus
    WAVE Eminus = $nameMinus

    Eminus = -Eplus
End

//============================================================
// 6. DIAGNOSTICS & PLOTTING HELPERS
//============================================================

// Plot F(E) vs E at a fixed B for debugging roots

//==============================================================================
// Solve_AllBranches_SNS_dGSJ_FullScan
//
// Purpose:
//   Solve subgap bound-state energies and/or track branches.
//
// Inputs:
//   B_T : [T]
//   L : input
//   W : input
//   Delta : [eV]
//   vF : input
//   lambdaL : [m]
//   T : input
//   mMax : input
//   nameE2D : input
//   nameM : input
//   nameS : input
//   [tol, maxIters] :
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Solve_AllBranches_SNS_dGSJ_FullScan(B_T, L, W, Delta, vF, lambdaL, T, mMax,nameE2D, nameM, nameS, [tol, maxIters])
    Wave   B_T
    Variable L, W, Delta, vF, lambdaL, T
    Variable mMax
    String nameE2D, nameM, nameS
    Variable tol, maxIters

    Variable nB = numpnts(B_T)
    if (nB <= 0)
        Abort "Solve_AllBranches_SNS_dGSJ_FullScan: B_T has no points."
    endif
    if (mMax < 0)
        Abort "Solve_AllBranches_SNS_dGSJ_FullScan: mMax must be >= 0."
    endif

    // Root-finder controls
    Variable tolLocal = (ParamIsDefault(tol)      || tol      <= 0) ? 1e-9 : tol
    Variable maxI     = (ParamIsDefault(maxIters) || maxIters <= 0) ? 200  : maxIters

    // Total number of (m,sSign) combinations:
    // m = -mMax..+mMax  (2*mMax+1 values)
    // sSign = -1,+1     (factor 2)
    Variable nBranchesMax = (2*mMax + 1)*2

    // ---------- Allocate output waves ----------

    Make/O/D/N=(nB, nBranchesMax) $nameE2D
    Wave E2D = $nameE2D
    E2D = NaN

    Make/O/D/N=(nBranchesMax) $nameM
    Make/O/D/N=(nBranchesMax) $nameS
    Wave mBranch = $nameM
    Wave sBranch = $nameS

    // Fill (m,sSign) mapping for each branch index j
    Variable j, m, sSign
    j = 0
    for (m = -mMax; m <= mMax; m += 1)
        for (sSign = -1; sSign <= 1; sSign += 2)
            if (j >= nBranchesMax)
                Abort "Solve_AllBranches_SNS_dGSJ_FullScan: internal index overflow."
            endif
            mBranch[j] = m
            sBranch[j] = sSign
            j += 1
        endfor
    endfor

    // ---------- Main solve loop over all B and all branches ----------

    Variable iB, Bval, Eroot

    for (iB = 0; iB < nB; iB += 1)

        Bval = B_T[iB]

        // For each (m,sSign) pair
        for (j = 0; j < nBranchesMax; j += 1)

            m     = mBranch[j]
            sSign = sBranch[j]

            // Attempt to find a subgap root for this (B, m, sSign)
            Eroot = SolveBranchAtB_dGSJ(Bval, L, W, Delta, vF, lambdaL, T, \
                                        m, sSign, tolLocal, maxI)

            // SolveBranchAtB_dGSJ returns NaN if no acceptable root
            E2D[iB][j] = Eroot

        endfor
    endfor

    // Return total number of branch slots used (all (m,sSign) pairs)
    return nBranchesMax
End

// Weight for trajectory at angle theta when screening current is along y.
// DeltaInd : induced gap at B=0 (eV)
// vF       : Fermi velocity (m/s)
// lambdaL  : London penetration depth (m)
// B        : magnetic field (T)
