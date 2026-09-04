#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma ModuleName=SNS_UI
#include "SNS_Core"
#include "SNS_Solver"
#include "SNS_Maps"
#include "SNS_GeometryFromMask"

//==============================================================================
// SNS_UI.ipf
//
// Plotting + diagnostics helpers for SNS ABS / DOS workflows.
//
// Keep this thin; anything that affects numerical results should live in
// SNS_Solver.ipf or SNS_Maps.ipf.
//==============================================================================


//==============================================================================
// Plot_AllBranches_SNS
//
// Purpose:
//   Plot diagnostic or result waves.
//
// Inputs:
//   nameE2D : input
//   nameB : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Plot_AllBranches_SNS(nameE2D, nameB)
    String nameE2D, nameB

    Wave E2D = $nameE2D    // 2D wave: [iB][jBranch]
    Wave B_T = $nameB      // 1D wave: magnetic field

    Variable nB        = DimSize(E2D, 0)
    Variable nBranches = DimSize(E2D, 1)
    if (nB != numpnts(B_T))
        Abort "Plot_AllBranches_SNS: size mismatch between E2D and B_T."
    endif

    Variable j

    // Make a new graph and plot all columns vs B_T
    // Display 
    for (j = 0; j < nBranches; j += 1)
        AppendToGraph E2D[][j] vs B_T
    endfor

    ModifyGraph mode=4        // lines
    ModifyGraph marker=0
    ModifyGraph rgb=(0,0,0)   // all black; tweak if you want colors
End

// Move waves with names "prefix + chXXX" into a subfolder and
// rename them to just "chXXX" (e.g. "E_allBranches_ch000" -> :E_allBranches:ch000)

//==============================================================================
// MoveBranchWavesToSubfolder
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   prefix : input
//   destFolder : input
//
// Outputs:
//   return : numeric
//   destFolder : output wave(s)
//
//==============================================================================
Function MoveBranchWavesToSubfolder(prefix, destFolder)
    String prefix, destFolder

    String pattern = prefix + "ch*"
    String wList   = WaveList(pattern, ";", "")
    Variable nWaves = ItemsInList(wList)
    Variable i
    String srcName, shortName

    for (i = 0; i < nWaves; i += 1)
        srcName = StringFromList(i, wList)

        // strip the prefix, keep "ch000"
        shortName = srcName[ strlen(prefix), Inf ]

        // correct subfolder reference: start with ":"
        Duplicate/O $srcName, $(":"+destFolder+":"+shortName)
        KillWaves/Z $srcName
    endfor

    return 0
End


//============================================================
// Estimate energy-axis offset from experimental dI/dV(E,B) map
// by fitting a Gaussian peak in the summed spectrum.
//
// expEB   : 2D wave, dim0 = E (meV), dim1 = B (e.g. mT)
// Returns : E_offset (same units as E axis, typically meV)
//
// Method:
//   1) Select B in [-epsB, +epsB].
//   2) Sum over B → Spec(E).
//   3) Find peak location (V_maxloc).
//   4) Restrict to E in [Epk-epsE, Epk+epsE].
//   5) Fit Spec(E) with Gaussian; return fitted center x0.
//============================================================

//==============================================================================
// AppendAllBranchesToGraph
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   E2D : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function AppendAllBranchesToGraph(E2D)
    Wave E2D

    // B-axis from settings (simulation in T)
    Wave B_T = root:SNS_Settings:B_T

    String gname = WinName(0, 1, 1)
    if (strlen(gname) == 0)
        Abort "AppendAllBranchesToGraph: no graph window found."
    endif

    Variable nBranches = DimSize(E2D, 1)
    Variable nB        = DimSize(E2D, 0)

    //--------------------------------------------------------
    // 1) Inspect the first image in the graph to get axis units
    //--------------------------------------------------------
    String imgList = ImageNameList(gname, ";")
    String xUnitImg = ""
    String yUnitImg = ""

    if (ItemsInList(imgList) > 0)
        String imgTrace = StringFromList(0, imgList)
        Wave imgW = ImageNameToWaveRef(gname, imgTrace)
        // dim 0 -> x-axis (energy), dim 1 -> y-axis (B) for 2D image
        xUnitImg = WaveUnits(imgW, 0)
        yUnitImg = WaveUnits(imgW, 1)
    endif

    // Simulation units (assumed; override if you actually set them)
    String simBunit = WaveUnits(B_T, 0)
    if (strlen(simBunit) == 0)
        simBunit = "T"
    endif

    String simEunit = WaveUnits(E2D, -1)     // data units of E2D
    if (strlen(simEunit) == 0)
        simEunit = "eV"
    endif

    //--------------------------------------------------------
    // 2) Decide scale factors to match the image
    //--------------------------------------------------------
    Variable scaleB = 1      // B_T  (T -> mT)
    Variable scaleE = 1      // E2D  (eV -> meV)

    // B axis: simulation T -> image mT
    if (CmpStr(simBunit, "T", 1) == 0 && CmpStr(yUnitImg, "mT", 1) == 0)
        scaleB = 1e3
    endif

    // E axis: simulation eV -> image meV
    if (CmpStr(simEunit, "eV", 1) == 0 && (CmpStr(xUnitImg, "meV", 1) == 0 || CmpStr(xUnitImg, "mV", 1) == 0 ))
        scaleE = 1e3
    endif

    //--------------------------------------------------------
    // 3) Build plot waves (scaled copies) in E2D's data folder
    //--------------------------------------------------------
    String df = GetWavesDataFolder(E2D, 1)
    String nameE2D = NameOfWave(E2D)

    // B-plot wave
    Make/O/D/N=(nB) $(df + "B_T_plot")
    Wave B_plot = $(df + "B_T_plot")
    B_plot = B_T * scaleB
    // match image B-unit if we have it
    if (strlen(yUnitImg) > 0)
        SetScale/P x, DimOffset(B_T, 0)*scaleB, DimDelta(B_T, 0)*scaleB, yUnitImg, B_plot
    endif

    // E-plot wave (2D copy, scaled)
    Make/O/D/N=(nB, nBranches) $(df + nameE2D + "_plot")
    Wave E2D_plot = $(df + nameE2D + "_plot")
    E2D_plot = E2D * scaleE
    // set data units for energies (y-axis label)
    if (strlen(xUnitImg) > 0)
        SetScale d, 0, 0, xUnitImg, E2D_plot     // data-units = xUnitImg (meV)
    else
        SetScale d, 0, 0, simEunit, E2D_plot
    endif

    //--------------------------------------------------------
    // 4) Append branches that actually contain data
    //--------------------------------------------------------
    Variable j, i
    Variable hasData

    for (j = 0; j < nBranches; j += 1)
        hasData = 0
        for (i = 0; i < nB; i += 1)
            if (numtype(E2D_plot[i][j]) == 0)
                hasData = 1
                break
            endif
        endfor

        if (hasData)
            // Correct order is: Ywave vs Xwave
            AppendToGraph /W=$gname B_plot vs E2D_plot[*][j]
        endif
    endfor

    ModifyGraph/W=$gname mode=2
End

//==============================================================================
// SNS_AppendSpecialBranches
//
// Purpose:
//   Optionally display an image with SNS_DisplayWithScales(...), then append:
//
//      L_N cluster max
//      L_N cluster median
//      W_eff cluster max
//      W_eff cluster median
//
//   Uses existing variable names from SNS_RayTrace_Hist_LN_Weff.
//   No variable renaming.
//
// Display convention:
//      L_N   : red
//      W_eff : black
//      max   : solid/thick
//      median: dashed/thin
//
// Important:
//   dirPath points to the folder containing RayTraceHist / RayTraceHist_3D.
//   branchPath optionally points to E_allBranches / E_allBranches_3D.
//   If branchPath is omitted, the original behavior is used:
//        branch folder = ":E_allBranches"
//   i.e. relative to the current data folder.
//
// Legend:
//   Uses Igor's built-in Legend with \s(traceName), so the legend displays the
//   actual plotted line color/style dynamically.
//==============================================================================
Function SNS_AppendSpecialBranches(dirPath, [is3D, makeRayTrace, branchPath, img, cmap, filetype, w_displaySize_pt, Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy, H_nm, maxPath_nm])
    String dirPath
    Variable is3D, makeRayTrace
    String branchPath
    Wave/Z img
    String cmap
    Variable filetype
    Wave/Z w_displaySize_pt
    Variable Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy, H_nm, maxPath_nm

    String oldDF = GetDataFolder(1)

    String d = dirPath
    Variable nd = strlen(d)
    if (nd > 0)
        if (!CmpStr(d[nd-1], ":"))
            d = d[0, nd-2]
        endif
    endif

    if (!DataFolderExists(d))
        Abort "SNS_AppendSpecialBranches: directory does not exist: " + d
    endif

    // ---------------------------
    // Defaults
    // ---------------------------
    if (ParamIsDefault(H_nm))
        H_nm = 0
    endif

    if (ParamIsDefault(is3D))
        is3D = (H_nm > 0)
    endif

    if (ParamIsDefault(makeRayTrace))
        makeRayTrace = 0
    endif

    if (ParamIsDefault(cmap))
        cmap = ""
    endif

    if (ParamIsDefault(filetype))
        filetype = 0
    endif

    if (ParamIsDefault(Bangle_deg))
        Bangle_deg = 225
    endif

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    if (ParamIsDefault(BTK_barrier))
        BTK_barrier = params.BTK_barrier
    endif
    if (ParamIsDefault(STSx))
        STSx = -4
    endif
    if (ParamIsDefault(STSy))
        STSy = -50
    endif
    if (ParamIsDefault(Vortexx))
        Vortexx = STSx
    endif
    if (ParamIsDefault(Vortexy))
        Vortexy = STSy
    endif
    if (ParamIsDefault(maxPath_nm))
        maxPath_nm = 1000
    endif

    // ---------------------------
    // Optional ray tracing
    // ---------------------------
    if (makeRayTrace)
        if (is3D && H_nm <= 0)
            SetDataFolder $oldDF
            Abort "SNS_AppendSpecialBranches: is3D=1 with makeRayTrace=1 requires H_nm > 0."
        endif

        SNS_RayTrace_Hist_LN_Weff(d, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, Vortexx=Vortexx, Vortexy=Vortexy, H_nm=H_nm, maxPath_nm=maxPath_nm, doDisplay=0)
    endif

    // ---------------------------
    // Optional map display
    // ---------------------------
    String gname = ""

    if (!ParamIsDefault(img) && WaveExists(img))
        if (!ParamIsDefault(w_displaySize_pt) && WaveExists(w_displaySize_pt))
            gname = SNS_DisplayWithScales(img, cmap=cmap, filetype=filetype, w_displaySize_pt=w_displaySize_pt)
        else
            gname = SNS_DisplayWithScales(img, cmap=cmap, filetype=filetype)
        endif

        if (strlen(gname) == 0)
            SetDataFolder $oldDF
            Abort "SNS_AppendSpecialBranches: SNS_DisplayWithScales failed."
        endif

        DoWindow/F $gname
    else
        gname = WinName(0, 1, 1)
        if (strlen(gname) == 0)
            SetDataFolder $oldDF
            Abort "SNS_AppendSpecialBranches: no active graph and no image supplied."
        endif
    endif

    // ---------------------------
    // Resolve histogram folder
    // ---------------------------
    String histFolder
    if (is3D)
        histFolder = d + ":RayTraceHist_3D"
    else
        histFolder = d + ":RayTraceHist"
    endif

    if (!DataFolderExists(histFolder))
        SetDataFolder $oldDF
        Abort "SNS_AppendSpecialBranches: missing RayTraceHist folder: " + histFolder
    endif

    // ---------------------------
    // Resolve branch folder
    //
    // Default preserves original behavior:
    //   :E_allBranches relative to current data folder
    // ---------------------------
    String branchFolder

    if (ParamIsDefault(branchPath) || strlen(branchPath) == 0)
        if (is3D)
            branchFolder = ":E_allBranches_3D"
        else
            branchFolder = ":E_allBranches"
        endif
    else
        branchFolder = branchPath
        Variable nb = strlen(branchFolder)
        if (nb > 0)
            if (!CmpStr(branchFolder[nb-1], ":"))
                branchFolder = branchFolder[0, nb-2]
            endif
        endif
    endif

    if (!DataFolderExists(branchFolder))
        SetDataFolder $oldDF
        Abort "SNS_AppendSpecialBranches: missing branch folder: " + branchFolder
    endif

    // ---------------------------
    // Read selected trace indices
    // ---------------------------
    Variable idx_L_max, idx_L_med, idx_W_max, idx_W_med

    idx_L_max = round(SNS_RayDiagValue(histFolder + ":", "L", "tailCluster", "max_trace", is3D=is3D))
    idx_L_med = round(SNS_RayDiagValue(histFolder + ":", "L", "tailCluster", "median_trace", is3D=is3D))
    idx_W_max = round(SNS_RayDiagValue(histFolder + ":", "W", "tailCluster", "max_trace", is3D=is3D))
    idx_W_med = round(SNS_RayDiagValue(histFolder + ":", "W", "tailCluster", "median_trace", is3D=is3D))

    if (numtype(idx_L_max) != 0 || numtype(idx_L_med) != 0 || numtype(idx_W_max) != 0 || numtype(idx_W_med) != 0)
        SetDataFolder $oldDF
        Abort "SNS_AppendSpecialBranches: missing or invalid L/W tailCluster trace diagnostics."
    endif

    // ---------------------------
    // Append and style four branches
    // ---------------------------
    String chWaveName
    String tracesBefore, tracesAfter, tr
    Variable nBefore, nAfter
    Variable i

    String leg_L_max = ""
    String leg_L_med = ""
    String leg_W_max = ""
    String leg_W_med = ""

    // --- L_N cluster max: red, solid/thick
    tracesBefore = TraceNameList(gname, ";", 1)

    sprintf chWaveName, "%s:ch%03d", branchFolder, idx_L_max
    Wave/Z ch_L_max = $chWaveName
    if (!WaveExists(ch_L_max))
        SetDataFolder $oldDF
        Abort "SNS_AppendSpecialBranches: missing L_N cluster-max branch wave: " + chWaveName
    endif

    AppendAllBranchesToGraph(ch_L_max)

    tracesAfter = TraceNameList(gname, ";", 1)
    nBefore = ItemsInList(tracesBefore)
    nAfter  = ItemsInList(tracesAfter)

    if (nAfter > nBefore)
        leg_L_max = StringFromList(nBefore, tracesAfter)
    endif

    for (i = nBefore; i < nAfter; i += 1)
        tr = StringFromList(i, tracesAfter)
        ModifyGraph/W=$gname rgb($tr)=(65535,0,0), lsize($tr)=2, lstyle($tr)=0
    endfor

    // --- L_N cluster median: red, dashed/thin
    tracesBefore = TraceNameList(gname, ";", 1)

    sprintf chWaveName, "%s:ch%03d", branchFolder, idx_L_med
    Wave/Z ch_L_med = $chWaveName
    if (!WaveExists(ch_L_med))
        SetDataFolder $oldDF
        Abort "SNS_AppendSpecialBranches: missing L_N cluster-median branch wave: " + chWaveName
    endif

    AppendAllBranchesToGraph(ch_L_med)

    tracesAfter = TraceNameList(gname, ";", 1)
    nBefore = ItemsInList(tracesBefore)
    nAfter  = ItemsInList(tracesAfter)

    if (nAfter > nBefore)
        leg_L_med = StringFromList(nBefore, tracesAfter)
    endif

    for (i = nBefore; i < nAfter; i += 1)
        tr = StringFromList(i, tracesAfter)
        ModifyGraph/W=$gname rgb($tr)=(65535,0,0), lsize($tr)=1, lstyle($tr)=3
    endfor

    // --- W_eff cluster max: black, solid/thick
    tracesBefore = TraceNameList(gname, ";", 1)

    sprintf chWaveName, "%s:ch%03d", branchFolder, idx_W_max
    Wave/Z ch_W_max = $chWaveName
    if (!WaveExists(ch_W_max))
        SetDataFolder $oldDF
        Abort "SNS_AppendSpecialBranches: missing W_eff cluster-max branch wave: " + chWaveName
    endif

    AppendAllBranchesToGraph(ch_W_max)

    tracesAfter = TraceNameList(gname, ";", 1)
    nBefore = ItemsInList(tracesBefore)
    nAfter  = ItemsInList(tracesAfter)

    if (nAfter > nBefore)
        leg_W_max = StringFromList(nBefore, tracesAfter)
    endif

    for (i = nBefore; i < nAfter; i += 1)
        tr = StringFromList(i, tracesAfter)
        ModifyGraph/W=$gname rgb($tr)=(0,0,0), lsize($tr)=2, lstyle($tr)=0
    endfor

    // --- W_eff cluster median: black, dashed/thin
    tracesBefore = TraceNameList(gname, ";", 1)

    sprintf chWaveName, "%s:ch%03d", branchFolder, idx_W_med
    Wave/Z ch_W_med = $chWaveName
    if (!WaveExists(ch_W_med))
        SetDataFolder $oldDF
        Abort "SNS_AppendSpecialBranches: missing W_eff cluster-median branch wave: " + chWaveName
    endif

    AppendAllBranchesToGraph(ch_W_med)

    tracesAfter = TraceNameList(gname, ";", 1)
    nBefore = ItemsInList(tracesBefore)
    nAfter  = ItemsInList(tracesAfter)

    if (nAfter > nBefore)
        leg_W_med = StringFromList(nBefore, tracesAfter)
    endif

    for (i = nBefore; i < nAfter; i += 1)
        tr = StringFromList(i, tracesAfter)
        ModifyGraph/W=$gname rgb($tr)=(0,0,0), lsize($tr)=1, lstyle($tr)=3
    endfor

    // ---------------------------
    // Built-in Igor legend
    // Uses actual plotted traces, so line color/style are shown dynamically.
    // ---------------------------
    String branchLegend = ""

    if (strlen(leg_L_max) > 0)
        branchLegend += "\\s(" + leg_L_max + ") L\\BN\\M cluster max, ch" + num2str(idx_L_max) + "\r"
    endif

    if (strlen(leg_L_med) > 0)
        branchLegend += "\\s(" + leg_L_med + ") L\\BN\\M cluster median, ch" + num2str(idx_L_med) + "\r"
    endif

    if (strlen(leg_W_max) > 0)
        branchLegend += "\\s(" + leg_W_max + ") w\\Beff\\M cluster max, ch" + num2str(idx_W_max) + "\r"
    endif

    if (strlen(leg_W_med) > 0)
        branchLegend += "\\s(" + leg_W_med + ") w\\Beff\\M cluster median, ch" + num2str(idx_W_med)
    endif

    if (strlen(branchLegend) > 0)
        Legend/C/N=SNSBranchLegend/J/F=0/A=RB/X=1/Y=1/W=$gname branchLegend
    endif

//    Print "SNS_AppendSpecialBranches:"
//    Print "  graph                   = " + gname
//    Print "  histFolder              = " + histFolder
//    Print "  branchFolder            = " + branchFolder
//    Print "  L_N cluster max trace   = " + num2str(idx_L_max)
//    Print "  L_N cluster median trace= " + num2str(idx_L_med)
//    Print "  W_eff cluster max trace = " + num2str(idx_W_max)
//    Print "  W_eff cluster median    = " + num2str(idx_W_med)

    SetDataFolder $oldDF

    return 0
End
//============================================================
// 2. GEOMETRY & PHYSICS HELPERS
//============================================================

// Trajectory angle θ = atan2(W, L)

//==============================================================================
// Diag_Plot_F_vs_E
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
//   nE : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Diag_Plot_F_vs_E(B, L, W, Delta, vF, lambdaL, T, nE)
    Variable B, L, W, Delta, vF, lambdaL, T
    Variable nE

    Variable tau  = SNS_tau_eVs(L, vF)
    Variable beta = SNS_beta(B, L, W, lambdaL)

    Make/FREE/D/N=4 pw
    FillParams_dGSJ(pw, Delta, tau, beta, T)

    Variable eps = 1e-9
    Variable lo  = -Delta + eps
    Variable hi  =  Delta - eps

    Make/O/D/N=(nE) F_E, E_axis
    SetScale/I x, lo, hi, E_axis
    E_axis = x

    F_E = PhaseEq_dGSJ(pw, E_axis)

    Display F_E vs E_axis
    ModifyGraph mode=4
End

//============================================================
// SNS: field-dependent Lorentzian broadening Gamma(B)
//============================================================

// Dimensionless field: bDim = abs(Bval) / B0
// Gamma0 : base broadening (Broadening argument)
// Bval   : magnetic field [T]
// B0     : characteristic field [T]

//==============================================================================
// SegmentedWeight_vsY
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   theta : input
//   B : [T]
//   DeltaInd : [eV]
//   vF : input
//   lambdaL : [m]
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SegmentedWeight_vsY(theta, B, DeltaInd, vF, lambdaL)
    Variable theta, B, DeltaInd, vF, lambdaL

    // 1) Doppler energy scale (up to constants): δE_max ∝ vF * lambdaL * B
    Variable deltaEmax = abs(vF * lambdaL * B)   // treat as an energy scale in eV units

    // If maximal Doppler < DeltaInd → no segmentation yet → use plain SNS model
    if (deltaEmax <= DeltaInd)
        return 1        // all angles equally allowed, no extra weighting
    endif

    // 2) Directional Doppler for current || y
    Variable deltaE = deltaEmax * sin(theta)

    // If this direction is still gapped in the homogeneous SC, suppress it
    if (abs(deltaE) <= DeltaInd)
        return 0
    endif

    // 3) For directions on the Bogoliubov arcs, use a DOS-like weight
    Variable denom2 = deltaE^2 - DeltaInd^2
    if (denom2 <= 0)
        return 0
    endif

    // Local DOS at E≈0 for Doppler-shifted s-wave:
    // N(0,δE) ∝ |δE| / sqrt(δE^2 - DeltaInd^2)
    return abs(deltaE) / sqrt(denom2)
End

//============================================================
// SS'S helpers (finite-momentum pairing in S')
//============================================================
//
// Geometry assumption (IMPORTANT):
//   - S–S' interfaces are normal to transport direction x (junction length L along x).
//   - The screening superflow in S' is parallel to the interface (along y).
//   - A ballistic trajectory is labeled by angle theta measured from x.
//   - Therefore the Doppler projection is proportional to v_y = vF*sin(theta).
//
// Simplified Doppler parameterization used here:
//   δ_theta(B) [eV] = vF * (lambdaL * B) * sin(theta)
//   where lambdaL is the London penetration depth of S' (small n_s => larger lambdaL).
//
// Notes:
//   - The threshold for propagation is |δ_theta| > DeltaP.
//
//============================================================


// Doppler energy δ_theta(B) in eV for superflow parallel to interface (along y)
// theta: trajectory angle from transport direction (x), so v_y = vF*sin(theta)
// δ_theta = vF * B * lambdaL * sin(theta)   (eV; signed)


//==========================================================================================
// [LEGACY_SSpS]
// The S–S'–S / SSpS solver code was moved to: SNS_Legacy_SSpS.ipf
// Include it only if you need legacy SSpS functionality.
//==========================================================================================
//#include "SNS_Legacy_SSpS.ipf"

//==============================================================================
// SNS_Test_BetaExtraEffect
//
// Purpose:
//   Internal helper.
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_Test_BetaExtraEffect()
    STRUCT SNS_Params p
    SNS_LoadParams(p)

    // Simple B-grid with a single field point
    Make/FREE/D/N=1 B_T
    B_T[0] = 0      // test at B=0T so only betaExtra matters

    // Artificial channel: pick some reasonable numbers
    Variable L = 500e-9    // [m]
    Variable W = 200e-9    // [m]
    Variable Tch = 1       // transparency

    String nameE2D_0   = "E2D_beta0"
    String nameM_0     = "m_beta0"
    String nameS_0     = "s_beta0"
    String nameE2D_pi  = "E2D_betapi"
    String nameM_pi    = "m_betapi"
    String nameS_pi    = "s_betapi"

    // betaExtra = 0
    Variable betaExtra0 = 0
    Variable nBr0 = Solve_AllBranches_SNS_dGSJ_betaExtra(B_T, L, W, \
        p.Delta, p.vF, p.lambdaL, p.T_K, betaExtra0, \
        nameE2D_0, nameM_0, nameS_0)

    // betaExtra = pi
    Variable betaExtraPi = pi
    Variable nBrPi = Solve_AllBranches_SNS_dGSJ_betaExtra(B_T, L, W, \
        p.Delta, p.vF, p.lambdaL, p.T_K, betaExtraPi, \
        nameE2D_pi, nameM_pi, nameS_pi)

    Wave E2D_0  = $nameE2D_0
    Wave E2D_pi = $nameE2D_pi

    Print "nBr0 = ", nBr0, "   nBrPi = ", nBrPi

    Print "E(B=0, branches) for betaExtra = 0:"
    Variable j
    for (j = 0; j < DimSize(E2D_0, 1); j += 1)
        Print "  j = ", j, "  E0 = ", E2D_0[0][j]
    endfor

    Print "E(B=0, branches) for betaExtra = pi:"
    for (j = 0; j < DimSize(E2D_pi, 1); j += 1)
        Print "  j = ", j, "  E0 = ", E2D_pi[0][j]
    endfor


    return 0
End


//============================================================
// Debug: print vortex phase for the first few real channels
// in a given geometry folder.
//============================================================

//==============================================================================
// SNS_Debug_VortexPhaseForChannels
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   dataFolder : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_Debug_VortexPhaseForChannels(dataFolder)
    String dataFolder

    String savedDF = GetDataFolder(1)
    SetDataFolder dataFolder

    Wave Hit1x_List = Hit1x_List
    Wave Hit1y_List = Hit1y_List
    Wave Hit2x_List = Hit2x_List
    Wave Hit2y_List = Hit2y_List
    Wave Ax_vortex  = Ax_vortex
    Wave Ay_vortex  = Ay_vortex

    Variable Nch = numpnts(Hit1x_List)
    if (Nch <= 0)
        Print "SNS_Debug_VortexPhaseForChannels: no channels."
        SetDataFolder savedDF
        return -1
    endif

    Variable j, betaChord
    Variable nSteps = 200

    Print "=== Vortex phase per channel (betaChord / 2π) ==="
    for (j = 0; j < Nch && j < 10; j += 1)   // first 10 channels
        SNS_VortexPhaseOnChord(Ax_vortex, Ay_vortex, \
            Hit1x_List[j], Hit1y_List[j], Hit2x_List[j], Hit2y_List[j], \
            nSteps, betaChord)
        Print "ch", j, ": betaChord = ", betaChord, "  beta/(2*pi) = ", betaChord/(2*pi)
    endfor

    SetDataFolder savedDF
    return 0
End

//==============================================================================
// Trajectory Selection UI for DOS(E,B) maps
// Igor Pro 9, WindowHook-based
//==============================================================================

#pragma IgorVersion=9.00

// ---- user-tunable defaults ----
Constant kSNS_TrajSelTopN    = 20
Constant kSNS_TrajSelFracCut = 0.01

// Internal state folder
Function SNS_TrajSel_EnsureFolder()
    NewDataFolder/O root:SNS_UI
    NewDataFolder/O root:SNS_UI:TrajSel
    return 0
End

// Lorentzian kernel used in DOS builder
Function SNS_TrajSel_Lorentz(dE, Gamma)
    Variable dE, Gamma
    return (Gamma/pi) / (dE*dE + Gamma*Gamma)
End

//------------------------------------------------------------------------------
// Helper: find a wave in a base data folder OR any direct subfolder
// baseDF must end with ":" (GetWavesDataFolder(w,1) gives that)
// Returns full path string or "" if not found.
//------------------------------------------------------------------------------
Function/S SNS_FindWaveInBaseOrSubfolders(baseDF, wName)
    String baseDF, wName

    // try base folder
    String p = baseDF + wName
    Wave/Z w0 = $p
    if (WaveExists(w0))
        return p
    endif

    // try direct subfolders
    Variable i, n = CountObjects(baseDF, 4)   // 4 = data folders
    for (i = 0; i < n; i += 1)
        String sub = GetIndexedObjName(baseDF, 4, i)
        if (strlen(sub) == 0)
            continue
        endif
        String subDF = baseDF + sub + ":"
        p = subDF + wName
        Wave/Z w1 = $p
        if (WaveExists(w1))
            return p
        endif
    endfor

    return ""
End

// Attach selection to a DOS-map window.
// - dos2D is the displayed DOS wave (NE x NB, scaled x=E[eV], y=B[T])
// - geomGraph is the window name of the geometry graph where rays should be drawn
// - hit lists are in the same data folder as your channel build output
Function SNS_UI_EnableTrajectorySelection(dos2D, geomGraph, Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, wChan)
    Wave dos2D
    String geomGraph
    Wave Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, wChan

    SNS_TrajSel_EnsureFolder()

    String win = WinName(0, 1, 1)
    if (strlen(win) == 0)
        Abort "SNS_UI_EnableTrajectorySelection: no top graph."
    endif

    // Store references as window userdata
    SetWindow $win, hook(SNS_TrajSelHook)=SNS_TrajSelHook

    // DOS wave path + DOS base folder (where E_allBranches_* are stored, possibly in subfolders)
    SetWindow $win, userdata(SNS_dos2D)=GetWavesDataFolder(dos2D,2)
    SetWindow $win, userdata(SNS_dosDF)=GetWavesDataFolder(dos2D,1)

    SetWindow $win, userdata(SNS_geomGraph)=geomGraph

    SetWindow $win, userdata(SNS_hit1x)=GetWavesDataFolder(Hit1x_List,2)
    SetWindow $win, userdata(SNS_hit1y)=GetWavesDataFolder(Hit1y_List,2)
    SetWindow $win, userdata(SNS_hit2x)=GetWavesDataFolder(Hit2x_List,2)
    SetWindow $win, userdata(SNS_hit2y)=GetWavesDataFolder(Hit2y_List,2)

    SetWindow $win, userdata(SNS_wChan)=GetWavesDataFolder(wChan,2)

    // Also store channel-count sanity
    SetWindow $win, userdata(SNS_Nch)=num2str(numpnts(wChan))

    Print "SNS Trajectory Selection enabled for window: ", win
    return 0
End


// Rank contributions at (E*,B*) and store ranked list in root:SNS_UI:TrajSel
// IMPORTANT: E_allBranches_ch%03d waves are searched in dosDF and its direct subfolders.
Function SNS_TrajSel_RankAtPoint(Estar, Bstar, dos2D, wChan, dosDF)
    Variable Estar, Bstar
    Wave dos2D, wChan
    String dosDF

    SNS_TrajSel_EnsureFolder()

    Wave/Z B_T = root:SNS_Settings:B_T
    if (!WaveExists(B_T))
        Print "SNS_TrajSel_RankAtPoint: missing root:SNS_Settings:B_T"
        return -1
    endif

    NVAR/Z SNS_Broadening = root:SNS_Settings:Broadening
    if (!NVAR_Exists(SNS_Broadening))
        Print "SNS_TrajSel_RankAtPoint: missing root:SNS_Settings:Broadening"
        return -1
    endif

    Variable nB = numpnts(B_T)
    Variable Nch = numpnts(wChan)
    if (Nch <= 0)
        Print "SNS_TrajSel_RankAtPoint: wChan empty."
        return -1
    endif

    Variable iB = round(x2pnt(B_T, Bstar))
    iB = max(0, min(nB-1, iB))

    Variable iE = round(x2pnt(dos2D, Estar))
    iE = max(0, min(DimSize(dos2D,0)-1, iE))
    Variable DOSpix = dos2D[iE][iB]
    if (numtype(DOSpix) != 0 || DOSpix <= 0)
        DOSpix = NaN
    endif

    // --- keep in sync with SNS_Maps.ipf ---
    Variable SNS_B0_T        = 0.1
    Variable SNS_GammaAlpha1 = 1
    Variable SNS_GammaAlpha2 = 0.01
    Variable SNS_GammaMin    = SNS_Broadening

    Variable SNS_DeltaBmax_T   = 0.5
    Variable SNS_DeltaFracDrop = 0.2

    Variable Bval = B_T[iB]
    Variable bDimGamma = (SNS_B0_T > 0) ? abs(Bval)/SNS_B0_T : 0
    Variable gammaScale = 1 + SNS_GammaAlpha1*bDimGamma + SNS_GammaAlpha2*bDimGamma*bDimGamma
    Variable GammaB = max(SNS_Broadening * gammaScale, SNS_GammaMin)

    Variable bDimDelta = (SNS_DeltaBmax_T > 0) ? abs(Bval)/SNS_DeltaBmax_T : 0
    Variable deltaScale = 1 - SNS_DeltaFracDrop*bDimDelta*bDimDelta
    if (deltaScale < 0)
        deltaScale = 0
    endif

    Variable sumW = sum(wChan)
    if (sumW <= 0)
        Print "SNS_TrajSel_RankAtPoint: sum(wChan) <= 0."
        return -1
    endif

    Make/FREE/D/N=0 scoreAll
    Make/FREE/I/N=0 chAll, brAll

    Variable ch, nBr, br
    String ePath
    Variable foundAny = 0

    for (ch = 0; ch < Nch; ch += 1)

        if (wChan[ch] <= 0)
            continue
        endif

        // *** CORRECT folder structure ***
        ePath = SNS_EBranchesPath_ForChannel(dosDF, ch)
        if (strlen(ePath) == 0)
            continue
        endif

        Wave E_all = $ePath
        foundAny = 1

        // B-index sanity: E_all dim0 must cover iB
        if (DimSize(E_all, 0) <= iB)
            continue
        endif

        nBr = DimSize(E_all, 1)
        for (br = 0; br < nBr; br += 1)

            Variable E0 = E_all[iB][br]
            if (numtype(E0) != 0)
                continue
            endif

            Variable E0s = E0 * deltaScale
            Variable S = (wChan[ch]/sumW) * SNS_TrajSel_Lorentz(Estar - E0s, GammaB)

            Redimension/N=(numpnts(scoreAll)+1) scoreAll
            Redimension/N=(numpnts(chAll)+1)    chAll, brAll
            scoreAll[numpnts(scoreAll)-1] = S
            chAll[numpnts(chAll)-1]       = ch
            brAll[numpnts(brAll)-1]       = br
        endfor
    endfor

    if (!foundAny)
        Print "SNS_TrajSel: no E_allBranches waves found under ", dosDF + "E_allBranches:ch###:"
        return -1
    endif

    Variable nTot = numpnts(scoreAll)
    if (nTot == 0)
        Print "SNS_TrajSel: E waves found, but no finite contributors at this point (iB=", iB, ")."
        return -1
    endif

    Sort/R scoreAll, scoreAll, chAll, brAll

// Keep TopN only (disable DOSpix fraction cutoff for now)
Variable keep = min(kSNS_TrajSelTopN, nTot)


    Duplicate/O/R=[0,keep-1] scoreAll, root:SNS_UI:TrajSel:score
    Duplicate/O/R=[0,keep-1] chAll,    root:SNS_UI:TrajSel:chIdx
    Duplicate/O/R=[0,keep-1] brAll,    root:SNS_UI:TrajSel:brIdx

    Variable/G root:SNS_UI:TrajSel:curSel = 0
    Variable/G root:SNS_UI:TrajSel:lastE  = Estar
    Variable/G root:SNS_UI:TrajSel:lastB  = B_T[iB]

    return 0
End



// Draw/refresh highlighted channel ray on the geometry graph.
// Uses one persistent trace name so we can replace it cleanly.
Function SNS_TrajSel_DrawRay(geomGraph, Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, ch)
    String geomGraph
    Wave Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
    Variable ch

    if (strlen(geomGraph) == 0)
        return 0
    endif
    if (WinType(geomGraph) == 0)
        return 0
    endif

    String xwName="SNS_raySel_X"
    String ywName="SNS_raySel_Y"
    Make/O/D/N=2 $xwName, $ywName
    Wave rayX=$xwName
    Wave rayY=$ywName

    rayX[0]=Hit1x_List[ch]; rayY[0]=Hit1y_List[ch]
    rayX[1]=Hit2x_List[ch]; rayY[1]=Hit2y_List[ch]

    // Remove previous instance of this trace if present
    RemoveFromGraph/W=$geomGraph/Z $ywName

    AppendToGraph/W=$geomGraph rayY vs rayX

    // Style: thick black
    ModifyGraph/W=$geomGraph lsize($ywName)=3, rgb($ywName)=(0,0,0)

    return 0
End


// Move selection index (+1/-1) and redraw ray + print info
Function SNS_TrajSel_Cycle(delta, geomGraph, Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List)
    Variable delta
    String geomGraph
    Wave Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List

    NVAR/Z curSel = root:SNS_UI:TrajSel:curSel
    Wave/Z chIdx  = root:SNS_UI:TrajSel:chIdx
    Wave/Z brIdx  = root:SNS_UI:TrajSel:brIdx
    Wave/Z score  = root:SNS_UI:TrajSel:score

    if (!NVAR_Exists(curSel) || !WaveExists(chIdx) || !WaveExists(brIdx) || !WaveExists(score))
        Print "SNS_TrajSel_Cycle: missing state (curSel/chIdx/brIdx/score)"
        return 0
    endif

    Variable n = numpnts(chIdx)
    Print "SNS_TrajSel_Cycle: n=", n, " curSel(before)=", curSel, " delta=", delta

    if (n <= 1)
        // nothing to cycle
        return 0
    endif

    // robust wrap using modulo arithmetic
    Variable idx = round(curSel) + (delta > 0 ? 1 : -1)
    idx = mod(idx, n)
    if (idx < 0)
        idx += n
    endif
    curSel = idx

    Variable ch = chIdx[curSel]
    Variable br = brIdx[curSel]
    Variable sc = score[curSel]

    SNS_TrajSel_DrawRay(geomGraph, Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, ch)
    Print "SNS_TrajSel: sel=", curSel, " ch=", ch, " br=", br, " score=", sc

    return 0
End



// Window hook: click selects, wheel cycles.
// Safe against missing userdata (won't silently error out).
Function SNS_TrajSelHook(s)
    STRUCT WMWinHookStruct &s

    // ------------------------------------------------------------
    // 22 = mouseWheel, 23 = spinUpdate
    // Cycle only; NEVER re-rank here.
    // ------------------------------------------------------------
    if (s.eventCode == 22 || s.eventCode == 23)

        Print "wheel eventCode=", s.eventCode, " wheelDy=", s.wheelDy, " wheelDx=", s.wheelDx

        Wave/Z chIdx = root:SNS_UI:TrajSel:chIdx
        if (!WaveExists(chIdx))
            Print "wheel: chIdx missing"
            return 0
        endif
        if (numpnts(chIdx) <= 1)
            Print "wheel: nSel<=1 (nSel=", numpnts(chIdx), ")"
            return 0
        endif

        Variable dy = s.wheelDy
        if (numtype(dy) != 0 || dy == 0)
            Print "wheel: dy invalid=", dy
            return 0
        endif

        // Resolve userdata safely (avoid Wave $("") runtime errors)
        String geomGraph = GetUserData(s.winName, "", "SNS_geomGraph")
        String p1 = GetUserData(s.winName, "", "SNS_hit1x")
        String p2 = GetUserData(s.winName, "", "SNS_hit1y")
        String p3 = GetUserData(s.winName, "", "SNS_hit2x")
        String p4 = GetUserData(s.winName, "", "SNS_hit2y")

        if (strlen(geomGraph) == 0 || strlen(p1) == 0 || strlen(p2) == 0 || strlen(p3) == 0 || strlen(p4) == 0)
            Print "wheel: missing geomGraph/hit userdata (hook on wrong window?)"
            return 0
        endif

        Wave/Z Hit1x_List = $p1
        Wave/Z Hit1y_List = $p2
        Wave/Z Hit2x_List = $p3
        Wave/Z Hit2y_List = $p4
        if (!WaveExists(Hit1x_List) || !WaveExists(Hit1y_List) || !WaveExists(Hit2x_List) || !WaveExists(Hit2y_List))
            Print "wheel: hit waves not found (paths bad)"
            return 0
        endif

        Variable delta = (dy > 0) ? +1 : -1
        SNS_TrajSel_Cycle(delta, geomGraph, Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List)
        return 1
    endif


    // ------------------------------------------------------------
    // 5 = mouseUp   (rank only here)
    // ------------------------------------------------------------
    if (s.eventCode == 5)

        String dosPath   = GetUserData(s.winName, "", "SNS_dos2D")
        String dosDF     = GetUserData(s.winName, "", "SNS_dosDF")
        String wChanPath = GetUserData(s.winName, "", "SNS_wChan")
        if (strlen(dosPath) == 0 || strlen(dosDF) == 0 || strlen(wChanPath) == 0)
            Print "click: missing dos2D/dosDF/wChan userdata (hook on wrong window?)"
            return 0
        endif

        Wave/Z dos2D = $dosPath
        Wave/Z wChan = $wChanPath
        if (!WaveExists(dos2D) || !WaveExists(wChan))
            Print "click: dos2D or wChan wave missing"
            return 0
        endif

        Variable Estar = AxisValFromPixel(s.winName, "bottom", s.mouseLoc.h)
        Variable Bstar = AxisValFromPixel(s.winName, "left",   s.mouseLoc.v)

        // rank
        SNS_TrajSel_RankAtPoint(Estar, Bstar, dos2D, wChan, dosDF)

        // draw best selection (safe userdata resolution)
        String geomGraph2 = GetUserData(s.winName, "", "SNS_geomGraph")
        String q1 = GetUserData(s.winName, "", "SNS_hit1x")
        String q2 = GetUserData(s.winName, "", "SNS_hit1y")
        String q3 = GetUserData(s.winName, "", "SNS_hit2x")
        String q4 = GetUserData(s.winName, "", "SNS_hit2y")
        if (strlen(geomGraph2) == 0 || strlen(q1) == 0 || strlen(q2) == 0 || strlen(q3) == 0 || strlen(q4) == 0)
            Print "click: missing geomGraph/hit userdata"
            return 1
        endif

        Wave/Z Hit1x_List2 = $q1
        Wave/Z Hit1y_List2 = $q2
        Wave/Z Hit2x_List2 = $q3
        Wave/Z Hit2y_List2 = $q4
        if (!WaveExists(Hit1x_List2) || !WaveExists(Hit1y_List2) || !WaveExists(Hit2x_List2) || !WaveExists(Hit2y_List2))
            Print "click: hit waves not found (paths bad)"
            return 1
        endif

        NVAR/Z curSel = root:SNS_UI:TrajSel:curSel
        Wave/Z chIdx2 = root:SNS_UI:TrajSel:chIdx
        if (NVAR_Exists(curSel) && WaveExists(chIdx2) && numpnts(chIdx2) > 0)
            SNS_TrajSel_DrawRay(geomGraph2, Hit1x_List2, Hit1y_List2, Hit2x_List2, Hit2y_List2, chIdx2[curSel])
        endif

        return 1
    endif

    return 0
End

//==============================================================================
// Disable trajectory selection hook on a given DOS window
//==============================================================================
Function SNS_UI_DisableTrajectorySelection(dosWin)
    String dosWin

    if (WinType(dosWin) == 0)
        Print "SNS_UI_DisableTrajectorySelection: window not found: ", dosWin
        return -1
    endif

    SetWindow $dosWin, hook(SNS_TrajSelHook)=$""

    Print "SNS Trajectory Selection disabled on window: ", dosWin
    return 0
End

//------------------------------------------------------------------------------
// Return full path to the E_allBranches wave for channel ch.
// Folder layout: dosDF:E_allBranches: contains waves named ch000, ch001, ...
// Returns "" if not found.
//------------------------------------------------------------------------------
Function/S SNS_EBranchesPath_ForChannel(dosDF, ch)
    String dosDF
    Variable ch

    String p
    sprintf p, "%sE_allBranches:ch%03d", dosDF, ch

    Wave/Z w = $p
    if (WaveExists(w))
        return p
    endif

    return ""
End

//==============================================================================
// Trajectory Selection UI for LDOS(E,r) line maps
// Requires per-position cache from SNS_LineDOS_BFixed_FromMask.
//==============================================================================

//------------------------------------------------------------------------------
// Read line-cache folder from LDOS wave note.
//------------------------------------------------------------------------------
Function/S SNS_LineTrajSel_GetCacheFolderFromLDOS(ldos2D)
    Wave ldos2D

    String nt = note(ldos2D)
    String cache = StringByKey("SNS_LineCacheFolder", nt, "=", ";")

    if (strlen(cache) == 0)
        return ""
    endif

    return SNS_LineDOS_TrailingColon(cache)
End


//------------------------------------------------------------------------------
// Return folder for line index ir:
//   root:SNS_LineTmp:LineDOSCache_<tag>:r0000:
//------------------------------------------------------------------------------
Function/S SNS_LineTrajSel_PosFolder(lineCacheFolder, ir)
    String lineCacheFolder
    Variable ir

    String baseDF = SNS_LineDOS_TrailingColon(lineCacheFolder)
    String outPosDF
    sprintf outPosDF, "%sr%04d:", baseDF, ir

    return outPosDF
End


//------------------------------------------------------------------------------
// Enable trajectory selection on an LDOS(E,r) graph.
//
// Usage:
//   Display; AppendImage LDOS_Conv_B40mT_0deg
//   SNS_UI_EnableLineTrajectorySelection(LDOS_Conv_B40mT_0deg, "GeomGraph0")
//------------------------------------------------------------------------------
Function SNS_UI_EnableLineTrajectorySelection(ldos2D, geomGraph, [lineCacheFolder])
    Wave ldos2D
    String geomGraph
    String lineCacheFolder

    SNS_TrajSel_EnsureFolder()

    String win = WinName(0, 1, 1)
    if (strlen(win) == 0)
        Abort "SNS_UI_EnableLineTrajectorySelection: no top graph."
    endif

    String cacheDF
    if (ParamIsDefault(lineCacheFolder))
        cacheDF = SNS_LineTrajSel_GetCacheFolderFromLDOS(ldos2D)
    else
        cacheDF = SNS_LineDOS_TrailingColon(lineCacheFolder)
    endif

    if (strlen(cacheDF) == 0)
        Abort "SNS_UI_EnableLineTrajectorySelection: could not determine SNS_LineCacheFolder from LDOS wave note."
    endif

    cacheDF = SNS_LineDOS_TrailingColon(cacheDF)

    DFREF dfrCache = $cacheDF
    if (!DataFolderRefStatus(dfrCache))
        Abort "SNS_UI_EnableLineTrajectorySelection: line cache folder does not exist: " + cacheDF
    endif

    SetWindow $win, hook(SNS_LineTrajSelHook)=SNS_LineTrajSelHook

    SetWindow $win, userdata(SNS_ldos2D)=GetWavesDataFolder(ldos2D, 2)
    SetWindow $win, userdata(SNS_lineCacheFolder)=cacheDF
    SetWindow $win, userdata(SNS_geomGraph)=geomGraph
    SetWindow $win, userdata(SNS_activePosFolder)=""

    Print "SNS Line Trajectory Selection enabled for window: ", win
    Print "  LDOS wave: ", GetWavesDataFolder(ldos2D, 2)
    Print "  cache:     ", cacheDF

    return 0
End


//------------------------------------------------------------------------------
// Rank contributing trajectories at LDOS point (E*, ir).
//
// Fixed-B line DOS: branch-energy B index is always iB = 0.
//
// Selection logic:
//   A branch contributes only if
//
//       abs(Estar - E0) <= GammaB
//
//   where GammaB is:
//     - root:SNS_Settings:Broadening for raw LDOS maps
//     - effective thermal/modulation convolution width for LDOS_Conv_* maps
//------------------------------------------------------------------------------ 
Function SNS_LineTrajSel_RankAtPoint(Estar, ir, ldos2D, inputPosFolder)
    Variable Estar, ir
    Wave ldos2D
    String inputPosFolder

    SNS_TrajSel_EnsureFolder()

    String rankPosDF = SNS_LineDOS_TrailingColon(inputPosFolder)

    Wave/Z wChan = $(rankPosDF + "wChan")
    if (!WaveExists(wChan))
        Print "SNS_LineTrajSel_RankAtPoint: missing wChan in ", rankPosDF
        return -1
    endif

    Variable GammaB = SNS_TrajSel_GetEffectiveEnergyWidth(ldos2D)
    if (numtype(GammaB) != 0 || GammaB <= 0)
        Print "SNS_LineTrajSel_RankAtPoint: invalid effective energy width."
        return -1
    endif

    Variable Nch = numpnts(wChan)
    if (Nch <= 0)
        Print "SNS_LineTrajSel_RankAtPoint: wChan empty."
        return -1
    endif

    Variable sumW = sum(wChan)
    if (sumW <= 0)
        Print "SNS_LineTrajSel_RankAtPoint: sum(wChan) <= 0."
        return -1
    endif

    Variable iB = 0

    Make/FREE/D/N=0 scoreAll
    Make/FREE/D/N=0 dEAll
    Make/FREE/D/N=0 E0All
    Make/FREE/I/N=0 chAll
    Make/FREE/I/N=0 brAll

    Variable ch, br, nBr
    String ePath
    Variable foundAny = 0

    for (ch = 0; ch < Nch; ch += 1)

        if (wChan[ch] <= 0)
            continue
        endif

        sprintf ePath, "%sE_allBranches:ch%03d", rankPosDF, ch
        Wave/Z E_all = $ePath
        if (!WaveExists(E_all))
            continue
        endif

        foundAny = 1

        if (DimSize(E_all, 0) <= iB)
            continue
        endif

        nBr = DimSize(E_all, 1)

        for (br = 0; br < nBr; br += 1)

            Variable E0 = E_all[iB][br]
            if (numtype(E0) != 0)
                continue
            endif

            Variable dE = Estar - E0

            // Hard spectral-window selection.
            if (abs(dE) > GammaB)
                continue
            endif

            // Rank accepted states by their Lorentzian contribution.
            Variable S = (wChan[ch] / sumW) * SNS_TrajSel_Lorentz(dE, GammaB)

            Redimension/N=(numpnts(scoreAll)+1) scoreAll
            Redimension/N=(numpnts(dEAll)+1)    dEAll
            Redimension/N=(numpnts(E0All)+1)    E0All
            Redimension/N=(numpnts(chAll)+1)    chAll
            Redimension/N=(numpnts(brAll)+1)    brAll

            scoreAll[numpnts(scoreAll)-1] = S
            dEAll[numpnts(dEAll)-1]       = dE
            E0All[numpnts(E0All)-1]       = E0
            chAll[numpnts(chAll)-1]       = ch
            brAll[numpnts(brAll)-1]       = br
        endfor
    endfor

    if (!foundAny)
        Print "SNS_LineTrajSel: no E_allBranches waves found in ", rankPosDF
        return -1
    endif

    Variable nTot = numpnts(scoreAll)
    if (nTot == 0)
        Print "SNS_LineTrajSel: no states within |E-E0| <= ", GammaB, " at E=", Estar
        return -1
    endif

    Sort/R scoreAll, scoreAll, dEAll, E0All, chAll, brAll

    Variable keep = min(kSNS_TrajSelTopN, nTot)

    Duplicate/O/R=[0, keep-1] scoreAll, root:SNS_UI:TrajSel:score
    Duplicate/O/R=[0, keep-1] dEAll,    root:SNS_UI:TrajSel:dE
    Duplicate/O/R=[0, keep-1] E0All,    root:SNS_UI:TrajSel:E0
    Duplicate/O/R=[0, keep-1] chAll,    root:SNS_UI:TrajSel:chIdx
    Duplicate/O/R=[0, keep-1] brAll,    root:SNS_UI:TrajSel:brIdx

    Variable/G root:SNS_UI:TrajSel:curSel = 0
    Variable/G root:SNS_UI:TrajSel:lastE  = Estar
    Variable/G root:SNS_UI:TrajSel:lastIr = ir
    Variable/G root:SNS_UI:TrajSel:lastWidth = GammaB

    Print "SNS_LineTrajSel: selected ", keep, " of ", nTot, " states within |E-E0| <= ", GammaB
    Print "  E=", Estar, " r-index=", ir, " posFolder=", rankPosDF

    return 0
End

//------------------------------------------------------------------------------
// Hook for LDOS(E,r): click ranks at selected E,r; mouse wheel cycles.
//------------------------------------------------------------------------------
Function SNS_LineTrajSelHook(s)
    STRUCT WMWinHookStruct &s

    // Mouse wheel / spinUpdate: cycle currently ranked trajectories.
    if (s.eventCode == 22 || s.eventCode == 23)

        Wave/Z chIdxWheel = root:SNS_UI:TrajSel:chIdx
        if (!WaveExists(chIdxWheel) || numpnts(chIdxWheel) <= 1)
            return 0
        endif

        Variable dy = s.wheelDy
        if (numtype(dy) != 0 || dy == 0)
            return 0
        endif

        String geomGraphWheel = GetUserData(s.winName, "", "SNS_geomGraph")
        String activePosDF = GetUserData(s.winName, "", "SNS_activePosFolder")

        if (strlen(geomGraphWheel) == 0 || strlen(activePosDF) == 0)
            Print "SNS_LineTrajSel wheel: missing geomGraph or active position folder."
            return 0
        endif

        activePosDF = SNS_LineDOS_TrailingColon(activePosDF)

        Wave/Z Hit1x_List_Wheel = $(activePosDF + "Hit1x_List")
        Wave/Z Hit1y_List_Wheel = $(activePosDF + "Hit1y_List")
        Wave/Z Hit2x_List_Wheel = $(activePosDF + "Hit2x_List")
        Wave/Z Hit2y_List_Wheel = $(activePosDF + "Hit2y_List")

        if (!WaveExists(Hit1x_List_Wheel) || !WaveExists(Hit1y_List_Wheel) || \
            !WaveExists(Hit2x_List_Wheel) || !WaveExists(Hit2y_List_Wheel))
            Print "SNS_LineTrajSel wheel: missing hit-list waves in ", activePosDF
            return 0
        endif

        Variable delta = (dy > 0) ? +1 : -1

        SNS_TrajSel_Cycle(delta, geomGraphWheel, \
            Hit1x_List_Wheel, Hit1y_List_Wheel, Hit2x_List_Wheel, Hit2y_List_Wheel)

        return 1
    endif


    // Mouse up: select LDOS point and rank trajectories.
    if (s.eventCode == 5)

        String ldosPath = GetUserData(s.winName, "", "SNS_ldos2D")
        String cacheDF  = GetUserData(s.winName, "", "SNS_lineCacheFolder")

        if (strlen(ldosPath) == 0 || strlen(cacheDF) == 0)
            Print "SNS_LineTrajSel click: missing LDOS/cache userdata."
            return 0
        endif

        Wave/Z ldos2D = $ldosPath
        if (!WaveExists(ldos2D))
            Print "SNS_LineTrajSel click: LDOS wave missing."
            return 0
        endif

        Variable Estar = AxisValFromPixel(s.winName, "bottom", s.mouseLoc.h)
        Variable Rstar = AxisValFromPixel(s.winName, "left",   s.mouseLoc.v)

        Variable r0 = DimOffset(ldos2D, 1)
        Variable dr = DimDelta(ldos2D, 1)

        if (dr == 0)
            Print "SNS_LineTrajSel click: LDOS y-axis has zero spacing."
            return 0
        endif

        Variable irClick = round((Rstar - r0) / dr)
        irClick = max(0, min(DimSize(ldos2D, 1)-1, irClick))

        String clickPosDF = SNS_LineTrajSel_PosFolder(cacheDF, irClick)

        DFREF dfrClickPos = $clickPosDF
        if (!DataFolderRefStatus(dfrClickPos))
            Print "SNS_LineTrajSel click: position cache not found: ", clickPosDF
            return 0
        endif

        SetWindow $s.winName, userdata(SNS_activePosFolder)=clickPosDF

        Variable errClick = SNS_LineTrajSel_RankAtPoint(Estar, irClick, ldos2D, clickPosDF)
        if (errClick != 0)
            return 1
        endif

        String geomGraphClick = GetUserData(s.winName, "", "SNS_geomGraph")

        Wave/Z Hit1x_List_Click = $(clickPosDF + "Hit1x_List")
        Wave/Z Hit1y_List_Click = $(clickPosDF + "Hit1y_List")
        Wave/Z Hit2x_List_Click = $(clickPosDF + "Hit2x_List")
        Wave/Z Hit2y_List_Click = $(clickPosDF + "Hit2y_List")

        NVAR/Z curSelClick = root:SNS_UI:TrajSel:curSel
        Wave/Z chIdxClick = root:SNS_UI:TrajSel:chIdx

        if (strlen(geomGraphClick) > 0 && \
            WaveExists(Hit1x_List_Click) && WaveExists(Hit1y_List_Click) && \
            WaveExists(Hit2x_List_Click) && WaveExists(Hit2y_List_Click) && \
            NVAR_Exists(curSelClick) && WaveExists(chIdxClick) && numpnts(chIdxClick) > 0)

            SNS_TrajSel_DrawRay(geomGraphClick, \
                Hit1x_List_Click, Hit1y_List_Click, Hit2x_List_Click, Hit2y_List_Click, \
                chIdxClick[curSelClick])
        endif

        Print "SNS_LineTrajSel: E=", Estar, " r-index=", irClick, " posFolder=", clickPosDF

        return 1
    endif

    return 0
End


//------------------------------------------------------------------------------
// Disable LDOS trajectory-selection hook.
//------------------------------------------------------------------------------
Function SNS_UI_DisableLineTrajectorySelection(ldosWin)
    String ldosWin

    if (WinType(ldosWin) == 0)
        Print "SNS_UI_DisableLineTrajectorySelection: window not found: ", ldosWin
        return -1
    endif

    SetWindow $ldosWin, hook(SNS_LineTrajSelHook)=$""

    Print "SNS Line Trajectory Selection disabled on window: ", ldosWin
    return 0
End

Function SNS_TrajSel_GetEffectiveEnergyWidth(ldos2D)
    Wave ldos2D

    NVAR/Z SNS_Broadening = root:SNS_Settings:Broadening
    if (!NVAR_Exists(SNS_Broadening))
        Print "SNS_TrajSel_GetEffectiveEnergyWidth: missing root:SNS_Settings:Broadening"
        return NaN
    endif

    String wName = NameOfWave(ldos2D)

    // Raw LDOS map: use Lorentzian Gamma
    if (strsearch(wName, "Conv", 0) < 0)
        return SNS_Broadening
    endif

    // Convolved LDOS map: use effective experimental-resolution window
    NVAR/Z T_K   = root:SNS_Settings:T_K
    NVAR/Z V_mod = root:SNS_Settings:V_mod

    if (!NVAR_Exists(T_K) || !NVAR_Exists(V_mod))
        Print "SNS_TrajSel_GetEffectiveEnergyWidth: missing T_K or V_mod."
        return NaN
    endif

    Variable NE = DimSize(ldos2D, 0)
    Variable E0 = DimOffset(ldos2D, 0)
    Variable dE = DimDelta(ldos2D, 0)

    Make/FREE/D/N=(NE) E_axis
    E_axis = E0 + p*dE

    SNS_MakeThermalKernel(E_axis, T_K, "SNS_K_T_sel")
    SNS_MakeModulationKernel(E_axis, V_mod, "SNS_K_Mod_sel")

    Wave SNS_K_T_sel   = SNS_K_T_sel
    Wave SNS_K_Mod_sel = SNS_K_Mod_sel

    Variable epsK  = 1e-4
    Variable NpadT = SNS_KernelHalfWidth(SNS_K_T_sel,   epsK)
    Variable NpadM = SNS_KernelHalfWidth(SNS_K_Mod_sel, epsK)

    Variable width = max(NpadT, NpadM) * abs(dE)

    KillWaves/Z SNS_K_T_sel, SNS_K_Mod_sel

    if (width <= 0)
        width = SNS_Broadening
    endif

    return width
End
