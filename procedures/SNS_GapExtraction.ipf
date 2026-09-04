#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma DefaultTab={3,20,4}

//==============================================================================
// SNS_GapExtraction.ipf
//
// Purpose:
//   Stand-alone SNS helpers for extracting E/B cuts, locating simple spectral
//   extrema, and storing gap-closing marker values for spectroscopy analysis.
//
// Code level:
//   These functions operate only on Igor waves and their dimension scaling.
//   They do not assume a specific SNS Hamiltonian or LDOS generator.
//
// Physics level:
//   The caller decides whether the selected peaks correspond to gap edges,
//   zero-bias crossings, or another spectral feature.
//
// Provenance:
//   Ported from STMtools/SymmetrizeAndFirstGapClosing.ipf and renamed with
//   SNS_ prefixes so the SNS release notebooks do not depend on STMtools.
//==============================================================================

Function SNS_ExtractESliceAt(W, E0, [outName, energyDim])
    Wave    W
    Variable E0
    String  outName
    Variable energyDim

    if (WaveDims(W) != 2)
        Abort "SNS_ExtractESliceAtMake: W must be 2D."
    endif

    if (ParamIsDefault(energyDim))
        energyDim = 0
    endif
    if (energyDim != 0 && energyDim != 1)
        Abort "SNS_ExtractESliceAtMake: energyDim must be 0 or 1."
    endif

    Variable otherDim = 1 - energyDim

    Variable nE = DimSize(W, energyDim)
    Variable nO = DimSize(W, otherDim)

    if (nE < 2)
        Abort "SNS_ExtractESliceAtMake: energy dimension has <2 points."
    endif
    if (nO < 1)
        Abort "SNS_ExtractESliceAtMake: other dimension has <1 point."
    endif

    // map E0 to fractional index along energyDim (assumes linear scaling)
    Variable off = DimOffset(W, energyDim)
    Variable dlt = DimDelta(W, energyDim)
    if (dlt == 0)
        Abort "SNS_ExtractESliceAtMake: energyDim delta is 0."
    endif

    Variable pFrac = (E0 - off)/dlt

    if (pFrac < 0 || pFrac > (nE-1))
        Abort "SNS_ExtractESliceAtMake: E0 outside energyDim range."
    endif

    Variable iExact = round(pFrac)
    Variable tolIdx = 1e-6

    Variable iLo, iHi, t
    if (abs(pFrac - iExact) < tolIdx)
        iLo = iExact
        iHi = iExact
        t   = 0
    else
        iLo = floor(pFrac)
        iHi = iLo + 1
        if (iHi >= nE)
            Abort "SNS_ExtractESliceAtMake: cannot interpolate at upper edge."
        endif
        t = pFrac - iLo
    endif

    if (ParamIsDefault(outName) || strlen(outName) == 0)
        outName = UniqueName("Ecut", 1, 0)
    endif
    outName = CleanupName(outName, 0)

    // Make a TRUE 1D output wave of length nO (all B points)
    Make/O/N=(nO) $outName
    Wave outW = $outName

    Variable p
    if (energyDim == 0)
        // rows are energy, columns are B
        for (p=0; p<nO; p+=1)
            outW[p] = (1-t)*W[iLo][p] + t*W[iHi][p]
        endfor
    else
        // columns are energy, rows are B
        for (p=0; p<nO; p+=1)
            outW[p] = (1-t)*W[p][iLo] + t*W[p][iHi]
        endfor
    endif

    // Put the "other" dimension scale (B) onto x-axis of output
    SetScale/P x, DimOffset(W,otherDim), DimDelta(W,otherDim), WaveUnits(W,otherDim), outW
    Note outW, "Ecut at E="+num2str(E0)+"; energyDim="+num2istr(energyDim)+"; iLo="+num2istr(iLo)+"; iHi="+num2istr(iHi)+"; t="+num2str(t)
End





// Extract a 1D cut at a given B0 from a 2D wave W,
// using the dimension scaling (no axis wave needed).
// Creates the output wave in the SAME data folder as W.
//
// W:   2D wave (typically rows=E, cols=B)
// B0:  field value to extract (e.g. 0)
// bDim (optional): which dimension is B (0 or 1). Default 1 (columns).
//
// Output wave is vs the OTHER dimension (typically energy).

Function SNS_ExtractBSliceAt(W, B0, [outName, bDim])
    Wave    W
    Variable B0
    String  outName
    Variable bDim

    if (WaveDims(W) != 2)
        Abort "SNS_ExtractBSliceAtMake: W must be a 2D wave."
    endif
    if (ParamIsDefault(bDim))
        bDim = 1
    endif
    if (bDim != 0 && bDim != 1)
        Abort "SNS_ExtractBSliceAtMake: bDim must be 0 or 1."
    endif

    Variable eDim = 1 - bDim
    Variable nB   = DimSize(W, bDim)
    Variable nE   = DimSize(W, eDim)

    // map B0 to fractional index along bDim (assumes linear scaling)
    Variable off = DimOffset(W, bDim)
    Variable dlt = DimDelta(W, bDim)
    if (dlt == 0)
        Abort "SNS_ExtractBSliceAtMake: B-axis delta is 0."
    endif

    Variable pFrac = (B0 - off)/dlt
    if (pFrac < 0 || pFrac > (nB-1))
        Abort "SNS_ExtractBSliceAtMake: B0 outside B-axis range."
    endif

    Variable iExact = round(pFrac)
    Variable tolIdx = 1e-6

    Variable iLo, iHi, t
    if (abs(pFrac - iExact) < tolIdx)
        iLo = iExact
        iHi = iExact
        t   = 0
    else
        iLo = floor(pFrac)
        iHi = iLo + 1
        if (iHi >= nB)
            Abort "SNS_ExtractBSliceAtMake: cannot interpolate at upper edge."
        endif
        t = pFrac - iLo
    endif

    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR   = GetWavesDataFolderDFR(W)
    SetDataFolder wDFR

    if (ParamIsDefault(outName) || strlen(outName)==0)
        outName = UniqueName("Bcut", 1, 0)
    endif
    outName = CleanupName(outName, 0)

    Make/O/N=(nE) $outName
    Wave outW = $outName

    Variable p
    if (bDim == 1)
        // cols are B, rows are E -> outW[p=row] = W[row][col]
        for (p=0; p<nE; p+=1)
            outW[p] = (1-t)*W[p][iLo] + t*W[p][iHi]
        endfor
    else
        // rows are B, cols are E -> outW[p=col] = W[row][col]
        for (p=0; p<nE; p+=1)
            outW[p] = (1-t)*W[iLo][p] + t*W[iHi][p]
        endfor
    endif

    // Put the ENERGY axis scaling onto x of output
    SetScale/P x, DimOffset(W,eDim), DimDelta(W,eDim), WaveUnits(W,eDim), outW
    Note outW, "Extracted slice at B="+num2str(B0)+"; bDim="+num2istr(bDim)+"; iLo="+num2istr(iLo)+"; iHi="+num2istr(iHi)+"; t="+num2str(t)

    SetDataFolder oldDFR
End




// ABS(first-derivative) peak finder with optional pre-processing (smooth/interp with selectable order).
// NO return value. Outputs in SAME folder as wIn.
//
// Pre-processing (before derivative):
//   preInterpFactor = 1  -> no interpolation
//   preInterpFactor > 1  -> upsample by interpolation (linear or cubic Hermite)
//   preSmoothPts > 0     -> smooth either before or after interpolation (controlled by preSmoothAfterInterp)
//
// Post-processing (after derivative):
//   smoothPts > 0        -> smooth the abs(first-derivative) wave before peak picking
//
// Optional params:
//   outAbsDerivName    : output abs(1st-derivative) wave name (default: "<wIn>_absD")
//   peakXName          : x positions of peaks (default: "<wIn>_absD_peakX")
//   peakYName          : peak heights (default: "<wIn>_absD_peakY")
//   peakCountVarName   : global var count (default: "<wIn>_absD_nPeaks")
//   thresh             : minimum peak height in absD to accept (default 0)
//   minDistPts         : minimum separation in POINTS of the absD wave (default 1)
//   smoothPts          : smoothing of absD after derivative (default 0)
//   preSmoothPts       : smoothing of input before derivative (default 0)
//   preInterpFactor    : upsampling factor before derivative (default 1)
//   preInterpMethod    : 0=linear, 1=cubic Hermite (default 1)
//   preSmoothAfterInterp: 0=smooth->interp, 1=interp->smooth (default 0)

Function SNS_AbsDerivFindPeaks1D(wIn, [outAbsDerivName, peakXName, peakYName, peakCountVarName, thresh, minDistPts, smoothPts, preSmoothPts, preInterpFactor, preInterpMethod, preSmoothAfterInterp])
    Wave    wIn
    String  outAbsDerivName, peakXName, peakYName, peakCountVarName
    Variable thresh, minDistPts, smoothPts, preSmoothPts, preInterpFactor, preInterpMethod, preSmoothAfterInterp

    // --- checks ---
    if (WaveDims(wIn) != 1)
        Abort "SNS_AbsDerivFindPeaks1D: input must be 1D."
    endif
    Variable n0 = numpnts(wIn)
    if (n0 < 3)
        Abort "SNS_AbsDerivFindPeaks1D: input must have at least 3 points."
    endif

    // defaults
    if (ParamIsDefault(thresh))
        thresh = 0
    endif
    if (ParamIsDefault(minDistPts))
        minDistPts = 1
    endif
    if (ParamIsDefault(smoothPts))
        smoothPts = 0
    endif
    if (ParamIsDefault(preSmoothPts))
        preSmoothPts = 0
    endif
    if (ParamIsDefault(preInterpFactor))
        preInterpFactor = 1
    endif
    if (ParamIsDefault(preInterpMethod))
        preInterpMethod = 1
    endif
    if (ParamIsDefault(preSmoothAfterInterp))
        preSmoothAfterInterp = 0
    endif

    minDistPts           = max(1, round(minDistPts))
    smoothPts            = max(0, round(smoothPts))
    preSmoothPts         = max(0, round(preSmoothPts))
    preInterpFactor      = max(1, round(preInterpFactor))
    preInterpMethod      = (preInterpMethod != 0)       // 0=linear, 1=cubic
    preSmoothAfterInterp = (preSmoothAfterInterp != 0)  // 0 or 1

    // x scaling (assumes linear). If unscaled, dx=1.
    Variable xOff = DimOffset(wIn, 0)
    Variable dx   = DimDelta(wIn, 0)
    if (dx == 0)
        dx = 1
    endif
    String xUnits = WaveUnits(wIn, 0)
    String dUnits = WaveUnits(wIn, -1)

    // --- switch to folder containing wIn ---
    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR   = GetWavesDataFolderDFR(wIn)
    SetDataFolder wDFR

    // --- default output names ---
    String base = CleanupName(NameOfWave(wIn), 0)

    if (ParamIsDefault(outAbsDerivName) || strlen(outAbsDerivName)==0)
        outAbsDerivName = base + "_absD"
    endif
    if (ParamIsDefault(peakXName) || strlen(peakXName)==0)
        peakXName = base + "_absD_peakX"
    endif
    if (ParamIsDefault(peakYName) || strlen(peakYName)==0)
        peakYName = base + "_absD_peakY"
    endif
    if (ParamIsDefault(peakCountVarName) || strlen(peakCountVarName)==0)
        peakCountVarName = base + "_absD_nPeaks"
    endif

    outAbsDerivName  = CleanupName(outAbsDerivName, 0)
    peakXName        = CleanupName(peakXName, 0)
    peakYName        = CleanupName(peakYName, 0)
    peakCountVarName = CleanupName(peakCountVarName, 0)

    // --------------------------
    // Build preprocessed input (FREE) -> wP, with x step dxP
    // --------------------------
    Variable nP, dxP, xOffP
    xOffP = xOff
    dxP   = dx

    Make/FREE/N=(n0) wBase
    wBase[] = wIn[p]

    Make/FREE/N=1 wP   // will be resized

    if (preInterpFactor <= 1)
        // no interpolation; just optional smoothing on original grid
        Redimension/N=(n0) wP
        wP[] = wBase[p]
        if (preSmoothPts > 0)
            Smooth preSmoothPts, wP
        endif
        nP  = n0
        dxP = dx
    else
        // interpolation involved
        if (!preSmoothAfterInterp)
            // smooth -> interp
            if (preSmoothPts > 0)
                Smooth preSmoothPts, wBase
            endif
        endif

        nP  = (n0 - 1)*preInterpFactor + 1
        dxP = dx / preInterpFactor

        Redimension/N=(nP) wP

        // slopes for cubic Hermite (computed on wBase)
        Make/FREE/N=(n0) m
        m[0]    = (wBase[1]    - wBase[0])    / dx
        m[n0-1] = (wBase[n0-1] - wBase[n0-2]) / dx
        Variable i
        for (i=1; i<n0-1; i+=1)
            m[i] = (wBase[i+1] - wBase[i-1]) / (2*dx)
        endfor

        Variable seg, sub, outIdx
        Variable t, t2, t3
        Variable y0, y1, m0, m1, y

        for (seg=0; seg<n0-1; seg+=1)
            y0 = wBase[seg]
            y1 = wBase[seg+1]
            m0 = m[seg]
            m1 = m[seg+1]

            for (sub=0; sub<preInterpFactor; sub+=1)
                t = sub / preInterpFactor
                outIdx = seg*preInterpFactor + sub

                if (!preInterpMethod)
                    // linear
                    y = (1-t)*y0 + t*y1
                else
                    // cubic Hermite
                    t2 = t*t
                    t3 = t2*t
                    y = (2*t3 - 3*t2 + 1)*y0 + (t3 - 2*t2 + t)*(dx*m0) + (-2*t3 + 3*t2)*y1 + (t3 - t2)*(dx*m1)
                endif
                wP[outIdx] = y
            endfor
        endfor
        wP[nP-1] = wBase[n0-1]

        if (preSmoothAfterInterp)
            // interp -> smooth (smooth happens on finer grid)
            if (preSmoothPts > 0)
                Smooth preSmoothPts, wP
            endif
        endif
    endif

    // --------------------------
    // SAVE the preprocessed wave used for derivative as "<wIn>_absD_preInt"
    // --------------------------
    String preName = CleanupName(base + "_absD_preInt", 0)
    Make/O/N=(nP) $preName
    WAVE wPre = $preName
    wPre[] = wP[p]
    SetScale/P x, xOffP, dxP, xUnits, wPre
    SetScale d, 0, 0, dUnits, wPre
    Note wPre, "Preprocessed input for absD; preSmoothPts="+num2istr(preSmoothPts)+ \
               "; preInterpFactor="+num2istr(preInterpFactor)+ \
               "; preInterpMethod="+SelectString(preInterpMethod,"linear","cubicHermite")+ \
               "; preOrder="+SelectString(preSmoothAfterInterp,"smooth->interp","interp->smooth")

    // --------------------------
    // Compute abs(first derivative) on wP grid
    // --------------------------
    Make/O/N=(nP) $outAbsDerivName
    WAVE wAbsD = $outAbsDerivName

    wAbsD[0]     = abs((wP[1]      - wP[0])      / dxP)
    wAbsD[nP-1]  = abs((wP[nP-1]   - wP[nP-2])   / dxP)

    for (i=1; i<nP-1; i+=1)
        wAbsD[i] = abs((wP[i+1] - wP[i-1]) / (2*dxP))
    endfor

    SetScale/P x, xOffP, dxP, xUnits, wAbsD
    SetScale d, 0, 0, dUnits, wAbsD

    if (smoothPts > 0)
        Smooth smoothPts, wAbsD
    endif

    Note wAbsD, "abs(d/dx) of "+NameOfWave(wIn)+ \
                "; preSmoothPts="+num2istr(preSmoothPts)+ \
                "; preInterpFactor="+num2istr(preInterpFactor)+ \
                "; preInterpMethod="+SelectString(preInterpMethod,"linear","cubicHermite")+ \
                "; preOrder="+SelectString(preSmoothAfterInterp,"smooth->interp","interp->smooth")+ \
                "; dxUsed="+num2str(dxP)+ \
                "; postSmoothPts="+num2istr(smoothPts)+ \
                "; thresh="+num2str(thresh)+ \
                "; minDistPts="+num2istr(minDistPts)

    // --------------------------
    // Peak finding on wAbsD
    // --------------------------
    Make/FREE/N=(nP) peakMask
    peakMask = 0

    Variable halfWin = minDistPts
    Variable left, right, j, maxIdx
    Variable maxVal

    for (i=1; i<nP-1; i+=1)
        if (wAbsD[i] >= thresh)
            if ((wAbsD[i] > wAbsD[i-1]) && (wAbsD[i] >= wAbsD[i+1]))
                left  = max(0, i - halfWin)
                right = min(nP-1, i + halfWin)

                maxVal = -1e308
                maxIdx = -1
                for (j=left; j<=right; j+=1)
                    if (wAbsD[j] > maxVal)
                        maxVal = wAbsD[j]
                        maxIdx = j
                    endif
                endfor

                if (maxIdx == i)
                    peakMask[i] = 1
                endif
            endif
        endif
    endfor

    Variable cnt = 0
    for (i=0; i<nP; i+=1)
        if (peakMask[i] != 0)
            cnt += 1
        endif
    endfor

    Make/O/N=(cnt) $peakXName, $peakYName
    WAVE wPX = $peakXName
    WAVE wPY = $peakYName

    Variable k = 0
    for (i=0; i<nP; i+=1)
        if (peakMask[i] != 0)
            wPX[k] = pnt2x(wAbsD, i)
            wPY[k] = wAbsD[i]
            k += 1
        endif
    endfor

    Variable/G $peakCountVarName = cnt

    Note wPX, "Peak x positions of "+NameOfWave(wAbsD)
    Note wPY, "Peak heights of "+NameOfWave(wAbsD)

    SetDataFolder oldDFR
End



//==============================================================================
// AbsDerivFindPeaksLineSTS
//
// ABS(first-derivative) peak finder for a 2D LineSTS wave.
//
// Optional:
//   outBaseName : replaces <wIn> in all output names
//
// Outputs:
//   <base>_absD
//   <base>_absD_preInt
//   <base>_absD_peakX
//   <base>_absD_peakY
//   <base>_absD_nPeaks
//
// Assumed default layout:
//   wIn[energy/bias][line position]
//   derivDim = 0
//==============================================================================

Function SNS_ABSSecondDerivFindPeaks1D(wIn, [outD2Name, peakXName, peakYName, peakCountVarName, thresh, minDistPts, smoothPts, preSmoothPts, preInterpFactor, preInterpMethod, preSmoothAfterInterp])
    Wave    wIn
    String  outD2Name, peakXName, peakYName, peakCountVarName
    Variable thresh, minDistPts, smoothPts, preSmoothPts, preInterpFactor, preInterpMethod, preSmoothAfterInterp

    // --- checks ---
    if (WaveDims(wIn) != 1)
        Abort "SNS_ABSSecondDerivFindPeaks1D: input must be 1D."
    endif
    Variable n0 = numpnts(wIn)
    if (n0 < 3)
        Abort "SNS_ABSSecondDerivFindPeaks1D: input must have at least 3 points."
    endif

    // defaults
    if (ParamIsDefault(thresh))
        thresh = 0
    endif
    if (ParamIsDefault(minDistPts))
        minDistPts = 1
    endif
    if (ParamIsDefault(smoothPts))
        smoothPts = 0
    endif
    if (ParamIsDefault(preSmoothPts))
        preSmoothPts = 0
    endif
    if (ParamIsDefault(preInterpFactor))
        preInterpFactor = 1
    endif
    if (ParamIsDefault(preInterpMethod))
        preInterpMethod = 1
    endif
    if (ParamIsDefault(preSmoothAfterInterp))
        preSmoothAfterInterp = 0
    endif

    minDistPts           = max(1, round(minDistPts))
    smoothPts            = max(0, round(smoothPts))
    preSmoothPts         = max(0, round(preSmoothPts))
    preInterpFactor      = max(1, round(preInterpFactor))
    preInterpMethod      = (preInterpMethod != 0)       // 0=linear, 1=cubic
    preSmoothAfterInterp = (preSmoothAfterInterp != 0)  // 0 or 1

    // x scaling (assumes linear). If unscaled, dx=1.
    Variable xOff = DimOffset(wIn, 0)
    Variable dx   = DimDelta(wIn, 0)
    if (dx == 0)
        dx = 1
    endif
    String xUnits = WaveUnits(wIn, 0)
    String dUnits = WaveUnits(wIn, -1)

    // --- switch to folder containing wIn ---
    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR   = GetWavesDataFolderDFR(wIn)
    SetDataFolder wDFR

    // --- default output names ---
    String base = CleanupName(NameOfWave(wIn), 0)

    if (ParamIsDefault(outD2Name) || strlen(outD2Name)==0)
        outD2Name = base + "_absD2"
    endif
    if (ParamIsDefault(peakXName) || strlen(peakXName)==0)
        peakXName = base + "_absD2_peakX"
    endif
    if (ParamIsDefault(peakYName) || strlen(peakYName)==0)
        peakYName = base + "_absD2_peakY"
    endif
    if (ParamIsDefault(peakCountVarName) || strlen(peakCountVarName)==0)
        peakCountVarName = base + "_absD2_nPeaks"
    endif

    outD2Name        = CleanupName(outD2Name, 0)
    peakXName        = CleanupName(peakXName, 0)
    peakYName        = CleanupName(peakYName, 0)
    peakCountVarName = CleanupName(peakCountVarName, 0)

    // --------------------------
    // Build preprocessed input (FREE) -> wP, with x step dxP
    // --------------------------
    Variable nP, dxP, xOffP
    xOffP = xOff
    dxP   = dx

    Make/FREE/N=(n0) wBase
    wBase[] = wIn[p]

    Make/FREE/N=1 wP   // will be resized

    if (preInterpFactor <= 1)
        // no interpolation; just optional smoothing on original grid
        Redimension/N=(n0) wP
        wP[] = wBase[p]
        if (preSmoothPts > 0)
            Smooth preSmoothPts, wP
        endif
        nP  = n0
        dxP = dx
    else
        // interpolation involved
        if (!preSmoothAfterInterp)
            // smooth -> interp
            if (preSmoothPts > 0)
                Smooth preSmoothPts, wBase
            endif
        endif

        nP  = (n0 - 1)*preInterpFactor + 1
        dxP = dx / preInterpFactor

        Redimension/N=(nP) wP

        // slopes for cubic Hermite (computed on wBase)
        Make/FREE/N=(n0) m
        m[0]    = (wBase[1]    - wBase[0])    / dx
        m[n0-1] = (wBase[n0-1] - wBase[n0-2]) / dx
        Variable i
        for (i=1; i<n0-1; i+=1)
            m[i] = (wBase[i+1] - wBase[i-1]) / (2*dx)
        endfor

        Variable seg, sub, outIdx
        Variable t, t2, t3
        Variable y0, y1, m0, m1, y

        for (seg=0; seg<n0-1; seg+=1)
            y0 = wBase[seg]
            y1 = wBase[seg+1]
            m0 = m[seg]
            m1 = m[seg+1]

            for (sub=0; sub<preInterpFactor; sub+=1)
                t = sub / preInterpFactor
                outIdx = seg*preInterpFactor + sub

                if (!preInterpMethod)
                    y = (1-t)*y0 + t*y1
                else
                    t2 = t*t
                    t3 = t2*t
                    y = (2*t3 - 3*t2 + 1)*y0 + (t3 - 2*t2 + t)*(dx*m0) + (-2*t3 + 3*t2)*y1 + (t3 - t2)*(dx*m1)
                endif
                wP[outIdx] = y
            endfor
        endfor
        wP[nP-1] = wBase[n0-1]

        if (preSmoothAfterInterp)
            // interp -> smooth (smooth happens on finer grid)
            if (preSmoothPts > 0)
                Smooth preSmoothPts, wP
            endif
        endif
    endif

    // --------------------------
    // SAVE the preprocessed wave used for derivative as "<wIn>_absD2_preInt"
    // --------------------------
    String preName = CleanupName(base + "_absD2_preInt", 0)
    Make/O/N=(nP) $preName
    WAVE wPre = $preName
    wPre[] = wP[p]
    SetScale/P x, xOffP, dxP, xUnits, wPre
    SetScale d, 0, 0, dUnits, wPre
    Note wPre, "Preprocessed input for absD2; preSmoothPts="+num2istr(preSmoothPts)+ \
               "; preInterpFactor="+num2istr(preInterpFactor)+ \
               "; preInterpMethod="+SelectString(preInterpMethod,"linear","cubicHermite")+ \
               "; preOrder="+SelectString(preSmoothAfterInterp,"smooth->interp","interp->smooth")

    // --------------------------
    // Compute abs(second derivative) on wP grid
    // --------------------------
    Make/O/N=(nP) $outD2Name
    WAVE wD2 = $outD2Name

    Variable dx2 = dxP*dxP

    wD2[0]    = abs((wP[2]     - 2*wP[1]     + wP[0])     / dx2)
    wD2[nP-1] = abs((wP[nP-1]  - 2*wP[nP-2]  + wP[nP-3])  / dx2)

    for (i=1; i<nP-1; i+=1)
        wD2[i] = abs((wP[i+1] - 2*wP[i] + wP[i-1]) / dx2)
    endfor

    SetScale/P x, xOffP, dxP, xUnits, wD2
    SetScale d, 0, 0, dUnits, wD2

    if (smoothPts > 0)
        Smooth smoothPts, wD2
    endif

    Note wD2, "abs(d2/dx2) of "+NameOfWave(wIn)+ \
              "; preSmoothPts="+num2istr(preSmoothPts)+ \
              "; preInterpFactor="+num2istr(preInterpFactor)+ \
              "; preInterpMethod="+SelectString(preInterpMethod,"linear","cubicHermite")+ \
              "; preOrder="+SelectString(preSmoothAfterInterp,"smooth->interp","interp->smooth")+ \
              "; dxUsed="+num2str(dxP)+ \
              "; postSmoothPts="+num2istr(smoothPts)

    // --------------------------
    // Peak finding on wD2
    // --------------------------
    Make/FREE/N=(nP) peakMask
    peakMask = 0

    Variable halfWin = minDistPts
    Variable left, right, j, maxIdx
    Variable maxVal

    for (i=1; i<nP-1; i+=1)
        if (wD2[i] >= thresh)
            if ((wD2[i] > wD2[i-1]) && (wD2[i] >= wD2[i+1]))
                left  = max(0, i - halfWin)
                right = min(nP-1, i + halfWin)

                maxVal = -1e308
                maxIdx = -1
                for (j=left; j<=right; j+=1)
                    if (wD2[j] > maxVal)
                        maxVal = wD2[j]
                        maxIdx = j
                    endif
                endfor

                if (maxIdx == i)
                    peakMask[i] = 1
                endif
            endif
        endif
    endfor

    Variable cnt = 0
    for (i=0; i<nP; i+=1)
        if (peakMask[i] != 0)
            cnt += 1
        endif
    endfor

    Make/O/N=(cnt) $peakXName, $peakYName
    WAVE wPX = $peakXName
    WAVE wPY = $peakYName

    Variable k = 0
    for (i=0; i<nP; i+=1)
        if (peakMask[i] != 0)
            wPX[k] = pnt2x(wD2, i)
            wPY[k] = wD2[i]
            k += 1
        endif
    endfor

    Variable/G $peakCountVarName = cnt

    SetDataFolder oldDFR
End




//==============================================================================
// ABSSecondDerivFindPeaksLineSTS
//
// ABS(second-derivative) peak finder for a 2D LineSTS wave.
//
// Optional:
//   outBaseName : replaces <wIn> in all output names
//
// Outputs:
//   <base>_absD2
//   <base>_absD2_preInt
//   <base>_absD2_peakX
//   <base>_absD2_peakY
//   <base>_absD2_nPeaks
//
// Assumed default layout:
//   wIn[energy/bias][line position]
//   derivDim = 0
//==============================================================================

Function SNS_FindMaxima1D(wIn, [outName, peakXName, peakYName, peakCountVarName, thresh, minDistPts, smoothPts, interpFactor, interpMethod, smoothAfterInterp, xMin, xMax])
    Wave    wIn
    String  outName, peakXName, peakYName, peakCountVarName
    Variable thresh, minDistPts, smoothPts, interpFactor, interpMethod, smoothAfterInterp
    Variable xMin, xMax

    if (WaveDims(wIn) != 1)
        Abort "SNS_FindMaxima1D_Make: input must be 1D."
    endif
    Variable n0 = numpnts(wIn)
    if (n0 < 3)
        Abort "SNS_FindMaxima1D_Make: input must have at least 3 points."
    endif

    // remember whether user provided xMin/xMax (ParamIsDefault may be affected if we later assign)
    Variable haveXMin = !ParamIsDefault(xMin)
    Variable haveXMax = !ParamIsDefault(xMax)

    if (ParamIsDefault(thresh))
        thresh = -1e308
    endif
    if (ParamIsDefault(minDistPts))
        minDistPts = 1
    endif
    if (ParamIsDefault(smoothPts))
        smoothPts = 0
    endif
    if (ParamIsDefault(interpFactor))
        interpFactor = 1
    endif
    if (ParamIsDefault(interpMethod))
        interpMethod = 1
    endif
    if (ParamIsDefault(smoothAfterInterp))
        smoothAfterInterp = 0
    endif

    minDistPts        = max(1, round(minDistPts))
    smoothPts         = max(0, round(smoothPts))
    interpFactor      = max(1, round(interpFactor))
    interpMethod      = (interpMethod != 0)         // 0=linear, 1=cubic
    smoothAfterInterp = (smoothAfterInterp != 0)    // 0 or 1

    // x scaling (assumes linear). If unscaled, dx=1.
    Variable xOff = DimOffset(wIn, 0)
    Variable dx   = DimDelta(wIn, 0)
    if (dx == 0)
        dx = 1
    endif
    String xUnits = WaveUnits(wIn, 0)
    String dUnits = WaveUnits(wIn, -1)

    // Switch to folder containing wIn
    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR   = GetWavesDataFolderDFR(wIn)
    SetDataFolder wDFR

    String base = CleanupName(NameOfWave(wIn), 0)

    // Always default to "<input>_proc" (no info in the name)
    if (ParamIsDefault(outName) || strlen(outName) == 0)
        outName = base + "_proc"
    endif
    outName = CleanupName(outName, 0)

    if (ParamIsDefault(peakXName) || strlen(peakXName) == 0)
        peakXName = outName + "_maxX"
    endif
    if (ParamIsDefault(peakYName) || strlen(peakYName) == 0)
        peakYName = outName + "_maxY"
    endif
    if (ParamIsDefault(peakCountVarName) || strlen(peakCountVarName) == 0)
        peakCountVarName = outName + "_nMax"
    endif

    peakXName        = CleanupName(peakXName, 0)
    peakYName        = CleanupName(peakYName, 0)
    peakCountVarName = CleanupName(peakCountVarName, 0)

    // --------------------------
    // Build processed wave outName
    // --------------------------
    // We'll use a FREE base array for interpolation source:
    Make/FREE/N=(n0) wBase
    wBase[] = wIn[p]

    // Case 1: No interpolation -> outName is (optionally) smoothed copy on original grid
    if (interpFactor <= 1)

        Duplicate/O wIn, $outName
        WAVE wOut = $outName

        if (smoothPts > 0)
            Smooth smoothPts, wOut
        endif

    else
        // Case 2: Interpolation involved
        if (!smoothAfterInterp)
            // smooth first on original grid (affects interpolation result)
            if (smoothPts > 0)
                Smooth smoothPts, wBase
            endif
        endif

        Variable newN = (n0 - 1)*interpFactor + 1
        Make/O/N=(newN) $outName
        WAVE wOut = $outName

        SetScale/P x, xOff, dx/interpFactor, xUnits, wOut
        SetScale d, 0, 0, dUnits, wOut

        // slopes for cubic Hermite
        Make/FREE/N=(n0) m
        m[0]    = (wBase[1]     - wBase[0])     / dx
        m[n0-1] = (wBase[n0-1]  - wBase[n0-2])  / dx
        Variable i
        for (i=1; i<n0-1; i+=1)
            m[i] = (wBase[i+1] - wBase[i-1]) / (2*dx)
        endfor

        Variable seg, sub, outIdx
        Variable t, t2, t3
        Variable y0, y1, m0, m1, y

        for (seg=0; seg<n0-1; seg+=1)
            y0 = wBase[seg]
            y1 = wBase[seg+1]
            m0 = m[seg]
            m1 = m[seg+1]

            for (sub=0; sub<interpFactor; sub+=1)
                t = sub / interpFactor
                outIdx = seg*interpFactor + sub

                if (!interpMethod)
                    // linear
                    y = (1-t)*y0 + t*y1
                else
                    // cubic Hermite
                    t2 = t*t
                    t3 = t2*t
                    y = (2*t3 - 3*t2 + 1)*y0 + (t3 - 2*t2 + t)*(dx*m0) + (-2*t3 + 3*t2)*y1 + (t3 - t2)*(dx*m1)
                endif

                wOut[outIdx] = y
            endfor
        endfor
        wOut[newN-1] = wBase[n0-1]

        if (smoothAfterInterp)
            // smooth after interpolation (smoothing happens on the finer grid)
            if (smoothPts > 0)
                Smooth smoothPts, wOut
            endif
        endif
    endif

    // Add a note to the processed wave
    WAVE wOut2 = $outName
    Note wOut2, "Processed from "+NameOfWave(wIn)+"; smoothPts="+num2istr(smoothPts)+"; interpFactor="+num2istr(interpFactor)+"; interpMethod="+SelectString(interpMethod,"linear","cubicHermite")+"; smoothAfterInterp="+num2istr(smoothAfterInterp)

    // --------------------------
    // Peak finding on processed wave (outName), restricted to x-range if requested
    // --------------------------
    Variable n = numpnts(wOut2)
    Make/FREE/N=(n) peakMask
    peakMask = 0

    // Determine x-range defaults based on processed wave scaling
    Variable x0 = pnt2x(wOut2, 0)
    Variable x1 = pnt2x(wOut2, n-1)
    Variable xWaveMin = min(x0, x1)
    Variable xWaveMax = max(x0, x1)

    if (!haveXMin)
        xMin = xWaveMin
    endif
    if (!haveXMax)
        xMax = xWaveMax
    endif

    // Ensure order and clamp to wave range
    Variable tmp
    if (xMin > xMax)
        tmp = xMin; xMin = xMax; xMax = tmp
    endif
    xMin = max(xMin, xWaveMin)
    xMax = min(xMax, xWaveMax)

    // Convert x-range to point-range (inclusive). Use ceil for start, floor for end.
    Variable pA = x2pnt(wOut2, xMin)
    Variable pB = x2pnt(wOut2, xMax)
    Variable pLo = min(pA, pB)
    Variable pHi = max(pA, pB)

    // ceil(pLo) implemented via floor(pLo + 1 - eps) to avoid needing ceil()
    Variable eps = 1e-9
    Variable iStart = floor(pLo + 1 - eps)
    Variable iEnd   = floor(pHi)

    // Keep inside safe bounds for neighbor checks (i-1, i+1)
    iStart = max(1, iStart)
    iEnd   = min(n-2, iEnd)

    Variable halfWin = minDistPts
    Variable left, right, j, maxIdx
    Variable maxVal

    for (i=iStart; i<=iEnd; i+=1)
        if (wOut2[i] >= thresh)
            if ((wOut2[i] > wOut2[i-1]) && (wOut2[i] >= wOut2[i+1]))
                // Restrict the competition window to the requested i-range
                left  = max(iStart, i - halfWin)
                right = min(iEnd,   i + halfWin)

                maxVal = -1e308
                maxIdx = -1
                for (j=left; j<=right; j+=1)
                    if (wOut2[j] > maxVal)
                        maxVal = wOut2[j]
                        maxIdx = j
                    endif
                endfor

                if (maxIdx == i)
                    peakMask[i] = 1
                endif
            endif
        endif
    endfor

    Variable cnt = 0
    for (i=0; i<n; i+=1)
        if (peakMask[i] != 0)
            cnt += 1
        endif
    endfor

    Make/O/N=(cnt) $peakXName, $peakYName
    WAVE wPX = $peakXName
    WAVE wPY = $peakYName

    Variable k = 0
    for (i=0; i<n; i+=1)
        if (peakMask[i] != 0)
            wPX[k] = pnt2x(wOut2, i)
            wPY[k] = wOut2[i]
            k += 1
        endif
    endfor

    Variable/G $peakCountVarName = cnt

    SetDataFolder oldDFR
End



//==============================================================================
// FindMaxLineSTS
//
// Direct local-maxima finder for a 2D line STS.
// Uses SNS_FindMaxima1D() on each line spectrum.
//
// Output layout matches the absD-style peak finder:
//
//   <base>_max
//   <base>_max_preInt
//   <base>_max_peakX
//   <base>_max_peakY
//   <base>_max_nPeaks
//
// Dimensions:
//   <base>_max_peakX  : [peak index][line index], data = peak x/energy
//   <base>_max_peakY  : [peak index][line index], data = peak value
//   <base>_max_nPeaks : [line index], data = number of peaks
//
// energyDim:
//   0 : energy/x is dim 0, line/r is dim 1
//   1 : energy/x is dim 1, line/r is dim 0
//
// Notes:
//   preInt stores the spectrum after preSmoothPts, before interpolation.
//   max stores the processed/interpolated spectrum used for peak finding.
//==============================================================================

Function SNS_SaveNthGapClosingFromMaxX(wMaxX, n, [avgN, xRes, outName])
    Wave wMaxX
    Variable n
    Variable avgN, xRes
    String outName

    if (WaveDims(wMaxX) != 1)
        Abort "SNS_SaveNthGapClosingFromMaxX: wMaxX must be 1D."
    endif

    n = round(n)
    if (n == 0)
        Abort "SNS_SaveNthGapClosingFromMaxX: n must be nonzero (n>0 for + side, n<0 for - side)."
    endif

    if (ParamIsDefault(avgN))
        avgN = 0
    endif
    avgN = (avgN != 0)   // force 0/1

    if (ParamIsDefault(xRes))
        xRes = 10
    endif
    xRes = abs(xRes)

    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR   = GetWavesDataFolderDFR(wMaxX)
    SetDataFolder wDFR

    Variable Nall = numpnts(wMaxX)
    Variable i, x
    Variable k = abs(n)
    Variable outVal, errVal

    if (!avgN)
        // -----------------------
        // OLD behavior (one sign)
        // -----------------------
        Make/FREE/N=(Nall) tmp
        Variable m = 0
        Variable wantNeg = (n < 0)

        for (i=0; i<Nall; i+=1)
            x = wMaxX[i]
            if (numtype(x) != 0)       // skip NaN/Inf
                continue
            endif

            if (!wantNeg)
                if (x > 0)
                    tmp[m] = x
                    m += 1
                endif
            else
                if (x < 0)
                    tmp[m] = x
                    m += 1
                endif
            endif
        endfor

        if (m < k)
            SetDataFolder oldDFR
            Abort "SNS_SaveNthGapClosingFromMaxX: not enough values of the requested sign."
        endif

        Redimension/N=(m) tmp
        Sort tmp, tmp   // ascending

        Variable sel
        if (!wantNeg)
            sel = tmp[k-1]      // k-th smallest positive
        else
            sel = tmp[m-k]      // k-th largest negative (closest to 0 first)
        endif

        outVal = abs(sel)
        errVal = xRes

    else
        // -----------------------------------------
        // avgN behavior: average +k and -k closings
        // -----------------------------------------
        Make/FREE/N=(Nall) tmpPos, tmpNeg
        Variable mPos = 0
        Variable mNeg = 0

        for (i=0; i<Nall; i+=1)
            x = wMaxX[i]
            if (numtype(x) != 0)       // skip NaN/Inf
                continue
            endif

            if (x > 0)
                tmpPos[mPos] = x
                mPos += 1
            elseif (x < 0)
                tmpNeg[mNeg] = x
                mNeg += 1
            endif
        endfor

        if (mPos < k || mNeg < k)
            SetDataFolder oldDFR
            Abort "SNS_SaveNthGapClosingFromMaxX: avgN=1 requires at least k positive AND k negative values."
        endif

        Redimension/N=(mPos) tmpPos
        Redimension/N=(mNeg) tmpNeg

        Sort tmpPos, tmpPos   // ascending positives
        Sort tmpNeg, tmpNeg   // ascending negatives

        Variable posSel = tmpPos[k-1]
        Variable negSel = tmpNeg[mNeg-k]

        outVal = 0.5*(abs(posSel) + abs(negSel))

        Variable errAsym = 0.5*abs(abs(posSel) - abs(negSel))
        errVal = max(errAsym, xRes)
    endif

    // ----- output variable names -----
    String prefix = ""
    if (!ParamIsDefault(outName))
        if (strlen(outName) > 0)
            prefix = CleanupName(outName, 0) + "_"
        endif
    endif

    String varName    = CleanupName(prefix + "No_" + num2istr(k) + "_GapClosing", 0)
    String errVarName = CleanupName(prefix + "No_" + num2istr(k) + "_GapClosingErr", 0)

    Variable/G $varName    = outVal
    Variable/G $errVarName = errVal

    SetDataFolder oldDFR
End









// ============================================================
// Helper (1D): Find symmetry center index by minimizing
//   score(c) = weighted mean_{Δ} | w[c+Δ] - w[c-Δ] |
// weights favor small Δ.
//
// Output in resW (length 3):
//   resW[0] = best center index (>=0), or -1 if failed
//   resW[1] = best score (smaller = more symmetric; units of w)
//   resW[2] = number of Δ-pairs used at best center
// ============================================================

Function/WAVE SNS_MakeE0BoundMarkersFromPeakWaves(wLowPeakX, wHighPeakX, [outName, prefix, lowPeakOrder, highPeakOrder])
    Wave wLowPeakX, wHighPeakX
    String outName, prefix
    Variable lowPeakOrder, highPeakOrder

    Variable i, x
    Variable nNegLow=0, nPosLow=0, nNegHigh=0, nPosHigh=0

    if (WaveDims(wLowPeakX) != 1 || WaveDims(wHighPeakX) != 1)
        Abort "SNS_MakeE0BoundMarkersFromPeakWaves: both inputs must be 1D waves."
    endif

    if (ParamIsDefault(outName) || strlen(outName) == 0)
        outName = "E0BoundX"
    endif
    outName = CleanupName(outName, 0)

    if (ParamIsDefault(prefix) || strlen(prefix) == 0)
        prefix = "E0"
    endif
    prefix = CleanupName(prefix, 0)

    if (ParamIsDefault(lowPeakOrder))
        lowPeakOrder = 1
    endif
    if (ParamIsDefault(highPeakOrder))
        highPeakOrder = 1
    endif

    lowPeakOrder  = round(lowPeakOrder)
    highPeakOrder = round(highPeakOrder)

    if (lowPeakOrder < 1)
        Abort "SNS_MakeE0BoundMarkersFromPeakWaves: lowPeakOrder must be >= 1."
    endif
    if (highPeakOrder < 1)
        Abort "SNS_MakeE0BoundMarkersFromPeakWaves: highPeakOrder must be >= 1."
    endif

    // Work in folder of the first wave
    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR   = GetWavesDataFolderDFR(wLowPeakX)
    SetDataFolder wDFR

    // Count entries needed for temporary waves
    for (i=0; i<numpnts(wLowPeakX); i+=1)
        x = wLowPeakX[i]
        if (numtype(x) != 0)
            continue
        endif
        if (x < 0)
            nNegLow += 1
        elseif (x > 0)
            nPosLow += 1
        endif
    endfor

    for (i=0; i<numpnts(wHighPeakX); i+=1)
        x = wHighPeakX[i]
        if (numtype(x) != 0)
            continue
        endif
        if (x < 0)
            nNegHigh += 1
        elseif (x > 0)
            nPosHigh += 1
        endif
    endfor

    if (nNegLow < lowPeakOrder || nPosLow < lowPeakOrder)
        SetDataFolder oldDFR
        Abort "SNS_MakeE0BoundMarkersFromPeakWaves: not enough lower-bound peaks on both sides."
    endif

    if (nNegHigh < highPeakOrder || nPosHigh < highPeakOrder)
        SetDataFolder oldDFR
        Abort "SNS_MakeE0BoundMarkersFromPeakWaves: not enough upper-bound peaks on both sides."
    endif

    // Create and fill temporary waves
    Make/FREE/N=(nNegLow)  low_neg
    Make/FREE/N=(nPosLow)  low_pos
    Make/FREE/N=(nNegHigh) high_neg
    Make/FREE/N=(nPosHigh) high_pos

    Variable iNegLow=0, iPosLow=0, iNegHigh=0, iPosHigh=0

    for (i=0; i<numpnts(wLowPeakX); i+=1)
        x = wLowPeakX[i]
        if (numtype(x) != 0)
            continue
        endif
        if (x < 0)
            low_neg[iNegLow] = x
            iNegLow += 1
        elseif (x > 0)
            low_pos[iPosLow] = x
            iPosLow += 1
        endif
    endfor

    for (i=0; i<numpnts(wHighPeakX); i+=1)
        x = wHighPeakX[i]
        if (numtype(x) != 0)
            continue
        endif
        if (x < 0)
            high_neg[iNegHigh] = x
            iNegHigh += 1
        elseif (x > 0)
            high_pos[iPosHigh] = x
            iPosHigh += 1
        endif
    endfor

    // Sort ascending
    Sort low_neg,   low_neg
    Sort low_pos,   low_pos
    Sort high_neg,  high_neg
    Sort high_pos,  high_pos

    // Pick peaks counted from closest to zero outward
    Variable E0LowNeg  = low_neg[nNegLow - lowPeakOrder]
    Variable E0LowPos  = low_pos[lowPeakOrder - 1]
    Variable E0HighNeg = high_neg[nNegHigh - highPeakOrder]
    Variable E0HighPos = high_pos[highPeakOrder - 1]

    // Save variables
    Variable/G $(prefix + "LowNeg")  = E0LowNeg
    Variable/G $(prefix + "LowPos")  = E0LowPos
    Variable/G $(prefix + "HighNeg") = E0HighNeg
    Variable/G $(prefix + "HighPos") = E0HighPos

    // Output wave
    Make/O/N=4 $outName
    Wave wOut = $outName

    wOut[0] = E0HighNeg
    wOut[1] = E0LowNeg
    wOut[2] = E0LowPos
    wOut[3] = E0HighPos

    Note/K wOut, "E0 bounds from peak waves; lowPeakOrder="+num2istr(lowPeakOrder)+"; highPeakOrder="+num2istr(highPeakOrder)

    SetDataFolder oldDFR
    return wOut
End





//==============================================================================
// SNS_MakeNearestPeakPairAroundZero
//
// From a 1D wave of x-positions, find the two values closest to 0:
//   - largest negative value
//   - smallest positive value
//
// Creates a 2-point output wave in the SAME data folder as wX:
//   [0] = nearest negative peak
//   [1] = nearest positive peak
//
// Optional:
//   outName : output wave name (default: "<wX>_nearestZeroPair")
//
// Returns:
//   reference to the output wave
//==============================================================================
Function/WAVE SNS_MakeNearestPeakPairAroundZero(wX, [outName])
    Wave wX
    String outName

    Variable i, x
    Variable bestNeg = -1e308
    Variable bestPos =  1e308
    Variable haveNeg = 0
    Variable havePos = 0

    if (WaveDims(wX) != 1)
        Abort "SNS_MakeNearestPeakPairAroundZero: input must be 1D."
    endif

    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR   = GetWavesDataFolderDFR(wX)
    SetDataFolder wDFR

    if (ParamIsDefault(outName) || strlen(outName) == 0)
        outName = CleanupName(NameOfWave(wX) + "_nearestZeroPair", 0)
    endif

    for (i=0; i<numpnts(wX); i+=1)
        x = wX[i]
        if (numtype(x) != 0)
            continue
        endif

        if (x < 0)
            if (x > bestNeg)
                bestNeg = x
                haveNeg = 1
            endif
        elseif (x > 0)
            if (x < bestPos)
                bestPos = x
                havePos = 1
            endif
        endif
    endfor

    if (!haveNeg || !havePos)
        SetDataFolder oldDFR
        Abort "SNS_MakeNearestPeakPairAroundZero: need at least one negative and one positive value."
    endif

    Make/O/N=2 $outName
    Wave wOut = $outName
    wOut[0] = bestNeg
    wOut[1] = bestPos

    Note/K wOut, "Nearest negative and positive peak around zero from " + NameOfWave(wX)

    SetDataFolder oldDFR
    return wOut
End




//==============================================================================
// MakeSecondDerivative1D
//
// Returns the signed second derivative d^2y/dx^2 of a 1D wave.
// Optional smoothing is applied before differentiation.
//
// Inputs:
//   wIn
//
// Optional:
//   smoothPts : smoothing points before differentiation (default 0)
//   outName   : output wave name (default: <input>_dd)
//==============================================================================


Function/WAVE SNS_MakeE0UpperBoundPair(wLDOSMaxX, wAbsD2PeakX, [outName, forceD2Fallback])
    Wave wLDOSMaxX, wAbsD2PeakX
    String outName
    Variable forceD2Fallback

    Variable i, x
    Variable bestNeg, bestPos
    Variable haveNeg, havePos

    Variable ldosNeg, ldosPos
    Variable haveLDOSNeg, haveLDOSPos

    String srcName

    if (WaveDims(wLDOSMaxX) != 1 || WaveDims(wAbsD2PeakX) != 1)
        Abort "SNS_MakeE0UpperBoundPair: both inputs must be 1D waves."
    endif

    if (ParamIsDefault(forceD2Fallback))
        forceD2Fallback = 0
    endif
    forceD2Fallback = (forceD2Fallback != 0)

    if (ParamIsDefault(outName) || strlen(outName) == 0)
        outName = "E0HighX"
    endif
    outName = CleanupName(outName, 0)

    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR   = GetWavesDataFolderDFR(wLDOSMaxX)
    SetDataFolder wDFR

    // ------------------------------------------------------------
    // Find nearest LDOS maxima around zero.
    // These are either used directly, or as reference for forced D2.
    // ------------------------------------------------------------
    ldosNeg = -1e308
    ldosPos =  1e308
    haveLDOSNeg = 0
    haveLDOSPos = 0

    for (i=0; i<numpnts(wLDOSMaxX); i+=1)
        x = wLDOSMaxX[i]
        if (numtype(x) != 0)
            continue
        endif

        if (x < 0)
            if (x > ldosNeg)
                ldosNeg = x
                haveLDOSNeg = 1
            endif
        elseif (x > 0)
            if (x < ldosPos)
                ldosPos = x
                haveLDOSPos = 1
            endif
        endif
    endfor

    // ------------------------------------------------------------
    // Default behavior:
    // use LDOS maxima directly if both signs are present.
    // ------------------------------------------------------------
    if (!forceD2Fallback && haveLDOSNeg && haveLDOSPos)
        bestNeg = ldosNeg
        bestPos = ldosPos
        haveNeg = 1
        havePos = 1
        srcName = "LDOSMax"

    else
        // --------------------------------------------------------
        // Forced D2 fallback:
        // use D2 peaks neighboring the LDOS maxima from the inside.
        // --------------------------------------------------------
        if (forceD2Fallback && haveLDOSNeg && haveLDOSPos)
            bestNeg =  1e308     // choose smallest x with x > ldosNeg
            bestPos = -1e308     // choose largest  x with x < ldosPos
            haveNeg = 0
            havePos = 0
            srcName = "AbsD2InsideLDOS"

            for (i=0; i<numpnts(wAbsD2PeakX); i+=1)
                x = wAbsD2PeakX[i]
                if (numtype(x) != 0)
                    continue
                endif

                // Negative side:
                // next larger D2 peak than the negative LDOS maximum.
                if (x < 0 && x > ldosNeg)
                    if (x < bestNeg)
                        bestNeg = x
                        haveNeg = 1
                    endif
                endif

                // Positive side:
                // next smaller D2 peak than the positive LDOS maximum.
                if (x > 0 && x < ldosPos)
                    if (x > bestPos)
                        bestPos = x
                        havePos = 1
                    endif
                endif
            endfor

            if (!haveNeg || !havePos)
                SetDataFolder oldDFR
                Abort "SNS_MakeE0UpperBoundPair: forceD2Fallback=1, but could not find D2 peaks inside both LDOS maxima."
            endif

        else
            // ----------------------------------------------------
            // Automatic fallback:
            // LDOS pair is missing, so use nearest D2 peaks around zero.
            // This keeps the previous fallback behavior.
            // ----------------------------------------------------
            bestNeg = -1e308
            bestPos =  1e308
            haveNeg = 0
            havePos = 0
            srcName = "AbsD2"

            for (i=0; i<numpnts(wAbsD2PeakX); i+=1)
                x = wAbsD2PeakX[i]
                if (numtype(x) != 0)
                    continue
                endif

                if (x < 0)
                    if (x > bestNeg)
                        bestNeg = x
                        haveNeg = 1
                    endif
                elseif (x > 0)
                    if (x < bestPos)
                        bestPos = x
                        havePos = 1
                    endif
                endif
            endfor
        endif
    endif

    if (!haveNeg || !havePos)
        SetDataFolder oldDFR
        Abort "SNS_MakeE0UpperBoundPair: could not find both negative and positive upper-bound peaks."
    endif

    Make/O/N=2 $outName
    Wave wOut = $outName
    wOut[0] = bestNeg
    wOut[1] = bestPos

    Note/K wOut, "E0 upper-bound pair from "+srcName+"; forceD2Fallback="+num2istr(forceD2Fallback)
    if (stringmatch(srcName, "AbsD2InsideLDOS"))
        Note wOut, "; LDOS reference pair = ["+num2str(ldosNeg)+", "+num2str(ldosPos)+"]"
    endif

    SetDataFolder oldDFR
    return wOut
End


// Return the average absolute position of the two peaks with the largest
// finite amplitudes. NaN is returned when fewer than two peaks are available.
Static Function/D SNS_GE_TwoStrongestPeakAverage(wPeakX, wPeakY)
    Wave wPeakX, wPeakY

    Variable best1 = -1e308
    Variable best2 = -1e308
    Variable idx1 = -1
    Variable idx2 = -1
    Variable i

    for (i = 0; i < min(numpnts(wPeakX), numpnts(wPeakY)); i += 1)
        if ((numtype(wPeakX[i]) == 0) && (numtype(wPeakY[i]) == 0))
            if (wPeakY[i] > best1)
                best2 = best1
                idx2 = idx1
                best1 = wPeakY[i]
                idx1 = i
            elseif (wPeakY[i] > best2)
                best2 = wPeakY[i]
                idx2 = i
            endif
        endif
    endfor

    if ((idx1 < 0) || (idx2 < 0))
        return NaN
    endif
    return 0.5 * (abs(wPeakX[idx1]) + abs(wPeakX[idx2]))
End


//==============================================================================
// SNS_ExtractLineSTSTwoPeakAverages
//
// Extract two line-dependent energy estimates from a 2D
// line-spectroscopy wave whose energy axis is dimension 0. For every spatial
// column, the function selects the two strongest peaks from both the absolute
// first derivative and the direct local-maxima analysis and stores the average
// absolute peak position. Processing delegates to the validated SNS 1D peak
// helpers, with explicit controls matching the historical line-STS workflow.
//==============================================================================
Function SNS_ExtractLineSTSTwoPeakAverages(wIn, [outAbsDName, outMaxName, outPositionName, absDSmoothPts, absDPreSmoothPts, absDInterpFactor, maxSmoothPts, maxPreSmoothPts, maxInterpFactor, maxThresh, xMin, xMax])
    Wave wIn
    String outAbsDName, outMaxName, outPositionName
    Variable absDSmoothPts, absDPreSmoothPts, absDInterpFactor
    Variable maxSmoothPts, maxPreSmoothPts, maxInterpFactor, maxThresh
    Variable xMin, xMax

    if (WaveDims(wIn) != 2)
        Abort "SNS_ExtractLineSTSTwoPeakAverages: input wave must be 2D."
    endif
    if (DimSize(wIn, 0) < 3 || DimSize(wIn, 1) < 1)
        Abort "SNS_ExtractLineSTSTwoPeakAverages: invalid input dimensions."
    endif

    if (ParamIsDefault(outAbsDName) || strlen(outAbsDName) == 0)
        outAbsDName = "E0_AbsD"
    endif
    if (ParamIsDefault(outMaxName) || strlen(outMaxName) == 0)
        outMaxName = "E0_Maxima"
    endif
    if (ParamIsDefault(outPositionName) || strlen(outPositionName) == 0)
        outPositionName = "E0_Position"
    endif
    if (ParamIsDefault(absDSmoothPts))
        absDSmoothPts = 1
    endif
    if (ParamIsDefault(absDPreSmoothPts))
        absDPreSmoothPts = 1
    endif
    if (ParamIsDefault(absDInterpFactor))
        absDInterpFactor = 20
    endif
    if (ParamIsDefault(maxSmoothPts))
        maxSmoothPts = 1
    endif
    if (ParamIsDefault(maxPreSmoothPts))
        maxPreSmoothPts = 0
    endif
    if (ParamIsDefault(maxInterpFactor))
        maxInterpFactor = 5
    endif
    if (ParamIsDefault(maxThresh))
        maxThresh = 0.8
    endif
    if (ParamIsDefault(xMin))
        xMin = -1.1
    endif
    if (ParamIsDefault(xMax))
        xMax = 1.1
    endif

    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR = GetWavesDataFolderDFR(wIn)
    SetDataFolder wDFR

    Variable nEnergy = DimSize(wIn, 0)
    Variable nLine = DimSize(wIn, 1)
    Make/O/D/N=(nLine) $outAbsDName, $outMaxName, $outPositionName
    Wave outAbsD = $outAbsDName
    Wave outMax = $outMaxName
    Wave outPosition = $outPositionName
    outAbsD = NaN
    outMax = NaN
    outPosition = DimOffset(wIn, 1) + p * DimDelta(wIn, 1)
    SetScale/P x, DimOffset(wIn, 1), DimDelta(wIn, 1), WaveUnits(wIn, 1), outAbsD, outMax, outPosition
    SetScale d, 0, 0, WaveUnits(wIn, 0), outAbsD, outMax
    SetScale d, 0, 0, WaveUnits(wIn, 1), outPosition

    String tmpSpecName = "SNS_GE_tmpSpec"
    String tmpMaxInputName = "SNS_GE_tmpMaxInput"
    String tmpAbsDName = "SNS_GE_tmpAbsD"
    String tmpAbsDXName = "SNS_GE_tmpAbsDX"
    String tmpAbsDYName = "SNS_GE_tmpAbsDY"
    String tmpAbsDNName = "SNS_GE_tmpAbsDN"
    String tmpProcName = "SNS_GE_tmpProc"
    String tmpMaxXName = "SNS_GE_tmpMaxX"
    String tmpMaxYName = "SNS_GE_tmpMaxY"
    String tmpMaxNName = "SNS_GE_tmpMaxN"

    Make/O/D/N=(nEnergy) $tmpSpecName, $tmpMaxInputName
    Wave tmpSpec = $tmpSpecName
    Wave tmpMaxInput = $tmpMaxInputName
    SetScale/P x, DimOffset(wIn, 0), DimDelta(wIn, 0), WaveUnits(wIn, 0), tmpSpec, tmpMaxInput
    SetScale d, 0, 0, WaveUnits(wIn, -1), tmpSpec, tmpMaxInput

    Variable line
    for (line = 0; line < nLine; line += 1)
        tmpSpec[] = wIn[p][line]
        SNS_AbsDerivFindPeaks1D(tmpSpec, outAbsDerivName=tmpAbsDName, peakXName=tmpAbsDXName, peakYName=tmpAbsDYName, peakCountVarName=tmpAbsDNName, smoothPts=absDSmoothPts, preSmoothPts=absDPreSmoothPts, preInterpFactor=absDInterpFactor)
        Wave tmpAbsDX = $tmpAbsDXName
        Wave tmpAbsDY = $tmpAbsDYName
        outAbsD[line] = SNS_GE_TwoStrongestPeakAverage(tmpAbsDX, tmpAbsDY)

        tmpMaxInput[] = tmpSpec[p]
        if (maxPreSmoothPts > 0)
            Smooth maxPreSmoothPts, tmpMaxInput
        endif
        SNS_FindMaxima1D(tmpMaxInput, outName=tmpProcName, peakXName=tmpMaxXName, peakYName=tmpMaxYName, peakCountVarName=tmpMaxNName, thresh=maxThresh, minDistPts=1, smoothPts=maxSmoothPts, interpFactor=maxInterpFactor, interpMethod=1, smoothAfterInterp=0, xMin=xMin, xMax=xMax)
        Wave tmpMaxX = $tmpMaxXName
        Wave tmpMaxY = $tmpMaxYName
        outMax[line] = SNS_GE_TwoStrongestPeakAverage(tmpMaxX, tmpMaxY)
    endfor

    Note/K outAbsD
    Note outAbsD, "Two-strongest-peak absolute-derivative estimate;source=" + NameOfWave(wIn)
    Note/K outMax
    Note outMax, "Two-strongest-peak direct-maxima estimate;source=" + NameOfWave(wIn)
    Note/K outPosition
    Note outPosition, "Explicit line-position coordinate;source=" + NameOfWave(wIn)

    KillWaves/Z $tmpSpecName, $tmpMaxInputName, $tmpAbsDName, $tmpAbsDXName, $tmpAbsDYName, $tmpProcName, $tmpMaxXName, $tmpMaxYName
    KillVariables/Z $tmpAbsDNName, $tmpMaxNName
    SetDataFolder oldDFR
    return nLine
End
