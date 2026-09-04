#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma DefaultTab={3,20,4}
#include <MatrixToXYZ>

//==============================================================================
// SNS_Helpers.ipf
//
// Generic SNS project helpers for Igor data-folder, path, string, normalization,
// acquisition-parameter export, and notebook workflow operations.
//
// These routines are code-level infrastructure. They do not implement SNS
// physics, ray tracing, LDOS, or model logic.
//==============================================================================

// Canonical acquisition settings read from SIDAM-created settings folders.
// Values remain in SI units here; panel-export helpers perform display-unit
// conversion only when populating explicitly requested CSV columns.
Structure SNS_ExperimentalSettings
    Variable bias_V
    Variable current_A
    Variable vmod_V
    Variable frequency_Hz
    Variable x_m
    Variable y_m
    String lockInStatus
    String biasSpectroscopyLockInRun
    Variable lockInActive
EndStructure


//==============================================================================
// SNS_DetermineActualCoordinates
//
// Code Purpose:
//   Rotate/flip a 2D STM image into a common absolute scan frame and shift its
//   x/y scaling so the image center is at the supplied absolute coordinates.
//
// Physics Role:
//   None. Coordinate-preparation helper for source STM images.
//
// Required Inputs:
//   wIn       : 2D image wave.
//   slowScan  : "up", "u", "down", or "d".
//   angleDeg  : scan angle [deg].
//   xPos_m    : absolute x center [m].
//   yPos_m    : absolute y center [m].
//
// Generated Outputs:
//   Creates/overwrites <inputWaveName>_xy in the current data folder and
//   returns that wave.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools DetermineActualCoordinates.
//==============================================================================
Function/WAVE SNS_DetermineActualCoordinates(wIn, slowScan, angleDeg, xPos_m, yPos_m)
    Wave wIn
    String slowScan
    Variable angleDeg, xPos_m, yPos_m

    String outName = NameOfWave(wIn) + "_xy"

    Variable isUp     = (cmpstr(slowScan, "up") == 0) || (cmpstr(slowScan, "u") == 0)
    Variable isDown   = (cmpstr(slowScan, "down") == 0) || (cmpstr(slowScan, "d") == 0)
    Variable is90     = (angleDeg == 90)
    Variable isSquare = (DimSize(wIn, 0) == DimSize(wIn, 1))

    if (!isUp && !isDown)
        Abort "SNS_DetermineActualCoordinates: slowScan must be 'up'/'u' or 'down'/'d'."
    endif

    Duplicate/O/FREE wIn, wTmp

    if (isUp)
        if (is90)
            wTmp = wIn[DimSize(wIn,1) - q - 1][p]
        else
            ImageRotate/Q/A=(-angleDeg)/O wTmp
        endif
    else
        if (is90 && isSquare)
            ImageRotate/Q/V/O wTmp
            wTmp = wIn[q][p]
            ImageRotate/Q/V/O wTmp
        else
            ImageRotate/Q/V/O wTmp
            ImageRotate/Q/A=(angleDeg)/O wTmp
            ImageRotate/Q/V/O wTmp
        endif
    endif

    Variable xCenter = DimOffset(wTmp, 0) + 0.5 * (DimSize(wTmp, 0) - 1) * DimDelta(wTmp, 0)
    Variable yCenter = DimOffset(wTmp, 1) + 0.5 * (DimSize(wTmp, 1) - 1) * DimDelta(wTmp, 1)

    Variable dx = xPos_m * 1e9 - xCenter
    Variable dy = yPos_m * 1e9 - yCenter

    SetScale/P x, DimOffset(wTmp,0) + dx, DimDelta(wTmp,0), WaveUnits(wTmp,0), wTmp
    SetScale/P y, DimOffset(wTmp,1) + dy, DimDelta(wTmp,1), WaveUnits(wTmp,1), wTmp

    Duplicate/O wTmp, $outName
    Wave wOut = $outName

    return wOut
End


//==============================================================================
// SNS_MakeWMask
//
// Code Purpose:
//   Create a binary or smoothed mask from a 2D topography/image wave.
//
// Physics Role:
//   None directly. The resulting mask is later consumed by SNS geometry and
//   ray/channel extraction routines.
//
// Required Inputs:
//   src  : 2D source wave, usually a coordinate-corrected topography.
//
// Optional Inputs:
//   zThresh  : threshold value in src units; default -0.5.
//   gaussN   : MatrixFilter Gaussian kernel size; default 3. Set <=1 to skip.
//
// Generated Outputs:
//   Creates/overwrites w_mask in the current data folder and returns it.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools MakeWMask.
//==============================================================================
Function/WAVE SNS_MakeWMask(src, [zThresh, gaussN])
    Wave src
    Variable zThresh, gaussN

    if (ParamIsDefault(zThresh))
        zThresh = -0.5
    endif
    if (ParamIsDefault(gaussN))
        gaussN = 3
    endif

    Duplicate/O src w_mask
    Wave w_mask

    w_mask = numtype(src[p][q]) == 2 ? 0 : ((src[p][q] >= zThresh) ? 1 : 0)

    if (gaussN > 1)
        MatrixFilter/N=(gaussN) gauss w_mask
    endif

    return w_mask
End


//==============================================================================
// SNS_DuplicateCutArea
//
// Code Purpose:
//   Duplicate a rectangular scaled-coordinate ROI from a 2D wave.
//
// Physics Role:
//   None. Display/export preparation helper.
//
// Required Inputs:
//   topo  : source 2D wave.
//
// Optional Inputs:
//   x1, x2, y1, y2  : scaled-coordinate ROI bounds. Defaults to full wave.
//
// Generated Outputs:
//   Creates/overwrites <sourceName>_cut in the current data folder and returns
//   the output wave.
//
// Side Effects:
//   Appends CutSourceWave and CutROI metadata to the output wave note.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools DuplicateCutArea.
//==============================================================================
Function/WAVE SNS_DuplicateCutArea(topo, [x1, x2, y1, y2])
    Wave topo
    Variable x1, x2, y1, y2

    Variable nx = DimSize(topo, 0)
    Variable ny = DimSize(topo, 1)

    Variable xOff = DimOffset(topo, 0)
    Variable xDel = DimDelta(topo, 0)
    Variable yOff = DimOffset(topo, 1)
    Variable yDel = DimDelta(topo, 1)

    Variable xEnd = xOff + xDel*(nx-1)
    Variable yEnd = yOff + yDel*(ny-1)

    if (ParamIsDefault(x1))
        x1 = xOff
    endif
    if (ParamIsDefault(x2))
        x2 = xEnd
    endif
    if (ParamIsDefault(y1))
        y1 = yOff
    endif
    if (ParamIsDefault(y2))
        y2 = yEnd
    endif

    Variable xStart = min(x1, x2)
    Variable xStop  = max(x1, x2)
    Variable yStart = min(y1, y2)
    Variable yStop  = max(y1, y2)

    String outName = CleanupName(NameOfWave(topo) + "_cut", 0)
    Duplicate/O/R=(xStart, xStop)(yStart, yStop) topo, $outName
    Wave out = $outName

    String oldNote = note(out)
    String addNote = ""

    addNote += "CutSourceWave=" + GetWavesDataFolder(topo, 2) + ";"
    addNote += "CutROI=(" + num2str(xStart) + "," + num2str(xStop) + "," + num2str(yStart) + "," + num2str(yStop) + ");"

    if (strlen(oldNote) > 0)
        if (cmpstr(oldNote[strlen(oldNote)-1, strlen(oldNote)-1], ";") != 0)
            oldNote += ";"
        endif
        Note/K out, oldNote + addNote
    else
        Note/K out, addNote
    endif

    return out
End


//==============================================================================
// SNS_MakeScaledAxisWave
//
// Code Purpose:
//   Create an explicit 1D axis wave from the x scaling of another 1D wave.
//
// Physics Role:
//   None. Data/export helper for plots and CSV exports that need explicit axis
//   values rather than implicit Igor wave scaling.
//
// Required Inputs:
//   w  : 1D source wave whose x scaling defines the output values.
//
// Optional Inputs:
//   outName  : output wave name. Default is <sourceName>_x.
//
// Generated Outputs:
//   Creates/overwrites outName in the current data folder and returns it.
//   The output data values are pnt2x(w,p), and the output data units are copied
//   from the source wave x-axis units.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools MakeXWaveFromWave.
//==============================================================================
Function/WAVE SNS_MakeScaledAxisWave(w, [outName])
    Wave w
    String outName

    Variable n = numpnts(w)
    if (n <= 0)
        Abort "SNS_MakeScaledAxisWave: input wave has no points."
    endif

    if (ParamIsDefault(outName) || strlen(outName) == 0)
        outName = CleanupName(NameOfWave(w) + "_x", 0)
    endif

    Make/O/D/N=(n) $outName
    Wave xWave = $outName

    xWave = pnt2x(w, p)

    SetScale d, 0, 0, WaveUnits(w, 0), xWave
    SetScale/P x, DimOffset(w, 0), DimDelta(w, 0), WaveUnits(w, 0), xWave

    return xWave
End


//==============================================================================
// SNS_ExtractColumnAtPositionNM
//
// Extract one energy spectrum from a 2D wave at a physical position along
// dimension 1. Linear interpolation is the default; nearest=1 selects the
// closest stored column. The output retains the dimension-0 scale and units.
//==============================================================================
Function/WAVE SNS_ExtractColumnAtPositionNM(w2D, pos_nm, [outName, nearest])
    Wave w2D
    Variable pos_nm, nearest
    String outName

    if (WaveDims(w2D) != 2)
        Abort "SNS_ExtractColumnAtPositionNM: input wave must be 2D."
    endif
    if (ParamIsDefault(outName) || strlen(outName) == 0)
        outName = CleanupName(NameOfWave(w2D) + "_cut", 0)
    endif
    if (ParamIsDefault(nearest))
        nearest = 0
    endif

    Variable q = (pos_nm - DimOffset(w2D, 1)) / DimDelta(w2D, 1)
    if (q < 0 || q > DimSize(w2D, 1) - 1)
        Abort "SNS_ExtractColumnAtPositionNM: requested position is outside dimension 1."
    endif

    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR = GetWavesDataFolderDFR(w2D)
    SetDataFolder wDFR

    Make/O/D/N=(DimSize(w2D, 0)) $outName
    Wave out = $outName
    SetScale/P x, DimOffset(w2D, 0), DimDelta(w2D, 0), WaveUnits(w2D, 0), out
    SetScale d, 0, 0, WaveUnits(w2D, -1), out

    if (nearest)
        out[] = w2D[p][round(q)]
    else
        Variable q0 = floor(q)
        Variable q1 = min(q0 + 1, DimSize(w2D, 1) - 1)
        Variable alpha = q - q0
        out[] = (1 - alpha) * w2D[p][q0] + alpha * w2D[p][q1]
    endif

    SetDataFolder oldDFR
    return out
End


//==============================================================================
// SNS_DuplicateRotatedArea
//
// Code Purpose:
//   Duplicate an optional rectangular ROI from a 2D wave and rotate it around
//   the ROI center.
//
// Physics Role:
//   None. Display/export preparation helper.
//
// Required Inputs:
//   topo  : source 2D wave.
//
// Optional Inputs:
//   x1, x2, y1, y2  : scaled-coordinate ROI bounds. Defaults to full wave.
//   angleDeg        : rotation angle [deg], counterclockwise in wave x/y frame.
//   expand          : 0 keep source ROI size; 1 expand output to fit rotation.
//
// Generated Outputs:
//   Creates/overwrites <sourceName>_rot in the current data folder and returns
//   the output wave.
//
// Side Effects:
//   Writes SourceWave, ROI, rotAngleDeg, and expand metadata to the wave note.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools DuplicateRotatedArea.
//==============================================================================
Function/WAVE SNS_DuplicateRotatedArea(topo, [x1, x2, y1, y2, angleDeg, expand])
    Wave topo
    Variable x1, x2, y1, y2, angleDeg, expand

    Variable nx = DimSize(topo, 0)
    Variable ny = DimSize(topo, 1)

    Variable xOff = DimOffset(topo, 0)
    Variable xDel = DimDelta(topo, 0)
    Variable yOff = DimOffset(topo, 1)
    Variable yDel = DimDelta(topo, 1)

    Variable xEnd = xOff + xDel*(nx-1)
    Variable yEnd = yOff + yDel*(ny-1)

    if (ParamIsDefault(x1))
        x1 = xOff
    endif
    if (ParamIsDefault(x2))
        x2 = xEnd
    endif
    if (ParamIsDefault(y1))
        y1 = yOff
    endif
    if (ParamIsDefault(y2))
        y2 = yEnd
    endif
    if (ParamIsDefault(angleDeg))
        angleDeg = 0
    endif
    if (ParamIsDefault(expand))
        expand = 0
    endif

    Variable xStart = min(x1, x2)
    Variable xStop  = max(x1, x2)
    Variable yStart = min(y1, y2)
    Variable yStop  = max(y1, y2)

    Duplicate/FREE/R=(xStart, xStop)(yStart, yStop) topo, roi

    Variable rnx = DimSize(roi, 0)
    Variable rny = DimSize(roi, 1)

    Variable rxOff = DimOffset(roi, 0)
    Variable rxDel = DimDelta(roi, 0)
    Variable ryOff = DimOffset(roi, 1)
    Variable ryDel = DimDelta(roi, 1)

    Variable rxEnd = rxOff + rxDel*(rnx-1)
    Variable ryEnd = ryOff + ryDel*(rny-1)

    Variable rxLo = min(rxOff, rxEnd), rxHi = max(rxOff, rxEnd)
    Variable ryLo = min(ryOff, ryEnd), ryHi = max(ryOff, ryEnd)

    Variable xC = (rxLo + rxHi)/2
    Variable yC = (ryLo + ryHi)/2

    Variable theta = angleDeg*pi/180
    Variable c = cos(theta)
    Variable s = sin(theta)

    Variable outNX = rnx, outNY = rny
    Variable outXOff = rxOff, outXDel = rxDel
    Variable outYOff = ryOff, outYDel = ryDel

    if (expand)
        Variable hx = (rxHi - rxLo)/2
        Variable hy = (ryHi - ryLo)/2

        Variable x1r = xC + ( +hx)*c - ( +hy)*s
        Variable y1r = yC + ( +hx)*s + ( +hy)*c

        Variable x2r = xC + ( +hx)*c - ( -hy)*s
        Variable y2r = yC + ( +hx)*s + ( -hy)*c

        Variable x3r = xC + ( -hx)*c - ( +hy)*s
        Variable y3r = yC + ( -hx)*s + ( +hy)*c

        Variable x4r = xC + ( -hx)*c - ( -hy)*s
        Variable y4r = yC + ( -hx)*s + ( -hy)*c

        Variable xMin = min(min(x1r,x2r), min(x3r,x4r))
        Variable xMax = max(max(x1r,x2r), max(x3r,x4r))
        Variable yMin = min(min(y1r,y2r), min(y3r,y4r))
        Variable yMax = max(max(y1r,y2r), max(y3r,y4r))

        Variable dx = abs(outXDel)
        Variable dy = abs(outYDel)
        if (dx <= 0)
            dx = 1
        endif
        if (dy <= 0)
            dy = 1
        endif

        outNX = ceil((xMax - xMin)/dx) + 1
        outNY = ceil((yMax - yMin)/dy) + 1
        if (outNX < 2)
            outNX = 2
        endif
        if (outNY < 2)
            outNY = 2
        endif

        if (outXDel >= 0)
            outXOff = xMin
        else
            outXOff = xMax
        endif
        if (outYDel >= 0)
            outYOff = yMin
        else
            outYOff = yMax
        endif
    endif

    String outName = CleanupName(NameOfWave(topo) + "_rot", 0)
    Make/O/N=(outNX,outNY) $outName
    Wave out = $outName

    SetScale/P x, outXOff, outXDel, WaveUnits(roi, 0), out
    SetScale/P y, outYOff, outYDel, WaveUnits(roi, 1), out
    SetScale d, 0, 0, WaveUnits(roi, -1), out

    if ((abs(angleDeg) < 1e-12) && (expand == 0))
        out = roi
    else
        out = Interp2D(roi, xC + (x-xC)*c + (y-yC)*s,  yC - (x-xC)*s + (y-yC)*c)
    endif

    String noteStr
    noteStr  = "SourceWave=" + GetWavesDataFolder(topo, 2) + ";"
    noteStr += "ROI=(" + num2str(xStart) + "," + num2str(xStop) + "," + num2str(yStart) + "," + num2str(yStop) + ");"
    noteStr += "rotAngleDeg=" + num2str(angleDeg) + ";"
    noteStr += "expand=" + num2str(expand) + ";"
    Note/K out, noteStr

    return out
End


//==============================================================================


//==============================================================================
// SNS_WaterfallFromSubFolders
//
// Code Purpose:
//   Build a 2D waterfall wave from matching STS subfolders in the current data folder.
//
// Physics Role:
//   None directly. Data-import assembly helper for measured STS sweeps.
//
// Required Inputs:
//   sname  : substring used to select source subfolders and output wave name.
//
// Optional Inputs:
//   Phi, T, ColumnRange  : parse angle, parse temperature, or set output columns.
//
// Generated Outputs:
//   Creates/overwrites wave sname in the current data folder.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools WaterfallFromSubFolders.
//==============================================================================
Function SNS_WaterfallFromSubFolders(string sname, [variable Phi, variable T, variable ColumnRange])
// Set field Phi=1 for in-plane field rotation. Parsing STS name for "Pxxx_", with xxx= Angle phi in degree 
// Set field T=1 for Temperature dependent measurement. Parsing STS name for "_xxxxmK", with xxxx= temperature in mK

    DFREF dfr = GetDataFolderDFR()
    String SubFolderName, FolderName, ListOfdIdVWaves, NameOfCurrentWave, stemp, stemp2
    String subList
    Variable i, j, vtemp, v_Bfieldrange
    Variable NumWaveInFolder
    Variable NumDataFolders
    Variable numAll

    NumWaveInFolder = 0
    FolderName = GetDataFolder(1)

    // --- build filtered list of subfolders containing sname ---
    subList = ""
    numAll = CountObjectsDFR(dfr, 4)          // 4 = data folder
    for (i = 0; i < numAll; i += 1)
        SubFolderName = GetIndexedObjNameDFR(dfr, 4, i)
        if (strsearch(SubFolderName, sname, 0) >= 0)
            subList = AddListItem(SubFolderName, subList, ";", Inf)
        endif
    endfor

    NumDataFolders = ItemsInList(subList)
    if (NumDataFolders == 0)
        Abort "SNS_WaterfallFromSubFolders: no subfolders contain \"" + sname + "\""
    endif

    Make/N=(NumDataFolders)/FREE wFieldValues

    // ---------- extract field / angle / temperature values ----------
    if (!ParamIsDefault(Phi))
        for (i = 0; i < NumDataFolders; i += 1)
            SubFolderName = StringFromList(i, subList, ";")
            vtemp = strsearch(SubFolderName, "_P", 0)
            stemp2 = SubFolderName[vtemp+2, vtemp+5]
            stemp = ReplaceString("P", stemp2, "")
            stemp2 = ReplaceString("_", stemp, "")
            vtemp = str2num(stemp2)
            wFieldValues[i] = vtemp
        endfor

    elseif (!ParamIsDefault(T))
        for (i = 0; i < NumDataFolders; i += 1)
            SubFolderName = StringFromList(i, subList, ";")
            vtemp = strsearch(SubFolderName, "mK", 0)
            stemp = SubFolderName[vtemp-5, vtemp-1]
            stemp2 = ReplaceString("m", stemp, "-")
            stemp = ReplaceString("K", stemp2, "")
            stemp2 = ReplaceString("_", stemp, "")
            vtemp = str2num(stemp2)
            wFieldValues[i] = vtemp
        endfor

    else
        for (i = 0; i < NumDataFolders; i += 1)
            SubFolderName = StringFromList(i, subList, ";")
            vtemp = strsearch(SubFolderName, "mT", 0)
            stemp = SubFolderName[vtemp-5, vtemp-1]
            stemp2 = ReplaceString("m", stemp, "-")
            stemp = ReplaceString("B", stemp2, "")
            stemp2 = ReplaceString("_", stemp, "")
            stemp = ReplaceString("y", stemp2, "")
            vtemp = str2num(stemp)
            wFieldValues[i] = vtemp
        endfor
    endif

    Sort wFieldValues, wFieldValues

    // ---------- determine v_Bfieldrange ----------
    if (ParamIsDefault(ColumnRange))
        v_Bfieldrange = NumDataFolders
    else
        v_Bfieldrange = ColumnRange
    endif

    vtemp = (-wFieldValues[0] + wFieldValues[DimSize(wFieldValues, 0)-1]) / (v_Bfieldrange - 1)

    do
        if (vtemp > round(vtemp))
            v_Bfieldrange += 1
            vtemp = (-wFieldValues[0] + wFieldValues[DimSize(wFieldValues, 0)-1]) / (v_Bfieldrange - 1)
        elseif (vtemp < round(vtemp))
            v_Bfieldrange -= 1
            vtemp = (-wFieldValues[0] + wFieldValues[DimSize(wFieldValues, 0)-1]) / (v_Bfieldrange - 1)
        endif
    while (!vtemp == trunc(vtemp))

    // ---------- build waterfall ----------
    if (!ParamIsDefault(Phi))

        for (i = 0; i < NumDataFolders; i += 1)
            SubFolderName = StringFromList(i, subList, ";")

            if (i == 0)
                SetDataFolder ":'" + SubFolderName + "':"
                ListOfdIdVWaves = WaveList("*LI_Demod_1_X*", ";", "")
                ListOfdIdVWaves = ListMatch(ListOfdIdVWaves, "*avg*")
                vtemp = DimSize($StringFromList(0, ListOfdIdVWaves, ";"), 0)
                SetDataFolder FolderName

                Make/FREE/D/O/N=(vtemp, v_Bfieldrange) wout = NaN
                NameOfCurrentWave = ":'" + SubFolderName + "':'"
                NameOfCurrentWave += StringFromList(0, ListOfdIdVWaves, ";")
                NameOfCurrentWave += "'"

                CopyScales/P $NameOfCurrentWave, wout
                SetScale/I y, wFieldValues[0], wFieldValues[DimSize(wFieldValues, 0)-1], "deg", wout

                vtemp = strsearch(SubFolderName, "_P", 0)
                stemp2 = SubFolderName[vtemp+2, vtemp+5]
                stemp = ReplaceString("P", stemp2, "")
                stemp2 = ReplaceString("_", stemp, "")
                vtemp = str2num(stemp2)

                Wave w = $NameOfCurrentWave
                wout[][scaleToIndex(wout, vtemp, 1)] = w[p][0]

            else
                SetDataFolder ":'" + SubFolderName + "':"
                ListOfdIdVWaves = WaveList("*LI_Demod_1_X*", ";", "")
                ListOfdIdVWaves = ListMatch(ListOfdIdVWaves, "*avg*")
                vtemp = DimSize($StringFromList(0, ListOfdIdVWaves, ";"), 0)
                SetDataFolder FolderName

                NameOfCurrentWave = ":'" + SubFolderName + "':'"
                NameOfCurrentWave += StringFromList(0, ListOfdIdVWaves, ";")
                NameOfCurrentWave += "'"

                vtemp = strsearch(SubFolderName, "_P", 0)
                stemp2 = SubFolderName[vtemp+2, vtemp+5]
                stemp = ReplaceString("P", stemp2, "")
                stemp2 = ReplaceString("_", stemp, "")
                vtemp = str2num(stemp2)

                Wave w2 = $NameOfCurrentWave
                wout[][scaleToIndex(wout, vtemp, 1)] = w2[p][0]
            endif
        endfor

    elseif (!ParamIsDefault(T))

        for (i = 0; i < NumDataFolders; i += 1)
            SubFolderName = StringFromList(i, subList, ";")

            if (i == 0)
                SetDataFolder ":'" + SubFolderName + "':"
                ListOfdIdVWaves = WaveList("*LI_Demod_1_X*", ";", "")
                ListOfdIdVWaves = ListMatch(ListOfdIdVWaves, "*avg*")
                vtemp = DimSize($StringFromList(0, ListOfdIdVWaves, ";"), 0)
                SetDataFolder FolderName

                Make/FREE/D/O/N=(vtemp, v_Bfieldrange) wout = NaN
                NameOfCurrentWave = ":'" + SubFolderName + "':'"
                NameOfCurrentWave += StringFromList(0, ListOfdIdVWaves, ";")
                NameOfCurrentWave += "'"

                CopyScales/P $NameOfCurrentWave, wout
                SetScale/I y, wFieldValues[0], wFieldValues[DimSize(wFieldValues, 0)-1], "mK", wout

                vtemp = strsearch(SubFolderName, "mK", 0)
                stemp = SubFolderName[vtemp-5, vtemp-1]
                stemp2 = ReplaceString("m", stemp, "-")
                stemp = ReplaceString("K", stemp2, "")
                stemp2 = ReplaceString("_", stemp, "")
                vtemp = str2num(stemp2)

                Wave w3 = $NameOfCurrentWave
                wout[][scaleToIndex(wout, vtemp, 1)] = w3[p][0]

            else
                SetDataFolder ":'" + SubFolderName + "':"
                ListOfdIdVWaves = WaveList("*LI_Demod_1_X*", ";", "")
                ListOfdIdVWaves = ListMatch(ListOfdIdVWaves, "*avg*")
                vtemp = DimSize($StringFromList(0, ListOfdIdVWaves, ";"), 0)
                SetDataFolder FolderName

                NameOfCurrentWave = ":'" + SubFolderName + "':'"
                NameOfCurrentWave += StringFromList(0, ListOfdIdVWaves, ";")
                NameOfCurrentWave += "'"

                vtemp = strsearch(SubFolderName, "mK", 0)
                stemp = SubFolderName[vtemp-5, vtemp-1]
                stemp2 = ReplaceString("m", stemp, "-")
                stemp = ReplaceString("K", stemp2, "")
                stemp2 = ReplaceString("_", stemp, "")
                vtemp = str2num(stemp2)

                Wave w4 = $NameOfCurrentWave
                wout[][scaleToIndex(wout, vtemp, 1)] = w4[p][0]
            endif
        endfor

    else

        for (i = 0; i < NumDataFolders; i += 1)
            SubFolderName = StringFromList(i, subList, ";")

            if (i == 0)
                SetDataFolder ":'" + SubFolderName + "':"
                ListOfdIdVWaves = WaveList("*LI_Demod_1_X*", ";", "")
                ListOfdIdVWaves = ListMatch(ListOfdIdVWaves, "*avg*")
                vtemp = DimSize($StringFromList(0, ListOfdIdVWaves, ";"), 0)
                SetDataFolder FolderName

                Make/FREE/D/O/N=(vtemp, v_Bfieldrange) wout = NaN
                NameOfCurrentWave = ":'" + SubFolderName + "':'"
                NameOfCurrentWave += StringFromList(0, ListOfdIdVWaves, ";")
                NameOfCurrentWave += "'"

                CopyScales/P $NameOfCurrentWave, wout
                SetScale/I y, wFieldValues[0], wFieldValues[DimSize(wFieldValues, 0)-1], "mT", wout

                vtemp = strsearch(SubFolderName, "mT", 0)
                stemp = SubFolderName[vtemp-5, vtemp-1]
                stemp2 = ReplaceString("m", stemp, "-")
                stemp = ReplaceString("B", stemp2, "")
                stemp2 = ReplaceString("_", stemp, "")
                stemp = ReplaceString("y", stemp2, "")
                vtemp = str2num(stemp)

                Wave w5 = $NameOfCurrentWave
                wout[][scaleToIndex(wout, vtemp, 1)] = w5[p][0]

            else
                SetDataFolder ":'" + SubFolderName + "':"
                ListOfdIdVWaves = WaveList("*LI_Demod_1_X*", ";", "")
                ListOfdIdVWaves = ListMatch(ListOfdIdVWaves, "*avg*")
                vtemp = DimSize($StringFromList(0, ListOfdIdVWaves, ";"), 0)
                SetDataFolder FolderName

                NameOfCurrentWave = ":'" + SubFolderName + "':'"
                NameOfCurrentWave += StringFromList(0, ListOfdIdVWaves, ";")
                NameOfCurrentWave += "'"

                vtemp = strsearch(SubFolderName, "mT", 0)
                stemp = SubFolderName[vtemp-5, vtemp-1]
                stemp2 = ReplaceString("m", stemp, "-")
                stemp = ReplaceString("B", stemp2, "")
                stemp2 = ReplaceString("_", stemp, "")
                stemp = ReplaceString("y", stemp2, "")
                vtemp = str2num(stemp)

                Wave w6 = $NameOfCurrentWave
                wout[][scaleToIndex(wout, vtemp, 1)] = w6[p][0]
            endif
        endfor
    endif

    Duplicate/O wout, $sname
End

//==============================================================================
// SNS_NormalizeGridXY_PixelByPixel
//
// Normalize each spectrum in a 3D GridSTS wave by its value at refIndex.
// Invalid or near-zero reference values produce NaN spectra.
//==============================================================================
Function SNS_NormalizeGridXY_PixelByPixel(grid, [refIndex, outName, eps])
    Wave grid
    Variable refIndex, eps
    String outName

    if (WaveDims(grid) != 3)
        Abort "SNS_NormalizeGridXY_PixelByPixel: grid must be three-dimensional."
    endif
    if (ParamIsDefault(refIndex))
        refIndex = 1
    endif
    if (ParamIsDefault(outName) || strlen(outName) == 0)
        outName = NameOfWave(grid) + "_nPixRef"
    endif
    if (ParamIsDefault(eps))
        eps = 1e-30
    endif

    Variable nx = DimSize(grid, 0)
    Variable ny = DimSize(grid, 1)
    Variable nz = DimSize(grid, 2)
    if (refIndex < 0 || refIndex >= nz)
        Abort "SNS_NormalizeGridXY_PixelByPixel: refIndex outside layer range."
    endif

    Duplicate/O/FREE grid, normalized
    Variable ix, iy, iz, refVal
    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)
            refVal = grid[ix][iy][refIndex]
            if (numtype(refVal) == 0 && abs(refVal) > eps)
                for (iz = 0; iz < nz; iz += 1)
                    normalized[ix][iy][iz] = grid[ix][iy][iz] / refVal
                endfor
            else
                for (iz = 0; iz < nz; iz += 1)
                    normalized[ix][iy][iz] = NaN
                endfor
            endif
        endfor
    endfor

    SetScale d, 0, 0, "", normalized
    Duplicate/O normalized, $outName
    return 0
End

//==============================================================================
// SNS_NormalizeWaterFall
//
// Code Purpose:
//   Normalize each spectrum in a waterfall by a fitted linear background outside a central range.
//
// Physics Role:
//   None directly. Experimental-data normalization helper.
//
// Required Inputs:
//   waterfall  : 1D or 2D conductance/bias wave.
//
// Optional Inputs:
//   v_maskrange  : central excluded range in the wave x-axis units; default 2.
//
// Generated Outputs:
//   Creates/overwrites <waterfall>_n.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools NormalizeWaterFall.
//==============================================================================
Function SNS_NormalizeWaterFall(wave waterfall [, variable v_maskrange])

if(paramisDefault(v_maskrange))
	v_maskrange=2 // in the wave x-axis units
endif

KillWaves/Z w_sigma

duplicate/O/FREE/R=[][0] waterfall w1
Redimension/N=-1 w1
duplicate/FREE/O w1 W_mask

make/N=2/D/O/FREE mask={x2pnt(w1,-v_maskrange),x2pnt(w1,v_maskrange)}
sort mask, mask

if(mask[0]<0)
	mask[0]=0
endif
if(mask[1]>(DimSize(w1,0)-1))
	mask[1]=(DimSize(w1,0)-1)
endif	
		if(v_maskrange>0)
			W_mask=1
			W_mask[mask[0],mask[1]]=0	
		else
			W_mask=0
			W_mask[mask[0],mask[1]]=1
		endif
	
			
string s_out=nameofWave(waterfall)+"_n"

duplicate/O/FREE waterfall wout
make/N=2/D/O/FREE w_fit_coef

variable i,v_biasoffset

if(DimSize(waterfall,1)==0)
		duplicate/O/FREE waterfall w1

		CurveFit/X=1/Q/w=2 line kwCWave=w_fit_coef w1 /M=W_mask
		w1/=w_fit_coef[0]+w_fit_coef[1]*x;
		wout[][i]=w1[p];
else
	for(i=0;i<DimSize(waterfall,1);i++)
		duplicate/O/FREE/R=[][i] waterfall w1
		Redimension/N=-1 w1
		if(numtype(w1[0][0])==0)
			CurveFit/X=1/Q/w=2 line kwCWave=w_fit_coef w1 /M=W_mask
			w1/=w_fit_coef[0]+w_fit_coef[1]*x;
			wout[][i]=w1[p];
		endif
	endfor
endif
SetScale d 0,0,"",wout

duplicate/O wout $s_out
KillWaves/Z w_sigma

END

//==============================================================================
// SNS_ScaleWaveMaxToReferenceMax
//
// Code Purpose:
//   Multiply a target wave so its finite global maximum equals the finite
//   global maximum of a reference wave, and return the applied scale factor.
//
// Physics Role:
//   None. Data-derived display/source-wave amplitude normalization.
//
// Required Inputs:
//   targetWave    : wave to scale in place.
//   referenceWave : wave whose finite global maximum defines the target value.
//
// Returns:
//   Applied scale factor, referenceMaximum/targetMaximum.
//
// Failure:
//   Aborts if either wave has no finite points, either maximum is non-finite,
//   or either maximum is not positive.
//==============================================================================
Function/D SNS_ScaleWaveMaxToReferenceMax(targetWave, referenceWave)
    Wave targetWave, referenceWave

    WaveStats/Q targetWave
    Variable targetMaximum = V_max
    Variable targetFinitePoints = V_npnts

    WaveStats/Q referenceWave
    Variable referenceMaximum = V_max
    Variable referenceFinitePoints = V_npnts

    if (targetFinitePoints <= 0 || numtype(targetMaximum) != 0 || targetMaximum <= 0)
        Abort "SNS_ScaleWaveMaxToReferenceMax: target wave has no positive finite maximum."
    endif
    if (referenceFinitePoints <= 0 || numtype(referenceMaximum) != 0 || referenceMaximum <= 0)
        Abort "SNS_ScaleWaveMaxToReferenceMax: reference wave has no positive finite maximum."
    endif

    Variable scaleFactor = referenceMaximum/targetMaximum
    targetWave *= scaleFactor
    return scaleFactor
End

//==============================================================================
// SNS_DetermineSCGap
//
// Code Purpose:
//   Estimate gap and bias-offset traces from derivative extrema and local Lorentzian fits.
//
// Physics Role:
//   Experimental gap-position diagnostic; does not modify SNS model parameters.
//
// Required Inputs:
//   waterfall  : normalized 2D waterfall wave.
//
// Optional Inputs:
//   v_smoothrange  : interpolation smoothing range; default 10/DimSize(dim0).
//
// Generated Outputs:
//   <waterfall>_SCgap, <waterfall>_SCgap_fromfit, <waterfall>_biasoffset_fromfit, <waterfall>_biasoffset.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools DetermineSCGap.
//==============================================================================
Function SNS_DetermineSCGap(wave waterfall [variable v_smoothrange])

KillWaves/Z SNS_tmp_SCGap_SS, SNS_tmp_SCGap_DIF
KillWaves/Z fit_SNS_tmp_SCGap_DIF, w_sigma

make/O/FREE/N=(DimSize(waterfall,1)) wout
make/O/FREE/N=(DimSize(waterfall,1)) wout_offset
make/O/FREE/N=(DimSize(waterfall,1)) wout_fromfit
make/O/FREE/N=(DimSize(waterfall,1)) wout_offsetfromfit
SetScale/P x DimOffset(waterfall,1),DimDelta(waterfall,1),waveunits(waterfall,1),wout
SetScale/P x DimOffset(waterfall,1),DimDelta(waterfall,1),waveunits(waterfall,1),wout_offset
SetScale/P x DimOffset(waterfall,1),DimDelta(waterfall,1),waveunits(waterfall,1),wout_fromfit
SetScale/P x DimOffset(waterfall,1),DimDelta(waterfall,1),waveunits(waterfall,1),wout_offsetfromfit

make/O/FREE/N=(4) w_coef_lor
string s_out=nameofWave(waterfall)+"_SCgap"

variable i,v_SCgap


if(paramisDefault(v_smoothrange))
	v_smoothrange=10/DimSize(waterfall,0) // in mV
endif

for(i=0;i<DimSize(waterfall,1);i++)
	duplicate/O/FREE/R=[][i] waterfall w1
	Redimension/N=-1 w1
	Interpolate2/T=3/N=(DimSize(w1,0)*10)/F=(v_smoothrange)/Y=SNS_tmp_SCGap_SS w1;
	Differentiate SNS_tmp_SCGap_SS/D=SNS_tmp_SCGap_DIF;DelayUpdate
	SNS_tmp_SCGap_DIF=x<0&& SNS_tmp_SCGap_DIF(x) > 0 ? Nan : SNS_tmp_SCGap_DIF
	SNS_tmp_SCGap_DIF=x>0&& SNS_tmp_SCGap_DIF(x) < 0 ? Nan : SNS_tmp_SCGap_DIF
	wavestats/Q SNS_tmp_SCGap_DIF	
	wout[i]=(v_maxloc-v_minloc)/2
	wout_offset[i]=v_maxloc-(v_maxloc-v_minloc)/2
	duplicate/O/FREE/R=[][i] waterfall w1
	Redimension/N=-1 w1
	Interpolate2/T=3/N=(DimSize(w1,0)*10)/F=(v_smoothrange)/Y=SNS_tmp_SCGap_SS w1;
	Differentiate SNS_tmp_SCGap_SS/D=SNS_tmp_SCGap_DIF;
	CurveFit/Q/w=2 lor, kwCWave=w_coef_lor, SNS_tmp_SCGap_DIF[x2pnt(SNS_tmp_SCGap_DIF,v_maxloc)-10,x2pnt(SNS_tmp_SCGap_DIF,v_maxloc)+10]
	v_SCgap=w_coef_lor[2];
	CurveFit/Q/w=2 lor, kwCWave=w_coef_lor, SNS_tmp_SCGap_DIF[x2pnt(SNS_tmp_SCGap_DIF,v_minloc)-10,x2pnt(SNS_tmp_SCGap_DIF,v_minloc)+10]
	v_SCgap-=w_coef_lor[2];
	v_SCgap/=2;
	wout_offsetfromfit[i]=w_coef_lor[2]+v_SCgap
	wout_fromfit[i]=v_SCgap	
endfor

duplicate/O wout $s_out
s_out=nameofWave(waterfall)+"_SCgap_fromfit"
duplicate/O wout_fromfit $s_out
s_out=nameofWave(waterfall)+"_biasoffset_fromfit"
duplicate/O wout_offsetfromfit $s_out
s_out=nameofWave(waterfall)+"_biasoffset"
duplicate/O wout_offset $s_out


KillWaves/Z SNS_tmp_SCGap_SS, SNS_tmp_SCGap_DIF
KillWaves/Z fit_SNS_tmp_SCGap_DIF, w_sigma



END

//==============================================================================
// SNS_ConvertConductance_nS_to_S
//
// Code Purpose:
//   Convert conductance wave data units from nS to S in place.
//
// Physics Role:
//   Unit hygiene for experimental conductance data.
//
// Required Inputs:
//   w  : conductance wave with data unit nS or already S.
//
// Generated Outputs:
//   Modifies w in place and updates its data scale unit to S.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools ConvertConductance_nS_to_S.
//==============================================================================
Function SNS_ConvertConductance_nS_to_S(w)
	Wave w

	String oldUnit = WaveUnits(w, -1)   // data unit ("d" scale)

	if (cmpstr(oldUnit, "S") == 0)
		return 0
	endif

	if (cmpstr(oldUnit, "nS") != 0)
		Abort "SNS_ConvertConductance_nS_to_S: expected data unit 'nS', but found '"+oldUnit+"'. Conversion aborted."
	endif

	// Convert values
	w /= 1e9

	// Recompute min/max after conversion
	WaveStats/Q w

	// Update data scaling and unit
	SetScale d, V_min, V_max, "S", w

	return 1
End

//==============================================================================
// SNS_MakeAxisCorrDuplicate
//
// Code Purpose:
//   Duplicate a 1D/2D wave and shift one or both axis offsets by supplied symmetry offsets.
//
// Physics Role:
//   None directly. Axis-registration helper for experimental data.
//
// Required Inputs:
//   wIn, offset  : source wave and primary axis offset correction.
//
// Optional Inputs:
//   wDim, offsetY  : primary dimension and optional correction of the other dimension.
//
// Generated Outputs:
//   Creates/overwrites <wIn>_corr in the source wave folder.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools MakeAxisCorrDuplicate.
//==============================================================================
Function SNS_MakeAxisCorrDuplicate(wIn, offset, [wDim, offsetY])
    Wave wIn
    Variable offset
    Variable wDim, offsetY

    Variable dims = WaveDims(wIn)
    if (dims != 1 && dims != 2)
        Abort "SNS_MakeAxisCorrDuplicate: wIn must be 1D or 2D."
    endif

    if (ParamIsDefault(wDim))
        wDim = 0
    endif
    wDim = round(wDim)

    if (dims == 1)
        if (wDim != 0)
            Abort "SNS_MakeAxisCorrDuplicate: for 1D waves, wDim must be 0."
        endif
    else
        if (wDim != 0 && wDim != 1)
            Abort "SNS_MakeAxisCorrDuplicate: for 2D waves, wDim must be 0 or 1."
        endif
    endif

    if (ParamIsDefault(offsetY))
        offsetY = 0
    endif

    // Work in folder of wIn
    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR   = GetWavesDataFolderDFR(wIn)
    SetDataFolder wDFR

    String base = CleanupName(NameOfWave(wIn), 0)
    String outName = CleanupName(base + "_corr", 0)

    Duplicate/O wIn, $outName
    WAVE wOut = $outName

    // Helper function-like block to apply one axis shift
    // (Igor has no local functions, so we do it inline twice)
    Variable oldOff, dlt, newOff
    String u

    // ---- Apply primary correction (offset) if nonzero ----
    if (offset != 0)
        oldOff = DimOffset(wOut, wDim)
        dlt    = DimDelta(wOut, wDim)
        u      = WaveUnits(wOut, wDim)
        newOff = oldOff - offset

        if (wDim == 0)
            SetScale/P x, newOff, dlt, u, wOut
        else
            SetScale/P y, newOff, dlt, u, wOut
        endif
    endif

    // ---- Apply secondary correction (offsetY) to the other dim if requested ----
    if (dims == 2 && offsetY != 0)
        Variable otherDim = 1 - wDim

        oldOff = DimOffset(wOut, otherDim)
        dlt    = DimDelta(wOut, otherDim)
        u      = WaveUnits(wOut, otherDim)
        newOff = oldOff - offsetY

        if (otherDim == 0)
            SetScale/P x, newOff, dlt, u, wOut
        else
            SetScale/P y, newOff, dlt, u, wOut
        endif
    endif

    // Add a note describing what happened
    String noteStr = "Axis corrected copy of "+NameOfWave(wIn)+": "
    noteStr += "wDim="+num2istr(wDim)+"; offset="+num2str(offset)

    if (dims == 2)
        noteStr += "; offsetY="+num2str(offsetY)
    endif

    noteStr += "; rule: newOffset = oldOffset - offset (per dim)."
    Note wOut, noteStr

    SetDataFolder oldDFR
End

//==============================================================================
// SNS_FlipB
//
// Code Purpose:
//   Reverse a wave along its magnetic-field dimension and negate the field-axis
//   coordinates while preserving point spacing and units.
//
// Required Inputs:
//   w  : 1D wave with B on dim 0, or 2D wave with B on dim 1.
//
// Generated Outputs:
//   Creates/overwrites <w>_Bflip in the source wave folder.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools FlipB.
//==============================================================================
Function SNS_FlipB(w)
    Wave w

    Variable dims = WaveDims(w)
    if (dims != 1 && dims != 2)
        Abort "SNS_FlipB: input must be a one- or two-dimensional wave."
    endif

    String outPath = GetWavesDataFolder(w, 1) + NameOfWave(w) + "_Bflip"
    Duplicate/O w, $outPath
    Wave wFlip = $outPath

    Variable n, b0, db

    if (dims == 1)
        n = DimSize(w, 0)
        b0 = DimOffset(w, 0)
        db = DimDelta(w, 0)

        wFlip[] = w[n - 1 - p]
        SetScale/P x, -b0 - db*(n - 1), db, WaveUnits(w, 0), wFlip
    else
        n = DimSize(w, 1)
        b0 = DimOffset(w, 1)
        db = DimDelta(w, 1)

        wFlip[][] = w[p][n - 1 - q]
        SetScale/P y, -b0 - db*(n - 1), db, WaveUnits(w, 1), wFlip
    endif
End

//==============================================================================
// SNS__FindSymCenterScan
//
// Code Purpose:
//   Internal helper for SNS_FindOffsetBySymmetry. Finds the best 1D symmetry center.
//==============================================================================
Static Function SNS__FindSymCenterScan(w1, resW)
    Wave w1
    Wave resW

    Redimension/N=3 resW
    resW = NaN
    resW[0] = -1

    if (WaveDims(w1) != 1)
        return 0
    endif

    Variable n = numpnts(w1)
    if (n < 5)
        return 0
    endif

    Variable dx = DimDelta(w1, 0)
    if (dx == 0)
        dx = 1
    endif
    Variable dxAbs = abs(dx)

    Variable rangeCoord = abs(pnt2x(w1, n-1) - pnt2x(w1, 0))
    Variable Bw = 0.25*rangeCoord
    if (Bw <= 0)
        Bw = dxAbs
    endif

    Variable span = n - 1
    Variable minDelta = max(2, round(0.1*span))
    minDelta = min(minDelta, floor(span/2))
    if (minDelta < 1)
        minDelta = 1
    endif

    Variable bestC = -1
    Variable bestScore = 1e308
    Variable bestPairs = 0

    Variable c, delta, deltaMax
    Variable sum, wsum, score, w, dd
    Variable pairsUsed
    Variable a, b

    for (c = minDelta; c <= (n-1-minDelta); c += 1)

        deltaMax = min(c, (n-1) - c)
        if (deltaMax < minDelta)
            continue
        endif

        sum = 0
        wsum = 0
        pairsUsed = 0

        for (delta = 1; delta <= deltaMax; delta += 1)

            a = w1[c + delta]
            b = w1[c - delta]
            if ((numtype(a) != 0) || (numtype(b) != 0))
                continue
            endif

            dd = delta * dxAbs
            w = 1 / (1 + (dd/Bw)*(dd/Bw))     // emphasize near center

            sum  += w * abs(a - b)
            wsum += w
            pairsUsed += 1
        endfor

        if (pairsUsed < minDelta || wsum <= 0)
            continue
        endif

        score = sum / wsum
        if (score < bestScore)
            bestScore = score
            bestC = c
            bestPairs = pairsUsed
        endif
    endfor

    if (bestC >= 0)
        resW[0] = bestC
        resW[1] = bestScore
        resW[2] = bestPairs
    endif

    return 0
End

//==============================================================================
// SNS_FindOffsetBySymmetry
//
// Code Purpose:
//   Estimate x/y symmetry centers of a 1D/2D wave within an optional ROI.
//
// Physics Role:
//   None directly. Alignment helper for experimental STS maps.
//
// Required Inputs:
//   wIn  : 1D or 2D source wave.
//
// Optional Inputs:
//   xNegCut, xPosCut, yNegCut, yPosCut  : scaled-coordinate ROI limits.
//
// Generated Outputs:
//   <wave>_symROI and scalar center/shift/quality diagnostics next to wIn.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools FindOffsetBySymmetry.
//==============================================================================
Function SNS_FindOffsetBySymmetry(wIn, [xNegCut, xPosCut, yNegCut, yPosCut])
    Wave wIn
    Variable xNegCut, xPosCut, yNegCut, yPosCut

    Variable dims = WaveDims(wIn)
    if (dims != 1 && dims != 2)
        Abort "SNS_FindOffsetBySymmetry: wIn must be 1D or 2D."
    endif

    // Work in folder of wIn
    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR   = GetWavesDataFolderDFR(wIn)
    SetDataFolder wDFR

    String base   = CleanupName(NameOfWave(wIn), 0)
    String roiName = CleanupName(base + "_symROI", 0)

    // -------- locals declared ONCE (Igor has no block scope) --------
    Variable nX, nY
    Variable xOff, dx, xFirst, xLast, xWaveMin, xWaveMax, xCutMin, xCutMax
    Variable yOff, dy, yFirst, yLast, yWaveMin, yWaveMax, yCutMin, yCutMax
    Variable pA, pB, pLo, pHi, eps, tmp
    Variable iMinX, iMaxX, jMinY, jMaxY
    Variable xUseMin, xUseMax, yUseMin, yUseMax
    Variable xSymVal, ySymVal, xScoreVal, yScoreVal
    Variable nUsedX, nUsedY

    // --------------------------
    // X cut (dim0) -> indices iMinX..iMaxX
    // --------------------------
    nX = DimSize(wIn, 0)
    if (nX < 5)
        SetDataFolder oldDFR
        Abort "SNS_FindOffsetBySymmetry: need at least 5 points in dim0."
    endif

    xOff = DimOffset(wIn, 0)
    dx   = DimDelta(wIn, 0)
    if (dx == 0)
        dx = 1
    endif

    xFirst = xOff
    xLast  = xOff + (nX-1)*dx
    xWaveMin = min(xFirst, xLast)
    xWaveMax = max(xFirst, xLast)

    xCutMin = xWaveMin
    xCutMax = xWaveMax
    if (!ParamIsDefault(xNegCut))
        xCutMin = xNegCut
    endif
    if (!ParamIsDefault(xPosCut))
        xCutMax = xPosCut
    endif
    if (xCutMin > xCutMax)
        tmp = xCutMin; xCutMin = xCutMax; xCutMax = tmp
    endif
    xCutMin = max(xCutMin, xWaveMin)
    xCutMax = min(xCutMax, xWaveMax)

    pA = (xCutMin - xOff)/dx
    pB = (xCutMax - xOff)/dx
    pLo = min(pA, pB)
    pHi = max(pA, pB)
    eps = 1e-9

    iMinX = floor(pLo + 1 - eps)   // ceil(pLo)
    iMaxX = floor(pHi)             // floor(pHi)

    iMinX = max(0, iMinX)
    iMaxX = min(nX-1, iMaxX)

    if (iMaxX - iMinX < 4)
        SetDataFolder oldDFR
        Abort "SNS_FindOffsetBySymmetry: X cut too small (need >=5 points)."
    endif

    xUseMin = min(xOff + iMinX*dx, xOff + iMaxX*dx)
    xUseMax = max(xOff + iMinX*dx, xOff + iMaxX*dx)

    // --------------------------
    // Y cut (dim1) -> indices jMinY..jMaxY (2D only)
    // --------------------------
    if (dims == 2)
        nY = DimSize(wIn, 1)
        if (nY < 5)
            SetDataFolder oldDFR
            Abort "SNS_FindOffsetBySymmetry: need at least 5 points in dim1."
        endif

        yOff = DimOffset(wIn, 1)
        dy   = DimDelta(wIn, 1)
        if (dy == 0)
            dy = 1
        endif

        yFirst = yOff
        yLast  = yOff + (nY-1)*dy
        yWaveMin = min(yFirst, yLast)
        yWaveMax = max(yFirst, yLast)

        yCutMin = yWaveMin
        yCutMax = yWaveMax
        if (!ParamIsDefault(yNegCut))
            yCutMin = yNegCut
        endif
        if (!ParamIsDefault(yPosCut))
            yCutMax = yPosCut
        endif
        if (yCutMin > yCutMax)
            tmp = yCutMin; yCutMin = yCutMax; yCutMax = tmp
        endif
        yCutMin = max(yCutMin, yWaveMin)
        yCutMax = min(yCutMax, yWaveMax)

        pA = (yCutMin - yOff)/dy
        pB = (yCutMax - yOff)/dy
        pLo = min(pA, pB)
        pHi = max(pA, pB)

        jMinY = floor(pLo + 1 - eps)
        jMaxY = floor(pHi)

        jMinY = max(0, jMinY)
        jMaxY = min(nY-1, jMaxY)

        if (jMaxY - jMinY < 4)
            SetDataFolder oldDFR
            Abort "SNS_FindOffsetBySymmetry: Y cut too small (need >=5 points)."
        endif

        yUseMin = min(yOff + jMinY*dy, yOff + jMaxY*dy)
        yUseMax = max(yOff + jMinY*dy, yOff + jMaxY*dy)
    else
        // placeholders for 1D
        yUseMin = NaN
        yUseMax = NaN
    endif

    // --------------------------
    // Create ROI wave for visibility
    // --------------------------
    if (dims == 1)
        Duplicate/O/R=[iMinX, iMaxX] wIn, $roiName
    else
        Duplicate/O/R=[iMinX, iMaxX][jMinY, jMaxY] wIn, $roiName
    endif
    Wave wROI = $roiName

    // --------------------------
    // 1D: compute only X center on ROI
    // --------------------------
    if (dims == 1)
        Make/FREE/N=3 res
        SNS__FindSymCenterScan(wROI, res)

        if (res[0] < 0)
            SetDataFolder oldDFR
            Abort "SNS_FindOffsetBySymmetry: could not find symmetry center in 1D wave."
        endif

        xSymVal   = pnt2x(wROI, res[0])
        xScoreVal = res[1]
        nUsedX    = 1

        Variable/G $(base + "_dimXSym")      = xSymVal
        Variable/G $(base + "_dimXShiftSym") = -xSymVal
        Variable/G $(base + "_dimXSymScore") = xScoreVal
        Variable/G $(base + "_dimXSymPairs") = nUsedX
        Variable/G $(base + "_dimXUseMin")   = xUseMin
        Variable/G $(base + "_dimXUseMax")   = xUseMax

        // ensure no Y outputs remain for 1D case
        String v
        v = base + "_dimYSym";      KillVariables/Z $v
        v = base + "_dimYShiftSym"; KillVariables/Z $v
        v = base + "_dimYSymScore"; KillVariables/Z $v
        v = base + "_dimYSymPairs"; KillVariables/Z $v
        v = base + "_dimYUseMin";   KillVariables/Z $v
        v = base + "_dimYUseMax";   KillVariables/Z $v

        SetDataFolder oldDFR
        return 0
    endif

    // --------------------------
    // 2D: slice-by-slice averaging
    // --------------------------
    Variable nXR = DimSize(wROI, 0)
    Variable nYR = DimSize(wROI, 1)

    // subsample number of slices for speed (keeps <= ~200 slices)
    Variable maxSlices = 200
    Variable stepRows = max(1, floor(nXR / maxSlices))   // for dimY center: run over rows
    Variable stepCols = max(1, floor(nYR / maxSlices))   // for dimX center: run over cols

    // slice buffers
    Make/FREE/N=3 resS
    Make/FREE/N=(nYR) wSliceY
    Make/FREE/N=(nXR) wSliceX

    // set scales once
    SetScale/P x, DimOffset(wROI, 1), DimDelta(wROI, 1), WaveUnits(wROI, 1), wSliceY
    SetScale/P x, DimOffset(wROI, 0), DimDelta(wROI, 0), WaveUnits(wROI, 0), wSliceX

    // stats + accumulation
    Variable row, col, val
    Variable s, ss, cnt, mean, var, std, wgt
    Variable sumWY, sumY, sumScoreY
    Variable sumWX, sumX, sumScoreX

    // ---- dimYSym: for each ROW slice (fixed dim0), find symmetry center along dim1 ----
    sumWY = 0; sumY = 0; sumScoreY = 0
    nUsedY = 0

    for (row = 0; row < nXR; row += stepRows)

        s = 0; ss = 0; cnt = 0
        for (col = 0; col < nYR; col += 1)
            val = wROI[row][col]
            wSliceY[col] = val
            if (numtype(val) == 0)
                s  += val
                ss += val*val
                cnt += 1
            endif
        endfor

        if (cnt < 5)
            continue
        endif

        mean = s / cnt
        var  = ss / cnt - mean*mean
        if (var < 0)
            var = 0
        endif
        std = sqrt(var)

        // skip nearly-flat slices (don't constrain symmetry)
        if (std <= 1e-9*max(abs(mean), 1e-12))
            continue
        endif

        SNS__FindSymCenterScan(wSliceY, resS)
        if (resS[0] < 0)
            continue
        endif

        ySymVal = pnt2x(wSliceY, resS[0])
        wgt = std * max(resS[2], 1)

        sumWY     += wgt
        sumY      += wgt * ySymVal
        sumScoreY += wgt * resS[1]
        nUsedY    += 1
    endfor

    if (sumWY <= 0 || nUsedY < 1)
        SetDataFolder oldDFR
        Abort "SNS_FindOffsetBySymmetry: could not determine dimY symmetry center (no usable row slices)."
    endif

    ySymVal   = sumY / sumWY
    yScoreVal = sumScoreY / sumWY

    // ---- dimXSym: for each COL slice (fixed dim1), find symmetry center along dim0 ----
    sumWX = 0; sumX = 0; sumScoreX = 0
    nUsedX = 0

    for (col = 0; col < nYR; col += stepCols)

        s = 0; ss = 0; cnt = 0
        for (row = 0; row < nXR; row += 1)
            val = wROI[row][col]
            wSliceX[row] = val
            if (numtype(val) == 0)
                s  += val
                ss += val*val
                cnt += 1
            endif
        endfor

        if (cnt < 5)
            continue
        endif

        mean = s / cnt
        var  = ss / cnt - mean*mean
        if (var < 0)
            var = 0
        endif
        std = sqrt(var)

        if (std <= 1e-9*max(abs(mean), 1e-12))
            continue
        endif

        SNS__FindSymCenterScan(wSliceX, resS)
        if (resS[0] < 0)
            continue
        endif

        xSymVal = pnt2x(wSliceX, resS[0])
        wgt = std * max(resS[2], 1)

        sumWX     += wgt
        sumX      += wgt * xSymVal
        sumScoreX += wgt * resS[1]
        nUsedX    += 1
    endfor

    if (sumWX <= 0 || nUsedX < 1)
        SetDataFolder oldDFR
        Abort "SNS_FindOffsetBySymmetry: could not determine dimX symmetry center (no usable col slices)."
    endif

    xSymVal   = sumX / sumWX
    xScoreVal = sumScoreX / sumWX

    // --------------------------
    // Write outputs (2D: X and Y)
    // --------------------------
    Variable/G $(base + "_dimXSym")      = xSymVal
    Variable/G $(base + "_dimXShiftSym") = -xSymVal
    Variable/G $(base + "_dimXSymScore") = xScoreVal
    Variable/G $(base + "_dimXSymPairs") = nUsedX
    Variable/G $(base + "_dimXUseMin")   = xUseMin
    Variable/G $(base + "_dimXUseMax")   = xUseMax

    Variable/G $(base + "_dimYSym")      = ySymVal
    Variable/G $(base + "_dimYShiftSym") = -ySymVal
    Variable/G $(base + "_dimYSymScore") = yScoreVal
    Variable/G $(base + "_dimYSymPairs") = nUsedY
    Variable/G $(base + "_dimYUseMin")   = yUseMin
    Variable/G $(base + "_dimYUseMax")   = yUseMax

    SetDataFolder oldDFR
    return 0
End

//==============================================================================
// SNS_Interp2D_Duplicate
//
// Code Purpose:
//   Duplicate a 2D wave and interpolate one selected dimension by an integer factor.
//
// Physics Role:
//   None directly. Display and alignment preparation helper.
//
// Required Inputs:
//   wIn  : source 2D wave.
//
// Optional Inputs:
//   outName, interpFactor, dim, method  : output name, upsampling factor, dimension, interpolation method.
//
// Generated Outputs:
//   Creates/overwrites the interpolated output wave in the source wave folder.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools Interp2D_Duplicate.
//==============================================================================
Function SNS_Interp2D_Duplicate(wIn, [outName, interpFactor, dim, method])
    Wave wIn
    String outName
    Variable interpFactor, dim, method

    if (WaveDims(wIn) != 2)
        Abort "SNS_Interp2D_Duplicate: wIn must be a 2D wave."
    endif

    if (ParamIsDefault(interpFactor))
        interpFactor = 2
    endif
    interpFactor = max(1, round(interpFactor))

    if (ParamIsDefault(dim))
        dim = 0
    endif
    dim = round(dim)
    if (dim != 0 && dim != 1)
        Abort "SNS_Interp2D_Duplicate: dim must be 0 or 1."
    endif

    if (ParamIsDefault(method))
        method = 1
    endif
    method = (method != 0)   // 0=linear, 1=cubic

    // Folder handling
    DFREF oldDFR = GetDataFolderDFR()
    DFREF wDFR   = GetWavesDataFolderDFR(wIn)
    SetDataFolder wDFR

    String base = CleanupName(NameOfWave(wIn), 0)
    if (ParamIsDefault(outName) || strlen(outName) == 0)
        outName = base + "_interp"
    endif
    outName = CleanupName(outName, 0)

    // If factor=1, simple duplicate and return
    if (interpFactor <= 1)
        Duplicate/O wIn, $outName
        WAVE wOut0 = $outName
        Note wOut0, "Interp2D: factor=1 (no interpolation); source="+NameOfWave(wIn)
        SetDataFolder oldDFR
        return 0
    endif

    Variable n0 = DimSize(wIn, 0)
    Variable n1 = DimSize(wIn, 1)

    Variable newN0 = n0
    Variable newN1 = n1
    if (dim == 0)
        newN0 = (n0 - 1)*interpFactor + 1
    else
        newN1 = (n1 - 1)*interpFactor + 1
    endif

    // Create output wave
    Make/O/N=(newN0, newN1) $outName
    WAVE wOut = $outName

    // Copy scaling/units and adjust interpolated dim delta
    Variable off0 = DimOffset(wIn, 0)
    Variable dlt0 = DimDelta(wIn, 0)
    if (dlt0 == 0)
        dlt0 = 1
    endif
    String u0 = WaveUnits(wIn, 0)

    Variable off1 = DimOffset(wIn, 1)
    Variable dlt1 = DimDelta(wIn, 1)
    if (dlt1 == 0)
        dlt1 = 1
    endif
    String u1 = WaveUnits(wIn, 1)

    // data units
    String dU = WaveUnits(wIn, -1)

    if (dim == 0)
        SetScale/P x, off0, dlt0/interpFactor, u0, wOut
        SetScale/P y, off1, dlt1,             u1, wOut
    else
        SetScale/P x, off0, dlt0,             u0, wOut
        SetScale/P y, off1, dlt1/interpFactor, u1, wOut
    endif
    SetScale d, 0, 0, dU, wOut

    // Interpolate
    Variable i, j, seg, sub, outIdx
    Variable t, t2, t3
    Variable y0, y1, m0, m1, y

    if (dim == 0)
        // Interpolate along rows (dim0) for each column j
        Make/FREE/N=(n0) colBase, m

        for (j=0; j<n1; j+=1)

            // extract column
            for (i=0; i<n0; i+=1)
                colBase[i] = wIn[i][j]
            endfor

            // compute slopes for cubic Hermite if needed
            if (method)
                m[0]    = (colBase[1]    - colBase[0])    / dlt0
                m[n0-1] = (colBase[n0-1] - colBase[n0-2]) / dlt0
                for (i=1; i<n0-1; i+=1)
                    m[i] = (colBase[i+1] - colBase[i-1]) / (2*dlt0)
                endfor
            endif

            // fill output for this column
            for (seg=0; seg<n0-1; seg+=1)
                y0 = colBase[seg]
                y1 = colBase[seg+1]

                if (method)
                    m0 = m[seg]
                    m1 = m[seg+1]
                endif

                for (sub=0; sub<interpFactor; sub+=1)
                    t = sub / interpFactor
                    outIdx = seg*interpFactor + sub

                    if (!method)
                        y = (1-t)*y0 + t*y1
                    else
                        t2 = t*t
                        t3 = t2*t
                        y  = (2*t3 - 3*t2 + 1)*y0 + (t3 - 2*t2 + t)*(dlt0*m0) \
                           + (-2*t3 + 3*t2)*y1 + (t3 - t2)*(dlt0*m1)
                    endif

                    wOut[outIdx][j] = y
                endfor
            endfor

            wOut[newN0-1][j] = colBase[n0-1]
        endfor

    else
        // Interpolate along columns (dim1) for each row i
        Make/FREE/N=(n1) rowBase, m

        for (i=0; i<n0; i+=1)

            // extract row
            for (j=0; j<n1; j+=1)
                rowBase[j] = wIn[i][j]
            endfor

            // slopes for cubic Hermite if needed
            if (method)
                m[0]    = (rowBase[1]    - rowBase[0])    / dlt1
                m[n1-1] = (rowBase[n1-1] - rowBase[n1-2]) / dlt1
                for (j=1; j<n1-1; j+=1)
                    m[j] = (rowBase[j+1] - rowBase[j-1]) / (2*dlt1)
                endfor
            endif

            // fill output for this row
            for (seg=0; seg<n1-1; seg+=1)
                y0 = rowBase[seg]
                y1 = rowBase[seg+1]

                if (method)
                    m0 = m[seg]
                    m1 = m[seg+1]
                endif

                for (sub=0; sub<interpFactor; sub+=1)
                    t = sub / interpFactor
                    outIdx = seg*interpFactor + sub

                    if (!method)
                        y = (1-t)*y0 + t*y1
                    else
                        t2 = t*t
                        t3 = t2*t
                        y  = (2*t3 - 3*t2 + 1)*y0 + (t3 - 2*t2 + t)*(dlt1*m0) \
                           + (-2*t3 + 3*t2)*y1 + (t3 - t2)*(dlt1*m1)
                    endif

                    wOut[i][outIdx] = y
                endfor
            endfor

            wOut[i][newN1-1] = rowBase[n1-1]
        endfor
    endif

    Note wOut, "Interp2D: source="+NameOfWave(wIn)+"; dim="+num2istr(dim)+ \
               "; factor="+num2istr(interpFactor)+"; method="+SelectString(method,"linear","cubicHermite")

    SetDataFolder oldDFR
    return 0
End

// SNS_MoveDataFoldersByString
//
// Code Purpose:
//   Move data folders in the current data folder whose names match a string
//   pattern into a target data folder.
//
// Physics Role:
//   None. Generic notebook/project organization helper.
//
// Call Context:
//   Public notebook helper. Used after importing SIDAM/STM source data to move
//   loaded folders into the project data-folder tree.
//
// Required Inputs:
//   s_tomove             : Igor string-match pattern for source folder names.
//   s_destinationfolder  : target Igor data-folder path.
//
// Generated Outputs:
//   Matching data folders are moved into s_destinationfolder.
//
// Side Effects:
//   Modifies the current data-folder tree. If MoveDataFolder reports that the
//   source still exists, the source folder is killed as in the STMtools helper.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools MoveDataFoldersByString.
//==============================================================================
Function SNS_MoveDataFoldersByString(s_tomove, s_destinationfolder)
    String s_tomove, s_destinationfolder

    if (strlen(s_destinationfolder) == 0)
        Abort "SNS_MoveDataFoldersByString: destination folder path must not be empty."
    endif

    if (CmpStr(s_destinationfolder[strlen(s_destinationfolder) - 1], ":") != 0)
        s_destinationfolder += ":"
    endif

    if (!DataFolderExists(s_destinationfolder))
        Abort "SNS_MoveDataFoldersByString: destination folder does not exist: " + s_destinationfolder
    endif

    String objName
    String moveList = ""
    Variable index, nMove
    DFREF dfr = GetDataFolderDFR()

    index = 0
    do
        objName = GetIndexedObjNameDFR(dfr, 4, index)
        if (strlen(objName) == 0)
            break
        endif

        if (StringMatch(objName, s_tomove))
            moveList = AddListItem(objName, moveList, ";", Inf)
        endif
        index += 1
    while (1)

    nMove = ItemsInList(moveList)
    for (index = 0; index < nMove; index += 1)
        objName = StringFromList(index, moveList)
        MoveDataFolder/O=2/Z $objName, $s_destinationfolder
        if (V_flag > 0)
            KillDataFolder/Z $objName
        endif
    endfor
End


//==============================================================================
// SNS_DeleteWavesByString
//
// Code Purpose:
//   Delete waves in the current data folder whose names match an Igor
//   StringMatch pattern.
//
// Physics Role:
//   None. Notebook cleanup helper to remove transient processing waves.
//
// Required Inputs:
//   s_todelete  : Igor StringMatch pattern for waves to delete. If list=1,
//                 semicolon-separated list of StringMatch patterns. If keep=1,
//                 these are the only wave-name patterns to keep.
//
// Optional Inputs:
//   list  : 0/default means s_todelete is one pattern. 1 means s_todelete is
//           a semicolon-separated pattern list.
//   keep  : 0/default deletes matching waves. 1 deletes waves that do not
//           match any listed pattern.
//
// Side Effects:
//   Kills matching waves in the current data folder.
//
// Compatibility:
//   SNS-prefixed equivalent of the legacy DeleteWavesByString notebook helper.
//==============================================================================
Function SNS_DeleteWavesByString(s_todelete, [list, keep])
    String s_todelete
    Variable list
    Variable keep

    String objName
    String deleteList = ""
    String pattern
    Variable matched
    Variable index, nDelete
    Variable patIndex, nPatterns
    DFREF dfr = GetDataFolderDFR()

    if (ParamIsDefault(list))
        list = 0
    endif
    if (ParamIsDefault(keep))
        keep = 0
    endif
    nPatterns = SelectNumber(list, 1, ItemsInList(s_todelete))

    index = 0
    do
        objName = GetIndexedObjNameDFR(dfr, 1, index)
        if (strlen(objName) == 0)
            break
        endif

        matched = 0
        for (patIndex = 0; patIndex < nPatterns; patIndex += 1)
            pattern = SelectString(list, s_todelete, StringFromList(patIndex, s_todelete))
            if (StringMatch(objName, pattern))
                matched = 1
                break
            endif
        endfor
        if ((!keep && matched) || (keep && !matched))
            deleteList = AddListItem(objName, deleteList, ";", Inf)
        endif
        index += 1
    while (1)

    nDelete = ItemsInList(deleteList)
    for (index = 0; index < nDelete; index += 1)
        objName = StringFromList(index, deleteList)
        KillWaves/Z $objName
    endfor
End


//==============================================================================
// SNS_DeleteDataFoldersByString
//
// Code Purpose:
//   Delete data folders in the current data folder whose names match a string
//   pattern.
//
// Physics Role:
//   None. Generic notebook/project organization helper.
//
// Required Inputs:
//   s_todelete  : Igor string-match pattern for folders to delete.
//
// Side Effects:
//   Kills matching data folders in the current data folder.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools DeleteDataFoldersByString.
//==============================================================================
Function SNS_DeleteDataFoldersByString(s_todelete)
    String s_todelete

    String objName
    String deleteList = ""
    Variable index, nDelete
    DFREF dfr = GetDataFolderDFR()

    index = 0
    do
        objName = GetIndexedObjNameDFR(dfr, 4, index)
        if (strlen(objName) == 0)
            break
        endif

        if (StringMatch(objName, s_todelete))
            deleteList = AddListItem(objName, deleteList, ";", Inf)
        endif
        index += 1
    while (1)

    nDelete = ItemsInList(deleteList)
    for (index = 0; index < nDelete; index += 1)
        objName = StringFromList(index, deleteList)
        KillDataFolder/Z $objName
    endfor
End


//==============================================================================
// SNS_RemovePrefixFromSIDAMFolder
//
// Code Purpose:
//   Remove a prefix from a loaded SIDAM data folder and from all waves inside it.
//
// Physics Role:
//   None. Import cleanup helper for source-data folders.
//
// Call Context:
//   Public notebook helper. Used when source files include figure-specific
//   prefixes but the analysis should keep the original internal names.
//
// Required Inputs:
//   dfPath  : loaded Igor data-folder path. Final ":" is optional.
//   prefix  : prefix to remove from the folder name and matching wave names.
//
// Optional Parameters:
//   dryRun   : 0 default, perform renaming; 1 report planned renaming only.
//   verbose  : 0 default, no output; 1 print renamed/planned objects.
//
// Generated Outputs:
//   Renamed waves and, when applicable, a renamed data folder.
//
// Side Effects:
//   Modifies wave names and possibly the data-folder name. Does not overwrite
//   existing target waves or folders.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools RemovePrefixFromSIDAMFolder.
//==============================================================================
Function SNS_RemovePrefixFromSIDAMFolder(dfPath, prefix, [dryRun, verbose])
    String dfPath, prefix
    Variable dryRun, verbose

    if (ParamIsDefault(dryRun))
        dryRun = 0
    endif

    if (ParamIsDefault(verbose))
        verbose = 0
    endif

    if (strlen(prefix) == 0)
        Abort "SNS_RemovePrefixFromSIDAMFolder: prefix must not be empty."
    endif

    if (strlen(dfPath) == 0)
        Abort "SNS_RemovePrefixFromSIDAMFolder: dfPath must not be empty."
    endif

    if (CmpStr(dfPath[strlen(dfPath) - 1], ":") != 0)
        dfPath += ":"
    endif

    if (!DataFolderExists(dfPath))
        Abort "SNS_RemovePrefixFromSIDAMFolder: data folder does not exist: " + dfPath
    endif

    DFREF oldDFR = GetDataFolderDFR()
    SetDataFolder $dfPath

    String wList = WaveList(prefix + "*", ";", "")
    Variable nWaves = ItemsInList(wList)

    Variable i
    String oldName, newName

    for (i = 0; i < nWaves; i += 1)
        oldName = StringFromList(i, wList)
        newName = oldName[strlen(prefix), strlen(oldName) - 1]

        WAVE/Z wTarget = $newName
        if (WaveExists(wTarget))
            SetDataFolder oldDFR
            Abort "SNS_RemovePrefixFromSIDAMFolder: target wave already exists: " + dfPath + newName
        endif
    endfor

    for (i = 0; i < nWaves; i += 1)
        oldName = StringFromList(i, wList)
        newName = oldName[strlen(prefix), strlen(oldName) - 1]

        if (dryRun)
            if (verbose)
                Print "Would rename wave: " + dfPath + oldName + "  ->  " + newName
            endif
        else
            Rename $oldName, $newName

            if (verbose)
                Print "Renamed wave: " + oldName + "  ->  " + newName
            endif
        endif
    endfor

    SetDataFolder oldDFR

    String dfNoColon = dfPath[0, strlen(dfPath) - 2]

    Variable lastColon = strsearch(dfNoColon, ":", strlen(dfNoColon), 1)
    if (lastColon < 0)
        Abort "SNS_RemovePrefixFromSIDAMFolder: invalid data-folder path: " + dfPath
    endif

    String parentPath = dfNoColon[0, lastColon]
    String folderName = dfNoColon[lastColon + 1, strlen(dfNoColon) - 1]

    if (strsearch(folderName, prefix, 0) != 0)
        if (verbose)
            Print "Folder name does not start with prefix, not renamed: " + folderName
        endif

        return 0
    endif

    String cleanFolderName = folderName[strlen(prefix), strlen(folderName) - 1]
    String cleanFolderPath = parentPath + cleanFolderName + ":"

    if (DataFolderExists(cleanFolderPath))
        Abort "SNS_RemovePrefixFromSIDAMFolder: target data folder already exists: " + cleanFolderPath
    endif

    if (dryRun)
        if (verbose)
            Print "Would rename data folder: " + dfPath + "  ->  " + cleanFolderPath
        endif
    else
        RenameDataFolder $dfNoColon, $cleanFolderName

        if (verbose)
            Print "Renamed data folder: " + dfPath + "  ->  " + cleanFolderPath
        endif
    endif

    return 0
End


//==============================================================================
// Canonical experimental settings and panel-parameter export
//==============================================================================

Static Function SNS__ResetExperimentalSettings(settings)
    STRUCT SNS_ExperimentalSettings &settings

    settings.bias_V = NaN
    settings.current_A = NaN
    settings.vmod_V = NaN
    settings.frequency_Hz = NaN
    settings.x_m = NaN
    settings.y_m = NaN
    settings.lockInStatus = ""
    settings.biasSpectroscopyLockInRun = ""
    settings.lockInActive = 0
End


// Lock-in parameters apply when the controller left the lock-in on throughout
// acquisition, or when a bias-spectroscopy sweep enabled it independently.
Static Function SNS__LockInWasActive(lockInStatus, biasSpectroscopyLockInRun)
    String lockInStatus, biasSpectroscopyLockInRun

    return CmpStr(LowerStr(lockInStatus), "on") == 0 || CmpStr(LowerStr(biasSpectroscopyLockInRun), "true") == 0
End


// Read every canonical acquisition field currently supported by this package.
// Missing numeric fields remain NaN and missing string fields remain empty so
// image and spectroscopy folders share one API.
Function SNS_ReadExperimentalSettings(dataFolderPath, settings)
    String dataFolderPath
    STRUCT SNS_ExperimentalSettings &settings

    SNS__ResetExperimentalSettings(settings)

    if (strlen(dataFolderPath) == 0)
        Abort "SNS_ReadExperimentalSettings: dataFolderPath is empty."
    endif
    if (CmpStr(dataFolderPath[strlen(dataFolderPath) - 1], ":") != 0)
        dataFolderPath += ":"
    endif
    if (!DataFolderExists(dataFolderPath))
        Abort "SNS_ReadExperimentalSettings: data folder does not exist: " + dataFolderPath
    endif

    String settingsPath = dataFolderPath + "settings:"
    if (!DataFolderExists(settingsPath))
        Abort "SNS_ReadExperimentalSettings: settings folder does not exist: " + settingsPath
    endif

    NVAR/Z bias = $(settingsPath + "Bias:Bias__V_")
    NVAR/Z current = $(settingsPath + "Current:Current__A_")
    NVAR/Z vmod = $(settingsPath + "'Lock-in':Amplitude")
    NVAR/Z lockInFrequencyNumeric = $(settingsPath + "'Lock-in':Frequency__Hz_")
    SVAR/Z lockInFrequencyText = $(settingsPath + "'Lock-in':Frequency__Hz_")
    NVAR/Z xpos = $(settingsPath + "X__m_")
    NVAR/Z ypos = $(settingsPath + "Y__m_")
    SVAR/Z lockInStatus = $(settingsPath + "'Lock-in':Lock_in_status")
    SVAR/Z biasSpectroscopyLockInRun = $(settingsPath + "'Bias Spectroscopy':Lock_In_run")

    Variable nFound = 0
    if (NVAR_Exists(bias))
        settings.bias_V = bias
        nFound += 1
    endif
    if (NVAR_Exists(current))
        settings.current_A = current
        nFound += 1
    endif
    if (NVAR_Exists(vmod))
        settings.vmod_V = vmod
        nFound += 1
    endif
    if (NVAR_Exists(lockInFrequencyNumeric))
        settings.frequency_Hz = lockInFrequencyNumeric
        nFound += 1
    elseif (SVAR_Exists(lockInFrequencyText))
        settings.frequency_Hz = str2num(lockInFrequencyText)
        nFound += 1
    endif
    if (NVAR_Exists(xpos))
        settings.x_m = xpos
        nFound += 1
    endif
    if (NVAR_Exists(ypos))
        settings.y_m = ypos
        nFound += 1
    endif
    if (SVAR_Exists(lockInStatus))
        settings.lockInStatus = lockInStatus
        nFound += 1
    endif
    if (SVAR_Exists(biasSpectroscopyLockInRun))
        settings.biasSpectroscopyLockInRun = biasSpectroscopyLockInRun
        nFound += 1
    endif
    settings.lockInActive = SNS__LockInWasActive(settings.lockInStatus, settings.biasSpectroscopyLockInRun)

    return nFound
End


Static Function/S SNS__PanelNormalizeFolder(tableFolder)
    String tableFolder

    if (strlen(tableFolder) == 0)
        Abort "SNS panel parameters: tableFolder is empty."
    endif
    if (CmpStr(tableFolder[strlen(tableFolder) - 1], ":") != 0)
        tableFolder += ":"
    endif
    return tableFolder
End


Static Function/S SNS__PanelKnownParameters()
    return "DisplayedFieldAxis;SelectedChannel;DOSScaleReference;CalculationStartLocal;CalculationFinishLocal;Bias_mV;Current_nA;Vmod_uV;Frequency_Hz;STS_X_nm;STS_Y_nm;Trajectory_STS_X_nm;Trajectory_STS_Y_nm;Bangle_deg;Bfixed_T;BxNominal_mT;BxOffset_mT;E0Offset_meV;E0Correction_meV;GridLayerEnergy_meV;GridReferenceEnergy_meV;LineStart_X_nm;LineStart_Y_nm;LineEnd_X_nm;LineEnd_Y_nm;LinePoints;SpatialStepFactor;SparseSpatialStep;GridSimulationPoints;GridSamplingStep;CalculationElapsed_s;CalculationSecondsPerInput;CalculationValidPoints;m_eff;E_F_eV;DeltaEff_meV;h_eff_nm;IntrinsicBroadening_uEV;Temperature_K;LowBroadeningTemperature_K;ModulationBroadening_uEV;BTK_barrier;EnergyGridPoints;Bmin_T;Bmax_T;FieldGridPoints;DOSScaleFactor;LowBroadeningDOSScaleFactor;VortexX_nm;VortexY_nm;VortexFlux;ScreeningModel;"
End


Static Function SNS__PanelParameterIsText(parameterName)
    String parameterName

    return CmpStr(parameterName, "Panel") == 0 || CmpStr(parameterName, "PanelType") == 0 || CmpStr(parameterName, "Input") == 0 || CmpStr(parameterName, "DisplayedFieldAxis") == 0 || CmpStr(parameterName, "SelectedChannel") == 0 || CmpStr(parameterName, "DOSScaleReference") == 0 || CmpStr(parameterName, "CalculationStartLocal") == 0 || CmpStr(parameterName, "CalculationFinishLocal") == 0
End


Static Function SNS__PanelParameterIsKnown(parameterName)
    String parameterName

    if (CmpStr(parameterName, "Panel") == 0 || CmpStr(parameterName, "PanelType") == 0 || CmpStr(parameterName, "Input") == 0)
        return 1
    endif
    return WhichListItem(parameterName, SNS__PanelKnownParameters(), ";", 0, 0) >= 0
End


Static Function SNS__PanelFindRow(tableFolder, panel)
    String tableFolder, panel

    tableFolder = SNS__PanelNormalizeFolder(tableFolder)
    Wave/T/Z panelWave = $(tableFolder + "Panel")
    if (!WaveExists(panelWave))
        Abort "SNS panel parameters: table is not initialized: " + tableFolder
    endif

    Variable row
    for (row = 0; row < numpnts(panelWave); row += 1)
        if (CmpStr(panelWave[row], panel) == 0)
            return row
        endif
    endfor
    Abort "SNS panel parameters: panel is not present in the initialized table: " + panel
End


Static Function SNS__PanelSetValueIfPresent(tableFolder, panel, parameterName, value)
    String tableFolder, panel, parameterName
    Variable value

    tableFolder = SNS__PanelNormalizeFolder(tableFolder)
    Wave/Z parameterWave = $(tableFolder + parameterName)
    if (!WaveExists(parameterWave))
        return 0
    endif
    parameterWave[SNS__PanelFindRow(tableFolder, panel)] = value
    return 1
End


Static Function SNS__PanelSetTextIfPresent(tableFolder, panel, parameterName, value)
    String tableFolder, panel, parameterName, value

    tableFolder = SNS__PanelNormalizeFolder(tableFolder)
    Wave/T/Z parameterWave = $(tableFolder + parameterName)
    if (!WaveExists(parameterWave))
        return 0
    endif
    parameterWave[SNS__PanelFindRow(tableFolder, panel)] = value
    return 1
End


// Create a panel table whose exported columns are fixed by parameterList.
// Panel, PanelType, and Input are mandatory and are prepended automatically.
Function SNS_InitPanelParameterExport(tableFolder, panelList, parameterList)
    String tableFolder, panelList, parameterList

    tableFolder = SNS__PanelNormalizeFolder(tableFolder)
    String tableNoColon = tableFolder[0, strlen(tableFolder) - 2]
    Variable nRows = ItemsInList(panelList)
    if (nRows <= 0)
        Abort "SNS_InitPanelParameterExport: panelList is empty."
    endif

    String requested = ""
    String parameterName
    Variable i, j
    for (i = 0; i < ItemsInList(parameterList); i += 1)
        parameterName = SNS__CSVTrim(StringFromList(i, parameterList))
        if (strlen(parameterName) == 0)
            continue
        endif
        if (!SNS__PanelParameterIsKnown(parameterName))
            Abort "SNS_InitPanelParameterExport: unknown parameter: " + parameterName
        endif
        if (CmpStr(parameterName, "Panel") == 0 || CmpStr(parameterName, "PanelType") == 0 || CmpStr(parameterName, "Input") == 0)
            Abort "SNS_InitPanelParameterExport: mandatory parameter must not be repeated: " + parameterName
        endif
        if (WhichListItem(parameterName, requested, ";", 0, 0) >= 0)
            Abort "SNS_InitPanelParameterExport: duplicate parameter: " + parameterName
        endif
        requested += parameterName + ";"
    endfor

    for (i = 0; i < nRows; i += 1)
        String panel = SNS__CSVTrim(StringFromList(i, panelList))
        if (strlen(panel) == 0)
            Abort "SNS_InitPanelParameterExport: panel name is empty."
        endif
        for (j = i + 1; j < nRows; j += 1)
            if (CmpStr(panel, SNS__CSVTrim(StringFromList(j, panelList))) == 0)
                Abort "SNS_InitPanelParameterExport: duplicate panel: " + panel
            endif
        endfor
    endfor

    NewDataFolder/O $tableNoColon
    String fullList = "Panel;PanelType;Input;" + requested
    String allKnown = "Panel;PanelType;Input;" + SNS__PanelKnownParameters()
    for (i = 0; i < ItemsInList(allKnown); i += 1)
        parameterName = StringFromList(i, allKnown)
        String staleWavePath = tableFolder + parameterName
        KillWaves/Z $staleWavePath
    endfor
    for (i = 0; i < ItemsInList(fullList); i += 1)
        parameterName = StringFromList(i, fullList)
        String wavePath = tableFolder + parameterName
        if (SNS__PanelParameterIsText(parameterName))
            Make/O/T/N=(nRows) $wavePath
            Wave/T textWave = $wavePath
            textWave = ""
            Note/K textWave
            Note textWave, "SNS panel-parameter column;Parameter=" + parameterName
        else
            Make/O/D/N=(nRows) $wavePath
            Wave numericWave = $wavePath
            numericWave = NaN
            Note/K numericWave
            Note numericWave, "SNS panel-parameter column;Parameter=" + parameterName
        endif
    endfor

    Wave/T panelWave = $(tableFolder + "Panel")
    for (i = 0; i < nRows; i += 1)
        panelWave[i] = SNS__CSVTrim(StringFromList(i, panelList))
    endfor

    String orderPath = tableFolder + "SNS_ParameterOrder"
    String requestedPath = tableFolder + "SNS_RequestedParameters"
    String/G $orderPath = fullList
    String/G $requestedPath = requested
    return nRows
End


// Populate one initialized panel row. Bias/current come from the setpoint
// folder. Vmod is exported only when that acquisition's general lock-in status
// was on or its bias-spectroscopy run switch was true. STS and trajectory
// positions may come from distinct folders.
Function SNS_AddExperimentalPanelParameters(tableFolder, panel, panelType, inputName, setpointFolder, stsFolder, trajectorySTSFolder, bangleDeg, fieldAxis)
    String tableFolder, panel, panelType, inputName
    String setpointFolder, stsFolder, trajectorySTSFolder, fieldAxis
    Variable bangleDeg

    STRUCT SNS_ExperimentalSettings setpointSettings
    STRUCT SNS_ExperimentalSettings stsSettings
    STRUCT SNS_ExperimentalSettings trajectorySettings
    Variable applicableVmod = NaN
    Variable applicableLockInFrequency = NaN
    SNS__ResetExperimentalSettings(setpointSettings)
    SNS__ResetExperimentalSettings(stsSettings)
    SNS__ResetExperimentalSettings(trajectorySettings)

    if (strlen(setpointFolder) > 0)
        SNS_ReadExperimentalSettings(setpointFolder, setpointSettings)
    endif
    if (strlen(stsFolder) > 0)
        SNS_ReadExperimentalSettings(stsFolder, stsSettings)
    else
        stsSettings.x_m = setpointSettings.x_m
        stsSettings.y_m = setpointSettings.y_m
    endif
    if (strlen(trajectorySTSFolder) > 0)
        SNS_ReadExperimentalSettings(trajectorySTSFolder, trajectorySettings)
    endif

    SNS__PanelSetTextIfPresent(tableFolder, panel, "PanelType", panelType)
    SNS__PanelSetTextIfPresent(tableFolder, panel, "Input", inputName)
    SNS__PanelSetTextIfPresent(tableFolder, panel, "DisplayedFieldAxis", fieldAxis)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "Bias_mV", setpointSettings.bias_V * 1e3)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "Current_nA", setpointSettings.current_A * 1e9)
    if (setpointSettings.lockInActive)
        applicableVmod = setpointSettings.vmod_V * 1e6
        applicableLockInFrequency = setpointSettings.frequency_Hz
    endif
    SNS__PanelSetValueIfPresent(tableFolder, panel, "Vmod_uV", applicableVmod)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "Frequency_Hz", applicableLockInFrequency)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "STS_X_nm", stsSettings.x_m * 1e9)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "STS_Y_nm", stsSettings.y_m * 1e9)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "Trajectory_STS_X_nm", trajectorySettings.x_m * 1e9)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "Trajectory_STS_Y_nm", trajectorySettings.y_m * 1e9)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "Bangle_deg", bangleDeg)
    return 0
End


// Populate initialized simulation columns from the active SNS settings.
Function SNS_AddSimulationPanelParameters(tableFolder, panel, selectedChannel, dosScaleFactor)
    String tableFolder, panel, selectedChannel
    Variable dosScaleFactor

    NVAR/Z mEff = root:SNS_Settings:m_eff
    NVAR/Z EF = root:SNS_Settings:E_F
    NVAR/Z Delta = root:SNS_Settings:Delta
    NVAR/Z hEff = root:SNS_Settings:lambdaL
    NVAR/Z broadening = root:SNS_Settings:Broadening
    NVAR/Z temperature = root:SNS_Settings:T_K
    NVAR/Z modulation = root:SNS_Settings:V_mod
    NVAR/Z barrier = root:SNS_Settings:BTK_barrier
    NVAR/Z nEnergy = root:SNS_Settings:NE
    NVAR/Z bMin = root:SNS_Settings:Bmin
    NVAR/Z bMax = root:SNS_Settings:Bmax
    NVAR/Z nField = root:SNS_Settings:NB

    SNS__PanelSetTextIfPresent(tableFolder, panel, "SelectedChannel", selectedChannel)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "m_eff", NVAR_Exists(mEff) ? mEff : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "E_F_eV", NVAR_Exists(EF) ? EF : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "DeltaEff_meV", NVAR_Exists(Delta) ? Delta * 1e3 : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "h_eff_nm", NVAR_Exists(hEff) ? hEff * 1e9 : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "IntrinsicBroadening_uEV", NVAR_Exists(broadening) ? broadening * 1e6 : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "Temperature_K", NVAR_Exists(temperature) ? temperature : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "ModulationBroadening_uEV", NVAR_Exists(modulation) ? modulation * 1e6 : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "BTK_barrier", NVAR_Exists(barrier) ? barrier : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "EnergyGridPoints", NVAR_Exists(nEnergy) ? nEnergy : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "Bmin_T", NVAR_Exists(bMin) ? bMin : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "Bmax_T", NVAR_Exists(bMax) ? bMax : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "FieldGridPoints", NVAR_Exists(nField) ? nField : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "DOSScaleFactor", dosScaleFactor)
    return 0
End


// Add the fixed field and sampled line geometry retained by a line-DOS
// simulation folder. All values are stored in the units used by the notebook:
// tesla, degrees, nanometres, point count, and the dimensionless step factor.
Function SNS_AddFixedLineSimulationParameters(tableFolder, panel, simulationFolder, spatialStepFactor)
    String tableFolder, panel, simulationFolder
    Variable spatialStepFactor

    if (strlen(simulationFolder) == 0)
        Abort "SNS_AddFixedLineSimulationParameters: simulationFolder is empty."
    endif
    if (CmpStr(simulationFolder[strlen(simulationFolder) - 1], ":") != 0)
        simulationFolder += ":"
    endif
    if (!DataFolderExists(simulationFolder))
        Abort "SNS_AddFixedLineSimulationParameters: data folder does not exist: " + simulationFolder
    endif

    NVAR/Z field = $(simulationFolder + "V_B_T")
    NVAR/Z angle = $(simulationFolder + "V_Bangle_deg")
    NVAR/Z xStart = $(simulationFolder + "V_xStart")
    NVAR/Z yStart = $(simulationFolder + "V_yStart")
    NVAR/Z xEnd = $(simulationFolder + "V_xEnd")
    NVAR/Z yEnd = $(simulationFolder + "V_yEnd")
    NVAR/Z nLine = $(simulationFolder + "V_NLinePts")

    if (!NVAR_Exists(field) || !NVAR_Exists(xStart) || !NVAR_Exists(yStart) || !NVAR_Exists(xEnd) || !NVAR_Exists(yEnd) || !NVAR_Exists(nLine))
        Abort "SNS_AddFixedLineSimulationParameters: required line-simulation variables are missing from " + simulationFolder
    endif

    SNS__PanelSetValueIfPresent(tableFolder, panel, "Bangle_deg", NVAR_Exists(angle) ? angle : NaN)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "Bfixed_T", field)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "LineStart_X_nm", xStart)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "LineStart_Y_nm", yStart)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "LineEnd_X_nm", xEnd)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "LineEnd_Y_nm", yEnd)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "LinePoints", nLine)
    SNS__PanelSetValueIfPresent(tableFolder, panel, "SpatialStepFactor", spatialStepFactor)
    return 0
End


Function SNS_SetPanelParameterValue(tableFolder, panel, parameterName, value)
    String tableFolder, panel, parameterName
    Variable value

    if (SNS__PanelParameterIsText(parameterName))
        Abort "SNS_SetPanelParameterValue: parameter is textual: " + parameterName
    endif
    if (!SNS__PanelSetValueIfPresent(tableFolder, panel, parameterName, value))
        Abort "SNS_SetPanelParameterValue: parameter was not selected during initialization: " + parameterName
    endif
    return 0
End


Function SNS_SetPanelParameterText(tableFolder, panel, parameterName, value)
    String tableFolder, panel, parameterName, value

    if (!SNS__PanelParameterIsText(parameterName))
        Abort "SNS_SetPanelParameterText: parameter is not textual: " + parameterName
    endif
    if (!SNS__PanelSetTextIfPresent(tableFolder, panel, parameterName, value))
        Abort "SNS_SetPanelParameterText: parameter was not selected during initialization: " + parameterName
    endif
    return 0
End


Static Function/S SNS__PanelCSVEscape(value)
    String value

    Variable quoteValue = strsearch(value, ",", 0) >= 0 || strsearch(value, "\"", 0) >= 0 || strsearch(value, "\r", 0) >= 0 || strsearch(value, "\n", 0) >= 0
    value = ReplaceString("\"", value, "\"\"")
    if (quoteValue)
        return "\"" + value + "\""
    endif
    return value
End


// Write initialized columns in their declared order. Numeric NaN values are
// emitted as empty cells; numeric values use locale-independent decimal points.
Function SNS_ExportPanelParametersCSV(tableFolder, pathName, fileName)
    String tableFolder, pathName, fileName

    tableFolder = SNS__PanelNormalizeFolder(tableFolder)
    String orderPath = tableFolder + "SNS_ParameterOrder"
    SVAR/Z parameterOrder = $orderPath
    if (!SVAR_Exists(parameterOrder))
        Abort "SNS_ExportPanelParametersCSV: table is not initialized: " + tableFolder
    endif

    Wave/T/Z panelWave = $(tableFolder + "Panel")
    if (!WaveExists(panelWave))
        Abort "SNS_ExportPanelParametersCSV: Panel wave is missing."
    endif

    Variable nCols = ItemsInList(parameterOrder)
    Variable nRows = numpnts(panelWave)
    Variable refNum, row, col
    Open/P=$pathName refNum as fileName

    String parameterName = StringFromList(0, parameterOrder)
    String line = SNS__PanelCSVEscape(parameterName)
    for (col = 1; col < nCols; col += 1)
        parameterName = StringFromList(col, parameterOrder)
        line += "," + SNS__PanelCSVEscape(parameterName)
    endfor
    fprintf refNum, "%s\r\n", line

    String valueString
    Variable value
    for (row = 0; row < nRows; row += 1)
        line = ""
        for (col = 0; col < nCols; col += 1)
            parameterName = StringFromList(col, parameterOrder)
            if (SNS__PanelParameterIsText(parameterName))
                Wave/T textWave = $(tableFolder + parameterName)
                valueString = SNS__PanelCSVEscape(textWave[row])
            else
                Wave numericWave = $(tableFolder + parameterName)
                value = numericWave[row]
                if (numtype(value) == 0)
                    sprintf valueString, "%.15g", value
                else
                    valueString = ""
                endif
            endif
            if (col == 0)
                line = valueString
            else
                line += "," + valueString
            endif
        endfor
        fprintf refNum, "%s\r\n", line
    endfor

    Close refNum
    return 0
End


//==============================================================================
// SNS_EnsureRepositoryPath
//
// Code Purpose:
//   Return the existing SNSRepositoryPath location, or prompt the user to
//   select the upload repository root when that named path is unavailable.
//
// Physics Role:
//   None. Shared path helper for self-contained analysis notebooks.
//
// Returns:
//   The upload repository root path, including its trailing colon.
//
// Failure:
//   Aborts when the repository path is unavailable after the folder prompt.
//==============================================================================
Function/S SNS_EnsureRepositoryPath()
    String repositoryRoot = ""

    PathInfo SNSRepositoryPath
    if (!V_flag)
        NewPath/Q/O/M="Select the upload repository root" SNSRepositoryPath
        PathInfo SNSRepositoryPath
    endif

    if (!V_flag)
        Abort "SNS_EnsureRepositoryPath: no upload repository folder was selected."
    endif

    repositoryRoot = S_path
    return repositoryRoot
End


//==============================================================================
// SNS_CacheCurrentPXPPathGlobals
//
// Code Purpose:
//   Cache the current experiment's home path and file name in root globals.
//
// Physics Role:
//   None. Path helper for notebook export commands.
//
// Generated Outputs:
//   root:PXP_HomePath  : folder containing the current experiment, or "".
//   root:PXP_FileName  : current experiment file name from IgorInfo(1).
//   root:PXP_FullPath  : best-effort full path.
//
// Returns:
//   The cached home path string.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools CacheCurrentPXPPathGlobals.
//==============================================================================
Function/S SNS_CacheCurrentPXPPathGlobals()
    String homePath = ""
    String pxpName = ""
    String fullPath = ""

    PathInfo home
    if (V_flag)
        homePath = S_path
    else
        homePath = ""
    endif

    pxpName = IgorInfo(1)

    if (strlen(homePath) > 0 && strlen(pxpName) > 0)
        fullPath = homePath + pxpName
    endif

    String/G root:PXP_HomePath = homePath
    String/G root:PXP_FileName = pxpName
    String/G root:PXP_FullPath = fullPath

    return homePath
End


//==============================================================================
// SNS_GetRotAngleFromWaveNote
//
// Code Purpose:
//   Read the rotation angle stored in a wave note by SNS_DuplicateRotatedArea /
//   SNS_DuplicateCutArea style workflows.
//
// Physics Role:
//   None. Display-coordinate helper for overlaying trajectories on rotated
//   image panels.
//
// Required Inputs:
//   w : wave whose note may contain rotAngleDeg=<value>.
//
// Returns:
//   rotAngleDeg as a numeric value, or NaN if the key is missing or invalid.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools GetRotAngleFromWaveNote.
//==============================================================================
Function/D SNS_GetRotAngleFromWaveNote(w)
    Wave w

    String s = StringByKey("rotAngleDeg", note(w), "=", ";")
    if (strlen(s) == 0)
        return NaN
    endif

    Variable ang = str2num(s)
    if (numtype(ang) != 0)
        return NaN
    endif

    return ang
End


//==============================================================================
// SNS_ExportTopoPLYForBlender
//
// Code Purpose:
//   Export a clean 2D topography wave as an ASCII PLY triangle mesh for
//   Blender or other 3D software. Each pixel becomes one vertex and each grid
//   cell is split into two triangles.
//
// Inputs:
//   topo           : 2D height-field wave.
//   outFile        : fully specified output-file path.
//   xSize, ySize   : optional full lateral sizes; wave scaling is used when
//                    omitted.
//   x0, y0         : optional origins when xSize/ySize are supplied.
//   zScale         : optional exported-height multiplier; default 1.
//   centerXY       : center the mesh laterally when nonzero; default 1.
//   flipY, flipZ   : optional axis reversals; default 0.
//   quiet          : suppress the run note when nonzero; default 0.
//
// Returns:
//   The output-file path.
//
// Compatibility:
//   SNS-prefixed standalone replacement for STMtools
//   ExportTopoPLYForBlender. Geometry, scaling, and face-winding behavior are
//   preserved; the caller now supplies the repository output path directly.
//==============================================================================
Function/S SNS_ExportTopoPLYForBlender(topo, outFile, [xSize, ySize, x0, y0, zScale, centerXY, flipY, flipZ, quiet])
    Wave topo
    String outFile
    Variable xSize, ySize, x0, y0, zScale, centerXY, flipY, flipZ, quiet

    Variable nx = DimSize(topo, 0)
    Variable ny = DimSize(topo, 1)
    Variable dx, dy, xStart, yStart
    Variable xCtr, yCtr
    Variable i, j, jj
    Variable xVal, yVal, zVal
    Variable refNum
    Variable nVerts, nFaces
    Variable v0, v1, v2, v3
    Variable reverseWinding

    if (WaveDims(topo) != 2)
        Abort "SNS_ExportTopoPLYForBlender: topo must be a 2D matrix wave."
    endif
    if (nx < 2 || ny < 2)
        Abort "SNS_ExportTopoPLYForBlender: topo must be at least 2x2."
    endif

    if (ParamIsDefault(zScale))
        zScale = 1
    endif
    if (ParamIsDefault(centerXY))
        centerXY = 1
    endif
    if (ParamIsDefault(flipY))
        flipY = 0
    endif
    if (ParamIsDefault(flipZ))
        flipZ = 0
    endif
    if (ParamIsDefault(quiet))
        quiet = 0
    endif

    if (!ParamIsDefault(xSize))
        xStart = ParamIsDefault(x0) ? 0 : x0
        dx = xSize / (nx - 1)
    else
        xStart = DimOffset(topo, 0)
        dx = DimDelta(topo, 0)
        if (dx == 0)
            dx = 1
        endif
    endif

    if (!ParamIsDefault(ySize))
        yStart = ParamIsDefault(y0) ? 0 : y0
        dy = ySize / (ny - 1)
    else
        yStart = DimOffset(topo, 1)
        dy = DimDelta(topo, 1)
        if (dy == 0)
            dy = 1
        endif
    endif

    xCtr = xStart + 0.5 * (nx - 1) * dx
    yCtr = yStart + 0.5 * (ny - 1) * dy
    nVerts = nx * ny
    nFaces = 2 * (nx - 1) * (ny - 1)
    reverseWinding = mod((flipY != 0) + (flipZ != 0), 2)

    Open refNum as outFile

    fprintf refNum, "ply\n"
    fprintf refNum, "format ascii 1.0\n"
    fprintf refNum, "comment Exported from Igor Pro 9\n"
    fprintf refNum, "element vertex %d\n", nVerts
    fprintf refNum, "property float x\n"
    fprintf refNum, "property float y\n"
    fprintf refNum, "property float z\n"
    fprintf refNum, "element face %d\n", nFaces
    fprintf refNum, "property list uchar uint vertex_indices\n"
    fprintf refNum, "end_header\n"

    for (j = 0; j < ny; j += 1)
        jj = flipY ? (ny - 1 - j) : j
        yVal = yStart + jj * dy
        if (centerXY)
            yVal -= yCtr
        endif

        for (i = 0; i < nx; i += 1)
            xVal = xStart + i * dx
            if (centerXY)
                xVal -= xCtr
            endif

            zVal = topo[i][j] * zScale
            if (flipZ)
                zVal = -zVal
            endif

            fprintf refNum, "%.9g %.9g %.9g\n", xVal, yVal, zVal
        endfor
    endfor

    for (j = 0; j < ny - 1; j += 1)
        for (i = 0; i < nx - 1; i += 1)
            v0 = j * nx + i
            v1 = v0 + 1
            v3 = (j + 1) * nx + i
            v2 = v3 + 1

            if (!reverseWinding)
                fprintf refNum, "3 %d %d %d\n", v0, v1, v2
                fprintf refNum, "3 %d %d %d\n", v0, v2, v3
            else
                fprintf refNum, "3 %d %d %d\n", v0, v2, v1
                fprintf refNum, "3 %d %d %d\n", v0, v3, v2
            endif
        endfor
    endfor

    Close refNum

    if (!quiet)
        SNS_Log("SNS_ExportTopoPLYForBlender: wrote " + outFile + "; vertices=" + num2str(nVerts) + "; faces=" + num2str(nFaces))
    endif

    return outFile
End


//==============================================================================
// SNS_SavePlotAs
//
// Code Purpose:
//   Export an Igor graph window using SavePICT, optionally drawing temporary
//   white cover bands around the plot area for clean vector clipping.
//
// Physics Role:
//   None. Figure-export helper.
//
// Required Inputs:
//   None. Defaults to the front graph when graphName is omitted.
//
// Optional Inputs:
//   filename       : output file name or full path. Extension is inferred for
//                    PDF/SVG when missing.
//   savePath       : Igor path string or filesystem folder used when filename
//                    is not a full path.
//   graphName      : graph window to export.
//   datatype       : SavePICT export type; -8 PDF, -9 SVG.
//   clipPlotArea   : nonzero adds temporary cover bands for vector exports.
//   clipSideOutPt  : fallback outside thickness for left/right/bottom [pt].
//   clipTopOutPt   : outside thickness for top cover band [pt].
//   clipInPt       : overlap into plot area [pt].
//   clipBotOutPt   : outside thickness for bottom cover band [pt].
//   clipLeftOutPt  : outside thickness for left cover band [pt].
//   clipRightOutPt : outside thickness for right cover band [pt].
//
// Returns:
//   Full exported path string, or "" if the graph window was not found.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools SavePlotAs.
//==============================================================================
Function/S SNS_SavePlotAs([filename, savePath, graphName, datatype, clipPlotArea, clipSideOutPt, clipTopOutPt, clipInPt, clipBotOutPt, clipLeftOutPt, clipRightOutPt])
    String filename, savePath, graphName
    Variable datatype, clipPlotArea, clipSideOutPt, clipTopOutPt, clipInPt, clipBotOutPt, clipLeftOutPt, clipRightOutPt

    if (ParamIsDefault(filename))
        filename = ""
    endif
    if (ParamIsDefault(savePath))
        savePath = ""
    endif
    if (ParamIsDefault(graphName))
        graphName = ""
    endif
    if (ParamIsDefault(datatype))
        datatype = -8
    endif

    Variable doClip = 0
    if (!ParamIsDefault(clipPlotArea))
        doClip = (clipPlotArea != 0)
    endif

    Variable sideOut = 1.0
    if (!ParamIsDefault(clipSideOutPt))
        if (numtype(clipSideOutPt) == 0)
            sideOut = clipSideOutPt
        endif
    endif
    sideOut = max(0.2, min(sideOut, 4.0))

    Variable leftOut = sideOut
    if (!ParamIsDefault(clipLeftOutPt))
        if (numtype(clipLeftOutPt) == 0)
            leftOut = clipLeftOutPt
        endif
    endif
    leftOut = max(0.2, min(leftOut, 4.0))

    Variable rightOut = sideOut
    if (!ParamIsDefault(clipRightOutPt))
        if (numtype(clipRightOutPt) == 0)
            rightOut = clipRightOutPt
        endif
    endif
    rightOut = max(0.2, min(rightOut, 4.0))

    Variable topOut = 3.0
    if (!ParamIsDefault(clipTopOutPt))
        if (numtype(clipTopOutPt) == 0)
            topOut = clipTopOutPt
        endif
    endif
    topOut = max(0.2, min(topOut, 12.0))

    Variable botOut = sideOut
    if (!ParamIsDefault(clipBotOutPt))
        if (numtype(clipBotOutPt) == 0)
            botOut = clipBotOutPt
        endif
    endif
    botOut = max(0.2, min(botOut, 12.0))

    Variable inPt = 0.4
    if (!ParamIsDefault(clipInPt))
        if (numtype(clipInPt) == 0)
            inPt = clipInPt
        endif
    endif
    inPt = max(0.0, min(inPt, 2.0))

    String win = graphName
    if (strlen(win) == 0)
        win = WinName(0, 1, 1)
    endif

    DoWindow $win
    if (V_flag == 0)
        Print "SNS_SavePlotAs: Graph window not found: " + win
        return ""
    endif

    Variable hasPath = 0
    if (strlen(filename) > 0)
        if (strsearch(filename, "\\", 0) >= 0 || strsearch(filename, "/", 0) >= 0 || strsearch(filename, ":", 0) >= 0)
            hasPath = 1
        endif
    endif

    String folder = ""
    String leaf = ""

    if (hasPath)
        Variable i, L = strlen(filename)
        Variable lastSep = -1
        for (i = L - 1; i >= 0; i -= 1)
            if (CmpStr(filename[i], ":") == 0 || CmpStr(filename[i], "\\") == 0 || CmpStr(filename[i], "/") == 0)
                lastSep = i
                break
            endif
        endfor

        if (lastSep >= 0)
            folder = filename[0, lastSep]
            if (lastSep < L - 1)
                leaf = filename[lastSep + 1, L - 1]
            else
                leaf = ""
            endif
        else
            folder = ""
            leaf = filename
        endif

        if (strlen(leaf) == 0)
            leaf = win
        endif
    else
        leaf = filename
        if (strlen(leaf) == 0)
            leaf = win
        endif

        if (strlen(savePath) > 0)
            folder = savePath
        else
            folder = SpecialDirPath("Desktop", 0, 0, 0)
        endif
    endif

    leaf = ReplaceString(":", leaf, "_")
    leaf = ReplaceString("/", leaf, "_")
    leaf = ReplaceString("\\", leaf, "_")

    if (strsearch(leaf, ".", 0) < 0)
        if (datatype == -8)
            leaf += ".pdf"
        elseif (datatype == -9)
            leaf += ".svg"
        endif
    endif

    if (strlen(folder) == 0)
        folder = SpecialDirPath("Desktop", 0, 0, 0)
    endif

    Variable hasBack = (strsearch(folder, "\\", 0) >= 0)
    Variable hasFwd = (strsearch(folder, "/", 0) >= 0)

    if (hasBack || hasFwd)
        String sep = SelectString(hasBack, "/", "\\")
        if (CmpStr(folder[strlen(folder) - 1], sep) != 0)
            folder += sep
        endif
    else
        if (CmpStr(folder[strlen(folder) - 1], ":") != 0)
            folder += ":"
        endif
    endif

    Variable isVector = (datatype == -8 || datatype == -9)
    String axisState = ""

    if (isVector && doClip)
        axisState = SNS__SaveAxisOnTopState(win)
        SNS__SetAllAxesOnTop(win, 1)
        SNS__AddSolidCoverBandsAxRel(win, leftOut, rightOut, topOut, botOut, inPt)
        DoUpdate/W=$win
    endif

    NewPath/Q/O/C PathToFolderForSaving, folder
    SavePICT/O/P=PathToFolderForSaving/E=(datatype)/WIN=$win as leaf

    if (isVector && doClip)
        SNS__RemoveSolidCoverBandsAxRel(win)
        SNS__RestoreAxisOnTopState(win, axisState)
        DoUpdate/W=$win
    endif

    return folder + leaf
End


Static Function/S SNS__GetAxisRecreation(axisInfoStr)
    String axisInfoStr

    String rec
    rec = StringByKey("RECREATION", axisInfoStr, ":", ";")
    if (strlen(rec) == 0)
        rec = StringByKey("RECREATION", axisInfoStr, "=", ";")
    endif
    return rec
End


Static Function/S SNS__SaveAxisOnTopState(win)
    String win

    String axes = AxisList(win)
    Variable n = ItemsInList(axes, ";")
    Variable i
    String state = ""

    for (i = 0; i < n; i += 1)
        String ax = StringFromList(i, axes, ";")
        if (strlen(ax) == 0)
            continue
        endif

        String info = AxisInfo(win, ax)
        String rec = SNS__GetAxisRecreation(info)

        Variable v = NumberByKey("axisOnTop(" + ax + ")", rec, "=")
        if (numtype(v) != 0)
            v = NumberByKey("axisOnTop(" + ax + ")", rec, ":")
        endif
        if (numtype(v) != 0)
            v = 0
        endif

        state += ax + ":" + num2str(v) + ";"
    endfor

    return state
End


Static Function SNS__SetAllAxesOnTop(win, v)
    String win
    Variable v

    String axes = AxisList(win)
    Variable n = ItemsInList(axes, ";")
    Variable i

    for (i = 0; i < n; i += 1)
        String ax = StringFromList(i, axes, ";")
        if (strlen(ax) == 0)
            continue
        endif
        Execute/Q ("ModifyGraph/W=" + PossiblyQuoteName(win) + " axisOnTop(" + ax + ")=" + num2str(v))
    endfor
End


Static Function SNS__RestoreAxisOnTopState(win, state)
    String win, state

    if (strlen(state) == 0)
        return 0
    endif

    Variable i, n = ItemsInList(state, ";")
    for (i = 0; i < n; i += 1)
        String item = StringFromList(i, state, ";")
        if (strlen(item) == 0)
            continue
        endif

        String ax = StringFromList(0, item, ":")
        Variable v = str2num(StringFromList(1, item, ":"))
        if (strlen(ax) == 0 || numtype(v) != 0)
            continue
        endif

        Execute/Q ("ModifyGraph/W=" + PossiblyQuoteName(win) + " axisOnTop(" + ax + ")=" + num2str(v))
    endfor
End


Static Function SNS__AddSolidCoverBandsAxRel(win, leftOutPt, rightOutPt, topOutPt, botOutPt, inPt)
    String win
    Variable leftOutPt, rightOutPt, topOutPt, botOutPt, inPt

    SNS__RemoveSolidCoverBandsAxRel(win)

    GetWindow $win axSize
    Variable axW = V_right - V_left
    Variable axH = V_bottom - V_top
    if (axW <= 0 || axH <= 0)
        return 0
    endif

    Variable outLeftX = leftOutPt / axW
    Variable outRightX = rightOutPt / axW
    Variable outTopY = topOutPt / axH
    Variable outBotY = botOutPt / axH
    Variable inX = inPt / axW
    Variable inY = inPt / axH

    Variable r = 65535, g = 65535, b = 65535

    SetDrawLayer/W=$win ProgAxes
    SetDrawEnv/W=$win push
    SetDrawEnv/W=$win xcoord=axrel, ycoord=axrel
    SetDrawEnv/W=$win gstart, gname=SNS_ExportClip

    SetDrawEnv/W=$win fillpat=1, fillfgc=(r,g,b), linethick=0, linefgc=(r,g,b,0)
    DrawRect/W=$win (-outLeftX), (-outTopY), (inX), (1 + outBotY)

    SetDrawEnv/W=$win fillpat=1, fillfgc=(r,g,b), linethick=0, linefgc=(r,g,b,0)
    DrawRect/W=$win (1 - inX), (-outTopY), (1 + outRightX), (1 + outBotY)

    SetDrawEnv/W=$win fillpat=1, fillfgc=(r,g,b), linethick=0, linefgc=(r,g,b,0)
    DrawRect/W=$win (-outLeftX), (-outTopY), (1 + outRightX), (inY)

    SetDrawEnv/W=$win fillpat=1, fillfgc=(r,g,b), linethick=0, linefgc=(r,g,b,0)
    DrawRect/W=$win (-outLeftX), (1 - inY), (1 + outRightX), (1 + outBotY)

    SetDrawEnv/W=$win gstop
    SetDrawEnv/W=$win pop
End


Static Function SNS__RemoveSolidCoverBandsAxRel(win)
    String win

    DrawAction/W=$win/L=ProgAxes getgroup=SNS_ExportClip
    if (V_flag)
        DrawAction/W=$win/L=ProgAxes getgroup=SNS_ExportClip, delete
    endif
End


//==============================================================================
// SNS_MatrixToXYZTriplet_ExportCSV
//
// Code Purpose:
//   Convert a 2D matrix wave to an XYZ triplet wave and export it as CSV.
//
// Physics Role:
//   None. Source-data export helper.
//
// Required Inputs:
//   matrixWave : 2D input wave. Its x/y wave scaling defines exported X/Y.
//   outputName : output wave name and CSV file stem.
//   xHeader    : CSV header for X column.
//   yHeader    : CSV header for Y column.
//   zHeader    : CSV header for Z column.
//
// Optional Inputs:
//   folderPath : Igor path string to export folder. If omitted, Igor asks for
//                a destination folder.
//
// Generated Outputs:
//   Creates/overwrites outputName as a 3-column XYZ wave in the current data
//   folder and writes outputName.csv to disk.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools MatrixToXYZTriplet_ExportCSV.
//   Requires Igor optional procedure include <MatrixToXYZ>.
//==============================================================================
Function SNS_MatrixToXYZTriplet_ExportCSV(matrixWave, outputName, xHeader, yHeader, zHeader, [folderPath])
    Wave matrixWave
    String outputName, xHeader, yHeader, zHeader
    String folderPath

    if (ParamIsDefault(folderPath))
        folderPath = ""
    endif

    if (WaveDims(matrixWave) != 2)
        Abort "SNS_MatrixToXYZTriplet_ExportCSV: " + NameOfWave(matrixWave) + " is not a two-dimensional wave."
    endif

    outputName = CleanupName(outputName, 1)
    if (strlen(outputName) == 0)
        Abort "SNS_MatrixToXYZTriplet_ExportCSV: outputName is not valid."
    endif

    fMatrixToXYZTriplet(matrixWave, outputName)
    Wave triplet = $outputName

    SetDimLabel 1, 0, X, triplet
    SetDimLabel 1, 1, Y, triplet
    SetDimLabel 1, 2, Z, triplet

    Note/K triplet
    Note triplet, "Column0=" + xHeader
    Note triplet, "Column1=" + yHeader
    Note triplet, "Column2=" + zHeader
    Note triplet, "SourceMatrix=" + NameOfWave(matrixWave)

    String pathName = "SNSMatrixXYZCSVPath"

    if (strlen(folderPath) == 0)
        NewPath/O/Q/M="Choose folder for CSV export" $pathName
    else
        if (CmpStr(folderPath[strlen(folderPath) - 1], ":") != 0)
            folderPath += ":"
        endif
        NewPath/O/C/Q $pathName, folderPath
    endif

    String fileName = outputName + ".csv"

    Variable refNum
    Open/P=$pathName refNum as fileName

    fprintf refNum, "%s,%s,%s\r", xHeader, yHeader, zHeader

    Variable nRows = DimSize(triplet, 0)
    Variable i
    for (i = 0; i < nRows; i += 1)
        fprintf refNum, "%.15g,%.15g,%.15g\r", triplet[i][0], triplet[i][1], triplet[i][2]
    endfor

    Close refNum
End


//==============================================================================
// SNS_TwoWaves_ExportCSV
//
// Code Purpose:
//   Combine two 1D waves into a two-column wave and export it as CSV.
//
// Physics Role:
//   None. Source-data export helper for branch overlays and line cuts.
//
// Required Inputs:
//   xWave      : 1D x-axis/source wave.
//   yWave      : 1D y-axis/source wave.
//   outputName : output wave name and CSV file stem.
//   xHeader    : CSV header for X column.
//   yHeader    : CSV header for Y column.
//
// Optional Inputs:
//   folderPath : Igor path string to export folder. If omitted, Igor asks for
//                a destination folder.
//
// Generated Outputs:
//   Creates/overwrites outputName as a 2-column wave in the current data folder
//   and writes outputName.csv to disk.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools TwoWaves_ExportCSV.
//==============================================================================
Function SNS_TwoWaves_ExportCSV(xWave, yWave, outputName, xHeader, yHeader, [folderPath])
    Wave xWave, yWave
    String outputName, xHeader, yHeader
    String folderPath

    if (ParamIsDefault(folderPath))
        folderPath = ""
    endif

    if (WaveDims(xWave) != 1)
        Abort "SNS_TwoWaves_ExportCSV: " + NameOfWave(xWave) + " is not a one-dimensional wave."
    endif

    if (WaveDims(yWave) != 1)
        Abort "SNS_TwoWaves_ExportCSV: " + NameOfWave(yWave) + " is not a one-dimensional wave."
    endif

    Variable nRows = numpnts(xWave)

    if (numpnts(yWave) != nRows)
        Abort "SNS_TwoWaves_ExportCSV: input waves must have the same number of points."
    endif

    outputName = CleanupName(outputName, 1)
    if (strlen(outputName) == 0)
        Abort "SNS_TwoWaves_ExportCSV: outputName is not valid."
    endif

    Make/O/D/N=(nRows, 2) $outputName
    Wave out = $outputName

    out[][0] = xWave[p]
    out[][1] = yWave[p]

    SetDimLabel 1, 0, X, out
    SetDimLabel 1, 1, Y, out

    Note/K out
    Note out, "Column0=" + xHeader
    Note out, "Column1=" + yHeader
    Note out, "SourceXWave=" + NameOfWave(xWave)
    Note out, "SourceYWave=" + NameOfWave(yWave)

    String pathName = "SNSTwoWavesCSVPath"

    if (strlen(folderPath) == 0)
        NewPath/O/Q/M="Choose folder for CSV export" $pathName
    else
        if (CmpStr(folderPath[strlen(folderPath) - 1], ":") != 0)
            folderPath += ":"
        endif
        NewPath/O/C/Q $pathName, folderPath
    endif

    String fileName = outputName + ".csv"

    Variable refNum
    Open/P=$pathName refNum as fileName

    fprintf refNum, "%s,%s\r", xHeader, yHeader

    Variable i
    for (i = 0; i < nRows; i += 1)
        fprintf refNum, "%.15g,%.15g\r", out[i][0], out[i][1]
    endfor

    Close refNum
End


//==============================================================================
// SNS__CSVTrim
//
// Code Purpose:
//   Remove leading and trailing spaces/tabs/newlines from a string.
//==============================================================================
Static Function/S SNS__CSVTrim(s)
    String s

    Variable n = strlen(s)
    if (n == 0)
        return ""
    endif

    Variable i0 = 0
    Variable i1 = n - 1
    Variable c

    for (i0 = 0; i0 < n; i0 += 1)
        c = char2num(s[i0])
        if (c != 32 && c != 9 && c != 10 && c != 13)
            break
        endif
    endfor

    for (i1 = n - 1; i1 >= i0; i1 -= 1)
        c = char2num(s[i1])
        if (c != 32 && c != 9 && c != 10 && c != 13)
            break
        endif
    endfor

    if (i1 < i0)
        return ""
    endif

    return s[i0, i1]
End


//==============================================================================
// SNS_MultiWaves_ExportCSV
//
// Code Purpose:
//   Combine multiple corresponding 1D waves into a multi-column wave and export
//   it as CSV.
//
// Physics Role:
//   None. Source-data export helper for summary data, fits, and line cuts.
//
// Required Inputs:
//   waveList   : semicolon-separated list of 1D wave names or full paths.
//   headerList : semicolon-separated list of CSV column headers.
//   outputName : output wave name and CSV file stem.
//
// Optional Inputs:
//   folderPath : Igor path string to export folder. If omitted, Igor asks for
//                a destination folder.
//
// Generated Outputs:
//   Creates/overwrites outputName as an N-column wave in the current data folder
//   and writes outputName.csv to disk.
//
// Overwrite Behavior:
//   Reusing outputName replaces the existing in-memory output wave via Make/O.
//   If outputName.csv already exists in the selected export folder, Open without
//   /R or /A opens it for writing and overwrites/truncates its previous contents.
//   The function does not append, create a versioned filename, or retain a backup.
//
// Compatibility:
//   SNS-prefixed replacement for STMtools MultiWaves_ExportCSV.
//==============================================================================
Function SNS_MultiWaves_ExportCSV(waveList, headerList, outputName, [folderPath])
    String waveList, headerList, outputName
    String folderPath

    if (ParamIsDefault(folderPath))
        folderPath = ""
    endif

    Variable nCols = ItemsInList(waveList)
    if (nCols <= 0)
        Abort "SNS_MultiWaves_ExportCSV: waveList is empty."
    endif

    if (ItemsInList(headerList) != nCols)
        Abort "SNS_MultiWaves_ExportCSV: number of headers must match number of waves."
    endif

    outputName = CleanupName(outputName, 1)
    if (strlen(outputName) == 0)
        Abort "SNS_MultiWaves_ExportCSV: outputName is not valid."
    endif

    String wName = SNS__CSVTrim(StringFromList(0, waveList))
    Wave/Z w0 = $wName
    if (!WaveExists(w0))
        Abort "SNS_MultiWaves_ExportCSV: wave not found: " + wName
    endif
    if (WaveDims(w0) != 1)
        Abort "SNS_MultiWaves_ExportCSV: " + NameOfWave(w0) + " is not a one-dimensional wave."
    endif

    Variable nRows = numpnts(w0)
    Variable c
    for (c = 1; c < nCols; c += 1)
        wName = SNS__CSVTrim(StringFromList(c, waveList))
        Wave/Z w = $wName
        if (!WaveExists(w))
            Abort "SNS_MultiWaves_ExportCSV: wave not found: " + wName
        endif
        if (WaveDims(w) != 1)
            Abort "SNS_MultiWaves_ExportCSV: " + NameOfWave(w) + " is not a one-dimensional wave."
        endif
        if (numpnts(w) != nRows)
            Abort "SNS_MultiWaves_ExportCSV: all input waves must have the same number of points."
        endif
    endfor

    Make/O/D/N=(nRows, nCols) $outputName
    Wave out = $outputName

    String header
    String dimLabel
    for (c = 0; c < nCols; c += 1)
        wName = SNS__CSVTrim(StringFromList(c, waveList))
        Wave wCol = $wName
        out[][c] = wCol[p]

        dimLabel = "C" + num2istr(c)
        SetDimLabel 1, c, $dimLabel, out
    endfor

    Note/K out
    Note out, "SNS_MultiWaves_ExportCSV"
    Note out, "Headers=" + headerList
    Note out, "Sources=" + waveList

    for (c = 0; c < nCols; c += 1)
        header = SNS__CSVTrim(StringFromList(c, headerList))
        wName  = SNS__CSVTrim(StringFromList(c, waveList))

        Note out, "Column" + num2istr(c) + "=" + header
        Note out, "Source" + num2istr(c) + "=" + wName
    endfor

    String pathName = "SNSMultiWavesCSVPath"
    if (strlen(folderPath) == 0)
        NewPath/O/Q/M="Choose folder for CSV export" $pathName
    else
        if (CmpStr(folderPath[strlen(folderPath) - 1], ":") != 0)
            folderPath += ":"
        endif
        NewPath/O/C/Q $pathName, folderPath
    endif

    String fileName = outputName + ".csv"
    Variable refNum
    Open/P=$pathName refNum as fileName

    String line = SNS__CSVTrim(StringFromList(0, headerList))
    String valStr
    for (c = 1; c < nCols; c += 1)
        line += "," + SNS__CSVTrim(StringFromList(c, headerList))
    endfor
    fprintf refNum, "%s\r", line

    Variable i
    for (i = 0; i < nRows; i += 1)
        sprintf line, "%.15g", out[i][0]
        for (c = 1; c < nCols; c += 1)
            sprintf valStr, "%.15g", out[i][c]
            line += "," + valStr
        endfor
        fprintf refNum, "%s\r", line
    endfor

    Close refNum
End
