#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma DefaultTab={3,20,4}

//==============================================================================
// SNS_DisplayHelpers
//
// Standalone display helpers for SNS.
// - No STMtools dependency.
// - Uses SIDAM if available; otherwise falls back to Igor's Display/AppendImage.
//
// Naming:
//   SNS_PrepareDisplayWave
//   SNS_PrepareFiniteDisplayWave
//   SNS_DisplayWithScales  (was DisplayWithScales)
//   SNS_DetermineActualCoordinatesOfGrid
//   SNS_ApplyDavosColorToImage
//   SNS_DrawImageArrow     (was DrawImageArrow)
//   SNS_DrawVerticalArrow  (was DrawVerticalArrow)
//   SNS_DrawLineArrowOnRotatedTopo
//   SNS_AppendGridRotatedTopoTrace
//   SNS_AppendLayoutGraphByPlotOrigin
//   SNS_DrawLayoutPanelLabel (was DrawLayoutPanelLabel)
//   SNS_MakeLDOS_E_MapFromSpectra
//   SNS_FormatValueWithUncertainty
//   SNS_FormatEnergyMeVWithUncertainty
//
//==============================================================================


// Return 1 if a function exists, 0 otherwise (Igor: Exists()=6 for functions).
static Function SNS_FuncExists(fname)
    String fname
    return (Exists(fname) == 6)
End


//------------------------------------------------------------------------------
// SNS_PrepareDisplayWave
//
// Create a uniformly scaled, display-only wave for a requested x/y window.
// The source is resampled in two explicit one-dimensional stages:
//   1) Every acquired y row is linearly interpolated onto the target x grid.
//   2) Every target y row is formed using one interpolation weight applied to
//      two complete x-resampled source rows.
// Values outside the measured source-center support are NaN. The source wave
// is never modified, and scientific exports must continue to use the source.
//
// Optional arguments:
//   displayWaveName
//       Output wave name in the source wave's data folder. Default is
//       <sourceName>_display.
//   nXDisplay, nYDisplay
//       Output point counts. Defaults are inferred by rounding the requested
//       span divided by the corresponding absolute source increment.
//
// Returns:
//   The resampled display wave. The source wave is never modified.
//------------------------------------------------------------------------------
Function/WAVE SNS_PrepareDisplayWave(sourceWave, xMin, xMax, yMin, yMax, [displayWaveName, nXDisplay, nYDisplay])
    Wave sourceWave
    Variable xMin, xMax, yMin, yMax
    String displayWaveName
    Variable nXDisplay, nYDisplay

    if (WaveDims(sourceWave) != 2)
        Abort "SNS_PrepareDisplayWave: source wave must be two-dimensional."
    endif
    if ((numtype(xMin) != 0) || (numtype(xMax) != 0) || (numtype(yMin) != 0) || (numtype(yMax) != 0))
        Abort "SNS_PrepareDisplayWave: requested bounds must be finite."
    endif
    if ((xMax <= xMin) || (yMax <= yMin))
        Abort "SNS_PrepareDisplayWave: each maximum bound must exceed its minimum bound."
    endif

    Variable nx = DimSize(sourceWave, 0)
    Variable ny = DimSize(sourceWave, 1)
    Variable dx = DimDelta(sourceWave, 0)
    Variable dy = DimDelta(sourceWave, 1)
    if ((nx <= 0) || (ny <= 0) || (numtype(dx) != 0) || (numtype(dy) != 0) || (dx == 0) || (dy == 0))
        Abort "SNS_PrepareDisplayWave: source wave has invalid dimension scaling."
    endif

    Variable xFirst = DimOffset(sourceWave, 0)
    Variable xLast = xFirst + dx*(nx-1)
    Variable yFirst = DimOffset(sourceWave, 1)
    Variable yLast = yFirst + dy*(ny-1)
    Variable xCenterMin = min(xFirst, xLast)
    Variable xCenterMax = max(xFirst, xLast)
    Variable yCenterMin = min(yFirst, yLast)
    Variable yCenterMax = max(yFirst, yLast)

    if ((xMax < xCenterMin) || (xMin > xCenterMax) || (yMax < yCenterMin) || (yMin > yCenterMax))
        Abort "SNS_PrepareDisplayWave: requested window does not overlap the source wave."
    endif

    if (ParamIsDefault(nXDisplay))
        nXDisplay = round((xMax-xMin)/abs(dx))
    endif
    if (ParamIsDefault(nYDisplay))
        nYDisplay = round((yMax-yMin)/abs(dy))
    endif
    if ((numtype(nXDisplay) != 0) || (numtype(nYDisplay) != 0))
        Abort "SNS_PrepareDisplayWave: output point counts must be positive."
    endif
    nXDisplay = round(nXDisplay)
    nYDisplay = round(nYDisplay)
    if ((nXDisplay < 1) || (nYDisplay < 1))
        Abort "SNS_PrepareDisplayWave: output point counts must be positive."
    endif

    if (ParamIsDefault(displayWaveName) || (strlen(displayWaveName) == 0))
        displayWaveName = CleanupName(NameOfWave(sourceWave) + "_display", 0)
    else
        displayWaveName = CleanupName(displayWaveName, 0)
    endif
    String displayWavePath = GetWavesDataFolder(sourceWave, 1) + displayWaveName

    Make/O/D/N=(nXDisplay,nYDisplay) $displayWavePath = NaN
    Wave displayWave = $displayWavePath

    Variable dxDisplay = (xMax-xMin)/nXDisplay
    Variable dyDisplay = (yMax-yMin)/nYDisplay
    Variable xOffsetDisplay = xMin + dxDisplay/2
    Variable yOffsetDisplay = yMin + dyDisplay/2
    SetScale/P x, xOffsetDisplay, dxDisplay, WaveUnits(sourceWave, 0), displayWave
    SetScale/P y, yOffsetDisplay, dyDisplay, WaveUnits(sourceWave, 1), displayWave
    SetScale d, 0, 0, WaveUnits(sourceWave, -1), displayWave

    Make/FREE/D/N=(nXDisplay,ny) energyResampled = NaN
    Variable interpolationTol = 1e-9
    Variable iTarget, sourceIndex, sourceIndex0, alpha
    Variable xTarget, yTarget

    for (iTarget = 0; iTarget < nXDisplay; iTarget += 1)
        xTarget = xOffsetDisplay + iTarget*dxDisplay
        sourceIndex = (xTarget-xFirst)/dx
        if ((sourceIndex >= -interpolationTol) && (sourceIndex <= nx-1+interpolationTol))
            sourceIndex = max(0, min(nx-1, sourceIndex))
            sourceIndex0 = floor(sourceIndex)
            alpha = sourceIndex-sourceIndex0
            if ((sourceIndex0 >= nx-1) || (alpha <= interpolationTol))
                energyResampled[iTarget][] = sourceWave[sourceIndex0][q]
            elseif ((1-alpha) <= interpolationTol)
                energyResampled[iTarget][] = sourceWave[sourceIndex0+1][q]
            else
                energyResampled[iTarget][] = (1-alpha)*sourceWave[sourceIndex0][q] + alpha*sourceWave[sourceIndex0+1][q]
            endif
        endif
    endfor

    for (iTarget = 0; iTarget < nYDisplay; iTarget += 1)
        yTarget = yOffsetDisplay + iTarget*dyDisplay
        sourceIndex = (yTarget-yFirst)/dy
        if ((sourceIndex >= -interpolationTol) && (sourceIndex <= ny-1+interpolationTol))
            sourceIndex = max(0, min(ny-1, sourceIndex))
            sourceIndex0 = floor(sourceIndex)
            alpha = sourceIndex-sourceIndex0
            if ((sourceIndex0 >= ny-1) || (alpha <= interpolationTol))
                displayWave[][iTarget] = energyResampled[p][sourceIndex0]
            elseif ((1-alpha) <= interpolationTol)
                displayWave[][iTarget] = energyResampled[p][sourceIndex0+1]
            else
                displayWave[][iTarget] = (1-alpha)*energyResampled[p][sourceIndex0] + alpha*energyResampled[p][sourceIndex0+1]
            endif
        endif
    endfor

    String oldNote = note(sourceWave)
    if ((strlen(oldNote) > 0) && (cmpstr(oldNote[strlen(oldNote)-1, strlen(oldNote)-1], ";") != 0))
        oldNote += ";"
    endif
    String displayNote = "DisplaySource=" + GetWavesDataFolder(sourceWave, 2) + ";"
    displayNote += "DisplayRequested=(" + num2str(xMin) + "," + num2str(xMax) + "," + num2str(yMin) + "," + num2str(yMax) + ");"
    displayNote += "DisplaySourceCenters=(" + num2str(xCenterMin) + "," + num2str(xCenterMax) + "," + num2str(yCenterMin) + "," + num2str(yCenterMax) + ");"
    displayNote += "DisplaySize=(" + num2istr(nXDisplay) + "," + num2istr(nYDisplay) + ");"
    displayNote += "DisplayInterpolation=separable-linear;DisplayFieldInterpolation=row-wise;DisplayOutsideSupport=NaN;"
    Note/K displayWave, oldNote + displayNote

    return displayWave
End


//------------------------------------------------------------------------------
// SNS_PrepareFiniteDisplayWave
//
// Create a display-only crop containing the largest axis-aligned rectangle of
// finite source pixels. The source wave is never modified. This is intended
// for rotated images whose NaN perimeter should remain in scientific exports
// but should not be shown in a figure panel.
//
// Optional argument:
//   displayWaveName
//       Output wave name in the source wave's data folder. Default is
//       <sourceName>_display.
//
// Returns:
//   The cropped display wave.
//------------------------------------------------------------------------------
Function/WAVE SNS_PrepareFiniteDisplayWave(sourceWave, [displayWaveName])
    Wave sourceWave
    String displayWaveName

    if (WaveDims(sourceWave) != 2)
        Abort "SNS_PrepareFiniteDisplayWave: source wave must be two-dimensional."
    endif

    Variable nx = DimSize(sourceWave, 0)
    Variable ny = DimSize(sourceWave, 1)
    if ((nx <= 0) || (ny <= 0))
        Abort "SNS_PrepareFiniteDisplayWave: source wave is empty."
    endif

    Make/FREE/D/N=(nx) finiteHeights = 0
    Make/FREE/D/N=(nx) stackIndex = 0

    Variable ix, iy, stackTop, currentHeight
    Variable rectHeight, rectWidth, rectArea, rectLeft
    Variable bestArea = 0
    Variable bestX0 = 0
    Variable bestX1 = -1
    Variable bestY0 = 0
    Variable bestY1 = -1

    for (iy = 0; iy < ny; iy += 1)
        for (ix = 0; ix < nx; ix += 1)
            if (numtype(sourceWave[ix][iy]) == 0)
                finiteHeights[ix] += 1
            else
                finiteHeights[ix] = 0
            endif
        endfor

        stackTop = -1
        for (ix = 0; ix <= nx; ix += 1)
            if (ix < nx)
                currentHeight = finiteHeights[ix]
            else
                currentHeight = 0
            endif
            do
                if (stackTop < 0)
                    break
                endif
                if (finiteHeights[stackIndex[stackTop]] <= currentHeight)
                    break
                endif
                rectHeight = finiteHeights[stackIndex[stackTop]]
                stackTop -= 1
                if (stackTop >= 0)
                    rectLeft = stackIndex[stackTop] + 1
                else
                    rectLeft = 0
                endif
                rectWidth = ix - rectLeft
                rectArea = rectHeight * rectWidth
                if (rectArea > bestArea)
                    bestArea = rectArea
                    bestX0 = rectLeft
                    bestX1 = ix - 1
                    bestY0 = iy - rectHeight + 1
                    bestY1 = iy
                endif
            while (1)

            if (ix < nx)
                stackTop += 1
                stackIndex[stackTop] = ix
            endif
        endfor
    endfor

    if (bestArea <= 0)
        Abort "SNS_PrepareFiniteDisplayWave: source wave contains no finite pixels."
    endif

    if (ParamIsDefault(displayWaveName) || (strlen(displayWaveName) == 0))
        displayWaveName = CleanupName(NameOfWave(sourceWave) + "_display", 0)
    else
        displayWaveName = CleanupName(displayWaveName, 0)
    endif
    String displayWavePath = GetWavesDataFolder(sourceWave, 1) + displayWaveName

    Duplicate/O/R=[bestX0,bestX1][bestY0,bestY1] sourceWave, $displayWavePath
    Wave displayWave = $displayWavePath
    SetScale/P x, DimOffset(sourceWave, 0) + bestX0*DimDelta(sourceWave, 0), DimDelta(sourceWave, 0), WaveUnits(sourceWave, 0), displayWave
    SetScale/P y, DimOffset(sourceWave, 1) + bestY0*DimDelta(sourceWave, 1), DimDelta(sourceWave, 1), WaveUnits(sourceWave, 1), displayWave
    SetScale d, 0, 0, WaveUnits(sourceWave, -1), displayWave

    String oldNote = note(sourceWave)
    if ((strlen(oldNote) > 0) && (cmpstr(oldNote[strlen(oldNote)-1, strlen(oldNote)-1], ";") != 0))
        oldNote += ";"
    endif
    String displayNote = "DisplaySource=" + GetWavesDataFolder(sourceWave, 2) + ";"
    displayNote += "DisplayFiniteRectangleIndices=(" + num2istr(bestX0) + "," + num2istr(bestX1) + "," + num2istr(bestY0) + "," + num2istr(bestY1) + ");"
    displayNote += "DisplaySize=(" + num2istr(bestX1-bestX0+1) + "," + num2istr(bestY1-bestY0+1) + ");"
    displayNote += "DisplayCrop=largest-axis-aligned-finite-rectangle;"
    Note/K displayWave, oldNote + displayNote

    return displayWave
End


// Resolve the SIDAM color-table wave path for native ModifyImage.
// Bind this path to a Wave reference before calling ModifyImage; do not paste
// the raw full path into an Execute-built ctab={...} command.
static Function/S SNS_NativeCtabPath(cmap)
    String cmap

    String key = LowerStr(cmap)
    String srcPath = ""

#if Exists("SIDAMColor")
    strswitch(key)
        case "vik":
            srcPath = SIDAM_DF_CTAB + "SciColMaps:'1_Diverging':vik"
            break
        case "grayc":
            srcPath = SIDAM_DF_CTAB + "SciColMaps:'0_Sequential':grayC"
            break
        case "davos":
            srcPath = SIDAM_DF_CTAB + "SciColMaps:'0_Sequential':davos"
            break
        default:
            return "Grays"
    endswitch

    return srcPath
#else
    return "Grays"
#endif
End


// Apply the displayed z range through Igor's native image machinery.
static Function SNS_ApplyNativeImageRange(g, img, zminVal, zmaxVal, ctabName)
    String g
    Wave img
    Variable zminVal, zmaxVal
    String ctabName

    if ((strlen(g) <= 0) || (numtype(zminVal) != 0) || (numtype(zmaxVal) != 0))
        return 0
    endif
    if (zmaxVal <= zminVal)
        return 0
    endif

    String imgName = NameOfWave(img)
    String imgList = ImageNameList(g, ";")
    if ((WhichListItem(imgName, imgList) < 0) && (ItemsInList(imgList) > 0))
        imgName = StringFromList(0, imgList)
    endif
    if (strlen(imgName) <= 0)
        return 0
    endif
    if (strlen(ctabName) <= 0)
        ctabName = "Grays"
    endif

    Wave/Z ctabWave = $ctabName
    if (WaveExists(ctabWave))
        ModifyImage /W=$g $imgName ctab={zminVal,zmaxVal,ctabWave,0}
        return 1
    endif

    String cmd
    cmd = "ModifyImage /W=" + g + " " + imgName
    cmd += " ctab={" + num2str(zminVal) + "," + num2str(zmaxVal) + "," + ctabName + ",0}"
    Execute cmd
    return 1
End


//------------------------------------------------------------------------------
// SNS_DisplayWithScales
//
// Minimal image display helper.
//
// Uses SIDAMDisplay/SIDAM* if available, otherwise falls back to
// Display + AppendImage.
//
// Optional arguments:
//   cmap
//       SIDAM colormap:
//           "vik"
//           "grayc"
//           "davos"
//
//   filetype
//       0 (default) : automatic color-scale label from wave z-units
//       1           : dI/dV Grid display
//                     force label = "dI/dV (\\U)"
//		 2			 : dI/dV or LDOS along line.
//       3           : dI/dV(E,B) / DOS(E,B) display.
//                     x-axis = energy or voltage, y-axis = magnetic field.
//                     Uses SIDAM's 3-sigma upper range for sequential maps,
//                     symmetric 3-sigma range for cmap="vik", and unit-aware
//                     default E/B axis windows.
//
//   w_displaySize_pt
//       Optional 2-point wave controlling graph display size in points.
//       w_displaySize_pt[0] = graph width in points
//       w_displaySize_pt[1] = graph height in points
//
//   showScales
//       1 (default) : draw the SIDAM scale bar and color scale.
//       0           : suppress both the scale bar and color scale.
//       2           : draw only the color scale.
//       3           : draw only the SIDAM scale bar.
//
//       Example:
//           Make/O/D/N=2 w_displaySize_pt = {512,512}
//           SNS_DisplayWithScales(myImage, cmap="davos", filetype=1, w_displaySize_pt=w_displaySize_pt)
//
// Returns:
//   Graph window name.
//
// Notes:
//   - Color scale is placed at the top-left.
//   - Graphs created by this helper close without a confirmation dialog.
//   - cmap="vik" forces a symmetric zero-centered color range using a
//     3-sigma robust estimate, with an extrema fallback for flat/sparse data.
//   - If w_displaySize_pt is supplied and valid, it overrides the automatic
//     size heuristic based on input wave dimensions.
//------------------------------------------------------------------------------
Function/S SNS_DisplayWithScales(img, [cmap, filetype, w_displaySize_pt, showScales])
    Wave img
    String cmap
    Variable filetype
    Wave/Z w_displaySize_pt
    Variable showScales

    if (ParamIsDefault(filetype))
        filetype = 0
    endif
    if (ParamIsDefault(showScales))
        showScales = 1
    endif
    if (numtype(showScales) != 0)
        showScales = 1
    endif

    // Preserve the historical Boolean behavior for values other than the
    // explicitly defined color-scale-only and scale-bar-only modes.
    Variable v_showScaleBar = (showScales != 0)
    Variable v_showColorScale = (showScales != 0)
    if (showScales == 2)
        v_showScaleBar = 0
    elseif (showScales == 3)
        v_showColorScale = 0
    endif

    Variable v_useManualSize = 0
    Variable v_width_pt
    Variable v_height_pt
    Variable v_isVik = 0
    Variable v_vikAbsMax = NaN
    Variable v_autoZmax = NaN
    String v_nativeCtab = "Grays"
    if (!ParamIsDefault(cmap) && SNS_FuncExists("SIDAMColor"))
        v_nativeCtab = SNS_NativeCtabPath(cmap)
    endif

    WaveStats/Q img
    v_autoZmax = 3 * V_sdev
    if ((numtype(v_autoZmax) != 0) || (v_autoZmax <= 0))
        v_autoZmax = V_max
    endif
    if ((numtype(v_autoZmax) != 0) || (v_autoZmax <= 0))
        v_autoZmax = max(abs(V_min), abs(V_max))
    endif

    if (!ParamIsDefault(cmap) && !CmpStr(LowerStr(cmap), "vik"))
        v_isVik = 1
        v_vikAbsMax = v_autoZmax
        if ((numtype(v_vikAbsMax) != 0) || (v_vikAbsMax <= 0))
            v_vikAbsMax = max(abs(V_min), abs(V_max))
        endif
    elseif (filetype == 3)
        WaveStats/Q img
        v_autoZmax = min(V_max, max(V_avg + 3 * V_sdev, 3 * V_sdev))
        if ((numtype(v_autoZmax) != 0) || (v_autoZmax <= 0))
            v_autoZmax = V_max
        endif
    endif

    if (!ParamIsDefault(w_displaySize_pt) && WaveExists(w_displaySize_pt))
        if (numpnts(w_displaySize_pt) >= 2)
            if (w_displaySize_pt[0] > 0 && w_displaySize_pt[1] > 0)
                v_useManualSize = 1
                v_width_pt = w_displaySize_pt[0]
                v_height_pt = w_displaySize_pt[1]
            endif
        endif
    endif

    String g = ""

#if Exists("SIDAMDisplay")
    if (SNS_FuncExists("SIDAMDisplay"))
        g = SIDAMDisplay(img)
        if (strlen(g) == 0)
            Abort "SNS_DisplayWithScales: SIDAMDisplay failed."
        endif
        // SIDAM helpers can create/retarget the graph internally. The active
        // top graph is the reliable handle immediately after SIDAMDisplay.
        String gTop = WinName(0, 1, 1)
        if (strlen(gTop) > 0)
            g = gTop
        endif

        if (v_useManualSize)
            ModifyGraph /W=$g width=v_width_pt, height=v_height_pt
        elseif (!CmpStr(WaveUnits(img,0), WaveUnits(img,1), 2) && CmpStr(WaveUnits(img,-1), "eV", 2))
            ModifyGraph /W=$g width=DimSize(img,0), height=DimSize(img,1)/DimDelta(img,0)*DimDelta(img,1)
        elseif (!CmpStr(WaveUnits(img,0), WaveUnits(img,1), 2) && !CmpStr(WaveUnits(img,-1), "eV", 2))
            ModifyGraph /W=$g width=DimSize(img,0)*10, height=DimSize(img,1)*10
        else
            ModifyGraph /W=$g margin(left)=44, margin(bottom)=36, margin(top)=8, margin(right)=8
            ModifyGraph /W=$g tick=0, noLabel=0, axThick=1
            SetAxis /A /W=$g left
            SetAxis /A /W=$g bottom
            ModifyGraph /W=$g width=200, height=200
            ModifyGraph /W=$g tick=2
        endif

#if Exists("SIDAMLayerAnnotation")
        if (DimSize(img,2) > 0 && SNS_FuncExists("SIDAMLayerAnnotation"))
            SIDAMLayerAnnotation("${value}", grfName=g, imgName=NameOfWave(img), digit=1)
        endif
#endif

#if Exists("SIDAMScalebar")
        if (v_showScaleBar && SNS_FuncExists("SIDAMScalebar"))
            if (filetype == 2)
                SIDAMScalebar(grfName=g, anchor="")
            elseif (filetype == 3)
                SIDAMScalebar(grfName=g, anchor="")
            elseif (filetype != 3)
                SIDAMScalebar(grfName=g, anchor="LB", fgRGBA={0,0,0,65535}, bgRGBA={65535,65535,65535,39321}, prefix=0)
            endif
        endif
#endif

#if Exists("SIDAMColor")
        if (!ParamIsDefault(cmap) && SNS_FuncExists("SIDAMColor"))
            strswitch(LowerStr(cmap))
                case "vik":
                    SIDAMColor(grfName=g, ctable=SIDAM_DF_CTAB+"SciColMaps:'1_Diverging':vik")
                    break
                case "grayc":
                    SIDAMColor(grfName=g, ctable=SIDAM_DF_CTAB+"SciColMaps:'0_Sequential':grayC")
                    break
                case "davos":
                    SIDAMColor(grfName=g, ctable=SIDAM_DF_CTAB+"SciColMaps:'0_Sequential':davos")
                    break
            endswitch
            gTop = WinName(0, 1, 1)
            if (strlen(gTop) > 0)
                g = gTop
            endif
        endif
#endif

        if (v_isVik && (numtype(v_vikAbsMax) == 0) && (v_vikAbsMax > 0))
            SNS_ApplyNativeImageRange(g, img, -v_vikAbsMax, v_vikAbsMax, v_nativeCtab)
        elseif (filetype == 3)
            SNS_ApplyNativeImageRange(g, img, 0, v_autoZmax, v_nativeCtab)
#if Exists("SIDAMRange")
        elseif (filetype == 2 && SNS_FuncExists("SIDAMRange"))
            SIDAMRange(grfName=g, zmin=0, zmaxmode=0, zmax=NaN)
        elseif (SNS_FuncExists("SIDAMRange"))
            SIDAMRange(grfName=g, zminmode=0, zmin=NaN, zmaxmode=0, zmax=NaN)
#endif
        endif

        if (v_showColorScale)
            ColorScale /C /N=textSNS /W=$g /B=(65535,65535,65535,65535*.6) /F=0 /A=LT frame=0.50, image=$NameOfWave(img), nticks=3, tickLen=2.00, tickUnit=1
            ColorScale /C /N=textSNS /W=$g lblMargin=0
            ColorScale /C /N=textSNS /W=$g /A=LT /X=0.00 /Y=0.00 vert=1, width=5, height=100
            ColorScale /C /N=textSNS /W=$g tickThick=0.50

            String zunit = WaveUnits(img,-1)

            if (filetype == 1 || filetype == 2 || filetype == 3)
                ColorScale /C /N=textSNS /W=$g "dI/dV (\\U)"
            elseif (!CmpStr(zunit,"nA"))
                ColorScale /C /N=textSNS /W=$g "I (\\U)"
            elseif (!CmpStr(zunit,"nm") || !CmpStr(zunit,"Å") || !CmpStr(zunit,"Å") || !CmpStr(zunit,"m"))
                ColorScale /C /N=textSNS /W=$g "Δz (\\U)"
            elseif (!CmpStr(zunit,"nS") || !CmpStr(zunit,"S"))
                ColorScale /C /N=textSNS /W=$g "dI/dV (\\U)"
            elseif (!CmpStr(zunit,"eV"))
                ColorScale /C /N=textSNS /W=$g "E (\\U)"
            else
                ColorScale /C /N=textSNS /W=$g "Intensity"
            endif
        endif

        if (filetype == 2)
            Label /W=$g left "r (\\U)"
            Label /W=$g bottom "Bias (\\U)"
            ModifyGraph /W=$g tick=2
            ModifyGraph /W=$g width=300, height=600
            ModifyGraph /W=$g width=0, height=0
            ModifyGraph /W=$g tick=2
        endif

        if (filetype == 3)
            Label /W=$g left "B (\\U)"
            Label /W=$g bottom "E (\\U)"
            if (!CmpStr(WaveUnits(img,0), "mV", 2) || !CmpStr(WaveUnits(img,0), "meV", 2))
                SetAxis /W=$g bottom -0.92, 0.92
            elseif (!CmpStr(WaveUnits(img,0), "V", 2) || !CmpStr(WaveUnits(img,0), "eV", 2))
                SetAxis /W=$g bottom -0.00092, 0.00092
            endif
            if (!CmpStr(WaveUnits(img,1), "mT", 2))
                SetAxis /W=$g left -505, 505
            elseif (!CmpStr(WaveUnits(img,1), "T", 2))
                SetAxis /W=$g left -505e-3, 505e-3
            endif
            ModifyGraph /W=$g tick=2, noLabel=0
            ModifyGraph /W=$g margin(left)=44, margin(bottom)=36
        endif

        return g
    endif
#endif

    Display/K=1
    g = WinName(0, 1, 1)
    AppendImage img

    ModifyGraph /W=$g margin(left)=44, margin(bottom)=36, margin(top)=8, margin(right)=8
    SetAxis /A /W=$g left
    SetAxis /A /W=$g bottom
    if (v_isVik && (numtype(v_vikAbsMax) == 0) && (v_vikAbsMax > 0))
        SNS_ApplyNativeImageRange(g, img, -v_vikAbsMax, v_vikAbsMax, v_nativeCtab)
    elseif (filetype == 3)
        SNS_ApplyNativeImageRange(g, img, 0, v_autoZmax, v_nativeCtab)
    endif

    if (v_useManualSize)
        ModifyGraph /W=$g width=v_width_pt, height=v_height_pt
    endif

    if (v_showColorScale)
        ColorScale /C /N=textSNS /W=$g /B=(65535,65535,65535,65535*.6) /F=0 /A=LT /X=0.00 /Y=0.00 vert=1, width=5, height=100, frame=0.5, nticks=3, tickLen=2, tickThick=0.5

        if (filetype == 1 || filetype == 2 || filetype == 3)
            ColorScale /C /N=textSNS /W=$g "dI/dV (\\U)"
        else
            ColorScale /C /N=textSNS /W=$g "Intensity"
        endif
    endif

    if (filetype == 2)
        Label /W=$g left "r (\\U)"
        Label /W=$g bottom "Bias (\\U)"
        ModifyGraph /W=$g tick=2
        ModifyGraph /W=$g width=300, height=600
        ModifyGraph /W=$g width=0, height=0
        ModifyGraph /W=$g tick=2
    endif

    if (filetype == 3)
        Label /W=$g left "B (\\U)"
        Label /W=$g bottom "E (\\U)"
        if (!CmpStr(WaveUnits(img,0), "mV", 2) || !CmpStr(WaveUnits(img,0), "meV", 2))
            SetAxis /W=$g bottom -0.92, 0.92
        elseif (!CmpStr(WaveUnits(img,0), "V", 2) || !CmpStr(WaveUnits(img,0), "eV", 2))
            SetAxis /W=$g bottom -0.00092, 0.00092
        endif
        if (!CmpStr(WaveUnits(img,1), "mT", 2))
            SetAxis /W=$g left -505, 505
        elseif (!CmpStr(WaveUnits(img,1), "T", 2))
            SetAxis /W=$g left -505e-3, 505e-3
        endif
        ModifyGraph /W=$g tick=2, noLabel=0
        ModifyGraph /W=$g margin(left)=44, margin(bottom)=36
    endif

    return g
End

//==============================================================================
// SNS_DrawImageArrow
//
// Purpose:
//   Draw an arrow and label on the top graph in image-axis coordinates.
//
//   The arrow is defined by:
//      arrowTail = requested/default position
//      arrowHead = arrowTail + rotated arrow vector
//
//   The transparent background box is built as an oriented polygon around the
//   arrow and label. The DrawPoly placement anchor is treated as the lower-right
//   corner of the polygon. The anchor is shifted perpendicular to the arrow by
//   half the label offset, so the box is centered around both the arrow and label.
//
// Inputs:
//   img      : input image wave whose x/y scaling defines drawing coordinates.
//   angleDeg : arrow angle in degrees relative to +x axis.
//   labelStr : label text.
//   pos_x    : optional arrow-tail x position in axis units.
//   pos_y    : optional arrow-tail y position in axis units.
//   scale    : optional size scale for arrow, label, and box. Default: 1.
//   color    : optional color keyword:
//              "black", "white", "red", "green", "blue".
//
// Outputs:
//   return : numeric
//            0 on successful completion.
//
//==============================================================================

Function SNS_DrawImageArrow(img, angleDeg, labelStr [, pos_x, pos_y, scale, color])
    Wave img
    Variable angleDeg
    String labelStr
    Variable pos_x, pos_y, scale
    String color

    if (ParamIsDefault(scale))
        scale = 1
    endif
    if (scale <= 0 || numtype(scale) != 0)
        scale = 1
    endif

    if (ParamIsDefault(color))
        color = "black"
    endif

    String graphWindowName = WinName(0, 1, 1)
    if (strlen(graphWindowName) == 0)
        Abort "SNS_DrawImageArrow: No graph window found."
    endif

    // ------------------------------
    // Image scaling and displayed axis range
    // ------------------------------
    Variable imageOffsetX_axis = DimOffset(img, 0)
    Variable imageSpacingX_axis = DimDelta(img, 0)
    Variable imagePointCountX = DimSize(img, 0)

    Variable imageOffsetY_axis = DimOffset(img, 1)
    Variable imageSpacingY_axis = DimDelta(img, 1)
    Variable imagePointCountY = DimSize(img, 1)

    Variable imageFirstX_axis = imageOffsetX_axis
    Variable imageLastX_axis = imageOffsetX_axis + imageSpacingX_axis * (imagePointCountX - 1)

    Variable imageFirstY_axis = imageOffsetY_axis
    Variable imageLastY_axis = imageOffsetY_axis + imageSpacingY_axis * (imagePointCountY - 1)

    Variable imageMinX_axis = min(imageFirstX_axis, imageLastX_axis)
    Variable imageMaxX_axis = max(imageFirstX_axis, imageLastX_axis)

    Variable imageMinY_axis = min(imageFirstY_axis, imageLastY_axis)
    Variable imageMaxY_axis = max(imageFirstY_axis, imageLastY_axis)

    Variable imageRangeX_axis = imageMaxX_axis - imageMinX_axis
    Variable imageRangeY_axis = imageMaxY_axis - imageMinY_axis
    Variable imageShortRange_axis = min(imageRangeX_axis, imageRangeY_axis)

    if (imageRangeX_axis <= 0 || imageRangeY_axis <= 0)
        Abort "SNS_DrawImageArrow: invalid image scaling."
    endif

    // ------------------------------
    // Arrow tail
    // ------------------------------
    Variable arrowTailX_axis
    Variable arrowTailY_axis

    if (ParamIsDefault(pos_x) || ParamIsDefault(pos_y) || numtype(pos_x) != 0 || numtype(pos_y) != 0)
        Variable defaultMarginFraction = 0.10
        arrowTailX_axis = imageMaxX_axis - defaultMarginFraction * imageRangeX_axis
        arrowTailY_axis = imageMaxY_axis - defaultMarginFraction * imageRangeY_axis
    else
        arrowTailX_axis = pos_x
        arrowTailY_axis = pos_y
    endif

    // ------------------------------
    // Annotation color
    // ------------------------------
    Variable annotationRed = 0
    Variable annotationGreen = 0
    Variable annotationBlue = 0

    String colorLower = LowerStr(color)

    if (cmpstr(colorLower, "white") == 0)
        annotationRed = 65535
        annotationGreen = 65535
        annotationBlue = 65535
    elseif (cmpstr(colorLower, "red") == 0)
        annotationRed = 65535
        annotationGreen = 0
        annotationBlue = 0
    elseif (cmpstr(colorLower, "green") == 0)
        annotationRed = 0
        annotationGreen = 45000
        annotationBlue = 0
    elseif (cmpstr(colorLower, "blue") == 0)
        annotationRed = 0
        annotationGreen = 0
        annotationBlue = 65535
    else
        annotationRed = 0
        annotationGreen = 0
        annotationBlue = 0
    endif

    Variable boxFillRed = 65535
    Variable boxFillGreen = 65535
    Variable boxFillBlue = 65535
    Variable boxFillAlpha = 39321

    // ------------------------------
    // Arrow head from requested angle
    // ------------------------------
    Variable arrowAngle_rad = angleDeg * pi / 180

    Variable requestedArrowDirectionX_unit = cos(arrowAngle_rad)
    Variable requestedArrowDirectionY_unit = sin(arrowAngle_rad)

    Variable arrowLength_axis = 0.075 * imageShortRange_axis * scale

    Variable arrowHeadX_axis = arrowTailX_axis + arrowLength_axis * requestedArrowDirectionX_unit
    Variable arrowHeadY_axis = arrowTailY_axis + arrowLength_axis * requestedArrowDirectionY_unit

    // ------------------------------
    // Actual arrow vector after rotation
    // ------------------------------
    Variable arrowVectorX_axis = arrowHeadX_axis - arrowTailX_axis
    Variable arrowVectorY_axis = arrowHeadY_axis - arrowTailY_axis

    Variable arrowVectorLength_axis = sqrt(arrowVectorX_axis^2 + arrowVectorY_axis^2)
    if (arrowVectorLength_axis <= 0)
        Abort "SNS_DrawImageArrow: zero arrow length."
    endif

    Variable arrowUnitX_axis = arrowVectorX_axis / arrowVectorLength_axis
    Variable arrowUnitY_axis = arrowVectorY_axis / arrowVectorLength_axis

    Variable arrowPerpendicularUnitX_axis = -arrowUnitY_axis
    Variable arrowPerpendicularUnitY_axis =  arrowUnitX_axis

    Variable arrowCenterX_axis = 0.5 * (arrowTailX_axis + arrowHeadX_axis)
    Variable arrowCenterY_axis = 0.5 * (arrowTailY_axis + arrowHeadY_axis)

    // ------------------------------
    // Label position
    // ------------------------------
    Variable labelOffsetPerpendicular_axis = 0.28 * arrowVectorLength_axis

    Variable labelX_axis = arrowCenterX_axis + labelOffsetPerpendicular_axis * arrowPerpendicularUnitX_axis
    Variable labelY_axis = arrowCenterY_axis + labelOffsetPerpendicular_axis * arrowPerpendicularUnitY_axis

    Variable labelFontSize = 12 * scale

    // ------------------------------
    // Box geometry relative to lower-right placement anchor
    // ------------------------------
    Variable boxBehindTail_axis = 0.05 * arrowVectorLength_axis
    Variable boxPastHead_axis = 1.05 * arrowVectorLength_axis
    Variable boxTotalLengthAlongArrow_axis = boxBehindTail_axis + boxPastHead_axis

    Variable boxHalfWidthPerpendicular_axis = 0.5 * (2 * abs(labelOffsetPerpendicular_axis) + 0.28 * arrowVectorLength_axis)
    Variable boxFullWidthPerpendicular_axis = 2 * boxHalfWidthPerpendicular_axis

    // Lower-right corner of the oriented box, shifted toward label by half
    // the label offset so arrow and label are centered together.
    Variable boxCenterShiftPerpendicular_axis = 0.5 * labelOffsetPerpendicular_axis

    Variable polygonAnchorX_axis = arrowTailX_axis + boxPastHead_axis * arrowUnitX_axis + (boxCenterShiftPerpendicular_axis - boxHalfWidthPerpendicular_axis) * arrowPerpendicularUnitX_axis
    Variable polygonAnchorY_axis = arrowTailY_axis + boxPastHead_axis * arrowUnitY_axis + (boxCenterShiftPerpendicular_axis - boxHalfWidthPerpendicular_axis) * arrowPerpendicularUnitY_axis

    // ------------------------------
    // Polygon waves with local offsets from lower-right corner
    // ------------------------------
    NewDataFolder/O root:Packages
    NewDataFolder/O root:Packages:SNS_TempDraw

    NVAR/Z drawObjectCounter = root:Packages:SNS_TempDraw:v_DrawObjectCounter
    if (!NVAR_Exists(drawObjectCounter))
        Variable/G root:Packages:SNS_TempDraw:v_DrawObjectCounter = 0
        NVAR drawObjectCounter = root:Packages:SNS_TempDraw:v_DrawObjectCounter
    endif

    drawObjectCounter += 1

    String polygonXWavePath = "root:Packages:SNS_TempDraw:wArrowBoxPolygonX_" + num2istr(drawObjectCounter)
    String polygonYWavePath = "root:Packages:SNS_TempDraw:wArrowBoxPolygonY_" + num2istr(drawObjectCounter)

    Make/O/D/N=5 $polygonXWavePath
    Make/O/D/N=5 $polygonYWavePath

    Wave arrowBoxPolygonX = $polygonXWavePath
    Wave arrowBoxPolygonY = $polygonYWavePath

    // Local offsets relative to lower-right polygon anchor.
    // Corner order:
    //   0 lower-right
    //   1 lower-left
    //   2 upper-left
    //   3 upper-right
    //   4 lower-right
    arrowBoxPolygonX[0] = 0
    arrowBoxPolygonY[0] = 0

    arrowBoxPolygonX[1] = -boxTotalLengthAlongArrow_axis * arrowUnitX_axis
    arrowBoxPolygonY[1] = -boxTotalLengthAlongArrow_axis * arrowUnitY_axis

    arrowBoxPolygonX[2] = -boxTotalLengthAlongArrow_axis * arrowUnitX_axis + boxFullWidthPerpendicular_axis * arrowPerpendicularUnitX_axis
    arrowBoxPolygonY[2] = -boxTotalLengthAlongArrow_axis * arrowUnitY_axis + boxFullWidthPerpendicular_axis * arrowPerpendicularUnitY_axis

    arrowBoxPolygonX[3] = boxFullWidthPerpendicular_axis * arrowPerpendicularUnitX_axis
    arrowBoxPolygonY[3] = boxFullWidthPerpendicular_axis * arrowPerpendicularUnitY_axis

    arrowBoxPolygonX[4] = arrowBoxPolygonX[0]
    arrowBoxPolygonY[4] = arrowBoxPolygonY[0]

    // ------------------------------
    // Draw
    // ------------------------------
    SetDrawLayer/W=$graphWindowName UserFront

    // Background box around arrow and label.
    // DrawPoly origin: shifted lower-right polygon anchor.
    SetDrawEnv/W=$graphWindowName xcoord=bottom, ycoord=left, fillfgc=(boxFillRed,boxFillGreen,boxFillBlue,boxFillAlpha), linefgc=(boxFillRed,boxFillGreen,boxFillBlue,boxFillAlpha), linethick=0
    DrawPoly/W=$graphWindowName polygonAnchorX_axis, polygonAnchorY_axis, 1, 1, arrowBoxPolygonX, arrowBoxPolygonY

    // Arrow.
    SetDrawEnv/W=$graphWindowName xcoord=bottom, ycoord=left, arrow=1, arrowlen=(10*scale), arrowsharp=0, linethick=(2*scale), linefgc=(annotationRed,annotationGreen,annotationBlue)
    DrawLine/W=$graphWindowName arrowTailX_axis, arrowTailY_axis, arrowHeadX_axis, arrowHeadY_axis

    // Text.
    SetDrawEnv/W=$graphWindowName xcoord=bottom, ycoord=left, textxjust=1, textyjust=1, fsize=labelFontSize, textrgb=(annotationRed,annotationGreen,annotationBlue)
    DrawText/W=$graphWindowName labelX_axis, labelY_axis, labelStr

    return 0
End

//------------------------------------------------------------------------------
// SNS_TagSetpoint
//
// Usage:
//   SNS_TagSetpoint()                // uses top graph
//   SNS_TagSetpoint(graphname="...") // explicit graph
//------------------------------------------------------------------------------
Function SNS_TagSetpoint([graphname])
    String graphname

    // Default to top graph if none provided
    if (ParamIsDefault(graphname))
        graphname = WinName(0,1,1)
    endif

    String oldfoldername = GetDataFolder(1)

    String imagename = StringFromList(0, ImageNameList(graphname,";"), ";")
    Wave image = ImageNameToWaveRef(graphname, imagename)

    String imagefoldername = GetWavesDataFolder(image, 2)
    imagefoldername = imagefoldername[0, strsearch(imagefoldername, ":", Inf, 1)]

    STRUCT SNS_ExperimentalSettings settings
    SNS_ReadExperimentalSettings(imagefoldername, settings)
    if (numtype(settings.current_A) != 0 || numtype(settings.bias_V) != 0)
        Abort "SNS_TagSetpoint: current or bias is unavailable in " + imagefoldername
    endif

    String spText
    if (abs(settings.current_A) < 1e-10)
        spText = "\\K(0,0,0)\\JR I\\B0\\M=" + num2str(round(settings.current_A*1e12)) + \
                 " pA, V\\B0\\M= " + num2str(round(settings.bias_V*1e3)) + " mV "
    else
        spText = "\\K(0,0,0)\\JR I\\B0\\M=" + num2str(round(settings.current_A*1e10)/10) + \
                 " nA, V\\B0\\M= " + num2str(round(settings.bias_V*1e3)) + " mV "
    endif

    // SNS style: semi-transparent white background; text color set via \K escape code
    TextBox/C/N=textSetpoint0/B=(65535,65535,65535,39321)/X=0/Y=0/F=0/Z=1/A=RB spText

    SetDataFolder oldfoldername
End

//------------------------------------------------------------------------------
// SNS_RemoveBackground
//
// Coarse-to-fine background removal for 2D STM images using histogram terrace peaks.
//
// Steps:
//   1) Coarse plane subtraction (P=1) using full-image ROI.
//   2) Histogram of coarse-corrected image; detect terrace peaks (local maxima).
//      Choose lowest or highest terrace peak (useHigh).
//   3) Build ROI mask selecting pixels within ±winZ (absolute units of img z)
//      around the chosen peak (on the coarse-corrected image).
//   4) Final ImageRemoveBackground on original image using that ROI and polyOrder.
//
// ROI convention for ImageRemoveBackground:
//   background pixels MUST be set to 1; all others can be !=1 (we use 64).
//
// Optional inputs:
//   useHigh    : 0 -> lowest terrace peak; 1 -> highest terrace peak, default 0
//   polyOrder  : polynomial order for final correction, default 1
//   nbins      : histogram bins, default 256
//   winZ       : absolute half-width around terrace peak (same units as img z).
//                Default: 2 Å if units look like Å/Angstrom, else 0.2 nm if "nm",
//                else 2 (raw units).
//   outName    : output wave name, default "SNS_BGRemoved"
//
// Returns:
//   Wave reference to the output wave.
//------------------------------------------------------------------------------
Function/WAVE SNS_RemoveBackground(img, [useHigh, polyOrder, nbins, winZ, doCoarse, outName])
    Wave img
    Variable useHigh, polyOrder, nbins, winZ, doCoarse
    String outName

    // ---- defaults ----
    if (ParamIsDefault(useHigh))
        useHigh = 0
    endif
    if (ParamIsDefault(polyOrder))
        polyOrder = 1
    endif
    if (ParamIsDefault(nbins))
        nbins = 256
    endif
    if (ParamIsDefault(doCoarse))
        doCoarse = 0
    endif
    if (ParamIsDefault(outName))    	
        outName = NameOfWave(img) + "_mbgnd";
    endif
    nbins = max(32, round(nbins))

    if (DimSize(img,0) <= 1 || DimSize(img,1) <= 1)
        Abort "SNS_RemoveBackground: img must be a 2D wave."
    endif

    // default winZ from z-units if not provided
    if (ParamIsDefault(winZ))
        String zU = WaveUnits(img,-1)
        if (CmpStr(zU,"Å",2) == 0 || CmpStr(zU,"Å",2) == 0 || CmpStr(zU,"A",2) == 0 || (strsearch(LowerStr(zU),"ang",0) >= 0))
            winZ = 2          // 2 Å
        elseif (CmpStr(LowerStr(zU),"nm",2) == 0)
            winZ = 0.2        // 0.2 nm = 2 Å
        else
            winZ = 2          // fallback: user should override
        endif
    endif
    if (winZ <= 0 || numtype(winZ) != 0)
        Abort "SNS_RemoveBackground: winZ must be > 0."
    endif

    String oldDF = GetDataFolder(1)
    SetDataFolder GetWavesDataFolder(img, 1)

    Variable nx = DimSize(img,0)
    Variable ny = DimSize(img,1)

    // ---- single precision working copy of source ----
    Duplicate/O img, SNS__tmpSrc
    Redimension/S SNS__tmpSrc

    // ---- choose image for histogram/ROI building ----
    // If doCoarse=1: build ROI from coarse-corrected image
    // If doCoarse=0: build ROI from raw image (single precision copy)
    if (doCoarse != 0)

        // coarse removal: plane fit to FULL image
        Make/O/B/U/N=(nx,ny) SNS__roiAll
        SNS__roiAll = 1
        String roiAllName = NameOfWave(SNS__roiAll)

        ImageRemoveBackground /R=$roiAllName /P=1 SNS__tmpSrc
        Wave M_RemovedBackground
        Duplicate/O M_RemovedBackground, SNS__tmpCoarse

        KillWaves/Z SNS__roiAll

    else
        // no coarse step: use raw single-precision copy for histogram/ROI
        Duplicate/O SNS__tmpSrc, SNS__tmpCoarse
    endif

    // ---- histogram of ROI-source image (SNS__tmpCoarse) ----
    WaveStats/Q SNS__tmpCoarse
    Variable vmin = V_min
    Variable vmax = V_max
    if (numtype(vmin) != 0 || numtype(vmax) != 0 || vmax <= vmin)
        KillWaves/Z SNS__tmpSrc, SNS__tmpCoarse
        SetDataFolder oldDF
        Abort "SNS_RemoveBackground: invalid data range for histogram."
    endif

    Make/O/D/N=(nbins) SNS__hist
    SetScale/I x, vmin, vmax, SNS__hist
    Histogram SNS__tmpCoarse, SNS__hist

    Make/O/B/U/N=(nbins) SNS__isPeak
    SNS__isPeak = 0

    Variable i
    for (i=1; i<nbins-1; i+=1)
        if (SNS__hist[i] > SNS__hist[i-1] && SNS__hist[i] >= SNS__hist[i+1] && SNS__hist[i] > 0)
            SNS__isPeak[i] = 1
        endif
    endfor

    // pick lowest or highest terrace peak
    Variable peakIdx = -1
    if (useHigh != 0)
        for (i=nbins-2; i>=1; i-=1)
            if (SNS__isPeak[i])
                peakIdx = i
                break
            endif
        endfor
    else
        for (i=1; i<=nbins-2; i+=1)
            if (SNS__isPeak[i])
                peakIdx = i
                break
            endif
        endfor
    endif

    // fallback: global maximum bin
    if (peakIdx < 0)
        Variable maxCount = -1
        for (i=0; i<nbins; i+=1)
            if (SNS__hist[i] > maxCount)
                maxCount = SNS__hist[i]
                peakIdx = i
            endif
        endfor
    endif

    // peak center in absolute units (same as img z)
    Variable binW = (vmax - vmin) / (nbins - 1)
    Variable peakCenter = vmin + peakIdx*binW
    Variable lo = peakCenter - winZ
    Variable hi = peakCenter + winZ

    // ---- ROI based on SNS__tmpCoarse values within ±winZ ----
    Make/O/B/U/N=(nx,ny) SNS__roiBG
    SNS__roiBG = 64
    SNS__roiBG = SelectNumber((SNS__tmpCoarse >= lo) && (SNS__tmpCoarse <= hi), 64, 1)
    String roiBGName = NameOfWave(SNS__roiBG)

    // ---- final background removal on original image ----
    Duplicate/O img, SNS__tmpFinalSrc
    Redimension/S SNS__tmpFinalSrc

    ImageRemoveBackground /R=$roiBGName /P=(polyOrder) SNS__tmpFinalSrc
	wave M_RemovedBackground
    Duplicate/O M_RemovedBackground, $outName
    Wave out = $outName

    // ---- cleanup ----
    KillWaves/Z SNS__tmpSrc, SNS__tmpCoarse, SNS__tmpFinalSrc
    KillWaves/Z SNS__hist, SNS__isPeak
    // keep SNS__roiBG by default for inspection; uncomment to kill:
    //KillWaves/Z SNS__roiBG

    SetDataFolder oldDF
    return out
End

//------------------------------------------------------------------------------
// SNS_LDOSTripletToMap
//
// Convert SNS LDOS triplet output into a 2D LDOS map.
//
// Inputs:
//   w_XY
//       XY coordinate wave
//
//   w_LDOS
//       LDOS intensity wave
//
//   s_xyzName
//       Name of intermediate XYZ wave
//
//   s_mapName
//       Name of output LDOS map
//
//   w_pos_nm
//       Position wave used to infer square output dimensions
//
// Optional:
//   v_tol
//       XYZTripletToMatrixEx tolerance (default 1e-05)
//
// Notes:
//   Interpolation is performed directly with ImageInterpolate using the same
//   bounds and Voronoi settings as WaveMetrics' XYZTripletToMatrixEx macro.
//------------------------------------------------------------------------------
Function SNS_LDOSTripletToMap(w_XY, w_LDOS, s_xyzName, s_mapName, w_pos_nm, [v_tol])
    Wave w_XY
    Wave w_LDOS
    Wave w_pos_nm
    String s_xyzName
    String s_mapName
    Variable v_tol

    if (ParamIsDefault(v_tol))
        v_tol = 1e-5
    endif

    Variable v_nGrid = sqrt(DimSize(w_pos_nm, 0))

    Concatenate/O {w_XY, w_LDOS}, $s_xyzName
    Wave w_xyz = $s_xyzName
    SNS_XYZTripletToMatrix(w_xyz, s_mapName, v_nGrid, v_nGrid, v_tol)

    SetScale/P x, DimOffset($s_mapName,0), DimDelta($s_mapName,0), "nm", $s_mapName
    SetScale/P y, DimOffset($s_mapName,1), DimDelta($s_mapName,1), "nm", $s_mapName
End

// SNS_XYZTripletToMatrix
//
// Function-safe equivalent of WaveMetrics' XYZTripletToMatrixEx macro.
// Calling ImageInterpolate directly avoids Execute parsing of macro string
// parameters while preserving the macro's bounds and Voronoi interpolation.
//------------------------------------------------------------------------------
Function SNS_XYZTripletToMatrix(wTriplet, mapName, nx, ny, tol)
    Wave wTriplet
    String mapName
    Variable nx, ny, tol

    nx = round(nx)
    ny = round(ny)
    if (WaveDims(wTriplet) != 2 || DimSize(wTriplet,1) < 3)
        Abort "SNS_XYZTripletToMatrix: input must be a two-dimensional XYZ triplet wave."
    endif
    if (nx < 2 || ny < 2)
        Abort "SNS_XYZTripletToMatrix: nx and ny must each be at least 2."
    endif

    Variable nTriplet = DimSize(wTriplet,0)
    ImageStats/M=1/G={0,nTriplet-1,0,0} wTriplet
    Variable xmin = V_min
    Variable xmax = V_max
    ImageStats/M=1/G={0,nTriplet-1,1,1} wTriplet
    Variable ymin = V_min
    Variable ymax = V_max

    Variable dx = (xmax-xmin)/(nx-1)
    Variable dy = (ymax-ymin)/(ny-1)
    if (numtype(dx) != 0 || numtype(dy) != 0 || dx <= 0 || dy <= 0)
        Abort "SNS_XYZTripletToMatrix: invalid XYZ coordinate range."
    endif

    // Match WaveMetrics XYZTripletToMatrixEx: extend the upper bounds so
    // ImageInterpolate produces exactly nx by ny output points.
    xmax += dx + dx/2
    ymax += dy + dy/2

    KillWaves/Z $mapName
    ImageInterpolate/DEST=$mapName/PFTL=(tol)/E=(NaN)/S={(xmin),(dx),(xmax),(ymin),(dy),(ymax)} Voronoi, wTriplet
End


//==============================================================================
// SNS_RotateXYCoordinates
//
// Rotate coordinate pairs without resampling their associated data. The input
// may be [point][x/y] or [x/y][point]. If xCenter and yCenter are omitted, the
// center of the finite coordinate bounding box is used.
//==============================================================================
Function/WAVE SNS_RotateXYCoordinates(coords, angleDeg, [outName, xCenter, yCenter])
    Wave coords
    Variable angleDeg, xCenter, yCenter
    String outName

    if (WaveDims(coords) != 2)
        Abort "SNS_RotateXYCoordinates: coordinates must be a two-dimensional wave."
    endif

    Variable pointsByRow
    Variable nPoints
    if (DimSize(coords, 1) == 2)
        pointsByRow = 1
        nPoints = DimSize(coords, 0)
    elseif (DimSize(coords, 0) == 2)
        pointsByRow = 0
        nPoints = DimSize(coords, 1)
    else
        Abort "SNS_RotateXYCoordinates: expected [point][x/y] or [x/y][point] coordinates."
    endif

    Variable centerXDefault = ParamIsDefault(xCenter)
    Variable centerYDefault = ParamIsDefault(yCenter)
    if (centerXDefault != centerYDefault)
        Abort "SNS_RotateXYCoordinates: xCenter and yCenter must be supplied together."
    endif

    Variable i, xValue, yValue, nFinite = 0
    Variable xMin, xMax, yMin, yMax
    for (i = 0; i < nPoints; i += 1)
        if (pointsByRow)
            xValue = coords[i][0]
            yValue = coords[i][1]
        else
            xValue = coords[0][i]
            yValue = coords[1][i]
        endif
        if ((numtype(xValue) != 0) || (numtype(yValue) != 0))
            continue
        endif
        if (nFinite == 0)
            xMin = xValue
            xMax = xValue
            yMin = yValue
            yMax = yValue
        else
            xMin = min(xMin, xValue)
            xMax = max(xMax, xValue)
            yMin = min(yMin, yValue)
            yMax = max(yMax, yValue)
        endif
        nFinite += 1
    endfor
    if (nFinite <= 0)
        Abort "SNS_RotateXYCoordinates: no finite coordinate pairs."
    endif

    if (centerXDefault)
        xCenter = 0.5*(xMin + xMax)
        yCenter = 0.5*(yMin + yMax)
    endif

    if (ParamIsDefault(outName) || (strlen(outName) == 0))
        outName = CleanupName(NameOfWave(coords) + "_rot", 0)
    else
        outName = CleanupName(outName, 0)
    endif

    DFREF oldDFR = GetDataFolderDFR()
    DFREF coordDFR = GetWavesDataFolderDFR(coords)
    SetDataFolder coordDFR
    Duplicate/O coords, $outName
    Wave rotated = $outName

    Variable angleRad = angleDeg*pi/180
    Variable c = cos(angleRad)
    Variable s = sin(angleRad)
    for (i = 0; i < nPoints; i += 1)
        if (pointsByRow)
            xValue = coords[i][0]
            yValue = coords[i][1]
        else
            xValue = coords[0][i]
            yValue = coords[1][i]
        endif
        if ((numtype(xValue) != 0) || (numtype(yValue) != 0))
            if (pointsByRow)
                rotated[i][0] = NaN
                rotated[i][1] = NaN
            else
                rotated[0][i] = NaN
                rotated[1][i] = NaN
            endif
            continue
        endif
        if (pointsByRow)
            rotated[i][0] = xCenter + (xValue-xCenter)*c - (yValue-yCenter)*s
            rotated[i][1] = yCenter + (xValue-xCenter)*s + (yValue-yCenter)*c
        else
            rotated[0][i] = xCenter + (xValue-xCenter)*c - (yValue-yCenter)*s
            rotated[1][i] = yCenter + (xValue-xCenter)*s + (yValue-yCenter)*c
        endif
    endfor

    Note/K rotated
    Note rotated, "SNS_CoordinateTransform=rotation;SNS_RotationAngle_deg=" + num2str(angleDeg) + ";SNS_RotationCenterX_nm=" + num2str(xCenter) + ";SNS_RotationCenterY_nm=" + num2str(yCenter) + ";"
    SetDataFolder oldDFR
    return rotated
End


static Function SNS__BuildGridXYZTriplet(img, coords, tripName, layer)
    Wave img, coords
    String tripName
    Variable layer

    Variable nx = DimSize(img, 0)
    Variable ny = DimSize(img, 1)
    Variable nTotal = nx * ny
    if (WaveDims(coords) != 2 || DimSize(coords, 0) != 2 || DimSize(coords, 1) != nTotal)
        Abort "SNS_DetermineActualCoordinatesOfGrid: coordinates must be 2 by nx*ny."
    endif

    Make/O/D/N=(nTotal, 3) $tripName
    Wave triplet = $tripName
    Variable ip, iq, n = 0
    for (iq = 0; iq < ny; iq += 1)
        for (ip = 0; ip < nx; ip += 1)
            triplet[n][0] = coords[0][n]
            triplet[n][1] = coords[1][n]
            triplet[n][2] = WaveDims(img) == 2 ? img[ip][iq] : img[ip][iq][layer]
            n += 1
        endfor
    endfor
End


// Reconstruct a regular absolute-XY GridSTS wave from retained coordinates.
// The 2D or 3D result is named <source>_xy.
Function SNS_DetermineActualCoordinatesOfGrid(img, coords)
    Wave img, coords

    Variable dims = WaveDims(img)
    if (dims != 2 && dims != 3)
        Abort "SNS_DetermineActualCoordinatesOfGrid: source must be 2D or 3D."
    endif

    Variable nx = DimSize(img, 0)
    Variable ny = DimSize(img, 1)
    Variable nz = dims == 3 ? DimSize(img, 2) : 1
    String destName = NameOfWave(img) + "_xy"
    String tripName = "SNS_grid_xyztriplet"
    String sliceName = NameOfWave(img) + "_xy_tmp"
    Variable layer

    KillWaves/Z $destName
    if (dims == 2)
        SNS__BuildGridXYZTriplet(img, coords, tripName, 0)
        Wave triplet2D = $tripName
        SNS_XYZTripletToMatrix(triplet2D, destName, nx, ny, 1e-5)
        KillWaves/Z $tripName
        return 0
    endif

    Make/O/D/N=(nx, ny, nz) $destName
    Wave dest = $destName
    SetScale/P z, DimOffset(img, 2), DimDelta(img, 2), WaveUnits(img, 2), dest
    SetScale d, 0, 0, WaveUnits(img, -1), dest
    for (layer = 0; layer < nz; layer += 1)
        SNS__BuildGridXYZTriplet(img, coords, tripName, layer)
        KillWaves/Z $sliceName
        Wave tripletLayer = $tripName
        SNS_XYZTripletToMatrix(tripletLayer, sliceName, nx, ny, 1e-5)
        Wave slice = $sliceName
        if (layer == 0)
            SetScale/P x, DimOffset(slice, 0), DimDelta(slice, 0), "nm", dest
            SetScale/P y, DimOffset(slice, 1), DimDelta(slice, 1), "nm", dest
        endif
        dest[][][layer] = slice[p][q]
        KillWaves/Z $tripName, $sliceName
    endfor
    return 0
End


Function SNS_ApplyDavosColorToImage(graphName, imageName, [rev])
    String graphName, imageName
    Variable rev

    if (ParamIsDefault(rev))
        rev = 0
    endif
    if (WinType(graphName) != 1)
        Abort "SNS_ApplyDavosColorToImage: graph not found: " + graphName
    endif
#if Exists("SIDAMColor")
    if (rev)
        SIDAMColor(grfName=graphName, ctable=SIDAM_DF_CTAB+"SciColMaps:'0_Sequential':davos", rev=1)
    else
        SIDAMColor(grfName=graphName, ctable=SIDAM_DF_CTAB+"SciColMaps:'0_Sequential':davos")
    endif
#else
    ModifyImage/W=$graphName $imageName ctab={*,*,Grays,rev}
#endif
    return 0
End

//==============================================================================
// SNS_MakeLDOS_E_MapFromSpectra
//
// Convert one [energy][point] LDOS wave from
// SNS_LDOSmap_BFixed_FromMask into a [x][y][energy] map stack. A Voronoi map
// from output pixels to input-point indices is constructed once and reused for
// every energy slice.
//
// Inputs:
//   xyWaveName       : XY_* wave, [point][x/y] in nm.
//   ldosWaveName     : LDOS_E_raw_* or LDOS_E_conv_* wave,
//                      [energy][point].
//   validWaveName    : corresponding ValidPt_* wave.
//
// Optional inputs:
//   outName          : output stack name. Default is <ldosWaveName>_map.
//   xyzTol           : Voronoi interpolation tolerance. Default is 1e-05.
//   nx, ny           : output grid dimensions. If omitted, infer a square
//                      grid from the total number of input positions.
//
// Returns:
//   The [x][y][energy] map stack. Pixels outside the spatial support of a
//   rotated or irregular point set and pixels mapped to invalid simulation
//   positions remain NaN.
//==============================================================================
Function/WAVE SNS_MakeLDOS_E_MapFromSpectra(xyWaveName, ldosWaveName, validWaveName, [outName, xyzTol, nx, ny])
    String xyWaveName, ldosWaveName, validWaveName
    String outName
    Variable xyzTol, nx, ny

    if (ParamIsDefault(xyzTol))
        xyzTol = 1e-5
    endif

    Wave/Z xyW = $xyWaveName
    Wave/Z ldosW = $ldosWaveName
    Wave/Z validW = $validWaveName
    if (!WaveExists(xyW))
        Abort "SNS_MakeLDOS_E_MapFromSpectra: missing XY wave: " + xyWaveName
    endif
    if (!WaveExists(ldosW))
        Abort "SNS_MakeLDOS_E_MapFromSpectra: missing LDOS wave: " + ldosWaveName
    endif
    if (!WaveExists(validW))
        Abort "SNS_MakeLDOS_E_MapFromSpectra: missing ValidPt wave: " + validWaveName
    endif

    Variable nE = DimSize(ldosW, 0)
    Variable nPts = DimSize(ldosW, 1)
    if ((WaveDims(ldosW) != 2) || (nE <= 0) || (nPts <= 0))
        Abort "SNS_MakeLDOS_E_MapFromSpectra: LDOS wave must be [energy][point]."
    endif
    if ((WaveDims(xyW) != 2) || (DimSize(xyW, 0) != nPts) || (DimSize(xyW, 1) < 2))
        Abort "SNS_MakeLDOS_E_MapFromSpectra: XY wave must be [point][x/y] and match the LDOS point count."
    endif
    if ((WaveDims(validW) != 1) || (DimSize(validW, 0) != nPts))
        Abort "SNS_MakeLDOS_E_MapFromSpectra: ValidPt wave must match the LDOS point count."
    endif

    if (ParamIsDefault(nx) || ParamIsDefault(ny) || (nx <= 0) || (ny <= 0))
        Variable nGrid = round(sqrt(nPts))
        if (nGrid*nGrid != nPts)
            Abort "SNS_MakeLDOS_E_MapFromSpectra: input positions do not form a square grid. Supply nx and ny."
        endif
        nx = nGrid
        ny = nGrid
    endif
    nx = round(nx)
    ny = round(ny)
    if (nx*ny != nPts)
        Abort "SNS_MakeLDOS_E_MapFromSpectra: nx*ny does not match the input position count."
    endif

    Variable ip, nValid = 0, nFinitePositions = 0
    for (ip = 0; ip < nPts; ip += 1)
        if ((numtype(xyW[ip][0]) == 0) && (numtype(xyW[ip][1]) == 0))
            nFinitePositions += 1
            if (validW[ip] != 0)
                nValid += 1
            endif
        endif
    endfor
    if (nValid <= 0)
        Abort "SNS_MakeLDOS_E_MapFromSpectra: no valid positions."
    endif

    Make/O/D/N=(nFinitePositions,3) SNS_tmp_LDOS_point_index_xyz
    Wave pointIndexXYZ = SNS_tmp_LDOS_point_index_xyz
    Variable j = 0
    for (ip = 0; ip < nPts; ip += 1)
        if ((numtype(xyW[ip][0]) != 0) || (numtype(xyW[ip][1]) != 0))
            continue
        endif
        pointIndexXYZ[j][0] = xyW[ip][0]
        pointIndexXYZ[j][1] = xyW[ip][1]
        pointIndexXYZ[j][2] = ip
        j += 1
    endfor

    String pointIndexMapName = "SNS_tmp_LDOS_point_index_map"
    SNS_XYZTripletToMatrix(pointIndexXYZ, pointIndexMapName, nx, ny, xyzTol)
    Wave pointIndexMap = $pointIndexMapName

    Variable ix, iy, mappedPoint, nMappedPixels = 0
    for (iy = 0; iy < DimSize(pointIndexMap, 1); iy += 1)
        for (ix = 0; ix < DimSize(pointIndexMap, 0); ix += 1)
            if (numtype(pointIndexMap[ix][iy]) != 0)
                continue
            endif
            mappedPoint = round(pointIndexMap[ix][iy])
            if ((mappedPoint < 0) || (mappedPoint >= nPts))
                Abort "SNS_MakeLDOS_E_MapFromSpectra: point-index interpolation produced an invalid source index."
            endif
            pointIndexMap[ix][iy] = mappedPoint
            nMappedPixels += 1
        endfor
    endfor
    if (nMappedPixels <= 0)
        Abort "SNS_MakeLDOS_E_MapFromSpectra: point-index interpolation produced no finite pixels."
    endif

    if (ParamIsDefault(outName) || (strlen(outName) == 0))
        outName = CleanupName(NameOfWave(ldosW) + "_map", 0)
    else
        outName = CleanupName(outName, 0)
    endif

    Make/O/D/N=(DimSize(pointIndexMap,0),DimSize(pointIndexMap,1),nE) $outName
    Wave mapW = $outName
    SetScale/P x, DimOffset(pointIndexMap,0), DimDelta(pointIndexMap,0), "nm", mapW
    SetScale/P y, DimOffset(pointIndexMap,1), DimDelta(pointIndexMap,1), "nm", mapW
    SetScale/P z, DimOffset(ldosW,0), DimDelta(ldosW,0), WaveUnits(ldosW,0), mapW
    SetScale d, 0, 0, WaveUnits(ldosW,-1), mapW
    mapW = NaN

    Variable iE
    for (iE = 0; iE < nE; iE += 1)
        for (iy = 0; iy < DimSize(pointIndexMap, 1); iy += 1)
            for (ix = 0; ix < DimSize(pointIndexMap, 0); ix += 1)
                if (numtype(pointIndexMap[ix][iy]) != 0)
                    continue
                endif
                mappedPoint = round(pointIndexMap[ix][iy])
                if (validW[mappedPoint] == 0)
                    continue
                endif
                mapW[ix][iy][iE] = ldosW[iE][mappedPoint]
            endfor
        endfor
    endfor

    String mapNote = note(ldosW)
    if ((strlen(mapNote) > 0) && (cmpstr(mapNote[strlen(mapNote)-1, strlen(mapNote)-1], ";") != 0))
        mapNote += ";"
    endif
    mapNote += "SNS_MapSource=" + GetWavesDataFolder(ldosW, 2) + ";SNS_MapDims=x,y,energy;SNS_MapInterpolation=Voronoi-point-index-map;SNS_OutsideSpatialSupport=NaN;SNS_MappedSpatialPixels=" + num2istr(nMappedPixels) + ";"
    Note/K mapW, mapNote

    KillWaves/Z SNS_tmp_LDOS_point_index_xyz, $pointIndexMapName
    return mapW
End


//==============================================================================
// SNS_MakeLDOS_E_Maps_FromSpectra
//
// Purpose:
//   Convert point-list LDOS spectra from SNS_LDOSmap_BFixed_FromMask into
//   3D image stacks suitable for E-slice display.
//
// Inputs:
//   xyWaveName       : XY_* wave, [point][x/y] in nm.
//   rawLDOSName      : LDOS_E_raw_* wave, [energy][point].
//   convLDOSName     : LDOS_E_conv_* wave, [energy][point].
//   validWaveName    : ValidPt_* wave.
//
// Optional Inputs:
//   outRawName       : output raw stack name.
//                      Default: "LDOS_E_map".
//   outConvName      : output broadened stack name.
//                      Default: "LDOS_conv_E_map".
//   xyzTol           : XYZTripletToMatrixEx tolerance.
//                      Default: 1e-05.
//   nx, ny           : output grid dimensions.
//                      If omitted, infer a square grid from the total number
//                      of input positions.
//
// Outputs:
//   outRawName       : raw LDOS map stack, [x][y][energy].
//   outConvName      : broadened LDOS map stack, [x][y][energy].
//   LDOS_E_axis_eV   : energy axis copied from rawLDOSName x scaling.
//
// Notes:
//   This helper keeps XY coordinates separate from the stored LDOS(E, point)
//   waves. It builds temporary XYZ triplets and interpolates each energy slice
//   with SNS_XYZTripletToMatrix.
//==============================================================================
Function SNS_MakeLDOS_E_Maps_FromSpectra(xyWaveName, rawLDOSName, convLDOSName, validWaveName, [outRawName, outConvName, xyzTol, nx, ny])
    String xyWaveName, rawLDOSName, convLDOSName, validWaveName
    String outRawName, outConvName
    Variable xyzTol, nx, ny

    if (ParamIsDefault(outRawName))
        outRawName = "LDOS_E_map"
    endif
    if (ParamIsDefault(outConvName))
        outConvName = "LDOS_conv_E_map"
    endif
    if (ParamIsDefault(xyzTol))
        xyzTol = 1e-05
    endif

    Wave/Z xyW = $xyWaveName
    Wave/Z rawW = $rawLDOSName
    Wave/Z convW = $convLDOSName

    if (!WaveExists(xyW))
        Abort "SNS_MakeLDOS_E_Maps_FromSpectra: missing XY wave: " + xyWaveName
    endif
    if (!WaveExists(rawW))
        Abort "SNS_MakeLDOS_E_Maps_FromSpectra: missing raw LDOS wave: " + rawLDOSName
    endif
    if (!WaveExists(convW))
        Abort "SNS_MakeLDOS_E_Maps_FromSpectra: missing convolved LDOS wave: " + convLDOSName
    endif

    Variable nE = DimSize(rawW, 0)
    Variable nPts = DimSize(rawW, 1)

    if ((nE <= 0) || (nPts <= 0))
        Abort "SNS_MakeLDOS_E_Maps_FromSpectra: LDOS waves must be [energy][point]."
    endif
    if ((DimSize(convW, 0) != nE) || (DimSize(convW, 1) != nPts))
        Abort "SNS_MakeLDOS_E_Maps_FromSpectra: raw and convolved LDOS dimensions differ."
    endif
    if ((DimSize(xyW, 0) != nPts) || (DimSize(xyW, 1) < 2))
        Abort "SNS_MakeLDOS_E_Maps_FromSpectra: XY wave must be [point][x/y] and match LDOS point count."
    endif

    Wave/Z validW = $validWaveName
    if (!WaveExists(validW))
        Abort "SNS_MakeLDOS_E_Maps_FromSpectra: missing ValidPt wave: " + validWaveName
    endif
    if (DimSize(validW, 0) != nPts)
        Abort "SNS_MakeLDOS_E_Maps_FromSpectra: ValidPt wave length does not match LDOS point count."
    endif
    Variable ip, k
    Variable nValid = 0
    for (ip = 0; ip < nPts; ip += 1)
        if (validW[ip] == 0)
            continue
        endif
        if ((numtype(xyW[ip][0]) == 0) && (numtype(xyW[ip][1]) == 0))
            nValid += 1
        endif
    endfor

    if (nValid <= 0)
        Abort "SNS_MakeLDOS_E_Maps_FromSpectra: no valid points."
    endif

    if (ParamIsDefault(nx) || ParamIsDefault(ny) || (nx <= 0) || (ny <= 0))
        Variable nGrid = round(sqrt(nPts))
        if (nGrid * nGrid != nPts)
            Abort "SNS_MakeLDOS_E_Maps_FromSpectra: input positions do not form a square grid. Supply nx and ny."
        endif
        nx = nGrid
        ny = nGrid
    endif

    if (round(nx) * round(ny) != nPts)
        Abort "SNS_MakeLDOS_E_Maps_FromSpectra: nx*ny does not match input position count."
    endif
    nx = round(nx)
    ny = round(ny)

    Make/O/D/N=(nValid, 3) SNS_tmpXYZ_raw
    Make/O/D/N=(nValid, 3) SNS_tmpXYZ_conv
    Wave tmpRawXYZ = SNS_tmpXYZ_raw
    Wave tmpConvXYZ = SNS_tmpXYZ_conv

    String tmpRawMapName = "SNS_tmp_LDOS_raw_map"
    String tmpConvMapName = "SNS_tmp_LDOS_conv_map"

    Variable j
    for (k = 0; k < nE; k += 1)
        j = 0
        for (ip = 0; ip < nPts; ip += 1)
            if (validW[ip] == 0)
                continue
            endif
            if ((numtype(xyW[ip][0]) != 0) || (numtype(xyW[ip][1]) != 0))
                continue
            endif

            tmpRawXYZ[j][0] = xyW[ip][0]
            tmpRawXYZ[j][1] = xyW[ip][1]
            tmpRawXYZ[j][2] = rawW[k][ip]

            tmpConvXYZ[j][0] = xyW[ip][0]
            tmpConvXYZ[j][1] = xyW[ip][1]
            tmpConvXYZ[j][2] = convW[k][ip]

            j += 1
        endfor

        Wave tmpXYZRaw = SNS_tmpXYZ_raw
        Wave tmpXYZConv = SNS_tmpXYZ_conv
        SNS_XYZTripletToMatrix(tmpXYZRaw, tmpRawMapName, nx, ny, xyzTol)
        SNS_XYZTripletToMatrix(tmpXYZConv, tmpConvMapName, nx, ny, xyzTol)

        Wave tmpRawMapW = $tmpRawMapName
        Wave tmpConvMapW = $tmpConvMapName

        if (k == 0)
            Make/O/D/N=(DimSize(tmpRawMapW, 0), DimSize(tmpRawMapW, 1), nE) $outRawName
            Make/O/D/N=(DimSize(tmpConvMapW, 0), DimSize(tmpConvMapW, 1), nE) $outConvName
            Wave LDOSrawInitW = $outRawName
            Wave LDOSconvInitW = $outConvName

            LDOSrawInitW = NaN
            LDOSconvInitW = NaN

            SetScale/P x, DimOffset(tmpRawMapW, 0), DimDelta(tmpRawMapW, 0), "nm", LDOSrawInitW
            SetScale/P y, DimOffset(tmpRawMapW, 1), DimDelta(tmpRawMapW, 1), "nm", LDOSrawInitW
            SetScale/P z, DimOffset(rawW, 0), DimDelta(rawW, 0), WaveUnits(rawW, 0), LDOSrawInitW

            SetScale/P x, DimOffset(tmpConvMapW, 0), DimDelta(tmpConvMapW, 0), "nm", LDOSconvInitW
            SetScale/P y, DimOffset(tmpConvMapW, 1), DimDelta(tmpConvMapW, 1), "nm", LDOSconvInitW
            SetScale/P z, DimOffset(convW, 0), DimDelta(convW, 0), WaveUnits(convW, 0), LDOSconvInitW
        endif

        Wave LDOSrawFillW = $outRawName
        Wave LDOSconvFillW = $outConvName
        LDOSrawFillW[][][k] = tmpRawMapW[p][q]
        LDOSconvFillW[][][k] = tmpConvMapW[p][q]
    endfor

    Make/O/D/N=(nE) LDOS_E_axis_eV
    Wave eAxisW = LDOS_E_axis_eV
    eAxisW = DimOffset(rawW, 0) + p * DimDelta(rawW, 0)
    SetScale/P x, DimOffset(rawW, 0), DimDelta(rawW, 0), WaveUnits(rawW, 0), eAxisW

    KillWaves/Z SNS_tmpXYZ_raw, SNS_tmpXYZ_conv
    KillWaves/Z $tmpRawMapName, $tmpConvMapName

    Wave LDOSrawDoneW = $outRawName
    Wave LDOSconvDoneW = $outConvName
    String rawNote = note(rawW)
    if ((strlen(rawNote) > 0) && (cmpstr(rawNote[strlen(rawNote)-1, strlen(rawNote)-1], ";") != 0))
        rawNote += ";"
    endif
    rawNote += "SNS_MapSource=" + GetWavesDataFolder(rawW, 2) + ";SNS_MapDims=x,y,energy;SNS_MapInterpolation=Voronoi-per-energy-slice;"
    Note/K LDOSrawDoneW, rawNote

    String convNote = note(convW)
    if ((strlen(convNote) > 0) && (cmpstr(convNote[strlen(convNote)-1, strlen(convNote)-1], ";") != 0))
        convNote += ";"
    endif
    convNote += "SNS_MapSource=" + GetWavesDataFolder(convW, 2) + ";SNS_MapDims=x,y,energy;SNS_MapInterpolation=Voronoi-per-energy-slice;"
    Note/K LDOSconvDoneW, convNote

    return 0
End

//==============================================================================
// SNS_Plot_ExpSim_Cuts
//
// Purpose:
//   Plot diagnostic cuts comparing experimental dI/dV(E,B) and simulated DOS(E,B).
//
//   The function extracts and plots:
//
//     1) Zero-bias / zero-energy line cut vs B:
//          expEB(E=0, B) and simEB(E=0, B)
//
//     2) Zero-field line cut vs E:
//          expEB(E, B=0) and simEB(E, B=0)
//
//   In the zero-field E-cut plot, the function indicates superconducting gap
//   positions as vertical dashed lines.
//
//   Gap annotation logic:
//     - If root:v_cfg_Delta2D_eV exists:
//           annotate ±Δ\BSS using root:v_cfg_Delta2D_eV.
//     - If root:v_cfg_Delta3D_eV exists:
//           additionally annotate ±Δ\Bbulk using root:v_cfg_Delta3D_eV.
//     - If root:v_cfg_Delta2D_eV does not exist, fall back to
//           root:SNS_Settings:Delta and annotate ±Δ.
//
// Inputs:
//   expEB : input experimental dI/dV(E,B) wave.
//           Dimension 0 is energy, dimension 1 is magnetic field.
//
//   simEB : input simulated DOS(E,B) wave.
//           Dimension 0 is energy, dimension 1 is magnetic field.
//
// Outputs:
//   return : numeric
//            0 on successful completion.
//
//==============================================================================

Function SNS_Plot_ExpSim_Cuts(expEB, simEB)
    Wave expEB, simEB

    if (DimSize(expEB,0) <= 1 || DimSize(expEB,1) <= 0)
        Abort "SNS_Plot_ExpSim_Cuts: expEB invalid."
    endif
    if (DimSize(simEB,0) <= 1 || DimSize(simEB,1) <= 0)
        Abort "SNS_Plot_ExpSim_Cuts: simEB invalid."
    endif

    //========================================================
    // a) Zero-bias / zero-energy vs B  (E = 0)
    //     -> select scaled x-range (0,0) and all B
    //========================================================

    Duplicate/O/R=(0)() expEB, ZBC
    Duplicate/O/R=(0)() simEB, DOS0

    Redimension/N=(DimSize(expEB,1)) ZBC
    Redimension/N=(DimSize(simEB,1)) DOS0

    SetScale/P x DimOffset(expEB,1), DimDelta(expEB,1), WaveUnits(expEB,1), ZBC
    SetScale/P x DimOffset(simEB,1)*1e3, DimDelta(simEB,1)*1e3, WaveUnits(expEB,1), DOS0

    Display/N=ZBC_Cut_Graph ZBC
    AppendToGraph DOS0
    ModifyGraph mode(DOS0)=0, rgb(DOS0)=(1,16019,65535)
    ModifyGraph tick=2, mirror=2, standoff=0
    Label left "ZBC (arb. u.)"
    Label bottom "B (\\U)"
    SetAxis left 0,*
    ModifyGraph lblMargin(left)=4, lblMargin(bottom)=4
    ModifyGraph tickUnit(left)=1, tickUnit(bottom)=1
    ModifyGraph zero(left)=1, zero(bottom)=1
    MoveWindow/W=ZBC_Cut_Graph 50,80,480,420

    //========================================================
    // b) Zero-field vs E  (B = 0)
    //     -> select all E and scaled y-range (0,0)
    //========================================================

    Duplicate/O/R=()(0) expEB, dIdV_B0
    Duplicate/O/R=()(0) simEB, DOS_B0

    SetScale/P x DimOffset(simEB,0)*1e3, DimDelta(simEB,0)*1e3, WaveUnits(expEB,0), DOS_B0

    Display/N=E_Cut_Graph dIdV_B0
    AppendToGraph DOS_B0
    ModifyGraph mode(DOS_B0)=0, rgb(DOS_B0)=(1,16019,65535)
    ModifyGraph tick=2, mirror=2, standoff=0
    Label left "LDOS (arb. u.)"
    Label bottom "E (\\U)"
    SetAxis left 0,*
    ModifyGraph lblMargin(left)=4, lblMargin(bottom)=4
    ModifyGraph tickUnit(left)=1, tickUnit(bottom)=1
    ModifyGraph zero(left)=1, zero(bottom)=1
    MoveWindow/W=E_Cut_Graph 500,80,930,420

    //========================================================
    // Indicate gap scales as dashed vertical lines
    //========================================================

    String gname2 = WinName(0,1,1)
    if (strlen(gname2) > 0)

        SetDrawLayer/W=$gname2 UserFront

        NVAR/Z DeltaSS = root:v_cfg_Delta2D_eV
        NVAR/Z DeltaBulk = root:v_cfg_Delta3D_eV
        NVAR/Z DeltaFallback = root:SNS_Settings:Delta

        Variable deltaSS_meV, deltaBulk_meV, deltaFallback_meV
        Variable offSS, offBulk, offFallback

        // ------------------------------
        // Delta_SS from root:v_cfg_Delta2D_eV
        // ------------------------------
        if (NVAR_Exists(DeltaSS))

            deltaSS_meV = DeltaSS * 1e3
            offSS = abs(deltaSS_meV) * 0.10

            SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, dash=1, linethick=1
            DrawLine/W=$gname2 deltaSS_meV, 0, deltaSS_meV, 1
            SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, textyjust=1
            DrawText/W=$gname2 deltaSS_meV + offSS, 0.58, "Δ\\BSS"

            SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, dash=1, linethick=1
            DrawLine/W=$gname2 -deltaSS_meV, 0, -deltaSS_meV, 1
            SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, textxjust=2, textyjust=1
            DrawText/W=$gname2 -deltaSS_meV - offSS, 0.58, "-Δ\\BSS"

        else

            // Fallback to old behavior if explicit 2D gap config is absent.
            if (NVAR_Exists(DeltaFallback))

                deltaFallback_meV = DeltaFallback * 1e3
                offFallback = abs(deltaFallback_meV) * 0.10

                SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, dash=1, linethick=1
                DrawLine/W=$gname2 deltaFallback_meV, 0, deltaFallback_meV, 1
                SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, textyjust=1
                DrawText/W=$gname2 deltaFallback_meV + offFallback, 0.50, "Δ"

                SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, dash=1, linethick=1
                DrawLine/W=$gname2 -deltaFallback_meV, 0, -deltaFallback_meV, 1
                SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, textxjust=2, textyjust=1
                DrawText/W=$gname2 -deltaFallback_meV - offFallback, 0.50, "-Δ"

            endif
        endif

        // ------------------------------
        // Delta_bulk from root:v_cfg_Delta3D_eV
        // ------------------------------
        if (NVAR_Exists(DeltaBulk))

            deltaBulk_meV = DeltaBulk * 1e3
            offBulk = abs(deltaBulk_meV) * 0.10

            SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, dash=1, linethick=1
            DrawLine/W=$gname2 deltaBulk_meV, 0, deltaBulk_meV, 1
            SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, textyjust=1
            DrawText/W=$gname2 deltaBulk_meV + offBulk, 0.42, "Δ\\Bbulk"

            SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, dash=1, linethick=1
            DrawLine/W=$gname2 -deltaBulk_meV, 0, -deltaBulk_meV, 1
            SetDrawEnv/W=$gname2 xcoord=bottom, ycoord=axrel, textxjust=2, textyjust=1
            DrawText/W=$gname2 -deltaBulk_meV - offBulk, 0.42, "-Δ\\Bbulk"

        endif
    endif

    return 0
End

//==============================================================================
// SNS_TestDrawRotatedBox
//
// Purpose:
//   Draw a test rotated polygon box on the top graph using DrawPoly.
//
// Inputs:
//   x0       : box center x position in graph-axis units.
//   y0       : box center y position in graph-axis units.
//   angleDeg : box angle in degrees.
//   L        : box length along angle direction.
//   W        : box width perpendicular to angle direction.
//
// Outputs:
//   return : numeric
//            0 on successful completion.
//
//==============================================================================

Function SNS_TestDrawRotatedBox(x0, y0, angleDeg, L, W)
    Variable x0, y0, angleDeg, L, W

    String win = WinName(0, 1, 1)
    if (strlen(win) == 0)
        Abort "SNS_TestDrawRotatedBox: No graph window found."
    endif

    NewDataFolder/O root:Packages
    NewDataFolder/O root:Packages:SNS_TempDraw

    Make/O/D/N=5 root:Packages:SNS_TempDraw:wTestBoxX
    Make/O/D/N=5 root:Packages:SNS_TempDraw:wTestBoxY

    Wave wX = root:Packages:SNS_TempDraw:wTestBoxX
    Wave wY = root:Packages:SNS_TempDraw:wTestBoxY

    Variable ang = angleDeg*pi/180

    Variable ux = cos(ang)
    Variable uy = sin(ang)
    Variable nx = -sin(ang)
    Variable ny = cos(ang)

    Variable hL = 0.5*L
    Variable hW = 0.5*W

    wX[0] = -hL*ux - hW*nx
    wY[0] = -hL*uy - hW*ny

    wX[1] =  hL*ux - hW*nx
    wY[1] =  hL*uy - hW*ny

    wX[2] =  hL*ux + hW*nx
    wY[2] =  hL*uy + hW*ny

    wX[3] = -hL*ux + hW*nx
    wY[3] = -hL*uy + hW*ny

    wX[4] = wX[0]
    wY[4] = wY[0]

    SetDrawLayer/W=$win UserFront
    SetDrawEnv/W=$win xcoord=bottom, ycoord=left, fillfgc=(65535,0,0,20000), linefgc=(65535,0,0), linethick=2
    DrawPoly/W=$win x0, y0, 1, 1, wX, wY

    return 0
End

//==============================================================================
// Display2DEG_OnsetOnLineSTS
//
// Purpose:
//   Display a large-bias LineSTS dI/dV map and overlay valid fitted 2DEG onset
//   energies as black cross markers.
//
//   The function detects which dIdV axis is the energy axis from wave units:
//
//      V, mV, eV, meV
//
//   If energy is dim 0:
//      x = onset energy
//      y = line position
//
//   If energy is dim 1:
//      x = line position
//      y = onset energy
//
//   Invalid NaN onset values are removed before plotting.
//   The label is placed next to the middle valid marker.
//
// Inputs:
//   dIdV      : 2D LineSTS wave to display.
//   onsetWave : 1D fitted onset wave from Fit2DEG_OnsetAndBroadening_2D.
//   labelText : optional label text. Default: "E\\BSS\\M".
//
// Outputs:
//   return : numeric
//            0 on successful completion.
//
//==============================================================================

Function Display2DEG_OnsetOnLineSTS(dIdV, onsetWave, [labelText])
    Wave dIdV
    Wave onsetWave
    String labelText

    if (ParamIsDefault(labelText))
        labelText = "E\\BSS\\M"
    endif

    if (DimSize(dIdV, 0) <= 0 || DimSize(dIdV, 1) <= 0)
        Abort "Display2DEG_OnsetOnLineSTS: dIdV must be a 2D wave."
    endif

    if (numpnts(onsetWave) <= 0)
        Abort "Display2DEG_OnsetOnLineSTS: onsetWave is empty."
    endif

    // ------------------------------
    // Detect energy axis from units
    // ------------------------------
    String dim0Unit = WaveUnits(dIdV, 0)
    String dim1Unit = WaveUnits(dIdV, 1)

    String dim0UnitLower = LowerStr(dim0Unit)
    String dim1UnitLower = LowerStr(dim1Unit)

    Variable dim0IsEnergy = 0
    Variable dim1IsEnergy = 0

    if (cmpstr(dim0UnitLower, "v") == 0 || cmpstr(dim0UnitLower, "mv") == 0 || cmpstr(dim0UnitLower, "ev") == 0 || cmpstr(dim0UnitLower, "mev") == 0)
        dim0IsEnergy = 1
    endif

    if (cmpstr(dim1UnitLower, "v") == 0 || cmpstr(dim1UnitLower, "mv") == 0 || cmpstr(dim1UnitLower, "ev") == 0 || cmpstr(dim1UnitLower, "mev") == 0)
        dim1IsEnergy = 1
    endif

    if (dim0IsEnergy && dim1IsEnergy)
        Abort "Display2DEG_OnsetOnLineSTS: both dIdV axes look like energy axes."
    endif

    if (!dim0IsEnergy && !dim1IsEnergy)
        Abort "Display2DEG_OnsetOnLineSTS: no energy axis found. Expected units V, mV, eV, or meV."
    endif

    Variable energyDim
    Variable lineDim

    if (dim0IsEnergy)
        energyDim = 0
        lineDim = 1
    else
        energyDim = 1
        lineDim = 0
    endif

    Variable linePointCount = DimSize(dIdV, lineDim)

    if (numpnts(onsetWave) != linePointCount)
        Abort "Display2DEG_OnsetOnLineSTS: onsetWave length does not match non-energy dimension of dIdV."
    endif

    // ------------------------------
    // Build overlay waves
    // ------------------------------
    NewDataFolder/O root:Packages
    NewDataFolder/O root:Packages:Onset2DEG_Display

    Duplicate/O onsetWave, root:Packages:Onset2DEG_Display:w_OnsetEnergy_OnsetOverlay
    Make/O/D/N=(linePointCount) root:Packages:Onset2DEG_Display:w_LinePos_OnsetOverlay

    Wave w_OnsetEnergy_OnsetOverlay = root:Packages:Onset2DEG_Display:w_OnsetEnergy_OnsetOverlay
    Wave w_LinePos_OnsetOverlay = root:Packages:Onset2DEG_Display:w_LinePos_OnsetOverlay

    w_LinePos_OnsetOverlay = DimOffset(dIdV, lineDim) + DimDelta(dIdV, lineDim) * p

    // Remove invalid onset entries without loops.
    w_LinePos_OnsetOverlay = numtype(w_OnsetEnergy_OnsetOverlay[p]) == 0 ? w_LinePos_OnsetOverlay[p] : NaN
    WaveTransform zapNaNs w_LinePos_OnsetOverlay
    WaveTransform zapNaNs w_OnsetEnergy_OnsetOverlay

    if (numpnts(w_OnsetEnergy_OnsetOverlay) <= 0)
        Abort "Display2DEG_OnsetOnLineSTS: no valid onset values after removing NaNs."
    endif

    // ------------------------------
    // Display and overlay
    // ------------------------------
    SNS_DisplayWithScales(dIdV, cmap="davos")

    String graphName = WinName(0, 1, 1)
    Variable labelIndex = floor((numpnts(w_OnsetEnergy_OnsetOverlay) - 1) / 2)

    if (energyDim == 0)

        // dIdV energy axis is graph bottom axis.
        // Overlay as y-position versus x-energy.
        AppendToGraph/W=$graphName w_LinePos_OnsetOverlay vs w_OnsetEnergy_OnsetOverlay

        ModifyGraph/W=$graphName mode(w_LinePos_OnsetOverlay)=3
        ModifyGraph/W=$graphName marker(w_LinePos_OnsetOverlay)=0
        ModifyGraph/W=$graphName rgb(w_LinePos_OnsetOverlay)=(0,0,0)
        ModifyGraph/W=$graphName msize(w_LinePos_OnsetOverlay)=4
        ModifyGraph/W=$graphName mrkThick(w_LinePos_OnsetOverlay)=1.5

        SetDrawLayer/W=$graphName UserFront
        SetDrawEnv/W=$graphName xcoord=bottom, ycoord=left, textrgb=(0,0,0), fsize=12, textxjust=0, textyjust=1
        DrawText/W=$graphName w_OnsetEnergy_OnsetOverlay[labelIndex], w_LinePos_OnsetOverlay[labelIndex], labelText

    else

        // dIdV energy axis is graph left axis.
        // Overlay as y-energy versus x-position.
        AppendToGraph/W=$graphName w_OnsetEnergy_OnsetOverlay vs w_LinePos_OnsetOverlay

        ModifyGraph/W=$graphName mode(w_OnsetEnergy_OnsetOverlay)=3
        ModifyGraph/W=$graphName marker(w_OnsetEnergy_OnsetOverlay)=0
        ModifyGraph/W=$graphName rgb(w_OnsetEnergy_OnsetOverlay)=(0,0,0)
        ModifyGraph/W=$graphName msize(w_OnsetEnergy_OnsetOverlay)=4
        ModifyGraph/W=$graphName mrkThick(w_OnsetEnergy_OnsetOverlay)=1.5

        SetDrawLayer/W=$graphName UserFront
        SetDrawEnv/W=$graphName xcoord=bottom, ycoord=left, textrgb=(0,0,0), fsize=12, textxjust=0, textyjust=1
        DrawText/W=$graphName w_LinePos_OnsetOverlay[labelIndex], w_OnsetEnergy_OnsetOverlay[labelIndex], labelText

    endif
    
    ModifyGraph tick=2;DelayUpdate
	Label left "Position (nm)";DelayUpdate
	Label bottom "Bias (mV)"
    return 0
End


//==============================================================================
// SNS_DisplayTopoWithLineSTS
//
// Purpose:
//   Display a topography image and overlay up to five LineSTS trajectories.
//
//   Each LineSTS folder must contain:
//
//      lineSTSFolder:pos:stspos
//
//   where stspos is a 2-row, N-column wave:
//      stspos[0][i] = x position
//      stspos[1][i] = y position
//
// Inputs:
//   topo    : topography image wave.
//   lineDF1 : first LineSTS folder.
//   lineDF2 : optional second LineSTS folder.
//   lineDF3 : optional third LineSTS folder.
//   lineDF4 : optional fourth LineSTS folder.
//   lineDF5 : optional fifth LineSTS folder.
//
// Outputs:
//   return : numeric
//            0 on successful completion.
//
//==============================================================================

Function SNS_DisplayTopoWithLineSTS(topo, lineDF1, [lineDF2, lineDF3, lineDF4, lineDF5])
    Wave topo
    String lineDF1
    String lineDF2
    String lineDF3
    String lineDF4
    String lineDF5

    String winTopo
    String lineFolderList
    String lineFolder
    String posWavePath
    String lineLabel

    Variable lineIndex
    Variable nLines
    Variable lastPointIndex

    Variable lineRed
    Variable lineGreen
    Variable lineBlue

    Wave/Z linePosWave

    // ------------------------------
    // Build folder list
    // ------------------------------
    lineFolderList = lineDF1 + ";"

    if (!ParamIsDefault(lineDF2))
        lineFolderList += lineDF2 + ";"
    endif
    if (!ParamIsDefault(lineDF3))
        lineFolderList += lineDF3 + ";"
    endif
    if (!ParamIsDefault(lineDF4))
        lineFolderList += lineDF4 + ";"
    endif
    if (!ParamIsDefault(lineDF5))
        lineFolderList += lineDF5 + ";"
    endif

    nLines = ItemsInList(lineFolderList, ";")

    if (nLines <= 0)
        Abort "SNS_DisplayTopoWithLineSTS: no LineSTS folders supplied."
    endif

    // ------------------------------
    // Display topography
    // ------------------------------
    SNS_DisplayWithScales(topo, cmap="gray")
    winTopo = WinName(0, 1, 1)

    if (strlen(winTopo) == 0)
        Abort "SNS_DisplayTopoWithLineSTS: no graph window found."
    endif

    SetDrawLayer/W=$winTopo UserFront

    // ------------------------------
    // Draw LineSTS trajectories
    // ------------------------------
    for (lineIndex = 0; lineIndex < nLines; lineIndex += 1)

        lineFolder = StringFromList(lineIndex, lineFolderList, ";")

        if (strlen(lineFolder) == 0)
            continue
        endif

        if (cmpstr(lineFolder[strlen(lineFolder)-1, strlen(lineFolder)-1], ":") != 0)
            lineFolder += ":"
        endif

        posWavePath = lineFolder + "pos:stspos"

        Wave/Z linePosWave = $posWavePath

        if (!WaveExists(linePosWave))
            Print "SNS_DisplayTopoWithLineSTS: missing ", posWavePath
            continue
        endif

        if (DimSize(linePosWave, 0) < 2 || DimSize(linePosWave, 1) < 2)
            Print "SNS_DisplayTopoWithLineSTS: invalid stspos wave in ", lineFolder
            continue
        endif

        lastPointIndex = DimSize(linePosWave, 1) - 1

		// Color cycle for up to five lines.
		// line 1: black, line 2: blue, line 3: red, line 4: green, line 5: magenta
		if (lineIndex == 0)
		    lineRed = 0;     lineGreen = 0;     lineBlue = 0
		elseif (lineIndex == 1)
		    lineRed = 0;     lineGreen = 0;     lineBlue = 65535
		elseif (lineIndex == 2)
		    lineRed = 65535; lineGreen = 0;     lineBlue = 0
		elseif (lineIndex == 3)
		    lineRed = 0;     lineGreen = 45000; lineBlue = 0
		else
		    lineRed = 65535; lineGreen = 0;     lineBlue = 65535
		endif

        SetDrawEnv/W=$winTopo xcoord=bottom, ycoord=left, linefgc=(lineRed,lineGreen,lineBlue), linethick=2, arrow=1

        DrawLine/W=$winTopo linePosWave[0][0], linePosWave[1][0], linePosWave[0][lastPointIndex], linePosWave[1][lastPointIndex]

        lineLabel = "LineSTS " + num2str(lineIndex + 1)

        SetDrawEnv/W=$winTopo xcoord=bottom, ycoord=left, textrgb=(lineRed,lineGreen,lineBlue), fsize=12, textxjust=0, textyjust=1

        DrawText/W=$winTopo linePosWave[0][lastPointIndex], linePosWave[1][lastPointIndex], lineLabel

    endfor

    return 0
End

//==============================================================================
// SNS_Display_ResolutionMatrix
//
// Unified display helper for resolution matrices.
//
// matrixKind:
//   0 = ZeroField  -> displays ZeroFieldRes_LDOS_ERes_n
//   1 = ZeroEnergy -> displays ZeroEnergyRes_ZBC_BRes, not normalized
//
// y-axis / resolution labels are identical for both:
//   860 mK
//   400 mK
//   10 mK / 1 um
//   10 um
//   100 um
//
// For ZeroField:
//   x-axis is energy, restricted to -deltaEff_eV ... +deltaEff_eV.
//
// For ZeroEnergy:
//   x-axis is magnetic field B, using the full stored B range.
//==============================================================================

Function SNS_Display_ResolutionMatrix(dataFolder, [matrixKind, resFolderName, outPrefix, winName, deltaEff_eV])
	String dataFolder
	Variable matrixKind
	String resFolderName, outPrefix, winName
	Variable deltaEff_eV

	if (ParamIsDefault(matrixKind))
		matrixKind = 0
	endif

	matrixKind = round(matrixKind)

	if ((matrixKind != 0) && (matrixKind != 1))
		Abort "SNS_Display_ResolutionMatrix: matrixKind must be 0=ZeroField or 1=ZeroEnergy."
	endif

	if (ParamIsDefault(resFolderName))
		if (matrixKind == 0)
			resFolderName = "ZeroField_ResolutionMatrix"
		else
			resFolderName = "ZeroEnergy_ResolutionMatrix"
		endif
	endif

	if (ParamIsDefault(outPrefix))
		if (matrixKind == 0)
			outPrefix = "ZeroFieldRes"
		else
			outPrefix = "ZeroEnergyRes"
		endif
	endif

	if (ParamIsDefault(winName))
		if (matrixKind == 0)
			winName = "ZeroFieldRes_LDOS_n"
		else
			winName = "ZeroEnergyRes_ZBC"
		endif
	endif

	String savedDF = GetDataFolder(1)
	String resDF = dataFolder + ":" + resFolderName

	if (!DataFolderExists(resDF + ":"))
		Abort "SNS_Display_ResolutionMatrix: resolution folder does not exist."
	endif

	String matrixName
	if (matrixKind == 0)
		matrixName = outPrefix + "_LDOS_ERes_n"
	else
		matrixName = outPrefix + "_ZBC_BRes"
	endif

	WAVE/Z M     = $(resDF + ":" + matrixName)
	WAVE/Z T_K   = $(resDF + ":" + outPrefix + "_T_K")
	WAVE/Z Lphi  = $(resDF + ":" + outPrefix + "_Lphi_um")

	if (!WaveExists(M))
		Abort "SNS_Display_ResolutionMatrix: missing matrix wave: " + matrixName
	endif
	if (!WaveExists(T_K) || !WaveExists(Lphi))
		Abort "SNS_Display_ResolutionMatrix: missing T_K or Lphi_um column waves."
	endif

	SetDataFolder $resDF

	SNS_DisplayWithScales(M, cmap="davos")
	DoWindow/C $winName

	// -------------------------------------------------------------------------
	// x-axis handling.
	// -------------------------------------------------------------------------
	if (matrixKind == 0)

		if (ParamIsDefault(deltaEff_eV))
			if (Exists("root:SNS_Settings:delta_eff") == 2)
				NVAR delta_eff_var = root:SNS_Settings:delta_eff
				deltaEff_eV = delta_eff_var
			elseif (Exists("root:SNS_Settings:Delta_eff") == 2)
				NVAR Delta_eff_var = root:SNS_Settings:Delta_eff
				deltaEff_eV = Delta_eff_var
			elseif (Exists("root:SNS_Settings:Delta") == 2)
				NVAR Delta_var = root:SNS_Settings:Delta
				deltaEff_eV = Delta_var
			else
				SetDataFolder $savedDF
				Abort "SNS_Display_ResolutionMatrix: supply deltaEff_eV or define root:SNS_Settings:delta_eff / Delta_eff / Delta."
			endif
		endif

		if (numtype(deltaEff_eV) != 0 || deltaEff_eV <= 0)
			SetDataFolder $savedDF
			Abort "SNS_Display_ResolutionMatrix: invalid deltaEff_eV."
		endif

		SetAxis bottom, -deltaEff_eV, deltaEff_eV
	endif

	// For matrixKind=1, keep full B range from the wave scaling.

	// -------------------------------------------------------------------------
	// Common y-axis / resolution ticks.
	// -------------------------------------------------------------------------
	Variable c860mK = SNS__FindClosestColumn_T_Lphi(T_K, Lphi, 0.86, 1, 1)
	Variable c400mK = SNS__FindClosestColumn_T_Lphi(T_K, Lphi, 0.4, 1, 1)
	Variable c10mK1 = SNS__FindClosestColumn_T_Lphi(T_K, Lphi, 0.01, 1, 1)
	Variable c10um  = SNS__FindClosestColumn_T_Lphi(T_K, Lphi, 0.01, 10, 0)
	Variable c100um = SNS__FindClosestColumn_T_Lphi(T_K, Lphi, 0.01, 100, 0)

	// Crop and reverse vertical axis:
	// top = 860 mK column, bottom = last resolution column.
	Variable nCols = DimSize(M, 1)
	Variable yTop = max(0, c860mK - 0.5)
	Variable yBot = nCols - 0.5
	SetAxis left, yBot, yTop

	Make/O/D/N=5 $(outPrefix + "_ytick_pos")
	Make/O/T/N=5 $(outPrefix + "_ytick_label")

	WAVE ytick = $(outPrefix + "_ytick_pos")
	WAVE/T ylab = $(outPrefix + "_ytick_label")

	ytick[0] = c860mK
	ytick[1] = c400mK
	ytick[2] = c10mK1
	ytick[3] = c10um
	ytick[4] = c100um

	ylab[0] = "860 mK"
	ylab[1] = "400 mK"
	ylab[2] = "10 mK / 1 um"
	ylab[3] = "10 um"
	ylab[4] = "100 um"

	ModifyGraph userticks(left)={ytick, ylab}
	ModifyGraph margin(left)=100

	Label left "Resolution";DelayUpdate

	if (matrixKind == 0)
		Label bottom "Energy (\\U)"
	else
		Label bottom "B (T)"
	endif

	// Separator between:
	//   above: fixed Lphi = LphiMin_um, varying T
	//   below: fixed T = Tmin_K, varying Lphi
	Variable cSep = c10mK1 + 0.5
	
	GetAxis/Q bottom
	Variable x0 = V_min
	Variable x1 = V_max
	Variable xMid = 0.5*(x0 + x1)
	
	SetDrawLayer UserFront
	
	// Dashed white separator line.
	SetDrawEnv xcoord=bottom, ycoord=left, dash=3, linethick=1, linefgc=(65535,65535,65535)
	DrawLine x0, cSep, x1, cSep
	
	// y-axis is reversed:
	//   visually above the line -> smaller y value
	//   visually below the line -> larger y value
	Variable dyText = 0.8
	
	// Above dashed line: temperature sweep at fixed Lphi.
	SetDrawEnv xcoord=bottom, ycoord=left, fsize=10, textxjust=1, textyjust=2, textrgb=(65535,65535,65535)
	DrawText xMid, cSep + dyText, "↓ fixed T = 10 mK; varying L\\BΦ\\M"
	
	// Below dashed line: Lphi sweep at fixed temperature.
	SetDrawEnv xcoord=bottom, ycoord=left, fsize=10, textxjust=1, textyjust=0, textrgb=(65535,65535,65535)
	DrawText xMid, cSep - dyText, "↑ fixed L\\BΦ\\M = 1 μm; varying T"
	
		// -------------------------------------------------------------------------
	// Second window: selected resolution traces.
	//   left axis : 860 mK
	//   right axis: 10 mK / 1 um, 10 um, 100 um
	// -------------------------------------------------------------------------
	Variable nX = DimSize(M, 0)
	Variable xOff = DimOffset(M, 0)
	Variable xDel = DimDelta(M, 0)
	String xUnit = WaveUnits(M, 0)

	String traceWinName = winName + "_traces"

	Make/O/D/N=(nX) $(outPrefix + "_trace_860mK")
	Make/O/D/N=(nX) $(outPrefix + "_trace_10mK_1um")
	Make/O/D/N=(nX) $(outPrefix + "_trace_10um")
	Make/O/D/N=(nX) $(outPrefix + "_trace_100um")

	WAVE tr860   = $(outPrefix + "_trace_860mK")
	WAVE tr10mK  = $(outPrefix + "_trace_10mK_1um")
	WAVE tr10um  = $(outPrefix + "_trace_10um")
	WAVE tr100um = $(outPrefix + "_trace_100um")

	tr860   = M[p][c860mK]
	tr10mK  = M[p][c10mK1]
	tr10um  = M[p][c10um]
	tr100um = M[p][c100um]

	SetScale/P x, xOff, xDel, xUnit, tr860
	SetScale/P x, xOff, xDel, xUnit, tr10mK
	SetScale/P x, xOff, xDel, xUnit, tr10um
	SetScale/P x, xOff, xDel, xUnit, tr100um

	DoWindow/K $traceWinName
	Display/N=$traceWinName tr860
	AppendToGraph/L/W=$traceWinName tr10mK
	AppendToGraph/R/W=$traceWinName tr10um
	AppendToGraph/R/W=$traceWinName tr100um

	ModifyGraph/W=$traceWinName lsize=2
	ModifyGraph/W=$traceWinName rgb($NameOfWave(tr860))=(0,0,0)
	ModifyGraph/W=$traceWinName rgb($NameOfWave(tr10mK))=(0,29298,45746)
	ModifyGraph/W=$traceWinName rgb($NameOfWave(tr10um))=(54741,24158,0)
	ModifyGraph/W=$traceWinName rgb($NameOfWave(tr100um))=(52428,31097,42919)

	Label/W=$traceWinName left "860 mK & 10 mK trace"
	Label/W=$traceWinName right "phase coherence-length traces"

	if (matrixKind == 0)
		Label/W=$traceWinName bottom "Energy (\\U)"
	else
		Label/W=$traceWinName bottom "B (T)"
	endif

	Legend/W=$traceWinName/C/N=ResolutionTraceLegend/J/F=0/A=RT

	SetDataFolder $savedDF
End

//==============================================================================
// SNS__FindClosestColumn_T_Lphi
//==============================================================================

Function SNS__FindClosestColumn_T_Lphi(T_K, Lphi_um, targetT_K, targetLphi_um, fixedLphiMode)
	WAVE T_K, Lphi_um
	Variable targetT_K, targetLphi_um, fixedLphiMode

	Variable n = numpnts(T_K)
	Variable i
	Variable best = 0
	Variable score, bestScore = Inf

	for (i = 0; i < n; i += 1)
		if (fixedLphiMode)
			score = abs(T_K[i] - targetT_K) + 0.01*abs(Lphi_um[i] - targetLphi_um)
		else
			score = abs(Lphi_um[i] - targetLphi_um) + 10*abs(T_K[i] - targetT_K)
		endif

		if (score < bestScore)
			bestScore = score
			best = i
		endif
	endfor

	return best
End

//==============================================================================
// SNS_DrawTrajectory
//
// Draw one SNS trajectory on the currently active graph.
// Each call creates a NEW trajectory trace.
//
// Rotation center priority:
//   1) rotCenterX / rotCenterY, if explicitly given
//   2) center of rotCenterImage, if given
//   3) center of first image in active graph
//   4) trajectory midpoint
//==============================================================================
Function SNS_DrawTrajectory(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, idx, [rotAngle, rotCenterImage, rotCenterX, rotCenterY, lineSize, r, g, b])
    Wave Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
    Variable idx
    Wave rotCenterImage
    Variable rotAngle, rotCenterX, rotCenterY, lineSize, r, g, b

    String graphName, imageList, imageName
    String trajXName, trajYName
    Variable nCh
    Variable x1, y1, x2, y2
    Variable x0, y0, theta, deg2rad = Pi/180
    Variable xTmp, yTmp
    Variable haveCenterX = 0
    Variable haveCenterY = 0
    Variable nCopy = 0

    if (ParamIsDefault(rotAngle))
        rotAngle = 0
    endif
    if (ParamIsDefault(lineSize))
        lineSize = 1
    endif
    if (ParamIsDefault(r))
        r = 0
    endif
    if (ParamIsDefault(g))
        g = 0
    endif
    if (ParamIsDefault(b))
        b = 0
    endif

    graphName = WinName(0, 1, 1)
    if (strlen(graphName) == 0)
        Abort "SNS_DrawTrajectory: no active graph window."
    endif
    if (WinType(graphName) != 1)
        Abort "SNS_DrawTrajectory: active window is not a graph."
    endif

    nCh = DimSize(Hit1x_List, 0)
    if (idx < 0 || idx >= nCh)
        Abort "SNS_DrawTrajectory: idx out of range."
    endif

    x1 = Hit1x_List[idx]
    y1 = Hit1y_List[idx]
    x2 = Hit2x_List[idx]
    y2 = Hit2y_List[idx]

    if (numtype(x1) != 0 || numtype(y1) != 0 || numtype(x2) != 0 || numtype(y2) != 0)
        Abort "SNS_DrawTrajectory: selected trajectory contains invalid endpoint coordinates."
    endif

    if (abs(rotAngle) > 0)

        if (!ParamIsDefault(rotCenterX))
            x0 = rotCenterX
            haveCenterX = 1
        endif
        if (!ParamIsDefault(rotCenterY))
            y0 = rotCenterY
            haveCenterY = 1
        endif

        if ((!haveCenterX || !haveCenterY) && !ParamIsDefault(rotCenterImage))
            if (!WaveExists(rotCenterImage))
                Abort "SNS_DrawTrajectory: rotCenterImage does not exist."
            endif
            if (!haveCenterX)
                x0 = DimOffset(rotCenterImage, 0) + 0.5 * (DimSize(rotCenterImage, 0) - 1) * DimDelta(rotCenterImage, 0)
                haveCenterX = 1
            endif
            if (!haveCenterY)
                y0 = DimOffset(rotCenterImage, 1) + 0.5 * (DimSize(rotCenterImage, 1) - 1) * DimDelta(rotCenterImage, 1)
                haveCenterY = 1
            endif
        endif

        if (!haveCenterX || !haveCenterY)
            imageList = ImageNameList(graphName, ";")
            imageName = StringFromList(0, imageList, ";")

            if (strlen(imageName) > 0)
                Wave/Z image = ImageNameToWaveRef(graphName, imageName)
                if (WaveExists(image))
                    if (!haveCenterX)
                        x0 = DimOffset(image, 0) + 0.5 * (DimSize(image, 0) - 1) * DimDelta(image, 0)
                        haveCenterX = 1
                    endif
                    if (!haveCenterY)
                        y0 = DimOffset(image, 1) + 0.5 * (DimSize(image, 1) - 1) * DimDelta(image, 1)
                        haveCenterY = 1
                    endif
                endif
            endif
        endif

        if (!haveCenterX)
            x0 = 0.5 * (x1 + x2)
        endif
        if (!haveCenterY)
            y0 = 0.5 * (y1 + y2)
        endif

        theta = rotAngle * deg2rad

        xTmp = x1 - x0
        yTmp = y1 - y0
        x1 = x0 + cos(theta) * xTmp - sin(theta) * yTmp
        y1 = y0 + sin(theta) * xTmp + cos(theta) * yTmp

        xTmp = x2 - x0
        yTmp = y2 - y0
        x2 = x0 + cos(theta) * xTmp - sin(theta) * yTmp
        y2 = y0 + sin(theta) * xTmp + cos(theta) * yTmp
    endif

    do
        sprintf trajXName, "trajX_%d_%d", idx, nCopy
        sprintf trajYName, "trajY_%d_%d", idx, nCopy
        Wave/Z testWave = $trajYName
        if (!WaveExists(testWave))
            break
        endif
        nCopy += 1
    while (1)

    Make/O/D/N=2 $trajXName, $trajYName
    Wave trajX = $trajXName
    Wave trajY = $trajYName

    trajX[0] = x1
    trajY[0] = y1
    trajX[1] = x2
    trajY[1] = y2

    AppendToGraph/W=$graphName trajY vs trajX
    ModifyGraph/W=$graphName lsize($trajYName)=lineSize
    ModifyGraph/W=$graphName lstyle($trajYName)=0
    ModifyGraph/W=$graphName rgb($trajYName)=(r,g,b)

    return 1
End


//==============================================================================
// SNS_DrawTrajectoryLMaxFromFolder
//
// Code Purpose:
//   Draw the previously extracted longest ray trajectory onto the front graph.
//
// Physics Role:
//   Displays the ray selected by SNS_ExtractModesForFolder as the maximum-L
//   trajectory. This is a visualization of already-computed channel geometry;
//   it does not recompute ray tracing or LDOS.
//
// Required Inputs:
//   dfPath : data folder containing v_Longest and Hit*_List_nm channel waves.
//
// Optional Inputs:
//   rotAngle  : display rotation angle [deg] applied to ray endpoints.
//   topoImage : source image used to determine the rotation center.
//   lineSize  : trajectory line thickness.
//   r,g,b     : trajectory color components in Igor 0..65535 range.
//   verbose   : nonzero prints selected ray diagnostics.
//
// Returns:
//   Selected ray index.
//
// Compatibility:
//   Prefers unit-suffixed Hit*_List_nm / L_N_List_nm / W_eff_List_nm waves.
//   Falls back to legacy unsuffixed names for old experiments.
//==============================================================================
Function SNS_DrawTrajectoryLMaxFromFolder(dfPath, [rotAngle, topoImage, lineSize, r, g, b, verbose])
    String dfPath
    Variable rotAngle, lineSize, r, g, b, verbose
    Wave topoImage

    String oldDF = GetDataFolder(1)
    SetDataFolder $dfPath

    NVAR/Z v_Longest
    Wave/Z Hit1x_List_nm
    Wave/Z Hit1y_List_nm
    Wave/Z Hit2x_List_nm
    Wave/Z Hit2y_List_nm
    Wave/Z Hit1x_List
    Wave/Z Hit1y_List
    Wave/Z Hit2x_List
    Wave/Z Hit2y_List
    Wave/Z W_eff_List_nm
    Wave/Z L_N_List_nm

    if (!NVAR_Exists(v_Longest))
        SetDataFolder $oldDF
        Abort "SNS_DrawTrajectoryLMaxFromFolder: v_Longest not found. Run SNS_ExtractModesForFolder first."
    endif

    if (WaveExists(Hit1x_List_nm) && WaveExists(Hit1y_List_nm) && WaveExists(Hit2x_List_nm) && WaveExists(Hit2y_List_nm))
        Wave Hit1x_use = Hit1x_List_nm
        Wave Hit1y_use = Hit1y_List_nm
        Wave Hit2x_use = Hit2x_List_nm
        Wave Hit2y_use = Hit2y_List_nm
    elseif (WaveExists(Hit1x_List) && WaveExists(Hit1y_List) && WaveExists(Hit2x_List) && WaveExists(Hit2y_List))
        Wave Hit1x_use = Hit1x_List
        Wave Hit1y_use = Hit1y_List
        Wave Hit2x_use = Hit2x_List
        Wave Hit2y_use = Hit2y_List
    else
        SetDataFolder $oldDF
        Abort "SNS_DrawTrajectoryLMaxFromFolder: missing Hit*_List_nm channel waves."
    endif

    if (ParamIsDefault(rotAngle))
        rotAngle = 0
    endif
    if (ParamIsDefault(lineSize))
        lineSize = 1
    endif
    if (ParamIsDefault(r))
        r = 65535
    endif
    if (ParamIsDefault(g))
        g = 0
    endif
    if (ParamIsDefault(b))
        b = 0
    endif
    if (ParamIsDefault(verbose))
        verbose = 0
    endif

    Variable idx = round(v_Longest)

    if (verbose && WaveExists(W_eff_List_nm) && WaveExists(L_N_List_nm))
        Print "Drawing longest trajectory:"
        Print "  idx = ", idx
        Print "  W_eff = ", W_eff_List_nm[idx], " nm"
        Print "  L_N   = ", L_N_List_nm[idx], " nm"
    endif

    if (ParamIsDefault(topoImage))
        SNS_DrawTrajectory(Hit1x_use, Hit1y_use, Hit2x_use, Hit2y_use, idx, rotAngle=rotAngle, lineSize=lineSize, r=r, g=g, b=b)
    else
        SNS_DrawTrajectory(Hit1x_use, Hit1y_use, Hit2x_use, Hit2y_use, idx, rotAngle=rotAngle, rotCenterImage=topoImage, lineSize=lineSize, r=r, g=g, b=b)
    endif

    SetDataFolder $oldDF
    return idx
End


//==============================================================================
// SNS_DrawTrajectoryWMaxFromFolder
//
// Code Purpose:
//   Draw the previously extracted W-max / gap-closing ray trajectory onto the
//   front graph.
//
// Physics Role:
//   Displays the ray selected by SNS_ExtractModesForFolder as v_Gap0. This is
//   a visualization of already-computed channel geometry; it does not recompute
//   ray tracing or LDOS.
//
// Required Inputs:
//   dfPath : data folder containing v_Gap0 and Hit*_List_nm channel waves.
//
// Optional Inputs:
//   rotAngle  : display rotation angle [deg] applied to ray endpoints.
//   topoImage : source image used to determine the rotation center.
//   lineSize  : trajectory line thickness.
//   r,g,b     : trajectory color components in Igor 0..65535 range.
//   verbose   : nonzero prints selected ray diagnostics.
//
// Returns:
//   Selected ray index.
//
// Compatibility:
//   Prefers unit-suffixed Hit*_List_nm / L_N_List_nm / W_eff_List_nm waves.
//   Falls back to legacy unsuffixed names for old experiments.
//==============================================================================
Function SNS_DrawTrajectoryWMaxFromFolder(dfPath, [rotAngle, topoImage, lineSize, r, g, b, verbose])
    String dfPath
    Variable rotAngle, lineSize, r, g, b, verbose
    Wave topoImage

    String oldDF = GetDataFolder(1)
    SetDataFolder $dfPath

    NVAR/Z v_Gap0
    Wave/Z Hit1x_List_nm
    Wave/Z Hit1y_List_nm
    Wave/Z Hit2x_List_nm
    Wave/Z Hit2y_List_nm
    Wave/Z Hit1x_List
    Wave/Z Hit1y_List
    Wave/Z Hit2x_List
    Wave/Z Hit2y_List
    Wave/Z W_eff_List_nm
    Wave/Z L_N_List_nm

    if (!NVAR_Exists(v_Gap0))
        SetDataFolder $oldDF
        Abort "SNS_DrawTrajectoryWMaxFromFolder: v_Gap0 not found. Run SNS_ExtractModesForFolder first."
    endif

    if (WaveExists(Hit1x_List_nm) && WaveExists(Hit1y_List_nm) && WaveExists(Hit2x_List_nm) && WaveExists(Hit2y_List_nm))
        Wave Hit1x_use = Hit1x_List_nm
        Wave Hit1y_use = Hit1y_List_nm
        Wave Hit2x_use = Hit2x_List_nm
        Wave Hit2y_use = Hit2y_List_nm
    elseif (WaveExists(Hit1x_List) && WaveExists(Hit1y_List) && WaveExists(Hit2x_List) && WaveExists(Hit2y_List))
        Wave Hit1x_use = Hit1x_List
        Wave Hit1y_use = Hit1y_List
        Wave Hit2x_use = Hit2x_List
        Wave Hit2y_use = Hit2y_List
    else
        SetDataFolder $oldDF
        Abort "SNS_DrawTrajectoryWMaxFromFolder: missing Hit*_List_nm channel waves."
    endif

    if (ParamIsDefault(rotAngle))
        rotAngle = 0
    endif
    if (ParamIsDefault(lineSize))
        lineSize = 1
    endif
    if (ParamIsDefault(r))
        r = 0
    endif
    if (ParamIsDefault(g))
        g = 0
    endif
    if (ParamIsDefault(b))
        b = 0
    endif
    if (ParamIsDefault(verbose))
        verbose = 0
    endif

    Variable idx = round(v_Gap0)

    if (verbose && WaveExists(W_eff_List_nm) && WaveExists(L_N_List_nm))
        Print "Drawing gap / Wmax trajectory:"
        Print "  idx = ", idx
        Print "  W_eff = ", W_eff_List_nm[idx], " nm"
        Print "  L_N   = ", L_N_List_nm[idx], " nm"
    endif

    if (ParamIsDefault(topoImage))
        SNS_DrawTrajectory(Hit1x_use, Hit1y_use, Hit2x_use, Hit2y_use, idx, rotAngle=rotAngle, lineSize=lineSize, r=r, g=g, b=b)
    else
        SNS_DrawTrajectory(Hit1x_use, Hit1y_use, Hit2x_use, Hit2y_use, idx, rotAngle=rotAngle, rotCenterImage=topoImage, lineSize=lineSize, r=r, g=g, b=b)
    endif

    SetDataFolder $oldDF
    return idx
End

//==============================================================================
// SNS_TagSTS
//
// Purpose:
//   Mark one STS/tip position on the first image in a graph.
//
// Code level:
//   Rotates the supplied x/y coordinate around the displayed image center (or an
//   optional centerTopo wave), appends a cross marker as real graph traces, and
//   optionally draws a text label.
//
// Physics level:
//   The marker identifies the measurement position used by the notebook. The
//   function does not infer or validate the physical meaning of that position.
//
// Inputs:
//   x, y              Position in image-axis units, usually nm.
//
// Optional inputs:
//   rotAngle          Rotation angle in degrees, matching SNS_DuplicateRotatedArea.
//   graphName         Target graph; defaults to the top graph.
//   labelText         Optional Igor formatted text label.
//   centerTopo        Optional wave defining the rotation center.
//   verbose           Print coordinate diagnostics if nonzero.
//   labelPos          Label direction: N, NE, E, SE, S, SW, W, NW.
//   labelDistPct      Label offset as percent of graph width/height.
//   crossSizePt       Cross size control.
//   labelSizePt       Label text size.
//
// Outputs:
//   Creates SNS_TagSTS_crossX_* and SNS_TagSTS_crossY_* marker waves in the
//   current data folder and appends them to graphName.
//
// Provenance:
//   Ported from STMtools/TagSTSpos.ipf TagSTS(...), with SNS_ naming.
//==============================================================================
Function SNS_TagSTS(x, y, [rotAngle, graphName, labelText, centerTopo, verbose, labelPos, labelDistPct, crossSizePt, labelSizePt])
    Variable x, y
    Variable rotAngle
    String graphName, labelText, labelPos
    Wave centerTopo
    Variable verbose
    Variable labelDistPct, crossSizePt, labelSizePt

    String imageList, imageName
    String lp
    String crossXName, crossYName

    Variable nx, ny
    Variable xOff, yOff, xDel, yDel
    Variable xEnd, yEnd
    Variable xLo, xHi, yLo, yHi
    Variable xC, yC
    Variable theta, c, s
    Variable xTmp, yTmp
    Variable xPlot, yPlot

    Variable xAxis0, xAxis1, yAxis0, yAxis1
    Variable dxLabel, dyLabel
    Variable sx, sy
    Variable xLabel, yLabel
    Variable xJust, yJust

    Variable dxImage, dyImage
    Variable dCross
    Variable crossThick
    Variable nCross

    // ---------- defaults ----------
    if (ParamIsDefault(graphName) || strlen(graphName) == 0)
        graphName = WinName(0, 1, 1)
    endif

    if (ParamIsDefault(rotAngle))
        rotAngle = 0
    endif

    if (ParamIsDefault(verbose))
        verbose = 0
    endif

    if (ParamIsDefault(labelDistPct))
        labelDistPct = 4
    endif

    if (ParamIsDefault(crossSizePt))
        crossSizePt = 10
    endif

    if (ParamIsDefault(labelSizePt))
        labelSizePt = 10
    endif

    if (ParamIsDefault(labelPos) || strlen(labelPos) == 0)
        lp = "ne"
    else
        lp = LowerStr(labelPos)
    endif

    if (strlen(graphName) == 0 || WinType(graphName) != 1)
        Abort "SNS_TagSTS: graphName is not an existing graph."
    endif

    // ---------- first image in target graph ----------
    imageList = ImageNameList(graphName, ";")
    imageName = StringFromList(0, imageList, ";")

    if (strlen(imageName) == 0)
        Abort "SNS_TagSTS: graph contains no image."
    endif

    Wave/Z image = ImageNameToWaveRef(graphName, imageName)
    if (!WaveExists(image))
        Abort "SNS_TagSTS: could not resolve image wave."
    endif

    // ---------- choose wave for rotation center ----------
    if (!ParamIsDefault(centerTopo))
        if (!WaveExists(centerTopo))
            Abort "SNS_TagSTS: centerTopo does not exist."
        endif

        nx = DimSize(centerTopo, 0)
        ny = DimSize(centerTopo, 1)

        xOff = DimOffset(centerTopo, 0)
        yOff = DimOffset(centerTopo, 1)
        xDel = DimDelta(centerTopo, 0)
        yDel = DimDelta(centerTopo, 1)
    else
        nx = DimSize(image, 0)
        ny = DimSize(image, 1)

        xOff = DimOffset(image, 0)
        yOff = DimOffset(image, 1)
        xDel = DimDelta(image, 0)
        yDel = DimDelta(image, 1)
    endif

    // ---------- rotation center ----------
    xEnd = xOff + xDel * (nx - 1)
    yEnd = yOff + yDel * (ny - 1)

    xLo = min(xOff, xEnd)
    xHi = max(xOff, xEnd)
    yLo = min(yOff, yEnd)
    yHi = max(yOff, yEnd)

    xC = 0.5 * (xLo + xHi)
    yC = 0.5 * (yLo + yHi)

    // ---------- rotate point ----------
    xPlot = x
    yPlot = y

    if (abs(rotAngle) > 1e-12)
        theta = rotAngle * Pi / 180
        c = cos(theta)
        s = sin(theta)

        xTmp = x - xC
        yTmp = y - yC

        xPlot = xC + c * xTmp - s * yTmp
        yPlot = yC + s * xTmp + c * yTmp
    endif

    // ---------- cross size ----------
    dxImage = abs(DimDelta(image, 0))
    dyImage = abs(DimDelta(image, 1))

    dCross = 4 * min(dxImage, dyImage) * crossSizePt / 10
    if (dCross <= 0)
        dCross = 3 * crossSizePt / 10
    endif

    crossThick = 1.5 * crossSizePt / 10
    if (crossThick < 0.5)
        crossThick = 0.5
    endif

    // ---------- normalize label position ----------
    if (StringMatch(lp, "north"))
        lp = "n"
    elseif (StringMatch(lp, "northeast"))
        lp = "ne"
    elseif (StringMatch(lp, "east"))
        lp = "e"
    elseif (StringMatch(lp, "southeast"))
        lp = "se"
    elseif (StringMatch(lp, "south"))
        lp = "s"
    elseif (StringMatch(lp, "southwest"))
        lp = "sw"
    elseif (StringMatch(lp, "west"))
        lp = "w"
    elseif (StringMatch(lp, "northwest"))
        lp = "nw"
    endif

    sx = 0
    sy = 0

    if (StringMatch(lp, "n"))
        sy = 1
    elseif (StringMatch(lp, "ne"))
        sx = 1
        sy = 1
    elseif (StringMatch(lp, "e"))
        sx = 1
    elseif (StringMatch(lp, "se"))
        sx = 1
        sy = -1
    elseif (StringMatch(lp, "s"))
        sy = -1
    elseif (StringMatch(lp, "sw"))
        sx = -1
        sy = -1
    elseif (StringMatch(lp, "w"))
        sx = -1
    elseif (StringMatch(lp, "nw"))
        sx = -1
        sy = 1
    else
        sx = 1
        sy = 1
        lp = "ne"
    endif

    // ---------- label position ----------
    GetAxis/W=$graphName/Q bottom
    xAxis0 = V_min
    xAxis1 = V_max

    GetAxis/W=$graphName/Q left
    yAxis0 = V_min
    yAxis1 = V_max

    dxLabel = labelDistPct / 100 * abs(xAxis1 - xAxis0)
    dyLabel = labelDistPct / 100 * abs(yAxis1 - yAxis0)

    xLabel = xPlot + sx * dxLabel
    yLabel = yPlot + sy * dyLabel

    // ---------- text justification ----------
    if (sx > 0)
        xJust = 0
    elseif (sx < 0)
        xJust = 2
    else
        xJust = 1
    endif

    if (sy > 0)
        yJust = 0
    elseif (sy < 0)
        yJust = 2
    else
        yJust = 1
    endif

    // ---------- unique cross wave names ----------
    nCross = 0

    do
        sprintf crossXName, "SNS_TagSTS_crossX_%d", nCross
        sprintf crossYName, "SNS_TagSTS_crossY_%d", nCross
        Wave/Z testCross = $crossYName
        if (!WaveExists(testCross))
            break
        endif
        nCross += 1
    while (1)

    // ---------- draw marker as real graph trace ----------
    Make/O/D/N=5 $crossXName, $crossYName

    Wave crossX = $crossXName
    Wave crossY = $crossYName

    crossX[0] = xPlot - dCross
    crossY[0] = yPlot

    crossX[1] = xPlot + dCross
    crossY[1] = yPlot

    crossX[2] = NaN
    crossY[2] = NaN

    crossX[3] = xPlot
    crossY[3] = yPlot - dCross

    crossX[4] = xPlot
    crossY[4] = yPlot + dCross

    AppendToGraph/W=$graphName crossY vs crossX
    ModifyGraph/W=$graphName mode($crossYName)=0
    ModifyGraph/W=$graphName lsize($crossYName)=crossThick
    ModifyGraph/W=$graphName rgb($crossYName)=(0,0,0)

    // ---------- optional label ----------
    if (!ParamIsDefault(labelText) && strlen(labelText) > 0)
        SetDrawLayer/W=$graphName UserFront
        SetDrawEnv/W=$graphName xcoord=bottom, ycoord=left, textrgb=(0,0,0), fsize=labelSizePt, textxjust=xJust, textyjust=yJust
        DrawText/W=$graphName xLabel, yLabel, labelText
    endif

    // ---------- optional debug output ----------
    if (verbose)
        Print "SNS_TagSTS: input x/y =", x, y
        Print "SNS_TagSTS: plotted x/y =", xPlot, yPlot
        Print "SNS_TagSTS: rotation center =", xC, yC
        Print "SNS_TagSTS: rotAngle =", rotAngle
        Print "SNS_TagSTS: label position =", lp
        Print "SNS_TagSTS: label x/y =", xLabel, yLabel
        Print "SNS_TagSTS: dCross =", dCross
    endif
End

//==============================================================================
// Fig-panel annotation helpers
//
// Purpose:
//   Stand-alone graph annotation utilities for analysis and presentation.
//
// Code level:
//   These helpers draw lines, bands, arrows, and text in graph coordinates.
//   They are display-only utilities and do not modify source data waves.
//
// Physics level:
//   Any interpretation of the marked positions is supplied by the caller.
//
// Provenance:
//   Ported from STMtools/DrawInGraph.ipf with SNS_ public names and SNS__
//   private support-helper names.
//==============================================================================
Function SNS__ActivateGraphTarget(graphName)
    String graphName

    Variable hashPos
    String hostName

    if (strlen(graphName) == 0)
        Abort "SNS__ActivateGraphTarget: empty graphName."
    endif

    hashPos = strsearch(graphName, "#", 0)

    if (hashPos >= 0)
        hostName = graphName[0, hashPos - 1]
        DoWindow/F $hostName
        SetActiveSubwindow $graphName
    else
        DoWindow/F $graphName
    endif

    return 0
End


//==============================================================================
// SNS__SetDrawEnvLineDynamic
//
// Dynamic SetDrawEnv helper for line drawing.
// This avoids fragile direct syntax such as xcoord=$xCoord.
//==============================================================================
Function SNS__SetDrawEnvLineDynamic(graphName, xCoord, yCoord, dash, linethick, red, green, blue, opaque, doSave)
    String graphName, xCoord, yCoord
    Variable dash, linethick, red, green, blue, opaque, doSave

    String cmd
    String savePart

    if (doSave)
        savePart = ", save"
    else
        savePart = ""
    endif

    cmd = "SetDrawEnv/W=" + graphName
    cmd += " xcoord= " + xCoord
    cmd += ", ycoord= " + yCoord
    cmd += ", dash= " + num2str(dash)
    cmd += ", linethick= " + num2str(linethick)

    // Keep RGB three-component here for parser robustness.
    // opaque is kept as an argument for backward compatibility.
    cmd += ", linefgc= (" + num2str(red) + "," + num2str(green) + "," + num2str(blue) + ")"
    cmd += savePart

    Execute/Q/Z cmd

    if (V_flag != 0)
        Print "Failed SetDrawEnv command:"
        Print cmd
        Abort "SNS__SetDrawEnvLineDynamic: SetDrawEnv failed."
    endif

    return 0
End


//==============================================================================
// SNS__SetDrawEnvFillDynamic
//
// Dynamic SetDrawEnv helper for filled objects.
//==============================================================================
Function SNS__SetDrawEnvFillDynamic(graphName, xCoord, yCoord, red, green, blue)
    String graphName, xCoord, yCoord
    Variable red, green, blue

    String cmd

    cmd = "SetDrawEnv/W=" + graphName
    cmd += " xcoord= " + xCoord
    cmd += ", ycoord= " + yCoord
    cmd += ", fillpat= 1"
    cmd += ", fillfgc= (" + num2str(red) + "," + num2str(green) + "," + num2str(blue) + ")"
    cmd += ", linefgc= (" + num2str(red) + "," + num2str(green) + "," + num2str(blue) + ")"

    Execute/Q/Z cmd

    if (V_flag != 0)
        Print "Failed SetDrawEnv command:"
        Print cmd
        Abort "SNS__SetDrawEnvFillDynamic: SetDrawEnv failed."
    endif

    return 0
End


//==============================================================================
// SNS__SetDrawEnvTextDynamic
//
// Dynamic SetDrawEnv helper for text drawing with relative coordinates.
//==============================================================================
Function SNS__SetDrawEnvTextDynamic(graphName, xCoord, yCoord, fsize, fontName, xJust, yJust, red, green, blue)
    String graphName, xCoord, yCoord, fontName
    Variable fsize, xJust, yJust, red, green, blue

    String cmd

    cmd = "SetDrawEnv/W=" + graphName
    cmd += " xcoord= " + xCoord
    cmd += ", ycoord= " + yCoord
    cmd += ", fsize= " + num2str(fsize)
    cmd += ", fname=\"" + fontName + "\""
    cmd += ", textxjust= " + num2str(xJust)
    cmd += ", textyjust= " + num2str(yJust)
    cmd += ", textrgb= (" + num2str(red) + "," + num2str(green) + "," + num2str(blue) + ")"

    Execute/Q/Z cmd

    if (V_flag != 0)
        Print "Failed SetDrawEnv command:"
        Print cmd
        Abort "SNS__SetDrawEnvTextDynamic: SetDrawEnv failed."
    endif

    return 0
End


//==============================================================================
// SNS__NormalizeSubwindowTarget
//
// Convert an active subwindow name into a full graph target.
//
// Examples:
//   host="G0", sw="AbsD2"      -> "G0#AbsD2"
//   host="G0", sw="#AbsD2"     -> "G0#AbsD2"
//   host="G0", sw="G0#AbsD2"   -> "G0#AbsD2"
//==============================================================================
Function/S SNS__NormalizeSubwindowTarget(hostName, swName)
    String hostName, swName

    if (strlen(hostName) == 0 || strlen(swName) == 0)
        return ""
    endif

    if (StringMatch(swName, "#*"))
        return hostName + swName
    endif

    if (strsearch(swName, "#", 0) >= 0)
        return swName
    endif

    return hostName + "#" + swName
End


//==============================================================================
// SNS__GraphTargetExists
//
// Returns 1 if graphName is an existing graph or graph subwindow.
// Returns 0 otherwise.
//
// Handles:
//   "Graph0"
//   "Graph0#SubGraph"
//==============================================================================
Function SNS__GraphTargetExists(graphName)
    String graphName

    Variable hashPos
    String hostName
    String subName
    String childList

    if (strlen(graphName) == 0)
        return 0
    endif

    hashPos = strsearch(graphName, "#", 0)

    // ---------- normal top-level graph ----------
    if (hashPos < 0)
        DoWindow $graphName
        return V_flag
    endif

    // ---------- subwindow graph ----------
    hostName = graphName[0, hashPos - 1]
    subName  = graphName[hashPos + 1, strlen(graphName) - 1]

    DoWindow $hostName
    if (!V_flag)
        return 0
    endif

    childList = ChildWindowList(hostName)

    if (WhichListItem(subName, childList, ";", 0, 0) >= 0)
        return 1
    endif

    if (WhichListItem("#" + subName, childList, ";", 0, 0) >= 0)
        return 1
    endif

    if (WhichListItem(graphName, childList, ";", 0, 0) >= 0)
        return 1
    endif

    return 0
End


//==============================================================================
// SNS__GraphHasAxisQuiet
//
// Returns 1 if graphName exists and has axisName.
// Returns 0 otherwise.
//
// Important:
//   Validates the graph/subwindow before calling AxisList().
//==============================================================================
Function SNS__GraphHasAxisQuiet(graphName, axisName)
    String graphName, axisName

    String axList

    if (strlen(graphName) == 0 || strlen(axisName) == 0)
        return 0
    endif

    if (!SNS__GraphTargetExists(graphName))
        return 0
    endif

    axList = AxisList(graphName)

    if (WhichListItem(axisName, axList, ";", 0, 0) >= 0)
        return 1
    endif

    return 0
End


//==============================================================================
// SNS__DefaultGraphForDrawing
//
// Robust default drawing target.
//
// Logic:
//   1. Get front host graph.
//   2. Try active subwindow.
//   3. Try host itself.
//   4. Search child subwindows from newest to oldest.
//
// This preserves old calls like:
//   SNS_DrawVerticalDashedLineAtX(E0HighNeg, dash=1)
// also in HOST/subwindow stacks.
//==============================================================================
Function/S SNS__DefaultGraphForDrawing([xAxis, yAxis])
    String xAxis, yAxis

    String hostName
    String activeSW
    String graphName
    String childList
    String child
    Variable i, n

    if (ParamIsDefault(xAxis) || strlen(xAxis) == 0)
        xAxis = "bottom"
    endif

    if (ParamIsDefault(yAxis) || strlen(yAxis) == 0)
        yAxis = ""
    endif

    hostName = WinName(0, 1)
    if (strlen(hostName) == 0)
        return ""
    endif

    // ---------- 1) try active subwindow ----------
    GetWindow $hostName, activeSW
    activeSW = S_value

    if (strlen(activeSW) > 0)
        graphName = SNS__NormalizeSubwindowTarget(hostName, activeSW)

        if (SNS__GraphHasAxisQuiet(graphName, xAxis))
            if (strlen(yAxis) == 0 || SNS__GraphHasAxisQuiet(graphName, yAxis))
                return graphName
            endif
        endif
    endif

    // ---------- 2) try host itself ----------
    if (SNS__GraphHasAxisQuiet(hostName, xAxis))
        if (strlen(yAxis) == 0 || SNS__GraphHasAxisQuiet(hostName, yAxis))
            return hostName
        endif
    endif

    // ---------- 3) search child subwindows from newest to oldest ----------
    childList = ChildWindowList(hostName)
    n = ItemsInList(childList, ";")

    for (i = n - 1; i >= 0; i -= 1)
        child = StringFromList(i, childList, ";")

        if (strlen(child) == 0)
            continue
        endif

        graphName = SNS__NormalizeSubwindowTarget(hostName, child)

        if (SNS__GraphHasAxisQuiet(graphName, xAxis))
            if (strlen(yAxis) == 0 || SNS__GraphHasAxisQuiet(graphName, yAxis))
                return graphName
            endif
        endif
    endfor

    return hostName
End


//==============================================================================
// SNS__AssertGraphAxisExists
//
// Check whether a standard axis exists on a graph or graph subwindow.
//==============================================================================
Function SNS__AssertGraphAxisExists(graphName, axisName)
    String graphName, axisName

    if (strlen(graphName) == 0)
        Abort "SNS__AssertGraphAxisExists: empty graph name."
    endif

    axisName = LowerStr(axisName)

    if (!StringMatch(axisName, "bottom") && !StringMatch(axisName, "top") && !StringMatch(axisName, "left") && !StringMatch(axisName, "right"))
        Abort "SNS__AssertGraphAxisExists: axis must be 'bottom', 'top', 'left', or 'right'."
    endif

    if (!SNS__GraphTargetExists(graphName))
        Abort "SNS__AssertGraphAxisExists: graph target does not exist: " + graphName
    endif

    if (!SNS__GraphHasAxisQuiet(graphName, axisName))
        Abort "SNS__AssertGraphAxisExists: axis '" + axisName + "' does not exist on graph '" + graphName + "'."
    endif
End

Function SNS_DrawVerticalDashedLineAtX(xVal, [xAxis, yCoord, y0, y1, dash, linethick, red, green, blue, opaque, graphName])
    Variable xVal
    String xAxis, yCoord, graphName
    Variable y0, y1, dash, linethick, red, green, blue, opaque

    Variable xMin, xMax
    Variable xRel

    if (numtype(xVal) != 0)
        Abort "SNS_DrawVerticalDashedLineAtX: xVal must be a finite number."
    endif

    if (ParamIsDefault(xAxis) || strlen(xAxis) == 0)
        xAxis = "bottom"
    endif
    xAxis = LowerStr(xAxis)
    if (!StringMatch(xAxis, "bottom") && !StringMatch(xAxis, "top"))
        xAxis = "bottom"
    endif

    if (ParamIsDefault(yCoord) || strlen(yCoord) == 0)
        yCoord = "plot"
    endif
    yCoord = LowerStr(yCoord)
    if (StringMatch(yCoord, "plot"))
        yCoord = "prel"
    elseif (StringMatch(yCoord, "axis"))
        yCoord = "axrel"
    elseif (StringMatch(yCoord, "window"))
        yCoord = "rel"
    elseif (!StringMatch(yCoord, "prel") && !StringMatch(yCoord, "axrel") && !StringMatch(yCoord, "rel"))
        yCoord = "prel"
    endif

    if (ParamIsDefault(y0))
        y0 = 0
    endif
    if (ParamIsDefault(y1))
        y1 = 1
    endif
    if (ParamIsDefault(dash))
        dash = 3
    endif
    if (ParamIsDefault(linethick))
        linethick = 1
    endif
    if (ParamIsDefault(red))
        red = 0
    endif
    if (ParamIsDefault(green))
        green = 0
    endif
    if (ParamIsDefault(blue))
        blue = 0
    endif
    if (ParamIsDefault(opaque))
        opaque = 65535
    endif

    if (ParamIsDefault(graphName) || strlen(graphName) == 0)
        graphName = SNS__DefaultGraphForDrawing(xAxis=xAxis)
    endif
    if (strlen(graphName) == 0)
        Abort "SNS_DrawVerticalDashedLineAtX: no graph window found."
    endif

    SNS__ActivateGraphTarget(graphName)

    SNS__AssertGraphAxisExists(graphName, xAxis)

    GetAxis/W=$graphName/Q $xAxis
    xMin = V_min
    xMax = V_max

    if (xMax == xMin)
        Abort "SNS_DrawVerticalDashedLineAtX: selected x axis has zero range."
    endif

    xRel = (xVal - xMin) / (xMax - xMin)

    SetDrawLayer/W=$graphName UserFront
    SNS__SetDrawEnvLineDynamic(graphName, "prel", yCoord, dash, linethick, red, green, blue, opaque, 1)

    DrawLine/W=$graphName xRel, y0, xRel, y1
End


//==============================================================================
// SNS_DrawVerticalDashedLinesAtPeaks
//
// Draw vertical lines at x-positions in wPeakX.
//==============================================================================
Function SNS_DrawVerticalDashedLinesAtPeaks(wPeakX, [xAxis, clear, yCoord, y0, y1, dash, linethick, red, green, blue, graphName])
    Wave wPeakX
    String xAxis, yCoord, graphName
    Variable clear, y0, y1, dash, linethick, red, green, blue

    Variable xMin, xMax
    Variable x, xRel
    Variable i, n

    if (ParamIsDefault(xAxis) || strlen(xAxis) == 0)
        xAxis = "bottom"
    endif
    xAxis = LowerStr(xAxis)
    if (!StringMatch(xAxis, "bottom") && !StringMatch(xAxis, "top"))
        xAxis = "bottom"
    endif

    if (ParamIsDefault(clear))
        clear = 1
    endif

    if (ParamIsDefault(yCoord) || strlen(yCoord) == 0)
        yCoord = "plot"
    endif
    yCoord = LowerStr(yCoord)
    if (StringMatch(yCoord, "plot"))
        yCoord = "prel"
    elseif (StringMatch(yCoord, "axis"))
        yCoord = "axrel"
    elseif (StringMatch(yCoord, "window"))
        yCoord = "rel"
    elseif (!StringMatch(yCoord, "prel") && !StringMatch(yCoord, "axrel") && !StringMatch(yCoord, "rel"))
        yCoord = "prel"
    endif

    if (ParamIsDefault(y0))
        y0 = 0
    endif
    if (ParamIsDefault(y1))
        y1 = 1
    endif
    if (ParamIsDefault(dash))
        dash = 3
    endif
    if (ParamIsDefault(linethick))
        linethick = 1
    endif
    if (ParamIsDefault(red))
        red = 0
    endif
    if (ParamIsDefault(green))
        green = 0
    endif
    if (ParamIsDefault(blue))
        blue = 0
    endif

    if (ParamIsDefault(graphName) || strlen(graphName) == 0)
        graphName = SNS__DefaultGraphForDrawing(xAxis=xAxis)
    endif
    if (strlen(graphName) == 0)
        Abort "SNS_DrawVerticalDashedLinesAtPeaks: no graph window found."
    endif

    SNS__ActivateGraphTarget(graphName)

    SNS__AssertGraphAxisExists(graphName, xAxis)

    GetAxis/W=$graphName/Q $xAxis
    xMin = V_min
    xMax = V_max

    if (xMax == xMin)
        Abort "SNS_DrawVerticalDashedLinesAtPeaks: selected x axis has zero range."
    endif

    SetDrawLayer/W=$graphName UserFront

    if (clear)
        DrawAction/W=$graphName delete
    endif

    SNS__SetDrawEnvLineDynamic(graphName, "prel", yCoord, dash, linethick, red, green, blue, 65535, 1)

    n = numpnts(wPeakX)

    for (i = 0; i < n; i += 1)
        x = wPeakX[i]

        if (numtype(x) != 0)
            continue
        endif

        xRel = (x - xMin) / (xMax - xMin)
        DrawLine/W=$graphName xRel, y0, xRel, y1
    endfor
End


//==============================================================================
// SNS_DrawLineArrowOnRotatedTopo
//
// Draw the measured LineSTS direction on a rotated topography. Endpoints are
// read from the retained pos:stspos wave belonging to the LineSTS data folder,
// then rotated about the geometric center of the pre-rotation topography.
//==============================================================================
Function SNS_DrawLineArrowOnRotatedTopo(graphName, wLineSTS, wTopoBeforeRotation, rotAngleDeg, [red, green, blue, lineWidth, dash, arrowMode, arrowLen, arrowFat])
    String graphName
    Wave wLineSTS, wTopoBeforeRotation
    Variable rotAngleDeg, red, green, blue, lineWidth, dash, arrowMode, arrowLen, arrowFat

    if (ParamIsDefault(red))
        red = 0
    endif
    if (ParamIsDefault(green))
        green = 0
    endif
    if (ParamIsDefault(blue))
        blue = 0
    endif
    if (ParamIsDefault(lineWidth))
        lineWidth = 1
    endif
    if (ParamIsDefault(dash))
        dash = 0
    endif
    if (ParamIsDefault(arrowMode))
        arrowMode = 1
    endif

    String posPath = GetWavesDataFolder(wLineSTS, 1) + "pos:stspos"
    Wave/Z positions = $posPath
    if (!WaveExists(positions) || WaveDims(positions) != 2 || DimSize(positions, 0) < 2 || DimSize(positions, 1) < 2)
        Abort "SNS_DrawLineArrowOnRotatedTopo: missing or invalid pos:stspos wave for " + NameOfWave(wLineSTS)
    endif

    Variable x0 = positions[0][0]
    Variable y0 = positions[1][0]
    Variable lastPoint = DimSize(positions, 1) - 1
    Variable x1 = positions[0][lastPoint]
    Variable y1 = positions[1][lastPoint]
    Variable xc = DimOffset(wTopoBeforeRotation, 0) + 0.5 * DimDelta(wTopoBeforeRotation, 0) * (DimSize(wTopoBeforeRotation, 0) - 1)
    Variable yc = DimOffset(wTopoBeforeRotation, 1) + 0.5 * DimDelta(wTopoBeforeRotation, 1) * (DimSize(wTopoBeforeRotation, 1) - 1)
    Variable phi = rotAngleDeg * pi / 180
    Variable xr0 = xc + (x0 - xc) * cos(phi) - (y0 - yc) * sin(phi)
    Variable yr0 = yc + (x0 - xc) * sin(phi) + (y0 - yc) * cos(phi)
    Variable xr1 = xc + (x1 - xc) * cos(phi) - (y1 - yc) * sin(phi)
    Variable yr1 = yc + (x1 - xc) * sin(phi) + (y1 - yc) * cos(phi)

    DoWindow/F $graphName
    if (V_flag == 0)
        Abort "SNS_DrawLineArrowOnRotatedTopo: graph not found: " + graphName
    endif
    if (ParamIsDefault(arrowLen) && ParamIsDefault(arrowFat))
        SetDrawEnv/W=$graphName xcoord=bottom, ycoord=left, linefgc=(red,green,blue), linethick=lineWidth, dash=dash, arrow=arrowMode
    elseif (!ParamIsDefault(arrowLen) && ParamIsDefault(arrowFat))
        SetDrawEnv/W=$graphName xcoord=bottom, ycoord=left, linefgc=(red,green,blue), linethick=lineWidth, dash=dash, arrow=arrowMode, arrowlen=arrowLen
    elseif (ParamIsDefault(arrowLen) && !ParamIsDefault(arrowFat))
        SetDrawEnv/W=$graphName xcoord=bottom, ycoord=left, linefgc=(red,green,blue), linethick=lineWidth, dash=dash, arrow=arrowMode, arrowfat=arrowFat
    else
        SetDrawEnv/W=$graphName xcoord=bottom, ycoord=left, linefgc=(red,green,blue), linethick=lineWidth, dash=dash, arrow=arrowMode, arrowlen=arrowLen, arrowfat=arrowFat
    endif
    DrawLine/W=$graphName xr0, yr0, xr1, yr1
    return 0
End


//==============================================================================


// Draw the acquired GridSTS perimeter as an XY trace on a rotated topography.
// Corners come directly from the retained pos:stspos wave.
Function SNS_AppendGridRotatedTopoTrace(graphName, wGrid, wTopoBeforeRotation, rotAngleDeg, [outPrefix, red, green, blue, lineWidth, dash])
    String graphName, outPrefix
    Wave wGrid, wTopoBeforeRotation
    Variable rotAngleDeg, red, green, blue, lineWidth, dash

    if (ParamIsDefault(outPrefix) || strlen(outPrefix) == 0)
        outPrefix = "SNS_GridRotatedTopoTrace"
    endif
    if (ParamIsDefault(red))
        red = 0
    endif
    if (ParamIsDefault(green))
        green = 0
    endif
    if (ParamIsDefault(blue))
        blue = 0
    endif
    if (ParamIsDefault(lineWidth))
        lineWidth = 1
    endif
    if (ParamIsDefault(dash))
        dash = 0
    endif

    String posPath = GetWavesDataFolder(wGrid, 1) + "pos:stspos"
    Wave/Z positions = $posPath
    Variable nx = DimSize(wGrid, 0)
    Variable ny = DimSize(wGrid, 1)
    if (!WaveExists(positions) || WaveDims(positions) != 2 || DimSize(positions, 0) < 2 || DimSize(positions, 1) != nx*ny)
        Abort "SNS_AppendGridRotatedTopoTrace: missing or invalid pos:stspos wave for " + NameOfWave(wGrid)
    endif

    Make/O/D/N=5 $(outPrefix + "_x"), $(outPrefix + "_y")
    Wave xWave = $(outPrefix + "_x")
    Wave yWave = $(outPrefix + "_y")
    Make/FREE/D/N=4 cornerIndex = {0, nx-1, nx*ny-1, nx*(ny-1)}
    Variable xc = DimOffset(wTopoBeforeRotation, 0) + 0.5*DimDelta(wTopoBeforeRotation, 0)*(DimSize(wTopoBeforeRotation, 0)-1)
    Variable yc = DimOffset(wTopoBeforeRotation, 1) + 0.5*DimDelta(wTopoBeforeRotation, 1)*(DimSize(wTopoBeforeRotation, 1)-1)
    Variable phi = rotAngleDeg*pi/180
    Variable i, sourceIndex, dx, dy
    for (i = 0; i < 4; i += 1)
        sourceIndex = cornerIndex[i]
        dx = positions[0][sourceIndex] - xc
        dy = positions[1][sourceIndex] - yc
        xWave[i] = xc + dx*cos(phi) - dy*sin(phi)
        yWave[i] = yc + dx*sin(phi) + dy*cos(phi)
    endfor
    xWave[4] = xWave[0]
    yWave[4] = yWave[0]
    SetScale d, 0, 0, "nm", xWave
    SetScale d, 0, 0, "nm", yWave

    AppendToGraph/W=$graphName yWave vs xWave
    ModifyGraph/W=$graphName mode($NameOfWave(yWave))=0, rgb($NameOfWave(yWave))=(red,green,blue), lsize($NameOfWave(yWave))=lineWidth, lstyle($NameOfWave(yWave))=dash
    return 0
End


//==============================================================================


//------------------------------------------------------------------------------
// SNS_AppendLayoutGraphByPlotOrigin
//
// Append an absolute-sized graph to a layout so that the top-left corner of
// its plot area (the rectangle bounded by the axes) is at the requested page
// coordinate. The source graph retains its own absolute plot width and height.
// Axis labels, annotations, and explicit margins are carried around that fixed
// plot area. The appended graph object is transparent, so its white window
// background does not create a separate white rectangle in the layout.
//
// Coordinates are layout-page points, measured from the top-left of the page.
// The source graph must already have its final margins and decorations.
//------------------------------------------------------------------------------
Function SNS_AppendLayoutGraphByPlotOrigin(layoutName, graphName, plotLeft, plotTop)
    String layoutName, graphName
    Variable plotLeft, plotTop

    if (WinType(layoutName) != 3)
        Abort "SNS_AppendLayoutGraphByPlotOrigin: layout window not found."
    endif
    if (WinType(graphName) != 1)
        Abort "SNS_AppendLayoutGraphByPlotOrigin: graph window not found."
    endif

    DoUpdate/W=$graphName

    GetWindow $graphName gsize
    Variable graphLeft = V_left
    Variable graphTop = V_top
    Variable graphRight = V_right
    Variable graphBottom = V_bottom

    GetWindow $graphName psize
    Variable sourcePlotLeft = V_left
    Variable sourcePlotTop = V_top
    Variable sourcePlotRight = V_right
    Variable sourcePlotBottom = V_bottom

    if ((sourcePlotRight <= sourcePlotLeft) || (sourcePlotBottom <= sourcePlotTop))
        Abort "SNS_AppendLayoutGraphByPlotOrigin: source graph has an invalid plot area."
    endif

    Variable objectLeft = plotLeft - (sourcePlotLeft - graphLeft)
    Variable objectTop = plotTop - (sourcePlotTop - graphTop)
    Variable objectRight = objectLeft + (graphRight - graphLeft)
    Variable objectBottom = objectTop + (graphBottom - graphTop)

    AppendLayoutObject/F=0/T=1/W=$layoutName/R=(objectLeft,objectTop,objectRight,objectBottom) graph $graphName
    return 0
End


//==============================================================================
// SNS_DrawLayoutPanelLabel
//
// Draw one bold panel letter centered on a translucent gray box in a layout.
// xPt and yPt are page coordinates in points; dxPt and dyPt are optional
// point offsets. Uses root:globalLabelSize when it exists, otherwise 8 pt.
//==============================================================================
Function SNS_DrawLayoutPanelLabel(layoutName, xPt, yPt, labelText, [dxPt, dyPt])
    String layoutName, labelText
    Variable xPt, yPt, dxPt, dyPt

    if (WinType(layoutName) != 3)
        Abort "SNS_DrawLayoutPanelLabel: layout window not found."
    endif

    NVAR/Z gLabelSize = root:globalLabelSize

    Variable fSizePt
    if (NVAR_Exists(gLabelSize))
        fSizePt = gLabelSize
    else
        fSizePt = 8
    endif

    if (ParamIsDefault(dxPt))
        dxPt = 0
    endif
    if (ParamIsDefault(dyPt))
        dyPt = 0
    endif

    Variable x0 = xPt + dxPt
    Variable y0 = yPt + dyPt
    Variable boxWPt = fSizePt - 1
    Variable boxHPt = fSizePt + 2
    Variable boxGray = 61000
    Variable boxAlpha = round(0.60 * 65535)

    SetDrawLayer/W=$layoutName UserFront
    SetDrawEnv/W=$layoutName fillpat=1,fillfgc=(boxGray,boxGray,boxGray,boxAlpha),linefgc=(boxGray,boxGray,boxGray),linethick=0
    DrawRect/W=$layoutName x0-boxWPt/2,y0-boxHPt/2,x0+boxWPt/2,y0+boxHPt/2

    SetDrawEnv/W=$layoutName textrgb=(0,0,0),fsize=fSizePt,fstyle=1,textxjust=1,textyjust=1
    DrawText/W=$layoutName x0,y0,labelText
    return 0
End


//==============================================================================
// SNS_DrawHorizontalDashedLineAtY
//
// Draw a horizontal dashed line at y = yVal.
//==============================================================================
Function SNS_DrawHorizontalDashedLineAtY(yVal, [yAxis, xCoord, x0, x1, dash, linethick, red, green, blue, opaque, graphName])
    Variable yVal
    String yAxis, xCoord, graphName
    Variable x0, x1, dash, linethick, red, green, blue, opaque

    Variable yMin, yMax
    Variable yRel

    if (numtype(yVal) != 0)
        Abort "SNS_DrawHorizontalDashedLineAtY: yVal must be a finite number."
    endif

    if (ParamIsDefault(yAxis) || strlen(yAxis) == 0)
        yAxis = "left"
    endif
    yAxis = LowerStr(yAxis)
    if (!StringMatch(yAxis, "left") && !StringMatch(yAxis, "right"))
        yAxis = "left"
    endif

    if (ParamIsDefault(xCoord) || strlen(xCoord) == 0)
        xCoord = "plot"
    endif
    xCoord = LowerStr(xCoord)
    if (StringMatch(xCoord, "plot"))
        xCoord = "prel"
    elseif (StringMatch(xCoord, "axis"))
        xCoord = "axrel"
    elseif (StringMatch(xCoord, "window"))
        xCoord = "rel"
    elseif (!StringMatch(xCoord, "prel") && !StringMatch(xCoord, "axrel") && !StringMatch(xCoord, "rel"))
        xCoord = "prel"
    endif

    if (ParamIsDefault(x0))
        x0 = 0
    endif
    if (ParamIsDefault(x1))
        x1 = 1
    endif
    if (ParamIsDefault(dash))
        dash = 3
    endif
    if (ParamIsDefault(linethick))
        linethick = 1
    endif
    if (ParamIsDefault(red))
        red = 0
    endif
    if (ParamIsDefault(green))
        green = 0
    endif
    if (ParamIsDefault(blue))
        blue = 0
    endif
    if (ParamIsDefault(opaque))
        opaque = 65535
    endif

    if (ParamIsDefault(graphName) || strlen(graphName) == 0)
        graphName = SNS__DefaultGraphForDrawing(yAxis=yAxis)
    endif
    if (strlen(graphName) == 0)
        Abort "SNS_DrawHorizontalDashedLineAtY: no graph window found."
    endif

    SNS__ActivateGraphTarget(graphName)

    SNS__AssertGraphAxisExists(graphName, yAxis)

    GetAxis/W=$graphName/Q $yAxis
    yMin = V_min
    yMax = V_max

    if (yMax == yMin)
        Abort "SNS_DrawHorizontalDashedLineAtY: selected y axis has zero range."
    endif

    yRel = (yMax - yVal) / (yMax - yMin)

    SetDrawLayer/W=$graphName UserFront
    SNS__SetDrawEnvLineDynamic(graphName, xCoord, "prel", dash, linethick, red, green, blue, opaque, 1)

    DrawLine/W=$graphName x0, yRel, x1, yRel
End


//==============================================================================
// SNS_DrawHorizontalDashedLinesAtPeaks
//
// Draw horizontal lines at y-positions in wPeakY.
//==============================================================================
Function SNS_DrawHorizontalDashedLinesAtPeaks(wPeakY, [yAxis, clear, xCoord, x0, x1, dash, linethick, red, green, blue, graphName])
    Wave wPeakY
    String yAxis, xCoord, graphName
    Variable clear, x0, x1, dash, linethick, red, green, blue

    Variable yMin, yMax
    Variable y, yRel
    Variable i, n

    if (ParamIsDefault(yAxis) || strlen(yAxis) == 0)
        yAxis = "left"
    endif
    yAxis = LowerStr(yAxis)
    if (!StringMatch(yAxis, "left") && !StringMatch(yAxis, "right"))
        yAxis = "left"
    endif

    if (ParamIsDefault(clear))
        clear = 1
    endif

    if (ParamIsDefault(xCoord) || strlen(xCoord) == 0)
        xCoord = "plot"
    endif
    xCoord = LowerStr(xCoord)
    if (StringMatch(xCoord, "plot"))
        xCoord = "prel"
    elseif (StringMatch(xCoord, "axis"))
        xCoord = "axrel"
    elseif (StringMatch(xCoord, "window"))
        xCoord = "rel"
    elseif (!StringMatch(xCoord, "prel") && !StringMatch(xCoord, "axrel") && !StringMatch(xCoord, "rel"))
        xCoord = "prel"
    endif

    if (ParamIsDefault(x0))
        x0 = 0
    endif
    if (ParamIsDefault(x1))
        x1 = 1
    endif
    if (ParamIsDefault(dash))
        dash = 3
    endif
    if (ParamIsDefault(linethick))
        linethick = 1
    endif
    if (ParamIsDefault(red))
        red = 0
    endif
    if (ParamIsDefault(green))
        green = 0
    endif
    if (ParamIsDefault(blue))
        blue = 0
    endif

    if (ParamIsDefault(graphName) || strlen(graphName) == 0)
        graphName = SNS__DefaultGraphForDrawing(yAxis=yAxis)
    endif
    if (strlen(graphName) == 0)
        Abort "SNS_DrawHorizontalDashedLinesAtPeaks: no graph window found."
    endif

    SNS__ActivateGraphTarget(graphName)

    SNS__AssertGraphAxisExists(graphName, yAxis)

    GetAxis/W=$graphName/Q $yAxis
    yMin = V_min
    yMax = V_max

    if (yMax == yMin)
        Abort "SNS_DrawHorizontalDashedLinesAtPeaks: selected y axis has zero range."
    endif

    SetDrawLayer/W=$graphName UserFront

    if (clear)
        DrawAction/W=$graphName delete
    endif

    SNS__SetDrawEnvLineDynamic(graphName, xCoord, "prel", dash, linethick, red, green, blue, 65535, 1)

    n = numpnts(wPeakY)

    for (i = 0; i < n; i += 1)
        y = wPeakY[i]

        if (numtype(y) != 0)
            continue
        endif

        yRel = (yMax - y) / (yMax - yMin)
        DrawLine/W=$graphName x0, yRel, x1, yRel
    endfor
End


//==============================================================================
// SNS_DrawVerticalBand
//
// Draw a filled vertical band between x0 and x1.
// x-values are interpreted in the selected x-axis units.
// Use xAxis="bottom" or xAxis="top".
//==============================================================================
Function SNS_DrawVerticalBand(x0, x1, [graphName, xAxis, layer, yCoord, y0, y1, red, green, blue])
    Variable x0, x1
    String graphName, xAxis, layer, yCoord
    Variable y0, y1, red, green, blue

    Variable xMin, xMax
    Variable x0Rel, x1Rel, tmp

    if (ParamIsDefault(xAxis) || strlen(xAxis) == 0)
        xAxis = "bottom"
    endif

    if (ParamIsDefault(layer) || strlen(layer) == 0)
        layer = "UserBack"
    endif

    if (ParamIsDefault(yCoord) || strlen(yCoord) == 0)
        yCoord = "plot"
    endif
    yCoord = LowerStr(yCoord)
    if (StringMatch(yCoord, "plot"))
        yCoord = "prel"
    elseif (StringMatch(yCoord, "axis"))
        yCoord = "axrel"
    elseif (StringMatch(yCoord, "window"))
        yCoord = "rel"
    elseif (!StringMatch(yCoord, "prel") && !StringMatch(yCoord, "axrel") && !StringMatch(yCoord, "rel"))
        yCoord = "prel"
    endif

    if (ParamIsDefault(y0))
        y0 = 0
    endif
    if (ParamIsDefault(y1))
        y1 = 1
    endif
    if (ParamIsDefault(red))
        red = 61166
    endif
    if (ParamIsDefault(green))
        green = 61166
    endif
    if (ParamIsDefault(blue))
        blue = 61166
    endif

    if (ParamIsDefault(graphName) || strlen(graphName) == 0)
        graphName = SNS__DefaultGraphForDrawing(xAxis=xAxis)
    endif
    if (strlen(graphName) == 0)
        Abort "SNS_DrawVerticalBand: no graph window found."
    endif

    SNS__ActivateGraphTarget(graphName)

    SNS__AssertGraphAxisExists(graphName, xAxis)

    GetAxis/W=$graphName/Q $xAxis
    xMin = V_min
    xMax = V_max

    if (xMax == xMin)
        Abort "SNS_DrawVerticalBand: selected x-axis has zero range."
    endif

    // Convert selected-axis coordinates to plot-relative coordinates.
    // This also works for reversed axes.
    x0Rel = (x0 - xMin) / (xMax - xMin)
    x1Rel = (x1 - xMin) / (xMax - xMin)

    if (x1Rel < x0Rel)
        tmp = x0Rel
        x0Rel = x1Rel
        x1Rel = tmp
    endif

    SetDrawLayer/W=$graphName $layer
    SNS__SetDrawEnvFillDynamic(graphName, "prel", yCoord, red, green, blue)

    DrawRect/W=$graphName x0Rel, y0, x1Rel, y1
End



//==============================================================================
// SNS_DrawHorizontalArrow
//
// Draw a horizontal arrow using axis coordinates.
// Internally converts to plot-relative coordinates.
//==============================================================================
Function SNS_DrawHorizontalArrow(x0, x1, y, [graphName, layer, xAxis, yAxis, arrowMode, arrowLen, dash, lineThick, red, green, blue])
    Variable x0, x1, y
    String graphName, layer, xAxis, yAxis
    Variable arrowMode, arrowLen, dash, lineThick, red, green, blue

    Variable xMin, xMax, yMin, yMax
    Variable x0r, x1r, yr

    if (ParamIsDefault(layer) || strlen(layer) == 0)
        layer = "UserFront"
    endif

    if (ParamIsDefault(xAxis) || strlen(xAxis) == 0)
        xAxis = "bottom"
    endif
    xAxis = LowerStr(xAxis)
    if (!StringMatch(xAxis, "bottom") && !StringMatch(xAxis, "top"))
        Abort "SNS_DrawHorizontalArrow: xAxis must be 'bottom' or 'top'."
    endif

    if (ParamIsDefault(yAxis) || strlen(yAxis) == 0)
        yAxis = "left"
    endif
    yAxis = LowerStr(yAxis)
    if (!StringMatch(yAxis, "left") && !StringMatch(yAxis, "right"))
        Abort "SNS_DrawHorizontalArrow: yAxis must be 'left' or 'right'."
    endif

    if (ParamIsDefault(graphName) || strlen(graphName) == 0)
        graphName = SNS__DefaultGraphForDrawing(xAxis=xAxis, yAxis=yAxis)
    endif
    if (strlen(graphName) == 0)
        Abort "SNS_DrawHorizontalArrow: no graph window found."
    endif

    if (ParamIsDefault(arrowMode))
        arrowMode = 3
    endif
    if (ParamIsDefault(arrowLen))
        arrowLen = 8
    endif
    if (ParamIsDefault(dash))
        dash = 0
    endif
    if (ParamIsDefault(lineThick))
        lineThick = 1
    endif
    if (ParamIsDefault(red))
        red = 0
    endif
    if (ParamIsDefault(green))
        green = 0
    endif
    if (ParamIsDefault(blue))
        blue = 0
    endif

    SNS__ActivateGraphTarget(graphName)

    SNS__AssertGraphAxisExists(graphName, xAxis)
    SNS__AssertGraphAxisExists(graphName, yAxis)

    GetAxis/W=$graphName/Q $xAxis
    xMin = V_min
    xMax = V_max
    if (xMax == xMin)
        Abort "SNS_DrawHorizontalArrow: selected x axis has zero range."
    endif

    GetAxis/W=$graphName/Q $yAxis
    yMin = V_min
    yMax = V_max
    if (yMax == yMin)
        Abort "SNS_DrawHorizontalArrow: selected y axis has zero range."
    endif

    x0r = (x0 - xMin) / (xMax - xMin)
    x1r = (x1 - xMin) / (xMax - xMin)
    yr  = (yMax - y) / (yMax - yMin)

    SetDrawLayer/W=$graphName $layer
    SetDrawEnv/W=$graphName xcoord=prel, ycoord=prel, arrow=arrowMode, arrowlen=arrowLen, dash=dash, linethick=lineThick, linefgc=(red,green,blue)

    DrawLine/W=$graphName x0r, yr, x1r, yr
End


//==============================================================================
// SNS_DrawVerticalArrow
//
// Draw a vertical arrow using axis coordinates.
// Internally converts to plot-relative coordinates.
//==============================================================================
Function SNS_DrawVerticalArrow(x, y0, y1, [graphName, layer, xAxis, yAxis, arrowMode, arrowLen, dash, lineThick, red, green, blue])
    Variable x, y0, y1
    String graphName, layer, xAxis, yAxis
    Variable arrowMode, arrowLen, dash, lineThick, red, green, blue

    Variable xMin, xMax, yMin, yMax
    Variable xr, y0r, y1r

    if (ParamIsDefault(layer) || strlen(layer) == 0)
        layer = "UserFront"
    endif

    if (ParamIsDefault(xAxis) || strlen(xAxis) == 0)
        xAxis = "bottom"
    endif
    xAxis = LowerStr(xAxis)
    if (!StringMatch(xAxis, "bottom") && !StringMatch(xAxis, "top"))
        Abort "SNS_DrawVerticalArrow: xAxis must be 'bottom' or 'top'."
    endif

    if (ParamIsDefault(yAxis) || strlen(yAxis) == 0)
        yAxis = "left"
    endif
    yAxis = LowerStr(yAxis)
    if (!StringMatch(yAxis, "left") && !StringMatch(yAxis, "right"))
        Abort "SNS_DrawVerticalArrow: yAxis must be 'left' or 'right'."
    endif

    if (ParamIsDefault(graphName) || strlen(graphName) == 0)
        graphName = SNS__DefaultGraphForDrawing(xAxis=xAxis, yAxis=yAxis)
    endif
    if (strlen(graphName) == 0)
        Abort "SNS_DrawVerticalArrow: no graph window found."
    endif

    if (ParamIsDefault(arrowMode))
        arrowMode = 3
    endif
    if (ParamIsDefault(arrowLen))
        arrowLen = 8
    endif
    if (ParamIsDefault(dash))
        dash = 0
    endif
    if (ParamIsDefault(lineThick))
        lineThick = 1
    endif
    if (ParamIsDefault(red))
        red = 0
    endif
    if (ParamIsDefault(green))
        green = 0
    endif
    if (ParamIsDefault(blue))
        blue = 0
    endif

    SNS__ActivateGraphTarget(graphName)

    SNS__AssertGraphAxisExists(graphName, xAxis)
    SNS__AssertGraphAxisExists(graphName, yAxis)

    GetAxis/W=$graphName/Q $xAxis
    xMin = V_min
    xMax = V_max
    if (xMax == xMin)
        Abort "SNS_DrawVerticalArrow: selected x axis has zero range."
    endif

    GetAxis/W=$graphName/Q $yAxis
    yMin = V_min
    yMax = V_max
    if (yMax == yMin)
        Abort "SNS_DrawVerticalArrow: selected y axis has zero range."
    endif

    xr  = (x - xMin) / (xMax - xMin)
    y0r = (yMax - y0) / (yMax - yMin)
    y1r = (yMax - y1) / (yMax - yMin)

    SetDrawLayer/W=$graphName $layer
    SetDrawEnv/W=$graphName xcoord=prel, ycoord=prel, arrow=arrowMode, arrowlen=arrowLen, dash=dash, linethick=lineThick, linefgc=(red,green,blue)

    DrawLine/W=$graphName xr, y0r, xr, y1r
End


Function SNS_DrawGraphText(x, y, txt, [graphName, layer, xCoord, yCoord, fsize, fontName, xJust, yJust, red, green, blue])
    Variable x, y
    String txt
    String graphName, layer, xCoord, yCoord, fontName
    Variable fsize, xJust, yJust, red, green, blue

    Variable xMin, xMax, yMin, yMax
    Variable xr, yr
    Variable useAxisX, useAxisY
    Variable useRelX, useRelY

    if (ParamIsDefault(layer) || strlen(layer) == 0)
        layer = "UserFront"
    endif

    if (ParamIsDefault(xCoord) || strlen(xCoord) == 0)
        xCoord = "bottom"
    endif
    if (ParamIsDefault(yCoord) || strlen(yCoord) == 0)
        yCoord = "left"
    endif

    xCoord = LowerStr(xCoord)
    yCoord = LowerStr(yCoord)

    useAxisX = StringMatch(xCoord, "bottom") || StringMatch(xCoord, "top")
    useAxisY = StringMatch(yCoord, "left")   || StringMatch(yCoord, "right")

    useRelX = StringMatch(xCoord, "prel") || StringMatch(xCoord, "axrel") || StringMatch(xCoord, "rel")
    useRelY = StringMatch(yCoord, "prel") || StringMatch(yCoord, "axrel") || StringMatch(yCoord, "rel")

    if (ParamIsDefault(graphName) || strlen(graphName) == 0)
        if (useAxisX && useAxisY)
            graphName = SNS__DefaultGraphForDrawing(xAxis=xCoord, yAxis=yCoord)
        else
            graphName = SNS__DefaultGraphForDrawing()
        endif
    endif
    if (strlen(graphName) == 0)
        Abort "SNS_DrawGraphText: no graph window found."
    endif

    if (ParamIsDefault(fsize))
        fsize = 10
    endif
    if (ParamIsDefault(fontName) || strlen(fontName) == 0)
        fontName = "Arial"
    endif
    if (ParamIsDefault(xJust))
        xJust = 1
    endif
    if (ParamIsDefault(yJust))
        yJust = 1
    endif
    if (ParamIsDefault(red))
        red = 0
    endif
    if (ParamIsDefault(green))
        green = 0
    endif
    if (ParamIsDefault(blue))
        blue = 0
    endif

    SNS__ActivateGraphTarget(graphName)

    SetDrawLayer/W=$graphName $layer

    if (useAxisX && useAxisY)

        SNS__AssertGraphAxisExists(graphName, xCoord)
        SNS__AssertGraphAxisExists(graphName, yCoord)

        GetAxis/W=$graphName/Q $xCoord
        xMin = V_min
        xMax = V_max
        if (xMax == xMin)
            Abort "SNS_DrawGraphText: selected x axis has zero range."
        endif

        GetAxis/W=$graphName/Q $yCoord
        yMin = V_min
        yMax = V_max
        if (yMax == yMin)
            Abort "SNS_DrawGraphText: selected y axis has zero range."
        endif

        xr = (x - xMin) / (xMax - xMin)
        yr = (yMax - y) / (yMax - yMin)

        SNS__SetDrawEnvTextDynamic(graphName, "prel", "prel", fsize, fontName, xJust, yJust, red, green, blue)
        DrawText/W=$graphName xr, yr, txt
        return 0
    endif

    if (useRelX && useRelY)
        SNS__SetDrawEnvTextDynamic(graphName, xCoord, yCoord, fsize, fontName, xJust, yJust, red, green, blue)
        DrawText/W=$graphName x, y, txt
        return 0
    endif

    Abort "SNS_DrawGraphText: use either axis coordinates for both axes or relative coordinates for both axes."
End


//==============================================================================
// SNS_FormatValueWithUncertainty / SNS_FormatEnergyMeVWithUncertainty
//
// Format fitted values for graph annotations using one significant uncertainty
// digit. The uncertainty is rounded upward and the value is rounded to the same
// decimal place. The energy wrapper switches small meV values to microelectron-
// volts when both |value| and uncertainty are below 0.01 meV.
//==============================================================================
Function/S SNS_FormatValueWithUncertainty(value, uncertainty, unitText)
    Variable value, uncertainty
    String unitText

    if (numtype(value) != 0 || numtype(uncertainty) != 0 || uncertainty <= 0)
        return "not available"
    endif

    Variable uncertaintyAbs = abs(uncertainty)
    Variable exponent10 = floor(log(uncertaintyAbs))
    Variable placeValue = 10^exponent10
    Variable uncertaintyDigit = ceil(uncertaintyAbs/placeValue - 1e-12)

    if (uncertaintyDigit >= 10)
        exponent10 += 1
        placeValue *= 10
        uncertaintyDigit = 1
    endif

    Variable roundedValue = round(value/placeValue)*placeValue
    Variable roundedUncertainty = uncertaintyDigit*placeValue
    Variable decimalPlaces = max(0, -exponent10)
    String numberFormat, valueText, uncertaintyText, result

    sprintf numberFormat, "%%.%df", decimalPlaces
    sprintf valueText, numberFormat, roundedValue
    sprintf uncertaintyText, numberFormat, roundedUncertainty

    result = valueText + " ± " + uncertaintyText
    if (strlen(unitText) > 0)
        result += " " + unitText
    endif
    return result
End


Function/S SNS_FormatEnergyMeVWithUncertainty(value_meV, uncertainty_meV)
    Variable value_meV, uncertainty_meV

    if (numtype(value_meV) == 0 && numtype(uncertainty_meV) == 0)
        if (abs(value_meV) < 0.01 && abs(uncertainty_meV) < 0.01)
            return SNS_FormatValueWithUncertainty(1e3*value_meV, 1e3*uncertainty_meV, "µeV")
        endif
    endif
    return SNS_FormatValueWithUncertainty(value_meV, uncertainty_meV, "meV")
End
