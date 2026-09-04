#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma DefaultTab={3,20,4}

//==============================================================================
// SNS_IslandSummary
//
// Purpose:
//   Maintain the cross-island scalar results used by the main-text island
//   comparison. Rows are keyed by DataName, so rerunning an island/set updates
//   its existing row rather than adding a duplicate. The finalized
//   supplemental-summary CSV files can also be validated and aggregated into
//   the same table with SNS_ImportIslandSummaryCSVs, or with
//   SNS_ImportIslandSummaryCSVsFromPaths when the files reside in separate
//   figure folders.
//
// Data location:
//   root:IslandSummary
//
// Notes:
//   Geometry uncertainties default to 10% of the reported length scale.
//   E0Avg_err_meV is the half-range between the supplied E0 bounds.
//   Fig. 2 ODR preparation, repeated-length grouping, and B0(W)/E0(L) fit
//   models are implemented here with SNS_IS_* names; STMtools is not required.
//==============================================================================


//------------------------------------------------------------------------------
// Ensure one numeric summary column exists and matches the row count.
//------------------------------------------------------------------------------
Static Function SNS_IS_EnsureNumCol(colName, nRows)
    String colName
    Variable nRows

    String fullPath = "root:IslandSummary:" + colName
    Wave/Z w = $fullPath

    if (!WaveExists(w))
        Make/D/O/N=(nRows) $fullPath
        Wave wNew = $fullPath
        wNew = NaN
        return 0
    endif

    Variable oldN = numpnts(w)
    if (oldN != nRows)
        Redimension/N=(nRows) w
        if (nRows > oldN)
            w[oldN, nRows-1] = NaN
        endif
    endif

    return 0
End


//==============================================================================
// SNS_InitIslandSummary
//
// Purpose:
//   Create or repair the complete summary-table schema without deleting rows.
//==============================================================================
Function SNS_InitIslandSummary()
    NewDataFolder/O root:IslandSummary

    Wave/Z/T wDataName = root:IslandSummary:DataName
    if (!WaveExists(wDataName))
        Make/T/O/N=0 root:IslandSummary:DataName
    endif

    Wave/T wDataName2 = root:IslandSummary:DataName
    Variable nRows = numpnts(wDataName2)

    SNS_IS_EnsureNumCol("h_Cu", nRows)
    SNS_IS_EnsureNumCol("Area_nm2", nRows)
    SNS_IS_EnsureNumCol("Perim_nm", nRows)
    SNS_IS_EnsureNumCol("PoverA_nmInv", nRows)
    SNS_IS_EnsureNumCol("AoverP_nm", nRows)
    SNS_IS_EnsureNumCol("B0_mT", nRows)
    SNS_IS_EnsureNumCol("B0_err_mT", nRows)
    SNS_IS_EnsureNumCol("E0Low_meV", nRows)
    SNS_IS_EnsureNumCol("E0High_meV", nRows)
    SNS_IS_EnsureNumCol("E0Avg_meV", nRows)
    SNS_IS_EnsureNumCol("E0Avg_err_meV", nRows)
    SNS_IS_EnsureNumCol("GlobalMaxL_nm", nRows)
    SNS_IS_EnsureNumCol("GlobalMaxL_err_nm", nRows)
    SNS_IS_EnsureNumCol("ClusterMaxL_nm", nRows)
    SNS_IS_EnsureNumCol("ClusterMaxL_err_nm", nRows)
    SNS_IS_EnsureNumCol("ClusterMedianL_nm", nRows)
    SNS_IS_EnsureNumCol("ClusterMedianL_err_nm", nRows)
    SNS_IS_EnsureNumCol("GlobalMaxW_nm", nRows)
    SNS_IS_EnsureNumCol("GlobalMaxW_err_nm", nRows)
    SNS_IS_EnsureNumCol("ClusterMaxW_nm", nRows)
    SNS_IS_EnsureNumCol("ClusterMaxW_err_nm", nRows)
    SNS_IS_EnsureNumCol("ClusterMedianW_nm", nRows)
    SNS_IS_EnsureNumCol("ClusterMedianW_err_nm", nRows)
End


//------------------------------------------------------------------------------
// Import one validated supplemental-summary CSV into root:IslandSummary.
// Called by SNS_ImportIslandSummaryCSVs after the destination table is cleared.
//------------------------------------------------------------------------------
Static Function SNS_IS_ImportOneSummaryCSV(pathName, fileName)
    String pathName, fileName

    String numericCols = "h_Cu;Area_nm2;Perim_nm;PoverA_nmInv;AoverP_nm;B0_mT;B0_err_mT;E0Low_meV;E0High_meV;E0Avg_meV;E0Avg_err_meV;GlobalMaxL_nm;GlobalMaxL_err_nm;ClusterMaxL_nm;ClusterMaxL_err_nm;ClusterMedianL_nm;ClusterMedianL_err_nm;GlobalMaxW_nm;GlobalMaxW_err_nm;ClusterMaxW_nm;ClusterMaxW_err_nm;ClusterMedianW_nm;ClusterMedianW_err_nm"
    Variable nNumericCols = ItemsInList(numericCols)
    Variable i, j

    DFREF oldDFR = GetDataFolderDFR()
    KillDataFolder/Z root:SNSIslandSummaryImportTemp
    NewDataFolder/O root:SNSIslandSummaryImportTemp
    SetDataFolder root:SNSIslandSummaryImportTemp

    LoadWave/J/D/W/K=0/A/O/P=$pathName fileName
    Variable nLoadedWaves = V_flag
    SetDataFolder oldDFR

    if (nLoadedWaves != nNumericCols + 1)
        Abort "SNS_ImportIslandSummaryCSVs: expected 24 columns in " + fileName
    endif

    Wave/Z/T inDataName = root:SNSIslandSummaryImportTemp:DataName
    if (!WaveExists(inDataName))
        Abort "SNS_ImportIslandSummaryCSVs: missing text column DataName in " + fileName
    endif

    Variable nRows = numpnts(inDataName)
    if (nRows < 1)
        Abort "SNS_ImportIslandSummaryCSVs: no data rows in " + fileName
    endif

    for (i = 0; i < nNumericCols; i += 1)
        String colName = StringFromList(i, numericCols)
        Wave/Z inCol = $("root:SNSIslandSummaryImportTemp:" + colName)
        if (!WaveExists(inCol))
            Abort "SNS_ImportIslandSummaryCSVs: missing numeric column " + colName + " in " + fileName
        endif
        if (numpnts(inCol) != nRows)
            Abort "SNS_ImportIslandSummaryCSVs: column-length mismatch for " + colName + " in " + fileName
        endif
    endfor

    SNS_InitIslandSummary()
    Wave/T outDataName = root:IslandSummary:DataName
    Variable oldN = numpnts(outDataName)

    for (i = 0; i < nRows; i += 1)
        if (strlen(inDataName[i]) == 0)
            Abort "SNS_ImportIslandSummaryCSVs: blank DataName in " + fileName
        endif
        for (j = 0; j < oldN; j += 1)
            if (CmpStr(inDataName[i], outDataName[j]) == 0)
                Abort "SNS_ImportIslandSummaryCSVs: duplicate DataName " + inDataName[i]
            endif
        endfor
        for (j = i + 1; j < nRows; j += 1)
            if (CmpStr(inDataName[i], inDataName[j]) == 0)
                Abort "SNS_ImportIslandSummaryCSVs: duplicate DataName " + inDataName[i]
            endif
        endfor
    endfor

    InsertPoints oldN, nRows, outDataName
    outDataName[oldN, oldN+nRows-1] = inDataName[p-oldN]

    for (i = 0; i < nNumericCols; i += 1)
        colName = StringFromList(i, numericCols)
        Wave sourceCol = $("root:SNSIslandSummaryImportTemp:" + colName)
        Wave targetCol = $("root:IslandSummary:" + colName)
        InsertPoints oldN, nRows, targetCol
        targetCol[oldN, oldN+nRows-1] = sourceCol[p-oldN]
    endfor

    return nRows
End


//==============================================================================
// SNS_ImportIslandSummaryCSVs
//
// Purpose:
//   Replace root:IslandSummary with the rows from a semicolon-separated list
//   of supplemental-summary CSV files. CSV headers must match the complete
//   SNS_InitIslandSummary schema. Duplicate DataName values are rejected.
//
// Inputs:
//   pathName     : name of an existing Igor symbolic path
//   fileList     : semicolon-separated CSV filenames
//   expectedRows : required final row count; use a negative value to skip
//                  the final count check
//==============================================================================
Function SNS_ImportIslandSummaryCSVs(pathName, fileList, expectedRows)
    String pathName, fileList
    Variable expectedRows

    Variable nFiles = ItemsInList(fileList)
    if (nFiles < 1)
        Abort "SNS_ImportIslandSummaryCSVs: fileList is empty."
    endif

    SNS_ClearIslandSummary()
    SNS_InitIslandSummary()

    Variable i, totalRows = 0
    for (i = 0; i < nFiles; i += 1)
        String fileName = StringFromList(i, fileList)
        totalRows += SNS_IS_ImportOneSummaryCSV(pathName, fileName)
    endfor
    KillDataFolder/Z root:SNSIslandSummaryImportTemp

    if ((expectedRows >= 0) && (totalRows != expectedRows))
        Abort "SNS_ImportIslandSummaryCSVs: imported " + num2istr(totalRows) + " rows; expected " + num2istr(expectedRows)
    endif

    return totalRows
End


//==============================================================================
// SNS_ImportIslandSummaryCSVsFromPaths
//
// Purpose:
//   Replace root:IslandSummary with supplemental-summary CSV files that reside
//   in separate folders. Each semicolon-separated symbolic path is paired with
//   the filename at the same list position. Schema, duplicate, and final-row
//   validation are identical to SNS_ImportIslandSummaryCSVs.
//
// Inputs:
//   pathList     : semicolon-separated Igor symbolic path names
//   fileList     : semicolon-separated CSV filenames
//   expectedRows : required final row count; use a negative value to skip
//                  the final count check
//==============================================================================
Function SNS_ImportIslandSummaryCSVsFromPaths(pathList, fileList, expectedRows)
    String pathList, fileList
    Variable expectedRows

    Variable nPaths = ItemsInList(pathList)
    Variable nFiles = ItemsInList(fileList)
    if (nFiles < 1)
        Abort "SNS_ImportIslandSummaryCSVsFromPaths: fileList is empty."
    endif
    if (nPaths != nFiles)
        Abort "SNS_ImportIslandSummaryCSVsFromPaths: path and file counts differ."
    endif

    SNS_ClearIslandSummary()
    SNS_InitIslandSummary()

    Variable i, totalRows = 0
    for (i = 0; i < nFiles; i += 1)
        String pathName = StringFromList(i, pathList)
        String fileName = StringFromList(i, fileList)
        totalRows += SNS_IS_ImportOneSummaryCSV(pathName, fileName)
    endfor
    KillDataFolder/Z root:SNSIslandSummaryImportTemp

    if ((expectedRows >= 0) && (totalRows != expectedRows))
        Abort "SNS_ImportIslandSummaryCSVsFromPaths: imported " + num2istr(totalRows) + " rows; expected " + num2istr(expectedRows)
    endif

    return totalRows
End


//------------------------------------------------------------------------------
// Return 1 for a finite, strictly positive scalar.
//------------------------------------------------------------------------------
Static Function SNS_IS_IsPositiveFinite(value)
    Variable value

    return ((numtype(value) == 0) && (value > 0))
End


//------------------------------------------------------------------------------
// Read the Fermi velocity used by the island-summary E0(L) model.
//------------------------------------------------------------------------------
Static Function/D SNS_IS_GetFermiVelocity()
    NVAR/Z vF = root:SNS_Settings:vF
    if (!NVAR_Exists(vF) || !SNS_IS_IsPositiveFinite(vF))
        Abort "SNS_IS_E0vsL_FitFunc: root:SNS_Settings:vF is missing or invalid."
    endif

    return vF
End


//==============================================================================
// SNS_IS_PrepareODRFitWaves
//
// Purpose:
//   Preserve the established Fig. 2 ODR preparation without STMtools. Include
//   only finite points with positive x, y, and strictly positive uncertainties.
//   Excluded points receive unit uncertainties so FuncFit never encounters
//   zero or NaN weights before applying the mask.
//==============================================================================
Function SNS_IS_PrepareODRFitWaves(xW, yW, sxW, syW, sxFitW, syFitW, maskW, label)
    Wave xW, yW, sxW, syW, sxFitW, syFitW, maskW
    String label

    Variable n = numpnts(xW)
    if ((numpnts(yW) != n) || (numpnts(sxW) != n) || (numpnts(syW) != n))
        Abort "SNS_IS_PrepareODRFitWaves: inconsistent wave lengths for " + label
    endif

    Redimension/N=(n) sxFitW, syFitW, maskW

    Variable i
    Variable nFit = 0
    for (i = 0; i < n; i += 1)
        if ((numtype(xW[i]) == 0) && (numtype(yW[i]) == 0) && (numtype(sxW[i]) == 0) && (numtype(syW[i]) == 0) && (xW[i] > 0) && (yW[i] > 0) && (abs(sxW[i]) > 0) && (abs(syW[i]) > 0))
            maskW[i] = 1
            sxFitW[i] = abs(sxW[i])
            syFitW[i] = abs(syW[i])
            nFit += 1
        else
            maskW[i] = 0
            sxFitW[i] = 1
            syFitW[i] = 1
        endif
    endfor

    if (nFit < 2)
        Abort "SNS_IS_PrepareODRFitWaves: fewer than two valid ODR points for " + label
    endif

    return nFit
End


//==============================================================================
// SNS_IS_MakeSortedE0PlotWaves_ClusterMaxL
//
// Purpose:
//   Preserve the established Fig. 2 repeated-length grouping without
//   STMtools. Rows with the same ClusterMaxL are combined using mean E0,
//   outermost E0 bounds, maximum length error, and the larger asymmetric E0
//   error as the symmetric ODR uncertainty. Outputs are sorted by E0.
//==============================================================================
Function SNS_IS_MakeSortedE0PlotWaves_ClusterMaxL()
    Wave E0AvgW = root:IslandSummary:E0Avg_meV
    Wave E0LowW = root:IslandSummary:E0Low_meV
    Wave E0HighW = root:IslandSummary:E0High_meV
    Wave LW = root:IslandSummary:ClusterMaxL_nm
    Wave LErrW = root:IslandSummary:ClusterMaxL_err_nm

    Variable n = numpnts(E0AvgW)
    if ((numpnts(E0LowW) != n) || (numpnts(E0HighW) != n) || (numpnts(LW) != n) || (numpnts(LErrW) != n))
        Abort "SNS_IS_MakeSortedE0PlotWaves_ClusterMaxL: inconsistent summary-table wave lengths."
    endif

    Make/O/D/N=0 root:IslandSummary:SNS_IS_tmp_E0Avg_meV
    Make/O/D/N=0 root:IslandSummary:SNS_IS_tmp_E0Low_meV
    Make/O/D/N=0 root:IslandSummary:SNS_IS_tmp_E0High_meV
    Make/O/D/N=0 root:IslandSummary:SNS_IS_tmp_ClusterMaxL_nm
    Make/O/D/N=0 root:IslandSummary:SNS_IS_tmp_ClusterMaxL_err_nm

    Wave tmpE0Avg = root:IslandSummary:SNS_IS_tmp_E0Avg_meV
    Wave tmpE0Low = root:IslandSummary:SNS_IS_tmp_E0Low_meV
    Wave tmpE0High = root:IslandSummary:SNS_IS_tmp_E0High_meV
    Wave tmpL = root:IslandSummary:SNS_IS_tmp_ClusterMaxL_nm
    Wave tmpLErr = root:IslandSummary:SNS_IS_tmp_ClusterMaxL_err_nm

    Variable i, ok, nTmp = 0
    Variable lo, hi
    for (i = 0; i < n; i += 1)
        ok = 1
        ok = ok && (numtype(E0AvgW[i]) == 0)
        ok = ok && (numtype(E0LowW[i]) == 0)
        ok = ok && (numtype(E0HighW[i]) == 0)
        ok = ok && (numtype(LW[i]) == 0)
        ok = ok && (numtype(LErrW[i]) == 0)
        ok = ok && (E0AvgW[i] > 0)
        ok = ok && (LW[i] > 0)

        if (ok)
            InsertPoints nTmp, 1, tmpE0Avg, tmpE0Low, tmpE0High, tmpL, tmpLErr
            lo = min(E0LowW[i], E0HighW[i])
            hi = max(E0LowW[i], E0HighW[i])
            tmpE0Avg[nTmp] = E0AvgW[i]
            tmpE0Low[nTmp] = lo
            tmpE0High[nTmp] = hi
            tmpL[nTmp] = LW[i]
            tmpLErr[nTmp] = abs(LErrW[i])
            nTmp += 1
        endif
    endfor

    if (nTmp < 1)
        Abort "SNS_IS_MakeSortedE0PlotWaves_ClusterMaxL: no valid E0/ClusterMaxL rows found."
    endif

    Sort tmpL, tmpL, tmpLErr, tmpE0Avg, tmpE0Low, tmpE0High

    Make/O/D/N=0 root:IslandSummary:E0_sort_meV
    Make/O/D/N=0 root:IslandSummary:E0_sort_err_low_meV
    Make/O/D/N=0 root:IslandSummary:E0_sort_err_high_meV
    Make/O/D/N=0 root:IslandSummary:E0_sort_err_meV
    Make/O/D/N=0 root:IslandSummary:ClusterMaxL_sort_nm
    Make/O/D/N=0 root:IslandSummary:ClusterMaxL_sort_err_nm

    Wave E0Out = root:IslandSummary:E0_sort_meV
    Wave E0ErrLo = root:IslandSummary:E0_sort_err_low_meV
    Wave E0ErrHi = root:IslandSummary:E0_sort_err_high_meV
    Wave E0ErrFit = root:IslandSummary:E0_sort_err_meV
    Wave LOut = root:IslandSummary:ClusterMaxL_sort_nm
    Wave LErrOut = root:IslandSummary:ClusterMaxL_sort_err_nm

    Variable nOut = 0
    Variable j, groupCount
    Variable LVal, tol, sumE0, lowBound, highBound, maxLErr
    Variable meanE0, errLo, errHi

    i = 0
    do
        LVal = tmpL[i]
        tol = max(1e-9, abs(LVal)*1e-12)
        j = i
        groupCount = 0
        sumE0 = 0
        lowBound = tmpE0Low[i]
        highBound = tmpE0High[i]
        maxLErr = abs(tmpLErr[i])

        do
            if (abs(tmpL[j] - LVal) > tol)
                break
            endif
            sumE0 += tmpE0Avg[j]
            lowBound = min(lowBound, tmpE0Low[j])
            highBound = max(highBound, tmpE0High[j])
            maxLErr = max(maxLErr, abs(tmpLErr[j]))
            groupCount += 1
            j += 1
        while (j < nTmp)

        meanE0 = sumE0 / groupCount
        errLo = abs(meanE0 - lowBound)
        errHi = abs(highBound - meanE0)
        InsertPoints nOut, 1, E0Out, E0ErrLo, E0ErrHi, E0ErrFit, LOut, LErrOut
        E0Out[nOut] = meanE0
        E0ErrLo[nOut] = errLo
        E0ErrHi[nOut] = errHi
        E0ErrFit[nOut] = max(errLo, errHi)
        LOut[nOut] = LVal
        LErrOut[nOut] = maxLErr
        nOut += 1
        i = j
    while (i < nTmp)

    Sort E0Out, E0Out, E0ErrLo, E0ErrHi, E0ErrFit, LOut, LErrOut
    Duplicate/O LOut, root:IslandSummary:ClusterMaxL_shade_nm
    Duplicate/O LOut, root:IslandSummary:ClusterMaxL_xerr_nm
    KillWaves/Z tmpE0Avg, tmpE0Low, tmpE0High, tmpL, tmpLErr
End


//==============================================================================
// SNS_IS_B0vsW_FitFunc
//
// Model: B0 = Phi0/(2*lambda_eff*W).
// coefW[0] is lambda_eff [m], W_nm is [nm], return value is [mT].
//==============================================================================
Function SNS_IS_B0vsW_FitFunc(coefW, W_nm) : FitFunc
    Wave coefW
    Variable W_nm

    Variable lambdaEff_m = coefW[0]
    if (!SNS_IS_IsPositiveFinite(lambdaEff_m) || !SNS_IS_IsPositiveFinite(W_nm))
        return NaN
    endif

    Variable Phi0_Wb = 2.0678338484619295e-15
    Variable W_m = W_nm * 1e-9
    Variable B0_T = Phi0_Wb / (2 * lambdaEff_m * W_m)
    return B0_T / 1e-3
End


//------------------------------------------------------------------------------
// Zero-field n=0 ballistic SNS root equation and derivative.
//------------------------------------------------------------------------------
Static Function/D SNS_IS_E0Equation(L_m, Delta_eV, E_eV, vF)
    Variable L_m, Delta_eV, E_eV, vF

    Variable hbar_eVs = 6.582119569e-16
    return (2*E_eV*L_m)/(hbar_eVs*vF) - 2*acos(E_eV/Delta_eV)
End


Static Function/D SNS_IS_E0EquationDerivative(L_m, Delta_eV, E_eV, vF)
    Variable L_m, Delta_eV, E_eV, vF

    Variable hbar_eVs = 6.582119569e-16
    Variable u = E_eV/Delta_eV
    if ((u <= 0) || (u >= 1))
        return 1e30
    endif

    return (2*L_m)/(hbar_eVs*vF) + 2/(Delta_eV*sqrt(1-u*u))
End


//------------------------------------------------------------------------------
// Safeguarded Newton solution of the zero-field n=0 ballistic SNS equation.
//------------------------------------------------------------------------------
Static Function/D SNS_IS_SolveE0_eV(L_m, Delta_eV, vF)
    Variable L_m, Delta_eV, vF

    if (!SNS_IS_IsPositiveFinite(L_m) || !SNS_IS_IsPositiveFinite(Delta_eV) || !SNS_IS_IsPositiveFinite(vF))
        return NaN
    endif

    Variable a = Delta_eV * 1e-14
    Variable b = Delta_eV * (1-1e-14)
    Variable fa = SNS_IS_E0Equation(L_m, Delta_eV, a, vF)
    Variable fb = SNS_IS_E0Equation(L_m, Delta_eV, b, vF)
    if ((numtype(fa) != 0) || (numtype(fb) != 0) || (fa*fb > 0))
        return Delta_eV
    endif

    Variable hbar_eVs = 6.582119569e-16
    Variable E_eV = pi*hbar_eVs*vF/(2*L_m)
    if ((numtype(E_eV) != 0) || (E_eV <= a) || (E_eV >= b))
        E_eV = 0.5*(a+b)
    endif

    Variable iter, f, df, Enew
    for (iter = 0; iter < 60; iter += 1)
        f = SNS_IS_E0Equation(L_m, Delta_eV, E_eV, vF)
        if (numtype(f) != 0)
            return Delta_eV
        endif
        if (abs(f) < 1e-13)
            break
        endif
        if (f > 0)
            b = E_eV
        else
            a = E_eV
        endif
        df = SNS_IS_E0EquationDerivative(L_m, Delta_eV, E_eV, vF)
        if (!SNS_IS_IsPositiveFinite(df))
            Enew = 0.5*(a+b)
        else
            Enew = E_eV - f/df
            if ((numtype(Enew) != 0) || (Enew <= a) || (Enew >= b))
                Enew = 0.5*(a+b)
            endif
        endif
        E_eV = Enew
    endfor

    return E_eV
End


//==============================================================================
// SNS_IS_E0vsL_FitFunc
//
// Zero-field n=0 ballistic SNS model. coefW[0] is Delta [meV], x is L [nm],
// and the return value is E0 [meV].
//==============================================================================
Function SNS_IS_E0vsL_FitFunc(coefW, x) : FitFunc
    Wave coefW
    Variable x

    Variable Delta_meV = coefW[0]
    if (!SNS_IS_IsPositiveFinite(x) || !SNS_IS_IsPositiveFinite(Delta_meV))
        return NaN
    endif

    Variable vF = SNS_IS_GetFermiVelocity()
    Variable E0_eV = SNS_IS_SolveE0_eV(x*1e-9, Delta_meV*1e-3, vF)
    if (numtype(E0_eV) != 0)
        return NaN
    endif

    return E0_eV * 1e3
End


//==============================================================================
// SNS_AppendIslandSummaryResult
//
// Purpose:
//   Insert or update one island/set result in root:IslandSummary.
//
// Optional inputs:
//   geomRelErrFrac : fractional uncertainty assigned to geometry values (0.10)
//   h_Cu           : Cu-island height (nm)
//   B0ErrScale     : multiplier applied to the supplied B0 uncertainty (1)
//   area_nm2       : island area (nm^2)
//   perim_nm       : island perimeter (nm)
//==============================================================================
Function SNS_AppendIslandSummaryResult(dataName, B0_mT, B0_err_mT, E0Low_meV, E0High_meV, globalMaxL_nm, clusterMaxL_nm, clusterMedianL_nm, globalMaxW_nm, clusterMaxW_nm, clusterMedianW_nm, [geomRelErrFrac, h_Cu, B0ErrScale, area_nm2, perim_nm])
    String dataName
    Variable B0_mT, B0_err_mT
    Variable E0Low_meV, E0High_meV
    Variable globalMaxL_nm, clusterMaxL_nm, clusterMedianL_nm
    Variable globalMaxW_nm, clusterMaxW_nm, clusterMedianW_nm
    Variable geomRelErrFrac, h_Cu, B0ErrScale, area_nm2, perim_nm

    if (ParamIsDefault(geomRelErrFrac))
        geomRelErrFrac = 0.10
    endif
    geomRelErrFrac = abs(geomRelErrFrac)

    if (ParamIsDefault(B0ErrScale))
        B0ErrScale = 1
    endif
    B0ErrScale = abs(B0ErrScale)

    Variable hCuVal = NaN
    Variable areaVal = NaN
    Variable perimVal = NaN
    if (!ParamIsDefault(h_Cu))
        hCuVal = h_Cu
    endif
    if (!ParamIsDefault(area_nm2))
        areaVal = abs(area_nm2)
    endif
    if (!ParamIsDefault(perim_nm))
        perimVal = abs(perim_nm)
    endif

    Variable pOverAVal = NaN
    Variable aOverPVal = NaN
    if ((numtype(areaVal) == 0) && (numtype(perimVal) == 0) && (areaVal > 0) && (perimVal > 0))
        pOverAVal = perimVal / areaVal
        aOverPVal = areaVal / perimVal
    endif

    Variable E0Lo = min(E0Low_meV, E0High_meV)
    Variable E0Hi = max(E0Low_meV, E0High_meV)
    Variable B0Val = abs(B0_mT)
    Variable B0ErrVal = abs(B0_err_mT) * B0ErrScale
    Variable E0AvgVal = 0.5 * (E0Lo + E0Hi)
    Variable E0AvgErrVal = 0.5 * abs(E0Hi - E0Lo)

    Variable globalMaxLVal = abs(globalMaxL_nm)
    Variable clusterMaxLVal = abs(clusterMaxL_nm)
    Variable clusterMedianLVal = abs(clusterMedianL_nm)
    Variable globalMaxWVal = abs(globalMaxW_nm)
    Variable clusterMaxWVal = abs(clusterMaxW_nm)
    Variable clusterMedianWVal = abs(clusterMedianW_nm)

    SNS_InitIslandSummary()

    Wave/T wDataName = root:IslandSummary:DataName
    Wave whCu = root:IslandSummary:h_Cu
    Wave wArea = root:IslandSummary:Area_nm2
    Wave wPerim = root:IslandSummary:Perim_nm
    Wave wPoverA = root:IslandSummary:PoverA_nmInv
    Wave wAoverP = root:IslandSummary:AoverP_nm
    Wave wB0 = root:IslandSummary:B0_mT
    Wave wB0Err = root:IslandSummary:B0_err_mT
    Wave wE0Low = root:IslandSummary:E0Low_meV
    Wave wE0High = root:IslandSummary:E0High_meV
    Wave wE0Avg = root:IslandSummary:E0Avg_meV
    Wave wE0AvgErr = root:IslandSummary:E0Avg_err_meV
    Wave wGlobalMaxL = root:IslandSummary:GlobalMaxL_nm
    Wave wGlobalMaxLErr = root:IslandSummary:GlobalMaxL_err_nm
    Wave wClusterMaxL = root:IslandSummary:ClusterMaxL_nm
    Wave wClusterMaxLErr = root:IslandSummary:ClusterMaxL_err_nm
    Wave wClusterMedianL = root:IslandSummary:ClusterMedianL_nm
    Wave wClusterMedianLErr = root:IslandSummary:ClusterMedianL_err_nm
    Wave wGlobalMaxW = root:IslandSummary:GlobalMaxW_nm
    Wave wGlobalMaxWErr = root:IslandSummary:GlobalMaxW_err_nm
    Wave wClusterMaxW = root:IslandSummary:ClusterMaxW_nm
    Wave wClusterMaxWErr = root:IslandSummary:ClusterMaxW_err_nm
    Wave wClusterMedianW = root:IslandSummary:ClusterMedianW_nm
    Wave wClusterMedianWErr = root:IslandSummary:ClusterMedianW_err_nm

    Variable n = numpnts(wDataName)
    Variable i
    Variable idx = -1
    for (i = 0; i < n; i += 1)
        if (CmpStr(wDataName[i], dataName) == 0)
            idx = i
            break
        endif
    endfor

    if (idx < 0)
        idx = n
        InsertPoints idx, 1, wDataName, whCu, wArea, wPerim, wPoverA, wAoverP, wB0, wB0Err, wE0Low, wE0High, wE0Avg, wE0AvgErr, wGlobalMaxL, wGlobalMaxLErr, wClusterMaxL, wClusterMaxLErr, wClusterMedianL, wClusterMedianLErr, wGlobalMaxW, wGlobalMaxWErr, wClusterMaxW, wClusterMaxWErr, wClusterMedianW, wClusterMedianWErr
    endif

    wDataName[idx] = dataName
    whCu[idx] = hCuVal
    wArea[idx] = areaVal
    wPerim[idx] = perimVal
    wPoverA[idx] = pOverAVal
    wAoverP[idx] = aOverPVal
    wB0[idx] = B0Val
    wB0Err[idx] = B0ErrVal
    wE0Low[idx] = E0Lo
    wE0High[idx] = E0Hi
    wE0Avg[idx] = E0AvgVal
    wE0AvgErr[idx] = E0AvgErrVal
    wGlobalMaxL[idx] = globalMaxLVal
    wGlobalMaxLErr[idx] = geomRelErrFrac * globalMaxLVal
    wClusterMaxL[idx] = clusterMaxLVal
    wClusterMaxLErr[idx] = geomRelErrFrac * clusterMaxLVal
    wClusterMedianL[idx] = clusterMedianLVal
    wClusterMedianLErr[idx] = geomRelErrFrac * clusterMedianLVal
    wGlobalMaxW[idx] = globalMaxWVal
    wGlobalMaxWErr[idx] = geomRelErrFrac * globalMaxWVal
    wClusterMaxW[idx] = clusterMaxWVal
    wClusterMaxWErr[idx] = geomRelErrFrac * clusterMaxWVal
    wClusterMedianW[idx] = clusterMedianWVal
    wClusterMedianWErr[idx] = geomRelErrFrac * clusterMedianWVal
End


//==============================================================================
// SNS_ShowIslandSummaryTable
//
// Purpose:
//   Display the complete cross-island summary table.
//==============================================================================
Function SNS_ShowIslandSummaryTable()
    SNS_InitIslandSummary()

    Wave/T wDataName = root:IslandSummary:DataName
    Wave whCu = root:IslandSummary:h_Cu
    Wave wArea = root:IslandSummary:Area_nm2
    Wave wPerim = root:IslandSummary:Perim_nm
    Wave wPoverA = root:IslandSummary:PoverA_nmInv
    Wave wAoverP = root:IslandSummary:AoverP_nm
    Wave wB0 = root:IslandSummary:B0_mT
    Wave wB0Err = root:IslandSummary:B0_err_mT
    Wave wE0Low = root:IslandSummary:E0Low_meV
    Wave wE0High = root:IslandSummary:E0High_meV
    Wave wE0Avg = root:IslandSummary:E0Avg_meV
    Wave wE0AvgErr = root:IslandSummary:E0Avg_err_meV
    Wave wGlobalMaxL = root:IslandSummary:GlobalMaxL_nm
    Wave wGlobalMaxLErr = root:IslandSummary:GlobalMaxL_err_nm
    Wave wClusterMaxL = root:IslandSummary:ClusterMaxL_nm
    Wave wClusterMaxLErr = root:IslandSummary:ClusterMaxL_err_nm
    Wave wClusterMedianL = root:IslandSummary:ClusterMedianL_nm
    Wave wClusterMedianLErr = root:IslandSummary:ClusterMedianL_err_nm
    Wave wGlobalMaxW = root:IslandSummary:GlobalMaxW_nm
    Wave wGlobalMaxWErr = root:IslandSummary:GlobalMaxW_err_nm
    Wave wClusterMaxW = root:IslandSummary:ClusterMaxW_nm
    Wave wClusterMaxWErr = root:IslandSummary:ClusterMaxW_err_nm
    Wave wClusterMedianW = root:IslandSummary:ClusterMedianW_nm
    Wave wClusterMedianWErr = root:IslandSummary:ClusterMedianW_err_nm

    DoWindow/K SNS_IslandSummaryTable
    Edit/K=1 wDataName, whCu, wArea, wPerim, wPoverA, wAoverP, wB0, wB0Err, wE0Low, wE0High, wE0Avg, wE0AvgErr, wGlobalMaxL, wGlobalMaxLErr, wClusterMaxL, wClusterMaxLErr, wClusterMedianL, wClusterMedianLErr, wGlobalMaxW, wGlobalMaxWErr, wClusterMaxW, wClusterMaxWErr, wClusterMedianW, wClusterMedianWErr
    DoWindow/C SNS_IslandSummaryTable
End


//==============================================================================
// SNS_ClearIslandSummary
//
// Purpose:
//   Delete the complete cross-island summary table.
//==============================================================================
Function SNS_ClearIslandSummary()
    KillDataFolder/Z root:IslandSummary
End
