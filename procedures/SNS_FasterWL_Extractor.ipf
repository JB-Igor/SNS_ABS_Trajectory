#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

//==============================================================================
// LWMax-only helpers
// New standalone functions for fast extraction of
//   L_max(r)  and  W_max(r)
// in nm, without touching the original SNS DOS functions.
//==============================================================================

#pragma rtGlobals=3

//==============================================================================
// SNS_EstimateNphi_FromMask_LWMaxOnly
//
// Estimate angular channel count for LWMax-only workflow.
// Inputs:
//   Nmask       : 2D mask wave
//   lambdaF_nm  : Fermi wavelength [nm]
//
// Returns:
//   estimated Nphi
//==============================================================================
Function SNS_EstimateNphi_FromMask_LWMaxOnly(Nmask, lambdaF_nm)
    Wave Nmask
    Variable lambdaF_nm

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if (nx <= 1 || ny <= 1 || lambdaF_nm <= 0)
        return 0
    endif

    Variable dx = DimDelta(Nmask, 0)
    Variable dy = DimDelta(Nmask, 1)
    Variable cellArea = abs(dx * dy)

    Variable ix, iy
    Variable areaN = 0

    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)
            if (Nmask[ix][iy] > 0.5)
                areaN += cellArea
            endif
        endfor
    endfor

    if (areaN <= 0)
        return 0
    endif

    Variable Lchar = sqrt(areaN)
    Variable Nphi = round(Pi * Lchar / lambdaF_nm)

    if (Nphi < 1)
        Nphi = 1
    endif

    return Nphi
End


//==============================================================================
// SNS_BuildChannelsFromMask2D_LWMaxOnly_nm
//
// Fast geometry-only version of SNS_BuildChannelsFromMask2D().
// For a single point r0 inside Nmask, trace all 2D SNS chords exactly as in
// the original builder, but only keep
//
//   Lmax_nm = max(L_N)
//   Wmax_nm = max(W_eff)
//
// No channel waves are created.
//
// Inputs:
//   Nmask        : 2D N-region mask, axes in nm
//   r0x, r0y     : point inside mask [nm]
//   phiB         : in-plane field angle [rad]
//   stepFac      : ray marching refinement factor
//
// Outputs (by reference):
//   Lmax_nm      : maximum chord length [nm]
//   Wmax_nm      : maximum magnetic width [nm]
//
// Optional:
//   NphiOverride : if > 0, use this angular resolution instead of automatic
//                  estimate
//
// Returns:
//   0   success
//  -1   invalid mask
//  -2   no valid S-N-S chord through r0
//==============================================================================
Function SNS_BuildChannelsFromMask2D_LWMaxOnly_nm(Nmask, r0x, r0y, phiB, stepFac, Lmax_nm, Wmax_nm, [NphiOverride])

    Wave Nmask
    Variable r0x, r0y
    Variable phiB
    Variable stepFac
    Variable &Lmax_nm, &Wmax_nm
    Variable NphiOverride

    Lmax_nm = NaN
    Wmax_nm = NaN

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if (nx <= 1 || ny <= 1)
        return -1
    endif

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    Variable x1 = x0 + dx*(nx - 1)
    Variable y1 = y0 + dy*(ny - 1)

    Variable xMin = min(x0, x1)
    Variable xMax = max(x0, x1)
    Variable yMin = min(y0, y1)
    Variable yMax = max(y0, y1)

    // starting point must lie inside N
    Variable ix0 = round((r0x - x0)/dx)
    Variable iy0 = round((r0y - y0)/dy)

    if (ix0 < 0 || ix0 >= nx || iy0 < 0 || iy0 >= ny)
        return -2
    endif
    if (Nmask[ix0][iy0] <= 0.5)
        return -2
    endif

    // angular resolution
    Variable Nphi
    if (!ParamIsDefault(NphiOverride) && (NphiOverride > 0))
        Nphi = round(NphiOverride)
    else
        STRUCT SNS_Params params
        SNS_LoadParams(params)
        Variable lambdaF_nm = params.LambdaF * 1e9
        Nphi = SNS_EstimateNphi_FromMask_LWMaxOnly(Nmask, lambdaF_nm)
    endif

    if (Nphi < 1)
        Nphi = 1
    endif

    // marching step
    Variable baseStep = min(abs(dx), abs(dy))
    if (stepFac <= 0)
        stepFac = 0.5
    endif
    Variable ds = stepFac * baseStep

    // direction perpendicular to B
    Variable epx = -sin(phiB)
    Variable epy =  cos(phiB)

    Variable j, phi, vx, vy
    Variable x, y, x_prev, y_prev
    Variable ix, iy, ix_prev, iy_prev
    Variable hit1x, hit1y, hit2x, hit2y
    Variable have1, have2
    Variable dxChord, dyChord, L_N_nm, W_eff_nm
    Variable foundAny = 0

    Variable Lmax_local_nm = -Inf
    Variable Wmax_local_nm = -Inf

    for (j = 0; j < Nphi; j += 1)

        phi = (j + 0.5) * Pi / Nphi
        vx  = cos(phi)
        vy  = sin(phi)

        // -------- march along +v --------
        x = r0x
        y = r0y
        x_prev = x
        y_prev = y
        ix_prev = round((x_prev - x0)/dx)
        iy_prev = round((y_prev - y0)/dy)
        have1 = 0

        do
            x += ds * vx
            y += ds * vy

            if (x < xMin || x > xMax || y < yMin || y > yMax)
                break
            endif

            ix = round((x - x0)/dx)
            iy = round((y - y0)/dy)
            if (ix < 0 || ix >= nx || iy < 0 || iy >= ny)
                break
            endif

            if (Nmask[ix][iy] <= 0.5)
                hit1x = x_prev
                hit1y = y_prev
                have1 = 1
                break
            endif

            x_prev = x
            y_prev = y
            ix_prev = ix
            iy_prev = iy
        while (1)

        // -------- march along -v --------
        x = r0x
        y = r0y
        x_prev = x
        y_prev = y
        ix_prev = round((x_prev - x0)/dx)
        iy_prev = round((y_prev - y0)/dy)
        have2 = 0

        do
            x -= ds * vx
            y -= ds * vy

            if (x < xMin || x > xMax || y < yMin || y > yMax)
                break
            endif

            ix = round((x - x0)/dx)
            iy = round((y - y0)/dy)
            if (ix < 0 || ix >= nx || iy < 0 || iy >= ny)
                break
            endif

            if (Nmask[ix][iy] <= 0.5)
                hit2x = x_prev
                hit2y = y_prev
                have2 = 1
                break
            endif

            x_prev = x
            y_prev = y
            ix_prev = ix
            iy_prev = iy
        while (1)

        if (have1 && have2)
            dxChord = hit2x - hit1x
            dyChord = hit2y - hit1y

            L_N_nm   = sqrt(dxChord*dxChord + dyChord*dyChord)
            W_eff_nm = abs(dxChord*epx + dyChord*epy)

            if (L_N_nm > 0)
                foundAny = 1
                if (L_N_nm > Lmax_local_nm)
                    Lmax_local_nm = L_N_nm
                endif
                if (W_eff_nm > Wmax_local_nm)
                    Wmax_local_nm = W_eff_nm
                endif
            endif
        endif
    endfor

    if (!foundAny)
        return -2
    endif

    Lmax_nm = Lmax_local_nm
    Wmax_nm = Wmax_local_nm

    return 0
End


//==============================================================================
// SNS_LineLWMax_FromMask_LWMaxOnly_nm
//
// Build R_axis, L_max, W_max along a line using the lightweight
// SNS_BuildChannelsFromMask2D_LWMaxOnly_nm().
//
// Inputs:
//   Nmask              : 2D mask wave, axes in nm
//   phiB               : in-plane field direction [rad]
//   xStart, yStart     : line start [nm]
//   xEnd, yEnd         : line end   [nm]
//   stepFac            : marching refinement for the ray tracing
//
// Optional:
//   dr_nm              : line-point spacing [nm]
//                        default = lambdaF in nm from SNS settings
//   NphiOverride       : fixed angular resolution for faster preview runs
//
// Outputs in caller folder:
//   R_axis_LWmaxOnly_<phiDeg>deg   [nm]
//   W_max_LWmaxOnly_<phiDeg>deg    [nm]
//   L_max_LWmaxOnly_<phiDeg>deg    [nm]
//
// Returns:
//   Nr = number of sampled line points
//==============================================================================
Function SNS_LineLWMax_FromMask_LWMaxOnly_nm(Nmask, phiB, xStart, yStart, xEnd, yEnd, stepFac, [dr_nm, NphiOverride])

    Wave Nmask
    Variable phiB
    Variable xStart, yStart, xEnd, yEnd
    Variable stepFac
    Variable dr_nm
    Variable NphiOverride

    DFREF dfrCaller = GetDataFolderDFR()

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    if (ParamIsDefault(dr_nm) || (dr_nm <= 0))
        dr_nm = params.LambdaF * 1e9
    endif

    Variable dxLine = xEnd - xStart
    Variable dyLine = yEnd - yStart
    Variable Lline_nm = sqrt(dxLine*dxLine + dyLine*dyLine)

    Variable Nr
    if (Lline_nm <= 0)
        Nr = 1
    else
        Nr = round(Lline_nm / dr_nm) + 1
        if (Nr < 2)
            Nr = 2
        endif
    endif

    Variable drEff_nm
    if (Nr > 1)
        drEff_nm = Lline_nm / (Nr - 1)
    else
        drEff_nm = 1
    endif

    SetDataFolder dfrCaller

    String tag = "LWmaxOnly_" + num2str(round(phiB*180/pi)) + "deg"
    String nameRaxis = "R_axis_" + tag
    String nameWmax  = "W_max_"  + tag
    String nameLmax  = "L_max_"  + tag

    Make/O/D/N=(Nr) $nameRaxis
    Make/O/D/N=(Nr) $nameWmax
    Make/O/D/N=(Nr) $nameLmax

    Wave Rline      = $nameRaxis
    Wave W_max_Line = $nameWmax
    Wave L_max_Line = $nameLmax

    Rline = (Nr > 1) ? p*drEff_nm : 0
    W_max_Line = NaN
    L_max_Line = NaN

    if (Nr > 1)
        SetScale/P x, 0, drEff_nm, "nm", Rline
        SetScale/P x, 0, drEff_nm, "nm", W_max_Line
        SetScale/P x, 0, drEff_nm, "nm", L_max_Line
    else
        SetScale/P x, 0, 1, "nm", Rline
        SetScale/P x, 0, 1, "nm", W_max_Line
        SetScale/P x, 0, 1, "nm", L_max_Line
    endif

    Variable ir, t, r0x, r0y
    Variable Lmax_nm, Wmax_nm
    Variable err

    for (ir = 0; ir < Nr; ir += 1)

        if (Nr > 1)
            t = ir / (Nr - 1)
        else
            t = 0
        endif

        r0x = xStart + t*dxLine
        r0y = yStart + t*dyLine

        if (!ParamIsDefault(NphiOverride) && (NphiOverride > 0))
            err = SNS_BuildChannelsFromMask2D_LWMaxOnly_nm(Nmask, r0x, r0y, phiB, stepFac, Lmax_nm, Wmax_nm, NphiOverride=NphiOverride)
        else
            err = SNS_BuildChannelsFromMask2D_LWMaxOnly_nm(Nmask, r0x, r0y, phiB, stepFac, Lmax_nm, Wmax_nm)
        endif

        if (err == 0)
            L_max_Line[ir] = Lmax_nm
            W_max_Line[ir] = Wmax_nm
        else
            L_max_Line[ir] = NaN
            W_max_Line[ir] = NaN
        endif
    endfor

    return Nr
End


//==============================================================================
// SNS_LineLWMax_FromMask_Compat_LWMaxOnly_nm
//
// Compatibility wrapper keeping your usual old-style call:
//
//   SNS_LineLWMax_FromMask_Compat_LWMaxOnly_nm(
//       Nmask, phiB, xStart, yStart, xEnd, yEnd, stepFac, Zbarrier
//       [, dr_nm, NphiOverride] )
//
// Zbarrier is ignored in this lightweight workflow.
//==============================================================================
Function SNS_LineLWMax_FromMask_Compat_LWMaxOnly_nm(Nmask, phiB, xStart, yStart, xEnd, yEnd, stepFac, Zbarrier, [dr_nm, NphiOverride])

    Wave Nmask
    Variable phiB
    Variable xStart, yStart, xEnd, yEnd
    Variable stepFac, Zbarrier
    Variable dr_nm, NphiOverride

    if (!ParamIsDefault(dr_nm) && !ParamIsDefault(NphiOverride))
        return SNS_LineLWMax_FromMask_LWMaxOnly_nm(Nmask, phiB, xStart, yStart, xEnd, yEnd, stepFac, dr_nm=dr_nm, NphiOverride=NphiOverride)
    elseif (!ParamIsDefault(dr_nm))
        return SNS_LineLWMax_FromMask_LWMaxOnly_nm(Nmask, phiB, xStart, yStart, xEnd, yEnd, stepFac, dr_nm=dr_nm)
    elseif (!ParamIsDefault(NphiOverride))
        return SNS_LineLWMax_FromMask_LWMaxOnly_nm(Nmask, phiB, xStart, yStart, xEnd, yEnd, stepFac, NphiOverride=NphiOverride)
    else
        return SNS_LineLWMax_FromMask_LWMaxOnly_nm(Nmask, phiB, xStart, yStart, xEnd, yEnd, stepFac)
    endif
End