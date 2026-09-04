//==============================================================================
// SNS_GeometryFromMask.ipf
//
// Mask-derived SNS geometry utilities:
//   - Extract channel endpoints (Hit lists) from an N-region mask
//   - Build line sampling positions for LDOS cuts
//   - Estimate semiclassical channel count Nphi
//   - Optional vortex A-field and Delta vortex map helpers
//   - Folder-level visualization helper (uses SIDAM if available)
//
// NOTE: This file is part of the SNS (S-N-S) code path.
//       Legacy SSpS code lives in SNS_Legacy_SSpS.ipf.
//==============================================================================

#pragma rtGlobals=3
#pragma IgorVersion=9.00

// Apply the optional SIDAM styling used by ray-overlay displays. When SIDAM
// is unavailable, use an Igor-native reversed grayscale image with the same
// fixed display range.
static Function SNS_StyleRayDisplayImage(image, winImage)
    Wave image
    String winImage

#if Exists("SIDAMColor")
    SIDAMColor(imgList=NameOfWave(image), rev=1)
#if Exists("SIDAMRange")
    SIDAMRange(grfName=winImage, zmin=-0.5, zmax=0.5)
#else
    String imageName = NameOfWave(image)
    ModifyImage/W=$winImage $imageName ctab={-0.5,0.5,Grays,1}
#endif
#else
    String imageName = NameOfWave(image)
    ModifyImage/W=$winImage $imageName ctab={-0.5,0.5,Grays,1}
#endif

    return 0
End

//==============================================================================
// SNS_BuildChannelsFromMask2D
//
// Purpose
//   Build the SNS channel ensemble for a single STS position r0 from a 2D
//   N-region mask, by tracing straight S–N–S chords through r0. The resulting
//   channel list (geometry + transparency + weights) is used as input for the
//   ABS/DOS solver.
//
// Inputs
//   Nmask    : 2D wave, N-region mask in sample coordinates
//              1 = inside N, 0 = outside; x/y axes must be scaled in nm.
//
//   r0x      : STS x-position in Nmask coordinates [nm].
//   r0y      : STS y-position in Nmask coordinates [nm].
//
//   phiB     : In-plane magnetic-field direction with respect to +x [rad].
//              Channel directions φ run 0..π; the orbital lever arm W_eff is
//              taken perpendicular to B.
//
//   stepFac  : Dimensionless step refinement for ray marching.
//              Actual step in mask units:
//              ds = stepFac * min(|Δx|, |Δy|)
//              where Δx, Δy are the axis scaling of Nmask [nm].
//              Recommended: 0.3–0.5 (≈ 2–3 steps per grid cell).
//
//
//   folder   : Data-folder path (string).
//              If non-empty, all output waves are created there (NewDataFolder/O/S).
//              If empty, outputs are created in the current data folder.
//
// Return value
//   0   : success
//  -1   : invalid Nmask (nx ≤ 1 or ny ≤ 1)
//  -2   : no valid S–N–S chord through r0 (no pair of hits found)
//
// Created output waves (1D, length = nValid channels)
//   phiList     : Channel direction angle φ ∈ (0,π) [rad].
//
//   L_N_List    : N-region chord length |r₂ − r₁| in mask units [nm].
//   W_eff_List  : Magnetic lever arm relative to B:
//                   W_eff = |(r₂ − r₁) · ê⊥|  [nm],
//                 where ê⊥ is the in-plane unit vector ⟂ B.
//
//   wChan       : Channel weights / multiplicities. Downstream solvers
//                 normalize by sum(wChan) where needed.
//                 (Currently set to 1 for all valid channels before normalization.)
//
//   T_eff_List  : Effective normal-state transmission 0 ≤ Teff ≤ 1,
//                 Teff = min(TL, TR) from BTK formula above.
//
//   Hit1x_List  : x-coordinate of first ( “left” ) impact point on N boundary [nm].
//   Hit1y_List  : y-coordinate of first impact point [nm].
//   Hit2x_List  : x-coordinate of second ( “right” ) impact point on N boundary [nm].
//   Hit2y_List  : y-coordinate of second impact point [nm].
//
// Notes
//   • Nphi is estimated from Nmask area and λF via SNS_EstimateNphi_FromMask.
//   • All geometric outputs (L_N_List, W_eff_List, Hit*_*List) are in the same
//     nm coordinate system as Nmask’s axis scaling.
//   • This routine only constructs geometry + Teff + weights; ABS spectra and
//     DOS(E,B) are computed by downstream solver functions.
//
//==============================================================================

// Current canonical output names:
//   phiList_rad
//   L_N_List_nm, W_eff_List_nm
//   Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm
//   wChan, T_eff_List
// All geometry outputs from this function are in nm. Downstream solver-facing
// code is responsible for explicit local nm -> m conversion where needed.
// Optional:
//   angularOffsetFrac : fractional offset of the azimuthal Fermi-grid spacing.
//                       offset angle = angularOffsetFrac * (pi/Nphi).
//                       Default 0 preserves the historical midpoint grid.
//   angularAverageN   : interleaved angular-offset averaging control.
//                       Default 1 preserves historical behavior.
//                       If >= 1: number of grids to append.
//                       If 0 < value < 1: fractional sub-offset step.
//                         Example angularAverageN=0.1 appends 10 grids.
Function SNS_BuildChannelsFromMask2D(Nmask, r0x, r0y, phiB, stepFac, folder, [angularOffsetFrac, angularAverageN])
    Wave    Nmask
    Variable r0x        
    Variable r0y        
    Variable phiB         
    Variable stepFac       
    String  folder
    Variable angularOffsetFrac
    Variable angularAverageN

 // --- Load standard SNS settings (Igor Pro 9) ---
    STRUCT SNS_Params params
    SNS_LoadParams(params)

    // Override physical parameters with initialized values
    Variable lambdaF    = params.LambdaF      // [m]
    
    Variable Nphi = round(SNS_EstimateNphi_FromMask(Nmask,lambdaF*1e9)) // Angular resolution estimated from island area and Fermi-wavelength
    Variable angularOffsetFracLocal
    if (ParamIsDefault(angularOffsetFrac) || numtype(angularOffsetFrac) != 0)
        angularOffsetFracLocal = 0
    else
        angularOffsetFracLocal = angularOffsetFrac - floor(angularOffsetFrac)
        if (angularOffsetFracLocal < 0)
            angularOffsetFracLocal += 1
        endif
    endif
    Variable nAngularAvg
    if (ParamIsDefault(angularAverageN) || numtype(angularAverageN) != 0)
        nAngularAvg = 1
    elseif ((angularAverageN > 0) && (angularAverageN < 1))
        nAngularAvg = max(1, round(1/angularAverageN))
    else
        nAngularAvg = max(1, round(angularAverageN))
    endif
    Variable angularAverageStepFracLocal = 1 / nAngularAvg

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if (nx <= 1 || ny <= 1)
        return -1    // invalid mask
    endif

    Variable x0  = DimOffset(Nmask, 0)
    Variable dx  = DimDelta(Nmask, 0)
    Variable y0  = DimOffset(Nmask, 1)
    Variable dy  = DimDelta(Nmask, 1)

    Variable xMin = x0
    Variable xMax = x0 + dx*(nx - 1)
    Variable yMin = y0
    Variable yMax = y0 + dy*(ny - 1)

    Variable baseStep = min(abs(dx), abs(dy))
    if (stepFac <= 0)
        stepFac = 0.5
    endif
    Variable ds = stepFac*baseStep

    // perpendicular to B
    Variable epx = -sin(phiB)
    Variable epy =  cos(phiB)

    String oldDF = GetDataFolder(1)
    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    // temporary raw arrays (size Nphi * nAngularAvg)
    Variable NphiRaw = Nphi * nAngularAvg
    Make/FREE/D/N=(NphiRaw) phi_raw, L_N_raw, W_eff_raw, w_raw, T_raw
    Make/FREE/D/N=(NphiRaw) hit1x_raw, hit1y_raw, hit2x_raw, hit2y_raw

    phi_raw   = NaN
    L_N_raw   = NaN
    W_eff_raw = NaN
    w_raw     = 0
    T_raw     = NaN

    hit1x_raw = NaN
    hit1y_raw = NaN
    hit2x_raw = NaN
    hit2y_raw = NaN

    Variable ia, j, idxRaw, phi, vx, vy
    Variable x, y, x_prev, y_prev
    Variable ix, iy, ix_prev, iy_prev
    Variable hit1x, hit1y, hit2x, hit2y
    Variable ix1_in, iy1_in, ix2_in, iy2_in
    Variable have1, have2
    Variable offsetThis

    for (ia = 0; ia < nAngularAvg; ia += 1)

        offsetThis = angularOffsetFracLocal + ia * angularAverageStepFracLocal

        for (j = 0; j < Nphi; j += 1)

        idxRaw = ia*Nphi + j

        phi = mod((j + 0.5 + offsetThis)*Pi/Nphi, Pi)
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
            x += ds*vx
            y += ds*vy

            if (x < xMin || x > xMax || y < yMin || y > yMax)
                break
            endif

            ix = round((x - x0)/dx)
            iy = round((y - y0)/dy)
            if (ix < 0 || ix >= nx || iy < 0 || iy >= ny)
                break
            endif

            if (Nmask[ix][iy] <= 0.5)
                // crossed N → outside; last point & index were inside N
                hit1x   = x_prev
                hit1y   = y_prev
                ix1_in  = ix_prev
                iy1_in  = iy_prev
                have1   = 1
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
            x -= ds*vx
            y -= ds*vy

            if (x < xMin || x > xMax || y < yMin || y > yMax)
                break
            endif

            ix = round((x - x0)/dx)
            iy = round((y - y0)/dy)
            if (ix < 0 || ix >= nx || iy < 0 || iy >= ny)
                break
            endif

            if (Nmask[ix][iy] <= 0.5)
                hit2x   = x_prev
                hit2y   = y_prev
                ix2_in  = ix_prev
                iy2_in  = iy_prev
                have2   = 1
                break
            endif

            x_prev = x
            y_prev = y
            ix_prev = ix
            iy_prev = iy

        while (1)

        if (have1 && have2)
            Variable dxChord = hit2x - hit1x
            Variable dyChord = hit2y - hit1y
            Variable L_N     = sqrt(dxChord*dxChord + dyChord*dyChord)
            Variable W_eff   = abs(dxChord*epx + dyChord*epy)

            if (L_N > 0)
                // -------- interface normals from Nmask gradient --------
                Variable ixm, ixp, iym, iyp
                Variable gx1, gy1, gx2, gy2
                Variable n1x, n1y, n2x, n2y, nlen

                // left interface
                ixm = max(ix1_in - 1, 0)
                ixp = min(ix1_in + 1, nx - 1)
                iym = max(iy1_in - 1, 0)
                iyp = min(iy1_in + 1, ny - 1)

                gx1 = Nmask[ixp][iy1_in] - Nmask[ixm][iy1_in]
                gy1 = Nmask[ix1_in][iyp] - Nmask[ix1_in][iym]
                nlen = sqrt(gx1*gx1 + gy1*gy1)
                if (nlen > 0)
                    n1x = gx1/nlen
                    n1y = gy1/nlen
                else
                    n1x = 0
                    n1y = 0
                endif

                // right interface
                ixm = max(ix2_in - 1, 0)
                ixp = min(ix2_in + 1, nx - 1)
                iym = max(iy2_in - 1, 0)
                iyp = min(iy2_in + 1, ny - 1)

                gx2 = Nmask[ixp][iy2_in] - Nmask[ixm][iy2_in]
                gy2 = Nmask[ix2_in][iyp] - Nmask[ix2_in][iym]
                nlen = sqrt(gx2*gx2 + gy2*gy2)
                if (nlen > 0)
                    n2x = gx2/nlen
                    n2y = gy2/nlen
                else
                    n2x = 0
                    n2y = 0
                endif

                // -------- angular transparencies TL, TR --------
                Variable cosInc1 = abs(vx*n1x + vy*n1y)
                Variable cosInc2 = abs(vx*n2x + vy*n2y)
                cosInc1 = min(1, max(0, cosInc1))
                cosInc2 = min(1, max(0, cosInc2))

					Variable TL, TR
					TL = SNS_ChannelTransmissionFromCos(cosInc1)
					TR = SNS_ChannelTransmissionFromCos(cosInc2)

                Variable Teff = min(TL, TR)

                phi_raw[idxRaw]   = phi
                L_N_raw[idxRaw]   = L_N
                W_eff_raw[idxRaw] = W_eff
                w_raw[idxRaw]     = 1
                T_raw[idxRaw]     = Teff

                // store impact points for this channel
                hit1x_raw[idxRaw] = hit1x
                hit1y_raw[idxRaw] = hit1y
                hit2x_raw[idxRaw] = hit2x
                hit2y_raw[idxRaw] = hit2y
            endif
        endif

    endfor
    endfor

    // -------- compact to valid channels --------
    Variable nValid = 0
    for (j = 0; j < NphiRaw; j += 1)
        if (w_raw[j] > 0)
            nValid += 1
        endif
    endfor
    if (nValid == 0)
        SetDataFolder $oldDF
        return -2    // no S–N–S chord
    endif

    Make/O/D/N=(nValid) phiList_rad, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List
    Make/O/D/N=(nValid) Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm

    Variable k = 0
    for (j = 0; j < NphiRaw; j += 1)
        if (w_raw[j] > 0)
            phiList_rad[k]     = phi_raw[j]
            L_N_List_nm[k]     = L_N_raw[j]
            W_eff_List_nm[k]   = W_eff_raw[j]
            wChan[k]           = 1 // Semiclassical free-propagation prefactor would scale ~1/sqrt(L),
            T_eff_List[k]      = T_raw[j]

            Hit1x_List_nm[k]   = hit1x_raw[j]
            Hit1y_List_nm[k]   = hit1y_raw[j]
            Hit2x_List_nm[k]   = hit2x_raw[j]
            Hit2y_List_nm[k]   = hit2y_raw[j]

            k += 1
        endif
    endfor

    Variable/G v_Nphi_base = Nphi
    Variable/G v_angularOffsetFrac = angularOffsetFracLocal
    Variable/G v_angularAverageN = nAngularAvg
    Variable/G v_angularAverageStepFrac = angularAverageStepFracLocal

    // --- convert nm → m for solver compatibility ---
    // Geometry outputs remain in nm; solver-facing conversions happen downstream.
    SetDataFolder $oldDF
    return 0
End

//==============================================================================
// SNS_BuildChannelsFromRectModes2D
//
// Purpose
//   Tier-0 analytic mode selector for low-mode diagnostics.  The function
//   approximates the mask by an island-aligned effective rectangle: the outer
//   bounding rectangle in the mask principal-axis frame is uniformly shrunk so
//   its area matches the actual mask area.  The function then enumerates
//   hard-wall rectangle quantum numbers (n,m), keeps the modes whose |k_nm|
//   are closest to k_F, converts each mode to its two unoriented ray families,
//   and traces those ray angles through the actual mask.
//
//   This is not a replacement for SNS_BuildChannelsFromMask2D().  It is an
//   opt-in way to choose physically motivated ray angles for an approximately
//   rectangular island without treating a dense angular-offset grid as extra
//   physical modes.
//
// Inputs
//   Nmask   : 2D binary N-region mask scaled in nm.
//   r0x,r0y : STM tip position in mask coordinates [nm].
//   phiB    : in-plane magnetic-field direction [rad].
//   stepFac : ray-marching step refinement.
//   folder  : output data folder. Empty string means current folder.
//
// Optional
//   nModes              : number of rectangle (n,m) states closest to k_F.
//                         Default: SNS_EstimateNphi_FromMask(...).
//   useTipWeight        : 1 default. Multiply each mode by
//                         |psi_nm(r_tip)|^2 in the representative rectangle.
//   includeMirrorAngles : 1 default. Trace alpha and pi-alpha for each mode,
//                         splitting the mode weight across both branches.
//
// Outputs
//   Canonical channel waves plus rectangle-mode diagnostics:
//     RectMode_n, RectMode_m, RectMode_E_eV, RectMode_EminusEF_eV
//     RectMode_kx_1_per_nm, RectMode_ky_1_per_nm
//     RectMode_alpha_rad, RectMode_tipWeight, RectMode_angleBranch
//     RectMode_Lx_nm, RectMode_Ly_nm, RectMode_xTip_nm, RectMode_yTip_nm
//     RectMode_AreaMask_nm2, RectMode_AreaOuter_nm2, RectMode_AreaScale
//     RectMode_Axis1_unit, RectMode_Axis2_unit, RectMode_Center_nm
//     RectMode_bandModel
//
// Notes
//   Mode selection is by |k_nm - k_F| for both 2DEG and spinless TI mode.
//   RectMode_E_eV is diagnostic only:
//     bandModel 0: E = hbar^2 k^2 / 2m*
//     bandModel 1: E = hbar v_D k
//
// Returns
//   0  : success
//  -1  : invalid mask
//  -2  : no valid ray chords from selected rectangle modes
//==============================================================================
Function SNS_BuildChannelsFromRectModes2D(Nmask, r0x, r0y, phiB, stepFac, folder, [nModes, useTipWeight, includeMirrorAngles])
    Wave Nmask
    Variable r0x, r0y, phiB, stepFac
    String folder
    Variable nModes
    Variable useTipWeight
    Variable includeMirrorAngles

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    if (ParamIsDefault(useTipWeight))
        useTipWeight = 1
    endif
    if (ParamIsDefault(includeMirrorAngles))
        includeMirrorAngles = 1
    endif

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if (nx <= 1 || ny <= 1)
        return -1
    endif

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    Variable ix, iy
    Variable ixMin = nx, ixMax = -1, iyMin = ny, iyMax = -1
    Variable nInside = 0
    Variable cellArea_nm2 = abs(dx*dy)
    Variable xCenter, yCenter
    Variable sumX = 0, sumY = 0
    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)
            if (Nmask[ix][iy] > 0.5)
                xCenter = x0 + ix*dx
                yCenter = y0 + iy*dy
                sumX += xCenter
                sumY += yCenter
                nInside += 1
                ixMin = min(ixMin, ix)
                ixMax = max(ixMax, ix)
                iyMin = min(iyMin, iy)
                iyMax = max(iyMax, iy)
            endif
        endfor
    endfor
    if ((ixMax < ixMin) || (iyMax < iyMin))
        return -1
    endif
    if (nInside <= 0)
        return -1
    endif

    Variable areaMask_nm2 = nInside * cellArea_nm2
    Variable xCentroid_nm = sumX / nInside
    Variable yCentroid_nm = sumY / nInside

    Variable cxx = 0, cyy = 0, cxy = 0
    Variable xr, yr
    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)
            if (Nmask[ix][iy] > 0.5)
                xCenter = x0 + ix*dx
                yCenter = y0 + iy*dy
                xr = xCenter - xCentroid_nm
                yr = yCenter - yCentroid_nm
                cxx += xr*xr
                cyy += yr*yr
                cxy += xr*yr
            endif
        endfor
    endfor
    cxx /= nInside
    cyy /= nInside
    cxy /= nInside

    Variable theta = 0.5 * atan2(2*cxy, cxx - cyy)
    Variable e1x = cos(theta)
    Variable e1y = sin(theta)
    Variable e2x = -sin(theta)
    Variable e2y = cos(theta)

    Variable u, v
    Variable uMin = 1e300, uMax = -1e300, vMin = 1e300, vMax = -1e300
    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)
            if (Nmask[ix][iy] > 0.5)
                xCenter = x0 + ix*dx
                yCenter = y0 + iy*dy
                xr = xCenter - xCentroid_nm
                yr = yCenter - yCentroid_nm
                u = xr*e1x + yr*e1y
                v = xr*e2x + yr*e2y
                uMin = min(uMin, u)
                uMax = max(uMax, u)
                vMin = min(vMin, v)
                vMax = max(vMax, v)
            endif
        endfor
    endfor

    Variable padU_nm = 0.5*(abs(dx*e1x) + abs(dy*e1y))
    Variable padV_nm = 0.5*(abs(dx*e2x) + abs(dy*e2y))
    Variable LxOuter_nm = (uMax - uMin) + 2*padU_nm
    Variable LyOuter_nm = (vMax - vMin) + 2*padV_nm
    Variable areaOuter_nm2 = LxOuter_nm * LyOuter_nm
    if (areaOuter_nm2 <= 0)
        return -1
    endif

    Variable areaScale = sqrt(areaMask_nm2 / areaOuter_nm2)
    areaScale = min(1, areaScale)
    Variable Lx_nm = LxOuter_nm * areaScale
    Variable Ly_nm = LyOuter_nm * areaScale

    Variable uCenter = 0.5*(uMin + uMax)
    Variable vCenter = 0.5*(vMin + vMax)
    Variable xRectCenter_nm = xCentroid_nm + uCenter*e1x + vCenter*e2x
    Variable yRectCenter_nm = yCentroid_nm + uCenter*e1y + vCenter*e2y

    Variable uTip = (r0x - xRectCenter_nm)*e1x + (r0y - yRectCenter_nm)*e1y
    Variable vTip = (r0x - xRectCenter_nm)*e2x + (r0y - yRectCenter_nm)*e2y
    Variable xTipRect_nm = uTip + 0.5*Lx_nm
    Variable yTipRect_nm = vTip + 0.5*Ly_nm

    Variable lambdaF_nm = params.LambdaF * 1e9
    Variable kF_nm = 2*pi / lambdaF_nm
    if (ParamIsDefault(nModes) || numtype(nModes) != 0 || nModes <= 0)
        nModes = round(SNS_EstimateNphi_FromMask(Nmask, lambdaF_nm))
    endif
    nModes = max(1, round(nModes))

    Variable nMax = max(1, ceil(kF_nm*Lx_nm/pi) + 4)
    Variable mMax = max(1, ceil(kF_nm*Ly_nm/pi) + 4)
    Variable nCand = nMax*mMax

    Make/FREE/D/N=(nCand) candDist, candN, candM, candE, candKx, candKy, candAlpha, candTipWeight

    Variable n, m, j, ic = 0
    Variable kx_nm, ky_nm, k_nm, E_eV, tipW
    Variable mStar = params.m_eff * m_e_SI
    for (n = 1; n <= nMax; n += 1)
        kx_nm = n*pi/Lx_nm
        for (m = 1; m <= mMax; m += 1)
            ky_nm = m*pi/Ly_nm
            k_nm = sqrt(kx_nm*kx_nm + ky_nm*ky_nm)
            if (params.SNS_bandModel == 1)
                E_eV = params.SNS_hbarvD_eVA * (k_nm / 10)
            else
                E_eV = (HBAR_SI^2 * (k_nm*1e9)^2) / (2*mStar*q_e)
            endif

            tipW = 1
            if (useTipWeight)
                if ((xTipRect_nm >= 0) && (xTipRect_nm <= Lx_nm) && (yTipRect_nm >= 0) && (yTipRect_nm <= Ly_nm))
                    tipW = sin(n*pi*xTipRect_nm/Lx_nm)^2 * sin(m*pi*yTipRect_nm/Ly_nm)^2
                else
                    tipW = 0
                endif
            endif

            candDist[ic] = abs(k_nm - kF_nm)
            candN[ic] = n
            candM[ic] = m
            candE[ic] = E_eV
            candKx[ic] = kx_nm
            candKy[ic] = ky_nm
            candAlpha[ic] = atan2(ky_nm, kx_nm)
            candTipWeight[ic] = tipW
            ic += 1
        endfor
    endfor

    Sort candDist, candDist, candN, candM, candE, candKx, candKy, candAlpha, candTipWeight

    Variable nSelModes = min(nModes, nCand)
    Variable nAnglesPerMode = includeMirrorAngles ? 2 : 1
    Variable nRaw = nSelModes*nAnglesPerMode

    Make/FREE/D/N=(nRaw) phi_raw, L_N_raw, W_eff_raw, w_raw, T_raw
    Make/FREE/D/N=(nRaw) hit1x_raw, hit1y_raw, hit2x_raw, hit2y_raw
    Make/FREE/D/N=(nRaw) modeN_raw, modeM_raw, modeE_raw, modeKx_raw, modeKy_raw, modeAlpha_raw, modeTipW_raw, modeBranch_raw
    phi_raw = NaN
    L_N_raw = NaN
    W_eff_raw = NaN
    w_raw = 0
    T_raw = NaN

    Variable xMin = x0
    Variable xMax = x0 + dx*(nx - 1)
    Variable yMin = y0
    Variable yMax = y0 + dy*(ny - 1)
    Variable baseStep = min(abs(dx), abs(dy))
    if (stepFac <= 0)
        stepFac = 0.5
    endif
    Variable ds = stepFac*baseStep
    Variable epx = -sin(phiB)
    Variable epy = cos(phiB)

    Variable imode, ib, idxRaw
    Variable phi, vx, vy, x, y, x_prev, y_prev
    Variable ix_prev, iy_prev, ix1_in, iy1_in, ix2_in, iy2_in
    Variable hit1x, hit1y, hit2x, hit2y
    Variable have1, have2

    for (imode = 0; imode < nSelModes; imode += 1)
        for (ib = 0; ib < nAnglesPerMode; ib += 1)
            idxRaw = imode*nAnglesPerMode + ib
            phi = candAlpha[imode]
            if (ib == 1)
                phi = pi - phi
            endif
            vx = cos(phi)*e1x + sin(phi)*e2x
            vy = cos(phi)*e1y + sin(phi)*e2y
            phi = mod(atan2(vy, vx), pi)
            if (phi < 0)
                phi += pi
            endif
            vx = cos(phi)
            vy = sin(phi)

            x = r0x
            y = r0y
            x_prev = x
            y_prev = y
            ix_prev = round((x_prev - x0)/dx)
            iy_prev = round((y_prev - y0)/dy)
            have1 = 0

            do
                x += ds*vx
                y += ds*vy
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
                    ix1_in = ix_prev
                    iy1_in = iy_prev
                    have1 = 1
                    break
                endif
                x_prev = x
                y_prev = y
                ix_prev = ix
                iy_prev = iy
            while (1)

            x = r0x
            y = r0y
            x_prev = x
            y_prev = y
            ix_prev = round((x_prev - x0)/dx)
            iy_prev = round((y_prev - y0)/dy)
            have2 = 0

            do
                x -= ds*vx
                y -= ds*vy
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
                    ix2_in = ix_prev
                    iy2_in = iy_prev
                    have2 = 1
                    break
                endif
                x_prev = x
                y_prev = y
                ix_prev = ix
                iy_prev = iy
            while (1)

            if (have1 && have2)
                Variable dxChord = hit2x - hit1x
                Variable dyChord = hit2y - hit1y
                Variable L_N = sqrt(dxChord*dxChord + dyChord*dyChord)
                Variable W_eff = abs(dxChord*epx + dyChord*epy)

                if (L_N > 0)
                    Variable ixm, ixp, iym, iyp
                    Variable gx1, gy1, gx2, gy2
                    Variable n1x, n1y, n2x, n2y, nlen

                    ixm = max(ix1_in - 1, 0)
                    ixp = min(ix1_in + 1, nx - 1)
                    iym = max(iy1_in - 1, 0)
                    iyp = min(iy1_in + 1, ny - 1)
                    gx1 = Nmask[ixp][iy1_in] - Nmask[ixm][iy1_in]
                    gy1 = Nmask[ix1_in][iyp] - Nmask[ix1_in][iym]
                    nlen = sqrt(gx1*gx1 + gy1*gy1)
                    if (nlen > 0)
                        n1x = gx1/nlen
                        n1y = gy1/nlen
                    else
                        n1x = 0
                        n1y = 0
                    endif

                    ixm = max(ix2_in - 1, 0)
                    ixp = min(ix2_in + 1, nx - 1)
                    iym = max(iy2_in - 1, 0)
                    iyp = min(iy2_in + 1, ny - 1)
                    gx2 = Nmask[ixp][iy2_in] - Nmask[ixm][iy2_in]
                    gy2 = Nmask[ix2_in][iyp] - Nmask[ix2_in][iym]
                    nlen = sqrt(gx2*gx2 + gy2*gy2)
                    if (nlen > 0)
                        n2x = gx2/nlen
                        n2y = gy2/nlen
                    else
                        n2x = 0
                        n2y = 0
                    endif

                    Variable cosInc1 = abs(vx*n1x + vy*n1y)
                    Variable cosInc2 = abs(vx*n2x + vy*n2y)
                    cosInc1 = min(1, max(0, cosInc1))
                    cosInc2 = min(1, max(0, cosInc2))

                    Variable TL, TR, Teff
                    TL = SNS_ChannelTransmissionFromCos(cosInc1)
                    TR = SNS_ChannelTransmissionFromCos(cosInc2)
                    Teff = min(TL, TR)

                    phi_raw[idxRaw] = phi
                    L_N_raw[idxRaw] = L_N
                    W_eff_raw[idxRaw] = W_eff
                    w_raw[idxRaw] = candTipWeight[imode] / nAnglesPerMode
                    T_raw[idxRaw] = Teff
                    hit1x_raw[idxRaw] = hit1x
                    hit1y_raw[idxRaw] = hit1y
                    hit2x_raw[idxRaw] = hit2x
                    hit2y_raw[idxRaw] = hit2y

                    modeN_raw[idxRaw] = candN[imode]
                    modeM_raw[idxRaw] = candM[imode]
                    modeE_raw[idxRaw] = candE[imode]
                    modeKx_raw[idxRaw] = candKx[imode]
                    modeKy_raw[idxRaw] = candKy[imode]
                    modeAlpha_raw[idxRaw] = candAlpha[imode]
                    modeTipW_raw[idxRaw] = candTipWeight[imode]
                    modeBranch_raw[idxRaw] = ib
                endif
            endif
        endfor
    endfor

    Variable nValid = 0
    for (j = 0; j < nRaw; j += 1)
        if (numtype(L_N_raw[j]) == 0)
            nValid += 1
        endif
    endfor
    if (nValid == 0)
        return -2
    endif

    String oldDF = GetDataFolder(1)
    DFREF dfrMaskSource = GetWavesDataFolderDFR(Nmask)
    SetDataFolder dfrMaskSource
    String maskImageList = WaveList("*_Z_mbgnd_xy", ";", "DIMS:2")
    String maskImageName = ""
    if (strlen(maskImageList) > 0)
        maskImageName = StringFromList(0, maskImageList, ";")
    endif
    SetDataFolder $oldDF

    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif
    DFREF dfrRectOut = GetDataFolderDFR()

    if ((strlen(folder) > 0) || (CmpStr(NameOfWave(Nmask), "w_mask") != 0))
        Duplicate/O Nmask, w_mask
    endif
    if (strlen(maskImageName) > 0)
        String staleRectImageList = WaveList("*_Z_mbgnd_xy", ";", "DIMS:2")
        Variable staleRectImageIndex
        for (staleRectImageIndex = 0; staleRectImageIndex < ItemsInList(staleRectImageList, ";"); staleRectImageIndex += 1)
            KillWaves/Z $(StringFromList(staleRectImageIndex, staleRectImageList, ";"))
        endfor
        SetDataFolder dfrMaskSource
        Wave RectImageSource = $maskImageName
        SetDataFolder dfrRectOut
        Duplicate/O RectImageSource, $maskImageName
    endif

    Make/O/D/N=(nValid) phiList_rad, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List
    Make/O/D/N=(nValid) Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm
    Make/O/D/N=(nValid) RectMode_n, RectMode_m, RectMode_E_eV, RectMode_EminusEF_eV
    Make/O/D/N=(nValid) RectMode_kx_1_per_nm, RectMode_ky_1_per_nm
    Make/O/D/N=(nValid) RectMode_alpha_rad, RectMode_tipWeight, RectMode_angleBranch
    Make/O/D/N=1 RectMode_Lx_nm, RectMode_Ly_nm, RectMode_xTip_nm, RectMode_yTip_nm, RectMode_kF_1_per_nm, RectMode_bandModel
    Make/O/D/N=1 RectMode_AreaMask_nm2, RectMode_AreaOuter_nm2, RectMode_AreaScale
    Make/O/D/N=2 RectMode_Axis1_unit, RectMode_Axis2_unit, RectMode_Center_nm

    Variable k = 0
    for (j = 0; j < nRaw; j += 1)
        if (numtype(L_N_raw[j]) == 0)
            phiList_rad[k] = phi_raw[j]
            L_N_List_nm[k] = L_N_raw[j]
            W_eff_List_nm[k] = W_eff_raw[j]
            wChan[k] = w_raw[j]
            T_eff_List[k] = T_raw[j]
            Hit1x_List_nm[k] = hit1x_raw[j]
            Hit1y_List_nm[k] = hit1y_raw[j]
            Hit2x_List_nm[k] = hit2x_raw[j]
            Hit2y_List_nm[k] = hit2y_raw[j]
            RectMode_n[k] = modeN_raw[j]
            RectMode_m[k] = modeM_raw[j]
            RectMode_E_eV[k] = modeE_raw[j]
            RectMode_EminusEF_eV[k] = modeE_raw[j] - params.E_F
            RectMode_kx_1_per_nm[k] = modeKx_raw[j]
            RectMode_ky_1_per_nm[k] = modeKy_raw[j]
            RectMode_alpha_rad[k] = modeAlpha_raw[j]
            RectMode_tipWeight[k] = modeTipW_raw[j]
            RectMode_angleBranch[k] = modeBranch_raw[j]
            k += 1
        endif
    endfor

    RectMode_Lx_nm[0] = Lx_nm
    RectMode_Ly_nm[0] = Ly_nm
    RectMode_xTip_nm[0] = xTipRect_nm
    RectMode_yTip_nm[0] = yTipRect_nm
    RectMode_kF_1_per_nm[0] = kF_nm
    RectMode_bandModel[0] = params.SNS_bandModel
    RectMode_AreaMask_nm2[0] = areaMask_nm2
    RectMode_AreaOuter_nm2[0] = areaOuter_nm2
    RectMode_AreaScale[0] = areaScale
    RectMode_Axis1_unit[0] = e1x
    RectMode_Axis1_unit[1] = e1y
    RectMode_Axis2_unit[0] = e2x
    RectMode_Axis2_unit[1] = e2y
    RectMode_Center_nm[0] = xRectCenter_nm
    RectMode_Center_nm[1] = yRectCenter_nm

    Variable/G v_RectMode_NSelectedModes = nSelModes
    Variable/G v_RectMode_NChannels = nValid
    Variable/G v_RectMode_UseTipWeight = useTipWeight
    Variable/G v_RectMode_IncludeMirrorAngles = includeMirrorAngles
    Variable/G v_RectMode_AreaMask_nm2 = areaMask_nm2
    Variable/G v_RectMode_AreaOuter_nm2 = areaOuter_nm2
    Variable/G v_RectMode_AreaScale = areaScale
    Variable/G v_RectMode_BandModel = params.SNS_bandModel

    SetDataFolder $oldDF
    return 0
End


//==============================================================================
// SNS_BuildLinePositions
//
// Purpose:
//   Construct output waves used by the SNS workflow.
//
// Inputs:
//   xStart : [nm]
//   yStart : [nm]
//   xEnd : [nm]
//   yEnd : [nm]
//   lambdaF : input
//   nameX : input
//   nameY : input
//   nameR : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_BuildLinePositions(xStart, yStart, xEnd, yEnd, lambdaF, nameX, nameY, nameR)
    Variable xStart, yStart, xEnd, yEnd     // nm
    Variable lambdaF                        // m
    String nameX, nameY, nameR

    Variable Lnm   = sqrt( (xEnd-xStart)^2 + (yEnd-yStart)^2 )   // nm
    Variable Lm    = Lnm*1e-9
    Variable ds    = lambdaF        // m, you can choose lambdaF/2 if you want denser
    Variable Nr    = max( 1, round(Lm/ds) )

    Make/O/D/N=(Nr) $nameX, $nameY, $nameR
    Wave Xline = $nameX
    Wave Yline = $nameY
    Wave Rline = $nameR

    Variable i, t
    for (i=0; i<Nr; i+=1)
        t        = i/(Nr-1.0)       // 0..1
        Xline[i] = xStart + t*(xEnd-xStart)
        Yline[i] = yStart + t*(yEnd-yStart)
        Rline[i] = t*Lnm            // distance along line (nm)
    endfor

    return Nr
End


//==============================================================================
// SNS_PrepareHighResMaskForPython
//
// Purpose:
//   Build a cropped, high-resolution binary mask for external Python normal-mode
//   diagnostics.
//
//   The source STM/topography mask often covers a much larger field of view than
//   the actual N island. This helper:
//     1. finds the occupied source-mask bounding box in nm,
//     2. adds physical padding,
//     3. crops the source mask,
//     4. resamples the cropped mask to a requested target spacing,
//     5. thresholds the result back to a binary mask.
//
// Inputs:
//   Nmask        : 2D binary mask, axes scaled in nm.
//   outMaskName  : output high-resolution mask wave name.
//   cropMaskName : output cropped source-resolution mask wave name.
//
// Optional Inputs:
//   targetDx_nm  : requested high-resolution spacing [nm].
//                  Default: 0.22 nm.
//   padding_nm   : padding added to occupied-mask bounding box [nm].
//                  Default: 5 nm.
//   threshold    : occupied-mask threshold.
//                  Default: 0.5.
//
// Outputs:
//   outMaskName  : high-resolution binary mask.
//   cropMaskName : cropped source-resolution binary mask.
//
// Returns:
//   0 on success.
//
// Notes:
//   This belongs in procedure code, not notebooks. Igor notebook cells should
//   call this helper rather than running loops.
//==============================================================================
Function SNS_PrepareHighResMaskForPython(Nmask, outMaskName, cropMaskName, [targetDx_nm, padding_nm, threshold])
    Wave Nmask
    String outMaskName, cropMaskName
    Variable targetDx_nm, padding_nm, threshold

    if (ParamIsDefault(targetDx_nm))
        targetDx_nm = 0.22
    endif
    if (ParamIsDefault(padding_nm))
        padding_nm = 5
    endif
    if (ParamIsDefault(threshold))
        threshold = 0.5
    endif

    if (targetDx_nm <= 0)
        Abort "SNS_PrepareHighResMaskForPython: targetDx_nm must be positive."
    endif

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if (nx <= 1 || ny <= 1)
        Abort "SNS_PrepareHighResMaskForPython: Nmask must be a 2D wave."
    endif

    Variable ix, iy
    Variable xHere_nm, yHere_nm
    Variable foundMaskPixel = 0
    Variable xMinFound_nm = Inf
    Variable xMaxFound_nm = -Inf
    Variable yMinFound_nm = Inf
    Variable yMaxFound_nm = -Inf

    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)
            if (Nmask[ix][iy] > threshold)
                xHere_nm = DimOffset(Nmask, 0) + ix*DimDelta(Nmask, 0)
                yHere_nm = DimOffset(Nmask, 1) + iy*DimDelta(Nmask, 1)
                xMinFound_nm = min(xMinFound_nm, xHere_nm)
                xMaxFound_nm = max(xMaxFound_nm, xHere_nm)
                yMinFound_nm = min(yMinFound_nm, yHere_nm)
                yMaxFound_nm = max(yMaxFound_nm, yHere_nm)
                foundMaskPixel = 1
            endif
        endfor
    endfor

    if (!foundMaskPixel)
        Abort "SNS_PrepareHighResMaskForPython: Nmask has no occupied pixels."
    endif

    Variable xAxisA_nm = DimOffset(Nmask, 0)
    Variable xAxisB_nm = DimOffset(Nmask, 0) + (nx - 1)*DimDelta(Nmask, 0)
    Variable yAxisA_nm = DimOffset(Nmask, 1)
    Variable yAxisB_nm = DimOffset(Nmask, 1) + (ny - 1)*DimDelta(Nmask, 1)

    Variable xLo_nm = max(xMinFound_nm - padding_nm, min(xAxisA_nm, xAxisB_nm))
    Variable xHi_nm = min(xMaxFound_nm + padding_nm, max(xAxisA_nm, xAxisB_nm))
    Variable yLo_nm = max(yMinFound_nm - padding_nm, min(yAxisA_nm, yAxisB_nm))
    Variable yHi_nm = min(yMaxFound_nm + padding_nm, max(yAxisA_nm, yAxisB_nm))

    Duplicate/O/R=(xLo_nm,xHi_nm)(yLo_nm,yHi_nm) Nmask, $cropMaskName
    Wave cropMask = $cropMaskName

    Variable nOutX = max(2, ceil(abs(xHi_nm - xLo_nm)/targetDx_nm) + 1)
    Variable nOutY = max(2, ceil(abs(yHi_nm - yLo_nm)/targetDx_nm) + 1)

    Make/O/D/N=(nOutX, nOutY) $outMaskName
    Wave outMask = $outMaskName
    SetScale/I x, xLo_nm, xHi_nm, "nm", outMask
    SetScale/I y, yLo_nm, yHi_nm, "nm", outMask

    outMask = interp2D(cropMask, DimOffset(outMask,0) + p*DimDelta(outMask,0), DimOffset(outMask,1) + q*DimDelta(outMask,1))
    outMask = numtype(outMask[p][q]) == 2 ? 0 : outMask[p][q]
    outMask = outMask[p][q] > threshold ? 1 : 0

    return 0
End


//==============================================================================
// SNS_PrepareAxisAlignedHighResMaskForPython
//
// Purpose:
//   Build a cropped, high-resolution binary mask in the island principal-axis
//   frame for external Python normal-mode diagnostics.
//
//   The output mask x/y axes are not lab x/y. They are:
//       x' =  (x-xc) cos(theta) + (y-yc) sin(theta)
//       y' = -(x-xc) sin(theta) + (y-yc) cos(theta)
//   where theta is the occupied-mask PCA long-axis angle in the lab frame.
//
//   If xPoint_nm / yPoint_nm and Bangle_deg are supplied, this helper also
//   stores transformed notebook-safe scalar variables in the caller folder:
//       v_PythonMaskTheta_deg
//       v_PythonMaskCenterX_nm, v_PythonMaskCenterY_nm
//       v_PythonMaskSTSx_nm, v_PythonMaskSTSy_nm
//       v_PythonMaskBangle_deg
//
// Inputs:
//   Nmask        : 2D binary mask, axes scaled in lab nm.
//   outMaskName  : output high-resolution axis-aligned mask wave name.
//   cropMaskName : output source-resolution axis-aligned crop wave name.
//
// Optional Inputs:
//   targetDx_nm  : requested high-resolution spacing [nm].
//                  Default: 0.22 nm.
//   padding_nm   : padding added to occupied-mask principal-axis bounds [nm].
//                  Default: 5 nm.
//   threshold    : occupied-mask threshold.
//                  Default: 0.5.
//   xPoint_nm,
//   yPoint_nm    : lab-frame point to transform into the output frame.
//   Bangle_deg   : lab-frame B angle to transform into the output frame.
//
// Returns:
//   0 on success.
//
// Notes:
//   Use this for the Python sparse sidecar when a rotated experimental mask is
//   close to rectangular. The sidecar then solves an axis-aligned mask, and the
//   Igor ray candidates must be generated in this same transformed frame using
//   v_PythonMaskSTSx_nm / v_PythonMaskSTSy_nm / v_PythonMaskBangle_deg.
//==============================================================================
Function SNS_PrepareAxisAlignedHighResMaskForPython(Nmask, outMaskName, cropMaskName, [targetDx_nm, padding_nm, threshold, xPoint_nm, yPoint_nm, Bangle_deg])
    Wave Nmask
    String outMaskName, cropMaskName
    Variable targetDx_nm, padding_nm, threshold
    Variable xPoint_nm, yPoint_nm, Bangle_deg

    if (ParamIsDefault(targetDx_nm))
        targetDx_nm = 0.22
    endif
    if (ParamIsDefault(padding_nm))
        padding_nm = 5
    endif
    if (ParamIsDefault(threshold))
        threshold = 0.5
    endif

    if (targetDx_nm <= 0)
        Abort "SNS_PrepareAxisAlignedHighResMaskForPython: targetDx_nm must be positive."
    endif

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if (nx <= 1 || ny <= 1)
        Abort "SNS_PrepareAxisAlignedHighResMaskForPython: Nmask must be a 2D wave."
    endif

    Variable dx = DimDelta(Nmask, 0)
    Variable dy = DimDelta(Nmask, 1)
    Variable x0 = DimOffset(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable ix, iy
    Variable xLab, yLab
    Variable nInside = 0
    Variable sumX = 0
    Variable sumY = 0

    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)
            if (Nmask[ix][iy] > threshold)
                xLab = x0 + ix*dx
                yLab = y0 + iy*dy
                sumX += xLab
                sumY += yLab
                nInside += 1
            endif
        endfor
    endfor

    if (nInside <= 0)
        Abort "SNS_PrepareAxisAlignedHighResMaskForPython: Nmask has no occupied pixels."
    endif

    Variable xCenter = sumX / nInside
    Variable yCenter = sumY / nInside

    Variable cxx = 0
    Variable cyy = 0
    Variable cxy = 0
    Variable xr, yr
    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)
            if (Nmask[ix][iy] > threshold)
                xLab = x0 + ix*dx
                yLab = y0 + iy*dy
                xr = xLab - xCenter
                yr = yLab - yCenter
                cxx += xr*xr
                cyy += yr*yr
                cxy += xr*yr
            endif
        endfor
    endfor
    cxx /= nInside
    cyy /= nInside
    cxy /= nInside

    Variable theta = 0.5 * atan2(2*cxy, cxx - cyy)
    Variable c = cos(theta)
    Variable s = sin(theta)

    Variable u, v
    Variable uMin = Inf
    Variable uMax = -Inf
    Variable vMin = Inf
    Variable vMax = -Inf
    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)
            if (Nmask[ix][iy] > threshold)
                xLab = x0 + ix*dx
                yLab = y0 + iy*dy
                xr = xLab - xCenter
                yr = yLab - yCenter
                u = xr*c + yr*s
                v = -xr*s + yr*c
                uMin = min(uMin, u)
                uMax = max(uMax, u)
                vMin = min(vMin, v)
                vMax = max(vMax, v)
            endif
        endfor
    endfor

    uMin -= padding_nm
    uMax += padding_nm
    vMin -= padding_nm
    vMax += padding_nm

    Variable dCrop = min(abs(dx), abs(dy))
    if (dCrop <= 0 || numtype(dCrop) != 0)
        dCrop = targetDx_nm
    endif

    Variable nCropX = max(2, ceil(abs(uMax - uMin)/dCrop) + 1)
    Variable nCropY = max(2, ceil(abs(vMax - vMin)/dCrop) + 1)
    Make/O/D/N=(nCropX, nCropY) $cropMaskName
    Wave cropMask = $cropMaskName
    SetScale/I x, uMin, uMax, "nm", cropMask
    SetScale/I y, vMin, vMax, "nm", cropMask

    Variable nOutX = max(2, ceil(abs(uMax - uMin)/targetDx_nm) + 1)
    Variable nOutY = max(2, ceil(abs(vMax - vMin)/targetDx_nm) + 1)
    Make/O/D/N=(nOutX, nOutY) $outMaskName
    Wave outMask = $outMaskName
    SetScale/I x, uMin, uMax, "nm", outMask
    SetScale/I y, vMin, vMax, "nm", outMask

    Variable uu, vv, xSrc, ySrc
    cropMask = NaN
    for (ix = 0; ix < nCropX; ix += 1)
        for (iy = 0; iy < nCropY; iy += 1)
            uu = DimOffset(cropMask,0) + ix*DimDelta(cropMask,0)
            vv = DimOffset(cropMask,1) + iy*DimDelta(cropMask,1)
            xSrc = xCenter + uu*c - vv*s
            ySrc = yCenter + uu*s + vv*c
            cropMask[ix][iy] = Interp2D(Nmask, xSrc, ySrc)
        endfor
    endfor
    cropMask = numtype(cropMask[p][q]) == 2 ? 0 : cropMask[p][q]
    cropMask = cropMask[p][q] > threshold ? 1 : 0

    outMask = NaN
    for (ix = 0; ix < nOutX; ix += 1)
        for (iy = 0; iy < nOutY; iy += 1)
            uu = DimOffset(outMask,0) + ix*DimDelta(outMask,0)
            vv = DimOffset(outMask,1) + iy*DimDelta(outMask,1)
            xSrc = xCenter + uu*c - vv*s
            ySrc = yCenter + uu*s + vv*c
            outMask[ix][iy] = Interp2D(Nmask, xSrc, ySrc)
        endfor
    endfor
    outMask = numtype(outMask[p][q]) == 2 ? 0 : outMask[p][q]
    outMask = outMask[p][q] > threshold ? 1 : 0

    Variable/G v_PythonMaskTheta_deg = theta*180/pi
    Variable/G v_PythonMaskCenterX_nm = xCenter
    Variable/G v_PythonMaskCenterY_nm = yCenter

    if (!ParamIsDefault(xPoint_nm) && !ParamIsDefault(yPoint_nm))
        Variable/G v_PythonMaskSTSx_nm = (xPoint_nm - xCenter)*c + (yPoint_nm - yCenter)*s
        Variable/G v_PythonMaskSTSy_nm = -(xPoint_nm - xCenter)*s + (yPoint_nm - yCenter)*c
    endif

    if (!ParamIsDefault(Bangle_deg))
        Variable bAng = Bangle_deg - theta*180/pi
        do
            bAng += 180
        while (bAng < 0)
        do
            bAng -= 180
        while (bAng >= 180)
        Variable/G v_PythonMaskBangle_deg = bAng
    endif

    return 0
End


//==============================================================================
// SNS_EstimateNphi_FromMask
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   Nmask : input
//   lambdaF : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_EstimateNphi_FromMask(Nmask, lambdaF)
    Wave    Nmask
    Variable lambdaF

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if (nx <= 1 || ny <= 1 || lambdaF <= 0)
        return 0
    endif

    Variable dx = DimDelta(Nmask, 0)
    Variable dy = DimDelta(Nmask, 1)
    Variable cellArea = abs(dx * dy)

    // --- compute N-region area ---
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

    // characteristic linear size
    Variable Lchar = sqrt(areaN)

    // N_phi ≈ π Lchar / λ_F
    Variable Nphi = round(Pi * Lchar / lambdaF)

    // safety floor
    if (Nphi < 1)
        Nphi = 1
    endif

    return Nphi
End

//==============================================================================
// SNS_MaskAreaPerim_FromParticles
//
// SAFE version: does not modify Nmask.
//
// Return wave:
//   w[0] = area_nm2
//   w[1] = perim_nm
//
// Input:
//   Nmask: 2D binary mask, N-region > 0.5, outside <= 0.5, axes in nm.
//==============================================================================
Function/WAVE SNS_MaskAreaPerim_FromParticles(Nmask)
    Wave Nmask

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)

    if (WaveDims(Nmask) != 2 || nx <= 1 || ny <= 1)
        Make/O/D/N=2 w
        w[0] = NaN
        w[1] = NaN
        return w
    endif

    Variable dx = abs(DimDelta(Nmask, 0))
    Variable dy = abs(DimDelta(Nmask, 1))

    if (dx <= 0 || dy <= 0)
        Make/O/D/N=2 w
        w[0] = NaN
        w[1] = NaN
        return w
    endif

    Variable area_nm2 = 0
    Variable perim_nm = 0

    Variable ix, iy

    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)

            if (Nmask[ix][iy] <= 0.5)
                continue
            endif

            // area
            area_nm2 += dx * dy

            // perimeter contribution from exposed pixel edges
            if (ix == 0 || Nmask[ix-1][iy] <= 0.5)
                perim_nm += dy
            endif
            if (ix == nx-1 || Nmask[ix+1][iy] <= 0.5)
                perim_nm += dy
            endif
            if (iy == 0 || Nmask[ix][iy-1] <= 0.5)
                perim_nm += dx
            endif
            if (iy == ny-1 || Nmask[ix][iy+1] <= 0.5)
                perim_nm += dx
            endif

        endfor
    endfor

    Make/O/D/N=2 w
    w[0] = area_nm2
    w[1] = perim_nm

    return w
End

//==============================================================================
// SNS_PlotChannelRay
//
// Purpose:
//   Plot diagnostic or result waves.
//
// Inputs:
//   Hit1x_List : input
//   Hit1y_List : input
//   Hit2x_List : input
//   Hit2y_List : input
//   idx : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_PlotChannelRay(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, idx)
    Wave    Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
    Variable idx

    Variable nCh = DimSize(Hit1x_List, 0)
    if (idx < 0 || idx >= nCh)
        Abort "SNS_PlotChannelRay: idx out of range."
    endif

    // Build unique names for this ray
    String rayXname, rayYname
    sprintf rayXname, "rayX_%d", idx
    sprintf rayYname, "rayY_%d", idx

    // 2-point line wave in real coordinates (same units as mask)
    Make/O/D/N=2 $rayXname, $rayYname
    Wave rayX = $rayXname
    Wave rayY = $rayYname

    rayX[0] = Hit1x_List[idx]
    rayY[0] = Hit1y_List[idx]
    rayX[1] = Hit2x_List[idx]
    rayY[1] = Hit2y_List[idx]

    // get top graph
    String topGraph = WinName(0, 1, 1)
    if (strlen(topGraph) == 0)
        return 0
    endif

    // append ray
    AppendToGraph/W=$topGraph rayY vs rayX

    // style this specific trace
    String cmd
    sprintf cmd, "ModifyGraph/W=%s lsize(%s)=2,lstyle(%s)=0", topGraph, rayYname, rayYname
    Execute cmd

    return 0
End



//==============================================================================
// SNS_MakeVortexA
//
// Purpose:
//   Construct output waves used by the SNS workflow.
//
// Inputs:
//   Nmask : input
//   AxName : input
//   AyName : input
//   xV_nm : input
//   yV_nm : input
//   nFlux : input
//   rc_nm : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_MakeVortexA(Nmask, AxName, AyName, xV_nm, yV_nm, nFlux, rc_nm)
    Wave   Nmask
    String AxName, AyName
    Variable xV_nm, yV_nm
    Variable nFlux
    Variable rc_nm

    if (WaveDims(Nmask) != 2)
        Abort "SNS_MakeVortexA: Nmask must be 2D."
    endif

    // We treat dim 0 as x, dim 1 as y to stay consistent with your geometry.
    Variable nX = DimSize(Nmask, 0)   // x index
    Variable nY = DimSize(Nmask, 1)   // y index

    // --- Create outputs ---
    Make/O/D/N=(nX, nY) $AxName, $AyName
    Wave Ax = $AxName
    Wave Ay = $AyName
    Ax = 0
    Ay = 0

    // Copy scaling from Nmask: dim 0 -> x, dim 1 -> y
    SetScale/P x, DimOffset(Nmask,0), DimDelta(Nmask,0), waveunits(Nmask,0), Ax, Ay
    SetScale/P y, DimOffset(Nmask,1), DimDelta(Nmask,1), waveunits(Nmask,1), Ax, Ay

    // Axis scaling
    Variable x0_nm = DimOffset(Nmask,0)
    Variable y0_nm = DimOffset(Nmask,1)
    Variable dx_nm = DimDelta(Nmask,0)   // x step
    Variable dy_nm = DimDelta(Nmask,1)   // y step

    // --- Flux parameters ---
    Variable h    = 6.62607015e-34       // J·s
    Variable q_e  = 1.602176634e-19      // C
    Variable Phi0 = h/(2*q_e)            // [Wb]
    Variable PhiV = nFlux * Phi0         // total vortex flux

    Variable ix, iy
    Variable x_nm, y_nm, rx, ry, r, Aphi, fCore
    Variable rMin_nm = 1e-3

    for (ix = 0; ix < nX; ix += 1)
        x_nm = x0_nm + ix*dx_nm

        for (iy = 0; iy < nY; iy += 1)
            y_nm = y0_nm + iy*dy_nm

            // Respect N region: zero outside N
            if (Nmask[ix][iy] != 1)
                Ax[ix][iy] = 0
                Ay[ix][iy] = 0
                continue
            endif

            rx = x_nm - xV_nm
            ry = y_nm - yV_nm
            r  = sqrt(rx*rx + ry*ry)

            if (r < rMin_nm)
                Ax[ix][iy] = 0
                Ay[ix][iy] = 0
                continue
            endif

            // Core smoothing: f(r/rc) = r^2 / (r^2 + rc^2)
            fCore = (r*r) / (r*r + rc_nm*rc_nm)

            // A_φ(r) = Φ_v / (2π r) * fCore
            Aphi = (PhiV / (2*pi*r)) * fCore

            // e_φ = (-sinθ, cosθ) = (-ry/r, rx/r)
            Ax[ix][iy] = -Aphi * (ry/r)   // A_x
            Ay[ix][iy] =  Aphi * (rx/r)   // A_y
        endfor
    endfor
End



//==============================================================================
// SNS_MakeDeltaVortexMap
//
// Purpose:
//   Construct output waves used by the SNS workflow.
//
// Inputs:
//   Nmask : input
//   nameDelta : [eV]
//   xV_nm : input
//   yV_nm : input
//   rc_nm : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_MakeDeltaVortexMap(Nmask, nameDelta, xV_nm, yV_nm, rc_nm)
    Wave   Nmask
    String nameDelta
    Variable xV_nm, yV_nm, rc_nm

    if (WaveDims(Nmask) != 2)
        Abort "SNS_MakeDeltaVortexMap: Nmask must be 2D."
    endif

    Variable nX = DimSize(Nmask, 0)   // x index
    Variable nY = DimSize(Nmask, 1)   // y index
    
    // --- Load standard SNS settings ---
    STRUCT SNS_Params params
    SNS_LoadParams(params)

    Variable Delta0      = params.Delta        // [eV]

    Make/O/D/N=(nX, nY) $nameDelta
    Wave DeltaMap = $nameDelta
    DeltaMap = Delta0   // default: bulk gap

    // copy scaling from Nmask (dim0->x, dim1->y)
    SetScale/P x, DimOffset(Nmask,0), DimDelta(Nmask,0), WaveUnits(Nmask,0), DeltaMap
    SetScale/P y, DimOffset(Nmask,1), DimDelta(Nmask,1), WaveUnits(Nmask,1), DeltaMap

    Variable x0_nm = DimOffset(Nmask,0)
    Variable dx_nm = DimDelta(Nmask,0)
    Variable y0_nm = DimOffset(Nmask,1)
    Variable dy_nm = DimDelta(Nmask,1)

    Variable ix, iy, x_nm, y_nm, rx, ry, r, u, fCore

    for (ix = 0; ix < nX; ix += 1)
        x_nm = x0_nm + ix*dx_nm
        for (iy = 0; iy < nY; iy += 1)
            y_nm = y0_nm + iy*dy_nm

            // optional: only modify inside N
            if (Nmask[ix][iy] != 1)
                DeltaMap[ix][iy] = Delta0
                continue
            endif

            rx = x_nm - xV_nm
            ry = y_nm - yV_nm
            r  = sqrt(rx*rx + ry*ry)

            // smooth core: f(r/rc) = 1 - exp(-(r/rc)^2)
            if (rc_nm > 0)
                u = r/rc_nm
                fCore = 1 - exp(-(u*u))
            else
                fCore = 1
            endif

            DeltaMap[ix][iy] = Delta0 * fCore
        endfor
    endfor
End






//==============================================================================
// SNS_ExtractModesForFolder
//
// Purpose:
//   Extract geometry-derived channels and related mode properties.
//
//   The function builds the 2D SNS channel ensemble for a selected STS point,
//   identifies special ballistic modes, and optionally creates diagnostic plots.
//
//   Display logic:
//     1) Topography with:
//          - gap-closing trajectory
//          - longest trajectory
//          - STS marker
//          - B-field arrow
//
//     2) Textbox containing:
//          - Ballistic block:
//              ℓ, w, B0
//          - optional Usadel 1D block, if precomputed values exist:
//              Ldiff, xprobe, Lφ, B0, and optionally D
//
//     3) Mode-property plot showing:
//          - L_N_List
//          - W_eff_List
//          - T_eff_List
//
// Inputs:
//   dfPath : input folder containing w_mask and image data.
//
// Optional:
//   Bangle_deg  : B-field angle in degrees. Default: 225.
//   STSx        : STS x position in nm. Default: -4.
//   STSy        : STS y position in nm. Default: -50.
//   Vortexx     : vortex x position in nm. Default: STSx.
//   Vortexy     : vortex y position in nm. Default: STSy.
//   doDisplay   : 1 default : create diagnostic plots.
//                 0         : skip all display/UI logic.
//   angularOffsetFrac
//               : optional fractional offset of the azimuthal Fermi-grid
//                 spacing passed to SNS_BuildChannelsFromMask2D().
//                 offset angle = angularOffsetFrac * (pi/Nphi).
//                 Default preserves the historical midpoint grid.
//   angularAverageN
//               : optional number of interleaved angular-offset grids passed
//                 to SNS_BuildChannelsFromMask2D(). Default 1 preserves the
//                 historical single-grid channel list. Values >1 append
//                 additional sub-offset channels into the same output waves.
//
// Outputs:
//   return : numeric
//            0 on successful completion.
//
// Created/updated variables:
//   root:LengthMax_nm
//   root:WPerp_nm
//   root:B0_ballistic_T
//
//   dfPath:v_Gap0
//   dfPath:v_Longest
//   dfPath:v_B0_ballistic_T
//
// Notes:
//   Usadel information is only displayed if the needed variables already exist
//   in dfPath, for example after running:
//
//      SNS_Usadel1D_GeometryFromRayFolder(...)
//      SNS_Usadel1D_LeverArmFromRayFolder(...)
//
//==============================================================================

Function SNS_ExtractModesForFolder(dfPath, [Bangle_deg, STSx, STSy, Vortexx, Vortexy, doDisplay, angularOffsetFrac, angularAverageN])
    String dfPath
    Variable Bangle_deg, STSx, STSy, Vortexx, Vortexy, doDisplay
    Variable angularOffsetFrac
    Variable angularAverageN

    String oldDF = GetDataFolder(1)
    SetDataFolder $dfPath

    Variable r=65535, g=65535, b=65535, a=65535


    // ---------------- Defaults ----------------
    if (ParamIsDefault(Bangle_deg))
        Bangle_deg = 225
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
    if (ParamIsDefault(doDisplay))
        doDisplay = 1
    endif

    Variable v_Bangle = Bangle_deg * pi/180

    // ---------------- Required waves ----------------
    Wave/Z w_mask
    if (!WaveExists(w_mask))
        SetDataFolder $oldDF
        Abort "SNS_ExtractModesForFolder: missing w_mask."
    endif

    String imgList = WaveList("*_Z_mbgnd_xy",";","DIMS:2")
    if (strlen(imgList)==0)
        SetDataFolder $oldDF
        Abort "SNS_ExtractModesForFolder: no *_Z_mbgnd_xy image found in folder."
    endif

    String imgName = StringFromList(0, imgList, ";")
    Wave image = $imgName

    // Do NOT use KillVariables/Z/A here.
    // That would remove precomputed v_Usadel_* variables from dfPath.

    // ---------------- Build geometry ----------------
    if (ParamIsDefault(angularOffsetFrac) && ParamIsDefault(angularAverageN))
        SNS_BuildChannelsFromMask2D(w_mask, STSx, STSy, v_Bangle, 0.5, "")
    elseif (ParamIsDefault(angularAverageN))
        SNS_BuildChannelsFromMask2D(w_mask, STSx, STSy, v_Bangle, 0.5, "", angularOffsetFrac=angularOffsetFrac)
    elseif (ParamIsDefault(angularOffsetFrac))
        SNS_BuildChannelsFromMask2D(w_mask, STSx, STSy, v_Bangle, 0.5, "", angularAverageN=angularAverageN)
    else
        SNS_BuildChannelsFromMask2D(w_mask, STSx, STSy, v_Bangle, 0.5, "", angularOffsetFrac=angularOffsetFrac, angularAverageN=angularAverageN)
    endif

    Wave L_N_List_nm
    Wave W_eff_List_nm
    Wave Hit1x_List_nm
    Wave Hit1y_List_nm
    Wave Hit2x_List_nm
    Wave Hit2y_List_nm
    Wave T_eff_List

    SNS_MaskAreaPerim_FromParticles(w_mask)
    Wave w

    Make/O/N=1 w_area_nm2 = w[0]
    Make/O/N=1 w_perim_nm = w[1]

    // ---------------- Identify special modes ----------------
    Make/O/D/N=(DimSize(L_N_List_nm,0)) A_eff_nm2 = L_N_List_nm * W_eff_List_nm

    WaveStats/Q W_eff_List_nm
    Variable/G v_Gap0 = V_maxloc

    WaveStats/Q L_N_List_nm
    Variable/G v_Longest = V_maxloc

    // ---------------- Store selected values ----------------
    Variable/G root:LengthMax_nm = L_N_List_nm[v_Longest]
    Variable/G root:WPerp_nm     = W_eff_List_nm[v_Gap0]

    // ---------------- Ballistic B0 estimate ----------------
    Variable ballisticWperp_m = abs(W_eff_List_nm[v_Gap0]) * 1e-9
    Variable ballisticB0_T = NaN

    NVAR/Z lambdaL_m = root:SNS_Settings:lambdaL
    if (NVAR_Exists(lambdaL_m) && ballisticWperp_m > 0)
        ballisticB0_T = pi * HBAR_SI / (2 * q_e * lambdaL_m * ballisticWperp_m)
    endif

    Variable/G v_B0_ballistic_T = ballisticB0_T
    Variable/G root:B0_ballistic_T = ballisticB0_T

    if (doDisplay)

        // ---------------- Topography ----------------
        SNS_DisplayWithScales(image, cmap="grayC")
        String winImage = WinName(0, 1, 1)

        SNS_StyleRayDisplayImage(image, winImage)

        SNS_PlotChannelRay(Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, v_Gap0)
        ModifyGraph rgb($StringFromList(ItemsInList(TraceNameList("", ";", 1))-1, TraceNameList("", ";", 1)))=(0,0,0)

        SNS_PlotChannelRay(Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, v_Longest)

        // ---------------- Textbox summary ----------------
        String txtBallistic, txtUsadel, txtBox
        Make/FREE/D/N=(numpnts(L_N_List_nm)) L_N_List
        L_N_List = L_N_List_nm

        txtBallistic = "Ballistic\r" + \
            "ℓ = " + num2str(round(L_N_List[v_Longest])) + " nm\r" + \
            "w = " + num2str(round(W_eff_List_nm[v_Gap0])) + " nm\r"

        if (numtype(ballisticB0_T) == 0)
            txtBallistic += "B\\B0\\M = " + num2str(round(1000*ballisticB0_T)/1000) + " T"
        else
            txtBallistic += "B\\B0\\M = n/a"
        endif

        // ---------------- Optional Usadel summary ----------------
        txtUsadel = ""

        NVAR/Z v_Usadel_Ldiff_nm
        NVAR/Z v_Usadel_xProbe_nm
        NVAR/Z v_Usadel_Lphi_use_nm
        NVAR/Z v_Usadel_Bpi_T

        if (NVAR_Exists(v_Usadel_Ldiff_nm) && NVAR_Exists(v_Usadel_xProbe_nm) && NVAR_Exists(v_Usadel_Lphi_use_nm) && NVAR_Exists(v_Usadel_Bpi_T))

            txtUsadel = "Usadel 1D\r" + \
                "L\\Bdiff\\M = " + num2str(round(v_Usadel_Ldiff_nm)) + " nm\r" + \
                "x\\Bprobe\\M = " + num2str(round(v_Usadel_xProbe_nm)) + " nm\r" + \
                "L\\Bφ\\M = " + num2str(round(v_Usadel_Lphi_use_nm)) + " nm\r" + \
                "B\\B0\\M = " + num2str(round(1000*v_Usadel_Bpi_T)/1000) + " T"

            NVAR/Z v_Usadel_D_cm2_per_s = v_Usadel_D_cm2_per_s
            NVAR/Z v_Usadel_D_cm2_s     = v_Usadel_D_cm2_s
            NVAR/Z v_Usadel_D_m2_s      = v_Usadel_D_m2_s

            if (NVAR_Exists(v_Usadel_D_cm2_per_s))
                txtUsadel += "\rD = " + num2str(round(100*v_Usadel_D_cm2_per_s)/100) + " cm\\S2\\M/s"
            elseif (NVAR_Exists(v_Usadel_D_cm2_s))
                txtUsadel += "\rD = " + num2str(round(100*v_Usadel_D_cm2_s)/100) + " cm\\S2\\M/s"
            elseif (NVAR_Exists(v_Usadel_D_m2_s))
                txtUsadel += "\rD = " + num2str(round(100*v_Usadel_D_m2_s*1e4)/100) + " cm\\S2\\M/s"
            else
                NVAR/Z lmfp_nm = root:v_cfg_lmfp_nm
                NVAR/Z vF_mps = root:SNS_Settings:vF

                if (NVAR_Exists(lmfp_nm) && NVAR_Exists(vF_mps))
                    Variable diffusionConstant_m2_per_s = 0.5 * vF_mps * lmfp_nm * 1e-9
                    Variable diffusionConstant_cm2_per_s = diffusionConstant_m2_per_s * 1e4

                    Variable/G v_Usadel_D_m2_per_s = diffusionConstant_m2_per_s
                    Variable/G v_Usadel_D_cm2_per_s = diffusionConstant_cm2_per_s

                    txtUsadel += "\rD = " + num2str(round(100*diffusionConstant_cm2_per_s)/100) + " cm\\S2\\M/s"
                endif
            endif

        endif

        txtBox = txtBallistic

        if (strlen(txtUsadel) > 0)
            txtBox += "\r\r" + txtUsadel
        endif

        TextBox/C/N=TextMode/X=0/Y=0/F=0/B=(r,g,b,a*0.6)/A=RT txtBox

        // ---------------- STS marker / setpoint / B arrow ----------------
        Make/O/N=1 tmpSTSx=STSx, tmpSTSy=STSy
        AppendToGraph tmpSTSy vs tmpSTSx
        ModifyGraph mode(tmpSTSy)=3, rgb(tmpSTSy)=(65535,65535,65535)

        SNS_TAGSetpoint()

        Variable xMin = DimOffset(image, 0)
        Variable xMax = xMin + DimDelta(image, 0) * (DimSize(image, 0) - 1)

        Variable yMin = DimOffset(image, 1)
        Variable yMax = yMin + DimDelta(image, 1) * (DimSize(image, 1) - 1)

        Variable pos_x = xMax - 0.15 * (xMax - xMin)
        Variable pos_y = yMin + 0.15 * (yMax - yMin)

        SNS_DrawImageArrow(image, Bangle_deg, "B", pos_x=pos_x, pos_y=pos_y, scale=1.5, color="black")

        GetWindow $WinName(0, 1, 1), wsize
        Variable wPx = abs(V_right - V_left)
        Variable hPx = abs(V_bottom - V_top)

        Variable maxDim = max(wPx, hPx)
        if (maxDim > 512)
            Variable ex = 512 / maxDim
            ModifyGraph/W=$WinName(0, 1, 1) expand=ex
        else
            ModifyGraph/W=$WinName(0, 1, 1) expand=1
        endif

        // ---------------- Mode property plot ----------------
        Display L_N_List_nm, W_eff_List_nm
        String winModes = WinName(0, 1, 1)

        AppendToGraph/R T_eff_List

        ModifyGraph rgb(L_N_List_nm)=(0,0,0)
        ModifyGraph rgb(W_eff_List_nm)=(1,16019,65535)
        ModifyGraph mode(L_N_List_nm)=2
        ModifyGraph mode(W_eff_List_nm)=2
        ModifyGraph lsize=5
        ModifyGraph mode(T_eff_List)=3
        ModifyGraph marker(T_eff_List)=8

        SetAxis right 0,1
        SetAxis left 0,*

        ModifyGraph muloffset(L_N_List_nm)={0,0}
        ModifyGraph muloffset(W_eff_List_nm)={0,0}
        ModifyGraph muloffset(T_eff_List)={0,0}
        ModifyGraph tick=2
        ModifyGraph mirror(bottom)=2
        ModifyGraph ZisZ=1
        ModifyGraph standoff(left)=0
        ModifyGraph standoff(bottom)=0

        Label left   "Length (nm)"
        Label bottom "Trajectory Nr."
        Label right  "Transparency"

        Legend/C/N=text0/J/F=0/A=RT "Mode\r" + \
            "\\s(L_N_List_nm) length (L)\r" + \
            "\\s(W_eff_List_nm) projection ⊥B\\Bext\\M (W)\r" + \
            "\\s(T_eff_List) transparency (T)"

        ModifyGraph margin(left)=33
        ModifyGraph margin(bottom)=30
        ModifyGraph margin(right)=35
        ModifyGraph margin(top)=6
        ModifyGraph width=220
        ModifyGraph height=220

    endif

    // ---------------- Restore folder ----------------
    SetDataFolder $oldDF

    return 0
End

//==============================================================================
// SNS_BuildModeWeightedChannelsFromRayWeights
//
// Purpose:
//   Import the Python normal-mode-to-ray projection table and build a new
//   canonical channel folder from an existing dense candidate ray folder.
//
//   This is meant for the low-mode diagnostic workflow:
//     1. build a dense candidate ray basis in Igor,
//     2. solve scalar hard-wall normal modes externally,
//     3. import normal_mode_ray_weights.txt,
//     4. run the usual LDOS calculation on the mode-weighted channel folder.
//
// Inputs:
//   candidateFolder : Igor data folder containing canonical dense ray waves:
//                       L_N_List_nm, W_eff_List_nm, wChan, T_eff_List,
//                       Hit1x/y_List_nm, Hit2x/y_List_nm, phiList_rad, w_mask.
//
//   exportFolder    : Igor path string for the folder containing
//                     normal_mode_ray_weights.txt.
//
//   outFolder       : output data folder for the selected/weighted channels.
//
// Optional:
//   topRanks        : keep ranks < topRanks. Default: 6. Use <=0 to keep all.
//   minWeight       : discard projection rows below this raw weight. Default: 0.
//   normalizePerMode: 1 default. Renormalize kept ray weights separately for
//                     each eigenmode so that each mode contributes one unit of
//                     channel weight across its retained rays.
//
// Outputs:
//   outFolder canonical channel waves plus:
//     ModeSolver_mode
//     ModeSolver_rank
//     ModeSolver_ray_index
//     ModeSolver_weight_raw
//     ModeSolver_weight_used
//
// Returns:
//   Number of generated weighted channels.
//==============================================================================
Function SNS_BuildModeWeightedChannelsFromRayWeights(candidateFolder, exportFolder, outFolder, [topRanks, minWeight, normalizePerMode])
    String candidateFolder
    String exportFolder
    String outFolder
    Variable topRanks
    Variable minWeight
    Variable normalizePerMode

    String savedDF = GetDataFolder(1)

    if (ParamIsDefault(topRanks))
        topRanks = 6
    endif
    if (ParamIsDefault(minWeight))
        minWeight = 0
    endif
    if (ParamIsDefault(normalizePerMode))
        normalizePerMode = 1
    endif

    DFREF dfrCand = $candidateFolder

    Wave Lcand   = dfrCand:L_N_List_nm
    Wave Wcand   = dfrCand:W_eff_List_nm
    Wave Tcand   = dfrCand:T_eff_List
    Wave X1cand  = dfrCand:Hit1x_List_nm
    Wave Y1cand  = dfrCand:Hit1y_List_nm
    Wave X2cand  = dfrCand:Hit2x_List_nm
    Wave Y2cand  = dfrCand:Hit2y_List_nm
    Wave Phicand = dfrCand:phiList_rad
    Wave/Z MaskCand = dfrCand:w_mask

    SetDataFolder dfrCand
    String candImageList = WaveList("*_Z_mbgnd_xy", ";", "DIMS:2")
    String candImageName = ""
    if (strlen(candImageList) > 0)
        candImageName = StringFromList(0, candImageList, ";")
    endif
    SetDataFolder $savedDF

    Variable nCand = numpnts(Lcand)
    if ((nCand <= 0) || \
        (nCand != numpnts(Wcand)) || \
        (nCand != numpnts(Tcand)) || \
        (nCand != numpnts(X1cand)) || \
        (nCand != numpnts(Y1cand)) || \
        (nCand != numpnts(X2cand)) || \
        (nCand != numpnts(Y2cand)) || \
        (nCand != numpnts(Phicand)))

        SetDataFolder $savedDF
        Abort "SNS_BuildModeWeightedChannelsFromRayWeights: candidate channel waves are inconsistent."
    endif
    if (!WaveExists(MaskCand))
        SetDataFolder $savedDF
        Abort "SNS_BuildModeWeightedChannelsFromRayWeights: candidate folder missing w_mask."
    endif

    NewPath/O/Q SNSModeRayWeightPath, exportFolder

    Variable refNum
    String line
    Variable modeVal, rankVal, rayVal, phiVal, weightVal
    Variable nRead, nAccepted = 0
    Variable maxMode = -1

    Open/R/P=SNSModeRayWeightPath refNum as "normal_mode_ray_weights.txt"
    FReadLine refNum, line
    do
        FReadLine refNum, line
        if (strlen(line) == 0)
            break
        endif
        sscanf line, "%g %g %g %g %g", modeVal, rankVal, rayVal, phiVal, weightVal
        nRead = V_flag
        if (nRead == 5)
            if (((topRanks <= 0) || (rankVal < topRanks)) && (weightVal >= minWeight))
                nAccepted += 1
                maxMode = max(maxMode, round(modeVal))
            endif
        endif
    while (1)
    Close refNum

    if ((nAccepted <= 0) || (maxMode < 0))
        SetDataFolder $savedDF
        Abort "SNS_BuildModeWeightedChannelsFromRayWeights: no usable rows in normal_mode_ray_weights.txt."
    endif

    Make/FREE/D/N=(nAccepted) tmpMode, tmpRank, tmpRay, tmpPhi, tmpWeight
    Make/FREE/D/N=(maxMode+1) weightSumByMode
    weightSumByMode = 0

    Variable i = 0
    Open/R/P=SNSModeRayWeightPath refNum as "normal_mode_ray_weights.txt"
    FReadLine refNum, line
    do
        FReadLine refNum, line
        if (strlen(line) == 0)
            break
        endif
        sscanf line, "%g %g %g %g %g", modeVal, rankVal, rayVal, phiVal, weightVal
        nRead = V_flag
        if (nRead == 5)
            if (((topRanks <= 0) || (rankVal < topRanks)) && (weightVal >= minWeight))
                tmpMode[i] = round(modeVal)
                tmpRank[i] = round(rankVal)
                tmpRay[i] = round(rayVal)
                tmpPhi[i] = phiVal
                tmpWeight[i] = weightVal
                weightSumByMode[tmpMode[i]] += weightVal
                i += 1
            endif
        endif
    while (1)
    Close refNum

    NewDataFolder/O $outFolder
    DFREF dfrOut = $outFolder
    SetDataFolder dfrOut

    Duplicate/O MaskCand, w_mask
    if (strlen(candImageName) > 0)
        SetDataFolder dfrCand
        Wave ImageCand = $candImageName
        SetDataFolder dfrOut
        Duplicate/O ImageCand, $candImageName
    endif

    Make/O/D/N=(nAccepted) L_N_List_nm, W_eff_List_nm, wChan, T_eff_List
    Make/O/D/N=(nAccepted) Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm
    Make/O/D/N=(nAccepted) phiList_rad
    Make/O/D/N=(nAccepted) ModeSolver_mode, ModeSolver_rank, ModeSolver_ray_index
    Make/O/D/N=(nAccepted) ModeSolver_weight_raw, ModeSolver_weight_used

    Variable idx, m, wUse
    for (i = 0; i < nAccepted; i += 1)
        idx = round(tmpRay[i])
        if ((idx < 0) || (idx >= nCand))
            SetDataFolder $savedDF
            Abort "SNS_BuildModeWeightedChannelsFromRayWeights: ray_index=" + num2str(idx) + " outside candidate range [0," + num2str(nCand-1) + "]. Re-export ray_candidates_* and rerun mask_normal_modes.py before importing weights."
        endif

        m = round(tmpMode[i])
        wUse = tmpWeight[i]
        if (normalizePerMode && (weightSumByMode[m] > 0))
            wUse = tmpWeight[i] / weightSumByMode[m]
        endif

        L_N_List_nm[i]    = Lcand[idx]
        W_eff_List_nm[i]  = Wcand[idx]
        T_eff_List[i]     = Tcand[idx]
        Hit1x_List_nm[i]  = X1cand[idx]
        Hit1y_List_nm[i]  = Y1cand[idx]
        Hit2x_List_nm[i]  = X2cand[idx]
        Hit2y_List_nm[i]  = Y2cand[idx]
        phiList_rad[i]    = Phicand[idx]
        wChan[i]          = wUse

        ModeSolver_mode[i]        = m
        ModeSolver_rank[i]        = tmpRank[i]
        ModeSolver_ray_index[i]   = idx
        ModeSolver_weight_raw[i]  = tmpWeight[i]
        ModeSolver_weight_used[i] = wUse
    endfor

    Variable/G v_ModeSolver_NModes = maxMode + 1
    Variable/G v_ModeSolver_NWeightedChannels = nAccepted
    Variable/G v_ModeSolver_TopRanks = topRanks
    Variable/G v_ModeSolver_MinWeight = minWeight
    Variable/G v_ModeSolver_NormalizePerMode = normalizePerMode

    SetDataFolder $savedDF
    return nAccepted
End
//==============================================================================
// SNS_Usadel1D_DisplayFromRayFolder
//
// Purpose:
//   Display Usadel-geometry diagnostics for one ray-tracing folder.
//
//   The function assumes that these helpers have already been run:
//
//      SNS_Usadel1D_GeometryFromRayFolder(...)
//      SNS_Usadel1D_LeverArmFromRayFolder(...)
//
//   It creates:
//
//     1) topography image with an effective 1D Usadel overlay:
//          - black line  : Ldiff effective 1D wire
//          - white cross : STS position
//          - cyan dashed : Lphi lever arm, perpendicular to B
//
//        The STS position is placed at xProbe along the Ldiff line,
//        measured from the nearer effective end.
//
//     2) textbox with effective Usadel parameters
//
//     3) mode-property plot showing chord length, near-edge distance,
//        fractional probe position, and optionally Lphi.
//
// Inputs:
//   dfPath      : folder containing ray-tracing and Usadel output waves.
//   Bangle_deg  : optional magnetic-field angle in degrees.
//                 Required for drawing Lphi direction and B arrow.
//   doImage     : optional. 1 default: show image diagnostic.
//   doModePlot  : optional. 1 default: show mode-property plot.
//
// Outputs:
//   return : numeric
//            0 on successful completion.
//
//==============================================================================

Function SNS_Usadel1D_DisplayFromRayFolder(dfPath, [Bangle_deg, doImage, doModePlot])
    String dfPath
    Variable Bangle_deg, doImage, doModePlot

    String oldDF = GetDataFolder(1)
    SetDataFolder $dfPath

    if (ParamIsDefault(doImage))
        doImage = 1
    endif
    if (ParamIsDefault(doModePlot))
        doModePlot = 1
    endif

    Variable textBoxRed = 65535
    Variable textBoxGreen = 65535
    Variable textBoxBlue = 65535
    Variable textBoxAlpha = 65535

    // ------------------------------
    // Required Usadel outputs
    // ------------------------------
    NVAR/Z v_Usadel_Ldiff_nm
    NVAR/Z v_Usadel_xProbe_nm
    NVAR/Z v_Usadel_fProbe
    NVAR/Z v_Usadel_Lmean_nm
    NVAR/Z v_Usadel_Lrms_nm
    NVAR/Z v_Usadel_Lmin_nm
    NVAR/Z v_Usadel_Lmax_nm
    NVAR/Z v_Usadel_dNearMean_nm

    if (!NVAR_Exists(v_Usadel_Ldiff_nm) || !NVAR_Exists(v_Usadel_xProbe_nm))
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_DisplayFromRayFolder: run SNS_Usadel1D_GeometryFromRayFolder first."
    endif

    NVAR/Z v_Usadel_Lphi_use_nm
    NVAR/Z v_Usadel_Bpi_T

    Wave/Z w_UsadelChord_nm
    Wave/Z w_UsadelDnear_nm
    Wave/Z w_UsadelFracNear
    Wave/Z w_UsadelWeight
    Wave/Z w_UsadelLphi_nm

    Wave/Z Hit1x_List
    Wave/Z Hit1y_List
    Wave/Z Hit2x_List
    Wave/Z Hit2y_List
    Wave/Z tmpSTSx
    Wave/Z tmpSTSy

    // ------------------------------
    // 1) Topography display
    // ------------------------------
    if (doImage)

        String imgList = WaveList("*_Z_mbgnd_xy", ";", "DIMS:2")
        if (strlen(imgList) == 0)
            imgList = WaveList("*_Z_xy", ";", "DIMS:2")
        endif

        if (strlen(imgList) > 0)

            String imgName = StringFromList(0, imgList, ";")
            Wave image = $imgName

            SNS_DisplayWithScales(image, cmap="grayC")
            String winImage = WinName(0, 1, 1)

            SNS_StyleRayDisplayImage(image, winImage)

            // ------------------------------
            // Effective 1D Usadel overlay
            // ------------------------------
            if (WaveExists(tmpSTSx) && WaveExists(tmpSTSy) && WaveExists(Hit1x_List) && WaveExists(Hit1y_List) && WaveExists(Hit2x_List) && WaveExists(Hit2y_List) && WaveExists(w_UsadelChord_nm))

                Variable STSx_axis = tmpSTSx[0]
                Variable STSy_axis = tmpSTSy[0]

                // Choose representative direction from the ray whose chord is
                // closest to Ldiff. This gives only the orientation of the 1D wire.
                Make/O/D/N=(numpnts(w_UsadelChord_nm)) w_UsadelDiffToLdiff
                w_UsadelDiffToLdiff = abs(w_UsadelChord_nm - v_Usadel_Ldiff_nm)
                WaveStats/Q w_UsadelDiffToLdiff
                Variable idxLdiff = V_minloc
                KillWaves/Z w_UsadelDiffToLdiff

                Variable hit1Distance_axis = sqrt((Hit1x_List[idxLdiff] - STSx_axis)^2 + (Hit1y_List[idxLdiff] - STSy_axis)^2)
                Variable hit2Distance_axis = sqrt((Hit2x_List[idxLdiff] - STSx_axis)^2 + (Hit2y_List[idxLdiff] - STSy_axis)^2)

                Variable nearHitX_axis
                Variable nearHitY_axis
                Variable farHitX_axis
                Variable farHitY_axis

                if (hit1Distance_axis <= hit2Distance_axis)
                    nearHitX_axis = Hit1x_List[idxLdiff]
                    nearHitY_axis = Hit1y_List[idxLdiff]
                    farHitX_axis  = Hit2x_List[idxLdiff]
                    farHitY_axis  = Hit2y_List[idxLdiff]
                else
                    nearHitX_axis = Hit2x_List[idxLdiff]
                    nearHitY_axis = Hit2y_List[idxLdiff]
                    farHitX_axis  = Hit1x_List[idxLdiff]
                    farHitY_axis  = Hit1y_List[idxLdiff]
                endif

                Variable representativeDirectionX_axis = farHitX_axis - nearHitX_axis
                Variable representativeDirectionY_axis = farHitY_axis - nearHitY_axis
                Variable representativeDirectionLength_axis = sqrt(representativeDirectionX_axis^2 + representativeDirectionY_axis^2)

                if (representativeDirectionLength_axis > 0)

                    Variable wireUnitX_axis = representativeDirectionX_axis / representativeDirectionLength_axis
                    Variable wireUnitY_axis = representativeDirectionY_axis / representativeDirectionLength_axis

                    Variable wirePerpendicularUnitX_axis = -wireUnitY_axis
                    Variable wirePerpendicularUnitY_axis =  wireUnitX_axis

                    // Effective Ldiff segment:
                    // STS position is xProbe from the near effective end.
                    Variable effectiveNearEndX_axis = STSx_axis - v_Usadel_xProbe_nm * wireUnitX_axis
                    Variable effectiveNearEndY_axis = STSy_axis - v_Usadel_xProbe_nm * wireUnitY_axis

                    Variable effectiveFarEndX_axis = STSx_axis + (v_Usadel_Ldiff_nm - v_Usadel_xProbe_nm) * wireUnitX_axis
                    Variable effectiveFarEndY_axis = STSy_axis + (v_Usadel_Ldiff_nm - v_Usadel_xProbe_nm) * wireUnitY_axis

                    // Draw Ldiff line.
                    SetDrawLayer/W=$winImage UserFront
                    SetDrawEnv/W=$winImage xcoord=bottom, ycoord=left, linefgc=(0,0,0), linethick=2.5
                    DrawLine/W=$winImage effectiveNearEndX_axis, effectiveNearEndY_axis, effectiveFarEndX_axis, effectiveFarEndY_axis

                    // End ticks.
                    Variable imageRangeX_axis = abs(DimDelta(image,0)) * max(1, DimSize(image,0)-1)
                    Variable imageRangeY_axis = abs(DimDelta(image,1)) * max(1, DimSize(image,1)-1)
                    Variable imageShortRange_axis = min(imageRangeX_axis, imageRangeY_axis)
                    Variable endTickLength_axis = 0.025 * imageShortRange_axis

                    SetDrawEnv/W=$winImage xcoord=bottom, ycoord=left, linefgc=(0,0,0), linethick=2
                    DrawLine/W=$winImage effectiveNearEndX_axis - 0.5*endTickLength_axis*wirePerpendicularUnitX_axis, effectiveNearEndY_axis - 0.5*endTickLength_axis*wirePerpendicularUnitY_axis, effectiveNearEndX_axis + 0.5*endTickLength_axis*wirePerpendicularUnitX_axis, effectiveNearEndY_axis + 0.5*endTickLength_axis*wirePerpendicularUnitY_axis
                    DrawLine/W=$winImage effectiveFarEndX_axis - 0.5*endTickLength_axis*wirePerpendicularUnitX_axis, effectiveFarEndY_axis - 0.5*endTickLength_axis*wirePerpendicularUnitY_axis, effectiveFarEndX_axis + 0.5*endTickLength_axis*wirePerpendicularUnitX_axis, effectiveFarEndY_axis + 0.5*endTickLength_axis*wirePerpendicularUnitY_axis

                    // STS marker as plotted cross.
                    Make/O/D/N=1 w_UsadelDraw_STSx = STSx_axis
                    Make/O/D/N=1 w_UsadelDraw_STSy = STSy_axis
                    AppendToGraph/W=$winImage w_UsadelDraw_STSy vs w_UsadelDraw_STSx
                    ModifyGraph/W=$winImage mode(w_UsadelDraw_STSy)=3, marker(w_UsadelDraw_STSy)=0, rgb(w_UsadelDraw_STSy)=(65535,65535,65535), msize(w_UsadelDraw_STSy)=5, mrkThick(w_UsadelDraw_STSy)=1.5

                    // Lphi lever arm: perpendicular to B, centered on STS.
                    if (!ParamIsDefault(Bangle_deg) && NVAR_Exists(v_Usadel_Lphi_use_nm))

                        Variable Bangle_rad = Bangle_deg * pi / 180
                        Variable lphiUnitX_axis = -sin(Bangle_rad)
                        Variable lphiUnitY_axis =  cos(Bangle_rad)

                        Variable lphiHalfLength_axis = 0.5 * v_Usadel_Lphi_use_nm

                        Variable lphiStartX_axis = STSx_axis - lphiHalfLength_axis * lphiUnitX_axis
                        Variable lphiStartY_axis = STSy_axis - lphiHalfLength_axis * lphiUnitY_axis

                        Variable lphiEndX_axis = STSx_axis + lphiHalfLength_axis * lphiUnitX_axis
                        Variable lphiEndY_axis = STSy_axis + lphiHalfLength_axis * lphiUnitY_axis

                        SetDrawEnv/W=$winImage xcoord=bottom, ycoord=left, linefgc=(1,16019,65535), linethick=2.5, dash=1
                        DrawLine/W=$winImage lphiStartX_axis, lphiStartY_axis, lphiEndX_axis, lphiEndY_axis

                        SetDrawEnv/W=$winImage xcoord=bottom, ycoord=left, textrgb=(1,16019,65535), fsize=12, textxjust=0, textyjust=1
                        DrawText/W=$winImage lphiEndX_axis, lphiEndY_axis, "L\\Bφ\\M"
                    endif

                    // Ldiff label.
                    SetDrawEnv/W=$winImage xcoord=bottom, ycoord=left, textrgb=(0,0,0), fsize=12, textxjust=0, textyjust=1
                    DrawText/W=$winImage 0.5*(effectiveNearEndX_axis + effectiveFarEndX_axis), 0.5*(effectiveNearEndY_axis + effectiveFarEndY_axis), "L\\Bdiff\\M"
                endif
            endif

            // ------------------------------
            // Text summary
            // ------------------------------
            String txtBox
            txtBox = "Usadel 1D geometry\r" + \
                "L\\Bdiff\\M = " + num2str(round(v_Usadel_Ldiff_nm)) + " nm\r" + \
                "x\\Bprobe\\M = " + num2str(round(v_Usadel_xProbe_nm)) + " nm\r" + \
                "f\\Bprobe\\M = " + num2str(round(1000*v_Usadel_fProbe)/1000) + "\r" + \
                "d\\Bnear\\M = " + num2str(round(v_Usadel_dNearMean_nm)) + " nm"

            if (NVAR_Exists(v_Usadel_Lphi_use_nm))
                txtBox += "\rL\\Bφ\\M = " + num2str(round(v_Usadel_Lphi_use_nm)) + " nm"
            endif

            if (NVAR_Exists(v_Usadel_Bpi_T))
                txtBox += "\rB\\Bπ\\M = " + num2str(round(1000*v_Usadel_Bpi_T)/1000) + " T"
            endif

            TextBox/W=$winImage/C/N=TextUsadel/X=0/Y=0/F=0/B=(textBoxRed,textBoxGreen,textBoxBlue,textBoxAlpha*0.6)/A=RT txtBox

            // B direction arrow.
            if (!ParamIsDefault(Bangle_deg))
                Variable imageXmin_axis = DimOffset(image, 0)
                Variable imageXmax_axis = imageXmin_axis + DimDelta(image, 0) * (DimSize(image, 0) - 1)

                Variable imageYmin_axis = DimOffset(image, 1)
                Variable imageYmax_axis = imageYmin_axis + DimDelta(image, 1) * (DimSize(image, 1) - 1)

                Variable BArrowX_axis = imageXmax_axis - 0.15 * (imageXmax_axis - imageXmin_axis)
                Variable BArrowY_axis = imageYmin_axis + 0.15 * (imageYmax_axis - imageYmin_axis)

                SNS_DrawImageArrow(image, Bangle_deg, "B", pos_x=BArrowX_axis, pos_y=BArrowY_axis, scale=1.5, color="White")
            endif

            // Scale image window similarly to SNS_ExtractModesForFolder.
            GetWindow $winImage, wsize
            Variable winWidthPx = abs(V_right - V_left)
            Variable winHeightPx = abs(V_bottom - V_top)
            Variable winMaxDimPx = max(winWidthPx, winHeightPx)

            if (winMaxDimPx > 512)
                ModifyGraph/W=$winImage expand=(512 / winMaxDimPx)
            else
                ModifyGraph/W=$winImage expand=1
            endif
        endif
    endif

    // ------------------------------
    // 2) Mode-property display
    // ------------------------------
    if (doModePlot)

        if (WaveExists(w_UsadelChord_nm) && WaveExists(w_UsadelDnear_nm) && WaveExists(w_UsadelFracNear))

            Display w_UsadelChord_nm, w_UsadelDnear_nm
            String winModes = WinName(0, 1, 1)

            AppendToGraph/R w_UsadelFracNear

            ModifyGraph/W=$winModes rgb(w_UsadelChord_nm)=(0,0,0)
            ModifyGraph/W=$winModes rgb(w_UsadelDnear_nm)=(1,16019,65535)
            ModifyGraph/W=$winModes mode(w_UsadelChord_nm)=2
            ModifyGraph/W=$winModes mode(w_UsadelDnear_nm)=2
            ModifyGraph/W=$winModes lsize=4

            ModifyGraph/W=$winModes mode(w_UsadelFracNear)=3
            ModifyGraph/W=$winModes marker(w_UsadelFracNear)=8
            ModifyGraph/W=$winModes rgb(w_UsadelFracNear)=(65535,0,0)

            if (WaveExists(w_UsadelLphi_nm))
                AppendToGraph/W=$winModes w_UsadelLphi_nm
                ModifyGraph/W=$winModes rgb(w_UsadelLphi_nm)=(32768,0,65535)
                ModifyGraph/W=$winModes mode(w_UsadelLphi_nm)=2
                ModifyGraph/W=$winModes lsize(w_UsadelLphi_nm)=3
            endif

            SetAxis/W=$winModes left 0,*
            SetAxis/W=$winModes right 0,1

            ModifyGraph/W=$winModes tick=2
            ModifyGraph/W=$winModes mirror(bottom)=2
            ModifyGraph/W=$winModes ZisZ=1
            ModifyGraph/W=$winModes standoff(left)=0
            ModifyGraph/W=$winModes standoff(bottom)=0

            Label/W=$winModes left   "Length (nm)"
            Label/W=$winModes bottom "Trajectory Nr."
            Label/W=$winModes right  "Fraction"

            String legText
            legText = "Usadel\r" + \
                "\\s(w_UsadelChord_nm) chord length\r" + \
                "\\s(w_UsadelDnear_nm) near-side distance\r" + \
                "\\s(w_UsadelFracNear) x\\Bprobe\\M / L"

            if (WaveExists(w_UsadelLphi_nm))
                legText += "\r\\s(w_UsadelLphi_nm) L\\Bφ\\M"
            endif

            Legend/W=$winModes/C/N=text0/J/F=0/A=RT legText

            ModifyGraph/W=$winModes margin(left)=40
            ModifyGraph/W=$winModes margin(bottom)=30
            ModifyGraph/W=$winModes margin(right)=40
            ModifyGraph/W=$winModes margin(top)=6
            ModifyGraph/W=$winModes width=240
            ModifyGraph/W=$winModes height=220

        endif
    endif

    SetDataFolder $oldDF
    return 0
End

//==============================================================================
// SNS_Usadel1D_DrawEffectiveGeometry
//
// Purpose:
//   Draw an accurate effective 1D Usadel representation on the top image graph.
//
//   The function draws:
//     1) Ldiff line segment, with the STS point located at xProbe from
//        the near effective end.
//     2) Lphi line segment, perpendicular to B, centered on the STS position.
//
// Inputs:
//   dfPath     : folder containing Usadel output variables and ray waves.
//   image      : image wave shown in the active graph.
//   Bangle_deg : magnetic-field direction in degrees.
//
// Outputs:
//   return : numeric
//            0 on successful completion.
//
//==============================================================================

Function SNS_Usadel1D_DrawEffectiveGeometry(dfPath, image, Bangle_deg)
    String dfPath
    Wave image
    Variable Bangle_deg

    String oldDF = GetDataFolder(1)
    SetDataFolder $dfPath

    String graphWindowName = WinName(0, 1, 1)
    if (strlen(graphWindowName) == 0)
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_DrawEffectiveGeometry: no graph window found."
    endif

    NVAR/Z v_Usadel_Ldiff_nm
    NVAR/Z v_Usadel_xProbe_nm
    NVAR/Z v_Usadel_Lphi_use_nm

    if (!NVAR_Exists(v_Usadel_Ldiff_nm) || !NVAR_Exists(v_Usadel_xProbe_nm))
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_DrawEffectiveGeometry: missing v_Usadel_Ldiff_nm or v_Usadel_xProbe_nm."
    endif

    Wave/Z tmpSTSx
    Wave/Z tmpSTSy
    if (!WaveExists(tmpSTSx) || !WaveExists(tmpSTSy))
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_DrawEffectiveGeometry: missing tmpSTSx/tmpSTSy."
    endif

    Wave/Z Hit1x_List
    Wave/Z Hit1y_List
    Wave/Z Hit2x_List
    Wave/Z Hit2y_List
    Wave/Z w_UsadelChord_nm

    if (!WaveExists(Hit1x_List) || !WaveExists(Hit1y_List) || !WaveExists(Hit2x_List) || !WaveExists(Hit2y_List) || !WaveExists(w_UsadelChord_nm))
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_DrawEffectiveGeometry: missing ray or Usadel chord waves."
    endif

    Variable STSx_axis = tmpSTSx[0]
    Variable STSy_axis = tmpSTSy[0]

    // ------------------------------
    // Choose a ray direction representative of Ldiff
    // ------------------------------
    Make/O/D/N=(numpnts(w_UsadelChord_nm)) w_UsadelDiffToLdiff
    w_UsadelDiffToLdiff = abs(w_UsadelChord_nm - v_Usadel_Ldiff_nm)
    WaveStats/Q w_UsadelDiffToLdiff
    Variable idxLdiff = V_minloc
    KillWaves/Z w_UsadelDiffToLdiff

    // Direction from near side to far side through STS.
    Variable hit1Distance = sqrt((Hit1x_List[idxLdiff] - STSx_axis)^2 + (Hit1y_List[idxLdiff] - STSy_axis)^2)
    Variable hit2Distance = sqrt((Hit2x_List[idxLdiff] - STSx_axis)^2 + (Hit2y_List[idxLdiff] - STSy_axis)^2)

    Variable nearHitX_axis
    Variable nearHitY_axis
    Variable farHitX_axis
    Variable farHitY_axis

    if (hit1Distance <= hit2Distance)
        nearHitX_axis = Hit1x_List[idxLdiff]
        nearHitY_axis = Hit1y_List[idxLdiff]
        farHitX_axis  = Hit2x_List[idxLdiff]
        farHitY_axis  = Hit2y_List[idxLdiff]
    else
        nearHitX_axis = Hit2x_List[idxLdiff]
        nearHitY_axis = Hit2y_List[idxLdiff]
        farHitX_axis  = Hit1x_List[idxLdiff]
        farHitY_axis  = Hit1y_List[idxLdiff]
    endif

    Variable representativeDirectionX_axis = farHitX_axis - nearHitX_axis
    Variable representativeDirectionY_axis = farHitY_axis - nearHitY_axis
    Variable representativeDirectionLength_axis = sqrt(representativeDirectionX_axis^2 + representativeDirectionY_axis^2)

    if (representativeDirectionLength_axis <= 0)
        SetDataFolder $oldDF
        Abort "SNS_Usadel1D_DrawEffectiveGeometry: invalid representative direction."
    endif

    Variable wireUnitX_axis = representativeDirectionX_axis / representativeDirectionLength_axis
    Variable wireUnitY_axis = representativeDirectionY_axis / representativeDirectionLength_axis

    // ------------------------------
    // Effective Ldiff segment
    // STS is located xProbe from the near end.
    // ------------------------------
    Variable Ldiff_axis = v_Usadel_Ldiff_nm
    Variable xProbe_axis = v_Usadel_xProbe_nm

    Variable effectiveNearEndX_axis = STSx_axis - xProbe_axis * wireUnitX_axis
    Variable effectiveNearEndY_axis = STSy_axis - xProbe_axis * wireUnitY_axis

    Variable effectiveFarEndX_axis = STSx_axis + (Ldiff_axis - xProbe_axis) * wireUnitX_axis
    Variable effectiveFarEndY_axis = STSy_axis + (Ldiff_axis - xProbe_axis) * wireUnitY_axis

    // ------------------------------
    // Effective Lphi segment
    // Lphi is perpendicular to B and centered on STS.
    // ------------------------------
    Variable Bangle_rad = Bangle_deg * pi / 180
    Variable lphiUnitX_axis = -sin(Bangle_rad)
    Variable lphiUnitY_axis =  cos(Bangle_rad)

    Variable Lphi_axis = 0
    if (NVAR_Exists(v_Usadel_Lphi_use_nm))
        Lphi_axis = v_Usadel_Lphi_use_nm
    endif

    Variable lphiHalf_axis = 0.5 * Lphi_axis

    Variable lphiStartX_axis = STSx_axis - lphiHalf_axis * lphiUnitX_axis
    Variable lphiStartY_axis = STSy_axis - lphiHalf_axis * lphiUnitY_axis

    Variable lphiEndX_axis = STSx_axis + lphiHalf_axis * lphiUnitX_axis
    Variable lphiEndY_axis = STSy_axis + lphiHalf_axis * lphiUnitY_axis

    // ------------------------------
    // Draw
    // ------------------------------
    SetDrawLayer/W=$graphWindowName UserFront

    // Ldiff: black thick line
    SetDrawEnv/W=$graphWindowName xcoord=bottom, ycoord=left, linefgc=(0,0,0), linethick=2.5
    DrawLine/W=$graphWindowName effectiveNearEndX_axis, effectiveNearEndY_axis, effectiveFarEndX_axis, effectiveFarEndY_axis

    // Near/far end ticks
    Variable tickLength_axis = 0.025 * min(abs(DimDelta(image,0))*(DimSize(image,0)-1), abs(DimDelta(image,1))*(DimSize(image,1)-1))
    Variable wirePerpX_axis = -wireUnitY_axis
    Variable wirePerpY_axis =  wireUnitX_axis

    SetDrawEnv/W=$graphWindowName xcoord=bottom, ycoord=left, linefgc=(0,0,0), linethick=2
    DrawLine/W=$graphWindowName effectiveNearEndX_axis - 0.5*tickLength_axis*wirePerpX_axis, effectiveNearEndY_axis - 0.5*tickLength_axis*wirePerpY_axis, effectiveNearEndX_axis + 0.5*tickLength_axis*wirePerpX_axis, effectiveNearEndY_axis + 0.5*tickLength_axis*wirePerpY_axis
    DrawLine/W=$graphWindowName effectiveFarEndX_axis - 0.5*tickLength_axis*wirePerpX_axis, effectiveFarEndY_axis - 0.5*tickLength_axis*wirePerpY_axis, effectiveFarEndX_axis + 0.5*tickLength_axis*wirePerpX_axis, effectiveFarEndY_axis + 0.5*tickLength_axis*wirePerpY_axis

    // Lphi: cyan line perpendicular to B
    if (Lphi_axis > 0)
        SetDrawEnv/W=$graphWindowName xcoord=bottom, ycoord=left, linefgc=(1,16019,65535), linethick=2.5, dash=1
        DrawLine/W=$graphWindowName lphiStartX_axis, lphiStartY_axis, lphiEndX_axis, lphiEndY_axis
    endif

    // STS marker: white cross
    Make/O/D/N=1 w_UsadelDraw_STSx = STSx_axis
    Make/O/D/N=1 w_UsadelDraw_STSy = STSy_axis
    AppendToGraph/W=$graphWindowName w_UsadelDraw_STSy vs w_UsadelDraw_STSx
    ModifyGraph/W=$graphWindowName mode(w_UsadelDraw_STSy)=3, marker(w_UsadelDraw_STSy)=0, rgb(w_UsadelDraw_STSy)=(65535,65535,65535), msize(w_UsadelDraw_STSy)=5, mrkThick(w_UsadelDraw_STSy)=1.5

    // Labels
    SetDrawEnv/W=$graphWindowName xcoord=bottom, ycoord=left, textrgb=(0,0,0), fsize=12, textxjust=0, textyjust=1
    DrawText/W=$graphWindowName 0.5*(effectiveNearEndX_axis + effectiveFarEndX_axis), 0.5*(effectiveNearEndY_axis + effectiveFarEndY_axis), "L\\Bdiff\\M"

    if (Lphi_axis > 0)
        SetDrawEnv/W=$graphWindowName xcoord=bottom, ycoord=left, textrgb=(1,16019,65535), fsize=12, textxjust=0, textyjust=1
        DrawText/W=$graphWindowName lphiEndX_axis, lphiEndY_axis, "L\\Bφ\\M"
    endif

    SetDataFolder $oldDF
    return 0
End


//==============================================================================
// SNS_ExtractRectModesForFolder
//
// Purpose:
//   Build tier-0 analytic rectangle-mode ray channels for a mask folder and
//   display them with the same diagnostic view used by SNS_ExtractModesForFolder.
//
//   This wrapper exists because SNS_ExtractModesForFolder intentionally builds
//   the ordinary equally-spaced angular ray ladder through
//   SNS_BuildChannelsFromMask2D().  Calling it for a rectangle-mode channel
//   folder would display the wrong rays.  This function instead calls
//   SNS_BuildChannelsFromRectModes2D() and then plots representative rays from
//   the resulting rectangle-mode channel waves.
//
// Inputs:
//   dfPath     : data folder containing w_mask and a *_Z_mbgnd_xy image.
//
// Optional:
//   Bangle_deg          : in-plane magnetic-field angle [deg]. Default 225.
//   STSx, STSy          : spectroscopy position [nm]. Defaults -4, -50.
//   DoDisplay           : 1 default. Show image, selected rays, and L/W plot.
//   nModes              : number of rectangle modes closest to kF. Default
//                         SNS_BuildChannelsFromRectModes2D default.
//   useTipWeight        : passed to SNS_BuildChannelsFromRectModes2D.
//   includeMirrorAngles : passed to SNS_BuildChannelsFromRectModes2D.
//
// Outputs:
//   Canonical channel waves in dfPath, plus rectangle-mode diagnostics.
//
// Returns:
//   Return code from SNS_BuildChannelsFromRectModes2D().
//==============================================================================
Function SNS_ExtractRectModesForFolder(dfPath, [Bangle_deg, STSx, STSy, DoDisplay, nModes, useTipWeight, includeMirrorAngles])
    String dfPath
    Variable Bangle_deg, STSx, STSy, DoDisplay
    Variable nModes, useTipWeight, includeMirrorAngles

    String oldDF = GetDataFolder(1)
    SetDataFolder $dfPath

    Variable r=65535, g=65535, b=65535, a=65535

    if (ParamIsDefault(Bangle_deg))
        Bangle_deg = 225
    endif
    if (ParamIsDefault(STSx))
        STSx = -4
    endif
    if (ParamIsDefault(STSy))
        STSy = -50
    endif
    if (ParamIsDefault(DoDisplay))
        DoDisplay = 1
    endif

    Wave/Z w_mask
    if (!WaveExists(w_mask))
        SetDataFolder $oldDF
        Abort "SNS_ExtractRectModesForFolder: missing w_mask."
    endif

    String imgList = WaveList("*_Z_mbgnd_xy",";","DIMS:2")
    if (strlen(imgList)==0)
        SetDataFolder $oldDF
        Abort "SNS_ExtractRectModesForFolder: no *_Z_mbgnd_xy image found in folder."
    endif

    String imgName = StringFromList(0, imgList, ";")
    Wave image = $imgName

    Variable v_Bangle = Bangle_deg * pi/180
    Variable err

    if (ParamIsDefault(nModes) && ParamIsDefault(useTipWeight) && ParamIsDefault(includeMirrorAngles))
        err = SNS_BuildChannelsFromRectModes2D(w_mask, STSx, STSy, v_Bangle, 0.5, "")
    elseif (ParamIsDefault(useTipWeight) && ParamIsDefault(includeMirrorAngles))
        err = SNS_BuildChannelsFromRectModes2D(w_mask, STSx, STSy, v_Bangle, 0.5, "", nModes=nModes)
    elseif (ParamIsDefault(nModes) && ParamIsDefault(includeMirrorAngles))
        err = SNS_BuildChannelsFromRectModes2D(w_mask, STSx, STSy, v_Bangle, 0.5, "", useTipWeight=useTipWeight)
    elseif (ParamIsDefault(nModes) && ParamIsDefault(useTipWeight))
        err = SNS_BuildChannelsFromRectModes2D(w_mask, STSx, STSy, v_Bangle, 0.5, "", includeMirrorAngles=includeMirrorAngles)
    elseif (ParamIsDefault(includeMirrorAngles))
        err = SNS_BuildChannelsFromRectModes2D(w_mask, STSx, STSy, v_Bangle, 0.5, "", nModes=nModes, useTipWeight=useTipWeight)
    elseif (ParamIsDefault(useTipWeight))
        err = SNS_BuildChannelsFromRectModes2D(w_mask, STSx, STSy, v_Bangle, 0.5, "", nModes=nModes, includeMirrorAngles=includeMirrorAngles)
    elseif (ParamIsDefault(nModes))
        err = SNS_BuildChannelsFromRectModes2D(w_mask, STSx, STSy, v_Bangle, 0.5, "", useTipWeight=useTipWeight, includeMirrorAngles=includeMirrorAngles)
    else
        err = SNS_BuildChannelsFromRectModes2D(w_mask, STSx, STSy, v_Bangle, 0.5, "", nModes=nModes, useTipWeight=useTipWeight, includeMirrorAngles=includeMirrorAngles)
    endif

    if (err != 0)
        SetDataFolder $oldDF
        return err
    endif

    Wave L_N_List_nm
    Wave W_eff_List_nm
    Wave Hit1x_List_nm
    Wave Hit1y_List_nm
    Wave Hit2x_List_nm
    Wave Hit2y_List_nm
    Wave T_eff_List

    SNS_MaskAreaPerim_FromParticles(w_mask)
    Wave w

    Make/O/N=1 w_area_nm2 = w[0]
    Make/O/N=1 w_perim_nm = w[1]

    Make/O/D/N=(DimSize(L_N_List_nm,0)) A_eff_nm2 = L_N_List_nm * W_eff_List_nm

    WaveStats/Q W_eff_List_nm
    Variable/G v_Gap0 = V_maxloc

    WaveStats/Q L_N_List_nm
    Variable/G v_Longest = V_maxloc

    Variable/G root:LengthMax_nm = L_N_List_nm[v_Longest]
    Variable/G root:WPerp_nm     = W_eff_List_nm[v_Gap0]

    Variable ballisticWperp_m = abs(W_eff_List_nm[v_Gap0]) * 1e-9
    Variable ballisticB0_T = NaN

    NVAR/Z lambdaL_m = root:SNS_Settings:lambdaL
    if (NVAR_Exists(lambdaL_m) && ballisticWperp_m > 0)
        ballisticB0_T = pi * HBAR_SI / (2 * q_e * lambdaL_m * ballisticWperp_m)
    endif

    Variable/G v_B0_ballistic_T = ballisticB0_T
    Variable/G root:B0_ballistic_T = ballisticB0_T

    if (DoDisplay)
        SNS_DisplayWithScales(image, cmap="grayC")
        String winImage = WinName(0, 1, 1)

        SNS_StyleRayDisplayImage(image, winImage)

        SNS_PlotChannelRay(Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, v_Gap0)
        ModifyGraph rgb($StringFromList(ItemsInList(TraceNameList("", ";", 1))-1, TraceNameList("", ";", 1)))=(0,0,0)

        SNS_PlotChannelRay(Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, v_Longest)

        String txtBallistic, txtBox
        Make/FREE/D/N=(numpnts(L_N_List_nm)) L_N_List
        L_N_List = L_N_List_nm

        txtBallistic = "Rect modes\r" + \
            "ell = " + num2str(round(L_N_List[v_Longest])) + " nm\r" + \
            "w = " + num2str(round(W_eff_List_nm[v_Gap0])) + " nm\r"

        if (numtype(ballisticB0_T) == 0)
            txtBallistic += "B\\B0\\M = " + num2str(round(1000*ballisticB0_T)/1000) + " T"
        else
            txtBallistic += "B\\B0\\M = n/a"
        endif

        txtBox = txtBallistic
        TextBox/C/N=TextMode/X=0/Y=0/F=0/B=(r,g,b,a*0.6)/A=RT txtBox

        Make/O/N=1 tmpSTSx=STSx, tmpSTSy=STSy
        AppendToGraph tmpSTSy vs tmpSTSx
        ModifyGraph mode(tmpSTSy)=3, rgb(tmpSTSy)=(65535,65535,65535)

        SNS_TAGSetpoint()

        Variable xMin = DimOffset(image, 0)
        Variable xMax = xMin + DimDelta(image, 0) * (DimSize(image, 0) - 1)
        Variable yMin = DimOffset(image, 1)
        Variable yMax = yMin + DimDelta(image, 1) * (DimSize(image, 1) - 1)

        Variable pos_x = xMax - 0.15 * (xMax - xMin)
        Variable pos_y = yMin + 0.15 * (yMax - yMin)

        SNS_DrawImageArrow(image, Bangle_deg, "B", pos_x=pos_x, pos_y=pos_y, scale=1.5, color="black")

        GetWindow $WinName(0, 1, 1), wsize
        Variable wPx = abs(V_right - V_left)
        Variable hPx = abs(V_bottom - V_top)
        Variable maxDim = max(wPx, hPx)
        if (maxDim > 512)
            Variable ex = 512 / maxDim
            ModifyGraph/W=$WinName(0, 1, 1) expand=ex
        else
            ModifyGraph/W=$WinName(0, 1, 1) expand=1
        endif

        Display L_N_List_nm, W_eff_List_nm
        String winModes = WinName(0, 1, 1)

        AppendToGraph/W=$winModes/R T_eff_List

        ModifyGraph/W=$winModes rgb(L_N_List_nm)=(0,0,0)
        ModifyGraph/W=$winModes rgb(W_eff_List_nm)=(1,16019,65535)
        ModifyGraph/W=$winModes mode(L_N_List_nm)=2
        ModifyGraph/W=$winModes mode(W_eff_List_nm)=2
        ModifyGraph/W=$winModes lsize=5
        ModifyGraph/W=$winModes mode(T_eff_List)=3
        ModifyGraph/W=$winModes marker(T_eff_List)=8

        SetAxis/W=$winModes right 0,1
        SetAxis/W=$winModes left 0,*

        ModifyGraph/W=$winModes muloffset(L_N_List_nm)={0,0}
        ModifyGraph/W=$winModes muloffset(W_eff_List_nm)={0,0}
        ModifyGraph/W=$winModes muloffset(T_eff_List)={0,0}
        ModifyGraph/W=$winModes tick=2
        ModifyGraph/W=$winModes mirror(bottom)=2
        ModifyGraph/W=$winModes ZisZ=1
        ModifyGraph/W=$winModes standoff(left)=0
        ModifyGraph/W=$winModes standoff(bottom)=0

        Label/W=$winModes left   "Length (nm)"
        Label/W=$winModes bottom "Trajectory Nr."
        Label/W=$winModes right  "Transparency"

        Legend/W=$winModes/C/N=text0/J/F=0/A=RT "Rect modes\r" + \
            "\\s(L_N_List_nm) length (L)\r" + \
            "\\s(W_eff_List_nm) projection perp B\\Bext\\M (W)\r" + \
            "\\s(T_eff_List) transparency (T)"

        ModifyGraph/W=$winModes margin(left)=33
        ModifyGraph/W=$winModes margin(bottom)=30
        ModifyGraph/W=$winModes margin(right)=35
        ModifyGraph/W=$winModes margin(top)=6
        ModifyGraph/W=$winModes width=220
        ModifyGraph/W=$winModes height=220
    endif

    SetDataFolder $oldDF
    return err
End


//==============================================================================
// SNS_IntegrateQSegmentNearest2D
//
// Purpose:
//   Integrate Q·dl along one projected 2D segment using midpoint sampling.
//
// Units:
//   Qx,Qy : rad/nm
//   x,y   : nm
//   return: rad
//==============================================================================
Function SNS_IntegrateQSegmentNearest2D(Qx, Qy, x1, y1, x2, y2, nStep)
    Wave Qx, Qy
    Variable x1, y1, x2, y2
    Variable nStep

    if (nStep < 1)
        nStep = 1
    endif

    Variable nx = DimSize(Qx, 0)
    Variable ny = DimSize(Qx, 1)

    Variable x0 = DimOffset(Qx, 0)
    Variable dx = DimDelta(Qx, 0)
    Variable y0 = DimOffset(Qx, 1)
    Variable dy = DimDelta(Qx, 1)

    Variable dlx = (x2 - x1) / nStep
    Variable dly = (y2 - y1) / nStep

    Variable beta = 0
    Variable k, t, x, y, ix, iy
    Variable qxVal, qyVal

    for (k = 0; k < nStep; k += 1)

        t = (k + 0.5) / nStep
        x = x1 + t*(x2 - x1)
        y = y1 + t*(y2 - y1)

        ix = round((x - x0) / dx)
        iy = round((y - y0) / dy)

        if ((ix < 0) || (ix >= nx) || (iy < 0) || (iy >= ny))
            continue
        endif

        qxVal = Qx[ix][iy]
        qyVal = Qy[ix][iy]

        if (numtype(qxVal) || numtype(qyVal))
            continue
        endif

        beta += qxVal*dlx + qyVal*dly
    endfor

    return beta
End



//==============================================================================
// SNS_BuildChannelsFromMask3D
//
// Purpose:
//   Construct a 3D SNS channel ensemble for a single STS position on the
//   top N–I interface, z = H_nm, using the shared 3D Cu(111)-faceted tracer.
//
// Interface transparency:
//   The interface model is controlled globally through SNS_Settings and
//   SNS_ChannelTransmissionFromCos(). This function no longer exposes a
//   Zbarrier / BTK_barrier input.
//
// Magnetic-area handling:
//   The tracer returns:
//      res[10] = signed Aeff_3D_nm2
//      res[11] = W_geom_nm
//
//   The existing solver-facing W_eff_List is kept, but is now the
//   solver-equivalent magnetic width:
//
//      W_eff_List = |Aeff_3D| / lambdaL_legacy
//
//   in meters. Therefore existing DOS helpers that compute the B-dependent
//   orbital phase from lambdaL * W_eff_List continue to work.
//
//   The actual geometric projection is preserved separately as:
//
//      W_geom_List [m]
//
//   and the actual signed area is stored as:
//
//      Aeff_3D_List_m2 [m^2]
//
// Optional Q-field mode:
//   If QxPhase and QyPhase are supplied, the builder additionally computes
//
//       betaExtra_List[ch] = betaLeg2 - betaLeg1
//
//   from the actual reflected trajectory, as returned by the shared tracer.
//
// Inputs:
//   Nmask      : 2D mask defining the N footprint; axes scaled in nm.
//   r0x, r0y   : STS position [nm]
//   phiB       : in-plane magnetic-field angle [rad]
//   stepFac    : ray/path integration step factor
//   H_nm       : N-region height [nm]
//   maxPath_nm : maximum allowed 3D trajectory length [nm]
//   folder     : output folder for channel waves
//
// Optional Inputs:
//   QxPhase, QyPhase : Q-field waves for betaExtra integration
//   qNstep           : Q integration step subdivision
//   qCoreHandling    : vortex-core handling mode
//   qMinUsedFrac     : minimum accepted Q integration fraction
//   xV_nm, yV_nm     : vortex center [nm]
//   rCore_nm         : vortex-core radius [nm]
//   PhaseReCore,
//   PhaseImCore      : phase field used for qCoreHandling=4
//
// Return values:
//    0 : success
//   -1 : invalid mask
//   -2 : no valid channels
//
// Outputs:
//   phiList, thetaList
//   L_N_List              [m]
//   W_eff_List            [m], solver-equivalent magnetic width
//   W_geom_List           [m], geometric |(S2-S1) dot e_perp|
//   Aeff_3D_List_m2       [m^2], signed effective magnetic area
//   wChan, T_eff_List
//   Hit1x_List, Hit1y_List [nm]
//   Hit2x_List, Hit2y_List [nm]
//   betaExtra_List         [rad], if Q mode is active
//==============================================================================
// Optional:
//   angularOffsetFrac : fractional offset of the azimuthal Fermi-grid spacing.
//                       offset angle = angularOffsetFrac * (pi/Nphi).
//                       Only the azimuthal phi grid is shifted; theta is unchanged.
//                       Default 0 preserves the historical midpoint grid.
Function SNS_BuildChannelsFromMask3D(Nmask, r0x, r0y, phiB, stepFac, H_nm, maxPath_nm, folder, [QxPhase, QyPhase, qNstep, qCoreHandling, qMinUsedFrac, xV_nm, yV_nm, rCore_nm, PhaseReCore, PhaseImCore, angularOffsetFrac])

    Wave    Nmask
    Variable r0x, r0y
    Variable phiB
    Variable stepFac
    Variable H_nm
    Variable maxPath_nm
    String  folder

    Wave QxPhase, QyPhase
    Variable qNstep, qCoreHandling, qMinUsedFrac
    Variable xV_nm, yV_nm, rCore_nm
    Wave PhaseReCore, PhaseImCore
    Variable angularOffsetFrac

    Variable haveQ = 0
    if (!ParamIsDefault(QxPhase) && !ParamIsDefault(QyPhase))
        haveQ = 1
    elseif (!ParamIsDefault(QxPhase) || !ParamIsDefault(QyPhase))
        Abort "SNS_BuildChannelsFromMask3D: supply both QxPhase and QyPhase, or neither."
    endif

    if (ParamIsDefault(qNstep))
        qNstep = 1
    endif
    qNstep = max(1, round(qNstep))

    if (ParamIsDefault(qCoreHandling))
        qCoreHandling = 0
    endif
    qCoreHandling = round(qCoreHandling)

    if ((qCoreHandling < 0) || (qCoreHandling > 4))
        Abort "SNS_BuildChannelsFromMask3D: qCoreHandling must be 0, 1, 2, 3, or 4."
    endif

    if (ParamIsDefault(qMinUsedFrac))
        qMinUsedFrac = 0.95
    endif

    if ((qMinUsedFrac <= 0) || (qMinUsedFrac > 1))
        Abort "SNS_BuildChannelsFromMask3D: qMinUsedFrac must be in (0,1]."
    endif

    if ((qCoreHandling > 0) && (ParamIsDefault(xV_nm) || ParamIsDefault(yV_nm) || ParamIsDefault(rCore_nm)))
        Abort "SNS_BuildChannelsFromMask3D: qCoreHandling>0 requires xV_nm, yV_nm, and rCore_nm."
    endif

    if ((qCoreHandling > 0) && !(rCore_nm > 0))
        Abort "SNS_BuildChannelsFromMask3D: rCore_nm must be > 0."
    endif

    if (qCoreHandling == 4)
        if (ParamIsDefault(PhaseReCore) || ParamIsDefault(PhaseImCore))
            Abort "SNS_BuildChannelsFromMask3D: qCoreHandling=4 requires PhaseReCore and PhaseImCore."
        endif
    endif

    if (haveQ)
        if ((DimSize(QxPhase,0) != DimSize(QyPhase,0)) || (DimSize(QxPhase,1) != DimSize(QyPhase,1)))
            Abort "SNS_BuildChannelsFromMask3D: QxPhase/QyPhase dimensions differ."
        endif
    endif

	 STRUCT SNS_Params params
	 SNS_LoadParams(params)
	
	 Variable lambdaF = params.LambdaF
	
	 // SNS_TraceOneChannel3D_Cu111 still receives the configured barrier
	 // internally; public callers do not pass Zbarrier.
	 Variable Zbarrier = params.BTK_barrier

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if (nx <= 1 || ny <= 1)
        return -1
    endif

    Variable lambdaF_nm = lambdaF * 1e9
    Variable Ndir = SNS_EstimateNdir_FromMask3D(Nmask, lambdaF_nm, H_nm)
    if (Ndir <= 0)
        Ndir = 8
    endif
    Ndir = max(Ndir, 8)
    Ndir = min(Ndir, 1024)

    Variable Ntheta = max(4, round(sqrt(Ndir)))
    Variable Nphi   = max(8, round(Ndir / Ntheta))
    Variable NchanRaw = Ntheta * Nphi
    Variable angularOffsetFracLocal
    if (ParamIsDefault(angularOffsetFrac) || numtype(angularOffsetFrac) != 0)
        angularOffsetFracLocal = 0
    else
        angularOffsetFracLocal = angularOffsetFrac - floor(angularOffsetFrac)
        if (angularOffsetFracLocal < 0)
            angularOffsetFracLocal += 1
        endif
    endif

    Make/FREE/D/N=(NchanRaw) phi_raw, theta_raw
    Make/FREE/D/N=(NchanRaw) L_N_raw, W_eff_raw, W_geom_raw, Aeff_raw_nm2
    Make/FREE/D/N=(NchanRaw) w_raw, T_raw
    Make/FREE/D/N=(NchanRaw) hit1x_raw, hit1y_raw, hit2x_raw, hit2y_raw
    Make/FREE/D/N=(NchanRaw) betaExtra_raw

    phi_raw = NaN
    theta_raw = NaN
    L_N_raw = NaN
    W_eff_raw = NaN
    W_geom_raw = NaN
    Aeff_raw_nm2 = NaN
    w_raw = 0
    T_raw = NaN
    hit1x_raw = NaN
    hit1y_raw = NaN
    hit2x_raw = NaN
    hit2y_raw = NaN
    betaExtra_raw = NaN

    String oldDF = GetDataFolder(1)
    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    KillWaves/Z phiList, thetaList, L_N_List, W_eff_List, W_geom_List, Aeff_3D_List_m2
    KillWaves/Z wChan, T_eff_List
    KillWaves/Z Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
    KillWaves/Z betaExtra_List
    KillWaves/Z trace3D_res_tmp

    Variable it, ip, j
    j = 0

    for (it = 0; it < Ntheta; it += 1)

        Variable theta = (it + 0.5)*(pi/2)/Ntheta

        for (ip = 0; ip < Nphi; ip += 1)

            Variable phi = mod((ip + 0.5 + angularOffsetFracLocal)*pi/Nphi, pi)

            phi_raw[j] = phi
            theta_raw[j] = theta

            Variable err
            if (haveQ)
                err = SNS_TraceOneChannel3D_Cu111(Nmask, r0x, r0y, theta, phi, phiB, stepFac, H_nm, maxPath_nm, "trace3D_res_tmp", 0, "", "", "", QxPhase=QxPhase, QyPhase=QyPhase, qNstep=qNstep, qCoreHandling=qCoreHandling, qMinUsedFrac=qMinUsedFrac, xV_nm=xV_nm, yV_nm=yV_nm, rCore_nm=rCore_nm, PhaseReCore=PhaseReCore, PhaseImCore=PhaseImCore)
            else
                err = SNS_TraceOneChannel3D_Cu111(Nmask, r0x, r0y, theta, phi, phiB, stepFac, H_nm, maxPath_nm, "trace3D_res_tmp", 0, "", "", "")
            endif

            if (err == 0)
                Wave tr = trace3D_res_tmp

                if (tr[0] > 0 && tr[9] > 0)
                    L_N_raw[j] = tr[1]
                    W_eff_raw[j] = tr[2]
                    hit1x_raw[j] = tr[3]
                    hit1y_raw[j] = tr[4]
                    hit2x_raw[j] = tr[5]
                    hit2y_raw[j] = tr[6]
                    T_raw[j] = tr[7]
                    betaExtra_raw[j] = tr[8]
                    w_raw[j] = tr[9]

                    if (numpnts(tr) > 10)
                        Aeff_raw_nm2[j] = tr[10]
                    endif
                    if (numpnts(tr) > 11)
                        W_geom_raw[j] = tr[11]
                    else
                        W_geom_raw[j] = tr[2]
                    endif
                endif
            endif

            j += 1
        endfor
    endfor

    Variable nRaw = NchanRaw
    Variable nValid = 0
    Variable k

    for (j = 0; j < nRaw; j += 1)
        if (w_raw[j] > 0)
            nValid += 1
        endif
    endfor

    if (nValid == 0)
        KillWaves/Z trace3D_res_tmp
        SetDataFolder $oldDF
        return -2
    endif

    Make/O/D/N=(nValid) phiList_rad, thetaList_rad
    Make/O/D/N=(nValid) L_N_List_nm, W_eff_List_nm, W_geom_List_nm, Aeff_3D_List_nm2
    Make/O/D/N=(nValid) wChan, T_eff_List
    Make/O/D/N=(nValid) Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm

    if (haveQ)
        Make/O/D/N=(nValid) betaExtra_List
    endif

    k = 0
    for (j = 0; j < nRaw; j += 1)
        if (w_raw[j] > 0)

            phiList_rad[k] = phi_raw[j]
            thetaList_rad[k] = theta_raw[j]

            // Geometry outputs stay in nm; solver-facing conversions happen downstream.
            L_N_List_nm[k] = L_N_raw[j]

            // Solver-equivalent magnetic width.
            W_eff_List_nm[k] = W_eff_raw[j]

            // Diagnostic geometric endpoint projection.
            W_geom_List_nm[k] = W_geom_raw[j]

            // Signed effective magnetic area.
            Aeff_3D_List_nm2[k] = Aeff_raw_nm2[j]

            wChan[k] = w_raw[j]
            T_eff_List[k] = T_raw[j]

            Hit1x_List_nm[k] = hit1x_raw[j]
            Hit1y_List_nm[k] = hit1y_raw[j]
            Hit2x_List_nm[k] = hit2x_raw[j]
            Hit2y_List_nm[k] = hit2y_raw[j]

            if (haveQ)
                betaExtra_List[k] = betaExtra_raw[j]
            endif

            k += 1
        endif
    endfor

    KillWaves/Z trace3D_res_tmp

    SetDataFolder $oldDF
    return 0
End

//==============================================================================
// SNS_IntegrateQSegmentNearest2D_CorePhase
//
// Helper for SNS_BuildChannelsFromMask3D.
//
// Integrates Q · dl over one projected 2D path segment [nm].
//
// coreHandling:
//   0 = ordinary nearest-neighbor Q integration
//   1 = project samples inside rCore_nm to core boundary
//   2 = reject segment if any sample enters core
//   3 = skip core samples, no correction
//   4 = geometrically split segment at core boundary and fill missing segment
//       using PhaseReCore/PhaseImCore:
//           beta_core = phase(exit) - phase(entry), unwrapped to [-pi, pi]
//
// Returns:
//   integrated phase segment [rad], or NaN if rejected/failed.
//
//==============================================================================
Function SNS_IntegrateQSegmentNearest2D_CorePhase(Qx, Qy, xA, yA, xB, yB, nStep, coreHandling, xV_nm, yV_nm, rCore_nm, minUsedFrac, [PhaseReCore, PhaseImCore])

    Wave Qx, Qy
    Variable xA, yA, xB, yB
    Variable nStep, coreHandling
    Variable xV_nm, yV_nm, rCore_nm
    Variable minUsedFrac
    Wave PhaseReCore, PhaseImCore

    nStep = max(1, round(nStep))
    coreHandling = round(coreHandling)

    Variable nQx = DimSize(Qx, 0)
    Variable nQy = DimSize(Qx, 1)
    Variable qxOffset = DimOffset(Qx, 0)
    Variable qxDelta  = DimDelta(Qx, 0)
    Variable qyOffset = DimOffset(Qx, 1)
    Variable qyDelta  = DimDelta(Qx, 1)

    Variable dxTot = xB - xA
    Variable dyTot = yB - yA
    Variable lenTot = sqrt(dxTot*dxTot + dyTot*dyTot)
    Variable tiny = 1e-12

    if (!(lenTot > tiny))
        return NaN
    endif

    Variable betaInt = 0
    Variable nUsed = 0
    Variable nBad = 0
    Variable k, ix, iy
    Variable t, xPos, yPos, xSample, ySample
    Variable qxVal, qyVal
    Variable rx, ry, r2, r
    Variable rCore2 = rCore_nm*rCore_nm
    Variable dlx = dxTot/nStep
    Variable dly = dyTot/nStep

    // ---------- mode 4: explicit split and phase-field fill ----------
    if (coreHandling == 4)

        if (ParamIsDefault(PhaseReCore) || ParamIsDefault(PhaseImCore))
            return NaN
        endif

        Variable fx = xA - xV_nm
        Variable fy = yA - yV_nm
        Variable aa = dxTot*dxTot + dyTot*dyTot
        Variable bb = 2*(fx*dxTot + fy*dyTot)
        Variable cc = fx*fx + fy*fy - rCore2
        Variable disc = bb*bb - 4*aa*cc

        if (disc > 0)
            Variable sqDisc = sqrt(disc)
            Variable t1 = (-bb - sqDisc)/(2*aa)
            Variable t2 = (-bb + sqDisc)/(2*aa)
            Variable tEntry = max(t1, 0)
            Variable tExit  = min(t2, 1)

            if ((tExit > tEntry) && (tExit >= 0) && (tEntry <= 1))

                Variable xEntry = xA + tEntry*dxTot
                Variable yEntry = yA + tEntry*dyTot
                Variable xExit  = xA + tExit*dxTot
                Variable yExit  = yA + tExit*dyTot

                ix = round((xEntry - qxOffset)/qxDelta)
                iy = round((yEntry - qyOffset)/qyDelta)
                if ((ix < 0) || (ix >= nQx) || (iy < 0) || (iy >= nQy))
                    return NaN
                endif
                Variable reEntry = PhaseReCore[ix][iy]
                Variable imEntry = PhaseImCore[ix][iy]

                ix = round((xExit - qxOffset)/qxDelta)
                iy = round((yExit - qyOffset)/qyDelta)
                if ((ix < 0) || (ix >= nQx) || (iy < 0) || (iy >= nQy))
                    return NaN
                endif
                Variable reExit = PhaseReCore[ix][iy]
                Variable imExit = PhaseImCore[ix][iy]

                if (numtype(reEntry) || numtype(imEntry) || numtype(reExit) || numtype(imExit))
                    return NaN
                endif

                Variable betaCore = atan2(imExit, reExit) - atan2(imEntry, reEntry)
                if (betaCore > pi)
                    betaCore -= 2*pi
                elseif (betaCore < -pi)
                    betaCore += 2*pi
                endif

                betaInt = betaCore

                Variable iSeg, xSegA, ySegA, xSegB, ySegB
                Variable segLen, nStepSeg, kSeg, tSeg
                Variable dlxSeg, dlySeg
                Variable nSegTotal = 0

                for (iSeg = 0; iSeg < 2; iSeg += 1)

                    if (iSeg == 0)
                        xSegA = xA
                        ySegA = yA
                        xSegB = xEntry
                        ySegB = yEntry
                    else
                        xSegA = xExit
                        ySegA = yExit
                        xSegB = xB
                        ySegB = yB
                    endif

                    segLen = sqrt((xSegB-xSegA)^2 + (ySegB-ySegA)^2)
                    if (!(segLen > tiny))
                        continue
                    endif

                    nStepSeg = max(1, ceil(nStep*segLen/lenTot))
                    nSegTotal += nStepSeg
                    dlxSeg = (xSegB-xSegA)/nStepSeg
                    dlySeg = (ySegB-ySegA)/nStepSeg

                    for (kSeg = 0; kSeg < nStepSeg; kSeg += 1)

                        tSeg = (kSeg + 0.5)/nStepSeg
                        xSample = xSegA + tSeg*(xSegB-xSegA)
                        ySample = ySegA + tSeg*(ySegB-ySegA)

                        ix = round((xSample - qxOffset)/qxDelta)
                        iy = round((ySample - qyOffset)/qyDelta)

                        if ((ix < 0) || (ix >= nQx) || (iy < 0) || (iy >= nQy))
                            nBad += 1
                            continue
                        endif

                        qxVal = Qx[ix][iy]
                        qyVal = Qy[ix][iy]

                        if (numtype(qxVal) || numtype(qyVal))
                            nBad += 1
                            continue
                        endif

                        betaInt += qxVal*dlxSeg + qyVal*dlySeg
                        nUsed += 1
                    endfor
                endfor

                if (nSegTotal < 1)
                    return NaN
                endif
                if ((nUsed/nSegTotal) < minUsedFrac)
                    return NaN
                endif

                return betaInt
            endif
        endif
    endif

    // ---------- modes 0-3, and mode 4 with no core crossing ----------
    for (k = 0; k < nStep; k += 1)

        t = (k + 0.5)/nStep
        xPos = xA + t*dxTot
        yPos = yA + t*dyTot
        xSample = xPos
        ySample = yPos

        if (coreHandling != 0)

            rx = xPos - xV_nm
            ry = yPos - yV_nm
            r2 = rx*rx + ry*ry

            if (r2 < rCore2)

                if (coreHandling == 2)
                    return NaN
                endif

                if (coreHandling == 3)
                    continue
                endif

                if (coreHandling == 1)
                    r = sqrt(r2)
                    if (r > tiny)
                        xSample = xV_nm + rCore_nm*rx/r
                        ySample = yV_nm + rCore_nm*ry/r
                    else
                        r = sqrt(dlx*dlx + dly*dly)
                        if (r > tiny)
                            xSample = xV_nm + rCore_nm*dlx/r
                            ySample = yV_nm + rCore_nm*dly/r
                        else
                            nBad += 1
                            continue
                        endif
                    endif
                endif
            endif
        endif

        ix = round((xSample - qxOffset)/qxDelta)
        iy = round((ySample - qyOffset)/qyDelta)

        if ((ix < 0) || (ix >= nQx) || (iy < 0) || (iy >= nQy))
            nBad += 1
            continue
        endif

        qxVal = Qx[ix][iy]
        qyVal = Qy[ix][iy]

        if (numtype(qxVal) || numtype(qyVal))
            nBad += 1
            continue
        endif

        betaInt += qxVal*dlx + qyVal*dly
        nUsed += 1
    endfor

    if ((nUsed/nStep) < minUsedFrac)
        return NaN
    endif

    return betaInt
End

//==============================================================================
// SNS_ExtractModes3DForFolder
//
// Purpose:
//   Extract 3D SNS trajectory channels from a 2D top-view mask using a
//   Cu(111)-consistent faceted prism geometry.
//
//   The STM position is located at the top surface:
//
//      z = H_nm
//
//   The supplied 2D mask is interpreted as the measured island outline at the
//   top surface. Toward the substrate/S plane:
//
//      z = 0
//
//   the island expands laterally according to close-packed fcc Cu(111)
//   stacking:
//
//      lateral / vertical = 1/sqrt(2)
//
//   equivalently:
//
//      side-wall angle = 35.264 deg from vertical
//                      = 54.736 deg from the Cu(111) surface plane.
//
//   The actual 3D channel construction is performed by
//   SNS_BuildChannelsFromMask3D(...), which should use the same Cu(111)
//   faceted-prism inside-test and side-wall reflection logic as the plotting
//   retracer.
//
//   The function stores the 3D channel lists generated by
//   SNS_BuildChannelsFromMask3D and identifies:
//
//      v_Gap0    : channel with largest W_eff_List
//      v_Longest : channel with largest L_N_List
//
//   If doDisplay=1, the function displays:
//
//      1. top-view image with projected selected 3D channels
//      2. trajectory-property plot of L_N_List, W_eff_List, T_eff_List
//
//   Optional:
//      If plotSideMasks=1, also displays two side-view silhouettes of the
//      Cu(111)-faceted prism and overlays the same selected channels:
//
//          SNS_3DMaskSide_XZ : x_parallel_B vs z
//          SNS_3DMaskSide_YZ : y_perp_B     vs z
//
//      Here x is always defined parallel to the in-plane magnetic field B,
//      and y is perpendicular to B in the sample plane.
//
// Inputs:
//   dfPath      : data folder containing w_mask and *_Z_mbgnd_xy image.
//   Bangle_deg  : in-plane B-field angle in degrees. Default: 225.
//   BTK_barrier : BTK barrier parameter. Default: 0.1.
//   STSx        : STM spectroscopy x position in nm. Default: -4.
//   STSy        : STM spectroscopy y position in nm. Default: -50.
//   Vortexx     : compatibility parameter for vortex maps. Default: STSx.
//   Vortexy     : compatibility parameter for vortex maps. Default: STSy.
//   H_nm        : island height in nm. Default: 15.
//   maxPath_nm  : path tracing limit passed to SNS_BuildChannelsFromMask3D.
//                 Default: 1000.
//   doDisplay   : 1 default, create diagnostic displays.
//                 0, suppress display output.
//   plotSideMasks:
//                 0 default, no side-view mask plots.
//                 1, create x-z and y-z side-view Cu(111)-prism plots and
//                    overlay the longest and largest-W_eff channels.
//
// Outputs in dfPath:
//   From SNS_BuildChannelsFromMask3D:
//      phiList
//      thetaList
//      L_N_List
//      W_eff_List
//      Hit1x_List
//      Hit1y_List
//      Hit2x_List
//      Hit2y_List
//      T_eff_List
//      wChan
//
//   Additional:
//      A_eff
//      v_Gap0
//      v_Longest
//      w_area_nm2
//      w_perim_nm
//
// Return:
//   0
//
// Notes:
//   This function itself is primarily a folder/display wrapper. The Cu(111)
//   geometry is enforced by the channel builder and retracer helpers.
//   Existing calls remain backward compatible because plotSideMasks is optional
//   and defaults to 0.
//==============================================================================
Function SNS_ExtractModes3DForFolder(dfPath, [Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy, H_nm, maxPath_nm, doDisplay, plotSideMasks])
    String   dfPath
    Variable Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy
    Variable H_nm, maxPath_nm, doDisplay, plotSideMasks

    String oldDF = GetDataFolder(1)
    SetDataFolder $dfPath

    Variable r=65535, g=65535, b=65535, a=65535

    // ---------------- Defaults ----------------
    if (ParamIsDefault(Bangle_deg))
        Bangle_deg = 225
    endif
    if (ParamIsDefault(BTK_barrier))
        BTK_barrier = 0.1
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
    if (ParamIsDefault(H_nm))
        H_nm = 15
    endif
    if (ParamIsDefault(maxPath_nm))
        maxPath_nm = 1000
    endif
    if (ParamIsDefault(doDisplay))
        doDisplay = 1
    endif
    if (ParamIsDefault(plotSideMasks))
        plotSideMasks = 0
    endif

    Variable v_Bangle = Bangle_deg * pi/180

    // ---------------- Required waves ----------------
    Wave/Z w_mask = w_mask

    if (!WaveExists(w_mask))
        SetDataFolder $oldDF
        Abort "SNS_ExtractModes3DForFolder: w_mask not found in folder."
    endif

    String imgList = WaveList("*_Z_mbgnd_xy",";","DIMS:2")
    if (strlen(imgList)==0)
        SetDataFolder $oldDF
        Abort "SNS_ExtractModes3DForFolder: No *_Z_mbgnd_xy image found in folder."
    endif

    String imgName = StringFromList(0, imgList, ";")
    Wave image = $imgName

    // Do not use broad KillVariables/Z/A here. It can remove unrelated
    // data-folder variables. The outputs below are overwritten explicitly.

    // ---------------- Build geometry / compatibility maps ----------------
    SNS_MakeVortexA(w_mask, "Ax_Vortex", "Ay_Vortex", Vortexx, Vortexy, 2, 50)
    SNS_MakeDeltaVortexMap(w_mask, "Delta_map", Vortexx, Vortexy, 100)

    // Builds 3D Cu(111)-faceted-prism channels if SNS_BuildChannelsFromMask3D
    // has been updated with SNS_IsInsideCu111Prism / SNS_Cu111PrismNormalComponent.
    // Remove old channel outputs before rebuilding.
	 KillWaves/Z phiList_rad, thetaList_rad
	 KillWaves/Z L_N_List_nm, W_eff_List_nm, W_geom_List_nm, Aeff_3D_List_nm2
	 KillWaves/Z Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm
	 KillWaves/Z T_eff_List, wChan
	 KillVariables/Z v_Gap0, v_Longest
    SNS_BuildChannelsFromMask3D(w_mask, STSx, STSy, v_Bangle, 0.5, H_nm, maxPath_nm, "")

    Wave/Z phiList_rad
    Wave/Z thetaList_rad
    Wave/Z L_N_List_nm
    Wave/Z W_eff_List_nm
    Wave/Z Hit1x_List_nm
    Wave/Z Hit1y_List_nm
    Wave/Z Hit2x_List_nm
    Wave/Z Hit2y_List_nm
    Wave/Z T_eff_List
    wave/Z W_geom_List_nm

    if (!WaveExists(phiList_rad) || !WaveExists(thetaList_rad) || !WaveExists(L_N_List_nm) || !WaveExists(W_eff_List_nm))
        SetDataFolder $oldDF
        Abort "SNS_ExtractModes3DForFolder: channel builder did not create required channel waves."
    endif

    if (numpnts(L_N_List_nm) <= 0)
        SetDataFolder $oldDF
        Abort "SNS_ExtractModes3DForFolder: no valid 3D channels generated."
    endif

    // ---------------- Area/perimeter ----------------
    SNS_MaskAreaPerim_FromParticles(w_mask)
    Wave/Z w

    if (WaveExists(w) && numpnts(w) >= 2)
        Make/O/N=1 w_area_nm2 = w[0]
        Make/O/N=1 w_perim_nm = w[1]
    else
        Make/O/N=1 w_area_nm2 = NaN
        Make/O/N=1 w_perim_nm = NaN
    endif

	 // ---------------- Identify special modes ----------------
	 // W_eff_List_nm is the solver-equivalent magnetic width.
	 // W_geom_List_nm is the geometric endpoint projection.
	 // For intuitive ray plotting, choose v_Gap0 from W_geom_List_nm.
	
	 Wave/Z W_geom_List_nm
	 if (WaveExists(W_geom_List_nm))
	     WaveStats/Q W_geom_List_nm
	     Variable/G v_Gap0 = V_maxloc
	 else
	     // Fallback for old channel sets.
	     WaveStats/Q W_eff_List_nm
	     Variable/G v_Gap0 = V_maxloc
	 endif
	
	WaveStats/Q L_N_List_nm
	Variable/G v_Longest = V_maxloc

    if (doDisplay)

        // ---------------- Topography ----------------
        SNS_DisplayWithScales(image, cmap="grayC")
        String winImage = WinName(0, 1, 1)

        SNS_StyleRayDisplayImage(image, winImage)

        // ---- Plot selected full 3D trajectories projected into top view ----
        // These wrappers should now use the Cu(111)-faceted retracer internally.
        Variable errPlot

        errPlot = SNS_PlotChannel3D_Surface(w_mask, STSx, STSy,             \
                                            thetaList_rad, phiList_rad, v_Gap0,     \
                                            0.5, H_nm, maxPath_nm,          \
                                            winImage, "ray3D_gap")

        if (errPlot == 0)
            String rayGapY = "ray3D_gap_Y_" + num2str(v_Gap0)
            ModifyGraph rgb($rayGapY)=(0,0,0)
        endif

        errPlot = SNS_PlotChannel3D_Surface(w_mask, STSx, STSy,             \
                                            thetaList_rad, phiList_rad, v_Longest,  \
                                            0.5, H_nm, maxPath_nm,          \
                                            winImage, "ray3D_long")

        if (errPlot == 0)
            String rayLongY = "ray3D_long_Y_" + num2str(v_Longest)
            ModifyGraph rgb($rayLongY)=(65535,0,0)
        endif

        // ---------------- Textbox ----------------
        String txtLongest, txtGap, txtBox

			if (WaveExists(W_geom_List_nm))
			    txtLongest = "Longest Mode\r" + \
			        "L = " + num2str(round(L_N_List_nm[v_Longest])) + " nm\r" + \
			        "W_geom = " + num2str(round(W_geom_List_nm[v_Longest])) + " nm\r" + \
			        "W_eff = " + num2str(round(W_eff_List_nm[v_Longest])) + " nm"
			
			    txtGap = "Largest geometric W\r" + \
			        "L = " + num2str(round(L_N_List_nm[v_Gap0])) + " nm\r" + \
			        "W_geom = " + num2str(round(W_geom_List_nm[v_Gap0])) + " nm\r" + \
			        "W_eff = " + num2str(round(W_eff_List_nm[v_Gap0])) + " nm"
			else
			    txtLongest = "Longest Mode\r" + \
			        "L = " + num2str(round(L_N_List_nm[v_Longest])) + " nm\r" + \
			        "W_eff = " + num2str(round(W_eff_List_nm[v_Longest])) + " nm"
			
			    txtGap = "Largest W_eff\r" + \
			        "L = " + num2str(round(L_N_List_nm[v_Gap0])) + " nm\r" + \
			        "W_eff = " + num2str(round(W_eff_List_nm[v_Gap0])) + " nm"
			endif

        txtBox = txtLongest + "\r" + txtGap

        TextBox/C/N=TextMode/X=0/Y=0/F=0/B=(r,g,b,a*0.6)/A=RT txtBox

        Make/O/N=1 tmpSTSx=STSx, tmpSTSy=STSy
        AppendToGraph tmpSTSy vs tmpSTSx
        ModifyGraph mode(tmpSTSy)=3, rgb(tmpSTSy)=(65535,65535,65535)

        SNS_TAGSetpoint()

        // ---------------- B-field arrow ----------------
        Variable xMin = DimOffset(image, 0)
        Variable xMax = xMin + DimDelta(image, 0) * (DimSize(image, 0) - 1)

        Variable yMin = DimOffset(image, 1)
        Variable yMax = yMin + DimDelta(image, 1) * (DimSize(image, 1) - 1)

        Variable pos_x = xMax - 0.15 * (xMax - xMin)
        Variable pos_y = yMin + 0.15 * (yMax - yMin)

        SNS_DrawImageArrow(image, Bangle_deg, "B", pos_x=pos_x, pos_y=pos_y, scale=1.5, color="White")

        // ---------------- Clamp top-view window size ----------------
        GetWindow $WinName(0, 1, 1), wsize
        Variable wPx = abs(V_right - V_left)
        Variable hPx = abs(V_bottom - V_top)

        Variable maxDim = max(wPx, hPx)
        if (maxDim > 512)
            Variable ex = 512 / maxDim
            ModifyGraph/W=$WinName(0, 1, 1) expand=ex
        else
            ModifyGraph/W=$WinName(0, 1, 1) expand=1
        endif

        // ---------------- Optional side-view Cu(111) prism plots ----------------
			if (plotSideMasks)
			    SNS_Show3DMaskTopFieldFrame_WithChannels(w_mask, STSx, STSy, thetaList_rad, phiList_rad, \
			        v_Longest, v_Gap0, Bangle_deg, H_nm, maxPath_nm)
			
			    SNS_Show3DMaskSideViews_WithChannels(w_mask, STSx, STSy, thetaList_rad, phiList_rad, \
			        v_Longest, v_Gap0, Bangle_deg, H_nm, maxPath_nm)
			endif

        // ---------------- Mode property plot ----------------
        Display L_N_List_nm, W_eff_List_nm
        String winModes = WinName(0, 1, 1)

        if (WaveExists(T_eff_List))
            AppendToGraph/R T_eff_List
        endif

        ModifyGraph rgb(L_N_List_nm)=(0,0,0)
        ModifyGraph rgb(W_eff_List_nm)=(1,16019,65535)
        ModifyGraph mode(L_N_List_nm)=2
        ModifyGraph mode(W_eff_List_nm)=2
        ModifyGraph lsize=5

        if (WaveExists(T_eff_List))
            ModifyGraph mode(T_eff_List)=3
            ModifyGraph marker(T_eff_List)=8
            SetAxis right 0,1
            ModifyGraph muloffset(T_eff_List)={0,0}
        endif

        SetAxis left 0,*

        ModifyGraph muloffset(L_N_List_nm)={0,0}
        ModifyGraph muloffset(W_eff_List_nm)={0,0}
        ModifyGraph tick=2
        ModifyGraph mirror(bottom)=2
        ModifyGraph ZisZ=1
        ModifyGraph standoff(left)=0
        ModifyGraph standoff(bottom)=0

        Label left   "Length (nm)"
        Label bottom "Trajectory Nr."

        if (WaveExists(T_eff_List))
            Label right  "Transparency"
            Legend/C/N=text0/J/F=0/A=RT "Mode\r" + \
                "\\s(L_N_List_nm) length (L)\r" + \
                "\\s(W_eff_List_nm) projection ⊥B\\Bext\\M (W)\r" + \
                "\\s(T_eff_List) transparency (T)"
        else
            Legend/C/N=text0/J/F=0/A=RT "Mode\r" + \
                "\\s(L_N_List_nm) length (L)\r" + \
                "\\s(W_eff_List_nm) projection ⊥B\\Bext\\M (W)"
        endif

        ModifyGraph margin(left)=33
        ModifyGraph margin(bottom)=30
        ModifyGraph margin(right)=35
        ModifyGraph margin(top)=6
        ModifyGraph width=220
        ModifyGraph height=220

    endif

    // ---------------- Restore folder ----------------
    SetDataFolder $oldDF

    return 0
End


//==============================================================================
// SNS_EstimateNdir_FromMask3D
//
// 3D analogue of SNS_EstimateNphi_FromMask.
//
// - Nmask    : 2D mask (same as 2D case)
// - lambdaF  : Fermi wavelength [same units as mask axes, e.g. nm]
// - H        : N-region height (STM→S distance) in same units as mask axes
//
// Returns an estimate of the *total* number of angular directions needed
// to resolve the geometry in 3D.
//
// 2D: Lchar = sqrt(areaN)
// 3D: Lchar3D = (areaN * H)^(1/3)
// N_dir ≈ π Lchar3D / λF
//==============================================================================

Function SNS_EstimateNdir_FromMask3D(Nmask, lambdaF, H)
    Wave    Nmask
    Variable lambdaF
    Variable H          // height in same units as DimDelta (e.g. nm)

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if (nx <= 1 || ny <= 1 || lambdaF <= 0 || H <= 0)
        return 0
    endif

    Variable dx = DimDelta(Nmask, 0)
    Variable dy = DimDelta(Nmask, 1)
    Variable cellArea = abs(dx * dy)

    // --- compute N-region area ---
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

    // 3D characteristic linear size from volume
    Variable volumeN = areaN * H             // [length^3]
    Variable Lchar3D = volumeN^(1.0/3.0)     // [length]

    // total number of directions ≈ π Lchar3D / λF
    Variable Ndir = round(pi * Lchar3D / lambdaF)

    if (Ndir < 1)
        Ndir = 1
    endif

    return Ndir
End

//==============================================================================
// SNS_ProjectChannel3D_Surface
//
// Purpose:
//   Backward-compatible top-view projection wrapper.
//   Uses SNS_ProjectChannel3D_SurfaceXYZ, which itself calls the shared tracer.
//==============================================================================
Function SNS_ProjectChannel3D_Surface(Nmask, r0x, r0y, theta, phi, stepFac, H_nm, maxPath_nm, rayXName, rayYName)
    Wave    Nmask
    Variable r0x, r0y
    Variable theta, phi
    Variable stepFac, H_nm, maxPath_nm
    String  rayXName, rayYName

    String zName = rayXName + "_Ztmp"

    Variable err = SNS_ProjectChannel3D_SurfaceXYZ(Nmask, r0x, r0y, theta, phi, stepFac, H_nm, maxPath_nm, rayXName, rayYName, zName)

    KillWaves/Z $zName

    return err
End




//==============================================================================
// SNS_PlotChannel3D_Surface
//
// Plot full 3D channel footprint (both legs) on a 2D STM/topography image.
//
// Inputs:
//   Nmask         : same mask used by the 3D builder.
//   r0x, r0y      : STM position [nm].
//   thetaList     : polar angles per channel (output from SNS_BuildChannelsFromMask3D).
//   phiList       : azimuth angles per channel (output from SNS_BuildChannelsFromMask3D).
//   idx           : channel index to plot.
//   stepFac       : same as in builder.
//   H_nm          : N height [nm].
//   maxPath_nm    : max path per leg [nm].
//   graphName     : target graph window name (typically STM image window).
//   rayBaseName   : base name for generated ray waves.
//
// Behavior:
//   Appends ray trace to graph, with solid line.
//==============================================================================
Function SNS_PlotChannel3D_Surface(Nmask, r0x, r0y, thetaList, phiList, idx, stepFac, H_nm, maxPath_nm, graphName, rayBaseName)
    Wave    Nmask
    Variable r0x, r0y
    Wave    thetaList, phiList
    Variable idx
    Variable stepFac, H_nm, maxPath_nm
    String  graphName, rayBaseName

    if (idx < 0 || idx >= DimSize(thetaList, 0))
        return -1
    endif

    Variable theta = thetaList[idx]
    Variable phi   = phiList[idx]

    String rayXname = rayBaseName + "_X_" + num2str(idx)
    String rayYname = rayBaseName + "_Y_" + num2str(idx)

    Variable err = SNS_ProjectChannel3D_Surface(Nmask, r0x, r0y, theta, phi, \
                                                stepFac, H_nm, maxPath_nm,   \
                                                rayXname, rayYname)
    if (err != 0)
        return err
    endif

    Wave rayY = $rayYname
    Wave rayX = $rayXname

    AppendToGraph/W=$graphName rayY vs rayX
    ModifyGraph mode($NameOfWave(rayY))=0, lsize($NameOfWave(rayY))=2

    return 0
End





//==============================================================================
// SNS_BuildDOSChannelsFromMask2D
//
// Purpose
//   Build a GLOBAL 2D channel ensemble for DOS / free-energy calculations by
//   scanning anchor points over the full N mask, calling the existing local
//   SNS_BuildChannelsFromMask2D(...) at each anchor point, and collapsing
//   duplicate trajectories based on their endpoint pair.
//
// Strategy
//   1. Loop over a regular anchor grid inside Nmask.
//   2. At each valid anchor point, call SNS_BuildChannelsFromMask2D(...).
//   3. Collect all returned rays.
//   4. Canonicalize endpoint order so forward/backward copies are identical.
//   5. Quantize endpoints (and optionally length) using tolerances.
//   6. Sort by the quantized key and merge duplicates.
//   7. Store one representative per unique chord, with wChan = multiplicity.
//
// Inputs
//   Nmask          : 2D mask in nm coordinates; N if >0.5
//   phiB           : field angle [rad]
//   stepFac        : same marching step factor used by SNS_BuildChannelsFromMask2D
//   Zbarrier       : BTK barrier parameter
//   folder         : output data folder; "" means current folder
//
// Optional inputs
//   anchorStepX_px : anchor-point step in x, in mask pixels; default 1
//   anchorStepY_px : anchor-point step in y, in mask pixels; default 1
//   endpointTol_nm : endpoint grouping tolerance [nm]; default = min(|dx|,|dy|)
//   lengthTol_nm   : length grouping tolerance [nm];   default = endpointTol_nm
//
// Outputs (same naming convention as SNS_BuildChannelsFromMask2D)
//   phiList        : representative ray angle [rad]
//   L_N_List       : chord length [m]
//   W_eff_List     : magnetic width [m]
//   wChan          : multiplicity of merged copies
//   T_eff_List     : representative transparency
//   Hit1x_List     : endpoint 1 x [nm]
//   Hit1y_List     : endpoint 1 y [nm]
//   Hit2x_List     : endpoint 2 x [nm]
//   Hit2y_List     : endpoint 2 y [nm]
//
// Additional diagnostics
//   nCopies_List   : same as wChan but kept explicitly as integer-like wave
//   nAnchorsUsed   : number of anchor points inside Nmask that were sampled
//
// Notes
//   - Endpoint order is canonicalized lexicographically.
//   - Merged representatives are simple averages over the duplicate group.
//   - wChan carries the duplicate multiplicity and should be used as the
//     geometric weight in the downstream DOS builder.
//==============================================================================

Function SNS_BuildDOSChannelsFromMask2D(Nmask, phiB, stepFac, Zbarrier, folder, [anchorStepX_px, anchorStepY_px, endpointTol_nm, lengthTol_nm])

    Wave    Nmask
    Variable phiB
    Variable stepFac
    Variable Zbarrier
    String  folder
    Variable anchorStepX_px, anchorStepY_px, endpointTol_nm, lengthTol_nm

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if ((nx <= 1) || (ny <= 1))
        return -1
    endif

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    if (ParamIsDefault(anchorStepX_px) || (anchorStepX_px < 1))
        anchorStepX_px = 1
    endif
    if (ParamIsDefault(anchorStepY_px) || (anchorStepY_px < 1))
        anchorStepY_px = 1
    endif

    if (ParamIsDefault(endpointTol_nm) || (endpointTol_nm <= 0))
        endpointTol_nm = min(abs(dx), abs(dy))
    endif
    if (ParamIsDefault(lengthTol_nm) || (lengthTol_nm <= 0))
        lengthTol_nm = endpointTol_nm
    endif

    String oldDF = GetDataFolder(1)
    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    // Temporary folder reused for each local builder call
    String tmpFolder = "root:Packages:SNS:tmpBuildDOS2D"
    NewDataFolder/O root:Packages
    NewDataFolder/O root:Packages:SNS
    NewDataFolder/O/S $tmpFolder
    SetDataFolder $oldDF
    if (strlen(folder) > 0)
        SetDataFolder $folder
    endif

    // Start empty global raw lists
    Make/O/D/N=0 phiList_raw, L_N_List_raw, W_eff_List_raw, T_eff_List_raw
    Make/O/D/N=0 Hit1x_List_raw, Hit1y_List_raw, Hit2x_List_raw, Hit2y_List_raw
    Make/O/T/N=0 key_raw

    Variable ix0, iy0, xA, yA
    Variable nAnchorsUsedLocal  = 0
    Variable nTot = 0
    Variable nLoc, iLoc
    Variable x1, y1, x2, y2, xt, yt
    Variable Lnm
    Variable qx1, qy1, qx2, qy2, qL
    String keyStr

    for (ix0 = 0; ix0 < nx; ix0 += anchorStepX_px)
        for (iy0 = 0; iy0 < ny; iy0 += anchorStepY_px)

            if (Nmask[ix0][iy0] <= 0.5)
                continue
            endif

            xA = x0 + dx*ix0
            yA = y0 + dy*iy0

            // Build local channel fan into the reusable temporary folder
            Variable err = SNS_BuildChannelsFromMask2D(Nmask, xA, yA, phiB, stepFac, tmpFolder)
            if (err != 0)
                continue
            endif

            nAnchorsUsedLocal  += 1

            Wave Lloc  = $(tmpFolder + ":L_N_List_nm")
            Wave Wloc  = $(tmpFolder + ":W_eff_List_nm")
            Wave Ploc  = $(tmpFolder + ":phiList_rad")
            Wave Tloc  = $(tmpFolder + ":T_eff_List")
            Wave X1loc = $(tmpFolder + ":Hit1x_List_nm")
            Wave Y1loc = $(tmpFolder + ":Hit1y_List_nm")
            Wave X2loc = $(tmpFolder + ":Hit2x_List_nm")
            Wave Y2loc = $(tmpFolder + ":Hit2y_List_nm")

            nLoc = numpnts(Lloc)
            if (nLoc <= 0)
                continue
            endif

            Redimension/N=(nTot + nLoc) phiList_raw, L_N_List_raw, W_eff_List_raw, T_eff_List_raw
            Redimension/N=(nTot + nLoc) Hit1x_List_raw, Hit1y_List_raw, Hit2x_List_raw, Hit2y_List_raw
            Redimension/N=(nTot + nLoc) key_raw

            for (iLoc = 0; iLoc < nLoc; iLoc += 1)

                Lnm = Lloc[iLoc]

                x1 = X1loc[iLoc]
                y1 = Y1loc[iLoc]
                x2 = X2loc[iLoc]
                y2 = Y2loc[iLoc]

                // Canonical endpoint order: lexicographic by (x,y)
                if ((x2 < x1) || ((x2 == x1) && (y2 < y1)))
                    xt = x1; yt = y1
                    x1 = x2; y1 = y2
                    x2 = xt; y2 = yt
                endif

                qx1 = round(x1 / endpointTol_nm)
                qy1 = round(y1 / endpointTol_nm)
                qx2 = round(x2 / endpointTol_nm)
                qy2 = round(y2 / endpointTol_nm)
                qL  = round(Lnm / lengthTol_nm)

                sprintf keyStr, "%012.0f_%012.0f_%012.0f_%012.0f_%012.0f", qx1, qy1, qx2, qy2, qL

                phiList_raw[nTot + iLoc]  = Ploc[iLoc]
                L_N_List_raw[nTot + iLoc] = Lloc[iLoc]
                W_eff_List_raw[nTot + iLoc] = Wloc[iLoc]
                T_eff_List_raw[nTot + iLoc] = Tloc[iLoc]

                Hit1x_List_raw[nTot + iLoc] = x1
                Hit1y_List_raw[nTot + iLoc] = y1
                Hit2x_List_raw[nTot + iLoc] = x2
                Hit2y_List_raw[nTot + iLoc] = y2

                key_raw[nTot + iLoc] = keyStr
            endfor

            nTot += nLoc
        endfor
    endfor

    if (nTot <= 0)
        SetDataFolder $oldDF
        return -2
    endif

    // Sort by endpoint/length key so duplicates become adjacent
    Sort key_raw, key_raw, phiList_raw, L_N_List_raw, W_eff_List_raw, T_eff_List_raw, Hit1x_List_raw, Hit1y_List_raw, Hit2x_List_raw, Hit2y_List_raw

    // First pass: count unique groups
    Variable nUnique = 0
    Variable i = 0
    do
        nUnique += 1
        String keyHere = key_raw[i]
        do
            i += 1
        while ((i < nTot) && (cmpstr(key_raw[i], keyHere) == 0))
    while (i < nTot)

    // Allocate final merged waves
    Make/O/D/N=(nUnique) phiList_rad, L_N_List_nm, W_eff_List_nm, wChan, T_eff_List
    Make/O/D/N=(nUnique) Hit1x_List_nm, Hit1y_List_nm, Hit2x_List_nm, Hit2y_List_nm, nCopies_List

    // Second pass: merge duplicates by averaging representative quantities
    Variable u = 0
    Variable j0, nGrp
    Variable sumPhi, sumL, sumW, sumT, sumX1, sumY1, sumX2, sumY2
    i = 0
    do
        String keyNow = key_raw[i]
        j0 = i

        sumPhi = 0
        sumL   = 0
        sumW   = 0
        sumT   = 0
        sumX1  = 0
        sumY1  = 0
        sumX2  = 0
        sumY2  = 0
        nGrp   = 0

        do
            sumPhi += phiList_raw[i]
            sumL   += L_N_List_raw[i]
            sumW   += W_eff_List_raw[i]
            sumT   += T_eff_List_raw[i]
            sumX1  += Hit1x_List_raw[i]
            sumY1  += Hit1y_List_raw[i]
            sumX2  += Hit2x_List_raw[i]
            sumY2  += Hit2y_List_raw[i]

            nGrp += 1
            i += 1
        while ((i < nTot) && (cmpstr(key_raw[i], keyNow) == 0))

        phiList_rad[u]   = sumPhi / nGrp
        L_N_List_nm[u]   = sumL   / nGrp
        W_eff_List_nm[u] = sumW   / nGrp
        T_eff_List[u]    = sumT   / nGrp

        Hit1x_List_nm[u] = sumX1  / nGrp
        Hit1y_List_nm[u] = sumY1  / nGrp
        Hit2x_List_nm[u] = sumX2  / nGrp
        Hit2y_List_nm[u] = sumY2  / nGrp

        wChan[u]       = nGrp
        nCopies_List[u]= nGrp

        u += 1
    while (i < nTot)

    Variable/G nAnchorsUsed = nAnchorsUsedLocal 

    SetDataFolder $oldDF
    return 0
End

//==============================================================================
// SNS_BuildTrialPointListFromMask2D
//
// Purpose
//   Build a sparse list of candidate vortex positions from a 2D N-region mask.
//
//   Unlike a rectangular xTry_nm × yTry_nm grid, this function returns paired
//   trial coordinates (xTry_nm[k], yTry_nm[k]) such that every candidate point
//   lies inside the mask.
//
// Inputs
//   Nmask      : 2D wave, N-region mask in sample coordinates.
//                Values > 0.5 are treated as valid N-region pixels.
//                X and Y scaling must be in nm.
//
//   folder     : Output data folder for the generated trial-point waves.
//
// Optional Inputs
//   stepX_px   : Sampling stride in x, in mask pixels. Default = 5.
//   stepY_px   : Sampling stride in y, in mask pixels. Default = stepX_px.
//
// Outputs (in folder)
//   xTry_nm    : 1D wave of trial x positions [nm]
//   yTry_nm    : 1D wave of trial y positions [nm]
//   nTryPoints : Number of valid trial points
//
// Notes
//   • xTry_nm and yTry_nm are paired point lists of equal length.
//   • The kth trial point is (xTry_nm[k], yTry_nm[k]).
//   • Only mask points with Nmask > 0.5 are included.
//==============================================================================
Function SNS_BuildTrialPointListFromMask2D(Nmask, folder, [stepX_px, stepY_px])

    Wave Nmask
    String folder
    Variable stepX_px, stepY_px

    if (ParamIsDefault(stepX_px))
        stepX_px = 5
    endif

    if (ParamIsDefault(stepY_px))
        stepY_px = stepX_px
    endif

    if ((stepX_px < 1) || (stepY_px < 1))
        Abort "SNS_BuildTrialPointListFromMask2D: step sizes must be >= 1."
    endif

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)

    if ((nx <= 1) || (ny <= 1))
        Abort "SNS_BuildTrialPointListFromMask2D: invalid mask."
    endif

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    String oldDF = GetDataFolder(1)
    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    Make/O/D/N=0 xTry_nm, yTry_nm

    Variable ix, iy
    Variable xnm, ynm
    Variable nTry = 0

    for (ix = 0; ix < nx; ix += stepX_px)
        for (iy = 0; iy < ny; iy += stepY_px)

            if (Nmask[ix][iy] <= 0.5)
                continue
            endif

            xnm = x0 + ix*dx
            ynm = y0 + iy*dy

            InsertPoints nTry, 1, xTry_nm, yTry_nm
            xTry_nm[nTry] = xnm
            yTry_nm[nTry] = ynm
            nTry += 1

        endfor
    endfor

    Variable/G nTryPoints = nTry

    SetDataFolder $oldDF
    return 0
End


//==============================================================================
// SNS_ComputeExteriorPhaseEnergy2D
//
// Purpose
//   Compute the phase-field energy of a solved exterior phase map using the
//   complex phase representation phaseRe_ext = cos(phi_ext) and
//   phaseIm_ext = sin(phi_ext).
//
//   The implemented discrete energy is
//
//       Fphase_raw = 1/2 * sum_ext dA *
//                    [ (dphi_x/dx)^2 + (dphi_y/dy)^2 ]
//
//   where dphi_x and dphi_y are nearest-neighbor wrapped phase differences
//   obtained directly from the complex phase field.
//
// Inputs
//   extMask        : 2D wave, exterior domain mask
//                    extMask > 0.5  -> exterior (valid phase)
//                    extMask <= 0.5 -> hole / invalid region
//
//   phaseRe_ext    : 2D wave, cos(phi_ext)
//   phaseIm_ext    : 2D wave, sin(phi_ext)
//
// Returns
//   Fphase_raw     : raw phase energy [rad^2]
//
// Notes
//   • Uses only phaseRe_ext / phaseIm_ext, not phi_ext.
//   • Only exterior-exterior bonds contribute.
//   • This avoids branch-cut artifacts from raw angle fields.
//==============================================================================
Function SNS_ComputeExteriorPhaseEnergy2D(extMask, phaseRe_ext, phaseIm_ext)

    Wave extMask, phaseRe_ext, phaseIm_ext

    Variable nx = DimSize(extMask, 0)
    Variable ny = DimSize(extMask, 1)

    if ((nx <= 1) || (ny <= 1))
        Abort "SNS_ComputeExteriorPhaseEnergy2D: invalid extMask dimensions."
    endif

    Variable dx = abs(DimDelta(extMask, 0))
    Variable dy = abs(DimDelta(extMask, 1))

    if ((dx <= 0) || (dy <= 0))
        Abort "SNS_ComputeExteriorPhaseEnergy2D: invalid wave scaling."
    endif

    Variable dA = dx * dy
    Variable ix, iy
    Variable re1, im1, re2, im2
    Variable dphi, Fphase_raw = 0

    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)

            if (extMask[ix][iy] <= 0.5)
                continue
            endif

            re1 = phaseRe_ext[ix][iy]
            im1 = phaseIm_ext[ix][iy]

            // x-bond
            if ((ix < nx-1) && (extMask[ix+1][iy] > 0.5))

                re2 = phaseRe_ext[ix+1][iy]
                im2 = phaseIm_ext[ix+1][iy]

                // wrapped phase difference from z2 * conj(z1)
                dphi = atan2(im2*re1 - re2*im1, re2*re1 + im2*im1)
                Fphase_raw += 0.5 * dA * (dphi/dx)^2

            endif

            // y-bond
            if ((iy < ny-1) && (extMask[ix][iy+1] > 0.5))

                re2 = phaseRe_ext[ix][iy+1]
                im2 = phaseIm_ext[ix][iy+1]

                dphi = atan2(im2*re1 - re2*im1, re2*re1 + im2*im1)
                Fphase_raw += 0.5 * dA * (dphi/dy)^2

            endif

        endfor
    endfor

    return Fphase_raw
End


//==============================================================================
// SNS_FindPointVortexMinEnergy2D_PointList
//
// Purpose
//   Find the minimum-energy configuration of a point-like vortex by scanning
//   a sparse list of candidate vortex positions.
//
//   The function always computes the exterior phase-field contribution and,
//   optionally, the ABS contribution.
//
//   For each trial point (xTry_nm[k], yTry_nm[k]), the function:
//
//     1. Solves the exterior phase field treating Nmask as a hole via
//          SNS_SolveExteriorPhaseField2D_FromBoundaryMaps(...)
//     2. Computes the per-channel extra phase betaExtra_List from the solved
//        boundary phase via
//          SNS_ComputeBetaExtraFromExteriorPhase2D(...)
//     3. Evaluates
//
//          Cost_try = CostABS_try + CostPhase_try
//
//        with
//
//          CostABS_try   = -sum_{ch,br} wChan[ch] * |E_{ch,br}|
//                          only if useABS != 0
//          CostPhase_try = alphaPhase_eV * Fphase_raw
//
// Inputs
//   Bval_T         : Magnetic field [T] for the evaluation.
//   Nmask          : 2D N-region mask used as hole geometry for the exterior
//                    phase solve. Values > 0.5 are treated as hole / excluded
//                    region. X and Y scaling must be in nm.
//   L_N_List       : Channel N lengths [m].
//   W_eff_List     : Channel effective widths [m].
//   wChan          : Channel weights / multiplicities.
//   T_eff_List     : Per-channel effective transparency.
//   Hit1x_List     : Endpoint 1 x coordinate [nm].
//   Hit1y_List     : Endpoint 1 y coordinate [nm].
//   Hit2x_List     : Endpoint 2 x coordinate [nm].
//   Hit2y_List     : Endpoint 2 y coordinate [nm].
//   xTry_nm_in     : 1D wave of trial x positions [nm].
//   yTry_nm_in     : 1D wave of trial y positions [nm].
//   folder         : Output data folder.
//
// Optional Inputs
//   nFlux          : Vortex winding number.
//                    Default = root:SNS_Settings:SNS_nFlux
//   alphaPhase_eV  : Energy prefactor multiplying the raw exterior phase
//                    energy. Default = 0.
//   phaseNIterMax  : Maximum number of iterations for exterior phase solve.
//                    Default = 100.
//   phaseTol_rad   : Convergence tolerance for exterior phase solve [rad].
//                    Default = 1e-6.
//   phaseOmega     : SOR relaxation parameter for exterior phase solve.
//                    Default = 1.6.
//   phaseSearchRad_px : Search radius used when sampling phase near the
//                    boundary. Default = 3.
//   useABS         : If nonzero, include ABS branch solving and ABS energy.
//                    Default = 0.
//   nSmoothIter    : Number of smoothing iterations used to precompute the
//                    boundary maps. Default = 2.
//
// Outputs (in folder)
//   Cost_try          : Total cost = CostABS_try + CostPhase_try [eV]
//   CostABS_try       : ABS contribution [eV]
//   CostPhase_try     : Exterior phase contribution [eV]
//   xTry_nm           : Copy of input x trial points [nm]
//   yTry_nm           : Copy of input y trial points [nm]
//   xV_best_nm        : Best-fit vortex x position [nm]
//   yV_best_nm        : Best-fit vortex y position [nm]
//   cost_best         : Minimum total cost [eV]
//   idx_best          : Index of the minimum in Cost_try
//   betaExtra_best    : Per-channel betaExtra [rad] for the best-fit trial point
//
// Notes
//   • Boundary geometry maps are computed once and reused for all trial points.
//   • Lower cost is better.
//   • By default, only the phase-field energy is evaluated.
//==============================================================================
Function SNS_FindPointVortexMinEnergy2D_PointList(Bval_T, Nmask, L_N_List, W_eff_List, wChan, T_eff_List, Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, xTry_nm_in, yTry_nm_in, folder, [nFlux, alphaPhase_eV, phaseNIterMax, phaseTol_rad, phaseOmega, phaseSearchRad_px, useABS, nSmoothIter])

    Variable Bval_T
    Wave Nmask
    Wave L_N_List, W_eff_List, wChan, T_eff_List
    Wave Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
    Wave xTry_nm_in, yTry_nm_in
    String folder
    Variable nFlux, alphaPhase_eV, phaseNIterMax, phaseTol_rad, phaseOmega, phaseSearchRad_px, useABS, nSmoothIter

    if (ParamIsDefault(nFlux))
        NVAR SNS_nFlux = root:SNS_Settings:SNS_nFlux
        nFlux = SNS_nFlux
    endif
    if (ParamIsDefault(alphaPhase_eV))
        alphaPhase_eV = 0
    endif
    if (ParamIsDefault(phaseNIterMax))
        phaseNIterMax = 100
    endif
    if (ParamIsDefault(phaseTol_rad))
        phaseTol_rad = 1e-6
    endif
    if (ParamIsDefault(phaseOmega))
        phaseOmega = 1.6
    endif
    if (ParamIsDefault(phaseSearchRad_px))
        phaseSearchRad_px = 3
    endif
    if (ParamIsDefault(useABS))
        useABS = 0
    endif
    if (ParamIsDefault(nSmoothIter))
        nSmoothIter = 2
    endif

    Variable Nch = numpnts(L_N_List)
    if ((Nch <= 0) || (Nch != numpnts(W_eff_List)) || (Nch != numpnts(wChan)) || \
        (Nch != numpnts(T_eff_List)) || (Nch != numpnts(Hit1x_List)) || \
        (Nch != numpnts(Hit1y_List)) || (Nch != numpnts(Hit2x_List)) || \
        (Nch != numpnts(Hit2y_List)))
        Abort "SNS_FindPointVortexMinEnergy2D_PointList: inconsistent channel waves."
    endif

    Variable Ntry = numpnts(xTry_nm_in)
    if ((Ntry <= 0) || (Ntry != numpnts(yTry_nm_in)))
        Abort "SNS_FindPointVortexMinEnergy2D_PointList: inconsistent trial-point waves."
    endif

    Variable sumW = sum(wChan)
    if (sumW <= 0)
        Abort "SNS_FindPointVortexMinEnergy2D_PointList: sum(wChan) <= 0."
    endif

    STRUCT SNS_Params SNS_p
    SNS_LoadParams(SNS_p)

    Variable Delta   = SNS_p.Delta
    Variable vF      = SNS_p.vF
    Variable lambdaL = SNS_p.lambdaL

    String oldDF = GetDataFolder(1)
    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    Duplicate/O xTry_nm_in, xTry_nm
    Duplicate/O yTry_nm_in, yTry_nm
    Make/O/D/N=(Ntry) Cost_try, CostABS_try, CostPhase_try
    Cost_try = NaN
    CostABS_try = 0
    CostPhase_try = NaN

    Make/O/D/N=(Nch) betaExtra_best
    betaExtra_best = NaN

    Make/FREE/D/N=1 B_one
    B_one[0] = Bval_T

    String tmpGeomFolder   = "root:Packages:SNS:tmpVortexGeom2D"
    String tmpSolverFolder = "root:Packages:SNS:tmpVortexMin2D"
    String tmpPhaseFolder  = "root:Packages:SNS:tmpVortexPhase2D"

    NewDataFolder/O root:Packages
    NewDataFolder/O root:Packages:SNS
    NewDataFolder/O $tmpGeomFolder
    NewDataFolder/O $tmpSolverFolder
    NewDataFolder/O $tmpPhaseFolder

    // --- precompute boundary geometry maps once ---
    SNS_ComputeBoundaryCurvatureFromMask2D(Nmask, tmpGeomFolder, nSmoothIter=nSmoothIter)

    Wave BoundaryMask_pre = $(tmpGeomFolder + ":BoundaryMask")
    Wave NxMap_pre        = $(tmpGeomFolder + ":NxMap")
    Wave NyMap_pre        = $(tmpGeomFolder + ":NyMap")

    String nameE2D = tmpSolverFolder + ":E2D_tmp"
    String nameM   = tmpSolverFolder + ":M_tmp"
    String nameS   = tmpSolverFolder + ":S_tmp"

    Variable iTry, j, k
    Variable xV_nm, yV_nm
    Variable Lch, Wch, Teff, weight
    Variable betaExtra
    Variable nBr, Eval
    Variable costABS_here, costPhase_here, costHere
    Variable bestCost_local = NaN
    Variable bestX_local = NaN
    Variable bestY_local = NaN
    Variable bestIdx_local = NaN

    for (iTry = 0; iTry < Ntry; iTry += 1)

        xV_nm = xTry_nm[iTry]
        yV_nm = yTry_nm[iTry]

        // --- solve exterior phase field using precomputed boundary maps ---
        SNS_SolveExteriorPhaseField2D_FromBoundaryMaps(Nmask, BoundaryMask_pre, NxMap_pre, NyMap_pre, \
            xV_nm, yV_nm, tmpPhaseFolder, \
            nIterMax=phaseNIterMax, tol_rad=phaseTol_rad, omega=phaseOmega)

        Wave extMask_tmp     = $(tmpPhaseFolder + ":extMask")
        Wave phaseRe_ext_tmp = $(tmpPhaseFolder + ":phaseRe_ext")
        Wave phaseIm_ext_tmp = $(tmpPhaseFolder + ":phaseIm_ext")

        // --- compute per-channel betaExtra from solved exterior phase ---
        SNS_ComputeBetaExtraFromExteriorPhase2D(extMask_tmp, phaseRe_ext_tmp, phaseIm_ext_tmp, \
            Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, tmpPhaseFolder, \
            nFlux=nFlux, searchRad_px=phaseSearchRad_px)

        Wave betaExtra_List_tmp = $(tmpPhaseFolder + ":betaExtra_List")

        // --- ABS contribution (optional; off by default) ---
        costABS_here = 0

        if (useABS != 0)
            for (j = 0; j < Nch; j += 1)

                weight = wChan[j] / sumW
                if (weight <= 0)
                    continue
                endif

                Lch  = L_N_List[j] * 1e-9
                Wch  = W_eff_List[j] * 1e-9
                Teff = T_eff_List[j]

                if (numtype(Lch) || numtype(Wch) || numtype(Teff))
                    continue
                endif

                betaExtra = betaExtra_List_tmp[j]
                if (numtype(betaExtra) != 0)
                    continue
                endif

                nBr = Solve_AllBranches_SNS_dGSJ_betaExtra(B_one, Lch, Wch, Delta, vF, lambdaL, Teff, betaExtra, nameE2D, nameM, nameS)
                if (nBr <= 0)
                    continue
                endif

                Wave E2D_tmp = $nameE2D

                for (k = 0; k < nBr; k += 1)
                    Eval = E2D_tmp[0][k]
                    if (numtype(Eval) == 0)
                        costABS_here -= weight * abs(Eval)
                    endif
                endfor

            endfor
        endif

        // --- phase-field contribution from solved exterior phase ---
        costPhase_here = alphaPhase_eV * SNS_ComputeExteriorPhaseEnergy2D(extMask_tmp, phaseRe_ext_tmp, phaseIm_ext_tmp)

        costHere = costABS_here + costPhase_here

        CostABS_try[iTry]   = costABS_here
        CostPhase_try[iTry] = costPhase_here
        Cost_try[iTry]      = costHere

        if ((numtype(bestCost_local) != 0) || (costHere < bestCost_local))
            bestCost_local = costHere
            bestX_local = xV_nm
            bestY_local = yV_nm
            bestIdx_local = iTry

            betaExtra_best[] = betaExtra_List_tmp[p]
        endif

    endfor

    Variable/G xV_best_nm = bestX_local
    Variable/G yV_best_nm = bestY_local
    Variable/G cost_best  = bestCost_local
    Variable/G idx_best   = bestIdx_local

    KillDataFolder/Z $tmpGeomFolder
    KillDataFolder/Z $tmpSolverFolder
    KillDataFolder/Z $tmpPhaseFolder
    SetDataFolder $oldDF

    return 0
End








//==============================================================================
// SNS_SolveExteriorPhaseField2D_FromBoundaryMaps
//
// Purpose
//   Solve a boundary-conditioned vortex phase field on the exterior of a 2D
//   hole mask using precomputed boundary maps.
//
//   This is the optimized version of SNS_SolveExteriorPhaseField2D(...):
//   the boundary geometry (BoundaryMask, NxMap, NyMap) is supplied directly,
//   so it does not need to be recomputed for every trial vortex position.
//
// Inputs
//   Nmask          : 2D hole mask; 1 inside hole, 0 in exterior.
//   BoundaryMask   : Boundary mask from SNS_ComputeBoundaryCurvatureFromMask2D(...)
//   NxMap          : x-component of local unit normal
//   NyMap          : y-component of local unit normal
//   xV_nm          : Trial vortex x position [nm]
//   yV_nm          : Trial vortex y position [nm]
//   folder         : Output data folder
//
// Optional Inputs
//   nIterMax       : Maximum number of SOR iterations. Default = 500.
//   tol_rad        : Convergence tolerance on max |Δchi|. Default = 1e-6.
//   omega          : SOR relaxation parameter. Default = 1.6.
//
// Outputs (in folder)
//   extMask        : Exterior-domain mask
//   phi0_ext       : Free-space vortex phase [rad]
//   chi_ext        : Harmonic correction [rad]
//   phi_ext        : Total phase [rad]
//   phaseRe_ext    : cos(phi_ext)
//   phaseIm_ext    : sin(phi_ext)
//   phaseIterUsed  : Number of iterations used
//   phaseMaxDelta  : Final max |Δchi|
//
// Notes
//   • Proper boundary condition is imposed as
//         ∂n chi = - ∂n phi0
//   • Wrapped gradients of phi0 are computed from its complex representation.
//==============================================================================
Function SNS_SolveExteriorPhaseField2D_FromBoundaryMaps(Nmask, BoundaryMask, NxMap, NyMap, xV_nm, yV_nm, folder, [nIterMax, tol_rad, omega])

    Wave Nmask, BoundaryMask, NxMap, NyMap
    Variable xV_nm, yV_nm
    String folder
    Variable nIterMax, tol_rad, omega

    if (ParamIsDefault(nIterMax))
        nIterMax = 500
    endif
    if (ParamIsDefault(tol_rad))
        tol_rad = 1e-6
    endif
    if (ParamIsDefault(omega))
        omega = 1.6
    endif

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if ((nx <= 2) || (ny <= 2))
        Abort "SNS_SolveExteriorPhaseField2D_FromBoundaryMaps: invalid mask dimensions."
    endif

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    Variable absdx = abs(dx)
    Variable absdy = abs(dy)

    String oldDF = GetDataFolder(1)
    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    Make/O/D/N=(nx, ny) extMask, phi0_ext, chi_ext, phi_ext, phaseRe_ext, phaseIm_ext
    Make/FREE/D/N=(nx, ny) phaseRe0_ext, phaseIm0_ext

    Variable ix, iy, xnm, ynm

    extMask = (Nmask <= 0.5) ? 1 : 0

    for (ix = 0; ix < nx; ix += 1)
        xnm = x0 + ix*dx
        for (iy = 0; iy < ny; iy += 1)
            ynm = y0 + iy*dy
            phi0_ext[ix][iy] = atan2(ynm - yV_nm, xnm - xV_nm)
        endfor
    endfor

    phaseRe0_ext = cos(phi0_ext)
    phaseIm0_ext = sin(phi0_ext)

    chi_ext = (extMask > 0.5) ? 0 : NaN

    Variable iter, maxDelta, delta, chiOld, chiNew
    Variable nHole, denom, sumChiExt, sumDrive
    Variable ixm, ixp, iym, iyp
    Variable reC, imC, reP, imP, reM, imM
    Variable dphiXP, dphiXM, dphiYP, dphiYM
    Variable dphi0_dx, dphi0_dy, dphi0_dn
    Variable nxLoc, nyLoc, dsEff

    for (iter = 0; iter < nIterMax; iter += 1)

        maxDelta = 0

        for (ix = 1; ix < nx-1; ix += 1)
            for (iy = 1; iy < ny-1; iy += 1)

                if (extMask[ix][iy] <= 0.5)
                    continue
                endif

                reC = phaseRe0_ext[ix][iy]
                imC = phaseIm0_ext[ix][iy]

                reP = phaseRe0_ext[ix+1][iy]
                imP = phaseIm0_ext[ix+1][iy]
                reM = phaseRe0_ext[ix-1][iy]
                imM = phaseIm0_ext[ix-1][iy]

                dphiXP = atan2(imP*reC - reP*imC, reP*reC + imP*imC)
                dphiXM = atan2(imC*reM - reC*imM, reC*reM + imC*imM)
                dphi0_dx = (dphiXP + dphiXM) / (2*absdx)

                reP = phaseRe0_ext[ix][iy+1]
                imP = phaseIm0_ext[ix][iy+1]
                reM = phaseRe0_ext[ix][iy-1]
                imM = phaseIm0_ext[ix][iy-1]

                dphiYP = atan2(imP*reC - reP*imC, reP*reC + imP*imC)
                dphiYM = atan2(imC*reM - reC*imM, reC*reM + imC*imM)
                dphi0_dy = (dphiYP + dphiYM) / (2*absdy)

                nxLoc = NxMap[ix][iy]
                nyLoc = NyMap[ix][iy]

                dphi0_dn = dphi0_dx*nxLoc + dphi0_dy*nyLoc

                dsEff = abs(nxLoc)*absdx + abs(nyLoc)*absdy
                if (dsEff <= 0)
                    dsEff = min(absdx, absdy)
                endif

                sumChiExt = 0
                sumDrive  = 0
                nHole     = 0

                ixm = ix - 1
                if (extMask[ixm][iy] > 0.5)
                    sumChiExt += chi_ext[ixm][iy]
                else
                    nHole += 1
                    sumDrive += -dsEff * dphi0_dn
                endif

                ixp = ix + 1
                if (extMask[ixp][iy] > 0.5)
                    sumChiExt += chi_ext[ixp][iy]
                else
                    nHole += 1
                    sumDrive += -dsEff * dphi0_dn
                endif

                iym = iy - 1
                if (extMask[ix][iym] > 0.5)
                    sumChiExt += chi_ext[ix][iym]
                else
                    nHole += 1
                    sumDrive += -dsEff * dphi0_dn
                endif

                iyp = iy + 1
                if (extMask[ix][iyp] > 0.5)
                    sumChiExt += chi_ext[ix][iyp]
                else
                    nHole += 1
                    sumDrive += -dsEff * dphi0_dn
                endif

                denom = 4 - nHole
                if (denom <= 0)
                    continue
                endif

                chiOld = chi_ext[ix][iy]
                chiNew = (sumChiExt + sumDrive) / denom
                chiNew = (1 - omega)*chiOld + omega*chiNew

                delta = abs(chiNew - chiOld)
                if (delta > maxDelta)
                    maxDelta = delta
                endif

                chi_ext[ix][iy] = chiNew

            endfor
        endfor

        if (maxDelta < tol_rad)
            break
        endif
    endfor

    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)

            if (extMask[ix][iy] <= 0.5)
                chi_ext[ix][iy] = NaN
                phi_ext[ix][iy] = NaN
                phaseRe_ext[ix][iy] = NaN
                phaseIm_ext[ix][iy] = NaN
                continue
            endif

            if ((ix == 0) || (ix == nx-1) || (iy == 0) || (iy == ny-1))
                chi_ext[ix][iy] = 0
            endif

            phi_ext[ix][iy] = phi0_ext[ix][iy] + chi_ext[ix][iy]
            phaseRe_ext[ix][iy] = cos(phi_ext[ix][iy])
            phaseIm_ext[ix][iy] = sin(phi_ext[ix][iy])

        endfor
    endfor

    Variable/G phaseIterUsed = iter + 1
    Variable/G phaseMaxDelta = maxDelta

    SetScale/P x, x0, dx, "nm", extMask, phi0_ext, chi_ext, phi_ext, phaseRe_ext, phaseIm_ext
    SetScale/P y, y0, dy, "nm", extMask, phi0_ext, chi_ext, phi_ext, phaseRe_ext, phaseIm_ext

    SetDataFolder $oldDF
    return 0
End

//==============================================================================
// SNS_WrapToPi
//
// Purpose
//   Wrap an angle to the interval [-pi, pi].
//
// Inputs
//   phi      : Input angle [rad]
//
// Returns
//   phiWrap  : Wrapped angle [rad]
//==============================================================================
Function SNS_WrapToPi(phi)

    Variable phi
    Variable phiWrap = phi

    do
        if (phiWrap > pi)
            phiWrap -= 2*pi
        elseif (phiWrap < -pi)
            phiWrap += 2*pi
        else
            break
        endif
    while (1)

    return phiWrap
End


//==============================================================================
// SNS_SolveExteriorPhaseField2D
//
// Purpose
//   Solve a boundary-conditioned vortex phase field on the exterior of a 2D
//   hole/N-region mask using:
//
//       phaseCorr = phaseFree + chi
//
//   where:
//
//       phaseFree = nFlux * atan2(y-yV, x-xV)
//
//   is the free-space point-vortex phase, and chi is a harmonic correction
//   solved only on the exterior superconducting domain.
//
// Mask convention:
//       Nmask > 0.5  -> excluded N/hole region
//       Nmask <= 0.5 -> exterior superconducting solve domain
//
// Boundary condition at excluded N/hole boundary:
//       ∂n(phaseFree + chi) = leakAlpha * ∂n(phaseFree)
//
//   equivalently:
//
//       ∂n(chi) = -(1 - leakAlpha) * ∂n(phaseFree)
//
//   leakAlpha:
//       0 -> hard zero-normal-current boundary, original behavior
//       1 -> no boundary correction, free-space vortex field
//       0<leakAlpha<1 -> partial normal-current leakage
//
// Outer box edge:
//       Dirichlet chi = 0
//
// Optional Inputs:
//   nIterMax       : Maximum number of SOR iterations. Default = 500.
//   tol_rad        : Convergence tolerance on max |Δchi| [rad]. Default = 1e-6.
//   omega          : SOR relaxation parameter. Default = 1.3.
//   nSmoothIter    : Smoothing iterations for diagnostic boundary maps. Default = 2.
//   nFlux          : Vortex winding number. Default from SNS settings, otherwise +1.
//   leakAlpha      : Boundary leakage/transparency parameter. Default = 0.
//
// Outputs:
//   Exterior-only / masked fields:
//       extMask
//       phaseFree_ext
//       phaseReFree_ext, phaseImFree_ext
//       chi_ext
//       phaseCorr_ext
//       phaseReCorr_ext, phaseImCorr_ext
//
//   Full fields, defined everywhere:
//       fullMask
//       phaseFree
//       phaseReFree, phaseImFree
//       phaseCorr
//       phaseReCorr, phaseImCorr
//
//   Scalars:
//       phaseIterUsed
//       phaseMaxDelta
//       phaseLeakAlpha
//
//   Diagnostics from SNS_ComputeBoundaryCurvatureFromMask2D:
//       BoundaryMask, NxMap, NyMap, KappaMap
//==============================================================================
Function SNS_SolveExteriorPhaseField2D(Nmask, xV_nm, yV_nm, folder, [nIterMax, tol_rad, omega, nSmoothIter, nFlux, leakAlpha])

    Wave Nmask
    Variable xV_nm, yV_nm
    String folder
    Variable nIterMax, tol_rad, omega, nSmoothIter, nFlux
    Variable leakAlpha

    if (ParamIsDefault(nFlux))
        NVAR/Z n_flux_set = root:SNS_Settings:SNS_nFlux
        if (NVAR_Exists(n_flux_set))
            nFlux = n_flux_set
        else
            nFlux = 1
        endif
    endif

    if (ParamIsDefault(nIterMax))
        nIterMax = 500
    endif
    if (ParamIsDefault(tol_rad))
        tol_rad = 1e-6
    endif
    if (ParamIsDefault(omega))
        omega = 1.3
    endif
    if (ParamIsDefault(nSmoothIter))
        nSmoothIter = 2
    endif
    if (ParamIsDefault(leakAlpha))
        leakAlpha = 0
    endif

    if (nIterMax < 1)
        Abort "SNS_SolveExteriorPhaseField2D: nIterMax must be >= 1."
    endif
    if (tol_rad <= 0)
        Abort "SNS_SolveExteriorPhaseField2D: tol_rad must be > 0."
    endif
    if ((omega <= 0) || (omega >= 2))
        Abort "SNS_SolveExteriorPhaseField2D: omega must satisfy 0 < omega < 2."
    endif
    if (nSmoothIter < 0)
        Abort "SNS_SolveExteriorPhaseField2D: nSmoothIter must be >= 0."
    endif
    if ((leakAlpha < 0) || (leakAlpha > 1))
        Abort "SNS_SolveExteriorPhaseField2D: leakAlpha must satisfy 0 <= leakAlpha <= 1."
    endif

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if ((nx <= 2) || (ny <= 2))
        Abort "SNS_SolveExteriorPhaseField2D: invalid mask dimensions."
    endif

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    Variable absdx = abs(dx)
    Variable absdy = abs(dy)

    if ((absdx <= 0) || (absdy <= 0))
        Abort "SNS_SolveExteriorPhaseField2D: invalid mask scaling."
    endif

    String oldDF = GetDataFolder(1)

    SNS_MaskAreaPerim_FromParticles(Nmask)
    Wave w
    Make/O/N=1 w_area_nm2 = w[0]
    Make/O/N=1 w_perim_nm = w[1]

    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    Make/O/D/N=(nx, ny) extMask, fullMask
    Make/O/D/N=(nx, ny) phaseFree_ext, phaseFree
    Make/O/D/N=(nx, ny) phaseReFree_ext, phaseImFree_ext
    Make/O/D/N=(nx, ny) phaseReFree, phaseImFree
    Make/O/D/N=(nx, ny) chi_ext
    Make/O/D/N=(nx, ny) phaseCorr_ext, phaseReCorr_ext, phaseImCorr_ext
    Make/O/D/N=(nx, ny) phaseCorr, phaseReCorr, phaseImCorr

    Variable ix, iy
    Variable xnm, ynm

    extMask  = (Nmask <= 0.5) ? 1 : 0
    fullMask = 1

    SNS_ComputeBoundaryCurvatureFromMask2D(Nmask, "", nSmoothIter=nSmoothIter)

    // Free-space vortex phase on full grid.
    for (ix = 0; ix < nx; ix += 1)
        xnm = x0 + ix*dx
        for (iy = 0; iy < ny; iy += 1)
            ynm = y0 + iy*dy
            phaseFree[ix][iy] = nFlux * atan2(ynm - yV_nm, xnm - xV_nm)
        endfor
    endfor

    phaseFree_ext = (extMask > 0.5) ? phaseFree : NaN

    phaseReFree = cos(phaseFree)
    phaseImFree = sin(phaseFree)

    phaseReFree_ext = (extMask > 0.5) ? phaseReFree : NaN
    phaseImFree_ext = (extMask > 0.5) ? phaseImFree : NaN

    // Initial chi: solved only outside N; NaN inside N.
    chi_ext = (extMask > 0.5) ? 0 : NaN

    Variable iter, iterUsed
    Variable maxDelta, delta, chiOld, chiNew
    Variable nHole, denom
    Variable sumChiExt, sumDrive
    Variable ixm, ixp, iym, iyp

    Variable reC, imC, reP, imP, reM, imM
    Variable dFreeXP, dFreeXM, dFreeYP, dFreeYM
    Variable dPhaseFree_dx, dPhaseFree_dy

    Variable sumChi, nChi, meanChi, gaugeDelta

    Variable leakFac = 1 - leakAlpha

    iterUsed = 0

    for (iter = 0; iter < nIterMax; iter += 1)

        maxDelta = 0
        iterUsed = iter + 1

        for (ix = 1; ix < nx-1; ix += 1)
            for (iy = 1; iy < ny-1; iy += 1)

                if (extMask[ix][iy] <= 0.5)
                    continue
                endif

                reC = phaseReFree[ix][iy]
                imC = phaseImFree[ix][iy]

                reP = phaseReFree[ix+1][iy]
                imP = phaseImFree[ix+1][iy]
                reM = phaseReFree[ix-1][iy]
                imM = phaseImFree[ix-1][iy]

                dFreeXP = atan2(imP*reC - reP*imC, reP*reC + imP*imC)
                dFreeXM = atan2(imC*reM - reC*imM, reC*reM + imC*imM)
                dPhaseFree_dx = (dFreeXP + dFreeXM) / (2*absdx)

                reP = phaseReFree[ix][iy+1]
                imP = phaseImFree[ix][iy+1]
                reM = phaseReFree[ix][iy-1]
                imM = phaseImFree[ix][iy-1]

                dFreeYP = atan2(imP*reC - reP*imC, reP*reC + imP*imC)
                dFreeYM = atan2(imC*reM - reC*imM, reC*reM + imC*imM)
                dPhaseFree_dy = (dFreeYP + dFreeYM) / (2*absdy)

                sumChiExt = 0
                sumDrive  = 0
                nHole     = 0

                // left neighbor
                ixm = ix - 1
                if (extMask[ixm][iy] > 0.5)
                    sumChiExt += chi_ext[ixm][iy]
                else
                    nHole += 1
                    sumDrive += +leakFac * absdx * dPhaseFree_dx
                endif

                // right neighbor
                ixp = ix + 1
                if (extMask[ixp][iy] > 0.5)
                    sumChiExt += chi_ext[ixp][iy]
                else
                    nHole += 1
                    sumDrive += -leakFac * absdx * dPhaseFree_dx
                endif

                // down neighbor
                iym = iy - 1
                if (extMask[ix][iym] > 0.5)
                    sumChiExt += chi_ext[ix][iym]
                else
                    nHole += 1
                    sumDrive += +leakFac * absdy * dPhaseFree_dy
                endif

                // up neighbor
                iyp = iy + 1
                if (extMask[ix][iyp] > 0.5)
                    sumChiExt += chi_ext[ix][iyp]
                else
                    nHole += 1
                    sumDrive += -leakFac * absdy * dPhaseFree_dy
                endif

                denom = 4 - nHole
                if (denom <= 0)
                    continue
                endif

                chiOld = chi_ext[ix][iy]
                chiNew = (sumChiExt + sumDrive) / denom
                chiNew = (1 - omega)*chiOld + omega*chiNew

                delta = abs(chiNew - chiOld)
                if (delta > maxDelta)
                    maxDelta = delta
                endif

                chi_ext[ix][iy] = chiNew

            endfor
        endfor

        // Remove constant null-mode of chi.
        sumChi = 0
        nChi   = 0

        for (ix = 1; ix < nx-1; ix += 1)
            for (iy = 1; iy < ny-1; iy += 1)
                if ((extMask[ix][iy] > 0.5) && (numtype(chi_ext[ix][iy]) == 0))
                    sumChi += chi_ext[ix][iy]
                    nChi += 1
                endif
            endfor
        endfor

        if (nChi > 0)
            meanChi = sumChi / nChi
            gaugeDelta = abs(meanChi)

            for (ix = 1; ix < nx-1; ix += 1)
                for (iy = 1; iy < ny-1; iy += 1)
                    if ((extMask[ix][iy] > 0.5) && (numtype(chi_ext[ix][iy]) == 0))
                        chi_ext[ix][iy] -= meanChi
                    endif
                endfor
            endfor

            if (gaugeDelta > maxDelta)
                maxDelta = gaugeDelta
            endif
        endif

        if (maxDelta < tol_rad)
            break
        endif

    endfor

    // Build masked exterior field and full field.
    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)

            if ((extMask[ix][iy] > 0.5) && (numtype(chi_ext[ix][iy]) == 0))

                if ((ix == 0) || (ix == nx-1) || (iy == 0) || (iy == ny-1))
                    chi_ext[ix][iy] = 0
                endif

                phaseCorr_ext[ix][iy] = phaseFree[ix][iy] + chi_ext[ix][iy]

                phaseReCorr_ext[ix][iy] = cos(phaseCorr_ext[ix][iy])
                phaseImCorr_ext[ix][iy] = sin(phaseCorr_ext[ix][iy])

                phaseCorr[ix][iy] = phaseCorr_ext[ix][iy]

            else

                chi_ext[ix][iy] = NaN
                phaseCorr_ext[ix][iy] = NaN
                phaseReCorr_ext[ix][iy] = NaN
                phaseImCorr_ext[ix][iy] = NaN

                // No imposed zero phase inside N:
                // full diagnostic/sampling phase uses free vortex phase there.
                phaseCorr[ix][iy] = phaseFree[ix][iy]

            endif

            phaseReCorr[ix][iy] = cos(phaseCorr[ix][iy])
            phaseImCorr[ix][iy] = sin(phaseCorr[ix][iy])

        endfor
    endfor

    Variable/G phaseIterUsed  = iterUsed
    Variable/G phaseMaxDelta  = maxDelta
    Variable/G phaseLeakAlpha = leakAlpha

    SetScale/P x, x0, dx, "nm", extMask, fullMask, phaseFree_ext, phaseFree, phaseReFree_ext, phaseImFree_ext, phaseReFree, phaseImFree, chi_ext, phaseCorr_ext, phaseReCorr_ext, phaseImCorr_ext, phaseCorr, phaseReCorr, phaseImCorr
    SetScale/P y, y0, dy, "nm", extMask, fullMask, phaseFree_ext, phaseFree, phaseReFree_ext, phaseImFree_ext, phaseReFree, phaseImFree, chi_ext, phaseCorr_ext, phaseReCorr_ext, phaseImCorr_ext, phaseCorr, phaseReCorr, phaseImCorr

    SetDataFolder $oldDF
    return 0
End

//==============================================================================
// SNS_HarmonicFaceCoeff
//==============================================================================
Function SNS_HarmonicFaceCoeff(a1, a2)
    Variable a1, a2

    if ((a1 <= 0) || (a2 <= 0))
        return 0
    endif

    return 2*a1*a2/(a1 + a2)
End


//==============================================================================
// SNS_SolveStiffnessPhaseField
//
// Purpose:
//   Solve a geometry-aware vortex phase field for an N island embedded in a
//   superconducting region, using a variable 2D phase stiffness.
//
//   The solved field is:
//
//       PhaseCorr = PhaseFree + PhaseChi
//
//   with:
//
//       PhaseFree = nFlux * atan2(y-yV, x-xV)
//
//   and PhaseChi chosen such that:
//
//       div[ PhaseStiffness(x,y) * grad(PhaseCorr) ] = 0
//
//   away from the vortex singularity.
//
//   Phase stiffness vs superfluid density:
//       The coefficient PhaseStiffness is the phase-stiffness prefactor K(r)
//       in the phase-only free energy
//
//           F = 1/2 ∫ K(r) |grad(phi)|^2 d^2r.
//
//       Varying this functional gives
//
//           div[ K(r) grad(phi) ] = 0.
//
//       Microscopically, K is proportional to n_s/m*, where n_s is the
//       superfluid density and m* is the effective mass. If m* is spatially
//       uniform, PhaseStiffness is proportional to superfluid density. In this
//       phenomenological solver only relative stiffness matters, so etaN can be
//       interpreted either as relative phase stiffness or relative superfluid
//       density.
//
//   This is a phase-only stiffness model. It does not include a dynamical vector
//   potential, London/Pearl magnetic screening, or a separate in-plane orbital
//   phase term.
//
// Mask convention:
//   Default, maskNIsOne = 1:
//
//       Nmask > 0.5  -> weakened N island
//       Nmask <= 0.5 -> superconducting region
//
//   If your mask convention is inverted, call with:
//
//       maskNIsOne = 0
//
// Stiffness model:
//       PhaseStiffnessRaw = 1      in SC
//       PhaseStiffnessRaw = etaN   in N
//
//   etaN = 0      : hard excluded / zero-stiffness island
//   etaN << 1     : weak proximity stiffness in N
//   etaN = 1      : uniform stiffness; mask has no effect
//
//   The raw stiffness step is then smoothed over coherenceLength_nm to avoid
//   unphysical pixel-sharp stiffness changes. The smoothing is implemented by
//   repeated 3x3 Gaussian MatrixFilter passes:
//
//       MatrixFilter /N=3 gauss PhaseStiffness
//
//   A single 3x3 Gaussian pass has an effective sigma ~ 0.7 pixel. Repeating
//   the filter nSmoothPass times gives approximately:
//
//       sigma_px ~ sqrt(nSmoothPass / 2)
//
//   so the number of passes is estimated as:
//
//       nSmoothPass ~ 2 * (coherenceLength_nm / pixelSize_nm)^2
//
//   where pixelSize_nm is the average of |dx| and |dy|.
//
// Numerical method:
//   Finite-volume relaxation using harmonic face averages of PhaseStiffness.
//   This makes the island boundary enter through the face couplings rather than
//   through a hard phase value inside N.
//
// Inputs:
//   Nmask              : 2D mask wave with x/y scaling in nm.
//   xV_nm, yV_nm       : vortex position [nm].
//   folder             : output data folder.
//
// Optional inputs:
//   etaN               : N-region stiffness relative to SC. Default = 0.02.
//   nIterMax           : maximum SOR iterations. Default = 5000.
//   tol                : convergence tolerance for max |ΔPhaseChi|. Default = 1e-8.
//   omega              : SOR relaxation parameter. Default = 1.4.
//   nFlux              : vortex winding. Default from root:SNS_Settings:SNS_nFlux,
//                        otherwise +1.
//   maskNIsOne         : mask convention flag. Default = 1.
//   coherenceLength_nm : stiffness smoothing length [nm]. Default = 10.
//                        Set <=0 to disable smoothing.
//
// Outputs:
//   PhaseStiffnessRaw
//       Unsmoothed relative phase stiffness from the mask.
//
//   PhaseStiffness
//       Smoothed local relative phase stiffness used in the finite-volume solve.
//
//   PhaseFree
//       Uncorrected free-vortex phase field:
//           nFlux * atan2(y-yV, x-xV)
//
//   PhaseReFree, PhaseImFree
//       Complex representation of PhaseFree.
//
//   PhaseChi
//       Harmonic correction generated by the variable-stiffness geometry.
//
//   PhaseCorr
//       Geometry-corrected scalar phase field:
//           PhaseCorr = PhaseFree + PhaseChi
//
//   PhaseReCorr, PhaseImCorr
//       Complex representation of PhaseCorr.
//
//   QxFree, QyFree, QmagFree
//       Gradient of the free vortex phase field, in rad/nm.
//
//   QxCorr, QyCorr, QmagCorr
//       Gradient of the corrected phase field, in rad/nm.
//
//   JxFree, JyFree, JmagFree
//       Stiffness-weighted current-like field from the free phase.
//
//   JxCorr, JyCorr, JmagCorr
//       Stiffness-weighted current-like field from the corrected phase.
//
//   StiffIterUsed
//   StiffMaxDelta
//   StiffEtaN
//   StiffMaskNIsOne
//   StiffCoherenceLength_nm
//   StiffSmoothPasses
//
// ABS usage:
//   For the current scalar phase-field ABS workflow, use:
//
//       betaExtra[ch] = PhaseCorr(hit2) - PhaseCorr(hit1)
//
//   computed as a wrapped phase difference using PhaseReCorr/PhaseImCorr.
//
//   Do not also add the separate in-plane orbital B*w*h_eff phase term when
//   using this vortex phase as the ABS phase contribution.
//==============================================================================
Function SNS_SolveStiffnessPhaseField(Nmask, xV_nm, yV_nm, folder, [etaN, nIterMax, tol, omega, nFlux, maskNIsOne, coherenceLength_nm])

    Wave Nmask
    Variable xV_nm, yV_nm
    String folder
    Variable etaN, nIterMax, tol, omega, nFlux, maskNIsOne
    Variable coherenceLength_nm

    if (ParamIsDefault(etaN))
        etaN = 0.02
    endif
    if (ParamIsDefault(nIterMax))
        nIterMax = 5000
    endif
    if (ParamIsDefault(tol))
        tol = 1e-8
    endif
    if (ParamIsDefault(omega))
        omega = 1.4
    endif
    if (ParamIsDefault(maskNIsOne))
        maskNIsOne = 1
    endif
    if (ParamIsDefault(coherenceLength_nm))
        coherenceLength_nm = 10
    endif
    if (ParamIsDefault(nFlux))
        NVAR/Z n_flux_set = root:SNS_Settings:SNS_nFlux
        if (NVAR_Exists(n_flux_set))
            nFlux = n_flux_set
        else
            nFlux = 1
        endif
    endif

    if ((etaN < 0) || (etaN > 1))
        Abort "SNS_SolveStiffnessPhaseField: etaN must satisfy 0 <= etaN <= 1."
    endif
    if ((omega <= 0) || (omega >= 2))
        Abort "SNS_SolveStiffnessPhaseField: omega must satisfy 0 < omega < 2."
    endif

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    Variable absdx = abs(dx)
    Variable absdy = abs(dy)

    if ((nx <= 2) || (ny <= 2))
        Abort "SNS_SolveStiffnessPhaseField: invalid mask dimensions."
    endif
    if ((absdx <= 0) || (absdy <= 0))
        Abort "SNS_SolveStiffnessPhaseField: invalid mask scaling."
    endif

    String oldDF = GetDataFolder(1)

    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    Make/O/D/N=(nx, ny) PhaseStiffnessRaw, PhaseStiffness
    Make/O/D/N=(nx, ny) PhaseFree, PhaseReFree, PhaseImFree
    Make/O/D/N=(nx, ny) PhaseChi
    Make/O/D/N=(nx, ny) PhaseCorr, PhaseReCorr, PhaseImCorr
    Make/O/D/N=(nx, ny) QxFree, QyFree, QmagFree
    Make/O/D/N=(nx, ny) QxCorr, QyCorr, QmagCorr
    Make/O/D/N=(nx, ny) JxFree, JyFree, JmagFree
    Make/O/D/N=(nx, ny) JxCorr, JyCorr, JmagCorr

    Variable ix, iy
    Variable xnm, ynm
    Variable rx, ry

    // Stiffness map before coherence-length smoothing.
    if (maskNIsOne)
        PhaseStiffnessRaw = (Nmask[p][q] > 0.5) ? etaN : 1
    else
        PhaseStiffnessRaw = (Nmask[p][q] > 0.5) ? 1 : etaN
    endif

    Duplicate/O PhaseStiffnessRaw, PhaseStiffness

    // Smooth stiffness over coherenceLength_nm using repeated 3x3 Gaussian filtering.
    Variable pixelSize_nm = 0.5 * (absdx + absdy)
    Variable smoothSigma_px
    Variable nSmoothPass
    Variable iSmooth

    if (coherenceLength_nm > 0)
        smoothSigma_px = coherenceLength_nm / pixelSize_nm
        nSmoothPass = max(1, round(2 * smoothSigma_px * smoothSigma_px))

        for (iSmooth = 0; iSmooth < nSmoothPass; iSmooth += 1)
            MatrixFilter /N=3 gauss PhaseStiffness
        endfor

        // Guard against any numerical overshoot at boundaries.
        PhaseStiffness = max(etaN, min(1, PhaseStiffness))
    else
        nSmoothPass = 0
    endif

    // Free vortex phase.
    for (ix = 0; ix < nx; ix += 1)
        xnm = x0 + ix*dx
        for (iy = 0; iy < ny; iy += 1)
            ynm = y0 + iy*dy
            rx = xnm - xV_nm
            ry = ynm - yV_nm

            PhaseFree[ix][iy]   = nFlux * atan2(ry, rx)
            PhaseReFree[ix][iy] = cos(PhaseFree[ix][iy])
            PhaseImFree[ix][iy] = sin(PhaseFree[ix][iy])
        endfor
    endfor

    PhaseChi = 0

    Variable invDx2 = 1/(absdx*absdx)
    Variable invDy2 = 1/(absdy*absdy)

    Variable iter, iterUsed
    Variable maxDelta, delta
    Variable chiOld, chiNew
    Variable aC, aR, aL, aU, aD
    Variable denom

    Variable reC, imC, reR, imR, reL, imL, reU, imU, reD, imD
    Variable dFreeR, dFreeL, dFreeU, dFreeD

    iterUsed = 0

    for (iter = 0; iter < nIterMax; iter += 1)

        maxDelta = 0
        iterUsed = iter + 1

        // Outer boundary: PhaseChi = 0
        PhaseChi[0][]    = 0
        PhaseChi[nx-1][] = 0
        PhaseChi[][0]    = 0
        PhaseChi[][ny-1] = 0

        for (ix = 1; ix < nx-1; ix += 1)
            for (iy = 1; iy < ny-1; iy += 1)

                aC = PhaseStiffness[ix][iy]

                if (aC <= 0)
                    PhaseChi[ix][iy] = 0
                    continue
                endif

                aR = SNS_HarmonicFaceCoeff(aC, PhaseStiffness[ix+1][iy])
                aL = SNS_HarmonicFaceCoeff(aC, PhaseStiffness[ix-1][iy])
                aU = SNS_HarmonicFaceCoeff(aC, PhaseStiffness[ix][iy+1])
                aD = SNS_HarmonicFaceCoeff(aC, PhaseStiffness[ix][iy-1])

                denom = (aR + aL)*invDx2 + (aU + aD)*invDy2

                if (denom <= 0)
                    continue
                endif

                reC = PhaseReFree[ix][iy]
                imC = PhaseImFree[ix][iy]

                reR = PhaseReFree[ix+1][iy]
                imR = PhaseImFree[ix+1][iy]
                reL = PhaseReFree[ix-1][iy]
                imL = PhaseImFree[ix-1][iy]
                reU = PhaseReFree[ix][iy+1]
                imU = PhaseImFree[ix][iy+1]
                reD = PhaseReFree[ix][iy-1]
                imD = PhaseImFree[ix][iy-1]

                dFreeR = atan2(imR*reC - reR*imC, reR*reC + imR*imC)
                dFreeL = atan2(imL*reC - reL*imC, reL*reC + imL*imC)
                dFreeU = atan2(imU*reC - reU*imC, reU*reC + imU*imC)
                dFreeD = atan2(imD*reC - reD*imC, reD*reC + imD*imC)

                chiNew = (aR*invDx2*(PhaseChi[ix+1][iy] + dFreeR) \
                        + aL*invDx2*(PhaseChi[ix-1][iy] + dFreeL) \
                        + aU*invDy2*(PhaseChi[ix][iy+1] + dFreeU) \
                        + aD*invDy2*(PhaseChi[ix][iy-1] + dFreeD)) / denom

                chiOld = PhaseChi[ix][iy]
                chiNew = (1 - omega)*chiOld + omega*chiNew

                delta = abs(chiNew - chiOld)
                if (delta > maxDelta)
                    maxDelta = delta
                endif

                PhaseChi[ix][iy] = chiNew

            endfor
        endfor

        if (maxDelta < tol)
            break
        endif

    endfor

    PhaseCorr = PhaseFree + PhaseChi
    PhaseReCorr = cos(PhaseCorr)
    PhaseImCorr = sin(PhaseCorr)

    QxFree = 0
    QyFree = 0
    QmagFree = 0

    QxCorr = 0
    QyCorr = 0
    QmagCorr = 0

    Variable dXP, dXM, dYP, dYM

    // Free phase gradient.
    for (ix = 1; ix < nx-1; ix += 1)
        for (iy = 1; iy < ny-1; iy += 1)

            reC = PhaseReFree[ix][iy]
            imC = PhaseImFree[ix][iy]

            reR = PhaseReFree[ix+1][iy]
            imR = PhaseImFree[ix+1][iy]
            reL = PhaseReFree[ix-1][iy]
            imL = PhaseImFree[ix-1][iy]
            reU = PhaseReFree[ix][iy+1]
            imU = PhaseImFree[ix][iy+1]
            reD = PhaseReFree[ix][iy-1]
            imD = PhaseImFree[ix][iy-1]

            dXP = atan2(imR*reC - reR*imC, reR*reC + imR*imC)
            dXM = atan2(imC*reL - reC*imL, reC*reL + imC*imL)

            dYP = atan2(imU*reC - reU*imC, reU*reC + imU*imC)
            dYM = atan2(imC*reD - reC*imD, reC*reD + imC*imD)

            QxFree[ix][iy] = (dXP + dXM) / (2*absdx)
            QyFree[ix][iy] = (dYP + dYM) / (2*absdy)

        endfor
    endfor

    // Corrected phase gradient.
    for (ix = 1; ix < nx-1; ix += 1)
        for (iy = 1; iy < ny-1; iy += 1)

            reC = PhaseReCorr[ix][iy]
            imC = PhaseImCorr[ix][iy]

            reR = PhaseReCorr[ix+1][iy]
            imR = PhaseImCorr[ix+1][iy]
            reL = PhaseReCorr[ix-1][iy]
            imL = PhaseImCorr[ix-1][iy]
            reU = PhaseReCorr[ix][iy+1]
            imU = PhaseImCorr[ix][iy+1]
            reD = PhaseReCorr[ix][iy-1]
            imD = PhaseImCorr[ix][iy-1]

            dXP = atan2(imR*reC - reR*imC, reR*reC + imR*imC)
            dXM = atan2(imC*reL - reC*imL, reC*reL + imC*imL)

            dYP = atan2(imU*reC - reU*imC, reU*reC + imU*imC)
            dYM = atan2(imC*reD - reC*imD, reC*reD + imC*imD)

            QxCorr[ix][iy] = (dXP + dXM) / (2*absdx)
            QyCorr[ix][iy] = (dYP + dYM) / (2*absdy)

        endfor
    endfor

    QxFree[0][]    = QxFree[1][q]
    QxFree[nx-1][] = QxFree[nx-2][q]
    QxFree[][0]    = QxFree[p][1]
    QxFree[][ny-1] = QxFree[p][ny-2]

    QyFree[0][]    = QyFree[1][q]
    QyFree[nx-1][] = QyFree[nx-2][q]
    QyFree[][0]    = QyFree[p][1]
    QyFree[][ny-1] = QyFree[p][ny-2]

    QmagFree = sqrt(QxFree^2 + QyFree^2)

    QxCorr[0][]    = QxCorr[1][q]
    QxCorr[nx-1][] = QxCorr[nx-2][q]
    QxCorr[][0]    = QxCorr[p][1]
    QxCorr[][ny-1] = QxCorr[p][ny-2]

    QyCorr[0][]    = QyCorr[1][q]
    QyCorr[nx-1][] = QyCorr[nx-2][q]
    QyCorr[][0]    = QyCorr[p][1]
    QyCorr[][ny-1] = QyCorr[p][ny-2]

    QmagCorr = sqrt(QxCorr^2 + QyCorr^2)

    JxFree = PhaseStiffness * QxFree
    JyFree = PhaseStiffness * QyFree
    JmagFree = sqrt(JxFree^2 + JyFree^2)

    JxCorr = PhaseStiffness * QxCorr
    JyCorr = PhaseStiffness * QyCorr
    JmagCorr = sqrt(JxCorr^2 + JyCorr^2)

    Variable/G StiffIterUsed = iterUsed
    Variable/G StiffMaxDelta = maxDelta
    Variable/G StiffEtaN = etaN
    Variable/G StiffMaskNIsOne = maskNIsOne
    Variable/G StiffCoherenceLength_nm = coherenceLength_nm
    Variable/G StiffSmoothPasses = nSmoothPass

    SetScale/P x, x0, dx, "nm", PhaseStiffnessRaw, PhaseStiffness, PhaseFree, PhaseReFree, PhaseImFree, PhaseChi, PhaseCorr, PhaseReCorr, PhaseImCorr, QxFree, QyFree, QmagFree, QxCorr, QyCorr, QmagCorr, JxFree, JyFree, JmagFree, JxCorr, JyCorr, JmagCorr
    SetScale/P y, y0, dy, "nm", PhaseStiffnessRaw, PhaseStiffness, PhaseFree, PhaseReFree, PhaseImFree, PhaseChi, PhaseCorr, PhaseReCorr, PhaseImCorr, QxFree, QyFree, QmagFree, QxCorr, QyCorr, QmagCorr, JxFree, JyFree, JmagFree, JxCorr, JyCorr, JmagCorr

    SetDataFolder $oldDF
    return 0
End


//==============================================================================
// SNS_MakeScreenedFieldsFromStiffnessPhase
//
// Purpose:
//   Apply a radial screening envelope to the phase-gradient / current-like
//   fields produced by SNS_SolveStiffnessPhaseField(...).
//
//   This function takes the output folder of SNS_SolveStiffnessPhaseField(...)
//   as input and constructs three screened field families:
//
//   1) Screened free phase-gradient field:
//          QScreenFree = f_screen(r) * QFree
//
//   2) Screened free-phase current-like field:
//          JScreenFree = f_screen(r) * JFree
//                      = f_screen(r) * PhaseStiffness * QFree
//
//   3) Screened corrected current-like field:
//          JScreenCorr = f_screen(r) * JCorr
//                      = f_screen(r) * PhaseStiffness * QCorr
//
//   For backwards compatibility with existing ScreeningModel=3 code, this
//   function also creates:
//
//          QxScreen = JxScreenCorr
//          QyScreen = JyScreenCorr
//          QmagScreen = JmagScreenCorr
//
// Model estimates:
//   screenModel = 0:
//       f(r) = exp(-r / screenLength)
//       screenLength_nm ~ lambdaL_nm
//
//   screenModel = 1:
//       Bulk Abrikosov/London asymptotic proxy.
//
//       Naive form:
//
//           f(r) = sqrt(1 + r / screenLength) * exp(-r / screenLength)
//
//       has a radial cusp at r=0 because f depends on |r| with nonzero slope.
//       To avoid core-centered numerical artifacts, this implementation uses
//       the smoothed radius
//
//           rEff = sqrt(r^2 + coreSmooth_nm^2)
//
//       and normalizes the envelope so f(0)=1:
//
//           uEff = rEff / screenLength
//           u0   = coreSmooth_nm / screenLength
//
//           f(r) = [sqrt(1+uEff) * exp(-uEff)]
//                  / [sqrt(1+u0)   * exp(-u0)]
//
//       with
//
//           coreSmooth_nm = max(0.5*pixelSize_nm, 1e-6*screenLength_nm)
//
//       Simple estimate:
//           screenLength_nm ~ lambdaL_nm
//
//   screenModel = 2:
//       f(r) = 1 / (1 + r / screenLength)
//       screenLength_nm ~ Lambda_nm = 2*lambdaL_nm^2/dFilm_nm
//
//   screenModel = 3:
//       f(r) = 1 / (1 + r / screenLength)^2
//       screenLength_nm ~ Lambda_nm = 2*lambdaL_nm^2/dFilm_nm
//
// Inputs:
//   phaseFolder        : data folder containing output waves from
//                        SNS_SolveStiffnessPhaseField(...):
//                          QxFree, QyFree
//                          JxFree, JyFree
//                          QxCorr, QyCorr
//                          JxCorr, JyCorr
//   xV_nm, yV_nm       : vortex center [nm]
//   screenLength_nm    : phenomenological length scale controlling the radial
//                        screening envelope [nm]
//   outFolder          : output folder. If empty, writes into phaseFolder.
//
// Optional Inputs:
//   screenModel        : radial envelope model. Default = 0.
//
// Outputs:
//   ScreenEnvelope
//   QxScreenFree, QyScreenFree, QmagScreenFree
//   JxScreenFree, JyScreenFree, JmagScreenFree
//   JxScreenCorr, JyScreenCorr, JmagScreenCorr
//   QxScreen, QyScreen, QmagScreen
//   screenModelUsed
//   screenLengthUsed_nm
//============================================================================== 
Function SNS_MakeScreenedFieldsFromStiffnessPhase(phaseFolder, xV_nm, yV_nm, screenLength_nm, outFolder, [screenModel])

    String phaseFolder
    Variable xV_nm, yV_nm
    Variable screenLength_nm
    String outFolder
    Variable screenModel

    if (ParamIsDefault(screenModel))
        screenModel = 0
    endif

    screenModel = round(screenModel)

    if (screenLength_nm <= 0)
        Abort "SNS_MakeScreenedFieldsFromStiffnessPhase: screenLength_nm must be > 0."
    endif
    if ((screenModel < 0) || (screenModel > 3))
        Abort "SNS_MakeScreenedFieldsFromStiffnessPhase: screenModel must be 0, 1, 2, or 3."
    endif
    if (strlen(phaseFolder) <= 0)
        Abort "SNS_MakeScreenedFieldsFromStiffnessPhase: phaseFolder must not be empty."
    endif

    String oldDF = GetDataFolder(1)

    Wave/Z QxFree_in = $(phaseFolder + ":QxFree")
    Wave/Z QyFree_in = $(phaseFolder + ":QyFree")

    Wave/Z JxFree_in = $(phaseFolder + ":JxFree")
    Wave/Z JyFree_in = $(phaseFolder + ":JyFree")

    Wave/Z QxCorr_in = $(phaseFolder + ":QxCorr")
    Wave/Z QyCorr_in = $(phaseFolder + ":QyCorr")

    Wave/Z JxCorr_in = $(phaseFolder + ":JxCorr")
    Wave/Z JyCorr_in = $(phaseFolder + ":JyCorr")

    if (!WaveExists(QxFree_in) || !WaveExists(QyFree_in) || \
        !WaveExists(JxFree_in) || !WaveExists(JyFree_in) || \
        !WaveExists(QxCorr_in) || !WaveExists(QyCorr_in) || \
        !WaveExists(JxCorr_in) || !WaveExists(JyCorr_in))

        Abort "SNS_MakeScreenedFieldsFromStiffnessPhase: phaseFolder missing QxFree/QyFree/JxFree/JyFree/QxCorr/QyCorr/JxCorr/JyCorr."
    endif

    Variable nx = DimSize(QxFree_in, 0)
    Variable ny = DimSize(QxFree_in, 1)

    if ((DimSize(QyFree_in,0) != nx) || (DimSize(QyFree_in,1) != ny) || \
        (DimSize(JxFree_in,0) != nx) || (DimSize(JxFree_in,1) != ny) || \
        (DimSize(JyFree_in,0) != nx) || (DimSize(JyFree_in,1) != ny) || \
        (DimSize(QxCorr_in,0) != nx) || (DimSize(QxCorr_in,1) != ny) || \
        (DimSize(QyCorr_in,0) != nx) || (DimSize(QyCorr_in,1) != ny) || \
        (DimSize(JxCorr_in,0) != nx) || (DimSize(JxCorr_in,1) != ny) || \
        (DimSize(JyCorr_in,0) != nx) || (DimSize(JyCorr_in,1) != ny))

        Abort "SNS_MakeScreenedFieldsFromStiffnessPhase: input field dimensions differ."
    endif

    Variable x0 = DimOffset(QxFree_in, 0)
    Variable dx = DimDelta(QxFree_in, 0)
    Variable y0 = DimOffset(QxFree_in, 1)
    Variable dy = DimDelta(QxFree_in, 1)

    Variable pixelSize_nm = 0.5 * (abs(dx) + abs(dy))
    Variable coreSmooth_nm = max(0.5 * pixelSize_nm, 1e-6 * screenLength_nm)

    if (strlen(outFolder) > 0)
        NewDataFolder/O/S $outFolder
    else
        SetDataFolder $phaseFolder
    endif

    Make/O/D/N=(nx, ny) ScreenEnvelope

    Make/O/D/N=(nx, ny) QxScreenFree, QyScreenFree, QmagScreenFree
    Make/O/D/N=(nx, ny) JxScreenFree, JyScreenFree, JmagScreenFree
    Make/O/D/N=(nx, ny) JxScreenCorr, JyScreenCorr, JmagScreenCorr

    Make/O/D/N=(nx, ny) QxScreen, QyScreen, QmagScreen

    Variable ix, iy
    Variable xnm, ynm, r, rEff, u, uEff, u0, f, f0

    u0 = coreSmooth_nm / screenLength_nm
    f0 = sqrt(1 + u0) * exp(-u0)

    for (ix = 0; ix < nx; ix += 1)
        xnm = x0 + ix*dx

        for (iy = 0; iy < ny; iy += 1)
            ynm = y0 + iy*dy

            r = sqrt((xnm - xV_nm)^2 + (ynm - yV_nm)^2)
            u = r / screenLength_nm

            if (screenModel == 0)
                f = exp(-u)

            elseif (screenModel == 1)
                rEff = sqrt(r*r + coreSmooth_nm*coreSmooth_nm)
                uEff = rEff / screenLength_nm
                f = (sqrt(1 + uEff) * exp(-uEff)) / f0

            elseif (screenModel == 2)
                f = 1 / (1 + u)

            else
                f = 1 / ((1 + u) * (1 + u))
            endif

            ScreenEnvelope[ix][iy] = f

            QxScreenFree[ix][iy] = f * QxFree_in[ix][iy]
            QyScreenFree[ix][iy] = f * QyFree_in[ix][iy]
            QmagScreenFree[ix][iy] = sqrt(QxScreenFree[ix][iy]^2 + QyScreenFree[ix][iy]^2)

            JxScreenFree[ix][iy] = f * JxFree_in[ix][iy]
            JyScreenFree[ix][iy] = f * JyFree_in[ix][iy]
            JmagScreenFree[ix][iy] = sqrt(JxScreenFree[ix][iy]^2 + JyScreenFree[ix][iy]^2)

            JxScreenCorr[ix][iy] = f * JxCorr_in[ix][iy]
            JyScreenCorr[ix][iy] = f * JyCorr_in[ix][iy]
            JmagScreenCorr[ix][iy] = sqrt(JxScreenCorr[ix][iy]^2 + JyScreenCorr[ix][iy]^2)

            QxScreen[ix][iy] = JxScreenCorr[ix][iy]
            QyScreen[ix][iy] = JyScreenCorr[ix][iy]
            QmagScreen[ix][iy] = JmagScreenCorr[ix][iy]
        endfor
    endfor

    Variable/G screenModelUsed = screenModel
    Variable/G screenLengthUsed_nm = screenLength_nm
    Variable/G screenCoreSmoothUsed_nm = coreSmooth_nm

    SetScale/P x, x0, dx, "nm", ScreenEnvelope
    SetScale/P x, x0, dx, "nm", QxScreenFree, QyScreenFree, QmagScreenFree
    SetScale/P x, x0, dx, "nm", JxScreenFree, JyScreenFree, JmagScreenFree
    SetScale/P x, x0, dx, "nm", JxScreenCorr, JyScreenCorr, JmagScreenCorr
    SetScale/P x, x0, dx, "nm", QxScreen, QyScreen, QmagScreen

    SetScale/P y, y0, dy, "nm", ScreenEnvelope
    SetScale/P y, y0, dy, "nm", QxScreenFree, QyScreenFree, QmagScreenFree
    SetScale/P y, y0, dy, "nm", JxScreenFree, JyScreenFree, JmagScreenFree
    SetScale/P y, y0, dy, "nm", JxScreenCorr, JyScreenCorr, JmagScreenCorr
    SetScale/P y, y0, dy, "nm", QxScreen, QyScreen, QmagScreen

    SetDataFolder $oldDF
    return 0
End

//==============================================================================
// SNS_DiagnoseBoundaryNormalCurrent_FromNmask
//
// Nmask convention:
//   Nmask > 0.5  : N/hole/excluded region
//   Nmask <= 0.5 : exterior/S solved region
//
// Inputs:
//   Re0, Im0 : full free vortex complex phase
//   ReC, ImC : corrected complex phase
//
// Outputs:
//   jn_free_face
//   jn_corr_face
//   jn_ratio_face
//   nFaces_face
//==============================================================================
Function SNS_DiagnoseBoundaryNormalCurrent_FromNmask(Nmask, Re0, Im0, ReC, ImC)

    Wave Nmask
    Wave Re0, Im0, ReC, ImC

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)

    Variable dx = abs(DimDelta(Nmask, 0))
    Variable dy = abs(DimDelta(Nmask, 1))

    Make/O/D/N=(nx,ny) jn_free_face, jn_corr_face, jn_ratio_face, nFaces_face
    jn_free_face  = NaN
    jn_corr_face  = NaN
    jn_ratio_face = NaN
    nFaces_face   = 0

    Variable ix, iy
    Variable nFaces
    Variable sumFree, sumCorr
    Variable dFree, dCorr

    for (ix = 1; ix < nx-1; ix += 1)
        for (iy = 1; iy < ny-1; iy += 1)

            // central pixel must be exterior/S
            if (Nmask[ix][iy] > 0.5)
                continue
            endif

            if ((numtype(ReC[ix][iy]) != 0) || (numtype(ImC[ix][iy]) != 0))
                continue
            endif

            nFaces  = 0
            sumFree = 0
            sumCorr = 0

            // N/hole on left. Use exterior derivative one pixel away to the right.
            if ((Nmask[ix-1][iy] > 0.5) && (Nmask[ix+1][iy] <= 0.5))
                if ((numtype(ReC[ix+1][iy]) == 0) && (numtype(ImC[ix+1][iy]) == 0))

                    dFree = atan2(Im0[ix+1][iy]*Re0[ix][iy] - Re0[ix+1][iy]*Im0[ix][iy], \
                                  Re0[ix+1][iy]*Re0[ix][iy] + Im0[ix+1][iy]*Im0[ix][iy]) / dx

                    dCorr = atan2(ImC[ix+1][iy]*ReC[ix][iy] - ReC[ix+1][iy]*ImC[ix][iy], \
                                  ReC[ix+1][iy]*ReC[ix][iy] + ImC[ix+1][iy]*ImC[ix][iy]) / dx

                    sumFree += -dFree
                    sumCorr += -dCorr
                    nFaces += 1
                endif
            endif

            // N/hole on right. Use exterior derivative one pixel away to the left.
            if ((Nmask[ix+1][iy] > 0.5) && (Nmask[ix-1][iy] <= 0.5))
                if ((numtype(ReC[ix-1][iy]) == 0) && (numtype(ImC[ix-1][iy]) == 0))

                    dFree = atan2(Im0[ix][iy]*Re0[ix-1][iy] - Re0[ix][iy]*Im0[ix-1][iy], \
                                  Re0[ix][iy]*Re0[ix-1][iy] + Im0[ix][iy]*Im0[ix-1][iy]) / dx

                    dCorr = atan2(ImC[ix][iy]*ReC[ix-1][iy] - ReC[ix][iy]*ImC[ix-1][iy], \
                                  ReC[ix][iy]*ReC[ix-1][iy] + ImC[ix][iy]*ImC[ix-1][iy]) / dx

                    sumFree += dFree
                    sumCorr += dCorr
                    nFaces += 1
                endif
            endif

            // N/hole below. Use exterior derivative one pixel upward.
            if ((Nmask[ix][iy-1] > 0.5) && (Nmask[ix][iy+1] <= 0.5))
                if ((numtype(ReC[ix][iy+1]) == 0) && (numtype(ImC[ix][iy+1]) == 0))

                    dFree = atan2(Im0[ix][iy+1]*Re0[ix][iy] - Re0[ix][iy+1]*Im0[ix][iy], \
                                  Re0[ix][iy+1]*Re0[ix][iy] + Im0[ix][iy+1]*Im0[ix][iy]) / dy

                    dCorr = atan2(ImC[ix][iy+1]*ReC[ix][iy] - ReC[ix][iy+1]*ImC[ix][iy], \
                                  ReC[ix][iy+1]*ReC[ix][iy] + ImC[ix][iy+1]*ImC[ix][iy]) / dy

                    sumFree += -dFree
                    sumCorr += -dCorr
                    nFaces += 1
                endif
            endif

            // N/hole above. Use exterior derivative one pixel downward.
            if ((Nmask[ix][iy+1] > 0.5) && (Nmask[ix][iy-1] <= 0.5))
                if ((numtype(ReC[ix][iy-1]) == 0) && (numtype(ImC[ix][iy-1]) == 0))

                    dFree = atan2(Im0[ix][iy]*Re0[ix][iy-1] - Re0[ix][iy]*Im0[ix][iy-1], \
                                  Re0[ix][iy]*Re0[ix][iy-1] + Im0[ix][iy]*Im0[ix][iy-1]) / dy

                    dCorr = atan2(ImC[ix][iy]*ReC[ix][iy-1] - ReC[ix][iy]*ImC[ix][iy-1], \
                                  ReC[ix][iy]*ReC[ix][iy-1] + ImC[ix][iy]*ImC[ix][iy-1]) / dy

                    sumFree += dFree
                    sumCorr += dCorr
                    nFaces += 1
                endif
            endif

            if (nFaces > 0)
                jn_free_face[ix][iy] = sumFree / nFaces
                jn_corr_face[ix][iy] = sumCorr / nFaces
                nFaces_face[ix][iy]  = nFaces

                if (abs(jn_free_face[ix][iy]) > 1e-12)
                    jn_ratio_face[ix][iy] = abs(jn_corr_face[ix][iy]) / abs(jn_free_face[ix][iy])
                endif
            endif

        endfor
    endfor

    SetScale/P x, DimOffset(Nmask,0), DimDelta(Nmask,0), WaveUnits(Nmask,0), jn_free_face, jn_corr_face, jn_ratio_face, nFaces_face
    SetScale/P y, DimOffset(Nmask,1), DimDelta(Nmask,1), WaveUnits(Nmask,1), jn_free_face, jn_corr_face, jn_ratio_face, nFaces_face

    Print "Boundary faces found = ", sum(nFaces_face)

End


//==============================================================================
// SNS_ExtractFreeVortexPhaseOnBoundary
//
// Purpose:
//   Extract free-vortex phase (phi0_ext) restricted to the boundary and
//   output its complex representation.
//
// Inputs:
//   BoundaryMask   : 2D wave, 1 on boundary pixels, 0 elsewhere
//   phi0_ext       : 2D wave, free vortex phase [rad]
//
// Outputs (created in current data folder):
//   phaseRe0_boundary_ext : cos(phi0_ext) on boundary, NaN elsewhere
//   phaseIm0_boundary_ext : sin(phi0_ext) on boundary, NaN elsewhere
//
// Notes:
//   • Uses complex representation to avoid branch-cut issues of atan2.
//   • Keeps same scaling as input waves.
//==============================================================================
Function SNS_ExtractFreeVortexPhaseOnBoundary(BoundaryMask, phi0_ext)

    Wave BoundaryMask
    Wave phi0_ext

    Variable nx = DimSize(phi0_ext, 0)
    Variable ny = DimSize(phi0_ext, 1)

    if ((DimSize(BoundaryMask,0) != nx) || (DimSize(BoundaryMask,1) != ny))
        Abort "SNS_ExtractFreeVortexPhaseOnBoundary: size mismatch."
    endif

    Make/O/D/N=(nx, ny) phaseRe0_boundary_ext, phaseIm0_boundary_ext

    phaseRe0_boundary_ext = NaN
    phaseIm0_boundary_ext = NaN

    Variable ix, iy

    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)

            if (BoundaryMask[ix][iy] < 0.5)
                phaseRe0_boundary_ext[ix][iy] = cos(phi0_ext[ix][iy])
                phaseIm0_boundary_ext[ix][iy] = sin(phi0_ext[ix][iy])
            endif

        endfor
    endfor

    // preserve axis scaling
    SetScale/P x, DimOffset(phi0_ext,0), DimDelta(phi0_ext,0), WaveUnits(phi0_ext,0), phaseRe0_boundary_ext, phaseIm0_boundary_ext
    SetScale/P y, DimOffset(phi0_ext,1), DimDelta(phi0_ext,1), WaveUnits(phi0_ext,1), phaseRe0_boundary_ext, phaseIm0_boundary_ext

End



//==============================================================================
// SNS_SampleExteriorPhaseNearest2D and SNS_EnsureFreeVortexPhase2D
//
// Purpose
//   Sample the exterior phase field near a given coordinate by finding the
//   nearest valid exterior pixel and reconstructing the phase from
//   phaseRe_ext / phaseIm_ext.
//
//   This avoids interpolation across the hole boundary, which can lead to
//   invalid phases when one interpolation corner lies inside the hole.
//
// Inputs
//   extMask        : 2D wave, exterior domain mask
//                    extMask <= 0.5  -> exterior (valid phase)
//                    extMask > 0.5 -> hole
//
//   phaseRe_ext    : 2D wave, cos(phi_ext)
//   phaseIm_ext    : 2D wave, sin(phi_ext)
//
//   x_nm           : Sample x coordinate [nm]
//   y_nm           : Sample y coordinate [nm]
//
// Optional Inputs
//   searchRad_px   : Search radius in pixels around nearest pixel
//                    Default = 3
//
// Returns
//   phi_samp       : Sampled phase [rad]
//
// Notes
//   • Uses nearest valid exterior pixel within search radius.
//   • Phase reconstructed from cos/sin to avoid branch-cut issues.
//   • Returns NaN if no valid pixel found.
//==============================================================================
//==============================================================================
// SNS_EnsureFreeVortexPhase2D
//
// Return a folder containing PhaseFree, PhaseReFree, and PhaseImFree on the
// Nmask grid. Compatible waves in preferredFolder are reused. Otherwise the
// analytical free-vortex phase nFlux*atan2(y-yV,x-xV) is generated in
// generatedFolder without running the stiffness-phase solver.
//==============================================================================
Function/S SNS_EnsureFreeVortexPhase2D(Nmask, xV_nm, yV_nm, nFlux, preferredFolder, generatedFolder)
    Wave Nmask
    Variable xV_nm, yV_nm, nFlux
    String preferredFolder, generatedFolder

    if (WaveDims(Nmask) != 2)
        Abort "SNS_EnsureFreeVortexPhase2D: Nmask must be two-dimensional."
    endif
    if (numtype(nFlux) != 0 || round(nFlux) != nFlux)
        Abort "SNS_EnsureFreeVortexPhase2D: nFlux must be a finite integer."
    endif

    String preferredBase = preferredFolder
    if (strlen(preferredBase) > 0 && CmpStr(preferredBase[strlen(preferredBase)-1], ":") != 0)
        preferredBase += ":"
    endif

    if (strlen(preferredBase) > 0)
        Wave/Z preferredRe = $(preferredBase + "PhaseReFree")
        Wave/Z preferredIm = $(preferredBase + "PhaseImFree")
        if (WaveExists(preferredRe) && WaveExists(preferredIm))
            Variable dimensionsMatch = WaveDims(preferredRe) == 2 && WaveDims(preferredIm) == 2 && \
                DimSize(preferredRe,0) == DimSize(Nmask,0) && DimSize(preferredRe,1) == DimSize(Nmask,1) && \
                DimSize(preferredIm,0) == DimSize(Nmask,0) && DimSize(preferredIm,1) == DimSize(Nmask,1)
            Variable scalingMatches = dimensionsMatch && \
                DimOffset(preferredRe,0) == DimOffset(Nmask,0) && DimDelta(preferredRe,0) == DimDelta(Nmask,0) && \
                DimOffset(preferredRe,1) == DimOffset(Nmask,1) && DimDelta(preferredRe,1) == DimDelta(Nmask,1) && \
                DimOffset(preferredIm,0) == DimOffset(Nmask,0) && DimDelta(preferredIm,0) == DimDelta(Nmask,0) && \
                DimOffset(preferredIm,1) == DimOffset(Nmask,1) && DimDelta(preferredIm,1) == DimDelta(Nmask,1)
            if (scalingMatches)
                return preferredBase
            endif
        endif
    endif

    if (strlen(generatedFolder) == 0)
        Abort "SNS_EnsureFreeVortexPhase2D: generatedFolder must not be empty."
    endif
    String generatedBase = generatedFolder
    if (CmpStr(generatedBase[strlen(generatedBase)-1], ":") == 0)
        generatedBase = generatedBase[0, strlen(generatedBase)-2]
    endif

    String oldDF = GetDataFolder(1)
    NewDataFolder/O/S $generatedBase

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    Variable x0 = DimOffset(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dx = DimDelta(Nmask, 0)
    Variable dy = DimDelta(Nmask, 1)
    if (dx == 0 || dy == 0)
        SetDataFolder $oldDF
        Abort "SNS_EnsureFreeVortexPhase2D: Nmask has invalid spatial scaling."
    endif

    Make/O/D/N=(nx, ny) PhaseFree, PhaseReFree, PhaseImFree
    SetScale/P x, x0, dx, WaveUnits(Nmask,0), PhaseFree, PhaseReFree, PhaseImFree
    SetScale/P y, y0, dy, WaveUnits(Nmask,1), PhaseFree, PhaseReFree, PhaseImFree
    Variable ix, iy, xnm, ynm
    for (ix = 0; ix < nx; ix += 1)
        xnm = x0 + ix*dx
        for (iy = 0; iy < ny; iy += 1)
            ynm = y0 + iy*dy
            PhaseFree[ix][iy] = nFlux*atan2(ynm-yV_nm, xnm-xV_nm)
            PhaseReFree[ix][iy] = cos(PhaseFree[ix][iy])
            PhaseImFree[ix][iy] = sin(PhaseFree[ix][iy])
        endfor
    endfor
    Note/K PhaseFree, "SNS_PhaseModel=free-vortex-analytic;SNS_nFlux=" + num2str(nFlux) + ";SNS_xV_nm=" + num2str(xV_nm) + ";SNS_yV_nm=" + num2str(yV_nm) + ";"
    Note/K PhaseReFree, note(PhaseFree)
    Note/K PhaseImFree, note(PhaseFree)

    SetDataFolder $oldDF
    return generatedBase + ":"
End


Function SNS_SampleExteriorPhaseNearest2D(extMask, phaseRe_ext, phaseIm_ext, x_nm, y_nm, [searchRad_px])

    Wave extMask, phaseRe_ext, phaseIm_ext
    Variable x_nm, y_nm
    Variable searchRad_px

    if (ParamIsDefault(searchRad_px))
        searchRad_px = 3
    endif

    Variable nx = DimSize(extMask, 0)
    Variable ny = DimSize(extMask, 1)

    Variable x0 = DimOffset(extMask, 0)
    Variable dx = DimDelta(extMask, 0)

    Variable y0 = DimOffset(extMask, 1)
    Variable dy = DimDelta(extMask, 1)

    Variable ix0 = round((x_nm - x0)/dx)
    Variable iy0 = round((y_nm - y0)/dy)

    ix0 = max(0, min(nx-1, ix0))
    iy0 = max(0, min(ny-1, iy0))

    Variable ix, iy, rad
    Variable bestIx = -1
    Variable bestIy = -1
    Variable d2, d2min = Inf
    Variable xpix, ypix

    for (rad = 0; rad <= searchRad_px; rad += 1)

        for (ix = max(0, ix0-rad); ix <= min(nx-1, ix0+rad); ix += 1)
        for (iy = max(0, iy0-rad); iy <= min(ny-1, iy0+rad); iy += 1)

            if (extMask[ix][iy] > 0.5)
                continue
            endif

            xpix = x0 + ix*dx
            ypix = y0 + iy*dy

            d2 = (xpix - x_nm)^2 + (ypix - y_nm)^2

            if (d2 < d2min)
                d2min = d2
                bestIx = ix
                bestIy = iy
            endif

        endfor
        endfor

        if (bestIx >= 0)
            break
        endif

    endfor

    if (bestIx < 0)
        return NaN
    endif

    return atan2( \
        phaseIm_ext[bestIx][bestIy], \
        phaseRe_ext[bestIx][bestIy] \
    )

End

//==============================================================================
// SNS_ComputeBetaExtraFromExteriorPhase2D
//
// Purpose
//   Compute the per-channel extra phase betaExtra_List from a solved exterior
//   phase field, sampling the phase at the two N–S boundary hit points of each
//   trajectory.
//
//   This replaces the simple atan2 point-vortex phase and instead uses the
//   boundary-conditioned phase obtained from SNS_SolveExteriorPhaseField2D.
//
//   For each channel γ:
//
//       betaExtra = nFlux * wrap( phi(hit2) - phi(hit1) )
//
//   where phi is reconstructed from phaseRe_ext / phaseIm_ext using nearest
//   valid exterior pixels.
//
// Inputs
//   extMask        : 2D wave, exterior domain mask
//                    extMask > 0.5  -> exterior (valid phase)
//                    extMask <= 0.5 -> hole / invalid region
//
//   phaseRe_ext    : 2D wave, cos(phi_ext)
//   phaseIm_ext    : 2D wave, sin(phi_ext)
//
//   Hit1x_List     : 1D wave, endpoint 1 x coordinate [nm]
//   Hit1y_List     : 1D wave, endpoint 1 y coordinate [nm]
//   Hit2x_List     : 1D wave, endpoint 2 x coordinate [nm]
//   Hit2y_List     : 1D wave, endpoint 2 y coordinate [nm]
//
//   folder         : Output data folder
//
// Optional Inputs
//   nFlux          : Integer vortex winding number
//                    Default = 1
//
//   searchRad_px   : Search radius (pixels) used when sampling phase near the
//                    boundary to ensure a valid exterior pixel is used.
//                    Default = 3
//
// Outputs (in folder)
//   phi1_List      : Phase at endpoint 1 [rad]
//   phi2_List      : Phase at endpoint 2 [rad]
//   betaExtra_List : Wrapped phase difference [rad]
//
// Notes
//   • Phase is reconstructed from cos/sin to avoid branch-cut artifacts.
//   • Sampling uses nearest valid exterior pixel to avoid mixing with hole.
//   • betaExtra_List can be passed directly to the ABS solver.
//   • extMask must correspond to the same geometry used in the phase solve.
//==============================================================================
Function SNS_ComputeBetaExtraFromExteriorPhase2D(extMask, phaseRe_ext, phaseIm_ext, Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, folder, [nFlux, searchRad_px])

    Wave extMask, phaseRe_ext, phaseIm_ext
    Wave Hit1x_List, Hit1y_List
    Wave Hit2x_List, Hit2y_List
    String folder
    Variable nFlux, searchRad_px

    if (ParamIsDefault(nFlux))
        nFlux = 1
    endif

    if (ParamIsDefault(searchRad_px))
        searchRad_px = 3
    endif

    Variable Nch = numpnts(Hit1x_List)

    if ((Nch <= 0) || (Nch != numpnts(Hit1y_List)) || (Nch != numpnts(Hit2x_List)) || (Nch != numpnts(Hit2y_List)))
        Abort "SNS_ComputeBetaExtraFromExteriorPhase2D: inconsistent hit-point waves."
    endif

    String oldDF = GetDataFolder(1)
    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    Make/O/D/N=(Nch) phi1_List, phi2_List, betaExtra_List

    Variable j
    Variable dphi

    for (j = 0; j < Nch; j += 1)

        phi1_List[j] = SNS_SampleExteriorPhaseNearest2D(extMask, phaseRe_ext, phaseIm_ext, Hit1x_List[j], Hit1y_List[j], searchRad_px=searchRad_px)
        phi2_List[j] = SNS_SampleExteriorPhaseNearest2D(extMask, phaseRe_ext, phaseIm_ext, Hit2x_List[j], Hit2y_List[j], searchRad_px=searchRad_px)

        if (numtype(phi1_List[j]) || numtype(phi2_List[j]))
            betaExtra_List[j] = NaN
        else
            dphi = SNS_WrapToPi(phi2_List[j] - phi1_List[j])
            betaExtra_List[j] = dphi
        endif

    endfor

    SetDataFolder $oldDF

    return 0
End

//==============================================================================
// SNS_ComputeBetaExtraFromQField2D
//
// betaExtra_List[ch] = ∫ Q · dl along Hit1 -> Hit2
//
// Qx,Qy units: rad/nm
// Hit coordinates: nm
// Output betaExtra_List: rad
//
// Optional inputs:
//   nStep
//       Number of midpoint integration steps. Default = 256.
//
//   useInterp
//       0 = legacy nearest-neighbor sampling
//       1 = bilinear sampling using Igor built-in Interp2D()
//       Default = 1.
//
//   minUsedFrac
//       Minimum fraction of valid Q samples required for accepting a trajectory.
//       Default = 0.95.
//
//   coreHandling
//       0 = no special core handling.
//       1 = boundary-project samples inside rCore_nm around vortex center.
//       2 = reject any trajectory with samples inside rCore_nm.
//       3 = cut out core samples; no correction. Legacy sample-skip behavior.
//       4 = geometrically split trajectory at excluded core:
//              Hit1 -> core entry
//              core exit -> Hit2
//           Integrate Q only on these two outside segments, then fill the
//           removed core segment using the supplied corrected phase field:
//
//              beta_core = phase(exit) - phase(entry)
//
//           where phase = atan2(PhaseImCore, PhaseReCore), unwrapped to
//           nearest branch in [-pi, pi].
//
//   xV_nm, yV_nm
//       Vortex center [nm]. Required if coreHandling > 0.
//
//   rCore_nm
//       Core exclusion/projection radius [nm]. Required if coreHandling > 0.
//
//   vortexSign
//       Kept for backward compatibility. Not used by coreHandling=4.
//
//   maxCoreCorrectionAbs
//       Optional rejection threshold for abs(beta_core) in coreHandling=4.
//       Default = Inf.
//
//   PhaseReCore, PhaseImCore
//       Corrected phase-field components used to fill the missing core segment
//       in coreHandling=4. Required for coreHandling=4.
//
// Diagnostics written to folder:
//   betaExtra_List
//   nUsed_List
//   usedFrac_List
//   nNaN_List
//   nCoreProj_List
//   nCoreHit_List
//   nCoreCut_List
//   coreCutLength_List
//   coreVortexJump_List        // stores phase-field core correction in mode 4
//   coreSmoothCorrection_List  // always 0 here
//   coreCorrection_List
//   coreImpact_List
//   badBeta_List
//
// badBeta_List codes:
//   0 = accepted
//   1 = insufficient valid samples
//   2 = rejected by coreHandling=2
//   3 = invalid hit-point coordinates
//   4 = abs(coreCorrection) exceeded maxCoreCorrectionAbs
//   5 = failed to evaluate PhaseReCore/PhaseImCore at core boundary
//
// Notes:
//   • Modes 0–3 preserve previous behavior.
//   • Mode 4 no longer uses analytic vortex jumps.
//   • Mode 4 uses Q outside the excluded disk and corrected phase difference
//     across the disk, matching the logic of SNS_ComputeBetaExtraFromExteriorPhase2D.
//==============================================================================
Function SNS_ComputeBetaExtraFromQField2D(Qx, Qy, Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, folder, [nStep, useInterp, minUsedFrac, coreHandling, xV_nm, yV_nm, rCore_nm, vortexSign, maxCoreCorrectionAbs, PhaseReCore, PhaseImCore])

    Wave Qx, Qy
    Wave Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
    String folder
    Variable nStep, useInterp, minUsedFrac, coreHandling
    Variable xV_nm, yV_nm, rCore_nm
    Variable vortexSign, maxCoreCorrectionAbs
    Wave PhaseReCore, PhaseImCore

    if (ParamIsDefault(nStep))
        nStep = 256
    endif
    nStep = round(nStep)
    if (nStep < 1)
        Abort "SNS_ComputeBetaExtraFromQField2D: nStep must be >= 1."
    endif

    if (ParamIsDefault(useInterp))
        useInterp = 1
    endif
    useInterp = round(useInterp)
    if ((useInterp < 0) || (useInterp > 1))
        Abort "SNS_ComputeBetaExtraFromQField2D: useInterp must be 0 or 1."
    endif

    if (ParamIsDefault(minUsedFrac))
        minUsedFrac = 0.95
    endif
    if ((minUsedFrac <= 0) || (minUsedFrac > 1))
        Abort "SNS_ComputeBetaExtraFromQField2D: minUsedFrac must be in (0,1]."
    endif

    if (ParamIsDefault(coreHandling))
        coreHandling = 0
    endif
    coreHandling = round(coreHandling)
    if ((coreHandling < 0) || (coreHandling > 4))
        Abort "SNS_ComputeBetaExtraFromQField2D: coreHandling must be 0, 1, 2, 3, or 4."
    endif

    if (coreHandling != 0)
        if (ParamIsDefault(xV_nm) || ParamIsDefault(yV_nm) || ParamIsDefault(rCore_nm))
            Abort "SNS_ComputeBetaExtraFromQField2D: coreHandling>0 requires xV_nm, yV_nm, and rCore_nm."
        endif
        if (!(rCore_nm > 0))
            Abort "SNS_ComputeBetaExtraFromQField2D: rCore_nm must be > 0."
        endif
    endif

    if (coreHandling == 4)
        if (ParamIsDefault(PhaseReCore) || ParamIsDefault(PhaseImCore))
            Abort "SNS_ComputeBetaExtraFromQField2D: coreHandling=4 requires PhaseReCore and PhaseImCore."
        endif
        if ((DimSize(PhaseReCore,0) != DimSize(Qx,0)) || (DimSize(PhaseReCore,1) != DimSize(Qx,1)) || \
            (DimSize(PhaseImCore,0) != DimSize(Qx,0)) || (DimSize(PhaseImCore,1) != DimSize(Qx,1)))
            Abort "SNS_ComputeBetaExtraFromQField2D: PhaseReCore/PhaseImCore dimensions must match Qx/Qy."
        endif
    endif

    if (ParamIsDefault(vortexSign))
        vortexSign = 1
    endif
    if (vortexSign >= 0)
        vortexSign = 1
    else
        vortexSign = -1
    endif

    if (ParamIsDefault(maxCoreCorrectionAbs))
        maxCoreCorrectionAbs = Inf
    endif

    Variable Nch = numpnts(Hit1x_List)
    if ((Nch <= 0) || (Nch != numpnts(Hit1y_List)) || \
        (Nch != numpnts(Hit2x_List)) || (Nch != numpnts(Hit2y_List)))
        Abort "SNS_ComputeBetaExtraFromQField2D: inconsistent hit-point waves."
    endif

    Variable nx = DimSize(Qx, 0)
    Variable ny = DimSize(Qx, 1)

    if ((DimSize(Qy, 0) != nx) || (DimSize(Qy, 1) != ny))
        Abort "SNS_ComputeBetaExtraFromQField2D: Qx/Qy dimensions differ."
    endif

    Variable x0 = DimOffset(Qx, 0)
    Variable dx = DimDelta(Qx, 0)
    Variable y0 = DimOffset(Qx, 1)
    Variable dy = DimDelta(Qx, 1)

    if ((dx == 0) || (dy == 0))
        Abort "SNS_ComputeBetaExtraFromQField2D: Qx wave scaling has zero step."
    endif

    String oldDF = GetDataFolder(1)
    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    Make/O/D/N=(Nch) betaExtra_List
    Make/O/D/N=(Nch) nUsed_List, usedFrac_List, nNaN_List
    Make/O/D/N=(Nch) nCoreProj_List, nCoreHit_List, nCoreCut_List
    Make/O/D/N=(Nch) coreCutLength_List, coreVortexJump_List
    Make/O/D/N=(Nch) coreSmoothCorrection_List, coreCorrection_List
    Make/O/D/N=(Nch) coreImpact_List, badBeta_List

    betaExtra_List            = NaN
    nUsed_List                = 0
    usedFrac_List             = NaN
    nNaN_List                 = 0
    nCoreProj_List            = 0
    nCoreHit_List             = 0
    nCoreCut_List             = 0
    coreCutLength_List        = 0
    coreVortexJump_List       = 0
    coreSmoothCorrection_List = 0
    coreCorrection_List       = 0
    coreImpact_List           = NaN
    badBeta_List              = 0

    Variable ch, k
    Variable x1, y1, x2, y2
    Variable dxTot, dyTot, lenTot
    Variable dlx, dly, t
    Variable x, y, xs, ys
    Variable ix, iy
    Variable qxVal, qyVal
    Variable beta, nUsed, nNaN, nCoreProj, nCoreHit
    Variable rCore2, tiny
    Variable rx, ry, r2, r
    Variable rejectCore
    Variable denomUsed, usedFrac

    Variable fx, fy
    Variable aa, bb, cc, disc, sqDisc
    Variable tA, tB, tEntry, tExit, tLo, tHi
    Variable hasCoreCut
    Variable xEntry, yEntry, xExit, yExit
    Variable cutDx, cutDy, cutLen
    Variable impact
    Variable betaCore
    Variable seg, segLen, nStepSeg, kSeg
    Variable xa, ya, xb, yb
    Variable dlxSeg, dlySeg, tSeg
    Variable nSegTotal
    Variable reEntry, imEntry, reExit, imExit
    Variable phEntry, phExit

    tiny = 1e-12
    rCore2 = 0
    if (coreHandling != 0)
        rCore2 = rCore_nm*rCore_nm
    endif

    for (ch = 0; ch < Nch; ch += 1)

        x1 = Hit1x_List[ch]
        y1 = Hit1y_List[ch]
        x2 = Hit2x_List[ch]
        y2 = Hit2y_List[ch]

        if (numtype(x1) || numtype(y1) || numtype(x2) || numtype(y2))
            badBeta_List[ch] = 3
            continue
        endif

        dxTot = x2 - x1
        dyTot = y2 - y1
        lenTot = sqrt(dxTot*dxTot + dyTot*dyTot)

        if (!(lenTot > tiny))
            badBeta_List[ch] = 3
            continue
        endif

        dlx = dxTot / nStep
        dly = dyTot / nStep

        beta = 0
        nUsed = 0
        nNaN = 0
        nCoreProj = 0
        nCoreHit = 0
        rejectCore = 0

        hasCoreCut = 0
        tEntry = NaN
        tExit = NaN
        betaCore = 0
        cutLen = 0
        impact = NaN

        // ---------- determine straight-line intersection with core disk ----------
        if (coreHandling != 0)

            fx = x1 - xV_nm
            fy = y1 - yV_nm

            aa = dxTot*dxTot + dyTot*dyTot
            bb = 2*(fx*dxTot + fy*dyTot)
            cc = fx*fx + fy*fy - rCore2

            impact = abs(fx*dyTot - fy*dxTot) / lenTot
            coreImpact_List[ch] = impact

            disc = bb*bb - 4*aa*cc

            if (disc >= 0)
                sqDisc = sqrt(disc)
                tA = (-bb - sqDisc) / (2*aa)
                tB = (-bb + sqDisc) / (2*aa)

                tLo = max(tA, 0)
                tHi = min(tB, 1)

                if ((tHi > tLo) && (tHi >= 0) && (tLo <= 1))
                    hasCoreCut = 1
                    tEntry = tLo
                    tExit  = tHi

                    xEntry = x1 + tEntry*dxTot
                    yEntry = y1 + tEntry*dyTot
                    xExit  = x1 + tExit*dxTot
                    yExit  = y1 + tExit*dyTot

                    cutDx = xExit - xEntry
                    cutDy = yExit - yEntry
                    cutLen = sqrt(cutDx*cutDx + cutDy*cutDy)

                    nCoreCut_List[ch] = 1
                    coreCutLength_List[ch] = cutLen
                    nCoreHit = max(1, round(nStep*(tExit - tEntry)))
                endif
            endif
        endif

        // ---------- coreHandling=4: split trajectory + fill cut using corrected phase ----------
        if ((coreHandling == 4) && hasCoreCut)

            if (useInterp)
                reEntry = Interp2D(PhaseReCore, xEntry, yEntry)
                imEntry = Interp2D(PhaseImCore, xEntry, yEntry)
                reExit  = Interp2D(PhaseReCore, xExit, yExit)
                imExit  = Interp2D(PhaseImCore, xExit, yExit)
            else
                ix = round((xEntry - x0) / dx)
                iy = round((yEntry - y0) / dy)

                if ((ix < 0) || (ix >= nx) || (iy < 0) || (iy >= ny))
                    betaExtra_List[ch] = NaN
                    badBeta_List[ch] = 5
                    continue
                endif

                reEntry = PhaseReCore[ix][iy]
                imEntry = PhaseImCore[ix][iy]

                ix = round((xExit - x0) / dx)
                iy = round((yExit - y0) / dy)

                if ((ix < 0) || (ix >= nx) || (iy < 0) || (iy >= ny))
                    betaExtra_List[ch] = NaN
                    badBeta_List[ch] = 5
                    continue
                endif

                reExit = PhaseReCore[ix][iy]
                imExit = PhaseImCore[ix][iy]
            endif

            if (numtype(reEntry) || numtype(imEntry) || numtype(reExit) || numtype(imExit))
                betaExtra_List[ch] = NaN
                badBeta_List[ch] = 5
                continue
            endif

            phEntry = atan2(imEntry, reEntry)
            phExit  = atan2(imExit,  reExit)

            betaCore = phExit - phEntry

            if (betaCore > pi)
                betaCore -= 2*pi
            elseif (betaCore < -pi)
                betaCore += 2*pi
            endif

            coreVortexJump_List[ch]       = betaCore
            coreSmoothCorrection_List[ch] = 0
            coreCorrection_List[ch]       = betaCore

            if (abs(betaCore) > maxCoreCorrectionAbs)
                betaExtra_List[ch] = NaN
                badBeta_List[ch] = 4
                continue
            endif

            nSegTotal = 0

            // Segment 0: Hit1 -> entry
            // Segment 1: exit -> Hit2
            for (seg = 0; seg < 2; seg += 1)

                if (seg == 0)
                    xa = x1
                    ya = y1
                    xb = xEntry
                    yb = yEntry
                else
                    xa = xExit
                    ya = yExit
                    xb = x2
                    yb = y2
                endif

                segLen = sqrt((xb-xa)^2 + (yb-ya)^2)

                if (!(segLen > tiny))
                    continue
                endif

                nStepSeg = max(1, ceil(nStep * segLen / lenTot))
                nSegTotal += nStepSeg

                dlxSeg = (xb - xa) / nStepSeg
                dlySeg = (yb - ya) / nStepSeg

                for (kSeg = 0; kSeg < nStepSeg; kSeg += 1)

                    tSeg = (kSeg + 0.5) / nStepSeg
                    xs = xa + tSeg*(xb - xa)
                    ys = ya + tSeg*(yb - ya)

                    if (useInterp)
                        qxVal = Interp2D(Qx, xs, ys)
                        qyVal = Interp2D(Qy, xs, ys)
                    else
                        ix = round((xs - x0) / dx)
                        iy = round((ys - y0) / dy)

                        if ((ix < 0) || (ix >= nx) || (iy < 0) || (iy >= ny))
                            nNaN += 1
                            continue
                        endif

                        qxVal = Qx[ix][iy]
                        qyVal = Qy[ix][iy]
                    endif

                    if (numtype(qxVal) || numtype(qyVal))
                        nNaN += 1
                        continue
                    endif

                    beta += qxVal*dlxSeg + qyVal*dlySeg
                    nUsed += 1

                endfor
            endfor

            nUsed_List[ch]     = nUsed
            nNaN_List[ch]      = nNaN
            nCoreProj_List[ch] = 0
            nCoreHit_List[ch]  = nCoreHit

            if (nSegTotal < 1)
                betaExtra_List[ch] = NaN
                usedFrac_List[ch] = 0
                badBeta_List[ch] = 1
                continue
            endif

            usedFrac = nUsed / nSegTotal
            usedFrac_List[ch] = usedFrac

            if (usedFrac < minUsedFrac)
                betaExtra_List[ch] = NaN
                badBeta_List[ch] = 1
                continue
            endif

            betaExtra_List[ch] = beta + betaCore
            badBeta_List[ch] = 0
            continue

        endif

        // ---------- modes 0-3, and mode 4 with no core intersection ----------
        for (k = 0; k < nStep; k += 1)

            t = (k + 0.5) / nStep
            x = x1 + t*dxTot
            y = y1 + t*dyTot

            xs = x
            ys = y

            if (coreHandling != 0)

                rx = x - xV_nm
                ry = y - yV_nm
                r2 = rx*rx + ry*ry

                if (r2 < rCore2)

                    nCoreHit += 1

                    if (coreHandling == 2)
                        rejectCore = 1
                        continue
                    endif

                    if (coreHandling == 3)
                        continue
                    endif

                    if (coreHandling == 1)

                        r = sqrt(r2)

                        if (r > tiny)
                            xs = xV_nm + rCore_nm*rx/r
                            ys = yV_nm + rCore_nm*ry/r
                        else
                            r = sqrt(dlx*dlx + dly*dly)
                            if (r > tiny)
                                xs = xV_nm + rCore_nm*dlx/r
                                ys = yV_nm + rCore_nm*dly/r
                            else
                                nNaN += 1
                                continue
                            endif
                        endif

                        nCoreProj += 1
                    endif
                endif
            endif

            if (useInterp)
                qxVal = Interp2D(Qx, xs, ys)
                qyVal = Interp2D(Qy, xs, ys)
            else
                ix = round((xs - x0) / dx)
                iy = round((ys - y0) / dy)

                if ((ix < 0) || (ix >= nx) || (iy < 0) || (iy >= ny))
                    nNaN += 1
                    continue
                endif

                qxVal = Qx[ix][iy]
                qyVal = Qy[ix][iy]
            endif

            if (numtype(qxVal) || numtype(qyVal))
                nNaN += 1
                continue
            endif

            beta += qxVal*dlx + qyVal*dly
            nUsed += 1

        endfor

        nUsed_List[ch]     = nUsed
        nNaN_List[ch]      = nNaN
        nCoreProj_List[ch] = nCoreProj
        nCoreHit_List[ch]  = nCoreHit

        if (rejectCore)
            betaExtra_List[ch] = NaN
            usedFrac_List[ch] = nUsed / nStep
            badBeta_List[ch] = 2
            continue
        endif

        denomUsed = nStep
        if ((coreHandling == 3) && (nCoreHit > 0))
            denomUsed = nStep - nCoreHit
        endif

        if (denomUsed < 1)
            betaExtra_List[ch] = NaN
            usedFrac_List[ch] = 0
            badBeta_List[ch] = 1
            continue
        endif

        usedFrac = nUsed / denomUsed
        usedFrac_List[ch] = usedFrac

        if (usedFrac < minUsedFrac)
            betaExtra_List[ch] = NaN
            badBeta_List[ch] = 1
            continue
        endif

        betaExtra_List[ch] = beta
        badBeta_List[ch] = 0

    endfor

    SetDataFolder $oldDF
    return Nch
End
//==============================================================================
// SNS_ComputeBoundaryCurvatureFromMask2D
//
// Purpose
//   Compute a local boundary-curvature map from a 2D binary N-region mask.
//
//   The method treats the mask as a smoothed level-set field M(x,y), computes
//   the local unit normal
//
//       n = grad(M) / |grad(M)|
//
//   and then evaluates the boundary curvature as
//
//       kappa = div(n)
//
//   on the mask grid.
//
//   The output curvature is meaningful near the N-boundary and is intended to
//   be sampled at boundary hit points such as Hit1x_List / Hit2x_List.
//
// Inputs
//   Nmask          : 2D wave, N-region mask in sample coordinates
//                    1 = inside N, 0 = outside; x/y axes must be scaled in nm.
//
//   folder         : Output data folder.
//
// Optional Inputs
//   nSmoothIter    : Number of 3×3 box-smoothing iterations applied to Nmask
//                    before computing gradients. Default = 2.
//
//   gradFloor      : Floor for |grad(M)| used when normalizing.
//                    Default = 1e-6.
//
// Outputs (in folder)
//   MaskSmooth     : Smoothed mask field
//   BoundaryMask   : Boundary-pixel mask (1 on boundary, 0 elsewhere)
//   NxMap          : x-component of local unit normal
//   NyMap          : y-component of local unit normal
//   KappaMap       : Curvature map [1/nm]
//   KappaBoundary  : Curvature restricted to BoundaryMask; NaN away from boundary
//
// Notes
//   • Positive / negative sign of curvature depends on mask convention and
//     normal orientation. For most uses, abs(KappaBoundary) is the relevant
//     kink/curvature strength.
//   • Because Nmask is pixelated, smoothing is important; nSmoothIter = 1–3 is
//     a sensible range.
//   • This helper computes curvature on the grid; it does not explicitly trace
//     and parameterize the boundary contour.
//==============================================================================
Function SNS_ComputeBoundaryCurvatureFromMask2D(Nmask, folder, [nSmoothIter, gradFloor])

    Wave Nmask
    String folder
    Variable nSmoothIter, gradFloor

    if (ParamIsDefault(nSmoothIter))
        nSmoothIter = 2
    endif
    if (ParamIsDefault(gradFloor))
        gradFloor = 1e-6
    endif

    if (nSmoothIter < 0)
        Abort "SNS_ComputeBoundaryCurvatureFromMask2D: nSmoothIter must be >= 0."
    endif
    if (gradFloor <= 0)
        Abort "SNS_ComputeBoundaryCurvatureFromMask2D: gradFloor must be > 0."
    endif

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)

    if ((nx <= 2) || (ny <= 2))
        Abort "SNS_ComputeBoundaryCurvatureFromMask2D: invalid mask dimensions."
    endif

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    if ((dx == 0) || (dy == 0))
        Abort "SNS_ComputeBoundaryCurvatureFromMask2D: invalid mask scaling."
    endif

    Variable absdx = abs(dx)
    Variable absdy = abs(dy)

    String oldDF = GetDataFolder(1)
    if (strlen(folder) > 0)
        NewDataFolder/O/S $folder
    endif

    Make/O/D/N=(nx, ny) MaskSmooth, BoundaryMask, NxMap, NyMap, KappaMap, KappaBoundary
    Duplicate/FREE Nmask, tmpA
    Duplicate/FREE Nmask, tmpB

    Variable ix, iy, it
    Variable ixm, ixp, iym, iyp
    Variable gx, gy, gmag
    Variable dnx_dx, dny_dy

    // ------------------------------------------------------------
    // 1) Smooth binary mask with repeated 3x3 box averaging
    // ------------------------------------------------------------
    tmpA = Nmask

    for (it = 0; it < nSmoothIter; it += 1)

        // copy borders unchanged
        tmpB = tmpA

        for (ix = 1; ix < nx-1; ix += 1)
            for (iy = 1; iy < ny-1; iy += 1)

                tmpB[ix][iy] = ( \
                    tmpA[ix-1][iy-1] + tmpA[ix][iy-1] + tmpA[ix+1][iy-1] + \
                    tmpA[ix-1][iy]   + tmpA[ix][iy]   + tmpA[ix+1][iy]   + \
                    tmpA[ix-1][iy+1] + tmpA[ix][iy+1] + tmpA[ix+1][iy+1] ) / 9

            endfor
        endfor

        tmpA = tmpB
    endfor

    MaskSmooth = tmpA

    // ------------------------------------------------------------
    // 2) Boundary mask from original binary Nmask
    // ------------------------------------------------------------
    BoundaryMask = 0

    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)

            if (Nmask[ix][iy] <= 0.5)
                continue
            endif

            if ((ix == 0) || (ix == nx-1) || (iy == 0) || (iy == ny-1))
                BoundaryMask[ix][iy] = 1
            else
                if ((Nmask[ix-1][iy] <= 0.5) || (Nmask[ix+1][iy] <= 0.5) || \
                    (Nmask[ix][iy-1] <= 0.5) || (Nmask[ix][iy+1] <= 0.5))
                    BoundaryMask[ix][iy] = 1
                endif
            endif

        endfor
    endfor

    // ------------------------------------------------------------
    // 3) Unit normal field from smoothed mask
    // ------------------------------------------------------------
    NxMap = 0
    NyMap = 0

    for (ix = 1; ix < nx-1; ix += 1)
        for (iy = 1; iy < ny-1; iy += 1)

            gx = (MaskSmooth[ix+1][iy] - MaskSmooth[ix-1][iy]) / (2*absdx)
            gy = (MaskSmooth[ix][iy+1] - MaskSmooth[ix][iy-1]) / (2*absdy)

            gmag = sqrt(gx*gx + gy*gy)
            if (gmag < gradFloor)
                gmag = gradFloor
            endif

            NxMap[ix][iy] = gx / gmag
            NyMap[ix][iy] = gy / gmag

        endfor
    endfor

    // copy nearest interior values to borders
    for (ix = 0; ix < nx; ix += 1)
        NxMap[ix][0]     = NxMap[ix][1]
        NyMap[ix][0]     = NyMap[ix][1]
        NxMap[ix][ny-1]  = NxMap[ix][ny-2]
        NyMap[ix][ny-1]  = NyMap[ix][ny-2]
    endfor
    for (iy = 0; iy < ny; iy += 1)
        NxMap[0][iy]     = NxMap[1][iy]
        NyMap[0][iy]     = NyMap[1][iy]
        NxMap[nx-1][iy]  = NxMap[nx-2][iy]
        NyMap[nx-1][iy]  = NyMap[nx-2][iy]
    endfor

    // ------------------------------------------------------------
    // 4) Curvature map: kappa = div(n)
    // ------------------------------------------------------------
    KappaMap = NaN

    for (ix = 1; ix < nx-1; ix += 1)
        for (iy = 1; iy < ny-1; iy += 1)

            dnx_dx = (NxMap[ix+1][iy] - NxMap[ix-1][iy]) / (2*absdx)
            dny_dy = (NyMap[ix][iy+1] - NyMap[ix][iy-1]) / (2*absdy)

            KappaMap[ix][iy] = dnx_dx + dny_dy

        endfor
    endfor

    // ------------------------------------------------------------
    // 5) Boundary-restricted curvature output
    // ------------------------------------------------------------
    KappaBoundary = NaN
    for (ix = 0; ix < nx; ix += 1)
        for (iy = 0; iy < ny; iy += 1)
            if (BoundaryMask[ix][iy] > 0.5)
                KappaBoundary[ix][iy] = KappaMap[ix][iy]
            endif
        endfor
    endfor

    SetScale/P x, x0, dx, "nm", MaskSmooth, BoundaryMask, NxMap, NyMap, KappaMap, KappaBoundary
    SetScale/P y, y0, dy, "nm", MaskSmooth, BoundaryMask, NxMap, NyMap, KappaMap, KappaBoundary

    SetDataFolder $oldDF
    return 0
End






//==============================================================================
// SNS_Cut_EnsureTrailingColon
//==============================================================================
Function/S SNS_Cut_EnsureTrailingColon(dfPath)
    String dfPath

    if (strlen(dfPath) == 0)
        return dfPath
    endif

    if (cmpstr(dfPath[strlen(dfPath)-1, strlen(dfPath)-1], ":") != 0)
        dfPath += ":"
    endif

    return dfPath
End


//==============================================================================
// SNS_Cut_PointInsideMask2D
//
// Return 1 if (x_nm,y_nm) lies inside mask > 0.5, otherwise 0.
//==============================================================================
Function SNS_Cut_PointInsideMask2D(mask, x_nm, y_nm)
    Wave mask
    Variable x_nm, y_nm

    if (WaveDims(mask) != 2)
        return 0
    endif

    Variable nx = DimSize(mask, 0)
    Variable ny = DimSize(mask, 1)

    Variable x0 = DimOffset(mask, 0)
    Variable dx = DimDelta(mask, 0)
    Variable y0 = DimOffset(mask, 1)
    Variable dy = DimDelta(mask, 1)

    if (dx == 0 || dy == 0)
        return 0
    endif

    Variable ix = round((x_nm - x0)/dx)
    Variable iy = round((y_nm - y0)/dy)

    if (ix < 0 || ix >= nx || iy < 0 || iy >= ny)
        return 0
    endif

    return (mask[ix][iy] > 0.5)
End







//==============================================================================
// SNS_MakeLocalLondonQFromStiffnessPhase
//
// Purpose:
//   Construct a local-London gauge-invariant vortex field from the output of
//   SNS_SolveStiffnessPhaseField(...).
//
//   Starting point from phase solver:
//
//       QCorr = grad(PhaseCorr)
//
//   This function solves, component-wise, for a dimensionless vector potential
//   ADim = (2e/hbar) A in units rad/nm:
//
//       laplacian(ADim) - [K(r)/lambdaV^2] ADim
//           = - [K(r)/lambdaV^2] QCorr
//
//   Then constructs:
//
//       QLondonCorr = QCorr - ADim
//
//   and current-like field:
//
//       JLondonCorr = K(r) * QLondonCorr
//
//   This is a local London approximation. It is physically more meaningful
//   than multiplying QCorr by a radial envelope, but it is still NOT a full
//   Pearl/nonlocal Maxwell-London solve.
//
// Inputs:
//   phaseFolder      : folder containing PhaseStiffness, QxCorr, QyCorr
//   lambdaV_nm       : vortex London screening length [nm]
//   outFolder        : output folder. If empty, writes into phaseFolder.
//
// Optional:
//   nIterMax         : max SOR iterations. Default 5000.
//   tol              : convergence tolerance. Default 1e-8.
//   omega            : SOR parameter. Default 1.4.
//
// Outputs:
//   AxDimLondon, AyDimLondon          [rad/nm]
//   QxLondonCorr, QyLondonCorr        [rad/nm]
//   QmagLondonCorr                    [rad/nm]
//   JxLondonCorr, JyLondonCorr
//   JmagLondonCorr
//   LondonIterUsed
//   LondonMaxDelta
//   LondonLambdaV_nm
//
// Boundary condition:
//   Dirichlet ADim = QCorr on the outer grid boundary,
//   so QLondonCorr -> 0 at the computational boundary.
//   This is reasonable only if the box boundary is far enough from the vortex.
//==============================================================================
Function SNS_MakeLocalLondonQFromStiffnessPhase(phaseFolder, lambdaV_nm, outFolder, [nIterMax, tol, omega])
    String phaseFolder
    Variable lambdaV_nm
    String outFolder
    Variable nIterMax, tol, omega

    if (ParamIsDefault(nIterMax))
        nIterMax = 5000
    endif
    if (ParamIsDefault(tol))
        tol = 1e-8
    endif
    if (ParamIsDefault(omega))
        omega = 1.4
    endif

    if (strlen(phaseFolder) <= 0)
        Abort "SNS_MakeLocalLondonQFromStiffnessPhase: phaseFolder must not be empty."
    endif
    if (lambdaV_nm <= 0)
        Abort "SNS_MakeLocalLondonQFromStiffnessPhase: lambdaV_nm must be > 0."
    endif
    if ((omega <= 0) || (omega >= 2))
        Abort "SNS_MakeLocalLondonQFromStiffnessPhase: omega must satisfy 0 < omega < 2."
    endif

    String oldDF = GetDataFolder(1)

    Wave/Z K_in  = $(phaseFolder + ":PhaseStiffness")
    Wave/Z Qx_in = $(phaseFolder + ":QxCorr")
    Wave/Z Qy_in = $(phaseFolder + ":QyCorr")

    if (!WaveExists(K_in) || !WaveExists(Qx_in) || !WaveExists(Qy_in))
        Abort "SNS_MakeLocalLondonQFromStiffnessPhase: missing PhaseStiffness/QxCorr/QyCorr."
    endif

    Variable nx = DimSize(Qx_in, 0)
    Variable ny = DimSize(Qx_in, 1)

    if ((DimSize(Qy_in,0) != nx) || (DimSize(Qy_in,1) != ny) || \
        (DimSize(K_in,0)  != nx) || (DimSize(K_in,1)  != ny))
        Abort "SNS_MakeLocalLondonQFromStiffnessPhase: input dimensions differ."
    endif

    Variable x0 = DimOffset(Qx_in, 0)
    Variable dx = DimDelta(Qx_in, 0)
    Variable y0 = DimOffset(Qx_in, 1)
    Variable dy = DimDelta(Qx_in, 1)

    Variable absdx = abs(dx)
    Variable absdy = abs(dy)

    if ((nx <= 2) || (ny <= 2) || (absdx <= 0) || (absdy <= 0))
        Abort "SNS_MakeLocalLondonQFromStiffnessPhase: invalid grid."
    endif

    if (strlen(outFolder) > 0)
        NewDataFolder/O/S $outFolder
    else
        SetDataFolder $phaseFolder
    endif

    Make/O/D/N=(nx, ny) AxDimLondon, AyDimLondon
    Make/O/D/N=(nx, ny) QxLondonCorr, QyLondonCorr, QmagLondonCorr
    Make/O/D/N=(nx, ny) JxLondonCorr, JyLondonCorr, JmagLondonCorr

    // Initial guess: partial screening.
    AxDimLondon = 0
    AyDimLondon = 0

    // Boundary condition: ADim = QCorr at outer boundary.
    AxDimLondon[0][]    = Qx_in[0][q]
    AxDimLondon[nx-1][] = Qx_in[nx-1][q]
    AxDimLondon[][0]    = Qx_in[p][0]
    AxDimLondon[][ny-1] = Qx_in[p][ny-1]

    AyDimLondon[0][]    = Qy_in[0][q]
    AyDimLondon[nx-1][] = Qy_in[nx-1][q]
    AyDimLondon[][0]    = Qy_in[p][0]
    AyDimLondon[][ny-1] = Qy_in[p][ny-1]

    Variable invDx2 = 1/(absdx*absdx)
    Variable invDy2 = 1/(absdy*absdy)
    Variable invLam2 = 1/(lambdaV_nm*lambdaV_nm)

    Variable iter, iterUsed
    Variable ix, iy
    Variable maxDelta, delta
    Variable mass, denom
    Variable oldVal, newVal

    iterUsed = 0

    for (iter = 0; iter < nIterMax; iter += 1)

        maxDelta = 0
        iterUsed = iter + 1

        // Re-impose boundary each iteration.
        AxDimLondon[0][]    = Qx_in[0][q]
        AxDimLondon[nx-1][] = Qx_in[nx-1][q]
        AxDimLondon[][0]    = Qx_in[p][0]
        AxDimLondon[][ny-1] = Qx_in[p][ny-1]

        AyDimLondon[0][]    = Qy_in[0][q]
        AyDimLondon[nx-1][] = Qy_in[nx-1][q]
        AyDimLondon[][0]    = Qy_in[p][0]
        AyDimLondon[][ny-1] = Qy_in[p][ny-1]

        for (ix = 1; ix < nx-1; ix += 1)
            for (iy = 1; iy < ny-1; iy += 1)

                mass = K_in[ix][iy] * invLam2
                if (mass < 0 || numtype(mass))
                    mass = 0
                endif

                denom = 2*invDx2 + 2*invDy2 + mass
                if (denom <= 0)
                    continue
                endif

                // --- Ax solve ---
                newVal = (invDx2*(AxDimLondon[ix+1][iy] + AxDimLondon[ix-1][iy]) \
                        + invDy2*(AxDimLondon[ix][iy+1] + AxDimLondon[ix][iy-1]) \
                        + mass*Qx_in[ix][iy]) / denom

                oldVal = AxDimLondon[ix][iy]
                newVal = (1 - omega)*oldVal + omega*newVal
                delta = abs(newVal - oldVal)
                if (delta > maxDelta)
                    maxDelta = delta
                endif
                AxDimLondon[ix][iy] = newVal

                // --- Ay solve ---
                newVal = (invDx2*(AyDimLondon[ix+1][iy] + AyDimLondon[ix-1][iy]) \
                        + invDy2*(AyDimLondon[ix][iy+1] + AyDimLondon[ix][iy-1]) \
                        + mass*Qy_in[ix][iy]) / denom

                oldVal = AyDimLondon[ix][iy]
                newVal = (1 - omega)*oldVal + omega*newVal
                delta = abs(newVal - oldVal)
                if (delta > maxDelta)
                    maxDelta = delta
                endif
                AyDimLondon[ix][iy] = newVal

            endfor
        endfor

        if (maxDelta < tol)
            break
        endif
    endfor

    QxLondonCorr = Qx_in - AxDimLondon
    QyLondonCorr = Qy_in - AyDimLondon
    QmagLondonCorr = sqrt(QxLondonCorr^2 + QyLondonCorr^2)

    JxLondonCorr = K_in * QxLondonCorr
    JyLondonCorr = K_in * QyLondonCorr
    JmagLondonCorr = sqrt(JxLondonCorr^2 + JyLondonCorr^2)

    Variable/G LondonIterUsed = iterUsed
    Variable/G LondonMaxDelta = maxDelta
    Variable/G LondonLambdaV_nm = lambdaV_nm

    SetScale/P x, x0, dx, "nm", AxDimLondon, AyDimLondon, QxLondonCorr, QyLondonCorr, QmagLondonCorr, JxLondonCorr, JyLondonCorr, JmagLondonCorr
    SetScale/P y, y0, dy, "nm", AxDimLondon, AyDimLondon, QxLondonCorr, QyLondonCorr, QmagLondonCorr, JxLondonCorr, JyLondonCorr, JmagLondonCorr

    SetDataFolder $oldDF
    return 0
End

//==============================================================================
// SNS_Usadel1D_GeometryFromRayFolder
//
// Purpose:
//   Convert the ballistic ray ensemble at one STS point into effective 1D
//   Usadel geometry parameters.
//
//   Ballistic-to-Usadel mapping:
//     Each ray/channel defines two distances from the STS point to the two
//     superconducting hit points:
//
//        d1_i = distance(STS, Hit1_i)
//        d2_i = distance(STS, Hit2_i)
//        L_i  = d1_i + d2_i
//
//     The effective 1D diffusion length is taken as the RMS chord length:
//
//        Ldiff = sqrt(<L_i^2>)
//
//     The probe position is mapped by preserving the average fractional
//     nearest-side position:
//
//        f_probe = < min(d1_i,d2_i) / L_i >
//        x_probe = f_probe * Ldiff
//
//     With optional transparency weighting:
//
//        <A_i>_w = sum_i w_i A_i / sum_i w_i,
//        w_i = max(T_eff_List[i], 0)
//
//   Physical meaning:
//     Ldiff sets the 1D Usadel diffusion length / Thouless scale.
//     x_probe sets where the local DOS is sampled on that compressed wire.
//     dNearMean is kept only as a diagnostic of the original 2D ray geometry.
//
// Inputs:
//   dfPath     : folder containing ray-tracing output.
//                Required waves:
//                  L_N_List
//                  Hit1x_List, Hit1y_List
//                  Hit2x_List, Hit2y_List
//
// Optional inputs:
//   STSx, STSy : STS coordinate in same coordinate system as Hit*_List.
//                If omitted, tmpSTSx/tmpSTSy are read from dfPath.
//   useWeights : optional. Default 0.
//                0 = equal channel weights.
//                1 = weight by max(T_eff_List,0), if T_eff_List exists.
//   doPrint    : optional. Default 0.
//                1 = print summary.
//
// Outputs in dfPath:
//   v_Usadel_Ldiff_nm      : effective 1D diffusion length = Lrms.
//   v_Usadel_xProbe_nm     : effective probe coordinate on 1D wire.
//   v_Usadel_fProbe        : fractional nearest-side probe position.
//   v_Usadel_Lmean_nm      : mean chord length.
//   v_Usadel_Lrms_nm       : RMS chord length.
//   v_Usadel_Lmin_nm       : minimum chord length.
//   v_Usadel_Lmax_nm       : maximum chord length.
//   v_Usadel_dNearMean_nm  : mean nearest-boundary distance in original geometry.
//   v_Usadel_hitToNm       : conversion factor from hit-coordinate units to nm.
//   v_Usadel_useWeights
//   v_Usadel_nChannels
//
//   w_UsadelChord_nm       : L_i = d1_i + d2_i.
//   w_UsadelD1_nm          : STS-to-hit1 distance.
//   w_UsadelD2_nm          : STS-to-hit2 distance.
//   w_UsadelDnear_nm       : min(d1_i,d2_i).
//   w_UsadelFracNear       : min(d1_i,d2_i)/L_i.
//   w_UsadelWeight         : weights used in averages.
//
// Notes:
//   L_N_List is used only to infer the nm conversion for hit coordinates.
//   The fractional probe position uses d1+d2 from endpoint geometry, because
//   it must be consistent with the actual STS position along each chord.
//
//==============================================================================

Function SNS_Usadel1D_GeometryFromRayFolder(dfPath, [STSx, STSy, useWeights, doPrint])
	String dfPath
	Variable STSx, STSy, useWeights, doPrint

	String oldDF = GetDataFolder(1)
	SetDataFolder $dfPath

	if (ParamIsDefault(useWeights))
		useWeights = 0
	endif
	if (ParamIsDefault(doPrint))
		doPrint = 0
	endif

	// Required ray-tracing waves.
	Wave/Z L_N_List_nm
	Wave/Z Hit1x_List_nm
	Wave/Z Hit1y_List_nm
	Wave/Z Hit2x_List_nm
	Wave/Z Hit2y_List_nm

	if (!WaveExists(L_N_List_nm) || !WaveExists(Hit1x_List_nm) || !WaveExists(Hit1y_List_nm) || !WaveExists(Hit2x_List_nm) || !WaveExists(Hit2y_List_nm))
		SetDataFolder $oldDF
		Abort "SNS_Usadel1D_GeometryFromRayFolder: missing L_N_List_nm or Hit*_List_nm waves."
	endif

	Variable nCh = DimSize(L_N_List_nm, 0)
	if (nCh <= 0)
		SetDataFolder $oldDF
		Abort "SNS_Usadel1D_GeometryFromRayFolder: empty channel list."
	endif

	// STS coordinate.
	if (ParamIsDefault(STSx) || ParamIsDefault(STSy))
		Wave/Z tmpSTSx
		Wave/Z tmpSTSy

		if (WaveExists(tmpSTSx) && WaveExists(tmpSTSy))
			STSx = tmpSTSx[0]
			STSy = tmpSTSy[0]
		else
			SetDataFolder $oldDF
			Abort "SNS_Usadel1D_GeometryFromRayFolder: provide STSx/STSy or store tmpSTSx/tmpSTSy in folder."
		endif
	endif

	// Build raw geometric chord from hit coordinates.
	Make/O/D/N=(nCh) w_UsadelChordRaw
	w_UsadelChordRaw = sqrt((Hit2x_List_nm - Hit1x_List_nm)^2 + (Hit2y_List_nm - Hit1y_List_nm)^2)

	// L_N_List_nm is already in nm.
	Make/O/D/N=(nCh) w_UsadelChordFromL_nm
	w_UsadelChordFromL_nm = L_N_List_nm

	// Hit coordinates are in nm.
	WaveStats/Q w_UsadelChordRaw
	Variable rawAvg = V_avg

	WaveStats/Q w_UsadelChordFromL_nm
	Variable Lavg_nm = V_avg

	Variable hitToNm = 1

	Make/O/D/N=(nCh) w_UsadelD1_nm
	Make/O/D/N=(nCh) w_UsadelD2_nm
	Make/O/D/N=(nCh) w_UsadelDnear_nm
	Make/O/D/N=(nCh) w_UsadelChord_nm
	Make/O/D/N=(nCh) w_UsadelFracNear
	Make/O/D/N=(nCh) w_UsadelWeight

	w_UsadelD1_nm = hitToNm * sqrt((STSx - Hit1x_List_nm)^2 + (STSy - Hit1y_List_nm)^2)
	w_UsadelD2_nm = hitToNm * sqrt((STSx - Hit2x_List_nm)^2 + (STSy - Hit2y_List_nm)^2)

	// Prefer actual one-sided distances for fractional position.
	w_UsadelChord_nm = w_UsadelD1_nm + w_UsadelD2_nm

	// If numerical mismatch is large, L_N_List remains a useful diagnostic,
	// but the point-position fractions require d1+d2.
	w_UsadelDnear_nm = min(w_UsadelD1_nm, w_UsadelD2_nm)
	w_UsadelFracNear = w_UsadelDnear_nm / w_UsadelChord_nm

	w_UsadelWeight = 1

	if (useWeights)
		Wave/Z T_eff_List
		if (WaveExists(T_eff_List))
			w_UsadelWeight = max(T_eff_List, 0)
		endif
	endif

	// Weighted stats.
	Variable i
	Variable wSum = 0
	Variable Lsum = 0
	Variable L2sum = 0
	Variable fsum = 0
	Variable dNearSum = 0

	Variable wi, Li, fi, di

	for (i=0; i<nCh; i+=1)

		wi = w_UsadelWeight[i]
		Li = w_UsadelChord_nm[i]
		fi = w_UsadelFracNear[i]
		di = w_UsadelDnear_nm[i]

		if (numtype(wi) == 0 && numtype(Li) == 0 && Li > 0 && numtype(fi) == 0)
			wSum += wi
			Lsum += wi * Li
			L2sum += wi * Li * Li
			fsum += wi * fi
			dNearSum += wi * di
		endif
	endfor

	if (wSum <= 0)
		SetDataFolder $oldDF
		Abort "SNS_Usadel1D_GeometryFromRayFolder: invalid weights or chord lengths."
	endif

	Variable Lmean_nm = Lsum / wSum
	Variable Lrms_nm  = sqrt(L2sum / wSum)
	Variable fProbe   = fsum / wSum
	Variable dNearMean_nm = dNearSum / wSum

	Variable xProbe_nm = fProbe * Lrms_nm

	WaveStats/Q w_UsadelChord_nm
	Variable Lmin_nm = V_min
	Variable Lmax_nm = V_max

	Variable/G v_Usadel_Ldiff_nm = Lrms_nm
	Variable/G v_Usadel_xProbe_nm = xProbe_nm
	Variable/G v_Usadel_fProbe = fProbe
	Variable/G v_Usadel_Lmean_nm = Lmean_nm
	Variable/G v_Usadel_Lrms_nm = Lrms_nm
	Variable/G v_Usadel_Lmin_nm = Lmin_nm
	Variable/G v_Usadel_Lmax_nm = Lmax_nm
	Variable/G v_Usadel_dNearMean_nm = dNearMean_nm
	Variable/G v_Usadel_hitToNm = hitToNm
	Variable/G v_Usadel_useWeights = useWeights
	Variable/G v_Usadel_nChannels = nCh

	KillWaves/Z w_UsadelChordRaw, w_UsadelChordFromL_nm

	if (doPrint)
		Print "SNS_Usadel1D_GeometryFromRayFolder:"
		Print "  nChannels      = ", nCh
		Print "  hitToNm        = ", hitToNm
		Print "  Lmean_nm       = ", Lmean_nm
		Print "  Lrms_nm        = ", Lrms_nm
		Print "  xProbe_nm      = ", xProbe_nm
		Print "  fProbe         = ", fProbe
		Print "  dNearMean_nm   = ", dNearMean_nm
	endif

	SetDataFolder $oldDF
	return 0
End

//==============================================================================
// SNS_Usadel1D_LeverArmFromRayFolder
//
// Purpose:
//   Estimate the magnetic phase lever arm Lphi_nm for the effective 1D Usadel
//   model from a ray-tracing folder.
//
//   Ballistic-to-Usadel mapping:
//     Each ballistic channel has a transverse magnetic lever arm
//
//        Lphi_i = abs(W_eff_List[i])
//
//     The effective 1D Usadel lever arm is taken from the channel statistics:
//
//        Lphi_mean = <Lphi_i>
//        Lphi_rms  = sqrt(<Lphi_i^2>)
//        Lphi_max  = max(Lphi_i)
//
//     with optional transparency weighting:
//
//        <A_i>_w = sum_i w_i A_i / sum_i w_i,
//        w_i = max(T_eff_List[i], 0)
//
//     The value used for the Usadel phase bias is:
//
//        Lphi_use = Lphi_rms   if useMax = 0
//        Lphi_use = Lphi_max   if useMax = 1
//
//     and the field-to-phase conversion is:
//
//        phi(B) = (2e/hbar) * B * lambdaL * Lphi_use
//
//   Uses abs(W_eff_List) as the per-channel phase lever arm. In the ballistic
//   ray geometry, W_eff_List is the effective S-S separation projected
//   perpendicular to B.
//
// Inputs:
//   dfPath      : ray-tracing folder containing W_eff_List.
//   useWeights  : optional. Default 0.
//                 0 = average channels equally.
//                 1 = weight channels by max(T_eff_List, 0), if T_eff_List exists.
//   useMax      : optional. Default 0.
//                 0 = use Lphi_rms as v_Usadel_Lphi_use_nm.
//                 1 = use Lphi_max as v_Usadel_Lphi_use_nm.
//   doPrint     : optional. Default 0.
//                 1 = print mean/rms/max/used Lphi and Bpi.
//
// Outputs in dfPath:
//   w_UsadelLphi_nm       : abs(W_eff_List) in nm.
//   w_UsadelLphiWeight    : per-channel weights used for averaging.
//   v_Usadel_Lphi_mean_nm : weighted/equal mean Lphi.
//   v_Usadel_Lphi_rms_nm  : weighted/equal RMS Lphi.
//   v_Usadel_Lphi_max_nm  : maximum Lphi.
//   v_Usadel_Lphi_use_nm  : selected lever arm, RMS or max.
//   v_Usadel_Bpi_T        : field where phi = pi using Lphi_use_nm.
//   v_Usadel_Lphi_useMax
//   v_Usadel_Lphi_useWeights
//
//==============================================================================

Function SNS_Usadel1D_LeverArmFromRayFolder(dfPath, [useWeights, useMax, doPrint])
	String dfPath
	Variable useWeights, useMax, doPrint

	String oldDF = GetDataFolder(1)
	SetDataFolder $dfPath

	if (ParamIsDefault(useWeights))
		useWeights = 0
	endif
	if (ParamIsDefault(useMax))
		useMax = 0
	endif
	if (ParamIsDefault(doPrint))
		doPrint = 0
	endif

	Wave/Z W_eff_List_nm
	if (!WaveExists(W_eff_List_nm))
		SetDataFolder $oldDF
		Abort "SNS_Usadel1D_LeverArmFromRayFolder: missing W_eff_List_nm."
	endif

	Variable nCh = DimSize(W_eff_List_nm, 0)
	if (nCh <= 0)
		SetDataFolder $oldDF
		Abort "SNS_Usadel1D_LeverArmFromRayFolder: empty W_eff_List_nm."
	endif

	Make/O/D/N=(nCh) w_UsadelLphi_nm = abs(W_eff_List_nm)
	Make/O/D/N=(nCh) w_UsadelLphiWeight = 1

	if (useWeights)
		Wave/Z T_eff_List
		if (WaveExists(T_eff_List))
			w_UsadelLphiWeight = max(T_eff_List, 0)
		endif
	endif

	Variable i, wi, Li
	Variable wSum = 0
	Variable Lsum = 0
	Variable L2sum = 0
	Variable Lmax = -Inf

	for (i=0; i<nCh; i+=1)
		wi = w_UsadelLphiWeight[i]
		Li = w_UsadelLphi_nm[i]

		if (numtype(wi) == 0 && numtype(Li) == 0 && wi >= 0)
			wSum += wi
			Lsum += wi * Li
			L2sum += wi * Li * Li
			Lmax = max(Lmax, Li)
		endif
	endfor

	if (wSum <= 0)
		SetDataFolder $oldDF
		Abort "SNS_Usadel1D_LeverArmFromRayFolder: invalid weights."
	endif

	Variable Lmean_nm = Lsum / wSum
	Variable Lrms_nm = sqrt(L2sum / wSum)

	Variable Luse_nm
	if (useMax)
		Luse_nm = Lmax
	else
		Luse_nm = Lrms_nm
	endif

	NVAR lambdaL_m = root:SNS_Settings:lambdaL
	Variable Bpi_T = pi * HBAR_SI / (2*q_e * lambdaL_m * (Luse_nm*1e-9))

	Variable/G v_Usadel_Lphi_mean_nm = Lmean_nm
	Variable/G v_Usadel_Lphi_rms_nm = Lrms_nm
	Variable/G v_Usadel_Lphi_max_nm = Lmax
	Variable/G v_Usadel_Lphi_use_nm = Luse_nm
	Variable/G v_Usadel_Bpi_T = Bpi_T
	Variable/G v_Usadel_Lphi_useMax = useMax
	Variable/G v_Usadel_Lphi_useWeights = useWeights

	if (doPrint)
		Print "SNS_Usadel1D_LeverArmFromRayFolder:"
		Print "  Lphi_mean_nm = ", Lmean_nm
		Print "  Lphi_rms_nm  = ", Lrms_nm
		Print "  Lphi_max_nm  = ", Lmax
		Print "  Lphi_use_nm  = ", Luse_nm
		Print "  Bpi_T        = ", Bpi_T
	endif

	SetDataFolder $oldDF
	return 0
End

//==============================================================================
// SNS_WarnLargeRayDisplay
//
// Warn before appending a very large number of ray traces to an Igor graph.
// This is intentionally display-only: it does not change channel lists,
// histograms, or LDOS calculations.  The threshold is on plotted graph traces,
// not on total source channels, so raySparsity is respected.
//==============================================================================
Function SNS_WarnLargeRayDisplay(context, nRequested, nPlotted, raySparsity)
	String context
	Variable nRequested, nPlotted, raySparsity

	Variable warnTraceCount = 2000

	if (nPlotted <= warnTraceCount)
		return 0
	endif

	String msg
	msg = context + ": requested ray display will append " + num2str(round(nPlotted)) + " graph traces"
	msg += " from " + num2str(round(nRequested)) + " selected/source rays"
	msg += " (raySparsity=" + num2str(round(raySparsity)) + ")."
	msg += " This can make Igor slow or unresponsive. Increase raySparsity or use a fan/summary display for quick inspection."

	SNS_Log(msg, level="WARN")

	DoAlert 1, msg + "\r\rContinue drawing these rays?"
	if (V_flag != 1)
		return 1
	endif

	return 0
End

//==============================================================================
// SNS_PlotAllChannelRays_Faint
//
// Purpose:
//   Plot ray-traced SNS channels as faint background trajectories on the
//   current graph.
//
//   Rays are selected by fixed sparsity, not by target number:
//
//      raySparsity = 1   : plot every ray
//      raySparsity = 10  : plot every 10th ray
//
//   The first plotted ray is always rayStartIndex. This makes the displayed
//   subset reproducible and independent of the total number of rays.
//
// Inputs:
//   Hit1x_List, Hit1y_List : first S-hit coordinates for each channel
//   Hit2x_List, Hit2y_List : second S-hit coordinates for each channel
//
// Optional:
//   raySparsity       : plot every raySparsity-th ray.
//                       Default: 1.
//   rayStartIndex     : first ray index to plot.
//                       Default: 0.
//   faintColor        : "gray", "black", "white", "red", "blue", "green",
//                       "cyan", "magenta", "yellow".
//                       Default: "gray".
//   colorTransparency : alpha value in [0,1].
//                       0 = transparent, 1 = opaque.
//                       Default: 0.20.
//   lineSize          : line thickness.
//                       Default: 0.25.
//
// Outputs:
//   return : numeric
//            0 on successful completion.
//==============================================================================
Function SNS_PlotAllChannelRays_Faint(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, [raySparsity, rayStartIndex, faintColor, colorTransparency, lineSize])
    Wave Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
    Variable raySparsity, rayStartIndex, colorTransparency, lineSize
    String faintColor

    // ---------------- Defaults ----------------
    if (ParamIsDefault(raySparsity))
        raySparsity = 1
    endif
    if (ParamIsDefault(rayStartIndex))
        rayStartIndex = 0
    endif
    if (ParamIsDefault(faintColor))
        faintColor = "gray"
    endif
    if (ParamIsDefault(colorTransparency))
        colorTransparency = 0.20
    endif
    if (ParamIsDefault(lineSize))
        lineSize = 0.25
    endif

    raySparsity = max(1, round(raySparsity))
    rayStartIndex = max(0, round(rayStartIndex))

    colorTransparency = max(0, min(1, colorTransparency))
    Variable aa = round(65535 * colorTransparency)

    // ---------------- Choose faint-ray base color ----------------
    Variable rr, gg, bb

    strswitch(LowerStr(faintColor))
        case "black":
            rr = 0
            gg = 0
            bb = 0
            break

        case "white":
            rr = 65535
            gg = 65535
            bb = 65535
            break

        case "red":
            rr = 65535
            gg = 0
            bb = 0
            break

        case "blue":
            rr = 0
            gg = 0
            bb = 65535
            break

        case "green":
            rr = 0
            gg = 45000
            bb = 0
            break

        case "cyan":
            rr = 0
            gg = 52000
            bb = 65535
            break

        case "magenta":
            rr = 65535
            gg = 0
            bb = 65535
            break

        case "yellow":
            rr = 65535
            gg = 52000
            bb = 0
            break

        case "gray":
        default:
            rr = 54000
            gg = 54000
            bb = 54000
            break
    endswitch

    // ---------------- Validate ray list ----------------
    Variable nRay = numpnts(Hit1x_List)

    if (nRay <= 0)
        return 0
    endif

    if (numpnts(Hit1y_List) != nRay || numpnts(Hit2x_List) != nRay || numpnts(Hit2y_List) != nRay)
        Abort "SNS_PlotAllChannelRays_Faint: endpoint lists have inconsistent lengths."
    endif

    if (rayStartIndex >= nRay)
        return 0
    endif

    Variable nRayToDraw = floor((nRay - rayStartIndex - 1) / raySparsity) + 1
    nRayToDraw = max(0, nRayToDraw)

    if (SNS_WarnLargeRayDisplay("SNS_PlotAllChannelRays_Faint", nRay, nRayToDraw, raySparsity))
        return 0
    endif

    // ---------------- Append faint ray traces ----------------
    Variable j
    String idxStr, xName, yName
    String traceList, lastTrace

    for (j = rayStartIndex; j < nRay; j += raySparsity)

        sprintf idxStr, "%04d", j

        xName = "tmpAllRayX_" + idxStr
        yName = "tmpAllRayY_" + idxStr

        Make/O/D/N=2 $xName, $yName
        Wave rayX = $xName
        Wave rayY = $yName

        rayX[0] = Hit1x_List[j]
        rayY[0] = Hit1y_List[j]
        rayX[1] = Hit2x_List[j]
        rayY[1] = Hit2y_List[j]

        AppendToGraph rayY vs rayX

        traceList = TraceNameList("", ";", 1)
        lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

        ModifyGraph rgb($lastTrace)=(rr, gg, bb, aa)
        ModifyGraph lsize($lastTrace)=lineSize
    endfor

    return 0
End

//==============================================================================
// SNS_MakeRotatedImageForDisplay
//
// Rotates image coordinates so that display x-axis points along xAxisAngle_deg
// in the original image coordinate system.
//
// New coordinates:
//   x' =  (x-xc) cosθ + (y-yc) sinθ
//   y' = -(x-xc) sinθ + (y-yc) cosθ
//
// The output image is sampled back from the original image using Interp2D.
//==============================================================================
Function SNS_MakeRotatedImageForDisplay(src, outName, xAxisAngle_deg)
	Wave src
	String outName
	Variable xAxisAngle_deg

	Variable nx = DimSize(src, 0)
	Variable ny = DimSize(src, 1)
	Variable dx = DimDelta(src, 0)
	Variable dy = DimDelta(src, 1)
	Variable x0 = DimOffset(src, 0)
	Variable y0 = DimOffset(src, 1)

	Variable x1 = x0 + dx*(nx-1)
	Variable y1 = y0 + dy*(ny-1)

	Variable xMin = min(x0, x1)
	Variable xMax = max(x0, x1)
	Variable yMin = min(y0, y1)
	Variable yMax = max(y0, y1)

	Variable xc = 0.5*(xMin+xMax)
	Variable yc = 0.5*(yMin+yMax)

	Variable th = xAxisAngle_deg*pi/180
	Variable c = cos(th)
	Variable s = sin(th)

	// Rotate four image corners to get output bounds.
	Make/FREE/D/N=4 xpCorner, ypCorner

	xpCorner[0] = (xMin-xc)*c + (yMin-yc)*s
	ypCorner[0] = -(xMin-xc)*s + (yMin-yc)*c

	xpCorner[1] = (xMax-xc)*c + (yMin-yc)*s
	ypCorner[1] = -(xMax-xc)*s + (yMin-yc)*c

	xpCorner[2] = (xMin-xc)*c + (yMax-yc)*s
	ypCorner[2] = -(xMin-xc)*s + (yMax-yc)*c

	xpCorner[3] = (xMax-xc)*c + (yMax-yc)*s
	ypCorner[3] = -(xMax-xc)*s + (yMax-yc)*c

	WaveStats/Q xpCorner
	Variable xpMin = V_min
	Variable xpMax = V_max

	WaveStats/Q ypCorner
	Variable ypMin = V_min
	Variable ypMax = V_max

	Variable dOut = min(abs(dx), abs(dy))
	if (dOut <= 0 || numtype(dOut) != 0)
		dOut = 1
	endif

	Variable nxOut = max(2, ceil((xpMax-xpMin)/dOut) + 1)
	Variable nyOut = max(2, ceil((ypMax-ypMin)/dOut) + 1)

	Make/O/D/N=(nxOut, nyOut) $outName
	Wave out = $outName

	SetScale/P x, xpMin, dOut, WaveUnits(src, 0), out
	SetScale/P y, ypMin, dOut, WaveUnits(src, 1), out

	Variable ix, iy, xp, yp, xx, yy

	for (ix = 0; ix < nxOut; ix += 1)
		xp = xpMin + ix*dOut

		for (iy = 0; iy < nyOut; iy += 1)
			yp = ypMin + iy*dOut

			// inverse transform
			xx = xc + xp*c - yp*s
			yy = yc + xp*s + yp*c

			if (xx < xMin || xx > xMax || yy < yMin || yy > yMax)
				out[ix][iy] = NaN
			else
				out[ix][iy] = Interp2D(src, xx, yy)
			endif
		endfor
	endfor

	return 0
End


//==============================================================================
// SNS_PlotChannelRay_Rotated
//
// Plots one trajectory in rotated display coordinates.
//==============================================================================
Function SNS_PlotChannelRay_Rotated(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, idx, xAxisAngle_deg, xCenter, yCenter)
	Wave Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
	Variable idx, xAxisAngle_deg, xCenter, yCenter

	NVAR/Z counter = root:SNS_Settings:v_RotRayCounter
	if (!NVAR_Exists(counter))
		Variable/G root:SNS_Settings:v_RotRayCounter = 0
		NVAR counter = root:SNS_Settings:v_RotRayCounter
	endif

	counter += 1

	String xName = "tmpRotRayX_" + num2istr(counter)
	String yName = "tmpRotRayY_" + num2istr(counter)

	Make/O/D/N=2 $xName, $yName
	Wave rx = $xName
	Wave ry = $yName

	Variable th = xAxisAngle_deg*pi/180
	Variable c = cos(th)
	Variable s = sin(th)

	Variable x, y

	x = Hit1x_List[idx]
	y = Hit1y_List[idx]
	rx[0] = (x-xCenter)*c + (y-yCenter)*s
	ry[0] = -(x-xCenter)*s + (y-yCenter)*c

	x = Hit2x_List[idx]
	y = Hit2y_List[idx]
	rx[1] = (x-xCenter)*c + (y-yCenter)*s
	ry[1] = -(x-xCenter)*s + (y-yCenter)*c

	AppendToGraph ry vs rx

	return 0
End

//==============================================================================
// SNS_ExtractModesForFolder_AllRays
//
// Purpose:
//   Figure-oriented version of SNS_ExtractModesForFolder.
//
//   It reuses the existing geometry/mode extraction logic, but changes the
//   topography display so that:
//      - ray-traced SNS channels are drawn faintly with fixed sparsity,
//      - the gap-closing / largest-W_eff trajectory is highlighted,
//      - the longest trajectory is highlighted,
//      - STS position, B-field arrow, textbox, and mode-property plot are kept.
//
//   Faint-ray selection:
//      raySparsity = 1   : plot every ray
//      raySparsity = 10  : plot every 10th ray
//
//   The first faint ray is always rayStartIndex. This makes the displayed
//   subset reproducible and independent of the total number of rays.
//
// Added:
//   Textbox now also reports:
//      E0(ell_cluster,max)
//   using the full zero-field ABS equation via:
//      SNS_ABSZeroFieldEnergy_meV_FromEll(ell_nm)
//
// Inputs:
//   dfPath : input folder containing w_mask and image data.
//
// Optional:
//   Bangle_deg        : B-field angle in degrees. Default: 225.
//   BTK_barrier       : BTK barrier parameter.
//   STSx, STSy        : STS position in nm.
//   Vortexx/y         : kept for interface compatibility.
//   raySparsity       : plot every raySparsity-th faint ray.
//                       Default: 1.
//   rayStartIndex     : first faint-ray index to plot.
//                       Default: 0.
//   faintColor        : base color string for faint rays.
//                       Supported values: "gray", "black", "white", "red",
//                       "blue", "green", "cyan", "magenta", "yellow".
//                       Default: "gray".
//   colorTransparency : alpha/transparency control for faint rays.
//                       Range is clipped to [0,1].
//                       0 = fully transparent.
//                       1 = fully opaque.
//                       Default: 0.20.
//   faintLineSize     : line size for faint rays. Default: 0.25.
//   useExistingChannels
//                    : 0 default, rebuild channels with SNS_ExtractModesForFolder.
//                      1, display canonical channel waves already in dfPath.
//
// Outputs:
//   Same geometry variables as SNS_ExtractModesForFolder.
//
// Requires:
//   SNS_PlotAllChannelRays_Faint(...)
//==============================================================================
Function SNS_ExtractModesForFolder_AllRays(dfPath, [Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy, raySparsity, rayStartIndex, faintColor, colorTransparency, faintLineSize, useExistingChannels])
    String dfPath
    Variable Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy
    Variable raySparsity, rayStartIndex, colorTransparency, faintLineSize
    Variable useExistingChannels
    String faintColor

    String oldDF = GetDataFolder(1)

    Variable r=65535, g=65535, b=65535, a=65535
    String traceList, lastTrace

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    // ---------------- Defaults ----------------
    if (ParamIsDefault(Bangle_deg))
        Bangle_deg = 225
    endif
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
    if (ParamIsDefault(raySparsity))
        raySparsity = 1
    endif
    if (ParamIsDefault(rayStartIndex))
        rayStartIndex = 0
    endif
    if (ParamIsDefault(faintColor))
        faintColor = "gray"
    endif
    if (ParamIsDefault(colorTransparency))
        colorTransparency = 0.20
    endif
    if (ParamIsDefault(faintLineSize))
        faintLineSize = 0.25
    endif
    if (ParamIsDefault(useExistingChannels))
        useExistingChannels = 0
    endif

    raySparsity = max(1, round(raySparsity))
    rayStartIndex = max(0, round(rayStartIndex))
    useExistingChannels = max(0, min(1, round(useExistingChannels)))

    // -------------------------------------------------------------------------
    // 1. Reuse existing extraction logic, but suppress its display.
    // -------------------------------------------------------------------------
    if (!useExistingChannels)
        SNS_ExtractModesForFolder(dfPath, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, Vortexx=Vortexx, Vortexy=Vortexy, doDisplay=0)
    endif

    SetDataFolder $dfPath

    // ---------------- Required waves ----------------
    Wave/Z w_mask
    if (!WaveExists(w_mask))
        SetDataFolder $oldDF
        Abort "SNS_ExtractModesForFolder_AllRays: missing w_mask."
    endif

    String imgList = WaveList("*_Z_mbgnd_xy",";","DIMS:2")
    String imgName = ""
    if (strlen(imgList)==0)
        Wave/Z maskImageFallback = w_mask
        if (WaveExists(maskImageFallback))
            imgName = "w_mask"
        else
            SetDataFolder $oldDF
            Abort "SNS_ExtractModesForFolder_AllRays: no *_Z_mbgnd_xy image or w_mask found in folder."
        endif
    else
        imgName = StringFromList(0, imgList, ";")
    endif

    Wave image = $imgName

    Wave/Z L_N_List = L_N_List_nm
    Wave/Z W_eff_List = W_eff_List_nm
    Wave/Z Hit1x_List = Hit1x_List_nm
    Wave/Z Hit1y_List = Hit1y_List_nm
    Wave/Z Hit2x_List = Hit2x_List_nm
    Wave/Z Hit2y_List = Hit2y_List_nm
    Wave/Z T_eff_List = T_eff_List

    if (!WaveExists(L_N_List) || !WaveExists(W_eff_List) || !WaveExists(Hit1x_List) || !WaveExists(Hit1y_List) || !WaveExists(Hit2x_List) || !WaveExists(Hit2y_List))
        SetDataFolder $oldDF
        Abort "SNS_ExtractModesForFolder_AllRays: missing canonical channel waves."
    endif

    Variable nRayAll = numpnts(L_N_List)
    if ((nRayAll <= 0) || (nRayAll != numpnts(W_eff_List)) || (nRayAll != numpnts(Hit1x_List)) || (nRayAll != numpnts(Hit1y_List)) || (nRayAll != numpnts(Hit2x_List)) || (nRayAll != numpnts(Hit2y_List)))
        SetDataFolder $oldDF
        Abort "SNS_ExtractModesForFolder_AllRays: inconsistent channel wave lengths."
    endif

    Variable v_Gap0_local, v_Longest_local
    NVAR/Z v_heff_local = root:SNS_Settings:lambdaL // [m] Historically still directly referring to london penetration depth.
    
    NVAR/Z v_Gap0
    if (NVAR_Exists(v_Gap0))
        v_Gap0_local = v_Gap0
    else
        WaveStats/Q W_eff_List
        v_Gap0_local = V_maxloc
    endif
    NVAR/Z v_Longest
    if (NVAR_Exists(v_Longest))
        v_Longest_local = v_Longest
    else
        WaveStats/Q L_N_List
        v_Longest_local = V_maxloc
    endif
    v_Gap0_local = round(v_Gap0_local)
    if (numtype(v_Gap0_local) != 0 || v_Gap0_local < 0 || v_Gap0_local >= nRayAll)
        WaveStats/Q W_eff_List
        v_Gap0_local = round(V_maxloc)
    endif
    v_Longest_local = round(v_Longest_local)
    if (numtype(v_Longest_local) != 0 || v_Longest_local < 0 || v_Longest_local >= nRayAll)
        WaveStats/Q L_N_List
        v_Longest_local = round(V_maxloc)
    endif
    v_Gap0_local = max(0, min(nRayAll - 1, v_Gap0_local))
    v_Longest_local = max(0, min(nRayAll - 1, v_Longest_local))
    NVAR/Z v_B0_ballistic_T
    Variable v_B0_ballistic_T_local = NaN
    if (NVAR_Exists(v_B0_ballistic_T))
        v_B0_ballistic_T_local = v_B0_ballistic_T    
    endif
    

    // -------------------------------------------------------------------------
    // Zero-field ABS energy from full equation for L tailCluster,max.
    //
    // The tailCluster diagnostic is resolved through SNS_RayDiagValue(...),
    // which maps canonical names first and legacy names second without falling
    // back to global maxima.
    // -------------------------------------------------------------------------
    String rayHistDF = dfPath + ":RayTraceHist:"
    Variable LclusterMax_nm = SNS_RayDiagValue(rayHistDF, "L", "tailCluster", "max")

    Variable E0_Lcluster_meV = NaN
    if (numtype(LclusterMax_nm) == 0 && LclusterMax_nm > 0)
        E0_Lcluster_meV = SNS_ABSZeroFieldEnergy_meV_FromEll(LclusterMax_nm)
    endif   

    Variable/G v_E0_Lcluster_meV = E0_Lcluster_meV
    Variable/G v_LclusterMax_forE0_nm = LclusterMax_nm
    
    // -------------------------------------------------------------------------
    // Zero-bias "gap closing" condition from W tailCluster,max.
    // -------------------------------------------------------------------------
    Variable v_Weff_max_nm = SNS_RayDiagValue(rayHistDF, "W", "tailCluster", "max")

    if (numtype(v_Weff_max_nm) == 0 && v_Weff_max_nm > 0 && numtype(v_heff_local) == 0 && v_heff_local > 0)
        Variable W_tailCluster_m = v_Weff_max_nm * 1e-9
        v_B0_ballistic_T_local = phi_0/(2*v_heff_local*W_tailCluster_m)
    endif

    // -------------------------------------------------------------------------
    // 2. Topography + faint sparse rays + highlighted rays
    // -------------------------------------------------------------------------
    SNS_DisplayWithScales(image, cmap="grayC")
    String winImage = WinName(0, 1, 1)
    String/G s_AllRays_winImage = winImage

    SNS_StyleRayDisplayImage(image, winImage)

    // Draw sparse faint ray subset first, so highlighted trajectories sit on top.
    SNS_PlotAllChannelRays_Faint(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, raySparsity=raySparsity, rayStartIndex=rayStartIndex, faintColor=faintColor, colorTransparency=colorTransparency, lineSize=faintLineSize)

    // ---------------- Highlight gap-closing / largest-W_eff mode ----------------
    SNS_PlotChannelRay(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, v_Gap0_local)

    traceList = TraceNameList("", ";", 1)
    lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

    ModifyGraph rgb($lastTrace)=(0,0,0)
    ModifyGraph lsize($lastTrace)=2.5

    // ---------------- Highlight longest trajectory ----------------
    SNS_PlotChannelRay(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, v_Longest_local)

    traceList = TraceNameList("", ";", 1)
    lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

    ModifyGraph rgb($lastTrace)=(65535,0,0)
    ModifyGraph lsize($lastTrace)=2.5

    // -------------------------------------------------------------------------
    // 3. Textbox summary
    // -------------------------------------------------------------------------
    String txtBallistic, txtUsadel, txtBox

    txtBallistic = "Ballistic\r" + \
        "ℓ = " + num2str(round(L_N_List[v_Longest_local])) + " nm\r" + \
        "w = " + num2str(round(W_eff_List[v_Gap0_local])) + " nm\r"

//    txtBallistic = "Ballistic\r" + \
//        "L = " + num2str(round(L_N_List[v_Longest_local])) + " nm\r" + \
//        "w = " + num2str(round(W_eff_List[v_Gap0_local])) + " nm\r"

    if (numtype(v_B0_ballistic_T_local) == 0)
        txtBallistic += "B\\B0\\M = " + num2str(round(1000*v_B0_ballistic_T_local)/1000) + " T"
    else
        txtBallistic += "B\\B0\\M = n/a"
    endif

    if (numtype(E0_Lcluster_meV) == 0)
        txtBallistic += "\rE\\B0\\M(ℓ\\Bcl,max\\M) = " + num2str(round(1000*E0_Lcluster_meV)/1000) + " meV"
        txtBallistic += "\rℓ\\Bcl,max\\M = " + num2str(round(LclusterMax_nm)) + " nm"
    else
        txtBallistic += "\rE\\B0\\M(ℓ\\Bcl,max\\M) = n/a"
    endif

    // ---------------- Optional Usadel summary ----------------
    txtUsadel = ""

    NVAR/Z v_Usadel_Ldiff_nm
    NVAR/Z v_Usadel_xProbe_nm
    NVAR/Z v_Usadel_Lphi_use_nm
    NVAR/Z v_Usadel_Bpi_T

    if (NVAR_Exists(v_Usadel_Ldiff_nm) && NVAR_Exists(v_Usadel_xProbe_nm) && NVAR_Exists(v_Usadel_Lphi_use_nm) && NVAR_Exists(v_Usadel_Bpi_T))

        txtUsadel = "Usadel 1D\r" + \
            "L\\Bdiff\\M = " + num2str(round(v_Usadel_Ldiff_nm)) + " nm\r" + \
            "x\\Bprobe\\M = " + num2str(round(v_Usadel_xProbe_nm)) + " nm\r" + \
            "L\\Bφ\\M = " + num2str(round(v_Usadel_Lphi_use_nm)) + " nm\r" + \
            "B\\B0\\M = " + num2str(round(1000*v_Usadel_Bpi_T)/1000) + " T"

        NVAR/Z v_Usadel_D_cm2_per_s
        NVAR/Z v_Usadel_D_cm2_s
        NVAR/Z v_Usadel_D_m2_s

        if (NVAR_Exists(v_Usadel_D_cm2_per_s))
            txtUsadel += "\rD = " + num2str(round(100*v_Usadel_D_cm2_per_s)/100) + " cm\\S2\\M/s"
        elseif (NVAR_Exists(v_Usadel_D_cm2_s))
            txtUsadel += "\rD = " + num2str(round(100*v_Usadel_D_cm2_s)/100) + " cm\\S2\\M/s"
        elseif (NVAR_Exists(v_Usadel_D_m2_s))
            txtUsadel += "\rD = " + num2str(round(100*v_Usadel_D_m2_s*1e4)/100) + " cm\\S2\\M/s"
        else
            NVAR/Z lmfp_nm = root:v_cfg_lmfp_nm
            NVAR/Z vF_mps = root:SNS_Settings:vF

            if (NVAR_Exists(lmfp_nm) && NVAR_Exists(vF_mps))
                Variable diffusionConstant_m2_per_s = 0.5 * vF_mps * lmfp_nm * 1e-9
                Variable diffusionConstant_cm2_per_s = diffusionConstant_m2_per_s * 1e4

                Variable/G v_Usadel_D_m2_per_s = diffusionConstant_m2_per_s
                Variable/G v_Usadel_D_cm2_per_s = diffusionConstant_cm2_per_s

                txtUsadel += "\rD = " + num2str(round(100*diffusionConstant_cm2_per_s)/100) + " cm\\S2\\M/s"
            endif
        endif
    endif

    txtBox = txtBallistic

    if (strlen(txtUsadel) > 0)
        txtBox += "\r\r" + txtUsadel
    endif

    TextBox/C/N=TextMode/X=0/Y=0/F=0/B=(r,g,b,a*0.6)/A=RT txtBox

    // -------------------------------------------------------------------------
    // 4. STS marker / setpoint / B arrow
    // -------------------------------------------------------------------------
    Make/O/N=1 tmpSTSx=STSx, tmpSTSy=STSy

    AppendToGraph tmpSTSy vs tmpSTSx

    traceList = TraceNameList("", ";", 1)
    lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

    ModifyGraph mode($lastTrace)=3
    ModifyGraph marker($lastTrace)=19
    ModifyGraph msize($lastTrace)=5
    ModifyGraph rgb($lastTrace)=(65535,65535,65535)

    SNS_TAGSetpoint()

    Variable xMin = DimOffset(image, 0)
    Variable xMax = xMin + DimDelta(image, 0) * (DimSize(image, 0) - 1)

    Variable yMin = DimOffset(image, 1)
    Variable yMax = yMin + DimDelta(image, 1) * (DimSize(image, 1) - 1)

    Variable pos_x = xMax - 0.15 * (xMax - xMin)
    Variable pos_y = yMin + 0.15 * (yMax - yMin)

    SNS_DrawImageArrow(image, Bangle_deg, "B", pos_x=pos_x, pos_y=pos_y, scale=1.5, color="black")

    GetWindow $WinName(0, 1, 1), wsize
    Variable wPx = abs(V_right - V_left)
    Variable hPx = abs(V_bottom - V_top)

    Variable maxDim = max(wPx, hPx)
    if (maxDim > 512)
        Variable ex = 512 / maxDim
        ModifyGraph/W=$WinName(0, 1, 1) expand=ex
    else
        ModifyGraph/W=$WinName(0, 1, 1) expand=1
    endif

    // -------------------------------------------------------------------------
    // 5. Mode property plot
    // -------------------------------------------------------------------------
    Display L_N_List, W_eff_List
    String winModes = WinName(0, 1, 1)
    String/G s_AllRays_winModes = winModes

    if (WaveExists(T_eff_List))
        AppendToGraph/R T_eff_List
    endif

    ModifyGraph/W=$winModes rgb(L_N_List_nm)=(0,0,0)
    ModifyGraph/W=$winModes rgb(W_eff_List_nm)=(1,16019,65535)
    ModifyGraph/W=$winModes mode(L_N_List_nm)=2
    ModifyGraph/W=$winModes mode(W_eff_List_nm)=2
    ModifyGraph/W=$winModes lsize=5
    if (WaveExists(T_eff_List))
        ModifyGraph/W=$winModes mode(T_eff_List)=3
        ModifyGraph/W=$winModes marker(T_eff_List)=8
    endif

    if (WaveExists(T_eff_List))
        SetAxis/W=$winModes right 0,1
    endif
    SetAxis/W=$winModes left 0,*

    ModifyGraph/W=$winModes muloffset(L_N_List_nm)={0,1}
    ModifyGraph/W=$winModes muloffset(W_eff_List_nm)={0,1}
    if (WaveExists(T_eff_List))
        ModifyGraph/W=$winModes muloffset(T_eff_List)={0,0}
    endif
    ModifyGraph/W=$winModes tick=2
    ModifyGraph/W=$winModes mirror(bottom)=2
    ModifyGraph/W=$winModes ZisZ=1
    ModifyGraph/W=$winModes standoff(left)=0
    ModifyGraph/W=$winModes standoff(bottom)=0

    Label/W=$winModes left   "Length (nm)"
    Label/W=$winModes bottom "Trajectory Nr."
    if (WaveExists(T_eff_List))
        Label/W=$winModes right  "Transparency"
    endif

    Legend/W=$winModes/C/N=text0/J/F=0/A=RT "Mode\r" + \
        "\\s(L_N_List_nm) length (L)\r" + \
        "\\s(W_eff_List_nm) projection ⊥B\\Bext\\M (W)\r" + \
        "\\s(T_eff_List) transparency (T)"

    ModifyGraph/W=$winModes margin(left)=33
    ModifyGraph/W=$winModes margin(bottom)=30
    ModifyGraph/W=$winModes margin(right)=35
    ModifyGraph/W=$winModes margin(top)=6
    ModifyGraph/W=$winModes width=220
    ModifyGraph/W=$winModes height=220

    // ---------------- Restore folder ----------------
    SetDataFolder $oldDF

    return 0
End


//==============================================================================
// SNS_MakeWGapRayContributions
//
// For a found w_gap, compute each trajectory's experimental-resolution-broadened
// zero-bias contribution:
//
//   A_i = weight_i * K_exp[ Sphi_i*pi*(w_i/w_gap - 1) ]
//
// Also creates:
//   outBase + "_contribution"
//   outBase + "_family_mask"     A_i >= 0.5*Amax
//   outBase + "_alpha"           alphaMax at Amax, alphaMin at HWHM
//
// W_nm and L_nm must be matched trajectory waves in nm.
//==============================================================================
Function SNS_MakeWGapRayContributions(W_nm, L_nm, wGap_nm, outBase, [weightWave, Emax_meV, dE_meV, alphaMin, alphaMax])
	Wave W_nm, L_nm
	Variable wGap_nm
	String outBase
	Wave weightWave
	Variable Emax_meV, dE_meV, alphaMin, alphaMax

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	if (ParamIsDefault(Emax_meV))
		Emax_meV = 3
	endif
	if (ParamIsDefault(dE_meV))
		dE_meV = 0.001
	endif
	if (ParamIsDefault(alphaMin))
		alphaMin = 0.10
	endif
	if (ParamIsDefault(alphaMax))
		alphaMax = 1.00
	endif

	Variable n = numpnts(W_nm)
	if (n <= 0 || numpnts(L_nm) != n)
		Abort "SNS_MakeWGapRayContributions: W_nm and L_nm must have same nonzero length."
	endif

	Variable hasWeights = !ParamIsDefault(weightWave) && WaveExists(weightWave) && (numpnts(weightWave) == n)

	if (numtype(wGap_nm) != 0 || wGap_nm <= 0)
		Abort "SNS_MakeWGapRayContributions: invalid wGap_nm."
	endif

	alphaMin = max(0, min(1, alphaMin))
	alphaMax = max(0, min(1, alphaMax))

	// --- Build experimental resolution kernel K_exp(E) ---
	Variable nE = 2*ceil(Emax_meV/dE_meV) + 1
	Variable dE_eV = dE_meV * 1e-3
	Variable E0_eV = -0.5*(nE-1)*dE_eV

	Make/O/D/N=(nE) root:SNS_Settings:tmp_Wgap_E_axis_eV
	Wave E_axis = root:SNS_Settings:tmp_Wgap_E_axis_eV
	E_axis = E0_eV + p*dE_eV
	SetScale/P x, E0_eV, dE_eV, "eV", E_axis

	SNS_MakeThermalKernel(E_axis, params.T_K, "root:SNS_Settings:tmp_Wgap_KT")
	SNS_MakeModulationKernel(E_axis, params.V_mod, "root:SNS_Settings:tmp_Wgap_KM")

	Wave KT = root:SNS_Settings:tmp_Wgap_KT
	Wave KM = root:SNS_Settings:tmp_Wgap_KM

	Duplicate/O KT, root:SNS_Settings:tmp_Wgap_Kexp
	Wave Kexp = root:SNS_Settings:tmp_Wgap_Kexp

	Convolve/A KM, Kexp

	WaveStats/Q Kexp
	if (V_max > 0)
		Kexp /= V_max
	endif

	// --- Output waves ---
	Make/O/D/N=(n) $(outBase + "_contribution")
	Make/O/D/N=(n) $(outBase + "_family_mask")
	Make/O/D/N=(n) $(outBase + "_alpha")

	Wave contrib = $(outBase + "_contribution")
	Wave mask = $(outBase + "_family_mask")
	Wave alpha = $(outBase + "_alpha")

	contrib = 0
	mask = 0
	alpha = 0

	Variable i, wi, elli, Sphi, Ei_meV, Ei_eV, wt

	for (i = 0; i < n; i += 1)
		wi = abs(W_nm[i])
		elli = L_nm[i]

		if (numtype(wi) != 0 || numtype(elli) != 0 || wi <= 0 || elli <= 0)
			continue
		endif

		Sphi = SNS_PhaseSlope_meVPerRad_FromParams(elli)
		if (numtype(Sphi) != 0 || Sphi <= 0)
			continue
		endif

		if (hasWeights)
			wt = weightWave[i]
			if (numtype(wt) != 0 || wt < 0)
				wt = 0
			endif
		else
			wt = 1
		endif

		Ei_meV = Sphi*pi*(wi/wGap_nm - 1)
		Ei_eV = Ei_meV * 1e-3

		if (abs(Ei_meV) <= Emax_meV)
			contrib[i] = wt * interp(Ei_eV, E_axis, Kexp)
		endif
	endfor

	WaveStats/Q contrib
	Variable Amax = V_max
	Variable Ahwhm = 0.5*Amax

	if (Amax > 0 && numtype(Amax) == 0)
		for (i = 0; i < n; i += 1)
			if (contrib[i] >= Ahwhm)
				mask[i] = 1

				if (Amax > Ahwhm)
					alpha[i] = alphaMin + (alphaMax-alphaMin)*(contrib[i]-Ahwhm)/(Amax-Ahwhm)
				else
					alpha[i] = alphaMax
				endif

				alpha[i] = max(alphaMin, min(alphaMax, alpha[i]))
			endif
		endfor
	endif

	Variable/G $(outBase + "_wGap_nm") = wGap_nm
	Variable/G $(outBase + "_Amax") = Amax
	Variable/G $(outBase + "_Ahwhm") = Ahwhm
	Variable/G $(outBase + "_alphaMin") = alphaMin
	Variable/G $(outBase + "_alphaMax") = alphaMax
	Variable/G $(outBase + "_T_K") = params.T_K
	Variable/G $(outBase + "_Vmod_V") = params.V_mod
	Variable/G $(outBase + "_Delta_eff_meV") = 1e3*params.Delta

	KillWaves/Z root:SNS_Settings:tmp_Wgap_E_axis_eV
	KillWaves/Z root:SNS_Settings:tmp_Wgap_KT
	KillWaves/Z root:SNS_Settings:tmp_Wgap_KM
	KillWaves/Z root:SNS_Settings:tmp_Wgap_Kexp

	return 0
End

//==============================================================================
// SNS_ExtractModesForFolder_BGapFamilyRays
//
// Field-space analogue of SNS_ExtractModesForFolder_WGapFamilyRays.
//
// Defines the gap-defining family from:
//
//   rho(B) = sum_i a_i K_exp[ Sphi_i * ( q(B)*w_i - pi ) ]
//
// The peak gives B_gap directly. The displayed w_gap is then:
//   w_gap = Phi0/(2*h_eff*B_gap)
//
// Optional display rotation:
//   displayXAxisAngle_deg sets the displayed x-axis direction in the original
//   image coordinate system.
//      displayXAxisAngle_deg = Bangle_deg      -> x-axis parallel to B
//      displayXAxisAngle_deg = Bangle_deg + 90 -> x-axis perpendicular to B
//
// Requires:
//   SNS_MakeRotatedImageForDisplay(...)
//   SNS_PlotChannelRay_Rotated(...)
//==============================================================================
Function SNS_ExtractModesForFolder_BGapFamilyRays(dfPath, [Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy, hEff_nm, BMin_T, BMax_T, dB_T, peakMinRelHeight, peakBox, alphaMin, alphaMax, alphaGamma, familyColor, lineSize, useWeights, displayXAxisAngle_deg])
	String dfPath
	Variable Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy
	Variable hEff_nm, BMin_T, BMax_T, dB_T
	Variable peakMinRelHeight, peakBox
	Variable alphaMin, alphaMax, alphaGamma
	Variable lineSize, useWeights
	Variable displayXAxisAngle_deg
	String familyColor

	String oldDF = GetDataFolder(1)

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	if (ParamIsDefault(Bangle_deg))
		Bangle_deg = 225
	endif
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
	if (ParamIsDefault(hEff_nm))
		NVAR/Z hEffCfg = root:v_cfg_h_eff_nm
		if (NVAR_Exists(hEffCfg))
			hEff_nm = hEffCfg
		else
			Abort "SNS_ExtractModesForFolder_BGapFamilyRays: pass hEff_nm or define root:v_cfg_h_eff_nm."
		endif
	endif
	if (ParamIsDefault(dB_T))
		dB_T = 0.001
	endif
	if (ParamIsDefault(peakMinRelHeight))
		peakMinRelHeight = 0.05
	endif
	if (ParamIsDefault(peakBox))
		peakBox = 3
	endif
	if (ParamIsDefault(alphaMin))
		alphaMin = 0.10
	endif
	if (ParamIsDefault(alphaMax))
		alphaMax = 1.00
	endif
	if (ParamIsDefault(alphaGamma))
		alphaGamma = 3.0
	endif
	if (ParamIsDefault(familyColor))
		familyColor = "cyan"
	endif
	if (ParamIsDefault(lineSize))
		lineSize = 1.0
	endif
	if (ParamIsDefault(useWeights))
		useWeights = 0
	endif

	Variable rotateDisplay = 0
	if (!ParamIsDefault(displayXAxisAngle_deg) && numtype(displayXAxisAngle_deg) == 0)
		rotateDisplay = 1
	endif

	dB_T = max(1e-5, dB_T)
	peakBox = max(1, round(peakBox))
	alphaMin = max(0, min(1, alphaMin))
	alphaMax = max(0, min(1, alphaMax))
	alphaGamma = max(0.1, alphaGamma)

	// -------------------------------------------------------------------------
	// 1. Extract rays, no display.
	// -------------------------------------------------------------------------
	SNS_ExtractModesForFolder(dfPath, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, Vortexx=Vortexx, Vortexy=Vortexy, doDisplay=0)

	SetDataFolder $dfPath

	Wave/Z w_mask
	if (!WaveExists(w_mask))
		SetDataFolder $oldDF
		Abort "SNS_ExtractModesForFolder_BGapFamilyRays: missing w_mask."
	endif

	String imgList = WaveList("*_Z_mbgnd_xy", ";", "DIMS:2")
	if (strlen(imgList) == 0)
		SetDataFolder $oldDF
		Abort "SNS_ExtractModesForFolder_BGapFamilyRays: no *_Z_mbgnd_xy image found."
	endif

	String imgName = StringFromList(0, imgList, ";")
	Wave image = $imgName

	Wave L_N_List
	Wave W_eff_List
	Wave Hit1x_List
	Wave Hit1y_List
	Wave Hit2x_List
	Wave Hit2y_List
	Wave T_eff_List

	Variable nRay = numpnts(W_eff_List)

	Make/O/D/N=(nRay) tmp_W_eff_List_nm_forBGap
	Make/O/D/N=(nRay) tmp_L_N_List_nm_forBGap

	tmp_W_eff_List_nm_forBGap = abs(W_eff_List) * 1e9
	tmp_L_N_List_nm_forBGap = L_N_List * 1e9

	Wave W_nm = tmp_W_eff_List_nm_forBGap
	Wave L_nm = tmp_L_N_List_nm_forBGap

	// -------------------------------------------------------------------------
	// 2. Field density and B-gap peak.
	// -------------------------------------------------------------------------
	if (ParamIsDefault(BMin_T) && ParamIsDefault(BMax_T))
		if (useWeights)
			SNS_MakeBResolutionDensity(W_nm, L_nm, "B_gap_density_exp_T", hEff_nm, weightWave=T_eff_List, dB_T=dB_T, doDisplay=0)
		else
			SNS_MakeBResolutionDensity(W_nm, L_nm, "B_gap_density_exp_T", hEff_nm, dB_T=dB_T, doDisplay=0)
		endif
	elseif (ParamIsDefault(BMax_T))
		if (useWeights)
			SNS_MakeBResolutionDensity(W_nm, L_nm, "B_gap_density_exp_T", hEff_nm, weightWave=T_eff_List, BMin_T=BMin_T, dB_T=dB_T, doDisplay=0)
		else
			SNS_MakeBResolutionDensity(W_nm, L_nm, "B_gap_density_exp_T", hEff_nm, BMin_T=BMin_T, dB_T=dB_T, doDisplay=0)
		endif
	else
		if (ParamIsDefault(BMin_T))
			BMin_T = 0
		endif

		if (useWeights)
			SNS_MakeBResolutionDensity(W_nm, L_nm, "B_gap_density_exp_T", hEff_nm, weightWave=T_eff_List, BMin_T=BMin_T, BMax_T=BMax_T, dB_T=dB_T, doDisplay=0)
		else
			SNS_MakeBResolutionDensity(W_nm, L_nm, "B_gap_density_exp_T", hEff_nm, BMin_T=BMin_T, BMax_T=BMax_T, dB_T=dB_T, doDisplay=0)
		endif
	endif

	Wave rho = B_gap_density_exp_T

	Variable BGap_T
	BGap_T = SNS_FindFirstBPeakInResolutionDensity(rho, minRelHeight=peakMinRelHeight, box=peakBox, doDisplay=0)

	if (numtype(BGap_T) != 0)
		SetDataFolder $oldDF
		Abort "SNS_ExtractModesForFolder_BGapFamilyRays: failed to find B_gap peak."
	endif

	Variable wFromB_nm = SNS_Wnm_FromB_T(BGap_T, hEff_nm)

	// -------------------------------------------------------------------------
	// 3. Ray family at B_gap.
	// -------------------------------------------------------------------------
	if (useWeights)
		SNS_MakeBGapRayContributions(W_nm, L_nm, BGap_T, hEff_nm, "B_gap_family", weightWave=T_eff_List, alphaMin=alphaMin, alphaMax=alphaMax, alphaGamma=alphaGamma)
	else
		SNS_MakeBGapRayContributions(W_nm, L_nm, BGap_T, hEff_nm, "B_gap_family", alphaMin=alphaMin, alphaMax=alphaMax, alphaGamma=alphaGamma)
	endif

	Wave contrib = B_gap_family_contribution
	Wave mask = B_gap_family_family_mask
	Wave alpha = B_gap_family_alpha

	WaveStats/Q contrib
	Variable Amax = V_max

	Variable i
	Variable bestTrace = NaN
	Variable bestA = -Inf

	for (i = 0; i < nRay; i += 1)
		if (contrib[i] > bestA)
			bestA = contrib[i]
			bestTrace = i
		endif
	endfor

	// -------------------------------------------------------------------------
	// 4. Topography display, optionally rotated.
	// -------------------------------------------------------------------------
	Variable xOrig0 = DimOffset(image, 0)
	Variable xOrig1 = xOrig0 + DimDelta(image, 0) * (DimSize(image, 0) - 1)
	Variable yOrig0 = DimOffset(image, 1)
	Variable yOrig1 = yOrig0 + DimDelta(image, 1) * (DimSize(image, 1) - 1)

	Variable xCenterRot = 0.5 * (min(xOrig0, xOrig1) + max(xOrig0, xOrig1))
	Variable yCenterRot = 0.5 * (min(yOrig0, yOrig1) + max(yOrig0, yOrig1))

	String displayImgName = NameOfWave(image)

	if (rotateDisplay)
		SNS_MakeRotatedImageForDisplay(image, "tmp_BGapFamily_rot_image", displayXAxisAngle_deg)
		displayImgName = "tmp_BGapFamily_rot_image"
	endif

	Wave displayImage = $displayImgName

	SNS_DisplayWithScales(displayImage, cmap="grayC")
	String winImage = WinName(0, 1, 1)

	SNS_StyleRayDisplayImage(displayImage, winImage)

	DoWindow/F $winImage

	Variable rCol, gCol, bCol

	if (StringMatch(familyColor, "black"))
		rCol = 0; gCol = 0; bCol = 0
	elseif (StringMatch(familyColor, "white"))
		rCol = 65535; gCol = 65535; bCol = 65535
	elseif (StringMatch(familyColor, "red"))
		rCol = 65535; gCol = 0; bCol = 0
	elseif (StringMatch(familyColor, "blue"))
		rCol = 0; gCol = 0; bCol = 65535
	elseif (StringMatch(familyColor, "green"))
		rCol = 0; gCol = 45000; bCol = 0
	elseif (StringMatch(familyColor, "magenta"))
		rCol = 65535; gCol = 0; bCol = 65535
	elseif (StringMatch(familyColor, "yellow"))
		rCol = 65535; gCol = 50000; bCol = 0
	else
		rCol = 1; gCol = 16019; bCol = 65535
	endif

	String traceList, lastTrace
	Variable nFamily = 0
	Variable aRGBA

	for (i = 0; i < nRay; i += 1)
		if (mask[i] < 0.5)
			continue
		endif

		DoWindow/F $winImage

		if (rotateDisplay)
			SNS_PlotChannelRay_Rotated(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, i, displayXAxisAngle_deg, xCenterRot, yCenterRot)
		else
			SNS_PlotChannelRay(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, i)
		endif

		traceList = TraceNameList(winImage, ";", 1)
		lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

		aRGBA = round(65535 * alpha[i])
		aRGBA = max(0, min(65535, aRGBA))

		ModifyGraph/W=$winImage rgb($lastTrace)=(rCol,gCol,bCol,aRGBA)
		ModifyGraph/W=$winImage lsize($lastTrace)=lineSize

		nFamily += 1
	endfor

	// Max-contributing ray on top.
	if (numtype(bestTrace) == 0)
		DoWindow/F $winImage

		if (rotateDisplay)
			SNS_PlotChannelRay_Rotated(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, bestTrace, displayXAxisAngle_deg, xCenterRot, yCenterRot)
		else
			SNS_PlotChannelRay(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, bestTrace)
		endif

		traceList = TraceNameList(winImage, ";", 1)
		lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

		ModifyGraph/W=$winImage rgb($lastTrace)=(0,0,0,65535)
		ModifyGraph/W=$winImage lsize($lastTrace)=2.5
	endif

	// -------------------------------------------------------------------------
	// STS marker / B arrow.
	// -------------------------------------------------------------------------
	Make/O/N=1 tmpSTSx, tmpSTSy

	if (rotateDisplay)
		Variable thSTS = displayXAxisAngle_deg * pi / 180
		tmpSTSx[0] = (STSx - xCenterRot) * cos(thSTS) + (STSy - yCenterRot) * sin(thSTS)
		tmpSTSy[0] = -(STSx - xCenterRot) * sin(thSTS) + (STSy - yCenterRot) * cos(thSTS)
	else
		tmpSTSx[0] = STSx
		tmpSTSy[0] = STSy
	endif

	AppendToGraph/W=$winImage tmpSTSy vs tmpSTSx

	traceList = TraceNameList(winImage, ";", 1)
	lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

	ModifyGraph/W=$winImage mode($lastTrace)=3
	ModifyGraph/W=$winImage marker($lastTrace)=19
	ModifyGraph/W=$winImage msize($lastTrace)=5
	ModifyGraph/W=$winImage rgb($lastTrace)=(65535,65535,65535,65535)

	DoWindow/F $winImage
	SNS_TAGSetpoint()

	Variable xMin = DimOffset(displayImage, 0)
	Variable xMax = xMin + DimDelta(displayImage, 0) * (DimSize(displayImage, 0) - 1)
	Variable yMin = DimOffset(displayImage, 1)
	Variable yMax = yMin + DimDelta(displayImage, 1) * (DimSize(displayImage, 1) - 1)

	Variable pos_x = xMax - 0.15 * (xMax - xMin)
	Variable pos_y = yMin + 0.15 * (yMax - yMin)

	Variable displayBangle_deg = Bangle_deg
	if (rotateDisplay)
		displayBangle_deg = Bangle_deg - displayXAxisAngle_deg
	endif

	SNS_DrawImageArrow(displayImage, displayBangle_deg, "B", pos_x=pos_x, pos_y=pos_y, scale=1.5, color="black")

	String txt
	txt = "Gap-defining field family\r"
	txt += "B\\Bgap\\M = " + num2str(round(1000*BGap_T)/1000) + " T\r"
	txt += "w(B\\Bgap\\M) = " + num2str(round(10*wFromB_nm)/10) + " nm\r"
	txt += "h\\Beff\\M = " + num2str(round(10*hEff_nm)/10) + " nm\r"
	txt += "N\\Bfamily\\M = " + num2str(nFamily) + " / " + num2str(nRay) + "\r"
	txt += "criterion: A\\Bi\\M ≥ 0.5 A\\Bmax\\M\r"
	txt += "alpha γ = " + num2str(round(10*alphaGamma)/10)

	if (rotateDisplay)
		txt += "\rdisplay x-axis = " + num2str(round(10 * displayXAxisAngle_deg) / 10) + "°"
	endif

	TextBox/C/N=BGapFamilyInfo/F=0/A=RT/X=0/Y=0/B=(65535,65535,65535,39321) txt

	GetWindow $winImage, wsize
	Variable wPx = abs(V_right - V_left)
	Variable hPx = abs(V_bottom - V_top)
	Variable maxDim = max(wPx, hPx)

	if (maxDim > 512)
		Variable ex = 512 / maxDim
		ModifyGraph/W=$winImage expand=ex
	else
		ModifyGraph/W=$winImage expand=1
	endif

	// -------------------------------------------------------------------------
	// 5. Field-density diagnostic plot.
	// -------------------------------------------------------------------------
	Display/K=1 rho
	String winRho = WinName(0, 1, 1)

	ModifyGraph/W=$winRho mode=0
	ModifyGraph/W=$winRho tick=2, mirror=2, standoff=0
	Label/W=$winRho bottom "B (T)"
	Label/W=$winRho left "resolution-broadened zero-bias weight"

	SetDrawLayer UserFront
	SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(1,16019,65535), linethick=2, dash=3
	DrawLine BGap_T, 0, BGap_T, 1.05

	ModifyGraph/W=$winRho width=260
	ModifyGraph/W=$winRho height=220

	KillWaves/Z tmp_W_eff_List_nm_forBGap
	KillWaves/Z tmp_L_N_List_nm_forBGap

	SetDataFolder $oldDF

	return BGap_T
End

//==============================================================================
// SNS_ExtractModesForFolder_WGapFamilyRays
//
// Like SNS_ExtractModesForFolder_AllRays, but the highlighted ray family is
// defined from the experimental-resolution-broadened zero-bias weight:
//
//   1. extract all rays,
//   2. compute rho_exp(w_c),
//   3. find high-w peak w_gap,
//   4. compute each ray's contribution at w_gap,
//   5. plot all rays with A_i >= 0.5 Amax,
//      alpha = alphaMax at Amax and alphaMin at HWHM,
//      with exponential/power-law fade controlled by alphaGamma.
//==============================================================================
Function SNS_ExtractModesForFolder_WGapFamilyRays(dfPath, [Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy, densityDW_nm, peakMinRelHeight, peakBox, alphaMin, alphaMax, alphaGamma, familyColor, lineSize, useWeights])
	String dfPath
	Variable Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy
	Variable densityDW_nm, peakMinRelHeight, peakBox
	Variable alphaMin, alphaMax, alphaGamma
	Variable lineSize, useWeights
	String familyColor

	String oldDF = GetDataFolder(1)

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	// ---------------- Defaults ----------------
	if (ParamIsDefault(Bangle_deg))
		Bangle_deg = 225
	endif
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
	if (ParamIsDefault(densityDW_nm))
		densityDW_nm = 1
	endif
	if (ParamIsDefault(peakMinRelHeight))
		peakMinRelHeight = 0.05
	endif
	if (ParamIsDefault(peakBox))
		peakBox = 3
	endif
	if (ParamIsDefault(alphaMin))
		alphaMin = 0.10
	endif
	if (ParamIsDefault(alphaMax))
		alphaMax = 1.00
	endif
	if (ParamIsDefault(alphaGamma))
		alphaGamma = 3.0
	endif
	if (ParamIsDefault(familyColor))
		familyColor = "cyan"
	endif
	if (ParamIsDefault(lineSize))
		lineSize = 1.0
	endif
	if (ParamIsDefault(useWeights))
		useWeights = 0
	endif

	densityDW_nm = max(0.05, densityDW_nm)
	peakBox = max(1, round(peakBox))
	alphaMin = max(0, min(1, alphaMin))
	alphaMax = max(0, min(1, alphaMax))
	alphaGamma = max(0.1, alphaGamma)

	// -------------------------------------------------------------------------
	// 1. Reuse existing extraction logic, no display.
	// -------------------------------------------------------------------------
	SNS_ExtractModesForFolder(dfPath, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, Vortexx=Vortexx, Vortexy=Vortexy, doDisplay=0)

	SetDataFolder $dfPath

	Wave/Z w_mask
	if (!WaveExists(w_mask))
		SetDataFolder $oldDF
		Abort "SNS_ExtractModesForFolder_WGapFamilyRays: missing w_mask."
	endif

	String imgList = WaveList("*_Z_mbgnd_xy", ";", "DIMS:2")
	if (strlen(imgList) == 0)
		SetDataFolder $oldDF
		Abort "SNS_ExtractModesForFolder_WGapFamilyRays: no *_Z_mbgnd_xy image found."
	endif

	String imgName = StringFromList(0, imgList, ";")
	Wave image = $imgName

	Wave L_N_List
	Wave W_eff_List
	Wave Hit1x_List
	Wave Hit1y_List
	Wave Hit2x_List
	Wave Hit2y_List
	Wave T_eff_List

	Variable nRay = numpnts(W_eff_List)

	Make/O/D/N=(nRay) tmp_W_eff_List_nm_forWGap
	Make/O/D/N=(nRay) tmp_L_N_List_nm_forWGap

	tmp_W_eff_List_nm_forWGap = abs(W_eff_List) * 1e9
	tmp_L_N_List_nm_forWGap = L_N_List * 1e9

	Wave W_nm = tmp_W_eff_List_nm_forWGap
	Wave L_nm = tmp_L_N_List_nm_forWGap

	// -------------------------------------------------------------------------
	// 2. Compute rho_exp(w_c) and find high-w peak.
	// -------------------------------------------------------------------------
	if (useWeights)
		SNS_MakeWResolutionDensity(W_nm, "W_eff_density_exp_nm", L_nm=L_nm, weightWave=T_eff_List, dw_nm=densityDW_nm, doDisplay=0)
	else
		SNS_MakeWResolutionDensity(W_nm, "W_eff_density_exp_nm", L_nm=L_nm, dw_nm=densityDW_nm, doDisplay=0)
	endif

	Wave rho = W_eff_density_exp_nm

	Variable wGap_nm
	wGap_nm = SNS_FindHighWPeakInResolutionDensity(rho, minRelHeight=peakMinRelHeight, box=peakBox, doDisplay=0,edgeSkipPts=1)

	if (numtype(wGap_nm) != 0)
		SetDataFolder $oldDF
		Abort "SNS_ExtractModesForFolder_WGapFamilyRays: failed to find w_gap peak."
	endif

	// -------------------------------------------------------------------------
	// 3. Compute HWHM-supported ray family at w_gap.
	// -------------------------------------------------------------------------
	if (useWeights)
		SNS_MakeWGapRayContributions(W_nm, L_nm, wGap_nm, "W_gap_family", weightWave=T_eff_List, alphaMin=alphaMin, alphaMax=alphaMax)
	else
		SNS_MakeWGapRayContributions(W_nm, L_nm, wGap_nm, "W_gap_family", alphaMin=alphaMin, alphaMax=alphaMax)
	endif

	Wave contrib = W_gap_family_contribution
	Wave mask = W_gap_family_family_mask
	Wave alpha = W_gap_family_alpha

	WaveStats/Q contrib
	Variable Amax = V_max
	Variable Ahwhm = 0.5 * Amax

	// Recompute alpha with exponential/power-law fade.
	Variable i
	Variable alphaFrac

	alpha = 0
	mask = 0

	if (Amax > 0 && numtype(Amax) == 0)
		for (i = 0; i < nRay; i += 1)
			if (contrib[i] >= Ahwhm)
				mask[i] = 1

				if (Amax > Ahwhm)
					alphaFrac = (contrib[i] - Ahwhm) / (Amax - Ahwhm)
					alphaFrac = max(0, min(1, alphaFrac))
					alpha[i] = alphaMin + (alphaMax - alphaMin) * alphaFrac^alphaGamma
				else
					alpha[i] = alphaMax
				endif

				alpha[i] = max(alphaMin, min(alphaMax, alpha[i]))
			endif
		endfor
	endif

	Variable/G W_gap_family_alphaGamma = alphaGamma

	// Identify max-contributing ray.
	Variable bestTrace = NaN
	Variable bestA = -Inf

	for (i = 0; i < nRay; i += 1)
		if (contrib[i] > bestA)
			bestA = contrib[i]
			bestTrace = i
		endif
	endfor

	// -------------------------------------------------------------------------
	// 4. Topography display.
	// -------------------------------------------------------------------------
	SNS_DisplayWithScales(image, cmap="grayC")
	String winImage = WinName(0, 1, 1)

	SNS_StyleRayDisplayImage(image, winImage)

	// Color choice.
	Variable rCol, gCol, bCol

	if (StringMatch(familyColor, "black"))
		rCol = 0; gCol = 0; bCol = 0
	elseif (StringMatch(familyColor, "white"))
		rCol = 65535; gCol = 65535; bCol = 65535
	elseif (StringMatch(familyColor, "red"))
		rCol = 65535; gCol = 0; bCol = 0
	elseif (StringMatch(familyColor, "blue"))
		rCol = 0; gCol = 0; bCol = 65535
	elseif (StringMatch(familyColor, "green"))
		rCol = 0; gCol = 45000; bCol = 0
	elseif (StringMatch(familyColor, "magenta"))
		rCol = 65535; gCol = 0; bCol = 65535
	elseif (StringMatch(familyColor, "yellow"))
		rCol = 65535; gCol = 50000; bCol = 0
	else
		// cyan default
		rCol = 1; gCol = 16019; bCol = 65535
	endif

	// -------------------------------------------------------------------------
	// 5. Plot all HWHM-supported rays with alpha fade.
	// -------------------------------------------------------------------------
	String traceList, lastTrace
	Variable nFamily = 0
	Variable aRGBA

	for (i = 0; i < nRay; i += 1)
		if (mask[i] < 0.5)
			continue
		endif

		SNS_PlotChannelRay(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, i)

		traceList = TraceNameList(winImage, ";", 1)
		lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

		aRGBA = round(65535 * alpha[i])
		aRGBA = max(0, min(65535, aRGBA))

		ModifyGraph/W=$winImage rgb($lastTrace)=(rCol,gCol,bCol,aRGBA)
		ModifyGraph/W=$winImage lsize($lastTrace)=lineSize

		nFamily += 1
	endfor

	// Highlight max-contributing ray on top.
	if (numtype(bestTrace) == 0)
		SNS_PlotChannelRay(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, bestTrace)

		traceList = TraceNameList(winImage, ";", 1)
		lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

		ModifyGraph/W=$winImage rgb($lastTrace)=(0,0,0,65535)
		ModifyGraph/W=$winImage lsize($lastTrace)=2.5
	endif

	// -------------------------------------------------------------------------
	// 6. STS marker and B arrow.
	// -------------------------------------------------------------------------
	Make/O/N=1 tmpSTSx=STSx, tmpSTSy=STSy

	AppendToGraph/W=$winImage tmpSTSy vs tmpSTSx

	traceList = TraceNameList(winImage, ";", 1)
	lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

	ModifyGraph/W=$winImage mode($lastTrace)=3
	ModifyGraph/W=$winImage marker($lastTrace)=19
	ModifyGraph/W=$winImage msize($lastTrace)=5
	ModifyGraph/W=$winImage rgb($lastTrace)=(65535,65535,65535,65535)

	DoWindow/F $winImage
	SNS_TAGSetpoint()

	Variable xMin = DimOffset(image, 0)
	Variable xMax = xMin + DimDelta(image, 0) * (DimSize(image, 0) - 1)
	Variable yMin = DimOffset(image, 1)
	Variable yMax = yMin + DimDelta(image, 1) * (DimSize(image, 1) - 1)

	Variable pos_x = xMax - 0.15 * (xMax - xMin)
	Variable pos_y = yMin + 0.15 * (yMax - yMin)

	SNS_DrawImageArrow(image, Bangle_deg, "B", pos_x=pos_x, pos_y=pos_y, scale=1.5, color="black")

	// -------------------------------------------------------------------------
	// 7. Textbox.
	// -------------------------------------------------------------------------
	String txt
	txt = "Gap-defining family\r"
	txt += "w\\Bgap\\M = " + num2str(round(10*wGap_nm)/10) + " nm\r"
	txt += "N\\Bfamily\\M = " + num2str(nFamily) + " / " + num2str(nRay) + "\r"
	txt += "criterion: A\\Bi\\M ≥ 0.5 A\\Bmax\\M\r"
	txt += "alpha γ = " + num2str(round(10*alphaGamma)/10) + "\r"
	txt += "A\\Bmax\\M = " + num2str(round(1000*Amax)/1000)

	TextBox/C/N=WGapFamilyInfo/F=0/A=RT/X=0/Y=0/B=(65535,65535,65535,39321) txt

	GetWindow $winImage, wsize
	Variable wPx = abs(V_right - V_left)
	Variable hPx = abs(V_bottom - V_top)
	Variable maxDim = max(wPx, hPx)

	if (maxDim > 512)
		Variable ex = 512 / maxDim
		ModifyGraph/W=$winImage expand=ex
	else
		ModifyGraph/W=$winImage expand=1
	endif

	// -------------------------------------------------------------------------
	// 8. Density diagnostic plot.
	// -------------------------------------------------------------------------
	Display/K=1 rho
	String winRho = WinName(0, 1, 1)

	ModifyGraph/W=$winRho mode=0
	ModifyGraph/W=$winRho tick=2, mirror=2, standoff=0
	Label/W=$winRho bottom "candidate w\\Bc\\M (nm)"
	Label/W=$winRho left "resolution-broadened zero-bias weight"

	SetDrawLayer UserFront
	SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(1,16019,65535), linethick=2, dash=3
	DrawLine wGap_nm, 0, wGap_nm, 1.05

	ModifyGraph/W=$winRho width=260
	ModifyGraph/W=$winRho height=220

	// Keep useful outputs, remove only temporary converted waves.
	KillWaves/Z tmp_W_eff_List_nm_forWGap
	KillWaves/Z tmp_L_N_List_nm_forWGap

	SetDataFolder $oldDF

	return wGap_nm
End


//==============================================================================
// SNS_EllFromABSZeroFieldEnergy_nm
//
// Inverts E0(ell) using the existing full zero-field ABS helper:
//     E0 = SNS_ABSZeroFieldEnergy_meV_FromEll(ell_nm)
//
// Assumes E0 decreases monotonically with ell.
//==============================================================================
Function SNS_EllFromABSZeroFieldEnergy_nm(Etarget_meV, [Lmin_nm, Lmax_nm])
	Variable Etarget_meV, Lmin_nm, Lmax_nm

	if (ParamIsDefault(Lmin_nm))
		Lmin_nm = 1
	endif
	if (ParamIsDefault(Lmax_nm))
		Lmax_nm = 50000
	endif

	Lmin_nm = max(0.1, Lmin_nm)
	Lmax_nm = max(Lmin_nm + 1, Lmax_nm)

	if (numtype(Etarget_meV) != 0)
		return NaN
	endif
	if (Etarget_meV <= 0)
		return Lmax_nm
	endif

	Variable Emin, Emax, Emid, Lmid
	Variable iter

	Emax = SNS_ABSZeroFieldEnergy_meV_FromEll(Lmin_nm)
	Emin = SNS_ABSZeroFieldEnergy_meV_FromEll(Lmax_nm)

	if (numtype(Emax) != 0 || numtype(Emin) != 0)
		return NaN
	endif

	// Target above the maximum accessible energy: clamp to lower length.
	if (Etarget_meV >= Emax)
		return Lmin_nm
	endif

	// Expand upper length until E(Lmax) falls below target.
	do
		Emin = SNS_ABSZeroFieldEnergy_meV_FromEll(Lmax_nm)

		if (numtype(Emin) != 0)
			return NaN
		endif
		if (Emin <= Etarget_meV)
			break
		endif

		Lmax_nm *= 2

		if (Lmax_nm > 1e7)
			return Lmax_nm
		endif
	while (1)

	// Bisection.
	for (iter = 0; iter < 80; iter += 1)
		Lmid = 0.5 * (Lmin_nm + Lmax_nm)
		Emid = SNS_ABSZeroFieldEnergy_meV_FromEll(Lmid)

		if (numtype(Emid) != 0)
			return NaN
		endif

		if (Emid > Etarget_meV)
			Lmin_nm = Lmid
		else
			Lmax_nm = Lmid
		endif
	endfor

	return 0.5 * (Lmin_nm + Lmax_nm)
End


//==============================================================================
// SNS_LResolutionWindowFromCluster_nm
//
// Full-equation experimental L-window around Lcluster.
// Returns DeltaL = Lhigh - Llow.
// Also stores last low/high in root:SNS_Settings.
//==============================================================================
Function SNS_LResolutionWindowFromCluster_nm(Lcluster_nm, EresHWHM_meV, [Lmin_nm, Lmax_nm])
	Variable Lcluster_nm, EresHWHM_meV, Lmin_nm, Lmax_nm

	if (ParamIsDefault(Lmin_nm))
		Lmin_nm = 1
	endif
	if (ParamIsDefault(Lmax_nm))
		Lmax_nm = 50000
	endif

	if (numtype(Lcluster_nm) != 0 || Lcluster_nm <= 0)
		return NaN
	endif
	if (numtype(EresHWHM_meV) != 0 || EresHWHM_meV <= 0)
		return NaN
	endif

	Variable Ecluster = SNS_ABSZeroFieldEnergy_meV_FromEll(Lcluster_nm)
	if (numtype(Ecluster) != 0)
		return NaN
	endif

	Variable Ehigh = Ecluster + EresHWHM_meV
	Variable Elow = Ecluster - EresHWHM_meV

	if (Elow < 0)
		Elow = 0
	endif

	Variable Llow = SNS_EllFromABSZeroFieldEnergy_nm(Ehigh, Lmin_nm=Lmin_nm, Lmax_nm=Lmax_nm)
	Variable Lhigh = SNS_EllFromABSZeroFieldEnergy_nm(Elow, Lmin_nm=Lmin_nm, Lmax_nm=Lmax_nm)

	if (numtype(Llow) != 0 || numtype(Lhigh) != 0)
		return NaN
	endif

	Variable/G root:SNS_Settings:v_last_Lres_Lcluster_nm = Lcluster_nm
	Variable/G root:SNS_Settings:v_last_Lres_Ecluster_meV = Ecluster
	Variable/G root:SNS_Settings:v_last_Lres_EresHWHM_meV = EresHWHM_meV
	Variable/G root:SNS_Settings:v_last_Lres_Llow_nm = Llow
	Variable/G root:SNS_Settings:v_last_Lres_Lhigh_nm = Lhigh
	Variable/G root:SNS_Settings:v_last_Lres_DeltaL_nm = Lhigh - Llow

	return Lhigh - Llow
End


//==============================================================================
// SNS_NumVarByPathOrDefault
//==============================================================================
Function SNS_NumVarByPathOrDefault(varPath, defaultVal)
    String varPath
    Variable defaultVal

    if (Exists(varPath) == 0)
        return defaultVal
    endif

    NVAR/Z val = $varPath
    if (NVAR_Exists(val))
        return val
    endif

    return defaultVal
End


//==============================================================================
// SNS_RayDiagValue
//
// Purpose:
//   Compatibility resolver for ray-family scalar diagnostics.
//
//   This helper centralizes the transition from legacy names such as
//   v_LN_extreme_max_nm / v_Weff_extreme_max_nm to the canonical naming
//   convention documented in:
//
//      SNS_ABS_trajectory_review/RAY_NAMING_CONVENTION.md
//
//   It intentionally uses a small explicit map. It must not fall back across
//   semantic boundaries; for example, a tailCluster request never falls back to
//   a global maximum.
//
// Inputs:
//   histDF    : data-folder path containing RayTraceHist/RayTraceHist_3D output.
//   quantity  : "L", "W", or "Wgeom".
//   selection : "global", "tailCluster", "minbinCluster", or "histCluster".
//   stat      : requested statistic, e.g. "max", "median", "max_trace".
//
// Optional:
//   is3D      : 0 default. If 1, prefer the canonical *_3D_* names and legacy
//               RayTraceHist_3D names.
//   defaultVal: NaN default. Returned when neither canonical nor legacy name
//               exists.
//
// Returns:
//   Numeric diagnostic value or defaultVal.
//==============================================================================
Function SNS_RayDiagValue(histDF, quantity, selection, stat, [is3D, defaultVal])
    String histDF, quantity, selection, stat
    Variable is3D, defaultVal

    if (ParamIsDefault(is3D))
        is3D = 0
    endif
    is3D = round(is3D) != 0

    if (ParamIsDefault(defaultVal))
        defaultVal = NaN
    endif

    String canonical = ""
    String legacy = ""
    String qCanon = ""
    String qLegacy = ""
    String statCanon = ""
    String statLegacy = ""

    Variable isL = !CmpStr(quantity, "L", 2)
    Variable isW = !CmpStr(quantity, "W", 2)
    Variable isWgeom = !CmpStr(quantity, "Wgeom", 2) || !CmpStr(quantity, "W_geom", 2)

    if (!isL && !isW && !isWgeom)
        return defaultVal
    endif

    if (isL)
        if (is3D)
            qCanon = "L_3D"
            qLegacy = "LN_3D"
        else
            qCanon = "L"
            qLegacy = "LN"
        endif
    elseif (isW)
        if (is3D)
            qCanon = "W_3D"
            qLegacy = "Weff_3D"
        else
            qCanon = "W"
            qLegacy = "Weff"
        endif
    else
        if (is3D)
            qCanon = "Wgeom_3D"
            qLegacy = "Wgeom_3D"
        else
            qCanon = "Wgeom"
            qLegacy = "Wgeom"
        endif
    endif

    // -------------------------------------------------------------------------
    // Global diagnostics.
    // -------------------------------------------------------------------------
    if (!CmpStr(selection, "global", 2))
        if (!CmpStr(stat, "max", 2))
            canonical = "v_" + qCanon + "_global_max_nm"
            legacy = "v_" + qLegacy + "_max_nm"
        endif

    // -------------------------------------------------------------------------
    // Supported upper-tail cluster diagnostics.
    // -------------------------------------------------------------------------
    elseif (!CmpStr(selection, "tailCluster", 2))
        if (!CmpStr(stat, "low", 2))
            statCanon = "low_nm"
            statLegacy = "low_nm"
        elseif (!CmpStr(stat, "high", 2))
            statCanon = "high_nm"
            statLegacy = "high_nm"
        elseif (!CmpStr(stat, "n", 2))
            statCanon = "n"
            statLegacy = "n"
        elseif (!CmpStr(stat, "max", 2))
            statCanon = "max_nm"
            statLegacy = "max_nm"
        elseif (!CmpStr(stat, "median", 2))
            statCanon = "median_nm"
            statLegacy = "median_nm"
        elseif (!CmpStr(stat, "sigma", 2))
            statCanon = "sigma_nm"
            statLegacy = "sigma_nm"
        elseif (!CmpStr(stat, "max_trace", 2))
            statCanon = "max_trace"
            statLegacy = "max_trace"
        elseif (!CmpStr(stat, "median_trace", 2))
            statCanon = "median_trace"
            statLegacy = "median_trace"
        elseif (!CmpStr(stat, "start", 2))
            statCanon = "start"
            statLegacy = "cluster_start"
        elseif (!CmpStr(stat, "size", 2))
            statCanon = "size"
            statLegacy = "cluster_size"
        elseif (!CmpStr(stat, "minSize", 2))
            statCanon = "minSize"
            statLegacy = "minClusterSize"
        endif

        if (strlen(statCanon) > 0)
            canonical = "v_" + qCanon + "_tailCluster_" + statCanon
            legacy = "v_" + qLegacy + "_extreme_" + statLegacy
        endif

    // -------------------------------------------------------------------------
    // Smallest occupied L-bin cluster diagnostics.
    // Currently only L has an identified minbinCluster implementation.
    // -------------------------------------------------------------------------
    elseif (!CmpStr(selection, "minbinCluster", 2))
        if (!isL)
            return defaultVal
        endif

        if (!CmpStr(stat, "low", 2))
            statCanon = "low_nm"
            statLegacy = "low_nm"
        elseif (!CmpStr(stat, "high", 2))
            statCanon = "high_nm"
            statLegacy = "high_nm"
        elseif (!CmpStr(stat, "n", 2))
            statCanon = "n"
            statLegacy = "n"
        elseif (!CmpStr(stat, "max", 2))
            statCanon = "max_nm"
            statLegacy = "max_nm"
        elseif (!CmpStr(stat, "median", 2))
            statCanon = "median_nm"
            statLegacy = "median_nm"
        elseif (!CmpStr(stat, "sigma", 2))
            statCanon = "sigma_nm"
            statLegacy = "sigma_nm"
        elseif (!CmpStr(stat, "max_trace", 2))
            statCanon = "max_trace"
            statLegacy = "max_trace"
        elseif (!CmpStr(stat, "median_trace", 2))
            statCanon = "median_trace"
            statLegacy = "median_trace"
        elseif (!CmpStr(stat, "start", 2))
            statCanon = "start"
            statLegacy = "cluster_start"
        elseif (!CmpStr(stat, "size", 2))
            statCanon = "size"
            statLegacy = "cluster_size"
        elseif (!CmpStr(stat, "minSize", 2))
            statCanon = "minSize"
            statLegacy = "minClusterSize"
        endif

        if (strlen(statCanon) > 0)
            canonical = "v_" + qCanon + "_minbinCluster_" + statCanon
            legacy = "v_" + qLegacy + "_minbin_" + statLegacy
        endif

    // -------------------------------------------------------------------------
    // Histogram-display cluster diagnostics.
    // Legacy v_HistCluster_* names are quantity-ambiguous and overwritten per
    // call; they are supported only as an immediate-call compatibility fallback.
    // -------------------------------------------------------------------------
    elseif (!CmpStr(selection, "histCluster", 2))
        if (!CmpStr(stat, "tailLow", 2) || !CmpStr(stat, "low", 2))
            canonical = "v_" + qCanon + "_histCluster_tailLow_nm"
            legacy = "v_HistCluster_tailLow_nm"
        elseif (!CmpStr(stat, "tailHigh", 2) || !CmpStr(stat, "high", 2))
            canonical = "v_" + qCanon + "_histCluster_tailHigh_nm"
            legacy = "v_HistCluster_tailHigh_nm"
        elseif (!CmpStr(stat, "nTailMembers", 2))
            canonical = "v_" + qCanon + "_histCluster_nTailMembers"
            legacy = "v_HistCluster_nTailMembers"
        elseif (!CmpStr(stat, "peakBin", 2))
            canonical = "v_" + qCanon + "_histCluster_peakBin"
            legacy = "v_HistCluster_peakBin"
        elseif (!CmpStr(stat, "peakBinHeight", 2))
            canonical = "v_" + qCanon + "_histCluster_peakBinHeight"
            legacy = "v_HistCluster_peakBinHeight"
        elseif (!CmpStr(stat, "binLow", 2))
            canonical = "v_" + qCanon + "_histCluster_binLow_nm"
            legacy = "v_HistCluster_binLow_nm"
        elseif (!CmpStr(stat, "binHigh", 2))
            canonical = "v_" + qCanon + "_histCluster_binHigh_nm"
            legacy = "v_HistCluster_binHigh_nm"
        elseif (!CmpStr(stat, "binWidth", 2))
            canonical = "v_" + qCanon + "_histCluster_binWidth_nm"
            legacy = "v_HistCluster_binWidth_nm"
        elseif (!CmpStr(stat, "max", 2))
            canonical = "v_" + qCanon + "_histCluster_max_nm"
            legacy = "v_HistCluster_clusterMax_nm"
        elseif (!CmpStr(stat, "max_trace", 2))
            canonical = "v_" + qCanon + "_histCluster_max_trace"
            legacy = "v_HistCluster_clusterMaxTrace"
        elseif (!CmpStr(stat, "nSelected", 2))
            canonical = "v_" + qCanon + "_histCluster_nSelected"
            legacy = ""
        elseif (!CmpStr(stat, "LclusterForW", 2))
            canonical = "v_" + qCanon + "_histCluster_LclusterForW_nm"
            legacy = "v_HistCluster_LclusterForW_nm"
        elseif (!CmpStr(stat, "SphiForW", 2))
            canonical = "v_" + qCanon + "_histCluster_SphiForW_meVPerRad"
            legacy = "v_HistCluster_SphiForW_meVPerRad"
        endif
    endif

    if (strlen(canonical) == 0)
        return defaultVal
    endif

    String histBase = histDF
    if (strlen(histBase) > 0 && strsearch(histBase, ":", strlen(histBase)-1) < 0)
        histBase += ":"
    endif

    String canonicalPath = histBase + canonical
    String legacyPath = ""
    Variable val = SNS_NumVarByPathOrDefault(canonicalPath, NaN)
    if (numtype(val) == 0)
        return val
    endif

    if (strlen(legacy) > 0)
        legacyPath = histBase + legacy
        val = SNS_NumVarByPathOrDefault(legacyPath, NaN)
        if (numtype(val) == 0)
            SNS_RayDiagLegacyHit(canonicalPath, legacyPath)
            return val
        endif
    endif

    // Convenience fallback: allow callers to pass the channel folder itself
    // instead of its RayTraceHist child. This still does not cross semantic
    // boundaries; it only redirects the same requested diagnostic to the folder
    // where SNS_RayTrace_Hist_LN_Weff(...) stores it.
    if (strsearch(histBase, ":RayTraceHist:", 0, 2) < 0)
        canonicalPath = histBase + "RayTraceHist:" + canonical
        val = SNS_NumVarByPathOrDefault(canonicalPath, NaN)
        if (numtype(val) == 0)
            return val
        endif

        if (strlen(legacy) > 0)
            legacyPath = histBase + "RayTraceHist:" + legacy
            val = SNS_NumVarByPathOrDefault(legacyPath, NaN)
            if (numtype(val) == 0)
                SNS_RayDiagLegacyHit(canonicalPath, legacyPath)
                return val
            endif
        endif
    endif

    return defaultVal
End

Static Function SNS_RayDiagLegacyHit(canonicalPath, legacyPath)
    String canonicalPath, legacyPath

    NewDataFolder/O root:SNS_Log

    NVAR/Z nLegacy = root:SNS_Log:v_RayDiagLegacyFallbackCount
    if (!NVAR_Exists(nLegacy))
        Variable/G root:SNS_Log:v_RayDiagLegacyFallbackCount = 0
    endif
    NVAR nLegacyCount = root:SNS_Log:v_RayDiagLegacyFallbackCount

    String/G root:SNS_Log:s_RayDiagLegacyLastCanonical = canonicalPath
    String/G root:SNS_Log:s_RayDiagLegacyLastLegacy = legacyPath

    nLegacyCount += 1
    if (nLegacyCount == 1)
        SNS_Log("SNS_RayDiagValue: legacy ray diagnostic fallback used. Further hits are counted in root:SNS_Log:v_RayDiagLegacyFallbackCount. first legacy=" + legacyPath + " ; expected canonical=" + canonicalPath, level="WARN")
        // If detailed legacy tracking becomes useful, replace the scalar
        // counter above with a small text/count wave keyed by legacyPath.
    endif

    return 0
End


//==============================================================================
// SNS_DisplayHistClusterFamilyRays
//
// Wrapper around SNS_RayTrace_Hist_LN_Weff(...).
//
// Selection logic:
//   1. Run histogram using supplied binL_nm / binW_nm.
//   2. Apply the metadata-defined upper-tail bounds from
//      v_*_extreme_low/high.
//   3. Inside that upper-tail interval, find the most populated supplied bin.
//   4. Define that bin as the max cluster.
//   5. Draw exactly the rays inside that selected max-cluster bin.
//   6. Solid black ray = maximum-value member inside that selected bin.
//
// Notes:
//   - useExpLBin/useExpWBin are ignored; kept only for call compatibility.
//   - useExistingChannels=1 skips SNS_ExtractModes* and histograms the
//     canonical channel waves already present in dfPath.
//   - useFanDisplay=0: draw selected rays as graph traces.
//   - useFanDisplay=1: draw fan fill, black max ray, and STS marker as draw objects.
//   - colorPreset overrides familyColor:
//       0 = black
//       1 = Nature blue #0072B2 -> Igor (0,29298,45746)
//       3 = Nature vermillion #D55E00 -> Igor (54741,24158,0)
//   - For W/W_geom with absW=1, all binning/selection uses abs(W).
//==============================================================================
Function SNS_DisplayHistClusterFamilyRays(dfPath, [quantity, Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy, H_nm, maxPath_nm, nBins, binL_nm, binW_nm, histLMin_nm, histLMax_nm, histWMin_nm, histWMax_nm, absW, normalize, useWeights, useExpLBin, useExpWBin, EresHWHM_meV, familyColor, colorPreset, lineSize, maxLineSize, alphaMin, alphaMax, alphaGamma, alphaBeta, raySparsity, rayStartIndex, displayXAxisAngle_deg, useFanDisplay, fanNSteps, useExistingChannels])
	String dfPath
	String quantity
	Variable Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy
	Variable H_nm, maxPath_nm
	Variable nBins, binL_nm, binW_nm, histLMin_nm, histLMax_nm, histWMin_nm, histWMax_nm
	Variable absW, normalize, useWeights
	Variable useExpLBin, useExpWBin, EresHWHM_meV
	String familyColor
	Variable colorPreset
	Variable lineSize, maxLineSize
	Variable alphaMin, alphaMax, alphaGamma, alphaBeta
	Variable raySparsity, rayStartIndex
	Variable displayXAxisAngle_deg
	Variable useFanDisplay, fanNSteps
	Variable useExistingChannels

	String oldDF = GetDataFolder(1)

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	// ---------------- Defaults ----------------
	if (ParamIsDefault(quantity))
		quantity = "L"
	endif
	if (ParamIsDefault(Bangle_deg))
		Bangle_deg = 225
	endif
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
	if (ParamIsDefault(H_nm))
		H_nm = 0
	endif
	if (ParamIsDefault(maxPath_nm))
		maxPath_nm = 1000
	endif
	if (ParamIsDefault(nBins))
		nBins = 50
	endif
	if (ParamIsDefault(absW))
		absW = 1
	endif
	if (ParamIsDefault(normalize))
		normalize = 1
	endif
	if (ParamIsDefault(useWeights))
		useWeights = (H_nm > 0)
	endif
	if (ParamIsDefault(useExpLBin))
		useExpLBin = 0
	endif
	if (ParamIsDefault(useExpWBin))
		useExpWBin = 0
	endif
	if (ParamIsDefault(EresHWHM_meV))
		EresHWHM_meV = 0.136
	endif
	if (ParamIsDefault(lineSize))
		lineSize = 1.0
	endif
	if (ParamIsDefault(maxLineSize))
		maxLineSize = 2.5
	endif
	if (ParamIsDefault(alphaMin))
		alphaMin = 0.04
	endif
	if (ParamIsDefault(alphaMax))
		alphaMax = 0.75
	endif
	if (ParamIsDefault(alphaGamma))
		alphaGamma = 1.0
	endif
	if (ParamIsDefault(alphaBeta))
		alphaBeta = 8.0
	endif
	if (ParamIsDefault(raySparsity))
		raySparsity = 1
	endif
	if (ParamIsDefault(rayStartIndex))
		rayStartIndex = 0
	endif
	if (ParamIsDefault(useFanDisplay))
		useFanDisplay = 0
	endif
	if (ParamIsDefault(fanNSteps))
		fanNSteps = 48
	endif
	if (ParamIsDefault(useExistingChannels))
		useExistingChannels = 0
	endif
	if (ParamIsDefault(colorPreset))
		colorPreset = NaN
	else
		colorPreset = round(colorPreset)
	endif

	if (ParamIsDefault(histLMin_nm))
		histLMin_nm = NaN
	endif
	if (ParamIsDefault(histLMax_nm))
		histLMax_nm = NaN
	endif
	if (ParamIsDefault(histWMin_nm))
		histWMin_nm = NaN
	endif
	if (ParamIsDefault(histWMax_nm))
		histWMax_nm = NaN
	endif

	Variable rotateDisplay = 0
	if (!ParamIsDefault(displayXAxisAngle_deg) && numtype(displayXAxisAngle_deg) == 0)
		rotateDisplay = 1
	endif

	alphaMin = max(0, min(1, alphaMin))
	alphaMax = max(0, min(1, alphaMax))
	alphaGamma = max(0.1, alphaGamma)
	alphaBeta = max(0.1, alphaBeta)
	raySparsity = max(1, round(raySparsity))
	rayStartIndex = max(0, round(rayStartIndex))
	EresHWHM_meV = max(1e-6, EresHWHM_meV)
	useFanDisplay = round(useFanDisplay)
	useFanDisplay = max(0, min(1, useFanDisplay))
	fanNSteps = max(8, round(fanNSteps))
	useExistingChannels = round(useExistingChannels)
	useExistingChannels = max(0, min(1, useExistingChannels))

	Variable use3D = (H_nm > 0)

	Variable isL = StringMatch(quantity, "L") || StringMatch(quantity, "L_N") || StringMatch(quantity, "ell")
	Variable isW = StringMatch(quantity, "W") || StringMatch(quantity, "W_eff") || StringMatch(quantity, "w")
	Variable isWGeom = StringMatch(quantity, "W_geom") || StringMatch(quantity, "Wgeom") || StringMatch(quantity, "w_geom")

	if (!isL && !isW && !isWGeom)
		SetDataFolder $oldDF
		Abort "SNS_DisplayHistClusterFamilyRays: quantity must be \"L\", \"W\", or \"W_geom\"."
	endif
	if (isWGeom && !use3D)
		SetDataFolder $oldDF
		Abort "SNS_DisplayHistClusterFamilyRays: W_geom requires H_nm > 0 / 3D mode."
	endif

	if (ParamIsDefault(familyColor))
		if (isL)
			familyColor = "red"
		elseif (isW)
			familyColor = "cyan"
		else
			familyColor = "blue"
		endif
	endif

	// -------------------------------------------------------------------------
	// 1. Run histogram function.
	// -------------------------------------------------------------------------
	Variable passBinL = !ParamIsDefault(binL_nm)
	Variable passBinW = !ParamIsDefault(binW_nm)

	if (passBinL && passBinW)
		SNS_RayTrace_Hist_LN_Weff(dfPath, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, Vortexx=Vortexx, Vortexy=Vortexy, H_nm=H_nm, maxPath_nm=maxPath_nm, nBins=nBins, binL_nm=binL_nm, binW_nm=binW_nm, histLMin_nm=histLMin_nm, histLMax_nm=histLMax_nm, histWMin_nm=histWMin_nm, histWMax_nm=histWMax_nm, absW=absW, normalize=normalize, useWeights=useWeights, doDisplay=0, useExistingChannels=useExistingChannels)
	elseif (passBinL)
		SNS_RayTrace_Hist_LN_Weff(dfPath, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, Vortexx=Vortexx, Vortexy=Vortexy, H_nm=H_nm, maxPath_nm=maxPath_nm, nBins=nBins, binL_nm=binL_nm, histLMin_nm=histLMin_nm, histLMax_nm=histLMax_nm, histWMin_nm=histWMin_nm, histWMax_nm=histWMax_nm, absW=absW, normalize=normalize, useWeights=useWeights, doDisplay=0, useExistingChannels=useExistingChannels)
	elseif (passBinW)
		SNS_RayTrace_Hist_LN_Weff(dfPath, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, Vortexx=Vortexx, Vortexy=Vortexy, H_nm=H_nm, maxPath_nm=maxPath_nm, nBins=nBins, binW_nm=binW_nm, histLMin_nm=histLMin_nm, histLMax_nm=histLMax_nm, histWMin_nm=histWMin_nm, histWMax_nm=histWMax_nm, absW=absW, normalize=normalize, useWeights=useWeights, doDisplay=0, useExistingChannels=useExistingChannels)
	else
		SNS_RayTrace_Hist_LN_Weff(dfPath, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, Vortexx=Vortexx, Vortexy=Vortexy, H_nm=H_nm, maxPath_nm=maxPath_nm, nBins=nBins, histLMin_nm=histLMin_nm, histLMax_nm=histLMax_nm, histWMin_nm=histWMin_nm, histWMax_nm=histWMax_nm, absW=absW, normalize=normalize, useWeights=useWeights, doDisplay=0, useExistingChannels=useExistingChannels)
	endif

	// -------------------------------------------------------------------------
	// 2. Locate source waves and histogram output folder.
	// -------------------------------------------------------------------------
	SetDataFolder $dfPath

	String imgList = WaveList("*_Z_mbgnd_xy", ";", "DIMS:2")
	String imgName = ""
	if (strlen(imgList) == 0)
		Wave/Z maskImageFallback = w_mask
		if (WaveExists(maskImageFallback))
			imgName = "w_mask"
		else
			SetDataFolder $oldDF
			Abort "SNS_DisplayHistClusterFamilyRays: no *_Z_mbgnd_xy image or w_mask found."
		endif
	else
		imgName = StringFromList(0, imgList, ";")
	endif

	Wave image = $imgName

	String outDF
	if (use3D)
		outDF = dfPath + ":RayTraceHist_3D"
	else
		outDF = dfPath + ":RayTraceHist"
	endif

	if (use3D)
		Wave/Z Hit1x_List = $(outDF + ":Hit1x_List_3D_nm")
		Wave/Z Hit1y_List = $(outDF + ":Hit1y_List_3D_nm")
		Wave/Z Hit2x_List = $(outDF + ":Hit2x_List_3D_nm")
		Wave/Z Hit2y_List = $(outDF + ":Hit2y_List_3D_nm")
	else
		Wave/Z Hit1x_List = $(outDF + ":Hit1x_List_nm")
		Wave/Z Hit1y_List = $(outDF + ":Hit1y_List_nm")
		Wave/Z Hit2x_List = $(outDF + ":Hit2x_List_nm")
		Wave/Z Hit2y_List = $(outDF + ":Hit2y_List_nm")
	endif

	if (!WaveExists(Hit1x_List) || !WaveExists(Hit1y_List) || !WaveExists(Hit2x_List) || !WaveExists(Hit2y_List))
		SetDataFolder $oldDF
		Abort "SNS_DisplayHistClusterFamilyRays: missing histogram Hit* ray-coordinate waves."
	endif

	String valuePath, LPath
	String labelText
	String diagQuantity
	Variable histMin_nm, histMax_nm, dHist_nm, nRay
	Variable tailLow_nm, tailHigh_nm

	if (isL)
		diagQuantity = "L"
		if (use3D)
			valuePath = outDF + ":L_N_List_3D_nm"
			LPath = outDF + ":L_N_List_3D_nm"
			histMin_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_3D_LminHist_nm", NaN)
			histMax_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_3D_LmaxHist_nm", NaN)
			dHist_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_3D_dLHist_nm", NaN)
		else
			valuePath = outDF + ":L_N_List_nm"
			LPath = outDF + ":L_N_List_nm"
			histMin_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_LminHist_nm", NaN)
			histMax_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_LmaxHist_nm", NaN)
			dHist_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_dLHist_nm", NaN)
		endif
		labelText = "ℓ"
	elseif (isW)
		diagQuantity = "W"
		if (use3D)
			valuePath = outDF + ":W_eff_List_3D_nm"
			LPath = outDF + ":L_N_List_3D_nm"
			histMin_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_3D_WminHist_nm", NaN)
			histMax_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_3D_WmaxHist_nm", NaN)
			dHist_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_3D_dWHist_nm", NaN)
		else
			valuePath = outDF + ":W_eff_List_nm"
			LPath = outDF + ":L_N_List_nm"
			histMin_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_WminHist_nm", NaN)
			histMax_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_WmaxHist_nm", NaN)
			dHist_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_dWHist_nm", NaN)
		endif
		labelText = "w"
	else
		diagQuantity = "Wgeom"
		valuePath = outDF + ":W_geom_List_3D_nm"
		LPath = outDF + ":L_N_List_3D_nm"
		histMin_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_3D_WminHist_nm", NaN)
		histMax_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_3D_WmaxHist_nm", NaN)
		dHist_nm = NumVarOrDefault(outDF + ":v_RayTraceHist_3D_dWHist_nm", NaN)
		labelText = "w\\Bgeom\\M"
	endif

	tailLow_nm = SNS_RayDiagValue(outDF + ":", diagQuantity, "tailCluster", "low", is3D=use3D)
	tailHigh_nm = SNS_RayDiagValue(outDF + ":", diagQuantity, "tailCluster", "high", is3D=use3D)

	Wave/Z value_nm = $valuePath
	Wave/Z L_for_slope_nm = $LPath

	if (!WaveExists(value_nm))
		SetDataFolder $oldDF
		Abort "SNS_DisplayHistClusterFamilyRays: missing value wave in histogram output folder."
	endif
	if (!WaveExists(L_for_slope_nm))
		SetDataFolder $oldDF
		Abort "SNS_DisplayHistClusterFamilyRays: missing L_N wave for W-resolution estimate."
	endif

	nRay = numpnts(value_nm)

	if (numtype(histMin_nm) != 0 || numtype(histMax_nm) != 0 || numtype(dHist_nm) != 0 || dHist_nm <= 0)
		SetDataFolder $oldDF
		Abort "SNS_DisplayHistClusterFamilyRays: invalid histogram metadata."
	endif
	if (numtype(tailLow_nm) != 0 || numtype(tailHigh_nm) != 0 || tailHigh_nm <= tailLow_nm)
		SetDataFolder $oldDF
		Abort "SNS_DisplayHistClusterFamilyRays: missing upper-tail cluster low/high metadata."
	endif

	// -------------------------------------------------------------------------
	// 3. Inside the metadata-defined upper-tail interval, find the most
	//    populated supplied bin.
	// -------------------------------------------------------------------------
	Variable nBinsUse = max(1, round((histMax_nm - histMin_nm) / dHist_nm))
	Make/FREE/D/N=(nBinsUse) tmpBinScore
	tmpBinScore = 0

	Variable i, k, xVal, iBinTmp
	Variable nTailMembers = 0

	for (i = 0; i < nRay; i += 1)
		xVal = value_nm[i]
		if ((isW || isWGeom) && absW)
			xVal = abs(xVal)
		endif

		if (numtype(xVal) != 0)
			continue
		endif

		if (xVal < tailLow_nm || xVal > tailHigh_nm)
			continue
		endif

		iBinTmp = floor((xVal - histMin_nm) / dHist_nm)
		iBinTmp = max(0, min(nBinsUse - 1, iBinTmp))

		tmpBinScore[iBinTmp] += 1
		nTailMembers += 1
	endfor

	WaveStats/Q tmpBinScore
	Variable peakHeight = V_max

	if (numtype(peakHeight) != 0 || peakHeight <= 0)
		Print "SNS_DisplayHistClusterFamilyRays debug:"
		Print "  histMin/max/d = ", histMin_nm, histMax_nm, dHist_nm
		Print "  tailLow/high  = ", tailLow_nm, tailHigh_nm
		Print "  nTailMembers  = ", nTailMembers
		SetDataFolder $oldDF
		Abort "SNS_DisplayHistClusterFamilyRays: no populated supplied bin inside the upper-tail interval."
	endif

	Variable peakBin = 0
	Variable bestBinCenter_nm = -1e300
	Variable binCenter_nm

	for (k = 0; k < nBinsUse; k += 1)
		if (tmpBinScore[k] >= peakHeight - 1e-12)
			binCenter_nm = histMin_nm + (k + 0.5)*dHist_nm
			if (binCenter_nm > bestBinCenter_nm)
				bestBinCenter_nm = binCenter_nm
				peakBin = k
			endif
		endif
	endfor

	Variable clusterLow_nm  = histMin_nm + peakBin*dHist_nm
	Variable clusterHigh_nm = clusterLow_nm + dHist_nm
	Variable clusterWidth_nm = clusterHigh_nm - clusterLow_nm

	Variable clusterMax_nm = -1e300
	Variable maxTrace = NaN

	for (i = 0; i < nRay; i += 1)
		xVal = value_nm[i]
		if ((isW || isWGeom) && absW)
			xVal = abs(xVal)
		endif

		if (numtype(xVal) != 0)
			continue
		endif

		if (xVal < tailLow_nm || xVal > tailHigh_nm)
			continue
		endif

		if (xVal >= clusterLow_nm && xVal <= clusterHigh_nm)
			if (xVal > clusterMax_nm)
				clusterMax_nm = xVal
				maxTrace = i
			endif
		endif
	endfor

	if (numtype(maxTrace) != 0 || maxTrace < 0 || maxTrace >= nRay)
		SetDataFolder $oldDF
		Abort "SNS_DisplayHistClusterFamilyRays: no max member found inside selected supplied-bin cluster."
	endif

	// -------------------------------------------------------------------------
	// 4. Display membership is always selected max-cluster bin.
	// -------------------------------------------------------------------------
	Variable binLow_nm = clusterLow_nm
	Variable binHigh_nm = clusterHigh_nm
	Variable binWidth_nm = clusterWidth_nm

	Variable Lcluster_forW_nm = NaN
	Variable Sphi_cluster_meVPerRad = NaN

	if (isW || isWGeom)
		Lcluster_forW_nm = L_for_slope_nm[maxTrace]
		Sphi_cluster_meVPerRad = SNS_PhaseSlope_meVPerRad_FromParams(Lcluster_forW_nm)
	endif

	Variable maxDist_nm = max(abs(clusterMax_nm - binLow_nm), abs(binHigh_nm - clusterMax_nm))
	if (maxDist_nm <= 0 || numtype(maxDist_nm) != 0)
		maxDist_nm = binWidth_nm
	endif
	if (maxDist_nm <= 0 || numtype(maxDist_nm) != 0)
		maxDist_nm = 1
	endif

	// -------------------------------------------------------------------------
	// 5. Build mask/alpha waves in histogram output folder.
	// -------------------------------------------------------------------------
	SetDataFolder $outDF

	Make/O/D/N=(nRay) HistCluster_family_mask
	Make/O/D/N=(nRay) HistCluster_family_alpha

	Wave mask = HistCluster_family_mask
	Wave alpha = HistCluster_family_alpha

	mask = 0
	alpha = 0

	Variable deficit, expMin
	expMin = exp(-alphaBeta)

	Variable nFamilySelected = 0
	Variable nFamilyPlotted = 0

	for (i = 0; i < nRay; i += 1)
		xVal = value_nm[i]
		if ((isW || isWGeom) && absW)
			xVal = abs(xVal)
		endif

		if (numtype(xVal) != 0)
			continue
		endif

		if (xVal >= clusterLow_nm && xVal <= clusterHigh_nm)
			mask[i] = 1
			nFamilySelected += 1

			deficit = abs(xVal - clusterMax_nm) / maxDist_nm
			deficit = max(0, min(1, deficit))

			alpha[i] = alphaMin + (alphaMax - alphaMin) * (exp(-alphaBeta * deficit^alphaGamma) - expMin) / (1 - expMin)
			alpha[i] = max(alphaMin, min(alphaMax, alpha[i]))
		endif
	endfor

	Variable/G v_HistCluster_tailLow_nm = tailLow_nm
	Variable/G v_HistCluster_tailHigh_nm = tailHigh_nm
	Variable/G v_HistCluster_nTailMembers = nTailMembers
	Variable/G v_HistCluster_peakBin = peakBin
	Variable/G v_HistCluster_peakBinHeight = peakHeight
	Variable/G v_HistCluster_clusterLow_nm = clusterLow_nm
	Variable/G v_HistCluster_clusterHigh_nm = clusterHigh_nm
	Variable/G v_HistCluster_binLow_nm = binLow_nm
	Variable/G v_HistCluster_binHigh_nm = binHigh_nm
	Variable/G v_HistCluster_binWidth_nm = binWidth_nm
	Variable/G v_HistCluster_clusterMax_nm = clusterMax_nm
	Variable/G v_HistCluster_clusterMaxTrace = maxTrace
	Variable/G v_HistCluster_absW = absW
	Variable/G v_HistCluster_LclusterForW_nm = Lcluster_forW_nm
	Variable/G v_HistCluster_SphiForW_meVPerRad = Sphi_cluster_meVPerRad
	Variable/G v_HistCluster_displayXAxisAngle_deg = displayXAxisAngle_deg
	Variable/G v_HistCluster_rotateDisplay = rotateDisplay
	Variable/G v_HistCluster_useFanDisplay = useFanDisplay
	Variable/G v_HistCluster_colorPreset = colorPreset

	if (!useFanDisplay)
		Variable nFamilyToDraw = 0

		for (i = 0; i < nRay; i += 1)
			if (mask[i] < 0.5)
				continue
			endif
			if (i < rayStartIndex)
				continue
			endif
			if (mod(i - rayStartIndex, raySparsity) != 0)
				continue
			endif

			nFamilyToDraw += 1
		endfor

		// The selected cluster-max ray is drawn once more on top.
		nFamilyToDraw += 1

		if (SNS_WarnLargeRayDisplay("SNS_DisplayHistClusterFamilyRays", nFamilySelected, nFamilyToDraw, raySparsity))
			SetDataFolder $oldDF
			return 0
		endif
	endif

	SetDataFolder $dfPath

	// -------------------------------------------------------------------------
	// 6. Display topography, optionally rotated.
	// -------------------------------------------------------------------------
	Variable xOrig0 = DimOffset(image, 0)
	Variable xOrig1 = xOrig0 + DimDelta(image, 0) * (DimSize(image, 0) - 1)
	Variable yOrig0 = DimOffset(image, 1)
	Variable yOrig1 = yOrig0 + DimDelta(image, 1) * (DimSize(image, 1) - 1)

	Variable xCenterRot = 0.5 * (min(xOrig0, xOrig1) + max(xOrig0, xOrig1))
	Variable yCenterRot = 0.5 * (min(yOrig0, yOrig1) + max(yOrig0, yOrig1))

	String displayImgName = NameOfWave(image)

	if (rotateDisplay)
		SNS_MakeRotatedImageForDisplay(image, "tmp_HistCluster_rot_image", displayXAxisAngle_deg)
		displayImgName = "tmp_HistCluster_rot_image"
	endif

	Wave displayImage = $displayImgName

	SNS_DisplayWithScales(displayImage, cmap="grayC")
	String winImage = WinName(0, 1, 1)

	SNS_StyleRayDisplayImage(displayImage, winImage)

	DoWindow/F $winImage

	// -------------------------------------------------------------------------
	// 7. Plot selected max-cluster-bin rays / fan.
	// -------------------------------------------------------------------------
	Variable rCol, gCol, bCol
	Variable fanAlphaRGBA = round(65535 * alphaMax)

	// Optional color presets:
	//   0 = black
	//   1 = Nature blue #0072B2
	//   3 = Nature vermillion #D55E00
	if (numtype(colorPreset) == 0)
		if (colorPreset == 0)
			rCol = 0
			gCol = 0
			bCol = 0
			fanAlphaRGBA = 49151
		elseif (colorPreset == 1)
			rCol = 0
			gCol = 29298
			bCol = 45746
			fanAlphaRGBA = 49151
		elseif (colorPreset == 3)
			rCol = 54741
			gCol = 24158
			bCol = 0
			fanAlphaRGBA = 49151
		else
			rCol = 65535
			gCol = 0
			bCol = 0
		endif
	else
		if (StringMatch(familyColor, "black"))
			rCol = 0; gCol = 0; bCol = 0
		elseif (StringMatch(familyColor, "white"))
			rCol = 65535; gCol = 65535; bCol = 65535
		elseif (StringMatch(familyColor, "cyan"))
			rCol = 1; gCol = 16019; bCol = 65535
		elseif (StringMatch(familyColor, "blue"))
			rCol = 0; gCol = 0; bCol = 65535
		elseif (StringMatch(familyColor, "green"))
			rCol = 0; gCol = 45000; bCol = 0
		elseif (StringMatch(familyColor, "magenta"))
			rCol = 65535; gCol = 0; bCol = 65535
		elseif (StringMatch(familyColor, "yellow"))
			rCol = 65535; gCol = 50000; bCol = 0
		else
			rCol = 65535; gCol = 0; bCol = 0
		endif
	endif

	String traceList, lastTrace
	Variable aRGBA

	Variable thRot = displayXAxisAngle_deg*pi/180
	Variable stsX_disp, stsY_disp

	if (rotateDisplay)
		stsX_disp = (STSx - xCenterRot)*cos(thRot) + (STSy - yCenterRot)*sin(thRot)
		stsY_disp = -(STSx - xCenterRot)*sin(thRot) + (STSy - yCenterRot)*cos(thRot)
	else
		stsX_disp = STSx
		stsY_disp = STSy
	endif

	if (useFanDisplay)

		SNS__DrawHistClusterEndpointWalkFan(winImage, Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, mask, maxTrace, STSx, STSy, displayXAxisAngle_deg, rotateDisplay, xCenterRot, yCenterRot, rCol, gCol, bCol, alphaMax, fanAlphaRGBA, lineSize*1.4, maxLineSize, 2.5)

		nFamilyPlotted = 1

	else

		for (i = 0; i < nRay; i += 1)
			if (mask[i] < 0.5)
				continue
			endif

			if (i < rayStartIndex)
				continue
			endif
			if (mod(i - rayStartIndex, raySparsity) != 0)
				continue
			endif

			DoWindow/F $winImage

			if (rotateDisplay)
				SNS_PlotChannelRay_Rotated(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, i, displayXAxisAngle_deg, xCenterRot, yCenterRot)
			else
				SNS_PlotChannelRay(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, i)
			endif

			traceList = TraceNameList(winImage, ";", 1)
			lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

			aRGBA = round(65535 * alpha[i])
			aRGBA = max(0, min(65535, aRGBA))

			ModifyGraph/W=$winImage rgb($lastTrace)=(rCol,gCol,bCol,aRGBA)
			ModifyGraph/W=$winImage lsize($lastTrace)=lineSize

			nFamilyPlotted += 1
		endfor

		// Solid cluster-max ray on top as graph trace.
		DoWindow/F $winImage

		if (rotateDisplay)
			SNS_PlotChannelRay_Rotated(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, maxTrace, displayXAxisAngle_deg, xCenterRot, yCenterRot)
		else
			SNS_PlotChannelRay(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, maxTrace)
		endif

		traceList = TraceNameList(winImage, ";", 1)
		lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

		ModifyGraph/W=$winImage rgb($lastTrace)=(0,0,0,65535)
		ModifyGraph/W=$winImage lsize($lastTrace)=maxLineSize

	endif

	// -------------------------------------------------------------------------
	// 8. STS marker and B arrow.
	// -------------------------------------------------------------------------
	if (!useFanDisplay)
		Make/O/N=1 tmpSTSx, tmpSTSy
		tmpSTSx[0] = stsX_disp
		tmpSTSy[0] = stsY_disp

		AppendToGraph/W=$winImage tmpSTSy vs tmpSTSx

		traceList = TraceNameList(winImage, ";", 1)
		lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

		ModifyGraph/W=$winImage mode($lastTrace)=3
		ModifyGraph/W=$winImage marker($lastTrace)=19
		ModifyGraph/W=$winImage msize($lastTrace)=5
		ModifyGraph/W=$winImage rgb($lastTrace)=(65535,65535,65535,65535)

		DoWindow/F $winImage
		SNS_TAGSetpoint()
	endif

	Variable xMin = DimOffset(displayImage, 0)
	Variable xMax = xMin + DimDelta(displayImage, 0) * (DimSize(displayImage, 0) - 1)
	Variable yMin = DimOffset(displayImage, 1)
	Variable yMax = yMin + DimDelta(displayImage, 1) * (DimSize(displayImage, 1) - 1)

	Variable pos_x = xMax - 0.15 * (xMax - xMin)
	Variable pos_y = yMin + 0.15 * (yMax - yMin)

	Variable displayBangle_deg = Bangle_deg
	if (rotateDisplay)
		displayBangle_deg = Bangle_deg - displayXAxisAngle_deg
	endif

	SNS_DrawImageArrow(displayImage, displayBangle_deg, "B", pos_x=pos_x, pos_y=pos_y, scale=1.5, color="black")

	// -------------------------------------------------------------------------
	// 9. Textbox and window scale.
	// -------------------------------------------------------------------------
	String txt
	txt = "Histogram-cluster display\r"
	txt += labelText + "\\Bcluster\\M = " + num2str(round(10 * clusterMax_nm) / 10) + " nm\r"
	txt += "tail = [" + num2str(round(10*tailLow_nm)/10) + ", " + num2str(round(10*tailHigh_nm)/10) + "] nm\r"
	txt += "display bin = [" + num2str(round(10*clusterLow_nm)/10) + ", " + num2str(round(10*clusterHigh_nm)/10) + "] nm\r"
	txt += "peak bin = " + num2str(peakBin) + ", N = " + num2str(round(10*peakHeight)/10) + "\r"
	txt += "Δ" + labelText + "\\Bbin\\M = " + num2str(round(10 * binWidth_nm) / 10) + " nm\r"
	txt += "max trace = " + num2str(maxTrace) + "\r"
	txt += "N\\Btail\\M = " + num2str(nTailMembers) + " / " + num2str(nRay) + "\r"
	txt += "N\\Bbin\\M = " + num2str(nFamilySelected) + " / " + num2str(nRay) + "\r"
	txt += "N\\Bshown\\M = " + num2str(nFamilyPlotted) + "\r"
	txt += "alpha β = " + num2str(round(10 * alphaBeta) / 10) + ", γ = " + num2str(round(10 * alphaGamma) / 10)

	if ((isW || isWGeom) && numtype(Sphi_cluster_meVPerRad) == 0)
		txt += "\rS\\Bφ\\M = " + num2str(round(1000 * Sphi_cluster_meVPerRad) / 1000) + " meV/rad"
	endif
	if (rotateDisplay)
		txt += "\rdisplay x-axis = " + num2str(round(10 * displayXAxisAngle_deg) / 10) + "°"
	endif
	if (useFanDisplay)
		txt += "\rfan display = 1"
	endif
	if (numtype(colorPreset) == 0)
		txt += "\rcolor preset = " + num2str(colorPreset)
	endif

	TextBox/C/N=HistClusterInfo/F=0/A=RT/X=0/Y=0/B=(65535,65535,65535,39321) txt

	GetWindow $winImage, wsize
	Variable wPx = abs(V_right - V_left)
	Variable hPx = abs(V_bottom - V_top)
	Variable maxDim = max(wPx, hPx)

	if (maxDim > 512)
		Variable ex = 512 / maxDim
		ModifyGraph/W=$winImage expand=ex
	else
		ModifyGraph/W=$winImage expand=1
	endif

	SetDataFolder $oldDF

	return maxTrace
End
//==============================================================================
// SNS__DrawHistClusterEndpointWalkFan
//
// Fan/envelope display.
//
// This version sorts selected rays by a 180-degree folded midpoint angle.
// This treats θ and θ+π as the same ray direction, which is appropriate for
// chord-like SNS trajectories.
//
// Steps:
//   1. collect selected rays
//   2. rotate endpoints to display coordinates
//   3. compute midpoint angle around STS
//   4. fold angle into [0, pi)
//   5. sort by folded angle
//   6. cut at largest folded-angle gap
//   7. assign A/B endpoint branches continuously
//   8. draw two STS-anchored polygons
//   9. draw black max ray and STS marker last
//
// IMPORTANT:
//   Pass unrotated STSx, STSy.
//==============================================================================
Function SNS__DrawHistClusterEndpointWalkFan(winName, Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, mask, maxTrace, STSx, STSy, displayXAxisAngle_deg, rotateDisplay, xCenterRot, yCenterRot, rCol, gCol, bCol, alphaFill, fanAlphaRGBA, edgeLineSize, maxLineSize, stsMarkerRadius_nm)
	String winName
	WAVE Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
	WAVE mask
	Variable maxTrace
	Variable STSx, STSy
	Variable displayXAxisAngle_deg, rotateDisplay, xCenterRot, yCenterRot
	Variable rCol, gCol, bCol
	Variable alphaFill, fanAlphaRGBA, edgeLineSize, maxLineSize, stsMarkerRadius_nm

	Variable nRay = numpnts(mask)
	Variable nSel = 0
	Variable i

	for (i = 0; i < nRay; i += 1)
		if (mask[i] >= 0.5)
			nSel += 1
		endif
	endfor

	if (nSel < 2)
		return 0
	endif

	Variable thRot = displayXAxisAngle_deg*pi/180
	Variable stsX_disp, stsY_disp

	if (rotateDisplay)
		stsX_disp = (STSx - xCenterRot)*cos(thRot) + (STSy - yCenterRot)*sin(thRot)
		stsY_disp = -(STSx - xCenterRot)*sin(thRot) + (STSy - yCenterRot)*cos(thRot)
	else
		stsX_disp = STSx
		stsY_disp = STSy
	endif

	// -------------------------------------------------------------------------
	// 1. Collect selected rays and sort key = midpoint angle folded to [0, pi).
	// -------------------------------------------------------------------------
	Make/O/D/N=(nSel) tmp_HistFan_sortAngle
	Make/O/D/N=(nSel) tmp_HistFan_index
	Make/O/D/N=(nSel) tmp_HistFan_H1x
	Make/O/D/N=(nSel) tmp_HistFan_H1y
	Make/O/D/N=(nSel) tmp_HistFan_H2x
	Make/O/D/N=(nSel) tmp_HistFan_H2y

	WAVE sortAngle = tmp_HistFan_sortAngle
	WAVE selIndex = tmp_HistFan_index
	WAVE H1x = tmp_HistFan_H1x
	WAVE H1y = tmp_HistFan_H1y
	WAVE H2x = tmp_HistFan_H2x
	WAVE H2y = tmp_HistFan_H2y

	Variable c = 0
	Variable h1xD, h1yD, h2xD, h2yD
	Variable mx, my, theta

	for (i = 0; i < nRay; i += 1)
		if (mask[i] < 0.5)
			continue
		endif

		SNS__GetRotatedRayEndpoints(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, i, displayXAxisAngle_deg, rotateDisplay, xCenterRot, yCenterRot, h1xD, h1yD, h2xD, h2yD)

		H1x[c] = h1xD
		H1y[c] = h1yD
		H2x[c] = h2xD
		H2y[c] = h2yD
		selIndex[c] = i

		mx = 0.5*(h1xD + h2xD)
		my = 0.5*(h1yD + h2yD)

		theta = atan2(my - stsY_disp, mx - stsX_disp)

		// Fold polar angle into nematic/chord angle: [0, pi).
		theta = mod(theta, pi)
		if (theta < 0)
			theta += pi
		endif

		sortAngle[c] = theta

		c += 1
	endfor

	if (c < 2)
		return 0
	endif

	Redimension/N=(c) sortAngle, selIndex, H1x, H1y, H2x, H2y
	nSel = c

	Sort sortAngle, sortAngle, selIndex, H1x, H1y, H2x, H2y

	// -------------------------------------------------------------------------
	// 2. Remove branch-cut jump in the folded [0, pi) domain.
	// -------------------------------------------------------------------------
	Variable maxGap = -1
	Variable gap, cutIndex = 0

	for (i = 0; i < nSel-1; i += 1)
		gap = sortAngle[i+1] - sortAngle[i]
		if (gap > maxGap)
			maxGap = gap
			cutIndex = i
		endif
	endfor

	// Folded-angle wrap gap is pi-periodic, not 2*pi-periodic.
	gap = sortAngle[0] + pi - sortAngle[nSel-1]
	if (gap > maxGap)
		maxGap = gap
		cutIndex = nSel - 1
	endif

	Variable startIndex = mod(cutIndex + 1, nSel)

	Make/O/D/N=(nSel) tmp_HistFan_ordIndex
	Make/O/D/N=(nSel) tmp_HistFan_ordH1x
	Make/O/D/N=(nSel) tmp_HistFan_ordH1y
	Make/O/D/N=(nSel) tmp_HistFan_ordH2x
	Make/O/D/N=(nSel) tmp_HistFan_ordH2y
	Make/O/D/N=(nSel) tmp_HistFan_ordAngle

	WAVE ordIndex = tmp_HistFan_ordIndex
	WAVE ordH1x = tmp_HistFan_ordH1x
	WAVE ordH1y = tmp_HistFan_ordH1y
	WAVE ordH2x = tmp_HistFan_ordH2x
	WAVE ordH2y = tmp_HistFan_ordH2y
	WAVE ordAngle = tmp_HistFan_ordAngle

	Variable src

	for (i = 0; i < nSel; i += 1)
		src = mod(startIndex + i, nSel)

		ordIndex[i] = selIndex[src]
		ordH1x[i] = H1x[src]
		ordH1y[i] = H1y[src]
		ordH2x[i] = H2x[src]
		ordH2y[i] = H2y[src]
		ordAngle[i] = sortAngle[src]
	endfor

	// -------------------------------------------------------------------------
	// 3. Assign endpoint branches continuously along ordered rays.
	// -------------------------------------------------------------------------
	Make/O/D/N=(nSel) tmp_HistFan_Ax
	Make/O/D/N=(nSel) tmp_HistFan_Ay
	Make/O/D/N=(nSel) tmp_HistFan_Bx
	Make/O/D/N=(nSel) tmp_HistFan_By

	WAVE Ax = tmp_HistFan_Ax
	WAVE Ay = tmp_HistFan_Ay
	WAVE Bx = tmp_HistFan_Bx
	WAVE By = tmp_HistFan_By

	Variable prevAx, prevAy, prevBx, prevBy
	Variable costSame, costSwap
	Variable firstTrace = ordIndex[0]
	Variable lastTrace = ordIndex[nSel-1]

	for (i = 0; i < nSel; i += 1)

		if (i == 0)
			Ax[i] = ordH1x[i]
			Ay[i] = ordH1y[i]
			Bx[i] = ordH2x[i]
			By[i] = ordH2y[i]
		else
			costSame = (ordH1x[i] - prevAx)^2 + (ordH1y[i] - prevAy)^2 + (ordH2x[i] - prevBx)^2 + (ordH2y[i] - prevBy)^2
			costSwap = (ordH2x[i] - prevAx)^2 + (ordH2y[i] - prevAy)^2 + (ordH1x[i] - prevBx)^2 + (ordH1y[i] - prevBy)^2

			if (costSame <= costSwap)
				Ax[i] = ordH1x[i]
				Ay[i] = ordH1y[i]
				Bx[i] = ordH2x[i]
				By[i] = ordH2y[i]
			else
				Ax[i] = ordH2x[i]
				Ay[i] = ordH2y[i]
				Bx[i] = ordH1x[i]
				By[i] = ordH1y[i]
			endif
		endif

		prevAx = Ax[i]
		prevAy = Ay[i]
		prevBx = Bx[i]
		prevBy = By[i]
	endfor

	// -------------------------------------------------------------------------
	// 4. Build STS-relative polygon waves.
	// -------------------------------------------------------------------------
	Make/O/D/N=(nSel+2) tmp_HistFan_polyAx
	Make/O/D/N=(nSel+2) tmp_HistFan_polyAy
	Make/O/D/N=(nSel+2) tmp_HistFan_polyBx
	Make/O/D/N=(nSel+2) tmp_HistFan_polyBy

	WAVE polyAx = tmp_HistFan_polyAx
	WAVE polyAy = tmp_HistFan_polyAy
	WAVE polyBx = tmp_HistFan_polyBx
	WAVE polyBy = tmp_HistFan_polyBy

	polyAx[0] = 0
	polyAy[0] = 0
	polyBx[0] = 0
	polyBy[0] = 0

	for (i = 0; i < nSel; i += 1)
		polyAx[i+1] = Ax[i] - stsX_disp
		polyAy[i+1] = Ay[i] - stsY_disp
		polyBx[i+1] = Bx[i] - stsX_disp
		polyBy[i+1] = By[i] - stsY_disp
	endfor

	polyAx[nSel+1] = 0
	polyAy[nSel+1] = 0
	polyBx[nSel+1] = 0
	polyBy[nSel+1] = 0

	alphaFill = max(0, min(1, alphaFill))

	Variable aRGBA
	if (numtype(fanAlphaRGBA) == 0)
		aRGBA = round(fanAlphaRGBA)
	else
		aRGBA = round(65535 * alphaFill)
	endif
	aRGBA = max(0, min(65535, aRGBA))

	// -------------------------------------------------------------------------
	// 5. Draw fan polygons, then black max ray, then STS marker.
	// -------------------------------------------------------------------------
	SetDrawLayer/W=$winName UserFront

	SetDrawEnv/W=$winName xcoord=bottom, ycoord=left, fillpat=3, fillbgc=(rCol,gCol,bCol,aRGBA), fillfgc=(rCol,gCol,bCol,aRGBA), linefgc=(rCol,gCol,bCol,0), linethick=0
	DrawPoly/W=$winName stsX_disp, stsY_disp, 1, 1, polyAx, polyAy

	SetDrawEnv/W=$winName xcoord=bottom, ycoord=left, fillpat=3, fillbgc=(rCol,gCol,bCol,aRGBA), fillfgc=(rCol,gCol,bCol,aRGBA), linefgc=(rCol,gCol,bCol,0), linethick=0
	DrawPoly/W=$winName stsX_disp, stsY_disp, 1, 1, polyBx, polyBy

	// Draw black max ray last relative to fan polygons.
	if (numtype(maxTrace) == 0 && maxTrace >= 0 && maxTrace < nRay)
		SNS__GetRotatedRayEndpoints(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, maxTrace, displayXAxisAngle_deg, rotateDisplay, xCenterRot, yCenterRot, h1xD, h1yD, h2xD, h2yD)

		SetDrawEnv/W=$winName xcoord=bottom, ycoord=left, linefgc=(0,0,0,65535), linethick=maxLineSize
		DrawLine/W=$winName h1xD, h1yD, h2xD, h2yD
	endif

	// Draw STS marker as final draw object in this layer.
	if (numtype(stsMarkerRadius_nm) != 0 || stsMarkerRadius_nm <= 0)
		stsMarkerRadius_nm = 2.5
	endif

	SetDrawEnv/W=$winName xcoord=bottom, ycoord=left, fillfgc=(65535,65535,65535,65535), linefgc=(0,0,0,65535), linethick=max(1, maxLineSize*0.5)
	DrawOval/W=$winName stsX_disp-stsMarkerRadius_nm, stsY_disp-stsMarkerRadius_nm, stsX_disp+stsMarkerRadius_nm, stsY_disp+stsMarkerRadius_nm

	Variable/G v_HistFan_firstTrace = firstTrace
	Variable/G v_HistFan_lastTrace = lastTrace
	Variable/G v_HistFan_nSel = nSel
	Variable/G v_HistFan_STSx_disp = stsX_disp
	Variable/G v_HistFan_STSy_disp = stsY_disp
	Variable/G v_HistFan_maxTrace = maxTrace
	Variable/G v_HistFan_alphaRGBA = aRGBA
	Variable/G v_HistFan_angleCutIndex = cutIndex
	Variable/G v_HistFan_angleStartIndex = startIndex
	Variable/G v_HistFan_angleMaxGap = maxGap
	Variable/G v_HistFan_anglePeriod = pi

	return 0
End
//==============================================================================
// SNS__GetRotatedRayEndpoints
//==============================================================================

Function SNS__GetRotatedRayEndpoints(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, traceIndex, displayXAxisAngle_deg, rotateDisplay, xCenterRot, yCenterRot, h1xD, h1yD, h2xD, h2yD)
	WAVE Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
	Variable traceIndex, displayXAxisAngle_deg, rotateDisplay, xCenterRot, yCenterRot
	Variable &h1xD, &h1yD, &h2xD, &h2yD

	Variable thRot = displayXAxisAngle_deg*pi/180
	Variable h1x = Hit1x_List[traceIndex]
	Variable h1y = Hit1y_List[traceIndex]
	Variable h2x = Hit2x_List[traceIndex]
	Variable h2y = Hit2y_List[traceIndex]

	if (rotateDisplay)
		h1xD = (h1x - xCenterRot)*cos(thRot) + (h1y - yCenterRot)*sin(thRot)
		h1yD = -(h1x - xCenterRot)*sin(thRot) + (h1y - yCenterRot)*cos(thRot)

		h2xD = (h2x - xCenterRot)*cos(thRot) + (h2y - yCenterRot)*sin(thRot)
		h2yD = -(h2x - xCenterRot)*sin(thRot) + (h2y - yCenterRot)*cos(thRot)
	else
		h1xD = h1x
		h1yD = h1y
		h2xD = h2x
		h2yD = h2y
	endif

	return 0
End
//==============================================================================
// SNS_ExtractModesForFolder_LGapFamilyRays
//
// Longest-trajectory analogue of SNS_ExtractModesForFolder_WGapFamilyRays.
// Defines the family from the experimental-resolution-broadened zero-field
// ABS energy density rho_exp(ell_c).
//
// Ray alpha is mapped from contribution deficit:
//
//   d_i = (Amax - A_i)/(Amax - Ahwhm)
//
//   alpha_i = alphaMin + (alphaMax-alphaMin)
//             * [exp(-alphaBeta*d_i^alphaGamma)-exp(-alphaBeta)]
//               / [1-exp(-alphaBeta)]
//
// This gives alphaMax at Amax and alphaMin at the HWHM cutoff,
// with rapid drop-off for small deviations when alphaBeta is large.
//
// Optional display rotation:
//   displayXAxisAngle_deg sets the displayed x-axis direction in the original
//   image coordinate system.
//      displayXAxisAngle_deg = Bangle_deg      -> x-axis parallel to B
//      displayXAxisAngle_deg = Bangle_deg + 90 -> x-axis perpendicular to B
//
// Requires:
//   SNS_MakeRotatedImageForDisplay(...)
//   SNS_PlotChannelRay_Rotated(...)
//==============================================================================
Function SNS_ExtractModesForFolder_LGapFamilyRays(dfPath, [Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy, densityDL_nm, peakMinRelHeight, peakBox, alphaMin, alphaMax, alphaGamma, alphaBeta, familyColor, lineSize, useWeights, displayXAxisAngle_deg])
	String dfPath
	Variable Bangle_deg, BTK_barrier, STSx, STSy, Vortexx, Vortexy
	Variable densityDL_nm, peakMinRelHeight, peakBox
	Variable alphaMin, alphaMax, alphaGamma, alphaBeta
	Variable lineSize, useWeights
	Variable displayXAxisAngle_deg
	String familyColor

	String oldDF = GetDataFolder(1)

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	if (ParamIsDefault(Bangle_deg))
		Bangle_deg = 225
	endif
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
	if (ParamIsDefault(densityDL_nm))
		densityDL_nm = 1
	endif
	if (ParamIsDefault(peakMinRelHeight))
		peakMinRelHeight = 0.05
	endif
	if (ParamIsDefault(peakBox))
		peakBox = 3
	endif
	if (ParamIsDefault(alphaMin))
		alphaMin = 0.05
	endif
	if (ParamIsDefault(alphaMax))
		alphaMax = 1.00
	endif
	if (ParamIsDefault(alphaGamma))
		alphaGamma = 1.0
	endif
	if (ParamIsDefault(alphaBeta))
		alphaBeta = 8.0
	endif
	if (ParamIsDefault(familyColor))
		familyColor = "red"
	endif
	if (ParamIsDefault(lineSize))
		lineSize = 1.0
	endif
	if (ParamIsDefault(useWeights))
		useWeights = 0
	endif

	Variable rotateDisplay = 0
	if (!ParamIsDefault(displayXAxisAngle_deg) && numtype(displayXAxisAngle_deg) == 0)
		rotateDisplay = 1
	endif

	densityDL_nm = max(0.05, densityDL_nm)
	peakBox = max(1, round(peakBox))
	alphaMin = max(0, min(1, alphaMin))
	alphaMax = max(0, min(1, alphaMax))
	alphaGamma = max(0.1, alphaGamma)
	alphaBeta = max(0.1, alphaBeta)

	SNS_ExtractModesForFolder(dfPath, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, Vortexx=Vortexx, Vortexy=Vortexy, doDisplay=0)

	SetDataFolder $dfPath

	Wave/Z w_mask
	if (!WaveExists(w_mask))
		SetDataFolder $oldDF
		Abort "SNS_ExtractModesForFolder_LGapFamilyRays: missing w_mask."
	endif

	String imgList = WaveList("*_Z_mbgnd_xy", ";", "DIMS:2")
	if (strlen(imgList) == 0)
		SetDataFolder $oldDF
		Abort "SNS_ExtractModesForFolder_LGapFamilyRays: no *_Z_mbgnd_xy image found."
	endif

	String imgName = StringFromList(0, imgList, ";")
	Wave image = $imgName

	Wave L_N_List
	Wave W_eff_List
	Wave Hit1x_List
	Wave Hit1y_List
	Wave Hit2x_List
	Wave Hit2y_List
	Wave T_eff_List

	Variable nRay = numpnts(L_N_List)

	Make/O/D/N=(nRay) tmp_L_N_List_nm_forLGap
	tmp_L_N_List_nm_forLGap = L_N_List * 1e9
	Wave L_nm = tmp_L_N_List_nm_forLGap

	if (useWeights)
		SNS_MakeLResolutionDensity(L_nm, "L_N_density_exp_nm", weightWave=T_eff_List, dl_nm=densityDL_nm, doDisplay=0)
	else
		SNS_MakeLResolutionDensity(L_nm, "L_N_density_exp_nm", dl_nm=densityDL_nm, doDisplay=0)
	endif

	Wave rho = L_N_density_exp_nm

	Variable ellGap_nm
	ellGap_nm = SNS_FindHighLPeakInResolutionDensity(rho, minRelHeight=peakMinRelHeight, box=peakBox, doDisplay=0)

	if (numtype(ellGap_nm) != 0)
		SetDataFolder $oldDF
		Abort "SNS_ExtractModesForFolder_LGapFamilyRays: failed to find ell_gap peak."
	endif

	if (useWeights)
		SNS_MakeLGapRayContributions(L_nm, ellGap_nm, "L_gap_family", weightWave=T_eff_List, alphaMin=alphaMin, alphaMax=alphaMax, alphaGamma=alphaGamma)
	else
		SNS_MakeLGapRayContributions(L_nm, ellGap_nm, "L_gap_family", alphaMin=alphaMin, alphaMax=alphaMax, alphaGamma=alphaGamma)
	endif

	Wave contrib = L_gap_family_contribution
	Wave mask = L_gap_family_family_mask
	Wave alpha = L_gap_family_alpha

	WaveStats/Q contrib
	Variable Amax = V_max
	Variable Ahwhm = 0.5 * Amax

	// -------------------------------------------------------------------------
	// Recompute alpha with deficit-exponential mapping.
	// -------------------------------------------------------------------------
	Variable i
	Variable alphaDeficit, expMin

	alpha = 0
	mask = 0
	expMin = exp(-alphaBeta)

	if (Amax > 0 && numtype(Amax) == 0)
		for (i = 0; i < nRay; i += 1)
			if (contrib[i] >= Ahwhm)
				mask[i] = 1

				if (Amax > Ahwhm)
					alphaDeficit = (Amax - contrib[i]) / (Amax - Ahwhm)
					alphaDeficit = max(0, min(1, alphaDeficit))

					alpha[i] = alphaMin + (alphaMax - alphaMin) * (exp(-alphaBeta * alphaDeficit^alphaGamma) - expMin) / (1 - expMin)
				else
					alpha[i] = alphaMax
				endif

				alpha[i] = max(alphaMin, min(alphaMax, alpha[i]))
			endif
		endfor
	endif

	Variable/G L_gap_family_alphaGamma = alphaGamma
	Variable/G L_gap_family_alphaBeta = alphaBeta

	Variable bestTrace = NaN
	Variable bestA = -Inf

	for (i = 0; i < nRay; i += 1)
		if (contrib[i] > bestA)
			bestA = contrib[i]
			bestTrace = i
		endif
	endfor

	// -------------------------------------------------------------------------
	// Display topography, optionally rotated.
	// -------------------------------------------------------------------------
	Variable xOrig0 = DimOffset(image, 0)
	Variable xOrig1 = xOrig0 + DimDelta(image, 0) * (DimSize(image, 0) - 1)
	Variable yOrig0 = DimOffset(image, 1)
	Variable yOrig1 = yOrig0 + DimDelta(image, 1) * (DimSize(image, 1) - 1)

	Variable xCenterRot = 0.5 * (min(xOrig0, xOrig1) + max(xOrig0, xOrig1))
	Variable yCenterRot = 0.5 * (min(yOrig0, yOrig1) + max(yOrig0, yOrig1))

	String displayImgName = NameOfWave(image)

	if (rotateDisplay)
		SNS_MakeRotatedImageForDisplay(image, "tmp_LGapFamily_rot_image", displayXAxisAngle_deg)
		displayImgName = "tmp_LGapFamily_rot_image"
	endif

	Wave displayImage = $displayImgName

	SNS_DisplayWithScales(displayImage, cmap="grayC")
	String winImage = WinName(0, 1, 1)

	SNS_StyleRayDisplayImage(displayImage, winImage)

	DoWindow/F $winImage

	Variable rCol, gCol, bCol

	if (StringMatch(familyColor, "black"))
		rCol = 0; gCol = 0; bCol = 0
	elseif (StringMatch(familyColor, "white"))
		rCol = 65535; gCol = 65535; bCol = 65535
	elseif (StringMatch(familyColor, "cyan"))
		rCol = 1; gCol = 16019; bCol = 65535
	elseif (StringMatch(familyColor, "blue"))
		rCol = 0; gCol = 0; bCol = 65535
	elseif (StringMatch(familyColor, "green"))
		rCol = 0; gCol = 45000; bCol = 0
	elseif (StringMatch(familyColor, "magenta"))
		rCol = 65535; gCol = 0; bCol = 65535
	elseif (StringMatch(familyColor, "yellow"))
		rCol = 65535; gCol = 50000; bCol = 0
	else
		rCol = 65535; gCol = 0; bCol = 0
	endif

	String traceList, lastTrace
	Variable nFamily = 0
	Variable aRGBA

	for (i = 0; i < nRay; i += 1)
		if (mask[i] < 0.5)
			continue
		endif

		DoWindow/F $winImage

		if (rotateDisplay)
			SNS_PlotChannelRay_Rotated(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, i, displayXAxisAngle_deg, xCenterRot, yCenterRot)
		else
			SNS_PlotChannelRay(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, i)
		endif

		traceList = TraceNameList(winImage, ";", 1)
		lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

		aRGBA = round(65535 * alpha[i])
		aRGBA = max(0, min(65535, aRGBA))

		ModifyGraph/W=$winImage rgb($lastTrace)=(rCol,gCol,bCol,aRGBA)
		ModifyGraph/W=$winImage lsize($lastTrace)=lineSize

		nFamily += 1
	endfor

	if (numtype(bestTrace) == 0)
		DoWindow/F $winImage

		if (rotateDisplay)
			SNS_PlotChannelRay_Rotated(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, bestTrace, displayXAxisAngle_deg, xCenterRot, yCenterRot)
		else
			SNS_PlotChannelRay(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, bestTrace)
		endif

		traceList = TraceNameList(winImage, ";", 1)
		lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

		ModifyGraph/W=$winImage rgb($lastTrace)=(0,0,0,65535)
		ModifyGraph/W=$winImage lsize($lastTrace)=2.5
	endif

	// -------------------------------------------------------------------------
	// STS marker and B arrow.
	// -------------------------------------------------------------------------
	Make/O/N=1 tmpSTSx, tmpSTSy

	if (rotateDisplay)
		Variable thSTS = displayXAxisAngle_deg * pi / 180
		tmpSTSx[0] = (STSx - xCenterRot) * cos(thSTS) + (STSy - yCenterRot) * sin(thSTS)
		tmpSTSy[0] = -(STSx - xCenterRot) * sin(thSTS) + (STSy - yCenterRot) * cos(thSTS)
	else
		tmpSTSx[0] = STSx
		tmpSTSy[0] = STSy
	endif

	AppendToGraph/W=$winImage tmpSTSy vs tmpSTSx

	traceList = TraceNameList(winImage, ";", 1)
	lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

	ModifyGraph/W=$winImage mode($lastTrace)=3
	ModifyGraph/W=$winImage marker($lastTrace)=19
	ModifyGraph/W=$winImage msize($lastTrace)=5
	ModifyGraph/W=$winImage rgb($lastTrace)=(65535,65535,65535,65535)

	DoWindow/F $winImage
	SNS_TAGSetpoint()

	Variable xMin = DimOffset(displayImage, 0)
	Variable xMax = xMin + DimDelta(displayImage, 0) * (DimSize(displayImage, 0) - 1)
	Variable yMin = DimOffset(displayImage, 1)
	Variable yMax = yMin + DimDelta(displayImage, 1) * (DimSize(displayImage, 1) - 1)

	Variable pos_x = xMax - 0.15 * (xMax - xMin)
	Variable pos_y = yMin + 0.15 * (yMax - yMin)

	Variable displayBangle_deg = Bangle_deg
	if (rotateDisplay)
		displayBangle_deg = Bangle_deg - displayXAxisAngle_deg
	endif

	SNS_DrawImageArrow(displayImage, displayBangle_deg, "B", pos_x=pos_x, pos_y=pos_y, scale=1.5, color="black")

	String txt
	txt = "Gap-defining ℓ family\r"
	txt += "ℓ\\Bgap\\M = " + num2str(round(10 * ellGap_nm) / 10) + " nm\r"
	txt += "N\\Bfamily\\M = " + num2str(nFamily) + " / " + num2str(nRay) + "\r"
	txt += "criterion: A\\Bi\\M ≥ 0.5 A\\Bmax\\M\r"
	txt += "alpha β = " + num2str(round(10 * alphaBeta) / 10) + ", γ = " + num2str(round(10 * alphaGamma) / 10) + "\r"
	txt += "A\\Bmax\\M = " + num2str(round(1000 * Amax) / 1000)

	if (rotateDisplay)
		txt += "\rdisplay x-axis = " + num2str(round(10 * displayXAxisAngle_deg) / 10) + "°"
	endif

	TextBox/C/N=LGapFamilyInfo/F=0/A=RT/X=0/Y=0/B=(65535,65535,65535,39321) txt

	GetWindow $winImage, wsize
	Variable wPx = abs(V_right - V_left)
	Variable hPx = abs(V_bottom - V_top)
	Variable maxDim = max(wPx, hPx)

	if (maxDim > 512)
		Variable ex = 512 / maxDim
		ModifyGraph/W=$winImage expand=ex
	else
		ModifyGraph/W=$winImage expand=1
	endif

	Display/K=1 rho
	String winRho = WinName(0, 1, 1)

	ModifyGraph/W=$winRho mode=0
	ModifyGraph/W=$winRho tick=2, mirror=2, standoff=0
	Label/W=$winRho bottom "candidate ℓ\\Bc\\M (nm)"
	Label/W=$winRho left "resolution-broadened zero-field weight"

	SetDrawLayer UserFront
	SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(65535,0,0), linethick=2, dash=3
	DrawLine ellGap_nm, 0, ellGap_nm, 1.05

	ModifyGraph/W=$winRho width=260
	ModifyGraph/W=$winRho height=220

	KillWaves/Z tmp_L_N_List_nm_forLGap

	SetDataFolder $oldDF

	return ellGap_nm
End
//==============================================================================
// SNS_DrawHitRaysOnTopWindow
//
// Purpose:
//   Draw ray-traced SNS channels on the current/topmost graph window using
//   supplied hit-coordinate waves.
//
//   This is a lightweight display helper. It does not build channels, does not
//   identify special modes, and does not create a new graph. It only overlays
//   two-point XY traces on the active graph.
//
// Inputs:
//   Hit1x_List, Hit1y_List : first S-hit coordinates for each channel
//   Hit2x_List, Hit2y_List : second S-hit coordinates for each channel
//
// Optional:
//   maxRays : maximum number of rays to draw.
//             -1 or omitted means draw all rays.
//             If maxRays > 0 and fewer rays are desired, the function
//             subsamples the channel list uniformly.
//
//   rayColor : base color string for rays.
//              Supported values:
//                 "gray"
//                 "black"
//                 "white"
//                 "red"
//                 "blue"
//                 "green"
//                 "cyan"
//                 "magenta"
//                 "yellow"
//              Default: "white".
//
//   colorTransparency : alpha/transparency control.
//                       Range is clipped to [0,1].
//                       0 = fully transparent.
//                       1 = fully opaque.
//                       Default: 0.20.
//
//   lineSize : line thickness.
//              Default: 0.25.
//
// Outputs:
//   return : numeric
//            0 on successful completion.
//
// Notes:
//   Coordinates must already match the coordinate scaling of the topmost graph.
//   The function appends each ray as a two-point XY trace and styles the last
//   appended trace using TraceNameList(...).
//==============================================================================
Function SNS_DrawHitRaysOnTopWindow(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, [maxRays, rayColor, colorTransparency, lineSize])
    Wave Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
    Variable maxRays, colorTransparency, lineSize
    String rayColor

    // ---------------- Defaults ----------------
    if (ParamIsDefault(maxRays))
        maxRays = -1
    endif
    if (ParamIsDefault(rayColor))
        rayColor = "white"
    endif
    if (ParamIsDefault(colorTransparency))
        colorTransparency = 0.20
    endif
    if (ParamIsDefault(lineSize))
        lineSize = 0.25
    endif

    // ---------------- Require active graph ----------------
    String gName = WinName(0, 1, 1)
    if (strlen(gName) == 0)
        Abort "SNS_DrawHitRaysOnTopWindow: no graph window is active."
    endif

    // ---------------- Clamp alpha ----------------
    colorTransparency = max(0, min(1, colorTransparency))
    Variable aa = round(65535 * colorTransparency)

    // ---------------- Choose base color ----------------
    Variable rr, gg, bb

    strswitch(LowerStr(rayColor))
        case "black":
            rr = 0
            gg = 0
            bb = 0
            break

        case "white":
            rr = 65535
            gg = 65535
            bb = 65535
            break

        case "red":
            rr = 65535
            gg = 0
            bb = 0
            break

        case "blue":
            rr = 0
            gg = 0
            bb = 65535
            break

        case "green":
            rr = 0
            gg = 45000
            bb = 0
            break

        case "cyan":
            rr = 0
            gg = 52000
            bb = 65535
            break

        case "magenta":
            rr = 65535
            gg = 0
            bb = 65535
            break

        case "yellow":
            rr = 65535
            gg = 52000
            bb = 0
            break

        case "gray":
        default:
            rr = 54000
            gg = 54000
            bb = 54000
            break
    endswitch

    // ---------------- Validate ray list ----------------
    Variable nRay = numpnts(Hit1x_List)

    if (nRay <= 0)
        return 0
    endif

    if (numpnts(Hit1y_List) != nRay || numpnts(Hit2x_List) != nRay || numpnts(Hit2y_List) != nRay)
        Abort "SNS_DrawHitRaysOnTopWindow: endpoint lists have inconsistent lengths."
    endif

    // ---------------- Optional uniform subsampling ----------------
    Variable step = 1

    if (maxRays > 0 && nRay > maxRays)
        step = ceil(nRay / maxRays)
    endif

    // ---------------- Append ray traces ----------------
    Variable j
    String idxStr, xName, yName
    String traceList, lastTrace

    for (j = 0; j < nRay; j += step)

        sprintf idxStr, "%04d", j

        xName = "tmpTopRayX_" + idxStr
        yName = "tmpTopRayY_" + idxStr

        Make/O/D/N=2 $xName, $yName
        Wave rayX = $xName
        Wave rayY = $yName

        rayX[0] = Hit1x_List[j]
        rayY[0] = Hit1y_List[j]
        rayX[1] = Hit2x_List[j]
        rayY[1] = Hit2y_List[j]

        AppendToGraph/W=$gName rayY vs rayX

        traceList = TraceNameList(gName, ";", 1)
        lastTrace = StringFromList(ItemsInList(traceList)-1, traceList, ";")

        ModifyGraph/W=$gName rgb($lastTrace)=(rr, gg, bb, aa)
        ModifyGraph/W=$gName lsize($lastTrace)=lineSize
    endfor

    return 0
End

//==============================================================================
// SNS_RotateHitCoordinatesAroundSTS
//
// Purpose:
//   Rotate ray-hit coordinates in an absolute coordinate system around the
//   selected STS position.
//
//   This is intended for the situation where:
//     - Hit1x/y and Hit2x/y are in the same absolute coordinates as STSx, STSy,
//     - an image/island is viewed at a different rotation angle,
//     - the STS position should remain fixed,
//     - the hit positions should move with the rotated island/image.
//
//   Mathematically, each point r is transformed as:
//
//       r_rot = r_STS + R(angleDeg) * (r - r_STS)
//
//   Therefore:
//     - STSx, STSy are unchanged,
//     - all hit coordinates rotate around STSx, STSy,
//     - distances from the STS point are preserved,
//     - relative ray geometry is preserved.
//
// Inputs:
//   Hit1x_List, Hit1y_List : first S-hit coordinates
//   Hit2x_List, Hit2y_List : second S-hit coordinates
//   angleDeg               : rotation angle in degrees
//   STSx, STSy             : rotation center in same coordinates as hits
//
// Optional:
//   outPrefix : prefix for output waves.
//               Default: "rotSTS"
//
// Outputs:
//   Creates/overwrites:
//      outPrefix + "_Hit1x"
//      outPrefix + "_Hit1y"
//      outPrefix + "_Hit2x"
//      outPrefix + "_Hit2y"
//
//   return : numeric
//            0 on successful completion.
//
// Sign convention:
//   Positive angle follows the standard mathematical convention:
//
//      x' = x0 + cos(a)*(x-x0) - sin(a)*(y-y0)
//      y' = y0 + sin(a)*(x-x0) + cos(a)*(y-y0)
//
//   If this appears opposite to Igor's image rotation convention in your
//   display, use -angleDeg.
//
//==============================================================================
Function SNS_RotateHitCoordinatesAroundSTS(Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List, angleDeg, STSx, STSy, [outPrefix])
    Wave Hit1x_List, Hit1y_List, Hit2x_List, Hit2y_List
    Variable angleDeg, STSx, STSy
    String outPrefix

    // ---------------- Defaults ----------------
    if (ParamIsDefault(outPrefix))
        outPrefix = "rotSTS"
    endif

    // ---------------- Validate input ----------------
    Variable nRay = numpnts(Hit1x_List)

    if (nRay <= 0)
        return 0
    endif

    if (numpnts(Hit1y_List) != nRay || numpnts(Hit2x_List) != nRay || numpnts(Hit2y_List) != nRay)
        Abort "SNS_RotateHitCoordinatesAroundSTS: endpoint lists have inconsistent lengths."
    endif

    // ---------------- Create output waves ----------------
    String h1xName = outPrefix + "_Hit1x"
    String h1yName = outPrefix + "_Hit1y"
    String h2xName = outPrefix + "_Hit2x"
    String h2yName = outPrefix + "_Hit2y"

    Make/O/D/N=(nRay) $h1xName, $h1yName, $h2xName, $h2yName

    Wave outH1x = $h1xName
    Wave outH1y = $h1yName
    Wave outH2x = $h2xName
    Wave outH2y = $h2yName

    // ---------------- Rotation ----------------
    Variable ang = angleDeg * pi / 180
    Variable ca = cos(ang)
    Variable sa = sin(ang)

    outH1x = STSx + ca * (Hit1x_List[p] - STSx) - sa * (Hit1y_List[p] - STSy)
    outH1y = STSy + sa * (Hit1x_List[p] - STSx) + ca * (Hit1y_List[p] - STSy)

    outH2x = STSx + ca * (Hit2x_List[p] - STSx) - sa * (Hit2y_List[p] - STSy)
    outH2y = STSy + sa * (Hit2x_List[p] - STSx) + ca * (Hit2y_List[p] - STSy)

    return 0
End

//==============================================================================
// SNS_RayTrace_Hist_LN_Weff
//
// Purpose:
//   Run 2D/3D ray tracing for a single STM position and create ray-trace
//   histograms of L_N, W_eff, and, for 3D, W_geom.
//
//   The extractor calls are intentionally kept as the geometry source:
//
//      2D: SNS_ExtractModesForFolder(...)
//      3D: SNS_ExtractModes3DForFolder(...)
//
// Interface transparency:
//   This function no longer accepts BTK_barrier / Zbarrier.
//   Interface transparency is controlled globally through SNS_Settings and is
//   applied inside the 2D/3D channel construction used by the extractors.
//
// Output folders:
//   H_nm <= 0:
//      <dfPath>:RayTraceHist
//
//   H_nm > 0:
//      <dfPath>:RayTraceHist_3D
//
// Histogram binning:
//   Weighted and unweighted histograms use the same manual binning path.
//   The only difference is the per-ray weight:
//
//      wt = wChan[i]    if valid weights are used
//      wt = 1           otherwise
//
// Extreme-tail estimator:
//   The upper-tail range is [0.9*global_max, global_max].
//   The function searches for a supported contiguous cluster in angular
//   ray-family order. If phiList_rad exists, diagnostics use phi modulo pi,
//   sorted and rotated so the branch cut lies in the largest angular gap.
//   This preserves continuum ray-tracing behavior for imported weighted
//   folders whose list order is not angular.
//
//   For the selected cluster it stores:
//      *_extreme_max_nm
//      *_extreme_median_nm
//      *_extreme_sigma_nm
//      *_extreme_max_trace
//      *_extreme_median_trace
//      *_extreme_cluster_start
//      *_extreme_cluster_size
//      *_extreme_minClusterSize
//
// Plotting:
//   Histogram plots use only waves in RayTraceHist / RayTraceHist_3D.
//   useExistingChannels=1 skips ray tracing and histograms the existing
//   canonical channel waves already present in dfPath. This is intended for
//   rectmode / mode-weighted sidecar folders.
//   The max marker is drawn by default.
//   The median marker is drawn only if plotMedian=1.
//
// Returns:
//   v_LN_extreme_max_trace for the selected output mode.
//
// Notes:
//   This is a ray-diagnostic / histogram helper. Long term it belongs in a
//   dedicated SNS_RayDiagnostics.ipf rather than SNS_SpatialMaps.ipf.
//==============================================================================
Function SNS_RayTrace_Hist_LN_Weff(dfPath, [Bangle_deg, STSx, STSy, Vortexx, Vortexy, H_nm, maxPath_nm, nBins, binL_nm, binW_nm, binT, histLMin_nm, histLMax_nm, histWMin_nm, histWMax_nm, histTMin, histTMax, absW, normalize, useWeights, padCityscape, showFitInfo, plotMedian, doDisplay, useExistingChannels])
    String dfPath
    Variable Bangle_deg, STSx, STSy, Vortexx, Vortexy
    Variable H_nm, maxPath_nm
    Variable nBins, binL_nm, binW_nm, binT, histLMin_nm, histLMax_nm, histWMin_nm, histWMax_nm, histTMin, histTMax
    Variable absW, normalize, useWeights, padCityscape, showFitInfo, plotMedian, doDisplay
    Variable useExistingChannels

    String oldDF = GetDataFolder(1)

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    // ---------------- Defaults ----------------
    if (ParamIsDefault(Bangle_deg))
        Bangle_deg = 225
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
    if (ParamIsDefault(H_nm))
        H_nm = 0
    endif
    if (ParamIsDefault(maxPath_nm))
        maxPath_nm = 1000
    endif

    Variable use3D = (H_nm > 0)

    if (ParamIsDefault(nBins))
        nBins = 50
    endif

    // Fixed-bin controls for consistent comparison histograms.
    // If binL_nm/binW_nm are omitted, the legacy nBins behavior is used.
    // T_eff is dimensionless, so it defaults to common [0, 1] bins.
    Variable useFixedBinL = !ParamIsDefault(binL_nm)
    Variable useFixedBinW = !ParamIsDefault(binW_nm)

    if (ParamIsDefault(binL_nm))
        binL_nm = NaN
    else
        binL_nm = abs(binL_nm)
        if (binL_nm <= 0 || numtype(binL_nm) != 0)
            SetDataFolder $oldDF
            Abort "SNS_RayTrace_Hist_LN_Weff: binL_nm must be positive."
        endif
    endif

    if (ParamIsDefault(binW_nm))
        binW_nm = NaN
    else
        binW_nm = abs(binW_nm)
        if (binW_nm <= 0 || numtype(binW_nm) != 0)
            SetDataFolder $oldDF
            Abort "SNS_RayTrace_Hist_LN_Weff: binW_nm must be positive."
        endif
    endif

    if (ParamIsDefault(binT))
        binT = 0.05
    else
        binT = abs(binT)
        if (binT <= 0 || numtype(binT) != 0)
            SetDataFolder $oldDF
            Abort "SNS_RayTrace_Hist_LN_Weff: binT must be positive."
        endif
    endif

    if (ParamIsDefault(histLMin_nm))
        histLMin_nm = NaN
    endif
    if (ParamIsDefault(histLMax_nm))
        histLMax_nm = NaN
    endif
    if (ParamIsDefault(histWMin_nm))
        histWMin_nm = NaN
    endif
    if (ParamIsDefault(histWMax_nm))
        histWMax_nm = NaN
    endif
    if (ParamIsDefault(histTMin))
        histTMin = 0
    endif
    if (ParamIsDefault(histTMax))
        histTMax = 1
    endif

    if (ParamIsDefault(absW))
        absW = 1
    endif
    if (ParamIsDefault(normalize))
        normalize = 1
    endif

    // Default weighting policy:
    //   2D : unweighted legacy behavior.
    //   3D : weighted by wChan if available.
    if (ParamIsDefault(useWeights))
        useWeights = use3D
    endif

    if (ParamIsDefault(padCityscape))
        padCityscape = 0
    endif
    if (ParamIsDefault(showFitInfo))
        showFitInfo = 1
    endif
    if (ParamIsDefault(plotMedian))
        plotMedian = 0
    endif
    if (ParamIsDefault(doDisplay))
        doDisplay = 1
    endif
    if (ParamIsDefault(useExistingChannels))
        useExistingChannels = 0
    endif

    nBins = max(2, round(nBins))
    useExistingChannels = round(useExistingChannels)
    useExistingChannels = max(0, min(1, useExistingChannels))

    // -------------------------------------------------------------------------
    // 1. Run ray tracing / mode extraction without extractor displays.
    // -------------------------------------------------------------------------
    if (!useExistingChannels)
        if (use3D)
            SNS_ExtractModes3DForFolder(dfPath, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, Vortexx=Vortexx, Vortexy=Vortexy, H_nm=H_nm, maxPath_nm=maxPath_nm, doDisplay=0)
        else
            SNS_ExtractModesForFolder(dfPath, Bangle_deg=Bangle_deg, STSx=STSx, STSy=STSy, Vortexx=Vortexx, Vortexy=Vortexy, doDisplay=0)
        endif
    endif

    SetDataFolder $dfPath

    Wave/Z src_L_N_List = L_N_List_nm
    Wave/Z src_W_eff_List = W_eff_List_nm
    Wave/Z src_W_geom_List = W_geom_List_nm
    Wave/Z src_wChan = wChan
    Wave/Z src_T_eff_List = T_eff_List
    Wave/Z src_phiList = phiList_rad
    Wave/Z src_Hit1x_List = Hit1x_List_nm
    Wave/Z src_Hit1y_List = Hit1y_List_nm
    Wave/Z src_Hit2x_List = Hit2x_List_nm
    Wave/Z src_Hit2y_List = Hit2y_List_nm

    if (!WaveExists(src_L_N_List))
        SetDataFolder $oldDF
        Abort "SNS_RayTrace_Hist_LN_Weff: missing L_N_List_nm after ray tracing."
    endif

    if (!WaveExists(src_W_eff_List))
        SetDataFolder $oldDF
        Abort "SNS_RayTrace_Hist_LN_Weff: missing W_eff_List_nm after ray tracing."
    endif

    Variable nRay = numpnts(src_L_N_List)

    if (nRay <= 0)
        SetDataFolder $oldDF
        Abort "SNS_RayTrace_Hist_LN_Weff: no valid ray-traced channels."
    endif

    if (numpnts(src_W_eff_List) != nRay)
        SetDataFolder $oldDF
        Abort "SNS_RayTrace_Hist_LN_Weff: L_N_List_nm and W_eff_List_nm have inconsistent lengths."
    endif

    if (!WaveExists(src_Hit1x_List) || !WaveExists(src_Hit1y_List) || !WaveExists(src_Hit2x_List) || !WaveExists(src_Hit2y_List))
        SetDataFolder $oldDF
        Abort "SNS_RayTrace_Hist_LN_Weff: missing endpoint Hit*List_nm waves after ray tracing."
    endif

    if ((numpnts(src_Hit1x_List) != nRay) || (numpnts(src_Hit1y_List) != nRay) || (numpnts(src_Hit2x_List) != nRay) || (numpnts(src_Hit2y_List) != nRay))
        SetDataFolder $oldDF
        Abort "SNS_RayTrace_Hist_LN_Weff: endpoint Hit*List_nm waves have inconsistent lengths."
    endif

    Variable hasGeomW = use3D && WaveExists(src_W_geom_List) && (numpnts(src_W_geom_List) == nRay)

    if (use3D && !hasGeomW)
        SNS_Log("SNS_RayTrace_Hist_LN_Weff: 3D mode requested, but valid W_geom_List_nm not found. Only W_eff histogram will be produced.", level="WARN")
    endif

    Variable hasWeights = useWeights && WaveExists(src_wChan) && (numpnts(src_wChan) == nRay)
    Variable hasT = WaveExists(src_T_eff_List) && (numpnts(src_T_eff_List) == nRay)

    if (useWeights && !hasWeights)
        SNS_Log("SNS_RayTrace_Hist_LN_Weff: useWeights=1 requested, but valid wChan not found. Falling back to unweighted histogram.", level="WARN")
    endif

    // -------------------------------------------------------------------------
    // 2. Create clean output folder and copy extractor products into it.
    // -------------------------------------------------------------------------
    if (use3D)
        KillDataFolder/Z RayTraceHist_3D
        NewDataFolder/O RayTraceHist_3D

        Duplicate/O src_L_N_List, :RayTraceHist_3D:L_N_List_3D_nm
        Duplicate/O src_W_eff_List, :RayTraceHist_3D:W_eff_List_3D_nm
        Duplicate/O src_Hit1x_List, :RayTraceHist_3D:Hit1x_List_3D_nm
        Duplicate/O src_Hit1y_List, :RayTraceHist_3D:Hit1y_List_3D_nm
        Duplicate/O src_Hit2x_List, :RayTraceHist_3D:Hit2x_List_3D_nm
        Duplicate/O src_Hit2y_List, :RayTraceHist_3D:Hit2y_List_3D_nm

        Duplicate/O src_L_N_List, :RayTraceHist_3D:tmp_L_N_work
        Duplicate/O src_W_eff_List, :RayTraceHist_3D:tmp_W_eff_work

        if (hasGeomW)
            Duplicate/O src_W_geom_List, :RayTraceHist_3D:W_geom_List_3D_nm
            Duplicate/O src_W_geom_List, :RayTraceHist_3D:tmp_W_geom_work
        endif

        if (WaveExists(src_wChan) && (numpnts(src_wChan) == nRay))
            Duplicate/O src_wChan, :RayTraceHist_3D:wChan_3D
            Duplicate/O src_wChan, :RayTraceHist_3D:tmp_wChan_work
        endif

        if (hasT)
            Duplicate/O src_T_eff_List, :RayTraceHist_3D:T_eff_List_3D
            Duplicate/O src_T_eff_List, :RayTraceHist_3D:tmp_T_eff_work
        endif

        if (WaveExists(src_phiList) && (numpnts(src_phiList) == nRay))
            Duplicate/O src_phiList, :RayTraceHist_3D:phiList_3D_rad
        endif
    else
        KillDataFolder/Z RayTraceHist
        NewDataFolder/O RayTraceHist

        Duplicate/O src_L_N_List, :RayTraceHist:L_N_List_nm
        Duplicate/O src_W_eff_List, :RayTraceHist:W_eff_List_nm
        Duplicate/O src_Hit1x_List, :RayTraceHist:Hit1x_List_nm
        Duplicate/O src_Hit1y_List, :RayTraceHist:Hit1y_List_nm
        Duplicate/O src_Hit2x_List, :RayTraceHist:Hit2x_List_nm
        Duplicate/O src_Hit2y_List, :RayTraceHist:Hit2y_List_nm

        Duplicate/O src_L_N_List, :RayTraceHist:tmp_L_N_work
        Duplicate/O src_W_eff_List, :RayTraceHist:tmp_W_eff_work

        if (WaveExists(src_wChan) && (numpnts(src_wChan) == nRay))
            Duplicate/O src_wChan, :RayTraceHist:wChan
            Duplicate/O src_wChan, :RayTraceHist:tmp_wChan_work
        endif

        if (hasT)
            Duplicate/O src_T_eff_List, :RayTraceHist:T_eff_List
            Duplicate/O src_T_eff_List, :RayTraceHist:tmp_T_eff_work
        endif

        if (WaveExists(src_phiList) && (numpnts(src_phiList) == nRay))
            Duplicate/O src_phiList, :RayTraceHist:phiList_rad
        endif
    endif

    // -------------------------------------------------------------------------
    // 3. Clean generic extractor/function leftovers from dfPath.
    //    The image data, settings folder, w_mask, and output folder are preserved.
    // -------------------------------------------------------------------------
    SetDataFolder $dfPath

    if (!useExistingChannels)
        KillWaves/Z L_N_List, W_eff_List, W_geom_List, wChan, phiList, thetaList
        KillWaves/Z L_N_List_nm, W_eff_List_nm, W_geom_List_nm, phiList_rad, thetaList_rad
        KillWaves/Z Hist_L_N_nm, Hist_W_eff_nm, Hist_W_geom_nm
        KillWaves/Z Hist_T_eff, Hist_T_eff_weighted
        KillWaves/Z Hist_L_N_weighted_nm, Hist_W_eff_weighted_nm, Hist_W_geom_weighted_nm
        KillWaves/Z Hist_L_N_plot_nm, Hist_W_eff_plot_nm, Hist_W_geom_plot_nm, Hist_T_eff_plot
        KillWaves/Z Hist_L_N_weighted_plot_nm, Hist_W_eff_weighted_plot_nm, Hist_W_geom_weighted_plot_nm, Hist_T_eff_weighted_plot
        KillWaves/Z L_N_extreme_nm, W_eff_extreme_nm, W_geom_extreme_nm
        KillWaves/Z tmp_L_cluster_vals, tmp_WEff_cluster_vals, tmp_WGeom_cluster_vals
        KillWaves/Z LN_center_x_nm, LN_center_y, Weff_center_x_nm, Weff_center_y, Wgeom_center_x_nm, Wgeom_center_y

        KillVariables/Z v_RayTraceHist_is3D, v_RayTraceHist_H_nm, v_RayTraceHist_maxPath_nm, v_RayTraceHist_padCityscape
        KillVariables/Z v_RayTraceHist_hasGeomW, v_RayTraceHist_useWeights
        KillVariables/Z v_LN_max_nm, v_LN_extreme_low_nm, v_LN_extreme_high_nm, v_LN_extreme_n
        KillVariables/Z v_LN_extreme_center_nm, v_LN_extreme_sigma_nm, v_LN_extreme_trace
        KillVariables/Z v_LN_extreme_cluster_start, v_LN_extreme_cluster_size, v_LN_extreme_minClusterSize
        KillVariables/Z v_Weff_max_nm, v_Weff_extreme_low_nm, v_Weff_extreme_high_nm, v_Weff_extreme_n
        KillVariables/Z v_Weff_extreme_center_nm, v_Weff_extreme_sigma_nm, v_Weff_extreme_trace
        KillVariables/Z v_Weff_extreme_cluster_start, v_Weff_extreme_cluster_size, v_Weff_extreme_minClusterSize
        KillVariables/Z v_Wgeom_max_nm, v_Wgeom_extreme_low_nm, v_Wgeom_extreme_high_nm, v_Wgeom_extreme_n
        KillVariables/Z v_Wgeom_extreme_center_nm, v_Wgeom_extreme_sigma_nm, v_Wgeom_extreme_trace
        KillVariables/Z v_Wgeom_extreme_cluster_start, v_Wgeom_extreme_cluster_size, v_Wgeom_extreme_minClusterSize

        KillStrings/Z s_RayTraceHist_HistL, s_RayTraceHist_HistW, s_RayTraceHist_HistW_eff, s_RayTraceHist_HistW_geom
        KillStrings/Z s_RayTraceHist_HistT, s_RayTraceHist_HistT_eff
        KillStrings/Z s_RayTraceHist_PlotL, s_RayTraceHist_PlotW, s_RayTraceHist_PlotW_eff, s_RayTraceHist_PlotW_geom
        KillStrings/Z s_RayTraceHist_PlotT, s_RayTraceHist_PlotT_eff
        KillStrings/Z s_RayTraceHist_WSource_eff, s_RayTraceHist_WSource_geom
    endif

    if (use3D)
        SetDataFolder :RayTraceHist_3D
    else
        SetDataFolder :RayTraceHist
    endif

    Wave L_work = tmp_L_N_work
    Wave W_eff_work = tmp_W_eff_work
    Wave/Z W_geom_work = tmp_W_geom_work
    Wave/Z Wt_work = tmp_wChan_work
    Wave/Z T_work = tmp_T_eff_work

    // -------------------------------------------------------------------------
    // 4. Metadata in the output folder.
    // -------------------------------------------------------------------------
    if (use3D)
        Variable/G v_RayTraceHist_3D_is3D = use3D
        Variable/G v_RayTraceHist_3D_H_nm = H_nm
        Variable/G v_RayTraceHist_3D_maxPath_nm = maxPath_nm
        Variable/G v_RayTraceHist_3D_Bangle_deg = Bangle_deg
        Variable/G v_RayTraceHist_3D_STSx_nm = STSx
        Variable/G v_RayTraceHist_3D_STSy_nm = STSy
        Variable/G v_RayTraceHist_3D_nRay = nRay
        Variable/G v_RayTraceHist_3D_hasGeomW = hasGeomW
        Variable/G v_RayTraceHist_3D_hasT = hasT
        Variable/G v_RayTraceHist_3D_useWeights = hasWeights
        Variable/G v_RayTraceHist_3D_normalize = normalize
        Variable/G v_RayTraceHist_3D_absW = absW
        Variable/G v_RayTraceHist_3D_nBins = nBins
        Variable/G v_RayTraceHist_3D_binL_nm = binL_nm
        Variable/G v_RayTraceHist_3D_binW_nm = binW_nm
        Variable/G v_RayTraceHist_3D_binT = binT
        Variable/G v_RayTraceHist_3D_histLMin_nm = histLMin_nm
        Variable/G v_RayTraceHist_3D_histLMax_nm = histLMax_nm
        Variable/G v_RayTraceHist_3D_histWMin_nm = histWMin_nm
        Variable/G v_RayTraceHist_3D_histWMax_nm = histWMax_nm
        Variable/G v_RayTraceHist_3D_histTMin = histTMin
        Variable/G v_RayTraceHist_3D_histTMax = histTMax
        Variable/G v_RayTraceHist_3D_padCityscape = padCityscape
        Variable/G v_RayTraceHist_3D_plotMedian = plotMedian
        String/G s_RayTraceHist_3D_SourceDF = dfPath
        Variable/G v_RayTraceHist_3D_useExistingChannels = useExistingChannels
    else
        Variable/G v_RayTraceHist_is3D = use3D
        Variable/G v_RayTraceHist_H_nm = H_nm
        Variable/G v_RayTraceHist_maxPath_nm = maxPath_nm
        Variable/G v_RayTraceHist_Bangle_deg = Bangle_deg
        Variable/G v_RayTraceHist_STSx_nm = STSx
        Variable/G v_RayTraceHist_STSy_nm = STSy
        Variable/G v_RayTraceHist_nRay = nRay
        Variable/G v_RayTraceHist_hasGeomW = hasGeomW
        Variable/G v_RayTraceHist_hasT = hasT
        Variable/G v_RayTraceHist_useWeights = hasWeights
        Variable/G v_RayTraceHist_normalize = normalize
        Variable/G v_RayTraceHist_absW = absW
        Variable/G v_RayTraceHist_nBins = nBins
        Variable/G v_RayTraceHist_binL_nm = binL_nm
        Variable/G v_RayTraceHist_binW_nm = binW_nm
        Variable/G v_RayTraceHist_binT = binT
        Variable/G v_RayTraceHist_histLMin_nm = histLMin_nm
        Variable/G v_RayTraceHist_histLMax_nm = histLMax_nm
        Variable/G v_RayTraceHist_histWMin_nm = histWMin_nm
        Variable/G v_RayTraceHist_histWMax_nm = histWMax_nm
        Variable/G v_RayTraceHist_histTMin = histTMin
        Variable/G v_RayTraceHist_histTMax = histTMax
        Variable/G v_RayTraceHist_padCityscape = padCityscape
        Variable/G v_RayTraceHist_plotMedian = plotMedian
        String/G s_RayTraceHist_SourceDF = dfPath
        Variable/G v_RayTraceHist_useExistingChannels = useExistingChannels
    endif

    Variable minClusterSize = max(5, round(0.01*nRay))

    // -------------------------------------------------------------------------
    // 5. Prepare nm lists.
    // -------------------------------------------------------------------------
    Make/O/D/N=(nRay) tmp_L_N_List_nm_work
    Make/O/D/N=(nRay) tmp_W_eff_List_nm_work

    tmp_L_N_List_nm_work = L_work

    if (absW)
        tmp_W_eff_List_nm_work = abs(W_eff_work)
    else
        tmp_W_eff_List_nm_work = W_eff_work
    endif

    if (hasGeomW)
        Make/O/D/N=(nRay) tmp_W_geom_List_nm_work

        if (absW)
            tmp_W_geom_List_nm_work = abs(W_geom_work)
        else
            tmp_W_geom_List_nm_work = W_geom_work
        endif
    endif

    if (use3D)
        Duplicate/O tmp_L_N_List_nm_work, L_N_List_3D_nm
        Duplicate/O tmp_W_eff_List_nm_work, W_eff_List_3D_nm

        if (hasGeomW)
            Duplicate/O tmp_W_geom_List_nm_work, W_geom_List_3D_nm
        endif
    else
        Duplicate/O tmp_L_N_List_nm_work, L_N_List_nm
        Duplicate/O tmp_W_eff_List_nm_work, W_eff_List_nm
    endif

    Wave L_nm = tmp_L_N_List_nm_work
    Wave W_eff_nm = tmp_W_eff_List_nm_work
    Wave/Z W_geom_nm = tmp_W_geom_List_nm_work
    Variable i, j

    // Diagnostic cluster searches must follow the same spirit as continuum ray
    // tracing: contiguous means contiguous in angular ray-family order. Imported
    // weighted folders can be grouped by mode/projection instead, so build a
    // temporary phi modulo pi ordering for cluster diagnostics only. Stored
    // channel waves are left untouched.
    Duplicate/O L_nm, tmp_L_cluster_nm
    Duplicate/O W_eff_nm, tmp_W_eff_cluster_nm
    Make/O/D/N=(nRay) tmp_cluster_orig_index
    tmp_cluster_orig_index = p

    if (hasGeomW)
        Duplicate/O W_geom_nm, tmp_W_geom_cluster_nm
    endif

    Wave/Z phi_sort_2D = phiList_rad
    Wave/Z phi_sort_3D = phiList_3D_rad
    Variable hasPhiSort = (use3D && WaveExists(phi_sort_3D) && numpnts(phi_sort_3D) == nRay) || (!use3D && WaveExists(phi_sort_2D) && numpnts(phi_sort_2D) == nRay)

    if (hasPhiSort)
        Make/O/D/N=(nRay) tmp_phi_cluster_sort
        if (use3D)
            tmp_phi_cluster_sort = phi_sort_3D - Pi*floor(phi_sort_3D/Pi)
        else
            tmp_phi_cluster_sort = phi_sort_2D - Pi*floor(phi_sort_2D/Pi)
        endif

        if (hasGeomW)
            Sort tmp_phi_cluster_sort, tmp_phi_cluster_sort, tmp_L_cluster_nm, tmp_W_eff_cluster_nm, tmp_W_geom_cluster_nm, tmp_cluster_orig_index
        else
            Sort tmp_phi_cluster_sort, tmp_phi_cluster_sort, tmp_L_cluster_nm, tmp_W_eff_cluster_nm, tmp_cluster_orig_index
        endif

        // Put the artificial 0/pi boundary into the largest empty angular gap.
        // This keeps circular projective-angle clusters from being split by the
        // chosen modulo-pi branch cut.
        if (nRay > 1)
            Variable maxPhiGap = tmp_phi_cluster_sort[0] + Pi - tmp_phi_cluster_sort[nRay-1]
            Variable cutIndex = 0
            Variable phiGap

            for (i = 0; i < nRay-1; i += 1)
                phiGap = tmp_phi_cluster_sort[i+1] - tmp_phi_cluster_sort[i]
                if (phiGap > maxPhiGap)
                    maxPhiGap = phiGap
                    cutIndex = i + 1
                endif
            endfor

            if (cutIndex > 0)
                Duplicate/O tmp_phi_cluster_sort, tmp_phi_cluster_sort_buf
                Duplicate/O tmp_L_cluster_nm, tmp_L_cluster_nm_buf
                Duplicate/O tmp_W_eff_cluster_nm, tmp_W_eff_cluster_nm_buf
                Duplicate/O tmp_cluster_orig_index, tmp_cluster_orig_index_buf

                if (hasGeomW)
                    Duplicate/O tmp_W_geom_cluster_nm, tmp_W_geom_cluster_nm_buf
                endif

                Variable srcIdx
                for (i = 0; i < nRay; i += 1)
                    srcIdx = mod(i + cutIndex, nRay)
                    tmp_phi_cluster_sort[i] = tmp_phi_cluster_sort_buf[srcIdx]
                    tmp_L_cluster_nm[i] = tmp_L_cluster_nm_buf[srcIdx]
                    tmp_W_eff_cluster_nm[i] = tmp_W_eff_cluster_nm_buf[srcIdx]
                    tmp_cluster_orig_index[i] = tmp_cluster_orig_index_buf[srcIdx]

                    if (hasGeomW)
                        tmp_W_geom_cluster_nm[i] = tmp_W_geom_cluster_nm_buf[srcIdx]
                    endif
                endfor
            endif
        endif
    else
        Make/O/D/N=(nRay) tmp_phi_cluster_sort
        tmp_phi_cluster_sort = p
    endif

    Wave L_cluster_nm = tmp_L_cluster_nm
    Wave W_eff_cluster_nm = tmp_W_eff_cluster_nm
    Wave/Z W_geom_cluster_nm = tmp_W_geom_cluster_nm
    Wave clusterOrigIndex = tmp_cluster_orig_index

    // -------------------------------------------------------------------------
    // 6. Manual histogram binning for weighted and unweighted cases.
    //    Legacy mode: one nBins value per quantity, each using its own min/max.
    //    Fixed-bin mode: binL_nm and/or binW_nm impose common bin widths.
    // -------------------------------------------------------------------------
    WaveStats/Q L_nm
    Variable LminData = V_min
    Variable LmaxData = V_max

    WaveStats/Q W_eff_nm
    Variable WEffMinData = V_min
    Variable WEffMaxData = V_max

    Variable WGeomMinData = NaN
    Variable WGeomMaxData = NaN

    if (hasGeomW)
        WaveStats/Q W_geom_nm
        WGeomMinData = V_min
        WGeomMaxData = V_max
    endif

    Variable LminHist, LmaxHist, dLHist, nBinsL
    Variable WEffMinHist, WEffMaxHist, dWEffHist, nBinsWEff
    Variable WGeomMinHist, WGeomMaxHist, dWGeomHist, nBinsWGeom
    Variable TMinHist, TMaxHist, dTHist, nBinsT

    if (useFixedBinL)
        if (numtype(histLMin_nm) == 0)
            LminHist = histLMin_nm
        else
            LminHist = 0
        endif

        if (numtype(histLMax_nm) == 0)
            LmaxHist = histLMax_nm
        else
            LmaxHist = LminHist + ceil((LmaxData - LminHist) / binL_nm) * binL_nm
        endif

        dLHist = binL_nm
        nBinsL = max(2, ceil((LmaxHist - LminHist) / dLHist))
        LmaxHist = LminHist + nBinsL * dLHist
    else
        LminHist = LminData
        LmaxHist = LmaxData
        if (LmaxHist <= LminHist)
            LmaxHist = LminHist + 1
        endif
        nBinsL = nBins
        dLHist = (LmaxHist - LminHist) / nBinsL
    endif

    if (useFixedBinW)
        if (numtype(histWMin_nm) == 0)
            WEffMinHist = histWMin_nm
        else
            WEffMinHist = 0
        endif

        if (numtype(histWMax_nm) == 0)
            WEffMaxHist = histWMax_nm
        else
            WEffMaxHist = WEffMinHist + ceil((WEffMaxData - WEffMinHist) / binW_nm) * binW_nm
        endif

        dWEffHist = binW_nm
        nBinsWEff = max(2, ceil((WEffMaxHist - WEffMinHist) / dWEffHist))
        WEffMaxHist = WEffMinHist + nBinsWEff * dWEffHist
    else
        WEffMinHist = WEffMinData
        WEffMaxHist = WEffMaxData
        if (WEffMaxHist <= WEffMinHist)
            WEffMaxHist = WEffMinHist + 1
        endif
        nBinsWEff = nBins
        dWEffHist = (WEffMaxHist - WEffMinHist) / nBinsWEff
    endif

    WGeomMinHist = NaN
    WGeomMaxHist = NaN
    dWGeomHist = NaN
    nBinsWGeom = 0

    if (hasGeomW)
        if (useFixedBinW)
            if (numtype(histWMin_nm) == 0)
                WGeomMinHist = histWMin_nm
            else
                WGeomMinHist = 0
            endif

            if (numtype(histWMax_nm) == 0)
                WGeomMaxHist = histWMax_nm
            else
                WGeomMaxHist = WGeomMinHist + ceil((WGeomMaxData - WGeomMinHist) / binW_nm) * binW_nm
            endif

            dWGeomHist = binW_nm
            nBinsWGeom = max(2, ceil((WGeomMaxHist - WGeomMinHist) / dWGeomHist))
            WGeomMaxHist = WGeomMinHist + nBinsWGeom * dWGeomHist
        else
            WGeomMinHist = WGeomMinData
            WGeomMaxHist = WGeomMaxData
            if (WGeomMaxHist <= WGeomMinHist)
                WGeomMaxHist = WGeomMinHist + 1
            endif
            nBinsWGeom = nBins
            dWGeomHist = (WGeomMaxHist - WGeomMinHist) / nBinsWGeom
        endif
    endif

    TMinHist = histTMin
    TMaxHist = histTMax
    if (numtype(TMinHist) != 0)
        TMinHist = 0
    endif
    if (numtype(TMaxHist) != 0)
        TMaxHist = 1
    endif
    if (TMaxHist <= TMinHist)
        TMaxHist = TMinHist + 1
    endif

    dTHist = binT
    nBinsT = max(2, ceil((TMaxHist - TMinHist) / dTHist))
    TMaxHist = TMinHist + nBinsT * dTHist

    Make/O/D/N=(nBinsL) tmp_Hist_L_work
    Make/O/D/N=(nBinsWEff) tmp_Hist_W_eff_work

    if (hasGeomW)
        Make/O/D/N=(nBinsWGeom) tmp_Hist_W_geom_work
    endif

    if (hasT)
        Make/O/D/N=(nBinsT) tmp_Hist_T_eff_work
    endif

    SetScale/P x, LminHist + 0.5*dLHist, dLHist, "nm", tmp_Hist_L_work
    SetScale/P x, WEffMinHist + 0.5*dWEffHist, dWEffHist, "nm", tmp_Hist_W_eff_work

    if (hasGeomW)
        SetScale/P x, WGeomMinHist + 0.5*dWGeomHist, dWGeomHist, "nm", tmp_Hist_W_geom_work
    endif

    if (hasT)
        SetScale/P x, TMinHist + 0.5*dTHist, dTHist, "", tmp_Hist_T_eff_work
    endif

    tmp_Hist_L_work = 0
    tmp_Hist_W_eff_work = 0

    if (hasGeomW)
        tmp_Hist_W_geom_work = 0
    endif

    if (hasT)
        tmp_Hist_T_eff_work = 0
    endif

    Variable iBin
    Variable wt
    Variable tVal

    for (i = 0; i < nRay; i += 1)

        if (hasWeights)
            wt = Wt_work[i]
            if (numtype(wt) != 0 || wt < 0)
                wt = 0
            endif
        else
            wt = 1
        endif

        if (L_nm[i] >= LminHist && L_nm[i] <= LmaxHist)
            iBin = floor((L_nm[i] - LminHist) / dLHist)
            iBin = max(0, min(nBinsL-1, iBin))
            tmp_Hist_L_work[iBin] += wt
        endif

        if (W_eff_nm[i] >= WEffMinHist && W_eff_nm[i] <= WEffMaxHist)
            iBin = floor((W_eff_nm[i] - WEffMinHist) / dWEffHist)
            iBin = max(0, min(nBinsWEff-1, iBin))
            tmp_Hist_W_eff_work[iBin] += wt
        endif

        if (hasGeomW)
            if (W_geom_nm[i] >= WGeomMinHist && W_geom_nm[i] <= WGeomMaxHist)
                iBin = floor((W_geom_nm[i] - WGeomMinHist) / dWGeomHist)
                iBin = max(0, min(nBinsWGeom-1, iBin))
                tmp_Hist_W_geom_work[iBin] += wt
            endif
        endif

        if (hasT)
            tVal = T_work[i]
            if (numtype(tVal) == 0 && tVal >= TMinHist && tVal <= TMaxHist)
                iBin = floor((tVal - TMinHist) / dTHist)
                iBin = max(0, min(nBinsT-1, iBin))
                tmp_Hist_T_eff_work[iBin] += wt
            endif
        endif
    endfor

    if (use3D)
        Variable/G v_RayTraceHist_3D_nBinsL = nBinsL
        Variable/G v_RayTraceHist_3D_nBinsW_eff = nBinsWEff
        Variable/G v_RayTraceHist_3D_nBinsW_geom = nBinsWGeom
        Variable/G v_RayTraceHist_3D_nBinsT_eff = nBinsT
        Variable/G v_RayTraceHist_3D_dLHist_nm = dLHist
        Variable/G v_RayTraceHist_3D_dWHist_nm = dWEffHist
        Variable/G v_RayTraceHist_3D_dTHist = dTHist
        Variable/G v_RayTraceHist_3D_LminHist_nm = LminHist
        Variable/G v_RayTraceHist_3D_LmaxHist_nm = LmaxHist
        Variable/G v_RayTraceHist_3D_WminHist_nm = WEffMinHist
        Variable/G v_RayTraceHist_3D_WmaxHist_nm = WEffMaxHist
        Variable/G v_RayTraceHist_3D_TminHist = TMinHist
        Variable/G v_RayTraceHist_3D_TmaxHist = TMaxHist
    else
        Variable/G v_RayTraceHist_nBinsL = nBinsL
        Variable/G v_RayTraceHist_nBinsW_eff = nBinsWEff
        Variable/G v_RayTraceHist_nBinsW_geom = nBinsWGeom
        Variable/G v_RayTraceHist_nBinsT_eff = nBinsT
        Variable/G v_RayTraceHist_dLHist_nm = dLHist
        Variable/G v_RayTraceHist_dWHist_nm = dWEffHist
        Variable/G v_RayTraceHist_dTHist = dTHist
        Variable/G v_RayTraceHist_LminHist_nm = LminHist
        Variable/G v_RayTraceHist_LmaxHist_nm = LmaxHist
        Variable/G v_RayTraceHist_WminHist_nm = WEffMinHist
        Variable/G v_RayTraceHist_WmaxHist_nm = WEffMaxHist
        Variable/G v_RayTraceHist_TminHist = TMinHist
        Variable/G v_RayTraceHist_TmaxHist = TMaxHist
    endif

    // -------------------------------------------------------------------------
    // 7. Optional normalization to probability per bin in percent.
    // -------------------------------------------------------------------------
    if (normalize)
        Variable sumL = sum(tmp_Hist_L_work)
        Variable sumWEff = sum(tmp_Hist_W_eff_work)

        if (sumL > 0)
            tmp_Hist_L_work = 100 * tmp_Hist_L_work / sumL
        endif

        if (sumWEff > 0)
            tmp_Hist_W_eff_work = 100 * tmp_Hist_W_eff_work / sumWEff
        endif

        if (hasGeomW)
            Variable sumWGeom = sum(tmp_Hist_W_geom_work)

            if (sumWGeom > 0)
                tmp_Hist_W_geom_work = 100 * tmp_Hist_W_geom_work / sumWGeom
            endif
        endif

        if (hasT)
            Variable sumT = sum(tmp_Hist_T_eff_work)

            if (sumT > 0)
                tmp_Hist_T_eff_work = 100 * tmp_Hist_T_eff_work / sumT
            endif
        endif
    endif

    String histLName, histWEffName, histWGeomName, histTName

    if (use3D)
        if (hasWeights)
            Duplicate/O tmp_Hist_L_work, Hist_L_N_3D_weighted_nm
            Duplicate/O tmp_Hist_W_eff_work, Hist_W_eff_3D_weighted_nm
            histLName = "Hist_L_N_3D_weighted_nm"
            histWEffName = "Hist_W_eff_3D_weighted_nm"

            if (hasT)
                Duplicate/O tmp_Hist_T_eff_work, Hist_T_eff_3D_weighted
                histTName = "Hist_T_eff_3D_weighted"
            else
                histTName = ""
            endif

            if (hasGeomW)
                Duplicate/O tmp_Hist_W_geom_work, Hist_W_geom_3D_weighted_nm
                histWGeomName = "Hist_W_geom_3D_weighted_nm"
            else
                histWGeomName = ""
            endif
        else
            Duplicate/O tmp_Hist_L_work, Hist_L_N_3D_nm
            Duplicate/O tmp_Hist_W_eff_work, Hist_W_eff_3D_nm
            histLName = "Hist_L_N_3D_nm"
            histWEffName = "Hist_W_eff_3D_nm"

            if (hasT)
                Duplicate/O tmp_Hist_T_eff_work, Hist_T_eff_3D
                histTName = "Hist_T_eff_3D"
            else
                histTName = ""
            endif

            if (hasGeomW)
                Duplicate/O tmp_Hist_W_geom_work, Hist_W_geom_3D_nm
                histWGeomName = "Hist_W_geom_3D_nm"
            else
                histWGeomName = ""
            endif
        endif
    else
        if (hasWeights)
            Duplicate/O tmp_Hist_L_work, Hist_L_N_weighted_nm
            Duplicate/O tmp_Hist_W_eff_work, Hist_W_eff_weighted_nm
            histLName = "Hist_L_N_weighted_nm"
            histWEffName = "Hist_W_eff_weighted_nm"
            histWGeomName = ""

            if (hasT)
                Duplicate/O tmp_Hist_T_eff_work, Hist_T_eff_weighted
                histTName = "Hist_T_eff_weighted"
            else
                histTName = ""
            endif
        else
            Duplicate/O tmp_Hist_L_work, Hist_L_N_nm
            Duplicate/O tmp_Hist_W_eff_work, Hist_W_eff_nm
            histLName = "Hist_L_N_nm"
            histWEffName = "Hist_W_eff_nm"
            histWGeomName = ""

            if (hasT)
                Duplicate/O tmp_Hist_T_eff_work, Hist_T_eff
                histTName = "Hist_T_eff"
            else
                histTName = ""
            endif
        endif
    endif

    if (use3D)
        String/G s_RayTraceHist_3D_HistL = histLName
        String/G s_RayTraceHist_3D_HistW = histWEffName
        String/G s_RayTraceHist_3D_HistW_eff = histWEffName
        String/G s_RayTraceHist_3D_HistW_geom = histWGeomName
        String/G s_RayTraceHist_3D_HistT = histTName
        String/G s_RayTraceHist_3D_HistT_eff = histTName
    else
        String/G s_RayTraceHist_HistL = histLName
        String/G s_RayTraceHist_HistW = histWEffName
        String/G s_RayTraceHist_HistW_eff = histWEffName
        String/G s_RayTraceHist_HistW_geom = histWGeomName
        String/G s_RayTraceHist_HistT = histTName
        String/G s_RayTraceHist_HistT_eff = histTName
    endif

    // -------------------------------------------------------------------------
    // 8. Supported upper-tail estimator for L_N.
    // -------------------------------------------------------------------------
    WaveStats/Q L_nm
    Variable LN_global_max_nm = V_max
    Variable LN_low_nm = 0.9 * LN_global_max_nm
    Variable LN_high_nm = LN_global_max_nm

    Variable LN_extreme_n = 0
    Variable LN_extreme_max_nm = NaN
    Variable LN_extreme_median_nm = NaN
    Variable LN_extreme_sigma_nm = NaN
    Variable LN_extreme_max_trace = NaN
    Variable LN_extreme_median_trace = NaN
    Variable LN_extreme_cluster_start = NaN
    Variable LN_extreme_cluster_size = NaN

    Variable curStart, curN
    Variable bestStartL = -1
    Variable bestNL = 0
    Variable bestMedianL = -Inf
    Variable testMedianL
    Variable nHalf

    Make/O/D/N=0 tmp_L_N_extreme_work_nm

    i = 0
    do
        if (i >= nRay)
            break
        endif

        if (L_cluster_nm[i] >= LN_low_nm && L_cluster_nm[i] <= LN_high_nm)
            curStart = i
            curN = 0

            do
                if (i >= nRay)
                    break
                endif
                if (!(L_cluster_nm[i] >= LN_low_nm && L_cluster_nm[i] <= LN_high_nm))
                    break
                endif
                curN += 1
                i += 1
            while (1)

            if (curN >= minClusterSize)
                Make/O/D/N=(curN) tmp_L_cluster_vals

                for (j = 0; j < curN; j += 1)
                    tmp_L_cluster_vals[j] = L_cluster_nm[curStart+j]
                endfor

                Sort tmp_L_cluster_vals, tmp_L_cluster_vals

                nHalf = floor(curN/2)
                if (mod(curN, 2) == 1)
                    testMedianL = tmp_L_cluster_vals[nHalf]
                else
                    testMedianL = 0.5*(tmp_L_cluster_vals[nHalf-1] + tmp_L_cluster_vals[nHalf])
                endif

                if ((curN > bestNL) || ((curN == bestNL) && (testMedianL > bestMedianL)))
                    bestNL = curN
                    bestStartL = curStart
                    bestMedianL = testMedianL
                endif
            endif
        else
            i += 1
        endif
    while (1)

    if (bestStartL >= 0)
        Make/O/D/N=(bestNL) tmp_L_N_extreme_work_nm

        for (j = 0; j < bestNL; j += 1)
            tmp_L_N_extreme_work_nm[j] = L_cluster_nm[bestStartL+j]
        endfor

        Duplicate/O tmp_L_N_extreme_work_nm, tmp_L_N_extreme_sorted
        Sort tmp_L_N_extreme_sorted, tmp_L_N_extreme_sorted

        nHalf = floor(bestNL/2)
        if (mod(bestNL, 2) == 1)
            LN_extreme_median_nm = tmp_L_N_extreme_sorted[nHalf]
        else
            LN_extreme_median_nm = 0.5*(tmp_L_N_extreme_sorted[nHalf-1] + tmp_L_N_extreme_sorted[nHalf])
        endif

        WaveStats/Q tmp_L_N_extreme_work_nm
        LN_extreme_max_nm = V_max

        if (bestNL > 1)
            LN_extreme_sigma_nm = V_sdev
        endif

        LN_extreme_n = bestNL
        LN_extreme_cluster_start = bestStartL
        LN_extreme_cluster_size = bestNL

        Variable bestDMaxL = Inf
        Variable bestDMedL = Inf
        Variable thisD

        for (i = 0; i < nRay; i += 1)
            thisD = abs(L_cluster_nm[i] - LN_extreme_max_nm)
            if (thisD < bestDMaxL)
                bestDMaxL = thisD
                LN_extreme_max_trace = clusterOrigIndex[i]
            endif

            thisD = abs(L_cluster_nm[i] - LN_extreme_median_nm)
            if (thisD < bestDMedL)
                bestDMedL = thisD
                LN_extreme_median_trace = clusterOrigIndex[i]
            endif
        endfor
    endif

    if (use3D)
        Duplicate/O tmp_L_N_extreme_work_nm, L_N_extreme_3D_nm

        Variable/G v_L_3D_global_max_nm = LN_global_max_nm
        Variable/G v_L_3D_tailCluster_low_nm = LN_low_nm
        Variable/G v_L_3D_tailCluster_high_nm = LN_high_nm
        Variable/G v_L_3D_tailCluster_n = LN_extreme_n
        Variable/G v_L_3D_tailCluster_max_nm = LN_extreme_max_nm
        Variable/G v_L_3D_tailCluster_median_nm = LN_extreme_median_nm
        Variable/G v_L_3D_tailCluster_sigma_nm = LN_extreme_sigma_nm
        Variable/G v_L_3D_tailCluster_max_trace = LN_extreme_max_trace
        Variable/G v_L_3D_tailCluster_median_trace = LN_extreme_median_trace
        Variable/G v_L_3D_tailCluster_start = LN_extreme_cluster_start
        Variable/G v_L_3D_tailCluster_size = LN_extreme_cluster_size
        Variable/G v_L_3D_tailCluster_minSize = minClusterSize
    else
        Duplicate/O tmp_L_N_extreme_work_nm, L_N_extreme_nm

        Variable/G v_L_global_max_nm = LN_global_max_nm
        Variable/G v_L_tailCluster_low_nm = LN_low_nm
        Variable/G v_L_tailCluster_high_nm = LN_high_nm
        Variable/G v_L_tailCluster_n = LN_extreme_n
        Variable/G v_L_tailCluster_max_nm = LN_extreme_max_nm
        Variable/G v_L_tailCluster_median_nm = LN_extreme_median_nm
        Variable/G v_L_tailCluster_sigma_nm = LN_extreme_sigma_nm
        Variable/G v_L_tailCluster_max_trace = LN_extreme_max_trace
        Variable/G v_L_tailCluster_median_trace = LN_extreme_median_trace
        Variable/G v_L_tailCluster_start = LN_extreme_cluster_start
        Variable/G v_L_tailCluster_size = LN_extreme_cluster_size
        Variable/G v_L_tailCluster_minSize = minClusterSize
    endif

    // -------------------------------------------------------------------------
    // 8b. Supported smallest-bin estimator for L_N.
    //
    // Uses the actual histogram binning above. Therefore, with binL_nm=20,
    // this finds the first occupied 20 nm ℓ-bin.
    //
    // Within the smallest occupied ℓ-bin:
    //   1. rays are selected by bin membership
    //   2. selected rays are split into contiguous angular-index clusters
    //   3. clusters with size < minClusterSize are rejected
    //   4. dominant cluster is largest population; tie-breaker = smaller median ℓ
    //   5. representative ℓ is the max member of that accepted small-ℓ cluster
    // -------------------------------------------------------------------------
    Variable LN_minbin_low_nm = NaN
    Variable LN_minbin_high_nm = NaN

    Variable LN_minbin_n = 0
    Variable LN_minbin_max_nm = NaN
    Variable LN_minbin_median_nm = NaN
    Variable LN_minbin_sigma_nm = NaN
    Variable LN_minbin_max_trace = NaN
    Variable LN_minbin_median_trace = NaN
    Variable LN_minbin_cluster_start = NaN
    Variable LN_minbin_cluster_size = NaN

    Variable iBinMinL = -1
    Variable binLo, binHi
    Variable inSmallBin
    Variable hasRayInBin

    Make/O/D/N=0 tmp_L_N_minbin_work_nm

    // Find the first occupied L histogram bin.
    for (iBin = 0; iBin < nBinsL; iBin += 1)
        binLo = LminHist + iBin*dLHist
        binHi = binLo + dLHist
        hasRayInBin = 0

        for (i = 0; i < nRay; i += 1)
            if (iBin == nBinsL-1)
                inSmallBin = (L_cluster_nm[i] >= binLo && L_cluster_nm[i] <= binHi)
            else
                inSmallBin = (L_cluster_nm[i] >= binLo && L_cluster_nm[i] < binHi)
            endif

            if (inSmallBin)
                hasRayInBin = 1
                break
            endif
        endfor

        if (hasRayInBin)
            iBinMinL = iBin
            LN_minbin_low_nm = binLo
            LN_minbin_high_nm = binHi
            break
        endif
    endfor

    Variable bestStartLmin = -1
    Variable bestNLmin = 0
    Variable bestMedianLmin = Inf
    Variable testMedianLmin

    if (iBinMinL >= 0)

        i = 0
        do
            if (i >= nRay)
                break
            endif

            if (iBinMinL == nBinsL-1)
                inSmallBin = (L_cluster_nm[i] >= LN_minbin_low_nm && L_cluster_nm[i] <= LN_minbin_high_nm)
            else
                inSmallBin = (L_cluster_nm[i] >= LN_minbin_low_nm && L_cluster_nm[i] < LN_minbin_high_nm)
            endif

            if (inSmallBin)
                curStart = i
                curN = 0

                do
                    if (i >= nRay)
                        break
                    endif

                    if (iBinMinL == nBinsL-1)
                        inSmallBin = (L_cluster_nm[i] >= LN_minbin_low_nm && L_cluster_nm[i] <= LN_minbin_high_nm)
                    else
                        inSmallBin = (L_cluster_nm[i] >= LN_minbin_low_nm && L_cluster_nm[i] < LN_minbin_high_nm)
                    endif

                    if (!inSmallBin)
                        break
                    endif

                    curN += 1
                    i += 1
                while (1)

                if (curN >= minClusterSize)
                    Make/O/D/N=(curN) tmp_L_minbin_vals

                    for (j = 0; j < curN; j += 1)
                        tmp_L_minbin_vals[j] = L_cluster_nm[curStart+j]
                    endfor

                    Sort tmp_L_minbin_vals, tmp_L_minbin_vals

                    nHalf = floor(curN/2)
                    if (mod(curN, 2) == 1)
                        testMedianLmin = tmp_L_minbin_vals[nHalf]
                    else
                        testMedianLmin = 0.5*(tmp_L_minbin_vals[nHalf-1] + tmp_L_minbin_vals[nHalf])
                    endif

                    if ((curN > bestNLmin) || ((curN == bestNLmin) && (testMedianLmin < bestMedianLmin)))
                        bestNLmin = curN
                        bestStartLmin = curStart
                        bestMedianLmin = testMedianLmin
                    endif
                endif
            else
                i += 1
            endif
        while (1)
    endif

    if (bestStartLmin >= 0)
        Make/O/D/N=(bestNLmin) tmp_L_N_minbin_work_nm

        for (j = 0; j < bestNLmin; j += 1)
            tmp_L_N_minbin_work_nm[j] = L_cluster_nm[bestStartLmin+j]
        endfor

        Duplicate/O tmp_L_N_minbin_work_nm, tmp_L_N_minbin_sorted
        Sort tmp_L_N_minbin_sorted, tmp_L_N_minbin_sorted

        nHalf = floor(bestNLmin/2)
        if (mod(bestNLmin, 2) == 1)
            LN_minbin_median_nm = tmp_L_N_minbin_sorted[nHalf]
        else
            LN_minbin_median_nm = 0.5*(tmp_L_N_minbin_sorted[nHalf-1] + tmp_L_N_minbin_sorted[nHalf])
        endif

        WaveStats/Q tmp_L_N_minbin_work_nm
        LN_minbin_max_nm = V_max

        if (bestNLmin > 1)
            LN_minbin_sigma_nm = V_sdev
        endif

        LN_minbin_n = bestNLmin
        LN_minbin_cluster_start = bestStartLmin
        LN_minbin_cluster_size = bestNLmin

        Variable bestDMaxLmin = Inf
        Variable bestDMedLmin = Inf

        for (i = 0; i < nRay; i += 1)
            thisD = abs(L_cluster_nm[i] - LN_minbin_max_nm)
            if (thisD < bestDMaxLmin)
                bestDMaxLmin = thisD
                LN_minbin_max_trace = clusterOrigIndex[i]
            endif

            thisD = abs(L_cluster_nm[i] - LN_minbin_median_nm)
            if (thisD < bestDMedLmin)
                bestDMedLmin = thisD
                LN_minbin_median_trace = clusterOrigIndex[i]
            endif
        endfor
    endif

    if (use3D)
        Duplicate/O tmp_L_N_minbin_work_nm, L_N_minbin_3D_nm

        Variable/G v_L_3D_minbinCluster_low_nm = LN_minbin_low_nm
        Variable/G v_L_3D_minbinCluster_high_nm = LN_minbin_high_nm
        Variable/G v_L_3D_minbinCluster_n = LN_minbin_n
        Variable/G v_L_3D_minbinCluster_max_nm = LN_minbin_max_nm
        Variable/G v_L_3D_minbinCluster_median_nm = LN_minbin_median_nm
        Variable/G v_L_3D_minbinCluster_sigma_nm = LN_minbin_sigma_nm
        Variable/G v_L_3D_minbinCluster_max_trace = LN_minbin_max_trace
        Variable/G v_L_3D_minbinCluster_median_trace = LN_minbin_median_trace
        Variable/G v_L_3D_minbinCluster_start = LN_minbin_cluster_start
        Variable/G v_L_3D_minbinCluster_size = LN_minbin_cluster_size
        Variable/G v_L_3D_minbinCluster_minSize = minClusterSize
    else
        Duplicate/O tmp_L_N_minbin_work_nm, L_N_minbin_nm

        Variable/G v_L_minbinCluster_low_nm = LN_minbin_low_nm
        Variable/G v_L_minbinCluster_high_nm = LN_minbin_high_nm
        Variable/G v_L_minbinCluster_n = LN_minbin_n
        Variable/G v_L_minbinCluster_max_nm = LN_minbin_max_nm
        Variable/G v_L_minbinCluster_median_nm = LN_minbin_median_nm
        Variable/G v_L_minbinCluster_sigma_nm = LN_minbin_sigma_nm
        Variable/G v_L_minbinCluster_max_trace = LN_minbin_max_trace
        Variable/G v_L_minbinCluster_median_trace = LN_minbin_median_trace
        Variable/G v_L_minbinCluster_start = LN_minbin_cluster_start
        Variable/G v_L_minbinCluster_size = LN_minbin_cluster_size
        Variable/G v_L_minbinCluster_minSize = minClusterSize
    endif

    // -------------------------------------------------------------------------
    // 9. Supported upper-tail estimator for W_eff.
    // -------------------------------------------------------------------------
    WaveStats/Q W_eff_nm
    Variable WEff_global_max_nm = V_max
    Variable WEff_low_nm = 0.9 * WEff_global_max_nm
    Variable WEff_high_nm = WEff_global_max_nm

    Variable WEff_extreme_n = 0
    Variable WEff_extreme_max_nm = NaN
    Variable WEff_extreme_median_nm = NaN
    Variable WEff_extreme_sigma_nm = NaN
    Variable WEff_extreme_max_trace = NaN
    Variable WEff_extreme_median_trace = NaN
    Variable WEff_extreme_cluster_start = NaN
    Variable WEff_extreme_cluster_size = NaN

    Variable bestStartWEff = -1
    Variable bestNWEff = 0
    Variable bestMedianWEff = -Inf
    Variable testMedianWEff

    Make/O/D/N=0 tmp_W_eff_extreme_work_nm

    i = 0
    do
        if (i >= nRay)
            break
        endif

        if (W_eff_cluster_nm[i] >= WEff_low_nm && W_eff_cluster_nm[i] <= WEff_high_nm)
            curStart = i
            curN = 0

            do
                if (i >= nRay)
                    break
                endif
                if (!(W_eff_cluster_nm[i] >= WEff_low_nm && W_eff_cluster_nm[i] <= WEff_high_nm))
                    break
                endif
                curN += 1
                i += 1
            while (1)

            if (curN >= minClusterSize)
                Make/O/D/N=(curN) tmp_WEff_cluster_vals

                for (j = 0; j < curN; j += 1)
                    tmp_WEff_cluster_vals[j] = W_eff_cluster_nm[curStart+j]
                endfor

                Sort tmp_WEff_cluster_vals, tmp_WEff_cluster_vals

                nHalf = floor(curN/2)
                if (mod(curN, 2) == 1)
                    testMedianWEff = tmp_WEff_cluster_vals[nHalf]
                else
                    testMedianWEff = 0.5*(tmp_WEff_cluster_vals[nHalf-1] + tmp_WEff_cluster_vals[nHalf])
                endif

                if ((curN > bestNWEff) || ((curN == bestNWEff) && (testMedianWEff > bestMedianWEff)))
                    bestNWEff = curN
                    bestStartWEff = curStart
                    bestMedianWEff = testMedianWEff
                endif
            endif
        else
            i += 1
        endif
    while (1)

    if (bestStartWEff >= 0)
        Make/O/D/N=(bestNWEff) tmp_W_eff_extreme_work_nm

        for (j = 0; j < bestNWEff; j += 1)
            tmp_W_eff_extreme_work_nm[j] = W_eff_cluster_nm[bestStartWEff+j]
        endfor

        Duplicate/O tmp_W_eff_extreme_work_nm, tmp_W_eff_extreme_sorted
        Sort tmp_W_eff_extreme_sorted, tmp_W_eff_extreme_sorted

        nHalf = floor(bestNWEff/2)
        if (mod(bestNWEff, 2) == 1)
            WEff_extreme_median_nm = tmp_W_eff_extreme_sorted[nHalf]
        else
            WEff_extreme_median_nm = 0.5*(tmp_W_eff_extreme_sorted[nHalf-1] + tmp_W_eff_extreme_sorted[nHalf])
        endif

        WaveStats/Q tmp_W_eff_extreme_work_nm
        WEff_extreme_max_nm = V_max

        if (bestNWEff > 1)
            WEff_extreme_sigma_nm = V_sdev
        endif

        WEff_extreme_n = bestNWEff
        WEff_extreme_cluster_start = bestStartWEff
        WEff_extreme_cluster_size = bestNWEff

        Variable bestDMaxWEff = Inf
        Variable bestDMedWEff = Inf

        for (i = 0; i < nRay; i += 1)
            thisD = abs(W_eff_cluster_nm[i] - WEff_extreme_max_nm)
            if (thisD < bestDMaxWEff)
                bestDMaxWEff = thisD
                WEff_extreme_max_trace = clusterOrigIndex[i]
            endif

            thisD = abs(W_eff_cluster_nm[i] - WEff_extreme_median_nm)
            if (thisD < bestDMedWEff)
                bestDMedWEff = thisD
                WEff_extreme_median_trace = clusterOrigIndex[i]
            endif
        endfor
    endif

    if (use3D)
        Duplicate/O tmp_W_eff_extreme_work_nm, W_eff_extreme_3D_nm

        Variable/G v_W_3D_global_max_nm = WEff_global_max_nm
        Variable/G v_W_3D_tailCluster_low_nm = WEff_low_nm
        Variable/G v_W_3D_tailCluster_high_nm = WEff_high_nm
        Variable/G v_W_3D_tailCluster_n = WEff_extreme_n
        Variable/G v_W_3D_tailCluster_max_nm = WEff_extreme_max_nm
        Variable/G v_W_3D_tailCluster_median_nm = WEff_extreme_median_nm
        Variable/G v_W_3D_tailCluster_sigma_nm = WEff_extreme_sigma_nm
        Variable/G v_W_3D_tailCluster_max_trace = WEff_extreme_max_trace
        Variable/G v_W_3D_tailCluster_median_trace = WEff_extreme_median_trace
        Variable/G v_W_3D_tailCluster_start = WEff_extreme_cluster_start
        Variable/G v_W_3D_tailCluster_size = WEff_extreme_cluster_size
        Variable/G v_W_3D_tailCluster_minSize = minClusterSize
    else
        Duplicate/O tmp_W_eff_extreme_work_nm, W_eff_extreme_nm

        Variable/G v_W_global_max_nm = WEff_global_max_nm
        Variable/G v_W_tailCluster_low_nm = WEff_low_nm
        Variable/G v_W_tailCluster_high_nm = WEff_high_nm
        Variable/G v_W_tailCluster_n = WEff_extreme_n
        Variable/G v_W_tailCluster_max_nm = WEff_extreme_max_nm
        Variable/G v_W_tailCluster_median_nm = WEff_extreme_median_nm
        Variable/G v_W_tailCluster_sigma_nm = WEff_extreme_sigma_nm
        Variable/G v_W_tailCluster_max_trace = WEff_extreme_max_trace
        Variable/G v_W_tailCluster_median_trace = WEff_extreme_median_trace
        Variable/G v_W_tailCluster_start = WEff_extreme_cluster_start
        Variable/G v_W_tailCluster_size = WEff_extreme_cluster_size
        Variable/G v_W_tailCluster_minSize = minClusterSize
    endif

    // -------------------------------------------------------------------------
    // 10. Supported upper-tail estimator for W_geom, if available.
    // -------------------------------------------------------------------------
    Variable WGeom_global_max_nm = NaN
    Variable WGeom_low_nm = NaN
    Variable WGeom_high_nm = NaN

    Variable WGeom_extreme_n = 0
    Variable WGeom_extreme_max_nm = NaN
    Variable WGeom_extreme_median_nm = NaN
    Variable WGeom_extreme_sigma_nm = NaN
    Variable WGeom_extreme_max_trace = NaN
    Variable WGeom_extreme_median_trace = NaN
    Variable WGeom_extreme_cluster_start = NaN
    Variable WGeom_extreme_cluster_size = NaN

    Make/O/D/N=0 tmp_W_geom_extreme_work_nm

    if (hasGeomW)

        WaveStats/Q W_geom_nm
        WGeom_global_max_nm = V_max
        WGeom_low_nm = 0.9 * WGeom_global_max_nm
        WGeom_high_nm = WGeom_global_max_nm

        Variable bestStartWGeom = -1
        Variable bestNWGeom = 0
        Variable bestMedianWGeom = -Inf
        Variable testMedianWGeom

        i = 0
        do
            if (i >= nRay)
                break
            endif

            if (W_geom_cluster_nm[i] >= WGeom_low_nm && W_geom_cluster_nm[i] <= WGeom_high_nm)
                curStart = i
                curN = 0

                do
                    if (i >= nRay)
                        break
                    endif
                    if (!(W_geom_cluster_nm[i] >= WGeom_low_nm && W_geom_cluster_nm[i] <= WGeom_high_nm))
                        break
                    endif
                    curN += 1
                    i += 1
                while (1)

                if (curN >= minClusterSize)
                    Make/O/D/N=(curN) tmp_WGeom_cluster_vals

                    for (j = 0; j < curN; j += 1)
                        tmp_WGeom_cluster_vals[j] = W_geom_cluster_nm[curStart+j]
                    endfor

                    Sort tmp_WGeom_cluster_vals, tmp_WGeom_cluster_vals

                    nHalf = floor(curN/2)
                    if (mod(curN, 2) == 1)
                        testMedianWGeom = tmp_WGeom_cluster_vals[nHalf]
                    else
                        testMedianWGeom = 0.5*(tmp_WGeom_cluster_vals[nHalf-1] + tmp_WGeom_cluster_vals[nHalf])
                    endif

                    if ((curN > bestNWGeom) || ((curN == bestNWGeom) && (testMedianWGeom > bestMedianWGeom)))
                        bestNWGeom = curN
                        bestStartWGeom = curStart
                        bestMedianWGeom = testMedianWGeom
                    endif
                endif
            else
                i += 1
            endif
        while (1)

        if (bestStartWGeom >= 0)
            Make/O/D/N=(bestNWGeom) tmp_W_geom_extreme_work_nm

            for (j = 0; j < bestNWGeom; j += 1)
                tmp_W_geom_extreme_work_nm[j] = W_geom_cluster_nm[bestStartWGeom+j]
            endfor

            Duplicate/O tmp_W_geom_extreme_work_nm, tmp_W_geom_extreme_sorted
            Sort tmp_W_geom_extreme_sorted, tmp_W_geom_extreme_sorted

            nHalf = floor(bestNWGeom/2)
            if (mod(bestNWGeom, 2) == 1)
                WGeom_extreme_median_nm = tmp_W_geom_extreme_sorted[nHalf]
            else
                WGeom_extreme_median_nm = 0.5*(tmp_W_geom_extreme_sorted[nHalf-1] + tmp_W_geom_extreme_sorted[nHalf])
            endif

            WaveStats/Q tmp_W_geom_extreme_work_nm
            WGeom_extreme_max_nm = V_max

            if (bestNWGeom > 1)
                WGeom_extreme_sigma_nm = V_sdev
            endif

            WGeom_extreme_n = bestNWGeom
            WGeom_extreme_cluster_start = bestStartWGeom
            WGeom_extreme_cluster_size = bestNWGeom

            Variable bestDMaxWGeom = Inf
            Variable bestDMedWGeom = Inf

            for (i = 0; i < nRay; i += 1)
                thisD = abs(W_geom_cluster_nm[i] - WGeom_extreme_max_nm)
                if (thisD < bestDMaxWGeom)
                    bestDMaxWGeom = thisD
                    WGeom_extreme_max_trace = clusterOrigIndex[i]
                endif

                thisD = abs(W_geom_cluster_nm[i] - WGeom_extreme_median_nm)
                if (thisD < bestDMedWGeom)
                    bestDMedWGeom = thisD
                    WGeom_extreme_median_trace = clusterOrigIndex[i]
                endif
            endfor
        endif

        Duplicate/O tmp_W_geom_extreme_work_nm, W_geom_extreme_3D_nm

        Variable/G v_Wgeom_3D_global_max_nm = WGeom_global_max_nm
        Variable/G v_Wgeom_3D_tailCluster_low_nm = WGeom_low_nm
        Variable/G v_Wgeom_3D_tailCluster_high_nm = WGeom_high_nm
        Variable/G v_Wgeom_3D_tailCluster_n = WGeom_extreme_n
        Variable/G v_Wgeom_3D_tailCluster_max_nm = WGeom_extreme_max_nm
        Variable/G v_Wgeom_3D_tailCluster_median_nm = WGeom_extreme_median_nm
        Variable/G v_Wgeom_3D_tailCluster_sigma_nm = WGeom_extreme_sigma_nm
        Variable/G v_Wgeom_3D_tailCluster_max_trace = WGeom_extreme_max_trace
        Variable/G v_Wgeom_3D_tailCluster_median_trace = WGeom_extreme_median_trace
        Variable/G v_Wgeom_3D_tailCluster_start = WGeom_extreme_cluster_start
        Variable/G v_Wgeom_3D_tailCluster_size = WGeom_extreme_cluster_size
        Variable/G v_Wgeom_3D_tailCluster_minSize = minClusterSize
    endif

    // -------------------------------------------------------------------------
    // 11. Optional cityscape padding for plot-only waves.
    // -------------------------------------------------------------------------
    String plotLName = ""
    String plotWEffName = ""
    String plotWGeomName = ""
    String plotTName = ""

    if (padCityscape)

        Variable dxL = DimDelta(tmp_Hist_L_work, 0)
        Variable x0L = DimOffset(tmp_Hist_L_work, 0)

        Variable dxWEff = DimDelta(tmp_Hist_W_eff_work, 0)
        Variable x0WEff = DimOffset(tmp_Hist_W_eff_work, 0)

        Make/O/D/N=(nBinsL+2) tmp_Hist_L_plot_work
        Make/O/D/N=(nBinsWEff+2) tmp_Hist_W_eff_plot_work

        tmp_Hist_L_plot_work[0] = 0
        tmp_Hist_L_plot_work[1, nBinsL] = tmp_Hist_L_work[p-1]
        tmp_Hist_L_plot_work[nBinsL+1] = 0

        tmp_Hist_W_eff_plot_work[0] = 0
        tmp_Hist_W_eff_plot_work[1, nBinsWEff] = tmp_Hist_W_eff_work[p-1]
        tmp_Hist_W_eff_plot_work[nBinsWEff+1] = 0

        SetScale/P x, x0L-dxL, dxL, WaveUnits(tmp_Hist_L_work, 0), tmp_Hist_L_plot_work
        SetScale/P x, x0WEff-dxWEff, dxWEff, WaveUnits(tmp_Hist_W_eff_work, 0), tmp_Hist_W_eff_plot_work

        if (use3D)
            if (hasWeights)
                Duplicate/O tmp_Hist_L_plot_work, Hist_L_N_3D_weighted_plot_nm
                Duplicate/O tmp_Hist_W_eff_plot_work, Hist_W_eff_3D_weighted_plot_nm
                plotLName = "Hist_L_N_3D_weighted_plot_nm"
                plotWEffName = "Hist_W_eff_3D_weighted_plot_nm"
            else
                Duplicate/O tmp_Hist_L_plot_work, Hist_L_N_3D_plot_nm
                Duplicate/O tmp_Hist_W_eff_plot_work, Hist_W_eff_3D_plot_nm
                plotLName = "Hist_L_N_3D_plot_nm"
                plotWEffName = "Hist_W_eff_3D_plot_nm"
            endif
        else
            if (hasWeights)
                Duplicate/O tmp_Hist_L_plot_work, Hist_L_N_weighted_plot_nm
                Duplicate/O tmp_Hist_W_eff_plot_work, Hist_W_eff_weighted_plot_nm
                plotLName = "Hist_L_N_weighted_plot_nm"
                plotWEffName = "Hist_W_eff_weighted_plot_nm"
            else
                Duplicate/O tmp_Hist_L_plot_work, Hist_L_N_plot_nm
                Duplicate/O tmp_Hist_W_eff_plot_work, Hist_W_eff_plot_nm
                plotLName = "Hist_L_N_plot_nm"
                plotWEffName = "Hist_W_eff_plot_nm"
            endif
        endif

        if (hasGeomW)
            Variable dxWGeom = DimDelta(tmp_Hist_W_geom_work, 0)
            Variable x0WGeom = DimOffset(tmp_Hist_W_geom_work, 0)

            Make/O/D/N=(nBinsWGeom+2) tmp_Hist_W_geom_plot_work
            tmp_Hist_W_geom_plot_work[0] = 0
            tmp_Hist_W_geom_plot_work[1, nBinsWGeom] = tmp_Hist_W_geom_work[p-1]
            tmp_Hist_W_geom_plot_work[nBinsWGeom+1] = 0

            SetScale/P x, x0WGeom-dxWGeom, dxWGeom, WaveUnits(tmp_Hist_W_geom_work, 0), tmp_Hist_W_geom_plot_work

            if (hasWeights)
                Duplicate/O tmp_Hist_W_geom_plot_work, Hist_W_geom_3D_weighted_plot_nm
                plotWGeomName = "Hist_W_geom_3D_weighted_plot_nm"
            else
                Duplicate/O tmp_Hist_W_geom_plot_work, Hist_W_geom_3D_plot_nm
                plotWGeomName = "Hist_W_geom_3D_plot_nm"
            endif
        endif

        if (hasT)
            Variable dxT = DimDelta(tmp_Hist_T_eff_work, 0)
            Variable x0T = DimOffset(tmp_Hist_T_eff_work, 0)

            Make/O/D/N=(nBinsT+2) tmp_Hist_T_eff_plot_work
            tmp_Hist_T_eff_plot_work[0] = 0
            tmp_Hist_T_eff_plot_work[1, nBinsT] = tmp_Hist_T_eff_work[p-1]
            tmp_Hist_T_eff_plot_work[nBinsT+1] = 0

            SetScale/P x, x0T-dxT, dxT, WaveUnits(tmp_Hist_T_eff_work, 0), tmp_Hist_T_eff_plot_work

            if (use3D)
                if (hasWeights)
                    Duplicate/O tmp_Hist_T_eff_plot_work, Hist_T_eff_3D_weighted_plot
                    plotTName = "Hist_T_eff_3D_weighted_plot"
                else
                    Duplicate/O tmp_Hist_T_eff_plot_work, Hist_T_eff_3D_plot
                    plotTName = "Hist_T_eff_3D_plot"
                endif
            else
                if (hasWeights)
                    Duplicate/O tmp_Hist_T_eff_plot_work, Hist_T_eff_weighted_plot
                    plotTName = "Hist_T_eff_weighted_plot"
                else
                    Duplicate/O tmp_Hist_T_eff_plot_work, Hist_T_eff_plot
                    plotTName = "Hist_T_eff_plot"
                endif
            endif
        endif

        if (use3D)
            String/G s_RayTraceHist_3D_PlotL = plotLName
            String/G s_RayTraceHist_3D_PlotW = plotWEffName
            String/G s_RayTraceHist_3D_PlotW_eff = plotWEffName
            String/G s_RayTraceHist_3D_PlotW_geom = plotWGeomName
            String/G s_RayTraceHist_3D_PlotT = plotTName
            String/G s_RayTraceHist_3D_PlotT_eff = plotTName
        else
            String/G s_RayTraceHist_PlotL = plotLName
            String/G s_RayTraceHist_PlotW = plotWEffName
            String/G s_RayTraceHist_PlotW_eff = plotWEffName
            String/G s_RayTraceHist_PlotW_geom = plotWGeomName
            String/G s_RayTraceHist_PlotT = plotTName
            String/G s_RayTraceHist_PlotT_eff = plotTName
        endif
    endif

    // -------------------------------------------------------------------------
    // 12. Optional display. Uses only waves in RayTraceHist / RayTraceHist_3D.
    // -------------------------------------------------------------------------
    if (doDisplay)

        String weightLabel = ""
        if (hasWeights)
            weightLabel = "weighted "
        endif

        String ensembleLabel
        if (use3D)
            ensembleLabel = "3D bulk trajectories\r"
        else
            ensembleLabel = "2D surface trajectories\r"
        endif

        // ---------------- ℓ_n histogram ----------------
        String histLTrace
        if (padCityscape)
            histLTrace = plotLName
            if (use3D)
                if (hasWeights)
                    Display/K=1 Hist_L_N_3D_weighted_plot_nm
                else
                    Display/K=1 Hist_L_N_3D_plot_nm
                endif
            else
                if (hasWeights)
                    Display/K=1 Hist_L_N_weighted_plot_nm
                else
                    Display/K=1 Hist_L_N_plot_nm
                endif
            endif
        else
            histLTrace = histLName
            if (use3D)
                if (hasWeights)
                    Display/K=1 Hist_L_N_3D_weighted_nm
                else
                    Display/K=1 Hist_L_N_3D_nm
                endif
            else
                if (hasWeights)
                    Display/K=1 Hist_L_N_weighted_nm
                else
                    Display/K=1 Hist_L_N_nm
                endif
            endif
        endif

        ModifyGraph mode=6
        ModifyGraph rgb=(0,0,0)
        ModifyGraph tick=2
        ModifyGraph mirror=2
        ModifyGraph standoff=0

        Label bottom "ℓ\\Bn\\M (nm)"

        if (normalize)
            Label left "Probability (%)"
        else
            if (hasWeights)
                Label left "Channel weight"
            else
                Label left "Channel count"
            endif
        endif

        if (showFitInfo)

            Wave Hist_L_display = $histLTrace
            WaveStats/Q Hist_L_display
            Variable yTopL = 1.08 * V_max
            if (yTopL <= 0 || numtype(yTopL) != 0)
                yTopL = 1
            endif

            if (numtype(LN_extreme_max_nm) == 0)
                SetDrawLayer UserFront
                SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(0,0,0), linethick=2, dash=3
                DrawLine LN_extreme_max_nm, 0, LN_extreme_max_nm, yTopL
            endif

            if (plotMedian && numtype(LN_extreme_median_nm) == 0)
                SetDrawLayer UserFront
                SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(0,0,0), linethick=1, dash=1
                DrawLine LN_extreme_median_nm, 0, LN_extreme_median_nm, yTopL
            endif

            String fitTxtL
            fitTxtL = "Representative large-ℓ\\Bn\\M cluster\r"
            fitTxtL += "ℓ\\Bn\\M: [0.9 max : max]\r"
            fitTxtL += "  global max = " + num2str(round(10*LN_global_max_nm)/10) + " nm\r"
            fitTxtL += "  cluster max = " + num2str(round(10*LN_extreme_max_nm)/10) + " nm\r"
            fitTxtL += "  cluster median = " + num2str(round(10*LN_extreme_median_nm)/10) + " nm\r"
            fitTxtL += "  σ = " + num2str(round(10*LN_extreme_sigma_nm)/10) + " nm\r"
            fitTxtL += "  N = " + num2str(LN_extreme_n) + " / " + num2str(nRay) + "\r"
            fitTxtL += "  N\\Bmin\\M = " + num2str(minClusterSize) + "\r"
            fitTxtL += "  max trace = " + num2str(LN_extreme_max_trace) + "\r"
            fitTxtL += "  median trace = " + num2str(LN_extreme_median_trace)

            TextBox/C/N=RepScaleInfo/F=0/A=LT/X=2/Y=2 fitTxtL
        endif

        Legend/C/N=text0/J/F=0/A=RT ensembleLabel + \
            "\\s(" + histLTrace + ") " + weightLabel + "ℓ\\Bn\\M"

        ModifyGraph margin(left)=45
        ModifyGraph margin(bottom)=35
        ModifyGraph margin(right)=10
        ModifyGraph margin(top)=10
        ModifyGraph width=300
        ModifyGraph height=230

        // ---------------- w_n / w_eff,n histogram ----------------
        String histWEffTrace
        if (padCityscape)
            histWEffTrace = plotWEffName
            if (use3D)
                if (hasWeights)
                    Display/K=1 Hist_W_eff_3D_weighted_plot_nm
                else
                    Display/K=1 Hist_W_eff_3D_plot_nm
                endif
            else
                if (hasWeights)
                    Display/K=1 Hist_W_eff_weighted_plot_nm
                else
                    Display/K=1 Hist_W_eff_plot_nm
                endif
            endif
        else
            histWEffTrace = histWEffName
            if (use3D)
                if (hasWeights)
                    Display/K=1 Hist_W_eff_3D_weighted_nm
                else
                    Display/K=1 Hist_W_eff_3D_nm
                endif
            else
                if (hasWeights)
                    Display/K=1 Hist_W_eff_weighted_nm
                else
                    Display/K=1 Hist_W_eff_nm
                endif
            endif
        endif

        ModifyGraph mode=6
        ModifyGraph rgb=(1,16019,65535)
        ModifyGraph tick=2
        ModifyGraph mirror=2
        ModifyGraph standoff=0

        if (use3D)
            Label bottom "w\\Beff,n\\M (nm)"
        else
            Label bottom "w\\Bn\\M (nm)"
        endif

        if (normalize)
            Label left "Probability (%)"
        else
            if (hasWeights)
                Label left "Channel weight"
            else
                Label left "Channel count"
            endif
        endif

        if (showFitInfo)

            Wave Hist_WEff_display = $histWEffTrace
            WaveStats/Q Hist_WEff_display
            Variable yTopWEff = 1.08 * V_max
            if (yTopWEff <= 0 || numtype(yTopWEff) != 0)
                yTopWEff = 1
            endif

            if (numtype(WEff_extreme_max_nm) == 0)
                SetDrawLayer UserFront
                SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(1,16019,65535), linethick=2, dash=3
                DrawLine WEff_extreme_max_nm, 0, WEff_extreme_max_nm, yTopWEff
            endif

            if (plotMedian && numtype(WEff_extreme_median_nm) == 0)
                SetDrawLayer UserFront
                SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(1,16019,65535), linethick=1, dash=1
                DrawLine WEff_extreme_median_nm, 0, WEff_extreme_median_nm, yTopWEff
            endif

            String fitTxtWEff
            if (use3D)
                fitTxtWEff = "Representative large-w\\Beff,n\\M cluster\r"
                fitTxtWEff += "w\\Beff,n\\M: [0.9 max : max]\r"
            else
                fitTxtWEff = "Representative large-w\\Bn\\M cluster\r"
                fitTxtWEff += "w\\Bn\\M: [0.9 max : max]\r"
            endif
            fitTxtWEff += "  global max = " + num2str(round(10*WEff_global_max_nm)/10) + " nm\r"
            fitTxtWEff += "  cluster max = " + num2str(round(10*WEff_extreme_max_nm)/10) + " nm\r"
            fitTxtWEff += "  cluster median = " + num2str(round(10*WEff_extreme_median_nm)/10) + " nm\r"
            fitTxtWEff += "  σ = " + num2str(round(10*WEff_extreme_sigma_nm)/10) + " nm\r"
            fitTxtWEff += "  N = " + num2str(WEff_extreme_n) + " / " + num2str(nRay) + "\r"
            fitTxtWEff += "  N\\Bmin\\M = " + num2str(minClusterSize) + "\r"
            fitTxtWEff += "  max trace = " + num2str(WEff_extreme_max_trace) + "\r"
            fitTxtWEff += "  median trace = " + num2str(WEff_extreme_median_trace)

            TextBox/C/N=RepScaleInfo/F=0/A=LT/X=2/Y=2 fitTxtWEff
        endif

        if (use3D)
            Legend/C/N=text0/J/F=0/A=RT ensembleLabel + \
                "\\s(" + histWEffTrace + ") " + weightLabel + "w\\Beff,n\\M"
        else
            Legend/C/N=text0/J/F=0/A=RT ensembleLabel + \
                "\\s(" + histWEffTrace + ") " + weightLabel + "w\\Bn\\M"
        endif

        ModifyGraph margin(left)=45
        ModifyGraph margin(bottom)=35
        ModifyGraph margin(right)=10
        ModifyGraph margin(top)=10
        ModifyGraph width=300
        ModifyGraph height=230

        // ---------------- w_geom,n histogram, only for 3D with W_geom ----------------
        if (hasGeomW)

            String histWGeomTrace
            if (padCityscape)
                histWGeomTrace = plotWGeomName
                if (hasWeights)
                    Display/K=1 Hist_W_geom_3D_weighted_plot_nm
                else
                    Display/K=1 Hist_W_geom_3D_plot_nm
                endif
            else
                histWGeomTrace = histWGeomName
                if (hasWeights)
                    Display/K=1 Hist_W_geom_3D_weighted_nm
                else
                    Display/K=1 Hist_W_geom_3D_nm
                endif
            endif

            ModifyGraph mode=6
            ModifyGraph rgb=(0,0,65535)
            ModifyGraph tick=2
            ModifyGraph mirror=2
            ModifyGraph standoff=0

            Label bottom "w\\Bgeom,n\\M (nm)"

            if (normalize)
                Label left "Probability (%)"
            else
                if (hasWeights)
                    Label left "Channel weight"
                else
                    Label left "Channel count"
                endif
            endif

            if (showFitInfo)

                Wave Hist_WGeom_display = $histWGeomTrace
                WaveStats/Q Hist_WGeom_display
                Variable yTopWGeom = 1.08 * V_max
                if (yTopWGeom <= 0 || numtype(yTopWGeom) != 0)
                    yTopWGeom = 1
                endif

                if (numtype(WGeom_extreme_max_nm) == 0)
                    SetDrawLayer UserFront
                    SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(0,0,65535), linethick=2, dash=3
                    DrawLine WGeom_extreme_max_nm, 0, WGeom_extreme_max_nm, yTopWGeom
                endif

                if (plotMedian && numtype(WGeom_extreme_median_nm) == 0)
                    SetDrawLayer UserFront
                    SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(0,0,65535), linethick=1, dash=1
                    DrawLine WGeom_extreme_median_nm, 0, WGeom_extreme_median_nm, yTopWGeom
                endif

                String fitTxtWGeom
                fitTxtWGeom = "Representative large-w\\Bgeom,n\\M cluster\r"
                fitTxtWGeom += "w\\Bgeom,n\\M: [0.9 max : max]\r"
                fitTxtWGeom += "  global max = " + num2str(round(10*WGeom_global_max_nm)/10) + " nm\r"
                fitTxtWGeom += "  cluster max = " + num2str(round(10*WGeom_extreme_max_nm)/10) + " nm\r"
                fitTxtWGeom += "  cluster median = " + num2str(round(10*WGeom_extreme_median_nm)/10) + " nm\r"
                fitTxtWGeom += "  σ = " + num2str(round(10*WGeom_extreme_sigma_nm)/10) + " nm\r"
                fitTxtWGeom += "  N = " + num2str(WGeom_extreme_n) + " / " + num2str(nRay) + "\r"
                fitTxtWGeom += "  N\\Bmin\\M = " + num2str(minClusterSize) + "\r"
                fitTxtWGeom += "  max trace = " + num2str(WGeom_extreme_max_trace) + "\r"
                fitTxtWGeom += "  median trace = " + num2str(WGeom_extreme_median_trace)

                TextBox/C/N=RepScaleInfo/F=0/A=LT/X=2/Y=2 fitTxtWGeom
            endif

            Legend/C/N=text0/J/F=0/A=RT ensembleLabel + \
                "\\s(" + histWGeomTrace + ") " + weightLabel + "w\\Bgeom,n\\M"

            ModifyGraph margin(left)=45
            ModifyGraph margin(bottom)=35
            ModifyGraph margin(right)=10
            ModifyGraph margin(top)=10
            ModifyGraph width=300
            ModifyGraph height=230
        endif

        // ---------------- T_eff,n histogram ----------------
        if (hasT)
            String histTTrace
            if (padCityscape)
                histTTrace = plotTName
                if (use3D)
                    if (hasWeights)
                        Display/K=1 Hist_T_eff_3D_weighted_plot
                    else
                        Display/K=1 Hist_T_eff_3D_plot
                    endif
                else
                    if (hasWeights)
                        Display/K=1 Hist_T_eff_weighted_plot
                    else
                        Display/K=1 Hist_T_eff_plot
                    endif
                endif
            else
                histTTrace = histTName
                if (use3D)
                    if (hasWeights)
                        Display/K=1 Hist_T_eff_3D_weighted
                    else
                        Display/K=1 Hist_T_eff_3D
                    endif
                else
                    if (hasWeights)
                        Display/K=1 Hist_T_eff_weighted
                    else
                        Display/K=1 Hist_T_eff
                    endif
                endif
            endif

            ModifyGraph mode=6
            ModifyGraph rgb=(54741,24158,0)
            ModifyGraph tick=2
            ModifyGraph mirror=2
            ModifyGraph standoff=0

            Label bottom "T\\Beff,n\\M"

            if (normalize)
                Label left "Probability (%)"
            else
                if (hasWeights)
                    Label left "Channel weight"
                else
                    Label left "Channel count"
                endif
            endif

            Legend/C/N=text0/J/F=0/A=RT ensembleLabel + \
                "\\s(" + histTTrace + ") " + weightLabel + "T\\Beff,n\\M"

            ModifyGraph margin(left)=45
            ModifyGraph margin(bottom)=35
            ModifyGraph margin(right)=10
            ModifyGraph margin(top)=10
            ModifyGraph width=300
            ModifyGraph height=230
        endif
    endif

    // -------------------------------------------------------------------------
    // 13. Remove internal work waves from the output folder.
    // -------------------------------------------------------------------------
    KillWaves/Z tmp_L_N_work, tmp_W_eff_work, tmp_W_geom_work, tmp_wChan_work, tmp_T_eff_work
    KillWaves/Z tmp_L_N_List_nm_work, tmp_W_eff_List_nm_work, tmp_W_geom_List_nm_work
    KillWaves/Z tmp_Hist_L_work, tmp_Hist_W_eff_work, tmp_Hist_W_geom_work, tmp_Hist_T_eff_work
    KillWaves/Z tmp_Hist_L_plot_work, tmp_Hist_W_eff_plot_work, tmp_Hist_W_geom_plot_work, tmp_Hist_T_eff_plot_work
    KillWaves/Z tmp_L_cluster_vals, tmp_WEff_cluster_vals, tmp_WGeom_cluster_vals
    KillWaves/Z tmp_L_N_extreme_work_nm, tmp_W_eff_extreme_work_nm, tmp_W_geom_extreme_work_nm
    KillWaves/Z tmp_L_N_extreme_sorted, tmp_W_eff_extreme_sorted, tmp_W_geom_extreme_sorted
    KillWaves/Z tmp_L_minbin_vals, tmp_L_N_minbin_work_nm, tmp_L_N_minbin_sorted

    Variable retTrace = LN_extreme_max_trace

    if (use3D)
        SNS_Log("SNS_RayTrace_Hist_LN_Weff: 3D outputs stored in " + dfPath + ":RayTraceHist_3D", level="INFO")
    else
        SNS_Log("SNS_RayTrace_Hist_LN_Weff: outputs stored in " + dfPath + ":RayTraceHist", level="INFO")
    endif

    SetDataFolder $oldDF

    return retTrace
End

//==============================================================================
// SNS_ExperimentalResolutionHWHM_meV
//
// Returns the half-width at half maximum [meV] of the combined
// thermal + lock-in modulation resolution kernel used by
// SNS_ApplyDOS_Broadening_TplusMod.
//
// This is the Gamma_exp to use in the w-bin criterion.
//==============================================================================
Function SNS_ExperimentalResolutionHWHM_meV([dE_meV, Emax_meV])
	Variable dE_meV, Emax_meV

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	if (ParamIsDefault(dE_meV))
		dE_meV = 0.001       // 1 µeV grid
	endif
	if (ParamIsDefault(Emax_meV))
		Emax_meV = 2.0
	endif

	Variable NE = 2*ceil(Emax_meV/dE_meV) + 1
	Variable dE_eV = dE_meV * 1e-3
	Variable E0_eV = -0.5*(NE-1)*dE_eV

	Make/O/D/N=(NE) root:SNS_Settings:tmp_res_E_axis_eV
	Wave E_axis = root:SNS_Settings:tmp_res_E_axis_eV
	E_axis = E0_eV + p*dE_eV
	SetScale/P x, E0_eV, dE_eV, "eV", E_axis

	SNS_MakeThermalKernel(E_axis, params.T_K, "root:SNS_Settings:tmp_res_KT")
	SNS_MakeModulationKernel(E_axis, params.V_mod, "root:SNS_Settings:tmp_res_KM")

	Wave KT = root:SNS_Settings:tmp_res_KT
	Wave KM = root:SNS_Settings:tmp_res_KM

	Duplicate/O KT, root:SNS_Settings:tmp_res_K
	Wave K = root:SNS_Settings:tmp_res_K

	// Combined experimental-resolution kernel:
	// same order as SNS_ApplyDOS_Broadening_TplusMod.
	Convolve/A KM, K

	WaveStats/Q K
	Variable halfMax = 0.5 * V_max

	Variable i0 = round((NE-1)/2)
	Variable i

	for (i = i0; i < NE-1; i += 1)
		if (K[i] >= halfMax && K[i+1] <= halfMax)
			Variable x1 = E_axis[i]
			Variable x2 = E_axis[i+1]
			Variable y1 = K[i]
			Variable y2 = K[i+1]

			Variable xHalf = x1 + (halfMax-y1)*(x2-x1)/(y2-y1)

			KillWaves/Z root:SNS_Settings:tmp_res_E_axis_eV
			KillWaves/Z root:SNS_Settings:tmp_res_KT
			KillWaves/Z root:SNS_Settings:tmp_res_KM
			KillWaves/Z root:SNS_Settings:tmp_res_K

			return abs(xHalf) * 1e3      // eV -> meV
		endif
	endfor

	KillWaves/Z root:SNS_Settings:tmp_res_E_axis_eV
	KillWaves/Z root:SNS_Settings:tmp_res_KT
	KillWaves/Z root:SNS_Settings:tmp_res_KM
	KillWaves/Z root:SNS_Settings:tmp_res_K

	return NaN
End

Function SNS_HbarVF_meVnm_FromParams()
	STRUCT SNS_Params params
	SNS_LoadParams(params)

	// hbar*vF [eV*s * m/s] = eV*m.
	// eV*m -> meV*nm gives factor 1e12.
	return HBAR_eVs * params.vF * 1e12
End


Function SNS_PhaseSlope_meVPerRad_FromParams(ell_nm)
	Variable ell_nm

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	Variable Delta_meV = 1e3 * params.Delta
	Variable hbarVF_meVnm = SNS_HbarVF_meVnm_FromParams()

	return 1 / (2/Delta_meV + 2*ell_nm/hbarVF_meVnm)
End


Function SNS_WBinWidthFromResolution_nm(w_nm, ell_nm, [GammaExp_meV])
	Variable w_nm, ell_nm
	Variable GammaExp_meV

	if (ParamIsDefault(GammaExp_meV))
		GammaExp_meV = SNS_ExperimentalResolutionHWHM_meV()
	endif

	Variable Sphi = SNS_PhaseSlope_meVPerRad_FromParams(ell_nm)

	if (numtype(GammaExp_meV) != 0 || GammaExp_meV <= 0)
		return NaN
	endif
	if (numtype(Sphi) != 0 || Sphi <= 0)
		return NaN
	endif

	return GammaExp_meV * w_nm / (pi * Sphi)
End

//==============================================================================
// SNS_MakeWResolutionDensity
//
// Computes unbinned experimental-resolution-broadened w-density:
//
//   rho(wc) = sum_i weight_i * K_exp[ Sphi_i*pi*(w_i/wc - 1) ]
//
// Input waves must be in nm.
// Typically use:
//   W_eff_List_nm
//   L_N_List_nm
//   optional wChan
//
// Output:
//   outName, x-scaled in nm, normalized to max=1 by default.
//
//==============================================================================



Function SNS_MakeWResolutionDensity(W_nm, outName, [L_nm, weightWave, wMin_nm, wMax_nm, dw_nm, Emax_meV, dE_meV, absW, normalize, doDisplay])
	Wave W_nm
	String outName
	Wave L_nm
	Wave weightWave
	Variable wMin_nm, wMax_nm, dw_nm, Emax_meV, dE_meV
	Variable absW, normalize, doDisplay

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	if (ParamIsDefault(absW))
		absW = 1
	endif
	if (ParamIsDefault(normalize))
		normalize = 1
	endif
	if (ParamIsDefault(doDisplay))
		doDisplay = 1
	endif
	if (ParamIsDefault(dw_nm))
		dw_nm = 1
	endif
	if (ParamIsDefault(Emax_meV))
		Emax_meV = 3
	endif
	if (ParamIsDefault(dE_meV))
		dE_meV = 0.001
	endif

	Variable n = numpnts(W_nm)
	if (n <= 0)
		Abort "SNS_MakeWResolutionDensity: empty W_nm wave."
	endif

	Variable hasL = !ParamIsDefault(L_nm) && WaveExists(L_nm) && (numpnts(L_nm) == n)
	Variable hasWeights = !ParamIsDefault(weightWave) && WaveExists(weightWave) && (numpnts(weightWave) == n)

	// Determine wc grid.
	WaveStats/Q W_nm
	Variable wDataMin = V_min
	Variable wDataMax = V_max

	if (absW)
		// If W_nm may contain signed values, find abs max manually.
		Variable i
		wDataMin = Inf
		wDataMax = -Inf
		for (i = 0; i < n; i += 1)
			if (numtype(W_nm[i]) == 0)
				wDataMin = min(wDataMin, abs(W_nm[i]))
				wDataMax = max(wDataMax, abs(W_nm[i]))
			endif
		endfor
	endif

	if (ParamIsDefault(wMin_nm))
		wMin_nm = max(1, floor(wDataMin))
	endif
	if (ParamIsDefault(wMax_nm))
		wMax_nm = ceil(wDataMax)
	endif

	if (dw_nm <= 0 || wMax_nm <= wMin_nm)
		Abort "SNS_MakeWResolutionDensity: invalid w grid."
	endif

	Variable nW = floor((wMax_nm - wMin_nm)/dw_nm) + 1

	// -------------------------------------------------------------------------
	// Build experimental resolution kernel K_exp(E), using the same ingredients
	// as SNS_ApplyDOS_Broadening_TplusMod: thermal + modulation convolution.
	// Energy axis is in eV because the existing kernel helpers use the DOS axis.
	// -------------------------------------------------------------------------
	Variable nE = 2*ceil(Emax_meV/dE_meV) + 1
	Variable dE_eV = dE_meV * 1e-3
	Variable E0_eV = -0.5*(nE-1)*dE_eV

	Make/O/D/N=(nE) root:SNS_Settings:tmp_Wdens_E_axis_eV
	Wave E_axis = root:SNS_Settings:tmp_Wdens_E_axis_eV
	E_axis = E0_eV + p*dE_eV
	SetScale/P x, E0_eV, dE_eV, "eV", E_axis

	SNS_MakeThermalKernel(E_axis, params.T_K, "root:SNS_Settings:tmp_Wdens_KT")
	SNS_MakeModulationKernel(E_axis, params.V_mod, "root:SNS_Settings:tmp_Wdens_KM")

	Wave KT = root:SNS_Settings:tmp_Wdens_KT
	Wave KM = root:SNS_Settings:tmp_Wdens_KM

	Duplicate/O KT, root:SNS_Settings:tmp_Wdens_Kexp
	Wave Kexp = root:SNS_Settings:tmp_Wdens_Kexp

	Convolve/A KM, Kexp

	WaveStats/Q Kexp
	if (V_max > 0)
		Kexp /= V_max
	endif

	// -------------------------------------------------------------------------
	// Compute rho_exp(wc)
	// -------------------------------------------------------------------------
	Make/O/D/N=(nW) $outName
	Wave rho = $outName
	SetScale/P x, wMin_nm, dw_nm, "nm", rho
	rho = 0

	Variable j, wc, wi, elli, Sphi, Ei_meV, Ei_eV, wt

	for (j = 0; j < nW; j += 1)
		wc = wMin_nm + j*dw_nm

		if (wc <= 0)
			rho[j] = NaN
			continue
		endif

		for (i = 0; i < n; i += 1)
			if (numtype(W_nm[i]) != 0)
				continue
			endif

			wi = W_nm[i]
			if (absW)
				wi = abs(wi)
			endif

			if (wi <= 0)
				continue
			endif

			if (hasL)
				elli = L_nm[i]
			else
				// Fallback: use wc as crude length scale only if no L wave supplied.
				// Better: always pass L_N_List_nm.
				elli = wc
			endif

			if (numtype(elli) != 0 || elli <= 0)
				continue
			endif

			Sphi = SNS_PhaseSlope_meVPerRad_FromParams(elli)
			if (numtype(Sphi) != 0 || Sphi <= 0)
				continue
			endif

			if (hasWeights)
				wt = weightWave[i]
				if (numtype(wt) != 0 || wt < 0)
					wt = 0
				endif
			else
				wt = 1
			endif

			// Energy mismatch at field q = pi/wc:
			// E_i = Sphi*pi*(wi/wc - 1)
			Ei_meV = Sphi*pi*(wi/wc - 1)
			Ei_eV = Ei_meV * 1e-3

			if (abs(Ei_meV) <= Emax_meV)
				rho[j] += wt * interp(Ei_eV, E_axis, Kexp)
			endif
		endfor
	endfor

	if (normalize)
		WaveStats/Q rho
		if (V_max > 0)
			rho /= V_max
		endif
	endif

	// Store useful metadata next to output wave.
	Variable/G $(outName + "_dw_nm") = dw_nm
	Variable/G $(outName + "_wMin_nm") = wMin_nm
	Variable/G $(outName + "_wMax_nm") = wMax_nm
	Variable/G $(outName + "_T_K") = params.T_K
	Variable/G $(outName + "_Vmod_V") = params.V_mod
	Variable/G $(outName + "_Delta_eff_meV") = 1e3*params.Delta

	KillWaves/Z root:SNS_Settings:tmp_Wdens_E_axis_eV
	KillWaves/Z root:SNS_Settings:tmp_Wdens_KT
	KillWaves/Z root:SNS_Settings:tmp_Wdens_KM
	KillWaves/Z root:SNS_Settings:tmp_Wdens_Kexp

	if (doDisplay)
		Display/K=1 rho
		ModifyGraph mode=0
		ModifyGraph tick=2, mirror=2, standoff=0
		Label bottom "candidate w\\Bc\\M (nm)"
		if (normalize)
			Label left "resolution-broadened zero-bias weight (norm.)"
		else
			Label left "resolution-broadened zero-bias weight"
		endif
	endif

	return 0
End

//==============================================================================
// SNS_QperNm_FromB_T
//
// q(B) in nm^-1 for the magnetic phase q*w.
// q(B) = 2*pi*h_eff*B/Phi0.
//==============================================================================
Function SNS_QperNm_FromB_T(B_T, hEff_nm)
	Variable B_T, hEff_nm

	Variable Phi0_Wb = 2.067833848e-15

	if (numtype(B_T) != 0 || numtype(hEff_nm) != 0 || hEff_nm <= 0)
		return NaN
	endif

	// hEff_nm * 1e-9 gives meters; q is 1/m; convert to 1/nm by *1e-9.
	return 2*pi*hEff_nm*B_T*1e-18/Phi0_Wb
End


//==============================================================================
// SNS_Wnm_FromB_T
//
// Converts field scale to effective transverse projection:
//   B = Phi0/(2*h_eff*w)
//==============================================================================
Function SNS_Wnm_FromB_T(B_T, hEff_nm)
	Variable B_T, hEff_nm

	Variable Phi0_Wb = 2.067833848e-15

	if (numtype(B_T) != 0 || numtype(hEff_nm) != 0 || B_T <= 0 || hEff_nm <= 0)
		return NaN
	endif

	return Phi0_Wb * 1e18 / (2*hEff_nm*B_T)
End


//==============================================================================
// SNS_MakeBResolutionDensity
//
// Computes experimental-resolution-broadened zero-bias field density:
//
//   rho(B) = sum_i a_i K_exp[ Sphi_i * ( q(B)*w_i - pi ) ]
//
// W_nm and L_nm must be matched trajectory waves in nm.
// Output wave x-scaling is B in T.
//==============================================================================
Function SNS_MakeBResolutionDensity(W_nm, L_nm, outName, hEff_nm, [weightWave, BMin_T, BMax_T, dB_T, Emax_meV, dE_meV, absW, normalize, doDisplay])
	Wave W_nm, L_nm
	String outName
	Variable hEff_nm
	Wave weightWave
	Variable BMin_T, BMax_T, dB_T, Emax_meV, dE_meV
	Variable absW, normalize, doDisplay

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	if (ParamIsDefault(absW))
		absW = 1
	endif
	if (ParamIsDefault(normalize))
		normalize = 1
	endif
	if (ParamIsDefault(doDisplay))
		doDisplay = 1
	endif
	if (ParamIsDefault(dB_T))
		dB_T = 0.001
	endif
	if (ParamIsDefault(Emax_meV))
		Emax_meV = 3
	endif
	if (ParamIsDefault(dE_meV))
		dE_meV = 0.001
	endif
	if (ParamIsDefault(BMin_T))
		BMin_T = 0
	endif

	Variable n = numpnts(W_nm)
	if (n <= 0 || numpnts(L_nm) != n)
		Abort "SNS_MakeBResolutionDensity: W_nm and L_nm must have same nonzero length."
	endif
	if (hEff_nm <= 0 || numtype(hEff_nm) != 0)
		Abort "SNS_MakeBResolutionDensity: invalid hEff_nm."
	endif

	Variable hasWeights = !ParamIsDefault(weightWave) && WaveExists(weightWave) && (numpnts(weightWave) == n)

	// Default BMax from smallest positive |w| first-zero condition.
	if (ParamIsDefault(BMax_T))
		Variable wMinPos = Inf
		Variable i, wi
		for (i = 0; i < n; i += 1)
			wi = W_nm[i]
			if (absW)
				wi = abs(wi)
			endif
			if (numtype(wi) == 0 && wi > 0)
				wMinPos = min(wMinPos, wi)
			endif
		endfor

		if (wMinPos < Inf)
			BMax_T = 2
		else
			Abort "SNS_MakeBResolutionDensity: could not infer BMax_T."
		endif
	endif

	if (dB_T <= 0 || BMax_T <= BMin_T)
		Abort "SNS_MakeBResolutionDensity: invalid B grid."
	endif

	Variable nB = floor((BMax_T - BMin_T)/dB_T) + 1

	// -------------------------------------------------------------------------
	// Build experimental resolution kernel K_exp(E)
	// -------------------------------------------------------------------------
	Variable nE = 2*ceil(Emax_meV/dE_meV) + 1
	Variable dE_eV = dE_meV * 1e-3
	Variable E0_eV = -0.5*(nE-1)*dE_eV

	Make/O/D/N=(nE) root:SNS_Settings:tmp_Bdens_E_axis_eV
	Wave E_axis = root:SNS_Settings:tmp_Bdens_E_axis_eV
	E_axis = E0_eV + p*dE_eV
	SetScale/P x, E0_eV, dE_eV, "eV", E_axis

	SNS_MakeThermalKernel(E_axis, params.T_K, "root:SNS_Settings:tmp_Bdens_KT")
	SNS_MakeModulationKernel(E_axis, params.V_mod, "root:SNS_Settings:tmp_Bdens_KM")

	Wave KT = root:SNS_Settings:tmp_Bdens_KT
	Wave KM = root:SNS_Settings:tmp_Bdens_KM

	Duplicate/O KT, root:SNS_Settings:tmp_Bdens_Kexp
	Wave Kexp = root:SNS_Settings:tmp_Bdens_Kexp

	Convolve/A KM, Kexp

	WaveStats/Q Kexp
	if (V_max > 0)
		Kexp /= V_max
	endif

	// -------------------------------------------------------------------------
	// Compute rho_exp(B)
	// -------------------------------------------------------------------------
	Make/O/D/N=(nB) $outName
	Wave rho = $outName
	SetScale/P x, BMin_T, dB_T, "T", rho
	rho = 0

	Variable j, B_T, q_nm, elli, Sphi, dEi_meV, dEi_eV, wt

	for (j = 0; j < nB; j += 1)
		B_T = BMin_T + j*dB_T
		q_nm = SNS_QperNm_FromB_T(B_T, hEff_nm)

		for (i = 0; i < n; i += 1)
			wi = W_nm[i]
			if (absW)
				wi = abs(wi)
			endif

			elli = L_nm[i]

			if (numtype(wi) != 0 || numtype(elli) != 0 || wi <= 0 || elli <= 0)
				continue
			endif

			Sphi = SNS_PhaseSlope_meVPerRad_FromParams(elli)
			if (numtype(Sphi) != 0 || Sphi <= 0)
				continue
			endif

			if (hasWeights)
				wt = weightWave[i]
				if (numtype(wt) != 0 || wt < 0)
					wt = 0
				endif
			else
				wt = 1
			endif

			dEi_meV = Sphi * (q_nm*wi - pi)
			dEi_eV = dEi_meV * 1e-3

			if (abs(dEi_meV) <= Emax_meV)
				rho[j] += wt * interp(dEi_eV, E_axis, Kexp)
			endif
		endfor
	endfor

	if (normalize)
		WaveStats/Q rho
		if (V_max > 0)
			rho /= V_max
		endif
	endif

	Variable/G $(outName + "_hEff_nm") = hEff_nm
	Variable/G $(outName + "_BMin_T") = BMin_T
	Variable/G $(outName + "_BMax_T") = BMax_T
	Variable/G $(outName + "_dB_T") = dB_T
	Variable/G $(outName + "_T_K") = params.T_K
	Variable/G $(outName + "_Vmod_V") = params.V_mod
	Variable/G $(outName + "_Delta_eff_meV") = 1e3*params.Delta

	KillWaves/Z root:SNS_Settings:tmp_Bdens_E_axis_eV
	KillWaves/Z root:SNS_Settings:tmp_Bdens_KT
	KillWaves/Z root:SNS_Settings:tmp_Bdens_KM
	KillWaves/Z root:SNS_Settings:tmp_Bdens_Kexp

	if (doDisplay)
		Display/K=1 rho
		ModifyGraph mode=0
		ModifyGraph tick=2, mirror=2, standoff=0
		Label bottom "B (T)"
		Label left "resolution-broadened zero-bias weight"
	endif

	return 0
End


//==============================================================================
// SNS_MakeBGapRayContributions
//
// At B_gap, compute each ray's zero-bias field contribution:
//
//   A_i = a_i K_exp[ Sphi_i * ( q(B_gap)*w_i - pi ) ]
//
// Creates:
//   outBase + "_contribution"
//   outBase + "_family_mask"
//   outBase + "_alpha"
//==============================================================================
Function SNS_MakeBGapRayContributions(W_nm, L_nm, BGap_T, hEff_nm, outBase, [weightWave, Emax_meV, dE_meV, absW, alphaMin, alphaMax, alphaGamma])
	Wave W_nm, L_nm
	Variable BGap_T, hEff_nm
	String outBase
	Wave weightWave
	Variable Emax_meV, dE_meV, absW, alphaMin, alphaMax, alphaGamma

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	if (ParamIsDefault(Emax_meV))
		Emax_meV = 3
	endif
	if (ParamIsDefault(dE_meV))
		dE_meV = 0.001
	endif
	if (ParamIsDefault(absW))
		absW = 1
	endif
	if (ParamIsDefault(alphaMin))
		alphaMin = 0.10
	endif
	if (ParamIsDefault(alphaMax))
		alphaMax = 1.00
	endif
	if (ParamIsDefault(alphaGamma))
		alphaGamma = 3.0
	endif

	alphaMin = max(0, min(1, alphaMin))
	alphaMax = max(0, min(1, alphaMax))
	alphaGamma = max(0.1, alphaGamma)

	Variable n = numpnts(W_nm)
	if (n <= 0 || numpnts(L_nm) != n)
		Abort "SNS_MakeBGapRayContributions: W_nm and L_nm must have same nonzero length."
	endif
	if (numtype(BGap_T) != 0 || BGap_T <= 0 || hEff_nm <= 0)
		Abort "SNS_MakeBGapRayContributions: invalid BGap_T or hEff_nm."
	endif

	Variable hasWeights = !ParamIsDefault(weightWave) && WaveExists(weightWave) && (numpnts(weightWave) == n)

	// Build experimental kernel.
	Variable nE = 2*ceil(Emax_meV/dE_meV) + 1
	Variable dE_eV = dE_meV * 1e-3
	Variable E0_eV = -0.5*(nE-1)*dE_eV

	Make/O/D/N=(nE) root:SNS_Settings:tmp_Bgap_E_axis_eV
	Wave E_axis = root:SNS_Settings:tmp_Bgap_E_axis_eV
	E_axis = E0_eV + p*dE_eV
	SetScale/P x, E0_eV, dE_eV, "eV", E_axis

	SNS_MakeThermalKernel(E_axis, params.T_K, "root:SNS_Settings:tmp_Bgap_KT")
	SNS_MakeModulationKernel(E_axis, params.V_mod, "root:SNS_Settings:tmp_Bgap_KM")

	Wave KT = root:SNS_Settings:tmp_Bgap_KT
	Wave KM = root:SNS_Settings:tmp_Bgap_KM

	Duplicate/O KT, root:SNS_Settings:tmp_Bgap_Kexp
	Wave Kexp = root:SNS_Settings:tmp_Bgap_Kexp

	Convolve/A KM, Kexp

	WaveStats/Q Kexp
	if (V_max > 0)
		Kexp /= V_max
	endif

	Make/O/D/N=(n) $(outBase + "_contribution")
	Make/O/D/N=(n) $(outBase + "_family_mask")
	Make/O/D/N=(n) $(outBase + "_alpha")

	Wave contrib = $(outBase + "_contribution")
	Wave mask = $(outBase + "_family_mask")
	Wave alpha = $(outBase + "_alpha")

	contrib = 0
	mask = 0
	alpha = 0

	Variable i, wi, elli, Sphi, q_nm, dEi_meV, dEi_eV, wt
	q_nm = SNS_QperNm_FromB_T(BGap_T, hEff_nm)

	for (i = 0; i < n; i += 1)
		wi = W_nm[i]
		if (absW)
			wi = abs(wi)
		endif
		elli = L_nm[i]

		if (numtype(wi) != 0 || numtype(elli) != 0 || wi <= 0 || elli <= 0)
			continue
		endif

		Sphi = SNS_PhaseSlope_meVPerRad_FromParams(elli)
		if (numtype(Sphi) != 0 || Sphi <= 0)
			continue
		endif

		if (hasWeights)
			wt = weightWave[i]
			if (numtype(wt) != 0 || wt < 0)
				wt = 0
			endif
		else
			wt = 1
		endif

		dEi_meV = Sphi * (q_nm*wi - pi)
		dEi_eV = dEi_meV * 1e-3

		if (abs(dEi_meV) <= Emax_meV)
			contrib[i] = wt * interp(dEi_eV, E_axis, Kexp)
		endif
	endfor

	WaveStats/Q contrib
	Variable Amax = V_max
	Variable Ahwhm = 0.5*Amax
	Variable alphaFrac

	if (Amax > 0 && numtype(Amax) == 0)
		for (i = 0; i < n; i += 1)
			if (contrib[i] >= Ahwhm)
				mask[i] = 1

				if (Amax > Ahwhm)
					alphaFrac = (contrib[i] - Ahwhm)/(Amax - Ahwhm)
					alphaFrac = max(0, min(1, alphaFrac))
					alpha[i] = alphaMin + (alphaMax-alphaMin)*alphaFrac^alphaGamma
				else
					alpha[i] = alphaMax
				endif

				alpha[i] = max(alphaMin, min(alphaMax, alpha[i]))
			endif
		endfor
	endif

	Variable/G $(outBase + "_BGap_T") = BGap_T
	Variable/G $(outBase + "_wFromB_nm") = SNS_Wnm_FromB_T(BGap_T, hEff_nm)
	Variable/G $(outBase + "_hEff_nm") = hEff_nm
	Variable/G $(outBase + "_Amax") = Amax
	Variable/G $(outBase + "_Ahwhm") = Ahwhm
	Variable/G $(outBase + "_alphaMin") = alphaMin
	Variable/G $(outBase + "_alphaMax") = alphaMax
	Variable/G $(outBase + "_alphaGamma") = alphaGamma

	KillWaves/Z root:SNS_Settings:tmp_Bgap_E_axis_eV
	KillWaves/Z root:SNS_Settings:tmp_Bgap_KT
	KillWaves/Z root:SNS_Settings:tmp_Bgap_KM
	KillWaves/Z root:SNS_Settings:tmp_Bgap_Kexp

	return 0
End


//==============================================================================
// SNS_FindFirstBPeakInResolutionDensity
//
// Finds the first finite-field peak of rho(B), searching from low B to high B.
// rho x-scaling must be B in T.
//==============================================================================
Function SNS_FindFirstBPeakInResolutionDensity(rho, [minRelHeight, box, BStart_T, BEnd_T, doDisplay])
	Wave rho
	Variable minRelHeight, box, BStart_T, BEnd_T, doDisplay

	if (ParamIsDefault(minRelHeight))
		minRelHeight = 0.05
	endif
	if (ParamIsDefault(box))
		box = 3
	endif
	if (ParamIsDefault(doDisplay))
		doDisplay = 1
	endif

	Variable n = numpnts(rho)
	if (n < 3)
		Abort "SNS_FindFirstBPeakInResolutionDensity: rho wave too short."
	endif

	Variable x0 = DimOffset(rho, 0)
	Variable dx = DimDelta(rho, 0)
	Variable x1 = x0 + (n - 1) * dx

	if (ParamIsDefault(BStart_T))
		BStart_T = x0 + dx      // skip exact zero-field point
	endif
	if (ParamIsDefault(BEnd_T))
		BEnd_T = x1
	endif

	WaveStats/Q/R=(BStart_T, BEnd_T) rho
	Variable localMax = V_max

	if (localMax <= 0 || numtype(localMax) != 0)
		return NaN
	endif

	Variable minLevel = minRelHeight * localMax

	// Search low field -> high field.
	FindPeak/Q/B=(box)/M=(minLevel)/R=(BStart_T, BEnd_T) rho

	if (V_flag != 0 || numtype(V_PeakLoc) != 0)
		SNS_Log("SNS_FindFirstBPeakInResolutionDensity: no first B peak found.", level="WARN")
		return NaN
	endif

	Variable BPeak_T = V_PeakLoc
	Variable peakHeight = V_PeakVal
	Variable edge1_T = V_LeadingEdgeLoc
	Variable edge2_T = V_TrailingEdgeLoc

	String base = NameOfWave(rho)
	Variable/G $(base + "_BPeak_T") = BPeak_T
	Variable/G $(base + "_peakHeight") = peakHeight
	Variable/G $(base + "_peakLow_T") = min(edge1_T, edge2_T)
	Variable/G $(base + "_peakHigh_T") = max(edge1_T, edge2_T)
	Variable/G $(base + "_peakWidth_T") = V_PeakWidth
	Variable/G $(base + "_minRelHeight") = minRelHeight
	Variable/G $(base + "_FindPeak_box") = box

	if (doDisplay)
		Display/K=1 rho
		ModifyGraph mode=0
		ModifyGraph tick=2, mirror=2, standoff=0
		Label bottom "B (T)"
		Label left "resolution-broadened zero-bias weight"

		SetDrawLayer UserFront
		SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(1,16019,65535), linethick=2, dash=3
		DrawLine BPeak_T, 0, BPeak_T, 1.05 * peakHeight
	endif

	return BPeak_T
End

Function SNS_FindHighWPeakInResolutionDensity(rho, [minRelHeight, box, edgeSkipPts, doDisplay])
	Wave rho
	Variable minRelHeight, box, edgeSkipPts, doDisplay

	if (ParamIsDefault(minRelHeight))
		minRelHeight = 0.05
	endif
	if (ParamIsDefault(box))
		box = 5
	endif
	if (ParamIsDefault(edgeSkipPts))
		edgeSkipPts = max(3, box)
	endif
	if (ParamIsDefault(doDisplay))
		doDisplay = 1
	endif

	Variable n = numpnts(rho)
	if (n < 3 + 2*edgeSkipPts)
		Abort "SNS_FindHighWPeakInResolutionDensity: rho wave too short."
	endif

	Variable x0 = DimOffset(rho, 0)
	Variable dx = DimDelta(rho, 0)

	WaveStats/Q rho
	Variable globalMax = V_max
	if (globalMax <= 0 || numtype(globalMax) != 0)
		return NaN
	endif

	Variable minLevel = minRelHeight * globalMax

	// Avoid boundary locking:
	// search high -> low, but start edgeSkipPts inside the high-w boundary.
	Variable pHigh = n - 1 - edgeSkipPts
	Variable pLow  = edgeSkipPts

	Variable xHigh = x0 + pHigh*dx
	Variable xLow  = x0 + pLow*dx

	FindPeak/Q/B=(box)/M=(minLevel)/R=(xHigh, xLow) rho

	if (V_flag != 0 || numtype(V_PeakLoc) != 0)
		SNS_Log("SNS_FindHighWPeakInResolutionDensity: no high-w peak found.", level="WARN")
		return NaN
	endif

	Variable wPeak_nm = V_PeakLoc
	Variable peakHeight = V_PeakVal

	Variable wEdge1_nm = V_LeadingEdgeLoc
	Variable wEdge2_nm = V_TrailingEdgeLoc
	Variable wPeakLow_nm = min(wEdge1_nm, wEdge2_nm)
	Variable wPeakHigh_nm = max(wEdge1_nm, wEdge2_nm)

	String base = NameOfWave(rho)
	Variable/G $(base + "_wPeak_nm") = wPeak_nm
	Variable/G $(base + "_peakHeight") = peakHeight
	Variable/G $(base + "_peakLow_nm") = wPeakLow_nm
	Variable/G $(base + "_peakHigh_nm") = wPeakHigh_nm
	Variable/G $(base + "_peakWidth_nm") = V_PeakWidth
	Variable/G $(base + "_minRelHeight") = minRelHeight
	Variable/G $(base + "_minLevel") = minLevel
	Variable/G $(base + "_FindPeak_box") = box
	Variable/G $(base + "_FindPeak_edgeSkipPts") = edgeSkipPts

	if (doDisplay)
		Display/K=1 rho
		ModifyGraph mode=0
		ModifyGraph tick=2, mirror=2, standoff=0
		Label bottom "candidate w\\Bc\\M (nm)"
		Label left "resolution-broadened zero-bias weight"

		SetDrawLayer UserFront
		SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(1,16019,65535), linethick=2, dash=3
		DrawLine wPeak_nm, 0, wPeak_nm, 1.05*peakHeight

		if (numtype(wPeakLow_nm) == 0 && numtype(wPeakHigh_nm) == 0)
			SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(0,0,0), linethick=1, dash=1
			DrawLine wPeakLow_nm, 0.5*peakHeight, wPeakHigh_nm, 0.5*peakHeight
		endif

		String txt
		txt = "high-w peak from FindPeak\r"
		txt += "w\\Bgap\\M = " + num2str(round(10*wPeak_nm)/10) + " nm\r"
		txt += "height = " + num2str(round(1000*peakHeight)/1000) + "\r"
		txt += "range = [" + num2str(round(10*wPeakLow_nm)/10) + ", " + num2str(round(10*wPeakHigh_nm)/10) + "] nm"
		TextBox/C/N=PeakInfo/F=0/A=LT/X=2/Y=2 txt
	endif

	return wPeak_nm
End





//==============================================================================
// SNS_ABSZeroFieldEnergy_meV_FromEll
//
// Zero-field ABS energy for one trajectory length ell_nm.
// Uses:
//   2 E ell/(hbar vF) - 2 acos(E/Delta_eff) = 0
//
// Returns E in meV.
//==============================================================================
Function SNS_ABSZeroFieldEnergy_meV_FromEll(ell_nm)
	Variable ell_nm

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	Variable Delta_meV = 1e3 * params.Delta
	Variable hbarVF_meVnm = HBAR_eVs * params.vF * 1e12

	if (numtype(ell_nm) != 0 || ell_nm <= 0)
		return NaN
	endif
	if (Delta_meV <= 0 || hbarVF_meVnm <= 0)
		return NaN
	endif

	Variable lo = 0
	Variable hi = Delta_meV * (1 - 1e-10)
	Variable mid, fmid
	Variable k

	// F(E) = 2Eell/hbarvF - 2acos(E/Delta)
	// F(0)<0, F(Delta)>0 for normal parameter range.
	for (k = 0; k < 80; k += 1)
		mid = 0.5 * (lo + hi)
		fmid = 2 * mid * ell_nm / hbarVF_meVnm - 2 * acos(mid / Delta_meV)

		if (fmid > 0)
			hi = mid
		else
			lo = mid
		endif
	endfor

	return 0.5 * (lo + hi)
End


//==============================================================================
// SNS_MakeLResolutionDensity
//
// Computes unbinned experimental-resolution-broadened ell-density:
//
//   rho(ell_c) = sum_i weight_i * K_exp[ E(ell_i) - E(ell_c) ]
//
// Input:
//   L_nm must be trajectory length wave in nm.
//
// Output:
//   outName, x-scaled in nm, normalized to max=1 by default.
//==============================================================================
Function SNS_MakeLResolutionDensity(L_nm, outName, [weightWave, lMin_nm, lMax_nm, dl_nm, Emax_meV, dE_meV, normalize, doDisplay])
	Wave L_nm
	String outName
	Wave weightWave
	Variable lMin_nm, lMax_nm, dl_nm, Emax_meV, dE_meV
	Variable normalize, doDisplay

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	if (ParamIsDefault(normalize))
		normalize = 1
	endif
	if (ParamIsDefault(doDisplay))
		doDisplay = 1
	endif
	if (ParamIsDefault(dl_nm))
		dl_nm = 1
	endif
	if (ParamIsDefault(Emax_meV))
		Emax_meV = 3
	endif
	if (ParamIsDefault(dE_meV))
		dE_meV = 0.001
	endif

	Variable n = numpnts(L_nm)
	if (n <= 0)
		Abort "SNS_MakeLResolutionDensity: empty L_nm wave."
	endif

	Variable hasWeights = !ParamIsDefault(weightWave) && WaveExists(weightWave) && (numpnts(weightWave) == n)

	WaveStats/Q L_nm
	Variable LDataMin = V_min
	Variable LDataMax = V_max

	if (ParamIsDefault(lMin_nm))
		lMin_nm = max(1, floor(LDataMin))
	endif
	if (ParamIsDefault(lMax_nm))
		lMax_nm = ceil(LDataMax)
	endif

	if (dl_nm <= 0 || lMax_nm <= lMin_nm)
		Abort "SNS_MakeLResolutionDensity: invalid ell grid."
	endif

	Variable nL = floor((lMax_nm - lMin_nm) / dl_nm) + 1

	// -------------------------------------------------------------------------
	// Build experimental resolution kernel K_exp(E)
	// -------------------------------------------------------------------------
	Variable nE = 2 * ceil(Emax_meV / dE_meV) + 1
	Variable dE_eV = dE_meV * 1e-3
	Variable E0_eV = -0.5 * (nE - 1) * dE_eV

	Make/O/D/N=(nE) root:SNS_Settings:tmp_Ldens_E_axis_eV
	Wave E_axis = root:SNS_Settings:tmp_Ldens_E_axis_eV
	E_axis = E0_eV + p * dE_eV
	SetScale/P x, E0_eV, dE_eV, "eV", E_axis

	SNS_MakeThermalKernel(E_axis, params.T_K, "root:SNS_Settings:tmp_Ldens_KT")
	SNS_MakeModulationKernel(E_axis, params.V_mod, "root:SNS_Settings:tmp_Ldens_KM")

	Wave KT = root:SNS_Settings:tmp_Ldens_KT
	Wave KM = root:SNS_Settings:tmp_Ldens_KM

	Duplicate/O KT, root:SNS_Settings:tmp_Ldens_Kexp
	Wave Kexp = root:SNS_Settings:tmp_Ldens_Kexp

	Convolve/A KM, Kexp

	WaveStats/Q Kexp
	if (V_max > 0)
		Kexp /= V_max
	endif

	// Precompute trajectory energies.
	Make/O/D/N=(n) root:SNS_Settings:tmp_Ldens_Ei_meV
	Wave Ei_meV = root:SNS_Settings:tmp_Ldens_Ei_meV

	Variable i
	for (i = 0; i < n; i += 1)
		Ei_meV[i] = SNS_ABSZeroFieldEnergy_meV_FromEll(L_nm[i])
	endfor

	// -------------------------------------------------------------------------
	// Compute rho_exp(ell_c)
	// -------------------------------------------------------------------------
	Make/O/D/N=(nL) $outName
	Wave rho = $outName
	SetScale/P x, lMin_nm, dl_nm, "nm", rho
	rho = 0

	Variable j, ellc, Ec_meV, dEi_meV, dEi_eV, wt

	for (j = 0; j < nL; j += 1)
		ellc = lMin_nm + j * dl_nm
		Ec_meV = SNS_ABSZeroFieldEnergy_meV_FromEll(ellc)

		if (numtype(Ec_meV) != 0)
			rho[j] = NaN
			continue
		endif

		for (i = 0; i < n; i += 1)
			if (numtype(Ei_meV[i]) != 0)
				continue
			endif

			if (hasWeights)
				wt = weightWave[i]
				if (numtype(wt) != 0 || wt < 0)
					wt = 0
				endif
			else
				wt = 1
			endif

			dEi_meV = Ei_meV[i] - Ec_meV
			dEi_eV = dEi_meV * 1e-3

			if (abs(dEi_meV) <= Emax_meV)
				rho[j] += wt * interp(dEi_eV, E_axis, Kexp)
			endif
		endfor
	endfor

	if (normalize)
		WaveStats/Q rho
		if (V_max > 0)
			rho /= V_max
		endif
	endif

	Variable/G $(outName + "_dl_nm") = dl_nm
	Variable/G $(outName + "_lMin_nm") = lMin_nm
	Variable/G $(outName + "_lMax_nm") = lMax_nm
	Variable/G $(outName + "_T_K") = params.T_K
	Variable/G $(outName + "_Vmod_V") = params.V_mod
	Variable/G $(outName + "_Delta_eff_meV") = 1e3 * params.Delta

	KillWaves/Z root:SNS_Settings:tmp_Ldens_E_axis_eV
	KillWaves/Z root:SNS_Settings:tmp_Ldens_KT
	KillWaves/Z root:SNS_Settings:tmp_Ldens_KM
	KillWaves/Z root:SNS_Settings:tmp_Ldens_Kexp
	KillWaves/Z root:SNS_Settings:tmp_Ldens_Ei_meV

	if (doDisplay)
		Display/K=1 rho
		ModifyGraph mode=0
		ModifyGraph tick=2, mirror=2, standoff=0
		Label bottom "candidate ℓ\\Bc\\M (nm)"
		Label left "resolution-broadened zero-field weight"
	endif

	return 0
End


//==============================================================================
// SNS_FindHighLPeakInResolutionDensity
//
// Finds the largest-ell peak of rho(ell_c), using Igor FindPeak.
// rho x-scaling must be ell_c in nm.
//==============================================================================
Function SNS_FindHighLPeakInResolutionDensity(rho, [minRelHeight, box, doDisplay])
	Wave rho
	Variable minRelHeight, box, doDisplay

	if (ParamIsDefault(minRelHeight))
		minRelHeight = 0.05
	endif
	if (ParamIsDefault(box))
		box = 3
	endif
	if (ParamIsDefault(doDisplay))
		doDisplay = 1
	endif

	Variable n = numpnts(rho)
	if (n < 3)
		Abort "SNS_FindHighLPeakInResolutionDensity: rho wave too short."
	endif

	Variable x0 = DimOffset(rho, 0)
	Variable dx = DimDelta(rho, 0)
	Variable xEnd = x0 + (n - 1) * dx

	WaveStats/Q rho
	Variable globalMax = V_max
	if (globalMax <= 0 || numtype(globalMax) != 0)
		return NaN
	endif

	Variable minLevel = minRelHeight * globalMax

	// Search high ell -> low ell.
	FindPeak/Q/B=(box)/M=(minLevel)/R=(xEnd, x0) rho

	if (V_flag != 0 || numtype(V_PeakLoc) != 0)
		SNS_Log("SNS_FindHighLPeakInResolutionDensity: no high-ell peak found.", level="WARN")
		return NaN
	endif

	Variable lPeak_nm = V_PeakLoc
	Variable peakHeight = V_PeakVal
	Variable edge1_nm = V_LeadingEdgeLoc
	Variable edge2_nm = V_TrailingEdgeLoc

	String base = NameOfWave(rho)
	Variable/G $(base + "_lPeak_nm") = lPeak_nm
	Variable/G $(base + "_peakHeight") = peakHeight
	Variable/G $(base + "_peakLow_nm") = min(edge1_nm, edge2_nm)
	Variable/G $(base + "_peakHigh_nm") = max(edge1_nm, edge2_nm)
	Variable/G $(base + "_peakWidth_nm") = V_PeakWidth
	Variable/G $(base + "_minRelHeight") = minRelHeight
	Variable/G $(base + "_FindPeak_box") = box

	if (doDisplay)
		Display/K=1 rho
		ModifyGraph mode=0
		ModifyGraph tick=2, mirror=2, standoff=0
		Label bottom "candidate ℓ\\Bc\\M (nm)"
		Label left "resolution-broadened zero-field weight"

		SetDrawLayer UserFront
		SetDrawEnv xcoord=bottom, ycoord=left, linefgc=(65535,0,0), linethick=2, dash=3
		DrawLine lPeak_nm, 0, lPeak_nm, 1.05 * peakHeight
	endif

	return lPeak_nm
End


//==============================================================================
// SNS_MakeLGapRayContributions
//
// At ell_gap, compute each ray's zero-field contribution:
//
//   A_i = weight_i * K_exp[ E(ell_i) - E(ell_gap) ]
//
// Creates:
//   outBase + "_contribution"
//   outBase + "_family_mask"
//   outBase + "_alpha"
//==============================================================================
Function SNS_MakeLGapRayContributions(L_nm, ellGap_nm, outBase, [weightWave, Emax_meV, dE_meV, alphaMin, alphaMax, alphaGamma])
	Wave L_nm
	Variable ellGap_nm
	String outBase
	Wave weightWave
	Variable Emax_meV, dE_meV, alphaMin, alphaMax, alphaGamma

	STRUCT SNS_Params params
	SNS_LoadParams(params)

	if (ParamIsDefault(Emax_meV))
		Emax_meV = 3
	endif
	if (ParamIsDefault(dE_meV))
		dE_meV = 0.001
	endif
	if (ParamIsDefault(alphaMin))
		alphaMin = 0.10
	endif
	if (ParamIsDefault(alphaMax))
		alphaMax = 1.00
	endif
	if (ParamIsDefault(alphaGamma))
		alphaGamma = 3.0
	endif

	alphaMin = max(0, min(1, alphaMin))
	alphaMax = max(0, min(1, alphaMax))
	alphaGamma = max(0.1, alphaGamma)

	Variable n = numpnts(L_nm)
	if (n <= 0)
		Abort "SNS_MakeLGapRayContributions: empty L_nm wave."
	endif

	Variable hasWeights = !ParamIsDefault(weightWave) && WaveExists(weightWave) && (numpnts(weightWave) == n)

	if (numtype(ellGap_nm) != 0 || ellGap_nm <= 0)
		Abort "SNS_MakeLGapRayContributions: invalid ellGap_nm."
	endif

	Variable EGap_meV = SNS_ABSZeroFieldEnergy_meV_FromEll(ellGap_nm)
	if (numtype(EGap_meV) != 0)
		Abort "SNS_MakeLGapRayContributions: could not compute E(ellGap)."
	endif

	// Build experimental kernel.
	Variable nE = 2 * ceil(Emax_meV / dE_meV) + 1
	Variable dE_eV = dE_meV * 1e-3
	Variable E0_eV = -0.5 * (nE - 1) * dE_eV

	Make/O/D/N=(nE) root:SNS_Settings:tmp_Lgap_E_axis_eV
	Wave E_axis = root:SNS_Settings:tmp_Lgap_E_axis_eV
	E_axis = E0_eV + p * dE_eV
	SetScale/P x, E0_eV, dE_eV, "eV", E_axis

	SNS_MakeThermalKernel(E_axis, params.T_K, "root:SNS_Settings:tmp_Lgap_KT")
	SNS_MakeModulationKernel(E_axis, params.V_mod, "root:SNS_Settings:tmp_Lgap_KM")

	Wave KT = root:SNS_Settings:tmp_Lgap_KT
	Wave KM = root:SNS_Settings:tmp_Lgap_KM

	Duplicate/O KT, root:SNS_Settings:tmp_Lgap_Kexp
	Wave Kexp = root:SNS_Settings:tmp_Lgap_Kexp

	Convolve/A KM, Kexp

	WaveStats/Q Kexp
	if (V_max > 0)
		Kexp /= V_max
	endif

	Make/O/D/N=(n) $(outBase + "_contribution")
	Make/O/D/N=(n) $(outBase + "_family_mask")
	Make/O/D/N=(n) $(outBase + "_alpha")

	Wave contrib = $(outBase + "_contribution")
	Wave mask = $(outBase + "_family_mask")
	Wave alpha = $(outBase + "_alpha")

	contrib = 0
	mask = 0
	alpha = 0

	Variable i, Ei_meV, dEi_meV, dEi_eV, wt

	for (i = 0; i < n; i += 1)
		Ei_meV = SNS_ABSZeroFieldEnergy_meV_FromEll(L_nm[i])
		if (numtype(Ei_meV) != 0)
			continue
		endif

		if (hasWeights)
			wt = weightWave[i]
			if (numtype(wt) != 0 || wt < 0)
				wt = 0
			endif
		else
			wt = 1
		endif

		dEi_meV = Ei_meV - EGap_meV
		dEi_eV = dEi_meV * 1e-3

		if (abs(dEi_meV) <= Emax_meV)
			contrib[i] = wt * interp(dEi_eV, E_axis, Kexp)
		endif
	endfor

	WaveStats/Q contrib
	Variable Amax = V_max
	Variable Ahwhm = 0.5 * Amax
	Variable alphaFrac

	if (Amax > 0 && numtype(Amax) == 0)
		for (i = 0; i < n; i += 1)
			if (contrib[i] >= Ahwhm)
				mask[i] = 1

				if (Amax > Ahwhm)
					alphaFrac = (contrib[i] - Ahwhm) / (Amax - Ahwhm)
					alphaFrac = max(0, min(1, alphaFrac))
					alpha[i] = alphaMin + (alphaMax - alphaMin) * alphaFrac^alphaGamma
				else
					alpha[i] = alphaMax
				endif

				alpha[i] = max(alphaMin, min(alphaMax, alpha[i]))
			endif
		endfor
	endif

	Variable/G $(outBase + "_ellGap_nm") = ellGap_nm
	Variable/G $(outBase + "_EGap_meV") = EGap_meV
	Variable/G $(outBase + "_Amax") = Amax
	Variable/G $(outBase + "_Ahwhm") = Ahwhm
	Variable/G $(outBase + "_alphaMin") = alphaMin
	Variable/G $(outBase + "_alphaMax") = alphaMax
	Variable/G $(outBase + "_alphaGamma") = alphaGamma

	KillWaves/Z root:SNS_Settings:tmp_Lgap_E_axis_eV
	KillWaves/Z root:SNS_Settings:tmp_Lgap_KT
	KillWaves/Z root:SNS_Settings:tmp_Lgap_KM
	KillWaves/Z root:SNS_Settings:tmp_Lgap_Kexp

	return 0
End




//==============================================================================
// SNS_Show3DMaskSideViews_WithChannels
//
// Purpose:
//   Create side-view silhouettes of the 3D island mask and overlay the same
//   selected 3D channels used in the top view.
//
// Side-view convention:
//   x = coordinate parallel to in-plane B
//   y = coordinate perpendicular to in-plane B
//   z = vertical coordinate, 0 at S plane, H_nm at top surface
//
// Display geometry:
//   Axis ranges are set by the island outline only, rounded outward to the next
//   integer nm. No extra padding is added.
//   The plotted z scale is exaggerated by zScaleExag = 10 relative to x/y.
//
// Overlay convention:
//   red  : longest trajectory
//   blue : largest geometric W trajectory
//==============================================================================
Function SNS_Show3DMaskSideViews_WithChannels(Nmask, r0x, r0y, thetaList, phiList, idxLongest, idxLargestW, Bangle_deg, H_nm, maxPath_nm)
    Wave Nmask
    Variable r0x, r0y
    Wave thetaList, phiList
    Variable idxLongest, idxLargestW
    Variable Bangle_deg, H_nm, maxPath_nm

    // -------------------------------------------------------------------------
    // 1. Create side-view mask silhouettes.
    // -------------------------------------------------------------------------
    SNS_Make3DMaskSideViews(Nmask, r0x, r0y, Bangle_deg, H_nm, "MaskSide_XZ", "MaskSide_YZ")

    Wave MaskSide_XZ
    Wave MaskSide_YZ

    // -------------------------------------------------------------------------
    // 2. Display x-z side view.
    // -------------------------------------------------------------------------
    SNS_AddSideMaskOutline(MaskSide_XZ, "SNS_3DMaskSide_XZ", "MaskSide_XZ", doNew=1)

    ModifyGraph/W=SNS_3DMaskSide_XZ tick=2, mirror=2, standoff=0
    Label/W=SNS_3DMaskSide_XZ bottom "x\\B∥B\\M (nm)"
    Label/W=SNS_3DMaskSide_XZ left   "z (nm)"

    // -------------------------------------------------------------------------
    // 3. Display y-z side view.
    // -------------------------------------------------------------------------
    SNS_AddSideMaskOutline(MaskSide_YZ, "SNS_3DMaskSide_YZ", "MaskSide_YZ", doNew=1)

    ModifyGraph/W=SNS_3DMaskSide_YZ tick=2, mirror=2, standoff=0
    Label/W=SNS_3DMaskSide_YZ bottom "y\\B⊥B\\M (nm)"
    Label/W=SNS_3DMaskSide_YZ left   "z (nm)"

    // -------------------------------------------------------------------------
    // 4. Overlay selected channels.
    // -------------------------------------------------------------------------
    Variable errPlot

    // Red: longest trajectory.
    errPlot = SNS_PlotChannel3D_OnSideMasks(Nmask, r0x, r0y, thetaList, phiList, \
        idxLongest, 0.5, H_nm, maxPath_nm, Bangle_deg, \
        "SNS_3DMaskSide_XZ", "SNS_3DMaskSide_YZ", "ray3D_long", 65535, 0, 0)

    // Blue: largest geometric W trajectory.
    errPlot = SNS_PlotChannel3D_OnSideMasks(Nmask, r0x, r0y, thetaList, phiList, \
        idxLargestW, 0.5, H_nm, maxPath_nm, Bangle_deg, \
        "SNS_3DMaskSide_XZ", "SNS_3DMaskSide_YZ", "ray3D_wmax", 0, 0, 65535)

    // -------------------------------------------------------------------------
    // 5. Axis ranges from side-view outlines only, rounded to integer nm.
    // -------------------------------------------------------------------------
    Variable xzMin = Inf
    Variable xzMax = -Inf
    Variable yzMin = Inf
    Variable yzMax = -Inf
    Variable zMin = Inf
    Variable zMax = -Inf

    Wave/Z xzOutlineX = MaskSide_XZ_outline_x
    Wave/Z xzOutlineZ = MaskSide_XZ_outline_z
    Wave/Z yzOutlineX = MaskSide_YZ_outline_x
    Wave/Z yzOutlineZ = MaskSide_YZ_outline_z

    Variable i

    if (WaveExists(xzOutlineX) && WaveExists(xzOutlineZ))
        for (i = 0; i < numpnts(xzOutlineX); i += 1)
            if (numtype(xzOutlineX[i]) == 0 && numtype(xzOutlineZ[i]) == 0)
                xzMin = min(xzMin, xzOutlineX[i])
                xzMax = max(xzMax, xzOutlineX[i])
                zMin  = min(zMin,  xzOutlineZ[i])
                zMax  = max(zMax,  xzOutlineZ[i])
            endif
        endfor
    endif

    if (WaveExists(yzOutlineX) && WaveExists(yzOutlineZ))
        for (i = 0; i < numpnts(yzOutlineX); i += 1)
            if (numtype(yzOutlineX[i]) == 0 && numtype(yzOutlineZ[i]) == 0)
                yzMin = min(yzMin, yzOutlineX[i])
                yzMax = max(yzMax, yzOutlineX[i])
                zMin  = min(zMin,  yzOutlineZ[i])
                zMax  = max(zMax,  yzOutlineZ[i])
            endif
        endfor
    endif

    if (numtype(xzMin) != 0)
        xzMin = -100
        xzMax = 100
    endif

    if (numtype(yzMin) != 0)
        yzMin = -50
        yzMax = 50
    endif

    if (numtype(zMin) != 0)
        zMin = 0
        zMax = H_nm
    endif

	// Add 10% padding in horizontal directions.
	// Keep z = 0 axis-aligned, add padding only above the island.
	Variable xzPad = 0.10 * (xzMax - xzMin)
	Variable yzPad = 0.10 * (yzMax - yzMin)
	Variable zPad  = 0.10 * (zMax - zMin)
	
	xzMin = floor(xzMin - xzPad)
	xzMax = ceil(xzMax + xzPad)
	
	yzMin = floor(yzMin - yzPad)
	yzMax = ceil(yzMax + yzPad)
	
	// Keep lower z boundary fixed at 0 if possible.
	zMin = 0
	zMax = ceil(zMax + zPad)
	
	Variable xzSpan = xzMax - xzMin
	Variable yzSpan = yzMax - yzMin
	Variable zSpan  = zMax  - zMin

    if (!(xzSpan > 0))
        xzSpan = 1
        xzMax = xzMin + 1
    endif

    if (!(yzSpan > 0))
        yzSpan = 1
        yzMax = yzMin + 1
    endif

    if (!(zSpan > 0))
        zSpan = 1
        zMax = zMin + 1
    endif

    SetAxis/W=SNS_3DMaskSide_XZ bottom, xzMin, xzMax
    SetAxis/W=SNS_3DMaskSide_XZ left,   zMin,  zMax

    SetAxis/W=SNS_3DMaskSide_YZ bottom, yzMin, yzMax
    SetAxis/W=SNS_3DMaskSide_YZ left,   zMin,  zMax

   // -------------------------------------------------------------------------
	// 6. Graph sizes with controlled z exaggeration.
	//
	// Important:
	//   Do NOT impose a minimum panel height/width here. A minimum size destroys
	//   the requested physical scale ratio and makes z appear exaggerated even
	//   when zScaleExag = 1.
	//
	// Visual scale:
	//   horizontal display scale ∝ 1 nm
	//   vertical display scale   ∝ zScaleExag * 1 nm
	//
	// Therefore:
	//   XZ width  / XZ height = xzSpan / (zSpan*zScaleExag)
	//   YZ width  / YZ height = yzSpan / (zSpan*zScaleExag)
	//
	// Use zScaleExag = 1 for true-angle debugging.
	// Use zScaleExag = 10 for the usual publication-style vertical exaggeration.
	// -------------------------------------------------------------------------
	Variable zScaleExag = 10
	// Variable zScaleExag = 1     // debug true x/y/z aspect
	
	Variable zPlotSpan = zSpan * zScaleExag
	
	Variable refSpan = max(max(xzSpan, yzSpan), zPlotSpan)
	Variable targetMaxPx = 320
	Variable pxPerNm = targetMaxPx / refSpan
	
	Variable xzWidthPx = round(xzSpan * pxPerNm)
	Variable yzWidthPx = round(yzSpan * pxPerNm)
	Variable zHeightPx = round(zPlotSpan * pxPerNm)
	
	// Only protect against pathological zero-size windows.
	xzWidthPx = max(10, xzWidthPx)
	yzWidthPx = max(10, yzWidthPx)
	zHeightPx = max(10, zHeightPx)
	
	ModifyGraph/W=SNS_3DMaskSide_XZ width=(xzWidthPx), height=(zHeightPx)
	ModifyGraph/W=SNS_3DMaskSide_YZ width=(yzWidthPx), height=(zHeightPx)
    // -------------------------------------------------------------------------
    // 7. Legends.
    // -------------------------------------------------------------------------
    DoWindow/F SNS_3DMaskSide_XZ
    Legend/C/N=text0/J/F=0/A=RT "Side view\r" + \
        "\\s(ray3D_long_Z_xz_" + num2str(idxLongest) + ") longest\r" + \
        "\\s(ray3D_wmax_Z_xz_" + num2str(idxLargestW) + ") largest W\\Bgeom\\M"

    DoWindow/F SNS_3DMaskSide_YZ
    Legend/C/N=text0/J/F=0/A=RT "Side view\r" + \
        "\\s(ray3D_long_Z_yz_" + num2str(idxLongest) + ") longest\r" + \
        "\\s(ray3D_wmax_Z_yz_" + num2str(idxLargestW) + ") largest W\\Bgeom\\M"

    return 0
End

//==============================================================================
// SNS_Make3DMaskSideViews
//
// Purpose:
//   Project a 2D N-region mask extruded from z=0 to z=H_nm into two side-view
//   silhouettes:
//
//      XZ view: x_parallel_B vs z
//      YZ view: y_perp_B     vs z
//
//   The side masks are binary silhouettes. They are not trajectory-density maps.
//==============================================================================
Function SNS_Make3DMaskSideViews(Nmask, r0x, r0y, Bangle_deg, H_nm, outXZName, outYZName)
    Wave Nmask
    Variable r0x, r0y
    Variable Bangle_deg, H_nm
    String outXZName, outYZName

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    Variable phiB = Bangle_deg*pi/180
    Variable bx = cos(phiB)
    Variable by = sin(phiB)
    Variable epx = -sin(phiB)
    Variable epy =  cos(phiB)

    Variable i, j
    Variable x, y, xB, yP

    Variable xBmin = Inf
    Variable xBmax = -Inf
    Variable yPmin = Inf
    Variable yPmax = -Inf

    // Determine projected coordinate ranges from inside-mask pixels.
    for (i = 0; i < nx; i += 1)
        x = x0 + i*dx
        for (j = 0; j < ny; j += 1)
            if (Nmask[i][j] > 0.5)
                y = y0 + j*dy

                xB = (x-r0x)*bx  + (y-r0y)*by
                yP = (x-r0x)*epx + (y-r0y)*epy

                xBmin = min(xBmin, xB)
                xBmax = max(xBmax, xB)
                yPmin = min(yPmin, yP)
                yPmax = max(yPmax, yP)
            endif
        endfor
    endfor

    if (numtype(xBmin) != 0 || numtype(yPmin) != 0)
        Abort "SNS_Make3DMaskSideViews: no inside pixels found in Nmask."
    endif

    Variable dxSide = min(abs(dx), abs(dy))
    Variable dzSide = dxSide

    Variable nXB = max(2, ceil((xBmax-xBmin)/dxSide) + 1)
    Variable nYP = max(2, ceil((yPmax-yPmin)/dxSide) + 1)
    Variable nZ  = max(2, ceil(H_nm/dzSide) + 1)
	 Variable dzZ = H_nm/(nZ-1)

    Make/O/D/N=(nXB,nZ) $outXZName
    Make/O/D/N=(nYP,nZ) $outYZName

    Wave sideXZ = $outXZName
    Wave sideYZ = $outYZName

    sideXZ = 0
    sideYZ = 0

    SetScale/P x, xBmin, dxSide, "nm", sideXZ
    SetScale/P y, 0, dzZ, "nm", sideXZ


    SetScale/P x, yPmin, dxSide, "nm", sideYZ
    SetScale/P y, 0, dzZ, "nm", sideYZ

    Variable ixB, iyP

    // Fill the full z-column for every projected inside-mask coordinate.
    for (i = 0; i < nx; i += 1)
        x = x0 + i*dx
        for (j = 0; j < ny; j += 1)
            if (Nmask[i][j] > 0.5)
                y = y0 + j*dy

                xB = (x-r0x)*bx  + (y-r0y)*by
                yP = (x-r0x)*epx + (y-r0y)*epy

                ixB = round((xB - xBmin)/dxSide)
                iyP = round((yP - yPmin)/dxSide)

                if (ixB >= 0 && ixB < nXB)
                    sideXZ[ixB][] = 1
                endif

                if (iyP >= 0 && iyP < nYP)
                    sideYZ[iyP][] = 1
                endif
            endif
        endfor
    endfor

    return 0
End

//==============================================================================
// SNS_PlotChannel3D_OnSideMasks
//
// Purpose:
//   Overlay one selected 3D channel on the x-z and y-z side-view mask windows.
//==============================================================================
Function SNS_PlotChannel3D_OnSideMasks(Nmask, r0x, r0y, thetaList, phiList, idx, stepFac, H_nm, maxPath_nm, Bangle_deg, winXZ, winYZ, rayBaseName, r, g, b)
    Wave Nmask
    Variable r0x, r0y
    Wave thetaList, phiList
    Variable idx
    Variable stepFac, H_nm, maxPath_nm, Bangle_deg
    String winXZ, winYZ, rayBaseName
    Variable r, g, b

    if (idx < 0 || idx >= numpnts(thetaList))
        return -1
    endif

    Variable theta = thetaList[idx]
    Variable phi   = phiList[idx]

    String rayXName = rayBaseName + "_X_" + num2str(idx)
    String rayYName = rayBaseName + "_Y_" + num2str(idx)
    String rayZName = rayBaseName + "_Z_" + num2str(idx)

    Variable err = SNS_ProjectChannel3D_SurfaceXYZ(Nmask, r0x, r0y, theta, phi, \
        stepFac, H_nm, maxPath_nm, rayXName, rayYName, rayZName)

    if (err != 0)
        return err
    endif

    Wave rayX = $rayXName
    Wave rayY = $rayYName
    Wave rayZ = $rayZName

    Variable n = numpnts(rayX)

    String rayXBName  = rayBaseName + "_XB_" + num2str(idx)
    String rayYPName  = rayBaseName + "_YP_" + num2str(idx)
    String rayZXZName = rayBaseName + "_Z_xz_" + num2str(idx)
    String rayZYZName = rayBaseName + "_Z_yz_" + num2str(idx)

    Make/O/D/N=(n) $rayXBName, $rayYPName, $rayZXZName, $rayZYZName

    Wave rayXB  = $rayXBName
    Wave rayYP  = $rayYPName
    Wave rayZXZ = $rayZXZName
    Wave rayZYZ = $rayZYZName

    Variable phiB = Bangle_deg*pi/180
    Variable bx = cos(phiB)
    Variable by = sin(phiB)
    Variable epx = -sin(phiB)
    Variable epy =  cos(phiB)

    rayXB = (rayX-r0x)*bx  + (rayY-r0y)*by
    rayYP = (rayX-r0x)*epx + (rayY-r0y)*epy

    rayZXZ = rayZ
    rayZYZ = rayZ

    AppendToGraph/W=$winXZ rayZXZ vs rayXB
    ModifyGraph/W=$winXZ mode($NameOfWave(rayZXZ))=0
    ModifyGraph/W=$winXZ lsize($NameOfWave(rayZXZ))=2
    ModifyGraph/W=$winXZ rgb($NameOfWave(rayZXZ))=(r,g,b)

    AppendToGraph/W=$winYZ rayZYZ vs rayYP
    ModifyGraph/W=$winYZ mode($NameOfWave(rayZYZ))=0
    ModifyGraph/W=$winYZ lsize($NameOfWave(rayZYZ))=2
    ModifyGraph/W=$winYZ rgb($NameOfWave(rayZYZ))=(r,g,b)

    return 0
End

//==============================================================================
// SNS_ProjectChannel3D_SurfaceXYZ
//
// Purpose:
//   Display/retrace wrapper. It does not contain independent ray physics.
//   It calls SNS_TraceOneChannel3D_Cu111 with savePath=1.
//
// Notes:
//   Uses legal Igor temporary wave name trace3D_plot_res_tmp.
//==============================================================================
Function SNS_ProjectChannel3D_SurfaceXYZ(Nmask, r0x, r0y, theta, phi, stepFac, H_nm, maxPath_nm, rayXName, rayYName, rayZName)
    Wave    Nmask
    Variable r0x, r0y
    Variable theta, phi
    Variable stepFac, H_nm, maxPath_nm
    String  rayXName, rayYName, rayZName

    // phiB does not affect tracing geometry, only W_eff in the result.
    // Use 0 here because this wrapper only needs rayX/rayY/rayZ.
    Variable phiB_dummy = 0
    Variable Zbarrier_dummy = 0

    Variable err = SNS_TraceOneChannel3D_Cu111(Nmask, r0x, r0y, theta, phi, phiB_dummy, stepFac, H_nm, maxPath_nm, "trace3D_plot_res_tmp", 1, rayXName, rayYName, rayZName)

    KillWaves/Z trace3D_plot_res_tmp

    return err
End

//==============================================================================
// SNS_AddSideMaskOutline
//
// Purpose:
//   Add or display a sloped side-view outline for an extruded Cu(111) island.
//
//   The top silhouette is defined by the projected 2D mask width.
//   The bottom silhouette is widened according to close-packed fcc Cu(111)
//   stacking:
//
//      d111 = a_Cu/sqrt(3)
//      lateral shift per layer = a_Cu/sqrt(6)
//      tan(alpha) = 1/sqrt(2)
//      alpha = 35.264 deg from vertical
//
//   Thus, for island height H:
//      lateral broadening per side = H/sqrt(2)
//
// Inputs:
//   wSide   : side-view mask wave.
//             x scaling = projected in-plane coordinate [nm]
//             y scaling = z [nm]
//   winName : graph window name.
//   baseName: base name for outline waves.
//   doNew   : 1 -> create new graph using Display
//             0 -> append to existing graph
//==============================================================================
Function SNS_AddSideMaskOutline(wSide, winName, baseName, [doNew])
    Wave wSide
    String winName, baseName
    Variable doNew

    if (ParamIsDefault(doNew))
        doNew = 0
    endif

    Variable nx = DimSize(wSide, 0)
    Variable nz = DimSize(wSide, 1)

    if (nx < 2 || nz < 2)
        return -1
    endif

    Variable xTopMin = DimOffset(wSide, 0)
    Variable dx      = DimDelta(wSide, 0)
    Variable xTopMax = xTopMin + dx*(nx-1)

    Variable zMin = DimOffset(wSide, 1)
    Variable dz   = DimDelta(wSide, 1)
    Variable zMax = zMin + dz*(nz-1)

	 Variable H_nm = zMax - zMin
	 Variable tanAlpha = 0
	 if (H_nm > 0)
	    tanAlpha = SNS_Cu111PrismSideExpand_nm(zMin, H_nm) / H_nm
	 endif
    Variable sideExpand_nm = H_nm * tanAlpha

    Variable xBotMin = xTopMin - sideExpand_nm
    Variable xBotMax = xTopMax + sideExpand_nm

    String xName = baseName + "_outline_x"
    String zName = baseName + "_outline_z"

    Make/O/D/N=5 $xName, $zName
    Wave xOutline = $xName
    Wave zOutline = $zName

    // Trapezoid: bottom wider than top.
    xOutline[0] = xBotMin
    zOutline[0] = zMin

    xOutline[1] = xBotMax
    zOutline[1] = zMin

    xOutline[2] = xTopMax
    zOutline[2] = zMax

    xOutline[3] = xTopMin
    zOutline[3] = zMax

    xOutline[4] = xBotMin
    zOutline[4] = zMin

    if (doNew)
        Display/K=1 zOutline vs xOutline
        DoWindow/C $winName
    else
        AppendToGraph/W=$winName zOutline vs xOutline
    endif

    ModifyGraph/W=$winName mode($NameOfWave(zOutline))=0
    ModifyGraph/W=$winName rgb($NameOfWave(zOutline))=(0,0,0)
    ModifyGraph/W=$winName lsize($NameOfWave(zOutline))=1.5

    return 0
End


//==============================================================================
// SNS_Cu111PrismSideExpand_nm
//
// Purpose:
//   Lateral expansion of a Cu(111) close-packed island side wall at height z.
//
// Geometry:
//   z = H_nm : measured top mask
//   z = 0    : bottom/S plane
//
//   fcc Cu(111) close-packed stacking gives:
//      lateral/vertical = 1/sqrt(2)
//      side-wall angle = 35.264 deg from vertical
//
// Return:
//   lateral expansion relative to the top mask [nm].
//==============================================================================
Function SNS_Cu111PrismSideExpand_nm(z_nm, H_nm)
    Variable z_nm, H_nm

    Variable tanAlpha = 1/sqrt(2)

    if (H_nm <= 0)
        return 0
    endif

    Variable zClamped = max(0, min(H_nm, z_nm))
    return (H_nm - zClamped) * tanAlpha
End


//==============================================================================
// SNS_IsInsideCu111Prism
//
// Purpose:
//   Test whether a point (x,y,z) is inside the Cu(111)-faceted prism.
//
//   The measured 2D mask is interpreted as the top mask at z=H_nm.
//   At lower z, the allowed in-plane region is the top mask expanded outward by
//   (H_nm-z)/sqrt(2).
//
// Return:
//   1 if inside, 0 if outside.
//==============================================================================
Function SNS_IsInsideCu111Prism(Nmask, x_nm, y_nm, z_nm, H_nm)
    Wave Nmask
    Variable x_nm, y_nm, z_nm, H_nm

    if (z_nm < -1e-9 || z_nm > H_nm + 1e-9)
        return 0
    endif

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    Variable pix = min(abs(dx), abs(dy))
    Variable expand_nm = SNS_Cu111PrismSideExpand_nm(z_nm, H_nm)

    Variable ix0 = round((x_nm - x0)/dx)
    Variable iy0 = round((y_nm - y0)/dy)

    // At the top surface, use the mask directly.
    if (expand_nm <= 0.5*pix)
        if (ix0 < 0 || ix0 >= nx || iy0 < 0 || iy0 >= ny)
            return 0
        endif
        return (Nmask[ix0][iy0] > 0.5)
    endif

    // Expanded lower slices: accept if the point lies within expand_nm of any
    // top-mask pixel.
    Variable rPix = ceil(expand_nm/pix) + 2

    Variable ixMin = max(0, ix0-rPix)
    Variable ixMax = min(nx-1, ix0+rPix)
    Variable iyMin = max(0, iy0-rPix)
    Variable iyMax = min(ny-1, iy0+rPix)

    if (ixMin > ixMax || iyMin > iyMax)
        return 0
    endif

    Variable i, j
    Variable xPix, yPix, rr
    Variable rCut = expand_nm + 0.75*pix
    Variable rCut2 = rCut*rCut

    for (i = ixMin; i <= ixMax; i += 1)
        xPix = x0 + i*dx
        for (j = iyMin; j <= iyMax; j += 1)
            if (Nmask[i][j] > 0.5)
                yPix = y0 + j*dy
                rr = (x_nm-xPix)^2 + (y_nm-yPix)^2
                if (rr <= rCut2)
                    return 1
                endif
            endif
        endfor
    endfor

    return 0
End


//==============================================================================
// SNS_Cu111PrismNormalComponent
//
// Purpose:
//   Return one component of the outward normal of the actual 3D prism boundary.
//
//   comp = 0 : nx
//   comp = 1 : ny
//   comp = 2 : nz
//
// Method:
//   This version does NOT estimate the normal from nearest top-mask pixels.
//   Instead, it numerically differentiates the actual inside/outside function:
//
//      SNS_IsInsideCu111Prism(Nmask, x, y, z, H_nm)
//
//   around the boundary hit point.
//
//   Since the inside indicator is 1 inside and 0 outside,
//   grad(I) points inward, so:
//
//      n_out = -grad(I)
//
//   This guarantees that the reflection normal is tied to the same geometry
//   used by the boundary-crossing test.
//
// Fallback:
//   If the binary finite-difference stencil fails, fall back to an analytic
//   side-wall normal using a nearest-boundary-pixel direction.
//==============================================================================
Function SNS_Cu111PrismNormalComponent(Nmask, x_nm, y_nm, z_nm, H_nm, comp)
    Wave Nmask
    Variable x_nm, y_nm, z_nm, H_nm, comp

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    Variable pix = min(abs(dx), abs(dy))
    if (!(pix > 0))
        pix = 1
    endif

    // -------------------------------------------------------------------------
    // 1. Preferred normal: finite difference of actual 3D inside-test.
    // -------------------------------------------------------------------------
    Variable nxOut = NaN
    Variable nyOut = NaN
    Variable nzOut = NaN

    Variable scale, epsXY, epsZ
    Variable ip, im
    Variable gx, gy, gz, nlen

    for (scale = 0; scale < 6; scale += 1)

        epsXY = pix * (0.5 + scale)
        epsZ  = epsXY

        // x derivative of inside indicator.
        ip = SNS_IsInsideCu111Prism(Nmask, x_nm + epsXY, y_nm, z_nm, H_nm)
        im = SNS_IsInsideCu111Prism(Nmask, x_nm - epsXY, y_nm, z_nm, H_nm)
        gx = (ip - im) / (2*epsXY)

        // y derivative of inside indicator.
        ip = SNS_IsInsideCu111Prism(Nmask, x_nm, y_nm + epsXY, z_nm, H_nm)
        im = SNS_IsInsideCu111Prism(Nmask, x_nm, y_nm - epsXY, z_nm, H_nm)
        gy = (ip - im) / (2*epsXY)

        // z derivative of inside indicator.
        ip = SNS_IsInsideCu111Prism(Nmask, x_nm, y_nm, z_nm + epsZ, H_nm)
        im = SNS_IsInsideCu111Prism(Nmask, x_nm, y_nm, z_nm - epsZ, H_nm)
        gz = (ip - im) / (2*epsZ)

        // inside indicator gradient points inward; outward normal is -grad.
        nxOut = -gx
        nyOut = -gy
        nzOut = -gz

        nlen = sqrt(nxOut*nxOut + nyOut*nyOut + nzOut*nzOut)

        if (nlen > 0)
            nxOut /= nlen
            nyOut /= nlen
            nzOut /= nlen

            if (comp == 0)
                return nxOut
            elseif (comp == 1)
                return nyOut
            else
                return nzOut
            endif
        endif
    endfor

    // -------------------------------------------------------------------------
    // 2. Fallback: nearest boundary pixel + analytic Cu/prism slope.
    // -------------------------------------------------------------------------
    Variable ix0 = round((x_nm - x0)/dx)
    Variable iy0 = round((y_nm - y0)/dy)

    Variable expand_nm = max(pix, SNS_Cu111PrismSideExpand_nm(z_nm, H_nm))
    Variable rPix = ceil(expand_nm/pix) + 8

    Variable ixMin = max(1, ix0-rPix)
    Variable ixMax = min(nx-2, ix0+rPix)
    Variable iyMin = max(1, iy0-rPix)
    Variable iyMax = min(ny-2, iy0+rPix)

    Variable bestR2 = Inf
    Variable bestIx = -1
    Variable bestIy = -1

    Variable i, j
    Variable xPix, yPix, rr
    Variable isBoundary

    for (i = ixMin; i <= ixMax; i += 1)
        xPix = x0 + i*dx

        for (j = iyMin; j <= iyMax; j += 1)

            if (Nmask[i][j] <= 0.5)
                continue
            endif

            isBoundary = 0
            if (Nmask[i+1][j] <= 0.5)
                isBoundary = 1
            endif
            if (Nmask[i-1][j] <= 0.5)
                isBoundary = 1
            endif
            if (Nmask[i][j+1] <= 0.5)
                isBoundary = 1
            endif
            if (Nmask[i][j-1] <= 0.5)
                isBoundary = 1
            endif

            if (!isBoundary)
                continue
            endif

            yPix = y0 + j*dy
            rr = (x_nm - xPix)^2 + (y_nm - yPix)^2

            if (rr < bestR2)
                bestR2 = rr
                bestIx = i
                bestIy = j
            endif
        endfor
    endfor

    if (bestIx >= 0)

        // Use radial direction from nearest boundary point to hit point.
        xPix = x0 + bestIx*dx
        yPix = y0 + bestIy*dy

        nxOut = x_nm - xPix
        nyOut = y_nm - yPix

        nlen = sqrt(nxOut*nxOut + nyOut*nyOut)

        if (nlen > 0)
            nxOut /= nlen
            nyOut /= nlen
        else
            nxOut = 1
            nyOut = 0
        endif

    else
        nxOut = 1
        nyOut = 0
    endif

    Variable tanAlpha = 0
    if (H_nm > 0)
        tanAlpha = SNS_Cu111PrismSideExpand_nm(0, H_nm) / H_nm
    endif

    nlen = sqrt(1 + tanAlpha*tanAlpha)

    nxOut /= nlen
    nyOut /= nlen
    nzOut = tanAlpha / nlen

    if (comp == 0)
        return nxOut
    elseif (comp == 1)
        return nyOut
    else
        return nzOut
    endif
End
//==============================================================================
// SNS_Show3DMaskTopFieldFrame_WithChannels
//
// Purpose:
//   Create an outline-only top-view mask in field-frame coordinates and overlay
//   the same selected 3D channels used in the side views.
//
// Coordinate convention:
//   x = coordinate parallel to in-plane B
//   y = coordinate perpendicular to in-plane B
//
// Display geometry:
//   Axis ranges are set by the island outline only, rounded outward to the next
//   integer nm. No extra padding is added.
//   The graph size follows the true x/y aspect ratio. For consistency with the
//   side views, the same reference scale includes z exaggerated by a factor 10.
//
// Overlay convention:
//   red  : longest trajectory
//   blue : largest geometric W trajectory
//==============================================================================
Function SNS_Show3DMaskTopFieldFrame_WithChannels(Nmask, r0x, r0y, thetaList, phiList, idxLongest, idxLargestW, Bangle_deg, H_nm, maxPath_nm)
    Wave Nmask
    Variable r0x, r0y
    Wave thetaList, phiList
    Variable idxLongest, idxLargestW
    Variable Bangle_deg, H_nm, maxPath_nm

    // -------------------------------------------------------------------------
    // 1. Display outline-only field-frame top view.
    // -------------------------------------------------------------------------
    SNS_AddTopMaskFieldFrameOutline(Nmask, r0x, r0y, Bangle_deg, "SNS_3DMaskTop_XY", "MaskTop_XY", doNew=1)

    ModifyGraph/W=SNS_3DMaskTop_XY tick=2, mirror=2, standoff=0
    Label/W=SNS_3DMaskTop_XY bottom "x\\B∥B\\M (nm)"
    Label/W=SNS_3DMaskTop_XY left   "y\\B⊥B\\M (nm)"

    // STS marker at origin.
    Make/O/D/N=1 TopField_STS_x = 0
    Make/O/D/N=1 TopField_STS_y = 0
    AppendToGraph/W=SNS_3DMaskTop_XY TopField_STS_y vs TopField_STS_x
    ModifyGraph/W=SNS_3DMaskTop_XY mode(TopField_STS_y)=3
    ModifyGraph/W=SNS_3DMaskTop_XY marker(TopField_STS_y)=19
    ModifyGraph/W=SNS_3DMaskTop_XY msize(TopField_STS_y)=3
    ModifyGraph/W=SNS_3DMaskTop_XY rgb(TopField_STS_y)=(65535,65535,65535)

    // -------------------------------------------------------------------------
    // 2. Overlay selected channels.
    // -------------------------------------------------------------------------
    Variable errPlot

    // Red: longest trajectory.
    errPlot = SNS_PlotChannel3D_OnTopFieldFrame(Nmask, r0x, r0y, thetaList, phiList, \
        idxLongest, 0.5, H_nm, maxPath_nm, Bangle_deg, \
        "SNS_3DMaskTop_XY", "ray3D_long_top", 65535, 0, 0)

    // Blue: largest geometric W trajectory.
    errPlot = SNS_PlotChannel3D_OnTopFieldFrame(Nmask, r0x, r0y, thetaList, phiList, \
        idxLargestW, 0.5, H_nm, maxPath_nm, Bangle_deg, \
        "SNS_3DMaskTop_XY", "ray3D_wmax_top", 0, 0, 65535)

    // -------------------------------------------------------------------------
    // 3. Axis range from outline only, rounded outward to integer nm.
    // -------------------------------------------------------------------------
    Variable xMin = Inf
    Variable xMax = -Inf
    Variable yMin = Inf
    Variable yMax = -Inf

    Wave/Z xOutline = MaskTop_XY_outline_x
    Wave/Z yOutline = MaskTop_XY_outline_y

    Variable i

    if (WaveExists(xOutline) && WaveExists(yOutline))
        for (i = 0; i < numpnts(xOutline); i += 1)
            if (numtype(xOutline[i]) == 0 && numtype(yOutline[i]) == 0)
                xMin = min(xMin, xOutline[i])
                xMax = max(xMax, xOutline[i])
                yMin = min(yMin, yOutline[i])
                yMax = max(yMax, yOutline[i])
            endif
        endfor
    endif

    if (numtype(xMin) != 0 || numtype(yMin) != 0)
        xMin = -100
        xMax = 100
        yMin = -50
        yMax = 50
    endif

 	 // Add 10% padding in x and y, then round outward to integer nm.
 	 Variable xPad = 0.10 * (xMax - xMin)
 	 Variable yPad = 0.10 * (yMax - yMin)
	
	 xMin = floor(xMin - xPad)
	 xMax = ceil(xMax + xPad)
	 yMin = floor(yMin - yPad)
	 yMax = ceil(yMax + yPad)
	
	 Variable xSpan = xMax - xMin
	 Variable ySpan = yMax - yMin

    if (!(xSpan > 0))
        xSpan = 1
        xMax = xMin + 1
    endif

    if (!(ySpan > 0))
        ySpan = 1
        yMax = yMin + 1
    endif

    SetAxis/W=SNS_3DMaskTop_XY bottom, xMin, xMax
    SetAxis/W=SNS_3DMaskTop_XY left,   yMin, yMax

    // -------------------------------------------------------------------------
    // 4. Graph size with correct relative scale.
    //
    // z-scale is exaggerated by 10x for consistency with side views.
    // The top view itself uses true x/y aspect ratio.
    // -------------------------------------------------------------------------
    Variable zScaleExag = 10
    Variable zPlotSpan = max(1, H_nm*zScaleExag)

    Variable refSpan = max(max(xSpan, ySpan), zPlotSpan)
    Variable targetMaxPx = 320
    Variable pxPerNm = targetMaxPx / refSpan

    Variable topWidthPx  = max(60, round(xSpan * pxPerNm))
    Variable topHeightPx = max(60, round(ySpan * pxPerNm))

    ModifyGraph/W=SNS_3DMaskTop_XY width=(topWidthPx), height=(topHeightPx)

    // -------------------------------------------------------------------------
    // 5. Legend.
    // -------------------------------------------------------------------------
    DoWindow/F SNS_3DMaskTop_XY
    Legend/C/N=text0/J/F=0/A=RT "Top view, field frame\r" + \
        "\\s(ray3D_long_top_YP_" + num2str(idxLongest) + ") longest\r" + \
        "\\s(ray3D_wmax_top_YP_" + num2str(idxLargestW) + ") largest W\\Bgeom\\M"

    return 0
End
//==============================================================================
// SNS_MakeTopMaskFieldFrame
//
// Purpose:
//   Rotate/project the 2D mask into field-frame top-view coordinates:
//
//      x = x_parallel_B
//      y = y_perp_B
//
//   The origin is the STM position.
//
// Output:
//      outName : binary mask image in field-frame coordinates.
//==============================================================================
Function SNS_MakeTopMaskFieldFrame(Nmask, r0x, r0y, Bangle_deg, outName)
    Wave Nmask
    Variable r0x, r0y
    Variable Bangle_deg
    String outName

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    Variable phiB = Bangle_deg*pi/180
    Variable bx = cos(phiB)
    Variable by = sin(phiB)
    Variable epx = -sin(phiB)
    Variable epy =  cos(phiB)

    Variable i, j
    Variable x, y, xB, yP

    Variable xBmin = Inf
    Variable xBmax = -Inf
    Variable yPmin = Inf
    Variable yPmax = -Inf

    // Determine field-frame range from inside-mask pixels.
    for (i = 0; i < nx; i += 1)
        x = x0 + i*dx
        for (j = 0; j < ny; j += 1)
            if (Nmask[i][j] > 0.5)
                y = y0 + j*dy

                xB = (x-r0x)*bx  + (y-r0y)*by
                yP = (x-r0x)*epx + (y-r0y)*epy

                xBmin = min(xBmin, xB)
                xBmax = max(xBmax, xB)
                yPmin = min(yPmin, yP)
                yPmax = max(yPmax, yP)
            endif
        endfor
    endfor

    if (numtype(xBmin) != 0 || numtype(yPmin) != 0)
        Abort "SNS_MakeTopMaskFieldFrame: no inside pixels found in Nmask."
    endif

    Variable dSide = min(abs(dx), abs(dy))

    Variable nXB = max(2, ceil((xBmax-xBmin)/dSide) + 1)
    Variable nYP = max(2, ceil((yPmax-yPmin)/dSide) + 1)

    Make/O/D/N=(nXB,nYP) $outName
    Wave topXY = $outName
    topXY = 0

    SetScale/P x, xBmin, dSide, "nm", topXY
    SetScale/P y, yPmin, dSide, "nm", topXY

    Variable ixB, iyP

    // Fill transformed top-view mask.
    for (i = 0; i < nx; i += 1)
        x = x0 + i*dx
        for (j = 0; j < ny; j += 1)
            if (Nmask[i][j] > 0.5)
                y = y0 + j*dy

                xB = (x-r0x)*bx  + (y-r0y)*by
                yP = (x-r0x)*epx + (y-r0y)*epy

                ixB = round((xB - xBmin)/dSide)
                iyP = round((yP - yPmin)/dSide)

                if (ixB >= 0 && ixB < nXB && iyP >= 0 && iyP < nYP)
                    topXY[ixB][iyP] = 1
                endif
            endif
        endfor
    endfor

    return 0
End

//==============================================================================
// SNS_PlotChannel3D_OnTopFieldFrame
//
// Purpose:
//   Overlay one selected 3D channel on the field-frame top-view mask.
//
//   Uses SNS_ProjectChannel3D_SurfaceXYZ so it is consistent with the actual
//   3D Cu(111)-faceted trajectory geometry.
//==============================================================================
Function SNS_PlotChannel3D_OnTopFieldFrame(Nmask, r0x, r0y, thetaList, phiList, idx, stepFac, H_nm, maxPath_nm, Bangle_deg, winXY, rayBaseName, r, g, b)
    Wave Nmask
    Variable r0x, r0y
    Wave thetaList, phiList
    Variable idx
    Variable stepFac, H_nm, maxPath_nm, Bangle_deg
    String winXY, rayBaseName
    Variable r, g, b

    if (idx < 0 || idx >= numpnts(thetaList))
        return -1
    endif

    Variable theta = thetaList[idx]
    Variable phi   = phiList[idx]

    String rayXName = rayBaseName + "_X_" + num2str(idx)
    String rayYName = rayBaseName + "_Y_" + num2str(idx)
    String rayZName = rayBaseName + "_Z_" + num2str(idx)

    Variable err = SNS_ProjectChannel3D_SurfaceXYZ(Nmask, r0x, r0y, theta, phi, \
        stepFac, H_nm, maxPath_nm, rayXName, rayYName, rayZName)

    if (err != 0)
        return err
    endif

    Wave rayX = $rayXName
    Wave rayY = $rayYName

    Variable n = numpnts(rayX)

    String rayXBName = rayBaseName + "_XB_" + num2str(idx)
    String rayYPName = rayBaseName + "_YP_" + num2str(idx)

    Make/O/D/N=(n) $rayXBName, $rayYPName

    Wave rayXB = $rayXBName
    Wave rayYP = $rayYPName

    Variable phiB = Bangle_deg*pi/180
    Variable bx = cos(phiB)
    Variable by = sin(phiB)
    Variable epx = -sin(phiB)
    Variable epy =  cos(phiB)

    rayXB = (rayX-r0x)*bx  + (rayY-r0y)*by
    rayYP = (rayX-r0x)*epx + (rayY-r0y)*epy

    AppendToGraph/W=$winXY rayYP vs rayXB
    ModifyGraph/W=$winXY mode($NameOfWave(rayYP))=0
    ModifyGraph/W=$winXY lsize($NameOfWave(rayYP))=2
    ModifyGraph/W=$winXY rgb($NameOfWave(rayYP))=(r,g,b)

    return 0
End

//==============================================================================
// SNS_AddTopMaskFieldFrameOutline
//
// Purpose:
//   Create/display an outline-only top-view mask in field-frame coordinates.
//
//   Coordinate convention:
//      x = coordinate parallel to in-plane B
//      y = coordinate perpendicular to in-plane B
//
//   The outline is extracted from the binary top-view mask after rotating into
//   the field frame. This is display-only.
//
// Inputs:
//   Nmask      : top-view binary mask.
//   r0x, r0y   : STS position [nm], used as field-frame origin.
//   Bangle_deg : in-plane field angle [deg].
//   winName    : graph window name.
//   baseName   : base name for generated outline waves.
//   doNew      : 1 create new graph, 0 append to existing graph.
//
// Output waves:
//   <baseName>_outline_x
//   <baseName>_outline_y
//==============================================================================
Function SNS_AddTopMaskFieldFrameOutline(Nmask, r0x, r0y, Bangle_deg, winName, baseName, [doNew])
    Wave Nmask
    Variable r0x, r0y, Bangle_deg
    String winName, baseName
    Variable doNew

    if (ParamIsDefault(doNew))
        doNew = 0
    endif

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)

    Variable x0 = DimOffset(Nmask, 0)
    Variable dx = DimDelta(Nmask, 0)
    Variable y0 = DimOffset(Nmask, 1)
    Variable dy = DimDelta(Nmask, 1)

    Variable phiB = Bangle_deg*pi/180
    Variable bx = cos(phiB)
    Variable by = sin(phiB)
    Variable epx = -sin(phiB)
    Variable epy =  cos(phiB)

    String xName = baseName + "_outline_x"
    String yName = baseName + "_outline_y"

    Make/O/D/N=0 $xName, $yName
    Wave xOut = $xName
    Wave yOut = $yName

    Variable i, j
    Variable isEdge, n
    Variable x, y, xB, yP

    for (i = 1; i < nx-1; i += 1)
        x = x0 + i*dx

        for (j = 1; j < ny-1; j += 1)
            if (Nmask[i][j] > 0.5)

                isEdge = (Nmask[i+1][j] <= 0.5) || (Nmask[i-1][j] <= 0.5) || \
                         (Nmask[i][j+1] <= 0.5) || (Nmask[i][j-1] <= 0.5)

                if (isEdge)
                    y = y0 + j*dy

                    xB = (x-r0x)*bx  + (y-r0y)*by
                    yP = (x-r0x)*epx + (y-r0y)*epy

                    n = numpnts(xOut)
                    Redimension/N=(n+1) xOut, yOut

                    xOut[n] = xB
                    yOut[n] = yP
                endif
            endif
        endfor
    endfor

    if (numpnts(xOut) <= 0)
        return -1
    endif

    if (doNew)
        Display/K=1 yOut vs xOut
        DoWindow/C $winName
    else
        AppendToGraph/W=$winName yOut vs xOut
    endif

    ModifyGraph/W=$winName mode($NameOfWave(yOut))=3
    ModifyGraph/W=$winName marker($NameOfWave(yOut))=19
    ModifyGraph/W=$winName msize($NameOfWave(yOut))=1
    ModifyGraph/W=$winName rgb($NameOfWave(yOut))=(0,0,0)

    return 0
End

//==============================================================================
// SNS_Cu111BoundaryFrac
//
// Purpose:
//   Given an inside point and an outside point, find the last inside fraction
//   along the segment using bisection.
//
// Return:
//   f in [0,1], where
//      r(f) = rIn + f*(rOut-rIn)
//   is still inside the Cu(111) prism.
//==============================================================================
Function SNS_Cu111BoundaryFrac(Nmask, xIn, yIn, zIn, xOut, yOut, zOut, H_nm)
    Wave Nmask
    Variable xIn, yIn, zIn, xOut, yOut, zOut, H_nm

    Variable lo = 0
    Variable hi = 1
    Variable mid, xm, ym, zm
    Variable k

    for (k = 0; k < 24; k += 1)
        mid = 0.5*(lo + hi)

        xm = xIn + mid*(xOut - xIn)
        ym = yIn + mid*(yOut - yIn)
        zm = zIn + mid*(zOut - zIn)

        if (SNS_IsInsideCu111Prism(Nmask, xm, ym, zm, H_nm))
            lo = mid
        else
            hi = mid
        endif
    endfor

    return lo
End

//==============================================================================
// SNS_TraceOneChannel3D_Cu111
//
// Purpose:
//   Trace one 3D SNS channel in the Cu(111)-faceted prism geometry.
//
//   This is the single shared ray engine used by:
//      SNS_BuildChannelsFromMask3D
//      SNS_ProjectChannel3D_Surface
//      SNS_ProjectChannel3D_SurfaceXYZ
//
// Interface transparency:
//   The barrier strength is read from SNS_Settings via SNS_Params.
//   This tracer no longer accepts a Zbarrier / BTK_barrier input.
//
//   The transmission is evaluated with the 3D bottom-interface incidence angle,
//   using cosInc = |vz| at the Cu-Nb interface:
//
//      T_3D = 4*cosInc^2 / (4*cosInc^2 + Z^2)
//
//   This is intentionally not identical to the 2D edge-interface helper.
//
// Geometry:
//   Nmask is the top mask at z = H_nm.
//   The bottom footprint expands laterally according to
//      SNS_Cu111PrismSideExpand_nm(z,H).
//
// Magnetic-area convention:
//   The legacy SNS setting params.lambdaL is currently the effective 2D height:
//
//      h_eff_2D = hSC_eff + H_nm
//
//   Therefore the 3D tracer uses:
//
//      hSC_eff_nm = max(0, params.lambdaL*1e9 - H_nm)
//
//   and computes the channel-dependent signed area
//
//      Aeff_nm2 = hSC_eff_nm * W_signed_nm + A_N_nm2
//
//   where
//
//      A_N_nm2 = integral_{S1->STM->S2} z d y_perp
//
//   is evaluated along the actual reflected 3D trajectory.
//
// Output:
//   Creates/overwrites resName wave:
//
//      res[0]  = success flag, 1 success, 0 fail
//      res[1]  = L_N_nm
//      res[2]  = W_eff_equiv_nm
//                solver-equivalent width such that
//                params.lambdaL * W_eff_equiv = |Aeff_3D|
//      res[3]  = Hit1x_nm
//      res[4]  = Hit1y_nm
//      res[5]  = Hit2x_nm
//      res[6]  = Hit2y_nm
//      res[7]  = T_eff
//      res[8]  = betaExtra, rad, NaN if no Q
//      res[9]  = w_raw
//      res[10] = Aeff_3D_nm2, signed
//      res[11] = W_geom_nm = |(S2-S1) dot e_perp|
//
//   If savePath=1, creates/overwrites rayXName, rayYName, rayZName containing:
//      leg 1, NaN separator, leg 2
//
// Return:
//   0 on success, <0 on failure.
//==============================================================================
Function SNS_TraceOneChannel3D_Cu111(Nmask, r0x, r0y, theta, phi, phiB, stepFac, H_nm, maxPath_nm, resName, savePath, rayXName, rayYName, rayZName, [QxPhase, QyPhase, qNstep, qCoreHandling, qMinUsedFrac, xV_nm, yV_nm, rCore_nm, PhaseReCore, PhaseImCore])
    Wave Nmask
    Variable r0x, r0y
    Variable theta, phi, phiB
    Variable stepFac, H_nm, maxPath_nm
    String resName
    Variable savePath
    String rayXName, rayYName, rayZName

    Wave QxPhase, QyPhase
    Variable qNstep, qCoreHandling, qMinUsedFrac
    Variable xV_nm, yV_nm, rCore_nm
    Wave PhaseReCore, PhaseImCore

    Make/O/D/N=12 $resName
    Wave res = $resName
    res = NaN
    res[0] = 0

    Variable haveQ = 0
    if (!ParamIsDefault(QxPhase) && !ParamIsDefault(QyPhase))
        haveQ = 1
    elseif (!ParamIsDefault(QxPhase) || !ParamIsDefault(QyPhase))
        return -100
    endif

    if (ParamIsDefault(qNstep))
        qNstep = 1
    endif
    qNstep = max(1, round(qNstep))

    if (ParamIsDefault(qCoreHandling))
        qCoreHandling = 0
    endif
    qCoreHandling = round(qCoreHandling)

    if (ParamIsDefault(qMinUsedFrac))
        qMinUsedFrac = 0.95
    endif

    Variable nx = DimSize(Nmask, 0)
    Variable ny = DimSize(Nmask, 1)
    if (nx <= 1 || ny <= 1)
        return -1
    endif

    Variable dx = DimDelta(Nmask, 0)
    Variable dy = DimDelta(Nmask, 1)

    Variable baseStep = min(abs(dx), abs(dy))
    if (stepFac <= 0)
        stepFac = 0.5
    endif
    Variable dsBase = stepFac*baseStep

    Variable maxBounces = 50
    Variable maxSteps = ceil(maxPath_nm/dsBase) + maxBounces + 20
    Variable maxPts = 2*maxSteps + 4

    if (savePath)
        Make/O/D/N=(maxPts) $rayXName, $rayYName, $rayZName
        Wave rayX = $rayXName
        Wave rayY = $rayYName
        Wave rayZ = $rayZName
        rayX = NaN
        rayY = NaN
        rayZ = NaN
    endif

    Variable p = 0

    Variable sinT = sin(theta)
    Variable cosT = cos(theta)

    Variable vx0 = cos(phi)*sinT
    Variable vy0 = sin(phi)*sinT
    Variable vz0 = -cosT

    Variable epx = -sin(phiB)
    Variable epy =  cos(phiB)

    STRUCT SNS_Params params
    SNS_LoadParams(params)
    
    Variable Zbarrier = params.BTK_barrier
    Variable lambdaLegacy_nm = params.lambdaL * 1e9
    Variable hSC_eff_nm = max(0, lambdaLegacy_nm - H_nm)

    // -------------------------------------------------------------------------
    // LEG 1: STM point -> S plane
    // -------------------------------------------------------------------------
    Variable vx = vx0
    Variable vy = vy0
    Variable vz = vz0

    Variable x = r0x
    Variable y = r0y
    Variable z = H_nm

    if (!SNS_IsInsideCu111Prism(Nmask, x, y, z, H_nm))
        return -2
    endif

    if (savePath)
        rayX[p] = x
        rayY[p] = y
        rayZ[p] = z
        p += 1
    endif

    Variable pathLen = 0
    Variable nBounce = 0
    Variable haveS1 = 0
    Variable Lhalf1_nm = NaN
    Variable xS1 = NaN
    Variable yS1 = NaN
    Variable betaLeg1 = 0
    Variable betaSeg
    Variable areaLeg1_nm2 = 0

    do
        Variable x_prev = x
        Variable y_prev = y
        Variable z_prev = z

        Variable ds = min(dsBase, maxPath_nm - pathLen)
        if (ds <= 0)
            break
        endif

        Variable x_new = x + ds*vx
        Variable y_new = y + ds*vy
        Variable z_new = z + ds*vz

        // Hit S plane.
        if ((z_new <= 0) && (z_prev > 0))
            Variable tS1 = z_prev / (z_prev - z_new)

            xS1 = x_prev + tS1*(x_new - x_prev)
            yS1 = y_prev + tS1*(y_new - y_prev)

            if (!SNS_IsInsideCu111Prism(Nmask, xS1, yS1, 0, H_nm))
                break
            endif

            Variable yp_prev = epx*x_prev + epy*y_prev
            Variable yp_hit  = epx*xS1   + epy*yS1
            areaLeg1_nm2 += 0.5*(z_prev + 0) * (yp_hit - yp_prev)

            if (haveQ)
                if (qCoreHandling == 0)
                    betaLeg1 += SNS_IntegrateQSegmentNearest2D(QxPhase, QyPhase, x_prev, y_prev, xS1, yS1, qNstep)
                else
                    betaSeg = SNS_IntegrateQSegmentNearest2D_CorePhase(QxPhase, QyPhase, x_prev, y_prev, xS1, yS1, qNstep, qCoreHandling, xV_nm, yV_nm, rCore_nm, qMinUsedFrac, PhaseReCore=PhaseReCore, PhaseImCore=PhaseImCore)
                    if (numtype(betaSeg))
                        break
                    endif
                    betaLeg1 += betaSeg
                endif
            endif

            Lhalf1_nm = pathLen + tS1*ds

            if (savePath && p < maxPts)
                rayX[p] = xS1
                rayY[p] = yS1
                rayZ[p] = 0
                p += 1
            endif

            haveS1 = 1
            break
        endif

        if (nBounce > maxBounces)
            break
        endif

        if (z_new > 0)
            Variable insidePrev = SNS_IsInsideCu111Prism(Nmask, x_prev, y_prev, z_prev, H_nm)
            Variable insideNow  = SNS_IsInsideCu111Prism(Nmask, x_new, y_new, z_new, H_nm)

            if (insidePrev && !insideNow)

                Variable fHit = SNS_Cu111BoundaryFrac(Nmask, x_prev, y_prev, z_prev, x_new, y_new, z_new, H_nm)

                Variable xHit = x_prev + fHit*(x_new - x_prev)
                Variable yHit = y_prev + fHit*(y_new - y_prev)
                Variable zHit = z_prev + fHit*(z_new - z_prev)

                yp_prev = epx*x_prev + epy*y_prev
                yp_hit  = epx*xHit   + epy*yHit
                areaLeg1_nm2 += 0.5*(z_prev + zHit) * (yp_hit - yp_prev)

                if (haveQ)
                    if (qCoreHandling == 0)
                        betaLeg1 += SNS_IntegrateQSegmentNearest2D(QxPhase, QyPhase, x_prev, y_prev, xHit, yHit, qNstep)
                    else
                        betaSeg = SNS_IntegrateQSegmentNearest2D_CorePhase(QxPhase, QyPhase, x_prev, y_prev, xHit, yHit, qNstep, qCoreHandling, xV_nm, yV_nm, rCore_nm, qMinUsedFrac, PhaseReCore=PhaseReCore, PhaseImCore=PhaseImCore)
                        if (numtype(betaSeg))
                            break
                        endif
                        betaLeg1 += betaSeg
                    endif
                endif

                pathLen += fHit*ds

                if (savePath && p < maxPts)
                    rayX[p] = xHit
                    rayY[p] = yHit
                    rayZ[p] = zHit
                    p += 1
                endif

                Variable nxN = SNS_Cu111PrismNormalComponent(Nmask, xHit, yHit, zHit, H_nm, 0)
                Variable nyN = SNS_Cu111PrismNormalComponent(Nmask, xHit, yHit, zHit, H_nm, 1)
                Variable nzN = SNS_Cu111PrismNormalComponent(Nmask, xHit, yHit, zHit, H_nm, 2)

                Variable dot = vx*nxN + vy*nyN + vz*nzN

                vx -= 2*dot*nxN
                vy -= 2*dot*nyN
                vz -= 2*dot*nzN

                // Restart just inside the wall.
                Variable eps = 1e-6*baseStep
                x = xHit - eps*nxN
                y = yHit - eps*nyN
                z = zHit - eps*nzN

                if (savePath && p < maxPts)
                    rayX[p] = x
                    rayY[p] = y
                    rayZ[p] = z
                    p += 1
                endif

                nBounce += 1
                continue

            else
                yp_prev = epx*x_prev + epy*y_prev
                Variable yp_new = epx*x_new + epy*y_new
                areaLeg1_nm2 += 0.5*(z_prev + z_new) * (yp_new - yp_prev)

                if (haveQ)
                    if (qCoreHandling == 0)
                        betaLeg1 += SNS_IntegrateQSegmentNearest2D(QxPhase, QyPhase, x_prev, y_prev, x_new, y_new, qNstep)
                    else
                        betaSeg = SNS_IntegrateQSegmentNearest2D_CorePhase(QxPhase, QyPhase, x_prev, y_prev, x_new, y_new, qNstep, qCoreHandling, xV_nm, yV_nm, rCore_nm, qMinUsedFrac, PhaseReCore=PhaseReCore, PhaseImCore=PhaseImCore)
                        if (numtype(betaSeg))
                            break
                        endif
                        betaLeg1 += betaSeg
                    endif
                endif

                x = x_new
                y = y_new
                z = z_new
                pathLen += ds

                if (savePath && p < maxPts)
                    rayX[p] = x
                    rayY[p] = y
                    rayZ[p] = z
                    p += 1
                endif
            endif
        endif

    while (1)

    if (!haveS1)
        return -3
    endif

    Variable cosInc1 = abs(vz)
    Variable cos2_1 = cosInc1*cosInc1
    Variable T_S1
    if (Zbarrier <= 0)
        T_S1 = 1
    else
        Variable denom1 = 4*cos2_1 + Zbarrier*Zbarrier
        T_S1 = (denom1 > 0) ? 4*cos2_1/denom1 : 0
    endif

    if (savePath && p < maxPts)
        rayX[p] = NaN
        rayY[p] = NaN
        rayZ[p] = NaN
        p += 1
    endif

    // -------------------------------------------------------------------------
    // LEG 2: opposite initial in-plane direction from STM point -> S plane
    // -------------------------------------------------------------------------
    Variable vx2 = -vx0
    Variable vy2 = -vy0
    Variable vz2 =  vz0

    Variable x2 = r0x
    Variable y2 = r0y
    Variable z2 = H_nm

    if (!SNS_IsInsideCu111Prism(Nmask, x2, y2, z2, H_nm))
        return -4
    endif

    if (savePath && p < maxPts)
        rayX[p] = x2
        rayY[p] = y2
        rayZ[p] = z2
        p += 1
    endif

    Variable pathLen2 = 0
    Variable nBounce2 = 0
    Variable haveS2 = 0
    Variable Lhalf2_nm = NaN
    Variable xS2 = NaN
    Variable yS2 = NaN
    Variable betaLeg2 = 0
    Variable areaLeg2_nm2 = 0

    do
        Variable x2_prev = x2
        Variable y2_prev = y2
        Variable z2_prev = z2

        Variable ds2 = min(dsBase, maxPath_nm - pathLen2)
        if (ds2 <= 0)
            break
        endif

        Variable x2_new = x2 + ds2*vx2
        Variable y2_new = y2 + ds2*vy2
        Variable z2_new = z2 + ds2*vz2

        // Hit S plane.
        if ((z2_new <= 0) && (z2_prev > 0))
            Variable tS2 = z2_prev / (z2_prev - z2_new)

            xS2 = x2_prev + tS2*(x2_new - x2_prev)
            yS2 = y2_prev + tS2*(y2_new - y2_prev)

            if (!SNS_IsInsideCu111Prism(Nmask, xS2, yS2, 0, H_nm))
                break
            endif

            Variable yp2_prev = epx*x2_prev + epy*y2_prev
            Variable yp2_hit  = epx*xS2     + epy*yS2
            areaLeg2_nm2 += 0.5*(z2_prev + 0) * (yp2_hit - yp2_prev)

            if (haveQ)
                if (qCoreHandling == 0)
                    betaLeg2 += SNS_IntegrateQSegmentNearest2D(QxPhase, QyPhase, x2_prev, y2_prev, xS2, yS2, qNstep)
                else
                    betaSeg = SNS_IntegrateQSegmentNearest2D_CorePhase(QxPhase, QyPhase, x2_prev, y2_prev, xS2, yS2, qNstep, qCoreHandling, xV_nm, yV_nm, rCore_nm, qMinUsedFrac, PhaseReCore=PhaseReCore, PhaseImCore=PhaseImCore)
                    if (numtype(betaSeg))
                        break
                    endif
                    betaLeg2 += betaSeg
                endif
            endif

            Lhalf2_nm = pathLen2 + tS2*ds2

            if (savePath && p < maxPts)
                rayX[p] = xS2
                rayY[p] = yS2
                rayZ[p] = 0
                p += 1
            endif

            haveS2 = 1
            break
        endif

        if (nBounce2 > maxBounces)
            break
        endif

        if (z2_new > 0)
            Variable insidePrev2 = SNS_IsInsideCu111Prism(Nmask, x2_prev, y2_prev, z2_prev, H_nm)
            Variable insideNow2  = SNS_IsInsideCu111Prism(Nmask, x2_new, y2_new, z2_new, H_nm)

            if (insidePrev2 && !insideNow2)

                Variable fHit2 = SNS_Cu111BoundaryFrac(Nmask, x2_prev, y2_prev, z2_prev, x2_new, y2_new, z2_new, H_nm)

                Variable xHit2 = x2_prev + fHit2*(x2_new - x2_prev)
                Variable yHit2 = y2_prev + fHit2*(y2_new - y2_prev)
                Variable zHit2 = z2_prev + fHit2*(z2_new - z2_prev)

                yp2_prev = epx*x2_prev + epy*y2_prev
                yp2_hit  = epx*xHit2   + epy*yHit2
                areaLeg2_nm2 += 0.5*(z2_prev + zHit2) * (yp2_hit - yp2_prev)

                if (haveQ)
                    if (qCoreHandling == 0)
                        betaLeg2 += SNS_IntegrateQSegmentNearest2D(QxPhase, QyPhase, x2_prev, y2_prev, xHit2, yHit2, qNstep)
                    else
                        betaSeg = SNS_IntegrateQSegmentNearest2D_CorePhase(QxPhase, QyPhase, x2_prev, y2_prev, xHit2, yHit2, qNstep, qCoreHandling, xV_nm, yV_nm, rCore_nm, qMinUsedFrac, PhaseReCore=PhaseReCore, PhaseImCore=PhaseImCore)
                        if (numtype(betaSeg))
                            break
                        endif
                        betaLeg2 += betaSeg
                    endif
                endif

                pathLen2 += fHit2*ds2

                if (savePath && p < maxPts)
                    rayX[p] = xHit2
                    rayY[p] = yHit2
                    rayZ[p] = zHit2
                    p += 1
                endif

                Variable nxN2 = SNS_Cu111PrismNormalComponent(Nmask, xHit2, yHit2, zHit2, H_nm, 0)
                Variable nyN2 = SNS_Cu111PrismNormalComponent(Nmask, xHit2, yHit2, zHit2, H_nm, 1)
                Variable nzN2 = SNS_Cu111PrismNormalComponent(Nmask, xHit2, yHit2, zHit2, H_nm, 2)

                Variable dot2 = vx2*nxN2 + vy2*nyN2 + vz2*nzN2

                vx2 -= 2*dot2*nxN2
                vy2 -= 2*dot2*nyN2
                vz2 -= 2*dot2*nzN2

                eps = 1e-6*baseStep
                x2 = xHit2 - eps*nxN2
                y2 = yHit2 - eps*nyN2
                z2 = zHit2 - eps*nzN2

                if (savePath && p < maxPts)
                    rayX[p] = x2
                    rayY[p] = y2
                    rayZ[p] = z2
                    p += 1
                endif

                nBounce2 += 1
                continue

            else
                yp2_prev = epx*x2_prev + epy*y2_prev
                Variable yp2_new = epx*x2_new + epy*y2_new
                areaLeg2_nm2 += 0.5*(z2_prev + z2_new) * (yp2_new - yp2_prev)

                if (haveQ)
                    if (qCoreHandling == 0)
                        betaLeg2 += SNS_IntegrateQSegmentNearest2D(QxPhase, QyPhase, x2_prev, y2_prev, x2_new, y2_new, qNstep)
                    else
                        betaSeg = SNS_IntegrateQSegmentNearest2D_CorePhase(QxPhase, QyPhase, x2_prev, y2_prev, x2_new, y2_new, qNstep, qCoreHandling, xV_nm, yV_nm, rCore_nm, qMinUsedFrac, PhaseReCore=PhaseReCore, PhaseImCore=PhaseImCore)
                        if (numtype(betaSeg))
                            break
                        endif
                        betaLeg2 += betaSeg
                    endif
                endif

                x2 = x2_new
                y2 = y2_new
                z2 = z2_new
                pathLen2 += ds2

                if (savePath && p < maxPts)
                    rayX[p] = x2
                    rayY[p] = y2
                    rayZ[p] = z2
                    p += 1
                endif
            endif
        endif

    while (1)

    if (!haveS2)
        return -5
    endif

    Variable cosInc2 = abs(vz2)
    Variable cos2_2 = cosInc2*cosInc2
    Variable T_S2
    if (Zbarrier <= 0)
        T_S2 = 1
    else
        Variable denom2 = 4*cos2_2 + Zbarrier*Zbarrier
        T_S2 = (denom2 > 0) ? 4*cos2_2/denom2 : 0
    endif

    Variable L_N_nm = Lhalf1_nm + Lhalf2_nm
    if (!(L_N_nm > 0))
        return -6
    endif

    Variable W_signed_nm = (xS2 - xS1)*epx + (yS2 - yS1)*epy
    Variable W_geom_nm = abs(W_signed_nm)

    // Area over physical path S1 -> STM -> S2.
    // Both legs were traced STM -> S, so use leg2 - leg1.
    Variable A_N_nm2 = areaLeg2_nm2 - areaLeg1_nm2
    Variable Aeff_nm2 = hSC_eff_nm * W_signed_nm + A_N_nm2

    Variable W_eff_equiv_nm
    if (lambdaLegacy_nm > 0)
        W_eff_equiv_nm = abs(Aeff_nm2 / lambdaLegacy_nm)
    else
        W_eff_equiv_nm = W_geom_nm
    endif

    Variable Teff = min(T_S1, T_S2)

    res[0]  = 1
    res[1]  = L_N_nm
    res[2]  = W_eff_equiv_nm
    res[3]  = xS1
    res[4]  = yS1
    res[5]  = xS2
    res[6]  = yS2
    res[7]  = Teff
    res[8]  = haveQ ? (betaLeg2 - betaLeg1) : NaN
    res[9]  = sinT
    res[10] = Aeff_nm2
    res[11] = W_geom_nm

    if (savePath)
        Redimension/N=(p) rayX, rayY, rayZ
    endif

    return 0
End


//==============================================================================
// SNS_CropMaskToActiveRect
//
// Purpose
//   Crop a sparse synthetic mask to the smallest rectangular index region containing
//   all active pixels plus a physical margin. This is intended for synthetic
//   notebooks where the structure is first defined in a large 2D field of view,
//   but Python sidecars should only see the compact active region plus enough
//   surrounding S/outside area for boundary detection.
//
// Created output
//   outName in the current data folder, with x/y scales copied from the source
//   coordinate system and shifted to the cropped origin.
//
// Notes
//   The crop is rectangular. If the requested rectangle exceeds the source
//   wave bounds, it is clipped to fit the available FOV.
//==============================================================================

Function SNS_CropMaskToActiveRect(maskW, outName, margin_nm, [threshold])
    Wave maskW
    String outName
    Variable margin_nm
    Variable threshold

    if (ParamIsDefault(threshold))
        threshold = 0.5
    endif

    Variable nx = DimSize(maskW, 0)
    Variable ny = DimSize(maskW, 1)
    if (nx <= 0 || ny <= 0)
        Abort "SNS_CropMaskToActiveRect: empty mask."
    endif

    Variable dx = DimDelta(maskW, 0)
    Variable dy = DimDelta(maskW, 1)
    Variable x0 = DimOffset(maskW, 0)
    Variable y0 = DimOffset(maskW, 1)
    Variable marginPts = max(0, ceil(margin_nm / max(1e-12, min(abs(dx), abs(dy)))))

    Variable i, j
    Variable iMin = nx
    Variable iMax = -1
    Variable jMin = ny
    Variable jMax = -1

    for (i = 0; i < nx; i += 1)
        for (j = 0; j < ny; j += 1)
            if (maskW[i][j] > threshold)
                iMin = min(iMin, i)
                iMax = max(iMax, i)
                jMin = min(jMin, j)
                jMax = max(jMax, j)
            endif
        endfor
    endfor

    if (iMax < iMin || jMax < jMin)
        Abort "SNS_CropMaskToActiveRect: no active pixels above threshold."
    endif

    iMin = max(0, iMin - marginPts)
    iMax = min(nx - 1, iMax + marginPts)
    jMin = max(0, jMin - marginPts)
    jMax = min(ny - 1, jMax + marginPts)

    Variable nCropX = iMax - iMin + 1
    Variable nCropY = jMax - jMin + 1

    Make/O/D/N=(nCropX, nCropY) $outName
    Wave outW = $outName
    outW = maskW[p + iMin][q + jMin]
    SetScale/P x, x0 + iMin * dx, dx, WaveUnits(maskW, 0), outW
    SetScale/P y, y0 + jMin * dy, dy, WaveUnits(maskW, 1), outW

    Variable/G v_CropMask_i0 = iMin
    Variable/G v_CropMask_j0 = jMin
    Variable/G v_CropMask_nx = nCropX
    Variable/G v_CropMask_ny = nCropY
    Variable/G v_CropMask_margin_nm = margin_nm
    Variable/G v_CropMask_x0_nm = x0 + iMin * dx
    Variable/G v_CropMask_y0_nm = y0 + jMin * dy
End


// Backward-compatible alias for notebooks from the brief square-crop iteration.
Function SNS_CropMaskToActiveSquare(maskW, outName, margin_nm, [threshold])
    Wave maskW
    String outName
    Variable margin_nm
    Variable threshold

    if (ParamIsDefault(threshold))
        SNS_CropMaskToActiveRect(maskW, outName, margin_nm)
    else
        SNS_CropMaskToActiveRect(maskW, outName, margin_nm, threshold=threshold)
    endif
End


//==============================================================================
// Export a 2D island mask for a Kwant / finite-difference BdG calculation.
//
// Output files:
//   <stem>_sites.tsv
//      id   ix   iy   x_nm   y_nm
//
//   <stem>_hoppings.tsv
//      id1  id2
//
//   <stem>_meta.tsv
//      key  value
//
// The export uses block downsampling. A block is kept if its active-pixel
// fraction exceeds fracMin. This is usually safer than taking only every nth
// pixel.
//==============================================================================

Function SNS_ExportMaskForKwant(maskW, pathName, stem, [threshold, block, fracMin, use8Connect])
    Wave maskW
    String pathName
    String stem
    Variable threshold, block, fracMin, use8Connect

    if (ParamIsDefault(threshold))
        threshold = 0.5
    endif
    if (ParamIsDefault(block))
        block = 1
    endif
    if (ParamIsDefault(fracMin))
        fracMin = 0.5
    endif
    if (ParamIsDefault(use8Connect))
        use8Connect = 0
    endif

    block = max(1, round(block))

    Variable nx = DimSize(maskW, 0)
    Variable ny = DimSize(maskW, 1)

    Variable dx = DimDelta(maskW, 0)
    Variable dy = DimDelta(maskW, 1)
    Variable x0 = DimOffset(maskW, 0)
    Variable y0 = DimOffset(maskW, 1)

    Variable nxB = ceil(nx / block)
    Variable nyB = ceil(ny / block)

    String siteMapName = UniqueName("kw_siteMap", 1, 0)
    Make/O/N=(nxB, nyB)/I $siteMapName
    Wave siteMap = $siteMapName
    siteMap = -1

    Variable iB, jB, id
    Variable i0, j0, activeFrac
    id = 0

    // Assign site IDs.
    for (iB = 0; iB < nxB; iB += 1)
        for (jB = 0; jB < nyB; jB += 1)
            i0 = iB * block
            j0 = jB * block

            activeFrac = SNS_BlockActiveFraction(maskW, i0, j0, block, threshold)

            if (activeFrac >= fracMin)
                siteMap[iB][jB] = id
                id += 1
            endif
        endfor
    endfor

    Variable nSites = id

    // -------------------------------------------------------------------------
    // Write sites.
    // -------------------------------------------------------------------------
    Variable ref
    String fSites = stem + "_sites.tsv"
    Open/P=$pathName ref as fSites

    fprintf ref, "id\tix\tiy\tx_nm\ty_nm\n"

    Variable x_nm, y_nm
    for (iB = 0; iB < nxB; iB += 1)
        for (jB = 0; jB < nyB; jB += 1)
            if (siteMap[iB][jB] >= 0)
                // Use block centre in the original scaled coordinate system.
                x_nm = x0 + (iB * block + 0.5 * (block - 1)) * dx
                y_nm = y0 + (jB * block + 0.5 * (block - 1)) * dy

                fprintf ref, "%d\t%d\t%d\t%.10g\t%.10g\n", siteMap[iB][jB], iB, jB, x_nm, y_nm
            endif
        endfor
    endfor

    Close ref

    // -------------------------------------------------------------------------
    // Write hoppings. Default is 4-neighbour connectivity.
    // Optional use8Connect adds diagonals, but I would keep it off initially.
    // -------------------------------------------------------------------------
    String fHop = stem + "_hoppings.tsv"
    Open/P=$pathName ref as fHop

    fprintf ref, "id1\tid2\n"

    for (iB = 0; iB < nxB; iB += 1)
        for (jB = 0; jB < nyB; jB += 1)
            if (siteMap[iB][jB] < 0)
                continue
            endif

            SNS_WriteHopIfValid(ref, siteMap, iB, jB, iB + 1, jB)
            SNS_WriteHopIfValid(ref, siteMap, iB, jB, iB, jB + 1)

            if (use8Connect)
                SNS_WriteHopIfValid(ref, siteMap, iB, jB, iB + 1, jB + 1)
                SNS_WriteHopIfValid(ref, siteMap, iB, jB, iB + 1, jB - 1)
            endif
        endfor
    endfor

    Close ref

    // -------------------------------------------------------------------------
    // Write metadata.
    // -------------------------------------------------------------------------
    String fMeta = stem + "_meta.tsv"
    Open/P=$pathName ref as fMeta

    fprintf ref, "key\tvalue\n"
    fprintf ref, "source_wave\t%s\n", NameOfWave(maskW)
    fprintf ref, "nx_original\t%d\n", nx
    fprintf ref, "ny_original\t%d\n", ny
    fprintf ref, "nx_export\t%d\n", nxB
    fprintf ref, "ny_export\t%d\n", nyB
    fprintf ref, "n_sites\t%d\n", nSites
    fprintf ref, "dx_original_nm\t%.10g\n", dx
    fprintf ref, "dy_original_nm\t%.10g\n", dy
    fprintf ref, "block\t%d\n", block
    fprintf ref, "dx_export_nm\t%.10g\n", dx * block
    fprintf ref, "dy_export_nm\t%.10g\n", dy * block
    fprintf ref, "x0_nm\t%.10g\n", x0
    fprintf ref, "y0_nm\t%.10g\n", y0
    fprintf ref, "threshold\t%.10g\n", threshold
    fprintf ref, "fracMin\t%.10g\n", fracMin
    fprintf ref, "use8Connect\t%d\n", use8Connect
    fprintf ref, "m_eff_over_me\t%.10g\n", 0.45
    fprintf ref, "EF_eV\t%.10g\n", 0.429
    fprintf ref, "Delta_eff_meV\t%.10g\n", 0.92
    fprintf ref, "h_eff_nm\t%.10g\n", 24
    fprintf ref, "lambda_F_nm\t%.10g\n", 2.8

    Close ref

    Print "SNS_ExportMaskForKwant: exported ", nSites, " sites to stem = ", stem
    Print "Files: ", fSites, ", ", fHop, ", ", fMeta

    KillWaves/Z siteMap
End


//------------------------------------------------------------------------------
// Fraction of active pixels inside a block.
//------------------------------------------------------------------------------

Function SNS_BlockActiveFraction(maskW, i0, j0, block, threshold)
    Wave maskW
    Variable i0, j0, block, threshold

    Variable nx = DimSize(maskW, 0)
    Variable ny = DimSize(maskW, 1)

    Variable ii, jj
    Variable count = 0
    Variable active = 0

    for (ii = i0; ii < min(i0 + block, nx); ii += 1)
        for (jj = j0; jj < min(j0 + block, ny); jj += 1)
            count += 1
            if (maskW[ii][jj] > threshold)
                active += 1
            endif
        endfor
    endfor

    if (count <= 0)
        return 0
    endif

    return active / count
End


//------------------------------------------------------------------------------
// Write one hopping if both sites exist.
//------------------------------------------------------------------------------

Function SNS_WriteHopIfValid(ref, siteMap, i1, j1, i2, j2)
    Variable ref
    Wave siteMap
    Variable i1, j1, i2, j2

    Variable nx = DimSize(siteMap, 0)
    Variable ny = DimSize(siteMap, 1)

    if (i2 < 0 || i2 >= nx || j2 < 0 || j2 >= ny)
        return 0
    endif

    if (siteMap[i1][j1] >= 0 && siteMap[i2][j2] >= 0)
        fprintf ref, "%d\t%d\n", siteMap[i1][j1], siteMap[i2][j2]
    endif

    return 1
End
