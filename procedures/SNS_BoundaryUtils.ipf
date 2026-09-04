#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma DefaultTab={3,20,4}

//==============================================================================
// SNS_BoundaryUtils.ipf
//
// Generic mask/boundary geometry utilities used by SNS ray tracing.
//
// Responsibilities:
//   - boundary-mask generation
//   - inside/outside checks
//   - chord / bounding-box helpers
//   - mask geometry helpers independent of ABS solving
//
// Notes:
//   Function names are intentionally left unchanged during this file rename.
//   This file should not solve ABS spectra, build DOS, or draw final figures.
//==============================================================================

//==============================================================================
// MakeBoundaryMaskNaNNeighbors
//
// Purpose:
//   Construct output waves used by the SNS workflow.
//
// Inputs:
//   src : input
//   outName : input
//   diag : input
//
// Outputs:
//   return : numeric
//   outName : output wave(s)
//
//==============================================================================
Function MakeBoundaryMaskNaNNeighbors(src, outName, diag)
    Wave src
    String outName
    Variable diag

    Variable nRows = DimSize(src, 0)
    Variable nCols = DimSize(src, 1)

    Make /O /U /B /N=(nRows, nCols) $outName
    Wave boundaryMask = $outName
    boundaryMask = 0

    Variable r, c, dr, dc, rr, cc
    Variable hasNaNNeighbor, v, isInterior

    for (r=0; r<nRows; r+=1)
        for (c=0; c<nCols; c+=1)
            v = src[r][c]
            if (numtype(v) == 2)          // this pixel is NaN => not part of island
                continue
            endif

            hasNaNNeighbor = 0

            // Check neighbors; treat out-of-bounds as NaN so border pixels become boundary
            // If diag==0, only check N,S,E,W. If diag==1, check all 8 neighbors.
            for (dr=-1; dr<=1 && !hasNaNNeighbor; dr+=1)
                for (dc=-1; dc<=1 && !hasNaNNeighbor; dc+=1)
                    if (dr==0 && dc==0)
                        continue
                    endif
                    if (!diag && abs(dr)+abs(dc) == 2)
                        continue        // skip diagonals in 4-neighborhood mode
                    endif

                    rr = r + dr
                    cc = c + dc

                    if (rr < 0 || rr >= nRows || cc < 0 || cc >= nCols)
                        hasNaNNeighbor = 1
                    else
                        if (numtype(src[rr][cc]) == 2)
                            hasNaNNeighbor = 1
                        endif
                    endif
                endfor
            endfor

            if (hasNaNNeighbor)
                boundaryMask[r][c] = 1
            endif
        endfor
    endfor

    // Optional: give the mask the same scaling as src
    SetScale /P x, DimOffset(src,0), DimDelta(src,0), "", boundaryMask
    SetScale /P y, DimOffset(src,1), DimDelta(src,1), "", boundaryMask

    // Convenience prints
    // Print "Boundary mask created:", NameOfWave(boundaryMask), " (1=boundary, 0=elsewhere)"
End


//==============================================================================
// FindDiameterFromMask
//
// Purpose:
//   Locate a feature or parameter from input data.
//
// Inputs:
//   boundaryMask : input
//   useScaling : input
//   outBase : input
//   keepD2 : input
//
// Outputs:
//   return : numeric
//   outBase : output wave(s)
//
//==============================================================================
Function FindDiameterFromMask(boundaryMask, useScaling, outBase, keepD2)
    Wave boundaryMask
    Variable useScaling
    String outBase
    Variable keepD2

    if (WaveDims(boundaryMask) != 2)
        Abort "FindDiameter_FromMask: input must be a 2D wave."
    endif

    Variable nRows = DimSize(boundaryMask, 0)
    Variable nCols = DimSize(boundaryMask, 1)

    // --- 1) Gather coordinates of boundary pixels (value == 1) ---
    Variable r, c, N=0
    for (r=0; r<nRows; r+=1)
        for (c=0; c<nCols; c+=1)
            if (boundaryMask[r][c] == 1)
                N += 1
            endif
        endfor
    endfor
    if (N < 2)
        Print "FindDiameter_FromMask: fewer than 2 boundary pixels."
        Make /O /N=2 $(outBase+"_X"), $(outBase+"_Y")
        Wave xPair = $(outBase+"_X"); Wave yPair = $(outBase+"_Y")
        xPair = NaN; yPair = NaN
        return 0
    endif

    Make /O /N=(N) $(outBase+"_AllX"), $(outBase+"_AllY")
    Wave allX = $(outBase+"_AllX")
    Wave allY = $(outBase+"_AllY")

    Variable xOff = DimOffset(boundaryMask, 1)
    Variable yOff = DimOffset(boundaryMask, 0)
    Variable xDel = DimDelta(boundaryMask, 1)
    Variable yDel = DimDelta(boundaryMask, 0)

    Variable k=0
    for (r=0; r<nRows; r+=1)
        for (c=0; c<nCols; c+=1)
            if (boundaryMask[r][c] == 1)
                if (useScaling)
                    allX[k] = xOff + c * xDel
                    allY[k] = yOff + r * yDel
                else
                    allX[k] = c
                    allY[k] = r
                endif
                k += 1
            endif
        endfor
    endfor

    // --- 2) Pairwise distance^2 matrix ---
    MatrixOP /FREE tmpDx = colrepeat(allX, N) - rowrepeat(allX, N)
    MatrixOP /FREE tmpDy = colrepeat(allY, N) - rowrepeat(allY, N)
    Duplicate /O tmpDx $(outBase+"_D2")
    Wave dist2 = $(outBase+"_D2")
    dist2 = tmpDx^2 + tmpDy^2

    Variable i
    for (i=0; i<N; i+=1)
        dist2[i][i] = NaN     // ignore self-pairs
    endfor

    // --- 3) Argmax via WaveStats ---
    WaveStats /Q dist2
    Variable bestD2 = V_max
    Variable iBest  = V_maxRowLoc
    Variable jBest  = V_maxColLoc

    // --- 4) Emit 2-point pair ---
    Make /O /N=2 $(outBase+"_X"), $(outBase+"_Y")
    Wave xPair = $(outBase+"_X"); Wave yPair = $(outBase+"_Y")
    xPair[0] = allX[iBest]; yPair[0] = allY[iBest]
    xPair[1] = allX[jBest]; yPair[1] = allY[jBest]

    if (useScaling)
        SetScale /P x, 0, 1, "", xPair, yPair
    endif

    Variable diameter = sqrt(bestD2)
    Print "Diameter =", diameter, "  (indices in boundary set:", iBest, ",", jBest, ")"

    if (!keepD2)
        KillWaves /Z dist2
    endif
    KillWaves /Z tmpDx, tmpDy

    return diameter
End



//==============================================================================
// HasAnyNaN
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   w : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function HasAnyNaN(w)
    Wave w
    Variable nr = DimSize(w,0), nc = DimSize(w,1)
    Variable r, c
    for (r=0; r<nr; r+=1)
        for (c=0; c<nc; c+=1)
            if (numtype(w[r][c]) == 2)
                return 1
            endif
        endfor
    endfor
    return 0
End


//==============================================================================
// InBounds
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   nr : input
//   nc : input
//   r : input
//   c : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function InBounds(nr, nc, r, c)
    Variable nr, nc, r, c
    return (r>=0 && r<nr && c>=0 && c<nc)
End


//==============================================================================
// IsInsideRC
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   islandMask : input
//   useNaNMode : input
//   r : input
//   c : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function IsInsideRC(islandMask, useNaNMode, r, c)
    Wave islandMask
    Variable useNaNMode, r, c
    Variable val = islandMask[r][c]
    if (useNaNMode)
        return (numtype(val) != 2)
    else
        return (val != 0)
    endif
End


//==============================================================================
// SegmentStaysInside
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   islandMask : input
//   useNaNMode : input
//   c0 : input
//   r0 : input
//   c1 : input
//   r1 : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SegmentStaysInside(islandMask, useNaNMode, c0, r0, c1, r1)
    Wave islandMask
    Variable useNaNMode, c0, r0, c1, r1

    Variable nr = DimSize(islandMask,0)
    Variable nc = DimSize(islandMask,1)

    Variable dc = c1 - c0
    Variable dr = r1 - r0
    Variable steps = max(abs(dc), abs(dr))
    Variable s, fx, fy, xi, yi

    // --- handle degenerate case (same pixel) ---
    if (steps <= 0)
        if (!InBounds(nr, nc, r0, c0))
            return 0
        endif
        if (IsInsideRC(islandMask, useNaNMode, r0, c0))
            return 1
        else
            return 0
        endif
    endif

    for (s=0; s<=steps; s+=1)
        fx = c0 + dc * (s/steps)
        fy = r0 + dr * (s/steps)
        xi = floor(fx)
        yi = floor(fy)

        // check the supercover set: (xi,yi) and neighbors if we crossed grid lines
        Variable checkC, checkR

        // (xi, yi)
        if (!InBounds(nr,nc,yi,xi) || !IsInsideRC(islandMask, useNaNMode, yi, xi))
            return 0
        endif

        // if we crossed a vertical grid line, also check (xi+1, yi)
        if (fx != xi)
            checkC = xi+1; checkR = yi
            if (!InBounds(nr,nc,checkR,checkC) || !IsInsideRC(islandMask, useNaNMode, checkR, checkC))
                return 0
            endif
        endif

        // if we crossed a horizontal grid line, also check (xi, yi+1)
        if (fy != yi)
            checkC = xi; checkR = yi+1
            if (!InBounds(nr,nc,checkR,checkC) || !IsInsideRC(islandMask, useNaNMode, checkR, checkC))
                return 0
            endif
        endif

        // if we crossed both (corner), also check (xi+1, yi+1)
        if (fx != xi && fy != yi)
            checkC = xi+1; checkR = yi+1
            if (!InBounds(nr,nc,checkR,checkC) || !IsInsideRC(islandMask, useNaNMode, checkR, checkC))
                return 0
            endif
        endif
    endfor

    return 1
End


//==============================================================================
// FindLongestInternalChord
//
// Purpose:
//   Locate a feature or parameter from input data.
//
// Inputs:
//   boundaryMask : input
//   islandMask : input
//   useScaling : input
//   outBase : input
//
// Outputs:
//   return : numeric
//   outBase : output wave(s)
//
//==============================================================================
Function FindLongestInternalChord(boundaryMask, islandMask, useScaling, outBase)
    Wave boundaryMask, islandMask
    Variable useScaling
    String outBase

    if (WaveDims(boundaryMask)!=2 || WaveDims(islandMask)!=2)
        Abort "FindLongestInternalChord: both inputs must be 2D."
    endif
    if (DimSize(boundaryMask,0)!=DimSize(islandMask,0) || DimSize(boundaryMask,1)!=DimSize(islandMask,1))
        Abort "FindLongestInternalChord: size mismatch between boundaryMask and islandMask."
    endif

    Variable nr = DimSize(boundaryMask,0)
    Variable nc = DimSize(boundaryMask,1)

    // Determine membership mode once
    Variable useNaNMode = HasAnyNaN(islandMask)

    // Gather boundary pixel indices
    Variable r, c, N=0
    for (r=0; r<nr; r+=1)
        for (c=0; c<nc; c+=1)
            if (boundaryMask[r][c] == 1)
                N += 1
            endif
        endfor
    endfor
    if (N < 2)
        Print "FindLongestInternalChord: fewer than 2 boundary pixels."
        Make /O /N=2 $(outBase+"_X"), $(outBase+"_Y")
        Wave xPair0 = $(outBase+"_X"); Wave yPair0 = $(outBase+"_Y")
        xPair0 = NaN; yPair0 = NaN
        return 0
    endif

    Make /O /N=(N) $(outBase+"_AllCol"), $(outBase+"_AllRow")
    Wave allC = $(outBase+"_AllCol")
    Wave allR = $(outBase+"_AllRow")

    Variable k=0
    for (r=0; r<nr; r+=1)
        for (c=0; c<nc; c+=1)
            if (boundaryMask[r][c] == 1)
                allC[k] = c
                allR[k] = r
                k += 1
            endif
        endfor
    endfor

    // Brute-force over boundary pairs, validate with supercover
    Variable bestD2 = -1, bi=-1, bj=-1
    Variable i, j, dc, dr, d2

    for (i=0; i<N; i+=1)
        for (j=i+1; j<N; j+=1)
            dc = allC[i] - allC[j]
            dr = allR[i] - allR[j]
            d2 = dc*dc + dr*dr
            if (d2 <= bestD2)
                continue
            endif
            if (SegmentStaysInside(islandMask, useNaNMode, allC[i], allR[i], allC[j], allR[j]))
                bestD2 = d2
                bi = i
                bj = j
            endif
        endfor
    endfor

    Make /O /N=2 $(outBase+"_X"), $(outBase+"_Y")
    Wave xPair = $(outBase+"_X"); Wave yPair = $(outBase+"_Y")

    Variable xOff = DimOffset(boundaryMask,1)
    Variable yOff = DimOffset(boundaryMask,0)
    Variable xDel = DimDelta(boundaryMask,1)
    Variable yDel = DimDelta(boundaryMask,0)

    if (bi >= 0)
        if (useScaling)
            xPair[0] = xOff + allC[bi]*xDel
            yPair[0] = yOff + allR[bi]*yDel
            xPair[1] = xOff + allC[bj]*xDel
            yPair[1] = yOff + allR[bj]*yDel
            SetScale /P x, 0, 1, "", xPair, yPair
        else
            xPair[0] = allC[bi];  yPair[0] = allR[bi]
            xPair[1] = allC[bj];  yPair[1] = allR[bj]
        endif

        Variable lengthPix = sqrt(bestD2)
        Variable lengthOut = useScaling ? sqrt( (xPair[1]-xPair[0])^2 + (yPair[1]-yPair[0])^2 ) : lengthPix
        Print "Longest internal chord =", lengthOut, SelectString(useScaling, " (pixels)", WaveUnits(islandMask, 0 ))
        return lengthOut
    else
        xPair = NaN; yPair = NaN
        Print "No straight internal chord found."
        return 0
    endif
End


//==============================================================================
// FindOBB_FromBoundary
//
// Purpose:
//   Locate a feature or parameter from input data.
//
// Inputs:
//   boundaryMask : input
//   useScaling : input
//   outBase : input
//
// Outputs:
//   return : numeric
//   outBase : output wave(s)
//
//==============================================================================
Function FindOBB_FromBoundary(boundaryMask, useScaling, outBase)
    Wave boundaryMask
    Variable useScaling
    String outBase

    Variable nRows = DimSize(boundaryMask, 0)
    Variable nCols = DimSize(boundaryMask, 1)
    if (nRows <= 0 || nCols <= 0)
        Abort "FindOBB_FromBoundary: boundaryMask must be a non-empty 2D wave."
    endif

    // --- use wave scaling for physical coordinates ---
    Variable xOff = DimOffset(boundaryMask, 1)
    Variable yOff = DimOffset(boundaryMask, 0)
    Variable xDel = DimDelta(boundaryMask, 1)
    Variable yDel = DimDelta(boundaryMask, 0)

    // ---- 1) collect boundary points in *axis units* ----
    Variable r, c, count
    count = 0
    for (r=0; r<nRows; r+=1)
        for (c=0; c<nCols; c+=1)
            if (boundaryMask[r][c] != 0)
                count += 1
            endif
        endfor
    endfor

    if (count == 0)
        Make /O /D /N=5 $(outBase+"_RectX"), $(outBase+"_RectY")
        Wave rx0 = $(outBase+"_RectX")
        Wave ry0 = $(outBase+"_RectY")
        rx0 = NaN
        ry0 = NaN
        Make /O /D /N=2 $(outBase+"_SideLengths")
        Wave sl0 = $(outBase+"_SideLengths")
        sl0 = 0
        Make /O /D /N=1 $(outBase+"_AngleDeg")
        Wave ang0 = $(outBase+"_AngleDeg")
        ang0[0] = 0
        Print "FindOBB_FromBoundary: no boundary pixels found."
        return NaN
    endif

    Make /O /D /N=(count) $(outBase+"_AllX"), $(outBase+"_AllY")
    Wave allX = $(outBase+"_AllX")
    Wave allY = $(outBase+"_AllY")

    Variable k
    k = 0
    for (r=0; r<nRows; r+=1)
        for (c=0; c<nCols; c+=1)
            if (boundaryMask[r][c] != 0)
                // convert (row, col) -> (x,y) using wave scaling
                allX[k] = xOff + c * xDel
                allY[k] = yOff + r * yDel
                k += 1
            endif
        endfor
    endfor

    // ---- 2) PCA: covariance of boundary points in axis units ----
    Variable N = count
    Variable i
    Variable meanX = 0
    Variable meanY = 0
    for (i=0; i<N; i+=1)
        meanX += allX[i]
        meanY += allY[i]
    endfor
    meanX /= N
    meanY /= N

    Variable cxx = 0
    Variable cyy = 0
    Variable cxy = 0
    for (i=0; i<N; i+=1)
        Variable dx = allX[i] - meanX
        Variable dy = allY[i] - meanY
        cxx += dx*dx
        cyy += dy*dy
        cxy += dx*dy
    endfor
    cxx /= N
    cyy /= N
    cxy /= N

    // principal axis angle theta (radians)
    // tan(2θ) = 2 cxy / (cxx - cyy)
    Variable theta
    if (abs(cxx - cyy) < 1e-20 && abs(cxy) < 1e-20)
        theta = 0
    else
        theta = 0.5 * atan2(2*cxy, cxx - cyy)
    endif

    Variable cosT = cos(theta)
    Variable sinT = sin(theta)

    // rotate into PCA frame: x' = cosT*(x-meanX) + sinT*(y-meanY)
    //                         y' = -sinT*(x-meanX) + cosT*(y-meanY)
    Variable minXp =  1e300
    Variable maxXp = -1e300
    Variable minYp =  1e300
    Variable maxYp = -1e300

    for (i=0; i<N; i+=1)
        Variable px = allX[i] - meanX
        Variable py = allY[i] - meanY

        Variable xp =  cosT*px + sinT*py
        Variable yp = -sinT*px + cosT*py

        if (xp < minXp)
            minXp = xp
        endif
        if (xp > maxXp)
            maxXp = xp
        endif
        if (yp < minYp)
            minYp = yp
        endif
        if (yp > maxYp)
            maxYp = yp
        endif
    endfor

    Variable width  = maxXp - minXp      // axis units (e.g. nm)
    Variable height = maxYp - minYp
    Variable area   = width * height

    // ---- 3) corners in PCA frame, then back to axis coordinates ----
    Make /O /D /N=5 $(outBase+"_RectX"), $(outBase+"_RectY")
    Wave rx = $(outBase+"_RectX")
    Wave ry = $(outBase+"_RectY")

    // corners in PCA (x',y'): (minXp,minYp), (minXp,maxYp), (maxXp,maxYp), (maxXp,minYp)
    Variable x0p = minXp, y0p = minYp
    Variable x1p = minXp, y1p = maxYp
    Variable x2p = maxXp, y2p = maxYp
    Variable x3p = maxXp, y3p = minYp

    // back-transform: (x,y) = (meanX,meanY) + x' u + y' v
    // with u = (cosT, sinT), v = (-sinT, cosT)
    rx[0] = meanX + cosT*x0p - sinT*y0p
    ry[0] = meanY + sinT*x0p + cosT*y0p

    rx[1] = meanX + cosT*x1p - sinT*y1p
    ry[1] = meanY + sinT*x1p + cosT*y1p

    rx[2] = meanX + cosT*x2p - sinT*y2p
    ry[2] = meanY + sinT*x2p + cosT*y2p

    rx[3] = meanX + cosT*x3p - sinT*y3p
    ry[3] = meanY + sinT*x3p + cosT*y3p

    rx[4] = rx[0]
    ry[4] = ry[0]

    // ---- 4) side lengths and angle (all in axis units now) ----
    Make /O /D /N=2 $(outBase+"_SideLengths")
    Wave sl = $(outBase+"_SideLengths")
    sl[0] = width
    sl[1] = height

    Make /O /D /N=1 $(outBase+"_AngleDeg")
    Wave ang = $(outBase+"_AngleDeg")
    ang[0] = theta * 180/pi

    if (useScaling)
        SetScale /P x, 0, 1, "", rx, ry, sl, ang
    endif

    Print "OBB (PCA-based, scaled) area =", area
    Print "Side lengths (axis units):", width, height
    Print "Angle (deg, major axis wrt +x):", ang[0]

    return area
End




//==============================================================================
// MakeIslandMask_FromBoundary
//
// Purpose:
//   Construct output waves used by the SNS workflow.
//
// Inputs:
//   boundaryMask : input
//   outName : input
//
// Outputs:
//   return : numeric
//   outName : output wave(s)
//
//==============================================================================
Function MakeIslandMask_FromBoundary(boundaryMask, outName)
    Wave boundaryMask
    String outName

    Variable nr = DimSize(boundaryMask, 0)
    Variable nc = DimSize(boundaryMask, 1)
    if (nr <= 0 || nc <= 0)
        Abort "MakeIslandMask_FromBoundary: boundaryMask must be non-empty 2D wave."
    endif

    Make /O /U /B /N=(nr, nc) outsideMask
    outsideMask = 0

    // queue for BFS (row, col)
    Make /O /I /N=(nr*nc) qR, qC
    Variable head = 0, tail = 0
    Variable r, c, rr, cc
    Variable dr, dc

    // enqueue all edge pixels that are NOT boundary (value 0)
    for (c=0; c<nc; c+=1)
        if (boundaryMask[0][c] == 0)
            outsideMask[0][c] = 1
            qR[tail] = 0
            qC[tail] = c
            tail += 1
        endif
        if (boundaryMask[nr-1][c] == 0 && outsideMask[nr-1][c] == 1-1)   // outsideMask==0
            outsideMask[nr-1][c] = 1
            qR[tail] = nr-1
            qC[tail] = c
            tail += 1
        endif
    endfor

    for (r=0; r<nr; r+=1)
        if (boundaryMask[r][0] == 0 && outsideMask[r][0] == 1-1)
            outsideMask[r][0] = 1
            qR[tail] = r
            qC[tail] = 0
            tail += 1
        endif
        if (boundaryMask[r][nc-1] == 0 && outsideMask[r][nc-1] == 1-1)
            outsideMask[r][nc-1] = 1
            qR[tail] = r
            qC[tail] = nc-1
            tail += 1
        endif
    endfor

    // --- BFS over 4-neighbors, respecting Igor's do-while syntax ---
    do
        if (head >= tail)
            break
        endif

        r = qR[head]
        c = qC[head]
        head += 1

        // 4-connected neighborhood
        for (dr=-1; dr<=1; dr+=1)
            for (dc=-1; dc<=1; dc+=1)
                if (abs(dr)+abs(dc) != 1)
                    continue
                endif

                rr = r + dr
                cc = c + dc
                if (rr < 0 || rr >= nr || cc < 0 || cc >= nc)
                    continue
                endif
                if (boundaryMask[rr][cc] != 0)
                    continue                    // boundary = barrier
                endif
                if (outsideMask[rr][cc] == 0)
                    outsideMask[rr][cc] = 1
                    qR[tail] = rr
                    qC[tail] = cc
                    tail += 1
                endif
            endfor
        endfor

    while (1)

    // --- build island mask: 1 = boundary OR interior; 0 = outside ---
    Make /O /U /B /N=(nr, nc) $outName
    Wave islandMask = $outName
    islandMask = 0

    for (r=0; r<nr; r+=1)
        for (c=0; c<nc; c+=1)
            if (boundaryMask[r][c] == 1)
                islandMask[r][c] = 1
            else
                // zero pixel: interior if not marked as outside
                if (outsideMask[r][c] == 0)
                    islandMask[r][c] = 1
                endif
            endif
        endfor
    endfor

    SetScale /P x, DimOffset(boundaryMask,0), DimDelta(boundaryMask,0), "", islandMask
    SetScale /P y, DimOffset(boundaryMask,1), DimDelta(boundaryMask,1), "", islandMask

    KillWaves /Z outsideMask, qR, qC
End



//==============================================================================
// FindMaxAxisAlignedInnerRect
//
// Purpose:
//   Locate a feature or parameter from input data.
//
// Inputs:
//   islandMask : input
//   useScaling : input
//   outBase : input
//
// Outputs:
//   return : numeric
//   outBase : output wave(s)
//
//==============================================================================
Function FindMaxAxisAlignedInnerRect(islandMask, useScaling, outBase)
    Wave islandMask
    Variable useScaling
    String outBase

    Variable nr = DimSize(islandMask, 0)
    Variable nc = DimSize(islandMask, 1)
    if (nr <= 0 || nc <= 0)
        Abort "FindMaxAxisAlignedInnerRect: islandMask must be non-empty 2D wave."
    endif

    // histH: nc columns; stackIdx: needs nc+1 capacity (sentinel column)
    Make /O /I /N=(nc)   histH
    Make /O /I /N=(nc+1) stackIdx
    histH = 0

    Variable bestAreaPix = 0
    Variable bestR0=-1, bestR1=-1, bestC0=-1, bestC1=-1
    Variable r, c
    Variable stackTop
    Variable curH, topIndex, h, left, width, area

    for (r=0; r<nr; r+=1)
        // update histogram for this row
        for (c=0; c<nc; c+=1)
            if (islandMask[r][c] != 0)
                histH[c] += 1
            else
                histH[c] = 0
            endif
        endfor

        stackTop = 0

        // c goes from 0..nc; at c==nc we use curH=0 as sentinel
        for (c=0; c<=nc; c+=1)
            if (c < nc)
                curH = histH[c]
            else
                curH = 0      // sentinel to flush stack
            endif

            // emulate "while (stackTop > 0 && curH < histH[stackIdx[stackTop-1]])"
            do
                if (stackTop <= 0)
                    break
                endif
                if (!(curH < histH[ stackIdx[stackTop-1] ]))
                    break
                endif

                topIndex = stackIdx[stackTop-1]
                stackTop -= 1
                h = histH[topIndex]
                if (h > 0)
                    if (stackTop == 0)
                        left = -1
                    else
                        left = stackIdx[stackTop-1]
                    endif
                    width = c - left - 1
                    area  = h * width
                    if (area > bestAreaPix)
                        bestAreaPix = area
                        bestR1 = r
                        bestR0 = r - h + 1
                        bestC1 = c - 1
                        bestC0 = left + 1
                    endif
                endif
            while (1)

            // push current index (0..nc); stackIdx has size nc+1 so this is safe
            stackIdx[stackTop] = c
            stackTop += 1
        endfor
    endfor

    Make /O /D /N=5 $(outBase+"_InnerRectX"), $(outBase+"_InnerRectY")
    Wave rx = $(outBase+"_InnerRectX")
    Wave ry = $(outBase+"_InnerRectY")
    Make /O /D /N=2 $(outBase+"_InnerSideLengths")
    Wave sl = $(outBase+"_InnerSideLengths")

    if (bestAreaPix <= 0 || bestR0 < 0)
        rx = NaN
        ry = NaN
        sl = 0
        if (useScaling)
            SetScale /P x, 0, 1, "", rx, ry, sl
        endif
        Print "FindMaxAxisAlignedInnerRect: no nonzero interior region found."
        return 0
    endif

    // --- convert best rectangle from indices to axis units ---
    Variable xOff = DimOffset(islandMask, 1)
    Variable yOff = DimOffset(islandMask, 0)
    Variable xDel = DimDelta(islandMask, 1)
    Variable yDel = DimDelta(islandMask, 0)

    Variable c0 = bestC0, c1 = bestC1
    Variable r0 = bestR0, r1 = bestR1

    Variable xC0 = xOff + c0 * xDel
    Variable xC1 = xOff + c1 * xDel
    Variable yC0 = yOff + r0 * yDel
    Variable yC1 = yOff + r1 * yDel

    Variable dx = abs(xDel)
    Variable dy = abs(yDel)

    Variable xMinC = min(xC0, xC1)
    Variable xMaxC = max(xC0, xC1)
    Variable yMinC = min(yC0, yC1)
    Variable yMaxC = max(yC0, yC1)

    Variable xLeft   = xMinC - 0.5*dx
    Variable xRight  = xMaxC + 0.5*dx
    Variable yTop    = yMinC - 0.5*dy
    Variable yBottom = yMaxC + 0.5*dy

    sl[0] = xRight - xLeft
    sl[1] = yBottom - yTop

    rx[0] = xLeft;   ry[0] = yTop
    rx[1] = xLeft;   ry[1] = yBottom
    rx[2] = xRight;  ry[2] = yBottom
    rx[3] = xRight;  ry[3] = yTop
    rx[4] = xLeft;   ry[4] = yTop

    if (useScaling)
        SetScale /P x, 0, 1, "", rx, ry, sl
    endif

    Variable areaAxis = sl[0] * sl[1]
    // Print "Largest inner axis-aligned rectangle area =", areaAxis
    return areaAxis
End



//==============================================================================
// FindMinOuterBoundingRect
//
// Purpose:
//   Locate a feature or parameter from input data.
//
// Inputs:
//   islandMask : input
//   useScaling : input
//   outBase : input
//
// Outputs:
//   return : numeric
//   outBase : output wave(s)
//
//==============================================================================
Function FindMinOuterBoundingRect(islandMask, useScaling, outBase)
    Wave islandMask
    Variable useScaling
    String outBase

    Variable nr = DimSize(islandMask, 0)
    Variable nc = DimSize(islandMask, 1)
    if (nr <= 0 || nc <= 0)
        Abort "FindMinOuterBoundingRect: islandMask must be non-empty 2D wave."
    endif

    Variable r, c
    Variable rMin = nr, rMax = -1, cMin = nc, cMax = -1

    for (r=0; r<nr; r+=1)
        for (c=0; c<nc; c+=1)
            if (islandMask[r][c] != 0)
                if (r < rMin)
                    rMin = r
                endif
                if (r > rMax)
                    rMax = r
                endif
                if (c < cMin)
                    cMin = c
                endif
                if (c > cMax)
                    cMax = c
                endif
            endif
        endfor
    endfor

    Make /O /D /N=5 $(outBase+"_OuterRectX"), $(outBase+"_OuterRectY")
    Wave rx = $(outBase+"_OuterRectX")
    Wave ry = $(outBase+"_OuterRectY")
    Make /O /D /N=2 $(outBase+"_OuterSideLengths")
    Wave sl = $(outBase+"_OuterSideLengths")

    if (rMax < rMin || cMax < cMin)
        rx = NaN
        ry = NaN
        sl = 0
        if (useScaling)
            SetScale /P x, 0, 1, "", rx, ry, sl
        endif
        Print "FindMinOuterBoundingRect: no island pixels found."
        return 0
    endif

    // expand by one pixel, clamp to image
    Variable r0 = max(0,      rMin - 1)
    Variable r1 = min(nr - 1, rMax + 1)
    Variable c0 = max(0,      cMin - 1)
    Variable c1 = min(nc - 1, cMax + 1)

    Variable xOff = DimOffset(islandMask, 1)
    Variable yOff = DimOffset(islandMask, 0)
    Variable xDel = DimDelta(islandMask, 1)
    Variable yDel = DimDelta(islandMask, 0)

    Variable xC0 = xOff + c0 * xDel
    Variable xC1 = xOff + c1 * xDel
    Variable yC0 = yOff + r0 * yDel
    Variable yC1 = yOff + r1 * yDel

    Variable dx = abs(xDel)
    Variable dy = abs(yDel)

    Variable xMinC = min(xC0, xC1)
    Variable xMaxC = max(xC0, xC1)
    Variable yMinC = min(yC0, yC1)
    Variable yMaxC = max(yC0, yC1)

    Variable xLeft   = xMinC - 0.5*dx
    Variable xRight  = xMaxC + 0.5*dx
    Variable yTop    = yMinC - 0.5*dy
    Variable yBottom = yMaxC + 0.5*dy

    sl[0] = xRight - xLeft
    sl[1] = yBottom - yTop

    rx[0] = xLeft;   ry[0] = yTop
    rx[1] = xLeft;   ry[1] = yBottom
    rx[2] = xRight;  ry[2] = yBottom
    rx[3] = xRight;  ry[3] = yTop
    rx[4] = xLeft;   ry[4] = yTop

    if (useScaling)
        SetScale /P x, 0, 1, "", rx, ry, sl
    endif

    Variable areaAxis = sl[0] * sl[1]
    //Print "Outer bounding rectangle area =", areaAxis
    return areaAxis
End

