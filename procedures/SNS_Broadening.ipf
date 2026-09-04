#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma DefaultTab={3,20,4}

//==============================================================================
// SNS_Broadening.ipf
//
// Experimental-resolution and post-processing broadening helpers.
//
// Responsibilities:
//   - thermal broadening kernel
//   - lock-in modulation broadening kernel
//   - kernel half-width estimate
//   - convolution of DOS(E,B) with temperature + modulation resolution
//
// Dependencies:
//   SNS_Core.ipf
//     - constants kB_eV_per_K if used by caller conventions
//
// Notes:
//   This file should not solve ABS spectra, build DOS from branches,
//   build ray channels, or draw plots.
//==============================================================================

//==============================================================================
// SNS_ThermalKernelScalar
//
// Positive thermal-resolution kernel, -df/dE. The explicit kB argument lets
// legacy fit functions retain their historical numerical constant while the
// repository broadening pipeline uses its higher-precision value.
//==============================================================================
Function SNS_ThermalKernelScalar(energy_eV, temperature_K, kB_eV_per_K)
    Variable energy_eV, temperature_K, kB_eV_per_K

    if (temperature_K <= 0 || kB_eV_per_K <= 0)
        return 0
    endif

    Variable kBT_eV = kB_eV_per_K*temperature_K
    Variable halfArgument = energy_eV/(2*kBT_eV)
    return 1/(4*kBT_eV*cosh(halfArgument)^2)
End


//==============================================================================
// SNS_MakeThermalKernel
//
// Purpose:
//   Make normalized thermal broadening kernel:
//       K_T(E) = [4 kB T cosh^2(E / 2 kB T)]^-1
//
// Inputs:
//   E_axis : energy axis wave [eV], uniformly spaced.
//   T_K    : effective temperature [K].
//
// Outputs:
//   wK     : output kernel wave, same length/scaling as E_axis.
//
// Notes:
//   Kernel is normalized so sum(wK)*dE = 1.
//==============================================================================
Function SNS_MakeThermalKernel(E_axis, T_K, nameKT)
    Wave   E_axis
    Variable T_K
    String nameKT

    Variable NE = DimSize(E_axis, 0)
    if (NE <= 1)
        Abort "SNS_MakeThermalKernel: E_axis too small."
    endif

    Make/O/D/N=(NE) $nameKT
    Wave KT = $nameKT

    Variable dE = DimDelta(E_axis, 0)
    if (dE == 0)
        Abort "SNS_MakeThermalKernel: zero energy spacing."
    endif

    if (T_K <= 0)
        // delta kernel: no thermal broadening
        KT = 0
        KT[round((NE - 1)/2)] = 1      // centered delta kernel
        return 0
    endif

    Variable kB_eV = 8.617333262e-5       // eV/K
    Variable i, E, val

    for (i = 0; i < NE; i += 1)
        E = E_axis[i]
        val = SNS_ThermalKernelScalar(E, T_K, kB_eV)
        KT[i] = val
    endfor

    // include dE so discrete sum approximates ∫K(E)dE = 1
    KT *= dE

    // renormalize for safety
    Variable sumKT = sum(KT)
    if (sumKT > 0)
        KT /= sumKT
    endif

    return 0
End


//==============================================================================
// SNS_MakeModulationKernel
//
// Purpose:
//   Make normalized lock-in modulation kernel for sinusoidal bias modulation.
//
//   For modulation amplitude Vmod_eV, approximate kernel:
//       K(E) = 2/(pi Vmod^2) sqrt(Vmod^2 - E^2)
//   for |E| <= Vmod, otherwise 0.
//
// Inputs:
//   E_axis   : energy axis wave [eV], uniformly spaced.
//   Vmod_eV  : lock-in modulation amplitude [eV].
//
// Outputs:
//   wK       : output kernel wave, same length/scaling as E_axis.
//
// Notes:
//   Kernel is normalized so sum(wK)*dE = 1.
//==============================================================================
Function SNS_MakeModulationKernel(E_axis, Eac_eV, nameKM)
    Wave   E_axis
    Variable Eac_eV
    String nameKM

    Variable NE = DimSize(E_axis, 0)
    if (NE <= 1)
        Abort "SNS_MakeModulationKernel: E_axis too small."
    endif

    Make/O/D/N=(NE) $nameKM
    Wave KM = $nameKM

    Variable dE = DimDelta(E_axis, 0)
    if (dE == 0)
        Abort "SNS_MakeModulationKernel: zero energy spacing."
    endif

    if (Eac_eV <= 0)
        // delta kernel: no modulation broadening
        KM = 0
        KM[round((NE - 1)/2)] = 1
        return 0
    endif

    Variable i, E, val
    Variable E2 = Eac_eV*Eac_eV

    for (i = 0; i < NE; i += 1)
        E = E_axis[i]
        if (abs(E) <= Eac_eV)
            // K_mod(E) = (2/(π Eac^2)) * sqrt(Eac^2 - E^2)
            val = (2/(Pi*E2)) * sqrt(E2 - E*E)
        else
            val = 0
        endif
        KM[i] = val
    endfor

    // include dE so discrete sum approximates ∫K(E)dE = 1
    KM *= dE

    // renormalize
    Variable sumKM = sum(KM)
    if (sumKM > 0)
        KM /= sumKM
    endif

    return 0
End

//==============================================================================
// SNS_KernelHalfWidth
//
// Purpose:
//   Estimate half-width in points for a centered kernel wave.
//
// Inputs:
//   K   : normalized kernel wave.
//   eps : relative cutoff compared with max(abs(K)).
//
// Outputs:
//   return : half-width in points.
//==============================================================================
Function SNS_KernelHalfWidth(K, eps)
    Wave K
    Variable eps

    Variable NE = DimSize(K, 0)
    if (NE <= 1)
        return 0
    endif

    if (eps <= 0 || numtype(eps))
        eps = 1e-4
    endif

    WaveStats/Q K
    Variable Kmax = max(abs(V_min), abs(V_max))
    if (Kmax <= 0 || numtype(Kmax))
        return 0
    endif

    Variable center = round((NE - 1)/2)
    Variable i, dist
    Variable Nhalf = 0

    for (i = 0; i < NE; i += 1)
        if (abs(K[i]) > eps*Kmax)
            dist = abs(i - center)
            if (dist > Nhalf)
                Nhalf = dist
            endif
        endif
    endfor

    return round(Nhalf)
End

//==============================================================================
// SNS_ApplyDOS_Broadening_TplusMod
//
// Purpose:
//   Apply experimental broadening to DOS(E,B):
//     1) thermal broadening
//     2) lock-in modulation broadening
//
// Inputs:
//   DOS_EB_in  : input DOS wave, dimensions E x B.
//   nameDOSout : output wave name.
//
// Reads from root:SNS_Settings:
//   T_K        : effective temperature [K]
//   V_mod      : lock-in modulation amplitude [eV]
//
// Outputs:
//   nameDOSout : broadened DOS wave, same dimensions/scaling as DOS_EB_in.
//
// Notes:
//   Rebuilds the energy axis from the x-scaling of DOS_EB_in.
//   Uses constant edge padding before convolution.
//==============================================================================
Function SNS_ApplyDOS_Broadening_TplusMod(DOS_EB_in, nameDOSout)
    Wave   DOS_EB_in
    String nameDOSout

    Variable NE = DimSize(DOS_EB_in, 0)
    Variable nB = DimSize(DOS_EB_in, 1)
    if (NE <= 1 || nB <= 0)
        Abort "SNS_ApplyDOS_Broadening_TplusMod: DOS_EB_in invalid."
    endif

    // --- settings from initialization ---
    NVAR T_K   = root:SNS_Settings:T_K
    NVAR V_mod = root:SNS_Settings:V_mod      // modulation amplitude [eV]

    if (numtype(T_K) != 0 || numtype(V_mod) != 0)
        Abort "SNS_ApplyDOS_Broadening_TplusMod: T_K or V_mod not initialized."
    endif

    // --- rebuild energy axis from DOS_EB_in scaling ---
    Variable E0   = DimOffset(DOS_EB_in, 0)
    Variable dE   = DimDelta(DOS_EB_in, 0)
    String   uE   = WaveUnits(DOS_EB_in, 0)

    Make/FREE/D/N=(NE) E_axis
    E_axis = E0 + p*dE

    // --- output DOS ---
    Make/O/D/N=(NE, nB) $nameDOSout
    Wave DOS_EB_out = $nameDOSout
    DOS_EB_out = 0

    // --- build kernels once (on the local E_axis grid) ---
    SNS_MakeThermalKernel(E_axis, T_K,   "SNS_K_T")
    SNS_MakeModulationKernel(E_axis, V_mod, "SNS_K_Mod")
    Wave SNS_K_T   = SNS_K_T
    Wave SNS_K_Mod = SNS_K_Mod

    // effective half-widths (in points); eps can be tuned
    Variable epsK  = 1e-4
    Variable NpadT = SNS_KernelHalfWidth(SNS_K_T,  epsK)
    Variable NpadM = SNS_KernelHalfWidth(SNS_K_Mod, epsK)
    Variable Npad  = max(NpadT, NpadM)
    if (Npad < 0)
        Npad = 0
    endif
    if (Npad > NE - 1)
        Npad = NE - 1
    endif

    // temp waves
    Make/FREE/D/N=(NE)           tmp
    Make/FREE/D/N=(NE + 2*Npad)  tmpPad

    Variable iB, i

    for (iB = 0; iB < nB; iB += 1)

        // tmp holds the column DOS(E,B) for this iB
        tmp[] = DOS_EB_in[p][iB]

        // constant padding into tmpPad
        // left pad: copy first value
        for (i = 0; i < Npad; i += 1)
            tmpPad[i] = tmp[0]
        endfor

        // central region: straight copy
        for (i = 0; i < NE; i += 1)
            tmpPad[Npad + i] = tmp[i]
        endfor

        // right pad: copy last value
        for (i = 0; i < Npad; i += 1)
            tmpPad[Npad + NE + i] = tmp[NE - 1]
        endfor

        // thermal broadening
        Convolve/A SNS_K_T, tmpPad
        // modulation broadening
        Convolve/A SNS_K_Mod, tmpPad

        // extract central region back into tmp
        tmp[] = tmpPad[p + Npad]

        // write back
        DOS_EB_out[][iB] = tmp[p]
    endfor

    // copy scales from input DOS (x from energy, y from B)
    SetScale/P x, DimOffset(DOS_EB_in,0), DimDelta(DOS_EB_in,0), uE,                DOS_EB_out
    SetScale/P y, DimOffset(DOS_EB_in,1), DimDelta(DOS_EB_in,1), WaveUnits(DOS_EB_in,1), DOS_EB_out
    
    KillWaves/Z SNS_K_T, SNS_K_Mod

    return 0
End
