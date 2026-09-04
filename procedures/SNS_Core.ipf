#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma ModuleName=SNS_Core

//==============================================================================
// SNS_Core.ipf
//
// Core constants, parameter bundle, default settings, and small utilities
// for the S–N–S ABS / DOS / LDOS codebase (Igor Pro 9).
//
// No geometry, solver, map-building, or UI code lives here.
//==============================================================================


//============================================================
// 0. CONSTANTS
//============================================================

Constant HBAR_SI		= 1.054571817e-34      // [J*s]
Constant q_e			= 1.602176634e-19      // [C]
Constant HBAR_eVs	= 6.582119569e-16      // [eV*s]
Constant m_e_SI	  	= 9.10938356e-31      // [kg]
Constant phi_0		= 2.067833848e-15 // [T*m^2]

//============================================================
// SNS parameter bundle
//============================================================
Structure SNS_Params
    Variable m_eff       // m*/m_e
    Variable E_F         // [eV]

    Variable vF          // [m/s]
    Variable LambdaF     // [m]

    Variable DOS2D_eV_Area   // [states/(eV m^2)] 2DEG
    Variable DOS3D_eV_Vol    // [states/(eV m^3)] 3DEG at E_F

    Variable Delta       // [eV]
    Variable lambdaL     // [m]
    Variable Broadening  // [eV]
    Variable T_K         // [K]
    Variable V_mod       // [eV]

    Variable NE          // energy grid points

    // B grid
    Variable Bmin        // [T]
    Variable Bmax        // [T]
    Variable NB          // field grid points

    // Ray tracing / interface transparency
    Variable BTK_barrier // BTK-like interface barrier parameter

    // Band model
    // 0 = normal parabolic 2DEG
    // 1 = TI Dirac cone, high-mu semiclassical limit
    Variable SNS_bandModel
    Variable SNS_hbarvD_eVA
    Variable SNS_muDirac_eV

    // Diffusive / Usadel transport settings
    Variable lmfp_bulk_nm          // impurity/bulk mean free path [nm]
    Variable lmfp_edge_prefactor   // edge mfp prefactor, l_edge = prefactor*A/P

    // Default impurity diagnostics
    Variable impurity_density_default_nm2 // [1/nm^2]
    Variable sigma_tr_default_nm          // near-unitary transport cross section [nm]
    Variable lmfp_bulk_default_nm         // derived near-unitary impurity mfp [nm]

    // Vortex / delta-map settings
    Variable SNS_useVortex
    Variable SNS_nFlux
    Variable SNS_nIntSteps
    Variable SNS_useDeltamap

    // Field-dependent broadening settings
    Variable SNS_GammaBase_eV
    Variable SNS_useGammaDopp
    Variable SNS_useGammaZeeman
    Variable SNS_useGammaPair
    Variable SNS_useGammaUser
    Variable SNS_useGammaCoh
    Variable SNS_GammaDoppAlpha
    Variable SNS_GammaZeemanScale
    Variable SNS_GammaPairScale
    Variable SNS_GammaPair_Bc_T
EndStructure


//==============================================================================
// SNS_InitDefaultSettings
//
// Purpose:
//   Initialize root:SNS_Settings and derived SNS parameters.
//
// Optional inputs:
//   mEff_in              : effective mass m*/m_e.
//   EF_in                : Fermi energy [eV].
//   Delta_in             : superconducting gap [eV].
//   lambdaL_in           : London penetration depth [m].
//   Broad_in             : base broadening [eV].
//   T_K_in               : effective temperature [K].
//   Vmod_in              : modulation broadening scale [eV].
//   NE_in                : number of energy points.
//   Bmin_in, Bmax_in     : magnetic-field range [T].
//   NB_in                : number of B points.
//   useVortex_in         : 0/1 point-vortex mode.
//   nFlux_in             : vortex flux count.
//   nIntSteps_in         : vortex integration steps.
//   useDeltaMap_in       : 0/1 delta-map mode.
//   lmfpBulk_in          : bulk/impurity mean free path [nm].
//                          If omitted, estimated for near-unitary 2D scatterers:
//                          sigma_tr = 4/kF = (2/pi)*lambdaF,
//                          n_imp = 29/(100 nm * 100 nm),
//                          lmfp_bulk = 1/(n_imp*sigma_tr).
//   lmfpEdgePrefactor_in : edge mfp prefactor, l_edge = prefactor*A/P.
//   BTKbarrier_in        : BTK-like interface barrier parameter.
//
// Notes:
//   New optional inputs are appended at the end for backward compatibility.
//   vF, LambdaF, DOS2D_eV_Area, and DOS3D_eV_Vol are derived from m_eff/E_F.
//
//==============================================================================
Function SNS_LoadSettingsConfig(pathName, fileName)
    String pathName, fileName

    if (strlen(pathName) == 0)
        Abort "SNS_LoadSettingsConfig: pathName must not be empty."
    endif
    if (strlen(fileName) == 0)
        Abort "SNS_LoadSettingsConfig: fileName must not be empty."
    endif

    NewDataFolder/O root:SNS_SettingsDefaults
    String oldDF = GetDataFolder(1)
    SetDataFolder root:SNS_SettingsDefaults
    KillVariables/A/Z
    KillStrings/A/Z
    SetDataFolder $oldDF

    Variable refNum
    String line, key, valueText, varPath
    Variable tabPos, value, valueIsNaNLiteral

    Open/R/P=$pathName refNum as fileName
    do
        FReadLine refNum, line
        if (strlen(line) == 0)
            break
        endif

        line = ReplaceString("\r", line, "")
        line = ReplaceString("\n", line, "")
        line = SNS__TrimString(line)

        if (strlen(line) == 0)
            continue
        endif
        if (CmpStr(line[0], "#") == 0)
            continue
        endif

        tabPos = strsearch(line, "\t", 0)
        if (tabPos < 0)
            Abort "SNS_LoadSettingsConfig: expected tab-separated key/value line: " + line
        endif

        key = SNS__TrimString(line[0, tabPos - 1])
        valueText = SNS__TrimString(line[tabPos + 1, strlen(line) - 1])

        if (strlen(key) == 0 || strlen(valueText) == 0)
            Abort "SNS_LoadSettingsConfig: empty key or value in line: " + line
        endif

        valueIsNaNLiteral = 0
        if (CmpStr(valueText, "pi") == 0)
            value = pi
        elseif (CmpStr(valueText, "NaN") == 0)
            value = NaN
            valueIsNaNLiteral = 1
        else
            value = str2num(valueText)
        endif

        if (numtype(value) == 2 && !valueIsNaNLiteral)
            Abort "SNS_LoadSettingsConfig: non-numeric value for key " + key + ": " + valueText
        endif

        varPath = "root:SNS_SettingsDefaults:" + key
        Variable/G $varPath = value
    while (1)
    Close refNum

    String/G root:SNS_SettingsDefaults:s_configPathName = pathName
    String/G root:SNS_SettingsDefaults:s_configFileName = fileName
End


Function SNS__SettingsDefaultValue(key)
    String key

    String varPath = "root:SNS_SettingsDefaults:" + key
    NVAR/Z value = $varPath
    if (!NVAR_Exists(value))
        Abort "SNS_InitDefaultSettings: missing settings default '" + key + "'. Load a complete settings config with SNS_LoadSettingsConfig(...)."
    endif

    return value
End


Function/S SNS__TrimString(s)
    String s

    Variable i0 = 0
    Variable i1 = strlen(s) - 1

    do
        if (i0 > i1)
            return ""
        endif
        if (CmpStr(s[i0], " ") != 0 && CmpStr(s[i0], "\t") != 0)
            break
        endif
        i0 += 1
    while (1)

    do
        if (i1 < i0)
            return ""
        endif
        if (CmpStr(s[i1], " ") != 0 && CmpStr(s[i1], "\t") != 0)
            break
        endif
        i1 -= 1
    while (1)

    return s[i0, i1]
End


Function SNS_LoadBuiltInSettingsDefaults_Legacy()
    NewDataFolder/O root:SNS_SettingsDefaults
    String oldDF = GetDataFolder(1)
    SetDataFolder root:SNS_SettingsDefaults
    KillVariables/A/Z
    KillStrings/A/Z
    SetDataFolder $oldDF

    Variable/G root:SNS_SettingsDefaults:m_eff = 0.45
    Variable/G root:SNS_SettingsDefaults:E_F_eV = 0.429
    Variable/G root:SNS_SettingsDefaults:Delta_eV = 0.92e-3
    Variable/G root:SNS_SettingsDefaults:h_eff_m = 24e-9
    Variable/G root:SNS_SettingsDefaults:Broadening_eV = 5e-6
    Variable/G root:SNS_SettingsDefaults:T_K = 0.86
    Variable/G root:SNS_SettingsDefaults:V_mod_eV = 50e-6
    Variable/G root:SNS_SettingsDefaults:NE = 1001
    Variable/G root:SNS_SettingsDefaults:Bmin_T = -0.5
    Variable/G root:SNS_SettingsDefaults:Bmax_T = 0.5
    Variable/G root:SNS_SettingsDefaults:NB = 1001
    Variable/G root:SNS_SettingsDefaults:BTK_barrier = 0.05
    Variable/G root:SNS_SettingsDefaults:bandModel = 0
    Variable/G root:SNS_SettingsDefaults:hbarvD_eVA = 3.0
    Variable/G root:SNS_SettingsDefaults:muDirac_eV = 0.2
    Variable/G root:SNS_SettingsDefaults:lmfp_edge_prefactor = pi
    Variable/G root:SNS_SettingsDefaults:useVortex = 0
    Variable/G root:SNS_SettingsDefaults:nFlux = 0
    Variable/G root:SNS_SettingsDefaults:nIntSteps = 64
    Variable/G root:SNS_SettingsDefaults:useDeltaMap = 0
    Variable/G root:SNS_SettingsDefaults:GammaBase_eV = NaN
    Variable/G root:SNS_SettingsDefaults:useGammaDopp = 0
    Variable/G root:SNS_SettingsDefaults:useGammaZeeman = 0
    Variable/G root:SNS_SettingsDefaults:useGammaPair = 1
    Variable/G root:SNS_SettingsDefaults:useGammaUser = 0
    Variable/G root:SNS_SettingsDefaults:useGammaCoh = 0
    Variable/G root:SNS_SettingsDefaults:GammaDoppAlpha = 0.01
    Variable/G root:SNS_SettingsDefaults:GammaZeemanScale = 0
    Variable/G root:SNS_SettingsDefaults:GammaPairScale_eV_T2 = 0

    String/G root:SNS_SettingsDefaults:s_configPathName = "<built-in legacy>"
    String/G root:SNS_SettingsDefaults:s_configFileName = "<built-in legacy>"
End


Function SNS_InitDefaultSettings([mEff_in, EF_in, Delta_in, lambdaL_in, Broad_in, T_K_in, Vmod_in, NE_in, Bmin_in, Bmax_in, NB_in, useVortex_in, nFlux_in, nIntSteps_in, useDeltaMap_in, lmfpBulk_in, lmfpEdgePrefactor_in, BTKbarrier_in,bandModel_in, hbarvD_eVA_in, muDirac_eV_in])
    Variable mEff_in, EF_in
    Variable Delta_in, lambdaL_in, Broad_in, T_K_in, Vmod_in, NE_in
    Variable Bmin_in, Bmax_in, NB_in
    Variable useVortex_in, nFlux_in, nIntSteps_in, useDeltaMap_in
    Variable lmfpBulk_in, lmfpEdgePrefactor_in, BTKbarrier_in
    Variable bandModel_in, hbarvD_eVA_in, muDirac_eV_in

    if (!DataFolderExists("root:SNS_SettingsDefaults:"))
        Abort "SNS_InitDefaultSettings: no settings defaults loaded. Run SNS_LoadSettingsConfig(...) first, or SNS_InitDefaultSettings_LegacyBuiltIn(...) for debug-only built-in defaults."
    endif

    NewDataFolder/O root:SNS_Settings
    
	 // ------------------------------
    // Band model
    // ------------------------------
    Variable/G root:SNS_Settings:SNS_bandModel = round(SNS__SettingsDefaultValue("bandModel"))
    Variable/G root:SNS_Settings:SNS_hbarvD_eVA = SNS__SettingsDefaultValue("hbarvD_eVA")
    Variable/G root:SNS_Settings:SNS_muDirac_eV = SNS__SettingsDefaultValue("muDirac_eV")

    NVAR SNS_bandModel  = root:SNS_Settings:SNS_bandModel
    NVAR SNS_hbarvD_eVA = root:SNS_Settings:SNS_hbarvD_eVA
    NVAR SNS_muDirac_eV = root:SNS_Settings:SNS_muDirac_eV

    if (!ParamIsDefault(bandModel_in))
        SNS_bandModel = round(bandModel_in)
    endif
    if (!ParamIsDefault(hbarvD_eVA_in))
        SNS_hbarvD_eVA = hbarvD_eVA_in
    endif
    if (!ParamIsDefault(muDirac_eV_in))
        SNS_muDirac_eV = muDirac_eV_in
    endif
	
    // ------------------------------
    // Electronic parameters
    // ------------------------------
    Variable/G root:SNS_Settings:m_eff = SNS__SettingsDefaultValue("m_eff")
    Variable/G root:SNS_Settings:E_F   = SNS__SettingsDefaultValue("E_F_eV")

    NVAR m_eff = root:SNS_Settings:m_eff
    NVAR E_F   = root:SNS_Settings:E_F

    if (!ParamIsDefault(mEff_in))
        m_eff = mEff_in
    endif
    if (!ParamIsDefault(EF_in))
        E_F = EF_in
    endif

    // ------------------------------
    // Ray tracing / interface transparency
    // ------------------------------
    Variable/G root:SNS_Settings:BTK_barrier = SNS__SettingsDefaultValue("BTK_barrier")
    NVAR BTK_barrier = root:SNS_Settings:BTK_barrier

    if (!ParamIsDefault(BTKbarrier_in))
        BTK_barrier = BTKbarrier_in
    endif

    // ------------------------------
    // Derived Fermi quantities
    // ------------------------------
    Variable m_star = m_eff * m_e_SI
    Variable EF_J   = E_F * q_e
    Variable kF     = sqrt(2 * m_star * EF_J) / HBAR_SI

    Variable/G root:SNS_Settings:vF
    Variable/G root:SNS_Settings:LambdaF

    NVAR vF      = root:SNS_Settings:vF
    NVAR LambdaF = root:SNS_Settings:LambdaF

    vF      = HBAR_SI * kF / m_star
    LambdaF = 2 * pi / kF
    
    // ------------------------------
    // Optional TI Dirac override
    // ------------------------------
    if (SNS_bandModel == 1)

        if (SNS_hbarvD_eVA <= 0 || numtype(SNS_hbarvD_eVA) != 0)
            Abort "SNS_InitDefaultSettings: TI mode requires SNS_hbarvD_eVA > 0."
        endif
        if (abs(SNS_muDirac_eV) <= 0 || numtype(SNS_muDirac_eV) != 0)
            Abort "SNS_InitDefaultSettings: TI mode requires nonzero SNS_muDirac_eV."
        endif

        vF = SNS_hbarvD_eVA / (HBAR_eVs * 1e10)
        LambdaF = (2*pi*SNS_hbarvD_eVA/abs(SNS_muDirac_eV)) * 1e-10
        E_F = SNS_muDirac_eV
    endif

    // ------------------------------
    // Diffusive / Usadel transport settings
    // ------------------------------
    Variable nImp_default_nm2 = 29 / (100.0 * 100.0)
    Variable sigmaTr_default_nm = (2.0 / pi) * (LambdaF * 1e9)
    Variable lmfpBulkDefault_nm = 1 / (nImp_default_nm2 * sigmaTr_default_nm)

    Variable/G root:SNS_Settings:impurity_density_default_nm2 = nImp_default_nm2
    Variable/G root:SNS_Settings:sigma_tr_default_nm = sigmaTr_default_nm
    Variable/G root:SNS_Settings:lmfp_bulk_default_nm = lmfpBulkDefault_nm

    Variable/G root:SNS_Settings:lmfp_bulk_nm = lmfpBulkDefault_nm
    Variable/G root:SNS_Settings:lmfp_edge_prefactor = SNS__SettingsDefaultValue("lmfp_edge_prefactor")

    NVAR lmfp_bulk_nm = root:SNS_Settings:lmfp_bulk_nm
    NVAR lmfp_edge_prefactor = root:SNS_Settings:lmfp_edge_prefactor

    if (!ParamIsDefault(lmfpBulk_in))
        lmfp_bulk_nm = lmfpBulk_in
    endif
    if (!ParamIsDefault(lmfpEdgePrefactor_in))
        lmfp_edge_prefactor = lmfpEdgePrefactor_in
    endif

    // ------------------------------
    // DOS for free 2DEG / 3DEG, including spin degeneracy
    // ------------------------------
    Variable DOS2D_SI = m_star / (pi * HBAR_SI^2)
    Variable DOS3D_SI = (1.0 / (2 * pi^2)) * (2 * m_star / (HBAR_SI^2))^1.5 * sqrt(EF_J)

    Variable/G root:SNS_Settings:DOS2D_eV_Area
    Variable/G root:SNS_Settings:DOS3D_eV_Vol

    NVAR DOS2D_eV_Area = root:SNS_Settings:DOS2D_eV_Area
    NVAR DOS3D_eV_Vol  = root:SNS_Settings:DOS3D_eV_Vol

    DOS2D_eV_Area = DOS2D_SI * q_e
    DOS3D_eV_Vol  = DOS3D_SI * q_e

    // ------------------------------
    // Physical / numerical parameters
    // ------------------------------
    Variable/G root:SNS_Settings:Delta      = SNS__SettingsDefaultValue("Delta_eV")
    Variable/G root:SNS_Settings:lambdaL    = SNS__SettingsDefaultValue("h_eff_m")
    Variable/G root:SNS_Settings:Broadening = SNS__SettingsDefaultValue("Broadening_eV")
    Variable/G root:SNS_Settings:T_K        = SNS__SettingsDefaultValue("T_K")
    Variable/G root:SNS_Settings:V_mod      = SNS__SettingsDefaultValue("V_mod_eV")
    Variable/G root:SNS_Settings:NE         = round(SNS__SettingsDefaultValue("NE"))

    NVAR Delta      = root:SNS_Settings:Delta
    NVAR lambdaL    = root:SNS_Settings:lambdaL
    NVAR Broadening = root:SNS_Settings:Broadening
    NVAR T_K        = root:SNS_Settings:T_K
    NVAR V_mod      = root:SNS_Settings:V_mod
    NVAR NE         = root:SNS_Settings:NE

    if (!ParamIsDefault(Delta_in))
        Delta = Delta_in
    endif
    if (!ParamIsDefault(lambdaL_in))
        lambdaL = lambdaL_in
    endif
    if (!ParamIsDefault(Broad_in))
        Broadening = Broad_in
    endif
    if (!ParamIsDefault(T_K_in))
        T_K = T_K_in
    endif
    if (!ParamIsDefault(Vmod_in))
        V_mod = Vmod_in
    endif
    if (!ParamIsDefault(NE_in))
        NE = round(NE_in)
    endif

	 // ------------------------------
	 // Gamma broadening model settings
	 // ------------------------------
	 Variable/G root:SNS_Settings:SNS_GammaBase_eV = SNS__SettingsDefaultValue("GammaBase_eV")
	
	 Variable/G root:SNS_Settings:SNS_useGammaDopp   = round(SNS__SettingsDefaultValue("useGammaDopp"))
	 Variable/G root:SNS_Settings:SNS_useGammaZeeman = round(SNS__SettingsDefaultValue("useGammaZeeman"))
	 Variable/G root:SNS_Settings:SNS_useGammaPair   = round(SNS__SettingsDefaultValue("useGammaPair"))
	 Variable/G root:SNS_Settings:SNS_useGammaUser   = round(SNS__SettingsDefaultValue("useGammaUser"))
	 Variable/G root:SNS_Settings:SNS_useGammaCoh    = round(SNS__SettingsDefaultValue("useGammaCoh"))
	
	 Variable/G root:SNS_Settings:SNS_GammaDoppAlpha   = SNS__SettingsDefaultValue("GammaDoppAlpha")
	 Variable/G root:SNS_Settings:SNS_GammaZeemanScale = SNS__SettingsDefaultValue("GammaZeemanScale")
	 Variable/G root:SNS_Settings:SNS_GammaPairScale = SNS__SettingsDefaultValue("GammaPairScale_eV_T2")

    // ------------------------------
    // B-grid
    // ------------------------------
    Variable/G root:SNS_Settings:Bmin = SNS__SettingsDefaultValue("Bmin_T")
    Variable/G root:SNS_Settings:Bmax = SNS__SettingsDefaultValue("Bmax_T")
    Variable/G root:SNS_Settings:NB   = round(SNS__SettingsDefaultValue("NB"))

    NVAR Bmin = root:SNS_Settings:Bmin
    NVAR Bmax = root:SNS_Settings:Bmax
    NVAR NB   = root:SNS_Settings:NB

    if (!ParamIsDefault(Bmin_in))
        Bmin = Bmin_in
    endif
    if (!ParamIsDefault(Bmax_in))
        Bmax = Bmax_in
    endif
    if (!ParamIsDefault(NB_in))
        NB = round(NB_in)
    endif

    Make/O/D/N=(NB) root:SNS_Settings:B_T
    Wave B_T = root:SNS_Settings:B_T
    SetScale/P x, Bmin, (Bmax - Bmin) / (NB - 1), "T", B_T
    B_T = x

    // ------------------------------
    // Vortex / delta-map settings
    // ------------------------------
    Variable/G root:SNS_Settings:SNS_useVortex   = round(SNS__SettingsDefaultValue("useVortex"))
    Variable/G root:SNS_Settings:SNS_nFlux       = SNS__SettingsDefaultValue("nFlux")
    Variable/G root:SNS_Settings:SNS_nIntSteps   = round(SNS__SettingsDefaultValue("nIntSteps"))
    Variable/G root:SNS_Settings:SNS_useDeltamap = round(SNS__SettingsDefaultValue("useDeltaMap"))

    NVAR SNS_useVortex   = root:SNS_Settings:SNS_useVortex
    NVAR SNS_nFlux       = root:SNS_Settings:SNS_nFlux
    NVAR SNS_nIntSteps   = root:SNS_Settings:SNS_nIntSteps
    NVAR SNS_useDeltamap = root:SNS_Settings:SNS_useDeltamap

    if (!ParamIsDefault(useVortex_in))
        SNS_useVortex = useVortex_in
    endif
    if (!ParamIsDefault(nFlux_in))
        SNS_nFlux = nFlux_in
    endif
    if (!ParamIsDefault(nIntSteps_in))
        SNS_nIntSteps = round(nIntSteps_in)
    endif
    if (!ParamIsDefault(useDeltaMap_in))
        SNS_useDeltamap = useDeltaMap_in
    endif
End


Function SNS_InitDefaultSettings_LegacyBuiltIn([mEff_in, EF_in, Delta_in, lambdaL_in, Broad_in, T_K_in, Vmod_in, NE_in, Bmin_in, Bmax_in, NB_in, useVortex_in, nFlux_in, nIntSteps_in, useDeltaMap_in, lmfpBulk_in, lmfpEdgePrefactor_in, BTKbarrier_in,bandModel_in, hbarvD_eVA_in, muDirac_eV_in])
    Variable mEff_in, EF_in
    Variable Delta_in, lambdaL_in, Broad_in, T_K_in, Vmod_in, NE_in
    Variable Bmin_in, Bmax_in, NB_in
    Variable useVortex_in, nFlux_in, nIntSteps_in, useDeltaMap_in
    Variable lmfpBulk_in, lmfpEdgePrefactor_in, BTKbarrier_in
    Variable bandModel_in, hbarvD_eVA_in, muDirac_eV_in

    NewDataFolder/O root:SNS_Settings
    
	 // ------------------------------
    // Band model
    // ------------------------------
    // 0 = parabolic 2DEG, default
    // 1 = TI Dirac cone, high-mu semiclassical limit
    Variable/G root:SNS_Settings:SNS_bandModel
    Variable/G root:SNS_Settings:SNS_hbarvD_eVA
    Variable/G root:SNS_Settings:SNS_muDirac_eV

    NVAR SNS_bandModel  = root:SNS_Settings:SNS_bandModel
    NVAR SNS_hbarvD_eVA = root:SNS_Settings:SNS_hbarvD_eVA
    NVAR SNS_muDirac_eV = root:SNS_Settings:SNS_muDirac_eV

    // Default MUST be normal 2DEG unless explicitly overridden.
    if (ParamIsDefault(bandModel_in))
        SNS_bandModel = 0
    else
        SNS_bandModel = round(bandModel_in)
    endif

    // Material defaults, but allow explicit override.
    if (ParamIsDefault(hbarvD_eVA_in))
        if (numtype(SNS_hbarvD_eVA) != 0)
            SNS_hbarvD_eVA = 3.0
        endif
    else
        SNS_hbarvD_eVA = hbarvD_eVA_in
    endif

    if (ParamIsDefault(muDirac_eV_in))
        if (numtype(SNS_muDirac_eV) != 0)
            SNS_muDirac_eV = 0.2
        endif
    else
        SNS_muDirac_eV = muDirac_eV_in
    endif
	
    // ------------------------------
    // Electronic parameters
    // ------------------------------
    Variable/G root:SNS_Settings:m_eff = 0.45
    Variable/G root:SNS_Settings:E_F   = 0.429

    NVAR m_eff = root:SNS_Settings:m_eff
    NVAR E_F   = root:SNS_Settings:E_F

    if (!ParamIsDefault(mEff_in))
        m_eff = mEff_in
    endif
    if (!ParamIsDefault(EF_in))
        E_F = EF_in
    endif

    // ------------------------------
    // Ray tracing / interface transparency
    // ------------------------------
    Variable/G root:SNS_Settings:BTK_barrier = 0.05
    NVAR BTK_barrier = root:SNS_Settings:BTK_barrier

    if (!ParamIsDefault(BTKbarrier_in))
        BTK_barrier = BTKbarrier_in
    endif

    // ------------------------------
    // Derived Fermi quantities
    // ------------------------------
    Variable m_star = m_eff * m_e_SI
    Variable EF_J   = E_F * q_e
    Variable kF     = sqrt(2 * m_star * EF_J) / HBAR_SI

    Variable/G root:SNS_Settings:vF
    Variable/G root:SNS_Settings:LambdaF

    NVAR vF      = root:SNS_Settings:vF
    NVAR LambdaF = root:SNS_Settings:LambdaF

    vF      = HBAR_SI * kF / m_star
    LambdaF = 2 * pi / kF
    
    // ------------------------------
    // Optional TI Dirac override
    // ------------------------------
    if (SNS_bandModel == 1)

        if (SNS_hbarvD_eVA <= 0 || numtype(SNS_hbarvD_eVA) != 0)
            Abort "SNS_InitDefaultSettings: TI mode requires SNS_hbarvD_eVA > 0."
        endif
        if (abs(SNS_muDirac_eV) <= 0 || numtype(SNS_muDirac_eV) != 0)
            Abort "SNS_InitDefaultSettings: TI mode requires nonzero SNS_muDirac_eV."
        endif

        // hbar*vD [eV A] = HBAR_eVs * vD[m/s] * 1e10
        vF = SNS_hbarvD_eVA / (HBAR_eVs * 1e10)

        // kF = |mu|/(hbar*vD), lambdaF = 2pi/kF.
        // hbar*vD is in eV A, so lambdaF is first in A, then converted to m.
        LambdaF = (2*pi*SNS_hbarvD_eVA/abs(SNS_muDirac_eV)) * 1e-10

        // In TI mode E_F means mu from Dirac point.
        E_F = SNS_muDirac_eV
    endif

    // ------------------------------
    // Diffusive / Usadel transport settings
    // ------------------------------
    // Default impurity-limited bulk mean free path assumes near-unitary
    // scatterers in a 2D metal:
    //
    //   sigma_tr = 4/kF = (2/pi) * lambdaF
    //   n_imp    = 29 / (100 nm * 100 nm)
    //   l_imp    = 1 / (n_imp * sigma_tr)
    //
    // This is distinct from the mean impurity spacing 1/sqrt(n_imp).

    Variable nImp_default_nm2 = 29 / (100.0 * 100.0)
    Variable sigmaTr_default_nm = (2.0 / pi) * (LambdaF * 1e9)
    Variable lmfpBulkDefault_nm = 1 / (nImp_default_nm2 * sigmaTr_default_nm)

    Variable/G root:SNS_Settings:impurity_density_default_nm2 = nImp_default_nm2
    Variable/G root:SNS_Settings:sigma_tr_default_nm = sigmaTr_default_nm
    Variable/G root:SNS_Settings:lmfp_bulk_default_nm = lmfpBulkDefault_nm

    Variable/G root:SNS_Settings:lmfp_bulk_nm = lmfpBulkDefault_nm
    Variable/G root:SNS_Settings:lmfp_edge_prefactor = pi

    NVAR lmfp_bulk_nm = root:SNS_Settings:lmfp_bulk_nm
    NVAR lmfp_edge_prefactor = root:SNS_Settings:lmfp_edge_prefactor

    if (!ParamIsDefault(lmfpBulk_in))
        lmfp_bulk_nm = lmfpBulk_in
    endif
    if (!ParamIsDefault(lmfpEdgePrefactor_in))
        lmfp_edge_prefactor = lmfpEdgePrefactor_in
    endif

    // ------------------------------
    // DOS for free 2DEG / 3DEG, including spin degeneracy
    // ------------------------------
    Variable DOS2D_SI = m_star / (pi * HBAR_SI^2)
    Variable DOS3D_SI = (1.0 / (2 * pi^2)) * (2 * m_star / (HBAR_SI^2))^1.5 * sqrt(EF_J)

    Variable/G root:SNS_Settings:DOS2D_eV_Area
    Variable/G root:SNS_Settings:DOS3D_eV_Vol

    NVAR DOS2D_eV_Area = root:SNS_Settings:DOS2D_eV_Area
    NVAR DOS3D_eV_Vol  = root:SNS_Settings:DOS3D_eV_Vol

    DOS2D_eV_Area = DOS2D_SI * q_e
    DOS3D_eV_Vol  = DOS3D_SI * q_e

    // ------------------------------
    // Physical / numerical parameters
    // ------------------------------
    Variable/G root:SNS_Settings:Delta      = 0.92e-3
    Variable/G root:SNS_Settings:lambdaL    = 24e-9 //this is really the effective height. lambdaL is a legacy name that needs to be fixed eventually.
    Variable/G root:SNS_Settings:Broadening = 5e-6
    Variable/G root:SNS_Settings:T_K        = 0.86
    Variable/G root:SNS_Settings:V_mod      = 50e-6
    Variable/G root:SNS_Settings:NE         = 1001

    NVAR Delta      = root:SNS_Settings:Delta
    NVAR lambdaL    = root:SNS_Settings:lambdaL
    NVAR Broadening = root:SNS_Settings:Broadening
    NVAR T_K        = root:SNS_Settings:T_K
    NVAR V_mod      = root:SNS_Settings:V_mod
    NVAR NE         = root:SNS_Settings:NE

    if (!ParamIsDefault(Delta_in))
        Delta = Delta_in
    endif
    if (!ParamIsDefault(lambdaL_in))
        lambdaL = lambdaL_in
    endif
    if (!ParamIsDefault(Broad_in))
        Broadening = Broad_in
    endif
    if (!ParamIsDefault(T_K_in))
        T_K = T_K_in
    endif
    if (!ParamIsDefault(Vmod_in))
        V_mod = Vmod_in
    endif
    if (!ParamIsDefault(NE_in))
        NE = round(NE_in)
    endif

	 // ------------------------------
	 // Gamma broadening model settings
	 // ------------------------------
	 Variable/G root:SNS_Settings:SNS_GammaBase_eV = NaN
	
	 Variable/G root:SNS_Settings:SNS_useGammaDopp   = 0
	 Variable/G root:SNS_Settings:SNS_useGammaZeeman = 0
	 Variable/G root:SNS_Settings:SNS_useGammaPair   = 1
	 Variable/G root:SNS_Settings:SNS_useGammaUser   = 0
	 Variable/G root:SNS_Settings:SNS_useGammaCoh    = 0
	
	 Variable/G root:SNS_Settings:SNS_GammaDoppAlpha   = 0.01
	 Variable/G root:SNS_Settings:SNS_GammaZeemanScale = 0
	
	 // Quadratic depairing coefficient [eV/T^2].
	 // Used as: Gamma_pair(B) = SNS_GammaPairScale * B^2.
	 // Keep default zero so the zero-field LDOS broadening remains unchanged
	 // unless explicitly fitted/enabled in the B-dependent fit.
	 Variable/G root:SNS_Settings:SNS_GammaPairScale = 0

    // ------------------------------
    // B-grid
    // ------------------------------
    Variable/G root:SNS_Settings:Bmin = -0.5
    Variable/G root:SNS_Settings:Bmax =  0.5
    Variable/G root:SNS_Settings:NB   = 1001

    NVAR Bmin = root:SNS_Settings:Bmin
    NVAR Bmax = root:SNS_Settings:Bmax
    NVAR NB   = root:SNS_Settings:NB

    if (!ParamIsDefault(Bmin_in))
        Bmin = Bmin_in
    endif
    if (!ParamIsDefault(Bmax_in))
        Bmax = Bmax_in
    endif
    if (!ParamIsDefault(NB_in))
        NB = round(NB_in)
    endif

    Make/O/D/N=(NB) root:SNS_Settings:B_T
    Wave B_T = root:SNS_Settings:B_T
    SetScale/P x, Bmin, (Bmax - Bmin) / (NB - 1), "T", B_T
    B_T = x

    // ------------------------------
    // Vortex / delta-map settings
    // ------------------------------
    Variable/G root:SNS_Settings:SNS_useVortex   = 0
    Variable/G root:SNS_Settings:SNS_nFlux       = 0
    Variable/G root:SNS_Settings:SNS_nIntSteps   = 64
    Variable/G root:SNS_Settings:SNS_useDeltamap = 0

    NVAR SNS_useVortex   = root:SNS_Settings:SNS_useVortex
    NVAR SNS_nFlux       = root:SNS_Settings:SNS_nFlux
    NVAR SNS_nIntSteps   = root:SNS_Settings:SNS_nIntSteps
    NVAR SNS_useDeltamap = root:SNS_Settings:SNS_useDeltamap

    if (!ParamIsDefault(useVortex_in))
        SNS_useVortex = useVortex_in
    endif
    if (!ParamIsDefault(nFlux_in))
        SNS_nFlux = nFlux_in
    endif
    if (!ParamIsDefault(nIntSteps_in))
        SNS_nIntSteps = round(nIntSteps_in)
    endif
    if (!ParamIsDefault(useDeltaMap_in))
        SNS_useDeltamap = useDeltaMap_in
    endif
End


//==============================================================================
// SNS_LoadParams
//
// Purpose:
//   Load root:SNS_Settings into SNS_Params.
//
// Notes:
//   Existing required legacy fields are loaded directly.
//   Newer fields use NVAR/Z fallback defaults for backward compatibility.
//
//==============================================================================
Function SNS_LoadParams(p)
    STRUCT SNS_Params &p

    NVAR m_eff      = root:SNS_Settings:m_eff
    NVAR E_F        = root:SNS_Settings:E_F
    NVAR vF         = root:SNS_Settings:vF
    NVAR LambdaF    = root:SNS_Settings:LambdaF
    NVAR Delta      = root:SNS_Settings:Delta
    NVAR lambdaL    = root:SNS_Settings:lambdaL
    NVAR Broadening = root:SNS_Settings:Broadening
    NVAR T_K        = root:SNS_Settings:T_K
    NVAR V_mod      = root:SNS_Settings:V_mod
    NVAR NE         = root:SNS_Settings:NE
    NVAR DOS2D_eV_Area = root:SNS_Settings:DOS2D_eV_Area
    NVAR DOS3D_eV_Vol  = root:SNS_Settings:DOS3D_eV_Vol

    p.m_eff         = m_eff
    p.E_F           = E_F
    p.vF            = vF
    p.LambdaF       = LambdaF
    p.DOS2D_eV_Area = DOS2D_eV_Area
    p.DOS3D_eV_Vol  = DOS3D_eV_Vol
    p.Delta         = Delta
    p.lambdaL       = lambdaL
    p.Broadening    = Broadening
    p.T_K           = T_K
    p.V_mod         = V_mod
    p.NE            = NE

    // ------------------------------
    // B-grid, backward-compatible
    // ------------------------------
    NVAR/Z Bmin = root:SNS_Settings:Bmin
    NVAR/Z Bmax = root:SNS_Settings:Bmax
    NVAR/Z NB   = root:SNS_Settings:NB

    p.Bmin = NVAR_Exists(Bmin) ? Bmin : -0.5
    p.Bmax = NVAR_Exists(Bmax) ? Bmax :  0.5
    p.NB   = NVAR_Exists(NB)   ? NB   : 1001

    // ------------------------------
    // Ray tracing / interface transparency
    // ------------------------------
    NVAR/Z BTK_barrier = root:SNS_Settings:BTK_barrier
    NVAR/Z legacyBTK_barrier = root:v_cfg_BTK_barrier

    if (NVAR_Exists(BTK_barrier))
        p.BTK_barrier = BTK_barrier
    elseif (NVAR_Exists(legacyBTK_barrier))
        p.BTK_barrier = legacyBTK_barrier
    else
        p.BTK_barrier = 0.05
    endif
    
    // ------------------------------
    // Band model
    // ------------------------------
    NVAR/Z SNS_bandModel  = root:SNS_Settings:SNS_bandModel
    NVAR/Z SNS_hbarvD_eVA = root:SNS_Settings:SNS_hbarvD_eVA
    NVAR/Z SNS_muDirac_eV = root:SNS_Settings:SNS_muDirac_eV

    p.SNS_bandModel  = NVAR_Exists(SNS_bandModel)  ? SNS_bandModel  : 0
    p.SNS_hbarvD_eVA = NVAR_Exists(SNS_hbarvD_eVA) ? SNS_hbarvD_eVA : 3.0
    p.SNS_muDirac_eV = NVAR_Exists(SNS_muDirac_eV) ? SNS_muDirac_eV : 0.2
    
    // ------------------------------
    // Diffusive / Usadel transport settings
    // ------------------------------
    NVAR/Z lmfp_bulk_nm = root:SNS_Settings:lmfp_bulk_nm
    NVAR/Z lmfp_edge_prefactor = root:SNS_Settings:lmfp_edge_prefactor

    p.lmfp_bulk_nm = NVAR_Exists(lmfp_bulk_nm) ? lmfp_bulk_nm : NaN
    p.lmfp_edge_prefactor = NVAR_Exists(lmfp_edge_prefactor) ? lmfp_edge_prefactor : pi

    // ------------------------------
    // Default impurity diagnostics
    // ------------------------------
    NVAR/Z impurity_density_default_nm2 = root:SNS_Settings:impurity_density_default_nm2
    NVAR/Z sigma_tr_default_nm = root:SNS_Settings:sigma_tr_default_nm
    NVAR/Z lmfp_bulk_default_nm = root:SNS_Settings:lmfp_bulk_default_nm

    p.impurity_density_default_nm2 = NVAR_Exists(impurity_density_default_nm2) ? impurity_density_default_nm2 : NaN
    p.sigma_tr_default_nm = NVAR_Exists(sigma_tr_default_nm) ? sigma_tr_default_nm : NaN
    p.lmfp_bulk_default_nm = NVAR_Exists(lmfp_bulk_default_nm) ? lmfp_bulk_default_nm : NaN

    // ------------------------------
    // Vortex / delta-map settings
    // ------------------------------
    NVAR/Z SNS_useVortex   = root:SNS_Settings:SNS_useVortex
    NVAR/Z SNS_nFlux       = root:SNS_Settings:SNS_nFlux
    NVAR/Z SNS_nIntSteps   = root:SNS_Settings:SNS_nIntSteps
    NVAR/Z SNS_useDeltamap = root:SNS_Settings:SNS_useDeltamap

    p.SNS_useVortex   = NVAR_Exists(SNS_useVortex)   ? SNS_useVortex   : 0
    p.SNS_nFlux       = NVAR_Exists(SNS_nFlux)       ? SNS_nFlux       : 0
    p.SNS_nIntSteps   = NVAR_Exists(SNS_nIntSteps)   ? SNS_nIntSteps   : 64
    p.SNS_useDeltamap = NVAR_Exists(SNS_useDeltamap) ? SNS_useDeltamap : 0

    // ------------------------------
    // Gamma broadening settings
    // ------------------------------
    NVAR/Z SNS_GammaBase_eV      = root:SNS_Settings:SNS_GammaBase_eV
    NVAR/Z SNS_useGammaDopp      = root:SNS_Settings:SNS_useGammaDopp
    NVAR/Z SNS_useGammaZeeman    = root:SNS_Settings:SNS_useGammaZeeman
    NVAR/Z SNS_useGammaPair      = root:SNS_Settings:SNS_useGammaPair
    NVAR/Z SNS_useGammaUser      = root:SNS_Settings:SNS_useGammaUser
    NVAR/Z SNS_useGammaCoh       = root:SNS_Settings:SNS_useGammaCoh
    NVAR/Z SNS_GammaDoppAlpha    = root:SNS_Settings:SNS_GammaDoppAlpha
    NVAR/Z SNS_GammaZeemanScale  = root:SNS_Settings:SNS_GammaZeemanScale
    NVAR/Z SNS_GammaPairScale    = root:SNS_Settings:SNS_GammaPairScale
    NVAR/Z SNS_GammaPair_Bc_T    = root:SNS_Settings:SNS_GammaPair_Bc_T

    p.SNS_GammaBase_eV      = NVAR_Exists(SNS_GammaBase_eV)     ? SNS_GammaBase_eV     : NaN
    p.SNS_useGammaDopp      = NVAR_Exists(SNS_useGammaDopp)     ? SNS_useGammaDopp     : 0
    p.SNS_useGammaZeeman    = NVAR_Exists(SNS_useGammaZeeman)   ? SNS_useGammaZeeman   : 0
    p.SNS_useGammaPair      = NVAR_Exists(SNS_useGammaPair)     ? SNS_useGammaPair     : 1
    p.SNS_useGammaUser      = NVAR_Exists(SNS_useGammaUser)     ? SNS_useGammaUser     : 0
    p.SNS_useGammaCoh       = NVAR_Exists(SNS_useGammaCoh)      ? SNS_useGammaCoh      : 0
    p.SNS_GammaDoppAlpha    = NVAR_Exists(SNS_GammaDoppAlpha)   ? SNS_GammaDoppAlpha   : 0.01
    p.SNS_GammaZeemanScale  = NVAR_Exists(SNS_GammaZeemanScale) ? SNS_GammaZeemanScale : 0
    p.SNS_GammaPairScale    = NVAR_Exists(SNS_GammaPairScale)   ? SNS_GammaPairScale   : 0.5 * Delta
    p.SNS_GammaPair_Bc_T    = NVAR_Exists(SNS_GammaPair_Bc_T)   ? SNS_GammaPair_Bc_T   : 0.5
End


//============================================================
// 1. GENERIC HELPERS
//============================================================

// Safe clamp

//==============================================================================
// Clamp
//
// Purpose:
//   Utility helper.
//
// Inputs:
//   val : input
//   lo : input
//   hi : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function Clamp(val, lo, hi)
    Variable val, lo, hi
    if (val < lo)
        return lo
    elseif (val > hi)
        return hi
    endif
    return val
End

// acos with argument clamped to [-1,1]

//==============================================================================
// AcosSafe
//
// Purpose:
//   Utility helper.
//
// Inputs:
//   u : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function AcosSafe(u)
    Variable u
    u = Clamp(u, -1, 1)
    return acos(u)
End

// Simple unwrap helper if you still need it for diagnostics

//==============================================================================
// UnwrapToPrev
//
// Purpose:
//   Internal helper.
//
// Inputs:
//   a : input
//   aPrev : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function UnwrapToPrev(a, aPrev)
    Variable a, aPrev

    Variable da = a - aPrev
    if (da > pi)
        a -= 2*pi
    elseif (da < -pi)
        a += 2*pi
    endif

    return a
End


// Plot all eigenvalue branches E_allBranches_full vs B_T

//==============================================================================
// SNS_Find_Boffset_From_ZBC
//
// Purpose:
//   Locate a feature or parameter from input data.
//
// Inputs:
//   expEB : input
//
// Outputs:
//   return : numeric
//
//==============================================================================
Function SNS_Find_Boffset_From_ZBC(expEB)
    Wave expEB

    Variable epsB = 25       // B units (e.g. mT) window around B=0
    Variable epsE = .05      // meV fit window half-width around the peak

    Variable E0 = 0          // default safe return

    if (DimSize(expEB,0) <= 1 || DimSize(expEB,1) <= 1)
        return 0
    endif

    //--------------------------------------------------------
    // 1) Sum over B in [-epsB, +epsB] -> Spec(E)
    //--------------------------------------------------------
    Duplicate/O/R=(-epsE,epsE)() expEB, EEwin
    SumDimension/D=0/DEST=Spec EEwin
    KillWaves/Z EEwin

    SetScale/P x, DimOffset(expEB,1), DimDelta(expEB,1), WaveUnits(expEB,1), Spec

    //--------------------------------------------------------
    // 2) Find dip location
    //--------------------------------------------------------
    WaveStats/Q Spec
    if (numtype(V_minloc) != 0)
        KillWaves/Z Spec
        return 0
    endif

    Duplicate/O/R=(V_minloc-epsB, V_minloc+epsB) Spec, SpecFit
    KillWaves/Z Spec

    if (DimSize(SpecFit,0) < 5)
        KillWaves/Z SpecFit
        return 0
    endif

    //--------------------------------------------------------
    // 3) Gaussian fit
    //    Igor "gauss" coefficients are typically:
    //      W_coef[0]=y0, W_coef[1]=A, W_coef[2]=x0, W_coef[3]=w
    //--------------------------------------------------------
    Duplicate/O SpecFit, SpecFit_cf

	// --- Constraints: enforce a DIP (negative amplitude) ---
	Make/T/O/N=2 T_Constraints
	T_Constraints[0] = "K1 < 0"   // amplitude must be negative
	T_Constraints[1] = "K3 > 0"   // width must be positive
	
	CurveFit/Q gauss, SpecFit_cf /C=T_Constraints
	Wave W_coef
	
	KillWaves/Z T_Constraints

    CurveFit/Q gauss, SpecFit_cf
    Wave W_coef

    if (!WaveExists(W_coef))
        KillWaves/Z SpecFit, SpecFit_cf
        return 0
    endif
    if (DimSize(W_coef,0) < 3)
        KillWaves/Z SpecFit, SpecFit_cf, W_coef
        return 0
    endif

    E0 = W_coef[2]
    if (numtype(E0) != 0)
        KillWaves/Z SpecFit, SpecFit_cf, W_coef
        return 0
    endif

    //--------------------------------------------------------
    // cleanup
    //--------------------------------------------------------
    KillWaves/Z SpecFit, SpecFit_cf, W_coef

    return E0
End
