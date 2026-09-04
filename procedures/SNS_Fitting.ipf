#pragma rtGlobals=3
#pragma IgorVersion=9.00

//==============================================================================
// SNS_Fitting.ipf
//
// Fit experimental LDOS data with ray-trace-based SNS spectra.
//
// Responsibilities:
//   - fit symmetrized LDOS(E) at a selected magnetic field with a mixed
//     2D surface + 3D bulk ray-trace DOS model
//   - fit ZBC(B) using the same mixed-DOS construction and fixed Stage-1
//     LDOS parameters
//   - copy and validate ray-trace channel waves from exported ray-trace folders
//   - fit large-bias surface-state onsets and average accepted line spectra
//   - provide the bulk-background and broadened-step fit functions used by the
//     bulk-versus-surface-state conductance decomposition
//
// Dependencies:
//   SNS_Core.ipf
//     - SNS_InitDefaultSettings
//     - root:SNS_Settings globals
//
//   SNS_DOS.ipf
//     - SNS_ComputeDOS_FromChannels
//
//   SNS_Broadening.ipf
//     - SNS_ApplyDOS_Broadening_TplusMod
//
// Notes:
//   Ray-trace folders are expected to use the canonical nm wave schema:
//     L_N_List_nm, W_eff_List_nm, Hit*List_nm for 2D
//     L_N_List_3D_nm, W_eff_List_3D_nm, Hit*List_3D_nm for 3D
//
//   Unit conversion to solver units is handled downstream by
//   SNS_ComputeDOS_FromChannels. This file should keep ray-trace geometry in nm.
//==============================================================================


//==============================================================================
// SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders
//
// Purpose:
//   Fit normalized experimental LDOS(E) with a mixed 2D surface + 3D bulk
//   ray-trace DOS spectrum at one selected magnetic field.
//
// Inputs:
//   LDOS_EB              : experimental LDOS wave, dimensions E x B
//   rayTraceFolder2D     : folder containing 2D ray-trace waves [nm schema]
//   rayTraceFolder3D     : folder containing 3D ray-trace waves [nm schema]
//   ratio2D0, ratio3D0   : starting relative 2D/3D spectral weights
//   DeltaBulk0_meV       : starting 3D/bulk gap [meV]
//   DeltaSurface0_meV    : starting 2D/surface gap [meV]
//
// Optional inputs:
//   B_mT or BIndex       : selected experimental field column
//   outPrefix            : output-wave prefix
//   deltaTol             : relative fit bound around starting gaps
//   ratioScaleTol        : relative fit bound around 2D/3D ratio scales
//   fitAbsMax_mV         : symmetric energy fit window half-width [mV]
//   xV_nm, yV_nm         : vortex position used by legacy endpoint phase [nm]
//   Broad0_meV           : starting intrinsic Lorentzian broadening [meV]
//   BroadMin/Max_meV     : broadening fit bounds [meV]
//   weightLowLDOS        : use low-LDOS weighted least squares
//   weightFloorRel       : relative floor for low-LDOS weights
//   doDisplay            : display fit graph
//
// Outputs:
//   Creates fit waves and scalar results using outPrefix.
//
// Fit parameters:
//   K0 : amplitude
//   K1 : Delta_bulk_fit_meV
//   K2 : Delta_surface_fit_meV
//   K3 : 2D ratio scale
//   K4 : 3D ratio scale
//   K5 : intrinsic Lorentzian broadening [meV]
//
// Notes:
//   The selected experimental field is also used for the model DOS.
//   Thus the W_eff / phase-bias contribution in SNS_ComputeDOS_FromChannels()
//   is active at nonzero B.
//==============================================================================
Function SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders(LDOS_EB, rayTraceFolder2D, rayTraceFolder3D, ratio2D0, ratio3D0, DeltaBulk0_meV, DeltaSurface0_meV, [B_mT, BIndex, outPrefix, deltaTol, ratioScaleTol, fitAbsMax_mV, xV_nm, yV_nm, Broad0_meV, BroadMin_meV, BroadMax_meV, weightLowLDOS, weightFloorRel, doDisplay])
	Wave LDOS_EB
	String rayTraceFolder2D, rayTraceFolder3D
	Variable ratio2D0, ratio3D0
	Variable DeltaBulk0_meV, DeltaSurface0_meV
	Variable B_mT, BIndex
	String outPrefix
	Variable deltaTol, ratioScaleTol, fitAbsMax_mV
	Variable xV_nm, yV_nm
	Variable Broad0_meV, BroadMin_meV, BroadMax_meV
	Variable weightLowLDOS, weightFloorRel
	Variable doDisplay

	String oldDF = GetDataFolder(1)

	SNS__AssertUniqueOptionalNames( \
		"B_mT;BIndex;outPrefix;deltaTol;ratioScaleTol;fitAbsMax_mV;" + \
		"xV_nm;yV_nm;Broad0_meV;BroadMin_meV;BroadMax_meV;" + \
		"weightLowLDOS;weightFloorRel;doDisplay")

	if (!DataFolderExists("root:SNS_Settings:"))
		SNS_InitDefaultSettings()
	endif

	NVAR settingsBroadening = root:SNS_Settings:Broadening

	if (ParamIsDefault(outPrefix))
		outPrefix = "Fit2D3D_RayTraceDOS"
	endif
	if (ParamIsDefault(deltaTol))
		deltaTol = 0.20
	endif
	if (ParamIsDefault(ratioScaleTol))
		ratioScaleTol = 0.25
	endif
	if (ParamIsDefault(fitAbsMax_mV))
		fitAbsMax_mV = 2.5
	endif
	if (ParamIsDefault(xV_nm))
		xV_nm = 0
	endif
	if (ParamIsDefault(yV_nm))
		yV_nm = 0
	endif
	if (ParamIsDefault(Broad0_meV))
		Broad0_meV = settingsBroadening * 1e3
	endif
	if (ParamIsDefault(BroadMin_meV))
		BroadMin_meV = max(1e-6, 0.25 * Broad0_meV)
	endif
	if (ParamIsDefault(BroadMax_meV))
		BroadMax_meV = 4 * Broad0_meV
	endif
	if (ParamIsDefault(weightLowLDOS))
		weightLowLDOS = 1
	endif
	if (ParamIsDefault(weightFloorRel))
		weightFloorRel = 0.05
	endif
	if (ParamIsDefault(doDisplay))
		doDisplay = 1
	endif

	if (numtype(DeltaBulk0_meV) != 0 || DeltaBulk0_meV <= 0)
		SetDataFolder $oldDF
		Abort "SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders: invalid DeltaBulk0_meV."
	endif
	if (numtype(DeltaSurface0_meV) != 0 || DeltaSurface0_meV <= 0)
		SetDataFolder $oldDF
		Abort "SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders: invalid DeltaSurface0_meV."
	endif
	if (DeltaBulk0_meV <= DeltaSurface0_meV)
		SetDataFolder $oldDF
		Abort "SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders: starting values must satisfy DeltaBulk0_meV > DeltaSurface0_meV."
	endif
	if (numtype(ratio2D0) != 0 || ratio2D0 < 0 || numtype(ratio3D0) != 0 || ratio3D0 < 0 || ratio2D0 + ratio3D0 <= 0)
		SetDataFolder $oldDF
		Abort "SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders: invalid 2D/3D ratio."
	endif
	if (numtype(deltaTol) != 0 || deltaTol <= 0 || deltaTol >= 1)
		SetDataFolder $oldDF
		Abort "SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders: invalid deltaTol."
	endif
	if (numtype(ratioScaleTol) != 0 || ratioScaleTol < 0 || ratioScaleTol >= 1)
		SetDataFolder $oldDF
		Abort "SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders: invalid ratioScaleTol."
	endif
	if (numtype(Broad0_meV) != 0 || Broad0_meV <= 0 || numtype(BroadMin_meV) != 0 || BroadMin_meV <= 0 || numtype(BroadMax_meV) != 0 || BroadMax_meV <= BroadMin_meV)
		SetDataFolder $oldDF
		Abort "SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders: invalid broadening start/bounds."
	endif
	if (numtype(weightFloorRel) != 0 || weightFloorRel <= 0)
		SetDataFolder $oldDF
		Abort "SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders: invalid weightFloorRel."
	endif

	// -------------------------------------------------------------------------
	// Setup fit folder and copy ray-trace waves.
	// -------------------------------------------------------------------------
	NewDataFolder/O root:SNS_RayTraceDOSFit

	SNS__CopyRayTraceDOSFitWavesFromFolder(rayTraceFolder2D, 0)
	SNS__CopyRayTraceDOSFitWavesFromFolder(rayTraceFolder3D, 1)

	SetDataFolder root:SNS_RayTraceDOSFit
	Variable/G g_ratio2D0 = ratio2D0
	Variable/G g_ratio3D0 = ratio3D0
	Variable/G g_xV_nm = xV_nm
	Variable/G g_yV_nm = yV_nm
	SetDataFolder $oldDF

	// -------------------------------------------------------------------------
	// Select experimental B column.
	// -------------------------------------------------------------------------
	Variable nEnergyRaw = DimSize(LDOS_EB, 0)
	Variable nFieldRaw = DimSize(LDOS_EB, 1)

	if (nEnergyRaw < 2 || nFieldRaw < 1)
		SetDataFolder $oldDF
		Abort "SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders: LDOS_EB must be E x B."
	endif

	Variable bAxisOffset = DimOffset(LDOS_EB, 1)
	Variable bAxisDelta = DimDelta(LDOS_EB, 1)
	String bAxisUnit = WaveUnits(LDOS_EB, 1)

	Variable bScaleToMilliT = 1
	if (StringMatch(bAxisUnit, "*T*") && !StringMatch(bAxisUnit, "*mT*"))
		bScaleToMilliT = 1000
	endif

	Variable selectedBIndex, ind
	Variable bestBdiff, bHereMilliT

	if (!ParamIsDefault(BIndex))
		selectedBIndex = round(BIndex)
	elseif (!ParamIsDefault(B_mT))
		bestBdiff = Inf
		selectedBIndex = 0
		for (ind = 0; ind < nFieldRaw; ind += 1)
			bHereMilliT = (bAxisOffset + ind*bAxisDelta) * bScaleToMilliT
			if (abs(bHereMilliT - B_mT) < bestBdiff)
				bestBdiff = abs(bHereMilliT - B_mT)
				selectedBIndex = ind
			endif
		endfor
	else
		bestBdiff = Inf
		selectedBIndex = 0
		for (ind = 0; ind < nFieldRaw; ind += 1)
			bHereMilliT = (bAxisOffset + ind*bAxisDelta) * bScaleToMilliT
			if (abs(bHereMilliT) < bestBdiff)
				bestBdiff = abs(bHereMilliT)
				selectedBIndex = ind
			endif
		endfor
	endif

	if (selectedBIndex < 0 || selectedBIndex >= nFieldRaw)
		SetDataFolder $oldDF
		Abort "SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders: selected B index out of range."
	endif

	Variable selectedB_mT = (bAxisOffset + selectedBIndex*bAxisDelta) * bScaleToMilliT

	// This is the field passed to SNS_ComputeDOS_FromChannels().
	SetDataFolder root:SNS_RayTraceDOSFit
	Make/O/D/N=1 bAxisFitT_fit
	bAxisFitT_fit[0] = selectedB_mT * 1e-3
	Variable/G g_selectedB_mT = selectedB_mT
	SetDataFolder $oldDF

	// -------------------------------------------------------------------------
	// Build symmetrized experimental fit wave.
	// -------------------------------------------------------------------------
	Variable eAxisOffset = DimOffset(LDOS_EB, 0)
	Variable eAxisDelta = DimDelta(LDOS_EB, 0)
	String eAxisUnit = WaveUnits(LDOS_EB, 0)

	Variable eScaleToMilliV = 1
	if (StringMatch(eAxisUnit, "*eV*") && !StringMatch(eAxisUnit, "*meV*") && !StringMatch(eAxisUnit, "*mV*"))
		eScaleToMilliV = 1000
	endif

	Make/O/D/N=(nEnergyRaw) root:SNS_RayTraceDOSFit:tmpErawMilliV_fit
	Make/O/D/N=(nEnergyRaw) root:SNS_RayTraceDOSFit:tmpYraw_fit

	Wave eRawMilliVWave = root:SNS_RayTraceDOSFit:tmpErawMilliV_fit
	Wave yRawWave = root:SNS_RayTraceDOSFit:tmpYraw_fit

	for (ind = 0; ind < nEnergyRaw; ind += 1)
		eRawMilliVWave[ind] = (eAxisOffset + ind*eAxisDelta) * eScaleToMilliV
		yRawWave[ind] = LDOS_EB[ind][selectedBIndex]
	endfor

	Variable eFirstMilliV = eRawMilliVWave[0]
	Variable eLastMilliV = eRawMilliVWave[nEnergyRaw-1]
	Variable eMinMilliV = min(eFirstMilliV, eLastMilliV)
	Variable eMaxMilliV = max(eFirstMilliV, eLastMilliV)

	Variable nFitPoints = 0
	Variable eHereMilliV

	for (ind = 0; ind < nEnergyRaw; ind += 1)
		eHereMilliV = eRawMilliVWave[ind]
		if (abs(eHereMilliV) <= fitAbsMax_mV && -eHereMilliV >= eMinMilliV && -eHereMilliV <= eMaxMilliV)
			nFitPoints += 1
		endif
	endfor

	if (nFitPoints < 5)
		SetDataFolder $oldDF
		Abort "SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders: too few fit points."
	endif

	Make/O/D/N=(nFitPoints) $(outPrefix + "_Efit_mV")
	Make/O/D/N=(nFitPoints) $(outPrefix + "_Yfit_sym")

	Wave eFitMilliVWave = $(outPrefix + "_Efit_mV")
	Wave yFitSymWave = $(outPrefix + "_Yfit_sym")

	Variable fitInd = 0
	for (ind = 0; ind < nEnergyRaw; ind += 1)
		eHereMilliV = eRawMilliVWave[ind]
		if (abs(eHereMilliV) <= fitAbsMax_mV && -eHereMilliV >= eMinMilliV && -eHereMilliV <= eMaxMilliV)
			eFitMilliVWave[fitInd] = eHereMilliV
			yFitSymWave[fitInd] = 0.5 * (interp(eHereMilliV, eRawMilliVWave, yRawWave) + interp(-eHereMilliV, eRawMilliVWave, yRawWave))
			fitInd += 1
		endif
	endfor

	KillWaves/Z eRawMilliVWave, yRawWave

	SetScale/P x, eFitMilliVWave[0], eFitMilliVWave[1]-eFitMilliVWave[0], "mV", eFitMilliVWave
	SetScale/P x, eFitMilliVWave[0], eFitMilliVWave[1]-eFitMilliVWave[0], "", yFitSymWave

	// -------------------------------------------------------------------------
	// Fit coefficients, constraints, and weights.
	// -------------------------------------------------------------------------
	Make/O/D/N=6 $(outPrefix + "_coef")
	Wave coefWave = $(outPrefix + "_coef")

	coefWave[0] = 1
	coefWave[1] = DeltaBulk0_meV
	coefWave[2] = DeltaSurface0_meV
	coefWave[3] = 1
	coefWave[4] = 1
	coefWave[5] = Broad0_meV

	Variable scaleLow = max(1e-6, 1 - ratioScaleTol)
	Variable scaleHigh = 1 + ratioScaleTol

	Make/O/T/N=12 $(outPrefix + "_constraints")
	Wave/T constraintTextWave = $(outPrefix + "_constraints")

	constraintTextWave[0] = "K0 > 0"
	constraintTextWave[1] = "K1 > " + num2str((1-deltaTol)*DeltaBulk0_meV)
	constraintTextWave[2] = "K1 < " + num2str((1+deltaTol)*DeltaBulk0_meV)
	constraintTextWave[3] = "K2 > " + num2str((1-deltaTol)*DeltaSurface0_meV)
	constraintTextWave[4] = "K2 < " + num2str((1+deltaTol)*DeltaSurface0_meV)
	constraintTextWave[5] = "K1 - K2 > 0"
	constraintTextWave[6] = "K3 > " + num2str(scaleLow)
	constraintTextWave[7] = "K3 < " + num2str(scaleHigh)
	constraintTextWave[8] = "K4 > " + num2str(scaleLow)
	constraintTextWave[9] = "K4 < " + num2str(scaleHigh)
	constraintTextWave[10] = "K5 > " + num2str(BroadMin_meV)
	constraintTextWave[11] = "K5 < " + num2str(BroadMax_meV)

	Make/O/D/N=(nFitPoints) $(outPrefix + "_fitWeights")
	Wave fitWeightWave = $(outPrefix + "_fitWeights")

	Variable weightFloor
	if (weightLowLDOS)
		WaveStats/Q yFitSymWave
		weightFloor = weightFloorRel * max(abs(V_max), abs(V_min))
		weightFloor = max(weightFloor, 1e-12)

		fitWeightWave = 1 / ((abs(yFitSymWave[p]) + weightFloor)^2)
	else
		weightFloor = NaN
		fitWeightWave = 1
	endif

	FuncFit/Q SNS__RayTrace2D3D_DOSFitFunc coefWave yFitSymWave /X=eFitMilliVWave /W=fitWeightWave /C=constraintTextWave

	Make/O/D/N=(nFitPoints) $(outPrefix + "_fit")
	Make/O/D/N=(nFitPoints) $(outPrefix + "_fit_2D_surface_raw")
	Make/O/D/N=(nFitPoints) $(outPrefix + "_fit_3D_bulk_raw")

	Wave yModelWave = $(outPrefix + "_fit")
	Wave ySurface2DWave = $(outPrefix + "_fit_2D_surface_raw")
	Wave yBulk3DWave = $(outPrefix + "_fit_3D_bulk_raw")

	SNS__RayTrace2D3D_BuildComponents(coefWave, eFitMilliVWave, ySurface2DWave, yBulk3DWave, yModelWave)

	Variable ampFit = coefWave[0]
	Variable deltaBulkFit_meV = coefWave[1]
	Variable deltaSurfaceFit_meV = coefWave[2]
	Variable scale2DFit = coefWave[3]
	Variable scale3DFit = coefWave[4]
	Variable broadFit_meV = coefWave[5]

	Variable rawWeight2D = scale2DFit * ratio2D0
	Variable rawWeight3D = scale3DFit * ratio3D0
	Variable rawWeightSum = rawWeight2D + rawWeight3D
	Variable frac2DFit = rawWeight2D/rawWeightSum
	Variable frac3DFit = rawWeight3D/rawWeightSum

	Variable/G $(outPrefix + "_B_mT") = selectedB_mT
	Variable/G $(outPrefix + "_B_index") = selectedBIndex
	Variable/G $(outPrefix + "_AmplitudeFit") = ampFit
	Variable/G $(outPrefix + "_DeltaBulkFit_meV") = deltaBulkFit_meV
	Variable/G $(outPrefix + "_DeltaSurfaceFit_meV") = deltaSurfaceFit_meV
	Variable/G $(outPrefix + "_BroadeningFit_meV") = broadFit_meV
	Variable/G $(outPrefix + "_BroadeningMin_meV") = BroadMin_meV
	Variable/G $(outPrefix + "_BroadeningMax_meV") = BroadMax_meV
	Variable/G $(outPrefix + "_s2D_fit") = scale2DFit
	Variable/G $(outPrefix + "_s3D_fit") = scale3DFit
	Variable/G $(outPrefix + "_ratio2D0") = ratio2D0
	Variable/G $(outPrefix + "_ratio3D0") = ratio3D0
	Variable/G $(outPrefix + "_frac2D_fit") = frac2DFit
	Variable/G $(outPrefix + "_frac3D_fit") = frac3DFit
	Variable/G $(outPrefix + "_fitAbsMax_mV") = fitAbsMax_mV
	Variable/G $(outPrefix + "_deltaTol") = deltaTol
	Variable/G $(outPrefix + "_ratioScaleTol") = ratioScaleTol
	Variable/G $(outPrefix + "_weightLowLDOS") = weightLowLDOS
	Variable/G $(outPrefix + "_weightFloorRel") = weightFloorRel
	Variable/G $(outPrefix + "_weightFloor") = weightFloor

	if (doDisplay)
		Display/K=1 yFitSymWave vs eFitMilliVWave
		String graphName = WinName(0, 1, 1)

		AppendToGraph/W=$graphName yModelWave vs eFitMilliVWave
		AppendToGraph/W=$graphName ySurface2DWave vs eFitMilliVWave
		AppendToGraph/W=$graphName yBulk3DWave vs eFitMilliVWave

		String yFitTraceName = NameOfWave(yFitSymWave)
		String yModelTraceName = NameOfWave(yModelWave)
		String ySurfaceTraceName = NameOfWave(ySurface2DWave)
		String yBulkTraceName = NameOfWave(yBulk3DWave)

		ModifyGraph/W=$graphName mode($yFitTraceName)=3
		ModifyGraph/W=$graphName marker($yFitTraceName)=8
		ModifyGraph/W=$graphName rgb($yFitTraceName)=(0,0,0)

		ModifyGraph/W=$graphName rgb($yModelTraceName)=(54741,24158,0)
		ModifyGraph/W=$graphName lsize($yModelTraceName)=2

		ModifyGraph/W=$graphName rgb($ySurfaceTraceName)=(0,0,65535)
		ModifyGraph/W=$graphName lstyle($ySurfaceTraceName)=2

		ModifyGraph/W=$graphName rgb($yBulkTraceName)=(0,40000,0)
		ModifyGraph/W=$graphName lstyle($yBulkTraceName)=2

		ModifyGraph/W=$graphName tick=2, mirror=2, standoff=0
		Label/W=$graphName bottom "E (mV)"
		Label/W=$graphName left "normalized LDOS"

		String txt
		txt = "2D(surface) + 3D(bulk) ray-trace DOS fit\r"
		txt += "weighted least squares, raw mix then resolution convolution\r"
		txt += "B = " + num2str(round(10*selectedB_mT)/10) + " mT\r"
		txt += "A = " + num2str(round(1000*ampFit)/1000) + "\r"
		txt += "Delta_bulk = " + num2str(round(1000*deltaBulkFit_meV)/1000) + " meV\r"
		txt += "Delta_surface = " + num2str(round(1000*deltaSurfaceFit_meV)/1000) + " meV\r"
		txt += "Broadening = " + num2str(round(1000*broadFit_meV)/1000) + " meV\r"
		txt += "s2D = " + num2str(round(1000*scale2DFit)/1000) + ", s3D = " + num2str(round(1000*scale3DFit)/1000) + "\r"
		txt += "f2D = " + num2str(round(1000*frac2DFit)/1000) + ", f3D = " + num2str(round(1000*frac3DFit)/1000) + "\r"
		txt += "weightLowLDOS = " + num2str(weightLowLDOS) + ", floorRel = " + num2str(weightFloorRel)

		TextBox/W=$graphName/C/N=RayTraceDOSFitInfo/F=0/A=RT/X=0/Y=0/B=(65535,65535,65535,39321) txt
	endif

	SetDataFolder $oldDF

	return deltaSurfaceFit_meV
End


//==============================================================================
// SNS__RayTrace2D3D_DOSFitFunc
//
// Purpose:
//   Igor FuncFit callback for the LDOS(E) mixed 2D+3D DOS model.
//
// Inputs:
//   coefWave    : fit coefficient wave [K0..K5]
//   xMilliVWave : fit energy axis [mV]
//
// Outputs:
//   yOutWave : model LDOS on xMilliVWave
//
// Notes:
//   The actual model construction is delegated to
//   SNS__RayTrace2D3D_BuildComponents.
//==============================================================================
Function SNS__RayTrace2D3D_DOSFitFunc(coefWave, yOutWave, xMilliVWave) : FitFunc
	Wave coefWave
	Wave yOutWave
	Wave xMilliVWave

	Make/O/D/N=(numpnts(xMilliVWave)) root:SNS_RayTraceDOSFit:tmpSurfaceRaw2D_fit
	Make/O/D/N=(numpnts(xMilliVWave)) root:SNS_RayTraceDOSFit:tmpBulkRaw3D_fit

	Wave ySurfaceRaw2DWave = root:SNS_RayTraceDOSFit:tmpSurfaceRaw2D_fit
	Wave yBulkRaw3DWave = root:SNS_RayTraceDOSFit:tmpBulkRaw3D_fit

	SNS__RayTrace2D3D_BuildComponents(coefWave, xMilliVWave, ySurfaceRaw2DWave, yBulkRaw3DWave, yOutWave)

	return 0
End


//==============================================================================
// SNS__RayTrace2D3D_BuildComponents
//
// Purpose:
//   Build the raw 2D surface DOS, raw 3D bulk DOS, and mixed broadened LDOS
//   for a supplied LDOS-fit coefficient wave.
//
// Inputs:
//   coefWave          : [A, Delta_bulk_meV, Delta_surface_meV,
//                       scale2D, scale3D, Broadening_meV]
//   xMilliVWave       : requested output energy axis [mV]
//
// Outputs:
//   ySurface2DOutWave : raw 2D surface component on xMilliVWave
//   yBulk3DOutWave    : raw 3D bulk component on xMilliVWave
//   yMixedBroadOutWave: mixed 2D+3D LDOS after thermal/modulation convolution
//
// Reads/writes:
//   Temporarily modifies root:SNS_Settings:Delta and Broadening, then restores
//   them before return.
//
// Notes:
//   Ray-trace channel lengths and endpoint coordinates are stored in nm in the
//   fit working folder. SNS_ComputeDOS_FromChannels handles conversion to the
//   solver's meter-based variables internally.
//==============================================================================
Function SNS__RayTrace2D3D_BuildComponents(coefWave, xMilliVWave, ySurface2DOutWave, yBulk3DOutWave, yMixedBroadOutWave)
	Wave coefWave
	Wave xMilliVWave
	Wave ySurface2DOutWave, yBulk3DOutWave, yMixedBroadOutWave

	String oldDF = GetDataFolder(1)
	SetDataFolder root:SNS_RayTraceDOSFit

	NVAR settingsDelta = root:SNS_Settings:Delta
	NVAR settingsBroadening = root:SNS_Settings:Broadening

	Variable deltaOriginal_eV = settingsDelta
	Variable broadOriginal_eV = settingsBroadening

	NVAR ratio2D0 = root:SNS_RayTraceDOSFit:g_ratio2D0
	NVAR ratio3D0 = root:SNS_RayTraceDOSFit:g_ratio3D0
	NVAR xV_nm = root:SNS_RayTraceDOSFit:g_xV_nm
	NVAR yV_nm = root:SNS_RayTraceDOSFit:g_yV_nm

	Wave bAxisFitTWave = root:SNS_RayTraceDOSFit:bAxisFitT_fit

	Wave lenList2DWave = root:SNS_RayTraceDOSFit:lenList2D_fit
	Wave weffList2DWave = root:SNS_RayTraceDOSFit:weffList2D_fit
	Wave wchanList2DWave = root:SNS_RayTraceDOSFit:wchanList2D_fit
	Wave teffList2DWave = root:SNS_RayTraceDOSFit:teffList2D_fit
	Wave rs1xList2DWave = root:SNS_RayTraceDOSFit:rs1xList2D_fit
	Wave rs1yList2DWave = root:SNS_RayTraceDOSFit:rs1yList2D_fit
	Wave rs2xList2DWave = root:SNS_RayTraceDOSFit:rs2xList2D_fit
	Wave rs2yList2DWave = root:SNS_RayTraceDOSFit:rs2yList2D_fit

	Wave lenList3DWave = root:SNS_RayTraceDOSFit:lenList3D_fit
	Wave weffList3DWave = root:SNS_RayTraceDOSFit:weffList3D_fit
	Wave wchanList3DWave = root:SNS_RayTraceDOSFit:wchanList3D_fit
	Wave teffList3DWave = root:SNS_RayTraceDOSFit:teffList3D_fit
	Wave rs1xList3DWave = root:SNS_RayTraceDOSFit:rs1xList3D_fit
	Wave rs1yList3DWave = root:SNS_RayTraceDOSFit:rs1yList3D_fit
	Wave rs2xList3DWave = root:SNS_RayTraceDOSFit:rs2xList3D_fit
	Wave rs2yList3DWave = root:SNS_RayTraceDOSFit:rs2yList3D_fit

	Variable amp = coefWave[0]
	Variable deltaBulkFit_meV = coefWave[1]
	Variable deltaSurfaceFit_meV = coefWave[2]
	Variable scale2D = coefWave[3]
	Variable scale3D = coefWave[4]
	Variable broadFit_meV = coefWave[5]

	if (amp <= 0 || deltaBulkFit_meV <= 0 || deltaSurfaceFit_meV <= 0 || scale2D <= 0 || scale3D <= 0 || broadFit_meV <= 0)
		ySurface2DOutWave = 0
		yBulk3DOutWave = 0
		yMixedBroadOutWave = 0
		SetDataFolder $oldDF
		return 0
	endif

	if (deltaBulkFit_meV <= deltaSurfaceFit_meV)
		ySurface2DOutWave = 0
		yBulk3DOutWave = 0
		yMixedBroadOutWave = 0
		SetDataFolder $oldDF
		return 0
	endif

	Variable rawWeight2D = scale2D * ratio2D0
	Variable rawWeight3D = scale3D * ratio3D0
	Variable rawWeightSum = rawWeight2D + rawWeight3D

	if (rawWeightSum <= 0)
		ySurface2DOutWave = 0
		yBulk3DOutWave = 0
		yMixedBroadOutWave = 0
		SetDataFolder $oldDF
		return 0
	endif

	Variable frac2D = rawWeight2D/rawWeightSum
	Variable frac3D = rawWeight3D/rawWeightSum

	Variable deltaBulkFit_eV = deltaBulkFit_meV * 1e-3
	Variable deltaSurfaceFit_eV = deltaSurfaceFit_meV * 1e-3
	Variable broadFit_eV = broadFit_meV * 1e-3

	// 2D surface raw DOS at selected field.
	settingsDelta = deltaSurfaceFit_eV
	settingsBroadening = broadFit_eV

	SNS_ComputeDOS_FromChannels( \
		bAxisFitTWave, lenList2DWave, weffList2DWave, wchanList2DWave, teffList2DWave, \
		rs1xList2DWave, rs1yList2DWave, rs2xList2DWave, rs2yList2DWave, \
		xV_nm, yV_nm, 1, "tmp_DOS2D_surface", "tmp_E2D_surface", \
		is3D=0)

	// 3D bulk raw DOS at selected field.
	settingsDelta = deltaBulkFit_eV
	settingsBroadening = broadFit_eV

	SNS_ComputeDOS_FromChannels( \
		bAxisFitTWave, lenList3DWave, weffList3DWave, wchanList3DWave, teffList3DWave, \
		rs1xList3DWave, rs1yList3DWave, rs2xList3DWave, rs2yList3DWave, \
		xV_nm, yV_nm, 1, "tmp_DOS3D_bulk", "tmp_E3D_bulk", \
		is3D=1)

	Wave energyAxis2DWave = root:SNS_RayTraceDOSFit:tmp_E2D_surface
	Wave energyAxis3DWave = root:SNS_RayTraceDOSFit:tmp_E3D_bulk
	Wave dosRaw2DWave = root:SNS_RayTraceDOSFit:tmp_DOS2D_surface
	Wave dosRaw3DWave = root:SNS_RayTraceDOSFit:tmp_DOS3D_bulk

	Variable nEnergy2D = numpnts(energyAxis2DWave)
	Variable nEnergy3D = numpnts(energyAxis3DWave)

	if (nEnergy2D < 2 || nEnergy3D < 2)
		ySurface2DOutWave = 0
		yBulk3DOutWave = 0
		yMixedBroadOutWave = 0
		settingsDelta = deltaOriginal_eV
		settingsBroadening = broadOriginal_eV
		SetDataFolder $oldDF
		return 0
	endif

	Make/O/D/N=(nEnergy2D) root:SNS_RayTraceDOSFit:tmpDosRaw2D_1D
	Make/O/D/N=(nEnergy3D) root:SNS_RayTraceDOSFit:tmpDosRaw3D_1D

	Wave dosRaw2D_1DWave = root:SNS_RayTraceDOSFit:tmpDosRaw2D_1D
	Wave dosRaw3D_1DWave = root:SNS_RayTraceDOSFit:tmpDosRaw3D_1D

	dosRaw2D_1DWave = dosRaw2DWave[p][0]
	dosRaw3D_1DWave = dosRaw3DWave[p][0]

	Make/O/D/N=(nEnergy2D) root:SNS_RayTraceDOSFit:tmpSurface2D_raw_onModelGrid
	Make/O/D/N=(nEnergy2D) root:SNS_RayTraceDOSFit:tmpBulk3D_raw_onModelGrid

	Wave surfaceRawModelGridWave = root:SNS_RayTraceDOSFit:tmpSurface2D_raw_onModelGrid
	Wave bulkRawModelGridWave = root:SNS_RayTraceDOSFit:tmpBulk3D_raw_onModelGrid

	surfaceRawModelGridWave = amp * frac2D * dosRaw2D_1DWave[p]
	bulkRawModelGridWave = amp * frac3D * interp(energyAxis2DWave[p], energyAxis3DWave, dosRaw3D_1DWave)

	Make/O/D/N=(nEnergy2D,1) root:SNS_RayTraceDOSFit:tmpMixedDOS_raw_EB
	Wave mixedRawEBWave = root:SNS_RayTraceDOSFit:tmpMixedDOS_raw_EB

	SetScale/P x, energyAxis2DWave[0], energyAxis2DWave[1]-energyAxis2DWave[0], "eV", mixedRawEBWave
	SetScale/P y, bAxisFitTWave[0], 1, "T", mixedRawEBWave

	mixedRawEBWave[][0] = surfaceRawModelGridWave[p] + bulkRawModelGridWave[p]

	// Apply experimental-resolution convolution to the mixed LDOS.
	settingsBroadening = broadFit_eV
	SNS_ApplyDOS_Broadening_TplusMod(mixedRawEBWave, "root:SNS_RayTraceDOSFit:tmpMixedDOS_broadened_EB")
	Wave mixedBroadEBWave = root:SNS_RayTraceDOSFit:tmpMixedDOS_broadened_EB

	// Restore global settings immediately.
	settingsDelta = deltaOriginal_eV
	settingsBroadening = broadOriginal_eV

	Make/O/D/N=(nEnergy2D) root:SNS_RayTraceDOSFit:tmpMixedDOS_broadened_1D
	Wave mixedBroad1DWave = root:SNS_RayTraceDOSFit:tmpMixedDOS_broadened_1D

	mixedBroad1DWave = mixedBroadEBWave[p][0]

	ySurface2DOutWave = interp(xMilliVWave[p]*1e-3, energyAxis2DWave, surfaceRawModelGridWave)
	yBulk3DOutWave = interp(xMilliVWave[p]*1e-3, energyAxis2DWave, bulkRawModelGridWave)
	yMixedBroadOutWave = interp(xMilliVWave[p]*1e-3, energyAxis2DWave, mixedBroad1DWave)

	SetDataFolder $oldDF

	return 0
End


//==============================================================================
// SNS__CopyRayTraceDOSFitWavesFromFolder
//
// Purpose:
//   Copy the canonical ray-trace channel waves from a 2D or 3D ray-trace folder
//   into root:SNS_RayTraceDOSFit for use by the fitting callbacks.
//
// Inputs:
//   rtFolder : source ray-trace folder
//   is3D     : 0 for 2D schema, nonzero for 3D schema
//
// Outputs:
//   Creates internal fit waves:
//     lenList2D/3D_fit, weffList2D/3D_fit, wchanList2D/3D_fit,
//     teffList2D/3D_fit, rs1/rs2 endpoint waves.
//
// Notes:
//   Source folders must use canonical nm endpoint names:
//     2D: Hit1x_List_nm, ...
//     3D: Hit1x_List_3D_nm, ...
//   Missing endpoint waves are treated as errors, not silently replaced by
//   zeros, because endpoint phase terms would otherwise become unphysical.
//==============================================================================
Function SNS__CopyRayTraceDOSFitWavesFromFolder(rtFolder, is3D)
	String rtFolder
	Variable is3D

	String oldDF = GetDataFolder(1)

	if (!DataFolderExists(rtFolder + ":"))
		Abort "SNS__CopyRayTraceDOSFitWavesFromFolder: folder does not exist: " + rtFolder
	endif

	NewDataFolder/O root:SNS_RayTraceDOSFit

	String strLenName
	String strWeffName
	String strWchanName
	String strTeffName
	if (is3D)
		strLenName = rtFolder + ":L_N_List_3D_nm"
		strWeffName = rtFolder + ":W_eff_List_3D_nm"
		strWchanName = rtFolder + ":wChan_3D"
		strTeffName = rtFolder + ":T_eff_List_3D"
	else
		strLenName = rtFolder + ":L_N_List_nm"
		strWeffName = rtFolder + ":W_eff_List_nm"
		strWchanName = rtFolder + ":wChan"
		strTeffName = rtFolder + ":T_eff_List"
	endif

	Wave/Z srcLenWave = $strLenName
	Wave/Z srcWeffWave = $strWeffName
	Wave/Z srcWchanWave = $strWchanName
	Wave/Z srcTeffWave = $strTeffName

	if (!WaveExists(srcLenWave))
		Abort "Missing ray-trace wave: " + strLenName
	endif
	if (!WaveExists(srcWeffWave))
		Abort "Missing ray-trace wave: " + strWeffName
	endif
	if (!WaveExists(srcWchanWave))
		Abort "Missing ray-trace wave: " + strWchanName
	endif

	Variable nChannels = numpnts(srcLenWave)

	if (is3D)
		Duplicate/O srcLenWave, root:SNS_RayTraceDOSFit:lenList3D_fit
		Duplicate/O srcWeffWave, root:SNS_RayTraceDOSFit:weffList3D_fit
		Duplicate/O srcWchanWave, root:SNS_RayTraceDOSFit:wchanList3D_fit

		if (WaveExists(srcTeffWave))
			Duplicate/O srcTeffWave, root:SNS_RayTraceDOSFit:teffList3D_fit
		else
			Make/O/D/N=(nChannels) root:SNS_RayTraceDOSFit:teffList3D_fit
			Wave teffOut3DWave = root:SNS_RayTraceDOSFit:teffList3D_fit
			teffOut3DWave = 1
		endif

		SNS__CopyRequiredRayTraceWave(rtFolder + ":Hit1x_List_3D_nm", "root:SNS_RayTraceDOSFit:rs1xList3D_fit", nChannels)
		SNS__CopyRequiredRayTraceWave(rtFolder + ":Hit1y_List_3D_nm", "root:SNS_RayTraceDOSFit:rs1yList3D_fit", nChannels)
		SNS__CopyRequiredRayTraceWave(rtFolder + ":Hit2x_List_3D_nm", "root:SNS_RayTraceDOSFit:rs2xList3D_fit", nChannels)
		SNS__CopyRequiredRayTraceWave(rtFolder + ":Hit2y_List_3D_nm", "root:SNS_RayTraceDOSFit:rs2yList3D_fit", nChannels)
	else
		Duplicate/O srcLenWave, root:SNS_RayTraceDOSFit:lenList2D_fit
		Duplicate/O srcWeffWave, root:SNS_RayTraceDOSFit:weffList2D_fit
		Duplicate/O srcWchanWave, root:SNS_RayTraceDOSFit:wchanList2D_fit

		if (WaveExists(srcTeffWave))
			Duplicate/O srcTeffWave, root:SNS_RayTraceDOSFit:teffList2D_fit
		else
			Make/O/D/N=(nChannels) root:SNS_RayTraceDOSFit:teffList2D_fit
			Wave teffOut2DWave = root:SNS_RayTraceDOSFit:teffList2D_fit
			teffOut2DWave = 1
		endif

		SNS__CopyRequiredRayTraceWave(rtFolder + ":Hit1x_List_nm", "root:SNS_RayTraceDOSFit:rs1xList2D_fit", nChannels)
		SNS__CopyRequiredRayTraceWave(rtFolder + ":Hit1y_List_nm", "root:SNS_RayTraceDOSFit:rs1yList2D_fit", nChannels)
		SNS__CopyRequiredRayTraceWave(rtFolder + ":Hit2x_List_nm", "root:SNS_RayTraceDOSFit:rs2xList2D_fit", nChannels)
		SNS__CopyRequiredRayTraceWave(rtFolder + ":Hit2y_List_nm", "root:SNS_RayTraceDOSFit:rs2yList2D_fit", nChannels)
	endif

	SetDataFolder $oldDF
End


//==============================================================================
// SNS_Conductance_Dynes
//
// Thermally broadened Dynes conductance for zero-field spectroscopy fits.
// This is the SNS-prefixed, repository-local form of the legacy
// Conductance_Dynes fit function. Energies and kBT are in meV.
//
// Coefficients:
//   w[0] temperature [K]
//   w[1] gap [meV]
//   w[2] Dynes broadening [meV]
//   w[3] amplitude
//   w[4] bias offset [meV]
//==============================================================================
Function SNS_Conductance_Dynes(w, x) : FitFunc
	Wave w
	Variable x

	Variable kB_meV_per_K = 8.62e-2
	Variable e1 = x - 5*kB_meV_per_K*w[0]
	Variable e2 = x + 5*kB_meV_per_K*w[0]
	Variable nIntegration = 1000
	Variable eStep = 10*kB_meV_per_K*w[0]/nIntegration
	Variable ePoint = e1
	Variable sumOdd = 0
	Variable sumEven = 0
	Variable sumTotal
	Variable i

	// SNS_ThermalKernelScalar returns 1/eV; 1e-3 converts it to 1/meV for
	// integration on this fit function's meV energy axis.
	sumTotal = SNS_DynesDOS(e1-w[4], w[1], w[2])*1e-3*SNS_ThermalKernelScalar((e1-x)*1e-3, w[0], 8.62e-5)
	sumTotal += SNS_DynesDOS(e2-w[4], w[1], w[2])*1e-3*SNS_ThermalKernelScalar((e2-x)*1e-3, w[0], 8.62e-5)

	for (i = 1; i < nIntegration/2; i += 1)
		ePoint += eStep
		sumOdd += 4*SNS_DynesDOS(ePoint-w[4], w[1], w[2])*1e-3*SNS_ThermalKernelScalar((ePoint-x)*1e-3, w[0], 8.62e-5)
		ePoint += eStep
		sumEven += 2*SNS_DynesDOS(ePoint-w[4], w[1], w[2])*1e-3*SNS_ThermalKernelScalar((ePoint-x)*1e-3, w[0], 8.62e-5)
	endfor

	sumTotal += sumEven + sumOdd
	return w[3]*sumTotal*eStep/3
End


//==============================================================================
// SNS_FitFunc_2DEGStep / SNS_FitFunc_Bulk3D
//
// Large-bias conductance components for spectroscopy analysis. Energy is in
// meV. The surface-state component is a Gaussian-broadened step; the bulk
// component follows the free-electron square-root background for EF = 7 eV.
//==============================================================================
Function SNS_FitFunc_2DEGStep(cw, x) : FitFunc
    Wave cw
    Variable x

    Variable a0 = cw[0]
    Variable a1 = cw[1]
    Variable amplitude = cw[2]
    Variable onset = cw[3]
    Variable sigma = cw[4]
    Variable arg = (x - onset) / (sqrt(2) * sigma)

    return a0 + a1*x + amplitude*0.5*(1 + erf(arg))
End


Function SNS_FitFunc_Bulk3D(cw, x) : FitFunc
    Wave cw
    Variable x

    Variable arg = 1 + x/7000
    if (arg < 0)
        return NaN
    endif

    return cw[0] * sqrt(arg)
End


//==============================================================================
// SNS_Fit2DEG_OnsetAndBroadening_2D
//
// Fit a broadened step along dim 0 of every spectrum in a two-dimensional
// line-spectroscopy wave. Rejected fits are retained as NaN diagnostics.
//==============================================================================
Function SNS_Fit2DEG_OnsetAndBroadening_2D(G2D, onsetGuess, sigmaGuess)
    Wave G2D
    Variable onsetGuess, sigmaGuess

    Variable biasPointCount = DimSize(G2D, 0)
    Variable positionPointCount = DimSize(G2D, 1)
    if (biasPointCount <= 5 || positionPointCount <= 0)
        Abort "SNS_Fit2DEG_OnsetAndBroadening_2D: input wave size invalid."
    endif

    Variable biasStart = DimOffset(G2D, 0)
    Variable biasStep = DimDelta(G2D, 0)
    String biasUnit = WaveUnits(G2D, 0)
    Variable biasEnd = biasStart + biasStep*(biasPointCount - 1)
    Variable biasMin = min(biasStart, biasEnd)
    Variable biasMax = max(biasStart, biasEnd)
    Variable biasRange = biasMax - biasMin
    Variable biasStepAbs = abs(biasStep)
    if (biasRange <= 0 || biasStepAbs <= 0)
        Abort "SNS_Fit2DEG_OnsetAndBroadening_2D: invalid bias scaling."
    endif

    Make/O/D/N=(biasPointCount) BiasAxis_2DEG
    BiasAxis_2DEG = biasStart + biasStep*p
    SetScale/P x, biasStart, biasStep, biasUnit, BiasAxis_2DEG

    String baseName = NameOfWave(G2D)
    String nameOnset = "Onset_2DEG_" + baseName
    String nameBroad = "Broad_2DEG_" + baseName
    String nameFitOK = "FitOK_2DEG_" + baseName
    String nameR2 = "R2_2DEG_" + baseName
    String nameRMSE = "RMSE_2DEG_" + baseName
    Make/O/D/N=(positionPointCount) $nameOnset, $nameBroad, $nameFitOK, $nameR2, $nameRMSE
    Wave onsetMap = $nameOnset
    Wave broadMap = $nameBroad
    Wave fitOKMap = $nameFitOK
    Wave r2Map = $nameR2
    Wave rmseMap = $nameRMSE
    onsetMap = NaN
    broadMap = NaN
    fitOKMap = 0
    r2Map = NaN
    rmseMap = NaN

    Variable positionStart = DimOffset(G2D, 1)
    Variable positionStep = DimDelta(G2D, 1)
    String positionUnit = WaveUnits(G2D, 1)
    SetScale/P x, positionStart, positionStep, positionUnit, onsetMap, broadMap, fitOKMap, r2Map, rmseMap

    Make/O/D/N=(biasPointCount) Trace_2DEG, Fit_2DEG, ResidualSq_2DEG, TraceVar_2DEG
    Make/O/D/N=5 Coef_2DEGStep

    Variable positionIndex, biasIndex, fitErrorCode
    Variable traceMean, residualMeanSq, traceVarMean, fitRMSE, fitR2
    Variable fitOnset, fitSigma, fitAccept, fitSigmaAbs
    Variable minSigmaAllowed = biasStepAbs
    Variable maxSigmaAllowed = 0.35*biasRange
    Variable traceFinitePoints, traceMin, traceMax, traceRange, initialAmplitude

    if (numtype(sigmaGuess) != 0 || sigmaGuess <= 0)
        sigmaGuess = 0.05*biasRange
    endif
    sigmaGuess = max(abs(sigmaGuess), minSigmaAllowed)

    for (positionIndex = 0; positionIndex < positionPointCount; positionIndex += 1)
        for (biasIndex = 0; biasIndex < biasPointCount; biasIndex += 1)
            Trace_2DEG[biasIndex] = G2D[biasIndex][positionIndex]
        endfor

        WaveStats/Q Trace_2DEG
        traceFinitePoints = V_npnts
        traceMean = V_avg
        traceMin = V_min
        traceMax = V_max
        traceRange = traceMax - traceMin
        if (traceFinitePoints < 6 || numtype(traceRange) != 0 || traceRange <= 0)
            continue
        endif

        initialAmplitude = Trace_2DEG[biasPointCount - 1] - Trace_2DEG[0]
        Coef_2DEGStep[0] = Trace_2DEG[0]
        Coef_2DEGStep[1] = 0
        Coef_2DEGStep[2] = initialAmplitude
        Coef_2DEGStep[3] = onsetGuess
        Coef_2DEGStep[4] = sigmaGuess

        fitErrorCode = 0
        try
            FuncFit/Q SNS_FitFunc_2DEGStep, Coef_2DEGStep, Trace_2DEG /X=BiasAxis_2DEG /D=Fit_2DEG; AbortOnRTE
        catch
            fitErrorCode = GetRTError(1)
        endtry
        if (fitErrorCode != 0)
            continue
        endif

        fitOnset = Coef_2DEGStep[3]
        fitSigma = Coef_2DEGStep[4]
        fitSigmaAbs = abs(fitSigma)
        ResidualSq_2DEG = (Trace_2DEG - Fit_2DEG)^2
        TraceVar_2DEG = (Trace_2DEG - traceMean)^2
        WaveStats/Q ResidualSq_2DEG
        residualMeanSq = V_avg
        fitRMSE = sqrt(residualMeanSq)
        WaveStats/Q TraceVar_2DEG
        traceVarMean = V_avg
        fitR2 = traceVarMean > 0 ? 1 - residualMeanSq/traceVarMean : NaN
        r2Map[positionIndex] = fitR2
        rmseMap[positionIndex] = fitRMSE

        fitAccept = 1
        if (numtype(fitOnset) != 0 || numtype(fitSigma) != 0)
            fitAccept = 0
        endif
        if (fitOnset < biasMin || fitOnset > biasMax)
            fitAccept = 0
        endif
        if (fitSigmaAbs <= minSigmaAllowed || fitSigmaAbs > maxSigmaAllowed)
            fitAccept = 0
        endif
        if (numtype(fitR2) != 0 || fitR2 < 0.6)
            fitAccept = 0
        endif
        if (fitAccept)
            onsetMap[positionIndex] = fitOnset
            broadMap[positionIndex] = fitSigmaAbs
            fitOKMap[positionIndex] = 1
        endif
    endfor

    return 0
End


//==============================================================================
// SNS_AvgSpectraFromValidFit
//
// Average dim-1 spectra whose matching onset entry is finite.
//==============================================================================
Function SNS_AvgSpectraFromValidFit(dIdV2D, onsetW, outWavePath, outCountVarPath)
    Wave dIdV2D, onsetW
    String outWavePath, outCountVarPath

    Variable nE = DimSize(dIdV2D, 0)
    Variable nSpec = DimSize(dIdV2D, 1)
    if (numpnts(onsetW) != nSpec)
        Abort "SNS_AvgSpectraFromValidFit: onset wave length does not match number of spectra."
    endif

    Make/O/D/N=(nE) $outWavePath
    Wave outW = $outWavePath
    SetScale/P x, DimOffset(dIdV2D,0), DimDelta(dIdV2D,0), WaveUnits(dIdV2D,0), outW

    Variable iE, iSpec, sumVal, nVal, nValidSpec = 0
    for (iSpec = 0; iSpec < nSpec; iSpec += 1)
        if (numtype(onsetW[iSpec]) == 0)
            nValidSpec += 1
        endif
    endfor

    for (iE = 0; iE < nE; iE += 1)
        sumVal = 0
        nVal = 0
        for (iSpec = 0; iSpec < nSpec; iSpec += 1)
            if (numtype(onsetW[iSpec]) == 0 && numtype(dIdV2D[iE][iSpec]) == 0)
                sumVal += dIdV2D[iE][iSpec]
                nVal += 1
            endif
        endfor
        outW[iE] = nVal > 0 ? sumVal/nVal : NaN
    endfor

    Variable/G $outCountVarPath = nValidSpec
    return nValidSpec
End


//==============================================================================
// SNS__CopyRequiredRayTraceWave
//
// Purpose:
//   Copy a required ray-trace wave after validating existence and channel count.
//
// Inputs:
//   srcName   : full source wave path
//   dstName   : full destination wave path
//   nChannels : required number of channels
//
// Outputs:
//   Duplicates srcName to dstName.
//==============================================================================
Function SNS__CopyRequiredRayTraceWave(srcName, dstName, nChannels)
	String srcName, dstName
	Variable nChannels

	Wave/Z srcWave = $srcName

	if (!WaveExists(srcWave))
		Abort "Missing ray-trace wave: " + srcName
	endif
	if (numpnts(srcWave) != nChannels)
		Abort "Ray-trace wave length mismatch: " + srcName
	endif

	Duplicate/O srcWave, $dstName
End


//==============================================================================
// SNS__AssertUniqueOptionalNames
//
// Purpose:
//   Defensive helper for detecting duplicate optional-parameter names in long
//   Igor function signatures.
//
// Inputs:
//   optionalNameList : semicolon-separated optional-parameter names
//
// Outputs:
//   Aborts if a duplicate name is found.
//==============================================================================
Function SNS__AssertUniqueOptionalNames(optionalNameList)
	String optionalNameList

	Variable nItems = ItemsInList(optionalNameList, ";")
	Variable ind1, ind2
	String name1, name2

	for (ind1 = 0; ind1 < nItems; ind1 += 1)
		name1 = StringFromList(ind1, optionalNameList, ";")
		name1 = LowerStr(ReplaceString(" ", name1, ""))

		if (strlen(name1) == 0)
			continue
		endif

		for (ind2 = ind1 + 1; ind2 < nItems; ind2 += 1)
			name2 = StringFromList(ind2, optionalNameList, ";")
			name2 = LowerStr(ReplaceString(" ", name2, ""))

			if (CmpStr(name1, name2) == 0)
				Abort "Duplicate optional parameter name in function signature: " + name1
			endif
		endfor
	endfor
End

//==============================================================================
// ZBC-vs-B fitting helpers
//
// Purpose:
//   Stage-2 fitting routines for zero-bias conductance versus magnetic field.
//
// Shared model:
//   These routines reuse the same 2D/3D gap handling and raw-mix-then-
//   broaden construction as SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders.
//   The only difference is that the model is evaluated at E = 0 for each
//   experimental B value.
//
// Internal dependencies:
//   SNS__AssertUniqueOptionalNames(...)
//   SNS__CopyRayTraceDOSFitWavesFromFolder(...)
//   SNS__RayTrace2D3D_BuildComponents(...)
//
// Recommended workflow:
//   fitMode = 0 : fit h_eff_nm only; gammaB2 fixed.
//   fitMode = 1 : fit gammaB2 only; h_eff fixed from previous run.
//   fitMode = 2 : fit both; diagnostic only, can be ill-conditioned.
//
// Notes:
//   h_eff is constrained to hEff0_nm*(1 +/- hEffBoundRel), default +/-25%.
//==============================================================================

//==============================================================================
// SNS_ZBCFit_StoreStage1Params
//
// Purpose:
//   Store Stage-1 LDOS fit results in a named data folder for later ZBC(B)
//   fitting.
//
// Inputs:
//   DeltaSurface0_meV : fitted 2D/surface gap [meV]
//   DeltaBulk0_meV    : fitted 3D/bulk gap [meV]
//   frac2D0, frac3D0  : fitted 2D/3D fractions
//   Amplitude0        : fitted LDOS amplitude
//   Broadening0_meV   : fitted intrinsic Lorentzian broadening [meV]
//
// Optional inputs:
//   outFolder         : destination data folder, default root:SNS_ZBCFit
//
// Outputs:
//   Creates scalar globals in outFolder.
//==============================================================================
Function SNS_ZBCFit_StoreStage1Params(DeltaSurface0_meV, DeltaBulk0_meV, frac2D0, frac3D0, Amplitude0, Broadening0_meV, [outFolder])
	Variable DeltaSurface0_meV, DeltaBulk0_meV
	Variable frac2D0, frac3D0
	Variable Amplitude0, Broadening0_meV
	String outFolder

	String oldDF = GetDataFolder(1)
	SNS__AssertUniqueOptionalNames("outFolder")

	if (ParamIsDefault(outFolder))
		outFolder = "root:SNS_ZBCFit"
	endif

	DeltaSurface0_meV = abs(DeltaSurface0_meV)
	DeltaBulk0_meV = abs(DeltaBulk0_meV)

	if (numtype(DeltaSurface0_meV) != 0 || DeltaSurface0_meV <= 0)
		Abort "SNS_ZBCFit_StoreStage1Params: invalid DeltaSurface0_meV."
	endif
	if (numtype(DeltaBulk0_meV) != 0 || DeltaBulk0_meV <= 0)
		Abort "SNS_ZBCFit_StoreStage1Params: invalid DeltaBulk0_meV."
	endif
	if (DeltaBulk0_meV <= DeltaSurface0_meV)
		Abort "SNS_ZBCFit_StoreStage1Params: require DeltaBulk0_meV > DeltaSurface0_meV."
	endif
	if (numtype(frac2D0) != 0 || frac2D0 < 0 || numtype(frac3D0) != 0 || frac3D0 < 0 || frac2D0 + frac3D0 <= 0)
		Abort "SNS_ZBCFit_StoreStage1Params: invalid fractions."
	endif
	if (numtype(Amplitude0) != 0 || Amplitude0 <= 0)
		Abort "SNS_ZBCFit_StoreStage1Params: invalid Amplitude0."
	endif
	if (numtype(Broadening0_meV) != 0 || Broadening0_meV <= 0)
		Abort "SNS_ZBCFit_StoreStage1Params: invalid Broadening0_meV."
	endif

	NewDataFolder/O/S $outFolder
	Variable fsum = frac2D0 + frac3D0

	Variable/G Delta_surface_0_meV = DeltaSurface0_meV
	Variable/G Delta_bulk_0_meV = DeltaBulk0_meV
	Variable/G frac2D_0 = frac2D0/fsum
	Variable/G frac3D_0 = frac3D0/fsum
	Variable/G Amplitude_0 = Amplitude0
	Variable/G Broadening_0_meV = Broadening0_meV

	SetDataFolder $oldDF
	return 0
End


//==============================================================================
// SNS_FitZBCvsB_2D3D_RayTrace
//
// Purpose:
//   Fit zero-bias conductance versus magnetic field using the same mixed
//   2D surface + 3D bulk ray-trace DOS model as the LDOS(E) fit.
//
// Inputs:
//   LDOS_EB              : experimental LDOS wave, dimensions E x B
//   rayTraceFolder2D     : folder containing 2D ray-trace waves [nm schema]
//   rayTraceFolder3D     : folder containing 3D ray-trace waves [nm schema]
//   DeltaSurface0_meV    : fixed Stage-1 surface gap [meV]
//   DeltaBulk0_meV       : fixed Stage-1 bulk gap [meV]
//   frac2D0, frac3D0     : fixed Stage-1 2D/3D fractions
//   Amplitude0           : fixed Stage-1 amplitude
//   Broadening0_meV      : fixed Stage-1 intrinsic broadening [meV]
//
// Optional inputs:
//   BMin_mT, BMax_mT     : B-axis fit window [mT]
//   hEff0_nm             : starting effective magnetic length [nm]
//   hEffBoundRel         : relative bound around hEff0_nm
//   gammaB2_0_meV_per_T2 : starting quadratic depairing coefficient [meV/T^2]
//   fitMode              : choose which Stage-2 parameters are fit
//   hEff_fixed_nm        : fixed h_eff for gamma-only fits [nm]
//   gammaB2_fixed...     : fixed gammaB2 for h-only fits [meV/T^2]
//   xV_nm, yV_nm         : vortex position used by legacy endpoint phase [nm]
//   outPrefix            : output-wave prefix
//   weight*              : weighting controls for ZBC(B)
//   doDisplay            : display fit graph
//
// Outputs:
//   Creates fit waves and scalar results in root:SNS_RayTraceDOSFit.
//
// Fit modes:
//   fitMode = 0 : fit h_eff_nm only; gammaB2 fixed.
//   fitMode = 1 : fit gammaB2 only; h_eff fixed.
//   fitMode = 2 : fit both; diagnostic only.
//
// Notes:
//   Reuses SNS__RayTrace2D3D_BuildComponents(), so the ZBC model follows the
//   same 2D surface / 3D bulk gap logic, fractional mixing, and experimental
//   resolution convolution as SNS_FitLDOS_2D3D_RayTraceDOS_FromFolders().
//
//   The base weight can emphasize low LDOS, exactly like the LDOS(E) fit:
//      W_base = 1/(abs(ZBC)+floor)^2
//   Additional multiplicative importance weights can emphasize:
//      |B| <= BFirstMax_mT
//      ZBC near the low/gapped branch.
//==============================================================================
Function SNS_FitZBCvsB_2D3D_RayTrace(LDOS_EB, rayTraceFolder2D, rayTraceFolder3D, DeltaSurface0_meV, DeltaBulk0_meV, frac2D0, frac3D0, Amplitude0, Broadening0_meV, [BMin_mT, BMax_mT, hEff0_nm, hEffBoundRel, gammaB2_0_meV_per_T2, fitMode, hEff_fixed_nm, gammaB2_fixed_meV_per_T2, xV_nm, yV_nm, outPrefix, weightLowLDOS, weightFloorRel, BFirstMax_mT, weightPreMax, lowZBCThresholdRel, weightLowZBC, weightCap, doDisplay])
	Wave LDOS_EB
	String rayTraceFolder2D, rayTraceFolder3D
	Variable DeltaSurface0_meV, DeltaBulk0_meV
	Variable frac2D0, frac3D0
	Variable Amplitude0, Broadening0_meV
	Variable BMin_mT, BMax_mT
	Variable hEff0_nm, hEffBoundRel, gammaB2_0_meV_per_T2
	Variable fitMode, hEff_fixed_nm, gammaB2_fixed_meV_per_T2
	Variable xV_nm, yV_nm
	String outPrefix
	Variable weightLowLDOS, weightFloorRel
	Variable BFirstMax_mT, weightPreMax, lowZBCThresholdRel, weightLowZBC, weightCap
	Variable doDisplay

	String oldDF = GetDataFolder(1)

	SNS__AssertUniqueOptionalNames( \
		"BMin_mT;BMax_mT;hEff0_nm;hEffBoundRel;gammaB2_0_meV_per_T2;fitMode;" + \
		"hEff_fixed_nm;gammaB2_fixed_meV_per_T2;xV_nm;yV_nm;outPrefix;" + \
		"weightLowLDOS;weightFloorRel;BFirstMax_mT;weightPreMax;" + \
		"lowZBCThresholdRel;weightLowZBC;weightCap;doDisplay")

	if (!DataFolderExists("root:SNS_Settings:"))
		SNS_InitDefaultSettings()
	endif

	DeltaSurface0_meV = abs(DeltaSurface0_meV)
	DeltaBulk0_meV = abs(DeltaBulk0_meV)

	if (ParamIsDefault(outPrefix))
		outPrefix = "ZBCvsB_2D3D_Fit"
	endif
	if (ParamIsDefault(hEff0_nm))
		hEff0_nm = 24
	endif
	if (ParamIsDefault(hEffBoundRel))
		hEffBoundRel = 0.25
	endif
	if (ParamIsDefault(gammaB2_0_meV_per_T2))
		gammaB2_0_meV_per_T2 = 0
	endif
	if (ParamIsDefault(fitMode))
		fitMode = 0
	endif
	fitMode = round(fitMode)
	if (ParamIsDefault(hEff_fixed_nm))
		hEff_fixed_nm = hEff0_nm
	endif
	if (ParamIsDefault(gammaB2_fixed_meV_per_T2))
		gammaB2_fixed_meV_per_T2 = gammaB2_0_meV_per_T2
	endif
	if (ParamIsDefault(xV_nm))
		xV_nm = 0
	endif
	if (ParamIsDefault(yV_nm))
		yV_nm = 0
	endif
	if (ParamIsDefault(weightLowLDOS))
		weightLowLDOS = 1
	endif
	if (ParamIsDefault(weightFloorRel))
		weightFloorRel = 0.05
	endif
	if (ParamIsDefault(weightPreMax))
		weightPreMax = 3
	endif
	if (ParamIsDefault(lowZBCThresholdRel))
		lowZBCThresholdRel = 0.25
	endif
	if (ParamIsDefault(weightLowZBC))
		weightLowZBC = 5
	endif
	if (ParamIsDefault(weightCap))
		weightCap = 15
	endif
	if (ParamIsDefault(doDisplay))
		doDisplay = 1
	endif

	if (fitMode < 0 || fitMode > 2)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: fitMode must be 0, 1, or 2."
	endif
	if (numtype(DeltaSurface0_meV) != 0 || DeltaSurface0_meV <= 0)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: invalid DeltaSurface0_meV."
	endif
	if (numtype(DeltaBulk0_meV) != 0 || DeltaBulk0_meV <= 0)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: invalid DeltaBulk0_meV."
	endif
	if (DeltaBulk0_meV <= DeltaSurface0_meV)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: require DeltaBulk0_meV > DeltaSurface0_meV."
	endif
	if (numtype(frac2D0) != 0 || frac2D0 < 0 || numtype(frac3D0) != 0 || frac3D0 < 0 || frac2D0 + frac3D0 <= 0)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: invalid fractions."
	endif
	if (numtype(Amplitude0) != 0 || Amplitude0 <= 0 || numtype(Broadening0_meV) != 0 || Broadening0_meV <= 0)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: invalid amplitude or broadening."
	endif
	if (numtype(hEff0_nm) != 0 || hEff0_nm <= 0 || numtype(hEff_fixed_nm) != 0 || hEff_fixed_nm <= 0)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: invalid h_eff start/fixed value."
	endif
	if (numtype(hEffBoundRel) != 0 || hEffBoundRel <= 0 || hEffBoundRel >= 1)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: invalid hEffBoundRel."
	endif
	if (numtype(gammaB2_0_meV_per_T2) != 0 || gammaB2_0_meV_per_T2 < 0 || numtype(gammaB2_fixed_meV_per_T2) != 0 || gammaB2_fixed_meV_per_T2 < 0)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: invalid gammaB2 start/fixed value."
	endif
	if (numtype(weightFloorRel) != 0 || weightFloorRel <= 0)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: invalid weightFloorRel."
	endif
	if (numtype(weightPreMax) != 0 || weightPreMax <= 0 || numtype(weightLowZBC) != 0 || weightLowZBC <= 0)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: invalid weight multipliers."
	endif
	if (numtype(weightCap) != 0 || weightCap < 1)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: invalid weightCap."
	endif
	if (numtype(lowZBCThresholdRel) != 0 || lowZBCThresholdRel < 0 || lowZBCThresholdRel > 1)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: invalid lowZBCThresholdRel."
	endif

	Variable nEnergyRaw = DimSize(LDOS_EB, 0)
	Variable nFieldRaw = DimSize(LDOS_EB, 1)
	if (nEnergyRaw < 2 || nFieldRaw < 2)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: LDOS_EB must be E x B."
	endif

	// Reuse the same working folder/helper as the LDOS(E) fit.
	NewDataFolder/O root:SNS_RayTraceDOSFit
	SNS__CopyRayTraceDOSFitWavesFromFolder(rayTraceFolder2D, 0)
	SNS__CopyRayTraceDOSFitWavesFromFolder(rayTraceFolder3D, 1)

	SetDataFolder root:SNS_RayTraceDOSFit
	Variable fracSum = frac2D0 + frac3D0
	Variable/G g_ratio2D0 = frac2D0/fracSum
	Variable/G g_ratio3D0 = frac3D0/fracSum
	Variable/G g_xV_nm = xV_nm
	Variable/G g_yV_nm = yV_nm
	Variable/G g_zbc_DeltaSurface0_meV = DeltaSurface0_meV
	Variable/G g_zbc_DeltaBulk0_meV = DeltaBulk0_meV
	Variable/G g_zbc_Amplitude0 = Amplitude0
	Variable/G g_zbc_Broadening0_meV = Broadening0_meV
	Variable/G g_zbc_fitMode = fitMode
	Variable/G g_zbc_hEff_fixed_nm = hEff_fixed_nm
	Variable/G g_zbc_gammaB2_fixed_meV_T2 = gammaB2_fixed_meV_per_T2
	SetDataFolder $oldDF

	// Axis unit conversion.
	Variable eAxisOffset = DimOffset(LDOS_EB, 0)
	Variable eAxisDelta = DimDelta(LDOS_EB, 0)
	String eAxisUnit = WaveUnits(LDOS_EB, 0)
	Variable eScaleToMilliV = 1
	if (StringMatch(eAxisUnit, "*eV*") && !StringMatch(eAxisUnit, "*meV*") && !StringMatch(eAxisUnit, "*mV*"))
		eScaleToMilliV = 1000
	endif

	Variable bAxisOffset = DimOffset(LDOS_EB, 1)
	Variable bAxisDelta = DimDelta(LDOS_EB, 1)
	String bAxisUnit = WaveUnits(LDOS_EB, 1)
	Variable bScaleToMilliT = 1
	if (StringMatch(bAxisUnit, "*T*") && !StringMatch(bAxisUnit, "*mT*"))
		bScaleToMilliT = 1000
	endif

	Make/O/D/N=(nEnergyRaw) root:SNS_RayTraceDOSFit:tmpZBC_E_mV
	Make/O/D/N=(nEnergyRaw) root:SNS_RayTraceDOSFit:tmpZBC_col
	Wave eMilliVWave = root:SNS_RayTraceDOSFit:tmpZBC_E_mV
	Wave colWave = root:SNS_RayTraceDOSFit:tmpZBC_col

	Variable ii, jj, bHere_mT
	for (ii = 0; ii < nEnergyRaw; ii += 1)
		eMilliVWave[ii] = (eAxisOffset + ii*eAxisDelta) * eScaleToMilliV
	endfor

	Variable nFit = 0
	for (jj = 0; jj < nFieldRaw; jj += 1)
		bHere_mT = (bAxisOffset + jj*bAxisDelta) * bScaleToMilliT
		if ((!ParamIsDefault(BMin_mT) && bHere_mT < BMin_mT) || (!ParamIsDefault(BMax_mT) && bHere_mT > BMax_mT))
			continue
		endif
		nFit += 1
	endfor

	if (nFit < 5)
		SetDataFolder $oldDF
		Abort "SNS_FitZBCvsB_2D3D_RayTrace: too few B points in fit window."
	endif

	SetDataFolder root:SNS_RayTraceDOSFit
	Make/O/D/N=(nFit) $(outPrefix + "_Bfit_mT")
	Make/O/D/N=(nFit) $(outPrefix + "_ZBC_exp")
	Make/O/D/N=(nFit) $(outPrefix + "_fitWeights")
	Make/O/D/N=(nFit) $(outPrefix + "_importance")
	Wave bFitMilliTWave = $(outPrefix + "_Bfit_mT")
	Wave zbcExpWave = $(outPrefix + "_ZBC_exp")
	Wave fitWeightWave = $(outPrefix + "_fitWeights")
	Wave importanceWave = $(outPrefix + "_importance")
	SetDataFolder $oldDF

	Variable kk = 0
	for (jj = 0; jj < nFieldRaw; jj += 1)
		bHere_mT = (bAxisOffset + jj*bAxisDelta) * bScaleToMilliT
		if ((!ParamIsDefault(BMin_mT) && bHere_mT < BMin_mT) || (!ParamIsDefault(BMax_mT) && bHere_mT > BMax_mT))
			continue
		endif

		for (ii = 0; ii < nEnergyRaw; ii += 1)
			colWave[ii] = LDOS_EB[ii][jj]
		endfor

		bFitMilliTWave[kk] = bHere_mT
		zbcExpWave[kk] = interp(0, eMilliVWave, colWave)
		kk += 1
	endfor

	// Build bounded weights. Unlike the LDOS(E) fit, do not let the
	// inverse-low-LDOS factor become arbitrarily large; otherwise h_eff can
	// collapse toward zero because the fit only tries to keep the central gap dark.
	WaveStats/Q zbcExpWave
	Variable zbcMin = V_min
	Variable zbcMax = V_max
	Variable zbcRange = max(zbcMax - zbcMin, 1e-12)
	Variable zbcScale = max(abs(zbcMax), abs(zbcMin))
	Variable weightFloor = max(weightFloorRel*zbcScale, 1e-12)
	Variable lowZBCCut = zbcMin + lowZBCThresholdRel*zbcRange

	Variable firstMaxUse_mT
	if (ParamIsDefault(BFirstMax_mT))
		firstMaxUse_mT = Inf
	else
		firstMaxUse_mT = abs(BFirstMax_mT)
	endif

	importanceWave = 1
	if (firstMaxUse_mT > 0)
		importanceWave = importanceWave[p] * (abs(bFitMilliTWave[p]) <= firstMaxUse_mT ? weightPreMax : 1)
	endif
	importanceWave = importanceWave[p] * (zbcExpWave[p] <= lowZBCCut ? weightLowZBC : 1)

	if (weightLowLDOS)
		// Soft low-ZBC emphasis: ranges from ~1 at high ZBC to ~2 at the minimum.
		// The explicit lowZBC multiplier above carries the main low-gap emphasis.
		fitWeightWave = importanceWave[p] * (1 + ((zbcMax - zbcExpWave[p])/(zbcRange + weightFloor))^2)
	else
		fitWeightWave = importanceWave[p]
	endif

	// Normalize and cap the relative weights so no small set of central points
	// can dominate the field-scale determination.
	WaveStats/Q fitWeightWave
	Variable minWeight = max(V_min, 1e-12)
	fitWeightWave = fitWeightWave[p] / minWeight
	fitWeightWave = min(fitWeightWave[p], weightCap)

	// h_eff bounds: default +/-25% around hEff0_nm.
	Variable hEffMin_nm = hEff0_nm * (1 - hEffBoundRel)
	Variable hEffMax_nm = hEff0_nm * (1 + hEffBoundRel)

	// Fit coefficients and constraints.
	SetDataFolder root:SNS_RayTraceDOSFit
	Make/O/D/N=2 $(outPrefix + "_coef")
	Wave coefWave = $(outPrefix + "_coef")
	SetDataFolder $oldDF

	coefWave[0] = hEff0_nm
	coefWave[1] = gammaB2_0_meV_per_T2

	if (fitMode == 0)
		SetDataFolder root:SNS_RayTraceDOSFit
		Make/O/D/N=1 $(outPrefix + "_coef_fit_hOnly")
		Wave coefHWave = $(outPrefix + "_coef_fit_hOnly")
		Make/O/T/N=2 $(outPrefix + "_constraints_hOnly")
		Wave/T constraintHWave = $(outPrefix + "_constraints_hOnly")
		SetDataFolder $oldDF

		coefHWave[0] = hEff0_nm
		constraintHWave[0] = "K0 > " + num2str(hEffMin_nm)
		constraintHWave[1] = "K0 < " + num2str(hEffMax_nm)
		FuncFit/Q SNS__ZBC2D3D_FitFunc_hOnly coefHWave zbcExpWave /X=bFitMilliTWave /W=fitWeightWave /C=constraintHWave
		coefWave[0] = coefHWave[0]
		coefWave[1] = gammaB2_fixed_meV_per_T2
	elseif (fitMode == 1)
		SetDataFolder root:SNS_RayTraceDOSFit
		Make/O/D/N=1 $(outPrefix + "_coef_fit_gammaOnly")
		Wave coefGWave = $(outPrefix + "_coef_fit_gammaOnly")
		Make/O/T/N=1 $(outPrefix + "_constraints_gammaOnly")
		Wave/T constraintGWave = $(outPrefix + "_constraints_gammaOnly")
		SetDataFolder $oldDF

		coefGWave[0] = gammaB2_0_meV_per_T2
		constraintGWave[0] = "K0 > -1e-12"
		FuncFit/Q SNS__ZBC2D3D_FitFunc_gammaOnly coefGWave zbcExpWave /X=bFitMilliTWave /W=fitWeightWave /C=constraintGWave
		coefWave[0] = hEff_fixed_nm
		coefWave[1] = coefGWave[0]
	else
		SetDataFolder root:SNS_RayTraceDOSFit
		Make/O/T/N=3 $(outPrefix + "_constraints")
		Wave/T constraintTextWave = $(outPrefix + "_constraints")
		SetDataFolder $oldDF
		constraintTextWave[0] = "K0 > " + num2str(hEffMin_nm)
		constraintTextWave[1] = "K0 < " + num2str(hEffMax_nm)
		constraintTextWave[2] = "K1 > -1e-12"
		FuncFit/Q SNS__ZBC2D3D_FitFunc coefWave zbcExpWave /X=bFitMilliTWave /W=fitWeightWave /C=constraintTextWave
	endif

	SetDataFolder root:SNS_RayTraceDOSFit
	Make/O/D/N=(nFit) $(outPrefix + "_ZBC_fit")
	Make/O/D/N=(nFit) $(outPrefix + "_Gamma_pair_meV")
	Make/O/D/N=(nFit) $(outPrefix + "_resid")
	Wave zbcFitWave = $(outPrefix + "_ZBC_fit")
	Wave gammaPairWave = $(outPrefix + "_Gamma_pair_meV")
	Wave residWave = $(outPrefix + "_resid")
	SetDataFolder $oldDF

	SNS__ZBC2D3D_BuildTrace(coefWave, bFitMilliTWave, zbcFitWave, gammaPairWave)
	residWave = zbcExpWave[p] - zbcFitWave[p]

	SetDataFolder root:SNS_RayTraceDOSFit
	Variable/G $(outPrefix + "_fitMode") = fitMode
	Variable/G $(outPrefix + "_h_eff_nm_fit") = coefWave[0]
	Variable/G $(outPrefix + "_gammaB2_meV_per_T2_fit") = coefWave[1]
	Variable/G $(outPrefix + "_h_eff_nm_fixed") = hEff_fixed_nm
	Variable/G $(outPrefix + "_gammaB2_meV_per_T2_fixed") = gammaB2_fixed_meV_per_T2
	Variable/G $(outPrefix + "_DeltaSurface0_meV") = DeltaSurface0_meV
	Variable/G $(outPrefix + "_DeltaBulk0_meV") = DeltaBulk0_meV
	Variable/G $(outPrefix + "_frac2D0") = frac2D0/fracSum
	Variable/G $(outPrefix + "_frac3D0") = frac3D0/fracSum
	Variable/G $(outPrefix + "_Amplitude0") = Amplitude0
	Variable/G $(outPrefix + "_Broadening0_meV") = Broadening0_meV
	Variable/G $(outPrefix + "_hEffBoundRel") = hEffBoundRel
	Variable/G $(outPrefix + "_hEffMin_nm") = hEffMin_nm
	Variable/G $(outPrefix + "_hEffMax_nm") = hEffMax_nm
	Variable/G $(outPrefix + "_weightLowLDOS") = weightLowLDOS
	Variable/G $(outPrefix + "_weightFloorRel") = weightFloorRel
	Variable/G $(outPrefix + "_weightFloor") = weightFloor
	Variable/G $(outPrefix + "_BFirstMax_mT") = firstMaxUse_mT
	Variable/G $(outPrefix + "_weightPreMax") = weightPreMax
	Variable/G $(outPrefix + "_lowZBCThresholdRel") = lowZBCThresholdRel
	Variable/G $(outPrefix + "_lowZBCCut") = lowZBCCut
	Variable/G $(outPrefix + "_weightLowZBC") = weightLowZBC
	Variable/G $(outPrefix + "_weightCap") = weightCap
	SetDataFolder $oldDF

	if (doDisplay)
		Display/K=1 zbcExpWave vs bFitMilliTWave
		String graphName = WinName(0, 1, 1)
		AppendToGraph/W=$graphName zbcFitWave vs bFitMilliTWave

		String expTrace = NameOfWave(zbcExpWave)
		String fitTrace = NameOfWave(zbcFitWave)

		ModifyGraph/W=$graphName mode($expTrace)=3, marker($expTrace)=8, rgb($expTrace)=(0,0,0)
		ModifyGraph/W=$graphName lsize($fitTrace)=2, rgb($fitTrace)=(54741,24158,0)
		ModifyGraph/W=$graphName tick=2, mirror=2, standoff=0
		Label/W=$graphName bottom "B (mT)"
		Label/W=$graphName left "ZBC"

		String modeTxt
		if (fitMode == 0)
			modeTxt = "fit h_eff only; gammaB2 fixed"
		elseif (fitMode == 1)
			modeTxt = "fit gammaB2 only; h_eff fixed"
		else
			modeTxt = "fit h_eff and gammaB2"
		endif

		String txt
		txt = "ZBC-vs-B fit, same gap handling as LDOS(E) fit\r"
		txt += modeTxt + "\r"
		txt += "h_eff = " + num2str(round(1000*coefWave[0])/1000) + " nm\r"
		txt += "gamma_B2 = " + num2str(round(1000*coefWave[1])/1000) + " meV/T^2\r"
		txt += "Delta_surf = " + num2str(round(1000*DeltaSurface0_meV)/1000) + " meV, Delta_bulk = " + num2str(round(1000*DeltaBulk0_meV)/1000) + " meV\r"
		txt += "Broadening0 = " + num2str(round(1000*Broadening0_meV)/1000) + " meV\r"
		txt += "h_eff bounds = [" + num2str(round(1000*hEffMin_nm)/1000) + ", " + num2str(round(1000*hEffMax_nm)/1000) + "] nm"
		txt += "weightLowLDOS = " + num2str(weightLowLDOS) + ", preMaxW = " + num2str(weightPreMax) + ", lowZBCW = " + num2str(weightLowZBC) + ", cap = " + num2str(weightCap)
		TextBox/W=$graphName/C/N=ZBCFitInfo/F=0/A=RT/X=0/Y=0/B=(65535,65535,65535,39321) txt
	endif

	KillWaves/Z eMilliVWave, colWave
	SetDataFolder $oldDF

	return coefWave[0]
End


//==============================================================================
// SNS__ZBC2D3D_FitFunc
//
// Purpose:
//   Igor FuncFit callback for the two-parameter ZBC(B) diagnostic fit.
//
// Inputs:
//   coefWave[0]  : h_eff_nm
//   coefWave[1]  : gammaB2_meV_per_T2
//   xBMilliTWave : B-axis fit points [mT]
//
// Outputs:
//   yOutWave : modeled ZBC(B)
//==============================================================================
Function SNS__ZBC2D3D_FitFunc(coefWave, yOutWave, xBMilliTWave) : FitFunc
	Wave coefWave
	Wave yOutWave
	Wave xBMilliTWave

	Make/O/D/N=(numpnts(xBMilliTWave)) root:SNS_RayTraceDOSFit:tmp_ZBC_gamma_pair_meV
	Wave gammaPairWave = root:SNS_RayTraceDOSFit:tmp_ZBC_gamma_pair_meV

	SNS__ZBC2D3D_BuildTrace(coefWave, xBMilliTWave, yOutWave, gammaPairWave)
	return 0
End


//==============================================================================
// SNS__ZBC2D3D_FitFunc_hOnly
//
// Purpose:
//   Igor FuncFit callback for the h_eff-only ZBC(B) fit.
//
// Inputs:
//   coefWave[0]  : h_eff_nm
//   xBMilliTWave : B-axis fit points [mT]
//
// Reads from root:SNS_RayTraceDOSFit:
//   g_zbc_gammaB2_fixed_meV_T2 : fixed gammaB2 coefficient [meV/T^2]
//
// Outputs:
//   yOutWave : modeled ZBC(B)
//==============================================================================
Function SNS__ZBC2D3D_FitFunc_hOnly(coefWave, yOutWave, xBMilliTWave) : FitFunc
	Wave coefWave
	Wave yOutWave
	Wave xBMilliTWave

	NVAR gammaFixed = root:SNS_RayTraceDOSFit:g_zbc_gammaB2_fixed_meV_T2

	Make/O/D/N=2 root:SNS_RayTraceDOSFit:tmp_ZBC_coef2_hOnly
	Wave coef2Wave = root:SNS_RayTraceDOSFit:tmp_ZBC_coef2_hOnly
	coef2Wave[0] = coefWave[0]
	coef2Wave[1] = gammaFixed

	Make/O/D/N=(numpnts(xBMilliTWave)) root:SNS_RayTraceDOSFit:tmp_ZBC_gamma_pair_meV
	Wave gammaPairWave = root:SNS_RayTraceDOSFit:tmp_ZBC_gamma_pair_meV

	SNS__ZBC2D3D_BuildTrace(coef2Wave, xBMilliTWave, yOutWave, gammaPairWave)
	return 0
End


//==============================================================================
// SNS__ZBC2D3D_FitFunc_gammaOnly
//
// Purpose:
//   Igor FuncFit callback for the gammaB2-only ZBC(B) fit.
//
// Inputs:
//   coefWave[0]  : gammaB2_meV_per_T2
//   xBMilliTWave : B-axis fit points [mT]
//
// Reads from root:SNS_RayTraceDOSFit:
//   g_zbc_hEff_fixed_nm : fixed effective magnetic length [nm]
//
// Outputs:
//   yOutWave : modeled ZBC(B)
//==============================================================================
Function SNS__ZBC2D3D_FitFunc_gammaOnly(coefWave, yOutWave, xBMilliTWave) : FitFunc
	Wave coefWave
	Wave yOutWave
	Wave xBMilliTWave

	NVAR hFixed = root:SNS_RayTraceDOSFit:g_zbc_hEff_fixed_nm

	Make/O/D/N=2 root:SNS_RayTraceDOSFit:tmp_ZBC_coef2_gammaOnly
	Wave coef2Wave = root:SNS_RayTraceDOSFit:tmp_ZBC_coef2_gammaOnly
	coef2Wave[0] = hFixed
	coef2Wave[1] = coefWave[0]

	Make/O/D/N=(numpnts(xBMilliTWave)) root:SNS_RayTraceDOSFit:tmp_ZBC_gamma_pair_meV
	Wave gammaPairWave = root:SNS_RayTraceDOSFit:tmp_ZBC_gamma_pair_meV

	SNS__ZBC2D3D_BuildTrace(coef2Wave, xBMilliTWave, yOutWave, gammaPairWave)
	return 0
End


//==============================================================================
// SNS__ZBC2D3D_BuildTrace
//
// Purpose:
//   Build a modeled ZBC(B) trace for a supplied Stage-2 coefficient wave.
//
// Inputs:
//   coefWave[0]      : h_eff_nm
//   coefWave[1]      : gammaB2_meV_per_T2
//   bMilliTWave      : B-axis fit points [mT]
//
// Outputs:
//   zbcOutWave       : modeled ZBC(B)
//   gammaPairOutWave : gammaB2 * B^2 contribution [meV]
//
// Reads/writes:
//   Temporarily modifies root:SNS_Settings:lambdaL, Broadening,
//   SNS_GammaBase_eV, SNS_useGammaPair, and SNS_GammaPairScale, then restores
//   them before return.
//
// Notes:
//   Reuses SNS__RayTrace2D3D_BuildComponents(), so the ZBC calculation uses the
//   same gap assignment and raw-mix-then-resolution-convolution path as the
//   LDOS(E) fit.
//
//   The 6-component LDOS coefficient wave is:
//     [Amplitude0, DeltaBulk0, DeltaSurface0, 1, 1, Broadening0]
//   with g_ratio2D0/g_ratio3D0 already set to the fixed Stage-1 fractions.
//==============================================================================
Function SNS__ZBC2D3D_BuildTrace(coefWave, bMilliTWave, zbcOutWave, gammaPairOutWave)
	Wave coefWave
	Wave bMilliTWave
	Wave zbcOutWave, gammaPairOutWave

	String oldDF = GetDataFolder(1)
	SetDataFolder root:SNS_RayTraceDOSFit

	NVAR settingsLambdaL = root:SNS_Settings:lambdaL
	NVAR settingsBroadening = root:SNS_Settings:Broadening
	NVAR settingsGammaBase = root:SNS_Settings:SNS_GammaBase_eV
	NVAR settingsUseGammaPair = root:SNS_Settings:SNS_useGammaPair
	NVAR settingsGammaPairScale = root:SNS_Settings:SNS_GammaPairScale

	Variable lambdaLOriginal_m = settingsLambdaL
	Variable broadOriginal_eV = settingsBroadening
	Variable gammaBaseOriginal_eV = settingsGammaBase
	Variable useGammaPairOriginal = settingsUseGammaPair
	Variable gammaPairScaleOriginal_eV_T2 = settingsGammaPairScale

	NVAR DeltaSurface0_meV = root:SNS_RayTraceDOSFit:g_zbc_DeltaSurface0_meV
	NVAR DeltaBulk0_meV = root:SNS_RayTraceDOSFit:g_zbc_DeltaBulk0_meV
	NVAR Amplitude0 = root:SNS_RayTraceDOSFit:g_zbc_Amplitude0
	NVAR Broadening0_meV = root:SNS_RayTraceDOSFit:g_zbc_Broadening0_meV

	Variable hEff_nm = coefWave[0]
	Variable gammaB2_meV_per_T2 = coefWave[1]

	if (hEff_nm <= 0 || gammaB2_meV_per_T2 < 0)
		zbcOutWave = 0
		gammaPairOutWave = NaN
		SetDataFolder $oldDF
		return 0
	endif

	Make/O/D/N=1 root:SNS_RayTraceDOSFit:tmpZBC_E0_mV
	Wave e0Wave = root:SNS_RayTraceDOSFit:tmpZBC_E0_mV
	e0Wave[0] = 0

	Make/O/D/N=6 root:SNS_RayTraceDOSFit:tmpZBC_LDOS_coef
	Wave ldosCoefWave = root:SNS_RayTraceDOSFit:tmpZBC_LDOS_coef
	ldosCoefWave[0] = Amplitude0
	ldosCoefWave[1] = DeltaBulk0_meV
	ldosCoefWave[2] = DeltaSurface0_meV
	ldosCoefWave[3] = 1
	ldosCoefWave[4] = 1
	ldosCoefWave[5] = Broadening0_meV

	Make/O/D/N=1 root:SNS_RayTraceDOSFit:tmpZBC_surface_component
	Make/O/D/N=1 root:SNS_RayTraceDOSFit:tmpZBC_bulk_component
	Make/O/D/N=1 root:SNS_RayTraceDOSFit:tmpZBC_mixed_component
	Wave surfaceWave = root:SNS_RayTraceDOSFit:tmpZBC_surface_component
	Wave bulkWave = root:SNS_RayTraceDOSFit:tmpZBC_bulk_component
	Wave mixedWave = root:SNS_RayTraceDOSFit:tmpZBC_mixed_component

	Make/O/D/N=1 root:SNS_RayTraceDOSFit:bAxisFitT_fit
	Wave bModelTWave = root:SNS_RayTraceDOSFit:bAxisFitT_fit

	// Apply Stage-2 physical settings globally for the reused DOS builder.
	settingsLambdaL = hEff_nm * 1e-9
	settingsBroadening = Broadening0_meV * 1e-3
	settingsGammaBase = Broadening0_meV * 1e-3
	settingsUseGammaPair = 1
	settingsGammaPairScale = gammaB2_meV_per_T2 * 1e-3     // [eV/T^2]

	Variable jj
	Variable bExp_T
	for (jj = 0; jj < numpnts(bMilliTWave); jj += 1)
		bExp_T = bMilliTWave[jj] * 1e-3
		bModelTWave[0] = bExp_T
		gammaPairOutWave[jj] = gammaB2_meV_per_T2 * bExp_T * bExp_T

		SNS__RayTrace2D3D_BuildComponents(ldosCoefWave, e0Wave, surfaceWave, bulkWave, mixedWave)
		zbcOutWave[jj] = mixedWave[0]
	endfor

	// Restore global settings.
	settingsLambdaL = lambdaLOriginal_m
	settingsBroadening = broadOriginal_eV
	settingsGammaBase = gammaBaseOriginal_eV
	settingsUseGammaPair = useGammaPairOriginal
	settingsGammaPairScale = gammaPairScaleOriginal_eV_T2

	SetDataFolder $oldDF
	return 0
End
