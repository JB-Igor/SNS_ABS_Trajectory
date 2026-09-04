#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3              // Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}      // Set default tab width in Igor Pro 9 and later

//==============================================================================
// SNS_Interface.ipf
//
// Interface / boundary transparency models for SNS channel construction.
//
// Dependencies:
//   SNS_Core.ipf
//     - Structure SNS_Params
//     - SNS_LoadParams
//
// Notes:
//   This file should not build channels and should not solve ABS spectra.
//   It only converts incidence/interface parameters into transmission factors.
//==============================================================================


//==============================================================================
// SNS_ChannelTransmissionFromCos
//
// Interface transparency from incidence cosine.
//
// Reads all model/interface parameters from SNS_LoadParams(params).
//
// band model:
//   0 = normal 2DEG BTK:
//       T = cos^2(theta) / [cos^2(theta) + Z^2]
//
//   1 = TI Dirac / Klein:
//       T = cos^2(theta) / [cos^2(theta) + sin^2(ZD) sin^2(theta)]
//==============================================================================
Function SNS_ChannelTransmissionFromCos(cosInc)
    Variable cosInc

    STRUCT SNS_Params params
    SNS_LoadParams(params)

    Variable Zbarrier = params.BTK_barrier
    Variable bandModel = params.SNS_bandModel

    cosInc = min(1, max(0, abs(cosInc)))
    Variable cos2 = cosInc*cosInc

    if (Zbarrier <= 0)
        return 1
    endif

    // ------------------------------
    // TI Dirac / Klein transparency
    // ------------------------------
    if (bandModel == 1)

        Variable sinZ = sin(Zbarrier)
        Variable sin2theta = 1 - cos2

        if (abs(sinZ) < 1e-12)
            return 1
        endif
        if (cos2 <= 0)
            return 0
        endif

        return cos2 / (cos2 + sinZ*sinZ*sin2theta)
    endif

    // ------------------------------
    // Normal 2DEG BTK transparency
    // ------------------------------
    if (cos2 <= 0)
        return 0
    endif

    return cos2 / (cos2 + Zbarrier*Zbarrier)
End
