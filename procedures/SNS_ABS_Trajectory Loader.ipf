#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma DefaultTab={3,20,4}

StrConstant ksSNS_ABS_TrajectoryVersion = "0.9.0"

//==============================================================================
// SNS_ABS_Trajectory Loader
//
// Purpose:
//   Load/unload the SNS_ABS_Trajectory package by inserting/removing includes.
//
//==============================================================================


//------------------------------------------------------------------------------
// Loaded-state helper
//------------------------------------------------------------------------------
Function SNS_ABS_Trajectory_IsLoaded()
    NVAR/Z flag = root:Packages:SNS_ABS_Trajectory:gLoaded
    return (NVAR_Exists(flag) && flag != 0)
End


//------------------------------------------------------------------------------
// Dynamic menu item string generator
// itemNumber = 1 -> Load item (hidden when loaded)
// itemNumber = 2 -> Unload item (hidden when not loaded)
//------------------------------------------------------------------------------
Function/S SNS_ABS_Trajectory_MenuItem(itemNumber)
    Variable itemNumber

    Variable loaded = SNS_ABS_Trajectory_IsLoaded()

    if (itemNumber == 1)
        if (loaded)
            return ""                      // hide "Load" when loaded
        else
            return "Load SNS_ABS_Trajectory"
        endif
    endif

    if (itemNumber == 2)
        if (loaded)
            return "Unload SNS_ABS_Trajectory"
        else
            return ""                      // hide "Unload" when not loaded
        endif
    endif

    return ""
End


//------------------------------------------------------------------------------
// Macros menu (dynamic so the item strings are re-evaluated)
//------------------------------------------------------------------------------
Menu "Macros", dynamic
    SNS_ABS_Trajectory_MenuItem(1), /Q, Load_SNS_ABS_Trajectory()
    SNS_ABS_Trajectory_MenuItem(2), /Q, Unload_SNS_ABS_Trajectory()
End


//==============================================================================
// Load_SNS_ABS_Trajectory
//==============================================================================
Function Load_SNS_ABS_Trajectory()

    // Ensure clean state
    Unload_SNS_ABS_Trajectory()

    Execute /P /Q /Z "INSERTINCLUDE \"SNS_Core\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_Interface\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_Logging\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_BoundaryUtils\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_Helpers\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_GapExtraction\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_DisplayHelpers\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_GeometryFromMask\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_IslandSummary\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_Solver\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_Broadening\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_DOS\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_SpatialMaps\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_Maps\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_Fitting\""
    Execute /P /Q /Z "INSERTINCLUDE \"SNS_UI\""
    

    // --- Optional legacy (off by default) ---
    //Execute /P /Q /Z "INSERTINCLUDE \"SNS_Legacy_SSpS\""

    Execute /P /Q /Z "COMPILEPROCEDURES "

    // Mark as loaded
    NewDataFolder/O root:Packages
    NewDataFolder/O root:Packages:SNS_ABS_Trajectory
    Variable/G root:Packages:SNS_ABS_Trajectory:gLoaded = 1
    String/G root:Packages:SNS_ABS_Trajectory:sReleaseVersion = ksSNS_ABS_TrajectoryVersion

    Print "SNS_ABS_Trajectory v" + ksSNS_ABS_TrajectoryVersion
End


//==============================================================================
// Unload_SNS_ABS_Trajectory
//==============================================================================
Function Unload_SNS_ABS_Trajectory()

    Execute /P /Q /Z "DELETEINCLUDE \"SNS_UI\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_Fitting\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_Maps\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_SpatialMaps\"" 
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_DOS\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_Broadening\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_Solver\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_IslandSummary\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_GeometryFromMask\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_DisplayHelpers\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_GapExtraction\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_Helpers\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_BoundaryUtils\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_Logging\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_Interface\""
    Execute /P /Q /Z "DELETEINCLUDE \"SNS_Core\""

    Execute /P /Q /Z "COMPILEPROCEDURES "

    // Mark as unloaded
    NVAR/Z flag = root:Packages:SNS_ABS_Trajectory:gLoaded
    if (NVAR_Exists(flag))
        flag = 0
    endif
End
