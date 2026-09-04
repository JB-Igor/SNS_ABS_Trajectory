#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3
#pragma ModuleName=SNS_Logging

//==============================================================================
// SNS_Logging.ipf
//
// Run-scoped logging for the SNS_ABS_Trajectory package.
//
// Concept:
//   - A "run" log is created once per top-level operation.
//   - All SNS_Log(...) calls append to the active run until a new run is started.
//
// Storage:
//   root:SNS_Log:run_YYYYMMDD_HHMMSS:SNS_LogWave   (text wave)
//
// State:
//   root:SNS_Log:activeDF     (string)  -> folder path of active run
//   root:SNS_Log:activeTag    (string)  -> human-readable tag
//   root:SNS_Log:activeNUsed  (var)     -> number of used lines in active log
//   root:SNS_Log:activeNAlloc (var)     -> allocated size of log wave
//
// Debug gating:
//   SNS_Log(..., level="DBG") only records when root:SNS_Settings:SNS_Debug != 0.
//==============================================================================


//============================================================
// 0. INTERNAL HELPERS
//============================================================

Static Function SNS_LogEnsureRoot()
    NewDataFolder/O root:SNS_Log
    return 0
End

Static Function/S SNS_LogSanitizeName(s)
    String s
    s = ReplaceString(":", s, "")
    s = ReplaceString("-", s, "")
    s = ReplaceString(" ", s, "_")
    s = ReplaceString("/", s, "_")
    s = ReplaceString("\\", s, "_")
    s = ReplaceString(".", s, "_")
    s = ReplaceString(",", s, "_")
    s = ReplaceString(";", s, "_")
    s = ReplaceString("[", s, "")
    s = ReplaceString("]", s, "")
    s = ReplaceString("(", s, "")
    s = ReplaceString(")", s, "")
    return s
End

// Ensure the active log wave exists; if not, auto-start.
Static Function SNS_LogEnsureActive()
    SNS_LogEnsureRoot()

    SVAR/Z activeDF = root:SNS_Log:activeDF
    if (!SVAR_Exists(activeDF))
        SNS_LogStart("auto")
        return 0
    endif

    // Ensure index vars exist (in case user deleted them)
    NVAR/Z nUsed  = root:SNS_Log:activeNUsed
    NVAR/Z nAlloc = root:SNS_Log:activeNAlloc
    if (!NVAR_Exists(nUsed))
        Variable/G root:SNS_Log:activeNUsed = 0
    endif
    if (!NVAR_Exists(nAlloc))
        Variable/G root:SNS_Log:activeNAlloc = 0
    endif

    return 0
End


//============================================================
// 1. PUBLIC API
//============================================================

//==============================================================================
// SNS_LogStart
//
// Create a new run log and set it active.
// Returns: data folder path string of the created run.
//==============================================================================
Function/S SNS_LogStart(runTag)
    String runTag

    SNS_LogEnsureRoot()

    // Timestamp safe for names
    String ts = Secs2Date(DateTime, -2) + "_" + Secs2Time(DateTime, 3)
    ts = SNS_LogSanitizeName(ts)

    String df = "root:SNS_Log:run_" + ts
    NewDataFolder/O $df

    // Preallocate log wave (avoid name collision with SNS_Log function!)
    Variable initAlloc = 200
    Make/O/T/N=(initAlloc) $(df + ":SNS_LogWave")

    // Track used/alloc for this run
    Variable/G root:SNS_Log:activeNUsed  = 0
    Variable/G root:SNS_Log:activeNAlloc = initAlloc

    // metadata
    Make/O/T/N=1 $(df + ":SNS_LogTagW")
    Wave/T tagW = $(df + ":SNS_LogTagW")
    tagW[0] = runTag
    Variable/G $(df + ":SNS_LogT0") = DateTime

    // set active
    String/G root:SNS_Log:activeDF  = df
    String/G root:SNS_Log:activeTag = runTag

    // Start line
    SNS_Log("Log started: " + runTag, level="INFO")
    return df
End


//==============================================================================
// SNS_Log
//
// Append a line to the active run log.
// Optional:
//   level : "DBG", "INFO", "WARN", "ERR" (default "INFO")
//   cap   : keep last cap lines (default 2000)
// Notes:
//   - "DBG" messages are only recorded when root:SNS_Settings:SNS_Debug != 0.
//==============================================================================
Function SNS_Log(msg, [level, cap])
    String msg
    String level
    Variable cap

    if (ParamIsDefault(level))
        level = "INFO"
    endif
    if (ParamIsDefault(cap))
        cap = 2000
    endif

    SNS_LogEnsureActive()

    // Debug gating
    if (CmpStr(level, "DBG") == 0)
        NVAR/Z SNS_Debug = root:SNS_Settings:SNS_Debug
        if (!NVAR_Exists(SNS_Debug) || SNS_Debug == 0)
            return 0
        endif
    endif

    SVAR activeDF = root:SNS_Log:activeDF
    Wave/T logW = $(activeDF + ":SNS_LogWave")

    NVAR nUsed  = root:SNS_Log:activeNUsed
    NVAR nAlloc = root:SNS_Log:activeNAlloc

    // If the log is already at the cap, drop a chunk before appending.
    // This avoids a grow/delete cycle on every single log line in noisy loops.
    if (cap > 0 && nUsed >= cap)
        Variable dropAtCap = max(1, ceil(cap / 10))
        dropAtCap = min(dropAtCap, nUsed)
        DeletePoints 0, dropAtCap, logW
        nUsed -= dropAtCap
        nAlloc = DimSize(logW, 0)
    endif

    // Ensure capacity
    if (nUsed >= nAlloc)
        Variable growBy = 200
        Redimension/N=(nAlloc + growBy) logW
        nAlloc += growBy
    endif

    // Write line
    String t = Secs2Time(DateTime, 3)
    logW[nUsed] = t + " [" + level + "] " + msg
    nUsed += 1

    // Cap length (drop oldest) only when needed
    if (nUsed > cap)
        Variable drop = nUsed - cap
        DeletePoints 0, drop, logW
        nUsed -= drop
        nAlloc = DimSize(logW, 0)
    endif

    return 0
End


//==============================================================================
// SNS_LogEnd
//
// Mark end time of active run and write a closing line.
// Shrinks log wave to used length for cleanliness.
//==============================================================================
Function SNS_LogEnd()
    SNS_LogEnsureActive()

    SVAR activeDF = root:SNS_Log:activeDF
    NVAR nUsed  = root:SNS_Log:activeNUsed
    NVAR nAlloc = root:SNS_Log:activeNAlloc

    Variable/G $(activeDF + ":SNS_LogT1") = DateTime
    SNS_Log("Log ended.", level="INFO")

    // shrink to used size
    Wave/T logW = $(activeDF + ":SNS_LogWave")
    if (nUsed < DimSize(logW,0))
        Redimension/N=(nUsed) logW
        nAlloc = nUsed
    endif

    return 0
End


//==============================================================================
// SNS_LogShowActive
//
// Open the active log wave in a table. Shrinks to used length first.
//==============================================================================
Function SNS_LogShowActive()
    SNS_LogEnsureActive()

    SVAR activeDF = root:SNS_Log:activeDF
    NVAR nUsed  = root:SNS_Log:activeNUsed
    NVAR nAlloc = root:SNS_Log:activeNAlloc

    Wave/T logW = $(activeDF + ":SNS_LogWave")
    if (nUsed < DimSize(logW,0))
        Redimension/N=(nUsed) logW
        nAlloc = nUsed
    endif

    Edit logW
    return 0
End


//==============================================================================
// SNS_LogClearActive
//
// Clear active log content (keeps run folder). Keeps allocation.
//==============================================================================
Function SNS_LogClearActive()
    SNS_LogEnsureActive()

    SVAR activeDF = root:SNS_Log:activeDF
    NVAR nUsed  = root:SNS_Log:activeNUsed
    NVAR nAlloc = root:SNS_Log:activeNAlloc

    Wave/T logW = $(activeDF + ":SNS_LogWave")
    nUsed = 0
    SNS_Log("Log cleared.", level="WARN")
    return 0
End


//==============================================================================
// SNS_LogExportActive
//
// Export the active run log as a text wave through an existing Igor path.
//
// Inputs:
//   pathName : existing Igor export path name, e.g. SNSFittingDiagExportPath.
//   fileName : output file name, e.g. "SNS_active_log.txt".
//
// Notes:
//   - This keeps notebooks from reaching into root:SNS_Log internals.
//   - The active log wave is not modified; a temporary copy is exported.
//==============================================================================
Function SNS_LogExportActive(pathName, fileName)
    String pathName, fileName

    SNS_LogEnsureActive()

    SVAR activeDF = root:SNS_Log:activeDF
    NVAR nUsed  = root:SNS_Log:activeNUsed

    Wave/T logW = $(activeDF + ":SNS_LogWave")

    Variable nExport = nUsed
    if (numtype(nExport) != 0 || nExport < 0)
        nExport = 0
    endif
    nExport = min(nExport, DimSize(logW,0))

    Make/O/T/N=(nExport) $(activeDF + ":SNS_LogExportWave")
    Wave/T exportW = $(activeDF + ":SNS_LogExportWave")
    if (nExport > 0)
        exportW[] = logW[p]
    endif

    Save/O/J/P=$pathName exportW as fileName
    KillWaves/Z $(activeDF + ":SNS_LogExportWave")
    return 0
End
