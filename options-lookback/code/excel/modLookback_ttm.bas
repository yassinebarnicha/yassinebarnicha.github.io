Attribute VB_Name = "modLookback"
Option Explicit

Private Const SH_INPUTS As String = "Inputs"
Private Const SH_OUTPUTS As String = "Outputs"
Private Const SH_CURVE As String = "CurveData"
Private Const OUT_FIRST_ROW As Long = 4

'========================
' Public macros (buttons)
'========================

Public Sub Run_MC_And_Analytic()
    Dim t0 As Variant, t1 As Variant, ttm As Double, optType As String
    Dim spot As Double, rate As Double, vol As Double
    Dim paths As Long, seed As Long
    Dim exePath As String
    Dim cmd As String, outText As String

    t0 = GetInputValue("Valuation date T0 (YYYY-MM-DD)")
    t1 = GetInputValue("Maturity date T (YYYY-MM-DD)")

    'Excel stores dates as serial day counts -> difference is days
    ttm = (CDbl(t1) - CDbl(t0)) / 365#
    If ttm <= 0# Then
        MsgBox "T must be after T0 (ttm must be positive).", vbExclamation
        Exit Sub
    End If

    optType = CStr(GetInputValue("Option type"))
    spot = CDbl(GetInputValue("Spot S0"))
    rate = CDbl(GetInputValue("Rate r (cc)"))
    vol = CDbl(GetInputValue("Volatility σ"))
    paths = CLng(GetInputValue("Paths (Monte Carlo)"))
    seed = CLng(GetInputValue("Seed"))
    exePath = CStr(GetInputValue("C++ executable path (optional)"))

    If Len(exePath) = 0 Then
        MsgBox "Please set the C++ executable path in Inputs.", vbExclamation
        Exit Sub
    End If

    cmd = Quote(exePath) & _
          " --ttm " & FormatNum(ttm) & _
          " --type " & optType & _
          " --spot " & FormatNum(spot) & _
          " --rate " & FormatNum(rate) & _
          " --vol " & FormatNum(vol) & _
          " --paths " & CStr(paths) & _
          " --seed " & CStr(seed) & _
          " --csv"

    outText = RunCommandCapture(cmd)
    If Len(outText) = 0 Then
        MsgBox "No output captured from executable. Check path/permissions.", vbExclamation
        Exit Sub
    End If

    WriteResultsFromCSV outText
    MsgBox "Done: Outputs updated.", vbInformation
End Sub

Public Sub Generate_Curves()
    Dim t0 As Variant, t1 As Variant, ttm As Double, optType As String
    Dim rate As Double, vol As Double
    Dim paths As Long, seed As Long
    Dim exePath As String
    Dim wsC As Worksheet
    Dim lastRow As Long, i As Long
    Dim spot As Double
    Dim cmd As String, outText As String

    t0 = GetInputValue("Valuation date T0 (YYYY-MM-DD)")
    t1 = GetInputValue("Maturity date T (YYYY-MM-DD)")
    ttm = (CDbl(t1) - CDbl(t0)) / 365#
    If ttm <= 0# Then
        MsgBox "T must be after T0 (ttm must be positive).", vbExclamation
        Exit Sub
    End If

    optType = CStr(GetInputValue("Option type"))
    rate = CDbl(GetInputValue("Rate r (cc)"))
    vol = CDbl(GetInputValue("Volatility σ"))
    paths = CLng(GetInputValue("Paths (Monte Carlo)"))
    seed = CLng(GetInputValue("Seed"))
    exePath = CStr(GetInputValue("C++ executable path (optional)"))

    If Len(exePath) = 0 Then
        MsgBox "Please set the C++ executable path in Inputs.", vbExclamation
        Exit Sub
    End If

    Set wsC = ThisWorkbook.Worksheets(SH_CURVE)
    lastRow = wsC.Cells(wsC.Rows.Count, "A").End(xlUp).Row

    Application.ScreenUpdating = False
    Application.Calculation = xlCalculationManual
    On Error GoTo CleanFail

    For i = 2 To lastRow
        spot = CDbl(wsC.Cells(i, "A").Value)

        cmd = Quote(exePath) & _
              " --ttm " & FormatNum(ttm) & _
              " --type " & optType & _
              " --spot " & FormatNum(spot) & _
              " --rate " & FormatNum(rate) & _
              " --vol " & FormatNum(vol) & _
              " --paths " & CStr(paths) & _
              " --seed " & CStr(seed) & _
              " --csv"

        outText = RunCommandCapture(cmd)
        If Len(outText) > 0 Then
            WriteCurvePointFromCSV wsC, i, outText
        End If
    Next i

CleanOK:
    Application.Calculation = xlCalculationAutomatic
    Application.ScreenUpdating = True
    MsgBox "Curves generated.", vbInformation
    Exit Sub

CleanFail:
    Application.Calculation = xlCalculationAutomatic
    Application.ScreenUpdating = True
    MsgBox "Error while generating curves: " & Err.Description, vbExclamation
End Sub

Public Sub Clear_Outputs()
    Dim wsO As Worksheet, wsC As Worksheet
    Set wsO = ThisWorkbook.Worksheets(SH_OUTPUTS)
    Set wsC = ThisWorkbook.Worksheets(SH_CURVE)

    wsO.Range("C4:F" & wsO.Cells(wsO.Rows.Count, "A").End(xlUp).Row).ClearContents
    wsC.Range("B2:E" & wsC.Cells(wsC.Rows.Count, "A").End(xlUp).Row).ClearContents

    MsgBox "Outputs cleared.", vbInformation
End Sub

'========================
' Command execution
'========================

Private Function RunCommandCapture(ByVal cmd As String) As String
    Dim tmpFile As String
    Dim fullCmd As String

    tmpFile = TempTextFile()
    fullCmd = "/bin/zsh -lc " & Quote(cmd & " > " & Quote(tmpFile) & " 2>&1")
    Call Shell(fullCmd, vbHide)
    RunCommandCapture = WaitAndReadFile(tmpFile, 5#)
End Function

Private Function WaitAndReadFile(ByVal path As String, ByVal timeoutSec As Double) As String
    Dim t0 As Double
    t0 = Timer
    Do While Timer - t0 < timeoutSec
        If FileExists(path) Then
            If FileLen(path) > 0 Then Exit Do
        End If
        DoEvents
    Loop
    If Not FileExists(path) Then
        WaitAndReadFile = ""
        Exit Function
    End If
    WaitAndReadFile = ReadAllText(path)
End Function

Private Function TempTextFile() As String
    Dim p As String
    p = Environ$("TMPDIR")
    If Len(p) = 0 Then p = "/tmp/"
    If Right$(p, 1) <> "/" Then p = p & "/"
    TempTextFile = p & "lookback_out_" & Format(Now, "yyyymmdd_hhnnss") & "_" & CStr(Int(Rnd() * 100000)) & ".txt"
End Function

Private Function ReadAllText(ByVal path As String) As String
    Dim f As Integer
    f = FreeFile
    Open path For Input As #f
    ReadAllText = Input$(LOF(f), #f)
    Close #f
End Function

Private Function FileExists(ByVal path As String) As Boolean
    On Error GoTo Nope
    FileExists = (Len(Dir$(path)) > 0)
    Exit Function
Nope:
    FileExists = False
End Function

Private Function Quote(ByVal s As String) As String
    Quote = """" & Replace(s, """", "\""") & """"
End Function

Private Function FormatNum(ByVal x As Double) As String
    FormatNum = Replace$(CStr(x), ",", ".")
End Function

'========================
' Reading Inputs
'========================

Private Function GetInputValue(ByVal paramName As String) As Variant
    Dim ws As Worksheet
    Dim c As Range
    Set ws = ThisWorkbook.Worksheets(SH_INPUTS)
    Set c = ws.Range("A:A").Find(What:=paramName, LookIn:=xlValues, LookAt:=xlWhole)
    If c Is Nothing Then Err.Raise vbObjectError + 100, , "Input parameter not found: " & paramName
    GetInputValue = ws.Cells(c.Row, "B").Value
End Function

'========================
' Parsing CSV output
'========================

Private Sub WriteResultsFromCSV(ByVal txt As String)
    Dim wsO As Worksheet
    Dim lines() As String
    Dim i As Long
    Dim line As String
    Dim cols() As String

    Set wsO = ThisWorkbook.Worksheets(SH_OUTPUTS)
    lines = Split(Replace$(txt, vbCrLf, vbLf), vbLf)

    wsO.Range("C4:F" & wsO.Cells(wsO.Rows.Count, "A").End(xlUp).Row).ClearContents

    For i = LBound(lines) To UBound(lines)
        line = Trim$(lines(i))
        If Len(line) = 0 Then GoTo ContinueLoop
        If Left$(line, 6) = "engine" Then GoTo ContinueLoop

        cols = Split(line, ",")
        If UBound(cols) < 5 Then GoTo ContinueLoop
        WriteOneMetric wsO, cols(0), cols(1), cols(2), cols(3), cols(4), cols(5)
ContinueLoop:
    Next i
End Sub

Private Sub WriteOneMetric(ByVal wsO As Worksheet, ByVal engine As String, ByVal metric As String, _
                           ByVal est As String, ByVal se As String, ByVal lo As String, ByVal hi As String)
    Dim row As Long
    row = FindOutputRow(wsO, metric, engine)
    If row = 0 Then Exit Sub

    wsO.Cells(row, "C").Value = CDblSafe(est)
    wsO.Cells(row, "D").Value = CDblSafe(se)
    wsO.Cells(row, "E").Value = CDblSafe(lo)
    wsO.Cells(row, "F").Value = CDblSafe(hi)
End Sub

Private Function FindOutputRow(ByVal wsO As Worksheet, ByVal metric As String, ByVal engine As String) As Long
    Dim lastRow As Long, r As Long
    Dim m As String, e As String

    lastRow = wsO.Cells(wsO.Rows.Count, "A").End(xlUp).Row
    For r = OUT_FIRST_ROW To lastRow
        m = LCase$(CStr(wsO.Cells(r, "A").Value))
        e = LCase$(CStr(wsO.Cells(r, "B").Value))
        If m = LCase$(MetricDisplay(metric)) And e = LCase$(EngineDisplay(engine)) Then
            FindOutputRow = r
            Exit Function
        End If
    Next r
    FindOutputRow = 0
End Function

Private Function MetricDisplay(ByVal metric As String) As String
    Select Case LCase$(metric)
        Case "price": MetricDisplay = "Price"
        Case "delta": MetricDisplay = "Delta"
        Case "gamma": MetricDisplay = "Gamma"
        Case "theta": MetricDisplay = "Theta"
        Case "rho":   MetricDisplay = "Rho"
        Case "vega":  MetricDisplay = "Vega"
        Case Else:    MetricDisplay = metric
    End Select
End Function

Private Function EngineDisplay(ByVal engine As String) As String
    Select Case LCase$(engine)
        Case "montecarlo": EngineDisplay = "MonteCarlo"
        Case "analytic":   EngineDisplay = "Analytic"
        Case Else:         EngineDisplay = engine
    End Select
End Function

Private Function CDblSafe(ByVal s As String) As Double
    On Error GoTo Fail
    CDblSafe = CDbl(Replace$(s, ",", "."))
    Exit Function
Fail:
    CDblSafe = 0#
End Function

Private Sub WriteCurvePointFromCSV(ByVal wsC As Worksheet, ByVal row As Long, ByVal txt As String)
    Dim lines() As String, i As Long
    Dim line As String, cols() As String
    Dim engine As String, metric As String
    Dim est As Double

    lines = Split(Replace$(txt, vbCrLf, vbLf), vbLf)
    For i = LBound(lines) To UBound(lines)
        line = Trim$(lines(i))
        If Len(line) = 0 Then GoTo ContinueLoop
        If Left$(line, 6) = "engine" Then GoTo ContinueLoop

        cols = Split(line, ",")
        If UBound(cols) < 2 Then GoTo ContinueLoop

        engine = LCase$(cols(0))
        metric = LCase$(cols(1))
        est = CDblSafe(cols(2))

        If metric = "price" Then
            If engine = "montecarlo" Then wsC.Cells(row, "B").Value = est
            If engine = "analytic" Then wsC.Cells(row, "C").Value = est
        ElseIf metric = "delta" Then
            If engine = "montecarlo" Then wsC.Cells(row, "D").Value = est
            If engine = "analytic" Then wsC.Cells(row, "E").Value = est
        End If

ContinueLoop:
    Next i
End Sub
