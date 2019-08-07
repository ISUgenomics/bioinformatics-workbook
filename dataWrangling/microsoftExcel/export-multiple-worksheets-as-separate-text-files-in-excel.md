---
title: "Introduction to Data Wrangling"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Export multiple worksheets as separate text files in Excel

If there are large number of worksheets (tabs) in your excel file that you need to export as a separate text file, follow these guidelines. Note that the worksheet label will be used as file name for the text file with the `.txt` extension.

Before proceeding, check if you have enabled the macros, i.e., if you don't see `DEVELOPER` tab in you empty spreadsheet,  click on `FILE`, `OPTIONS` and `Customize Ribbon`. You should see a check box on the right hand side, for the `DEVELOPER` tab, check it and click `OK`.

Click on `DEVELOPER` and then `Macros`, type in some name (eg. `import_text`), click `create`.

Paste the below code on the popped window:
```
' ---------------------- Directory Choosing Helper Functions -----------------------
' Excel and VBA do not provide any convenient directory chooser or file chooser
' dialogs, but these functions will provide a reference to a system DLL
' with the necessary capabilities
Private Type BROWSEINFO ' used by the function GetFolderName
    hOwner As Long
    pidlRoot As Long
    pszDisplayName As String
    lpszTitle As String
    ulFlags As Long
    lpfn As Long
    lParam As Long
    iImage As Long
End Type

Private Declare Function SHGetPathFromIDList Lib "shell32.dll" _
    Alias "SHGetPathFromIDListA" (ByVal pidl As Long, ByVal pszPath As String) As Long
Private Declare Function SHBrowseForFolder Lib "shell32.dll" _
    Alias "SHBrowseForFolderA" (lpBrowseInfo As BROWSEINFO) As Long

Function GetFolderName(Msg As String) As String
' returns the name of the folder selected by the user
Dim bInfo As BROWSEINFO, path As String, r As Long
Dim X As Long, pos As Integer
    bInfo.pidlRoot = 0& ' Root folder = Desktop
    If IsMissing(Msg) Then
        bInfo.lpszTitle = "Select a folder."
        ' the dialog title
    Else
        bInfo.lpszTitle = Msg ' the dialog title
    End If
    bInfo.ulFlags = &H1 ' Type of directory to return
    X = SHBrowseForFolder(bInfo) ' display the dialog
    ' Parse the result
    path = Space$(512)
    r = SHGetPathFromIDList(ByVal X, ByVal path)
    If r Then
        pos = InStr(path, Chr$(0))
        GetFolderName = Left(path, pos - 1)
    Else
        GetFolderName = ""
    End If
End Function
'---------------------- END Directory Chooser Helper Functions ----------------------

Public Sub DoTheExport()
Dim FName As Variant
Dim Sep As String
Dim wsSheet As Worksheet
Dim nFileNum As Integer
Dim csvPath As String


Sep = InputBox("Enter a single delimiter character (e.g., comma or semi-colon)", _
"Export To Text File")
'csvPath = InputBox("Enter the full path to export CSV files to: ")

csvPath = GetFolderName("Choose the folder to export CSV files to:")
If csvPath = "" Then
    MsgBox ("You didn't choose an export directory. Nothing will be exported.")
    Exit Sub
End If

For Each wsSheet In Worksheets
wsSheet.Activate
nFileNum = FreeFile
Open csvPath & "\" & _
  wsSheet.Name & ".csv" For Output As #nFileNum
ExportToTextFile CStr(nFileNum), Sep, False
Close nFileNum
Next wsSheet

End Sub



Public Sub ExportToTextFile(nFileNum As Integer, _
Sep As String, SelectionOnly As Boolean)

Dim WholeLine As String
Dim RowNdx As Long
Dim ColNdx As Integer
Dim StartRow As Long
Dim EndRow As Long
Dim StartCol As Integer
Dim EndCol As Integer
Dim CellValue As String

Application.ScreenUpdating = False
On Error GoTo EndMacro:

If SelectionOnly = True Then
With Selection
StartRow = .Cells(1).Row
StartCol = .Cells(1).Column
EndRow = .Cells(.Cells.Count).Row
EndCol = .Cells(.Cells.Count).Column
End With
Else
With ActiveSheet.UsedRange
StartRow = .Cells(1).Row
StartCol = .Cells(1).Column
EndRow = .Cells(.Cells.Count).Row
EndCol = .Cells(.Cells.Count).Column
End With
End If

For RowNdx = StartRow To EndRow
WholeLine = ""
For ColNdx = StartCol To EndCol
If Cells(RowNdx, ColNdx).Value = "" Then
CellValue = ""
Else
CellValue = Cells(RowNdx, ColNdx).Value
End If
WholeLine = WholeLine & CellValue & Sep
Next ColNdx
WholeLine = Left(WholeLine, Len(WholeLine) - Len(Sep))
Print #nFileNum, WholeLine
Next RowNdx

EndMacro:
On Error GoTo 0
Application.ScreenUpdating = True

End Sub
```
Run the macro, but clicking on `Developer` tab, and `Run Macro`

[Table of contents](../dataAcquisition/dAc_introduction.md)
