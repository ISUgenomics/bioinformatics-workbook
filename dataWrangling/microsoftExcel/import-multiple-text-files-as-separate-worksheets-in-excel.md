---
title: "Introduction to Data Wrangling"
layout: single
header:
  overlay_color: "444444"
  overlay_image: /assets/images/dna.jpg
---

# Import multiple text files as separate worksheets in Excel #
If there are large number of text files that you need to import as a separate worksheet, follow these guidelines. Note that the file name of the text file will be used to label the worksheet (tab), without the `.txt` extension.
Before proceeding, check if you have enabled the macros, i.e., if you don't see `DEVELOPER` tab in you empty spreadsheet,  click on `FILE`, `OPTIONS` and `Customize Ribbon`. You should see a check box on the right hand side, for the `DEVELOPER` tab, check it and click `OK`.
Click on `DEVELOPER` and then `Macros`, type in some name (eg. `import_text`), click `create`.
Paste the below code on the popped window:

```
Sub CombineTextFiles()
    Dim FilesToOpen
    Dim x As Integer
    Dim wkbAll As Workbook
    Dim wkbTemp As Workbook
    Dim sDelimiter As String

    On Error GoTo ErrHandler
    Application.ScreenUpdating = False

    sDelimiter = "|"

    FilesToOpen = Application.GetOpenFilename _
      (FileFilter:="Text Files (*.txt), *.txt", _
      MultiSelect:=True, Title:="Text Files to Open")

    If TypeName(FilesToOpen) = "Boolean" Then
        MsgBox "No Files were selected"
        GoTo ExitHandler
    End If

    x = 1
    Set wkbTemp = Workbooks.Open(FileName:=FilesToOpen(x))
    wkbTemp.Sheets(1).Copy
    Set wkbAll = ActiveWorkbook
    wkbTemp.Close (False)
    wkbAll.Worksheets(x).Columns("A:A").TextToColumns _
      Destination:=Range("A1"), DataType:=xlDelimited, _
      TextQualifier:=xlDoubleQuote, _
      ConsecutiveDelimiter:=False, _
      Tab:=False, Semicolon:=False, _
      Comma:=False, Space:=False, _
      Other:=True, OtherChar:="|"
    x = x + 1

    While x <= UBound(FilesToOpen)
        Set wkbTemp = Workbooks.Open(FileName:=FilesToOpen(x))
        With wkbAll
            wkbTemp.Sheets(1).Move After:=.Sheets(.Sheets.Count)
            .Worksheets(x).Columns("A:A").TextToColumns _
              Destination:=Range("A1"), DataType:=xlDelimited, _
              TextQualifier:=xlDoubleQuote, _
              ConsecutiveDelimiter:=False, _
              Tab:=False, Semicolon:=False, _
              Comma:=False, Space:=False, _
              Other:=True, OtherChar:=sDelimiter
        End With
        x = x + 1
    Wend

ExitHandler:
    Application.ScreenUpdating = True
    Set wkbAll = Nothing
    Set wkbTemp = Nothing
    Exit Sub

ErrHandler:
    MsgBox Err.Description
    Resume ExitHandler
End Sub
```

Now you are all set to import files. Click on `Macros` again, and then `Run`. A window should pop up asking location of the text files. Navigate to the files location and select all text files. Click `OK`. All text files will be imported to a new spreadsheet. It might take a while to complete, if there are large number of files, so be patient.

[Table of contents](../../dataAcquisition/dAc_introduction.md)
