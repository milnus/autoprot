<#

    .SYNOPSIS
    The autoprot pipeline allows for absolute quantification of proteins from raw mass spectrometry 
    (MS) files in an automated manner.

    .DESCRIPTION
    The autoprot pipeline covers data analysis from both DIA and DDA methods, where a fully open-source 
    option is avalaible for DIA methods. Raw data from labelled, label-free and standard-free approaches 
    can be analysed with the pipeline. The normalisation of peptide intensities into protein intensities 
    is performed with seven different algorithms to identify the optimal algorithm for the current experiment. 
    The incorporated algorithms are Top3, Topall, iBAQ, APEX, NSAF, LFAQ, and xTop.

    .PARAMETER mode
    [string] mandatory specify the acquisition mode as "DDA", "DIA" or "directDIA".

    .PARAMETER approach
    [string] mandatory specify the quantification approach as "label", "unlabel" or "free".

    .PARAMETER InputDir
    [directory] mandatory specify the input directory containing all input files with raw MS spectra. 
    The output directory will be located in the input directory after the run.

    .PARAMETER ExpName
    [string] mandatory specify the name of the experiment.

    .PARAMETER fasta
    [file] mandatory specify the FASTA file with the proteome sequences.

    .PARAMETER totalProt
    [file] mandatory specify the file with total protein amount for each sample. 
    Optional to include the cell volume to be used for each sample

    .PARAMETER SpecLib
    [file] (mandatory for "DIA" mode) specify the file with the spectral library for the "DIA" mode.

    .PARAMETER BGSfasta
    [file] (mandatory for "directDIA" mode with Spectronaut) specify the FASTA file in .BGSfasta format, 
    which is required for the "directDIA" mode using Spectronaut (not open source).

    .PARAMETER ISconc
    [file] (mandatory for "label" and "unlabel" approaches) specify the file with the absolute 
    concentrations of each standard peptide ("label" approach) or protein ("unlabel" approach).

    .INPUTS
    Input can only be provided through the corresponding arguments.

    .OUTPUTS
    The output directory will be located in the input directory after the run and 
    will contain seven files with a protein concentration table, one for each algorithm.

    .EXAMPLE
    autoprot.ps1 [args]

    .EXAMPLE
    autoprot.ps1 -mode "directDIA" -approach "free" -InputDir "$PSScriptRoot\..\Examples\Input" 
    -ExpName "test_run" -fasta "$PSScriptRoot\..\Examples\Input\" -totalProt "$PSScriptRoot\..\Examples\Input\" 
    -BGSfasta "$PSScriptRoot\..\Examples\Input\"

    .LINK
    https://github.com/biosustain/autoprot

#>

param(
    [Parameter(Mandatory=$true)][string] $mode,
    [Parameter(Mandatory=$true)][string] $approach,
    [Parameter(Mandatory=$true)][string] $InputDir,
    [Parameter(Mandatory=$true)][string] $ExpName,
    [Parameter(Mandatory=$true)][string] $fasta,
    [Parameter(Mandatory=$true)][string] $totalProt,
    [string] $SpecLib,
    [string] $BGSfasta,
    [string] $ISconc
)

$ErrorActionPreference = 'Stop'
if ($mode -ne "DDA" -and $mode -ne "DIA" -and $mode -ne "directDIA") {
    Write-Error -Message "-mode must be either `"DDA`", `"DIA`" or `"directDIA`""
}
if ($approach -ne "label" -and $approach -ne "unlabel" -and $approach -ne "free") {
    Write-Error -Message "-approach must be either `"label`", `"unlabel`" or `"free`""
}
if ($mode -eq "DIA") {
    if (!$SpecLib) {Write-Error -Message "File with spectral library is required"}
}
if ($mode -eq "directDIA") {
    if (!$BGSfasta -or $BGSfasta -notmatch ".bgsfasta") {Write-Error -Message "File with fasta in BGSfasta format required"}
}
if ($approach -eq "label" -or $approach -eq "unlabel") {
    if (!$ISconc) {Write-Error -Message "File with IS concentrations is required"}
}

$OutputDir = Join-Path $InputDir ((Get-Date -format 'yyyyMMdd_HHmmss') + "_" + $ExpName + "_" + $mode + "_" + $approach)
$intermediate = "$OutputDir\intermediate_results"
$DIAanalysis = "$PSScriptRoot\DIA_analysis"
$conversions = "$PSScriptRoot\Conversions"
$normalisation = "$PSScriptRoot\Normalisation"
$quantification = "$PSScriptRoot\Quantification"
New-Item -ItemType Directory -Path $OutputDir | Out-Null
New-Item -ItemType Directory -Path $intermediate | Out-Null

## DIA analysis using Spectronaut
if ($mode -eq "DIA" -or $mode -eq "directDIA") {
    $fileType = ".*\.raw"
    if ($mode -eq "DIA") {
        $settings = "$DIAanalysis\DIA_settings.prop"
        $SNargsList = "-d $InputDir -a $SpecLib -s $settings -o $intermediate -n $ExpName -f $fileType"
    }
    elseif ($mode -eq "directDIA") {
        $settings = "$DIAanalysis\directDIA_settings.prop"
        $SNargsList = "-direct -d $InputDir -fasta $BGSfasta -s $settings -o $intermediate -n $ExpName -f $fileType"
    }

    Start-Process -FilePath spectronaut -ArgumentList $SNargsList -Wait
    $SNoutputDir = (Get-ChildItem -Path $intermediate -Filter ("*" + $ExpName) -Recurse -Directory).Fullname
    $SNreport = Join-Path $intermediate ($ExpName + "_SNreport.tsv")
    Copy-Item -Path (Join-Path $SNoutputDir ($ExpName + "_Report_pipeline_report (Normal).xls")) -Destination $SNreport
    $samples = Import-Csv $SNreport -Delimiter "`t" | Select-Object -ExpandProperty run_id -Unique
    if ($approach -eq "free") {$INreport = $SNreport}
}

## DIA analysis using DIA-NN


## Conversion of DDA results in PeptideGroups (.CSV) format
if ($mode -eq "DDA") {
    & "$conversions\DDAconversion.ps1" -InputFilePath $DDA -name $ExpName -OutputDirPath $intermediate
    $PDreport = Join-Path $intermediate ($ExpName + "_PDreport.tsv")
    $samples = Import-Csv $PDreport -Delimiter "`t" | Select-Object -ExpandProperty run_id -Unique
    if ($approach -eq "free") {$INreport = $PDreport}
}

## IS extraction for label approach
if ($approach -eq "label") {
    if ($mode -eq "DIA" -or $mode -eq "directDIA") {
        $report = $SNreport
    }
    elseif ($mode -eq "DDA") {
        $report = $PDreport
    }

    & "$conversions\ISextraction.ps1" -InputFilePath $report -ISpepFilePath $ISconc -name $ExpName -samples $samples -OutputDirPath $intermediate
    $ISreport = Join-Path $intermediate ($ExpName + "_ISpep_int.csv")
    $NLreport = Join-Path $intermediate ($ExpName + "_NL.tsv")
    $INreport = $NLreport
}

## Normalisation using Top3, Topall, iBAQ, APEX, NSAF, xTop or LFAQ
$methods = "top","all","iBAQ","APEX","NSAF","LFAQ","xTop"

## using aLFQ - R package
& "$conversions\Report_to_OpenSWATH.ps1" -InputFilePath $INreport -name $ExpName -OutputDirPath $intermediate
$aLFQrscript = "$normalisation\PeptoProtInference.R"
$aLFQargsList = "$aLFQrscript $intermediate $ExpName $fasta"
Start-Process -FilePath Rscript -ArgumentList $aLFQargsList -Wait

## using xTop - Python package
& "$conversions\Report_to_xTopinput.ps1" -InputFilePath $INreport -name $ExpName -samples $samples -OutputDirPath $intermediate
$xTopInput = $ExpName + "_xTop.csv"
$xTopPYscript = Join-Path ($env:Path -split ";" | Where-Object {$_ -match "xtop"}) "xTop_pipeline.py"
$xTopargsList = "`"$xTopPYscript`" $xTopInput"
Start-Process -FilePath python -ArgumentList $xTopargsList -WorkingDirectory $intermediate -Wait
$xTopOutputDir = (Get-ChildItem -Path $intermediate -Filter "*export" -Recurse -Directory).Fullname
Copy-Item -LiteralPath (Join-Path $xTopOutputDir ("[" + $xTopInput + "] Intensity xTop.csv")) -Destination (Join-Path $intermediate ($ExpName + "_prot_int_xTop.csv"))

## using LFAQ - C++ executables and LFAQ.py wrapper
$LFAQintermediate = "$intermediate\LFAQintermediate"
New-Item -ItemType Directory -Path $LFAQintermediate | Out-Null
& "$conversions\convert_input_LFAQ.ps1" -InputFilePath $INreport -name $ExpName -samples $samples -OutputDirPath $LFAQintermediate
$LFAQPYscript = "$normalisation\lfaq.py"
$LFAQexe = $env:Path -split ";" | Where-Object {$_ -match "LFAQ"}
$randomConcFile = "$normalisation\randomConcFile.csv"
foreach ($sam in $samples) {
    $samFile = Join-Path $LFAQintermediate ($ExpName + "_" + $sam + ".csv")
    $LFAQargsList = "$LFAQPYscript $samFile $fasta $LFAQintermediate `"$LFAQexe`" --IdentificationFileType `"PeakView`" --IdentifierOfStandardProtein `"P`" --StandardProteinsFilePath $randomConcFile"
    Start-Process -FilePath python -ArgumentList $LFAQargsList -Wait
    Rename-Item -Path "$LFAQintermediate\ProteinResultsExperimentOnlyOne.txt" -NewName ("ProteinResults_" + $sam + ".txt")
}
& "$conversions\convert_output_LFAQ.ps1" -InputFilesPath $LFAQintermediate -name $ExpName -samples $samples
Copy-Item -Path (Join-Path $LFAQintermediate ($ExpName + "_prot_int_LFAQ.csv")) -Destination (Join-Path $intermediate ($ExpName + "_prot_int_LFAQ.csv"))

## Absolute quantification using label or label-free approach
$AQPYscript = "$quantification\AbsQuant.py"
if ($approach -eq "label") {
    $AQargsList = "$AQPYscript --label `"label`" --name $ExpName --inDir $intermediate --sam $samples --met $methods --tot $totalProt --Sint $ISreport --Sconc $ISconc"
}
elseif ($approach -eq "unlabel") {
    $AQargsList = "$AQPYscript --label `"unlabel`" --name $ExpName --inDir $intermediate --sam $samples --met $methods --tot $totalProt --Sconc $ISconc"
}
elseif ($approach -eq "free") {
    $AQargsList = "$AQPYscript --label `"free`" --name $ExpName --inDir $intermediate --sam $samples --met $methods --tot $totalProt --fasta $fasta"
}
Start-Process -FilePath python -ArgumentList $AQargsList -Wait
$AQoutputDir = "$intermediate\Absolute_quantification"
foreach ($m in $methods) {
    $outputFileName = $ExpName + "_prot_conc_" + $m + ".csv"
    Copy-Item -Path (Join-Path $AQoutputDir $outputFileName) -Destination (Join-Path $OutputDir $outputFileName)
}
