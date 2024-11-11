#$LogFile = ".\ExpDataSets\ConFem0.log"
#if (Test-Path $LogFile) {
#    Remove-Item $LogFile
$dataset=$args[0]
$lines = Get-Content $dataset

Measure-Command{
foreach ($line in $lines) {
	Write-Host "$line"
	"$Line" | .\ConFemAll\ConFemAll.exe
}
} 
