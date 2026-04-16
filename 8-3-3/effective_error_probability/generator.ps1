# param([double]$value)
Set-Location -Path $PSScriptRoot\effective_error_probability\build\bin

# Write-Host "Reencoding error = $value"

# & .\main.exe "1erasure" $value
# & .\main.exe "flag" $value

$values = @(1e-4, 3e-4, 5e-4, 1e-3, 2e-3)
# $values = @(1e-4, 3e-4, 5e-4, 1e-3)

foreach($value in $values){
    Write-Host "Reencoding error = $value"

    & .\main.exe "1erasure" $value
    # & .\main.exe "flag" $value
}