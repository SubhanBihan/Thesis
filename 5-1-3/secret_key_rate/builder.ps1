Set-Location -Path $PSScriptRoot

Remove-Item -Path ".\build\CMakeCache.txt" -Force

cmake -S . -B ./build -DCMAKE_BUILD_TYPE=Release
cmake --build ./build --config Release

$proceed = Read-Host "Proceed to generate binaries (output)? [y/n]"

if($proceed.ToUpper() -eq 'Y'){
    .\build\Release\generator.ps1
} else {
    exit 0
}
