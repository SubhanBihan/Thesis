Set-Location -Path $PSScriptRoot

Remove-Item -Path ".\build\CMakeCache.txt" -Force

cmake -S . -B .\build -G "MinGW Makefiles" -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release
cmake --build .\build --config Release

$proceed = Read-Host "Proceed to generate binaries (output)? [y/n]"

if($proceed.ToUpper() -eq 'Y'){
    Remove-Item -Path ".effective_error_probability\build\bin\binary_data" -Force -Recurse
    .\generator.ps1
    Write-Host "Output data generated successfully."
} else {
    exit 0
}

$proceed2 = Read-Host "Place/replace new data in SKR folder? [y/n] "

if($proceed2.ToUpper() -eq 'Y'){
    # Define source and destination directories
    $src = ".\binary_data"
    $dest = "..\..\..\secret_key_rate\build\effective_error_probability"

    # Check if the $dest variable is not null or empty
    if (-not [string]::IsNullOrEmpty($dest)) {
        # Construct the destination folder path
        $destFolder = Join-Path $dest (Split-Path $src -Leaf)

        # Check if the destination already has a folder with the same name
        if (Test-Path $destFolder) {
            # If the folder exists, remove it safely
            Remove-Item -Recurse -Force $destFolder
        }

        # Copy the source folder to the destination
        Copy-Item -Path $src -Destination $dest -Recurse

        Write-Host "Folder copy/pasted successfully."
    } else {
        # Error handling if $dest is null or empty
        Write-Host "Error: Destination path is null or empty."
        exit 1
    }
}