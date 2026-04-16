#!/bin/bash
cd "$(dirname "$0")" || exit   # Sets path to that of this directory (where this file is saved)
rm -f ./build/CMakeCache.txt

# If using WSL might have to write cmake.exe instead
cmake -S . -B ./build -DCMAKE_BUILD_TYPE=Release
cmake --build ./build --config Release

read -r -p "Proceed to generate binaries (output)? [y/n] " proceed

if [[ "$proceed" == "Y" ]] || [[ "$proceed" == "y" ]];
then
    bash ./build/Release/generator.sh
else
    exit 0
fi

