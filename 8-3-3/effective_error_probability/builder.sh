#!/bin/bash
cd "$(dirname "$0")" || exit

rm -f ./build/CMakeCache.txt

# If using WSL might have to write cmake.exe instead
cmake -S . -B ./build -DCMAKE_BUILD_TYPE=Release
cmake --build ./build --config Release

read -r -p "Proceed to generate binaries (output)? [y/n] " proceed

if [[ "$proceed" == "Y" ]] || [[ "$proceed" == "y" ]];
then
    cd ./build/bin/ || exit
    rm -rf ./binary_data
    bash generator.sh
    echo "Output data generated successfully."
else
    exit 0
fi

read -r -p "Place/replace new data in SKR folder? [y/n] " proceed2

if [[ "$proceed2" == "Y" ]] || [[ "$proceed2" == "y" ]];
then
    src="./binary_data"
    dest="../../../secret_key_rate/build/effective_error_probability"

    # Check if the destination already has a folder with the same name
    if [[ -d "${dest}/$(basename $src)" ]]; then
        # If the folder exists, remove it
        rm -rf "${dest:?}/$(basename $src)"
    fi

    # Copy the source folder to the destination
    cp -r "$src" "$dest"

    echo "Folder copy/pasted successfully."
fi