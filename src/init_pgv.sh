#!/bin/bash
set -euxo pipefail

datafiles_json=$1
datadir=$2
publicdir=$3
settings=$4
pgv_dir=$5
build=$6

echo "Checking if PGV instance exists..."
if [[ -d "$pgv_dir" ]]; then
  echo "PGV instance already exists in $pgv_dir"
  echo "Updating local repo..."
  cd "$pgv_dir"
  git pull origin main
else
  echo "Downloading pgv..."
  git clone "https://github.com/mskilab-org/pgv.git" "$pgv_dir"
fi

echo "Setting up PGV instance with PGVdb data..."
if [[ "$datafiles_json" != "$pgv_dir/public/datafiles.json" ]]; then
    cp -f "$datafiles_json" "$pgv_dir/public/datafiles.json"
fi

if [[ "$settings" != "$pgv_dir/public/settings.json" ]]; then
    cp -f "$settings" "$pgv_dir/public/settings.json"
fi

origin_directory=$datadir
target_directory="$pgv_dir/public/data"

# Create target directory structure if it doesn't exist
mkdir -p "$target_directory"

copy_directory() {
  local origin_dir="$1"
  local target_dir="$2"

  # Loop through files and directories in the origin directory
  for item in "$origin_dir"/*; do
    if [[ -d "$item" ]]; then
      # If it's a directory, create the corresponding directory in the target directory if it doesn't exist
      local subdirectory="${item##*/}"
      local target_subdirectory="$target_dir/$subdirectory"
      if [[ ! -d "$target_subdirectory" ]]; then
        mkdir -p "$target_subdirectory"
      fi
      
      # Recursively copy the subdirectory
      copy_directory "$item" "$target_subdirectory"
    elif [[ -f "$item" ]]; then
      # If it's a file, create a symlink in the target directory to the file in the origin directory if it doesn't exist
      local file="${item##*/}"
      local target_file="$target_dir/$file"
      if [[ ! -e "$target_file" ]]; then
        ln -s "$item" "$target_file"
      fi
    fi
  done
}

# Call the copy_directory function with the origin and target directory paths
copy_directory "$origin_directory" "$target_directory"

echo "Launching PGV..."
cd $pgv_dir

if [ "$build" = "TRUE" ]; then
    echo "Building PGV..."
    yarn install

    # Check if gene files are available and if not then download
    if [ ! -s public/genes/hg19.arrow ]; then
        echo 'Downloading hg19.arrow'
        wget -P public/genes https://mskilab.s3.amazonaws.com/pgv/hg19.arrow
    fi
    if [ ! -s public/genes/hg19_chr.arrow ]; then
        echo 'Downloading hg19_chr.arrow'
        wget -P public/genes https://mskilab.s3.amazonaws.com/pgv/hg19_chr.arrow
    fi
    if [ ! -s public/genes/hg38.arrow ]; then
        echo 'Downloading hg38.arrow'
        wget -P public/genes https://mskilab.s3.amazonaws.com/pgv/hg38.arrow
    fi
    if [ ! -s public/genes/hg38_chr.arrow ]; then
        echo 'Downloading hg38_chr.arrow'
        wget -P public/genes https://mskilab.s3.amazonaws.com/pgv/hg38_chr.arrow
    fi

    if [ ! -s public/datafiles.json ]; then
        echo 'No datafiles.json was found so using the DEMO data'
        echo "To use your own data, don't forget to update datafiles.json"
        cp -f public/datafiles0.json public/datafiles.json
    fi

    rm -rf ./build

    yarn build
else
    ./start.sh
fi
