#!/usr/bin/env bash
# Build architecture-neutral parser assemblies for the Thermo RAW helper.
# Older parser distributions contain x64-only managed assemblies; this pinned
# build produces compatible AnyCPU assemblies for all supported targets.
set -euo pipefail

if [[ "$#" -ne 2 ]]; then
    printf 'Usage: %s DOTNET_EXECUTABLE CACHE_DIRECTORY\n' "$0" >&2
    exit 2
fi
dotnet_command="$1"
cache_dir="$2"
parser_commit=cf548e4a3fadcc796291a01fdbd088af986a5d72
vendor_commit=47eeda0d54a9ca28ace8e6548b7c01b6843207c2
vendor_version=8.0.42
source_sha=dbe5ab5d79ee90f843b46e7e6518aa0b792cb430eb5c8c5db762d33ddc14597a
raw_sha=af259d466e24227ba6706872330a0736d6fc2f552ad9d46a78a4697f0f62e1a6
data_sha=0b1d2fff2f81532e0c56471b22bc0e2c0c6685f37aa43b0b134ad82ae75ae444
output_dir="$cache_dir/parser-anycpu"
cache_key="$parser_commit:$vendor_commit:$vendor_version"
marker="$output_dir/.radiant-thermo-dependencies"

if [[ -f "$marker" && "$(cat -- "$marker")" == "$cache_key"
    && -f "$output_dir/ThermoRawFileParser.dll"
    && -f "$output_dir/ThermoFisher.CommonCore.RawFileReader.dll"
    && -f "$output_dir/ThermoFisher.CommonCore.Data.dll"
    && -f "$output_dir/log4net.dll" ]]; then
    printf 'Using cached AnyCPU Thermo RAW dependencies: %s\n' "$output_dir"
    exit 0
fi

mkdir -p -- "$cache_dir"

download_verified() {
    local url="$1" destination="$2" expected_sha="$3"
    if [[ -f "$destination" ]]; then
        local existing_sha
        existing_sha="$(sha256sum -- "$destination")"
        if [[ "${existing_sha%% *}" == "$expected_sha" ]]; then
            return
        fi
    fi
    curl -fsSL --retry 3 "$url" -o "$destination.partial"
    printf '%s  %s\n' "$expected_sha" "$destination.partial" | sha256sum --check -
    mv -- "$destination.partial" "$destination"
}

archive="$cache_dir/parser-source-$parser_commit.tar.gz"
download_verified \
    "https://codeload.github.com/compomics/ThermoRawFileParser/tar.gz/$parser_commit" \
    "$archive" "$source_sha"
vendor_url="https://raw.githubusercontent.com/thermofisherlsms/RawFileReader/$vendor_commit/Libs/NetCore/Net8"
download_verified \
    "$vendor_url/thermofisher.commoncore.rawfilereader.$vendor_version.nupkg" \
    "$cache_dir/thermofisher.commoncore.rawfilereader.$vendor_version.nupkg" "$raw_sha"
download_verified \
    "$vendor_url/thermofisher.commoncore.data.$vendor_version.nupkg" \
    "$cache_dir/thermofisher.commoncore.data.$vendor_version.nupkg" "$data_sha"

parser_source_dir="$cache_dir/parser-source-$parser_commit"
mkdir -p -- "$parser_source_dir"
tar -xzf "$archive" --strip-components=1 --directory "$parser_source_dir"
# This vendor release supplies genuine AnyCPU assemblies. The original 8.0.6
# dependency is x64-only; neither its PE headers nor its signed binaries are edited.
sed -i \
    "s/Include=\"ThermoFisher.CommonCore.RawFileReader\" Version=\"8.0.6\"/Include=\"ThermoFisher.CommonCore.RawFileReader\" Version=\"$vendor_version\"/" \
    "$parser_source_dir/ThermoRawFileParser.csproj"
"$dotnet_command" publish "$parser_source_dir/ThermoRawFileParser.csproj" \
    --configuration Release --framework net8.0 --no-self-contained \
    --output "$output_dir" \
    --source "$cache_dir" --source https://api.nuget.org/v3/index.json \
    -p:PlatformTarget=AnyCPU -p:Platform=AnyCPU -p:UseAppHost=false

cat > "$output_dir/radiant-thermo-dependencies.json" <<EOF
{
  "parser_source_commit": "$parser_commit",
  "parser_source_sha256": "$source_sha",
  "parser_platform_target": "AnyCPU",
  "vendor_source_commit": "$vendor_commit",
  "rawfilereader_version": "$vendor_version",
  "rawfilereader_package_sha256": "$raw_sha",
  "data_package_sha256": "$data_sha"
}
EOF
printf '%s\n' "$cache_key" > "$marker"
printf 'Prepared AnyCPU Thermo RAW dependencies: %s\n' "$output_dir"
