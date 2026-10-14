#!/usr/bin/env bash
set -euo pipefail

source_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd -- "$source_dir/../.." && pwd)"
parser_dir="${THERMO_RAW_FILE_PARSER_DIR:-}"
output_dir="$repo_dir/build/bin/thermoraw"
sdk_dir="$repo_dir/build/thermoraw/dotnet"
sdk_version=8.0.426
host_arch="$(uname -m)"
case "$host_arch" in
    x86_64|amd64) runtime_id=linux-x64; sdk_arch=x64 ;;
    aarch64|arm64) runtime_id=linux-arm64; sdk_arch=arm64 ;;
    *) runtime_id=""; sdk_arch="" ;;
esac

for argument in "$@"; do
    case "$argument" in
        --parser-dir=*) parser_dir="${argument#*=}" ;;
        --output=*) output_dir="${argument#*=}" ;;
        --dotnet=*) sdk_dir="${argument#*=}" ;;
        --runtime=*) runtime_id="${argument#*=}" ;;
        --help|-h)
            printf 'Usage: %s --parser-dir=THERMORAWFILEPARSER_DIR [--output=DIR] [--dotnet=SDK_DIR] [--runtime=RID]\n' "$0"
            printf 'Default runtime: linux-x64 on x86_64; linux-arm64 on aarch64/arm64 (AWS Graviton).\n'
            printf 'Use --runtime=linux-arm64 to cross-publish from an x86_64 host.\n'
            printf 'ARM64 builds prepare pinned AnyCPU dependencies automatically; --parser-dir is required for x64.\n'
            exit 0 ;;
        *) printf 'Unknown argument: %s\n' "$argument" >&2; exit 2 ;;
    esac
done
if [[ -z "$runtime_id" ]]; then
    printf 'No default .NET runtime for host architecture %s; provide --runtime=RID.\n' "$host_arch" >&2
    exit 2
fi
if [[ "$runtime_id" != linux-arm64 ]]; then
    if [[ -z "$parser_dir" || ! -f "$parser_dir/ThermoRawFileParser.dll" ]]; then
        printf 'Provide --parser-dir=DIR for an unpacked .NET 8 ThermoRawFileParser distribution.\n' >&2
        exit 2
    fi
    parser_dir="$(cd -- "$parser_dir" && pwd)"
    for assembly in ThermoRawFileParser ThermoFisher.CommonCore.Data ThermoFisher.CommonCore.RawFileReader log4net; do
        if [[ ! -f "$parser_dir/$assembly.dll" ]]; then
            printf 'Required parser distribution assembly missing: %s\n' "$assembly.dll" >&2
            exit 2
        fi
    done
fi
if [[ ! -x "$sdk_dir/dotnet" ]]; then
    if command -v dotnet >/dev/null; then
        dotnet_command="$(command -v dotnet)"
    else
        mkdir -p -- "$sdk_dir" "$repo_dir/build/thermoraw"
        installer="$repo_dir/build/thermoraw/dotnet-install.sh"
        curl -fsSL https://dot.net/v1/dotnet-install.sh -o "$installer"
        sdk_arguments=(--version "$sdk_version" --install-dir "$sdk_dir" --no-path)
        if [[ -n "$sdk_arch" ]]; then
            sdk_arguments+=(--architecture "$sdk_arch")
        fi
        bash "$installer" "${sdk_arguments[@]}"
        dotnet_command="$sdk_dir/dotnet"
    fi
else
    dotnet_command="$sdk_dir/dotnet"
fi
export DOTNET_CLI_TELEMETRY_OPTOUT=1
export DOTNET_NOLOGO=1
if [[ "$runtime_id" == linux-arm64 ]]; then
    arm_dependencies_dir="$repo_dir/build/thermoraw/arm-dependencies"
    bash "$source_dir/prepare-arm64-dependencies.sh" "$dotnet_command" "$arm_dependencies_dir"
    parser_dir="$arm_dependencies_dir/parser-anycpu"
fi
"$dotnet_command" publish "$source_dir/RadiantThermoRawReader.csproj" \
    --configuration Release --runtime "$runtime_id" --self-contained true \
    --output "$output_dir" \
    "-p:ThermoRawFileParserDir=$parser_dir"
for license in LICENSE THERMO_LICENSE; do
    if [[ -f "$parser_dir/$license" ]]; then
        cp -- "$parser_dir/$license" "$output_dir/$license.ThermoRawFileParser"
    fi
done
if [[ -f "$parser_dir/radiant-arm64-dependencies.json" ]]; then
    cp -- "$parser_dir/radiant-arm64-dependencies.json" "$output_dir/"
fi
printf 'Direct Thermo RAW reader: %s/RadiantThermoRawReader\n' "$output_dir"
printf 'Target runtime: %s (build host: %s)\n' "$runtime_id" "$host_arch"
