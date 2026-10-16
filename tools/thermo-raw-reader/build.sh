#!/usr/bin/env bash
set -euo pipefail

source_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd -- "$source_dir/../.." && pwd)"
parser_dir="${THERMO_RAW_FILE_PARSER_DIR:-}"
output_dir="$repo_dir/build/bin/thermoraw"
sdk_dir=""
host_arch="$(uname -m)"
case "$host_arch" in
    x86_64|amd64) runtime_id=linux-x64 ;;
    aarch64|arm64) runtime_id=linux-arm64 ;;
    *) runtime_id="" ;;
esac

for argument in "$@"; do
    case "$argument" in
        --parser-dir=*) parser_dir="${argument#*=}" ;;
        --output=*) output_dir="${argument#*=}" ;;
        --dotnet=*) sdk_dir="${argument#*=}" ;;
        --runtime=*) runtime_id="${argument#*=}" ;;
        --help|-h)
            printf 'Usage: %s [--parser-dir=THERMORAWFILEPARSER_DIR] [--output=DIR] [--dotnet=SDK_DIR] [--runtime=RID]\n' "$0"
            printf 'Default runtime: linux-x64 on x86_64; linux-arm64 on aarch64/arm64 (AWS Graviton).\n'
            printf 'Use --runtime=linux-arm64 to cross-publish from an x86_64 host.\n'
            printf 'Pinned AnyCPU dependencies are prepared automatically; --parser-dir can select an existing distribution.\n'
            exit 0 ;;
        *) printf 'Unknown argument: %s\n' "$argument" >&2; exit 2 ;;
    esac
done
if [[ -z "$runtime_id" ]]; then
    printf 'No default .NET runtime for host architecture %s; provide --runtime=RID.\n' "$host_arch" >&2
    exit 2
fi
if [[ "$runtime_id" != linux-arm64 && -n "$parser_dir" ]]; then
    if [[ ! -f "$parser_dir/ThermoRawFileParser.dll" ]]; then
        printf 'The selected parser directory must contain ThermoRawFileParser.dll.\n' >&2
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
if [[ -n "$sdk_dir" ]]; then
    dotnet_command="$sdk_dir/dotnet"
    if [[ ! -x "$dotnet_command" ]]; then
        printf 'Configured .NET SDK directory does not contain an executable dotnet: %s\n' "$sdk_dir" >&2
        exit 2
    fi
else
    if ! dotnet_command="$(command -v dotnet)"; then
        printf 'A .NET 8 SDK is required. Run install-build-deps-ubuntu.sh first or use --dotnet=SDK_DIR.\n' >&2
        exit 2
    fi
fi
if ! "$dotnet_command" --list-sdks 2>/dev/null | awk '$1 ~ /^8\./ { found = 1 } END { exit !found }'; then
    printf 'The selected dotnet does not provide a .NET 8 SDK: %s\n' "$dotnet_command" >&2
    exit 2
fi
export DOTNET_CLI_TELEMETRY_OPTOUT=1
export DOTNET_NOLOGO=1
if [[ "$runtime_id" == linux-arm64 || -z "$parser_dir" ]]; then
    dependencies_dir="$repo_dir/build/thermoraw/dependencies"
    bash "$source_dir/prepare-dependencies.sh" "$dotnet_command" "$dependencies_dir"
    parser_dir="$dependencies_dir/parser-anycpu"
fi
"$dotnet_command" publish "$source_dir/RadiantThermoRawReader.csproj" \
    --configuration Release --runtime "$runtime_id" \
    --output "$output_dir" \
    -p:RuntimeFrameworkVersion=8.0.0 \
    -p:TargetLatestRuntimePatch=false \
    "-p:ThermoRawFileParserDir=$parser_dir"
for license in LICENSE THERMO_LICENSE; do
    if [[ -f "$parser_dir/$license" ]]; then
        cp -- "$parser_dir/$license" "$output_dir/$license.ThermoRawFileParser"
    fi
done
if [[ -f "$parser_dir/radiant-thermo-dependencies.json" ]]; then
    cp -- "$parser_dir/radiant-thermo-dependencies.json" "$output_dir/"
fi
printf 'Direct Thermo RAW reader: %s/RadiantThermoRawReader\n' "$output_dir"
printf 'Target runtime: %s (build host: %s)\n' "$runtime_id" "$host_arch"
