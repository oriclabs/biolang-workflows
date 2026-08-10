param(
    [Parameter(Mandatory = $true)]
    [string]$BioLangRepo,

    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$BioLangArguments
)

$source = (Resolve-Path -LiteralPath $BioLangRepo -ErrorAction Stop).Path
$executable = Join-Path $source "target/release/bl.exe"
$packages = Join-Path $source "packages"

if (-not (Test-Path -LiteralPath $executable -PathType Leaf)) {
    throw "release executable not found: $executable; run cargo build --release -p bl-cli"
}
if (-not (Test-Path -LiteralPath (Join-Path $packages "singlecell/biolang.toml") -PathType Leaf)) {
    throw "BioLang singlecell package not found under: $packages"
}
if ($BioLangArguments.Count -eq 0) {
    throw "provide BioLang arguments, for example: run workflows/single-cell/seurat_standard_workflow.bl"
}

$env:BIOLANG_PATH = $packages
if (-not (Test-Path Env:BIOLANG_GPU)) {
    $env:BIOLANG_GPU = "off"
}
& $executable @BioLangArguments
exit $LASTEXITCODE
