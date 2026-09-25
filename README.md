<p align="center"><img alt="Radiant DIA logo" src="./static/radiant-logo-lockup.png" style="width: 40%;"></p>

**Radiant DIA**™ is a search tool for bottom-up DIA-MS proteomics that offers:

* cutting-edge DIA processing algorithms, with excellent speed, sensitivity, and quantitative precision/accuracy
* extensive customizability via a parameters file
* multi-architecture support (Intel/ARM) for maximum flexibility and efficiency

<p align="center"><img alt="Radiant DIA logo" src="./static/radiant-hero-graphic.png" style="width: 40%; margin-bottom: 0px;"></p>

## Using Radiant DIA

> [!TIP]
> For desktop users, it's recommended to use the [Radiant-Fulcrum GUI](https://github.com/seerbio/radiant-fulcrum-gui/)
> or the [Radiant-Fulcrum Docker container](https://github.com/seerbio/radiant-fulcrum-container/) (CLI).
> 
> For pipeline developers, or if you're interested in running just the Radiant DIA search engine, continue reading.

Radiant DIA is available as a Docker container or a Linux binary package.

Docker images are available for both Intel (`x86`) and ARM / Apple Silicon (`aarch64`):

```
docker run --rm 718843040700.dkr.ecr.us-west-2.amazonaws.com/seer/radiant-dia:latest 
```

You can download Ubuntu/Debian binaries from the [Releases page](https://github.com/seerbio/radiant/releases).
To install and use, you can then run:

On Intel (`x86`):

```
apt-get install -y radiantdia_<version>+amd64.deb
/usr/local/bin/radiant/RadiantDIA --help
```

On ARM / Apple Silicon (`aarch64`):

```
apt-get install -y radiantdia_<version>+arm64.deb
/usr/local/bin/radiant/RadiantDIA --help
```

## Fragment competition

One boolean controls fragment competition before calibration fitting and after
main-pass LDA scoring, before neural-network candidate selection, training, and
inference:

```toml
[MS2Params]
competitionEnabled = true
```

The default is `true`, including when the option is omitted. Set it to `false`
to skip competition at both stages. Fragment competition replaces the previous
shared-evidence filter at both locations; disabling it does not restore that
filter. Fulcrum consumes Radiant's output without a Python competition step.

Competition links equal-charge, coeluting candidates with neutral masses within
5 ppm and at least four distinct shared supported fragments within 20 ppm.
Supported fragments require positive intensity and trace cosine at least 0.5.
The apex separation must be no greater than half the narrower peak width.
Each connected group retains the candidate with the greatest unshared
intensity × cosine-squared evidence. Groups without unshared evidence are
rejected, while singletons are retained. Unshared means absent from the other
group members' extracted fragment lists. Equal evidence prefers higher LDA
score, then decoys, peptide sequence, apex, mass, and input order. The algorithm
does not use NN scores, q-values, or protein identities to select candidates.
The old `ionsSharedToReject` setting does not configure fragment competition.

Candidates rejected in the main pass remain excluded from NN processing,
TIMS target-decoy pair completion, and final results. Existing cached results
are unchanged; use a rebuilt executable and fresh native search outputs.

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for information about building from sources
or contributing to development.
