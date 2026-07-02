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
docker run --rm seerbio/radiant-dia:latest 
```

You can download Ubuntu/Debian binaries from the [Releases page](https://github.com/seerbio/radiant/releases).
To install and use, you can then run:

On Intel (`x86`):

```
apt-get install -y radiantdia_<version>+amd64.deb
RadiantDIA --help
```

On ARM / Apple Silicon (`aarch64`):

```
apt-get install -y radiantdia_<version>+arm64.deb
RadiantDIA --help
```

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for information about building from sources
or contributing to development.
