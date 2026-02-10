# Radiant DIA

<p align="center"><img alt="Radiant DIA logo" src="./static/radiant-logo.png" style="width: 40%;"></p>

**Radiant DIA**™ is a search tool for bottom-up DIA-MS proteomics that offers:

* cutting-edge DIA processing algorithms, with excellent speed, sensitivity, and quantitative precision/accuracy
* extensive customizability via a parameters file
* multi-architecture support (Intel/ARM) for maximum flexibility and efficiency

Radiant DIA is usable as a Docker container or Linux binary.

## Running

Docker images are available for both Intel (`x86`) and ARM / Apple Silicon (`aarch64`):

```
docker run --rm 718843040700.dkr.ecr.us-west-2.amazonaws.com/seer/radiant-dia:latest 
```

It is also possible to install binaries on Ubuntu/Debian:

On Intel (`x86`):

```
aws s3 cp s3://seer-buildfiles/radiantdia_<version>+amd64.deb .
apt-get install -y radiantdia_<version>+amd64.deb
/usr/local/bin/radiant/RadiantDIA --help
```

On ARM / Apple Silicon (`aarch64`):

```
aws s3 cp s3://seer-buildfiles/radiantdia_<version>+arm64.deb .
apt-get install -y radiantdia_<version>+arm64.deb
/usr/local/bin/radiant/RadiantDIA --help
```

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for information about building from sources
or contributing to development.
