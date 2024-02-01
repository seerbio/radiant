# Pythia DIA

<p align="center"><img alt="Pythia logo" src="./static/pythia-logo.png" style="width: 60%;"></p>

**Pythia DIA** is a search tool for bottom-up DIA-MS proteomics that offers:

* near-parity with DIA-NN IDs, using cutting-edge DIA processing algorithms
* extensive customizability via a parameters file
* multi-architecture support (Intel/ARM) for maximum flexibility and efficiency

Pythia DIA is usable as a Docker container or Linux binary.

## Running

Docker images are available for both Intel (`x86`) and ARM / Apple Silicon (`aarch64):

```
docker run --rm 718843040700.dkr.ecr.us-west-2.amazonaws.com/seer/pythia-dia:latest 
```

It is also possible to install binaries on Ubuntu/Debian:

On Intel (`x86`):

```
aws s3 cp s3://seer-buildfiles/pythiadia_1.1.5+amd64.deb .
apt-get install -y pythiadia_1.1.5+amd64.deb
/usr/local/bin/PythiaDIACpp/PythiaDIA --help
    ```

On ARM / Apple Silicon (`aarch64`):

```
aws s3 cp s3://seer-buildfiles/pythiadia_1.1.5+arm64.deb .
apt-get install -y pythiadia_1.1.5+arm64.deb
/usr/local/bin/PythiaDIACpp/PythiaDIA --help
```

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for information about building PythiaDIA from sources
or contributing to its development.
