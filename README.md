*Pythia DIA* is a search tool for bottom-up DIA-MS proteomics that offers:

* near-parity with DIA-NN IDs, using cutting-edge DIA processing algorithms
* extensive customizability via a parameters file
* multi-architecture support (Intel/ARM) for maximum flexibility and efficiency

Pythia DIA is usable as a Docker container or Linux binary.

## Running

On Intel (`x86`):

* Docker:

    ```
    docker run --rm 718843040700.dkr.ecr.us-west-2.amazonaws.com/seer/pythia-dia:0.0.3-x86_64
    ```

* Binary:

    ```
    aws s3 cp s3://seer-buildfiles/pythiadia_0.0.3+amd64.deb .
    apt-get install -y pythiadia_0.0.3+amd64.deb
    /usr/local/bin/PythiaDIACpp/PythiaDIA --help
    ```

On ARM / Apple Silicon (`aarch64`):

* Docker:

    ```
    docker run --rm 718843040700.dkr.ecr.us-west-2.amazonaws.com/seer/pythia-dia:0.0.3-aarch64
    ```

* Binary:

    ```
    aws s3 cp s3://seer-buildfiles/pythiadia_0.0.3+arm64.deb .
    apt-get install -y pythiadia_0.0.3+arm64.deb
    /usr/local/bin/PythiaDIACpp/PythiaDIA --help
    ```

## Development

See [CONTRIBUTING.md]() for information about building PythiaDIA from sources
or contributing to its development.
