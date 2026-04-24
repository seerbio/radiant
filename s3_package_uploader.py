import boto3
import os
import shutil
import subprocess
import time

package_name = os.getenv("PACKAGE_NAME")

arch = subprocess.run("echo \"$(dpkg-architecture | grep 'DEB_BUILD_ARCH=' | cut -d = -f 2)\"", shell=True, capture_output=True, encoding="utf8").stdout.strip()

print(f"Got {package_name=}")
print(f"Got {arch=}")

session = boto3.Session(
    # Assume creds are configured for the env. or instance
    region_name="us-west-2"
)

bucket_name: str = 'seer-buildfiles'

s3 = session.resource('s3')
my_bucket = s3.Bucket(bucket_name)

# Locate the local package file

local_package: str = f"{package_name}+{arch}.deb"
local_package_path: str = os.path.abspath(local_package)

assert os.path.exists(local_package_path), f"Did not find local package {local_package}"
print(f"Found package at {local_package_path}")

# Construct the final filename for the uploaded package

upload_package_name: str = f"{package_name}+{arch}.deb"

print(f"Will upload {local_package} to s3://{bucket_name}/{upload_package_name}")

s3.meta.client.upload_file(local_package_path, bucket_name, upload_package_name)

local_output_dir = os.getenv("LOCAL_DEB_OUTPUT_DIR")
if local_output_dir:
    assert os.path.isdir(local_output_dir), (
        f"LOCAL_DEB_OUTPUT_DIR is not a directory: {local_output_dir}"
    )

    copied_package_path = os.path.join(local_output_dir, upload_package_name)
    print(f"Will copy {local_package_path} to {copied_package_path}")
    shutil.copy2(local_package_path, copied_package_path)
