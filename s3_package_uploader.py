import boto3
import os
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

local_package: str = f"{package_name}-{arch}.deb"

assert os.path.exists(local_package), f"Did not find local package {local_package}"
print(f"Found package at {local_package}")

# Construct the final filename for the uploaded package
# TODO: pick final filenaming conv.
#now: str = str(time.time())
now = "test"

upload_package_name: str = f"{now}_{package_name}-{arch}.deb"

print(f"Will upload {local_package} to s3://{bucket_name}/{upload_package_name}")

s3.meta.client.upload_file(f"./{local_package}", bucket_name, upload_package_name)
