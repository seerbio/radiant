import boto3
import os
import time

access_key = os.getenv("AWS_ACCESS_KEY_ID", "")
secret = os.getenv("AWS_SECRET_ACCESS_KEY", "")
token = os.getenv("AWS_SESSION_TOKEN", "")
package_name = os.getenv("PACKAGE_NAME")

print(access_key)
print(secret)
print(token)
print(package_name)

session = boto3.Session(
    aws_access_key_id=access_key,
    aws_secret_access_key=secret,
    aws_session_token=token,
    region_name="us-west-2"
)

bucket_name: str = 'seer-buildfiles'

s3 = session.resource('s3')
my_bucket = s3.Bucket(bucket_name)

now: str = str(time.time())
now = "test"

upload_package_name: str = f"{now}_{package_name}"
s3.meta.client.upload_file(f"./{package_name}", bucket_name, upload_package_name)
