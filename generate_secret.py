import os
from django.core.management.utils import get_random_secret_key

print("Generating a secret key and writing secrets.py.")
key = get_random_secret_key()

with open(os.path.join(os.getcwd(), "my_secrets", "secrets.py"), "w") as f:
    f.write("SECRET_KEY = \"%s\"\r\n" % key)