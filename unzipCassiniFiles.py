import os
import zipfile
import sys

path = sys.argv[1]
files = os.listdir(path)
for file in files:
    filePathAbs = os.path.join(path, file)
    baseName, ext = os.path.splitext(file)
    if ext.lower() == 'zip':
        fzip = zipfile.ZipFile(filePathAbs)
        dst = os.path.join(path, baseName)
        fzip.extractall(dst)
