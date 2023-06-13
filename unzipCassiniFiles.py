import os
import zipfile
import sys

lostFiles = []
lostFilesLogFileName = 'lost.log'
for path in sys.argv[1:]:
    print(path)
    files = os.listdir(path)
    for file in files:
        filePathAbs = os.path.join(path, file)
        baseName, ext = os.path.splitext(file)
        if ext.lower() == '.zip':
            if zipfile.is_zipfile(filePathAbs):
                print(filePathAbs)
                fzip = zipfile.ZipFile(filePathAbs)
                dst = os.path.join(path, baseName)
                fzip.extractall(dst)
                print(dst)
            else:
                lostFiles.append(file)
with open(lostFilesLogFileName, 'w') as f:
    print(*lostFiles, sep=',', file=f)
print("end")
