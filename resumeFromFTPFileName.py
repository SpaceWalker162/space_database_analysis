import os
import databaseTools as dbt
from ftplib import FTP_TLS
import otherTools as ot
import copy
import ast

ftpAddr ='cdaweb.gsfc.nasa.gov'
remoteDataDir = '/pub/data'
#downloadDataDir = 'F:\\data'
#databaseDir = '\\\\10.249.183.237\\data'
downloadDataDir = '/media/yufei/Elements/data'
databaseDirs = ['/media/yufei/server', '/media/yufei/Elements/data']
databaseDirs.append(downloadDataDir)
databaseDirs = set(databaseDirs)
verbose = True


#dataNameDict = {'ace': {'mag': {'level_2_cdaweb': 'mfi_h3'}, 'swepam': {'level_2_cdaweb': 'swe_h0'}}}

missionName = 'mms'
spacecraftNames = []
mmsNumbers = [1,2,3,4]
mmsNumbers = [1]
for mmsN in mmsNumbers:
    spacecraftNames.append(missionName+str(mmsN))
instrumentations = [['fpi', 'fast', 'l2', 'des-moms'], ['fpi', 'brst', 'l2', 'des-moms']]
#instrumentations = [['fpi', 'fast', 'l2', 'dis-moms', '2016']]
dataNameDict = {missionName: {}}
instrumentationsDict = ot.list2dict(instrumentations)
for spacecraftName in spacecraftNames:
    dataNameDict[missionName][spacecraftName] = instrumentationsDict
#print(dataNameDict)

databaseDataDict = {databaseDir: dataNameDict}

ftp = FTP_TLS(ftpAddr)
lgMess = ftp.login()
print(lgMess)
ftp.cwd(remoteDataDir)

#readFTPDirLog = 'readFTPDir.log'
readftpFileNameLog = 'readftpFileName.log'
#fFTPFile = open(readftpFileNameLog, 'w')
#fFTPDir = open(readFTPDirLog, 'w')
#logFileHandle = (fFTPDir, fFTPFile)
##logFileHandle = None
#fileInfoDict = dbt.readFTPFileInfoRecursively(ftp, ftpPath=dataNameDict, verbose=verbose, facts=['size'], logFileHandle=logFileHandle)
#fFTPDir.close()
#fFTPFile.close()

readftpFileNameLog = 'readftpFileName.log'
with open(readftpFileNameLog, 'r') as f:
    info = f.read().splitlines()
infoList = []
numberOfEntries = len(info)//2
for i in range(numberOfEntries):
    filePathList = info[i*2].split('/')
    infoDS = info[i*2+1]
    infoD = ast.literal_eval(infoDS)
    filePathList.append(infoD)
    infoList.append(filePathList)
fileInfoDict = ot.list2dict(infoList)

fileInfoDictLocal  = dbt.readFileInfoRecursively(path=databaseDataDict, verbose=verbose, facts=['size'])

toDownload = dbt.compareDictRecursively(fileInfoDict, fileInfoDictLocal[databaseDir])

toDownloadList_ = ot.dict2list(toDownload)
filePaths = [toD[:-2] for toD in toDownloadList_]
numberOfFiles = len(filePaths)
with open('toDownloadList.txt', 'w') as f:
    for filePath in filePaths:
        print(filePath, file=f)
for i in range(numberOfFiles):
    filePath = filePaths[i]
    l_ = [remoteDataDir]
    l_.extend(filePath)
    remoteFileName = '/'.join(l_)
    downloadToDir = os.path.join(downloadDataDir, *filePath[:-1])
    if not os.path.exists(downloadToDir):
         os.makedirs(downloadToDir)
    downloadedFileName = os.path.join(downloadDataDir, *filePath)
    _ = dbt.ftpDownload(ftp, remoteFileName, downloadedFileName, verbose=verbose)
    if verbose:
        print('{}/{}'.format(i, numberOfFiles))
print('good job!')
ftp.close()
