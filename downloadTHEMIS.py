import os
import databaseTools as dbt
from ftplib import FTP_TLS
import otherTools as ot
#import reconnecting_ftp
import ast

ftpAddr ='cdaweb.gsfc.nasa.gov'
remoteDataDir = '/pub/data'
downloadDataDir = '/media/yufei/Elements/data'
#downloadDataDir = '/mnt/data'
databaseDirs = ['/media/yufei/server']
#databaseDirs = ['\\\\10.249.183.237\\data']
databaseDirs = []
databaseDirs.append(downloadDataDir)
databaseDirs = set(databaseDirs)
verbose = True


#dataNameDict = {'ace': {'mag': {'level_2_cdaweb': 'mfi_h3'}, 'swepam': {'level_2_cdaweb': 'swe_h0'}}}

missionName = 'themis'
spacecraftNames = ['thc', 'thd']
instrumentations = [['l2', 'fgm']]
dataNameDict = {missionName: {}}
instrumentationsDict = ot.list2dict(instrumentations)
for spacecraftName in spacecraftNames:
    dataNameDict[missionName][spacecraftName] = instrumentationsDict
#

ftp = dbt.reconnectingFTP(ftpAddr)
#ftp = reconnecting_ftp.Client(ftpAddr, '', )
lgMess = ftp.login()
print(lgMess)
ftp.cwd(remoteDataDir)

readFTPDirLog = 'readFTPDir.log'
readftpFileNameLog = 'readftpFileName.log'
fFTPFile = open(readftpFileNameLog, 'w')
fFTPDir = open(readFTPDirLog, 'w')
logFileHandle = (fFTPDir, fFTPFile)
#logFileHandle = None
fileInfoDict = dbt.readFTPFileInfoRecursively(ftp, ftpPath=dataNameDict, verbose=verbose, facts=['size'], logFileHandle=logFileHandle)
fFTPDir.close()
fFTPFile.close()
ftp.quit()

def resumeFromFTPFileNameLog(readftpFileNameLog='readftpFileName.log'):
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
    return fileInfoDict

fileInfoDict = resumeFromFTPFileNameLog(readftpFileNameLog)

localFileDictTree = ot.DictTree()
for databaseDir in databaseDirs:
    databaseDataDict = {databaseDir: dataNameDict}
    fileInfoDictLocal = dbt.readFileInfoRecursively(path=databaseDataDict, verbose=verbose, facts=['size'])
    localFileDictTree = localFileDictTree.union(fileInfoDictLocal[databaseDir])

fileInfoDict = ot.DictTree(fileInfoDict)
toDownload = fileInfoDict.difference(localFileDictTree)

tolerateFailure = True
toDownloadList_ = toDownload.toList()
filePaths = [toD[:-2] for toD in toDownloadList_]
dbt.downloadFTPFromFileList(remoteFileNames=filePaths, ftpAddr=ftpAddr, ftpPath=remoteDataDir, downloadedFileNames=None, destPath=downloadDataDir, verbose=verbose, keepDownloading=True, tolerateFailure=tolerateFailure, blocksize=1024*32)
