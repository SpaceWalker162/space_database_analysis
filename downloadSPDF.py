import os
from ftplib import FTP_TLS
import copy
import logging
import space_database_analysis.otherTools as ot
import space_database_analysis.databaseTools as dbt
#import reconnecting_ftp

downloadDataDir = '/media/yufei/Elements1/data' #This is where you want to store the downloaded files
#downloadDataDir = '..\\data'

databaseDirs = ['/media/yufei/server'] #This is a list that contain all database. The files that exist in these databases will be omitted in downloading.
#databaseDirs = ['\\\\10.249.183.237\\data']
#databaseDirs = []

readFTP = False # if True, the program reads from ftp the list of the files you would like to present in your downloadDataDir and databaseDirs after it ends. This suit the first run of the program. After a run with this parameter set True, the list of files will be stored locally for later use, such as a second run when the first run is not successful. In this case, this parameter should be set to False and so the program will read the locally stored list to save the time spent on reading ftp.
ftpAddr ='cdaweb.gsfc.nasa.gov' #Address of the ftp
remoteDataDir = '/pub/data' #ftp database directory
verbose = True
downloadBlocksize = 1024*32

loggingHandlers = [logging.FileHandler('download.log')]
if verbose:
    loggingHandlers.append(logging.StreamHandler())
logging.basicConfig(level=logging.DEBUG, handlers=loggingHandlers)

#dataNameDict = {'ace': {'mag': {'level_2_cdaweb': 'mfi_h3'}, 'swepam': {'level_2_cdaweb': 'swe_h0'}}}

missionName = 'mms'
spacecraftNames = []
mmsNumbers = [1,2]
for mmsNumber in mmsNumbers:
    spacecraftNames.append(missionName+str(mmsNumber))
#for i in range(4):
#    spacecraftNames.append(missionName+str(i+1))
#instrumentations = [['fpi', 'fast', 'l2', 'dis-moms'], ['fpi', 'brst', 'l2', 'dis-moms']]
#instrumentations = [['mec', 'srvy', 'l2', 'epht89q'], ['mec', 'srvy', 'l2', 'epht89d']]
instrumentations = [['fgm', 'srvy', 'l2']]
#instrumentations = [['edp', 'slow', 'l2', 'dce']]
#instrumentations = [['fpi', 'fast', 'l2', 'dis-moms', '2016']]
dataNameDict = {missionName: {}}
instrumentationsDict = ot.list2dict(instrumentations)
for spacecraftName in spacecraftNames:
    dataNameDict[missionName][spacecraftName] = copy.deepcopy(instrumentationsDict)
#

#missionName = 'themis'
#spacecraftNames = ['thc', 'thd']
#instrumentations = [['l2', 'fgm']]
#dataNameDict = {missionName: {}}
#instrumentationsDict = ot.list2dict(instrumentations)
#for spacecraftName in spacecraftNames:
#    dataNameDict[missionName][spacecraftName] = copy.deepcopy(instrumentationsDict)

databaseDirs.append(downloadDataDir)
databaseDirs = set(databaseDirs)
if readFTP:
    fileInfoDict = dbt.readFTPFileInfoRecursively(ftpAddr=ftpAddr, ftpDir=remoteDataDir, ftpPath=dataNameDict, verbose=verbose, facts=['size'])
else:
    fileInfoDict = dbt.loadFileInfoDict()

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
print('number of files to download: {}'.format(len(filePaths)))
#dbt.downloadFTPFromFileList(remoteFileNames=filePaths, ftpAddr=ftpAddr, ftpPath=remoteDataDir, downloadedFileNames=None, destPath=downloadDataDir, verbose=verbose, keepDownloading=True, tolerateFailure=tolerateFailure, blocksize=downloadBlocksize)
ftpDownloader = dbt.FTPDownloadCommander(remoteFileNames=filePaths, ftpAddr=ftpAddr, ftpPath=remoteDataDir, downloadedFileNames=None, destPath=downloadDataDir, verbose=verbose, keepDownloading=True, blocksize=downloadBlocksize, timeout=20, workerNumber=2, monitorInterval=10)
ftpDownloader.processQueue()
print('End of the Program')
