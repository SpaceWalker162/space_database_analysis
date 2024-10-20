import os
from ftplib import FTP_TLS
import copy
import logging
import space_database_analysis.otherTools as ot
import space_database_analysis.databaseTools as dbt
#import reconnecting_ftp

def downloadSPDF(downloadDatabaseDir, databaseDirs, dataNameDict, fileNamesSource='internet', logFileDir='', protocol='ftp', workerNumber=2, delay_between_downloads=2):
    '''
    Purpose:
        download data files from CDAWeb through FTP
    Parameters:
        fileNamesSource: if 'internet', read the FTP server to find file names for downloading; if not 'internet', the function will read from a stored file in the logFileDir
    <<<<<<<<<<<< Example <<<<<<<<<<

    import os
    import copy
    import sys
    import space_database_analysis.otherTools as ot
    import space_database_analysis.databaseTools.downloadSPDF as downloadSPDF
    import socket

    if sys.platform == 'linux':
        if socket.gethostname() == 'yufei-OptiPlex-5040-2':
            downloadDataDir = '/home/yufei/Documents/database' #This is where you want to store the downloaded files
        else:
            downloadDataDir = '/media/yufei/Elements/database' #This is where you want to store the downloaded files
        databaseDirs = ['/home/yufei/Documents/remoteDatabase'] #This is a list that contain all database. The files that exist in these databases will be omitted in downloading.
    elif sys.platform == 'win32':
        downloadDataDir = '..\\database'
        databaseDirs = ['\\\\10.249.183.237\\pub']

    #databaseDirs = []

    readInternetForFileNames = True # if True, the program reads from ftp the list of the files you would like to present in your downloadDataDir and databaseDirs after it ends. This suit the first run of the program. After a run with this parameter set True, the list of files will be stored locally for later use, such as a second run when the first run is not successful. In this case, this parameter should be set to False and so the program will read the locally stored list to save the time spent on reading ftp.

    #dataNameDict = {'ace': {'mag': {'level_2_cdaweb': 'mfi_h3'}, 'swepam': {'level_2_cdaweb': 'swe_h0'}}}

    #missionName = 'mms'
    #spacecraftNames = []
    #mmsNumbers = [1]
    #for mmsNumber in mmsNumbers:
    #    spacecraftNames.append(missionName+str(mmsNumber))
    #for i in range(4):
    #    spacecraftNames.append(missionName+str(i+1))
    #instrumentations = [['fpi', 'fast', 'l2', 'dis-moms'], ['fpi', 'brst', 'l2', 'dis-moms']]
    #instrumentations = [['mec', 'srvy', 'l2', 'epht89q'], ['mec', 'srvy', 'l2', 'epht89d']]
    #instrumentations = [['fgm', 'srvy', 'l2']]
    #instrumentations = [['edp', 'slow', 'l2', 'dce', '2019'], ['edp', 'fast', 'l2', 'dce', '2019']]
    #instrumentations = [['fpi', 'fast', 'l2', 'dis-partmoms', '2015'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2016'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2017'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2018']]
    #instrumentations = [['fpi', 'fast', 'l2', 'dis-partmoms', '2015'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2016'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2017'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2018']]
    #dataNameDict = {missionName: {}}
    #directory = ot.Directory(lis=instrumentations)
    #directory.generate_dic_from_lis()
    #instrumentationsDict = directory.dic
    #instrumentationsDict = {'feeps': {'srvy': {'l2': {'ion': {'2015':{}, '2016':{}, '2017': {}, '2018': {}}}}}}
    #for spacecraftName in spacecraftNames:
    #    dataNameDict[missionName][spacecraftName] = copy.deepcopy(instrumentationsDict)
    #
    dataNameDict = {'wind': {'3dp': {'3dp_plsp': {}}}}

    #missionName = 'themis'
    #spacecraftNames = ['thc', 'thd']
    #instrumentations = [['l2', 'fgm']]
    #dataNameDict = {missionName: {}}
    #instrumentationsDict = ot.list2dict(instrumentations)
    #for spacecraftName in spacecraftNames:
    #    dataNameDict[missionName][spacecraftName] = copy.deepcopy(instrumentationsDict)

    logFileDir = os.path.expanduser('~/Documents/MyFiles/works/project_working/temp/downloadSPDF')
    protocol = 'http'
    if readInternetForFileNames:
        fileNamesSource = 'internet'
    else:
        fileNamesSource = 'from_log'
    downloadSPDF.downloadSPDF(downloadDataDir, databaseDirs, dataNameDict, logFileDir=logFileDir, fileNamesSource=fileNamesSource, protocol=protocol)
    print('end')
    >>>>>>>>>>>>> Example >>>>>>>>>>>>>


    '''
    host ='cdaweb.gsfc.nasa.gov' #Address of the ftp
    remoteDataDir = '/pub/data' #ftp database directory
    verbose = True
    downloadBlocksize = 1024*32
    downloadDataDir = os.path.join(downloadDatabaseDir, 'data')

    os.makedirs(logFileDir, exist_ok=True)
    loggingHandlers = [logging.FileHandler(os.path.join(logFileDir, 'download.log'))]
    formatter = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s",
                                  "%Y-%m-%d %H:%M:%S")
    for handler in loggingHandlers:
        handler.setFormatter(formatter)
    if verbose:
        loggingHandlers.append(logging.StreamHandler())
    logging.basicConfig(level=logging.DEBUG, handlers=loggingHandlers)

    databaseDirs.append(downloadDatabaseDir)
    databaseDirs = set(databaseDirs)
    databaseDataDirs = [os.path.join(databaseDir_, 'data') for databaseDir_ in databaseDirs]

    if protocol == 'ftp':
        facts = ['size']
    elif protocol == 'http':
        facts = ['size-h']

    if fileNamesSource == 'internet':
        logging.info('reading {} for files to download'.format(host))
        fileInfoDict = dbt.readFTPHTTPFileInfoRecursively(host=host, commonPath=remoteDataDir, path=dataNameDict, verbose=verbose, facts=facts, logFileDir=logFileDir, protocol=protocol)
    else:
        fileInfoDict = dbt.loadFileInfoDict(logFileDir=logFileDir)

    localFileDictTree = ot.DictTree()
    for databaseDataDir in databaseDataDirs:
        databaseDataDict = {databaseDataDir: dataNameDict}
        fileInfoDictLocal = dbt.readFileInfoRecursively(path=databaseDataDict, verbose=verbose, facts=facts)
        localFileDictTree = localFileDictTree.union(fileInfoDictLocal[databaseDataDir])
    fileInfoDict = ot.DictTree(fileInfoDict)
    toDownload = fileInfoDict.difference(localFileDictTree)

    toDownloadList_ = toDownload.toList()
    filePaths = []
    total_size = 0 # in Bytes
    for toD in toDownloadList_:
        filePaths.append(toD[:-2])
        if toD[-2] == 'size-h':
            total_size += ot.sizeof(toD[-1])
        elif toD[-2] == 'size':
            total_size += float(toD[-1])
    print('number of files to download: {}'.format(len(filePaths)))
    print('total size to download: {}'.format(ot.sizeof_fmt(total_size)))
    ftpDownloader = dbt.FileDownloadCommander(remoteFileNames=filePaths, host=host, downloadRootPath=remoteDataDir, downloadedFileNames=None, destPath=downloadDataDir, verbose=verbose, keepDownloading=True, blocksize=downloadBlocksize, timeout=120, workerNumber=workerNumber, monitorInterval=10, protocol=protocol, delay_between_downloads=delay_between_downloads)
    ftpDownloader.processQueue()
    print('End of the Program')
