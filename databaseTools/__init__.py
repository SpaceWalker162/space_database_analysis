__author__ = 'Yufei Zhou'
## This file contains functions and classes related to building a database

import tarfile
from ftplib import FTP_TLS
import numpy as np
import os
import json
import cdasws
from datetime import datetime, timedelta
import urllib3
from urllib.parse import urlparse, urljoin
from html.parser import HTMLParser
import time
import space_database_analysis.otherTools as ot
from pathlib import PurePath
import copy
import ast
#import progressbar
import functools
import threading
import socket
#import concurrent.futures
import logging
import ctypes
import queue
#from _thread import interrupt_main

_readftpFileNameLog ='readftpFileName.log'
_readFTPDirLog = 'readFTPDir.log'
## run this file under the destination directory 

dataTypeTransformationDict_PDS = {
        'TIME': 'CDF_EPOCH',
        'ASCII_REAL': 'CDF_FLOAT',
        }



class Database:

    def __init__(self, databasePaths=[]):
        '''
        Purposes:
            pass

        datasets_info stores detailed information for all datasets
        additional_datasets_info stores information needed by databaseUserTools to load data
        '''
        self.paths = databasePaths
        self.path = self.paths[0]
        self.datasets_info = {}
        self.additional_datasets_info = {}
        self.datasets_info_paths = []
        self.additional_datasets_info_paths = []
        self._metadata_path = 'metadata'
        self._database_dataset_info_file_name = os.path.join(self._metadata_path, 'datasets_info')
        self._database_additional_dataset_info_file_name = os.path.join(self._metadata_path, 'additional_datasets_info')
        for path in self.paths:
            datasets_info_file_path = os.path.join(path, self._database_dataset_info_file_name)
            additional_datasets_info_file_path = os.path.join(path, self._database_additional_dataset_info_file_name)
            if os.path.exists(datasets_info_file_path):
                self.datasets_info_paths.append(datasets_info_file_path)
            if os.path.exists(additional_datasets_info_file_path):
                self.additional_datasets_info_paths.append(additional_datasets_info_file_path)

    @classmethod
    def removeOutdatedFiles(cls, databasePaths, verbose=False):
        '''
        Purpose: SPDF from time to time updates the version of data files. Therefore in the database there could be multiple files being multiple versions of the same data. This function removes the files of older version if a newer version presents in the database.
        '''
        outdatedFilePaths = cls.getOutdatedFilePaths(databasePaths, verbose=verbose)
        print(outdatedFilePaths)
        for filePath in outdatedFilePaths:
            os.remove(filePath)
            logging.info('removed: {}'.format(filePath))

    @staticmethod
    def getOutdatedFilePaths(databasePaths, verbose=False):
        def getVersionFromFileName(fileName):
            '''
            Expected file name format example: bedatkjo-wekjfjk_jehwoi_jhfiwu10398_39823_v2.1.2.cdf
            '''
            fileNameWithoutExt, ext = os.path.splitext(fileName)
            fileNameComp = fileNameWithoutExt.split('_')
            fileNameBase = '_'.join(fileNameComp[:-1])
            if fileNameBase:
                pass
            else:
                return
            try:
                assert ext.lower() == '.cdf'
                assert fileNameComp[-1][0].lower() == 'v'
                fileVersion = [int(s_) for s_ in fileNameComp[-1][1:].split('.')]
            except:
                return
            return (fileNameBase, fileVersion, ext)
        ##
        def getDirDic(para):
            return tuple(para[:-4])

        outdatedFilePaths = []
        for databaseDir in databasePaths:
#            dataDir = os.path.join(databaseDir, 'data')
            dataDir = databaseDir
            fileInfoDict = readFileInfoRecursively(path=dataDir, verbose=verbose, facts=None)
            dirKeys = set(ot.doAtLeavesOfADict(dic=fileInfoDict, do=getDirDic))
            for ind, dirKey in enumerate(dirKeys):
                if dirKey:
                    dirPath = os.path.join(*dirKey)
                    dic = fileInfoDict.getSubDictTreeByKeys(dirKey)
                    fileNames = list(dic.keys())
                else:
                    continue
                while len(fileNames) > 0:
                    fileName = fileNames.pop()
                    ret = getVersionFromFileName(fileName)
                    if ret:
                        fileNameBase, fileVersion, fileExt = ret
                    else:
                        continue
                    for fileName_ in fileNames:
                        ret_ = getVersionFromFileName(fileName_)
                        if ret_:
                            fileNameBase_, fileVersion_, fileExt_ = ret_
                        else:
                            continue
                        if fileNameBase == fileNameBase_ and fileExt == fileExt_:
                            for vInd in range(len(fileVersion)):
                                if fileVersion[vInd] > fileVersion_[vInd]:
                                    fileNameOfOldVersion = fileName_
                                elif fileVersion[vInd] < fileVersion_[vInd]:
                                    fileNameOfOldVersion = fileName
                                else:
                                    fileNameOfOldVersion = None
                                if fileNameOfOldVersion:
                                    outdatedFilePaths.append(os.path.join(dataDir, dirPath, fileNameOfOldVersion))
                                    break
                            else:
                                print(fileName)
                                print(fileName_)
                                raise Exception('something went wrong')
        return outdatedFilePaths

    def define_add_info(self, datasetID):
        dID = datasetID
        dataset = {'Id': dID}
        if dID[:3] == 'MMS':
            dataset_path = os.path.join('mms', *dID.split('_')).lower()
            dataset.update({'dataset_path': dataset_path})
            if 'FPI_FAST' in dID:
                dataset.update({
                    'dataset_file_time_gap': '2 hour',
                    'file_naming_convention': '{datasetIDLower}_%Y%m%d%H%M%S',
                                })
            if 'EDP_FAST' in dID:
                dataset.update({'dataset_file_time_gap': '1 day'})
            if 'EDP_SLOW' in dID:
                dataset.update({'dataset_file_time_gap': '1 day'})
        elif dID[:4] == 'OMNI':
            dataset.update({
                'dataset_file_time_gap': '1 month',
                'file_naming_convention': '{datasetIDLower}_%Y%m%d',
                })
            if 'HRO' in dID:
                dataset_path = os.path.join('omni', 'omni_cdaweb', dID[5:].lower())
                dataset.update({'dataset_path': dataset_path})
        return dataset

    def define_dataset_file_naming_convention(self, datasetID, get_from_CDAWeb_metadata=False):
        '''
        Purpose: determine file naming convention for a dataset
        '''
        if get_from_CDAWeb_metadata: # this branch needs further development and is deprecated for the moment
            dataset = self.datasets_info[datasetID]
            for d_ in dataset['AdditionalMetadata']:
                if d_['Type'] == 'SPDF_SKT_CDF':
                    cdf_file_info_url = d_['value']
            o = urlparse(cdf_file_info_url)
            for path in self.paths:
                file_path = os.path.join(path, *o.path.split('/')[2:])
                if os.path.exists(file_path):
                    with open(file_path, 'r') as f:
                        dic = json.load(f)
                    cdfGlobalAttributes = dic[list(dic.keys())[0]].get('CDFglobalAttributes')
                    for item in cdfGlobalAttributes:
                        file_naming_convention_dic = item.get('File_naming_convention')
                        if file_naming_convention_dic:
                            file_naming_convention = file_naming_convention_dic[0]['0']
        return file_naming_convention

    def define_datasets_file_naming_convention(self, get_from_CDAWeb_metadata=False):
        '''
        Purpose: determine file naming conventions for all dataset
        '''
        # this function is under development
        self.loadDatasets_info()
        for datasetID in self.datasets_info.keys():
            dataset_naming_convention = self.define_dataset_file_naming_convention(datasetID, get_from_CDAWeb_metadata=get_from_CDAWeb_metadata)


    def make_additional_datasets_info(self, get_info_from_CDAWeb=False, dry_run=False):
        '''
        Purpose:
            Create a file storing support information about the datasets in the database
        Parameters:
            get_info_from_CDAWeb: if False, the function will get datasetID from the file defined by save_name

        '''
        datasets_dic = {}
        if get_info_from_CDAWeb:
            cdaswsObj = cdasws.CdasWs()
            datasets = cdaswsObj.get_datasets()
            for dataset in datasets:
                datasets_dic[dataset['Id']] = dataset
        else:
            for path in self.paths:
                save_path = os.path.join(path, self._database_dataset_info_file_name)
                if os.path.exists(save_path):
                    with open(save_path, 'r') as f:
                        datasets_dic = json.load(f)
                    break
        for path in self.paths:
            add_info_file_path = os.path.join(path, self._database_additional_dataset_info_file_name)
            if os.path.exists(add_info_file_path):
                with open(add_info_file_path, 'r') as f:
                    add_info = json.load(f)
                break
        else:
            add_info = {}
        for key, dataset in datasets_dic.items():
            dataset_add_info = self.define_add_info(key)
            if len(dataset_add_info) > 1:
                if key not in add_info:
                    add_info[key] = {}
                add_info[key].update(dataset_add_info)

        for path in self.paths:
            add_info_file_path = os.path.join(path, self._database_additional_dataset_info_file_name)
            with open(add_info_file_path, 'w') as f:
                json.dump(add_info, f)


    def make_datasets_info(self, get_info_from_CDAWeb=False, dry_run=False):
        '''
        Purpose:
            Create a file storing support information about the datasets in the database
        Parameters:
            get_info_from_CDAWeb: if False, the function will get info from a previously saved file.

        '''
        datasets_dic = {}
        if get_info_from_CDAWeb:
            cdaswsObj = cdasws.CdasWs()
            datasets = cdaswsObj.get_datasets()
            for dataset in datasets:
                datasets_dic[dataset['Id']] = dataset
        for path in self.paths:
            save_path = os.path.join(path, self._database_dataset_info_file_name)
            with open(save_path, 'w') as f:
                json.dump(datasets_dic, f)

    def _loadDataset_info(self, datasetID=None):
        '''
        This function is deprecated
        load dataset_info into self.datasets_info
        Parameters:
            datasetID: a string either being 'all' or t
        '''
        with open(self.datasets_info_paths[0], 'r') as f:
            self.datasets_info = json.load(f)

    def loadDatasets_info(self):
        if self.datasets_info_paths:
            with open(self.datasets_info_paths[0], 'r') as f:
                self.datasets_info = json.load(f)

    def load_additional_datasets_info(self):
        if self.additional_datasets_info_paths:
            with open(self.additional_datasets_info_paths[0], 'r') as f:
                self.additional_datasets_info = json.load(f)


class CDAWebHTMLParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.data_counts = 0
        self.objFact = {}

    def handle_starttag(self, tag, attrs):
        if tag == 'a' and attrs[0][0] == 'href':
            objName = attrs[0][1]
            if objName[-1] == '/':
                self.objName = attrs[0][1][:-1]
                self.objFact['type'] = 'dir'
            else:
                self.objName = attrs[0][1]
                self.objFact['type'] = 'file'

    def handle_endtag(self, tag):
        pass

    def handle_data(self, data):
        self.data_counts += 1
        if self.data_counts == 3:
            size = data.strip()
            if size == '-':
                self.objFact['type'] = 'dir'
            else:
                self.objFact['type'] = 'file'
#                if size[-3] == '.':
#                    size = str(ot.round_number(float(size[:-1]))) + size[-1]
                self.objFact['size-h'] = size


class HTTP(urllib3.PoolManager):
    def __init__(self, *argt, host=None, **argd):
        super().__init__(*argt, **argd)
        self.host = host
        self.current_path = PurePath('/')
        self.type = 'https'

    def cwd(self, path):
        self.current_path = PurePath(self.current_path, path)

    def mlsd(self, path, facts=None):
        url_ = self.type+'://'+self.host+str(PurePath(self.current_path, path))
        resp = self.request('GET', url_)
        objs = self._parse_cdaweb_html(resp.data.decode('ascii'))
        if facts is None:
            return [(obj[0], None) for obj in objs]
        else:
#            if 'size' in facts:
#                for objName, objFact in objs:
#                    if 'size' in objFact:
#                        sizeStr = objFact.get('size')
#                        size_unit_map = {'K': 1024, 'M': 1024**2}
#                        objFact['size'] = float(sizeStr[:-1]) * int(size_unit_map[sizeStr[-1]])
            return [(obj[0], {key: obj[1].get(key) for key in facts}) for obj in objs]

    @staticmethod
    def _parse_cdaweb_html(con):
        objs = []
        in_content = False
        for line in con.split('\n'):
            if "Parent Directory" in line:
                in_content = True
                continue
            if in_content:
                cdawebHTMLParser = CDAWebHTMLParser()
                cdawebHTMLParser.feed(line)
                if hasattr(cdawebHTMLParser, 'objName'):
                    objs.append((cdawebHTMLParser.objName, cdawebHTMLParser.objFact))
        return objs


def dirsInit(path, dictOfDirs):
    for supperDir, subDir in dictOfDirs.items():
        makeDirsFromDict(path, supperDir, subDir)


def makeDirsFromDict(path, key, value):
    path = os.path.join(path, key)
    if isinstance(value, str):
        subPath = os.path.join(path, value)
        if not os.path.exists(subPath):
            os.makedirs(subPath)
    elif isinstance(value, dict):
        for key_, value_ in value.items():
            makeDirsFromDict(path, key_, value_)
    else:
        raise Exception('Bad Dictionary')


def makeDictFromDirs(dictOfDirs, ls, path, cd=None, level=None):
    fileInfos = ls(path, ['type'])
    for fileName, fileInfo in fileInfos:
        if fileInfo['type'] == 'file':
            dictOfDirs[fileName] = 'file'
        elif fileInfo['type'] == 'dir':
            dictOfDirs[fileName] = {}
            path_ = '/'.join([path, fileName])
            makeDictFromDirs(dictOfDirs[fileName], ls, path_, cd, level)

def urlPrepare(start, end, dataset, action):
    website = 'https://csa.esac.esa.int/csa/aio/'
    info='?DATASET_ID='+dataset+'&START_DATE='+start.strftime("%Y-%m-%dT%H:%M:%SZ")+'&END_DATE='+end.strftime("%Y-%m-%dT%H:%M:%SZ")+'&DELIVERY_FORMAT=CDF&NON_BROWSER=1&CSACOOKIE=32125B6E3E58027D230F59657A5E015B6A1450722A465D302F0947606504526F3700546B2F5905330B07467527075C6A650B5C6C'
    if action == 'async':
        url = website + 'async-product-action' + info
    elif action == 'sync':
        url = website + 'product-action' + info
    return url

def syncDownload(start, end, dataset, fileName):
    'start and end in datetime form, dataset is DATASET_ID'
    url = urlPrepare(start, end, dataset, action='sync')
    fileSizeInM, timeCost, speed = download(url, fileName)
    return fileSizeInM, timeCost, speed

def asyncRequest(start, end, dataset): 
    'start and end in datetime form'
    url = urlPrepare(start, end, dataset, action='async')
    response = download(url, fileName=None)
    url = response.data.split()[16].decode('ascii')
    return url

def download(url, fileName='untitled'):
    downloadStart=datetime.now()
    http=urllib3.PoolManager()
    response = http.request('GET', url)
    timeCost=datetime.now()-downloadStart
    if fileName:
        with open(fileName, "wb") as file:
            file.write(response.data)
            fileSizeInM=os.stat(fileName).st_size/10**6
            speed=fileSizeInM/timeCost.total_seconds()
        return fileSizeInM, timeCost, speed
    else:
        return response

class ReconnectingFTP(FTP_TLS):
    def __init__(self, *argt, **argd):
        super().__init__(*argt, **argd)
        self.path = ''

    def login(self, *argt, **argd):
        ret = super().login(*argt, **argd)
        return ret

    def cwd(self, *argt, **argd):
        super().cwd(*argt, **argd)
        self.path = self.pwd()

    def reconnect(self, dt=20):
        self.close()
        reInterval = dt
        logging.info('''Something went wrong. Reconnecting FTP after {}s'''.format(reInterval))
        time.sleep(reInterval)
        self.connect()
        self.login()
        self.cwd(self.path)
        logging.info('''Reconnecting succeeded''')

    def connect(self, *argt, tryTimes=10, **argd):
        success = False
        while not success and tryTimes > 0:
            try:
                ret = super().connect(*argt, **argd)
                try:
                    self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
                    self.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPINTVL, 75)
                    self.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPIDLE, 60)
                except: logging.debug('socket setting failed')
                success = True
            except Exception as e:
                logging.info(e)
                tryTimes -= 1
        return ret

    def mlsd(self, *argt, tryTimes=10, **argd):
        success = False
        while not success and tryTimes > 0:
            try:
                objs = list(super().mlsd(*argt, **argd))
                success = True
            except Exception as e:
                logging.info(e)
                tryTimes -= 1
                logging.info('MLSD went wrong. Remaining try times: {}'.format(tryTimes))
                if tryTimes > 0:
                    self.reconnect()
        return objs

class FileDownloadCommander:
    def __init__(self, remoteFileNames, host=None, downloadRootPath=None, downloadedFileNames=None, destPath='.', verbose=False, keepDownloading=False, blocksize=32768, timeout=20, workerNumber=2, monitorInterval=10, protocol='ftp'):
        '''
        Parameters:
            remoteFileNames: a list of string, such as ['mms/mms1/.../fileName', ...], or a list of list of string
        '''
        self.remoteFileNames = remoteFileNames
        self.host = host
        self.downloadRootPath = downloadRootPath
        self.downloadedFileNames = downloadedFileNames
        self.destPath = destPath
        self.verbose = verbose
        self.keepDownloading = keepDownloading
        self.blocksize = blocksize
        self.timeout = timeout
        self.workerNumber = workerNumber
        self.monitorInterval = monitorInterval
        self.workers = []
        self.monitors = []
        self.preparePendingWorks()
        self.protocol = protocol

    def preparePendingWorks(self, source='remoteFileNames'):
        if source == 'remoteFileNames':
            self.pendingWorks = []
            rawWork = self.remoteFileNames
            for fInd, remoteFileName in enumerate(rawWork):
                if isinstance(remoteFileName, list):
                    remoteFilePath = remoteFileName
                    remoteFileName = str(PurePath(*remoteFileName))
                elif isinstance(remoteFileName, str):
                    remoteFilePath = remoteFileName.split('/')
                remoteFileFullName = str(PurePath(self.downloadRootPath, remoteFileName))
#                    remoteFilePath = remoteFileName.split('/')
#                    remoteFileFullPath = []
#                    remoteFileFullPath.append(self.downloadRootPath)
#                    remoteFileFullPath.extend(remoteFilePath)
#                    remoteFileFullName = '/'.join(remoteFileFullPath)
#                else:
#                    remoteFilePath = remoteFileName
#                    remoteFileFullPath = []
#                    remoteFileFullPath.append(self.downloadRootPath)
#                    remoteFileFullPath.extend(remoteFilePath)
#                    remoteFileFullName = '/'.join(remoteFileFullPath)
                if self.downloadedFileNames:
                    downloadedFileName = self.downloadedFileNames[fInd]
                    if isinstance(downloadedFileName, str):
                        downloadedFileFullName = os.path.join(self.destPath, downloadedFileName)
                    else:
                        downloadedFileFullName = os.path.join(self.destPath, *downloadedFileName)
                else:
                    downloadedFileFullName = os.path.join(self.destPath, *remoteFilePath)
                downloadToDir = os.path.split(downloadedFileFullName)[0]
                if not os.path.exists(downloadToDir):
                     os.makedirs(downloadToDir)
                self.pendingWorks.append((remoteFileFullName, downloadedFileFullName))
        elif source == 'failedWorks':
            self.pendingWorks = self.failedWorks
        self.worksN = len(self.pendingWorks)
        self.processedN = 0
        self.failedWorks = []

    def reportProgress(self):
        logging.info('progress: {}/{}'.format(self.processedN, self.worksN))
        logging.info('failed: {}'.format(len(self.failedWorks)))
        logging.info('active threads: {}'.format(threading.active_count()))

    def processQueue(self):
        def listenToWorkers():
            for i in reversed(range(len(self.workers))):
                worker = self.workers[i]
                if worker.is_alive():
                    if worker.finishedAWork.is_set():
                        self.processedN += 1
                        worker.finishedAWork.clear()
                        if len(self.pendingWorks) > 0:
                            work = self.pendingWorks.pop()
                            worker.pendingWorks.put(work)
                        else:
                            self.killWorkerAndMonitor(i)
                        self.reportProgress()
                else:
                    self.processedN += 1
                    self.failedWorks.append(worker.currentWork)
                    self.reportProgress()
                    self.killWorkerAndMonitor(i)

        def listenToMonitors():
            for i in reversed(range(len(self.monitors))):
                monitor = self.monitors[i]
                worker = self.workers[i]
                if monitor.badWorker.is_set():
                    self.processedN += 1
                    self.failedWorks.append(worker.currentWork)
                    self.killWorkerAndMonitor(i)

        while len(self.pendingWorks) > 0 or len(self.workers) > 0:
            if len(self.workers) < self.workerNumber:
                if len(self.pendingWorks) > 0:
                    self.addWorkerAndMonitor()
            if len(self.workers) > 0:
                listenToWorkers()
                listenToMonitors()
        if self.keepDownloading and len(self.failedWorks) > 0:
            self.preparePendingWorks(source='failedWorks')
            self.processQueue()

    def addWorkerAndMonitor(self):
        fInfo = queue.LifoQueue()
        worker = FileDownloader(host=self.host, fInfo=fInfo, verbose=self.verbose, blocksize=self.blocksize, timeout=self.timeout, protocol=self.protocol)
        work = self.pendingWorks.pop()
        worker.pendingWorks.put(work)
        worker.start()
        self.workers.append(worker)
        monitor = FileDownloadMonitor(fInfo=fInfo)
        monitor.start()
        self.monitors.append(monitor)

    def killWorkerAndMonitor(self, i):
        self.workers[i].kill()
        self.monitors[i].kill()
        self.workers[i].join()
        self.monitors[i].join()
        del self.workers[i]
        del self.monitors[i]
        logging.debug('A worker is fired')
        logging.info('active threads: {}'.format(threading.active_count()))


class FileDownloader(threading.Thread):
    def __init__(self, host=None, currentWork=None, fInfo=None, verbose=False, blocksize=32768, timeout=20, protocol='ftp'):
        self.host = host
        self.fInfo = fInfo
        self.fInd = -1
        self.currentWork = currentWork
        self.pendingWorks = queue.Queue()
        self.finishedAWork = threading.Event()
        self.verbose = verbose
        self.blocksize = blocksize
        self.timeout = timeout
        self.protocol = protocol
        super().__init__(target=self.process, daemon=True)

    def callback(self, cont):
        fInd, f = self.fInfo.get()
        f.write(cont)
        self.fInfo.put((fInd, f))

    def process(self):
        if self.protocol == 'ftp':
            logging.info('connecting ftp server at {}'.format(self.host))
            self.ftp = ReconnectingFTP(self.host, timeout=self.timeout)
            lgMess = self.ftp.login()
            logging.info(lgMess)
        while True:
            self.currentWork = self.pendingWorks.get()
            srcName, dstName = self.currentWork
            try:
                downloadStart = datetime.now()
                with open(dstName, 'wb') as f:
                    self.fInd += 1
                    self.fInfo.put((self.fInd, f))
                    if self.protocol == 'ftp':
                        logging.info('downloading {} from {}'.format(dstName, srcName))
                        self.ftp.retrbinary('RETR '+srcName, self.callback, blocksize=self.blocksize)
                    elif self.protocol == 'http':
                        url_ = urljoin('https://' + self.host, srcName)
                        logging.info('downloading {} from {}'.format(dstName, url_))
                        http = urllib3.PoolManager()
                        resp = http.request("GET", url_, preload_content=False)

                        for chunk in resp.stream(self.blocksize):
                            self.callback(chunk)
                        resp.release_conn()
                    fInd, f = self.fInfo.get()
            except KeyboardInterrupt:
                self._handle_exception_during_downloading(f)
                raise KeyboardInterrupt
            except Exception as e:
                logging.info(e)
                self._handle_exception_during_downloading(f)
                logging.info('failed in downloading {}'.format(srcName))
                raise
            timeCost = datetime.now()-downloadStart
            fileSizeInM = os.stat(dstName).st_size/10**6
            speed = fileSizeInM/timeCost.total_seconds()
            logging.info('{} downloaded at {}'.format(dstName, datetime.now()))
            logging.info(" size: {}M, time cost: {}, download speed: {:.3f}M/s" .format(fileSizeInM, timeCost, speed))
            self.finishedAWork.set()

    def _handle_exception_during_downloading(self, f):
        if hasattr(self, 'ftp'):
            self.ftp.abort()
            self.ftp.quit()
        dstName = f.name
        f.close()
        os.remove(dstName)
        logging.info("File Removed: {}".format(dstName))


    def get_id(self):
        # returns id of the respective thread
        if hasattr(self, '_thread_id'):
            return self._thread_id
        for id, thread in threading._active.items():
            if thread is self:
                return id

    def kill(self):
        thread_id = self.get_id()
        res = ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id,
              ctypes.py_object(SystemExit))
        if res > 1:
            ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id, 0)
            logging.info('Exception raise failure')
        else:
            logging.info('Exception raise success')
#        self.join()


class FileDownloadMonitor(threading.Thread):
    def __init__(self, fInfo, alarmingSpeed=8192*5, allowSlowSpeedN=3, monitorInterval=10):
        self.fInfo = fInfo
        self.alarmingSpeed = alarmingSpeed
        self.allowSlowSpeedN = allowSlowSpeedN
        self.monitorInterval = monitorInterval
        self.slowSpeedN = 0
        self.downloadedSize = 0
        self.fInd = -1
        self.badWorker = threading.Event()
        super().__init__(target=self.watch, daemon=True)

    def check(self, workerProgress):
        psChange = workerProgress - self.workerProgress
        self.workerProgress = workerProgress
        if psChange < self.alarmingSize/8:
            self.slowSpeedN += 1
        if self.slowSpeedN > self.allowSlowSpeedN:
            logging.info('Error: download speed slower than {}kB/s'.format(self.alarmingSpeed/8/1024))
            self.badWorker.set()

    def watch(self):
        fInd, f = self.fInfo.get()
        if fInd > self.fInd:
            self.fInd = fInd
            self.workerProgress = f.tell()
        self.fInfo.put((fInd, f))
        time.sleep(self.monitorInterval)
        self.alarmingSize = self.alarmingSpeed * self.monitorInterval
        while not self.badWorker.is_set():
            fInd, f = self.fInfo.get()
            if fInd > self.fInd:
                self.fInd = fInd
                self.workerProgress = f.tell()
            elif fInd == self.fInd:
                workerProgress = f.tell()
                self.check(workerProgress)
            self.fInfo.put((fInd, f))
            time.sleep(self.monitorInterval)

    def get_id(self):
        # returns id of the respective thread
        if hasattr(self, '_thread_id'):
            return self._thread_id
        for id, thread in threading._active.items():
            if thread is self:
                return id

    def kill(self):
        thread_id = self.get_id()
        res = ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id,
              ctypes.py_object(SystemExit))
        if res > 1:
            ctypes.pythonapi.PyThreadState_SetAsyncExc(thread_id, 0)
            logging.info('Exception raise failure')
        else:
            logging.info('Exception raise success')


class DownloadChecker:
    def __init__(self, filename, mode='wb', blocksize=8192, alarmingSpeed=8192*5, allowSlowSpeedN=3):
        self.filename = filename
        self.mode = mode
        self.blocksize = blocksize
        self.alarmingSpeed = alarmingSpeed
        self.allowSlowSpeedN = allowSlowSpeedN
        self.slowSpeedN = 0
        self.downloadedSize = 0
        self.lastCheckTime = datetime.now()
        self.lastCheckSize = 0
        self.f = open(self.filename, self.mode)
        self.filePos = self.f.tell()
        self.watch()
#        self.progressbar = progressbar.ProgressBar(max)

    def watch(self, reportInterval=20):
        alarmingSize = self.alarmingSpeed * reportInterval
        self.stop = threading.Event()
        def func_monitor(self):
            time.sleep(reportInterval)
            while not self.stop.is_set():
                pos = self.f.tell()
                posChange = pos - self.filePos
                self.filePos = pos
                if posChange < alarmingSize/8:
                    self.slowSpeedN +=1
                if self.slowSpeedN > self.allowSlowSpeedN:
#                    raise RuntimeError('Error: download speed slower than {}kB/s'.format(self.alarmingSpeed/8/1024))
                    print('Error: download speed slower than {}kB/s'.format(self.alarmingSpeed/8/1024))
                    sys.exit('Error: download speed slower than {}kB/s'.format(self.alarmingSpeed/8/1024))
#                    interrupt_main()
                time.sleep(reportInterval)

        monitor = threading.Thread(target=func_monitor, args=(self,), daemon=True)
        monitor.start()

    def check(self):
        dt = (datetime.now() - self.lastCheckTime).total_seconds()
        self.lastCheckTime = datetime.now()
        alarmingSize = self.alarmingSpeed * dt
        if self.downloadedSize - self.lastCheckSize < alarmingSize:
#            print('Warning: download speed slower than {}kB/s'.format(self.alarmingSpeed/1024))
            self.slowSpeedN += 1
        self.lastCheckSize = self.downloadedSize
        if self.slowSpeedN > self.allowSlowSpeedN:
            raise RuntimeError('Error: download speed slower than {}kB/s'.format(self.alarmingSpeed/8/1024))

    def write(self, cont):
        self.f.write(cont)
#        self.downloadedSize += self.blocksize
#        self.check()

    def close(self):
        self.stop.set()
        self.f.close()

    def __enter__(self):
#        self.f = open(self.filename, self.mode)
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.close()


def ftpDownload(ftp, remoteFileName, downloadedFileName, blocksize=8192, verbose=False):
    downloadStart = datetime.now()
    if verbose:
        print('downloading: ' + downloadedFileName)
    try:
        with DownloadChecker(filename=downloadedFileName, blocksize=blocksize) as f:
            ftp.retrbinary('RETR '+remoteFileName, f.write, blocksize=blocksize)
    except KeyboardInterrupt:
        os.remove(downloadedFileName)
        print("File Removed: {}".format(downloadedFileName))
        raise KeyboardInterrupt
    except EOFError:
        print('failed in downloading {}'.format(remoteFileName))
        raise
    timeCost = datetime.now()-downloadStart
    fileSizeInM = os.stat(downloadedFileName).st_size/10**6
    speed = fileSizeInM/timeCost.total_seconds()
    if verbose:
        print('downloaded at {}'.format(datetime.now()))
        print(" size: {}M, time cost: {}, download speed: {:.4f}M/s" .format(fileSizeInM, timeCost, speed))
    return fileSizeInM, timeCost, speed


def downloadSPDF(remotePath, localPath, host='cdaweb.gsfc.nasa.gov', years=None, monthsList=None):
    '''This function is deprecated. Use SyncSPDF instead'''
    ftp = FTP_TLS(host)
    ftp.login()
    ftp.cwd(remotePath)
    if years == None:
        years = ftp.nlst()
    numberOfFiles = 0
    numberOfDownloadedFiles = 0
    for i, year in enumerate(years):
        pathToYear = remotePath + '/' + year
        ftp.cwd(pathToYear)
        if monthsList is None:
            months = ftp.nlst()
        else:
            months = monthsList[i]
        for month in months:
            pathToMonth = pathToYear + '/' + month
            ftp.cwd(pathToMonth)
            files = ftp.nlst()
            localPath_ = os.path.join(localPath, year, month)
            if not os.path.exists(localPath_):
                os.makedirs(localPath_)
            for f in files:
                numberOfFiles +=1
                downloadedFileName = os.path.join(localPath_, f)
                if os.path.exists(downloadedFileName) and os.stat(downloadedFileName).st_size/10**3 > 1:
                    numberOfDownloadedFiles +=1
                else:
                    fileSizeInM, timeCost, speed = ftpDownload(ftp, f, downloadedFileName)
                    print('file name: ' + downloadedFileName)
                    print(' {}'.format(datetime.now()))
                    print(" size: {}M, time cost: {} download speed: {}M/s" .format(fileSizeInM, timeCost, speed))
                    numberOfDownloadedFiles +=1
    ftp.quit()
    return numberOfFiles, numberOfDownloadedFiles


def ftpAndLocalCWD(ftp, path):
    if isinstance(path, list):
        ftpPath = '/'.join(path)
        localPath = os.path.join(*path)
    elif isinstance(path, str):
        ftpPath = path
        localPath = path
    if not os.path.exists(localPath):
        os.makedirs(localPath)
    ftp.cwd(ftpPath)
    os.chdir(localPath)


def goToDataDir(ftp, remoteDataDir, localDataDir):
    ftp.cwd(remoteDataDir)
    os.chdir(localDataDir)


def syncSPDF(localDataDir, dataDict=None, host='cdaweb.gsfc.nasa.gov', dataPath=None, verbose=False, allowOverwrite=False, overwriteOption=None):
    '''
    Purpose:
        To download files from the spdf data server
    Parameters:
        localDataDir: can be either relative to the cwd or a absolute path
        dataDict: e.g. {'ace': {'mag': {'level_2_cdaweb': 'mfi_h3'}, 'swepam': {'level_2_cdaweb': 'swe_h0'}}}
        localDataDir: e.g. /mnt/data
        dataPath: a parameter in the form of a path that overwrite dataDict. For example, if mms/mms1/fpi/ is passed to this parameter, the function will download all files under this path from spdf. If path/to/file/fileName is passed to the parameter, the file will be downloaded.
    '''
    ftp = FTP_TLS(host)
    lgMess = ftp.login()
    print(lgMess)
    remoteDataDir = '/pub/data/'
    goToDataDir(ftp, remoteDataDir, localDataDir)
    localDataDir = os.getcwd()
    if dataPath is not None:
        pathInList = dataPath.split('/')
        if len(pathInList) == 1:
            ftpAndLocalCWD(ftp, pathInList)
            downloadFTPRecursively(ftp, verbose=verbose, allowOverwrite=allowOverwrite, overwriteOption=overwriteOption)
        elif len(pathInList) > 1:
            ftpAndLocalCWD(ftp, pathInList[:-1])
            files = ftp.mlsd(facts=['type'])
            for file in files:
                fileName = file[0]
                if pathInList[-1] == fileName:
                    fileType = file[1]['type']
                    if fileType == 'file':
                        downloadedFileName = fileName
                        fileSizeInM, timeCost, speed = ftpDownload(ftp, fileName, downloadedFileName, verbose=verbose, allowOverwrite=allowOverwrite, overwriteOption=overwriteOption)
                    elif fileType == 'dir':
                        ftp.cwd(file[0])
                        downloadFTPRecursively(ftp, verbose=verbose, allowOverwrite=allowOverwrite, overwriteOption=overwriteOption)
            else:
                raise Exception('not found file/dir')
        else:
            raise Exception('bad dataPath')
    elif dataDict is not None:
        datasets = ot.dict2list(dataDict)
        for dataset in datasets:
            goToDataDir(ftp, remoteDataDir, localDataDir)
            ftpAndLocalCWD(ftp, dataset)
            downloadFTPRecursively(ftp, verbose=verbose, allowOverwrite=allowOverwrite, overwriteOption=overwriteOption)
    ftp.quit()


def downloadFTPRecursively(ftp, depth=0, verbose=False, allowOverwrite=False, overwriteOption=None):
    '''
    Purpose:
        to download all files and directories recursively from the current workpath of ftp to the current workpath of the local machine
    Parameter:
        depth: the depth of current path from the first call of this function
    '''
    objs = ftp.mlsd(facts=['type', 'size'])
    for objName, objFact in objs:
        if objFact['type'] == 'dir':
            ftpAndLocalCWD(ftp, objName)
            downloadFTPRecursively(ftp, depth=depth+1, verbose=verbose)
        elif objFact['type'] == 'file':
            downloadedFileName = objName
            downloadQ = True
            if os.path.exists(downloadedFileName):
                existingFileSize = os.stat(downloadedFileName).st_size
                overwriteQ = False
                if allowOverwrite:
                    if overwriteOption is None:
                        overwriteQ = True
                    elif overwriteOption == 'sizeNotMatch':
                        if objFact['size'] != existingFileSize:
                            overwriteQ = True
                    elif overwriteOption == 'sizeSmall':
                        if existingFileSize/10**3 < 1:
                            overwriteQ = True
                if not overwriteQ:
                    downloadQ = False
            if downloadQ:
                _ = ftpDownload(ftp, objName, objName, verbose=verbose)
    ftpAndLocalCWD(ftp, '..')


def downloadFTPTree(ftp, ftpPath=None, localPath=None, verbose=False):
    '''
    Purpose:
        to download file and directories from ftpPath to localPath
    '''
    if ftpPath is not None:
        ftp.cwd(ftpPath)
    if localPath is not None:
        os.chdir(localPath)
    downloadFTPRecursively(ftp, verbose=verbose)


def readFTPHTTPFileInfoRecursively(client=None, host=None, commonPath=None, path='.', verbose=False, facts=None, logFileDir='', logFileHandle=None, protocol='ftp'):
    '''
    This function is a part of a collection of functions to download data from cdaweb.
    Parameters:
        commonPath: the common path leading from the host to the directory where 
        path: can be string, such as 'mms/mms1/fpi', or list of string, such as ['mms', 'mms1', 'fpi'], or list of list of string, such as [['mms', 'mms1', 'fpi'], ['mms', 'mms2', 'fpi']], or a dict, such as {'mms', {'mms1': 'fpi', 'mms2': 'fpi'}}
    '''
    logging.debug('in reading function')
    def func_base(para, **keywords):
        '''
        Parameters:
            **keywords: should contain ftp, verbose, facts, logFileDir, logFileHandle

        '''
        path_ = '/'.join(para[:-2])
        fileInfoDict = readFTPHTTPFileInfoRecursively(path=path_, **keywords)
        para[-1][para[-3]] = fileInfoDict
        return None

    if logFileHandle is None:
        readFTPDirLog = os.path.join(logFileDir, _readFTPDirLog)
        readftpFileNameLog = os.path.join(logFileDir, _readftpFileNameLog)
        fFTPFile = open(readftpFileNameLog, 'w')
        fFTPDir = open(readFTPDirLog, 'w')
        logFileHandle = (fFTPDir, fFTPFile)
        closeLogHandle = True
    else:
        closeLogHandle = False

    if protocol == 'ftp':
        if client is None:
            logging.info('connecting ftp server at {}'.format(host))
            client = ReconnectingFTP(host)
            lgMess = client.login()
            logging.info(lgMess)
            client.cwd(commonPath)
        func = functools.partial(func_base, client=client, verbose=verbose, facts=facts, logFileDir=logFileDir, logFileHandle=logFileHandle)
    elif protocol == 'http':
        if client is None:
            client = HTTP(host=host)
            client.cwd(commonPath)
        func = functools.partial(func_base, client=client, facts=facts, verbose=verbose, logFileDir=logFileDir, logFileHandle=logFileHandle)
    logging.debug('client defined according to protocol {}'.format(protocol))

    if isinstance(path, str): # this is the main branch
        fileInfoDict = {}
        if path[-1] == '/':
            path = path[:-1]
        if verbose:
            logging.info('reading {protocol} {path}'.format(protocol=protocol.upper(), path=path))
        fact0 = 'type'
        if facts is None:
            facts_ = [fact0]
        else:
            facts_ = facts.copy()
            facts_.append(fact0)
        facts_ = set(facts_)
        objs = client.mlsd(path=path, facts=facts_)
        for objName, objFact in objs:
            objAbsName = str(PurePath(path, objName))
#            objAbsName = '/'.join([path, objName])
            if objFact['type'] == 'dir':
                fileInfoDict_ = readFTPHTTPFileInfoRecursively(client=client, path=objAbsName, verbose=verbose, facts=facts, logFileDir=logFileDir, logFileHandle=logFileHandle)
                fileInfoDict[objName] = fileInfoDict_
            elif objFact['type'] == 'file':
                del objFact['type']
                if 'size' in facts:
                    objFact['size'] = int(objFact['size'])
                if 'size-h' in facts:
                    pass
                fileInfoDict[objName] = objFact
                if logFileHandle:
                    print(objAbsName, file=logFileHandle[1])
                    print(objFact, file=logFileHandle[1])
        if verbose:
            logging.info('reading {protocol} {path} done'.format(protocol=protocol, path=path))
        if logFileHandle:
            print(path, file=logFileHandle[0])
    elif isinstance(path, list):
        if isinstance(path[0], str):
            path = ot.list2dict([path])
            fileInfoDict = readFTPHTTPFileInfoRecursively(client=client, path=path, verbose=verbose, facts=facts, logFileHandle=logFileHandle)
        elif isinstance(path[0], list):
            path = ot.list2dict(path)
            fileInfoDict = readFTPHTTPFileInfoRecursively(client=client, path=path, verbose=verbose, facts=facts, logFileHandle=logFileHandle)
    elif isinstance(path, dict):
        fileInfoDict = copy.deepcopy(path)
        _ = ot.doAtLeavesOfADict(fileInfoDict, do=func)
    if closeLogHandle:
        fFTPDir.close()
        fFTPFile.close()
    return fileInfoDict


def loadFileInfoDict(logFileDir=''):
    readftpFileNameLog = os.path.join(logFileDir, _readftpFileNameLog)
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


def readFileInfoRecursively(path='.', verbose=False, facts='stats'):
    '''
    Parameters:
        path: can be string, such as 'mms/mms1/fpi', or list of string, such as ['mms', 'mms1', 'fpi'], or list of list of string, such as [['mms', 'mms1', 'fpi'], ['mms', 'mms2', 'fpi']], or a dict, such as {'mms', {'mms1': {'fpi': {}}, 'mms2': {'fpi': {}, 'fgm': {}}}}
    return:
        fileInfoDict: a dict in the form of {'mms': {'mms1': {'fpi': {filename: {'__info': infoOfTheFilefilename}, '__info': infoOfTheDirectorybrst}, '__info': infoOfTheDirectoryfpi}, '__info': infoOfTheDirectorymms1}, '__info': infoOfTheDirectorymms}
    '''
    def func_base(para, verbose, facts):
        path_ = os.path.join(*para[:-2])
        fileInfoDict = readFileInfoRecursively(path=path_, verbose=verbose, facts=facts)
        para[-1][para[-3]] = fileInfoDict
        return None
    func = functools.partial(func_base, verbose=verbose, facts=facts)

    if isinstance(path, str): # this is the major branch of the function
        fileInfoDict = ot.DictTree()
        if verbose:
            print('reading {}'.format(path))
        if os.path.exists(path):
            objs = os.scandir(path=path)
            for obj in objs:
                if obj.is_dir():
                    fileInfoDict_ = readFileInfoRecursively(path=obj.path, verbose=verbose, facts=facts)
                    fileInfoDict[obj.name] = fileInfoDict_
                elif obj.is_file():
                    stat = obj.stat()
                    if facts == 'stats':
                        objFact = stat
                    elif isinstance(facts, list):
                        objFact = {}
                        if 'size' in facts:
                            objFact['size'] = stat.st_size
                        if 'size-h' in facts:
                            ss_ = ot.sizeof_fmt(stat.st_size, fmt=False)[:-2].strip()
                            num = float(ss_[:-1])
                            if num >= 9.95:
                                sizeH = str(ot.round_number(num)) + ss_[-1]
                            else:
                                sizeH = str(ot.round_number(num, 1)) + ss_[-1]
                            objFact['size-h'] = sizeH
                    elif facts is None:
                        objFact = {'name': obj.name}
                    fileInfoDict[obj.name] = objFact
            if verbose:
                print('reading {} done'.format(path))
        else:
            if verbose:
                print('path does not exist: {}'.format(path))
    elif isinstance(path, list):
        if isinstance(path[0], str):
            path = ot.list2dict([path])
            fileInfoDict = readFileInfoRecursively(path=path, verbose=verbose, facts=facts)
        elif isinstance(path[0], list):
            path = ot.list2dict(path)
            fileInfoDict = readFileInfoRecursively(path=path, verbose=verbose, facts=facts)
    elif isinstance(path, dict):
        fileInfoDict = copy.deepcopy(path)
        _ = ot.doAtLeavesOfADict(fileInfoDict, do=func)
    else:
        raise Exception('Error: path input error')
    return fileInfoDict



def compareDictRecursively(d1, d2):
    d1Sole = {}
    for key, d1Value in d1.items():
        d2Value = d2.get(key)
        if d2Value:
            if d2Value == d1Value:
                continue
            elif isinstance(d2Value, dict):
                d1SoleKey = compareDictRecursively(d1Value, d2Value)
            else:
                d1SoleKey = d1Value
        else:
            d1SoleKey = d1Value
        d1Sole[key] = d1SoleKey
    return d1Sole


def selectDownloadFromFTP(fileNamesFTP, fileFactsFTP, fileNamesLocal, fileFactsLocal, ftpDir=None, localDir=None):
    '''
    Parameters:
        fileNamesFTP: a list of file names from ftp server. 
        ftpDir: a string of the ftp dir. The full name of a file from ftp server is composed of ftp directory and file name. If the parameter ftpDir is set as None, the parameter fileNamesFTP should be given by file name. If ftpDir is set as ftp directory, fileNamesFTP should also contain ftp directory so as to provide the full name.
   '''
    if ftpDir is not None:
        if ftpDir[-1] == '/':
            ftpDir = ftpDir[:-1]
        lenFTPDir = len(ftpDir) + 1
        fileInfosFTP = set([(fileName[lenFTPDir:], fileFact) for fileName, fileFact in zip(fileNamesFTP, fileFactsFTP)])
    else:
        fileInfosFTP = set([(fileName, fileFact) for fileName, fileFact in zip(fileNamesFTP, fileFactsFTP)])
    if localDir is not None:
        localDirList = os.path.split(localDir)
        lenLocalDir = len(localDirList)
        fileInfosLocal = []
        for fileName, fileFact in zip(fileNamesLocal, fileFactsLocal):
            fileNameList = os.path.split(fileName)
            fileName = os.path.join(*fileNameList[lenLocalDir:])
            fileInfo = (fileName, fileFact)
            fileInfosLocal.append(fileInfo)
        fileInfosLocal = set(fileInfosLocal)
    else:
        fileInfosLocal = set([(fileName, fileFact) for fileName, fileFact in zip(fileNamesLocal, fileFactsLocal)])
    fToDownloadInfo_ = fileInfosFTP.difference(fileInfosLocal)
    if ftpDir is None:
        fToDownloadInfo = fToDownloadInfo_ 
    else:
        fToDownloadInfo = [(ftpDir + '/' + element[0], element[1]) for element in fToDownloadInfo_]
    return fToDownloadInfo


def downloadFTPFromFileList(remoteFileNames, host=None, ftpPath=None, ftp=None, downloadedFileNames=None, destPath='.', verbose=False, tolerateFailure=True, keepDownloading=False, blocksize=32768, timeout=20):
    '''
    Parameters:
        remoteFileNames: a list of string, such as ['mms/mms1/.../fileName', ...], or a list of list of string
    '''
    if not keepDownloading:
        failedFiles = []
        if downloadedFileNames:
            failedDownloadedFileNames = []
        else:
            failedDownloadedFileNames = None
        returnedVar = [failedFiles, failedDownloadedFileNames]
        numberOfFiles = len(remoteFileNames)
        if ftp is None:
            closeFTP = True
            ftp = ReconnectingFTP(host, timeout=timeout)
            lgMess = ftp.login()
            print(lgMess)
            ftp.cwd(ftpPath)
        else: closeFTP = False
        for fInd, remoteFileName in enumerate(remoteFileNames):
            if isinstance(remoteFileName, str):
                remoteFilePath = remoteFileName.split('/')
            else:
                remoteFilePath = remoteFileName
                remoteFileName = '/'.join(remoteFilePath)
            if downloadedFileNames:
                downloadedFileName = downloadedFileNames[fInd]
                if isinstance(downloadedFileName, str):
                    downloadedFileFullName = os.path.join(destPath, downloadedFileName)
                else:
                    downloadedFileFullName = os.path.join(destPath, *downloadedFileName)
            else:
                downloadedFileFullName = os.path.join(destPath, *remoteFilePath)
            downloadToDir = os.path.split(downloadedFileFullName)[0]
            if not os.path.exists(downloadToDir):
                 os.makedirs(downloadToDir)
            try:
                _ = ftpDownload(ftp, remoteFileName, downloadedFileFullName, blocksize=blocksize, verbose=verbose)
            except Exception as e:
                print(e)
                if tolerateFailure:
                    failedFiles.append(remoteFileName)
                    if downloadedFileNames:
                        failedDownloadedFileNames.append(downloadedFileName)
                    ftp.reconnect()
                else: raise
            if verbose:
                print('progress: {}/{}'.format(fInd+1, numberOfFiles))
                print('failed: {}'.format(len(failedFiles)))
        if closeFTP:
            ftp.quit()
        return returnedVar
    else:
        downloadTurns = 0
        while len(remoteFileNames) > 0:
            downloadTurns +=1
            print('Start downloading turn: {}'.format(downloadTurns))
            remoteFileNames, downloadedFileNames = downloadFTPFromFileList(remoteFileNames=remoteFileNames, host=host, ftpPath=ftpPath, ftp=ftp, downloadedFileNames=downloadedFileNames, destPath=destPath, verbose=verbose, tolerateFailure=tolerateFailure, blocksize=blocksize, timeout=timeout)
            print('failed in downloading {} files'.format(len(remoteFileNames)))


def checkCompressedFiles(path, recordIn='./check.log', requestTo='./address/', fileType='r:gz'):
    'tarFileswithProblem is a list containing the file names of compressed files'
    tarFiles = os.listdir(path) 
    tarFilesWithProblem=[]
    smallCriterion = 3
    numberOfFailures=0
    for i, tarFileName in enumerate(tarFiles):
        tic = datetime.now()
        addressFile = tarFileName[:-6]+"txt"
        info=tarFileName.split('--')
        start=datetime(int(info[1]),int(info[2]),1,0,0,0)
        end = datetime(int(info[3]),int(info[4]),1,0,0,0)
        dataset=info[0]
        interval=end-start
        day=interval/interval.days
        timeTags=[[(start+i*day).strftime("%Y%m%d"),(start+(i+1)*day).strftime("%Y%m%d")] for i in range(0, interval.days)]
        lostTimeTags = []
        smallTimeTags = []
        tar = tarfile.open(path+tarFileName,fileType)
        tarInfos = tar.getmembers()
        for timeTag in timeTags:
            for tarInfo in tarInfos:
                if timeTag[0] in tarInfo.name and timeTag[1] in tarInfo.name and "5VPS" in tarInfo.name:
                    if tarInfo.size/10**3 < smallCriterion:
                        smallTimeTags.append(timeTag)
                    break
            else:
               lostTimeTags.append(timeTag)
        if lostTimeTags or smallTimeTags:
             numberOfFailures += 1
             tarFilesWithProblem.append(tarFileName)
             if recordIn:
                 with open(recordIn,'a') as file:
                     head = tarFileName+'\n'
                     file.write(head)
                     if lostTimeTags:
                         file.write(" "*4+"lost:\n")
                         for item in lostTimeTags:
                             file.write(" "*8+'%s\n' % item)
                     if smallTimeTags:
                         file.write(" "*4+"less than {}M:\n".format(smallCriterion))
                         for item in smallTimeTags:
                             file.write(" "*8+'%s\n' % item)
             if requestTo:
                 if not os.path.exists(requestTo):
                     os.makedirs(requestTo)
                 asyncRequest(start, end, dataset, requestTo+addressFile)
        print("checked: {}/{} lost: {}, time cost: {}".format(i+1, len(tarFiles), numberOfFailures, datetime.now()-tic))
    return tarFilesWithProblem


def isComplete(fileName, start, end, dataset, dataFormat='CDF', fileType='r:gz'):
    try:
        smallCriterion=0.1 # in k 
        interval=end-start
        day=interval/interval.days
        timeTags=[[(start+i*day).strftime("%Y%m%d"),(start+(i+1)*day).strftime("%Y%m%d")] for i in range(0, interval.days)]
        lostTimeTags = []
        smallTimeTags = []
        tar = tarfile.open(fileName, fileType)
        tarInfos = tar.getmembers()
        for timeTag in timeTags:
            for tarInfo in tarInfos:
                criteria = timeTag.copy()
                criteria.extend([dataset,dataFormat.lower()])
                if all(criterion in tarInfo.name for criterion in criteria):
                    if tarInfo.size/10**3 < smallCriterion:
                        smallTimeTags.append(timeTag)
                    break
            else:
               lostTimeTags.append(timeTag)
        if lostTimeTags or smallTimeTags:
            with open('error.log','a') as file:
                head = fileName+'\n'
                file.write(head)
                if lostTimeTags:
                    file.write(" "*4+"lost:\n")
                    for item in lostTimeTags:
                        file.write(" "*8+'%s\n' % item)
                if smallTimeTags:
                    file.write(" "*4+"less than {}M:\n".format(smallCriterion))
                    for item in smallTimeTags:
                        file.write(" "*8+'%s\n' % item)
            return False
        else:
            return True
    except:
        with open('error.log','a') as file:
            file.write(fileName+' exception\n')
        return False

def transformPDSdataToCDF(databaseDir, dataTransDict=None, stringCriteria=['FGM_KSM_1M'], infoFileExtension='.LBL', cdfFileNameFMT=None, recordPath=None):
    raise Exception('This function was moved to the submodule databaseDataFormat')
# start=datetime(2002,1,1)
# end=datetime(2002,7,1)
# dataset='C1_CP_CIS-CODIF_HS_H1_MOMENTS'
# isComplete('data/C1/C1_CP_CIS-CODIF_HS_H1_MOMENTS/compressed/C1_CP_CIS-CODIF_HS_H1_MOMENTS--2002--01--2002--07--cdf.tar.gz', start, end, dataset)
if __name__ == "__main__":
    pass
