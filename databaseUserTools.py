__author__ = 'Yufei Zhou'

import numpy as np
#import constants as con
import cdflib    # see github.com/MAVENSDC/cdflib
#import tarfile
from datetime import datetime
from datetime import timedelta
import os
import sys
import subprocess
from cycler import cycler
from itertools import combinations
import json
from pprint import pprint
from scipy.signal import butter, lfilter, freqz
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
import scipy.interpolate as interpolate
import scipy.special
from databaseTools import *
import logging
from wolframclient.language import wl
from wolframclient.evaluation import WolframLanguageSession
import otherTools as ot
import databaseTools as dbt
import dataAnalysisTools as dat

'''
<A, B> means either A or B
[A, B] means both A and B
'''


multispacecraftMissions = ['cluster', 'mms']
piledInstrumentation = ['5VPS', 'CIS']
instrumentationFileUnderYear = [('cluster', '5VPS'), ('cluster', 'CIS'), ('ace', 'swepam'), ('ace', 'mag')]
instrumentationFileUnderMonth = None 
instrumentationDiviedIntoMonth = ['FULL', 'mms']
missionInfo = {
        'mms': {
            'epochType': 'CDF_TIME_TT2000',
            'multispacecraftMission': True},
        'cluster': {
            'epochType': 'CDF_EPOCH',
            'multispacecraftMission': True},
        'ace': {
            'epochType': 'CDF_EPOCH',
            'multispacecraftMission': False}
        }


class Database:
    '''

    '''
    def __init__(self, path):
        self.path = path

    def loadInfo(self, fileDictInfoFile=None):
        if fileDictInfoFile is None:
            self.fileDict = dbt.readFileInfoRecursively(path=self.path)
            self.variableDict = aaaaa
        else:
            self.fileDict = json.load(fileDictInfoFile)
            self.variableDict = aaaaa

    def printContentTree(self, root=None, year=False, month=False, variables=False, file=sys.stdout):
        def cutBranches(contentDict, until='year'):
            if until == 'year':
                keyLen = 4
            elif until == 'month':
                keyLen = 2
            for key, item in contentDict.items():
                if len(key) == keyLen:
                    try:
                        int(key)
                        return None
                    except:
                        pass
                if isinstance(item, dict):
                    branch = cutBranches(item, until=until)
                    contentDict[key] = branch
            return contentDict

        if root is None:
            root = self.path
        contentDict = dbt.readFileInfoRecursively(path=root)
        if year:
            if month:
                pass
            else:
                datasets = cutBranches(contentDict, until='month')
        elif variables:
            datasets = cutBranches(contentDict, until='year')
        else:
            datasets = cutBranches(contentDict, until='year')
            contentDictTree = ot.DictTree(datasets)
        contentDictTree.print()


class Instrumentation:
    '''
    Purpose:
        A object
    properties:
        mission:
        spacecraft:
        instrumentation: a list in terms of the directory names at all levels below spacecraft and above year or files. For example, ['fgm', 'brst', 'l2']
        instrumentationPath: a string representing the path to the instrumentation after spacecraft, say, 'swepam/level_2_cdaweb/swe_h0'.
    methods:
    '''
    def __init__(self, mission=None, spacecraft=None, instrumentation=None, instrumentationPath=None, path=None, splitedPath=None, workDataDir=None, epochType='CDF_EPOCH'):
        self.workDataDir = workDataDir
        if isinstance(path, str):
            path = os.path.normpath(path)
            splitedPath = path.split(os.sep)
        if isinstance(splitedPath, list):
            path = os.path.join(*splitedPath)
            splitedPathCopy = splitedPath.copy()
            mission = splitedPathCopy[0]
            if missionInfo[mission]['multispacecraftMission']:
                del splitedPathCopy[0]
                spacecraft = splitedPathCopy.pop(0)
            else:
                spacecraft = splitedPathCopy.pop(0)
            instrumentation = splitedPathCopy
        self.mission = mission # a string
        self.spacecraft = spacecraft # a string
        if isinstance(instrumentation, str):
            self.instrumentation = [instrumentation]
            self.instrumentationPath = instrumentation
        elif isinstance(instrumentation, list):
            self.instrumentation = instrumentation
            self.instrumentationPath = os.path.join(*instrumentation)
        self.path = path
        self.splitedPath = splitedPath
        self.epochType = epochType


class DataFile(Instrumentation):
    '''
    '''
    def __init__(self, mission=None, spacecraft=None, instrumentation=None, instrumentationPath=None, pathToDataset=None, splitedPathToDataset=None, dateTime=None, fileName=None, instrumentationObj=None, cdfFile=None, workDataDir=None, fileSize='allSize', epochType='CDF_EPOCH', silence=True, filePath=None, defineCDFFileQ=True):
        '''
        Parameter:
            workDataDir: the path to the directory of the working database
            splitedPathToDataset: a list in terms of the path to the dataset. For example, ['mms', 'mms1', 'fgm', 'brst', 'l2']
            epoch: the epoch at which the destined file contain data
            fileSize: a string in the form of "<num" or ">num", where num is an integer whose unit is bytes. This parameter define either that the size of the file should be less than num bytes, or that its size greater than num bytes. The default value for this parameter is "allSize" which allow for all file size.
            silence: to print nothing

        '''
        if isinstance(instrumentationObj, Instrumentation):
            mission = instrumentationObj.mission
            spacecraft = instrumentationObj.spacecraft
            instrumentation = instrumentationObj.instrumentation
            instrumentationPath = instrumentationObj.instrumentationPath
            pathToDataset = instrumentationObj.path
            splitedPathToDataset = instrumentationObj.splitedPath
        super().__init__(mission=mission, spacecraft=spacecraft, instrumentation=instrumentation, instrumentationPath=instrumentationPath, path=pathToDataset, splitedPath=splitedPathToDataset, workDataDir=workDataDir, epochType=epochType)
        self.pathToDataset = self.path
        del self.path
#        if epoch is not None:
#            self.epoch = epoch
#            self.dateTime = datetime(*cdflib.cdfepoch.breakdown(epoch)[:6])
#        elif dateTime is not None:
        self.dateTime = dateTime
#            self.epoch = cdflib.cdfepoch.compute_epoch(ot.datetime2list(dateTime))
        self.fileName = fileName
        self.filePath = filePath
        self.cdfFile = cdfFile
#        for crit in instrumentationFileUnderYear:
#            if all([name in self.pathToDataset for name in crit]):
##                self.pathToFile = self.pathToDataset
#                self.pathToFile = os.path.join(self.pathToDataset, self.dateTime.strftime('%Y'))
#                break
#        else:
#            self.pathToFile = os.path.join(self.pathToDataset, self.dateTime.strftime('%Y'), self.dateTime.strftime('%m'))
        if not self.filePath:
            self.findFilePath(size=fileSize, silence=silence)
        if self.filePath:
            if defineCDFFileQ:
                self.defineCDFFile()


    def defineCDFFile(self):
        self.cdfFile = cdflib.CDF(self.filePath)

    def findFilePath(self, size='allSize', silence=True):
        '''
        Purpose:
            This function is to find the file stored in self.workDataDir.
        '''
        start = self.dateTime
        end = start + timedelta(seconds=1)
        absolutePathToDataset = os.path.join(self.workDataDir, self.pathToDataset)
        criteria = []
        timeTag = None
        logging.debug("mission: {}".format(self.mission))
        if self.mission == 'mms':
            searchMethod = 'general'
            criteria.extend(self.instrumentation.copy())
            criteria.append(start.strftime("%Y%m%d"))
            if 'fgm' in self.instrumentationPath:
                if 'brst' in self.instrumentationPath:
                    timeTag = start
            elif 'fpi' in self.instrumentationPath:
                if 'brst' in self.instrumentationPath:
                    timeTag = start
        elif self.mission == 'cluster':
            if 'cis-hia' in self.instrumentationPath:
                searchMethod = 'general'
                infoList_ = self.instrumentation[1].split('_')
                mode = infoList_[2][:-4]
                sensitivity = infoList_[3][0]
                dataName_ = infoList_[4]
                if dataName_ == 'phasespacedens':
                    dataNameAbb = 'psd'
                elif dataName_ == 'diffenergyflux':
                    dataNameAbb = 'pef'
                restructedStrCriterion = '_'.join(['cis-hia', sensitivity+'s', mode, 'ions', dataNameAbb])
                criteria.append(restructedStrCriterion)
                criteria.append('v2022')
                criteria.append(start.strftime("%Y%m%d"))
            else:
                searchMethod = 'ClusterCAA'
                interval = [start, end]
                logging.debug('ClusterCAA searching criteria (instrumentation)')
                logging.debug(self.instrumentation)
                criteria.extend(self.instrumentation.copy()[1:])
                criteria.append(start.strftime("%Y%m%d"))
        elif self.mission == 'ace':
            searchMethod = 'general'
            criteria.append(start.strftime("%Y%m%d"))
        else:
            searchMethod = 'general'
            raise Exception('mission not defined!')
        logging.debug("Looking for data files in: {}".format(absolutePathToDataset))
        if searchMethod == 'general':
            fileNames, absolutePathToFile = self.findFileNamesUnderDataset(absolutePathToDataset, func=findFileNames, strings=criteria, size=size, timeTag=timeTag)
#            = findFileNames(absolutePathToFile)
        elif searchMethod == 'ClusterCAA':
            fileNames, absolutePathToFile = self.findFileNamesUnderDataset(absolutePathToDataset, func=findFileNamesInInterval, interval=interval, strings=self.instrumentation)
#            fileNames = findFileNamesInInterval(absolutePathToFile, interval=interval, strings=self.instrumentation)
        numFiles = len(fileNames)
        if numFiles == 1:
            self.fileName = fileNames[0]
            self.filePath = os.path.join(absolutePathToFile, self.fileName)
            if not silence:
                print(self.fileName)
        elif numFiles == 0:
            logging.warning("No file was found.")
#            raise Exception('No file was found.')
        elif numFiles > 1:
            logging.debug("More than one files were found:" + ("\n{}"*numFiles).format(*fileNames))
#            raise Exception('More than one files were found.')


    def findFileNamesUnderDataset(self, absolutePathToDataset, func, **para):
        fileNames = func(absolutePathToDataset, **para)
        absolutePathToFile = absolutePathToDataset
        if len(fileNames) == 0:
            absolutePathToDatasetYear = os.path.join(absolutePathToDataset, self.dateTime.strftime('%Y'))
            absolutePathToFile = absolutePathToDatasetYear
            logging.debug("Not found. Looking for data files in: {}".format(absolutePathToDatasetYear))
            fileNames = func(absolutePathToDatasetYear, **para)
            if len(fileNames) == 0:
                absolutePathToDatasetYearMonth = os.path.join(absolutePathToDatasetYear, self.dateTime.strftime('%m'))
                absolutePathToFile = absolutePathToDatasetYearMonth
                logging.debug("Not found. Looking for data files in: {}".format(absolutePathToDatasetYearMonth))
                fileNames = func(absolutePathToDatasetYearMonth, **para)
        return fileNames, absolutePathToFile

class Spacecraft:
    '''
    Properties:
        mission:
        name:
        instrumentationAll:
        data: a dictionary where data is stored
        workDataDir: the main database
        workDataDirsBak: a list of paths. In case main database does not contain the wanted data, the program could look for data in these backup database.
    '''
    def __init__(self, mission=None, name=None, instrumentationAll=None, data=None, dataCleaned=None, workDataDir=None, workDataDirsBak=None, epochType=None):
        self.mission = mission
        self.name = name
        self.instrumentationAll = instrumentationAll
        if epochType is None:
           epochType = missionInfo[self.mission]['epochType']
        self.epochType = epochType
        if data is None:
            self.data = {}
        elif isinstance(data, dict):
            self.data = data
        else:
            raise Exception('wrong data input format')
        if dataCleaned is None:
            self.dataCleaned = {}
        elif isinstance(dataCleaned, dict):
            self.dataCleaned = dataCleaned
        else:
            raise Exception('wrong dataCleaned input format')
        self.workDataDir = workDataDir
        self.workDataDirsBak = workDataDirsBak

    def loadData(self, datetimeRange, instrumentation=None, instrumentationRetrivingName=None, variables=None, variableRetrivingNames=None, datasetAndVariables=None, instrumentationVariablesWithRetrivingName=None, cleanData=False, tStepPecentageCriterion=0.9, lowpassCutoff=None, inPlace=False, gapThreshold=None, minNumberOfPoints=None, returnShiftQ=False, copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak=True):
        '''
        Parameters:
            instrumentation: a list in terms of the directory names at all levels below spacecraft and above year or files. For example, ['fgm', 'brst', 'l2']
            datasetAndVariables: this parameter should not be used by user. 
            instrumentationVariablesWithRetrivingName: a dictionary in terms of the path to the datasets. Its leaf value is a list of variable names. The names are defined by the corresponding cdfFile. For example, {'Cluster': {'C1' : {'C1_CP_FGM_FULL': ['FGM', ('time_tags__C1_CP_FGM_FULL', 't'), ('B_vec_xyz_gse__C1_CP_FGM_FULL', 'B')]}}, 'mms': {'mms1': {'fgm': {'brst': {'l2': ['fgmBrst', ('Epoch', 't'), ('Btotal', 'BTotal')]}}}}}. Please note that 'Epoch' and 'Btotal' are improvised. It may also be a list of lists, with each sublist in the form of ['Cluster', 'C1', 'C1_CP_FGM_FULL', ['FGM', ('time_tags__C1_CP_FGM_FULL', 't'), ('B_vec_xyz_gse__C1_CP_FGM_FULL', 'B')]]. To retrieve data, for example, of 'B_vec_xyz_gse__C1_CP_FGM_FULL', use C1.data['FGM']['B']
        Note:
            Necessary parameters include: datetimeRange, <[<[instrumentation, variables], datasetsAndVariables>, instrumentationRetrivingName, variableRetrivingNames], instrumentationVariablesWithRetrivingName>
            To retrieve data, use Spacecraft.data[instrumentationName][variableRetrivingName]
        '''
        if instrumentationVariablesWithRetrivingName:
            # recursive branch of loadData
            if isinstance(instrumentationVariablesWithRetrivingName, dict):
                instrumentationVariablesWithRetrivingNameList = ot.dict2list(instrumentationVariablesWithRetrivingName)
                logging.debug(instrumentationVariablesWithRetrivingNameList)
            if isinstance(instrumentationVariablesWithRetrivingName, list):
                instrumentationVariablesWithRetrivingNameList = instrumentationVariablesWithRetrivingName
            for ins in instrumentationVariablesWithRetrivingNameList:
                dataset = ins[:-1]
                logging.debug(dataset)
                instrumentationRetrivingName = ins[-1][0]
                variables = [varPair[0] for varPair in ins[-1][1:]]
                variableRetrivingNames = [varPair[1] for varPair in ins[-1][1:]]
                datasetAndVariables = dataset
                datasetAndVariables.append(variables)
                self.loadData(datetimeRange=datetimeRange, instrumentationRetrivingName=instrumentationRetrivingName, variableRetrivingNames=variableRetrivingNames, datasetAndVariables=datasetAndVariables, cleanData=cleanData, tStepPecentageCriterion=tStepPecentageCriterion, lowpassCutoff=lowpassCutoff, inPlace=inPlace, gapThreshold=gapThreshold, minNumberOfPoints=minNumberOfPoints, returnShiftQ=returnShiftQ)
        else:
            # major branch of loadData
            variableRetrivingNames = variableRetrivingNames.copy()
#            if 't' not in variableRetrivingNames:
#                variables.insert(0, 0)
#                variableRetrivingNames.insert(0, 't')
            if datasetAndVariables is None:
                if isinstance(instrumentation, str):
                    instrumentation = [instrumentation]
                if self.mission in multispacecraftMissions:
                    datasetAndVariables = [self.mission, self.name, *instrumentation, variables]
                else:
                    datasetAndVariables = [self.name, *instrumentation, variables]
                variablesWithRetrivingNames = list(zip(variables, variableRetrivingNames))
            else:
                variablesWithRetrivingNames = list(zip(datasetAndVariables[-1], variableRetrivingNames))
            datasetsAndVariables = [datasetAndVariables]
            data_ = readData(self.workDataDir, datasetsAndVariables, datetimeRange=datetimeRange, epochType=self.epochType, workDataDirsBak=self.workDataDirsBak, copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak=copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak)
            logging.info(data_[0].keys())
            dataInDict = dict([(varRetName, data_[0][varName]) for varName, varRetName in variablesWithRetrivingNames])
            self.data.update({instrumentationRetrivingName: dataInDict})
            if cleanData:
                self.cleanData(instrumentation=instrumentationRetrivingName, tStepPecentageCriterion=tStepPecentageCriterion, lowpassCutoff=lowpassCutoff, inPlace=inPlace, gapThreshold=gapThreshold, minNumberOfPoints=minNumberOfPoints, returnShiftQ=returnShiftQ)

    def cleanData(self, instrumentation, tStepPecentageCriterion=0.9, lowpassCutoff=None, inPlace=False, tName='t', gapThreshold=None, minNumberOfPoints=None, returnShiftQ=None):
        '''
        Purpose: the time series of the dataset under an instrumentation may be irregular.
        Parameters:
            instrumentation: self-explanatory
            lowpassCutoff:
            inPlace: if True, spacecraft.data[instrumentation] is replaced by the cleaned data. Otherwise, the cleaned data is available at spacecraft.dataCleaned[instrumentation]
        '''
        self.dataCleaned.update({instrumentation: {}})
        t = self.data[instrumentation][tName]
        tLen = len(t)
        shapes = []
        variableNames = []
        variables = []
        for key, var in self.data[instrumentation].items():
            if key == tName:
                pass
            else:
                shape = var.shape
                if len(shape) > 1:
                    shapes.append(var.shape[1:])
                else:
                    shapes.append([1])
                variableNames.append(key)
                variables.append(var.reshape((tLen, -1)))
        data_ = np.concatenate(variables, axis=1)
        packedReturns = dat.dataFillAndLowPass(t, data_, axis=0, tStepPecentageCriterion=tStepPecentageCriterion, lowpassCutoff=lowpassCutoff, gapThreshold=gapThreshold, minNumberOfPoints=minNumberOfPoints, returnShiftQ=returnShiftQ)
        if returnShiftQ:
            tHomogeneous, data, shiftQs = packedReturns
            self.dataCleaned[instrumentation].update({'shiftQs': shiftQs})
        else:
            tHomogeneous, data = packedReturns
        self.dataCleaned[instrumentation].update({tName: tHomogeneous})
        tHomogeneousLen = len(tHomogeneous)
        colIndFrom = 0
        for i, varName in enumerate(variableNames):
            shape = shapes[i]
            nColumn = np.prod(shape)
            colIndTo = colIndFrom + nColumn
            var_ = data[:, colIndFrom:colIndTo]
            var = var_.reshape((tHomogeneousLen, *shape))
            self.dataCleaned[instrumentation].update({varName: var})
        if inPlace:
            self.data[instrumentation] = self.dataCleaned[instrumentation]

##
def extractFiles(databaseDir, workDataDir, datasets, interval, keepIfExist=True):
    '''
    Purpose:
        This function is to extract files from a database directory to a work data directory the datasets during a period of time.
    parameters:
        databaseDir: implied by its name
        workDataDir: implied by its name
        datasets: a dictionary in terms of the path to find the compressed file. For example, {'Cluster': {'C1': 'C1_CP_FGM_FULL', 'C2': 'C2_CP_FGM_FULL'}}
        interval: a list of two elements, each a datetime object. This parameter determines during which interval the user requests the data. Its first element is the start time, and its second element end time.
    Caveat:
        This function only works for Cluster at the moment
    '''
    def getInterval(dataset, dateTime):
        if 'FGM_FULL' in dataset:
            start = datetime(dateTime.year, dateTime.month, 1)
            if dateTime.month == 12:
                end = datetime(dateTime.year+1, 1, 1)
            else:
                end = datetime(dateTime.year, dateTime.month+1, 1)
        elif 'FGM_5VPS' in dataset:
            if dateTime.month % 2 == 0:
                monthShift = -1
            elif dateTime.month % 2 == 1:
                monthShift = 0
            start = datetime(dateTime.year, dateTime.month+monthShift, 1)
            if dateTime.month > 10:
                end = datetime(dateTime.year+1, 1, 1)
            else:
                end = datetime(dateTime.year, dateTime.month + 2 + monthShift, 1)
        elif 'CIS-HIA' in dataset:
            start = datetime(dateTime.year, 1, 1)
            end = datetime(dateTime.year+1, 1, 1)
        return [start, end]

    datasetsKeys = list(datasets.keys())
    if len(datasetsKeys) == 1 and datasetsKeys[0] == 'Cluster':
        pass
    else:
        raise Exception('This function only works for Cluster at the moment')
    pathlist_ = ot.dict2list(datasets)
    start_ = interval[0]
    end_ = interval[1]
    with open('extraction.err', 'w') as file:
        for l_ in pathlist_:
            path = os.path.join(*l_)
            print(path)
            dataset = l_[2]
            start = start_
            while start < end_:
                start, end = getInterval(dataset, start)
                pathDatabase = os.path.join(databaseDir, path, 'compressed')
                if 'FGM_FULL' in dataset:
                    pathWorkData = os.path.join(workDataDir, path, start.strftime("%Y"), start.strftime("%m"))
                if 'FGM_5VPS' in dataset:
                    pathWorkData = os.path.join(workDataDir, path, start.strftime("%Y"))
                elif 'CIS-HIA' in dataset:
                    pathWorkData = os.path.join(workDataDir, path, start.strftime("%Y"))
                if not os.path.exists(pathWorkData):
                    os.makedirs(pathWorkData)
        #        else:
        #            files = os.listdir(pathWorkData)
        #            for file in files:
        #                filePath = os.path.join(pathWorkData, file)
        #                if os.path.isfile(filePath):
        #                    os.unlink(filePath)
                fileName = '--'.join([dataset, start.strftime("%Y--%m"),
                                     end.strftime("%Y--%m"), 'cdf.tar.gz'])
                srcName = os.path.join(pathDatabase, fileName)
        #        tar = tarfile.open()
                if keepIfExist:
                    cmdArgs = ['tar', '-xkf', srcName, '--strip-components=2', '--one-top-level={:s}'.format(pathWorkData)]
    #                cmd = ('tar -xkf {:s} --strip-components 2 ' + '--one-top-level={:s}').format(srcName, pathWorkData)
                else:
                    cmdArgs = ['tar', '-xvf', srcName, '--strip-components=2', '--one-top-level={:s}'.format(pathWorkData)]
    #                cmd = ('tar -xvf {:s} --strip-components 2 ' + '--one-top-level={:s}').format(srcName, pathWorkData)
    #            try:
    #                os.popen(cmd)
    #            except Exception:
    #                with open('extraction.error', 'a') as f:
    #                    f.write('{}\n'.format(fileName))
                process = subprocess.Popen(cmdArgs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                print(stdout.decode())
                print(stderr.decode())
                print(stdout.decode(), file=file)
                print(stderr.decode(), file=file)
    #            for line in stdout.decode():
    #                print(line)
    #            for line in stderr:
    #                print(line.decode())
                start = end
##


def readCDFInfo(cdfFile):
    '''
    Purpose:
        This function is to show what variables are contained in the cdfFile.
    Parameters:
        cdfFile: a cdf file object defined by cdflib
    '''
    info = cdfFile.cdf_info()
    variables = info['zVariables']
    if len(variables) == 0:
        variables = info['rVariables']
    varInfo = []
    varInfoDict = {}
    for i, var in enumerate(variables):
        varInfoDict[var] = {}
        varInfoDict[var].update(cdfFile.varinq(variable=var))
        varInfoDict[var].update(cdfFile.varattsget(variable=var))
        varInfo.append(varInfoDict[var])
    return varInfo, varInfoDict

##
def findFileNames(path, strings=None, size='allSize', timeTag=None):
    '''
    Purpose:
        This function is to find in the directory defined by <path> the names of files which according to 
    Parameters:
        path: a string defining a directory
        strings: a list of strings that should be contained in the objective file name.
        size: a string in the form of "<num" or ">num", where num is an integer whose unit is bytes. This parameter define either that the size of the file should be less than num bytes, or that its size greater than num bytes. The default value for this parameter is "allSize" which allow for all file size.
        timeTag: a datetime object that specify the time accurate to second. This parameter is designed for mms brst data which is segmented over a few minutes.
    '''
    if os.path.exists(path):
        with os.scandir(path) as it:
            fileNamesMeetingSize = []
            foundFileNames = []
            for entry in it:
                if entry.is_file():
                    if size[0] == '<' and entry.stat().st_size < int(size[1:]):
                        fileNamesMeetingSize.append(entry.name)
                    elif size[0] == '>' and entry.stat().st_size > int(size[1:]):
                        fileNamesMeetingSize.append(entry.name)
                    elif size == 'allSize':
                        fileNamesMeetingSize.append(entry.name)
            logging.debug("number of files meeting size criterion: {}".format(len(fileNamesMeetingSize)))
            if strings:
                if not isinstance(strings[0], list):
                    logging.debug("string criteria for searching:")
                    logging.debug(strings)
                    initDeltaTime = timedelta(minutes=20)
                    deltaTime = initDeltaTime
                    while fileNamesMeetingSize:
                        if all(string in fileNamesMeetingSize[0] for string in strings):
                            if timeTag:
                                items = fileNamesMeetingSize[0].split('_')
                                for item in items:
                                    if item.isnumeric():
                                        year = int(item[:4])
                                        month = int(item[4:6])
                                        day = int(item[6:8])
                                        hour = int(item[8:10])
                                        minute = int(item[10:12])
                                        second = int(item[12:14])
                                        timeOfFile = datetime(year, month, day, hour, minute, second)
                                        deltaTime_ = timeTag - timeOfFile
                                        if deltaTime_ > timedelta(seconds=1):
                                            if deltaTime_ < deltaTime:
                                                if deltaTime == initDeltaTime:
                                                    deltaTime = deltaTime_
                                                    foundFileNames.append(fileNamesMeetingSize.pop(0))
                                                else:
                                                    deltaTime = deltaTime_
                                                    del foundFileNames[-1]
                                                    foundFileNames.append(fileNamesMeetingSize.pop(0))
                                            else: del fileNamesMeetingSize[0]
                                        else: del fileNamesMeetingSize[0]
                            else:
                                foundFileNames.append(fileNamesMeetingSize.pop(0))
                        else: del fileNamesMeetingSize[0]
                else:
                    while fileNamesMeetingSize:
                        for stringList in strings:
                            if all(string in fileNamesMeetingSize[0] for string in strings):
                                foundFileNames.append(fileNamesMeetingSize.pop(0))
                                break
                        else: del fileNamesMeetingSize[0]
    else:
        logging.debug("path does not exist: {}".format(path))
        foundFileNames = []
    return foundFileNames
##

def findFileNamesInInterval(path, interval=None, timeFMT='%Y%m%d', strings=None, size='allSize', findAll=False):
    '''
    Purpose:
        This function is dedicated to find file names in the form of *20020522*2002023*, such as Cluster files.
    Parameters:
        path: the path to the folder that contains the file
        interval: a list of two datetime object, start and end, that specifies the temporal interval of the data.
    '''
    foundFileNames = []
    criteria_ = []
    if strings:
        criteria_.extend(strings)
    if interval:
        start = interval[0]
        end = interval[1]
        unFoundTimeTags = []
        deltaTime = end - start
        numberOfDays = deltaTime.days
        if deltaTime.total_seconds() % (3600*24):
            numberOfDays += 1
        day = timedelta(days=1)
        timeTags = [[(start+i*day).strftime("%Y%m%d"),
                    (start+(i+1)*day).strftime("%Y%m%d")]
                    for i in range(numberOfDays)]
    while timeTags:
        criteria = criteria_.copy()
        criteria.extend(timeTags[0])
        foundFileName = findFileNames(path, strings=criteria, size=size)
        if foundFileName:
            foundFileNames.extend(foundFileName)
            del timeTags[0]
        else:
            unFoundTimeTags.append(timeTags.pop(0))
    if findAll and unFoundTimeTags:
        raise NameError('files not found: {}'.format(unFoundTimeTags))
    elif unFoundTimeTags:
        print('warning! files not found: {}'.format(unFoundTimeTags))
    return foundFileNames

##
def defineFiles(workDataDir, datasets, epoch, size='allSize', silence=True):
    '''
    Purpose:
        This function is to find the files stored in workDataDir. It returns a list of DataFile objects whose attribute DataFile.cdfFile is a cdfFile Object.
    Parameter:
        workDataDir: the path to the directory of the working database
        datasets: a dictionay in terms of the path to dataset. For example, {'Cluster': {'C1' : 'C1_CP_FGM_FULL', 'C2': 'C2_CP_FGM_FULL'}, 'mms': {'mms1': {'fgm': {'brst': 'l2'}}}}
        epoch: the epoch at which the destined files contain data
        size: a string in the form of "<num" or ">num", where num is an integer whose unit is bytes. This parameter define either that the size of the file should be less than num bytes, or that its size greater than num bytes. The default value for this parameter is "allSize" which allow for all file size.
        silence: to print nothing
    '''
    pathSplitedList = ot.dict2list(datasets) # in the form of [['Cluster', 'C1', 'C1_CP_FGM_FULL'], ['Cluster', ...]]
    dataFiles = []
    for splitedPathToDataset in pathSplitedList:
        dataFile = defineFile(workDataDir, splitedPathToDataset, epoch, size=size, silence=silence)
        dataFiles.append(dataFile)
    return dataFiles
##

def defineFile(workDataDir, splitedPathToDataset, epoch, size='AllSize', silence=True):
    '''
    Purpose:
        This function is to find the file stored in workDataDir. It returns a DataFile objects whose attribute DataFile.cdfFile is a cdfFile Object.
    Parameter:
        workDataDir: the path to the directory of the working database
        splitedPathToDataset: a list in terms of the path to the dataset. For example, ['mms', 'mms1', 'fgm', 'brst', 'l2']
        epoch: the epoch at which the destined file contain data
        size: a string in the form of "<num" or ">num", where num is an integer whose unit is bytes. This parameter define either that the size of the file should be less than num bytes, or that its size greater than num bytes. The default value for this parameter is "allSize" which allow for all file size.
        silence: to print nothing
    '''
    start = datetime(*cdflib.cdfepoch.breakdown(epoch)[:6])
    end = start + timedelta(seconds=2)
    interval = [start, end]
    dataFile = DataFile(splitedPathToDataset=splitedPathToDataset, dateTime=start)
    absolutePathToFile = os.path.join(workDataDir, dataFile.pathToFile)
    criteria = dataFile.instrumentation.copy()
    if dataFile.mission == 'mms':
        if 'fgm' in dataFile.instrumentationPath:
            if 'brst' in dataFile.instrumentationPath:
                criteria.append(start.strftime("%Y%m%d"))
                timeTag = start
            elif 'srvy' in dataFile.instrumentationPath:
                criteria.append(start.strftime("%Y%m%d"))
                timeTag = None
            fileNames = findFileNames(absolutePathToFile, strings=criteria, size=size, timeTag=timeTag)
            if len(fileNames) == 1:
                dataFile.fileName = fileNames[0]
            else:
                raise Exception('The number of found file names is not 1')
            if not silence:
                print(dataFile.fileName)
            dataFile.cdfFile = cdflib.CDF(os.path.join(absolutePathToFile, dataFile.fileName))
        elif 'fpi' in dataFile.instrumentationPath:
            if 'fast' in dataFile.instrumentationPath: # this branch need modefication. Its current version does not work
#                    raise Exception('The program fails to define mms/fpi/fast dataFile, because the program is not complete and need modification')
                criteria.append(start.strftime("%Y%m%d"))

                fileNames = findFileNames(absolutePathToFile, strings=criteria, size=size, timeTag=None)
                if len(fileNames) == 1:
                    dataFile.fileName = fileNames[0]
                else:
                    raise Exception('The number of found file names is not 1')
                dataFile.cdfFile = cdflib.CDF(os.path.join(absolutePathToFile, dataFile.fileName))
            elif 'brst' in dataFile.instrumentationPath:
                criteria.append(start.strftime("%Y%m%d"))
                timeTag = start
                fileNames = findFileNames(absolutePathToFile, strings=criteria, size=size, timeTag=timeTag)
                if len(fileNames) == 1:
                    dataFile.fileName = fileNames[0]
                else:
                    raise Exception('The number of found file names is not 1')
                dataFile.cdfFile = cdflib.CDF(os.path.join(absolutePathToFile, dataFile.fileName))
    elif dataFile.mission == 'Cluster':
        if 'FGM' in dataFile.instrumentationPath:
            criteria.append(start.strftime("%Y%m%d"))
            fileNames = findFileNamesInInterval(absolutePathToFile, interval=interval, strings=dataFile.instrumentation)
            if len(fileNames) == 1:
                dataFile.fileName = fileNames[0]
            else:
                raise Exception('The number of found file names is not 1')
            dataFile.cdfFile = cdflib.CDF(os.path.join(absolutePathToFile, dataFile.fileName))
    if not silence:
        print(dataFile.fileName)
    return dataFile


def readData(workDataDir, datasetsAndVariables, datetimeRange, epochType='CDF_EPOCH', workDataDirsBak=None, copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak=True):
    '''
    Purpose:
        This function is to load data.
    Parameters:
        workDataDir: the path to the directory of working data base
        datasetAndVariables: a dictionary in terms of the path to the datasets. Its leaf value is a list of variable names. The names are defined by the corresponding cdfFile. For example, {'Cluster': {'C1' : {'C1_CP_FGM_FULL': ['time_tags__C1_CP_FGM_FULL', 'B_vec_xyz_gse__C1_CP_FGM_FULL']}}, 'mms': {'mms1': {'fgm': {'brst': {'l2': ['Epoch', 'Btotal']}}}}}. Please note that Epoch and Btotal are improvised. It may also be a list of lists, with each sublist in the form of ['Cluster', 'C1', 'C1_CP_FGM_FULL', ['time_tags__C1_CP_FGM_FULL', 'B_vec_xyz_gse__C1_CP_FGM_FULL']].
        datetimeRange: a list of two elements which define the time interval during which the data is to be retrieved. [start, end]
        workDataDirsBak: a list of backup databases
        copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak: self-explanatory
    Return:
        variablesAllDataset: a list of dict. The list is for multiple datasets, the keys of a dict is for the variables in a dataset.
    '''
    start = datetimeRange[0]
    end = datetimeRange[1]
    if isinstance(datasetsAndVariables, dict):
        datasetsAndVariablesList = ot.dict2list(datasetsAndVariables)
    elif isinstance(datasetsAndVariables, list):
        datasetsAndVariablesList = datasetsAndVariables
    variablesAllDataset = []
    for datasetAndVariables in datasetsAndVariablesList:
        splitedPathToDataset = datasetAndVariables[:-1]
        variableNames = datasetAndVariables[-1]
        start_ = start
        variablesInADatasetOverDatetimeRange = [] # first index for data from different files over a epoch range, second index for variables
        variablesInADatasetIndependantOnTime = []
        while start_ < end:
#            start_NextDayDateTime = datetime(*cdflib.cdfepoch.breakdown(start_)[:6]) + timedelta(days=1)
            start_NextDayDateTime = start_ + timedelta(days=1)
            endOfTheDay = datetime(start_NextDayDateTime.year, start_NextDayDateTime.month, start_NextDayDateTime.day)
#            epochEndOfTheDay = cdflib.cdfepoch.compute_epoch(ot.datetime2list(endOfTheDay))
            if endOfTheDay >= end:
                datetimeRange = [start_, end]
            else:
                datetimeRange = [start_, endOfTheDay]
            dataFile = DataFile(workDataDir=workDataDir, splitedPathToDataset=splitedPathToDataset, dateTime=start_)
            if dataFile.filePath:
                pass
            elif workDataDirsBak:
                for workDataDirBak in workDataDirsBak:
                    dataFile = DataFile(workDataDir=workDataDirBak, splitedPathToDataset=splitedPathToDataset, dateTime=start_, defineCDFFileQ=False)
                    if dataFile.filePath:
                        if copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak:
                            relpath = os.path.relpath(dataFile.filePath, workDataDirBak)
                            destFilePath = os.path.join(workDataDir, relpath)
                            logging.info('copying {}'.format(destFilePath))
                            print('file not found in wordDataDir, but found in workDataDirBak: {}'.format(destFilePath))
                            print('now copying...')
                            destDirPath = os.path.dirname(destFilePath)
                            if not os.path.exists(destDirPath):
                                os.makedirs(destDirPath)
                            cmdArgs = ['cp', dataFile.filePath, destFilePath]
                            process = subprocess.Popen(cmdArgs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            process.wait()
                            logging.info('{} copied'.format(destFilePath))
                            print('file copied')
                            dataFile = DataFile(filePath=destFilePath)
                        break
            if not dataFile.filePath:
                raise Exception('data file not found')
            assert dataFile.cdfFile is not None
#            logging.info(splitedPathToDataset)
#            logging.info(dataFile.cdfFile)
            logging.info('reading data file: {}'.format(dataFile.filePath))
            logging.info('datetimeRange: {}'.format(datetimeRange))
            dataMajor, dataAux = readDataFromACdfFile(dataFile.cdfFile, variableNames, datetimeRange, epochType=epochType)
            logging.info('reading data file done: {}'.format(dataFile.filePath))
            variablesInADatasetOverDatetimeRange.append(dataMajor)
            variablesInADatasetIndependantOnTime.append(dataAux)
            start_ = endOfTheDay
        variables = {}
        for var in variablesInADatasetOverDatetimeRange[0].keys():
            variables[var] = np.concatenate([varsOfATimeRange[var] for varsOfATimeRange in variablesInADatasetOverDatetimeRange if varsOfATimeRange[var].size>0], axis=0)
        for fileInd in range(len(variablesInADatasetIndependantOnTime)-1):
            logging.debug(variablesInADatasetIndependantOnTime[fileInd].keys())
            for key in variablesInADatasetIndependantOnTime[fileInd].keys():
                assert np.all(variablesInADatasetIndependantOnTime[fileInd][key] == variablesInADatasetIndependantOnTime[fileInd+1][key])
        variables.update(variablesInADatasetIndependantOnTime[0])
#        for indOfVariable in range(len(variablesInADatasetOverDatetimeRange[0])):
#            variables.append(np.concatenate([variables[indOfVariable] for variables in variablesInADatasetOverDatetimeRange], axis=0))
        variablesAllDataset.append(variables)
    return variablesAllDataset

##
def readDataFromACdfFile(cdfFile, variables=None, datetimeRange=None, epochType='CDF_EPOCH'):
    varInfo, varInfoDict = readCDFInfo(cdfFile)
    epochDataInd = 0
    dataMajor = {}
    dataAux = {} # data not dependant on epoch
    if epochType == 'CDF_EPOCH':
        epochRange = [cdflib.cdfepoch.compute_epoch(ot.datetime2list(dateTime)) for dateTime in datetimeRange]
    elif epochType == 'CDF_EPOCH16':
        epochRange = [cdflib.cdfepoch.compute_epoch16(ot.datetime2list(dateTime)) for dateTime in datetimeRange]
    elif epochType == 'CDF_TIME_TT2000':
        epochRange = [cdflib.cdfepoch.compute_tt2000(ot.datetime2list(dateTime)) for dateTime in datetimeRange]
    recordRange = cdfFile.epochrange(epochDataInd, *epochRange)
    logging.info('record range:')
    logging.info(recordRange)
    for var in variables:
        majorData = True
        depend0 = varInfoDict[var].get('DEPEND_0', None)
        if depend0 is None:
            if varInfo[epochDataInd]['Variable'] == var:
                pass
            else:
                majorData = False
        else:
            assert depend0 == varInfo[epochDataInd]['Variable']
        if majorData:
            if recordRange is not None:
                dataMajor[var] = cdfFile.varget(var, startrec=recordRange[0], endrec=recordRange[1])
            else:
                dataMajor[var] = np.array([])
        else:
            dataAux[var] = cdfFile.varget(var)
    return dataMajor, dataAux


def readPDSData(fileName, dataFileExtension='.TAB', infoFileExtension='.xml', sep=None):
    infoFile = fileName + infoFileExtension
    xmlTree = ET.parse(infoFile)
    root = xmlTree.getroot()
    for l1 in root:
        if 'File_Area_Observational' in l1.tag:
            for l2 in l1:
                if 'Table_Character' in l2.tag:
                    for l3 in l2:
                        if 'Record_Character' in l3.tag:
                            record_character = l3
    dataDict = {}
    columnNames = []
    for child in record_character:
        if 'Field_Character' in child.tag:
            for field in child:
                tagName = field.tag.split('}')[-1]
                if tagName == 'name':
                    columnName = field.text
                    columnNames.append(columnName)
                    dataDict[columnName] = {}
                    break
            else:
                raise Exception("Didn't find name!")
            for field in child:
                tagName = field.tag.split('}')[-1]
                dataDict[columnName][tagName] = field.text
    if dataFileExtension.lower() == '.tab':
        dataFile = fileName + dataFileExtension
    with open(dataFile, 'r') as f:
        info = f.readlines()

    for lineInd, line in enumerate(info):
        if sep is None:
            lineInfo = info[lineInd].strip().split()
        else:
            lineInfo = info[lineInd].strip().split(sep=sep)
        for columnInd, columnInfo in enumerate(lineInfo):
            columnName = columnNames[columnInd]
            dataType = dataDict[columnName]['data_type']
            if dataType == 'ASCII_Date_Time_YMD_UTC':
                dateStr = columnInfo[:-1]
                yearmmdd, hhmmss = dateStr.split('T')
                year, month, day = yearmmdd.split('-')
                hour, minute, secondInfo = hhmmss.split(':')
                second, millisecond = secondInfo.split('.')
                data_ = cdflib.cdfepoch.compute_epoch(ot.datetime2list(datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), int(millisecond)*1000)))
            elif dataType == 'ASCII_String':
                data_ = columnInfo
            elif dataType in ['ASCII_Integer', 'ASCII_Real']:
                data_ = float(columnInfo)
            if lineInd == 0:
                dataDict[columnName]['data'] = []
            dataDict[columnName]['data'].append(data_)
    for key in dataDict.keys():
        dataDict[key]['data'] = np.array(dataDict[key]['data'])
    return dataDict, columnNames
