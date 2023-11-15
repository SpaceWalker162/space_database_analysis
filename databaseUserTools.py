__author__ = 'Yufei Zhou'

import numpy as np
import pandas as pd
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
import struct
from pprint import pprint
from scipy.signal import butter, lfilter, freqz
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
import scipy.interpolate as interpolate
import scipy.special
import xml.etree.ElementTree as ET
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


multispacecraftMissions = ['cluster', 'mms', 'themis']
piledInstrumentation = ['5VPS', 'CIS']
instrumentationFileUnderYear = [('cluster', '5VPS'), ('cluster', 'CIS'), ('ace', 'swepam'), ('ace', 'mag')]
instrumentationFileUnderMonth = None 
instrumentationDiviedIntoMonth = ['FULL', 'mms']
missionInfo = {
        'mms': {
            'epochType': 'CDF_TIME_TT2000',
            'multispacecraftMission': True,
            'dataset': [
                {
                'dataset_flags': ['fpi', 'fast'],
                'dataset_file_time_gap': timedelta(seconds=7200)
                },
                {
                'dataset_flags': ['edp', 'fast'],
                'dataset_file_time_gap': timedelta(days=1)
                },
                {
                'dataset_flags': ['edp', 'slow'],
                'dataset_file_time_gap': timedelta(days=1)
                }
                        ]
            },
        'cluster': {
            'epochType': 'CDF_EPOCH',
            'multispacecraftMission': True},
        'themis': {
            'epochType': 'CDF_EPOCH',
            'multispacecraftMission': True},
        'ace': {
            'epochType': 'CDF_EPOCH',
            'multispacecraftMission': False},
        'wind': {
            'epochType': 'CDF_EPOCH',
            'multispacecraftMission': False},
        'geotail': {
            'epochType': 'CDF_EPOCH',
            'multispacecraftMission': False},
        'cassini': {
            'epochType': 'CDF_EPOCH',
            'multispacecraftMission': False},
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
    def __init__(self, mission=None, spacecraft=None, instrumentation=None, instrumentationPath=None, pathToDataset=None, splitedPathToDataset=None, dateTime=None, timeRange=None, filePeriodRange=None, dataset_file_time_gap=None, fileName=None, instrumentationObj=None, cdfFile=None, workDataDir=None, fileSize='allSize', epochType='CDF_EPOCH', silence=True, filePath=None, filePaths=[], defineCDFFileQ=True):
        '''
        Parameter:
            workDataDir: the path to the directory of the working database
            splitedPathToDataset: a list in terms of the path to the dataset. For example, ['mms', 'mms1', 'fgm', 'brst', 'l2']
            dateTime: the epoch at which the destined file contain data
            timeRange: the time range in which all appropriate files shall be found. This is desigend for mms brst data.
            filePeriodRange: the start time to the end time of the file to be found.
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
        self.timeRange = timeRange
        self.filePeriodRange = filePeriodRange
        self.dataset_file_time_gap = dataset_file_time_gap
        if self.filePeriodRange is not None:
            self.dataset_file_time_gap = self.filePeriodRange[1] - self.filePeriodRange[0]
#            self.epoch = cdflib.cdfepoch.compute_epoch(ot.datetime2list(dateTime))
        self.fileName = fileName
        self.filePath = filePath
        self.filePaths = filePaths
        self.cdfFile = cdfFile
#        for crit in instrumentationFileUnderYear:
#            if all([name in self.pathToDataset for name in crit]):
##                self.pathToFile = self.pathToDataset
#                self.pathToFile = os.path.join(self.pathToDataset, self.dateTime.strftime('%Y'))
#                break
#        else:
#            self.pathToFile = os.path.join(self.pathToDataset, self.dateTime.strftime('%Y'), self.dateTime.strftime('%m'))
        if not self.filePath and not self.filePaths:
            self.findFilePath(size=fileSize, silence=silence)
        if self.filePath:
            if defineCDFFileQ:
                self.defineCDFFile()
        if self.filePaths:
            if defineCDFFileQ:
                self.defineCDFFiles()


    def defineCDFFile(self):
        self.cdfFile = cdflib.CDF(self.filePath)

    def defineCDFFiles(self):
        self.cdfFiles = []
        for path in self.filePaths:
            self.cdfFiles.append(cdflib.CDF(path))

    def findFilePath(self, size='allSize', silence=True):
        '''
        Purpose:
            This function is to find the file stored in self.workDataDir.
        '''
        if self.dateTime:
            start = self.dateTime
            end = start + timedelta(seconds=1)
        absolutePathToDataset = os.path.join(self.workDataDir, self.pathToDataset)
        criteria = []
        timeTag = None
        logging.debug("mission: {}".format(self.mission))
        if self.mission == 'mms':
            if 'brst' in self.instrumentationPath:
                searchMethod = 'allFilesInTimeRange'
            else:
                searchMethod = 'general'
                criteria.extend(self.instrumentation.copy())
                if self.dataset_file_time_gap is not None:
                    if self.dataset_file_time_gap < timedelta(days=1):
                        criteria.append(self.filePeriodRange[0].strftime("%Y%m%d%H"))
                    elif self.dataset_file_time_gap == timedelta(days=1):
                        criteria.append(start.strftime("%Y%m%d"))
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
        elif self.mission == 'cassini':
            searchMethod = 'general'
            if 'CO-E_SW_J_S-MAG-4-SUMM-1MINAVG-V2.0' in self.instrumentation:
                criteria.append(start.strftime("%Y"))
            else:
                criteria.append(start.strftime("%Y%m%d"))
        else:
            searchMethod = 'general'
            raise Exception('mission not defined!')
        logging.debug("Looking for data files in: {}".format(absolutePathToDataset))
        logging.debug("search method: {}".format(searchMethod))
        if searchMethod == 'general':
            fileNames, absolutePathToFile = self.findFileNamesUnderDataset(absolutePathToDataset, func=findFileNames, strings=criteria, size=size, timeTag=timeTag)
#            = findFileNames(absolutePathToFile)
        elif searchMethod == 'allFilesInTimeRange':
            getTimeFromNameFunc_ = getTimeFromName
            logging.debug('searchmethod:')
            logging.debug(self.timeRange)
            fileNames, absolutePathToFile = self.findFileNamesUnderDataset(absolutePathToDataset, func=findFileNamesInTimeRange, timeRange=self.timeRange, getTimeFromNameFunc=getTimeFromNameFunc_, strings=criteria, size=size)
        elif searchMethod == 'ClusterCAA':
            fileNames, absolutePathToFile = self.findFileNamesUnderDataset(absolutePathToDataset, func=findFileNamesInInterval, interval=interval, strings=self.instrumentation)
#            fileNames = findFileNamesInInterval(absolutePathToFile, interval=interval, strings=self.instrumentation)
        else:
            raise Exception('unknown method for searching files')
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
            if searchMethod == 'allFilesInTimeRange':
                paths = []
                for fileName in fileNames:
                    paths.append(os.path.join(absolutePathToFile, fileName))
                self.filePaths = paths
            else:
                logging.warning("More than one files were found:" + ("\n{}"*numFiles).format(*fileNames))

#            raise Exception('More than one files were found.')


    def findFileNamesUnderDataset(self, absolutePathToDataset, func, **para):
        fileNames = func(absolutePathToDataset, **para)
        absolutePathToFile = absolutePathToDataset
        if len(fileNames) == 0:
            if self.dateTime:
                absolutePathToDatasetYear = os.path.join(absolutePathToDataset, self.dateTime.strftime('%Y'))
            elif self.timeRange:
                if self.timeRange[0].year == self.timeRange[1].year:
                    absolutePathToDatasetYear = os.path.join(absolutePathToDataset, self.timeRange[0].strftime('%Y'))
                else:
                    raise Exception('time range across years')
            absolutePathToFile = absolutePathToDatasetYear
            logging.debug("Not found. Looking for data files in: {}".format(absolutePathToDatasetYear))
            fileNames = func(absolutePathToDatasetYear, **para)
            if len(fileNames) == 0:
                if self.dateTime:
                    absolutePathToDatasetYearMonth = os.path.join(absolutePathToDatasetYear, self.dateTime.strftime('%m'))
                elif self.timeRange:
                    if self.timeRange[0].month == self.timeRange[1].month:
                        absolutePathToDatasetYearMonth = os.path.join(absolutePathToDatasetYear, self.timeRange[0].strftime('%m'))
                    else:
                        raise Exception('time range across months')
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

    def loadData(self, datetimeRange=None, instrumentation=None, instrumentationRetrivingName=None, variables=None, variableRetrivingNames=None, datasetAndVariables=None, variablesWithRetrivingNames=None, instrumentationVariablesWithRetrivingName=None, cleanData=False, tStepPecentageCriterion=0.9, lowpassCutoff=None, inPlace=False, gapThreshold=None, minNumberOfPoints=None, returnShiftQ=False, copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak=True, fromFile=None):
        '''
        Parameters:
            datetimeRange: a list of two datetime objects, representing the start and end of the data to be loaded.
            instrumentation: a list in terms of the directory names at all levels below spacecraft and above year or files. For example, ['fgm', 'brst', 'l2']
            instrumentationRetrivingName: a string.
            datasetAndVariables: this parameter should not be used by user. 
            instrumentationVariablesWithRetrivingName: a dictionary in terms of the path to the datasets. Its leaf value is a list of variable names. The names are defined by the corresponding cdfFile. For example, {'Cluster': {'C1' : {'C1_CP_FGM_FULL': ['FGM', ('time_tags__C1_CP_FGM_FULL', 't'), ('B_vec_xyz_gse__C1_CP_FGM_FULL', 'B')]}}, 'mms': {'mms1': {'fgm': {'brst': {'l2': ['fgmBrst', ('Epoch', 't'), ('Btotal', 'BTotal')]}}}}}. Please note that 'Epoch' and 'Btotal' are improvised. This parameter may also be a list of lists, with each sublist in the form of ['Cluster', 'C1', 'C1_CP_FGM_FULL', ['FGM', ('time_tags__C1_CP_FGM_FULL', 't'), ('B_vec_xyz_gse__C1_CP_FGM_FULL', 'B')]]. To retrieve data, for example, of 'B_vec_xyz_gse__C1_CP_FGM_FULL', use C1.data['FGM']['B']
            fromFile: to load data from a file, the file name is given by this parameter.
        Note:
            Necessary parameters include: <fromFile, [datetimeRange, <[<[instrumentation, variables], datasetsAndVariables>, instrumentationRetrivingName, variableRetrivingNames], instrumentationVariablesWithRetrivingName>]>
            To retrieve data, use Spacecraft.data[instrumentationName][variableRetrivingName]
        '''
        if fromFile:
            cdfFile = cdflib.CDF(fromFile)
            if variablesWithRetrivingNames is not None:
                variableNames = variablesWithRetrivingNames
            else:
                variableNames = None
            dataFromACdfFile = readDataFromACdfFile(cdfFile, variableNames, epochType=self.epochType)
            self.data.update({instrumentationRetrivingName: dataFromACdfFile})
        else:
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
                if variablesWithRetrivingNames:
                    variables = [tu[0] for tu in variablesWithRetrivingNames]
                else:
                    variableRetrivingNames = variableRetrivingNames.copy()
                    if datasetAndVariables is not None:
                        variables = datasetAndVariables[-1]
                    variablesWithRetrivingNames = list(zip(variables, variableRetrivingNames))
                if datasetAndVariables is None:
                    if isinstance(instrumentation, str):
                        instrumentation = [instrumentation]
                    if self.mission in multispacecraftMissions:
                        datasetAndVariables = [self.mission, self.name, *instrumentation, variables]
                    else:
                        datasetAndVariables = [self.name, *instrumentation, variables]
    #            else:
    #                variablesWithRetrivingNames = list(zip(datasetAndVariables[-1], variableRetrivingNames))
                datasetsAndVariables = [datasetAndVariables]
                data_ = self.readData(self.workDataDir, datasetsAndVariables, datetimeRange=datetimeRange, epochType=self.epochType, workDataDirsBak=self.workDataDirsBak, copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak=copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak)
                logging.info(data_[0].keys())
                dataInDict = dict([(varRetName, data_[0][varName]) for varName, varRetName in variablesWithRetrivingNames])
                self.data.update({instrumentationRetrivingName: dataInDict})
                if cleanData:
                    self.cleanData(instrumentation=instrumentationRetrivingName, tStepPecentageCriterion=tStepPecentageCriterion, lowpassCutoff=lowpassCutoff, inPlace=inPlace, gapThreshold=gapThreshold, minNumberOfPoints=minNumberOfPoints, returnShiftQ=returnShiftQ)

    def readData(self, workDataDir, datasetsAndVariables, datetimeRange, epochType='CDF_EPOCH', workDataDirsBak=None, copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak=True):
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
            logging.info('finding files for dataset and variables: ')
            logging.info(datasetAndVariables)
            splitedPathToDataset = datasetAndVariables[:-1]
            variableNames = datasetAndVariables[-1]
            start_ = start
            variablesInADatasetOverDatetimeRange = [] # first index for data from different files over a epoch range, second index for variables
            variablesInADatasetIndependantOnTime = []
            if 'brst' in splitedPathToDataset:
                allowedTimeExtension = timedelta(seconds=600)
                dataFileTimeRange = [start-allowedTimeExtension, end]
                logging.debug('brst time range')
                logging.debug(dataFileTimeRange)
                dataFile = DataFile(workDataDir=workDataDir, splitedPathToDataset=splitedPathToDataset, timeRange=dataFileTimeRange)
                workDataFile = dataFile
                if workDataDirsBak:
                    for workDataDirBak in workDataDirsBak:
                        dataFileBak = DataFile(workDataDir=workDataDirBak, splitedPathToDataset=splitedPathToDataset, timeRange=dataFileTimeRange, defineCDFFileQ=False)
                        if dataFileBak.filePaths:
                            desiredDataPaths = []
                            for filePathBak in dataFileBak.filePaths:
                                workFilePath = os.path.join(os.path.relpath(filePathBak, workDataDirBak), workDataDir)
                                if workFilePath in workDataFile.filePaths:
                                    desiredDataPaths.append(workFilePath)
                                else:
                                    if copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak:
                                        destFilePath = workDataDirCopy(filePathBak, workDataDir, workDataDirBak)
                                        desiredDataPaths.append(destFilePath)
                                    else:
                                        desiredDataPaths.append(filePathBak)
                desiredDataPaths.sort()
                dataFile = DataFile(filePaths=desiredDataPaths)
                for fileInd, cdfFile in enumerate(dataFile.cdfFiles):
                    logging.info('reading data file: {}'.format(dataFile.filePaths[fileInd]))
                    logging.info('datetimeRange: {}'.format(datetimeRange))
                    dataMajor, dataAux = readDataFromACdfFile(cdfFile, variableNames, datetimeRange, epochType=epochType)
                    logging.info('reading data file done: {}'.format(dataFile.filePath))
                    variablesInADatasetOverDatetimeRange.append(dataMajor)
                    variablesInADatasetIndependantOnTime.append(dataAux)
            else:
                logging.debug('not in brst')
                while start_ < end:
                    for possible_dataset in missionInfo[self.mission]['dataset']:
                        if all([flag in splitedPathToDataset for flag in possible_dataset['dataset_flags']]):
                            dataset_file_time_gap = possible_dataset.get('dataset_file_time_gap')
                            break
                    else:
                        dataset_file_time_gap = timedelta(days=1)
                    if dataset_file_time_gap == timedelta(days=1):
                        beginOfTheFilePeriod = datetime(start_.year, start_.month, start_.day)
                        endOfTheFilePeriod = datetime(start_.year, start_.month, start_.day+1)
                    elif dataset_file_time_gap < timedelta(days=1):
                        beginOfTheDay = datetime(start_.year, start_.month, start_.day)
                        completePeriodBefore = (start_ - beginOfTheDay)//dataset_file_time_gap
                        beginOfTheFilePeriod = beginOfTheDay + completePeriodBefore * dataset_file_time_gap
                        endOfTheFilePeriod = beginOfTheDay + (completePeriodBefore+1)*dataset_file_time_gap 
                    if endOfTheFilePeriod >= end:
                        datetimeRange = [start_, end]
                    else:
                        datetimeRange = [start_, endOfTheFilePeriod]
                    dataFile = DataFile(workDataDir=workDataDir, splitedPathToDataset=splitedPathToDataset, dateTime=start_, filePeriodRange=[beginOfTheFilePeriod, endOfTheFilePeriod])
                    if dataFile.filePath:
                        pass
                    elif workDataDirsBak:
                        for workDataDirBak in workDataDirsBak:
                            dataFile = DataFile(workDataDir=workDataDirBak, splitedPathToDataset=splitedPathToDataset, dateTime=start_, filePeriodRange=[beginOfTheFilePeriod, endOfTheFilePeriod], defineCDFFileQ=False)
                            if dataFile.filePath:
                                if copyDataFileToWorkDataDirIfDataOnlyInWorkDataDirsBak:
                                    filePath = dataFile.filePath

                                    destFilePath = workDataDirCopy(filePath, workDataDir, workDataDirBak)
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
                    start_ = endOfTheFilePeriod
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

    def renameData(self, instrumentation, variablesAndNewNames=None):
        '''
        Purpose:
            rename data.
        Parameters:
            instrumentation: a string for the name of the instrumentation.
            variablesAndNewNames: example: [('Epoch', 't'), ('v_gse', 'v')]
        '''
        dic = self.data[instrumentation]
        for old_key, new_key in variablesAndNewNames:
            dic[new_key] = dic.pop(old_key)


class CelestialReferenceFrames:
    def __init__(self):
        self.vectors_in_ICRF = {}

    def ksmBasisInICRFBasis(self, t, sunPositionFilePath=None, saturnPositionFilePath=None):
        if 'sun_position' not in self.vectors_in_ICRF:
            dataDict, _ = readHoriaonsData(sunPositionFilePath)
            tSun = dataDict['Calendar Date (TDB)']['data']
            posCartesianSun = np.concatenate([dataDict['X']['data'][:, None], dataDict['Y']['data'][:, None], dataDict['Z']['data'][:, None]], axis=-1)
            self.vectors_in_ICRF.update({'sun_position': {'t': tSun, 'vector': posCartesianSun}})
        if 'saturn_position' not in self.vectors_in_ICRF:
            dataDict, _ = readHoriaonsData(saturnPositionFilePath)
            tSaturn = dataDict['Calendar Date (TDB)']['data']
            posCartesianSaturn = np.concatenate([dataDict['X']['data'][:, None], dataDict['Y']['data'][:, None], dataDict['Z']['data'][:, None]], axis=-1)
            self.vectors_in_ICRF.update({'saturn_position': {'t': tSaturn, 'vector': posCartesianSaturn}})
        t_range_available_sun = self.vectors_in_ICRF['sun_position']['t'][[0, -1]]
        t_range_available_saturn = self.vectors_in_ICRF['saturn_position']['t'][[0, -1]]
        t_ranges = np.concatenate([t_range_available_sun[None, :], t_range_available_saturn[None, :]], axis=0)
        t_range_available = [np.max(t_ranges[:, 0]), np.min(t_ranges[:, 1])]
        assert t_range_available[0] <= t[0] and t[1] <= t_range_available[1]
        posDataDict = {
                'tSaturn': self.vectors_in_ICRF['saturn_position']['t'],
                'posCartesianSaturn': self.vectors_in_ICRF['saturn_position']['vector'],
                'tSun': self.vectors_in_ICRF['sun_position']['t'],
                'posCartesianSun': self.vectors_in_ICRF['sun_position']['vector'],
                       }
        ksmBasisInICRFBasisData = dat.ksmBasisInICRFBasis(t, **posDataDict)
        return ksmBasisInICRFBasisData

##
def workDataDirCopy(filePath, workDataDir, workDataDirBak):
    relpath = os.path.relpath(filePath, workDataDirBak)
    destFilePath = os.path.join(workDataDir, relpath)
    logging.warning('file not found in wordDataDir, but found in workDataDirBak: {}'.format(destFilePath))
    logging.warning('now copying...')
    os.makedirs(os.path.dirname(destFilePath), exist_ok=True)
    cmdArgs = ['cp', filePath, destFilePath]
    process = subprocess.Popen(cmdArgs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    logging.info('{} copied'.format(destFilePath))
    logging.warning('file copied')
    return destFilePath


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
    variables = info.zVariables
    if len(variables) == 0:
        variables = info.rVariables
    varInfo = []
    varInfoDict = {}
    for i, var in enumerate(variables):
        varInfoDict[var] = {}
        varInfoDict[var]['varInfo'] = cdfFile.varinq(variable=var)
        varInfoDict[var]['varAtts'] = cdfFile.varattsget(variable=var)
        varInfo.append(varInfoDict[var])
    return varInfo, varInfoDict


def getTimeFromName(name):
    '''
    Purpose: get a datetime.datetime object from the name of a data file such as mms brst data file: mms1_fpi_brst_l2_dis-moms_20151007163101_v3.3.0.cdf
    Parameters:
        name: the name of the file
    '''
    timeString = name.split('_')[-2]
    fmt = "%Y%m%d%H%M%S"
    return datetime.strptime(timeString, fmt)

def findFileNamesInTimeRange(path, timeRange=None, getTimeFromNameFunc=None, strings=None, size='allSize', ext='.cdf'):
    '''
    Purpose:
        This function is to find in the directory defined by <path> the names of those files whose data are in the timeRange. This function is designed for mms brst data file whose name is irregular in epoch division and only contains the start of the epoch.
        from whose name its time information can be obtained by getTimeFromNameFunc
    Parameters:
        path: a string defining a directory
        timeRange: a list of two datetime.datetime objects
        getTimeFromNameFunc: a function whose input is a string of file name and output is a datetime.datetime object
        strings: a list of strings that should be contained in the objective file name.
        size: a string in the form of "<num" or ">num", where num is an integer whose unit is bytes. This parameter define either that the size of the file should be less than num bytes, or that its size greater than num bytes. The default value for this parameter is "allSize" which allow for all file size.
    '''
    logging.debug('find name in time range:')
    logging.debug(timeRange)
    if os.path.exists(path):
        with os.scandir(path) as it:
            fileNamesMeetingSize = []
            foundFileNames = []
            for entry in it:
                if entry.is_file():
                    _, fileExt = os.path.splitext(entry.name)
                    if ext and fileExt.lower() != ext:
                        continue
                    if size[0] == '<' and entry.stat().st_size < int(size[1:]):
                        fileNamesMeetingSize.append(entry.name)
                    elif size[0] == '>' and entry.stat().st_size > int(size[1:]):
                        fileNamesMeetingSize.append(entry.name)
                    elif size == 'allSize':
                        fileNamesMeetingSize.append(entry.name)
            logging.debug("number of files with extension {ext} and meeting size criterion: {nf}".format(ext=ext, nf=len(fileNamesMeetingSize)))
        fileNamesAfterSelection_ = fileNamesMeetingSize
        if strings:
            fileNamesAfterSelection = []
            if not isinstance(strings[0], list):
                logging.debug("string criteria for searching:")
                logging.debug(strings)
                while fileNamesAfterSelection_:
                    fileName_ = fileNamesAfterSelection_.pop(0)
                    if all(string in fileName_ for string in strings):
                        fileNamesAfterSelection.append(fileName_)
            else:
                while fileNamesAfterSelection_:
                    fileName_ = fileNamesAfterSelection_.pop(0)
                    for stringList in strings:
                        if all(string in fileName_ for string in strings):
                            fileNamesAfterSelection.append(fileName_)
                            break
            fileNamesAfterSelection_ = fileNamesAfterSelection
        fileNamesAfterSelection = []
        for fileName in fileNamesAfterSelection_:
            datetimeOfFile = getTimeFromNameFunc(fileName)
            if datetimeOfFile >= timeRange[0] and datetimeOfFile <= timeRange[1]:
                fileNamesAfterSelection.append(fileName)
        foundFileNames = fileNamesAfterSelection
    else:
        logging.debug("path does not exist: {}".format(path))
        foundFileNames = []
    return foundFileNames

##
def findFileNames(path, strings=None, size='allSize', timeTag=None, ext='.cdf'):
    '''
    Purpose:
        This function is to find in the directory defined by <path> the names of files which according to 
    Parameters:
        path: a string defining a directory
        strings: a list of strings that should be contained in the objective file name.
        size: a string in the form of "<num" or ">num", where num is an integer whose unit is bytes. This parameter define either that the size of the file should be less than num bytes, or that its size greater than num bytes. The default value for this parameter is "allSize" which allow for all file size.
        timeTag: a datetime object that specify the time accurate to second. This parameter is designed for mms brst data which is segmented over a few minutes.
    '''
    if timeTag:
        logging.debug('time tag')
        logging.debug(timeTag)
    logging.debug('size: {}'.format(size))
    if os.path.exists(path):
        with os.scandir(path) as it:
            fileNamesMeetingSize = []
            foundFileNames = []
            for entry in it:
                if entry.is_file():
                    _, fileExt = os.path.splitext(entry.name)
                    if ext and fileExt.lower() != ext:
                        continue
                    if size[0] == '<' and entry.stat().st_size < int(size[1:]):
                        fileNamesMeetingSize.append(entry.name)
                    elif size[0] == '>' and entry.stat().st_size > int(size[1:]):
                        fileNamesMeetingSize.append(entry.name)
                    elif size == 'allSize':
                        fileNamesMeetingSize.append(entry.name)
            logging.debug("number of files with extension {ext} and meeting size criterion: {nf}".format(ext=ext, nf=len(fileNamesMeetingSize)))
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
        logging.info('finding files for dataset and variables: ')
        logging.info(datasetAndVariables)
        splitedPathToDataset = datasetAndVariables[:-1]
        variableNames = datasetAndVariables[-1]
        start_ = start
        variablesInADatasetOverDatetimeRange = [] # first index for data from different files over a epoch range, second index for variables
        variablesInADatasetIndependantOnTime = []
        while start_ < end:
            start_NextDayDateTime = start_ + timedelta(days=1)
            endOfTheDay = datetime(start_NextDayDateTime.year, start_NextDayDateTime.month, start_NextDayDateTime.day)
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
    if variables:
        varInfo, varInfoDict = readCDFInfo(cdfFile)
        epochDataInd = 0 # in most cases, epoch is the first zVariable
        dataMajor = {}
        dataAux = {} # data not dependant on epoch
        timeRange = [ot.datetime2list(dateTime, epochType=epochType) for dateTime in datetimeRange]
        for var in variables:
            logging.debug('var: ' + var)
            majorData = True
            depend0 = varInfoDict[var]['varAtts'].get('DEPEND_0', None)
            logging.info('depend0: '+ str(depend0))
            if depend0 is None:
                if varInfo[epochDataInd]['varInfo'].Variable == var:
                    pass
                else:
                    majorData = False
            else:
                assert depend0 == varInfo[epochDataInd]['varInfo'].Variable
            logging.info('isMajorData: '+str(majorData))
            if majorData:
                if timeRange is not None:
                    dataMajor[var] = cdfFile.varget(var, starttime=timeRange[0], endtime=timeRange[1])
                    logging.debug('data type: {}'.format(str(type(dataMajor[var]))))
                    if dataMajor[var] is None:
                        dataMajor[var] = np.array([])
                else:
                    dataMajor[var] = np.array([])
                    raise Exception('time range is None')
            else:
                dataAux[var] = cdfFile.varget(var)
        return dataMajor, dataAux
    else:
        cdfInfo = cdfFile.cdf_info()
        dataFromACdfFile = {}
        for var in cdfInfo.zVariables:
            dataFromACdfFile[var] = cdfFile.varget(var)
        return dataFromACdfFile


def readPDSData(fileName, dataFileExtension='.TAB', infoFileExtension='.xml', sep=None):
    '''
    Purpose: read data file from PDS (https://pds-ppi.igpp.ucla.edu/)
    Parameters:
        fileName: path to file without extension of the file.
        dataFileExtension: this parameter is deprecated. The data file will be inferred from the info file.
        sep: if dataFileExtension.lower() == '.tab', the separator used in the tab file
    '''
    fileDir, filebaseName = os.path.split(fileName)
    infoFile = fileName + infoFileExtension
    columnNames = []
    dataDict = {}
    if infoFileExtension == '.xml':
        xmlTree = ET.parse(infoFile)
        root = xmlTree.getroot()
        for l1 in root:
            if 'File_Area_Observational' in l1.tag:
                for l2 in l1:
                    if 'File' in l2.tag:
                        file_ = l2
                        for child_ in file_:
                            if 'file_name' in child_.tag:
                                dataFileName = child_.text
                    if 'Table_Character' in l2.tag:
                        for l3 in l2:
                            if 'Record_Character' in l3.tag:
                                record_character = l3
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
    elif infoFileExtension == '.LBL':
        with open(infoFile, 'r') as f:
            info = f.readlines()
        lblObjectDict = {'objectCounts': 0}
        lblObjectDict_ = lblObjectDict
        for lineInd, line in enumerate(info):
            if '=' in line:
                split_parts_ = line.split('=')
                field_ = split_parts_[0]
                value_ = '='.join(split_parts_[1:])
                field, value = field_.strip().strip('"'), value_.strip().strip('"')
                if field == "END_OBJECT":
                    assert value == lblObjectDict_['objType']
                    lblObjectDict_ = lblObjectDict_['upperDict']
                elif field == "OBJECT":
                    newObj = {'objType': value,
                              'upperDict': lblObjectDict_,
                              'objectCounts': 0}
                    objectCounts = lblObjectDict_['objectCounts']
                    lblObjectDict_['objectCounts'] += 1
                    lblObjectDict_[str(objectCounts)] = newObj
                    lblObjectDict_ = newObj
                else:
                    lblObjectDict_[field.lower()] = value
        dataFileName = lblObjectDict['^table'].strip('"')
        fmtFileName = lblObjectDict['0'].get("^structure")
        if fmtFileName:
            fmtFileName = fmtFileName.strip('"')
            pathToFMTFile = os.path.join(fileDir, fmtFileName)
            with open(pathToFMTFile, 'r') as f:
                info = f.readlines()
            objectInd = -1
            objectDicts = []
            for lineInd, line in enumerate(info):
                if line.count('=') == 1:
                    field_, value_ = line.split('=')
                    field, value = field_.strip(), value.strip()
                    if field == 'OBJECT':
                        objectInd +=1
                        objectDict = {'fmtStartLine': lineInd}
                        objectDicts.append(objectDict)
                    elif field == 'END_OBJECT':
                        objectDict['fmtEndLine_inclusive'] = lineInd
                        if lineInd == len(info)-1:
                            objectInfo = info[objectDict['fmtStartLine']:]
                        else:
                            objectInfo = info[objectDict['fmtStartLine']:objectDict['fmtEndLine_inclusive']+1]
                        for objectLineInd, objectLine in enumerate(objectInfo):
                            if '=' in objectLine:
                                objectLineItems = objectLine.split('=')
                                if objectLineItems[0].isupper():
                                    objectField, objectValue = objectLineItems[0].strip(), '='.join(objectLineItems[1:]).strip()
                                    if objectField == 'DESCRIPTION':
                                        objectDict[objectField] = ''.join(objectInfo[objectLineInd:-1])
                                        break
                                    else:
                                        objectDict[objectField] = objectValue
            totalBytes = 0
            itemNums = []
            itemBytes = []
            dataNames = []
            fmt_i_list = []
            for objectDict in objectDicts:
                bytesNum = int(objectDict['BYTES'])
                totalBytes += bytesNum
                itemNum = int(objectDict.get('ITEMS', 1))
                itemNums.append(itemNum)
                itemBytes.append(objectDict.get('ITEMS_BYTES', objectDict['BYTES']))
                dataName = objectDict['NAME']
                dataNames.append(dataName)
                dataType_ = objectDict['DATA_TYPE']
                if '/*' in dataType_:
                    index_ = dataType_.find('/*')
                    dataType = dataType_[:index_].strip()
                    objectDict['DATA_TYPE'] = dataType
                else:
                    dataType = dataType_
                if dataType == 'PC_REAL':
                    fmt_i_list.extend('f'*itemNum)
                elif dataType == 'DATE':
                    fmt_i_list.extend('c'*bytesNum)
                    itemNums[-1] = bytesNum
                elif dataType == 'LSB_UNSIGNED_INTEGER':
                    if bytesNum == 1:
                        fmt_i_list.append('B')
                    elif bytesNum == 2:
                        fmt_i_list.append('H')
            fmt_i = '<'+''.join(fmt_i_list)
            fmtsz = struct.calcsize(fmt_i)
            assert fmtsz == totalBytes
        else:
            columnsDict = lblObjectDict['0']
            columnCounts = columnsDict['objectCounts']
            for columnInd in range(columnCounts):
                columnDict = columnsDict[str(columnInd)]
                columnName = columnDict['name']
                columnNames.append(columnName)
                toDelKeys = ['objType', 'upperDict', 'objectCounts']
                for key in toDelKeys:
                    columnDict.pop(key, None)
                dataDict[columnName] = columnDict
    dataFilePath = os.path.join(fileDir, dataFileName)
    _, dataFileExtension = os.path.splitext(dataFilePath)
    if dataFileExtension.lower() == '.tab':
        data_type = {}
        for columnInd, columnName in enumerate(columnNames):
            dataType = dataDict[columnName]['data_type']
            if dataType in ['ASCII_Date_Time_YMD_UTC', 'TIME']:
                timeName = columnName
                timeType = dataType
            elif dataType == 'ASCII_String':
                pass
            elif dataType in ['ASCII_Integer', 'ASCII_Real']:
                data_type[columnName] = np.float64
        if sep is None:
            sep = '\s+'
        data_ = pd.read_table(dataFilePath, sep=sep, names=columnNames, dtype=data_type)
        if timeType == 'ASCII_Date_Time_YMD_UTC':
            epoch = cdflib.cdfepoch.parse(list(data_[timeName].str[:-1]))
        elif timeType == 'TIME':
            timeString_ = data_[timeName].iloc[0]
            timeStringLen = len(timeString_)
            head_, tail_ = timeString_.split('T')
            number_of_dashes = list(head_).count('-')
            if number_of_dashes == 1:
                stringExample = '2007-173T00:00:00'
                t_strings = list(data_[timeName])
                epoch = np.zeros(len(t_strings))
                fmt = '%Y-%jT%H:%M:%S'
                for tInd, t_string in enumerate(t_strings):
                    epoch[tInd] = dat.datetime2epoch(datetime.strptime(t_string, fmt))
            elif number_of_dashes == 2:
                if timeStringLen < 24:
                    stringExample = '2007-06-23T00:00:00.000'
                    supp = stringExample[timeStringLen:]
                else:
                    raise Exception('time string too long: {}'.format(timeString_))
                print(list(data_[timeName] + supp)[:10])
                epoch = cdflib.cdfepoch.parse(list(data_[timeName] + supp))
        data_[timeName] = epoch
        for key in dataDict.keys():
            dataDict[key]['data'] = data_[key].to_numpy()
    elif dataFileExtension.lower() == '.dat':
        dataEntryList = []
        with open(dataFilePath, 'rb') as f:
            while True:
                entry = f.read(fmtsz)
                if not entry:
                    break
                data = struct.unpack(fmt_i, entry)
                dataEntryList.append(data)
        itemInd = 0
        for objectInd, objectDict in enumerate(objectDicts):
            itemNum = itemNums[objectInd]
            entryIndNextStart = itemInd + itemNum
            numberOfRecords = len(dataEntryList)
            dataList = []
            for entryInd, entry in enumerate(dataEntryList):
                dataEntry = entry[itemInd:entryIndNextStart]
                dataList.append(dataEntry)
            if objectDict['DATA_TYPE'] == 'DATE':
                dataArray = np.zeros(numberOfRecords)
                for entryInd in range(numberOfRecords):
                    dateStr = ''.join([ele.decode() for ele in dataList[entryInd]])
                    year, dayInfo = dateStr.split('-')
                    dayOfYear, timeInfo = dayInfo.split('T')
                    hour, minute, secondInfo = timeInfo.split(':')
                    second, millisecond = secondInfo.split('.')
                    dataArray[entryInd] = cdflib.cdfepoch.compute_epoch(ot.datetime2list(datetime(int(year), 1, 1, int(hour), int(minute), int(second), int(millisecond)*1000) + timedelta(days=int(dayOfYear)-1)))
            else:
                dataArray = np.array(dataList)
            objectDict['data'] = dataArray
            columnName = objectDict['NAME']
            columnNames.append(columnName)
            dataDict[columnName] = objectDict
            itemInd = entryIndNextStart

#    with open(dataFile, 'r') as f:
#        info = f.readlines()
#    for lineInd, line in enumerate(info):
#        print(lineInd)
#        if sep is None:
#            lineInfo = info[lineInd].strip().split()
#        else:
#            lineInfo = info[lineInd].strip().split(sep=sep)
#        for columnInd, columnInfo in enumerate(lineInfo):
#            columnName = columnNames[columnInd]
#            dataType = dataDict[columnName]['data_type']
#            if dataType == 'ASCII_Date_Time_YMD_UTC':
#                data_ = ascii_date_time_ymd_utc2epoch(columnInfo)
#            elif dataType == 'TIME':
#                data_ = cdflib.cdfepoch.parse(columnInfo+'00')
#
#            elif dataType == 'ASCII_String':
#                data_ = columnInfo
#            elif dataType in ['ASCII_Integer', 'ASCII_Real']:
#                data_ = float(columnInfo)
#            if lineInd == 0:
#                dataDict[columnName]['data'] = []
#            dataDict[columnName]['data'].append(data_)
#    for key in dataDict.keys():
#        dataDict[key]['data'] = np.array(dataDict[key]['data'])
    return dataDict, columnNames


def ascii_date_time_ymd_utc2epoch(datetimeStr=None):
    '''
    Purpose: transform datetime
    Parameters:
        datetimeStr: e.g. 1979-06-20T01:12:34.332Z
    Return:
        epoch: 
    '''
    yearmmdd, hhmmss = datetimeStr[:-1].split('T')
    year, month, day = yearmmdd.split('-')
    hour, minute, secondInfo = hhmmss.split(':')
    second, millisecond = secondInfo.split('.')
    epoch = cdflib.cdfepoch.compute_epoch(ot.datetime2list(datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), int(millisecond)*1000)))
    return epoch


def readHoriaonsData(dataFilePath):
    with open(dataFilePath, 'r') as f:
        info_ = f.readlines()

    info = []
    inTable = False
    for ind, line in enumerate(info_):
        if line.strip() == '$$EOE':
            inTable = False
        if inTable:
            info.append(line)
        if line.strip() == '$$SOE':
            inTable = True
            tableStartLine = ind
    columnNames_ = info_[tableStartLine-2].split(',')
    columnNames = [name.strip() for name in columnNames_]

    dataDict = {}
    for lineInd, line in enumerate(info):
        lineInfo = info[lineInd].strip(', \n\t').split(',')
        lineInfo = [infoItem.strip() for infoItem in lineInfo]
        for columnInd, columnInfo in enumerate(lineInfo):
            columnName = columnNames[columnInd]
            if columnName == 'Calendar Date (TDB)':
                columnInfo = columnInfo[:-5]
                fmt = 'A.D. %Y-%b-%d %H:%M:%S'
                data_ = dat.datetime2epoch(datetime.strptime(columnInfo, fmt))
            else:
                data_ = float(columnInfo)
            if lineInd == 0:
                dataDict[columnName] = {'data': []}
            dataDict[columnName]['data'].append(data_)
    for key in dataDict.keys():
        dataDict[key]['data'] = np.array(dataDict[key]['data'])
    return dataDict, columnNames

