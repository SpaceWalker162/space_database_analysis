__author__ = 'Yufei Zhou'
## This file contains functions and classes related to work with a local database. In case of insufficiency of the local database, access to remote database and internet may be needed

import numpy as np
import pandas as pd
#import constants as con
import cdflib    # see github.com/MAVENSDC/cdflib
#import tarfile
from datetime import datetime, timedelta, timezone
import shutil
import os
import sys
import subprocess
from cycler import cycler
from itertools import combinations
import json
import cdasws
import functools
from urllib.parse import urlparse
import urllib3
import contextlib
import struct
from pprint import pprint
from scipy.signal import butter, lfilter, freqz
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
import scipy.interpolate as interpolate
import scipy.special
import xml.etree.ElementTree as ET
import logging
from wolframclient.language import wl
from wolframclient.evaluation import WolframLanguageSession
import space_database_analysis.otherTools as ot
import space_database_analysis.databaseTools as dbt
import space_database_analysis.dataAnalysisTools as dat
#import space_database_analysis.examples.makeDatasetInfo
#import space_database_analysis.examples.cdasws_example
#from space_database_analysis.databaseTools import *

'''
<A, B> means either A or B
[A, B] means both A and B
'''

_mms_edp_bitmask = np.array([0, 1, 0, 1, 1, 1, 2, 0, 2, 1, 2, 0]) # see Section 11.5.1 Bitmasks and Quality Flgs, Calibration and Measurement Algorithms Document for MMS Project

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
        'omni': {
            'epochType': 'CDF_EPOCH',
            'multispacecraftMission': False},
        }

_epoch_types = ['CDF_EPOCH', 'CDF_EPOCH16', 'CDF_TIME_TT2000']

##
class Database(dbt.Database):
    '''

    '''
    def __init__(self, path):
        self.path = path

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


class FileList(list):
    '''
    file list for files in a database
    '''
    def __init__(self, filenames):
        super().__init__(self.classify_files(filenames))

    @staticmethod
    def _select_prefered_files(filenames):
        filenames_prefered = []
        for filename in filenames:
            # ext should be '.cdf'
            filename_base, ext = os.path.splitext(filename)
            if filename_base[-2:] == '.z':
                filenames_prefered.append(filename)
            else:
                if filename_base+'.z.cdf' in filenames:
                    pass
                else:
                    filenames_prefered.append(filename)
        return filenames_prefered

    @staticmethod
    def _compare_file_version(filename1, filename2):
        '''
        Parameters:
            filename1: str
            filename2: str
        '''
        if filename1 == filename2:
            return ot.Same
        filename1, _ = os.path.splitext(filename1)
        filename2, _ = os.path.splitext(filename2)
        filename1, ext1 = os.path.splitext(filename1)
        filename2, ext2 = os.path.splitext(filename2)
        if ext1 == ext2:
            raise Exception('filename mismatch')
        elif ext1 == '.z': # filename1 precedes filename2
            return ot.Precede
        elif ext2 == '.z': # filename2 precedes filename1
            return ot.Follow

    @classmethod
    def _sort_file_by_version(cls, filenames):
        sorted_files = [filenames[0]]
        for ind, filename in enumerate(filenames[1:]):
            direction_func = functools.partial(cls._compare_file_version, filename)
            pos = ot.bisect_solution(sorted_files, direction_func)
            sorted_files.insert(pos, filename)
        return sorted_files


    @classmethod
    def classify_files(cls, filenames):
        filenames = sorted(filenames)
        classified = []
        last_filenamebase = ''
        for ind in range(len(filenames)):
            filename = filenames[ind]
            filenamebase_, _ = os.path.splitext(filename)
            filenamebase_, z_ext = os.path.splitext(filenamebase_)
            components = filenamebase_.split('_')
            version = components[-1]
            filenamebase = '_'.join(components)
            if filenamebase == last_filenamebase:
                file_cat.append(filename)
            else:
                file_cat = []
                classified.append(file_cat)
                file_cat.append(filename)
            last_filenamebase = filenamebase
        for ind, file_class in enumerate(classified):
            classified[ind] = cls._sort_file_by_version(file_class)
        return classified

    def refine_files(self):
        for ind, file_class in enumerate(self):
            self[ind] = [file_class[0]]

    @classmethod
    def _compare_file_list(cls, file_list, file_name):
        if [file_name] in file_list:
            return ot.Precede
        for file_class in file_list:
            cla = cls.classify_files([file_class[0], file_name])
            if len(cla) == 1:
                if cla[0][0] == file_name:
                    return
                elif cla[0][0] == file_class[0]:
                    return ot.Precede

    def compare_file(self, file_name, update=False):
        if [file_name] in self:
            if update:
                return
            else:
                return ot.Precede
        for file_class in self:
            cla = self.classify_files([file_class[0], file_name])
            if len(cla) == 1:
                if cla[0][0] == file_name:
                    if update:
                        file_class.insert(0, file_name)
                        return True
                    else:
                        return
                elif cla[0][0] == file_class[0]:
                    if update:
                        return
                    else:
                        return ot.Precede
        if update:
            self.append([file_name])
            return True

##
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
    def __init__(self, mission=None, name=None, instrumentationAll=None, data=None, dataCleaned=None, workDataDir=None, workDataDirsBak=None, epoch_type_to_load_data=None):
        '''
        epochType: defines the type of epoch to which the epoch in cdf files will be transformed


        '''
        self.mission = mission
        self.name = name
        self.instrumentationAll = instrumentationAll
        self.epoch_type_to_load_data = epoch_type_to_load_data
#        if epochType is None:
#            pass
#           epochType = missionInfo[self.mission]['epochType']
#        self.epochType = epochType
        if data is None:
            self.data = {}
        elif isinstance(data, dict):
            self.data = data
        else:
            raise Exception('wrong data input format')
        self.metadata = {}
        if dataCleaned is None:
            self.dataCleaned = {}
        elif isinstance(dataCleaned, dict):
            self.dataCleaned = dataCleaned
        else:
            raise Exception('wrong dataCleaned input format')
        self.workDataDir = workDataDir
        self.workDataDirsBak = workDataDirsBak

    def loadData(self, datetimeRange=None, instrumentation=None, instrumentationRetrivingName=None, datasets_variables_with_retrieving_names=None, variables=None, variableRetrivingNames=None, datasetAndVariables=None, variablesWithRetrivingNames=None, instrumentationVariablesWithRetrivingName=None, cleanData=False, tStepPecentageCriterion=0.9, lowpassCutoff=None, inPlace=False, gapThreshold=None, minNumberOfPoints=None, returnShiftQ=False, copy_if_not_exist=True, search_online=False, fromFile=None, useMask=False, maskValue=-1.0*10**30, sparse_factor=None, allow_nonidentital_supporting_data=False):
        '''
        Parameters:
            datetimeRange: a list of two datetime objects, representing the start and end of the data to be loaded.
            instrumentation: a list in terms of the directory names at all levels below spacecraft and above year or files. For example, ['fgm', 'brst', 'l2']
            instrumentationRetrivingName: a string.
            datasetAndVariables: this parameter should not be used by user.
            instrumentationVariablesWithRetrivingName: a dictionary in terms of the path to the datasets. Its leaf value is a list of variable names. The names are defined by the corresponding cdfFile. For example, {'Cluster': {'C1' : {'C1_CP_FGM_FULL': ['FGM', ('time_tags__C1_CP_FGM_FULL', 't'), ('B_vec_xyz_gse__C1_CP_FGM_FULL', 'B')]}}, 'mms': {'mms1': {'fgm': {'brst': {'l2': ['fgmBrst', ('Epoch', 't'), ('Btotal', 'BTotal')]}}}}}. Please note that 'Epoch' and 'Btotal' are improvised. This parameter may also be a list of lists, with each sublist in the form of ['Cluster', 'C1', 'C1_CP_FGM_FULL', ['FGM', ('time_tags__C1_CP_FGM_FULL', 't'), ('B_vec_xyz_gse__C1_CP_FGM_FULL', 'B')]]. To retrieve data, for example, of 'B_vec_xyz_gse__C1_CP_FGM_FULL', use C1.data['FGM']['B']
            datasets_variables_with_retrieving_names: a dictionary of datasets. Its leaf value is a list of variable names. The names are defined by the corresponding cdfFile. For example, {'MMS1_FGM_SRVY_L2': ['FGM', ('Epoch', 't'), ('mms1_fgm_srvy_l2_bvec_gse', 'B')], 'MMS2_FPI_FAST_L2_DIS-MOMS': ['FPI', ('Epoch', 't'), ('mms2_fpi_fast_l2_dis-moms_density', 'n'), 'mms2_fpi_fast_l2_dis-moms_velocity_gse']}. Please note that 'Epoch' and 'Btotal' are improvised. To retrieve data 'B_vec_xyz_gse__C1_CP_FGM_FULL', use spacecraft.data['FGM']['B'], to retrieve 'mms2_fpi_fast_l2_dis-moms_velocity_gse', use spacecraft.data['FPI']['mms2_fpi_fast_l2_dis-moms_velocity_gse']
            fromFile: to load data from a file, the file name is given by this parameter.
            sparse_factor: When loading a large chunk of high resolution data it is sometimes ideal for the purpose of saving memory to take only one record, say, every 1000 records. If sparse_factor is None, the full data will be loaded. Otherwise it should be an integer such as 1000 to specify the step in loading data
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
            dataFromACdfFile, varInfoDict = readDataFromACdfFile(cdfFile, variableNames)
            self.data.update({instrumentationRetrivingName: dataFromACdfFile})
            self.metadata.update({instrumentationRetrivingName: varInfoDict})
        else:
            assert datasets_variables_with_retrieving_names
            datasets_info = loadDatasets_info(self.workDataDir, self.workDataDirsBak, copy_if_not_exist=copy_if_not_exist)
            for datasetID, item in datasets_variables_with_retrieving_names.items():
                datasetRetrievingName = item[0]
                variableNamesAndRetrievingNames_ = item[1:]
                variableNamesAndRetrievingNames = []
                for varRet in variableNamesAndRetrievingNames_:
                    if type(varRet) is str:
                        variableNamesAndRetrievingNames.append((varRet, varRet))
                    elif type(varRet) is tuple:
                        variableNamesAndRetrievingNames.append(varRet)
                variableNames = [varRet[0] for varRet in variableNamesAndRetrievingNames]
                dataset_info = datasets_info.get(datasetID)
                if dataset_info:
                    dataset = Dataset(dataset_info=dataset_info, databasePath=self.workDataDir, databaseBakPaths=self.workDataDirsBak)
                else:
                    dataset = Dataset(datasetID=datasetID, databasePath=self.workDataDir, databaseBakPaths=self.workDataDirsBak)
                dataset.load_data(variableNames=variableNames, datetimeRange=datetimeRange, copy_if_not_exist=copy_if_not_exist, search_online=search_online, sparse_factor=sparse_factor, allow_nonidentital_supporting_data=allow_nonidentital_supporting_data)
                if self.epoch_type_to_load_data:
                    for varName in variableNames:
                        data_type = dataset.varInfoDict[varName]['varInfo'].Data_Type_Description
                        if data_type in _epoch_types:
                            if not data_type == self.epoch_type_to_load_data:
                                para_ = {data_type: dataset.data[varName]}
                                dataset.data[varName] = dat.Epochs(**para_).get_data(fm=self.epoch_type_to_load_data)
                for varName, retName in variableNamesAndRetrievingNames:
                    try:
                        dataset.data[retName] = dataset.data.pop(varName) # 'pop' was 'get', to avoid multiple occurrence in the returned data dict we use pop instead.
                    except:
                        dataset.data[retName] = dataset.data.get(varName)
                    dataset.varInfoDict[retName] = dataset.varInfoDict.get(varName)
                self.data.update({datasetRetrievingName: dataset.data})
                self.metadata.update({datasetRetrievingName: dataset.varInfoDict})
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

    def maskData(self, instrumentations=None, value=-1.0*10**31):
        if instrumentations is not None:
            pass
        else:
            instrumentations = self.data.keys()
        for instrumentation in instrumentations:
            dic = self.data[instrumentation]
            for key, varData in dic.items():
                if np.issubdtype(varData.dtype, np.number):
                    varData = np.ma.masked_equal(varData, value)
                    varData = varData.filled(np.nan)
                else:
                    pass
                dic[key] = varData


class Dataset:
    def __init__(self, datasetID=None, dataset_info=None, databasePath=None, databaseBakPaths=None, copy_if_not_exist=True):
        '''
        Expected database directory:
            database:
                data:
                    mms
                    cluster
                software:
                    something
        '''
        self.databasePath = databasePath
        self.database = dbt.Database(databasePaths=[databasePath])
        self.databaseBakPaths = databaseBakPaths
        if databaseBakPaths is not None:
            self.databaseBak = dbt.Database(databasePaths=databaseBakPaths)
            self.dataBakPaths = [os.path.join(databaseBakPath, 'data') for databaseBakPath in self.databaseBakPaths]
        else:
            self.databaseBak = None
            self.dataBakPaths = None
        self.dataPath = os.path.join(self.databasePath, 'data')
        if dataset_info:
            self.datasetID = dataset_info['Id']
            self.dataset_info = dataset_info
        else:
            assert datasetID
            self.datasetID = datasetID
            datasets_info = loadDatasets_info(self.databasePath, self.databaseBakPaths, copy_if_not_exist=copy_if_not_exist)
            self.dataset_info = datasets_info.get(datasetID, {})

        dataset_file_time_gap_str = self.dataset_info.get('dataset_file_time_gap')
        if dataset_file_time_gap_str:
            if dataset_file_time_gap_str == 'irregular':
                dataset_file_time_gap = 'irregular'
#            logging.info(dataset_file_time_gap_str)
            else:
                num, unit = dataset_file_time_gap_str.split(' ')
                num = int(num)
                if unit == 'hour':
                    unit = timedelta(seconds=3600)
                elif unit == 'day':
                    unit = timedelta(days=1)
                elif unit == 'month':
                    unit = ['month']
                dataset_file_time_gap = num * unit
        else:
            dataset_file_time_gap = timedelta(days=1)
        self.dataset_file_time_gap = dataset_file_time_gap
        logging.debug('dataset file time gap: {}'.format(dataset_file_time_gap))

        self.datasetPathInDatabaseData = self.dataset_info.get('dataset_path')
        if not self.datasetPathInDatabaseData:
            self.datasetPathInDatabaseData = os.path.join(*self.datasetID.lower().split('_'))
        self.datasetAbsolutePath = os.path.join(self.dataPath, self.datasetPathInDatabaseData)
        if self.dataBakPaths is not None:
            self.datasetAbsoluteBakPaths = [os.path.join(dataBakPath, self.datasetPathInDatabaseData) for dataBakPath in self.dataBakPaths]
        else:
            self.datasetAbsoluteBakPaths = None

    def _get_file_time_limits(self, dateTime):
        if isinstance(self.dataset_file_time_gap, timedelta):
            if self.dataset_file_time_gap == timedelta(days=1):
                beginOfTheFilePeriod = datetime(dateTime.year, dateTime.month, dateTime.day, tzinfo=timezone.utc)
                endOfTheFilePeriod = beginOfTheFilePeriod + timedelta(days=1)
            elif self.dataset_file_time_gap < timedelta(days=1):
                beginOfTheDay = datetime(dateTime.year, dateTime.month, dateTime.day, tzinfo=timezone.utc)
                completePeriodBefore = (dateTime - beginOfTheDay)//self.dataset_file_time_gap
                beginOfTheFilePeriod = beginOfTheDay + completePeriodBefore * self.dataset_file_time_gap
                endOfTheFilePeriod = beginOfTheDay + (completePeriodBefore+1)*self.dataset_file_time_gap 
        elif isinstance(self.dataset_file_time_gap, list):
            assert self.dataset_file_time_gap[0] == 'month'
            beginOfTheFilePeriod = datetime(dateTime.year, dateTime.month, 1, tzinfo=timezone.utc)
            if dateTime.month == 12:
                endOfTheFilePeriod = datetime(dateTime.year+1, dateTime.month+len(self.dataset_file_time_gap)-12, 1, tzinfo=timezone.utc)
            else:
                endOfTheFilePeriod = datetime(dateTime.year, dateTime.month+len(self.dataset_file_time_gap), 1, tzinfo=timezone.utc)
        return beginOfTheFilePeriod, endOfTheFilePeriod

    def _define_search_criteria(self, beginOfTheFilePeriod=None, endOfTheFilePeriod=None, dateTime=None, size='allSize', **para):
        search_criteria = {}
        strings = []
        search_criteria.update({'strings': strings, 'size': size})
        file_naming_convention = self.dataset_info.get('file_naming_convention')
        if file_naming_convention:
            file_naming_convention = beginOfTheFilePeriod.strftime(file_naming_convention.format(datasetIDLower=self.datasetID.lower()))
            strings.append(file_naming_convention)
        else:
            strings.append(self.datasetID.lower())
            strings.append(beginOfTheFilePeriod.strftime("%Y%m%d"))
#        search_criteria.update({'searchMethod': searchMethod, 'strings': criteria, 'timeTag': timeTag, 'size': size})

#        if dateTime:
#            start = dateTime
#            end = start + timedelta(seconds=1)
#        elif beginOfTheFilePeriod:
#            start = beginOfTheFilePeriod
#            end = start + timedelta(seconds=1)
#        criteria = []
#        timeTag = None
#        logging.debug("mission: {}".format(self.mission))
#        if self.mission == 'mms':
#            if 'brst' in self.instrumentationPath:
#                searchMethod = 'allFilesInTimeRange'
#            else:
#                searchMethod = 'general'
#                criteria.extend(self.instrumentation.copy())
#                if self.dataset_file_time_gap is not None:
#                    if self.dataset_file_time_gap < timedelta(days=1):
#                        criteria.append(self.filePeriodRange[0].strftime("%Y%m%d%H"))
#                    elif self.dataset_file_time_gap == timedelta(days=1):
#                        criteria.append(start.strftime("%Y%m%d"))
#        elif self.mission == 'cluster':
#            if 'cis-hia' in self.instrumentationPath:
#                searchMethod = 'general'
#                infoList_ = self.instrumentation[1].split('_')
#                mode = infoList_[2][:-4]
#                sensitivity = infoList_[3][0]
#                dataName_ = infoList_[4]
#                if dataName_ == 'phasespacedens':
#                    dataNameAbb = 'psd'
#                elif dataName_ == 'diffenergyflux':
#                    dataNameAbb = 'pef'
#                restructedStrCriterion = '_'.join(['cis-hia', sensitivity+'s', mode, 'ions', dataNameAbb])
#                criteria.append(restructedStrCriterion)
#                criteria.append('v2022')
#                criteria.append(start.strftime("%Y%m%d"))
#            else:
#                searchMethod = 'ClusterCAA'
#                interval = [start, end]
#                search_criteria.update({'interval': interval})
#                logging.debug('ClusterCAA searching criteria (instrumentation)')
#                logging.debug(self.instrumentation)
#                criteria.extend(self.instrumentation.copy()[1:])
#                criteria.append(start.strftime("%Y%m%d"))
#        elif self.mission == 'ace':
#            searchMethod = 'general'
#            criteria.append(start.strftime("%Y%m%d"))
#        elif self.mission == 'cassini':
#            searchMethod = 'general'
#            if 'CO-E_SW_J_S-MAG-4-SUMM-1MINAVG-V2.0' in self.instrumentation:
#                criteria.append(start.strftime("%Y"))
#            else:
#                criteria.append(start.strftime("%Y%m%d"))
#        else:
#            searchMethod = 'general'
#            raise Exception('mission not defined!')

#        search_criteria.update({'searchMethod': searchMethod, 'strings': criteria, 'timeTag': timeTag, 'size': size})
        return search_criteria

    def _get_file_paths(self, absolutePathToDataset, dateTime=None, search_func=None, **search_criteria):
        '''
        find file path under dataset path
        Parameters:
            search_func: a function having input [directory, **search_criteria] to search for the file in the directory meeting criteria defined by search_criteria.
        '''
        fileNames = search_func(absolutePathToDataset, **search_criteria)
        absolutePathToFile = absolutePathToDataset
        if len(fileNames) == 0:
            absolutePathToDatasetYear = os.path.join(absolutePathToDataset, dateTime.strftime('%Y'))
            absolutePathToFile = absolutePathToDatasetYear
            logging.debug("Not found. Looking for data files in: {}".format(absolutePathToDatasetYear))
            fileNames = search_func(absolutePathToDatasetYear, **search_criteria)
            if len(fileNames) == 0:
                absolutePathToDatasetYearMonth = os.path.join(absolutePathToDatasetYear, dateTime.strftime('%m'))
                absolutePathToFile = absolutePathToDatasetYearMonth
                logging.debug("Not found. Looking for data files in: {}".format(absolutePathToDatasetYearMonth))
                fileNames = search_func(absolutePathToDatasetYearMonth, **search_criteria)
        filePaths = [os.path.join(absolutePathToFile, fileName) for fileName in fileNames]
        numFiles = len(filePaths)
#        filePath = None
        if numFiles == 0:
            logging.warning("No file was found under {} while looking for dataset: {}".format(absolutePathToDataset, self.datasetID))
#        elif numFiles == 1:
#            filePath = filePaths[0]
#            logging.info("file found: {}".format(filePath))
#        elif numFiles > 1:
#            logging.warning("More than one files were found in {}:" + ("\n{}"*numFiles).format(absolutePathToFile, *fileNames))
        return filePaths


    def _get_file_paths_from_multiple_sources(self, datetimeRange, dateTime=None, copy_if_not_exist=True, search_online=False, search_method='regular', **para):
        '''
        Parameters:
            search_method: when searching for files of a dataset whose interval are regularlly divided, let search_method='regular'. While the intervals are not regularlly distributed, such as those of mms brst data, use 'irregular'.

        '''
        if search_method == 'regular':
            search_func = findFileNames
            beginOfTheFilePeriod, endOfTheFilePeriod = self._get_file_time_limits(datetimeRange[0])
            search_criteria = self._define_search_criteria(beginOfTheFilePeriod, endOfTheFilePeriod, dateTime=None, size='allSize', **para)
            logging.debug('looking for file with string criteria: {}'.format(search_criteria['strings']))
        elif search_method == 'irregular':
            search_func = findFileNamesInTimeRange
            timedeltas = [-timedelta(seconds=3600*1), timedelta(seconds=3600*1)]
            search_criteria = {'timeRange': [datetime_+timedelta_ for datetime_, timedelta_ in zip(datetimeRange, timedeltas)], 'getTimeFromNameFunc': getTimeFromName, 'strings': [self.datasetID.lower()]}
            logging.debug('looking for file with string criteria: {} in datetime range {}/{}'.format(search_criteria['strings'], *search_criteria['timeRange']))

        filePaths = FileList(self._get_file_paths(self.datasetAbsolutePath, dateTime=datetimeRange[0], search_func=search_func, **search_criteria))
        filePaths.refine_files()
        filePaths = [file_[0] for file_ in filePaths]
        file_list = FileList([os.path.relpath(path_, self.datasetAbsolutePath) for path_ in filePaths])

        if self.dataBakPaths:
            for datasetAbsoluteBakPath in self.datasetAbsoluteBakPaths:
                fileBakPaths = self._get_file_paths(datasetAbsoluteBakPath, dateTime=datetimeRange[0], search_func=search_func, **search_criteria)
                for fileBakPath in fileBakPaths:
                    fileBakPathUnderDataset = os.path.relpath(fileBakPath, datasetAbsoluteBakPath)
                    destFilePath = os.path.join(self.datasetAbsolutePath, fileBakPathUnderDataset)
                    not_exist_Q = file_list.compare_file(fileBakPathUnderDataset, update=True)
                    if not_exist_Q:
                        if copy_if_not_exist:
                            logging.info('now copying {} to {} ...'.format(fileBakPath, destFilePath))
                            os.makedirs(os.path.dirname(destFilePath), exist_ok=True)
                            shutil.copyfile(fileBakPath, destFilePath)
                            try:
                                shutil.copystat(fileBakPath, destFilePath)
                            except: pass
#                            shutil.copy2(fileBakPath, destFilePath)
                            logging.info('file copied')
                            filePaths.append(destFilePath)
                        else:
                            filePaths.append(fileBakPath)

        if search_online:
            filePaths = self.get_online_files(datetimeRange=datetimeRange, search_criteria=search_criteria)

        return filePaths

    def _get_online_file(self, fileURL):
        o = urlparse(fileURL)
        filePathInDatabase = os.path.join(*o.path.split('/')[2:])
        existingFilePath = ''
        try:
            existingFilePath = self.database._get_file_path(filePathInDatabase)
        except FileNotFoundError:
            try:
                existingFilePath = self.databaseBak._get_file_path(filePathInDatabase)
            except FileNotFoundError: pass
        if existingFilePath:
            return existingFilePath
        else: # download the online file
            destFilePath = os.path.join(self.database.path, filePathInDatabase)
            os.makedirs(os.path.dirname(destFilePath), exist_ok=True)
            logging.info('downloading {}\n from {}'.format(destFilePath, fileURL))
            http = urllib3.PoolManager()
            blocksize = 32768
            with open(destFilePath, 'wb') as f:
                with contextlib.closing(http.request("GET", fileURL, preload_content=False)) as resp:
                    for chunk in resp.stream(blocksize):
                        f.write(chunk)
            logging.info('downloaded')
            return destFilePath


    def get_online_files(self, datetimeRange, **para):
        cdaswsObj = cdasws.CdasWs()
        logging.info('looking for files from CDAWeb with dataset ID {} for time period {} -- {} ...'.format(self.datasetID, *datetimeRange))
        try:
            status, files = cdaswsObj.get_original_files(self.datasetID, *datetimeRange)
            allFileNames = [file_['Name'] for file_ in files]
            logging.info('files found from CDAWeb for the temporal period: {} -- {}\n {}'.format(*datetimeRange, allFileNames))
            file = None
            search_criteria = para.get('search_criteria')
            if search_criteria:
                logging.debug('search online files with string criteria: {}'.format(search_criteria['strings']))
                downloadedFiles = []
                for file in files:
                    strings = search_criteria['strings']
                    if all(string in file['Name'] for string in strings):
                        downloadedFiles.append(self._get_online_file(file['Name']))
#                        return [self._get_online_file(file['Name'])]
            else:
                downloadedFiles = []
                for file in files:
                    downloadedFiles.append(self._get_online_file(file['Name']))
            return downloadedFiles
        except TypeError:
            pass


    def _load_data_from_files(self, filePaths, variableNames, datetimeRanges=None, sparse_factor=None):
        '''
        Parameters:
            variablesNames: a list of variable names. It can also be '!all' which will make the function to load all data.
        '''
        variablesInADatasetOverDatetimeRange = [] # first index for data from different files over a epoch range, second index for variables
        variablesInADatasetIndependantOnTime = []
        varInfoDict = {}
        for filePath, datetimeRange in zip(filePaths, datetimeRanges):
            if not filePath:
                logging.debug('data file not found')
            else:
                cdfFile = cdflib.CDF(filePath)
                logging.debug('reading data file: {}'.format(filePath))
                logging.debug('datetimeRange: {}'.format(datetimeRange))
                dataMajor, dataAux, varInfoDict = readDataFromACdfFile(cdfFile, variableNames, datetimeRange)
                logging.debug('reading data file done: {}'.format(filePath))
                if sparse_factor:
                    numberOfRecords = len(dataMajor[list(dataMajor.keys())[0]])
                    dataMajor = dat.mask_dict_of_ndarray(dataMajor, slice(0, numberOfRecords, sparse_factor), copy=True)
                variablesInADatasetOverDatetimeRange.append(dataMajor)
                variablesInADatasetIndependantOnTime.append(dataAux)
        return variablesInADatasetOverDatetimeRange, variablesInADatasetIndependantOnTime, varInfoDict

    def load_data(self, variableNames=None, datetimeRange=None, copy_if_not_exist=True, search_online=False, sparse_factor=None, allow_nonidentital_supporting_data=False):
        '''
        Parameter:
            sparse_factor: When loading a large chunk of high resolution data it is sometimes ideal for the purpose of saving memory to take only one record, say, every 1000 records. If sparse_factor is None, the full data will be loaded. Otherwise it should be an integer such as 1000 to specify the step in loading data
        '''
        self.data = {}
        self.varInfoDict = {}
        self.datetimeRange = datetimeRange
        self.variables_to_load = variableNames
        start, end = self.datetimeRange

        if isinstance(self.dataset_file_time_gap, str) and self.dataset_file_time_gap == 'irregular':
            filePaths = self._get_file_paths_from_multiple_sources(datetimeRange, copy_if_not_exist=copy_if_not_exist, search_online=search_online, search_method='irregular')
            datetimeRanges = [datetimeRange]*len(filePaths)
        else:
            start_ = start
            filePaths = []
            datetimeRanges = []
            while start_ < end:
                beginOfTheFilePeriod, endOfTheFilePeriod = self._get_file_time_limits(start_)
                logging.debug('begin/end of the file period: {}/{}'.format(beginOfTheFilePeriod, endOfTheFilePeriod))
                if endOfTheFilePeriod >= end:
                    datetimeRange_ = [start_, end]
                else:
                    datetimeRange_ = [start_, endOfTheFilePeriod]
                filePaths_ = self._get_file_paths_from_multiple_sources(datetimeRange=datetimeRange_, dateTime=None, copy_if_not_exist=copy_if_not_exist, search_online=search_online)
                if filePaths_:
                    filePaths.append(filePaths_[0])
                    datetimeRanges.append(datetimeRange_)
                start_ = endOfTheFilePeriod
        if not filePaths:
            logging.debug('data file not found')
        else:
            variablesInADatasetOverDatetimeRange, variablesInADatasetIndependantOnTime, varInfoDict = self._load_data_from_files(filePaths, variableNames, datetimeRanges=datetimeRanges, sparse_factor=sparse_factor)
            variables = {}
            if len(variablesInADatasetOverDatetimeRange) > 0:
                for var in variablesInADatasetOverDatetimeRange[0].keys():
                    vardata_to_concatenate = [varsOfATimeRange[var] for varsOfATimeRange in variablesInADatasetOverDatetimeRange if varsOfATimeRange[var].size>0]
                    if len(vardata_to_concatenate) > 0:
                        variables[var] = np.concatenate(vardata_to_concatenate, axis=0)
                    else:
                        variables[var] = np.array([])
                        logging.warning('no data found for {}'.format(var))
            else:
                logging.warning('no data over time found for this load of dataset')
            # Checking no temporal overlap
#            for var in variables.keys():
#                if varInfoDict[var]['varInfo'].Data_Type in [31, 32, 33]:
#                    # this variable is epoch
#                    t_diff = np.diff(variables[var])
            # Checking supporting data in all files are the same
            for fileInd in range(len(variablesInADatasetIndependantOnTime)-1):
                logging.debug(variablesInADatasetIndependantOnTime[fileInd].keys())
                for key in variablesInADatasetIndependantOnTime[fileInd].keys():
                    try:
                        assert np.all(variablesInADatasetIndependantOnTime[fileInd][key] == variablesInADatasetIndependantOnTime[fileInd+1][key])
                    except AssertionError as e:
                        logging.warning('Supporting data for {} not same in {} and {}'.format(key, *filePaths[fileInd:fileInd+2]))
                        if allow_nonidentital_supporting_data:
                            pass
                        else:
                            raise AssertionError

            variables.update(variablesInADatasetIndependantOnTime[0])
            self.data = variables
            self.varInfoDict = varInfoDict

    def load_data_example(self, **kwargs):
        '''
        '''
        datetimeRange = [datetime.fromisoformat(dtStr) for dtStr in self.dataset_info.get('example_time_interval')]
        self.load_data(variableNames=None, datetimeRange=datetimeRange, **kwargs)

    def findDataFiles(self):
        pass

    def get_data(var, datetimeRange=None):
        '''
        It might need to load data into the dataset multiple times. Each time the loaded data could be different. Therefore we need the capability to load data into the existing dataset. This function will return data of a variable in a time interval and will load data if that data is not in memory.
        Parameters:
            var: the name of the variable
            datetimeRange: if None, the data stored in the memory will be returned. Otherwise it should be a list of two datetime objects. If additional data need to be loaded into the memory to fulfill the datetime range, the function will do so and return.
        '''
        pass

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
def spacecraftsDataResampling(spacecrafts, source, dataName, timeName='t', resamplingT=None, standardSCInd=0):
    '''
    Parameters:
        resamplingT: an one-dimensional array of epoch to be used as the epoch for the resampled data.
        standardSCInd: an integer representing which spacecraft in the list of spacecrafts will provide its epoch as the standard one. If this is used, resamplingT will be
    '''
    dimensionOfData = spacecrafts[0].data[source][dataName].shape
    numberOfSpacecrafts = len(spacecrafts)
    if isinstance(standardSCInd, int):
        resamplingT = spacecrafts[standardSCInd].data[source][timeName]
        resamplingSCInds = np.delete(np.arange(numberOfSpacecrafts), standardSCInd)
    elif isinstance(resamplingT, np.ndarray):
        resamplingSCInds = np.arange(numberOfSpacecrafts)
    else:
        raise Exception('epoch for resampled data is not defined')
    dataAllSpacecrafts = np.zeros((len(resamplingT), numberOfSpacecrafts, *dimensionOfData[1:]))
    for scInd in resamplingSCInds:
        spacecraft = spacecrafts[scInd]
        t_ = spacecraft.data[source][timeName]
        data_ = spacecraft.data[source][dataName]
        dataAllSpacecrafts[:, scInd, ...] = dat.dataFillAndLowPass(t_, data_, resamplingT=resamplingT)
    if isinstance(standardSCInd, int):
        dataAllSpacecrafts[:, standardSCInd, ...] = spacecrafts[standardSCInd].data[source][dataName]
    return resamplingT, dataAllSpacecrafts




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
    if len(datasetsKeys) == 1 and datasetsKeys[0].lower() == 'cluster':
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
    try:
        info = cdfFile.cdf_info()
    except Exception as e:
        logging.warning("failed when reading info from CDF file {filename}".format(filename=cdfFile.file))
        raise e
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
    return datetime.strptime(timeString, fmt).replace(tzinfo=timezone.utc)

def findFileNamesInTimeRange(path, timeRange=None, getTimeFromNameFunc=None, strings=None, size='allSize', ext='.cdf', **para):
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
    def select_based_on_size(path, size='allSize', ext='.cdf'):
        with os.scandir(path) as it:
            fileNamesMeetingSize = []
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
        return fileNamesMeetingSize

    def select_based_on_name_strings(fileNames, strings):
        fileNamesAfterSelection = []
        if not isinstance(strings[0], list):
            logging.debug("string criteria for searching:")
            logging.debug(strings)
            while fileNames:
                fileName_ = fileNames.pop(0)
                if all(string in fileName_ for string in strings):
                    fileNamesAfterSelection.append(fileName_)
        else:
            while fileNames:
                fileName_ = fileNames.pop(0)
                for stringList in strings:
                    if all(string in fileName_ for string in strings):
                        fileNamesAfterSelection.append(fileName_)
                        break
        return fileNamesAfterSelection

    def select_based_on_name_time(fileNames, timeRange=None):
        fileNamesAfterSelection = []
        for fileName in fileNames:
            datetimeOfFile = getTimeFromNameFunc(fileName)
            if datetimeOfFile >= timeRange[0] and datetimeOfFile <= timeRange[1]:
                fileNamesAfterSelection.append(fileName)
        return fileNamesAfterSelection

    logging.debug('To find name in time range: {}/{} in {}'.format(*timeRange, path))
    if os.path.exists(path):
        fileNamesAfterSelection = select_based_on_size(path, size=size)
        if strings:
            if len(fileNamesAfterSelection) > 0:
                fileNamesAfterSelection = select_based_on_name_strings(fileNamesAfterSelection, strings=strings)
        if len(fileNamesAfterSelection) > 0:
            fileNamesAfterSelection = select_based_on_name_time(fileNamesAfterSelection, timeRange=timeRange)
        foundFileNames = fileNamesAfterSelection
    else:
        logging.debug("path does not exist: {}".format(path))
        foundFileNames = []
    logging.debug('files meeting all criteria: {}'.format(foundFileNames))
    return foundFileNames

##
def findFileNames(path, strings=None, size='allSize', timeTag=None, ext='.cdf', **para):
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
                                        timeOfFile = datetime(year, month, day, hour, minute, second, tzinfo=timezone.utc)
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
def findFileNamesInDatetimeRange(path, datetimeRange):
    '''
    This function is dedicated to find filenames in the form of mms1_fgm_brst_l2_20160904013434_v5.87.0.cdf 
    Parameters:
        path: the path to the folder that contains the file
        datetimeRange: a list of two datetime object, start and end, that specifies the temporal interval of the data.
    '''
    pass


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
def readDataFromACdfFile(cdfFile=None, variables=None, datetimeRange=None, cdf_file_path=None):
    if not cdfFile:
        cdfFile = cdflib.CDF(cdf_file_path)
    varInfo, varInfoDict = readCDFInfo(cdfFile)
    if variables:
        epochDataInd = 0 # in most cases, epoch is the first zVariable
        dataMajor = {}
        dataAux = {} # data not dependant on epoch
        for var in variables:
            fillval = varInfoDict[var]['varAtts'].get('FILLVAL', None)
            if fillval is None:
                pad = cdfFile.varinq(var).Pad
            else:
                pad = fillval
            logging.debug('var: ' + var)
            majorData = True
            depend0 = varInfoDict[var]['varAtts'].get('DEPEND_0', None)
            logging.debug('depend0: '+ str(depend0))
            if depend0 is None:
                if varInfoDict[var]['varInfo'].Data_Type_Description in _epoch_types:
                    epochDataName = var
                else:
                    majorData = False
            else:
#                assert depend0 == varInfo[epochDataInd]['varInfo'].Variable
                epochDataName = depend0
            logging.debug('isMajorData: '+str(majorData))
            if majorData:
                if datetimeRange is not None:
                    epochType = varInfoDict[epochDataName]['varInfo'].Data_Type_Description
                    timeRange = [ot.datetime2list(dateTime, epochType=epochType) for dateTime in datetimeRange]
                    varData = cdfFile.varget(var, starttime=timeRange[0], endtime=timeRange[1])
                    try:
                        if np.issubdtype(varData.dtype, np.number):
                            varData = np.ma.masked_equal(varData, pad)
                            varData = varData.filled(np.nan)
                    except: pass
                    dataMajor[var] = varData
                    logging.debug('data type: {}'.format(str(type(dataMajor[var]))))
                    if dataMajor[var] is None:
                        dataMajor[var] = np.array([])
                else:
                    dataMajor[var] = np.array([])
                    raise Exception('time range is None')
            else:
                dataAux[var] = cdfFile.varget(var)
        return dataMajor, dataAux, varInfoDict
    else:
        cdfInfo = cdfFile.cdf_info()
        dataFromACdfFile = {}
        for var in cdfInfo.zVariables:
            pad = cdfFile.varinq(var).Pad
            varData = cdfFile.varget(var)
            if np.issubdtype(varData.dtype, np.number):
                varData = np.ma.masked_equal(varData, pad)
                if np.issubdtype(varData.dtype, np.integer):
                    varData = varData.astype(np.float32)
                varData = varData.filled(np.nan)
            else:
                pass
            dataFromACdfFile[var] = varData
        return dataFromACdfFile, {}, varInfoDict

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
                    dataArray[entryInd] = cdflib.cdfepoch.compute_epoch(ot.datetime2list(datetime(int(year), 1, 1, int(hour), int(minute), int(second), int(millisecond)*1000, tzinfo=timezone.utc) + timedelta(days=int(dayOfYear)-1)))
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
    epoch = cdflib.cdfepoch.compute_epoch(ot.datetime2list(datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), int(millisecond)*1000, tzinfo=timezone.utc)))
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


def update_datasets_info(databasePath, databaseBakPaths):
    '''
    once the additional_datasets_info is updated on the server, user should run this function on the user side to get the new additional_datasets_info
    '''
    try:
        database = dbt.Database(databaseBakPaths)
        for path in databaseBakPaths:
            if path in database.additional_datasets_info_paths[0]:
                databasePathBak = path

        filePath = database.additional_datasets_info_paths[0]
        relpath = os.path.relpath(filePath, databasePathBak)
        destFilePath = os.path.join(databasePath, relpath)
        logging.debug('updating additional datasets info: rsync {} to {} ...'.format(filePath, destFilePath))
        os.makedirs(os.path.dirname(destFilePath), exist_ok=True)
        cmdArgs = ['rsync', '--update', filePath, destFilePath]
        process = subprocess.Popen(cmdArgs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        logging.info('dataset info updated')
    except: logging.warning('failed in updating dataset info')


def loadDatasets_info(databasePath, databaseBakPaths=None, copy_if_not_exist=True):
    database = dbt.Database([databasePath])
    database.load_additional_datasets_info()
    if database.additional_datasets_info:
        datasets_info = database.additional_datasets_info
    elif databaseBakPaths:
        database = dbt.Database(databaseBakPaths)
        if copy_if_not_exist:
            for path in databaseBakPaths:
                if path in database.additional_datasets_info_paths[0]:
                    databasePathBak = path
            destFilePath = workDataDirCopy(database.additional_datasets_info_paths[0], databasePath, databasePathBak)
            with open(destFilePath, 'r') as f:
                datasets_info = json.load(f)
        else:
            with open(database.additional_datasets_info_paths[0], 'r') as f:
                datasets_info = json.load(f)
    return datasets_info

def workDataDirCopy(filePath, workDataDir, workDataDirBak):
    relpath = os.path.relpath(filePath, workDataDirBak)
    destFilePath = os.path.join(workDataDir, relpath)
    logging.warning('file not found in wordDataDir {}, but found in workDataDirBak: {}'.format(workDataDir, workDataDirBak))
    logging.warning('now copying to {} ...'.format(destFilePath))
    os.makedirs(os.path.dirname(destFilePath), exist_ok=True)
    cmdArgs = ['cp', filePath, destFilePath]
    process = subprocess.Popen(cmdArgs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    logging.info('{} copied'.format(destFilePath))
    logging.warning('file copied')
    return destFilePath
