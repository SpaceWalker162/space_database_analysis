# this file contains obsolete functions of databaseUserTools.py

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
    def __init__(self, mission=None, spacecraft=None, instrumentation=None, instrumentationPath=None, pathToDataset=None, splitedPathToDataset=None, datasetID=None, dateTime=None, timeRange=None, filePeriodRange=None, dataset_file_time_gap=None, fileName=None, instrumentationObj=None, cdfFile=None, workDataDir=None, fileSize='allSize', epochType='CDF_EPOCH', silence=True, filePath=None, filePaths=[], defineCDFFileQ=True):
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
        absolutePathToDataset = os.path.join(self.workDataDir, self.pathToDataset)
        self.define_search_criteria(size=size)
        searchMethod = self.search_criteria['searchMethod']
        logging.debug("Looking for data files in: {}".format(absolutePathToDataset))
        logging.debug("search method: {}".format(searchMethod))
        if searchMethod == 'general':
            fileNames, absolutePathToFile = self.findFileNamesUnderDataset(absolutePathToDataset, func=findFileNames, **self.search_criteria)
#            = findFileNames(absolutePathToFile)
        elif searchMethod == 'allFilesInTimeRange':
            getTimeFromNameFunc_ = getTimeFromName
            logging.debug('searchmethod:')
            logging.debug(self.timeRange)
            fileNames, absolutePathToFile = self.findFileNamesUnderDataset(absolutePathToDataset, func=findFileNamesInTimeRange, timeRange=self.timeRange, getTimeFromNameFunc=getTimeFromNameFunc_, **self.search_criteria)
        elif searchMethod == 'ClusterCAA':
            fileNames, absolutePathToFile = self.findFileNamesUnderDataset(absolutePathToDataset, func=findFileNamesInInterval, interval=self.search_criteria['interval'], strings=self.instrumentation)
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
                logging.warning("More than one files were found in {}:" + ("\n{}"*numFiles).format(absolutePathToFile, *fileNames))

#            raise Exception('More than one files were found.')

    def define_search_criteria(self, size='allSize'):
        self.search_criteria = {}
        if self.dateTime:
            start = self.dateTime
            end = start + timedelta(seconds=1)
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
                self.search_criteria.update({'interval': interval})
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

        self.search_criteria.update({'searchMethod': searchMethod, 'strings': criteria, 'timeTag': timeTag, 'size': size})

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


def readData(workDataDir, datasetsAndVariables, datetimeRange, epochType='CDF_EPOCH', workDataDirsBak=None, saveDataFileInWorkDataDir=True):
    '''
    Purpose:
        This function is to load data.
    Parameters:
        workDataDir: the path to the directory of working data base
        datasetAndVariables: a dictionary in terms of the path to the datasets. Its leaf value is a list of variable names. The names are defined by the corresponding cdfFile. For example, {'Cluster': {'C1' : {'C1_CP_FGM_FULL': ['time_tags__C1_CP_FGM_FULL', 'B_vec_xyz_gse__C1_CP_FGM_FULL']}}, 'mms': {'mms1': {'fgm': {'brst': {'l2': ['Epoch', 'Btotal']}}}}}. Please note that Epoch and Btotal are improvised. It may also be a list of lists, with each sublist in the form of ['Cluster', 'C1', 'C1_CP_FGM_FULL', ['time_tags__C1_CP_FGM_FULL', 'B_vec_xyz_gse__C1_CP_FGM_FULL']].
        datetimeRange: a list of two elements which define the time interval during which the data is to be retrieved. [start, end]
        workDataDirsBak: a list of backup databases
        saveDataFileInWorkDataDir: self-explanatory
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
        logging.debug('finding files for dataset and variables: ')
        logging.debug(datasetAndVariables)
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
                        if saveDataFileInWorkDataDir:
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
            dataMajor, dataAux = readDataFromACdfFile(dataFile.cdfFile, variableNames, datetimeRange)
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

    def loadData(self, datetimeRange=None, instrumentation=None, instrumentationRetrivingName=None, datasets_variables_with_retrieving_names=None, variables=None, variableRetrivingNames=None, datasetAndVariables=None, variablesWithRetrivingNames=None, instrumentationVariablesWithRetrivingName=None, cleanData=False, tStepPecentageCriterion=0.9, lowpassCutoff=None, inPlace=False, gapThreshold=None, minNumberOfPoints=None, returnShiftQ=False, copy_if_not_exist=True, search_online=False, fromFile=None, useMask=False, maskValue=-1.0*10**30, sparse_factor=None):
        '''
        Parameters:
            datetimeRange: a list of two datetime objects, representing the start and end of the data to be loaded.
            instrumentation: a list in terms of the directory names at all levels below spacecraft and above year or files. For example, ['fgm', 'brst', 'l2']
            instrumentationRetrivingName: a string.
            datasetAndVariables: this parameter should not be used by user.
            instrumentationVariablesWithRetrivingName: a dictionary in terms of the path to the datasets. Its leaf value is a list of variable names. The names are defined by the corresponding cdfFile. For example, {'Cluster': {'C1' : {'C1_CP_FGM_FULL': ['FGM', ('time_tags__C1_CP_FGM_FULL', 't'), ('B_vec_xyz_gse__C1_CP_FGM_FULL', 'B')]}}, 'mms': {'mms1': {'fgm': {'brst': {'l2': ['fgmBrst', ('Epoch', 't'), ('Btotal', 'BTotal')]}}}}}. Please note that 'Epoch' and 'Btotal' are improvised. This parameter may also be a list of lists, with each sublist in the form of ['Cluster', 'C1', 'C1_CP_FGM_FULL', ['FGM', ('time_tags__C1_CP_FGM_FULL', 't'), ('B_vec_xyz_gse__C1_CP_FGM_FULL', 'B')]]. To retrieve data, for example, of 'B_vec_xyz_gse__C1_CP_FGM_FULL', use C1.data['FGM']['B']
            datasets_variables_with_retrieving_names: a dictionary of datasets. Its leaf value is a list of variable names. The names are defined by the corresponding cdfFile. For example, {'MMS1_FGM_SRVY_L2': ['FGM', ('Epoch', 't'), ('mms1_fgm_srvy_l2_bvec_gse', 'B')], 'MMS2_FPI_FAST_L2_DIS-MOMS': ['fgmBrst', ('Epoch', 't'), ('mms2_fpi_fast_l2_dis-moms_density', 'n')]}. Please note that 'Epoch' and 'Btotal' are improvised. To retrieve data 'B_vec_xyz_gse__C1_CP_FGM_FULL', use C1.data['FGM']['B']
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
            dataFromACdfFile = readDataFromACdfFile(cdfFile, variableNames)
            self.data.update({instrumentationRetrivingName: dataFromACdfFile})
        else:
            assert datasets_variables_with_retrieving_names
            datasets_info = loadDatasets_info(self.workDataDir, self.workDataDirsBak, copy_if_not_exist=copy_if_not_exist)
            for datasetID, item in datasets_variables_with_retrieving_names.items():
                datasetRetrievingName = item[0]
                variableNamesAndRetrievingNames = item[1:]
                variableNames = [varRet[0] for varRet in variableNamesAndRetrievingNames]
                dataset_info = datasets_info.get(datasetID)
                if dataset_info:
                    dataset = Dataset(dataset_info=dataset_info, databasePath=self.workDataDir, databaseBakPaths=self.workDataDirsBak)
                else:
                    dataset = Dataset(datasetID=datasetID, databasePath=self.workDataDir, databaseBakPaths=self.workDataDirsBak)
                dataset.load_data(variableNames=variableNames, datetimeRange=datetimeRange, copy_if_not_exist=copy_if_not_exist, search_online=search_online, sparse_factor=sparse_factor)
                if self.epoch_type_to_load_data:
                    for varName in variableNames:
                        data_type = dataset.varInfoDict[varName]['varInfo'].Data_Type_Description
                        if data_type in _epoch_types:
                            if not data_type == self.epoch_type_to_load_data:
                                para_ = {data_type: dataset.data[varName]}
                                dataset.data[varName] = dat.Epochs(**para_).get_data(fm=self.epoch_type_to_load_data)
                for varName, retName in variableNamesAndRetrievingNames:
                    dataset.data[retName] = dataset.data.get(varName)
                    dataset.varInfoDict[retName] = dataset.varInfoDict.get(varName)
                self.data.update({datasetRetrievingName: dataset.data})
                self.metadata.update({datasetRetrievingName: dataset.varInfoDict})
            if cleanData:
                self.cleanData(instrumentation=instrumentationRetrivingName, tStepPecentageCriterion=tStepPecentageCriterion, lowpassCutoff=lowpassCutoff, inPlace=inPlace, gapThreshold=gapThreshold, minNumberOfPoints=minNumberOfPoints, returnShiftQ=returnShiftQ)

    def readData(self, workDataDir, datasetsAndVariables, datetimeRange, workDataDirsBak=None, copy_if_not_exist=True, search_online=False):
        '''
        Purpose:
            This function is to load data.
        Parameters:
            workDataDir: the path to the directory of working data base
            datasetAndVariables: a dictionary in terms of the path to the datasets. Its leaf value is a list of variable names. The names are defined by the corresponding cdfFile. For example, {'Cluster': {'C1' : {'C1_CP_FGM_FULL': ['time_tags__C1_CP_FGM_FULL', 'B_vec_xyz_gse__C1_CP_FGM_FULL']}}, 'mms': {'mms1': {'fgm': {'brst': {'l2': ['Epoch', 'Btotal']}}}}}. Please note that Epoch and Btotal are improvised. It may also be a list of lists, with each sublist in the form of ['Cluster', 'C1', 'C1_CP_FGM_FULL', ['time_tags__C1_CP_FGM_FULL', 'B_vec_xyz_gse__C1_CP_FGM_FULL']].
            datetimeRange: a list of two elements which define the time interval during which the data is to be retrieved. [start, end]
            workDataDirsBak: a list of backup databases
            saveDataFileInWorkDataDir: self-explanatory
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

        # load dataset info
        datasets_info = loadDatasets_info(self.workDataDir, self.workDataDirsBak)
        for datasetAndVariables in datasetsAndVariablesList:
            logging.debug('finding files for dataset and variables: ')
            logging.debug(datasetAndVariables)
            splitedPathToDataset = datasetAndVariables[:-1]
            variableNames = datasetAndVariables[-1]
            start_ = start
            variablesInADatasetOverDatetimeRange = [] # first index for data from different files over a epoch range, second index for variables
            datasetID = '_'.join(splitedPathToDataset).upper()
            dataset_info = datasets_info[datasetID]
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
                                    if copy_if_not_exist:
                                        destFilePath = workDataDirCopy(filePathBak, workDataDir, workDataDirBak)
                                        desiredDataPaths.append(destFilePath)
                                    else:
                                        desiredDataPaths.append(filePathBak)
                desiredDataPaths.sort()
                dataFile = DataFile(filePaths=desiredDataPaths)
                for fileInd, cdfFile in enumerate(dataFile.cdfFiles):
                    logging.debug('reading data file: {}'.format(dataFile.filePaths[fileInd]))
                    logging.debug('datetimeRange: {}'.format(datetimeRange))
                    dataMajor, dataAux = readDataFromACdfFile(cdfFile, variableNames, datetimeRange)
                    logging.debug('reading data file done: {}'.format(dataFile.filePath))
                    variablesInADatasetOverDatetimeRange.append(dataMajor)
                    variablesInADatasetIndependantOnTime.append(dataAux)
            else:
                logging.debug('not in brst')
                dataset = Dataset(dataset_info=dataset_info, databasePath=self.workDataDir, databaseBakPaths=self.workDataDirsBak)
                dataset.load_data(variableNames=None, datetimeRange=None, copy_if_not_exist=copy_if_not_exist, search_online=search_online)

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
