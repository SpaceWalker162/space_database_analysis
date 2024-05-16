def transformPDSdataToCDF(databaseDir, dataTransDict=None, stringCriteria=['FGM_KSM_1M'], infoFileExtension='.LBL', cdfFileNameFMT=None, recordPath=None):
    dataTransDict = copy.deepcopy(dataTransDict)
    transformedFiles = []
    if not os.path.exists(recordPath):
        with open(recordPath, 'w') as f:
            pass
    else:
        with open(recordPath, 'r') as f:
            info_ = f.readlines()
        for line in info_:
            transformedFiles.append(line.strip())
    try:
        years = os.listdir(databaseDir)
        for year in years:
            try: yearNum = int(year)
            except: yearNum = None
            if yearNum is None:
                continue
            yearDataDir = os.path.join(databaseDir, year)
            logging.info(yearDataDir)
            fileNames = os.listdir(yearDataDir)
            for fileName in fileNames:
                filePath = os.path.join(yearDataDir, fileName)
                logging.info('filePath: {}'.format(filePath))
                stringCriteria.append(infoFileExtension)
                cdfDir = yearDataDir
                if all([string in fileName for string in stringCriteria]):
                    if filePath in transformedFiles:
                        continue
                    filebaseName, ext = os.path.splitext(filePath)
                    dataDict, _ = dut.readPDSData(fileName=filebaseName, dataFileExtension='.TAB', infoFileExtension=infoFileExtension, sep=None)
                    numberOfOutVariables = len(dataTransDict)
                    for vInd in range(numberOfOutVariables):
                        inDataNames = dataTransDict[str(vInd)]['inName']
                        if isinstance(inDataNames, str):
                            dataDict_ = dataDict[inDataNames]
                            data_ = dataDict_['data']
                            dataType = dataDict_['data_type']
                        elif isinstance(inDataNames, list):
                            data_List = [dataDict[inName]['data'][:, None] for inName in inDataNames]
                            data_ = np.concatenate(data_List, axis=-1)
                            dataType = dataDict[inDataNames[0]]['data_type']
                        else:
                            raise Exception('unknown input data')
                        dataTransDict[str(vInd)]['data'] = data_
                        dataTransDict[str(vInd)]['in_data_type'] = dataType
                        dataTransDict[str(vInd)]['out_data_type'] = dataTypeTransformationDict_PDS[dataType]

                    currentDate = datetime(int(year), 1, 1)
                    cdfFileName = currentDate.strftime(cdfFileNameFMT)
                    cdfFilePath = os.path.join(cdfDir, cdfFileName)
                    if os.path.exists(cdfFilePath):
                        os.remove(cdfFilePath)
                    cdfFile = cdflib.cdfwrite.CDF(cdfFilePath)
                    globalAttrs = {}
                    globalAttrs['file author'] = {0: 'Yufei Zhou'}
                    globalAttrs['data format source file'] = {0: fileName}
                    cdfFile.write_globalattrs(globalAttrs)
                    for vInd in range(len(dataTransDict)):
                        dataTransDict_ = dataTransDict[str(vInd)]
                        varData = dataTransDict_['data']
                        varDType = dataTransDict_['out_data_type']
                        varDType = cdflib.cdfwrite.CDF._datatype_token(varDType)
                        var_spec = {'Variable': dataTransDict_['outName'],
                                'Data_Type': varDType,
                                'Num_Elements': 1,
                                'Rec_Vary': True,
                                'Dim_Sizes': varData.shape[1:],
                                    }
                        var_attrs = {}
                        depend0 = dataTransDict_.get('DEPEND_0', None)
                        if depend0:
                            var_attrs['DEPEND_0'] = depend0
                        cdfFile.write_var(var_spec, var_attrs=var_attrs, var_data=varData)
                    cdfFile.close()
                    print(cdfFilePath)
                    transformedFiles.append(filePath)
    except Exception as e:
        with open(recordPath, 'w') as f:
            for transformedFile in transformedFiles:
                print(transformedFile, file=f)
        raise e
