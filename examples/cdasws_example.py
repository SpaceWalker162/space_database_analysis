# this file present example usage of the module cdasws
if __name__ == '__main__':
    import cdasws
    import os
    from datetime import datetime
    import numpy as np
    import space_database_analysis.dataAnalysisTools as dat

    cdaswsObj = cdasws.CdasWs()
    observatories = cdaswsObj.get_observatories()
    for observatory in observatories:
        if 'ac' in observatory['Name'].lower():
            print(observatory['Name'])
    cdaswsObj.get_instrument_types()
    cdaswsObj.get_instrument_types(observatory='AC')
    instruments = cdaswsObj.get_instruments(observatory='AC')
    instruments
    datasets = cdaswsObj.get_datasets(observatory='AC', instrument='MAG')
    datasets
    for dataset in datasets:
        print(dataset['Id'])
    dataset = datasets[-1]
    dataset = datasets[4]
    dataset
    from urllib.parse import urlparse
    from urllib.request import urlretrieve
    for d_ in dataset['AdditionalMetadata']:
        if d_['Type'] == 'SPDF_SKT_CDF':
            cdf_file_info_url = d_['value']
    o = urlparse(cdf_file_info_url)
    file_name = os.path.split(o.path)[-1]
    save_path = os.path.expanduser('~/Downloads')
    full_file_path = os.path.join(save_path, file_name)
    urlretrieve(cdf_file_info_url, full_file_path)

    datasetID = 'WI_H0_MFI'
    timeInterval = cdaswsObj.get_example_time_interval(datasetID)
    variable_names = cdaswsObj.get_variable_names(datasetID)
    timeInterval.basic_iso_format()
    timeInterval._start
    timeInterval._end
    timeInterval = cdaswsObj.get_example_time_interval(dataset['Id'])
    variable_names = cdaswsObj.get_variable_names(dataset['Id'])
    variables = cdaswsObj.get_variables(dataset['Id'])
    eventTime = dat.Epochs(dateTimeList=[datetime(2016, 1, 8, 10)], epochType='CDF_TIME_TT2000').epochs
    epochObj = dat.Epochs(epochs=eventTime + np.array([-1, 1])*10**9 * 60*60*24, epochType='CDF_TIME_TT2000')
    epochObj.dateTimeList
    status, data = cdaswsObj.get_data(dataset['Id'], variable_names, *epochObj.dateTimeList)
    status

    datetimeRange = [datetime(2015, 10, 7), datetime(2015, 10, 7, 8)]
    datetimeRange = [datetime(2002, 2, 7), datetime(2002, 2, 9)]
    status, files = cdaswsObj.get_original_files(dataset['Id'], *datetimeRange)
    status
    len(files)
    files[0]
    files
    print('done')
    len(data)
    data.keys()
    data['mms1_edp_dce_gse_fast_l2']
    data[0]
    ##
    for observatory in observatories:
        if 'omni' in observatory['Name'].lower():
            print(observatory)
    instruments = cdaswsObj.get_instrument_types(observatory='OMNI (1AU IP Data)')
    print(instruments)
    datasets = cdaswsObj.get_datasets(observatory='OMNI (1AU IP Data)', instrument=instruments[2]['Name'])
    datasets = cdaswsObj.get_datasets(observatory='OMNI', instrument=instruments[3]['Name'])
    datasets = cdaswsObj.get_datasets(observatory='OMNI (1AU IP Data)')
    datasets = cdaswsObj.get_datasets()
    datasets_all = datasets
    datasets_all[:2]
    type(datasets)
    datasets
    dataset = datasets[2]
    for dataset in datasets:
        print(dataset['Id'])
    status, files = cdaswsObj.get_original_files(dataset['Id'], *datetimeRange)
    ##
    database = dbt.Database([remoteWorkDataDir])
    datasets_info = database.loadDatasets_info()

