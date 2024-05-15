# this file present example usage of the module cdasws
import cdasws
import os
from datetime import datetime
import numpy as np
import space_database_analysis.dataAnalysisTools as dat

cdaswsObj = cdasws.CdasWs()
observatories = cdaswsObj.get_observatories()
for observatory in observatories:
    if 'mms' in observatory['Name'].lower():
        print(observatory['Name'])
cdaswsObj.get_instrument_types()
cdaswsObj.get_instrument_types(observatory='MMS1')
cdaswsObj.get_instruments(observatory='MMS1')
datasets = cdaswsObj.get_datasets(observatory='MMS1', instrument='SCM')
datasets
dataset = datasets[-1]
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

timeInterval = cdaswsObj.get_example_time_interval(dataset['Id'])
variable_names = cdaswsObj.get_variable_names(dataset['Id'])
variables = cdaswsObj.get_variables(dataset['Id'])
eventTime = dat.Epochs(dateTimeList=[datetime(2016, 1, 8, 10)], epochType='CDF_TIME_TT2000').epochs
epochObj = dat.Epochs(epochs=eventTime + np.array([-1, 1])*10**9 * 60*60*24, epochType='CDF_TIME_TT2000')
epochObj.dateTimeList
status, data = cdaswsObj.get_data(dataset['Id'], variable_names, *epochObj.dateTimeList)
status
status, files = cdaswsObj.get_original_files(datasetId, [beginOfTheFilePeriod, endOfTheFilePeriod])
len(files)
files[0]
files
print('done')
len(data)
data.keys()
data['mms1_edp_dce_gse_fast_l2']
data[0]
