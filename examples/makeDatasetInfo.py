import os
import space_database_analysis.databaseTools as dbt
import space_database_analysis.databaseUserTools as dut

if __name__ == '__main__':
    databasePath = os.path.expanduser('~/Documents/remoteDatabase')
    database = dbt.Database([databasePath])
    database.make_additional_datasets_info()
    print('end')
    ##
    dataset_dic = dut.loadDatasets_info(databasePath)
    dataset_dic.keys()
    for key, item in dataset_dic.items():
        if 'mms' in key.lower():
            ii = item

    ii
