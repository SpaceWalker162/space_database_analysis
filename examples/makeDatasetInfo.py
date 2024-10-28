import os
import space_database_analysis.databaseTools as dbt
import logging

if __name__ == '__main__':
    loggingHandlers = []
    loggingHandlers.append(logging.StreamHandler())
    formatter = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s",
                                  "%Y-%m-%d %H:%M:%S")
    loggingHandlers[0].setFormatter(formatter)
    logging.basicConfig(level=logging.WARNING, handlers=loggingHandlers)
    logging.getLogger().setLevel(logging.INFO)

    databasePath = os.path.expanduser('~/Documents/remoteDatabase')
    database = dbt.Database([databasePath])

    ## to initialize, uncomment the following lines and run
    database.make_additional_datasets_info()
    database.add_example_time_interval_to_add_dataset_info()
    
    ## to add user defined information to datasets in the database use the following codes
    add_dic = {}
#    for spacecraftInd in range(1, 5):
#        spacecraftName = 'C' + str(spacecraftInd)
#
#        datasetID = spacecraftName + '_CP_FGM_FULL'
#        dataset = {'Id': datasetID}
#        dataset_path = os.path.join('cluster', spacecraftName.lower(), 'fgm', 'mfield_3dvect_fullreso')
#        dataset.update({'dataset_path': dataset_path})
#        add_dic[datasetID] = dataset
#
#        datasetID = spacecraftName + '_CP_FGM_5VPS'
#        dataset = {'Id': datasetID}
#        dataset_path = os.path.join('cluster', spacecraftName.lower(), 'fgm', 'mfield_3dvect_5vectpersec')
#        dataset.update({'dataset_path': dataset_path})
#        add_dic[datasetID] = dataset

    for datasetID in ['WI_H0_MFI', 'WI_H0_SWE', 'WI_H1_SWE', 'WI_K0_SWE']:
        dataset = {'Id': datasetID}
        datasetID_com = datasetID.split('_')
        inst = datasetID_com[-1]
        dataset_path = os.path.join('wind', inst.lower(), '_'.join([inst.lower(), datasetID_com[1].lower()]))
        dataset.update({'dataset_path': dataset_path})
        add_dic[datasetID] = dataset

    database.load_additional_datasets_info()
    database.update_additional_datasets_info(add_dic)
    database.save_additinal_datasets_info()
