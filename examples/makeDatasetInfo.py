import os
import space_database_analysis.databaseTools as dbt

if __name__ == '__main__':
    databasePath = os.path.expanduser('~/Documents/remoteDatabase')
    database = dbt.Database([databasePath])

    ## to initialize, uncomment the following line and run
    database.make_additional_datasets_info()
    
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

    for name in ['MOM', 'SPEC']:
        datasetID = 'MVN_SWI_L2_ONBOARDSVY' + name
        dataset = {'Id': datasetID}
        dataset_path = os.path.join('maven', 'swia', 'l2', 'onboardsvy'+name.lower())
        dataset.update({'dataset_path': dataset_path})
        add_dic[datasetID] = dataset

    datasetID = 'MVN_MAG_L2_SUNSTATE-1SEC'
    dataset = {'Id': datasetID}
    dataset_path = os.path.join('maven', 'mag', 'l2', 'sunstate-1sec', 'cdfs')
    dataset.update({'dataset_path': dataset_path})
    add_dic[datasetID] = dataset

    database.load_additional_datasets_info()
    database.update_additional_datasets_info(add_dic)
#    database.save_additinal_datasets_info()
