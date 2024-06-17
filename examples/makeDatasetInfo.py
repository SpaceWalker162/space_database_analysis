import os
import space_database_analysis.databaseTools as dbt
import space_database_analysis.databaseUserTools as dut

if __name__ == '__main__':
    databasePath = os.path.expanduser('~/Documents/remoteDatabase')
    database = dbt.Database([databasePath])
    database.make_additional_datasets_info()
    ##
    databasePath_local = os.path.expanduser('~/Documents/database')
    databasePath_remote = os.path.expanduser('~/Documents/remoteDatabase')
