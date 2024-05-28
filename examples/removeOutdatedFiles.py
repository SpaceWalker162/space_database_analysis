import os
import logging
import space_database_analysis.databaseTools as dbt

remoteWorkDataDir = os.path.expanduser('~/Documents/remoteDatabase/data/mms/mms1/edp/slow/l2/dce/2015')
databasePaths = [remoteWorkDataDir]

logPath = os.path.expanduser("~/Documents/MyFiles/works/research_working/.logFiles")
scriptFileName = "removeOutdatedFiles"
logFileName = os.path.join(logPath, scriptFileName +'.log')
verbose = True
loggingHandlers = [logging.FileHandler(logFileName)]
formatter = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s",
                              "%Y-%m-%d %H:%M:%S")
loggingHandlers[0].setFormatter(formatter)
if verbose:
    loggingHandlers.append(logging.StreamHandler())
#logging.basicConfig(level=logging.DEBUG, handlers=loggingHandlers)
logging.basicConfig(level=logging.INFO, handlers=loggingHandlers)
#logging.getLogger().setLevel(logging.INFO)

dbt.Database.removeOutdatedFiles(databasePaths)
print('end')
