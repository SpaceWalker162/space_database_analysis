import os
import copy
import space_database_analysis.otherTools as ot
from space_database_analysis.databaseTools.downloadSPDF import downloadSPDF

downloadDataDir = '/media/yufei/Elements/data' #This is where you want to store the downloaded files
#downloadDataDir = '..\\data'

databaseDirs = ['/home/yufei/Documents/remoteDatabase'] #This is a list that contain all database. The files that exist in these databases will be omitted in downloading.
#databaseDirs = ['\\\\10.249.183.237\\data']
#databaseDirs = []

readFTP = True # if True, the program reads from ftp the list of the files you would like to present in your downloadDataDir and databaseDirs after it ends. This suit the first run of the program. After a run with this parameter set True, the list of files will be stored locally for later use, such as a second run when the first run is not successful. In this case, this parameter should be set to False and so the program will read the locally stored list to save the time spent on reading ftp.

#dataNameDict = {'ace': {'mag': {'level_2_cdaweb': 'mfi_h3'}, 'swepam': {'level_2_cdaweb': 'swe_h0'}}}

missionName = 'mms'
spacecraftNames = []
mmsNumbers = [1,2]
for mmsNumber in mmsNumbers:
    spacecraftNames.append(missionName+str(mmsNumber))
#for i in range(4):
#    spacecraftNames.append(missionName+str(i+1))
#instrumentations = [['fpi', 'fast', 'l2', 'dis-moms'], ['fpi', 'brst', 'l2', 'dis-moms']]
#instrumentations = [['mec', 'srvy', 'l2', 'epht89q'], ['mec', 'srvy', 'l2', 'epht89d']]
#instrumentations = [['fgm', 'srvy', 'l2']]
instrumentations = [['edp', 'slow', 'l2', 'dce', '2019'], ['edp', 'fast', 'l2', 'dce', '2019']]
#instrumentations = [['fpi', 'fast', 'l2', 'dis-moms', '2016']]
dataNameDict = {missionName: {}}
instrumentationsDict = ot.list2dict(instrumentations)
for spacecraftName in spacecraftNames:
    dataNameDict[missionName][spacecraftName] = copy.deepcopy(instrumentationsDict)
#

#missionName = 'themis'
#spacecraftNames = ['thc', 'thd']
#instrumentations = [['l2', 'fgm']]
#dataNameDict = {missionName: {}}
#instrumentationsDict = ot.list2dict(instrumentations)
#for spacecraftName in spacecraftNames:
#    dataNameDict[missionName][spacecraftName] = copy.deepcopy(instrumentationsDict)

logFileDir = os.path.expanduser('~/Documents/MyFiles/works/project_working/temp/downloadSPDF')
downloadSPDF(downloadDataDir, databaseDirs, dataNameDict, logFileDir=logFileDir, fileNamesSource='from_log')
