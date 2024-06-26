import os
import copy
import sys
import space_database_analysis.otherTools as ot
import space_database_analysis.databaseTools.downloadSPDF as downloadSPDF
import socket

if sys.platform == 'linux':
    if socket.gethostname() == 'yufei-OptiPlex-5040-2':
        downloadDatabaseDir = '/home/yufei/Documents/database' #This is where you want to store the downloaded files
    else:
        downloadDatabaseDir = '/media/yufei/Elements/database' #This is where you want to store the downloaded files
    databaseDirs = ['/home/yufei/Documents/remoteDatabase'] #This is a list that contain all database. The files that exist in these databases will be omitted in downloading.
elif sys.platform == 'win32':
    downloadDatabaseDir = '..\\database'
    databaseDirs = ['\\\\10.249.183.237\\pub']

#databaseDirs = []

readInternetForFileNames = True # if True, the program reads from ftp the list of the files you would like to present in your downloadDatabaseDir and databaseDirs after it ends. This suit the first run of the program. After a run with this parameter set True, the list of files will be stored locally for later use, such as a second run when the first run is not successful. In this case, this parameter should be set to False and so the program will read the locally stored list to save the time spent on reading ftp.

#dataNameDict = {'ace': {'mag': {'level_2_cdaweb': 'mfi_h3'}, 'swepam': {'level_2_cdaweb': 'swe_h0'}}}

#missionName = 'mms'
#spacecraftNames = []
#mmsNumbers = [1]
#for mmsNumber in mmsNumbers:
#    spacecraftNames.append(missionName+str(mmsNumber))
#for i in range(4):
#    spacecraftNames.append(missionName+str(i+1))
#instrumentations = [['fpi', 'fast', 'l2', 'dis-moms'], ['fpi', 'brst', 'l2', 'dis-moms']]
#instrumentations = [['mec', 'srvy', 'l2', 'epht89q'], ['mec', 'srvy', 'l2', 'epht89d']]
#instrumentations = [['fgm', 'srvy', 'l2']]
#instrumentations = [['edp', 'slow', 'l2', 'dce', '2019'], ['edp', 'fast', 'l2', 'dce', '2019']]
#instrumentations = [['fpi', 'fast', 'l2', 'dis-partmoms', '2015'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2016'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2017'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2018']]
#instrumentations = [['fpi', 'fast', 'l2', 'dis-partmoms', '2015'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2016'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2017'], ['fpi', 'fast', 'l2', 'dis-partmoms', '2018']]
#dataNameDict = {missionName: {}}
#directory = ot.Directory(lis=instrumentations)
#directory.generate_dic_from_lis()
#instrumentationsDict = directory.dic
#instrumentationsDict = {'feeps': {'srvy': {'l2': {'ion': {'2015':{}, '2016':{}, '2017': {}, '2018': {}}}}}}
#for spacecraftName in spacecraftNames:
#    dataNameDict[missionName][spacecraftName] = copy.deepcopy(instrumentationsDict)
#
dataNameDict = {'wind': {'3dp': {'3dp_plsp': {}}}}

#missionName = 'themis'
#spacecraftNames = ['thc', 'thd']
#instrumentations = [['l2', 'fgm']]
#dataNameDict = {missionName: {}}
#instrumentationsDict = ot.list2dict(instrumentations)
#for spacecraftName in spacecraftNames:
#    dataNameDict[missionName][spacecraftName] = copy.deepcopy(instrumentationsDict)

logFileDir = os.path.expanduser('~/Documents/MyFiles/works/project_working/temp/downloadSPDF')
protocol = 'http'
if readInternetForFileNames:
    fileNamesSource = 'internet'
else:
    fileNamesSource = 'from_log'
downloadSPDF.downloadSPDF(downloadDatabaseDir, databaseDirs, dataNameDict, logFileDir=logFileDir, fileNamesSource=fileNamesSource, protocol=protocol)
print('end')
