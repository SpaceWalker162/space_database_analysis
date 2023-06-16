
#https://pds-ppi.igpp.ucla.edu/search/view/?f=yes&id=pds://PPI/CO-E_J_S_SW-CAPS-3-CALIBRATED-V1.0/DATA/CALIBRATED
#
#https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/CO-E_J_S_SW-CAPS-3-CALIBRATED-V1.0/DATA/CALIBRATED/2004/001

import databaseTools as dbt
from datetime import datetime
from datetime import timedelta
import json
import time
import os
import numpy as np
import logging

verbose = True
loggingHandlers = [logging.FileHandler('download.log')]
if verbose:
    loggingHandlers.append(logging.StreamHandler())
logging.basicConfig(level=logging.DEBUG, handlers=loggingHandlers)

logDir = os.path.expanduser('~/Documents/MyFiles/works/project_working/.logFiles/')

databaseDir = ('/media/yufei/Elements1/data/')

urlBase = "https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/CO-E_J_S_SW-CAPS-3-CALIBRATED-V1.0/DATA/CALIBRATED/"
#
instrumentation = ['cassini', 'MAG', 'CO-E_SW_J_S-MAG-3-RDR-FULL-RES-V2.0', 'DATA']
urlBase = "https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/{ins:}/DATA/{year:}/{days:}"

urlBase = "https://pds-ppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/{ins:}/DATA/2004/001_031_JAN"


with open("../.logFiles/lost.log", 'r') as f:
    info = f.read().strip()
fileNames = info.split(",")
year = 2004
for fileName in fileNames:
    basename, ext = os.path.splitext(fileName)
    print(basename)
    print(int(basename)-1)
    currentDate = datetime(year, 1, 1) + timedelta(days=(int(basename)-1))
    yearStr = currentDate.strftime('%Y')
    daysInYear = currentDate.strftime('%j')
    path = os.path.join(databaseDir, 'cassini', 'CAPS', 'CO-E_J_S_SW-CAPS-3-CALIBRATED-V1.0', yearStr)
    if not os.path.exists(path):
        os.makedirs(path)
    downloadedFileName = os.path.join(path, daysInYear + '.zip')
    url = urlBase + currentDate.strftime('%Y/%j')
    logging.info("downloading: {}".format(url))
    fileSizeInM, timeCost, speed = dbt.download(url, fileName=downloadedFileName)
    logging.info("year: {}, date: {}\n size: {}M cost: {}, speed: {}"
          .format(year, daysInYear, fileSizeInM, timeCost, speed))
    currentDate = currentDate + timedelta(days=1)

#start = datetime(2004, 1, 1)
#end = datetime(2012, 1, 1) + timedelta(days=153)
##interval = end - start
##numberOfDays = interval.days
#currentDate = start
#while currentDate != end + timedelta(days=1):
##for days in range(numberOfDays):
##    currentDate = start + timedelta(days=days)
#    year = currentDate.strftime('%Y')
#    daysInYear = currentDate.strftime('%j')
#    path = os.path.join(databaseDir, 'Cassini', 'CO-E_J_S_SW-CAPS-3-CALIBRATED-V1.0', year)
#    if not os.path.exists(path):
#        os.makedirs(path)
#    downloadedFileName = os.path.join(path, daysInYear + '.zip')
#    url = urlBase + currentDate.strftime('%Y/%j')
#    logging.info("downloading: {}".format(url))
#    fileSizeInM, timeCost, speed = dbt.download(url, fileName=downloadedFileName)
#    logging.info("year: {}, date: {}\n size: {}M cost: {}, speed: {}"
#          .format(year, daysInYear, fileSizeInM, timeCost, speed))
#    currentDate = currentDate + timedelta(days=1)
