import numpy as np
#import constants as con
import cdflib    # see github.com/MAVENSDC/cdflib
#import tarfile
from datetime import datetime
from datetime import timedelta
import os
import sys
#import subprocess
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker as ticker
from matplotlib import rcParams
from itertools import combinations
import json
from pprint import pprint
import tensorflow.nn as nn
from scipy.signal import butter, lfilter, freqz
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
import scipy.interpolate as interpolate
import scipy.special 
from databaseTools import *
from wolframclient.language import wl
from wolframclient.evaluation import WolframLanguageSession


def extractFiles(databaseDir, workDataDir, spacecrafts, datasets, interval):
    '''
    parameters:
        databaseDir: implied by its name
        workDataDir: implied by its name
        spacecrafts: a list of spacecraft names
        datasets: a list of dataset names corresponding to the list of spacecraft names
    '''
    start_ = interval[0]
    start_ = datetime(start_.year, start_.month, 1)
    end_ = interval[1]
    if end_ == datetime(end_.year, end_.month, 1):
        pass
    elif end_.month == 12:
        end_ = datetime(end_.year+1, 1, 1)
    else:
        end_ = datetime(end_.year, end_.month+1, 1)
    interval = 1
    start = start_
    if start.month == 12:
        end = datetime(start.year+1, 1, 1)
    else:
        end = datetime(start.year, start.month+1, 1)
    while start.year != end_.year or start.month != end_.month:
        for spacecraft, dataset in zip(spacecrafts, datasets):
            pathDatabase = os.path.join(databaseDir, spacecraft,
                                        dataset, 'compressed')
            pathWorkData = os.path.join(workDataDir, spacecraft, dataset, start.strftime("%Y"), start.strftime("%m"))
            if not os.path.exists(pathWorkData):
                os.makedirs(pathWorkData)
    #        else:
    #            files = os.listdir(pathWorkData)
    #            for file in files:
    #                filePath = os.path.join(pathWorkData, file)
    #                if os.path.isfile(filePath):
    #                    os.unlink(filePath)
            fileName = '--'.join([dataset, start.strftime("%Y--%m"),
                                 end.strftime("%Y--%m"), 'cdf.tar.gz'])
            srcName = os.path.join(pathDatabase, fileName)
    #        tar = tarfile.open()
            cmd = ('tar -xkf {:s} --strip-components 2 ' +
                   '--one-top-level={:s}').format(srcName, pathWorkData)
            try:
                os.popen(cmd)
            except Exception:
                with open('extraction.error', 'a') as f:
                    f.write('{}\n'.format(fileName))
        start = end
        if start.month == 12:
            end = datetime(start.year+1, 1, 1)
        else:
            end = datetime(start.year, start.month+1, 1)


def readCDFInfo(cdfFile):
    info = cdfFile.cdf_info()
    variables = info['zVariables']
    if len(variables) == 0:
        variables = info['rVariables']
    varInfo = [None]*len(variables)
    for i, var in enumerate(variables):
        varInfo[i] = {}
        varInfo[i]['name'] = var
        varInfo[i].update(cdfFile.varattsget(variable=var))
        varInfo[i].update(cdfFile.varinq(variable=var))
    return varInfo

##
def findFileNames(path, strings=None, size='allSize', timeTag=None):
    with os.scandir(path) as it:
        fileNamesMeetingSize = []
        foundFileNames = []
        for entry in it:
            if entry.is_file:
                if size[0] == '<' and entry.stat().st_size < int(size[1:]):
                    fileNamesMeetingSize.append(entry.name)
                elif size[0] == '>' and entry.stat().st_size > int(size[1:]):
                    fileNamesMeetingSize.append(entry.name)
                elif size == 'allSize':
                    fileNamesMeetingSize.append(entry.name)
        if strings:
            if not isinstance(strings[0], list):
                initDeltaTime = timedelta(minutes=20)
                deltaTime = initDeltaTime
                while fileNamesMeetingSize:
                    if all(string in fileNamesMeetingSize[0] for string in strings):
                        if timeTag:
                            items = fileNamesMeetingSize[0].split('_')
                            for item in items:
                                if item.isnumeric():
                                    year = int(item[:4])
                                    month = int(item[4:6])
                                    day = int(item[6:8])
                                    hour = int(item[8:10])
                                    minute = int(item[10:12])
                                    second = int(item[12:14])
                                    timeOfFile = datetime(year, month, day, hour, minute, second)
                                    deltaTime_ = timeTag - timeOfFile
                                    if deltaTime_ > timedelta(seconds=1):
                                        if deltaTime_ < deltaTime:
                                            if deltaTime == initDeltaTime:
                                                deltaTime = deltaTime_
                                                foundFileNames.append(fileNamesMeetingSize.pop(0))
                                            else:
                                                deltaTime = deltaTime_
                                                del foundFileNames[-1]
                                                foundFileNames.append(fileNamesMeetingSize.pop(0))
                                        else: del fileNamesMeetingSize[0]
                                    else: del fileNamesMeetingSize[0]
                        else:
                            foundFileNames.append(fileNamesMeetingSize.pop(0))
                    else: del fileNamesMeetingSize[0]
            else:
                while fileNamesMeetingSize:
                    for stringList in strings:
                        if all(string in fileNamesMeetingSize[0] for string in strings):
                            foundFileNames.append(fileNamesMeetingSize.pop(0))
                            break
                    else: del fileNamesMeetingSize[0]
        return foundFileNames
##

def findFileNamesInInterval(path, interval=None, strings=None, size='allSize', findAll=False):
    foundFileNames = []
    criteria_ = []
    if strings:
        criteria_.extend(strings)
    if interval:
        start = interval[0]
        end = interval[1]
        unFoundTimeTags = []
        deltaTime = end - start
        numberOfDays = deltaTime.days
        if deltaTime.total_seconds() % (3600*24):
            numberOfDays += 1
        day = timedelta(days=1)
        timeTags = [[(start+i*day).strftime("%Y%m%d"), 
                    (start+(i+1)*day).strftime("%Y%m%d")]
                    for i in range(numberOfDays)]
    while timeTags:
        criteria = criteria_.copy()
        criteria.extend(timeTags[0])
        foundFileName = findFileNames(path, strings=criteria, size=size)
        if foundFileName:
            foundFileNames.extend(foundFileName)
            del timeTags[0]
        else:
            unFoundTimeTags.append(timeTags.pop(0))
    if findAll and unFoundTimeTags:
        raise NameError('files not found: {}'.format(unFoundTimeTags))
    elif unFoundTimeTags:
        print('warning! files not found: {}'.format(unFoundTimeTags))
    return foundFileNames

##
def defineFiles(workDataDir, spacecrafts, datasets, epoch, size='allSize', silence=True):
    start = datetime(*cdflib.cdfepoch.breakdown(epoch)[:6])
    end = start + timedelta(seconds=2)
    interval = [start, end]
    if any([dataset in datasets[0] for dataset in ['5VPS', 'CIS']]):
#        workDataDir = '/home/yufei/Documents/MyFiles/works/python/bowShock/data'
        paths = [os.path.join(workDataDir, spacecraft, dataset)
                 for spacecraft, dataset in zip(spacecrafts, datasets)]
    elif 'FULL' in datasets[0]:
        paths = [os.path.join(workDataDir, spacecraft, dataset, end.strftime('%Y'), end.strftime('%m')) for spacecraft, dataset in zip(spacecrafts, datasets)]
    elif 'mms' in spacecrafts[0]:  ## dataset in the form 'fgm/srvy/l2', and the directories under it is year/month
        paths = [os.path.join(workDataDir, spacecraft, *dataset.split('/'), end.strftime('%Y'), end.strftime('%m')) for spacecraft, dataset in zip(spacecrafts, datasets)]
#    fileSizeCriterion = '>'+int(10**6)
    fileNamesList = [None]*len(spacecrafts)
    cdfFileList = [None]*len(spacecrafts)
    for i in range(len(spacecrafts)):
        criteria = datasets[i].split('/')
        path = paths[i]
        if 'fgm' in datasets[i]:
            if 'brst' in datasets[i]:
                criteria.append(start.strftime("%Y%m%d"))
                timeTag = start
            elif 'srvy' in datasets[i]:
                criteria.append(start.strftime("%Y%m%d"))
                timeTag = None
            fileNamesList[i] = findFileNames(path, strings=criteria, size=size, timeTag=timeTag)
            if not silence:
                print(fileNamesList[i][0])
            cdfFileList[i] = cdflib.CDF(os.path.join(path, fileNamesList[i][0]))
        elif 'FGM' in datasets[i]:
            criteria.append(start.strftime("%Y%m%d"))
            fileNamesList[i] = findFileNamesInInterval(path, interval=interval, strings=datasets[i])
            if not silence:
                print(fileNamesList[i][0])
            cdfFileList[i] = cdflib.CDF(os.path.join(path, fileNamesList[i][0]))
        elif 'fpi' in datasets[i]:
            if 'fast' in datasets[i]:
                cdfFiles_ = []
                for fileName in fileNamesList[i]:
                    if not silence:
                        print(fileName)
                    cdfFiles_.append(cdflib.CDF(os.path.join(path, fileName)))
                cdfFileList[i] = cdfFiles_
            elif 'brst' in datasets[i]:
                criteria.append(start.strftime("%Y%m%d"))
                timeTag = start
                fileNamesList[i] = findFileNames(path, strings=criteria, size=size, timeTag=timeTag)
                if not silence:
                    print(fileNamesList[i][0])
                cdfFileList[i] = cdflib.CDF(os.path.join(path, fileNamesList[i][0]))
    return cdfFileList
##

def datetime2list(dateTime):
    return [dateTime.year, dateTime.month, dateTime.day,
            dateTime.hour, dateTime.minute, dateTime.second,
            dateTime.microsecond//10**3]


#def xyz2spherical(unitVectors):
#    if len(unitVectors.shape) == 1:
#        unitVectors = unitVectors[None, :]
#    unitVectors *= np.sign(unitVectors[:, 0])[:, None]
#    thetas = np.arccos(unitVectors[:, 2])
#    phis = np.arcsin(unitVectors[:, 1]/np.sin(thetas))
#    return np.array([thetas, phis]).squeeze().T

##
def xyz2spherical(unitVectors, halfSphere=True, twopi=True):
    '''twopi means that phi varies from zero to 2pi'''
    if len(unitVectors.shape) == 1:
        unitVectors = unitVectors[None, :]
    sign0 = np.sign(unitVectors[..., 0])[..., None]
    unitVectors = unitVectors*sign0
    thetas = np.arccos(unitVectors[..., 2])
    phis = np.arcsin(unitVectors[..., 1]/np.sin(thetas))
    if not halfSphere:
        sign0 = sign0.squeeze()
        thetas = sign0*thetas + np.pi/2*(1-sign0)
        signPhis = np.sign(phis)
        phis = (1+sign0)/2*(phis+(1-signPhis)*np.pi) + (1-sign0)/2*(-phis+np.pi+np.pi/2*signPhis)
        if not twopi:
            sign1 = np.sign(np.pi-phis)
            phis -= (1-sign1)*np.pi
    return np.stack([thetas, phis], axis=-1).squeeze()
##

def spherical2xyz(vectors):
    shape = vectors.shape
    if len(shape) == 1:
        vectors = vectors[None, :]
    if shape[-1] == 2:
        x = np.sin(vectors[..., 0])*np.cos(vectors[..., 1])
        y = np.sin(vectors[..., 0])*np.sin(vectors[..., 1])
        z = np.cos(vectors[..., 0])
    elif shape[-1] == 3:
        x = vectors[..., 0]*np.sin(vectors[..., 1])*np.cos(vectors[..., 2])
        y = vectors[..., 0]*np.sin(vectors[..., 1])*np.sin(vectors[..., 2])
        z = vectors[..., 0]*np.cos(vectors[..., 1])
    return np.stack([x, y, z], axis=-1).squeeze()


def normalized(array, axis=-1):
    norm_ = np.linalg.norm(array, axis=axis)
    dim = len(array.shape)
    if dim == 1:
        array = array / norm_
    else:
        array = array / norm_[..., None]
    return  array


def align(vectorListA, vectorListB):
    'align B with A'
    mainComponent_ = np.argmax(np.abs(vectorListA).flatten())
    mainVector = mainComponent_ // 3
    mainComponent = mainComponent_ % 3
    sign0 = np.sign(vectorListA[mainVector, mainComponent])
    leng = len(vectorListB.shape)
#    if leng == 1:
#        BAlignedWithA = vectorListB[None, :]
#    elif leng == 2:
#        BAlignedWithA = vectorListB
#    else:
#        raise Exception('attention!')
    for vectorB in vectorListB:
        if np.sign(vectorB[mainComponent]) == -sign0:
            vectorB *= -1
    return np.array([mainVector, mainComponent])


def transNormal(normalList):
    'make vectors in normalList point at the same direction'
    quality = True
    mainVectorComponent = align(normalList, normalList)
    if normalList[mainVectorComponent[0], 0] < 0:
        normalList[:] *= -1
    sign0 = np.sign(normalList[0])
    for normal in normalList:
        sign_ = np.sign(normal)
        if any(sign_ == -sign0):
            quality = False
    return quality


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff /nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5, axis=-1):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data, axis=axis)
    return y


def inclusive_range(*args):
    numargs = len(args)
    if numargs == 0:
        raise TypeError("you need to write at least a value")
    elif numargs == 1:
        stop = args[0]
        start = 0
        step = 1
    elif numargs == 2:
        (start, stop) = args
        step = 1
    elif numargs == 3:
        (start, stop, step) = args
    else:
        raise TypeError("Inclusive range was expected at most 3 arguments,got {}".format(numargs))
    i = start
    while i <= stop:
        yield i
        i += step


def includeEnd(len2List):
    return [len2List[0], len2List[1]+1]


def format_day(t, pos=None):
    return cdflib.cdfepoch.encode(t)[11:16]


def format_dayTT2000(t, pos=None):
    return cdflib.cdfepoch.encode(int(t))[11:16]


def format_hour(t, pos=None):
    return cdflib.cdfepoch.encode(t)[14:19]


def format_hourTT2000(t, pos=None):
    return cdflib.cdfepoch.encode(int(t))[14:19]


def format_min(t, pos=None):
    return cdflib.cdfepoch.encode(t)[17:21]


def format_minTT2000(t, pos=None):
    return cdflib.cdfepoch.encode(int(t))[17:21]


def format_hourMinTT2000(t, pos=None):
    return cdflib.cdfepoch.encode(int(t))[14:21]

dayFormatter = FuncFormatter(format_day)
dayFormatterTT2000 = FuncFormatter(format_dayTT2000)
hourFormatter = FuncFormatter(format_hour)
hourFormatterTT2000 = FuncFormatter(format_hourTT2000)
minFormatter = FuncFormatter(format_min)
minFormatterTT2000 = FuncFormatter(format_minTT2000)
hourMinFormatterTT2000 = FuncFormatter(format_hourMinTT2000)


def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def gauss_fit(x, y):
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
    return popt


def allocateToBins(arrayA, axis=0, numberOfBins=5, scale='linear'):
    bins = []
    if scale == 'log':
        arrayA = np.log(arrayA)
    binLimits = np.linspace(arrayA.min(), arrayA.max(), numberOfBins+1)
    for i in range(numberOfBins):
        binLeft = binLimits[i]
        binRight = binLimits[i+1]
        if i+1 < numberOfBins:
            bins.append(arrayA[(binLeft <= arrayA) & (arrayA < binRight)])
        else:
            bins.append(arrayA[(binLeft <= arrayA) & (arrayA <= binRight)])
    if scale == 'log':
        binLimits = np.exp(binLimits)
        for bin_ in bins:
            bin_ = np.exp(bin_)
    return bins, binLimits


def showIdentifiedEvents(workDataDir, spacecrafts, epoch, tBAllSpacecrafts, tBList, bMagList, xList, plotPlasma=True, allowNotPlotPlasma=False, bUpperLimit=60):
#    datasets = ['C1_CP_FGM_5VPS']
#    epochList = np.repeat(epoch, 4)
#    cdfFileList = defineFiles(workDataDir, spacecrafts, datasets, epochList[0])
#    tBList, bMagList, xList = readData(spacecrafts, datasets, cdfFileList, epochList, lowpass=None, variables=['tBList', 'bMagList', 'xList'])
    if plotPlasma:
        if 'C' in spacecrafts[0]:
            cdfFileCIS = defineFiles(workDataDir, ['C1'], ['C1_CP_CIS-HIA_ONBOARD_MOMENTS'], epoch)[0]
            tP = cdfFileCIS.varget(0)
            if len(tP) < 10:
                if allowNotPlotPlasma:
                    plotPlasma = False
                    print('warning! less than 10 points in CIS')
                else:
                    raise Exception('less than 10 points in CIS')
            density = cdfFileCIS.varget(4)
            velocity = np.empty([len(tP), 4])
            velocity[:, :3] = cdfFileCIS.varget(6)
            velocity[:, 3] = np.sqrt(np.sum(velocity[:, :3]**2, axis=1))
        elif 'mms' in spacecrafts[0]:
            plasmaDatasets = ['fpi/fast/l2/dis-moms']*4
            plasmaCDFFileList = defineFiles(workDataDir, spacecrafts, plasmaDatasets, epoch)
            tPList, densityList, velocityList = readData(spacecrafts, plasmaDatasets, plasmaCDFFileList, variables=['tPList', 'densityList', 'velocityList'])
            tP = tPList[0]
            density = densityList[0]
            velocity = np.empty([len(tP), 4])
            velocity[:, :3] = velocityList[0]
            velocity[:, 3] = np.sqrt(np.sum(velocity[:, :3]**2, axis=1))
    if 'C' in spacecrafts[0]:
        myFormatter = dayFormatter
    elif 'mms' in spacecrafts[0]:
        myFormatter = dayFormatterTT2000
    colors = ['k', 'b', 'r', 'c']
    linestyles = ['-', '--', '-.', ':']
    default_cycler = (cycler(color=colors) +
                      cycler(linestyle=linestyles))
    color = 'g'
    plt.rc('lines', linewidth=0.4)
    plt.rc('axes', prop_cycle=default_cycler)
    fig, axes = plt.subplots(nrows=3, sharex=True)
    axX, axPd, axB = axes
    axPv = axPd.twinx()  # axPv = axes of particle velocity
    lw = 0.4
    left = 0.15
    bottom = 0.1
#    height = 0.4
#    width = 0.7
#    wholeWidth = 0.7
    plt.subplots_adjust(left=left, bottom=bottom, hspace=0.05)
    for i in range(len(spacecrafts)):
        axB.plot(tBList[i], bMagList[i], lw=lw, label='C{}'.format(i+1))
    xyz = ['x', 'y', 'z']
    for k in range(3):
        axX.plot(tBList[0], xList[0][:, k]/6371, lw=lw, label=xyz[k])
    xlim = axB.get_xlim()
    axB.xaxis.set_major_formatter(myFormatter)
    axB.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axB.tick_params(which='both', direction='in', top=True, right=True)
    axB.set_ylim(0, bUpperLimit)
    axB.set_ylabel('B [nT]')
    axB.legend(frameon=False)
    if 'C' in spacecrafts[0]:
        beginTime = cdflib.cdfepoch.encode(xlim[0])[:-4]
    elif 'mms' in spacecrafts[0]:
        beginTime = cdflib.cdfepoch.encode(int(xlim[0]))[:-4]
    axB.set_xlabel('time [min:s] since {}'.format(beginTime))

    axX.set_ylabel('x [$R_E$]')
    axX.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axX.tick_params(which='both', direction='in', top=True, right=True)
    axX.tick_params(which='both', direction='in', top=True, right=True, labelbottom=False)
    axX.xaxis.set_major_formatter(myFormatter)
    axX.legend(frameon=False)
    if plotPlasma:
        axPd.plot(tP, density, lw=lw, label='density') #axPd.set_xlim(xlim)
        axPd.xaxis.set_major_formatter(myFormatter)
        axPd.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axPd.tick_params(which='both', direction='in', top=True)
        axPd.set_ylabel('density [$cm^{-3}$]')
        #axPd.legend(loc='upper left')
        axPv.plot(tP, velocity[:, 3], color=color, lw=lw, label='velocity')
        axPv.tick_params(which='both', direction='in', right=True)
        axPv.tick_params(axis='y', labelcolor=color)
        axPv.set_ylabel('velocity [km/s]')
    ylim = axB.get_ylim()
    if tBAllSpacecrafts:
        for i in range(len(spacecrafts)):
            for tB in tBAllSpacecrafts[i]:
                axB.plot([tB]*2, ylim + 0.2*np.diff(ylim)*np.array([-1, 1]), color=colors[i], ls=linestyles[i], lw=lw)
    return fig


def importKnownEvent(event, maxHalfWindowLength=2, spacecrafts=['C1', 'C2', 'C3', 'C4']):
    epochLimList = np.zeros((len(spacecrafts), 2))
    epochList = np.zeros((len(spacecrafts)))
    for i in range(len(spacecrafts)):
        dateAndTime = event[i].split()
        date = dateAndTime[0].split('-')
        microsecond = dateAndTime[1].split('.')
        time = microsecond[0].split(':')
        eventEpochVector = list(map(int, date + time + [microsecond[1]]))
        eventEpoch = cdflib.cdfepoch.compute_epoch(eventEpochVector)
        epochList[i] = eventEpoch
        epochLimList[i, :] = np.array([eventEpoch, eventEpoch]) + 1000* np.array([-maxHalfWindowLength, maxHalfWindowLength])
#    eventEpoch =  np.mean(epochList)
    return epochLimList, epochList


def multirowPlot(rowNumber):
    default_cycler = (cycler(color=['k', 'b', 'r', 'c']) +
                      cycler(linestyle=['-', '--', '-.', ':']))
    plt.rc('lines', linewidth=0.4)
    plt.rc('axes', prop_cycle=default_cycler)
    plt.close(bMagfig)
    bMagfig = plt.figure()
    lw = 0.4
    left = 0.15
    bottom = 0.1
    height = 0.4
    width = 0.7
    wholeWidth = 0.7
    axB = bMagfig.add_axes([left, bottom, width, height])
    axX = bMagfig.add_axes([left, bottom+height, width, height])


def showEvent(workDataDir, spacecrafts, datasets, epochList, lowpass=None, onlyIonAndBMag=False, direction=None):
    'plot CIS, bMagList, dB/dt, zoomIn bMagList near the event'
    datasets = ['C1_CP_FGM_5VPS', 'C2_CP_FGM_5VPS',
                'C3_CP_FGM_5VPS', 'C4_CP_FGM_5VPS']
    epoch = epochList[0]
    cdfFileList = defineFiles(workDataDir, spacecrafts, datasets, epoch)
    tBList, tStep, bMagList, bVectorsList, xList = readData(spacecrafts, datasets, cdfFileList, epochList, lowpass=lowpass, variables=['tBList', 'tStep', 'bMagList', 'bVectorsList', 'xList'], tStepAtStart=False, gsm=False, auxFile=None)
    # CIS data
    cdfFileCIS = defineFiles(workDataDir, ['C1'], ['C1_CP_CIS-HIA_ONBOARD_MOMENTS'], epoch)[0]
    tP = cdfFileCIS.varget(0)
    density = cdfFileCIS.varget(4)
    velocity = np.empty([len(tP), 4])
    velocity[:, :3] = cdfFileCIS.varget(6)
    velocity[:, 3] = np.sqrt(np.sum(velocity[:, :3]**2, axis=1))

    colors = ['k', 'b', 'r', 'c']
    linestyles = ['-', '--', '-.', ':']
    default_cycler = (cycler(color=colors) +
                      cycler(linestyle=linestyles))
    color = 'g'
    plt.rc('lines', linewidth=0.4)
    plt.rc('axes', prop_cycle=default_cycler)
    fig = plt.figure(figsize=(10, 6))
    lw = 0.4
    tIndicesList = [None]*len(spacecrafts)
    xlim = np.array([epochList.min(), epochList.max()]) + 20*60*1000*np.array([-1, 1])
    if onlyIonAndBMag:
        bMagWidth = 0.7
    else:
        bMagWidth = 0.3
    for i in range(len(spacecrafts)):
        tIndicesList[i] = np.arange(*includeEnd([np.argmin(np.abs(tBList[i]-epoch_)) for epoch_ in xlim]))
    axB = fig.add_axes([0.1, 0.1, bMagWidth, 0.5])
    axPd = fig.add_axes([0.1, 0.6, bMagWidth, 0.3])  # axPd = axes of particle density
    axPv = axPd.twinx()  # axPv = axes of particle velocity
    # plot
    for i in range(len(spacecrafts)):
        axB.plot(tBList[i][tIndicesList[i]], bMagList[i][tIndicesList[i]], lw=lw, label='C{}'.format(i+1))
        ylim = axB.get_ylim()
        axB.set_ylim(ylim)
        axB.plot([epochList[i]]*2, ylim + 0.2*np.diff(ylim)*np.array([-1, 1]), color=colors[i], ls=linestyles[i], lw=lw)
    axB.set_xlim(xlim)
    if np.diff(axB.get_xlim()) > 3600*1000:
        myFormatter = dayFormatter
    else: myFormatter = hourFormatter
    axB.xaxis.set_major_formatter(myFormatter)
    axB.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axB.tick_params(which='both', direction='in', top=True, right=True)
    axB.set_xlabel('time [min:s] since {}'.format(cdflib.cdfepoch.encode(xlim[0])[:-4]))
    axB.set_ylabel('B [nT]')
    axB.legend(loc='upper left', frameon=False)

    axPd.plot(tP, density, lw=lw, label='density')
    axPd.set_xlim(xlim)
    axPd.xaxis.set_major_formatter(myFormatter)
    axPd.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axPd.tick_params(which='both', direction='in', top=True, labelbottom=False)
    axPd.set_ylabel('density [$cm^{-3}$]')
    #axPd.legend(loc='upper left')
    axPv.plot(tP, velocity[:, 3], color=color, lw=lw, label='velocity')
    axPv.tick_params(axis='y', labelcolor=color)
    axPv.set_ylabel('velocity [km/s]')
    axPv.tick_params(axis='y', colors='g')
    axPv.yaxis.label.set_color('g')
    #axPv.legend(loc='upper right')
    if not onlyIonAndBMag:
        zoomInHalfWindow = 9
        zoomInTIndicesList = [None]*len(spacecrafts)
        zoomInXlim = np.array([epochList.min(), epochList.max()]) + 1000*np.array([-zoomInHalfWindow, zoomInHalfWindow])
        dBList = [None]*len(spacecrafts)
        for i in range(len(spacecrafts)):
            zoomInTIndices = np.arange(*includeEnd([np.argmin(np.abs(tBList[i]-epoch_)) for epoch_ in zoomInXlim]))
            zoomInTIndicesList[i] = zoomInTIndices
            dBList[i] = 1000*np.diff(bMagList[i][np.append(zoomInTIndices, zoomInTIndices.max()+1)]) / np.diff(tBList[i][np.append(zoomInTIndices, zoomInTIndices.max()+1)])
        axZoomInB = fig.add_axes([0.53, 0.1, 0.4, 0.5])
        axZoomInDB = fig.add_axes([0.53, 0.6, 0.4, 0.3])
        for i in range(len(spacecrafts)):
            axZoomInB.plot(tBList[i][zoomInTIndicesList[i]], bMagList[i][zoomInTIndicesList[i]], lw=lw, label='C{}'.format(i+1))
            axZoomInDB.plot(tBList[i][zoomInTIndicesList[i]], dBList[i], lw=lw, label='C{}'.format(i+1))
        #axZoomInB.set_xlim(zoomInXlim)
        #axZoomInB.set_ylim([0, tBList[0][zoomInXlim])
        axZoomInB.xaxis.set_ticks(np.linspace(*zoomInXlim, 4))
        axZoomInB.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axZoomInB.xaxis.set_major_formatter(minFormatter)
    #    axZoomInB.xaxis.set_minor_formatter(minFormatter)
        axZoomInB.tick_params(which='both', direction='in', top=True, right=True)
        axZoomInB.set_xlabel('time [s] since {}'.format(cdflib.cdfepoch.encode(zoomInXlim[0])[14:19]))
        axZoomInB.set_ylabel('B [nT]')
        axZoomInDB.xaxis.set_ticks(np.linspace(*zoomInXlim, 4))
        axZoomInDB.tick_params(which='both', direction='in', top=True, right=True, labelbottom=False)
        axZoomInDB.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axZoomInDB.xaxis.set_major_formatter(minFormatter)
        #axZoomInDB.get_xaxis().set_visible(False)
        axZoomInDB.set_ylabel('dB/dt [nT/s]')
# compute additional information
    xGSE = np.empty((len(spacecrafts), 3))
    intervalLength = 120
    numberOfPoints = int(intervalLength*1000/tStep)
    lowerLimit_ = int(intervalLength/2*1000/tStep)
    bUpstream = np.empty((len(spacecrafts), 3))
    for i in range(len(spacecrafts)):
        indexOfCenter_ = np.argmin(np.abs(tBList[i]-epochList[i]))
        xGSE[i] = xList[i][indexOfCenter_]/6371
        if direction == None:
            out = None
            if bMagList[i][indexOfCenter_ + lowerLimit_] < bMagList[i][indexOfCenter_ - lowerLimit_]:
                if out != None and out is False:
                    raise Exception('inconsistent B data')
                out = True
            elif bMagList[i][indexOfCenter_ + lowerLimit_] > bMagList[i][indexOfCenter_ - lowerLimit_]:
                if out != None and out is True:
                    raise Exception('inconsistent B data')
                out = False
        elif direction == 'inbound':
            out = False
        elif direction == 'outbound':
            out = True
        if out:
            rangeOfB = np.arange(indexOfCenter_ + lowerLimit_, indexOfCenter_ + lowerLimit_ + numberOfPoints)
        elif out is False:
            rangeOfB = np.arange(indexOfCenter_ - lowerLimit_ -numberOfPoints, indexOfCenter_ - lowerLimit_)
        bUpstream[i, :] = np.mean(bVectorsList[i][rangeOfB], axis=0)
    timingVelocity, timingNormal, timingShape = timingVelocityAndNormal(epochList/1000, xGSE, silence=False, getShape=True)
    timingVelocity = timingVelocity*6371
#    combs = combinations(range(4), 2)
#    number_ = np.math.factorial(4)//np.math.factorial(2)//np.math.factorial(2)
#    diffX = np.empty((number_, 3))
#    difft = np.empty(number_)
#    for i, comb in enumerate(list(combs)):
#        diffX[i] = x[comb[0]] - x[comb[1]]
#        difft[i] = (epochList[comb[0]] - epochList[comb[1]])/1000
#    m = 1/16 * difft[None, :] @ diffX @ RInverse
#    timingVelocity = 1/np.linalg.norm(m)*6371
#    timingNormal = m*timingVelocity/6371
    xGSECenter = np.mean(xGSE, axis=0)
    thetaBn = np.arccos(normalized(bUpstream) @ timingNormal)
    thetaBn = np.where(thetaBn > np.pi/2, np.pi-thetaBn, thetaBn) 
    threeItemsF = '[' +'{{:.{d}f}}, '*2 + '{{:.{d}f}}' + ']'
    fig.text(0.1, 0.95, 'constellation position:\n [{:.2f},{:.2f},{:.2f}] [$R_E$]'.format(*[xGSECenter[i] for i in range(3)]))
    fig.text(0.3, 0.95, 'normal:\n [{:.2f},{:.2f},{:.2f}]'.format(*[timingNormal[i] for i in range(3)]))
    fig.text(0.5, 0.95, 'velocity:\n {:.0f} [km/s]'.format(timingVelocity))
    fig.text(0.6, 0.95, r'$\theta_{{Bn}}$:' '\n' r'{:.0f}$^\circ$,{:.0f}$^\circ$,{:.0f}$^\circ$,{:.0f}$^\circ$'.format(*[thetaBn[i]/np.pi*180 for i in range(4)]))
    fig.text(0.82, 0.35, 'shape:\n '+threeItemsF.format(d=0).format(*[timingShape[0][i]*100 for i in range(3)])+' [$10^{{-2}}R_E$]\n '+'\n '.join([threeItemsF.format(d=2).format(*timingShape[1][:, i]) for i in range(3)]))
    return fig, xGSECenter, timingVelocity, timingNormal, timingShape, thetaBn

##
def gse2sm(vectors, tB, auxFile):
    vectorsInSM = np.zeros_like(vectors)
    t_aux = auxFile.varget(2)
    gse_gsm = auxFile.varget(30)
    dipole = auxFile.varget(31)
    gse_gsm = np.interp(tB, t_aux, gse_gsm).squeeze()
    dipole = np.interp(tB, t_aux, dipole).squeeze()
    vectorsInSM[:, 1] = np.cos(gse_gsm * np.pi/180) * vectors[:, 1] - np.sin(gse_gsm * np.pi/180) * vectors[:, 2]
    Bz1 = np.sin(gse_gsm * np.pi/180) * vectors[:, 1] + np.cos(gse_gsm * np.pi/180) * vectors[:, 2]
    vectorsInSM[:, 0] = np.cos(dipole * np.pi/180) * vectors[:, 0] - np.sin(dipole * np.pi/180) * Bz1
    vectorsInSM[:, 2] = np.sin(dipole * np.pi/180) * vectors[:, 0] + np.cos(dipole * np.pi/180) * Bz1
    return vectorsInSM
##

def gse2gsm(vectors, tB, auxFile):
    vectorsInGSM = np.zeros_like(vectors)
    t_aux = auxFile.varget(2)
    gse_gsm = auxFile.varget(30)/180*np.pi
    gse_gsm = np.interp(tB, t_aux, gse_gsm)
    cosTheta = np.cos(gse_gsm)
    sinTheta = np.sin(gse_gsm)
    vectorsInGSM[:, 0] = vectors[:, 0]
    vectorsInGSM[:, 1] = cosTheta * vectors[:, 1] - sinTheta * vectors[:, 2]
    vectorsInGSM[:, 2] = sinTheta * vectors[:, 1] + cosTheta * vectors[:, 2]
    return vectorsInGSM
##

##
def readData(spacecrafts, datasets, cdfFileList, epochList=None, lowpass=None, variables=['tBList', 'tStep', 'bVectorsList', 'xList'], tStepAtStart=False, coordinateSystem='gse', auxFile=None, evenlySpaced=False):
    tBList = [None]*len(spacecrafts)
    tPList = [None]*len(spacecrafts)
    sampleRateList = [None]*len(spacecrafts)
    bVectorsList = [None]*len(spacecrafts)
    bMagList = [None]*len(spacecrafts)
    velocityList = [None]*len(spacecrafts)
    densityList = [None]*len(spacecrafts)
    xList = [None]*len(spacecrafts)
    for i in range(len(spacecrafts)):
        if 'C' in spacecrafts[i]:
            if 'tBList' in variables:
                tBList[i] = cdfFileList[i].varget(0)
                if 'tStep' in variables:
                    if '5VPS' in datasets[0]:
                        tStep = 200
                    elif tStepAtStart:
                        tStep = (tBList[i][200] - tBList[i][0])/200
                    else:
                        tIndice = np.argmin(np.abs(tBList[i]-epochList[i]))
                        tStep = (tBList[i][tIndice+100] - tBList[i][tIndice-100])/200
            if 'bVectorsList' in variables:
                if lowpass:
                    if '5VPS' in datasets[0]:
                        bVectorsList[i] = butter_lowpass_filter(cdfFileList[i].varget(2), lowpass, 5, axis=0)
                    elif 'FULL' in datasets[0]:
                        bVectorsList[i] = butter_lowpass_filter(cdfFileList[i].varget(2), lowpass, 22, axis=0)
                else:
                    bVectorsList[i] = cdfFileList[i].varget(2)
                if coordinateSystem == 'gsm':
                    bVectorsList[i] = gse2gsm(bVectorsList[i], tBList[i], auxFile)
            if 'bMagList' in variables:
                if lowpass:
                    if '5VPS' in datasets[0]:
                        bMagList[i] = butter_lowpass_filter(cdfFileList[i].varget(3), lowpass, 5, axis=0)
                    elif 'FULL' in datasets[0]:
                        bMagList[i] = butter_lowpass_filter(cdfFileList[i].varget(3), lowpass, 22, axis=0)
                else:
                    bMagList[i] = cdfFileList[i].varget(3)
            if 'xList' in variables:
                xList[i] = cdfFileList[i].varget(4)
                if coordinateSystem == 'gsm':
                    xList[i] = gse2gsm(xList[i], tBList[i], auxFile)
        elif 'mms' in spacecrafts[i]:
            if 'fgm' in datasets[i]:
                if 'tBList' in variables:
                    tBList[i] = cdfFileList[i].varget(0)
                    if 'sampleRateList' in variables:
                        sampleRateList[i] = cdfFileList[0].varget(20)
                if 'bVectorsList' in variables:
                    if coordinateSystem == 'gse':
                        bVectors = cdfFileList[i].varget(1)[:, :3]
                    elif coordinateSystem == 'gsm':
                        bVectors = cdfFileList[i].varget(2)[:, :3]
                    if lowpass:
                        if 'srvy' in datasets[0]:
                            diffT = np.diff(tBList[i])
                            pos = np.where(diffT < 10**8)[0]
                            bVectors = np.delete(bVectors, pos[::2], axis=0)
                            bVectorsList[i] = butter_lowpass_filter(bVectors, lowpass, 8, axis=0)
                        elif 'brst' in datasets[0]:
                            bVectorsList[i] = butter_lowpass_filter(bVectors, lowpass, 128, axis=0)
                    else:
                        bVectorsList[i] = bVectors
                if 'bMagList' in variables:
                    if coordinateSystem == 'gse':
                        bMag = cdfFileList[i].varget(1)[:, 3]
                    elif coordinateSystem == 'gsm':
                        bMag = cdfFileList[i].varget(2)[:, 3]
                    if evenlySpaced:
                        diffT = np.diff(tBList[i])
                        pos = np.where(diffT < 10**8)[0]
                        bMag = np.delete(bMag, pos[::2], axis=0)
                        tBList[i] = np.delete(tBList[i], pos[::2])
                        tStep = 125001692
                    if lowpass and evenlySpaced:
                        if 'srvy' in datasets[0]:
                            bMagList[i] = butter_lowpass_filter(bMag, lowpass, 8, axis=0)
                    else:
                        bMagList[i] = bMag
                if 'xList' in variables:
                    tX = cdfFileList[i].varget(6)
                    if coordinateSystem == 'gse':
                        x = cdfFileList[i].varget(7)
                    elif coordinateSystem == 'gsm':
                        x = cdfFileList[i].varget(8)
                    cs = interpolate.CubicSpline(tX, x[:, :3])
                    xList[i] = cs(tBList[i])
            elif 'fpi' in datasets[i]:
                if coordinateSystem == 'gse':
                    if 'fast' in datasets[i]:
                        if 'tPList' in variables:
                            tP_ = []
                            for cdfFile in cdfFileList[i]:
                                tP_.append(cdfFile.varget(0))
                            tp_ = np.concatenate(tP_)
                            permutation = np.argsort(tp_)
                            tPList[i] = tp_[permutation]
                        if 'velocityList' in variables:
                            v_ = []
                            for cdfFile in cdfFileList[i]:
                                v_.append(cdfFile.varget(22) - cdfFile.varget(23))
                            velocityList[i] = np.concatenate(v_)[permutation]
                        if 'densityList' in variables:
                            d_ = []
                            for cdfFile in cdfFileList[i]:
                                d_.append(cdfFile.varget(16))
                            densityList[i] = np.concatenate(d_)[permutation]
                    elif 'brst' in datasets[i]:
                        if 'tPList' in variables:
                            tPList[i] = cdfFileList[i].varget(0)
                        if 'velocityList' in variables:
                            velocityList[i] = cdfFileList[i].varget(22) - cdfFileList[i].varget(23)
                        if 'densityList' in variables:
                            densityList[i] = cdfFileList[i].varget(16)
    if isinstance(variables, list):
        scope = locals()
        return [eval(variable, scope) for variable in variables]
##
##
def identifyBSCrossingByBMag(workDataDir, spacecrafts, datasets, epochs, intervalLength=120, compressionRatio=2, maxStd=0.1, maxOuterBMag=60, allowedExtraRange=0.5, displacements=[0], minSheathBMag=10, getArrayList=False, silence=False, show=False, allowNotPlotPlasma=False, evenlySpaced=True, lowpass=None):
    in_out = ['inbound', 'outbound']
    dictOfRawUnidentifiedEvents = {}
    for in_out_ in in_out:
        dictOfRawUnidentifiedEvents[in_out_] = {}
        for spacecraft in spacecrafts:
            dictOfRawUnidentifiedEvents[in_out_][spacecraft] = {}
    unidentifiedEventInd = np.zeros(4, dtype=int)
    unidentifiedEventIndIn_out = np.zeros((2, 4), dtype=int)
#    epoch = epochs[0]
    for epoch in epochs:
        if not silence:
            print(cdflib.cdfepoch.encode(epoch))
        epochList = np.repeat(epoch, 4)
        filesize = '>' + str(50*10**3) # in bytes
        try:
            cdfFileList = defineFiles(workDataDir, spacecrafts, datasets, epochList[0], size=filesize)
            tBList, tStep, bMagList, xList = readData(spacecrafts, datasets, cdfFileList,  epochList, lowpass=lowpass, variables=['tBList', 'tStep', 'bMagList', 'xList'], evenlySpaced=evenlySpaced)
            dataIsObtained = True
        except:
            dataIsObtained = False
            print("warning! didn't get data")
        if not dataIsObtained:
            continue
        else:
            for scInd in range(len(spacecrafts)):
                if 'C' in spacecrafts[scInd]:
                    average = int(intervalLength*1000 / tStep)
                elif 'mms' in spacecrafts[scInd]:
                    average = int(intervalLength*10**9 / tStep)
                for displacement in displacements:
                    bMagWithDisplacement = bMagList[scInd][int(average*displacement):]
                    tBWithDisplacement = tBList[scInd][int(average*displacement):]
                    xWithDisplacement = xList[scInd][int(average*displacement):]
                    lengthOfBMag = len(bMagWithDisplacement)
                    bMagListAfterDrop = bMagWithDisplacement[:-(lengthOfBMag % average)].reshape((-1, average))
                    tBListAfterDrop = tBWithDisplacement[:-(lengthOfBMag % average)].reshape((-1, average))
                    xListAfterDrop = xWithDisplacement[:-(lengthOfBMag % average)].transpose((1, 0)).reshape((3, -1, average)).transpose((1, 2, 0))
                    bMagListAveraged = np.mean(bMagListAfterDrop, axis=1)
                    bMagListStd = np.std(bMagListAfterDrop, axis=1) / bMagListAveraged
                    ratio = bMagListAveraged[2:] / bMagListAveraged[:-2]
                    bMagMean = (bMagListAveraged[2:] + bMagListAveraged[:-2])/2
                    numberOfEffectiveAverages = len(ratio)
                    magAndStdMask = (bMagListAveraged < maxOuterBMag) & (bMagListStd < maxStd)
                    sheathBMagMask = bMagListAveraged > minSheathBMag
                    tAverageInd = np.arange(len(ratio))
                    rangeMask = tBWithDisplacement[(tAverageInd + 3)*average - 1] - tBWithDisplacement[(tAverageInd + 1)*average] < tStep*(2*average + allowedExtraRange*average)
                    inMask = (ratio > compressionRatio) & magAndStdMask[:-2] & rangeMask & sheathBMagMask[2:]
                    outMask = (ratio < 1/compressionRatio) & magAndStdMask[2:] & rangeMask & sheathBMagMask[:-2]
                    masks = [inMask, outMask]
                    tAverageInd = [None]*2  # tAverageInd[0] for inbound motion [1] for outbound motion
                    tBs = [None]*2  # same as the line above (tAverageInd)
                    xs = [None]*2
                    forePostTraversalBAverages = [None]*2
                    identifiedEvents = [None]*2
                    for in_outInd, mask in enumerate(masks):
                        dictUnderDirection = dictOfRawUnidentifiedEvents[in_out[in_outInd]]
                        tAverageInd[in_outInd] = np.where(mask)[0]
                        bMags = bMagMean[tAverageInd[in_outInd]]
                        diffBMag = np.abs(bMagListAfterDrop[tAverageInd[in_outInd]+1]-bMags[:, None])
                        columns = np.argmin(diffBMag, axis=1)
                        meanMask = diffBMag[np.arange(len(columns)), columns]/bMags < 0.01
                        tAverageInd[in_outInd] = tAverageInd[in_outInd][meanMask]
                        columns = columns[meanMask]
                        forePostTraversalBAverages[in_outInd] = bMagListAveraged[tAverageInd[in_outInd][:, None]+np.array([0, 2])[None, :]]
                        tBs[in_outInd] = tBListAfterDrop[tAverageInd[in_outInd]+1, columns]
                        xs[in_outInd] = xListAfterDrop[tAverageInd[in_outInd]+1, columns, :]
                        for eventInd_ in range(len(tBs[in_outInd])):
                            tB = tBs[in_outInd][eventInd_]
                            x = xs[in_outInd][eventInd_]
                            forePostTraversalBAverage = forePostTraversalBAverages[in_outInd][eventInd_]
                            dictUnderDirection[spacecrafts[scInd]][str(unidentifiedEventIndIn_out[in_outInd, scInd])] = {'epoch': tB, 'xGSE': (x/6371).tolist(), 'forePostTraversalBAverage': forePostTraversalBAverage, 'eventInd': unidentifiedEventInd[scInd]}
                            unidentifiedEventInd[scInd] += 1
                            unidentifiedEventIndIn_out[in_outInd, scInd] +=1
                    if not silence:
                        print(spacecrafts[scInd])
                        print('number of events: {}'.format(unidentifiedEventInd[scInd]))
    dictOfUnidentifiedEvents = {}
    eventEpochsList = [[], []]
    eventXsList = [[], []]
    if show:
        getArrayList = True
    for in_outInd in range(2):
        direction = in_out[in_outInd]
        dictOfUnidentifiedEvents[direction] = {}
        for spacecraft in spacecrafts:
            if 'C' in spacecraft:
                repeatRange = intervalLength * 1000
            if 'mms' in spacecraft:
                repeatRange = intervalLength * 1000 * 10**6
            dictOfRepeatedEvents = dictOfRawUnidentifiedEvents[direction][spacecraft]
            if getArrayList:
                dictOfUnidentifiedEvents[direction][spacecraft], eventEpochs_, eventXs_ = dropRepeatedEventsInADict(dictOfRepeatedEvents, repeatRange, getArray=getArrayList, getDict=True)
                eventEpochsList[in_outInd].append(eventEpochs_)
                eventXsList[in_outInd].append(eventXs_)
            else:
                dictOfUnidentifiedEvents[direction][spacecraft], = dropRepeatedEventsInADict(dictOfRepeatedEvents, repeatRange, getArray=getArrayList, getDict=True)
    if show:
        tBAllSpacecrafts = [None]*len(spacecrafts)
        if len(epochs) != 1:
            print("Too many days")
        else:
            for scInd in range(len(spacecrafts)):
                tBAllSpacecrafts[scInd] = np.hstack([eventEpochsList[k][scInd] for k in range(2)])
#                    print(in_out[in_outInd])
#                    print(cdflib.cdfepoch.encode(tInd[in_outInd]))
#                plt.close('all')
#                plt.ion()
            figEvents = showIdentifiedEvents(workDataDir, spacecrafts, epoch, tBAllSpacecrafts, tBList, bMagList, xList, allowNotPlotPlasma=allowNotPlotPlasma)
            return dictOfUnidentifiedEvents, figEvents
    returnedVariables = [dictOfUnidentifiedEvents]
    if getArrayList:
        returnedVariables.append(eventEpochsList)
        returnedVariables.append(eventXsList)
    return returnedVariables
##

def dict2ArrayOfEvents(dictOfEvents):
    numberOfEvents = len(dictOfEvents)
    eventEpochs_ = []
    eventXs_ = []
    for eventInd in range(numberOfEvents):
        event = dictOfEvents[str(eventInd)]
        eventEpochs_.append(event['epoch'])
        eventXs_.append(event['xGSE'])
    return np.array(eventEpochs_), np.array(eventXs_)


def dropRepeatedEventsInADict(dictOfRepeatedEvents, repeatRange=0, getArray=True, getDict=False):
    eventEpochs_, eventXs_ = dict2ArrayOfEvents(dictOfRepeatedEvents)
    permutation_ = np.argsort(eventEpochs_)
    if len(permutation_) > 0:
        mask_ = np.append([True], np.diff(eventEpochs_[permutation_]) > repeatRange)
    else:
        mask_ = np.array([], dtype=bool)
    eventEpochsNoRepeated = eventEpochs_[permutation_][mask_]
    eventXsNoRepeated = eventXs_[permutation_][mask_]
    returnedVariables = []
    if getDict:
        dictOfEventsOfASpacecraft = {}
        for ind in range(len(eventEpochsNoRepeated)):
            dictOfEventsOfASpacecraft[str(ind)] = {'epoch': eventEpochsNoRepeated[ind], 'xGSE': eventXsNoRepeated[ind]}
        returnedVariables.append(dictOfEventsOfASpacecraft)
    if getArray:
        returnedVariables.append(eventEpochsNoRepeated)
        returnedVariables.append(eventXsNoRepeated)
    return returnedVariables


def confirmBSCrossingByCIS(workDataDir, spacecrafts, epochs, direction, compressionRatio=2):
    epoch = epochs[0]
    cdfFileCIS = defineFiles(workDataDir, ['C1'], ['C1_CP_CIS-HIA_ONBOARD_MOMENTS'], epoch)[0]
    tP = cdfFileCIS.varget(0)
    if len(tP) < 10:
        bowShockCrossing = False
    else:
        density = cdfFileCIS.varget(4)
        velocity = np.empty([len(tP), 4])
        velocity[:, :3] = cdfFileCIS.varget(6)
        velocity[:, 3] = np.sqrt(np.sum(velocity[:, :3]**2, axis=1))

        diff_ = np.abs(tP[:, None] - (epoch + 1000*60*np.array([-3, -1, 1, 3]))[None, :])
        limits_ = np.argmin(diff_, axis=0)
        if any(np.diff(limits_) == 0):
            bowShockCrossing = False
        else:
            rightRange = np.arange(*limits_[2:])
            leftRange = np.arange(*limits_[:2])
            ratioDensity = np.mean(density[rightRange]) / np.mean(density[leftRange])
            ratioVelocity = np.mean(velocity[rightRange, 3]) / np.mean(velocity[leftRange, 3])
            if direction == 'inbound':
                bowShockCrossing = (ratioDensity > compressionRatio) & (ratioVelocity < 1.5/compressionRatio)
            elif direction == 'outbound':
                bowShockCrossing = (ratioDensity < 1/compressionRatio) & (ratioVelocity > compressionRatio/1.5)
    return bowShockCrossing


def mvab(bVectors, average=1):
    m = np.empty([3, 3])
    if average > 1:
        bVectors = bVectors[:-(len(bVectors) % average), :]
        bVectors_ = np.zeros((len(bVectors)//average, 3))
        for i in range(average):
            bVectors_ += bVectors[i::average, :]
        bVectors = bVectors_/average
    for j in range(3):
        m[j, :] = np.mean(bVectors[:, j, None] *
                          bVectors[:, :], axis=0) -\
                  np.mean(bVectors[:, j, None], axis=0) *\
                  np.mean(bVectors[:, :], axis=0)
    eigenSystem_ = np.linalg.eig(m)
    permutation = np.argsort(eigenSystem_[0])
    eigenValues_ = eigenSystem_[0][permutation]
    ratio = eigenValues_[1]/eigenValues_[0]
    eigenSystem = (eigenValues_, eigenSystem_[1][:, permutation])
    return eigenSystem, ratio


def computeNormalsListAndRatiosList(tBList, tStep, epochList, epochLimList, bVectorsList, spacecrafts, datasets, lowpass=None, minHalfWindowLength=np.array([0.5]*4), dilution=1, tIndicesNumber=0, dilutionOfWindow=1, average=1):
    ratiosList = [None]*len(spacecrafts)
    normalsList = [None]*len(spacecrafts)
    halfWindowLengthsList = [None]*len(spacecrafts)
    tIndicesList = [None]*len(spacecrafts)
    numberOfPointsList = [None]*len(spacecrafts)
    for i in range(len(spacecrafts)):
        if 'FULL' in datasets[0]:
            if not (tIndicesNumber is None):
                diff_ = np.abs(tBList[i]-epochList[i])
                center_ = np.argmin(diff_)
                tIndices = np.arange(-tIndicesNumber, tIndicesNumber+1) + center_
                if diff_[center_] > 400:
                    raise Exception('''didn't find center''')
            else:
                tIndices = np.arange(*[np.argmin(np.abs(tBList[i]-epoch_)) for epoch_ in epochLimList[i]])
                tIndices = tIndices[::int(1/dilution)]
        elif '5VPS' in datasets[0]:
            if not (tIndicesNumber is None):
                tIndices = np.arange(-tIndicesNumber, tIndicesNumber+1) + np.argmin(np.abs(tBList[i]-epochList[i]))
            else:
                tIndices = np.arange(*[np.argmin(np.abs(tBList[i]-epoch_)) for epoch_ in epochLimList[i]])
        tIndicesList[i] = tIndices
        maxHalfWindowLength = (np.diff([np.argmin(np.abs(tBList[i]-epoch_)) for epoch_ in epochLimList[i, :]]))//2+3
        halfWindowLengths = np.arange(int(minHalfWindowLength[i]*1000/tStep), maxHalfWindowLength)
        halfWindowLengths = halfWindowLengths[::int(1/dilutionOfWindow)]
        halfWindowLengthsList[i] = halfWindowLengths.copy()
        ratios = np.zeros([len(halfWindowLengths), len(tIndices)])
        normals = np.zeros([len(halfWindowLengths), len(tIndices), 3])
        numberOfPoints = np.zeros([len(halfWindowLengths), len(tIndices)])
        for indexOfCenter, center_ in enumerate(tIndices):
            for indexOfWindowLength, halfWindowLength_ in enumerate(halfWindowLengths):
                leftLimit = center_-halfWindowLength_
                rightLimit = center_+halfWindowLength_
                tIndices_ = np.arange(leftLimit, rightLimit+1)
                eigenSystem_, ratio_ = mvab(bVectorsList[i][tIndices_], average=average)
                numberOfPoints[indexOfWindowLength, indexOfCenter] = len(tIndices_)
                ratios[indexOfWindowLength, indexOfCenter] = ratio_ * numberOfPoints[indexOfWindowLength, indexOfCenter]
                normals[indexOfWindowLength, indexOfCenter, :] = eigenSystem_[1][:, 0].copy()
    #            if ratioList[i] < ratio_:
    #                ratioList[i] = ratio_
    #                eigenSystemList[i] = eigenSystem_
    #                normalList[i, :] = eigenSystem_[1][:, permutation][:, 0].copy()
    #                centerList[i] = center_
    #                halfWindowLengthList[i] = halfWindowLength_
        ratiosList[i] = ratios.copy()
        normalsList[i] = normals.copy()
        numberOfPointsList[i] = numberOfPoints.copy()
    return ratiosList, numberOfPointsList, normalsList, halfWindowLengthsList, tIndicesList


def showRatioAndAngleDiff(spacecrafts, tBList, tStep, epochLimList, tIndicesList, halfWindowLengthsList, ratiosList, halfWindowLengthList, centerList, ratioList, diffSphAfterConv):
    # plot ratios to see which are the best choice of window length and center
    left = 0.1
    bottom = 0.2
    height = 0.37
    width = 0.18
    wholeWidth = 0.21
    figRatioAndAngleDiff = plt.figure(figsize=(16, 9))
    colorBarLabel = '$\lambda_2/\lambda_3$'
    for i in range(len(spacecrafts)):
        axRatios = figRatioAndAngleDiff.add_axes([left+wholeWidth*i, bottom, width, height])
        axColorBar = figRatioAndAngleDiff.add_axes([left+wholeWidth*i, bottom-0.05, width, 0.01])
        cont = axRatios.pcolormesh(tBList[i][tIndicesList[i]], halfWindowLengthsList[i]*tStep/1000, ratiosList[i]/halfWindowLengthsList[i][:, None]/2, cmap='plasma')
        axRatios.plot([tBList[i][centerList[i]]]*2, halfWindowLengthsList[i][[0, -1]]*tStep/1000, color='k')
        axRatios.plot(tBList[i][tIndicesList[i][[0, -1]]], halfWindowLengthList[[i, i]]*tStep/1000, color='k', lw=1, linestyle='--')
        axRatios.text(0.05, 0.85, 'C{}'.format(i+1), transform=axRatios.transAxes, color='k')
        axRatios.xaxis.set_major_formatter(minFormatter)
        cbar = figRatioAndAngleDiff.colorbar(cont, cax=axColorBar, ax=axRatios, orientation='horizontal')
        cbar.set_label(colorBarLabel+'={:.1f}'.format(ratioList[i]/halfWindowLengthList[i]/2), rotation=0)
    xlabel_ = cdflib.cdfepoch.encode([epochLimList.min(), epochLimList.max()])
    xlabel = 'window center [s] in {}-{}'.format(xlabel_[0][:-4], xlabel_[1][-9:-4])
    ylabel = 'half window length [s]'
    figRatioAndAngleDiff.text(0.45, 0.05, xlabel)
    figRatioAndAngleDiff.text(0.05, bottom+height-0.1, ylabel, rotation='vertical')
    axDirectionList = [None]*len(spacecrafts)
    axColorBar = figRatioAndAngleDiff.add_axes([left+wholeWidth*len(spacecrafts), bottom+height, 0.01, height])
    for i in range(len(spacecrafts)):
        axDirectionList[i] = figRatioAndAngleDiff.add_axes([left+wholeWidth*i, bottom+height, width, height])
        cont = axDirectionList[i].pcolormesh(
                tBList[i][tIndicesList[i]],
                halfWindowLengthsList[i]*tStep/1000,
                np.log(diffSphAfterConv[i]),
#                levels=8,
#                locator=ticker.LogLocator(base=2),
                cmap='gray')
        axDirectionList[i].plot([epochLimList[i, j] for j in range(2)], halfWindowLengthList[[i, i]]*tStep/1000, color='k', lw=1, linestyle='--')
        axDirectionList[i].plot([tBList[i][centerList[i]]]*2, halfWindowLengthsList[i][[0, -1]]*tStep/1000, color='k', lw=0.5)
        axDirectionList[i].text(0.05, 0.85, 'C{}'.format(i+1), transform=axDirectionList[i].transAxes, color='k')
        axDirectionList[i].get_xaxis().set_visible(False)
    cbar = figRatioAndAngleDiff.colorbar(cont, cax=axColorBar, ax=axDirectionList)
    cbar.set_label('~1/$\\theta$')
    return figRatioAndAngleDiff


getMiddle = lambda x: np.take(x, x.size//2)


def fourLongBars(spacecrafts, epochLimList, cdfFileList, datasets, zoomInHalfWindow=20, lowpass=None, sameT=False):
    bMagList = [None]*len(spacecrafts)
    tBList = [None]*len(spacecrafts)
    for i in range(len(spacecrafts)):
        tBList[i] = cdfFileList[i].varget(0)
        if lowpass:
            if '5VPS' in datasets[0]:
                bMagList[i] = butter_lowpass_filter(cdfFileList[i].varget(3), lowpass, 5, axis=0)
            elif 'FULL' in datasets[0]:
                bMagList[i] = butter_lowpass_filter(cdfFileList[i].varget(3), lowpass, 22, axis=0)
        else:
            bMagList[i] = cdfFileList[i].varget(3)
    #zoomInXlim = cdflib.cdfepoch.compute_epoch([[2002, 4, 7, 1, 26, 20, 0],
    #                                            [2002, 4, 7, 1, 27, 5, 0]])
    zoomInTIndicesList = [None]*len(spacecrafts)
    if sameT:
        zoomInXlimList = np.array([[epochList.min(), epochList.max()]]*len(spacecrafts)) + zoomInHalfWindow*1000*np.array([-1, 1])[None, :]
    else:
        zoomInXlimList = epochLimList + zoomInHalfWindow*1000*np.array([-1, 1])[None, :]
    for i in range(len(spacecrafts)):
        zoomInTIndices = np.arange(*includeEnd([np.argmin(np.abs(tBList[i]-epoch_)) for epoch_ in zoomInXlimList[i]]))
        zoomInTIndicesList[i] = zoomInTIndices
    left = 0.1
    bottom = 0.06
    height = 0.18
    width = 0.7
    if sameT:
        wholeHeight = height
    else:
        wholeHeight = height + 0.07
    colors = ["k", 'b', 'r', 'c']
    color = 'g'
    plt.rc('lines', linewidth=0.4)
    figFourLongBars = plt.figure(figsize=(10, 6))
    lw = 0.4
    axZoomInBList = [None]*len(spacecrafts)
    axZoomInDBList = [None]*len(spacecrafts)
    for i in range(len(spacecrafts)):
        axZoomInBList[i] = figFourLongBars.add_axes([left, bottom+i*wholeHeight, width, height])
        axZoomInBList[i].plot(tBList[i][zoomInTIndicesList[i]], bMagList[i][zoomInTIndicesList[i]], color=colors[i], lw=lw, label='C{}'.format(i+1))
    for i in range(len(spacecrafts)):
#        axZoomInBList[i].xaxis.set_ticks(np.linspace(*zoomInXlimList[i], 4))
        axZoomInBList[i].tick_params(which='both', direction='in', top=True, right=True)
        axZoomInBList[i].set_ylabel('B [nT]')
        axZoomInBList[i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axZoomInBList[i].xaxis.set_major_formatter(hourFormatter)
#        axZoomInBList[i].xaxis.set_minor_formatter(minFormatter)
        if sameT:
            if i !=0:
                axZoomInBList[i].tick_params(which='both', labelbottom=False)
            axZoomInBList[0].set_xlabel('time [s] since {}'.format(cdflib.cdfepoch.encode(zoomInXlimList[i][0])[:19]))
        else:
            axZoomInBList[i].set_xlabel('time [s] since {}'.format(cdflib.cdfepoch.encode(zoomInXlimList[i][0])[:19]))
#        axZoomInBList[i].legend(loc='upper left', frameon=False)
    return figFourLongBars

##
def showMVABDataset(spacecrafts, epochLimList, centerList, halfWindowLengthList, cdfFileList, tBList, xList, datasets, lowpass=None, tDtype=np.float64):
    bMagList, = readData(spacecrafts, datasets, cdfFileList, variables=['bMagList'], lowpass=lowpass)
    #zoomInXlim = cdflib.cdfepoch.compute_epoch([[2002, 4, 7, 1, 26, 20, 0],
    #                                            [2002, 4, 7, 1, 27, 5, 0]])
    if 'C' in spacecrafts[0]:
        zoomInHalfWindow = 5*10**3
        myFormatter = minFormatter
    if 'mms' in spacecrafts[0]:
        zoomInHalfWindow = 5*10**9
        myFormatter = minFormatterTT2000
    zoomInTIndicesList = [None]*len(spacecrafts)
    zoomInXlimList = epochLimList + zoomInHalfWindow*np.array([-1, 1])[None, :]
    print(zoomInXlimList)
    dBList = [None]*len(spacecrafts)
    for i in range(len(spacecrafts)):
        zoomInTIndices = np.arange(*includeEnd([np.argmin(np.abs(tBList[i]-epoch_)) for epoch_ in zoomInXlimList[i]]))
        zoomInTIndicesList[i] = zoomInTIndices
        dBList[i] = 1000*np.diff(bMagList[i][np.append(zoomInTIndices, zoomInTIndices.max()+1)]) / np.diff(tBList[i][np.append(zoomInTIndices, zoomInTIndices.max()+1)])
    left = 0.1
    bottom = 0.1
    height = 0.5
    width = 0.19
    wholeWidth = 0.21
    colors = ["k", 'b', 'r', 'c']
    color = 'g'
    plt.rc('lines', linewidth=0.4)
    figZoomIn = plt.figure(figsize=(10, 6))
    lw = 0.4
    axZoomInBList = [None]*len(spacecrafts)
    axZoomInDBList = [None]*len(spacecrafts)
    for i in range(len(spacecrafts)):
        axZoomInBList[i] = figZoomIn.add_axes([left+wholeWidth*i, bottom, width, height])
        axZoomInDBList[i] = figZoomIn.add_axes([left+wholeWidth*i, bottom+height, width, 0.3])
        axZoomInBList[i].plot(tBList[i][zoomInTIndicesList[i]], bMagList[i][zoomInTIndicesList[i]], color=colors[i], lw=lw, label='C{}'.format(i+1))
        axZoomInDBList[i].plot(tBList[i][zoomInTIndicesList[i]], dBList[i], color=colors[i], lw=lw, label='C{}'.format(i+1))
    ylimB_ = np.array([axZoomInBList[i].get_ylim() for i in range(len(spacecrafts))])
    ylimB_ = [ylimB_.flatten().min(), ylimB_.flatten().max()]
    ylimDB_ = axZoomInDBList[0].get_ylim()
    for i in range(len(spacecrafts)):
        axZoomInBList[i].plot([tBList[i][centerList[i]-halfWindowLengthList[i]]]*2, ylimB_, color=colors[i], ls='--', lw=lw)
        axZoomInBList[i].plot([tBList[i][centerList[i]+halfWindowLengthList[i]]]*2, ylimB_, color=colors[i], ls='--', lw=lw)
        axZoomInDBList[i].plot([tBList[i][centerList[i]-halfWindowLengthList[i]]]*2, ylimDB_, color=colors[i], ls='--', lw=lw)
        axZoomInDBList[i].plot([tBList[i][centerList[i]+halfWindowLengthList[i]]]*2, ylimDB_, color=colors[i], ls='--', lw=lw)
        #axZoomInBList[i].set_xlim(zoomInXlimList[i])
        #axZoomInBList[i].set_ylim([0, tBList[0][zoomInXlimList[i]])
        axZoomInBList[i].set_ylim(ylimB_)
        axZoomInDBList[i].set_ylim(ylimDB_)
        axZoomInBList[i].xaxis.set_ticks(np.linspace(*zoomInXlimList[i], 4))
        axZoomInBList[i].tick_params(which='both', direction='in', top=True, right=True, labelleft=False)
        if i == 0:
            axZoomInBList[i].tick_params(which='both', labelleft=True)
            axZoomInBList[i].set_ylabel('B [nT]')
        axZoomInBList[i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axZoomInBList[i].xaxis.set_major_formatter(myFormatter)
#        axZoomInBList[i].xaxis.set_minor_formatter(myFormatter)
        axZoomInBList[i].tick_params(which='both', direction='in', top=True, right=True)
        beginTime = cdflib.cdfepoch.encode(zoomInXlimList[i][0].astype(tDtype))[:19]
        axZoomInBList[i].set_xlabel('time [s]\n since {}'.format(beginTime))
#        axZoomInBList[i].legend(loc='upper left', frameon=False)
        axZoomInDBList[i].xaxis.set_ticks(np.linspace(*zoomInXlimList[i], 4))
        axZoomInDBList[i].tick_params(which='both', direction='in', top=True, right=True, labelbottom=False, labelleft=False)
        if i == 0:
            axZoomInDBList[i].tick_params(which='both', labelleft=True)
            axZoomInDBList[i].set_ylabel('dB/dt [nT/s]')
        axZoomInDBList[i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        axZoomInDBList[i].xaxis.set_major_formatter(myFormatter)
        axZoomInDBList[i].set_title('C{}'.format(i+1))
        #axZoomInDBList[i].get_xaxis().set_visible(False)
    return figZoomIn
#
#    plt.plot(tBList[0][zoomInTIndicesList[0]], bVectorsList[0][zoomInTIndicesList[0],0])
#    plt.plot(tBList[0][zoomInTIndicesList[0]], bVectorsList[0][zoomInTIndicesList[0],1])
#    plt.plot(tBList[0][zoomInTIndicesList[0]], bVectorsList[0][zoomInTIndicesList[0],2])
##

def chooseMVABInterval(spacecrafts, ratiosList, numberOfPointsList, normalsList, halfWindowLengthsList, tIndicesList, epochLimList, minimumRatio=3, percentOfAngle=0.3, show=False, tBList=None, cdfFileList=None, xList=None, percentOfRatio=None, tStep=None, datasets=None, lowpass=None):
    normalList = np.zeros((len(spacecrafts), 3))
    ratioList = np.zeros(len(spacecrafts))
    pointsNumberList = np.zeros(len(spacecrafts))
    centerList = np.zeros(len(spacecrafts), dtype=int)
    halfWindowLengthList = np.zeros(len(spacecrafts), dtype=int)
    diffSph = [None]*len(spacecrafts)
    diffSphAfterConv = [None]*len(spacecrafts)
    if len(normalsList[0].squeeze().shape) == 2:
        angleFilter = np.array([0.05, 0.1, 0.1, 0.2, 0.1, 0.05, 0.05])
        for i in range(len(spacecrafts)):
            diffSph[i] = np.arccos(np.sum(normalsList[i][:-1].squeeze()*normalsList[i][1:].squeeze(), axis=1))
            diffSphAfterConv[i] = nn.convolution(1/diffSph[i][None, : ,None], angleFilter[:, None, None], padding='SAME').numpy().squeeze()
            columnNumber_ = len(tIndicesList[i])
            diffSphAfterConv_ = diffSphAfterConv[i].flatten()
            permutation = np.argsort(diffSphAfterConv_)[::-1]
            number_ = int(len(permutation)*percentOfAngle)
            permutation2 = np.argsort((ratiosList[i]/numberOfPointsList[i]).flatten()[permutation[:number_]])[::-1]
            number_2 = int(number_*percentOfRatio)
            potential = ratiosList[i].flatten()[permutation[:number_][permutation2[:number_2]]]
            potentialArg = np.argsort(potential)[::-1]
            arg = potentialArg[0]
            indexOfCenter = permutation[permutation2[arg]] % columnNumber_
            indexOfWindowLength = permutation[permutation2[arg]] // columnNumber_
            halfWindowLength_ = halfWindowLengthsList[i][indexOfWindowLength]
            center_ = tIndicesList[i][indexOfCenter]
            ratioList[i] = potential[arg]
            pointsNumberList[i] = numberOfPointsList[i][indexOfWindowLength, indexOfCenter]
            normalList[i, :] = normalsList[i][indexOfWindowLength, indexOfCenter]
            centerList[i] = center_
            halfWindowLengthList[i] = halfWindowLength_
        return centerList, halfWindowLengthList, normalList, ratioList, pointsNumberList, None, None
    elif len(normalsList[0].squeeze().shape) == 3:
        directions = ['y', 'x']
        for i in range(len(spacecrafts)):
            diffSph[i] = np.concatenate(2*[0.5*np.ones_like(normalsList[i][..., 0, None])], axis=2)
            for indexOfDirection, direction in enumerate(directions):
                if indexOfDirection == 0:
                    tensorDown = normalsList[i][:-1, ...]
                    tensorUp = normalsList[i][1:, ...]
                    theta = np.arccos(np.abs(np.sum(tensorDown * tensorUp, axis=2)))/np.pi
                    diffSph[i][:-1, :, indexOfDirection] = theta
                if indexOfDirection == 1:
                    tensorLeft = normalsList[i][:, :-1, :]
                    tensorRight = normalsList[i][:, 1:, :]
                    theta = np.arccos(np.abs(np.sum(tensorLeft * tensorRight, axis=2)))/np.pi
                    diffSph[i][:, :-1, indexOfDirection] = theta
        angleFilter_ = np.array([[0.05, 0.15, 0.05],
                                 [0.1, 0.25, 0.1],
                                 [0.1, 0.25, 0.1],
                                 [0.2, 0.4, 0.2],
                                 [0.1, 0.25, 0.1],
                                 [0.05, 0.15, 0.05],
                                 [0.05, 0.15, 0.05]])
        angleFilter = np.repeat(angleFilter_[:, :, None], 2 , axis=2)
        for i in range(len(spacecrafts)):
            diffSphAfterConv[i] = nn.convolution(1/diffSph[i][None, ...], angleFilter[..., None], padding='SAME').numpy().squeeze()
            cornerDrop = []
            columnNumber_ = len(tIndicesList[i])
            leftLimit = epochLimList[i, 0]
            rightLimit = epochLimList[i, 1]
            for indexOfCenter, center_ in enumerate(tIndicesList[i]):
                for indexOfWindowLength, halfWindowLength_ in enumerate(halfWindowLengthsList[i]):
                    if tBList[i][center_ - halfWindowLength_]  < leftLimit or tBList[i][center_ + halfWindowLength_] > rightLimit:
                        cornerDrop.append(indexOfWindowLength*columnNumber_+indexOfCenter)
            diffSphAfterConv_ = diffSphAfterConv[i].flatten()
            diffSphAfterConv_[cornerDrop] = diffSphAfterConv_.min()
            permutation = np.argsort(diffSphAfterConv_)[::-1]
            number_ = int(len(permutation)*percentOfAngle)
    #        _, columnNumber_ = diffSphAfterConv[i].shape
            permutation2 = np.argsort((ratiosList[i]/halfWindowLengthsList[i][:, None]).flatten()[permutation[:number_]])[::-1]
            number_2 = int(number_*percentOfRatio)
            potential = ratiosList[i].flatten()[permutation[:number_][permutation2[:number_2]]]
            potentialArg = np.argsort(potential)[::-1]
            for arg in potentialArg:
                indexOfCenter = permutation[permutation2[arg]] % columnNumber_
                indexOfWindowLength = permutation[permutation2[arg]] // columnNumber_
                halfWindowLength_ = halfWindowLengthsList[i][indexOfWindowLength]
                center_ = tIndicesList[i][indexOfCenter]
                if 0 < indexOfCenter < columnNumber_ and all(indexOfWindowLength != np.array([0])):
                    ratioList[i] = potential[arg]
                    pointsNumberList[i] = numberOfPointsList[i][indexOfWindowLength, indexOfCenter]
                    normalList[i, :] = normalsList[i][indexOfWindowLength, indexOfCenter]
                    centerList[i] = center_
                    halfWindowLengthList[i] = halfWindowLength_
                    break
            else:
                raise Exception("didn't find best normal")

    #        potential = ratiosList[i].flatten()[permutation[:number_]]
    #        potentialArg = np.argsort(potential)
    #        for arg in potentialArg[::-1]:
    #            indexOfCenter = permutation[arg] % columnNumber_
    #            indexOfWindowLength = permutation[arg] // columnNumber_
    #            halfWindowLength_ = halfWindowLengthsList[i][indexOfWindowLength]
    #            # if potential[arg] > 2 and 0 < indexOfCenter < columnNumber_ and all(indexOfWindowLength !=np.array([0, 1])):
    #            if potential[arg] / (halfWindowLength_*2) > minimumRatio and 0 < indexOfCenter < columnNumber_ and all(indexOfWindowLength != np.array([0, 1])):
    #                ratioList[i] = potential[arg]
    #                normalList[i, :] = normalsList[i][indexOfWindowLength, indexOfCenter]
    #                centerList[i] = tIndicesList[i][indexOfCenter]
    #                halfWindowLengthList[i] = halfWindowLength_
    #                break
    #        else:
    #            raise Exception('ratio too small')
        if show:
            figRatioAndAngleDiff = showRatioAndAngleDiff(spacecrafts, tBList, tStep, epochLimList, tIndicesList, halfWindowLengthsList, ratiosList, halfWindowLengthList, centerList, ratioList, diffSphAfterConv)
            figEventLimit = showMVABDataset(spacecrafts, epochLimList, centerList, halfWindowLengthList, cdfFileList, tBList, xList, datasets, lowpass)
            return centerList, halfWindowLengthList, normalList, ratioList, pointsNumberList, figRatioAndAngleDiff, figEventLimit
        else:
            return centerList, halfWindowLengthList, normalList, ratioList, pointsNumberList, None, None


def computeGradBs(bVectorLists, xGSEs, order=1):
    '''see doi:10.1029/2002JA009612 Appendix B.
    bVectorList and xGSE in the form of [time index, spacecraft index, cartesian index]'''
    numberOfSpacecrafts = xGSEs.shape[1]
    x = xGSEs - np.mean(xGSEs, axis=1)[:, None, :]
    R = np.transpose(x, (0, 2, 1)) @ x / numberOfSpacecrafts
    RInverse = np.linalg.inv(R)
    G0 = np.transpose(bVectorLists, (0, 2, 1)) @ x @ RInverse / numberOfSpacecrafts
    if order == 1:
        LagrangianMultiplier = -np.trace(G0, axis1=1, axis2=2)/np.trace(RInverse, axis1=1, axis2=2)
        G = G0 + LagrangianMultiplier[:, None, None] * RInverse
    else:
        G = G0
    return G


def normalFromPressureGradient(bVectorLists, xGSEs):
    'see doi:10.1029/2002JA009612 equation (3)'
    gradBs = computeGradBs(bVectorLists, xGSEs)
    normals_ = np.mean(bVectorLists, axis=1)[:, None, :] @ gradBs
    normals = -normalized(normals_.squeeze())
    return normals


def mca(bVectorLists=None, xGSEs=None, gradBs=None, bVectorAtCenters=None):
    '''see doi:10.1029/2002JA009612 equation (1)'''
    if gradBs is None:
        gradBs = computeGradBs(bVectorLists, xGSEs)
    if bVectorAtCenters is None:
        bVectorAtCenters = np.mean(bVectorLists, axis=1)
    bMagAtCenters = np.linalg.norm(bVectorAtCenters, axis=-1)
    term1 = np.sum(bVectorAtCenters[:, None, :] * gradBs, axis=2) / bMagAtCenters[:, None]**2
    term2 = - np.sum(np.sum(bVectorAtCenters[:, None, :] * gradBs, axis=2)* bVectorAtCenters, axis=1)[:, None] * bVectorAtCenters / bMagAtCenters[:, None]**4
    curvatureVectors = term1 + term2
    curvatures = np.linalg.norm(curvatureVectors, axis=1)
    normals = curvatureVectors / curvatures[:, None]
    return curvatures, normals


def dipoleField(xGSE, M=-30438):
    'xGSE in km, B in nT or xGSE in m, B in T'
    x1 = xGSE[..., 0][..., None]
    x2 = xGSE[..., 1][..., None]
    x3 = xGSE[..., 2][..., None]
    r = np.linalg.norm(xGSE, axis=-1)[..., None]
    return M*np.concatenate([3*x1*x3, 3*x2*x3, (3*x3**2-r**2)], axis=-1) / r**5


def magnetosphericField(xGSE, M1=-30438, M2=-28*30438):
    'xGSE in km, B in nT or xGSE in m, B in T'
    x1 = xGSE[..., 0][..., None]
    x2 = xGSE[..., 1][..., None]
    x3 = xGSE[..., 2][..., None]
    x4 = xGSE[..., 0][..., None] - 40
    r1 = np.linalg.norm(xGSE, axis=-1)[..., None]
    xGSE2 = xGSE.copy()
    xGSE2[..., 0] = x4.squeeze()
    r2 = np.linalg.norm(xGSE2, axis=-1)[..., None]
    b1 = M1*np.concatenate([3*x1*x3, 3*x2*x3, (3*x3**2-r1**2)], axis=-1) / r1**5
    b2 = M2*np.concatenate([3*x4*x3, 3*x2*x3, (3*x3**2-r2**2)], axis=-1) / r2**5
    return b1+b2

##
def lundquistForceFreeField(x=None, R=None, B0=1, coordinateSystem='Cartesian'):
    if not (x is None):
        r = np.linalg.norm(x[..., 0:2], axis=-1)
        theta = np.arctan(x[..., 1]/x[..., 0])
        theta = (np.sign(x[..., 0]) + 1)/2*theta - (1 - np.sign(x[..., 0]))/2*np.pi*np.sign(x[..., 1])
    bTheta = B0*scipy.special.j1(r)
    bZ = B0*scipy.special.j0(r)
    if coordinateSystem =='Cartesian':
        b = np.stack([-bTheta*np.sin(theta), bTheta*np.cos(theta), bZ], axis=-1)
    return b
##

def harrisBField(xGSE, B0=1, h=1):
    x1 = xGSE[..., 0][..., None]
    x2 = xGSE[..., 1][..., None]
    x3 = xGSE[..., 2][..., None]
    return np.concatenate([B0*np.tanh(x3/h), np.zeros_like(x3), np.zeros_like(x3)], axis=-1)

##
def list2arrayAccordingToEpochs(epochs, tBList, listOfLists, allowFailure=False):
    numberOfPoints = len(epochs)
    numberOfSpacecrafts = len(tBList)
    centerLists = np.zeros((numberOfPoints, numberOfSpacecrafts), dtype=int)
    returnedVariables = [None]*len(listOfLists)
    tStep = epochs[1]-epochs[0]
    epochsNotMatched = []
    for k, lists in enumerate(listOfLists):
        returnedVariables[k] = np.zeros((numberOfPoints, numberOfSpacecrafts, *lists[0].shape[1:]), dtype=lists[0].dtype)
    for i in range(numberOfSpacecrafts):
        centerLists[0, i] = np.argmin(np.abs(tBList[i]-epochs[0]))
        for k in range(1, numberOfPoints):
            lastCenter = centerLists[k-1, i]
            rightLim = lastCenter+4
            diff = np.abs(tBList[i][lastCenter:rightLim] - epochs[k])
            forwardSteps = np.argmin(diff)
            if forwardSteps > 2:
                raise Exception('''step too big''')
            if diff[forwardSteps]/tStep > 0.5:
                message = '''lost a point:\n spacecraftInd: {}\n epochInd: {}\n epoch: {}\n minimum diff={} \n lastCenter: {}\n forwardSteps: {}\n gap: {}'''.format(i, k, epochs[k], diff[forwardSteps], lastCenter, forwardSteps, tBList[i][lastCenter+forwardSteps+1] - tBList[i][lastCenter+forwardSteps])
                if not allowFailure:
                    raise Exception(message)
                else:
                    epochsNotMatched.append(epochs[k])
                    print(message)
            centerLists[k, i] = lastCenter  + forwardSteps
        for k, array in enumerate(returnedVariables):
            array[:, i] = listOfLists[k][i][centerLists[:, i]]
    if allowFailure:
        returnedVariables.append(epochsNotMatched)
    return returnedVariables
##

def pressureGradientOfTimeSeriesB(tBList, bVectorsList, xList, epochs):
    bVectorLists, xGSEs = list2arrayAccordingToEpochs(epochs, tBList, [bVectorsList, xList])
#    numberOfPoints = len(epochs)
#    numberOfSpacecrafts = len(tBList)
#    centerLists = np.zeros((numberOfPoints, numberOfSpacecrafts), dtype=int)
#    bVectorLists = np.zeros((numberOfPoints, numberOfSpacecrafts, 3))
#    xGSEs = np.zeros((numberOfPoints, numberOfSpacecrafts, 3))
#    for i in range(numberOfSpacecrafts):
#        centerLists[0, i] = np.argmin(np.abs(tBList[i]-epochs[0]))
#        for k in range(1, numberOfPoints):
#            lastCenter = centerLists[k-1, i]
#            rightLim = lastCenter+5
#            forwardSteps = np.argmin(np.abs(tBList[i][lastCenter:rightLim] - epochs[k]))
#            if forwardSteps > 3:
#                raise Exception('''step too big''')
#            centerLists[k, i] = lastCenter  + forwardSteps
#        bVectorLists[:, i] = bVectorsList[i][centerLists[:, i]]
#        xGSEs[:, i] = xList[i][centerLists[:, i]]
    normals = normalFromPressureGradient(bVectorLists, xGSEs)
    return normals


def timing(normalList, xGSE, getShape=False):
    x = xGSE - np.mean(xGSE, axis=0)
    R = x.T @ x / 4
    volume = 8/3*np.sqrt(np.linalg.det(R))
    RInverse = np.linalg.inv(R)
    nablaN = np.empty((3, 3))
    combs = combinations(range(4), 2)
    number_ = np.math.factorial(4)//np.math.factorial(2)//np.math.factorial(2)
    diffX = np.empty((number_, 3))
    if len(normalList.shape) == 1:
        diffN = np.empty(number_)
    elif len(normalList.shape) == 2:
        diffN = np.empty((number_, 3))
    for i, comb in enumerate(list(combs)):
        diffX[i] = x[comb[0]] - x[comb[1]]
        diffN[i] = normalList[comb[0]] - normalList[comb[1]]
    if len(diffN.shape) == 1:
        diffNTranspose = diffN[None, :]
    elif len(diffN.shape) == 2:
        diffNTranspose = diffN.T
    nablaN = 1/16 * diffNTranspose @ diffX @ RInverse
    if getShape:
        eigenSystemOfR = np.linalg.eig(R)
        permutation = np.argsort(eigenSystemOfR[0])
        timingShape = (np.sqrt(eigenSystemOfR[0])[permutation], eigenSystemOfR[1][:, permutation])
        returnedVariables = [nablaN, timingShape]
        return returnedVariables
    else:
        return nablaN


def calculateTimingShape(xGSE):
    x = xGSE - np.mean(xGSE, axis=0)
    R = x.T @ x / 4
    eigenSystemOfR = np.linalg.eig(R)
    permutation = np.argsort(eigenSystemOfR[0])
    timingShape = (np.sqrt(eigenSystemOfR[0])[permutation], eigenSystemOfR[1][:, permutation])
    return timingShape


def timingVelocityAndNormal(t, xGSE, silence=False, getShape=False):
    m = timing(t, xGSE, getShape=getShape)
    if getShape:
        timingShape = m[1]
        m = m[0]
    timingVelocity = 1/np.linalg.norm(m)
    timingNormal = m*timingVelocity
    timingNormal = timingNormal.squeeze()
    if silence is False:
        print("timing velocity: {:.1f}km/s,\n normal: {}".format(timingVelocity*6371, timingNormal))
    returnedVariables = [timingVelocity, timingNormal]
    if getShape:
        returnedVariables.append(timingShape)
    return returnedVariables


def nfa(normalList, xGSE, projection=True, noProjection=True):
    x = xGSE - np.mean(xGSE, axis=0)
    normalCenter = normalized(np.mean(normalList, axis=0))
    meanDifferenceOfNormals = 2*np.mean(np.arccos(normalList @ normalCenter[:, None]))
    nablaN = timing(normalList, x)
    if projection:
        nablaNTilde = nablaN - nablaN @ normalCenter[:, np.newaxis] @ normalCenter[None, :]
    else:
        nablaNTilde = nablaN
    eigenSystem = np.linalg.eig(nablaNTilde)
    permutation_ = np.argsort(np.abs(eigenSystem[0]))[::-1]
    eigenValues = eigenSystem[0][permutation_]
    eigenVectors = eigenSystem[1][:, permutation_]
    return normalCenter, meanDifferenceOfNormals, eigenValues, eigenVectors, nablaNTilde


def showNFA(xGSE, normalList, eigenValues, eigenVectors, timingNormal):
    meanPosition = np.mean(xGSE, axis=0)
    x = xGSE - meanPosition
    fig = plt.figure(figsize=(6, 6))
    ax3d = fig.add_subplot(111, projection='3d')
    colors = ["k", 'b', 'r', 'c']
    for i, c in enumerate(colors):
        ax3d.scatter(*[xGSE[i, j] for j in range(3)], c=c)
    l = np.diff(ax3d.get_xlim())/10
    for i, c in enumerate(colors):
        ax3d.quiver(*[xGSE[i, j] for j in range(3)], *[normalList[i, j] for j in range(3)], colors=c, length=l, normalize=True)
    ax3d.scatter(*[meanPosition[j] for j in range(3)], c='m')
    normalCenter = normalized(np.mean(normalList, axis=0))
    ax3d.quiver(*[meanPosition[j] for j in range(3)], *[normalCenter[j] for j in range(3)], colors='m', linestyle='-.', linewidths=2, length=2*l, normalize=True)
    toEarth = -normalized(meanPosition)
    ax3d.quiver(*[meanPosition[j] for j in range(3)], *[toEarth[j] for j in range(3)], colors='y', linestyle='-', linewidths=2, length=3*l, normalize=True)
    ax3d.quiver(*[meanPosition[j] for j in range(3)], *[timingNormal[j] for j in range(3)], colors='g', linestyle='-', linewidths=2, length=3*l, normalize=True)
    lengths = [2*l, 4*l]
    for i in range(2):
        ax3d.quiver(*[meanPosition[j] for j in range(3)], *[eigenVectors[j, i] for j in range(3)], colors='m', length=lengths[i], linewidths=3, normalize=True)
    ax3d.set_xlabel('x [$R_E$]')
    ax3d.set_ylabel('y [$R_E$]')
    ax3d.set_zlabel('z [$R_E$]')
    extent = np.mean(np.linalg.norm(x, axis=1))
    angle = np.arccos(eigenVectors[:, 0] @ eigenVectors[:, 1])
    coeff_ = np.array([-np.cos(angle), 1])
    vectorOrthogonalToV1 = (normalized(coeff_[None, :] @ eigenVectors[:, :2].T)).squeeze()
    numberOfPointsInQuadrant = 9
    endPoints = np.zeros((numberOfPointsInQuadrant*4+1, 3, 2))
    resolutionInQuadrants = np.array([angle, np.pi-angle]*2)/numberOfPointsInQuadrant
    for i in range(4*numberOfPointsInQuadrant):
        quadrant = i // numberOfPointsInQuadrant
        theta = quadrant*np.pi/2 + i % numberOfPointsInQuadrant * resolutionInQuadrants[quadrant]
        normalCurvature = np.array([np.cos(theta), np.sin(theta)])**2 @ eigenValues[:2] 
        direction = np.cos(theta)*eigenVectors[:, 0] + np.sin(theta)*vectorOrthogonalToV1
        angle_ = extent * normalCurvature
        endPoints[i, :, 1] = extent * (-np.sin(angle_)*normalCenter + np.cos(angle_)*direction)
        endPoints[-1, :, 1] = endPoints[0, :, 1]
    endPoints = endPoints + meanPosition[None, :, None]
    ax3d.plot_surface(*[endPoints[:, i, :] for i in range(3)], color=(0, 0, 0, 0.2))
    return fig


def normalListCalibrate(normalList, ratioList, timingNormal, halfWindowLengthList, centerList, bVectorsList, average):
    criterion = np.cos(45/180*np.pi)
    for i in range(normalList.shape[0]):
        if np.abs(normalList[i] @ timingNormal) < criterion:
            center_ = centerList[i]
            halfWindowLength_ = halfWindowLengthList[i]
            leftLimit = center_-halfWindowLength_
            rightLimit = center_+halfWindowLength_
            tIndices_ = np.arange(leftLimit, rightLimit+1)
            eigenSystem_, ratio_ = mvab(bVectorsList[i][tIndices_], average=average)
            normalList[i] = eigenSystem_[1][:, 1].copy()
            ratioList[i] = 1/ratio_*len(tIndices_)
    quality_ = transNormal(normalList)


def presentCurvature(normalList, xGSE, ratioList, pointsNumberList, silence=True, projection=True):
    errorOfNormals = 1/np.sqrt(np.mean(ratioList))
    np.set_printoptions(precision=3)
    massOfCenter = np.mean(xGSE, axis=0)
    x = xGSE - massOfCenter
    normalCenter, meanDifferenceOfNormals, eigenValues, eigenVectors, nablaNTilde = nfa(normalList, xGSE, projection=projection)
    ratio = ratioList/pointsNumberList
    e1e2 = float((eigenVectors[:, 0] @ eigenVectors[:, 1]))
    try:
        e1ne2n = [float((eigenVectors[:, :2].T @ normalCenter[:, None])[i]) for i in range(2)]
        quality = True
    except:
        e1ne2n = [float('Nan')]*2
        quality = False
    dictOfEvent = {'quality': quality,
                   'normalList': normalList.tolist(),
                   'numberOfPoints': pointsNumberList.tolist(),
                   'ratio': ratio.tolist(),
                   'errorOfNormals': errorOfNormals,
                   'meanDifferenceOfNormals': meanDifferenceOfNormals,
                   'positionsInMassOfCenterRef': x.tolist(),
                   'massOfCenter': massOfCenter.tolist(),
                   'principalCurvature': eigenValues.tolist(),
                   'principalDirection': eigenVectors.tolist(),
                   'overallError': [e1e2, *e1ne2n]}
    if silence is False:
        print('quality: {}'.format(quality))
        print('normal at four points:\n {}'.format(normalList))
        print('number of points: {}\nlambda2/lambda3: {}\nerror of normals: {:.3f}\nmean difference: {:.3f}'.format(pointsNumberList, ratio, errorOfNormals, meanDifferenceOfNormals))
        print('constellation position:\n {} RE'.format(xGSE[0]))
        print('constellation seperation (in km):\n {}'.format(x*6371))
        print('principal curvature in 1/RE:\n {}'.format(eigenValues))
        print('principal direction:\n {}'.format(eigenVectors))
#        if '!' in  e1ne2n:
#            print('$e_1 \cdot e_2$ $e_1 \cdot n$ $e_2 \cdot n$\n'+'{:.2f} {} {}'.format(e1e2, *e1ne2n))
#        else:
        print('$e_1 \cdot e_2$ $e_1 \cdot n$ $e_2 \cdot n$\n'+'{:.2f} {:.2f} {:.2f}'.format(e1e2, *e1ne2n))
        if projection:
            print('nablaNTilde')
        else:
            print('nablaN')
        print(nablaNTilde)
    return eigenValues, eigenVectors, dictOfEvent

##
def mainCalculation(workDataDir, spacecrafts, datasets, event=None, eventInd=None, lowpass=None, minHalfWindowLength=None, percentOfAngle=0.07, percentOfRatio=0.4, tIndicesNumber=0, dilution=1, dilutionOfWindow=1, silence=False, projection=True, show=True, saveFig=False, directoryFigRatioAndAngleDiff='figRatioAndAngleDiff', directoryFigEventLimit='figEventLimit', directoryFigNFA='figNFA', wlSession=None, mvabRange=None, epochList=None, coordinateSystem='gse'):
    eventNew = {}
    if 'C' in spacecrafts[0]:
        tDtype = np.float64
    elif 'mms' in spacecrafts[0]:
        tDtype = np.int64
    if np.any(mvabRange):
        epochLimList = mvabRange
        if epochList is None:
            epochList = np.mean(epochLimList, axis=1, dtype=tDtype)
        cdfFileList = defineFiles(workDataDir, spacecrafts, datasets, epochList[0])
        tBList, bVectorsList, xList = readData(spacecrafts, datasets, cdfFileList,  epochList, lowpass, variables=['tBList', 'bVectorsList', 'xList'], coordinateSystem=coordinateSystem)
        normalList = np.zeros((4, 3))
        ratioList = np.zeros(4)
        halfWindowLengthList = np.zeros(4, dtype=int)
        centerList = np.zeros(4, dtype=int)
        pointsNumberList = np.zeros(4)
        epochLimListOfData = np.zeros_like(epochLimList, dtype=tDtype)
        for i in range(4):
            tIndices_ = np.arange(*includeEnd([np.argmin(np.abs(tBList[i]-epoch_)) for epoch_ in epochLimList[i]]))
            epochLimListOfData[i] = tBList[i][tIndices_[[0, -1]]]
            eigenSystem_, ratio_ = mvab(bVectorsList[i][tIndices_])
            pointsNumberList[i] = len(tIndices_)
            centerList[i] = int(getMiddle(tIndices_))
            halfWindowLengthList[i] = pointsNumberList[i]//2
            ratioList[i] = ratio_ * pointsNumberList[i]
            normalList[i] = eigenSystem_[1][:, 0]
        if not silence:
            print([cdflib.cdfepoch.encode(epochLimListOfData[i]) for i in range(4)])
        if show:
            figEventLimit = showMVABDataset(spacecrafts, epochLimList, centerList, halfWindowLengthList, cdfFileList, tBList, xList, datasets, lowpass, tDtype=tDtype)
        else:
            figEventLimit = figRatioAndAngleDiff = None
        #    pprint(eigenSystem_)
        #    print(ratio_)
        #    print(cdflib.cdfepoch.encode(epochLimList[i]))
#        quality_ = transNormal(normalList)
#        eigenValues, eigenVectors, eventInfoNew = presentCurvature(normalList, xGSE, ratioList, pointsNumberList, silence=silence, projection=projection)
    else:
        if eventInd is None:
            eventInd = event['eventInd']
        epochList = np.array(event['epochs']).squeeze()
        epochLimList = np.array(event['epochLimList']).squeeze()
        average = 1
        minimumRatio = 5
        cdfFileList = defineFiles(workDataDir, spacecrafts, datasets, epochList[0])
        tBList, tStep, bVectorsList, xList = readData(spacecrafts, datasets, cdfFileList,  epochList, lowpass)
        ratiosList, numberOfPointsList, normalsList, halfWindowLengthsList, tIndicesList = computeNormalsListAndRatiosList(tBList, tStep, epochList, epochLimList, bVectorsList, spacecrafts, datasets, lowpass=lowpass, minHalfWindowLength=minHalfWindowLength, dilution=dilution, tIndicesNumber=tIndicesNumber, dilutionOfWindow=dilutionOfWindow, average=average)
        centerList, halfWindowLengthList, normalList, ratioList, pointsNumberList, figRatioAndAngleDiff, figEventLimit = chooseMVABInterval(spacecrafts, ratiosList, numberOfPointsList, normalsList, halfWindowLengthsList, tIndicesList, epochLimList, minimumRatio, percentOfAngle, show=show, tBList=tBList, cdfFileList=cdfFileList, xList=xList, percentOfRatio=percentOfRatio, tStep=tStep, datasets=datasets)
    xGSE = np.empty((4, 3))
    t = np.empty(4)
    for i in range(4):
        xGSE[i] = xList[i][centerList[i]]/6371
        if 'C' in spacecrafts[i]:
            t[i] = tBList[i][centerList[i]]/1000
        if 'mms' in spacecrafts[i]:
            t[i] = tBList[i][centerList[i]]/10**9
    timingVelocity, timingNormal = timingVelocityAndNormal(t, xGSE, silence=silence)
    eventNew.update({'normal': timingNormal.tolist(),
              'velocity': timingVelocity})
    if 'C' in spacecrafts[0]:
        normalListCalibrate(normalList, ratioList, timingNormal, halfWindowLengthList, centerList, bVectorsList, average=average)
    elif 'mms' in spacecrafts[0]:
        quality_ = transNormal(normalList)
    eigenValues, eigenVectors, eventInfoNew = presentCurvature(normalList, xGSE, ratioList, pointsNumberList, silence=silence, projection=projection)
    eventNew.update(eventInfoNew)
    normalCenter = normalized(np.mean(normalList, axis=0))
    #_ = align(timingNormal[None, :], normalList)
    #timingNormal = timingNormal.squeeze()
    if not show:
        figNFA = None
    else:
        figNFA = showNFA(xGSE, normalList, eigenValues, eigenVectors, timingNormal)
        if saveFig:
            figBaseName = '{}-{}.pdf'.format(eventInd, cdflib.cdfepoch.encode(epochList[0])[:-7])
            if not os.path.exists(directoryFigRatioAndAngleDiff):
                os.makedirs(directoryFigRatioAndAngleDiff)
            if not os.path.exists(directoryFigEventLimit):
                os.makedirs(directoryFigEventLimit)
            if not os.path.exists(directoryFigNFA):
                os.makedirs(directoryFigNFA)
            pathFigRatioAndAngleDiff = os.path.join(directoryFigRatioAndAngleDiff, 'ratioAndAngleDiff-'+figBaseName)
            pathFigEventLimit = os.path.join(directoryFigEventLimit, 'eventLimit-'+figBaseName)
            pathFigNFA = os.path.join(directoryFigNFA, 'nfa-'+figBaseName)
            figRatioAndAngleDiff.savefig(pathFigRatioAndAngleDiff)
            figEventLimit.savefig(pathFigEventLimit)
            figNFA.savefig(pathFigNFA)
            eventNew.update({'path:figRatioAndAngleDiff': pathFigRatioAndAngleDiff, 'path:figEventLimit': pathFigEventLimit, 'path:figNFA': pathFigNFA})
    xGSECenter = np.mean(xGSE, axis=0)
    if wlSession:
        modelCurvatures, modelDirections = modelCurvaturesAndDirections(wlSession, xGSECenter)
        for i in range(4):
            modelCurvatures, modelDirections = modelCurvaturesAndDirections(wlSession, xGSE[i])
            normalList[i] = modelDirections[:, 2]
        normalCenter, meanDifferenceOfNormals, eigenValues, eigenVectors, nablaNTilde = nfa(normalList, xGSE, projection=projection)
        e1e2 = float((eigenVectors[:, 0] @ eigenVectors[:, 1]))
        try:
            e1ne2n = [float((eigenVectors[:, :2].T @ normalCenter[:, None])[i]) for i in range(2)]
            quality = True
        except:
            e1ne2n = [float('Nan')]*2
            quality = False
        eventNew.update({'J05Curvature': modelCurvatures.tolist(),
                         'J05Direction': modelDirections.tolist(),
                         'J05NFACurvature': eigenValues.tolist(),
                         'J05NFADirection': eigenVectors.tolist(),
                         'J05NFAError': [e1e2, *e1ne2n]})
    if np.any(mvabRange):
        return eventNew, figEventLimit, figNFA
    else:
        return eventNew, figRatioAndAngleDiff, figEventLimit, figNFA
##

def eventSummary(event, file=sys.stdout):
    MVABItems = ['normalList', 'numberOfPoints', 'ratio', 'positionsInMassOfCenterRef']
#    nfaItems = ['errorOfNormals', 'meanDifferenceOfNormals', 'principalCurvature', 'principalDirection', 'overallError']
    itemWidth = np.max([len(item) for item in MVABItems])
    firstColumnWidth = 15
    nameColumnWidth = 25
    width = firstColumnWidth + nameColumnWidth*4
    header_fmt = lambda x:('{{:{}}}'+'{{:^{}}}'*x).format(firstColumnWidth, *[nameColumnWidth]*x)
    content_fmt = lambda x:('{{:{}}}'+'{{:<{}}}'*x).format(firstColumnWidth, *[nameColumnWidth]*x)
    print('='*width, file=file)
    print('-'*width, file=file)
    print('eventInd: {}'.format(event['eventInd']), file=file)
    threeItemsF = '[' +'{{:.{d}f}}, '*2 + '{{:.{d}f}}' + ']'
    complexFormatter =  lambda c, d:'{{c.real:.{d}f}}+{{c.imag:.{d}f}}j'.format(d=d).format(c=c)
    threeItemsC = lambda x, d: '[' + ', '.join([complexFormatter(x[i], d) for i in range(3)]) + ']'
    print('Constellation position [RE]: '+threeItemsF.format(d=2).format(*event['massOfCenter']), file=file)
    print('Event time: {}'.format(cdflib.cdfepoch.encode(event['epochs'][0])), file=file)
    print('Timing velocity [km/s]: {:.1f}'.format(event['velocity']*6371*1000), file=file)
    print('Timing normal [in GSE]: '+threeItemsF.format(d=3).format(*event['normal']), file=file)
    print('='*width, file=file)
    print(header_fmt(4).format('', *['C{}'.format(i) for i in range(1,5)]), file=file)
    print('-'*width, file=file)
    for item in MVABItems:
        if isinstance(event[item][0], list):
            vectors = [threeItemsF.format(d=3).format(*event[item][i]) for i in range(4)]
            if item ==  'positionsInMassOfCenterRef':
                print(content_fmt(4).format('x, y, z [RE]', *vectors), file=file)
            else:
                print(content_fmt(4).format(item, *vectors), file=file)
        else:
            items_ = ['{:.1f}'.format(event[item][i]) for i in range(4)]
            print(content_fmt(4).format(item, *items_), file=file)
    print('='*width, file=file)
    print('{}: {:.3f},  {}: {:.3f}'.format('errorOfNormals', event['errorOfNormals'], 'meanDifferenceOfNormals', event['meanDifferenceOfNormals']), file=file)
    print('semiaxes: ' + threeItemsF.format(d=3).format(*event['timingShape'][0]), file=file)
    coordinates = ['x', 'y', 'z']
    for i, x in enumerate(coordinates):
        if i == 0:
            stringList = ['timingShape:'+x]
        else:
            stringList = [' '*len('timingShape:')+x]
        stringList.append(threeItemsF.format(d=3).format(*event['timingShape'][1][i]))
        print(' '.join(stringList), file=file)
    print('='*width, file=file)
    print(header_fmt(3).format('', 'nfa+data', 'nfa+J05', 'J05'), file=file)
    print('-'*width, file=file)
    methodsMeta = ['principal', 'J05NFA', 'J05']
    stringList = ['curvature']
    for methodInd, method in enumerate(methodsMeta):
        method = method + 'Curvature'
        if isinstance(event[method][0], float):
            stringList.append(threeItemsF.format(d=3).format(*event[method]))
        elif isinstance(event[method][0], complex):
            stringList.append(threeItemsC(x=event[method], d=3))
    print(content_fmt(3).format(*stringList), file=file)
    for i, x in enumerate(coordinates):
        methods = ['principalDirection', 'J05NFADirection', 'J05Direction']
        if i == 0:
            stringList = ['direction:'+x]
        else:
            stringList = [' '*len('direction:')+x]
        for methodInd, method in enumerate(methods):
            if isinstance(event[method][i][0], float):
                stringList.append(threeItemsF.format(d=3).format(*event[method][i]))
            elif isinstance(event[method][i][0], complex):
                stringList.append(threeItemsC(x=event[method][i], d=3))
        print(content_fmt(3).format(*stringList), file=file)
    print(content_fmt(3).format('error', threeItemsF.format(d=3).format(*event['overallError']), threeItemsF.format(d=3).format(*event['J05NFAError']), ' '), file=file)
    print('-'*width, file=file)
    print('='*width, file=file)
#    for item in nfaItems:
##        if isinstance(event[item], float):
##            print('{}: {:.3f}'.format(item, event[item]), file=file)
#        elif isinstance(event[item][0], float):
#            print('{}: '.format(item) + threeItemsF.format(d=3).format(*event[item]), file=file)
#        elif isinstance(event[item][0], complex):
#            print('{}: '.format(item) + threeItemsC(x=event[item], d=3), file=file)
#        elif item == 'principalDirection':
#            coordinates = ['x', 'y', 'z']
#            print('{}: '.format(item), file=file)
#            for i, x in enumerate(coordinates):
#                if isinstance(event[item][i][0], float):
#                    print(' '*(len(item)+1) + x + threeItemsF.format(d=3).format(*event[item][i]), file=file)
#                elif isinstance(event[item][i][0], complex):
#                    print(' '*(len(item)+1) + x + threeItemsC(x=event[item][i], d=3), file=file)
#    print('-'*width, file=file)
#    print('='*width, file=file)
#
##
def calculateAllEvents(workDataDir, spacecrafts, datasets, dictOfEvents, eventsInds=None, MVABRange=None, lowpass=None, tIndicesNumber=3, writeSummary=True, summaryFile=None, directoryForSaving=os.getcwd(), wlSession=None, allowFailure=False, show=True, saveFig=True):
    lost = []
    if saveFig or writeSummary:
        if not os.path.exists(directoryForSaving):
            os.makedirs(directoryForSaving)
        if saveFig:
            directoryFigRatioAndAngleDiff = os.path.join(directoryForSaving, 'figRatioAndAngleDiff')
            directoryFigEventLimit = os.path.join(directoryForSaving, 'figEventLimit')
            directoryFigNFA = os.path.join(directoryForSaving, 'figNFA')
    if (summaryFile is None) and writeSummary:
        localSummary = True
        summaryFile = open(os.path.join(directoryForSaving, 'summary.txt'), 'a')
    else:
        localSummary = False
    if eventsInds is None:
        eventsInds = list(dictOfEvents.keys())
    for eventInd_ in eventsInds:
        eventInd = str(eventInd_)
        print(eventInd)
        event = dictOfEvents[eventInd]
        epochList = np.array(event['epochs'])
        if MVABRange is None:
            maxHalfWindowLength = np.array([60]*4)
            epochLimList = epochList[:, None] + maxHalfWindowLength[:, None]*1000*np.array([-1, 1])[None, :]
        event.update({'epochLimList': epochLimList.tolist()})
        minHalfWindowLength = np.array([40]*4)
        try:
            eventNew, figRatioAndAngleDiff, figEventLimit, figNFA = mainCalculation(workDataDir, spacecrafts, datasets, event, eventInd=eventInd, lowpass=lowpass, minHalfWindowLength=minHalfWindowLength, percentOfAngle=0.07, percentOfRatio=0.4, tIndicesNumber=tIndicesNumber, dilution=1, dilutionOfWindow=0.5, silence=True, projection=True, show=show, saveFig=saveFig, directoryFigRatioAndAngleDiff=directoryFigRatioAndAngleDiff, directoryFigEventLimit=directoryFigEventLimit, directoryFigNFA=directoryFigNFA, wlSession=wlSession)
            plt.close('all')
            event.update(eventNew)
            if writeSummary:
                eventSummary(event, summaryFile)
        except:
            if allowFailure:
                lost.append(eventInd_)
            else:
                raise Exception
    if localSummary:
        summaryFile.close()
    return lost
##
##
def initModel(session, coeff = 'coeff = {a11 -> 0.45, a22 -> 1, a33 -> 0.8, a12 -> 0.18, a14 -> 46.6, a24 -> -2.2, a34 -> -0.6, a44 -> -618};'):
    '''J05 model by default'''
    initModelCMD = '''surface = 
  a11*(cont*x)^2 + a22*(cont*y)^2 + a33*(cont*z)^2 + 
   a12*(cont*x)*(cont*y) + a14*(cont*x) + a24*(cont*y) + 
   a34*(cont*z) + a44;
contraction = Solve[surface == 0, cont];
zsInXY = Solve[surface == 0, z];
r = {x, y, z} /. zsInXY;
rx = D[r, x];
rxx = D[rx, x];
rxy = D[rx, y];
ry = D[r, y];
ryy = D[ry, y];
EE = rx[[1]].rx[[1]];
FF = rx[[1]].ry[[1]];
GG = ry[[1]].ry[[1]];
n = FullSimplify[
    Normalize[
     Grad[surface, {x, y, z}]], {x, y, z, cont, a11, a22, a33, a12, 
      a14, a24, a34, a44} \[Element] Reals] /. zsInXY;
nx = D[n, x];
ny = D[n, y];
LL = Table[-rx[[i]].nx[[i]], {i, 1, 2}];
MM = Table[-rx[[i]].ny[[i]], {i, 1, 2}];
NN = Table[-ry[[i]].ny[[i]], {i, 1, 2}];
curvatures = 
  Solve[a k^2 + b k + c == 0, k] /. {a -> (EE GG - FF^2), 
    b -> -(LL GG - 2 MM FF + NN EE), c -> (LL NN - MM^2)};
abc = {a -> (EE MM - FF LL), b -> (EE NN - GG LL), 
   c -> (FF NN - GG MM)};
directionsDx = Solve[a dx^2 + b dx dy + c dy^2 == 0, dx];
directionsSym = 
  Table[(dx rx[[j]] + dy ry[[j]]) /. 
      directionsDx[[i]] /. {dy -> 1} /. 
    Table[abc[[k, 1]] -> abc[[k, 2, j]], {k, 1, 3}], {j, 1, 2}, {i, 1,
     2}];
w = {u, v};
rw = Table[(rx[[i]] (dx /. directionsDx[[j]]) + ry[[i]] dy) /. 
    dy -> 1 , {i, 1, 2}, {j, 1, 2}];
nw = Table[(nx[[i]] (dx /. directionsDx[[j]]) + ny[[i]] dy) /. 
    dy -> 1 , {i, 1, 2}, {j, 1, 2}];
rwrw = Table[
   rw[[i, j]].rw[[i, j]] /. 
    Table[abc[[k, 1]] -> abc[[k, 2, i]], {k, 1, 3}], {i, 1, 2}, {j, 1,
     2}];
curvaturesWithDirections = 
  Table[(nw[[i, j]].rw[[i, j]]/(rw[[i, j]].rw[[i, j]])) /. 
    Table[abc[[k, 1]] -> abc[[k, 2, i]], {k, 1, 3}], {i, 1, 2}, {j, 1,
     2}];'''
    session.evaluate(initModelCMD)
    session.evaluate(coeff)

def modelProjection(session, xGSE):
    xGSE_ = 'xGSE = {{x -> {}, y -> {}, z -> {}}};'.format(*xGSE)
    session.evaluate(xGSE_)
    session.evaluate('contValue = contraction/.coeff/.xGSE')
    numberOfContraction = 0
    for i in range(1, 3):
        if (session.evaluate('cont/.contValue[[{}]]'.format(i))) > 0:
            numberOfContraction +=1
            if numberOfContraction == 1:
                contractionInd = i
            elif numberOfContraction == 2:
                raise Exception('two contractions > 0')
    cont = session.evaluate('cont/.contValue[[{}]]'.format(contractionInd))
    return cont*xGSE


def modelCurvaturesAndDirections(session, xGSE):
    xGSE_ = 'xGSE = {{x -> {}, y -> {}, z -> {}}};'.format(*xGSE)
    session.evaluate(xGSE_)
    session.evaluate('contValue = contraction/.coeff/.xGSE')
    numberOfContraction = 0
    for i in range(1, 3):
        if (session.evaluate('cont/.contValue[[{}]]'.format(i))) > 0:
            numberOfContraction +=1
            if numberOfContraction == 1:
                contractionInd = i
            elif numberOfContraction == 2:
                raise Exception('two contractions > 0')
    session.evaluate('zsInXYValue = zsInXY/.coeff/.xGSE[[{1,2}]]'+'/.contValue[[{}]]'.format(contractionInd))
    for i in range(1, 3):
        if np.abs(session.evaluate('z/.zsInXYValue[[{}]]'.format(i)) - xGSE[2]) < 10**(-4):
            zInd = i
    curvatures_ = np.array(session.evaluate('curvaturesWithDirections[[{}]] /. coeff /. contValue[[{}]] /. xGSE'.format(zInd, contractionInd)))
    curvatures = np.append(curvatures_, 0)
    directions = np.zeros((3, 3))
    directions[:, :2] = np.array(session.evaluate('Normalize /@ (directionsSym[[{}]] /. coeff /. contValue[[{}]]) /. xGSE'.format(zInd, contractionInd))).T
    directions[:, 2] = np.array(session.evaluate('n[[{}]] /. coeff /. xGSE /. contValue[[{}]]'.format(zInd, contractionInd)))
    return curvatures, directions
##

def initGradModel(session, model='dipole', M=-30438, B0=1, M1=-30438, M2=-30438*28):
    '''in dipole model input xGSE in RE get B in nT in Earth's dipole field if M=-30438'''
    if model == 'dipole':
        initModelCMD1 = '''
        b={{3 x z M r^(-5), 3 y z M r^(-5), (3z^2-r^2) M r^(-5)}}/.{{M->{}, r->Sqrt[x^2 + y^2 + z^2]}};'''.format(M)
    elif model == 'magnetosphericField':
        initModelCMD1 = '''
        b1={{3 x z M r^(-5), 3 y z M r^(-5), (3z^2-r^2) M r^(-5)}}/.{{M->{}, r->Sqrt[x^2 + y^2 + z^2]}};
        b2={{3 (x-40) z M r^(-5), 3 y z M r^(-5), (3z^2-r^2) M r^(-5)}}/.{{M->{}, r->Sqrt[(x-40)^2 + y^2 + z^2]}};
        b=b1 + b2'''.format(M1, M2)
    elif model == 'lundquistForceFreeField':
        initModelCMD1 = '''
        bCylin = {{0, B0 BesselJ[1, r], B0 BesselJ[0, r]}}/.B0->{};
        b = TransformedField["Polar" -> "Cartesian", bCylin[[{{1, 2}}]], {{r, \[Theta]}} -> {{x, y}}];
        AppendTo[b, bCylin[[3]] /. r -> Sqrt[x^2 + y^2]]; '''.format(B0)
    initModelCMD2 = '''
        gradB=Grad[b, {x, y, z}];
        grad2B=Grad[Grad[b, {x, y, z}], {x, y, z}]'''
    session.evaluate(initModelCMD1)
    session.evaluate(initModelCMD2)


def gradAndSecondGradB(session, xGSE):
    session.evaluate('''xx={};yy={};zz={};'''.format(*xGSE))
    gradBCMD = '''gradB /. {M -> m, x -> xx, y -> yy, z -> zz} // N'''
    grad2BCMD = '''grad2B /. {M -> m, x -> xx, y -> yy, z -> zz} // N'''
    gradB = np.array(session.evaluate(gradBCMD))
    grad2B = np.array(session.evaluate(grad2BCMD))
    return gradB, grad2B
##

def initCurvatureAndTorsionModel(session, model='dipole', M=-30438, B0=1, M1=-30438, M2=-30438*28):
    '''quadratic gradient of B'''
    if model == 'dipole':
        initModelCMD1 = '''b={{3 x z M r^(-5), 3 y z M r^(-5), (3z^2-r^2) M r^(-5)}}/.{{M->{}, r->Sqrt[x^2 + y^2 + z^2]}};'''.format(M)
    elif model == 'magnetosphericField':
        initModelCMD1 = '''
        b1={{3 x z M r^(-5), 3 y z M r^(-5), (3z^2-r^2) M r^(-5)}}/.{{M->{}, r->Sqrt[x^2 + y^2 + z^2]}};
        b2={{3 (x-40) z M r^(-5), 3 y z M r^(-5), (3z^2-r^2) M r^(-5)}}/.{{M->{}, r->Sqrt[(x-40)^2 + y^2 + z^2]}};
        b=b1 + b2'''.format(M1, M2)
    elif model == 'lundquistForceFreeField':
        initModelCMD1 = '''
        bCylin = {{0, B0 BesselJ[1, r], B0 BesselJ[0, r]}}/.B0->{};
        b = TransformedField["Polar" -> "Cartesian", bCylin[[{{1, 2}}]], {{r, \[Theta]}} -> {{x, y}}];
        AppendTo[b, bCylin[[3]] /. r -> Sqrt[x^2 + y^2]]; '''.format(B0)
    if model == 'magnetosphericField':
        initModelCMD2 = '''bNormalized = b/Sqrt[b.b];
        curvature = Grad[bNormalized, {x, y, z}].bNormalized;
        curvatureNorm = Sqrt[curvature.curvature];
        binormal = Cross[bNormalized, curvature/curvatureNorm];
        torsion = 1/curvatureNorm Grad[curvature, {x, y, z}].bNormalized.binormal;'''
    else:
        initModelCMD2 = '''
        bNormalized = FullSimplify[Normalize[b], {x, y, z, B0} \[Element] Reals];
        curvature = Grad[bNormalized, {x, y, z}].bNormalized;
        binormal = FullSimplify[ Normalize[Cross[b, curvature]], {x, y, z, B0} \[Element] Reals];
        torsion = FullSimplify[ 1/Norm[curvature] Grad[ curvature, {x, y, z}].bNormalized.binormal, {x, y, z, B0} \[Element] Reals]'''
    session.evaluate(initModelCMD1)
    session.evaluate(initModelCMD2)


def curvatureAndTorsion(session, xGSE):
    session.evaluate('''xx={};yy={};zz={};'''.format(*xGSE))
    curvatureCMD = '''curvature /. {x -> xx, y -> yy, z -> zz} // N'''
    torsionCMD = '''torsion /. {x -> xx, y -> yy, z -> zz} // N'''
    curvature = np.array(session.evaluate(curvatureCMD))
    torsion = np.array(session.evaluate(torsionCMD))
    return curvature, torsion
##

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def bundleEvents(xGSEList, column=0, maximumRange=2):
    xGSE1 = xGSEList[:, column]
    diff1 = xGSE1[:, None] - xGSE1[None, :]
    rowInds, columnInds = np.where(np.abs(diff1) < maximumRange)
    xGSEBundles = []
    for i in range(len(xGSEList[:, column])):
        xGSEBundle = xGSEList[columnInds[rowInds==i]]
        if column+1 < xGSEList.shape[-1]:
            xGSEBundle = bundleEvents(xGSEBundle, column=column+1, maximumRange=maximumRange)
        xGSEBundles.append(xGSEBundle)
    return xGSEBundles


def vector2symmetricalMat(G3Reconstructed):
    length = int(np.sqrt(len(G3Reconstructed)*2))
    G3 = np.zeros((length, length))
    xs, ys = np.triu_indices(length)
    G3[xs, ys] = G3Reconstructed
    G3[ys, xs] = G3Reconstructed
    return G3
def symmetricalMat2vector(c):
    xs, ys = np.triu_indices(c.shape[-1])
    cReconstructed = c[xs, ys]
    return cReconstructed


def reconstruct(R4=None, c=None, sixCombs=[(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]):
    RReconstructed = np.zeros((6, 6))
    cReconstructed = np.zeros((3, 6))
    for i, firstTwoIndices in enumerate(sixCombs):
        if firstTwoIndices[0] == firstTwoIndices[1]:
            factor_ = 1
        else:
            factor_ = 2
        if not (c is None):
            cReconstructed[:, i] = c[:, firstTwoIndices[0], firstTwoIndices[1]] * factor_
        if not (R4 is None):
            for k, secondTwoIndices in enumerate(sixCombs):
                if secondTwoIndices[0] == secondTwoIndices[1]:
                    RReconstructed[i, k] = R4[firstTwoIndices[0], firstTwoIndices[1], secondTwoIndices[0], secondTwoIndices[1]] * factor_
                else:
                    RReconstructed[i, k] = 2*R4[firstTwoIndices[0], firstTwoIndices[1], secondTwoIndices[0], secondTwoIndices[1]] * factor_
    return RReconstructed, cReconstructed

##
def iterateOnce(b, x, R, RInverse, R3, G3, eigenSolve=True, eigenSystemOfR4=None, numberOfSpacecrafts=None):
    b0 = np.mean(b, axis=0) - 1/2*np.sum(np.sum(R[None, ...]*G3, axis=-1), axis=-1)
    G = (1/numberOfSpacecrafts*b.T @ x - 1/2*np.sum(np.sum(R3[None, ...]*G3[:, None, ...], axis=-1), axis=-1)) @ RInverse
    c = 2*(1/numberOfSpacecrafts * np.sum(b.transpose((1, 0))[..., None, None] * (x[:, None, :] * x[:, :, None])[None, ...], axis=1) - b0[:, None, None] * R[None, ...] - np.sum(R3[None, ...] * G[:, None, None, :], axis=-1))
    _, cReconstructed = reconstruct(c=c)
    if eigenSolve:
        transformMat = normalized(eigenSystemOfR4[1]).T  # transformMat is \xi_{\mu i} where \mu is index of the new base and i index of the old base
        cReconstructedInNewBase = np.sum(transformMat.T[None, :, :] * cReconstructed[:, :, None], axis=1)
        G3ReconstructedInNewBase = cReconstructedInNewBase / eigenSystemOfR4[0][None, :]
        G3Reconstructed = np.sum(G3ReconstructedInNewBase[:, :, None] * transformMat[None, :, :], axis=1)
    else:
        G3Reconstructed = np.linalg.solve(RReconstructed[None, ...], cReconstructed)
    G3 = np.zeros_like(G3)
    for i in range(3):
        G3[i] = vector2symmetricalMat(G3Reconstructed[i])
    return b0, G, G3

##
def multipointsCalculateGradAndGrad2(xGSEs, b, numberOfTurns=None, silence=False, eigenSolve=True, allowNotConverge=True, x=None, xCenter=None):
    converged = True
    if xCenter is None:
        xCenter = np.mean(xGSEs, axis=0)
    if x is None:
        x = xGSEs - xCenter[None, :]
    numberOfSpacecrafts = x.shape[0]
    R = np.transpose(x, (1, 0)) @ x / numberOfSpacecrafts
    eigenSystemOfR = np.linalg.eig(R)
    permutation = np.argsort(eigenSystemOfR[0])
    timingShape = (np.sqrt(eigenSystemOfR[0])[permutation], eigenSystemOfR[1][:, permutation])
#    print(timingShape[0])
    LDRatio = 2*timingShape[0][2]/np.linalg.norm(xCenter)
    RInverse = np.linalg.inv(R)
    R3 = 1/numberOfSpacecrafts * np.sum(x[:, :, None, None] * x[:, None, :, None] * x[:, None, None, :], axis=0)
    R4 = 1/numberOfSpacecrafts * np.sum(x[:, :, None, None, None] * x[:, None, :, None, None] * x[:, None, None, :, None] * x[:, None, None, None, :], axis=0)
    G3 = np.zeros((3, 3, 3))
    G3Last = np.ones((3, 3, 3))
#    R4Reconstructed = R4.reshape((9, 9))
#    eigenSystemOfR4 = np.linalg.eig(R4Reconstructed)
    RReconstructed, _ = reconstruct(R4=R4)
    eigenSystemOfR4 = np.linalg.eig(RReconstructed)
    if numberOfTurns:
        G3All = np.zeros((numberOfTurns, 3, 3, 3))
        G2All = np.zeros((numberOfTurns, 3, 3))
        b0All = np.zeros((numberOfTurns, 3))
        for turn in range(numberOfTurns):
            if not silence:
                print(turn)
            b0, G, G3 = iterateOnce(b, x, R, RInverse, R3, G3, eigenSolve=True, eigenSystemOfR4=eigenSystemOfR4, numberOfSpacecrafts=numberOfSpacecrafts)
            b0All[turn] = b0
            G2All[turn] = G
            G3All[turn] = G3
    elif numberOfTurns is None:
        turn = 0
        G3All = []
        while np.any(np.abs((G3-G3Last)/G3Last) > 0.00005):
#            if not silence:
#                print(turn)
            turn +=1
            if turn > 10**6:
                if allowNotConverge:
                    converged = False
                    print('not converge')
                    break
                else:
                    raise Exception('not converge')
            G3Last = G3
            b0, G, G3 = iterateOnce(b, x, R, RInverse, R3, G3, eigenSolve=True, eigenSystemOfR4=eigenSystemOfR4, numberOfSpacecrafts=numberOfSpacecrafts)
            G3All.append(G3)
        if not silence:
            print('L/D={}'.format(LDRatio))
            print('converge after {} steps'.format(turn))
        b0All = np.repeat(b0[None, ...], 2, axis=0)
        G2All = np.repeat(G[None, ...], 2, axis=0)
#        G3All = np.repeat(G3[None, ...], 2, axis=0)
        G3All = np.stack(G3All)
    return G3All, G2All, b0All, LDRatio, eigenSystemOfR4, converged
##

def levicivita(arg):
    if len(arg) == 3:
        i = arg[0]
        j = arg[1]
        k = arg[2]
        return (-i+j)*(-i+k)*(-j+k)/2


def makeLeviCivitaTensor(order=3):
    if order == 3:
        levicivitaTensor = np.zeros((3, 3, 3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    levicivitaTensor[i, j, k] = levicivita((i, j, k))
    return levicivitaTensor

##
def curvatureAndTorsionFormSecondGradient(b0LastTurnAllContraction, G2LastTurnAllContraction, G3LastTurnAllContraction, calculateTorsion=True):
    curvatures, normals = mca(gradBs=G2LastTurnAllContraction, bVectorAtCenters=b0LastTurnAllContraction)
    binormals = normalized(np.cross(normalized(b0LastTurnAllContraction), normals))
    returnedVariables = [curvatures, normals, binormals]
    if calculateTorsion:
        b0LastTurnAllContractionMag = np.linalg.norm(b0LastTurnAllContraction, axis=-1)
        torsions = 1/(curvatures*b0LastTurnAllContractionMag**3) * ((binormals[:, None, :] @ G2LastTurnAllContraction @ G2LastTurnAllContraction @ b0LastTurnAllContraction[:, :, None]).squeeze() + np.sum(np.sum(np.sum(b0LastTurnAllContraction[:, None, None, :] * G3LastTurnAllContraction, axis=-1) * b0LastTurnAllContraction[:, None, :], axis=-1) * binormals, axis=-1))
        returnedVariables.append(torsions)
    return returnedVariables

##
def errorVSLD(xCenter, xInCenterOfMass, numberOfContraction=20, model='lundquistForceFreeField', wlSession=None, errorDef='averaged', silence=False, numberOfTurns=None, returnG3CurTor=False, returnG2AndError=False):
    numberOfSpacecrafts = len(xInCenterOfMass)
    gradB, grad2B = gradAndSecondGradB(wlSession, xCenter)
    curvatureVector, torsion = curvatureAndTorsion(wlSession, xCenter)
    curvature = np.linalg.norm(curvatureVector)
    ## error vs L/D
    G3LastTurnAllContraction = np.zeros((numberOfContraction, 3, 3, 3))
    G2LastTurnAllContraction = np.zeros((numberOfContraction, 3, 3))
    b0LastTurnAllContraction = np.zeros((numberOfContraction, 3))
    R = xInCenterOfMass.T @ xInCenterOfMass / numberOfSpacecrafts
    eigenSystemOfR = np.linalg.eig(R)
    permutation = np.argsort(eigenSystemOfR[0])
    timingShape = (np.sqrt(eigenSystemOfR[0])[permutation], eigenSystemOfR[1][:, permutation])
    if model == 'dipole' or model == 'magnetosphericField':
        contractions_ = 10**np.linspace(-3, -0.5, numberOfContraction)/timingShape[0][2]
    elif model == 'lundquistForceFreeField':
        contractions_ = np.linspace(0.01, 0.6, numberOfContraction)/timingShape[0][2]
    contractions = np.linalg.norm(xCenter)*contractions_ / 2
    LDRatios = np.zeros(numberOfContraction)
    for k, contraction in enumerate(contractions):
        if not silence:
            print(k)
        x = xInCenterOfMass*contraction
        xGSEs = x + xCenter
        if model == 'dipole':
            b = dipoleField(xGSEs, M=-30438)
        elif model == 'lundquistForceFreeField':
            b = lundquistForceFreeField(x=xGSEs, B0=60)
        elif model == 'magnetosphericField':
            b = magnetosphericField(xGSE=xGSEs)
#        G3All, G2All, b0All, LDRatios[k], eigenSystemOfR4, converged = multipointsCalculateGradAndGrad2(xGSEs, b, numberOfTurns=1000, silence=True, eigenSolve=True)
        G3All, G2All, b0All, LDRatios[k], eigenSystemOfR4, converged = multipointsCalculateGradAndGrad2(xGSEs, b, numberOfTurns=numberOfTurns, silence=silence, eigenSolve=True)
        G3LastTurnAllContraction[k] = G3All[-1]
        G2LastTurnAllContraction[k] = G2All[-1]
        b0LastTurnAllContraction[k] = b0All[-1]
    print('end')
    if errorDef == 'averaged':
        error = 100 * np.abs((G3LastTurnAllContraction-grad2B[None, ...])/np.mean(np.abs(grad2B)))
        G2Error = 100 * np.abs((G2LastTurnAllContraction-gradB[None, ...])/np.mean(np.abs(gradB)))
    else:
        error = 100 * np.abs(G3LastTurnAllContraction/grad2B[None, ...] - 1)
        G2Error = 100 * np.abs(G2LastTurnAllContraction/gradB[None, ...] - 1)
    curvatures, normals, binormals, torsions = curvatureAndTorsionFormSecondGradient(b0LastTurnAllContraction, G2LastTurnAllContraction, G3LastTurnAllContraction, calculateTorsion=True)
    errorCurvatures = 100 * np.abs(curvatures/curvature - 1)
    errorTorsions = 100 * np.abs(torsions/torsion - 1)
    returnedVariables = [error, errorCurvatures, errorTorsions, LDRatios]
    if returnG3CurTor:
        returnedVariables.extend([G3LastTurnAllContraction, curvatures, torsions])
    if returnG2AndError:
        returnedVariables.extend([G2LastTurnAllContraction, G2Error])
    return returnedVariables

##
def calculateCurvaturesAndTorsions(xCenters, xInCenterOfMass, model='magnetosphericField', wlSession=None, numberOfTurns=None, silence=True):
    numberOfXCenters = len(xCenters)
    G3LastTurnAllXCenters = np.zeros((numberOfXCenters, 3, 3, 3))
    G2LastTurnAllXCenters = np.zeros((numberOfXCenters, 3, 3))
    b0LastTurnAllXCenters = np.zeros((numberOfXCenters, 3))
    curvatures = np.zeros(numberOfXCenters)
    torsions = np.zeros(numberOfXCenters)
    LDRatios = np.zeros(numberOfXCenters)
    for k in range(numberOfXCenters):
        if len(xInCenterOfMass.shape) == 2:
            x = xInCenterOfMass
        elif len(xInCenterOfMass.shape) > 2:
            x = xInCenterOfMass[k]
        xCenter = xCenters[k]
        print("xCenter:")
        print(xCenter)
        gradB, grad2B = gradAndSecondGradB(wlSession, xCenter)
        curvatureVector, torsion = curvatureAndTorsion(wlSession, xCenter)
        curvature = np.linalg.norm(curvatureVector)
        curvatures[k] = curvature
        torsions[k] = torsion
        xGSEs = x + xCenter
        if model == 'dipole':
            b = dipoleField(xGSEs, M=-30438)
        elif model == 'lundquistForceFreeField':
            b = lundquistForceFreeField(x=xGSEs, B0=60)
        elif model == 'magnetosphericField':
            b = magnetosphericField(xGSE=xGSEs)
    #        G3All, G2All, b0All, LDRatios[k], eigenSystemOfR4, converged = multipointsCalculateGradAndGrad2(xGSEs, b, numberOfTurns=1000, silence=True, eigenSolve=True)
        G3All, G2All, b0All, LDRatios[k], eigenSystemOfR4, converged = multipointsCalculateGradAndGrad2(xGSEs, b, numberOfTurns=numberOfTurns, silence=True, eigenSolve=True)
        G3LastTurnAllXCenters[k] = G3All[-1]
        G2LastTurnAllXCenters[k] = G2All[-1]
        b0LastTurnAllXCenters[k] = b0All[-1]
    curvatures10Points, normals, binormals, torsions10Points = curvatureAndTorsionFormSecondGradient(b0LastTurnAllXCenters, G2LastTurnAllXCenters, G3LastTurnAllXCenters, calculateTorsion=True)
    return curvatures, torsions, curvatures10Points, torsions10Points, LDRatios


##
def myLegend(ax, linestyles, colors, pos, kappa=False, tau=False):
    bottom, left, width, height, wholeWidth = pos
    xs, ys = np.triu_indices(3)
    ax.text(left-wholeWidth-0.05, bottom+height*6, '$i=$', transform=ax.transAxes)
    for indIn6 in range(6):
        for indIn3 in range(3):
            xRange = np.linspace(left+wholeWidth*indIn3, left+width*(indIn3+1), 2)
            yRange = np.repeat(bottom+height*indIn6, 2)
            ax.plot(xRange, yRange, ls=linestyles[indIn3], color=colors[indIn6], transform=ax.transAxes)
        row = xs[indIn6]
        column = ys[indIn6]
        label_ = '$B_{{i,{},{}}}$'.format(row+1, column+1)
        ax.text(left-wholeWidth-0.05, bottom+height*indIn6, label_, transform=ax.transAxes)
    for indIn3 in range(3):
        label_ = '${}$'.format(indIn3+1)
        ax.text(left+wholeWidth*indIn3, bottom+height*6, label_, transform=ax.transAxes)
    if kappa:
        ax.text(left+wholeWidth*3, bottom+height*5, '$\delta \kappa$', transform=ax.transAxes)
        ax.plot(np.linspace(left+wholeWidth*4, left+wholeWidth*4+width, 2), np.repeat(bottom+height*5+height/4, 2), ls='-.', color='k', transform=ax.transAxes)
    if tau:
        ax.text(left+wholeWidth*3, bottom+height*4, '$\delta \\tau$', transform=ax.transAxes)
        ax.plot(np.linspace(left+wholeWidth*4, left+wholeWidth*4+width, 2), np.repeat(bottom+height*4+height/4, 2), ls=':', color='k', transform=ax.transAxes)
