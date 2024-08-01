__author__ = 'Yufei Zhou'

import matplotlib as mpl
import numpy as np
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import matplotlib.axes._axes
import space_database_analysis.dataAnalysisTools as dat
import cdflib
import datetime as dt
import logging

plt.rcParams['font.size'] = 24
plt.rcParams['font.size'] = 18
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['axes.grid'] = 'True'


## plot related
def plot_time_series(self, *args, scalex=True, scaley=True, data=None, **kwargs):
    '''
    Time series may contain gap, this function replace the original method plot to better tackle the gaps in plot by not plotting them.
    Parameters:
        gap_threshold: two data point with a gap larger than this will not be connected in plot. If this parameter is not given, 5 * minimum gap in time series will be used.
        time_data: when the first parameter args[0] is not time, this parameter should be provided to break all data into blocks. This parameter can be used in the case of plotting trajectory of spacecraft that is not continuous
    '''
    assert len(args) == 2
    x = args[0]
    y = args[1]
    gap_threshold = kwargs.get('gap_threshold')
    time_data = kwargs.get('time_data')
    if time_data is None:
        tdiff = np.diff(x)
    else:
        tdiff = np.diff(time_data)
        _ = kwargs.pop('time_data')

    try:
        _ = kwargs.pop('gap_threshold')
    except: pass
    try:
        enable = kwargs.pop('enable')
    except: enable = True
    if enable:
        if gap_threshold is None:
            gap_threshold = 5*np.min(tdiff)
        break_points = np.nonzero(tdiff > gap_threshold)[0] + 1
        break_points = np.insert(break_points, 0, 0)
        break_points = np.append(break_points, len(x))
        plots_ = []
        for ind in range(len(break_points)-1):
            s_ = slice(*break_points[ind:ind+2])
            x_ = x[s_]
            y_ = y[s_]
            plot_ = self.plot(x_, y_, scalex=scalex, scaley=scaley, data=data, **kwargs)
            plots_.append(plot_)
        return plots_
    else:
        return self.plot(*args, scalex=scalex, scaley=scaley, data=data, **kwargs)

mpl.axes.Axes.plot_time_series = plot_time_series


def plotCartesianVectorTimeSeries(ax, t, data, norm=True, label=None, **kwargs):
    labelsF = ['${}_x$', '${}_y$', '${}_z$']
    colors = ['m', 'g', 'b']
    if label:
        labels = [labelF.format(label) for labelF in labelsF]
    for ind in range(3):
        ax.plot_time_series(t, data[:, ind], color=colors[ind], label=labels[ind], **kwargs)
    if norm:
        ax.plot_time_series(t, np.linalg.norm(data, axis=-1), color='k', label=label, **kwargs)
    ax.set_ylabel(label)
    ax.grid(True)
    return ax

## <<<< plot time format
'''
Usage:
    myFormatter = dat.dayFormatter
    ax.xaxis.set_major_formatter(myFormatter)
'''
def format_month(t, pos=None):
    return cdflib.cdfepoch.encode(t)[8:13]

def format_monthTT2000(t, pos=None):
    tStr = cdflib.cdfepoch.encode(int(t))
    return tStr[8:13] + tStr[14:16]

def format_day(t, pos=None):
    return cdflib.cdfepoch.encode(t)[11:16]

def format_dayTT2000(t, pos=None):
    return cdflib.cdfepoch.encode(int(t))[11:16]

def format_UT(t, pos=None):
    tStr = cdflib.cdfepoch.encode(t)
    return tStr[11:13] + tStr[14:19]

def format_UTTT2000(t, pos=None):
    tStr = cdflib.cdfepoch.encode_tt2000(t)
    return tStr[11:13] + tStr[14:19]

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

def format_doy(t, pos=None):
    '''day of year'''
    epochs = dat.Epochs(CDF_EPOCH=t)
    return epochs.get_data('datetime').strftime('%j')

def format_doyThour(t, pos=None):
    '''day of year'''
    epochs = dat.Epochs(CDF_EPOCH=t)
    return epochs.get_data('datetime').strftime('%jT%H')

monthFormatter = FuncFormatter(format_month)
monthFormatterTT2000 = FuncFormatter(format_monthTT2000)
dayFormatter = FuncFormatter(format_day)
dayFormatterTT2000 = FuncFormatter(format_dayTT2000)
utFormatter = FuncFormatter(format_UT)
utFormatterTT2000 = FuncFormatter(format_UTTT2000)
hourFormatter = FuncFormatter(format_hour)
hourFormatterTT2000 = FuncFormatter(format_hourTT2000)
minFormatter = FuncFormatter(format_min)
minFormatterTT2000 = FuncFormatter(format_minTT2000)
hourMinFormatterTT2000 = FuncFormatter(format_hourMinTT2000)
doyFormatter = FuncFormatter(format_doy)
doyThourFormatter = FuncFormatter(format_doyThour)

## plot time format >>>>

def plotTimeSeriesOfAngle(ax, t, angle, **para):
    deltaAngle = np.diff(angle)
    upArgs, = np.nonzero(deltaAngle < -180)
    upAngleStart = angle[upArgs]
    upAngleEnd = angle[upArgs+1] + 360
    tUpIntersection = (360 - upAngleStart)*(t[upArgs+1] - t[upArgs])/(upAngleEnd - upAngleStart) + t[upArgs]
    numberOfUpIntersection = len(tUpIntersection)
    logging.debug('number of up intersection: {}'.format(numberOfUpIntersection))
    if numberOfUpIntersection > 0:
        logging.debug('up intersection:')
        logging.debug(cdflib.cdfepoch.breakdown(tUpIntersection))
        data = np.stack([angle, t], axis=-1)
        dataUpBlocks = dat.makeBlocks(data, upArgs+1)
        for upBlockInd in range(numberOfUpIntersection+1):
            dataUpBlock = dataUpBlocks[upBlockInd]
            if upBlockInd < numberOfUpIntersection:
                dataUpBlock = np.append(dataUpBlock, np.array([360, tUpIntersection[upBlockInd]])[None, :], axis=0)
            if upBlockInd > 0:
                dataUpBlock = np.insert(dataUpBlock, 0, np.array([0, tUpIntersection[upBlockInd-1]])[None, :], axis=0)
            angle = dataUpBlock[:, 0]
            t = dataUpBlock[:, -1]
            deltaAngle = np.diff(angle, axis=0)
            downArgs, = np.nonzero(deltaAngle >= 180)
            logging.debug("current up intersection:")
            logging.debug(cdflib.cdfepoch.breakdown(t[0]))
            tOfInterest = dat.Epoch(dateTime=dt.datetime(2002, 2, 9, 0, 34, 0)).epoch
            if np.abs(tOfInterest - t[0]) < 5000:
                logging.debug('current block:')
                logging.debug(cdflib.cdfepoch.breakdown(dataUpBlock[:, -1]))
                logging.debug('angle:')
                logging.debug(dataUpBlock[:, 0])
                logging.debug('down args:')
                logging.debug(downArgs)
            if len(downArgs) > 0:
                downAngleStart = angle[downArgs]
                downAngleEnd = angle[downArgs+1] - 360
                tDownIntersection = -downAngleStart*(t[downArgs+1] - t[downArgs])/(downAngleEnd - downAngleStart) + t[downArgs]
                dataDownBlocks = dat.makeBlocks(dataUpBlock, downArgs+1)
                numberOfDownIntersection = len(tDownIntersection)
                logging.debug('number of down intersection: {}'.format(numberOfDownIntersection))
                for downBlockInd in range(numberOfDownIntersection+1):
                    dataDownBlock = dataDownBlocks[downBlockInd]
                    if downBlockInd < numberOfDownIntersection:
                        dataDownBlock = np.append(dataDownBlock, np.array([0, tDownIntersection[downBlockInd]])[None, :], axis=0)
                    if downBlockInd > 0:
                        dataDownBlock = np.insert(dataDownBlock, 0, np.array([360, tDownIntersection[downBlockInd-1]])[None, :], axis=0)
                    if np.abs(tOfInterest - t[0]) < 5000:
                        logging.debug('current down block:')
                        logging.debug(cdflib.cdfepoch.breakdown(dataDownBlock[:, -1]))
                        logging.debug('angle:')
                        logging.debug(dataDownBlock[:, 0])
                    logging.debug("current down intersection:")
                    logging.debug(cdflib.cdfepoch.breakdown(dataDownBlock[0, -1]))
                    plot_, = ax.plot(dataDownBlock[:, -1], dataDownBlock[:, 0], **para)
            else:
                plot_, = ax.plot(dataUpBlock[:, -1], dataUpBlock[:, 0], **para)
    else:
        plot_, = ax.plot(t, angle, **para)
    return plot_


def standardizePhaseSpaceDensity(f, theta, phi, energyTable):
    '''
    Purpose:
    Parameters:
        f: origional phase space density
        theta: theta table
        phi: phi table
    returns:
        phaseSpaceDensity: an array of four dimensions [time, energyTable, vthetaTable, vPhiTable]
        vThetaTable: sorted theta table from 0 to pi
        vPhiTable: sorted phi table from 0 to 2pi
        energyTable: sorted energy table
    '''
    argEnergy = np.argsort(energyTable, axis=-1)
    energyTable = energyTable[..., argEnergy]
    f_ = f[:, :, :, argEnergy]
    vThetaTable_ = dat.sphericalAngleTransform(np.radians(theta), coordinate='theta', standardOutRange=True, inRange=[-np.pi/2, np.pi/2])
    vPhiTable_ = dat.sphericalAngleTransform(np.radians(phi), coordinate='phi', standardOutRange=True, inRange=[-np.pi, np.pi])
    argTheta = np.argsort(vThetaTable_)
    argPhi = np.argsort(vPhiTable_)
    vThetaTable = vThetaTable_[argTheta]
    vPhiTable = vPhiTable_[argPhi]
    f_ = f_[:, argTheta]
    f_ = f_[:, :, argPhi]
    phaseSpaceDensity = np.moveaxis(f_, source=-1, destination=1)
    return phaseSpaceDensity, vThetaTable, vPhiTable, energyTable


def plotMultiplePhaseSpaceDensity(epochStart, t, f, theta, phi, energyTable, datasetName=None, plotTGap=10, integration=None):
    '''
    Purpose:
        The name of this function is self-explanatory
    Parameters:
        epochStart:
        plotTGap: = 10 # in second
    '''
#    nPoints = np.prod(f_.shape[1:])
    f_, vThetaTable, vPhiTable, energyTable = standardizePhaseSpaceDensity(f, theta, phi, energyTable)
    tIndOfStart = epochStart.epochRecord(t, tolerance=20)
    tSteps = np.diff(t)
    uni_, counts = np.unique(tSteps, return_counts=True)
    tStep = uni_[counts.argmax()]
    print(uni_)
    print(counts)
    tGapN = int(np.ceil(plotTGap*1000/tStep))
    if integration is None:
        fig, axes = plt.subplots(nrows=3, ncols=4, layout='constrained', subplot_kw={'projection': 'polar'})
    else:
        fig, axes = plt.subplots(nrows=3, ncols=4, layout='constrained')
    for axr in range(3):
        for axc in range(4):
            ax = axes[axr][axc]
            tInd = (axr*4+axc) * tGapN + tIndOfStart
            print(tInd)
            if tInd > len(f_)-1:
                break
            phaseSpaceDensity = f_[tInd]
            datasetEnergyTableDependantOnTime = ['ls_sw']
            for ind_ in range(len(datasetEnergyTableDependantOnTime)):
                datasetEnergyTableDependantOnTime.append(datasetEnergyTableDependantOnTime[ind_].upper())
            if  any([datasetName_ in datasetName for datasetName_ in datasetEnergyTableDependantOnTime]):
                energyTableDependantOnTime = True
            else:
                energyTableDependantOnTime = False
            if energyTableDependantOnTime:
                vTable = 13.841 * np.sqrt(energyTable[tInd]) # in km/s
            else:
                vTable = 13.841 * np.sqrt(energyTable) # in km/s
            if integration is None:
                plotPhaseSpaceDensityCut(ax, phaseSpaceDensity, vTable=vTable, vThetaTable=vThetaTable, vPhiTable=vPhiTable)
            else:
                plotPhaseSpaceDensity2D(ax, phaseSpaceDensity, vTable=vTable, vThetaTable=vThetaTable, vPhiTable=vPhiTable, integration=integration)
            ax.set_title(dt.datetime(*cdflib.cdfepoch.breakdown(t[tInd])[:6]))
    return fig


def plotPhaseSpaceDensity2D(ax, phaseSpaceDensity, vTable, vThetaTable, vPhiTable, integration=None, integrationStepsN=100):
    '''
    Purpose:
    Parameters:
        integration: if None, plot a cut of phase space density. Otherwise, it should be a character, <'x', 'y', 'z'>, which represents a dimension to be integrated out.
    '''
    interp = interpolatePhaseSpaceDensity(phaseSpaceDensity, vTable, vThetaTable, vPhiTable)
    vRange = vTable[-1] * np.array([-1, 1])
    vPlotGrid = np.linspace(*vRange, 100)
    vGrid = np.meshgrid(vPlotGrid, vPlotGrid, indexing='ij')
    vzGrid = np.zeros_like(vGrid[0])
    if integration is None:
        phaseSDInterp = interp(*vGrid, vzGrid)
        ax.pcolormesh(*vGrid, np.log10(phaseSDInterp), shading='auto')
    else:
        if integration == 'z':
            vGridAllIntegration = [np.repeat(grid[..., None], integrationStepsN, axis=-1) for grid in vGrid]
            vzGridAllIntegration = vzGrid[..., None] + np.linspace(*vRange, integrationStepsN)
            phaseSDInterp = np.sum(interp(*vGridAllIntegration, vzGridAllIntegration), axis=-1)
    ax.pcolormesh(*vGrid, np.log10(phaseSDInterp), shading='auto')
#    xyLim = [-2000, 2000]
#    xyTicks = np.arange(-2000, 2001, 1000)
#    ax.set_ylim(xyLim)
#    ax.set_xlim(xyLim)
#    ax.set_xticks(xyTicks, labels=[])
#    ax.set_yticks(xyTicks)
    ax.grid(True)
#    ax.axis('equal')


def plotPhaseSpaceDensityCut(ax, phaseSpaceDensity, vTable, vThetaTable, vPhiTable):
    vxyTable = vTable
    phaseSpaceDensityCut = np.log10(phaseSpaceDensity[:, 3:5].sum(axis=1))
    vMesh, vPhiMesh = np.meshgrid(vxyTable, vPhiTable, indexing='ij')
    ax.pcolormesh(vPhiMesh, vMesh, phaseSpaceDensityCut)
    ax.grid(True)

def interpolatePhaseSpaceDensity(phaseSpaceDensity, vTable, vThetaTable, vPhiTable):
    phaseSpaceDensity = phaseSpaceDensity.flatten()
    vPointsMeshgridSpherical = np.meshgrid(vTable, vThetaTable, vPhiTable, indexing='ij')
    vPointsCartesian = dat.spherical2cartesian(np.stack(vPointsMeshgridSpherical, axis=-1).reshape(-1, 3))
    interp = dat.interpolate.NearestNDInterpolator(vPointsCartesian, phaseSpaceDensity)
    return interp


def integratingPhaseSpaceDensity(phaseSpaceDensity, vTable, vThetaTable, vPhiTable, integration, integrationStepsN=100):
    '''
    Purpose:
    Parameters:
        phaseSpaceDensity: a ndarray of shape [time, vTable, vThetaTable, vPhiTable]
        integration: if None, plot a cut of phase space density. Otherwise, it should be a list of numbers, e.g. [1, 2], in which each number represents the dimension to be integrated out.
    '''
    vRange = vTable[-1] * np.array([-1, 1])
    vInterpGrid = np.linspace(*vRange, integrationStepsN)
    interpoationGrid = np.meshgrid(vInterpGrid, vInterpGrid, vInterpGrid, indexing='ij')
    numberOfRecords = phaseSpaceDensity.shape[0]
#    numberOfIntegrationDim = len(integration)
#    integratedShape = [integrationStepsN] *(3 - numberOfIntegrationDim)
    interpolatedData = []
#    psdIntegrated = np.zeros((numberOfRecords, *integratedShape))
    for tInd in range(numberOfRecords):
        print(tInd)
        phaseSpaceDensity_ = phaseSpaceDensity[tInd]
        interp = interpolatePhaseSpaceDensity(phaseSpaceDensity_, vTable, vThetaTable, vPhiTable)
        interpolatedData_ = interp(*interpoationGrid)
        interpolatedData.append(interpolatedData_)
    interpolatedData = np.stack(interpolatedData, axis=0)
    integration = np.sort(np.array(integration))[::-1]
    phaseSDIntegrated_ = interpolatedData
    for axis in integration:
        phaseSDIntegrated_ = phaseSDIntegrated_.sum(axis=axis+1)
    psdInterpInteg = phaseSDIntegrated_
    return psdInterpInteg, vInterpGrid

def plot_time_series_spectrogram(fig, axes, epochs, data, nperseg, noverlap, gap_threshold='6*', epoch_fmt='CDF_TIME_TT2000', vmin=None, vmax=None, cax_width=0.02):
    '''

    '''
    cax_gap = cax_width / 3
    t = epochs.get_data(epoch_fmt)
    if epoch_fmt == 'CDF_EPOCH':
        factor = 10**3
    elif epoch_fmt == 'CDF_TIME_TT2000':
        factor = 10**9
    frequencyTable, t_spectrogram_split, data_spectrogram_split = dat.spectrogram(t/factor, data, nperseg, noverlap, gap_threshold='6*', window='hamming')
    for t_ in t_spectrogram_split:
        t_[:] *= factor
    energy_spectrum = np.concatenate(data_spectrogram_split, axis=0)
    tTable = np.concatenate(t_spectrogram_split, axis=0)
    if isinstance(axes, list):
        pass
    else:
        axes = [axes]
        energy_spectrum = energy_spectrum[..., None]
    for componentInd in range(energy_spectrum.shape[-1]):
        ax = axes[componentInd]
        data = energy_spectrum[..., componentInd]
        if np.any(tTable):
            tMesh, energyMesh = np.meshgrid(tTable, frequencyTable, indexing='ij')
            fUni = np.unique(data)
            if vmin is None:
                vmin = fUni[fUni>0][0]
            if vmax is None:
                vmax = data.max()
            pcm_ = ax.pcolormesh(tMesh, energyMesh, data, norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax), shading='auto', cmap='jet')
    #            pcm_ = ax.pcolormesh(tMesh, energyMesh, data, shading='auto', cmap='jet')
    #        ax.set_yscale('log')
    #        y_major = pt.mpl.ticker.LogLocator(base = 10.0, numticks = 1)
    #        ax.yaxis.set_major_locator(y_major)
    #        y_minor = pt.mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
    #        ax.yaxis.set_minor_locator(y_minor)
    #        ax.set_yticks(10.0**np.array([-2, -1, 0]))
    #        ax.set_ylim(frequencyTable[[1, -1]])
            ax.set_ylabel('freq [Hz]')
            ax_pos = ax.get_position()
            cax_pos = [ax_pos.x1+cax_gap, ax_pos.y0, cax_width, ax_pos.y1-ax_pos.y0]
            cax = fig.add_axes(cax_pos)
            cbar = fig.colorbar(pcm_, cax=cax, orientation='vertical', location='right', extend='max', label=r'$\log_{10}$(nT$^2$/Hz)', ticks=mpl.ticker.LogLocator(base=10.0, numticks=3))
    #        cax.set_yticks(10**np.array([6, 7, 8, 9]))
            #cax.set_ylabel('DEF')

def plot_PSD_spectrogram_from_partial_numberdensity(fig, ax, epochs, energy_table, energy_delta, partial_numberdensity, epoch_fmt='CDF_TIME_TT2000', cax_width=0.02):
    '''
    this function work for mms fpi part-moments, partial-numberdensity

    '''
    cax_gap = cax_width / 3
    numberdensity_over_energy = np.zeros_like(partial_numberdensity)
    numberdensity_over_energy[:, :-1] = -(np.diff(partial_numberdensity, axis=-1))
    numberdensity_over_energy[:, -1] = partial_numberdensity[:, -1]
    f_xv = numberdensity_over_energy/(2*energy_delta*dat.eV)/(2*energy_table*dat.eV/dat.mass_proton)**(1/2) * dat.mass_proton/2 *10**18 # in s^3 km^-3
    data = f_xv
    tMesh, energyMesh = np.meshgrid(epochs.get_data(epoch_fmt), energy_table, indexing='ij')
    fUni = np.unique(data)
    vmin = fUni[fUni>0][0]
    pcm_ = ax.pcolormesh(tMesh, energyMesh, data, norm=mpl.colors.LogNorm(vmin=vmin, vmax=data.max()), shading='auto', cmap='jet')
    ax.set_yscale('log')
    y_major = mpl.ticker.LogLocator(base = 10.0, numticks = 6)
    ax.yaxis.set_major_locator(y_major)
    y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
    ax.yaxis.set_minor_locator(y_minor)
    ax_pos = ax.get_position()
    cax_pos = [ax_pos.x1+cax_gap, ax_pos.y0, cax_width, ax_pos.y1-ax_pos.y0]
    cax = fig.add_axes(cax_pos)
    cbar = fig.colorbar(pcm_, cax=cax, orientation='vertical', location='right', extend='max', label=r'f [$\mathrm{s}^3\mathrm{km}^{-6}$]', ticks=mpl.ticker.LogLocator(base=10.0, numticks=3))
    ax.set_ylabel('energy [eV]')
