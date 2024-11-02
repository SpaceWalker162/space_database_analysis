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
            kwargs:
                enable: bool, Defualt:True. If true, this function will plot the data into segments divided by gaps in the data. Otherwise, this function is identical to matplotlib.axes.Axes.plot
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
    y_major = mpl.ticker.MaxNLocator(nbins=4, min_n_ticks=3)
    self.yaxis.set_major_locator(y_major)
    if enable:
        if gap_threshold is None:
#            gap_threshold = 5*np.min(tdiff)
            values, counts = np.unique(tdiff, return_counts=True)
            gap_threshold = 5*values[np.argmax(counts)]
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
#        self.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
        return plots_
    else:
        plot_ = self.plot(*args, scalex=scalex, scaley=scaley, data=data, **kwargs)
#        self.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(5))
        return plot_
mpl.axes.Axes.plot_time_series = plot_time_series

def plot_multiple_time_series(self, t, data, labels, **kwargs):
    '''
    Parameters:
        t: a list of 1d array
        data: a list of 1d array
    '''
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color'][:len(data)]

    for ind in range(len(data)):
        self.plot_time_series(t[ind], data[ind], color=colors[ind], label=labels[ind], **kwargs)

    y_major = mpl.ticker.MaxNLocator(nbins=4, min_n_ticks=4)
    self.yaxis.set_major_locator(y_major)
    self.grid(True)

    labelHeight = 0.2
    vectorLabelX = 1.01
    vectorLabelPos = np.zeros((4, 2))
    vectorLabelPos[:, 0] = vectorLabelX
    vectorLabelPos[:, 1] = 0.85 - np.arange(0, 4) * labelHeight
    for ind in range(len(labels)):
        self.text(vectorLabelPos[ind, 0], vectorLabelPos[ind, 1], s=labels[ind], color=colors[ind], transform=self.transAxes)
mpl.axes.Axes.plot_multiple_time_series = plot_multiple_time_series

def plotCartesianVectorTimeSeries(self, t, data, norm=True, label=None, **kwargs):
    labelsF = ['${}_x$', '${}_y$', '${}_z$']
    colors = ['m', 'g', 'b']
    if label:
        labels = [labelF.format(label) for labelF in labelsF]
    for ind in range(3):
        self.plot_time_series(t, data[:, ind], color=colors[ind], label=labels[ind], **kwargs)
    if norm:
        self.plot_time_series(t, np.linalg.norm(data[..., :3], axis=-1), color='k', label=label, **kwargs)
    y_major = mpl.ticker.MaxNLocator(nbins=4, symmetric=True, min_n_ticks=4)
    self.yaxis.set_major_locator(y_major)
    self.set_ylabel(label)
    self.grid(True)

    labelHeight = 0.2
    vectorLabelX = 1.01
    vectorLabelPos = np.zeros((4, 2))
    vectorLabelPos[:, 0] = vectorLabelX
    vectorLabelPos[:, 1] = 0.85 - np.arange(0, 4) * labelHeight

    for ind in range(3):
        self.text(vectorLabelPos[ind, 0], vectorLabelPos[ind, 1], s=labels[ind], color=colors[ind], transform=self.transAxes)
    self.text(vectorLabelPos[3, 0], vectorLabelPos[3, 1], s='$'+label+'$', color='k', transform=self.transAxes)
mpl.axes.Axes.plotCartesianVectorTimeSeries = plotCartesianVectorTimeSeries


def plotVerticalLinesAcrossMultiplePanels(fig, axes, epochs, notations=None, color='k', coordinates='fig'):
    '''
    Parameters:
        coordinates: if "fig", the internal coordinate system for the plotted vertical lines is figure coordinate. In zoomming in using the interactive tool, the vertical lines will remain where they were. If "ax", the internal coordinate system is ax, in which case to make sure the newly added vertical lines will not change the value of ax.get_ylim(), the ylim will be fixed here. Therefore if "ax" is to be used, put this function at the end of the program.
    '''
    if isinstance(notations, str):
        if notations == 'alphabetical':
            notations = [chr(notation) for notation in np.arange(0, len(epochs)) + ord('a')]
        elif notations == 'numerical':
            notations = [str(ind) for ind in range(1, len(epochs)+1)]

    if coordinates == 'fig':
        figInv = fig.transFigure.inverted()
        for ind, epoch_ in enumerate(epochs):
            p1 = figInv.transform(axes[0].transData.transform(np.array([epoch_, axes[0].get_ylim()[1]])))
            p2 = figInv.transform(axes[-1].transData.transform(np.array([epoch_, axes[-1].get_ylim()[0]])))
            points = np.stack([p1, p2])
            fig.add_artist(mpl.lines.Line2D(*points.T, color=color, ls='--', lw=1.5))
    elif coordinates == 'ax':
        for ax in axes:
            ylim = ax.get_ylim()
            for ind, epoch_ in enumerate(epochs):
                ax.plot(np.repeat(epoch_, 2), ylim, color=color, ls='--', lw=1.5)
            ax.set_ylim(ylim)
    if notations is None:
        pass
    elif isinstance(notations, list):
        for ind, epoch_ in enumerate(epochs):
            notation = notations[ind]
            axes[0].text(x=epoch_, y=axes[0].get_ylim()[1], s=notation, transform=axes[0].transData, horizontalalignment='center', verticalalignment='bottom')


def plotSphericalCoordinatesOfVectorTimeSeries(ax, t , data, **para):
    '''
    Parameters:
        t: time
        data: time series of vector data in Cartesian coordinates, such as magnetic field vector array of shape [numberOfEpochs, 3]
    '''
    try:
        subscript = '_' + para.pop('subscript')
    except: subscript = ''
    B = data
    BSpherical = dat.cartesian2spherical(B)/np.pi*180
    _ = ax.plot_time_series(t, BSpherical[:, 1], color='k', ls='-', **para)
    ax.set_ylabel('$\\theta{}$ [deg]'.format(subscript))
    ax.set_ylim([0, 180])
    ax.set_yticks([0, 90, 180])

    axTwin = ax.twinx()
    #dat = reload(dat)
    plotPhi = plotTimeSeriesOfAngle(axTwin, t, BSpherical[:, 2], color='b', ls='-', label='$\\phi$', **para)
    #plotPhi, = axTwin.plot(t, vSpherical[:, 2], color='b', ls='-', label='$\\phi$')
    axTwin.spines["right"].set_color("b")
    axTwin.set_ylabel('$\\phi{}$ [deg]'.format(subscript))
    axTwin.set_ylim([0, 360])
    axTwin.set_yticks([0, 90, 180, 270, 360])
    axTwin.tick_params(which='both', direction='in')
    axTwin.yaxis.label.set_color(plotPhi[0][0].get_color())
    axTwin.tick_params(axis='y', colors=plotPhi[0][0].get_color())
    axTwin.grid(True)

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

def format_show_min_TT2000(t, pos=None):
    '''e.g. 30^\mathrm{m}'''
    return '$' + cdflib.cdfepoch.encode(int(t))[14:16] + '^\mathrm{m}$'

def format_show_hour_min_TT2000(t, pos=None):
    '''e.g. 13^\mathrm{h}30^\mathrm{m}'''
    return '$' + cdflib.cdfepoch.encode(int(t))[11:13] + '^\mathrm{h}'+ cdflib.cdfepoch.encode(int(t))[14:16] + '^\mathrm{m}$'

def format_show_min_sec_TT2000(t, pos=None):
    '''e.g. 30^\mathrm{m}42^\mathrm{s}'''
    return '$' + cdflib.cdfepoch.encode(int(t))[14:16] + '^\mathrm{m}$' +'$' + cdflib.cdfepoch.encode(int(t))[17:19] + '^\mathrm{s}$'


showMinFormatterTT2000 = FuncFormatter(format_show_min_TT2000)
showHourMinFormatterTT2000 = FuncFormatter(format_show_hour_min_TT2000)
showMinSecFormatterTT2000 = FuncFormatter(format_show_min_sec_TT2000)
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
                    logging.debug("current down intersection:")
                    logging.debug(cdflib.cdfepoch.breakdown(dataDownBlock[0, -1]))
                    plot_ = ax.plot_time_series(dataDownBlock[:, -1], dataDownBlock[:, 0], **para)
            else:
                plot_ = ax.plot_time_series(dataUpBlock[:, -1], dataUpBlock[:, 0], **para)
    else:
        plot_ = ax.plot_time_series(t, angle, **para)
    return plot_


def standardizePhaseSpaceDensity(f, energyTable, thetaTable, phiTable, order=[0, 3, 2, 1], mission='Cluster'):
    '''
    Purpose:
    Parameters:
        f: f is the origional phase space density with f[axis0, axis1, axis2, axis3] such that np.transpose(f, order) gives f[t, energyTable, thetaTable, phiTable]. MMS FPI dis-dist is [t, phiTable, thetaTable, energyTable], in which case order should be [0, 3, 2, 1]. Cluster HIA order should be [0, 3, 1, 2]
        thetaTable: theta table in degrees
        phiTable: phi table in degrees
        order: see f above.
    returns:
        phaseSpaceDensity: an array of four dimensions [time, energyTable, vthetaTableTable, vPhiTable]
        vThetaTable: sorted thetaTable table from 0 to pi
        vPhiTable: sorted phiTable table from 0 to 2pi
        energyTable: sorted energy table
    '''
    if mission == 'cluster':
        order = [0, 3, 1, 2]
        f = np.transpose(f, order)
        argEnergy = np.argsort(energyTable, axis=-1)
        energyTable = energyTable[..., argEnergy]
        f_ = f[:, argEnergy]
        vThetaTable_ = dat.sphericalAngleTransform(np.radians(thetaTable), coordinate='thetaTable', standardOutRange=True, inRange=[-np.pi/2, np.pi/2])
        vPhiTable_ = dat.sphericalAngleTransform(np.radians(phiTable), coordinate='phiTable', standardOutRange=True, inRange=[-np.pi, np.pi])
        argTheta = np.argsort(vThetaTable_)
        argPhi = np.argsort(vPhiTable_)
        vThetaTable = vThetaTable_[argTheta]
        vPhiTable = vPhiTable_[argPhi]
        f_ = f_[:, :, argTheta]
        f_ = f_[:, :, :, argPhi]
        phaseSpaceDensity = f_
        return phaseSpaceDensity, energyTable, vThetaTable, vPhiTable
    else:
        raise Exception('not supported')


def plotMultiplePhaseSpaceDensity(epochStart, t, f, theta, phi, energyTable, datasetName=None, plotTGap=10, integration=None):
    '''
    Purpose:
        plot phase space density observation at multiple times.
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


def plotPhaseSpaceDensity2D(ax, phaseSpaceDensity, vTable=None, vThetaTable=None, vPhiTable=None, phasePointsSpherical=None, integration=None, integrationStepsN=200, vPlotRange=None, vmin=None, vmax=None, fig=None, cax_width=0.02, gridPointsNumber=200, levels_max=None, levels_N=None, level_values=None):
    '''
    Purpose:
    Parameters:
        integration: if None, plot a cut of phase space density. Otherwise, it should be a character, <'x', 'y', 'z'>, which represents a dimension to be integrated out.
        levels_max, levels_N, level_values: provide only one of them specifying the levels. levels_max (int) will finds up to a max number of intervals with ticks at nice locations. levels_N (int) specify the exact number of levels evenly logarithmically distributed. level_values (array) specify the exact levels.
    '''
    interp = interpolatePhaseSpaceDensity(phaseSpaceDensity, vTable=vTable, vThetaTable=vThetaTable, vPhiTable=vPhiTable, phasePointsSpherical=phasePointsSpherical)
    if vTable is not None:
        maxv = vTable[-1]
    elif phasePointsSpherical is not None:
        maxv = np.max(phasePointsSpherical[:, 0])
    vRange = maxv * np.array([-1, 1])
    if vPlotRange is None:
        vPlotRange = vRange
    vPlotGrid = np.linspace(*vPlotRange, 100)
    vGrid = np.meshgrid(vPlotGrid, vPlotGrid, indexing='ij') # coordinates grid
    integrationGrid = np.zeros_like(vGrid[0])
    if integration is None:
        phaseSDInterp = interp(*vGrid, integrationGrid).astype(np.float64)
        ax.pcolormesh(*vGrid, np.log10(phaseSDInterp), shading='auto')
    else:
        integration_code = {'x': 1, 'y': 2, 'z': 3}
        if isinstance(integration, str):
            integration = integration_code[integration]
        vGridAllIntegration = [np.repeat(np.expand_dims(grid, axis=integration-1), integrationStepsN, axis=integration-1) for grid in vGrid]
        integrationGridAllIntegration = np.expand_dims(integrationGrid, axis=integration-1) + np.swapaxes(np.linspace(*vRange, integrationStepsN)[:, None, None], 0, integration-1)
        interpolation_grid = []
        for axis in range(1, 4):
            if integration == axis:
                interpolation_grid.append(integrationGridAllIntegration)
            else:
                interpolation_grid.append(vGridAllIntegration.pop(0))
        phaseSDInterp = np.sum(interp(*interpolation_grid).astype(np.float64) * np.diff(vRange)/integrationStepsN, axis=integration-1)
#    pcm_ = ax.pcolormesh(*vGrid, np.log10(phaseSDInterp), shading='auto')
    data = phaseSDInterp
    fUni = np.unique(data)
    if vmin is None:
        vmin = fUni[fUni>0][0]
    if vmax is None:
        vmax = data.max()
#    pcm_ = ax.pcolormesh(*vGrid, phaseSDInterp, norm=mpl.colors.LogNorm(vmin=vmin, vmax=data.max()), shading='auto')
#    pcm_ = ax.pcolormesh(*vGrid, phaseSDInterp, norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax), shading='auto', cmap='viridis')
    provided_level_option_count = 0
    for level_option in [levels_max, levels_N, level_values]:
        if level_option is None: pass
        else: provided_level_option_count +=1
    assert provided_level_option_count == 1
    if levels_max is not None:
        levels = levels_max
    elif levels_N is not None:
        levels = (vmax/vmin)**(np.arange(levels_N+1)/levels_N) * vmin
    elif level_values is not None:
        levels = level_values
    pcm_ = ax.contour(*vGrid, phaseSDInterp, levels=levels, norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax), colors='k')
    pcm_ = ax.contourf(*vGrid, phaseSDInterp, levels=levels, norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax), cmap='viridis')
    if fig is not None:
        cax_gap = cax_width / 3
        ax_pos = ax.get_position()
        cax_pos = [ax_pos.x1+cax_gap, ax_pos.y0, cax_width, ax_pos.y1-ax_pos.y0]
        cax = fig.add_axes(cax_pos)
        cbar = fig.colorbar(pcm_, cax=cax, orientation='vertical', extend='max', label=r'$f$', ticks=mpl.ticker.LogLocator(base=10.0, numticks=3))
#        cax.set_yticks(10.0**np.array([-8, -7, -6, -5, -4]))
#cax.set_ylabel('DEF')
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

def interpolatePhaseSpaceDensity(phaseSpaceDensity, vTable=None, vThetaTable=None, vPhiTable=None, phasePointsSpherical=None):
    phaseSpaceDensity = phaseSpaceDensity.flatten()
    if phasePointsSpherical is None:
        vPointsMeshgridSpherical = np.meshgrid(vTable, vThetaTable, vPhiTable, indexing='ij')
        phasePointsCartesian = dat.spherical2cartesian(np.stack(vPointsMeshgridSpherical, axis=-1).reshape(-1, 3))
    else:
        phasePointsCartesian = dat.spherical2cartesian(phasePointsSpherical)
    interp = dat.interpolate.NearestNDInterpolator(phasePointsCartesian, phaseSpaceDensity)
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
            cbar = fig.colorbar(pcm_, cax=cax, orientation='vertical', extend='max', label=r'$\log_{10}$(nT$^2$/Hz)', ticks=mpl.ticker.LogLocator(base=10.0, numticks=3))
#            cbar = fig.colorbar(pcm_, cax=cax, orientation='vertical', location='right', extend='max', label=r'$\log_{10}$(nT$^2$/Hz)', ticks=mpl.ticker.LogLocator(base=10.0, numticks=3))
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
    cbar = fig.colorbar(pcm_, cax=cax, orientation='vertical', extend='max', label=r'f [$\mathrm{s}^3\mathrm{km}^{-6}$]', ticks=mpl.ticker.LogLocator(base=10.0, numticks=3))
#    cbar = fig.colorbar(pcm_, cax=cax, orientation='vertical', location='right', extend='max', label=r'f [$\mathrm{s}^3\mathrm{km}^{-6}$]', ticks=mpl.ticker.LogLocator(base=10.0, numticks=3))
    ax.set_ylabel('E [eV]')

def plot_omnidirectional_differential_energy_flux_spectrogram(fig, ax, epochs, energy_table, omnidirectional_differential_energy_flux, epoch_fmt='CDF_TIME_TT2000', vmin=None, vmax=None, cax_width=0.02, ylabel='Energy', cbar_label='DEF'):
    '''
    this function work for mms fpi moms, mms1_dis_energyspectr_omni_fast

    '''
    cax_gap = cax_width / 3
    data = omnidirectional_differential_energy_flux
    tMesh, energyMesh = np.meshgrid(epochs.get_data(epoch_fmt), energy_table, indexing='ij')
    fUni = np.unique(data)
    if vmin is None:
        vmin = fUni[fUni>0][0]
    if vmax is None:
        vmax = data.max()
    pcm_ = ax.pcolormesh(tMesh, energyMesh, data, norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax), shading='auto', cmap='jet')
    ax.set_yscale('log')
    y_major = mpl.ticker.LogLocator(base = 10.0, numticks = 6)
    ax.yaxis.set_major_locator(y_major)
    y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
    ax.yaxis.set_minor_locator(y_minor)
    ax_pos = ax.get_position()
    cax_pos = [ax_pos.x1+cax_gap, ax_pos.y0, cax_width, ax_pos.y1-ax_pos.y0]
    cax = fig.add_axes(cax_pos)
#    cbar = fig.colorbar(pcm_, cax=cax, orientation='vertical', extend='max', label=r'DEF [$\mathrm{keV}/(\mathrm{cm}^2\cdot\mathrm{s}\cdot\mathrm{sr}\cdot\mathrm{keV})$]', ticks=mpl.ticker.LogLocator(base=10.0, numticks=3))
    cbar = fig.colorbar(pcm_, cax=cax, orientation='vertical', extend='max', label=cbar_label, ticks=mpl.ticker.LogLocator(base=10.0, numticks=3))
    ax.set_ylabel(ylabel)

## additional information
def add_indices_to_axes(axes, notations='alphabetical', pos=(0.05, 0.97)):
    if isinstance(notations, str):
        if notations == 'alphabetical':
            notations = [chr(notation) for notation in np.arange(0, len(axes)) + ord('a')]
        elif notations == 'numerical':
            notations = [str(ind) for ind in range(1, len(axes)+1)]

    for ax, notation in zip(axes, notations):
        ax.text(*pos, s=notation, transform=ax.transAxes, va='top', ha='left')

def add_summary_brackets_to_axes(fig, axes, text, gap_between_bracket_and_axes):
    figInv = fig.transFigure.inverted()
    bracket_width = 2*plt.rcParams['font.size'] # in display unit
    gap_between_text_and_bracket = 0.5*plt.rcParams['font.size']

    p1 = figInv.transform(axes[0].transAxes.transform(np.array([0, 1]))) - np.array([gap_between_bracket_and_axes, 0])
    p2 = p1 - figInv.transform(np.array([bracket_width, 0]))
    p4 = figInv.transform(axes[-1].transAxes.transform(np.array([0, 0]))) - np.array([gap_between_bracket_and_axes, 0])
    p3 = p4 - figInv.transform(np.array([bracket_width, 0]))
    points = np.stack([p1, p2, p3, p4])
    fig.add_artist(mpl.lines.Line2D(*points.T, color='k', ls='-'))
    text_pos = (p2 + p3)/2
    fig.text(x=text_pos[0] - figInv.transform(np.array([gap_between_text_and_bracket, 0]))[0], y=text_pos[1], s=text, transform=fig.transFigure, horizontalalignment='right', verticalalignment='center', rotation=90)

def find_best_round_gap_for_time_ticks(datetimeRange, ticks_N=6):
    possible_gaps_seconds = np.array([1, 2, 6, 10, 20, 30, 40])
    possible_gaps_minutes_in_seconds = possible_gaps_seconds * 60
    possible_gaps_hours_in_seconds = np.array([1, 2, 5, 6, 10, 12]) * 60 * 60
    possible_gaps = np.concatenate((possible_gaps_seconds, possible_gaps_minutes_in_seconds, possible_gaps_hours_in_seconds))
    round_gap_seconds = possible_gaps[np.argmin(np.abs(possible_gaps - (datetimeRange[1] - datetimeRange[0]).total_seconds()/ticks_N))]
    return dt.timedelta(seconds=int(round_gap_seconds))

def find_CDF_TIME_TT2000_xticks(datetimeRange, round_gap):
    '''
    Parameters:
        round_gap: timedelta object specifying the distance between two major ticks

    '''
    start_ = dat.datetime_ceil(datetimeRange[0], round_gap=round_gap)
    xTicks = []
    while datetimeRange[1] > start_:
        xTicks.append(dat.Epochs(datetime=start_).get_data('CDF_TIME_TT2000')[0])
        start_ += round_gap
    return xTicks

def set_CDF_TIME_TT2000_xticks(self, round_gap=None, ticks_N=6):
    xlim = self.get_xlim()
    datetimeRange = dat.Epochs(CDF_TIME_TT2000=xlim).get_data('datetime')
    if round_gap is None:
        round_gap = find_best_round_gap_for_time_ticks(datetimeRange, ticks_N=ticks_N)
    xTicks = find_CDF_TIME_TT2000_xticks(datetimeRange, round_gap)
    self.set_xticks(xTicks)
mpl.axes.Axes.set_CDF_TIME_TT2000_xticks = set_CDF_TIME_TT2000_xticks

def rainbow_text_interactive(x, y, strings, colors, orientation='horizontal', ax=None, **kwargs):
    """
    This function only works for interactive backend.
    Take a list of *strings* and *colors* and place them next to each
    other, with text strings[i] being shown in colors[i].

    Parameters
    ----------
    x, y : float
        Text position in data coordinates.
    strings : list of str
        The strings to draw.
    colors : list of color
        The colors to use.
    orientation : {'horizontal', 'vertical'}
    ax : Axes, optional
        The Axes to draw into. If None, the current axes will be used.
    **kwargs
        All other keyword arguments are passed to plt.text(), so you can
        set the font size, family, etc.
    """
    if ax is None:
        ax = plt.gca()
    t = ax.transData
    canvas = ax.figure.canvas

    assert orientation in ['horizontal', 'vertical']
    if orientation == 'vertical':
        kwargs.update(rotation=90, verticalalignment='bottom')

    for s, c in zip(strings, colors):
        text = ax.text(x, y, s, color=c, transform=t, **kwargs)

        # Need to draw to update the text position.
        text.draw(canvas.get_renderer())
        ex = text.get_window_extent()
        if orientation == 'horizontal':
            t = mpl.transforms.offset_copy(text.get_transform(), x=ex.width, units='dots')
        else:
            t = mpl.transforms.offset_copy(text.get_transform(), y=ex.height, units='dots')

def rainbow_text_output(x, y, strings, colors, orientation='horizontal', ax=None, **kwargs):
    ha = kwargs.get('ha', 'left')
    if ha == 'left':
        text = ax.text(x, y, strings[0], color=colors[0], **kwargs)
        # Subsequent words, positioned with annotate(), relative to the preceding one.
        for s, c in zip(strings[1:], colors[1:]):
            text = ax.annotate(s, xycoords=text, xy=(1, 0), verticalalignment="bottom", color=c)  # custom properties
    elif ha == 'right':
        text = ax.text(x, y, strings[-1], color=colors[-1], **kwargs)
        # Subsequent words, positioned with annotate(), relative to the preceding one.
        for s, c in zip(strings[-2::-1], colors[-2::-1]):
            text = ax.annotate(s, xycoords=text, xy=(0, 0), verticalalignment="bottom", ha=ha, color=c)  # custom properties

def plot_flag_bar(fig, ax, epochs, data, epoch_fmt='CDF_TIME_TT2000', vmin=None, vmax=None, cax_width=0.02, ylabel=None, cbar_tickAndLabels=None, cbar_label=None):
    '''
    Parameters:
        cbar_tickAndLabels: a dictionary e.g. {'ticks': [0, 1], 'labels': ['SW', 'SH']}

    '''
    cax_gap = cax_width / 3
    if isinstance(epochs, dat.Epochs):
        t = epochs.get_data(epoch_fmt)
    elif isinstance(epochs, np.ndarray):
        t = epochs
    fUni = np.unique(data)
    if vmin is None:
        vmin = fUni[fUni>0][0]
    if vmax is None:
        vmax = data.max()
    tMesh, valueMesh = np.meshgrid(np.append(t, t[-1]+t[-1]-t[-2]), np.array([0, 1]), indexing='ij')
    pcm_ = ax.pcolormesh(tMesh, valueMesh, data[:, None], vmin=vmin, vmax=vmax, shading='flat', cmap='jet')
    ax.set_yticks([])
    ax_pos = ax.get_position()
    cax_pos = [ax_pos.x1+cax_gap, ax_pos.y0, cax_width, ax_pos.y1-ax_pos.y0]
    cax = fig.add_axes(cax_pos)
    cbar = fig.colorbar(pcm_, cax=cax, orientation='vertical', extend='max', label=cbar_label)
    if cbar_tickAndLabels is not None:
        cbar.set_ticks(**cbar_tickAndLabels)
    ax.set_ylabel(ylabel)

def plot_with_arrowhead(self, x, y, plot_kwargs=None, arrow_kwargs=None):
    xy = np.vstack((x, y))
    self.plot(x, y, **plot_kwargs)
    arrowStart = xy[:, -1]
    arrowVector = xy[:, -1] - xy[:, -2]
    self.arrow(x=arrowStart[0], y=arrowStart[1], dx=arrowVector[0], dy=arrowVector[1], **arrow_kwargs)
mpl.axes.Axes.plot_with_arrowhead = plot_with_arrowhead
