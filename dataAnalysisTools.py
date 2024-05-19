__author__ = 'Yufei Zhou'

import numpy as np
import cdflib    # see github.com/MAVENSDC/cdflib
import space_database_analysis.otherTools as ot
import functools
import matplotlib as mpl
from datetime import datetime
from itertools import combinations
from scipy.signal import butter, lfilter, freqz
from scipy.optimize import fsolve
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import scipy.special
import logging
import spacepy.coordinates as sppcoo


'''
<A, B> means either A or B
[A, B] means both A and B
'''

k_B = 1.380649*10**(-23) # J/Klevin
au_in_meters = 1.4959787 * 10**11
radius_Earth_in_meters = 6.378137 * 10**6
radius_Jupiter_in_meters = 6.9911 * 10**7
radius_Saturn_in_meters = 5.8232 * 10**7
radius_Mars_in_meters = 3.3895 * 10**6
radius_Mercury_in_meters = 2.4397 * 10**6
solarNorthPoleOfRotationInICRFEquatorial = np.array([63.87, 286.13])/180*np.pi # latitude and longitude
saturnianNorthPoleOfRotationInICRFEquatorial = np.array([83.54, 40.58])/180*np.pi

mass_proton = 1.6726 * 10**(-27)
mass_electron = 9.11 * 10**(-31)

eV = 1.602176634 * 10**(-19)

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff /nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5, axis=-1):
    '''
    Parameters:
        data: a numpy.ndarray object
        cutoff: the cutoff frequency
        fs: the frequency of data
    '''
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data, axis=axis)
    return y


def moving_average(data, n=1, axis=0):
    if axis != 0:
        data = np.swapaxes(data, axis, 0)
    ret = np.cumsum(data, axis=0)
    ret[n:] = ret[n:] - ret[:-n]
    ret = ret[n-1:] / n
    if axis != 0:
        ret = np.swapaxes(ret, axis, 0)
    return ret


def mvab(bVectors, returnStatisticalError=False, errorEstimation='analytical'):
    '''
    see Sonnerup and Scheible, Minimum and Maximum Variance Analysis in Analysis Methods for Multi-Spacecraft Data, ESA Publications Division, 1998
    Purpose:
        To determine along which direction the component of a vector field vary slowly, abruptly, or intermediately
    Parameter:
        bVectors: a ndarray of dimensions of (..., numberOfPoints, 3)
        returnStatisticalError: is True, ...
        errorEstimation: possible options include: analytical, bootstrap
    Return:
        eigenSystem: (eigenValues, eigenVectors). eigenValues is a ndarray of dimension (3,) whose elements are from the smallest to the greatest. eigenVectors is a ndarray (3, 3). Its second index of 3 dimensions correspond are for three eigenvectors.
    '''
    bShape = bVectors.shape
    m = np.empty([*bShape[:-2], 3, 3])
    for j in range(3):
        m[..., j, :] = np.mean(bVectors[..., :, j, None] *
                          bVectors[..., :, :], axis=-2) -\
                  np.mean(bVectors[..., :, j, None], axis=-2) *\
                  np.mean(bVectors[..., :, :], axis=-2)
    eigenSystem_ = np.linalg.eig(m)
    permutation = np.argsort(eigenSystem_[0])
    eigenValues_ = np.take_along_axis(eigenSystem_[0], permutation, axis=-1)
#    eigenValues_ = eigenSystem_[0][permutation]
    ratio = eigenValues_[..., 1]/eigenValues_[..., 0]
    eigenSystem = eigenValues_, np.take_along_axis(eigenSystem_[1], permutation[..., None, :], axis=-1)
    returnedVariables = [eigenSystem, ratio]
    if returnStatisticalError:
        numberOfPoints = len(bVectors)
        if errorEstimation == 'analytical':
            para = {'ratio': ratio, 'numberOfPoints': numberOfPoints}
            statisticalError = mvabErrorAnalysis(method='analytical', para=para)
        elif errorEstimation == 'bootstrap':
            sampleN = 1000
            sampleInds = np.floor(np.random.rand(numberOfPoints, sampleN) * numberOfPoints).astype(int)
            sampleInds
            samples = bVectors[sampleInds].swapaxes(0, 1)
            eigenSystemBootstrap, ratioBootstrap = mvab(samples)
            print(eigenSystemBootstrap[1].shape)
            minimumDirectionsBootstrap = eigenSystemBootstrap[1][..., :, 0]
            q_ = transNormal(minimumDirectionsBootstrap)
            meanDirection = normalized(np.mean(minimumDirectionsBootstrap, axis=0))
            print(meanDirection)
            statisticalError = angleBetweenVectors(meanDirection[None, :], minimumDirectionsBootstrap)
        returnedVariables.append(statisticalError)
    return returnedVariables


def mvabErrorAnalysis(method='analytical', para=None):
    '''
    Parameters:
        method: if "analytical", see BengtSonnerup1998, if "bootstrap", also see BengtSonnerup1998
    Returns:
        statisticalError: error in degrees
    '''
    if method == 'analytical':
        ratio = para['ratio']
        numberOfPoints = para['numberOfPoints']
        statisticalError = 1/(ratio-1) * np.sqrt(ratio/(numberOfPoints-1))*180/np.pi
        return statisticalError


def nfa(normalList, pos, projection=True):
    x = pos - np.mean(pos, axis=0)
    normalCenter = normalized(np.mean(normalList, axis=0))
    meanDifferenceOfNormals = 2*np.mean(np.arccos(normalList @ normalCenter[:, None]))
    nablaN = timing(normalList, x)
    if projection:
        nablaNTilde = nablaN - nablaN @ normalCenter[:, np.newaxis] @ normalCenter[None, :]
    else:
        nablaNTilde = nablaN
    eigenSystem = np.linalg.eig(nablaNTilde)
    permutation_ = np.argsort(np.abs(eigenSystem[0]))[::-1]
    print(permutation_)
    eigenValues = eigenSystem[0][permutation_]
    eigenVectors = eigenSystem[1][:, permutation_]
    errors =  np.linalg.norm(eigenVectors * eigenVectors[:, [1,2,0]], axis=0)
    return normalCenter, meanDifferenceOfNormals, eigenValues, eigenVectors, nablaNTilde, errors


def timing(normalList, xGSE, getShape=False):
    x = xGSE - np.mean(xGSE, axis=0)
    R = x.T @ x / 4
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


def gradOfVectors(vectors, x, method='Harvey'):
    if method == 'Harvey':
        numberOfSpacecrafts = vectors.shape[-2]
        shapeOfVectors = vectors.shape
        vectors = vectors.reshape((-1, numberOfSpacecrafts, 3))
        xx = x.reshape((-1, numberOfSpacecrafts, 3))
        gradOfV = np.zeros((len(vectors), 3, 3))
        for j in range(len(vectors)):
            xGSE = xx[j]
            normalList = vectors[j]
            x = xGSE - np.mean(xGSE, axis=0)
            R = x.T @ x / numberOfSpacecrafts
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
            gradOfV[j] = nablaN
        shapeOfGradV = list(shapeOfVectors[:-2])
        shapeOfGradV.extend([3, 3])
        gradOfV = gradOfV.reshape(shapeOfGradV)
    elif method == 'Shen':
        gradOfV = grad(vectors, x)
    return gradOfV


def timingVelocityAndNormal(t, xGSE, silence=False, getShape=False):
    m = timing(t, xGSE, getShape=getShape)
    if getShape:
        timingShape = m[1]
        m = m[0][0]
    print(m)
    timingVelocity = 1/np.linalg.norm(m)
    print(timingVelocity)
    timingNormal = m*timingVelocity
    timingNormal = timingNormal.squeeze()
    if silence is False:
        print("timing velocity: {:.1f}km/s,\n timing normal: {}".format(timingVelocity*6371, timingNormal))
    returnedVariables = [timingVelocity, timingNormal]
    if getShape:
        returnedVariables.append(timingShape)
    return returnedVariables


def grad(vectorLists, x, divergence0=False):
    '''see doi:10.1029/2002JA009612 Appendix B.
    Parameters:
        vectorList and x in the form of [time index, point index, cartesian index]
    Returns:
        G: G_{ij} = v_{i,j}
    '''
    numberOfPoints = x.shape[1]
    xInCenterOfMass = x - np.mean(x, axis=1)[:, None, :]
    R = np.transpose(xInCenterOfMass, (0, 2, 1)) @ xInCenterOfMass / numberOfPoints
    RInverse = np.linalg.inv(R)
    G0 = np.transpose(vectorLists, (0, 2, 1)) @ xInCenterOfMass @ RInverse / numberOfPoints
    if divergence0:
        LagrangianMultiplier = -np.trace(G0, axis1=1, axis2=2)/np.trace(RInverse, axis1=1, axis2=2)
        G = G0 + LagrangianMultiplier[:, None, None] * RInverse
    else:
        G = G0
    return G


def normalized(array, axis=-1):
    norm_ = np.linalg.norm(array, axis=axis)
    dim = len(array.shape)
    if dim == 1:
        array = array / norm_
    else:
        array = array / np.expand_dims(norm_, axis=axis)
    return  array


def align(vectorListA, vectorListB):
    'align B with A'
    mainComponent_ = np.argmax(np.abs(vectorListA).flatten())
    mainVector = mainComponent_ // 3
    mainComponent = mainComponent_ % 3
    sign0 = np.sign(vectorListA[mainVector, mainComponent])
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


## theoretical models
def kroneckerDelta(i, j):
    return 1 - np.sign(np.abs(i-j))


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


def regularTetrahedron(a=1, alpha=0, beta=0, gamma=0):
    '''
    Parameters:
        a: the length of the sides of the regular tetrahedron.
        alpha, beta, gamma are the three Euler angles, see, e.g., Section 7.1.2 of Group Theory in Physics by Wu-Ki Tung (1985).
    '''
    def R1(psi):
        return np.array([[1, 0, 0], [0, np.cos(psi), -np.sin(psi)], [0, np.sin(psi), np.cos(psi)]])
    def R2(psi):
        return np.array([[np.cos(psi), 0, np.sin(psi)], [0, 1, 0], [-np.sin(psi), 0, np.cos(psi)]])
    def R3(psi):
        return np.array([[np.cos(psi), -np.sin(psi), 0], [np.sin(psi), np.cos(psi), 0], [0, 0, 1]])
    RT_ = np.array([[0, 0, np.sqrt(6)/3*a], [np.sqrt(1/3)*a, 0, 0], [-np.sqrt(1/3)/2*a, a/2, 0], [-np.sqrt(1/3)/2*a, -a/2, 0]])
    RT = (R3(gamma)@R2(beta)@R3(alpha)@RT_.T).T
    return RT


def magnetopauseSubsolarDistance(Bz=17, Dp=3):
    '''Shue et al, 1998 doi.org/10.1029/98JA01103'''
    r0 = (10.22 + 1.29*np.tanh(0.184*(Bz + 8.14)))*Dp**(-1/6.6)
    alpha = (0.58 - 0.007*Bz)*(1 + 0.024*np.log(Dp))
    return r0, alpha


def magnetopause(theta, subsolarDistance=None, alpha=None):
    '''Shue et al, 1998 doi.org/10.1029/98JA01103'''
    r = subsolarDistance*(2/(1 + np.cos(np.pi/2)))**alpha
    return r


def dipoleField(xGSE, M=-30438):
    'xGSE in RE, B in nT'
    x1 = xGSE[..., 0][..., None]
    x2 = xGSE[..., 1][..., None]
    x3 = xGSE[..., 2][..., None]
    r = np.linalg.norm(xGSE, axis=-1)[..., None]
    return M*np.concatenate([3*x1*x3, 3*x2*x3, (3*x3**2-r**2)], axis=-1) / r**5

def corotationalElectricField(pos, M=-30348, omega=np.array([0, 0, 2*np.pi/24/3600])):
    'pos in RE, B in nT'
    velocity = np.cross(omega, pos*radius_Earth_in_meters)
    magneticField = dipoleField(pos, M=M)*10**(-9)
    e = -np.cross(velocity, magneticField)
    return e

def magnetosphericField(xGSE, M1=-30438, model="mirror", M2=-28*30438, subsolarDistance=None):
    '''xGSE in RE, B in nT, subsolarDistance should be in RE'''
    if model == "mirror":
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
    elif model == "Legendre":
        B1 = 2500
        B2 = 2100
        x1 = xGSE[..., 0][..., None]
        x2 = xGSE[..., 1][..., None]
        x3 = xGSE[..., 2][..., None]
        r1 = np.linalg.norm(xGSE, axis=-1)[..., None]
        b1 = M1*np.concatenate([3*x1*x3, 3*x2*x3, (3*x3**2-r1**2)], axis=-1) / r1**5
        b2 = (1/subsolarDistance)**3 * np.concatenate([-B2/subsolarDistance*x3, np.zeros_like(x3), B1-B2/subsolarDistance*x1], axis=-1)
        return b1+b2

##
def lundquistForceFreeField(x, B0=1, xCoordinateSystem='Cartesian', coordinateSystem='Cartesian'):
    '''
        input: 
            x, a $... \times 3$ array for the observation points in Cartesian coordinates system
            B0, a scale for the strength of the magnetic field
    '''
    if xCoordinateSystem == 'Cartesian':
        xPolar = cartesian2polar(x[..., 0:2])
        r = xPolar[..., 0]
        theta = xPolar[..., 1]
    bTheta = B0*scipy.special.j1(r)
    bZ = B0*scipy.special.j0(r)
    if coordinateSystem =='Cartesian':
        b = np.stack([-bTheta*np.sin(theta), bTheta*np.cos(theta), bZ], axis=-1)
    return b
##

def movingField(field, v, x, **para):
    '''
        input:
            field, a field model whose first input is x in Cartesian coordinates system
            v, a vector of three or four components
            x, a $... \times 4$ array for the observation points in Cartesian coordinates system.
        Note:
            this function is not complete for v
    '''
    if v.shape[-1] == 4:
        v = v[..., 1:]
    if v.shape[-1] == 3:
        xForField = x[..., 1:] - v * x[..., 0, None]
    f = field(xForField, **para)
    return f


def harrisBField(xGSE, B0=1, h=1):
    x3 = xGSE[..., 2][..., None]
    return np.concatenate([B0*np.tanh(x3/h), np.zeros_like(x3), np.zeros_like(x3)], axis=-1)


def chargedSpherePotential(r=None, x=None, rho=1, epsilon=1, a=1, ret=['potential']):
    '''
    r = distance from the center
    a = radius of the sphere
    rho = charge density of in the sphere
    epsilon = permittivity
    ret = return. potential and electricField can be returned. electricField is returned in Cartesian coordinates.
    '''
    totalCharge = 4/3*np.pi*a**3*rho
    returnedVariables = []
    for item in ret:
        if item == 'potential':
            outerPotential = 1/(4*np.pi*epsilon)*totalCharge/r
            innerPotential = -1/6/epsilon*rho*r**2 + 1/(2*epsilon)*a**2*rho
            potential = (np.sign(r-a)+1)/2*outerPotential + (np.sign(a-r)+1)/2*innerPotential
            returnedVariables.append(potential)
        elif item == 'electricField':
            r = np.linalg.norm(x, axis=-1)[..., None]
            outerCoeff = a**3*rho/3/epsilon/r**3
            outerElectricField = outerCoeff*x
            innerCoeff = rho/3/epsilon
            innerElectricField = innerCoeff*x
            electricField = (np.sign(r-a)+1)/2*outerElectricField + (np.sign(a-r)+1)/2*innerElectricField
            returnedVariables.append(electricField)
    if len(returnedVariables) == 1:
        return potential
    else:
        return returnedVariables


def chargedBall2Potential(r, b=1, epsilon=1, a=1):
    '''
    r = distance from the center
    a = radius of the sphere
    rho = charge density of in the sphere
    b = rho * r^2
    epsilon = permittivity
    '''
    const = b / epsilon
    outerPotential = const * a / r
    innerPotential = const + const * np.log(a/r)
    return (np.sign(r-a)+1)/2*outerPotential + (np.sign(a-r)+1)/2*innerPotential



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

## plot time format >>>>

def datetime2epoch(dateTime, epochType='CDF_EPOCH'):
    if epochType == 'CDF_EPOCH':
        return cdflib.cdfepoch.compute_epoch(ot.datetime2list(dateTime))
    elif epochType == 'CDF_TIME_TT2000':
        return cdflib.cdfepoch.compute_tt2000(ot.datetime2list(dateTime))
    else:
        raise Exception("datetime object resolution to microsecond")

def epoch2datetime(epoch, epochType='CDF_EPOCH'):
    if epochType == 'CDF_EPOCH':
        return datetime(*cdflib.cdfepoch.breakdown_epoch(epoch)[:6])
    elif epochType == 'CDF_TIME_TT2000':
        return datetime(*cdflib.cdfepoch.breakdown_tt2000(epoch)[:6])

class Epoch:
    '''
    This class is deprecated. use Epochs instead.
    '''
    def __init__(self, dateTime=None, epoch=None, epochType='CDF_EPOCH'):
        self.epochType = epochType
        if dateTime:
            self.dateTime = dateTime
            self.epoch = datetime2epoch(self.dateTime)
        elif epoch:
            self.epoch = epoch
            self.dateTime = epoch2datetime(self.epoch)

    def epochRecord(self, t, tolerance=1):
        '''
        Purpose: to find the index of a epoch in a time series.
        Parameters:
            t: a ndarray of time series
            tolerance (in second): the maximal allowed difference between t[record] and self.epoch
        '''
        if self.epochType == 'CDF_EPOCH':
            epoch = self.epoch
            record = np.argmin(np.abs(t - epoch))
            if np.abs(t[record] - epoch) < tolerance*1000:
                return record
            else:
                raise Exception("record not found")

class Epochs:
    '''
    Purpose:
        To facilitate the transformation between datetime, epoch, and time series record index
    Properties:
        dateTimeList: a list
        epochs: a ndarray
    '''
    def __init__(self, dateTimeList=None, epochs=None, epochType='CDF_EPOCH'):
        '''
        Parameters:
            dateTimeList: a list nested to any degree of depth whose final element is datetime object. Input either this parameter or epochs to initialize the instance.
            epochs: a list nested to any degree of depth or a ndarray.
            epochType: CDF_EPOCH, CDF_TIME_TT2000
        '''
        datetime2epochWithEpochType = functools.partial(datetime2epoch, epochType=epochType)
        epoch2datetimeWithEpochType = functools.partial(epoch2datetime, epochType=epochType)
        self.epochType = epochType
        if dateTimeList is not None:
            self.dateTimeList = dateTimeList
            self.epochs = np.array(map_multi_dimensional_list(datetime2epochWithEpochType, self.dateTimeList))
        elif epochs is not None:
            if isinstance(epochs, list):
                epochs = np.array(epochs)
            self.epochs = epochs
            self.dateTimeList = map_multi_dimensional_list(epoch2datetimeWithEpochType, epochs.tolist())

    def to_epochs(self, epochType='CDF_TIME_TT2000'):
        '''
        transform to other format of epochs
        '''
        datetimes = cdflib.epochs.CDFepoch.to_datetime(self.epochs)
        if epochType == 'CDF_TIME_TT2000':
            epochs = cdflib.epochs.CDFepoch.compute_tt2000(datetimes)
        return epochs

    def epochRecords(self, ts, tolerance=1):
        '''
        Parameters:
            ts: a ndarray [..., timeIndex]
            tolerance (in second): the maximal allowed difference between t[record] and self.epoch
        '''
        epochs = self.epochs
        records = np.argmin(np.abs(ts - epochs[..., None]), axis=-1)
        if self.epochType == 'CDF_EPOCH':
            if np.all(np.abs(ts[..., records] - epochs) < tolerance*1000):
                return records
            else:
                raise Exception("records not found")
        if self.epochType == 'CDF_TIME_TT2000':
            if np.all(np.abs(ts[..., records] - epochs) < tolerance*10**9):
                return records
            else:
                raise Exception("records not found")


def map_multi_dimensional_list(func, l):
    if hasattr(l, '__iter__') and len(l) > 0:
        if type(l[0]) != list:
            return [func(v) for v in l]
        else:
            return [map_multi_dimensional_list(func, v) for v in l]
    else:
        return []


def dataFillAndLowPass(t, data, axis=0, resamplingT=None, tDistribution='evenlySpaced', tStepPecentageCriterion=0.9, lowpassCutoff=None, gapThreshold=None, minNumberOfPoints=None, returnShiftQ=False, badpointsMask=None):
    '''
    Purposes:
        To fill missed data in a time series and filter the time series using low pass filter
    Parameters:
        t: an one-dimensional array
        data: a data array with one axis corresponding to the time series
        axis: the time axis of the data
        resamplingT: if given, the data is resampled at times in this parameter
        tDistribution: possible options include:
                evenlySpaced: return times stamps adjusted to be evenly spaced and the associated data.
                original: return the original input t and the data such that the bad points are replaced with linear interpolation.

        parameters associated with "evenlySpaced":
            tStepPecentageCriterion: the portion of major temporal gap should be larger than this parameter, otherwise the program raise exception.
            lowpassCutoff: lowpass cutoff frequency.
            gapThreshold: a number. The t is divided into blocks of t such that the gap between consecutive blocks is larger than gapThreshold.
            minNumberOfPoints: the minimum number of points in every block.
        parameters associated with "original":
            badpointsMask: an vector of true and false representing at which time stamp the data is bad.
    Return:
        if resamplingT is not None:
            return resampledData
        if resamplingT is None:
            if tDistribution == 'evenlySpaced':
                return tHomogeneous, dataProcessed
            elif tDistribution == 'original':
                return nothing, modify the input data in place

    Note:
        t is the epoch for Cluster
    '''
    if resamplingT is None:
        if tDistribution == 'evenlySpaced':
            tDiff = np.diff(t)
            if gapThreshold:
                section = 1 + np.argwhere(tDiff > gapThreshold)
                tBlocks = np.split(t, section) # each element in tBlocks is a consecutive time series in which the maximal gap ls less than gapThreshold
                dataBlocks = np.aplit(data, section, axis=axis)
            else:
                tBlocks = [t]
                dataBlocks = [data]
            processedTBlocks = []
            processedDataBlocks = []
            shiftQs = []
            for tBInd in range(len(tBlocks)):
                tBlock = tBlocks[tBInd]
                dataBlock = dataBlocks[tBInd]
                tLen = len(tBlock)
                if minNumberOfPoints:
                    if minNumberOfPoints > tLen:
                        continue
                tBlockDiff = np.diff(tBlock)
                unique, counts = np.unique(tBlockDiff, return_counts=True)
                tStep = unique[counts.argmax()]
                if len(unique) > 1:
                    if counts.max() / counts.sum() > tStepPecentageCriterion:
                        remainders = np.mod(unique, tStep)
                        uniqueRe, countsRe = np.unique(remainders, return_counts=True)
                        if all(remainders == 0):
                            shiftQ = False
                        else:
                            shiftQ = True
                        tHomogeneous = np.arange(t[0], t[-1], tStep)
                        cs = interpolate.CubicSpline(tBlock, dataBlock, axis=axis)
                        dataBlockHomogeneous = cs(tHomogeneous)
                    else:
                        print('tStep: {}'.format(tStep))
                        print('unique:')
                        print(unique)
                        print('counts:')
                        print(counts)
                        raise Exception('The time tags are irregular')
                elif len(unique) == 1:
                    shiftQ = False
                    tHomogeneous = t
                    dataBlockHomogeneous = dataBlock
                else:
                    raise Exception('Attention ! Known Problem')
                shiftQs.append(shiftQ)
                processedTBlocks.append(tHomogeneous)
                if lowpassCutoff is not None:
                    fs = 1/tStep
                    dataBlockHomogeneous = butter_lowpass_filter(dataBlockHomogeneous, lowpassCutoff, fs, axis=axis)
                processedDataBlocks.append(dataBlockHomogeneous)
            tHomogeneous = np.concatenate(processedTBlocks)
            dataProcessed = np.concatenate(processedDataBlocks, axis=axis)
            returnedVariables = [tHomogeneous, dataProcessed]
            if returnShiftQ:
                returnedVariables.append(shiftQs)
            return returnedVariables
        elif tDistribution == 'original':
            if not axis == 0:
                data = np.swapaxes(data, 0, axis)
            "First, check head and tail of data. If the first continuous set of bad points start from the first point, fill them with the value at the first good point. If the last continuous set of bad points end with the last point, fill them with the value at the last good point."
            indicesOfBadPoints = np.nonzero(badpointsMask)[0]
            indicesBreakLogical = np.diff(indicesOfBadPoints) - 1
            indicesOfBreakpointsInBadPoints = np.nonzero(indicesBreakLogical)[0]
            if badpointsMask[-1]:
                startOfLastContinuousBadPoints = indicesOfBadPoints[indicesOfBreakpointsInBadPoints[-1]+1]
                data[startOfLastContinuousBadPoints:] = data[startOfLastContinuousBadPoints-1]
                badpointsMask[startOfLastContinuousBadPoints:] = False
            if badpointsMask[0]:
                endOfFirstContinuousBadPoints = indicesOfBadPoints[indicesOfBreakpointsInBadPoints[0]]
                data[:endOfFirstContinuousBadPoints] = data[endOfFirstContinuousBadPoints+1]
                badpointsMask[:endOfFirstContinuousBadPoints] = False
            "second, linearly interpolate."
            indicesOfBadPoints = np.nonzero(badpointsMask)[0]
            indicesBreakLogical = np.diff(indicesOfBadPoints) - 1
            indicesOfBreakpointsInBadPoints = np.nonzero(indicesBreakLogical)[0]
            indicesOfTheEndOfAllContinuousBadPoints = indicesOfBadPoints[np.append(indicesOfBreakpointsInBadPoints, len(indicesOfBadPoints)-1)]
            indicesOfTheStartOfAllContinuousBadPoints = indicesOfBadPoints[np.insert(indicesOfBreakpointsInBadPoints+1, 0, 0)]
            indicesOfTheStartAndEndOfContinuousBadPoints = np.concatenate((indicesOfTheStartOfAllContinuousBadPoints[:, None], indicesOfTheEndOfAllContinuousBadPoints[:, None]), axis=1)
            linearInterpolationRange = indicesOfTheStartAndEndOfContinuousBadPoints + np.array([-1, 1])
            if len(data.shape) > 1:
                slope = (np.diff(data[linearInterpolationRange, ...], axis=1).squeeze() / np.diff(t[linearInterpolationRange], axis=1)).squeeze()
            elif len(data.shape) == 1:
                slope = (np.diff(data[linearInterpolationRange, ...], axis=1) / np.diff(t[linearInterpolationRange], axis=1)).squeeze()
            print(slope.shape)
            for i in range(len(linearInterpolationRange)):
                data[linearInterpolationRange[i, 0]:linearInterpolationRange[i, 1]] = (data[linearInterpolationRange[i, 0]] + (t[linearInterpolationRange[i, 0]:linearInterpolationRange[i, 1]] - t[linearInterpolationRange[i, 0]])[:, None]*slope[i, None, ...]).squeeze()
    else:
        logging.debug('t shape: {}'.format(t.shape))
        logging.debug('data shape: {}'.format(data.shape))
        cs = interpolate.CubicSpline(t, data, axis=axis)
        resampledData = cs(resamplingT)
        return resampledData


def fillData(t, data):
    badpointsMask = data < -10**10
    numberOfPoints = len(data)
    if len(badpointsMask.shape) > 1:
        badpointsMask.reshape((numberOfPoints, -1))
        badpointsMask = np.any(badpointsMask, axis=1).squeeze()
    numberOfBadPoints = np.count_nonzero(badpointsMask)
    logging.info("bad points number {}/{}".format(numberOfBadPoints, numberOfPoints))
    if numberOfBadPoints > 0:
        dataFillAndLowPass(t, data, axis=0, tDistribution='original', badpointsMask=badpointsMask)


def volumetricAnalysis(pos):
    '''
    Purpose:
        To calculate the characteristic directions and lengths of a distribution of points
    Parameter:
        pos: a ndarray of dimension [..., numberOfPoints, numberOfDimensions]
    '''
    x = pos - np.mean(pos, axis=-2)[..., None, :]
    numberOfPoints = x.shape[-2]
    R = x.swapaxes(-1, -2) @ x / numberOfPoints
    eigenSystemOfR = np.linalg.eig(R)
    permutation = np.argsort(eigenSystemOfR[0], axis=-1)
    eigenValues = np.take_along_axis(eigenSystemOfR[0], permutation, axis=-1)
    eigenVectors = np.take_along_axis(eigenSystemOfR[1], permutation[..., None, :], axis=-1)
    shape = (np.sqrt(eigenValues), eigenVectors)
    return shape

def alfvenSpeed(BStrength, n):
    '''
    Paramters:
        BStrength: magnetic field magnitude in nT
        n: proton number density in cm^{-3}
    Return:
        alfvenSpeed: in km/s
    '''
    alfvenSpeed = 21.812 * BStrength / np.sqrt(n)
    return alfvenSpeed


## Coordinates
def cartesian2spherical(vectors):
    '''
    Parameters:
        vectors: ndarray of shape [..., 3]
    return:
        a ndarray of shape [..., 3], the last dimension stands for r, theta and phi
    '''
    r = np.linalg.norm(vectors, axis=-1)
    unitVectors = vectors / r[..., None]
    thetaPhi = unitVectorFromCartesian2Spherical(unitVectors, halfSphere=False)
    return np.concatenate((r[..., None], thetaPhi), axis=-1)


def vectorRThetaPhi2VectorCartesian(pos, vector):
    '''
    Purpose: v_r. v_theta, v_phi to v_x, v_y, v_z
    Parameters:
        pos: an array [..., 3], representing position of the vector. The last dimension represent three components in Cartesian coordinates.
        vector: an array [..., 3]. The last dimension represent three vector components v_r, v_theta, v_phi, which have identical units.
    Return:
        vectorCartesian: an array [..., 3]
    '''
    r = np.linalg.norm(pos, axis=-1)
    x, y, z = np.moveaxis(pos, -1, 0)
    rtp2xyzMat = np.zeros(pos.shape+(3,)) 
#    rtp2xyzMat[..., 0, 0] = x/r
#    rtp2xyzMat[..., 0, 1] = y/r
#    rtp2xyzMat[..., 0, 2] = z/r
#    rtp2xyzMat[..., 1, 0] = -y/r
#    rtp2xyzMat[..., 1, 1] = x/r
#    rtp2xyzMat[..., 2, 0] = x*z/r**2
#    rtp2xyzMat[..., 2, 1] = z*y/r**2
#    rtp2xyzMat[..., 2, 2] = (-x*y-x**2)/r**2
#    vectorCartesian = np.sum(vector[..., :, None] * rtp2xyzMat, axis=-2)
    rtp2xyzMat[..., 0, 0] = x/r
    rtp2xyzMat[..., 0, 1] = y/r
    rtp2xyzMat[..., 0, 2] = z/r
    rtp2xyzMat[..., 1, 0] = x*z/r**2
    rtp2xyzMat[..., 1, 1] = y*z/r**2
    rtp2xyzMat[..., 1, 2] = -(x**2+y**2)/r**2
    rtp2xyzMat[..., 2, 0] = -y/r
    rtp2xyzMat[..., 2, 1] = x/r
    rtp2xyzMat = normalized(rtp2xyzMat)
    vectorCartesian = cartesian1ToCartesian2(vecInC1=vector, c1BasisInC2Basis=rtp2xyzMat)
    return vectorCartesian


def unitVectorFromCartesian2SphericalOld(unitVectors, halfSphere=True, twopi=True):
    '''
    Parameters:
        unitVectors: ndarray of shape [..., 3]
        halfSphere: if True, the unit vectors in opposite directions are considered the same so as to return the same theta and phi, with phi ranging from -\pi/2 to \pi/2
        twopi: if true and if halfSphere is not true, phi varies from 0 to 2\pi, else from -\pi to -\pi
    return:
        a ndarray of shape [..., 2], the last dimension stands for theta and phi
    '''
    if len(unitVectors.shape) == 1:
        unitVectors = unitVectors[None, :]
    sign0 = np.sign(unitVectors[..., 0])[..., None]
    unitVectors = unitVectors*sign0
    print(unitVectors)
    thetas = np.arccos(unitVectors[..., 2])
    print(thetas)
    phis = np.arcsin(unitVectors[..., 1]/np.sin(thetas))
    print(phis)
    if not halfSphere:
        sign0 = sign0.squeeze()
        thetas = sign0*thetas + np.pi/2*(1-sign0)
        signPhis = np.sign(phis)
        phis = (1+sign0)/2*(phis+(1-signPhis)*np.pi) + (1-sign0)/2*(-phis+np.pi+np.pi/2*signPhis)
        if not twopi:
            sign1 = np.sign(np.pi-phis)
            phis -= (1-sign1)*np.pi
    return np.stack([thetas, phis], axis=-1).squeeze()

def unitVectorFromCartesian2Spherical(unitVectors, halfSphere=True, twopi=True):
    '''
    Parameters:
        unitVectors: ndarray of shape [..., 3]
        halfSphere: if True, the unit vectors in opposite directions are considered the same so as to return the same theta and phi, with phi ranging from -\pi/2 to \pi/2
        twopi: if true and if halfSphere is not true, phi varies from 0 to 2\pi, else from -\pi to -\pi
    return:
        a ndarray of shape [..., 2], the last dimension stands for theta and phi
    '''
    if len(unitVectors.shape) == 1:
        unitVectors = unitVectors[None, :]
    theta = np.arccos(unitVectors[..., 2])
    phi_x = np.arccos(unitVectors[..., 0]/np.sin(theta))
    change = np.sign(np.sign(unitVectors[..., 1]) +0.5)
    phi = np.pi*(1-change) + change*phi_x
    phi[np.abs(unitVectors[..., 2])==1] = np.NAN
    if halfSphere:
        flipMask = np.pi/2 <= phi & phi < np.pi*3/2
        theta[flipMask] = np.pi-theta[flipMask]
        phi[flipMask] -= - np.pi
        phiFourthQuadrantMask = phi >= np.pi*3/2
        phi[phiFourthQuadrantMask] -= 2*np.pi
    else:
        if twopi:
            pass
        else:
            phi[phi>=np.pi] -= 2*np.pi
    return np.stack([theta, phi], axis=-1).squeeze()

def sphericalAngleTransform(data, coordinate=None, standardOutRange=True, inRange=None, outRange=None):
    '''
    Purpose: transform spherical anagles, theta and phi
    Parameters:
        data:
        coordinate: 'theta 'or 'phi'
        standardOutRange: if True, output theta in [0, pi] and phi in [0, 2*pi]
        inRange: clarifies the input data range
        outRange: in case of standardOutRange==False, specifies the output data range
    '''
    if coordinate == 'theta':
        if not standardOutRange:
            if inRange == 'standardOutRange':
                if all(np.abs(np.array(outRange) - np.array([-1, 1])*np.pi/2) < 10**(-5)):
                    dataTransformed = np.pi/2 - data
            else:
                dataIntermediate = sphericalAngleTransform(data, coordinate=coordinate, standardOutRange=True, inRange=inRange)
                dataTransformed = sphericalAngleTransform(dataIntermediate, coordinate=coordinate, standardOutRange=False, inRange='standardOutRange', outRange=outRange)
        else:
            if all(np.abs(np.array(inRange) - np.array([-1, 1])*np.pi/2) < 10**(-5)):
                dataTransformed = np.pi/2 - data
            elif all(np.abs(np.array(inRange) - np.array([-1, 1])*90) < 10**(-5)):
                dataTransformed = np.pi/2 - data/180*np.pi
            else:
                raise Exception('bad input range')
    if coordinate == 'phi':
        if not standardOutRange:
            if inRange == 'standardOutRange':
                if all(np.abs(np.array(outRange) - np.array([-1, 1])*np.pi) < 10**(-5)):
                    dataTransformed = data - (np.sign(data-np.pi) + 1)*np.sign(data-np.pi)*np.pi
            else:
                dataIntermediate = sphericalAngleTransform(data, coordinate=coordinate, standardOutRange=True, inRange=inRange)
                dataTransformed = sphericalAngleTransform(dataIntermediate, coordinate=coordinate, standardOutRange=False, inRange='standardOutRange', outRange=outRange)
        else:
            if all(np.abs(np.array(inRange) - np.array([-1, 1])*np.pi) < 10**(-5)):
                dataTransformed = data + (np.sign(data) - 1)*np.sign(data)*np.pi
            elif all(np.abs(np.array(inRange) - np.array([0, 24])) < 10**(-5)):
                data_ = (data - 12)/12*np.pi
                dataTransformed = data_ + (np.sign(data_) - 1)*np.sign(data_)*np.pi
            else:
                raise Exception('bad input range')
    return dataTransformed


def normalizeAngle(angle):
    '''
    purpose: to make angle within [0,2pi)
    '''
    return np.mod(angle, np.pi*2)

def cartesian2polar(vectors):
    r = np.linalg.norm(vectors, axis=-1)
    theta_ = np.arctan(vectors[..., 1]/vectors[..., 0])
    theta = theta_ + (1 - np.sign(vectors[..., 0]))/2*np.pi*np.sign(vectors[..., 1])
    return np.stack((r, theta), axis=-1)


def polar2cartesian(vectors):
    '''
    Parameters:
        vectors: can be a ndarray either of shape [..., 2] for r and theta, or of shape [..., 1] only for theta
    '''
    shape = vectors.shape
    if len(shape) == 1:
        vectors = vectors[None, :]
    if shape[-1] == 1:
        x = np.cos(vectors)
        y = np.sin(vectors)
    elif shape[-1] == 2:
        x = vectors[..., 0]*np.cos(vectors[..., 1])
        y = vectors[..., 0]*np.sin(vectors[..., 1])
    return np.stack([x, y], axis=-1).squeeze()


def spherical2cartesian(vectors):
    '''
    Parameters:
        vectors: can be a ndarray either of shape [..., 2] for theta and phi, or of shape [..., 3] for r, theta, and phi
    '''
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


def cartesian1ToCartesian2(vecInC1, c1BasisInC2Basis=None, c2BasisInC1Basis=None):
    '''
    Purpose: transform vectors between cartesian coordinates
    Parameter:
        vecInC1: an array [..., 3]
        c1BasisInC2Basis: an array [..., 3, 3]
        c2BasisInC1Basis: an array [..., 3, 3]
    '''
    if c1BasisInC2Basis is not None:
        vecInC2 = np.sum(vecInC1[..., None] * c1BasisInC2Basis, axis=-2)
    elif c2BasisInC1Basis is not None:
        vecInC2 = np.sum(vecInC1[..., None] * np.linalg.inv(c2BasisInC1Basis), axis=-2)
    return vecInC2

def c1BasisInC2Basis(c1BasisInC3Basis, c2BasisInC3Basis):
    '''
    Purpose: as indicated by the function name and the parameter names. 
    Parameter:
        c1BasisInC3Basis: an array [..., 3, 3]
        c2BasisInC3Basis: an array [..., 3, 3]
    '''
    c3BasisInC2Basis = np.linalg.inv(c2BasisInC3Basis)
    return c1BasisInC3Basis @ c3BasisInC2Basis


def unit_quaternion_spherical_linear_interpolation(tNew, tOri, q):
    '''
    Purpose:
        to interpolate the unit quaternion for coordinate transformation of vectors. see the Calibration and Measurement Algorithms Document for MMS, G-7
    Parameters:
        tNew: the x-coordinates at which to evaluate the interpolated values
        tOri: 1-D array [n], the x-coordinates of the original data
        q: an array [n, ..., 4], original data
    '''
    tNewStartInd = 0
    qNew = np.zeros((len(tNew), *q.shape[1:]))
    for tInd in range(len(tOri)-1):
        tNewToSearch = tNew[tNewStartInd:]
        t1 = tOri[tInd]
        t2 = tOri[tInd+1]
        t12 = tOri[tInd:tInd+2]
        q12 = q[tInd:tInd+2]
        tDiff = (tNewToSearch[:, None] - t12[None, :])
        mask = np.prod(tDiff * np.array([1, -1])[None, :], axis=-1, dtype=np.float64) >= 0
        t = tNewToSearch[mask]
        u = (t - t1)/(t2-t1)
#        print('len t: ', len(t))
        qNew_ = quaternionMultiply((q12[0][None, :], quaternionPower(quaternionMultiply((sppcoo.quaternionConjugate(q12[0]), q12[1])), u)))
#        tt_ = quaternionPower(quaternionMultiply(sppcoo.quaternionConjugate(q12[0]), q12[1]), u)
#        qNew_ = quaternionMultiply(q12[0][None, :], tt_)
        tNewNextStartInd = tNewStartInd + len(t)
        qNew[tNewStartInd:tNewNextStartInd] = qNew_
        tNewStartInd = tNewNextStartInd
    return qNew


def angleBetweenVectors(v1, v2):
    '''
    Purpose:
        to calculate the angle in degrees between two vectors, v1 and v2
    Parameters:
        v1 and v2 are of same shape, ndarray(..., 3)
    '''
    v1Unit = v1 / np.linalg.norm(v1, axis=-1)[..., None]
    v2Unit = v2 / np.linalg.norm(v2, axis=-1)[..., None]
    angle = np.arccos(np.sum(v1Unit * v2Unit, axis=-1))/np.pi*180
    return angle


def binomial(a, b):
    return np.math.factorial(a)//np.math.factorial(b)//np.math.factorial(a-b)


def leastSquarePolynomialApproximation(j, x, d, omega=None, regularizationMethod=None, regPara=None, solver='direct', numberOfIterations=1000):
    '''
    This Function do something
    Parameters:
        j: the sampled data of dimension [..., M, s] where M is the number of samples, s is the number of the functions
        x: the sampling nodes of dimension [..., M, r] where r is the number of variables of the sampled function
        d: the total degree of the polynomial
        omega: the weight matrix
        retularizationMethod: as implied by the name.
        solver: 'iterative' use iterative solver, 'direct' multiply the inverse of the coefficient matrix
    Return:
        if solver == 'direct', returns include:
            c: of dimension [...,h, s] where h=binomial(d+r,r), is the coefficient of polynomial approximation
    for more info see in Notes: Mathematics, subsubsection: multivariate least square approximation on canonical basis.
    '''
    r = x.shape[-1]
    M = x.shape[-2]
    vandermondeMatrix = multivariateVandermondeMatrix(x, d)
    h = vandermondeMatrix.shape[-1]
    hBlock = np.int_(scipy.special.binom(r + np.arange(d+1) - 1, r - 1))
    if omega is None:
        omega = np.identity(M)
    R = np.swapaxes(vandermondeMatrix, -1, -2) @ omega @ vandermondeMatrix
    J = np.swapaxes(vandermondeMatrix, -1, -2) @ omega @ j
    if regularizationMethod == 'Tikhonov':
        if isinstance(regPara, float):
            R = R - M * regPara * np.identity(h)
        else:
            raise Exception('bad regularization parameter for the Tikhonov method')
    if solver == 'iterative':
        c_ = np.zeros(R.shape[:-1] + (J.shape[-1],))
        if numberOfIterations is not None:
            cAllIteration = np.zeros(c_.shape + (numberOfIterations,))
            for i in range(numberOfIterations):
                c_ = lsIterate(R, J, c_, hBlock, d)
                cAllIteration[..., i] = c_
        return cAllIteration
    elif solver == 'direct':
        c = np.linalg.solve(R, J)
        return c


def leastSquarePolynomialApproximationErrorEstimation(x, d, omega=None, dj=None):
    '''
    Purpose:
        to estimate error for the function leastSquarePolynomialApproximation.
    Parameters:
        x: the sampling nodes of dimension [..., M, r] where M is the number of samples, r is the number of variables of the sampled function
        d: the total degree of the polynomial
        omega: the weight matrix
        dj: dimension [..., s] where s is the number of functions
    Return:
        if solver == 'direct', returns include:
            c: of dimension [...,h, s] where h=binomial(d+r,r), is the coefficient of polynomial approximation
    for more info see in Notes: Mathematics, subsubsection: multivariate least square approximation on canonical basis.
    '''
    r = x.shape[-1]
    M = x.shape[-2]
    vandermondeMatrix = multivariateVandermondeMatrix(x, d)
    h = vandermondeMatrix.shape[-1]
    if omega is None:
        omega = np.identity(M)
    R = np.swapaxes(vandermondeMatrix, -1, -2) @ omega @ vandermondeMatrix
    RInverse = np.linalg.inv(R)
    dg = np.linalg.norm(RInverse @ np.swapaxes(vandermondeMatrix, -1, -2) @ omega, axis=-1)[..., None] * dj[..., None, :]
    return dg


def lsIterate(R, J, c, hBlock, d):
    h = c.shape[-1]
    for i in range(d+1):
        hp = np.sum(hBlock[:i])
        hd = hBlock[i]
        b = J[..., hp:hp+hd, :] - R[..., hp:hp+hd, :hp] @ c[..., :hp, :] - R[..., hp:hp+hd, hp+hd:] @ c[..., hp+hd:, :]
        A = R[..., hp:hp+hd, hp:hp+hd]
        x = np.linalg.solve(A, b)
        c[..., hp:hp+hd, :] = x
    return c

def multivariateVandermondeMatrix(x, degree):
    '''
    Parameters:
        x: the positions, an array of dimension [..., M, r] where M is the number of samples and r is the number of variables
        degree: the highest degree of the Vandermonde matrix
    '''
    r = x.shape[-1]
    powerList = []
    for i in range(degree + 1):
       powerList.extend(polynomialSpaceBasisDegree(i, r))
    power = np.array(powerList)
    vandermondeMatrix = np.prod(x[..., None, :]**power, axis=-1)
    return vandermondeMatrix

def polynomialSpaceBasisDegree(d, r):
    '''
    This Function produce the possible multi-indices of total degree d in r-dimensional space.
    '''
    comb = []
    rn = r - 1
    for i in range(d+1):
        if i == 0:
            comb.append([d] + [0]*rn)
            if rn == 0:
                break
        elif i > 0:
            if rn > 0:
                combn = polynomialSpaceBasisDegree(i, rn)
                for item in combn:
                    item.insert(0, d-i)
                comb.extend(combn)
    return comb


def evaulatePolynomial(d, c, x):
    vandermondeMatrix = multivariateVandermondeMatrix(x, d)
    p = np.sum(vandermondeMatrix * c[None, :], axis=-1)
    return p


def polynomialInterpolation(j, x, d, selection=None):
    '''
    This Function do something
    Parameters:
        j: the sampled data of dimension [..., M] where M is the number of samples
        x: the sampling nodes of dimension [..., M, r] where r is the number of variables of the sampled function
        d: the total degree of the polynomial
        selection: the method for selection if more points are available than needed
    '''
    r = x.shape[-1]
    M = x.shape[-2]
    h= np.int_(scipy.special.binom(r+ d, r))
    if M < h:
        raise Exception('Sampled data not enough!')
    if M > h:
        if selection is None:
            raise Exception('Too many sampled data!')
        elif selection == 'order':
            x = x[..., :h, :]
            j = j[..., :h]
        elif selection == 'leja':
            pass
    vandermondeMatrix = multivariateVandermondeMatrix(x, d)
    c = np.linalg.solve(vandermondeMatrix, j)
    return c


def makeBlocks(data, breakPoints, axis=0):
    '''
    Purpose: break data into blocks, each breakPoints is the start of a block.
    '''
    if axis != 0:
        data = data.swapaxes(axis, 0)
    blocks = []
    blocks.append(data[:breakPoints[0]])
    for breakPointsInd in range(len(breakPoints)-1):
        blocks.append(data[breakPoints[breakPointsInd]:breakPoints[breakPointsInd+1]])
    blocks.append(data[breakPoints[-1]:])
    if axis != 0:
        for blockInd in range(len(blocks)):
            blocks[blockInd] = blocks[blockInd].swapaxes(axis, 0)
    return blocks


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
        dataUpBlocks = makeBlocks(data, upArgs+1)
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
            tOfInterest = Epoch(dateTime=datetime(2002, 2, 9, 0, 34, 0)).epoch
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
                dataDownBlocks = makeBlocks(dataUpBlock, downArgs+1)
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
    vThetaTable_ = sphericalAngleTransform(np.radians(theta), coordinate='theta', standardOutRange=True, inRange=[-np.pi/2, np.pi/2])
    vPhiTable_ = sphericalAngleTransform(np.radians(phi), coordinate='phi', standardOutRange=True, inRange=[-np.pi, np.pi])
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
            ax.set_title(datetime(*cdflib.cdfepoch.breakdown(t[tInd])[:6]))
    return fig


def interpolatePhaseSpaceDensity(phaseSpaceDensity, vTable, vThetaTable, vPhiTable):
    phaseSpaceDensity = phaseSpaceDensity.flatten()
    vPointsMeshgridSpherical = np.meshgrid(vTable, vThetaTable, vPhiTable, indexing='ij')
    vPointsCartesian = spherical2cartesian(np.stack(vPointsMeshgridSpherical, axis=-1).reshape(-1, 3))
    interp = interpolate.NearestNDInterpolator(vPointsCartesian, phaseSpaceDensity)
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

def interp(x, xp, fp):
    '''
    Purpose: linear interpolation
    Parameters:
        x: the x-coordinates at which to evaluate the interpolated values
        xp: 1-D array [n], the x-coordinates of the original data
        fp: an array [n, ...], original data
    '''
    dataShape = fp.shape
    data = fp.reshape([dataShape[0], -1])
    data_interpolated = np.zeros([len(x), data.shape[-1]])
    for ind in range(data.shape[-1]):
        data_interpolated[:, ind] = np.interp(x, xp, data[:, ind])
    return data_interpolated

def rtnBasisInJSOBasisFromHorizonsData(t, **posDataDict):
    '''
    Parameters:
        posDataDict: e.g. = {
                'tJupiter': tJupiter,
                'vJupiter': vJupiter,
                'posCartesianJupiter': posCartesianJupiter,
                'tSun': tSun,
                'vSun': vSun,
                'posCartesianSun': posCartesianSun,
                'tVoyager': tVoyager,
                'posCartesianVoyager': posCartesianVoyager,
                       }
    '''
    solarEquatorInICRFEcliptic = np.array([7.25, 75.76])/180*np.pi
    solarPoleInICRFEclipticSpherical = np.array([solarEquatorInICRFEcliptic[0], np.mod(solarEquatorInICRFEcliptic[1]-np.pi/2, np.pi*2)])
    solarPoleInICRFEclipticCartesian = spherical2cartesian(solarPoleInICRFEclipticSpherical)
    vJupiter = interp(t, posDataDict['tJupiter'], posDataDict['vJupiter'])
    vSun = interp(t, posDataDict['tSun'], posDataDict['vSun'])
    posCartesianVoyager = interp(t, posDataDict['tVoyager'], posDataDict['posCartesianVoyager'])
    posCartesianJupiter = interp(t, posDataDict['tJupiter'], posDataDict['posCartesianJupiter'])
    posCartesianSun = interp(t, posDataDict['tSun'], posDataDict['posCartesianSun'])
    vJupiterOrbital = vJupiter - vSun
    vJupiterOrbitalNormalized = normalized(vJupiterOrbital)
    x_axis = normalized(posCartesianSun-posCartesianJupiter)
    z_axis = np.cross(vJupiterOrbitalNormalized, x_axis)
    y_axis = np.cross(z_axis, x_axis)
    jsoBasisInICRFBasis = np.concatenate([x_axis[..., None, :], y_axis[..., None, :], z_axis[..., None, :]], axis=-2)
    r_axis_rtn = normalized(posCartesianVoyager - posCartesianSun)
    t_axis_rtn = normalized(np.cross(solarPoleInICRFEclipticCartesian, r_axis_rtn))
    n_axis_rtn = normalized(np.cross(r_axis_rtn, t_axis_rtn))
    rtnBasisInICRFBasis = np.concatenate([r_axis_rtn[..., None, :], t_axis_rtn[..., None, :], n_axis_rtn[..., None, :]], axis=-2)
    rtnBasisInJSOBasis = c1BasisInC2Basis(rtnBasisInICRFBasis, jsoBasisInICRFBasis)
    return rtnBasisInJSOBasis
#
def jsoBasisInICRFBasisFromHorizonsData(t, **posDataDict):
    '''
    Parameters:
        posDataDict: e.g. = {
                'tJupiter': tJupiter,
                'vJupiter': vJupiter,
                'posCartesianJupiter': posCartesianJupiter,
                'tSun': tSun,
                'vSun': vSun,
                'posCartesianSun': posCartesianSun,
                'tVoyager': tVoyager,
                'posCartesianVoyager': posCartesianVoyager,
                       }
    '''
    vJupiter = interp(t, posDataDict['tJupiter'], posDataDict['vJupiter'])
    vSun = interp(t, posDataDict['tSun'], posDataDict['vSun'])
    posCartesianVoyager = interp(t, posDataDict['tVoyager'], posDataDict['posCartesianVoyager'])
    posCartesianJupiter = interp(t, posDataDict['tJupiter'], posDataDict['posCartesianJupiter'])
    posCartesianSun = interp(t, posDataDict['tSun'], posDataDict['posCartesianSun'])
    vJupiterOrbital = vJupiter - vSun
    vJupiterOrbitalNormalized = normalized(vJupiterOrbital)
    jupiterToVoyager = (posCartesianVoyager - posCartesianJupiter)*1000/(radius_Jupiter_in_meters)
    x_axis = normalized(posCartesianSun-posCartesianJupiter)
    z_axis = np.cross(vJupiterOrbitalNormalized, x_axis)
    y_axis = np.cross(z_axis, x_axis)
    jsoBasisInICRFBasis = np.concatenate([x_axis[..., None, :], y_axis[..., None, :], z_axis[..., None, :]], axis=-2)
    posVoyagerInJSO = cartesian1ToCartesian2(jupiterToVoyager, c2BasisInC1Basis=jsoBasisInICRFBasis)
    return posVoyagerInJSO

def ksmBasisInICRFBasis(t, tSaturn, posCartesianSaturn, tSun, posCartesianSun):
    '''
    see doi: 10.1029/2010JA016349 for KSM definition
    '''
    posCartesianSaturn = interp(t, tSaturn, posCartesianSaturn)
    posCartesianSun = interp(t, tSun, posCartesianSun)
    magneticAxisSpherical = np.array([np.pi/2, 0]) + np.array([-1, 1])*saturnianNorthPoleOfRotationInICRFEquatorial
    magneticAxisCartesian = spherical2cartesian(magneticAxisSpherical)
    shapeOfPoints = posCartesianSun.shape
    magneticAxisCartesian = np.repeat(magneticAxisCartesian[None, :], np.prod(shapeOfPoints[:-1]), axis=0).reshape(shapeOfPoints)
    planetToSun = posCartesianSun-posCartesianSaturn
    x_axis = normalized(planetToSun, axis=-1)
    y_axis = normalized(np.cross(magneticAxisCartesian, x_axis))
    z_axis = normalized(np.cross(x_axis, y_axis))
    ksmBasisInICRFBasis = np.concatenate([x_axis[..., None, :], y_axis[..., None, :], z_axis[..., None, :]], axis=-2)
    return ksmBasisInICRFBasis

def ssqBasisInICRFBasis(t, tSaturn, posCartesianSaturn, tSun, posCartesianSun):
    '''
    see doi: 10.1002/2018JA025214 for SSQ definition
    '''
    posCartesianSaturn = interp(t, tSaturn, posCartesianSaturn)
    posCartesianSun = interp(t, tSun, posCartesianSun)
    magneticAxisSpherical = np.array([np.pi/2, 0]) + np.array([-1, 1])*saturnianNorthPoleOfRotationInICRFEquatorial
    magneticAxisCartesian = spherical2cartesian(magneticAxisSpherical)
    shapeOfPoints = posCartesianSun.shape
    magneticAxisCartesian = np.repeat(magneticAxisCartesian[None, :], np.prod(shapeOfPoints[:-1]), axis=0).reshape(shapeOfPoints)
    planetToSun = posCartesianSun-posCartesianSaturn
    z_axis = normalized(magneticAxisCartesian)
    y_axis = normalized(np.cross(magneticAxisCartesian, normalized(planetToSun)))
    x_axis = normalized(np.cross(y_axis, z_axis))
    ssqBasisInICRFBasis = np.concatenate([x_axis[..., None, :], y_axis[..., None, :], z_axis[..., None, :]], axis=-2)
    return ssqBasisInICRFBasis

def plotCartesianVectorTimeSeries(ax, t, data, norm=True, label=None):
    labelsF = ['${}_x$', '${}_y$', '${}_z$']
    colors = ['m', 'g', 'b']
    if label:
        labels = [labelF.format(label) for labelF in labelsF]
    for ind in range(3):
        ax.plot(t, data[:, ind], color=colors[ind], label=labels[ind])
    if norm:
        ax.plot(t, np.linalg.norm(data, axis=-1), color='k', label=label)
    ax.set_ylabel(label)
    ax.grid(True)
    return ax

def tetrahedronQuality(pos, method='qGM'):
    '''
    see Section 13.3 of Analysis Methods for Multi-Spacecraft Data for the methods
    '''
    x = pos - np.mean(pos, axis=-2)[..., None, :]
    numberOfPoints = x.shape[-2]
    R = x.swapaxes(-1, -2) @ x / numberOfPoints
    trueVolume = np.sqrt(np.linalg.det(R))*8/3
    combs = combinations(range(numberOfPoints), 3)
    numberOfSurfaces = combinationCounts(numberOfPoints, 3)
    surfaces = np.zeros((*x.shape[:-2], numberOfSurfaces, 3, 3))
    for i, comb in enumerate(list(combs)):
        surfaces[..., i, :, :] = x[..., comb, :]
    areaOfSurfaces = surfaceOfThreePoints(surfaces)
    trueSurface = np.sum(areaOfSurfaces, axis=-1)

    combs = combinations(range(numberOfPoints), 2)
    numberOfSides = combinationCounts(numberOfPoints, 2)
    sides = np.zeros((*x.shape[:-2], numberOfSides, 2, 3))
    for i, comb in enumerate(list(combs)):
        sides[..., i, :, :] = x[..., comb, :]
    averageLength = np.mean(np.linalg.norm(np.diff(sides, axis=-2).squeeze(), axis=-1), axis=-1)
    idealSurface = averageLength**2 *np.sqrt(3)/4 * numberOfSurfaces
    idealVolume = averageLength**3 /6/np.sqrt(2)
    if method == 'qGM':
        qualityParameter = trueVolume/idealVolume + trueSurface/idealSurface +1
    elif method == 'ellipsoidalShape':
        eigenSystemOfR = np.linalg.eig(R)
        timingShape = eigenSystemOfR
#        permutation = np.argsort(eigenSystemOfR[0])
#        timingShape = (np.sqrt(eigenSystemOfR[0])[permutation], eigenSystemOfR[1][:, permutation])
        qualityParameter = timingShape
    return qualityParameter


def surfaceOfThreePoints(pos):
    vec = np.diff(pos, axis=-2)
    surface = 0.5*np.linalg.norm(np.cross(vec[..., 0, :], vec[..., 1, :]), axis=-1)
    return surface

def combinationCounts(totalN, combN):
    return np.math.factorial(totalN)//np.math.factorial(combN)//np.math.factorial(totalN-combN)


## functions below are for producing data used for a study. These will be transfer to some other standalone tool file
def dumpAxData(cdfFile, ax):
        ax_children = ax.get_children()
        for childInd, ax_child in enumerate(ax_children):
            varName = 'object_{}'.format(len(cdfFile.zvars))
            varData = None
            if type(ax_child) == mpl.lines.Line2D:
                varData = ax_child.get_xydata()
                var_attrs = {
                        'type': str(type(ax_child)),
                        'data_description': 'abscissa and ordinate of the curve',
                        'label': ax_child.get_label(),
                        }
            elif type(ax_child) == mpl.patches.FancyArrow:
                varData = np.array([ax_child._x, ax_child._y, ax_child._dx, ax_child._dy])
                var_attrs = {
                        'type': str(type(ax_child)),
                        'data_description': 'arrow x, y, dx, dy',
                        }
            elif type(ax_child) == mpl.text.Text:
                varData = np.array(ax_child.get_position())
                var_attrs = {
                        'type': str(type(ax_child)),
                        'data_description': 'position of the text',
                        'text': ax_child.get_text(),
                        }
            varDType = 'CDF_DOUBLE'
            if varData is not None:
                dumpPlotData(cdfFile=cdfFile, varData=varData, varName=varName, varDType=varDType, dim_size=[], num_elements=1, rec_vary=True, var_attrs=var_attrs)


def dumpPlotData(cdfFile, varData=None, varName=None, varDType=None, dim_size=[], num_elements=1, rec_vary=True, var_attrs={}):
    if not dim_size:
        dim_size = varData.shape[1:]
    varDType = cdflib.cdfwrite.CDF._datatype_token(varDType)
    var_spec = {
            'Variable': varName,
            'Data_Type': varDType,
            'Num_Elements': num_elements,
            'Rec_Vary': rec_vary,
            'Dim_Sizes': dim_size,
                }
    cdfFile.write_var(var_spec, var_attrs=var_attrs, var_data=varData)

def getAxesInfoDict(ax):
    axInfoDict = {
            'panel_title': ax.get_title(),
            'panel_x_label': ax.get_xlabel(),
            'panel_y_label': ax.get_ylabel(),
            }
    return axInfoDict

def get_twin(ax, axis):
    assert axis in ("x", "y")
    siblings = getattr(ax, f"get_shared_{axis}_axes")().get_siblings(ax)
    twins = []
    for sibling in siblings:
        if sibling.bbox.bounds == ax.bbox.bounds and sibling is not ax:
            twins.append(sibling)
    return twins

## quaternion
def quaternionPower(q, p, scalarPos='last'):
    '''
    Parameters:
        q: quaternions of dimension [..., 4]
        p: power
    '''
    if scalarPos == 'last':
        halfTheta = np.arccos(q[..., -1])
        qPoweredVector = q[..., :-1] / np.sin(halfTheta) * np.sin(halfTheta*p)[:, None]
        qPoweredScalar = np.cos(halfTheta*p)
        qPowered = np.concatenate((qPoweredVector, qPoweredScalar[..., None]), axis=-1)
    return qPowered

def quaternionMultiply(qs, scalarPos='last'):
    '''
    Purpose:
        see the function name
    Parameters:
        qs: a list of q
        q1: quaternions of dimension [..., 4]
        q2: quaternions of dimension [..., 4]
    '''
    qs = list(qs)
    if len(qs) == 2:
        q1, q2 = qs
        scalar = np.array([q1[..., -1] * q2[..., -1] - np.sum(q1[..., :-1] * q2[..., :-1], axis=-1)]).swapaxes(0, -1)
        vector = np.cross(q1[..., :-1], q2[..., :-1]) + q1[..., -1, None] * q2[..., :-1] + q2[..., -1, None] * q1[..., :-1]
        return np.concatenate((vector, scalar), axis=-1)
    elif len(qs) > 2:
       q2 = qs.pop()
       q1 = qs.pop()
       q = quaternionMultiply((q1, q2))
       while len(qs) > 0:
           q2 = q
           q1 = qs.pop()
           q = quaternionMultiply((q1, q2))
       return q
    else: raise Exception

#np.array([np.arange(6).reshape((2,3))]).swapaxes(0, -1)
