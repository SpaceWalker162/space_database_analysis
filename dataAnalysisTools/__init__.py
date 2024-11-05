__author__ = 'Yufei Zhou'

import numpy as np
import cdflib    # see github.com/MAVENSDC/cdflib
import space_database_analysis.otherTools as ot
import space_database_analysis.databaseUserTools as dut
from . import constants as const
import functools
#from datetime import datetime
import datetime as dt
from astropy.time import Time, TimeDelta
from astropy import units as u
from itertools import combinations
from scipy.signal import butter, lfilter, freqz
from scipy.optimize import fsolve
import scipy.interpolate as interpolate
import scipy.signal as signal
import scipy.special
import logging
import spacepy.coordinates as sppcoo


'''
<A, B> means either A or B
[A, B] means both A and B
'''

class mms:
    def __init__(self, spacecrafts=None, workDataDir=None, workDataDirsBak=None):
        if spacecrafts:
            self.spacecrafts = spacecrafts
        else:
            missionName = 'mms'
            spacecraftNames = ['mms'+str(i) for i in range(1,5)]
            spacecrafts = []
            for spacecraftName in spacecraftNames:
                spacecraft = dut.Spacecraft(mission=missionName, name=spacecraftName, workDataDir=workDataDir, workDataDirsBak=workDataDirsBak)
                spacecrafts.append(spacecraft)
            self.spacecrafts = spacecrafts
        self.data = {}

    def _load_edp_mec(self, datetimeRange, mode='fast', **kwargs):
        try:
            quality_from = kwargs.pop('quality_from')
        except: quality_from = None
        for spacecraft in self.spacecrafts:
            spacecraftName = spacecraft.name
            if quality_from is None:
                quality_from = 'quality'
            if quality_from == 'quality':
                quality_var_ret = (spacecraftName+'_edp_quality_'+mode+'_l2', 'quality')
            elif quality_from == 'bitmask':
                quality_var_ret = (spacecraftName+'_edp_bitmask_'+mode+'_l2', 'bitmask')
            datasets_variables_with_retrieving_names = {
                (spacecraftName+'_edp_'+mode+'_l2_dce').upper(): ['EDP', (spacecraftName+'_edp_epoch_'+mode+'_l2', 't'), (spacecraftName+'_edp_dce_gse_'+mode+'_l2', 'E'), (spacecraftName+'_edp_dce_dsl_'+mode+'_l2', 'E_dsl'), (spacecraftName+'_edp_dce_err_'+mode+'_l2', 'E_err'), quality_var_ret],
                (spacecraftName+'_mec_srvy_l2_epht89q').upper(): ['MEC89Q', ('Epoch', 't'), (spacecraftName+'_mec_r_gse', 'xGSE'), (spacecraftName+'_mec_v_gse', 'vGSE'), (spacecraftName+'_mec_quat_eci_to_dsl', 'quat_eci_to_dsl'), (spacecraftName+'_mec_quat_eci_to_gse', 'quat_eci_to_gse')],
                }
            spacecraft.loadData(datetimeRange=datetimeRange, datasets_variables_with_retrieving_names=datasets_variables_with_retrieving_names, **kwargs)
            source = 'EDP'
            if quality_from == 'bitmask':
                bitmask = ot.decimal2binaryArray(spacecraft.data[source]['bitmask'], order='<')
                expected_bitmask_size = len(dut._mms_edp_bitmask)
                if bitmask.shape[-1] < expected_bitmask_size:
                    bitmask_ = np.zeros((*bitmask.shape[:-1], expected_bitmask_size))
                    bitmask_[..., :bitmask.shape[-1]] = bitmask
                    bitmask = bitmask_
                quality_flags = np.min(np.where(bitmask, dut._mms_edp_bitmask, 3), axis=-1)
                quality_mask = quality_flags > 1
            elif quality_from == 'quality':
                quality_mask = spacecraft.data[source]['quality'] > 1
            spacecraft.data[source] = mask_dict_of_ndarray(spacecraft.data[source], quality_mask)
            quat_dsl_to_gse = sppcoo.quaternionMultiply(sppcoo.quaternionConjugate(spacecraft.data['MEC89Q']['quat_eci_to_dsl']), spacecraft.data['MEC89Q']['quat_eci_to_gse'])
            spacecraft.data['MEC89Q']['quat_dsl_to_gse'] = quat_dsl_to_gse

    def chargeCalculation(self, datetimeRange, mode='fast', **kwargs):
        self._load_edp_mec(datetimeRange=datetimeRange, mode=mode, **kwargs)
        data = self.chargeCalculationWithDataLoaded(self.spacecrafts)
        self.data.update(data)


    @staticmethod
    def chargeCalculationWithDataLoaded(spacecrafts, error_estimation=True):
        numberOfSpacecrafts = len(spacecrafts)
        source = 'EDP'
        t_list = []
        data_list = []
        for spacecraft in spacecrafts:
            t_ = spacecraft.data[source]['t']
            if error_estimation:
                data_ = np.concatenate((spacecraft.data[source]['E'], spacecraft.data[source]['E_err']), axis=-1)
            else:
                data_ = spacecraft.data[source]['E']
            mask_ = ~np.isnan(data_).any(axis=-1)
            data_list.append(data_[mask_])
            t_list.append(t_[mask_])
        resamplingT, synchronized_data = data_synchronization(data_list, t_list)
        EAllSpacecrafts = synchronized_data[..., :3]
        if error_estimation:
            EErrAllSpacecrafts = synchronized_data[..., 3:]

        source = 'MEC89Q'
        dataAllSpacecrafts = np.zeros((len(resamplingT), numberOfSpacecrafts, 3))
        for indOfSC in range(numberOfSpacecrafts):
            data_ = spacecrafts[indOfSC].data[source]['xGSE']
            t_ = spacecrafts[indOfSC].data[source]['t']
            dataAllSpacecrafts[:, indOfSC, :] = dataFillAndLowPass(t_, data_, resamplingT=resamplingT)
        xGSEAllSpacecrafts = dataAllSpacecrafts
        xGSEConstellationCenter = np.mean(xGSEAllSpacecrafts, axis=-2)
        xGSEAllSpacecraftsInConstellationCenter = xGSEAllSpacecrafts - xGSEConstellationCenter[..., None, :]
        print('approximating with input data shape {}...'.format(EAllSpacecrafts.shape))
        coeff = leastSquarePolynomialApproximation(j=EAllSpacecrafts, x=xGSEAllSpacecraftsInConstellationCenter, d=1, omega=None, regularizationMethod=None, regPara=None, solver='direct')
        print('approximated')
        timingShape = volumetricAnalysis(xGSEAllSpacecrafts)
        gradE = coeff[:, 1:, :]
        rho = np.trace(gradE, axis1=-1, axis2=-2) * 55.2636 # in unit of e/m^3
        data = {}
        data['rho'] = rho
        data['resamplingT'] = resamplingT
        data['timingShape'] = timingShape[0][:, -1]/timingShape[0][:, 0]
        data['x_gse'] = xGSEConstellationCenter
        if error_estimation:
            EErrAverageOverAllSpacecrafts = np.mean(EErrAllSpacecrafts, axis=-2)
            print('estimating errors...')
            gradEErr = leastSquarePolynomialApproximationErrorEstimation(x=xGSEAllSpacecraftsInConstellationCenter, d=1, omega=None, dj=EErrAverageOverAllSpacecrafts)
            print('estimated')
            rhoErr = np.linalg.norm(np.diagonal(gradEErr[:, 1:], axis1=-2, axis2=-1), axis=-1) * 55.2636 # in unit of e/m^3
            data['rhoErr'] = rhoErr
        return data

    def plot_PSD_spectrogram_from_partial_numberdensity(self, ax, mms_ind=1, datetimeRange=None):
        source = 'FPI_DIS-PARTMOMS'
        t =  self.spacecrafts[mms_ind-1].data[source]['t']
        if t is not None:
            energyTable_ = spacecraft.data[source]['E'][0]
            energy_delta = spacecraft.data[source]['E_delta'][0]
            energyTable = energyTable_ + energy_delta
            partial_numberdensity = spacecraft.data[source]['f(E)'] * 10**6 # in m^-3

class maven:

    @staticmethod
    def swia_calculate_moments_from_3d(datetimeRange, workDataDir, workDataDirsBak=None, copy_if_not_exist=False, search_online=False):
        spacecraft = dut.Spacecraft(mission='maven', name='maven', workDataDir=workDataDir, workDataDirsBak=workDataDirsBak)
        datasets_variables_with_retrieving_names = {
            'MVN_SWI_L2_COARSESVY3D': ['SWI-dist', ('epoch', 't'), ('atten_state', 'atten_state'), ('diff_en_fluxes', 'diff_en_fluxes'), ('energy_coarse', 'energy'), ('theta_coarse', 'theta'), ('theta_atten_coarse', 'theta_atten'), ('phi_coarse', 'phi')],
            }
        spacecraft.loadData(datetimeRange=datetimeRange, datasets_variables_with_retrieving_names=datasets_variables_with_retrieving_names, copy_if_not_exist=copy_if_not_exist, search_online=search_online)

        source = 'SWI-dist'
        t_3d = spacecraft.data[source]['t']
        atten_state = spacecraft.data[source]['atten_state']
        def_ion = np.swapaxes(spacecraft.data[source]['diff_en_fluxes'], 1, 3) / (u.cm**2*u.s)
        energy = spacecraft.data[source]['energy']
        theta = spacecraft.data[source]['theta']
        theta_atten_closed = spacecraft.data[source]['theta_atten']
        phi = spacecraft.data[source]['phi']

        theta_spherical = sphericalAngleTransform(theta, coordinate='theta', inRange=[-90, 90])
        phi_spherical = sphericalAngleTransform(phi, coordinate='phi', inRange=[0, 360])
        phi_half_window = normalizeAngle(np.diff(phi_spherical))/2
        phi_bound = np.zeros((phi_spherical.shape[0]+1, *phi_spherical.shape[1:]))
        phi_bound[1:-1] = phi_half_window +  phi_spherical[:-1]
        phi_bound[0] = phi_spherical[0] - phi_half_window[0]
        phi_bound[-1] = phi_spherical[-1] + phi_half_window[-1]
        phi_bound = normalizeAngle(phi_bound)
        phi_bin_width = normalizeAngle(np.diff(phi_bound))
        phi_bound = np.insert(np.cumsum(normalizeAngle(np.diff(phi_bound))), 0, 0) + phi_bound[0]

        theta_half_window = np.diff(theta_spherical, axis=0)/2
        theta_bound = np.zeros((theta_spherical.shape[0]+1, *theta_spherical.shape[1:]))
        theta_bound[1:-1] = theta_half_window +  theta_spherical[:-1]
        theta_bound[0] = theta_spherical[0] - theta_half_window[0]
        theta_bound[-1] = theta_spherical[-1] + theta_half_window[-1]

        solid_angle = np.moveaxis(np.abs(-(np.cos(theta_bound[1:]) - np.cos(theta_bound[:-1]))[:, None] * phi_bin_width[None, :, None]), -1, 0)

        energy_log = np.log(energy)
        energy_log_half_window = np.diff(energy_log)/2
        energy_log_bound = np.zeros((energy_log.shape[0]+1, *energy_log.shape[1:]))
        energy_log_bound[1:-1] = energy_log_half_window +  energy_log[:-1]
        energy_log_bound[0] = energy_log[0] - energy_log_half_window[0]
        energy_log_bound[-1] = energy_log[-1] + energy_log_half_window[-1]
        energy_bound = np.e**energy_log_bound
        energy_width = -np.diff(energy_bound)
        speed =  np.sqrt(energy*u.eV * 2 / const.m_p)
        v_spherical_coord = np.zeros((len(energy), len(theta), len(phi), 3))
        v_spherical_coord[:, :, :, 0] = speed.to('cm/s')[:, None, None]
        v_spherical_coord[:, :, :, 2] = phi[None, None, :]
        for thetaInd in range(len(theta)):
            v_spherical_coord[:, thetaInd, :, 1] = theta[thetaInd, :, None]
        v_cartesian_coord = spherical2cartesian(v_spherical_coord) * u.cm/u.s
        speed_cartesian_coord = np.linalg.norm(v_cartesian_coord, axis=-1)
        factor_n = (1/speed*energy_width/energy)[None, :, None, None] * solid_angle[None, :, :, :] * def_ion
        density = np.sum(factor_n, axis=(1, 2, 3)).to('cm^-3')

        angle_factors = np.swapaxes(angle_factors_in_integrating_1st_moment(theta_bound, phi_bound), axis1=0, axis2=1)[None, ...]
        factor_v = ((energy_width/energy)[None, :, None, None] * def_ion)[..., None] * angle_factors
        bulk_velocity = np.sum(factor_v, axis=(1, 2, 3))/density[:, None]

        temperature_tensor_minus_bulkv = []
        for ind in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]:
            temperature_tensor_minus_bulkv.append(density[:, None] * (bulk_velocity[:, ind[0]] * bulk_velocity[:, ind[1]])[:, None])
        temperature_tensor_minus_bulkv = np.concatenate(temperature_tensor_minus_bulkv, axis=-1)
        angle_factors = np.swapaxes(angle_factors_in_integrating_2nd_moment(theta_bound, phi_bound), axis1=0, axis2=1)
        factor_temperature = ((speed * energy_width/energy)[None, :, None, None] * def_ion)[..., None] * angle_factors[None, ...]
        temperature_tensor = (const.m_p * (np.sum(factor_temperature, axis=(1, 2, 3)) - temperature_tensor_minus_bulkv) /density[:, None]).to('eV')
        temperature = np.mean(temperature_tensor[:, [0, 3, 5]], axis=-1)
        return t_3d, density, bulk_velocity, temperature


def synchronize_data_from_spacecraft_for_gradient_calculation(spacecrafts, data_source='FGM', data_name='B', x_source='MEC89Q', x_name='x'):
    '''
    synchronize field data and position data
    Parameters:
        pass
    Return:
        resamplingT: 1d array the epochs for the resulted data
        data_array: of size [len(resamplingT), len(spacecrafts), original_data.shape[1:]]
        x_array: of size [len(resamplingT), len(spacecrafts), original_x.shape[1:]]
    '''
    t_data_list = [sc.data[data_source]['t'] for sc in spacecrafts]
    data_list = [sc.data[data_source][data_name] for sc in spacecrafts]
    t_x_list = [sc.data[x_source]['t'] for sc in spacecrafts]
    x_list = [sc.data[x_source][x_name] for sc in spacecrafts]
    resamplingT, resampled_data_list = data_synchronization(data_list+x_list, t_data_list+t_x_list, resamplingT=t_data_list[0])
    data_array = np.concatenate([data[:, None, ...] for data in resampled_data_list[:len(spacecrafts)]], axis=1)
    x_array = np.concatenate([data[:, None, ...] for data in resampled_data_list[len(spacecrafts):]], axis=1)
    return resamplingT, data_array, x_array


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
    '''
    Parameters:
        n: an integer, the window size of the moving average. If 1, this function returns the original data.
    '''
    if axis != 0:
        data = np.swapaxes(data, axis, 0)
    ret = np.cumsum(data, axis=0)
    ret[n:] = ret[n:] - ret[:-n]
    ret = ret[n-1:] / n
    if axis != 0:
        ret = np.swapaxes(ret, axis, 0)
    return ret

def data_synchronization(data_list, t_list, minimum_gap_ratio=1.5, minimum_gap_abs=None, resamplingT=None, return_data_list=True):
    '''
    Purpose:
        pass
    Parameters:
        return_data_list: if True, return the data_list with time synchronized. If false, return a data array of size [t_size, len(data_list), data_shape[1:]]
        resamplingT: if not given, the function will try to find the best resamplingT. 
    '''
    if resamplingT is None:
        resamplingT = define_synchronized_t(t_list=t_list, minimum_gap_ratio=minimum_gap_ratio, minimum_gap_abs=minimum_gap_abs)
    resampled_data_list = []
    for data_ind, data in enumerate(data_list):
        t = t_list[data_ind]
        resampled_data_list.append(dataFillAndLowPass(t, data, resamplingT=resamplingT))
    if return_data_list:
        return resamplingT, resampled_data_list
    else:
        data_array = np.zeros((len(resamplingT), len(data_list), *data_list[0].shape[1:]))
        for data_ind, data in enumerate(data_list):
            data_array[:, data_ind] = resampled_data_list[data_ind]
        return resamplingT, data_array

def define_synchronized_t(t_list, minimum_gap_ratio=1.5, minimum_gap_abs=None):
    stdInd = np.argmin(np.array([len(t) for t in t_list]))
    if not minimum_gap_abs:
        stdT = t_list[stdInd]
        minimum_gap_abs = np.min(np.diff(stdT))*minimum_gap_ratio
    for t_list_ind, t in enumerate(t_list):
        if t_list_ind == stdInd:
            continue
        pos = np.searchsorted(t, stdT)
        pos[pos == len(t)] = 0
        left = np.abs(stdT - t[pos-1])
        right = np.abs(t[pos] - stdT)
        diff_ = np.concatenate((left[:, None], right[:, None]), axis=-1)
        stdT = stdT[np.min(diff_, axis=-1) < minimum_gap_abs]
    return stdT

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

def electric_current_from_magetic_measurement(x, B, grad_method='Zhou', unit_str='mA/km^2'):
    '''
    Parameters:
        x: use astropy quantity, x and B should be are synchronized and be of [numberOfEvents, numberOfMeasurementPoints, 3]
        B: use astropy quantity
    Examples:
        data_list = []
        t_list = []
        std_source = 'FGM'
        dataName = 'B'
        x_source = 'MEC89Q'
        unit_str = 'mA/km^2'
        for spacecraft in spacecrafts:
            t_list.extend([spacecraft.data[std_source]['t'], spacecraft.data[x_source]['t']])
            data_list.extend([spacecraft.data[std_source][dataName], spacecraft.data[x_source]['x']])
        resamplingT = t_list[0]
        resamplingT, resampled_data_list = dat.data_synchronization(data_list, t_list, resamplingT=resamplingT)
        B = np.concatenate([d[:, None, ...] for d in resampled_data_list[::2]], axis=1)[..., :3] * u.nT
        x = np.concatenate([d[:, None, ...] for d in resampled_data_list[1::2]], axis=1) * u.km
        data = j_mag = dat.electric_current_from_magetic_measurement(x, B, grad_method='Zhou', unit_str='mA/km^2')
    '''
    #field_dict_list = [{'t': spacecraft.data[std_source]['t'], 'f': spacecraft.data[std_source][dataName]} for spacecraft in spacecrafts]
    #x_dict_list = [{'t': spacecraft.data[x_source]['t'], 'x': spacecraft.data[x_source][dataName]} for spacecraft in spacecrafts]
    #for scInd in range(len(spacecrafts)):
    #    field_dict = field_dict_list[scInd]
    #    x_dict = x_dict_list[scInd]
    #    x = dat.linear_interpolate_array(field_dict['t'], x_dict['t'], x_dict['x'])
    #    data_list.append(np.concatenate((field_dict['f'], x), axis=-1))
    #    t_list.append(field_dict['t'])
    return curl(B.to('nT').value, x.to('km').value, d=1) * (u.nT/u.km / const.mu0).to(unit_str)


def electric_current_from_plasma_measurement(electron, proton, unit_str='mA/km^2'):
    '''
    Purpose:
        to calculate electric current from \sum_{s=p,e} n_s q_s \vec{v}_s
    Parameters:
        electron: a dict containing three keys 't', 'n', and 'v'
        proton: a dict containing three keys 't', 'n', and 'v'
        unit_str: the units of returned data 
    '''
    data = np.concatenate((electron['v'], electron['n'][:, None]), axis=-1)
    data_t_p = linear_interpolate_array(proton['t'], electron['t'], data)
    v_e = data_t_p[:, :3]
    n_e = data_t_p[:, 3]
    return (proton['n'][:, None] * proton['v'] - n_e[:, None] * v_e) * (u.cm**(-3) * u.km/u.s * const.e.si).to(unit_str)

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


def datetime2epoch(dateTime, epochType='CDF_EPOCH'):
    if epochType == 'CDF_EPOCH':
        return cdflib.cdfepoch.compute_epoch(ot.datetime2list(dateTime))
    elif epochType == 'CDF_TIME_TT2000':
        return cdflib.cdfepoch.compute_tt2000(ot.datetime2list(dateTime))
    else:
        raise Exception("datetime object resolution to microsecond")

def epoch2datetime(epoch, epochType='CDF_EPOCH'):
    if epochType == 'CDF_EPOCH':
        return dt.datetime(*cdflib.cdfepoch.breakdown_epoch(epoch)[:6])
    elif epochType == 'CDF_TIME_TT2000':
        return dt.datetime(*cdflib.cdfepoch.breakdown_tt2000(epoch)[:6])

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

    def __init__(self, datetime=None, CDF_EPOCH=None, CDF_TIME_TT2000=None, astropy_Time=None, standard_fm='astropy_Time'):
        '''
        Parameters:
            datetime: a list nested to any degree of depth whose final element is datetime object. Input either this parameter or epochs to initialize the instance. This format refers to datetime.datetime
            CDF_EPOCH: a list nested to any degree of depth or a ndarray.
            CDF_TIME_TT2000: a list nested to any degree of depth or a ndarray.
            standard_fm: standard internal format for storage. 'astropy_Time ' is the default. 'CDF_TIME_TT2000_components' is deprecated.
        '''
        self.standard_fm = standard_fm
        self.data = {}
        if isinstance(datetime, dt.datetime):
            datetime = [datetime]
        self.data['datetime'] = datetime
        if CDF_EPOCH is not None:
            CDF_EPOCH = np.array(CDF_EPOCH)
            if len(CDF_EPOCH.shape) == 0:
                CDF_EPOCH = CDF_EPOCH[..., None]
        self.data['CDF_EPOCH'] = CDF_EPOCH
        if CDF_TIME_TT2000 is not None:
            CDF_TIME_TT2000 = np.array(CDF_TIME_TT2000)
            if len(CDF_TIME_TT2000.shape) == 0:
                CDF_TIME_TT2000 = CDF_TIME_TT2000[..., None]
        self.data['CDF_TIME_TT2000'] = CDF_TIME_TT2000
        self.data['astropy_Time'] = astropy_Time
        for data_format, data in self.data.items():
            if data is not None:
                if self.standard_fm == 'CDF_TIME_TT2000_components': # this option is deprecated
                    components = Epochs.breakdown(data=data, fm=data_format)
                    self.standard_data = components
                elif self.standard_fm == 'astropy_Time':
                    if self.standard_fm == data_format:
                        self.standard_data = self.data['astropy_Time']
                    else:

                        self.standard_data = Epochs.convert_to_standard_data(data, fm=data_format)
                break

#        datetime2epochWithEpochType = functools.partial(datetime2epoch, epochType=epochType)
#        epoch2datetimeWithEpochType = functools.partial(epoch2datetime, epochType=epochType)
#        self.epochType = epochType
#        if dateTimeList is not None:
#            self.dateTimeList = dateTimeList
#            self.epochs = np.array(map_multi_dimensional_list(datetime2epochWithEpochType, self.dateTimeList))
#        elif epochs is not None:
#            if isinstance(epochs, list):
#                epochs = np.array(epochs)
#            self.epochs = epochs
#            self.dateTimeList = map_multi_dimensional_list(epoch2datetimeWithEpochType, epochs.tolist())
    @classmethod
    def convert_to_standard_data(cls, data, fm):
        if fm == 'CDF_EPOCH':
            encoded_data = cdflib.epochs.CDFepoch.encode_epoch(data)
            return Time(encoded_data)
        elif fm == 'CDF_TIME_TT2000':
            return Time(2000, format='jyear') + TimeDelta(data/10**9*u.s)
        elif fm == 'datetime':
            return Time(data)

    @classmethod
    def breakdown(cls, data, fm):
        if fm == 'CDF_EPOCH':
            components = cdflib.epochs.CDFepoch.breakdown_epoch(data)
        elif fm == 'CDF_TIME_TT2000':
            components = cdflib.epochs.CDFepoch.breakdown_tt2000(data)
        elif fm == 'datetime':
            components = np.array(map_multi_dimensional_list(functools.partial(ot.datetime2list, epochType='CDF_TIME_TT2000'), data)).squeeze()
#            components = ot.datetime2list(data, epochType=epochType)
        return components

    @classmethod
    def compute(cls, components, fm): # from components to other format
        if fm == 'datetime':
            if len(components.shape) == 1:
                data = [dt.datetime(*components[:6])]
            elif len(components.shape) == 2:
                data = [dt.datetime(*components_[:6]) for components_ in components]
            else: raise Exception('wrong input of components')
        else:
            if fm == 'CDF_EPOCH':
                data = np.array(cdflib.epochs.CDFepoch.compute_epoch(components))
            elif fm == 'CDF_TIME_TT2000':
                data = cdflib.epochs.CDFepoch.compute_tt2000(components)
            if len(data.shape) == 0:
                data = data[..., None]
        return data

    @staticmethod
    def convert_from_astropy_Time(astropy_Time, fm):
        if fm == 'datetime':
            data = astropy_Time.utc.to_datetime(timezone=dt.timezone.utc)
        else:
            if fm == 'CDF_EPOCH':
                data = Epochs.compute(Epochs.breakdown(astropy_Time.utc.datetime, fm='datetime'), fm='CDF_EPOCH')
            elif fm == 'CDF_TIME_TT2000':
                data = (astropy_Time - Time(2000, format='jyear')).sec * 10**9
            if len(data.shape) == 0:
                data = data[..., None]
        return data

    def get_data(self, fm=None):
        if fm:
            if self.data[fm] is not None:
                pass
            else:
                if self.standard_fm == 'astropy_Time':
                    self.data[fm] = Epochs.convert_from_astropy_Time(self.standard_data, fm)
                elif self.standard_fm == 'CDF_TIME_TT2000_components':
                    self.data[fm] = self.compute(self.standard_data, fm)
        else:
            fm = self.standard_fm
        return self.data[fm]

    def epochRecords(self, ts, fm='CDF_TIME_TT2000', tolerance=1):
        '''
        Purposes:
            ts is a epoch series, this method finds the indices of the epoch series where the values are closest to the epochs in the Epochs object.
        Parameters:
            ts: a ndarray [..., timeIndex]
            tolerance (in second): the maximal allowed difference between t[record] and self.epoch
        '''
        epochs = self.get_data(fm=fm)
        records = np.argmin(np.abs(ts - epochs[..., None]), axis=-1)
        if fm == 'CDF_EPOCH':
            unit_conversion = 1000
        elif fm == 'CDF_TIME_TT2000':
            unit_conversion = 10**9
        tolerance_mask = np.abs(ts[..., records] - epochs) > tolerance*unit_conversion
        if np.any(tolerance_mask):
            records = records.astype(np.float64)
            records[tolerance_mask] = np.nan # nan representing out of tolerance range
        return records

    def findepochrange(self, fm='datetime', start_time=None, end_time=None):
        for fm_local, data in self.data.items():
            if data:
                ts = []
                if start_time:
                    ts.append(start_time)
                if end_time:
                    ts.append(end_time)
                para = {fm: ts}
                epochs = Epochs(**para)
                ts_local = epochs.get_data(fm_local)
                if end_time:
                    end_time = ts_local.pop()
                if start_time:
                    start_time = ts_local.pop()
                return cdflib.epochs.CDFepoch.findepochrange(data, start_time, end_time)


def map_multi_dimensional_list(func, l):
    if hasattr(l, '__iter__') and len(l) > 0:
        if type(l[0]) != list:
            return [func(v) for v in l]
        else:
            return [map_multi_dimensional_list(func, v) for v in l]
    else:
        return []


def dataFillAndLowPass(t, data, axis=0, resamplingT=None, tDistribution='evenlySpaced', tStepPecentageCriterion=0.9, lowpassCutoff=None, gapThreshold=None, minNumberOfPoints=2, returnShiftQ=False, badpointsMask=None):
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
            tStepPecentageCriterion: In the general case of irregually epoch records, there are multiple possible time steps between consecutive records. The portion of the most frequent temporal step in all occurance of all time steps should be larger than this parameter, otherwise the program raise exception.
            lowpassCutoff: lowpass cutoff frequency.
            gapThreshold: can be a number, in which case the t is divided into blocks of t such that the gap between consecutive blocks is larger than gapThreshold. or a string in the form of 'num*', such as '3*', which set the gapThreshold to be three times the minimum gap in the t.
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
            assert np.all(tDiff >= 0)
            if type(gapThreshold) is str and gapThreshold[-1] == '*':
                gapThreshold = float(gapThreshold[:-1])*np.min(tDiff)
            if gapThreshold:
                section = 1 + np.argwhere(tDiff > gapThreshold).squeeze()
                if len(section.shape) == 0:
                    section = np.array([section])
                tBlocks = np.split(t, section) # each element in tBlocks is a consecutive time series in which the maximal gap ls less than gapThreshold
                dataBlocks = np.split(data, section, axis=axis)
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
                        tHomogeneous = np.arange(tBlock[0], tBlock[-1], tStep)
#                        cs = interpolate.CubicSpline(tBlock, dataBlock, axis=axis)
#                        dataBlockHomogeneous = cs(tHomogeneous)
                        dataBlockHomogeneous = linear_interpolate_array(tHomogeneous, tBlock, dataBlock, axis=axis)
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
#        cs = interpolate.CubicSpline(t, data, axis=axis)
#        resampledData = cs(resamplingT)
        resampledData = linear_interpolate_array(resamplingT, t, data, axis=axis)
        return resampledData


def linear_interpolate_array(x_new, x_ori, array, axis=0):
    '''
    Parameters:
        x_new: the new x coordinates after interpolation
        x_ori: the x coordinate of array along axis=axis
        axis: the axis of the array for coordinate x
    '''
    if len(array.shape) == 1:
        return np.interp(x_new, x_ori, array)
    assert len(x_ori) == array.shape[axis]
    if axis != 0:
        array = array.swapaxes(axis, 0)
    shape_swaped = array.shape
    numberOfFunctions = np.prod(array.shape[1:])
    array = array.reshape((len(x_ori), numberOfFunctions))
    y = np.zeros((len(np.atleast_1d(x_new)), numberOfFunctions))
    for ind in range(array.shape[1]):
        y[:, ind] = np.interp(x_new, x_ori, array[:, ind])
    y = y.reshape((len(np.atleast_1d(x_new)), *shape_swaped[1:]))
    if axis != 0:
        y = y.swapaxes(axis, 0)
    return y


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
    alfvenSpeed = BStrength / np.sqrt(const.mu0 * const.m_p * n)
    return alfvenSpeed

def sonicSpeed(T, gamma=5/3):
    return (gamma * const.k_B * T / const.m_p)**(1/2)

def magnetosonicSpeed(BStrength, n, T):
    return np.sqrt(alfvenSpeed(BStrength, n)**2 + sonicSpeed(T)**2)

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
            elif np.all(np.array(inRange) == np.array([0, 360])):
                dataTransformed = data/180*np.pi
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
        to interpolate the unit quaternion for coordinate transformation of vectors. see the Calibration and Measurement Algorithms Document for MMS, Section 3, G-7
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

quaternionFromMatrix = sppcoo.quaternionFromMatrix


def angleBetweenVectors(v1, v2, ignore_direction=False):
    '''
    Purpose:
        to calculate the angle in degrees between two vectors, v1 and v2
    Parameters:
        v1 and v2 are of same shape, ndarray(..., 3)
        ignore_direction: v1 is deemed equivalent to -v1. Thus the angle is always less than 90 degrees.
    '''
    v1Unit = v1 / np.linalg.norm(v1, axis=-1)[..., None]
    v2Unit = v2 / np.linalg.norm(v2, axis=-1)[..., None]
    angle = np.arccos(np.sum(v1Unit * v2Unit, axis=-1))/np.pi*180
    if ignore_direction:
        angle = np.abs(180 * (np.sign(angle - 90 - 10**(-31)) + 1)/2 - angle)
    return angle


def binomial(a, b):
    return np.math.factorial(a)//np.math.factorial(b)//np.math.factorial(a-b)


def leastSquarePolynomialApproximation(j, x, d, omega=None, regularizationMethod=None, regPara=None, solver='direct', numberOfIterations=1000):
    '''
    This Function do something
    Parameters:
        j: the sampled data of dimension [..., M, s] where M is the number of samples, s is the number of functions
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

def linearGradient(j, x, d=1, omega=None, regularizationMethod=None, regPara=None, solver='direct', numberOfIterations=1000):
    xMean = np.mean(x, axis=-2)
    x = x - xMean[..., None, :]
    coeff = leastSquarePolynomialApproximation(j, x, d=d, omega=omega)
    return coeff[..., 1:4, :]

def curl(j, x, d=1, omega=None):
    gradients = linearGradient(j, x, d, omega=omega)
    leviCivitaTensor = makeLeviCivitaTensor()
    return np.sum(np.sum(gradients[..., None] * leviCivitaTensor[None, ...], axis=-3), axis=-2)

def div(j, x, d=1, omega=None):
    return np.trace(linearGradient(j, x, d, omega=omega), axis1=-2, axis2=-1)

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

def rotateVectorUsingQuaternion(vec, quat, scalarPos='last'):
    if isinstance(vec, u.Quantity):
        vec_unit = vec.unit
    if scalarPos == 'last':
        quat_inverse = np.copy(quat)
        quat_inverse[..., :3] *= -1
        vec_quat = np.zeros((*vec.shape[:-1], 4))
        vec_quat[..., :3] = vec
        new_vec = quaternionMultiply([quat, vec_quat, quat_inverse])[..., :3]
    if isinstance(vec, u.Quantity):
        new_vec *= vec_unit
    return new_vec

def mask_dict_of_ndarray(dic, mask, copy=False):
    dic_new = {}
    for key, item in dic.items():
        if copy:
            dic_new[key] = np.copy(item[mask])
        else:
            dic_new[key] = item[mask]
    return dic_new

def concatenate_dict_of_ndarray(dic_list, axis=None):
    dic_new = {}
    for key in dic_list[0].keys():
        if isinstance(axis, int):
            dic_new[key] = np.concatenate([dic[key] for dic in dic_list], axis=axis)
        elif axis is None:
            dic_new[key] = np.concatenate([dic[key][None, ...] for dic in dic_list], axis=0)
    return dic_new

def arg_split(data, gap_threshold=1):
    '''
    data: 1d ndarray
    gap_threshold: can be a number, in which case the data is divided into blocks such that the gap between consecutive blocks is larger than gap_threshold. or a string in the form of 'num*', such as '3*', which set the gap_threshold to be three times the minimum gap in the data.
    return:
        break_points: the split data is retrieved by
            data_split = []
            for ind in range(len(break_points)-1):
                s_ = slice(*break_points[ind:ind+2])
                data_split.append(data[s_])
    '''
    tdiff = np.diff(data)
    if type(gap_threshold) is str and gap_threshold[-1] == '*':
        gap_threshold = float(gap_threshold[:-1])*np.min(tdiff)
    break_points = np.nonzero(tdiff > gap_threshold)[0] + 1
    break_points = np.insert(break_points, 0, 0)
    break_points = np.append(break_points, len(data))
    return break_points


def data_split(t, data, gap_threshold=1):
    '''
    data: 2d ndarray with data[:, 0] be the time data acoording to which the data will be splited.
    gap_threshold: can be a number, in which case the data is divided into blocks such that the gap between consecutive blocks is larger than gap_threshold. or a string in the form of 'num*', such as '3*', which set the gap_threshold to be three times the minimum gap in the data.
    return:
    '''
    break_points = arg_split(t, gap_threshold=gap_threshold)
    t_split = []
    data_split = []
    for ind in range(len(break_points)-1):
        s_ = slice(*break_points[ind:ind+2])
        t_split.append(t[s_])
        data_split.append(data[s_])
    return t_split, data_split


def spectrogram(t, data, nperseg, noverlap, gap_threshold='6*', window='hamming', verbose=True, minNumberOfPoints=2):
    '''
    A wrapper of the scipy.signal.spectrogram that split the data into chuncks containing evenly spaced data points and perform scipy.signal.spectrogram on these chuncks individually.
    '''
    if t.size == 0:
        logging.warning('spectrogram input size 0')
        return [[]]*3
    t_, data_ = dataFillAndLowPass(t, data, axis=0, resamplingT=None, tDistribution='evenlySpaced', tStepPecentageCriterion=0.2, lowpassCutoff=None, gapThreshold=gap_threshold, minNumberOfPoints=minNumberOfPoints, returnShiftQ=False, badpointsMask=None)
    sample_spacing = (t_[1] - t_[0])
    t_split, data_split_ = data_split(t_, data_, gap_threshold=gap_threshold)
    t_spectrogram_split = []
    data_spectrogram_split = []
    for t_, data_ in zip(t_split, data_split_):
        if len(t_) >= nperseg:
            freq_spectrogram, t_spectrogram, sxx = signal.spectrogram(data_, fs=1/sample_spacing, nperseg=nperseg, noverlap=noverlap, axis=0, window=window, detrend='linear')
            sxx = np.moveaxis(sxx, -1, 0)
            t_spectrogram_split.append(t_spectrogram+t_[0])
            data_spectrogram_split.append(sxx)
    if verbose:
        print('number of split: {}'.format(len(t_spectrogram_split)))
        print('counts in each split: ', [len(t_) for t_ in t_spectrogram_split])
    return freq_spectrogram, t_spectrogram_split, data_spectrogram_split


def CalculateDynamicPressure(n, v, spacecraft='mms'):
    if spacecraft in ['mms', 'maven']:
        speed = np.linalg.norm(v, axis=-1)
        dynamicPressure = 1.67 * 10**(-6) * n * speed**2 # nPa
        return dynamicPressure

def datetime_floor(t, round_gap=None):
    '''
    Purpose: find the floor of a time in a day
    Parameters:
        round_gap: timedelta object
    '''
    day_start = dt.datetime(t.year, t.month, t.day, tzinfo=t.tzinfo)
    return day_start + dt.timedelta(seconds=(t - day_start).total_seconds()//round_gap.total_seconds() * round_gap.total_seconds())

def datetime_ceil(t, round_gap=None):
    '''
    Purpose: find the ceil of a time in a day
    Parameters:
        round_gap: timedelta object
    '''
    return datetime_floor(t, round_gap=round_gap) + round_gap


def angle_factor_in_integrating_1st_moment(theta_bound, phi_bound, axis=None, check_phi=True):
    r'''
    Purpose:
        this function is a part of a group functions to calculate moment from distribution function.
    Parameters:
        theta_bound: pass
        phi_bound: make sure np.diff(phi_bound) gives the expected difference between the upper and lower bound. For example, if np.array([3/2*np.pi, 7/4*np.pi, 1/4*np.pi, 1/2*np.pi]) is provided, np.diff may result in a bad diff between 7/4*np.pi and 1/4*np.pi. The intended difference may be 1/2*np.pi. Therefore, make sure convert this array to np.array([3/2*np.pi, 7/4*np.pi, 9/4*np.pi, 5/2*np.pi]) before using this function
        axis: 0, 1, or 2, to choose the angle factor used in calculating vx, vy, or vz.
    Return:
        angle_factor \gamma: defined according to v_i = \sum\sum \parens(\Delta E/E) \diff{J_E}{E} \gamma_i, where \diff{J_E}{E}=E\vec{v}f d^3v/dE is differential energy flux.
    '''
    if check_phi:
        if np.any(np.abs(np.diff(phi_bound)) > np.pi):
            raise Exception('check phi_bound to make sure np.diff(phi_bound) gives the expected difference between the upper and lower bound.')
    if axis == 0: # x-axis
        def theta_func(theta):
            return theta/2 - 1/4*np.sin(2*theta)
        def phi_func(phi):
            return np.sin(phi)
    elif axis == 1:
        def theta_func(theta):
            return theta/2 - 1/4*np.sin(2*theta)
        def phi_func(phi):
            return -np.cos(phi)
    elif axis == 2:
        def theta_func(theta):
            return -1/2*np.cos(theta)**2
        def phi_func(phi):
            return phi
    return (theta_func(theta_bound[1:]) - theta_func(theta_bound[:-1]))[..., None] * (phi_func(phi_bound[1:]) - phi_func(phi_bound[:-1]))[None, :]


def angle_factors_in_integrating_1st_moment(theta_bound, phi_bound, axis=['x', 'y', 'z'], check_phi=True):
    r'''
    Purpose:
        this function is a part of a group functions to calculate moment from distribution function.
    Parameters:
        theta_bound: pass
        phi_bound: make sure np.diff(phi_bound) gives the expected difference between the upper and lower bound. For example, if np.array([3/2*np.pi, 7/4*np.pi, 1/4*np.pi, 1/2*np.pi]) is provided, np.diff may result in a bad diff between 7/4*np.pi and 1/4*np.pi. The intended difference may be 1/2*np.pi. Therefore, make sure convert this array to np.array([3/2*np.pi, 7/4*np.pi, 9/4*np.pi, 5/2*np.pi]) before using this function
        axis: a list such as ['x', 'y', 'z'] or a str 'x' to choose the angle factor used in calculating vx, vy, or vz.
    Return:
        angle_factor \gamma: defined according to v_i = \sum\sum \parens(\Delta E/E) \diff{J_E}{E} \gamma_i, where \diff{J_E}{E}=E\vec{v}f d^3v/dE is differential energy flux.
    Examples:
        # to calculate bulk_velocity from maven swia coarsesvy3d data
        angle_factors = np.swapaxes(dat.angle_factors_in_integrating_1st_moment(theta_bound, phi_bound), axis1=0, axis2=1)[None, ...]
        factor_v = ((energy_width/energy)[None, :, None, None] * def_ion)[..., None] * angle_factors * u.cm**(-3) *u.cm/u.s
        bulk_velocity = np.sum(factor_v, axis=(1, 2, 3))/density[:, None]
    '''
    if check_phi:
        if np.any(np.abs(np.diff(phi_bound)) > np.pi):
            raise Exception('check phi_bound to make sure np.diff(phi_bound) gives the expected difference between the upper and lower bound.')

    if isinstance(axis, list):
        angle_factors = []
        for axis_ in axis:
            angle_factors.append(angle_factors_in_integrating_1st_moment(theta_bound, phi_bound, axis=axis_, check_phi=check_phi)[..., None])
        return np.concatenate(angle_factors, axis=-1)
    else: # main branch
        if axis == 'x': # x-axis
            def theta_func(theta):
                return theta/2 - 1/4*np.sin(2*theta)
            def phi_func(phi):
                return np.sin(phi)
        elif axis == 'y':
            def theta_func(theta):
                return theta/2 - 1/4*np.sin(2*theta)
            def phi_func(phi):
                return -np.cos(phi)
        elif axis == 'z':
            def theta_func(theta):
                return -1/2*np.cos(theta)**2
            def phi_func(phi):
                return phi
        theta_diff = np.diff(theta_bound, axis=0)
        phi_diff = np.diff(phi_bound, axis=0)
        sign = 1
        if np.all(theta_diff < 0): sign *= -1
        elif np.all(theta_diff > 0): pass
        else: raise Exception('Unexpected order of theta_bound')
        if np.all(phi_diff < 0): sign *= -1
        elif np.all(phi_diff > 0): pass
        else: raise Exception('Unexpected order of phi_bound')
        return sign * (theta_func(theta_bound[1:]) - theta_func(theta_bound[:-1]))[..., None] * (phi_func(phi_bound[1:]) - phi_func(phi_bound[:-1]))[None, :]

def angle_factors_in_integrating_2nd_moment(theta_bound, phi_bound, axis=['xx', 'xy', 'xz', 'yy', 'yz', 'zz'], check_phi=True):
    r'''
    Purpose:
        This function is a part of a group functions to calculate moment from distribution function. Check angle_factors_in_integrating_1st_moment for more information.
    Parameters:
        theta_bound:
        phi_bound: make sure np.diff(phi_bound) gives the expected difference between the upper and lower bound. For example, if np.array([3/2*np.pi, 7/4*np.pi, 1/4*np.pi, 1/2*np.pi]) is provided, np.diff may result in a bad diff between 7/4*np.pi and 1/4*np.pi. The intended difference may be 1/2*np.pi. Therefore, make sure convert this array to np.array([3/2*np.pi, 7/4*np.pi, 9/4*np.pi, 5/2*np.pi]) before using this function
    Return:
        angle_factor \gamma: defined according to v_iv_j = \sum\sum v\parens(\Delta E/E) \diff{J_E}{E} \gamma_{ij}, where \diff{J_E}{E}=E\vec{v}f d^3v/dE is differential energy flux.
    '''
    if check_phi:
        if np.any(np.abs(np.diff(phi_bound)) > np.pi):
            raise Exception('check phi_bound to make sure np.diff(phi_bound) gives the expected difference between the upper and lower bound.')

    if isinstance(axis, list):
        angle_factors = []
        for axis_ in axis:
            angle_factors.append(angle_factors_in_integrating_2nd_moment(theta_bound, phi_bound, axis=axis_, check_phi=check_phi)[..., None])
        return np.concatenate(angle_factors, axis=-1)
    else: # main branch
        if axis == 'xx':
            def theta_func(theta):
                return -(3*np.cos(theta))/4 + 1/12*np.cos(3*theta)
            def phi_func(phi):
                return phi/2 + 1/4*np.sin(2*phi)
        elif axis == 'xy':
            def theta_func(theta):
                return -(3*np.cos(theta))/4 + 1/12*np.cos(3*theta)
            def phi_func(phi):
                return -(1/2)*np.cos(phi)**2
        elif axis == 'xz':
            def theta_func(theta):
                return np.sin(theta)**3/3
            def phi_func(phi):
                return phi
        elif axis == 'yy':
            def theta_func(theta):
                return -(3*np.cos(theta))/4 + 1/12*np.cos(3*theta)
            def phi_func(phi):
                return phi/2 - 1/4*np.sin(2*phi)
        elif axis == 'yz':
            def theta_func(theta):
                return np.sin(theta)**3/3
            def phi_func(phi):
                return -np.cos(phi)
        elif axis == 'zz':
            def theta_func(theta):
                return -(1/3)*np.cos(theta)**3
            def phi_func(phi):
                return phi
        else: raise Exception('Unexpected axis {}'.format(axis))
        theta_diff = np.diff(theta_bound, axis=0)
        phi_diff = np.diff(phi_bound, axis=0)
        sign = 1
        if np.all(theta_diff < 0): sign *= -1
        elif np.all(theta_diff > 0): pass
        else: raise Exception('Unexpected order of theta_bound')
        if np.all(phi_diff < 0): sign *= -1
        elif np.all(phi_diff > 0): pass
        else: raise Exception('Unexpected order of phi_bound')
        return sign * (theta_func(theta_bound[1:]) - theta_func(theta_bound[:-1]))[..., None] * (phi_func(phi_bound[1:]) - phi_func(phi_bound[:-1]))[None, :]


def pearson_correlation_coefficient(x, y):
    '''
    Parameters:
        x: 1d array or 2d array with the first axis representing observations and the second representing different variables.
        y: same shape as x
    '''
    x_ = x - np.mean(x, axis=0)
    y_ = y - np.mean(y, axis=0)
    return np.mean(x_ * y_, axis=0) / (np.linalg.norm(x_, axis=0) * np.linalg.norm(y_, axis=0) / len(x_))

def shifted_pearson_correlation_coefficient(x, y):
    '''
    Purpose:
        this function is for finding the timing lag between data obtained by two spacecraft
    Parameters:
        x: the shorter one. [observation_number, variable_number]
        y: the longer one. [more_observation_number, variable_number]
    Return:
        rho: the correlation coefficient of length len(y) - len(x) + 1. rho[i] = pearson_correlation_coefficient(x, y[i:len(x)+i])
    '''
    obsN = len(x)
    obsN_y = len(y)
    rhoN = obsN_y-obsN+1
    x_ = np.repeat(x[..., None], rhoN, axis=-1).reshape((obsN, -1))
    y_ = np.swapaxes(np.lib.stride_tricks.sliding_window_view(y, obsN, axis=0), axis1=0, axis2=-1).reshape((obsN, -1))
    rho = pearson_correlation_coefficient(x_, y_).reshape((-1, rhoN))
    return rho
