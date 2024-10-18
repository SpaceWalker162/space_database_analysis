import numpy as np
from astropy.constants import *

solarNorthPoleOfRotationInICRFEquatorial = np.array([63.87, 286.13])/180*np.pi # latitude and longitude
saturnianNorthPoleOfRotationInICRFEquatorial = np.array([83.54, 40.58])/180*np.pi

radius_Mercury = Constant('R_Mercury', 'radius of Mercury', 2.4397 * 10**6, unit='m', uncertainty=0.0, reference='wiki')
radius_Earth = Constant('R_E', 'radius of Earth', 6.378137 * 10**6, unit='m', uncertainty=0.0, reference='wiki')
radius_Mars = Constant('R_Mars', 'radius of Mars', 3.3895 * 10**6, unit='m', uncertainty=0.0, reference='wiki')
radius_Jupiter = Constant('R_Jupiter', 'radius of Jupiter', 6.9911 * 10**7, unit='m', uncertainty=0.0, reference='wiki')
radius_Saturn = Constant('R_Saturn', 'radius of Saturn', 5.8232 * 10**7, unit='m', uncertainty=0.0, reference='wiki')

magnetic_permeability_in_vacuum = mu0
