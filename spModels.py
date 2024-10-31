import numpy as np
#import constants as con
import cdflib    # see github.com/MAVENSDC/cdflib
#import tarfile
from datetime import datetime
from datetime import timedelta
import os
import sys
import logging
from wolframclient.language import wl
import wolframclient
from wolframclient.evaluation import WolframLanguageSession
import space_database_analysis.dataAnalysisTools as dat


class BowShockPositionModel:
    def __init__(self, wlPath=None, modelName='J05', adiabaticIndex=5/3, solarWindParams=False, B=None, BStrength=None, n=None, rho=None, v=None, vStrength=None, T=None, P=None, PDyn=None):
        '''
        Parameters:
            wlPath: the path to wolframkernel, for example, '/usr/local/Wolfram/Mathematica/12.0/Executables/WolframKernel'
            modelName: J05 by default, see doi:10.1016/j.pss.2004.09.032
            adiabaticIndex: usually taken as 5/3, but suggested to take 2.15 in order to find a better description of the bow shock position by Zhuang and Russell (1981, doi:10.1029/JA086iA04p02191)
            B: interplanetary magnetic field, a ndarray of three components, in nT
            n: number density in cm^{-3}
            rho: mass density
            v: velocity of solar wind, a ndarray of three components
            T: temperature
            P: pressure
        '''
        self.modelName = modelName
        self.adiabaticIndex = adiabaticIndex
        self.B = B
        if B is not None:
            BStrength = np.linalg.norm(B, axis=-1)
        self.BStrength = BStrength
        self.n = n
        self.rho = rho
        self.v = v
        if v is not None:
            vStrength = np.linalg.norm(v, axis=-1)
        self.vStrength = vStrength
        self.T = T
        self.P = P
        if wlPath is not None:
            self.wlSession = WolframLanguageSession(wlPath)
            self.initMathematicaModel()
        if solarWindParams:
            if self.modelName == 'J05':
                BStrength = self.BStrength
                gamma = self.adiabaticIndex
                N = self.n
                V = self.vStrength
                C = 91.55
                D = 0.937 * (0.846 + 0.042 * BStrength)
                alfvenSpeed = dat.alfvenSpeed(BStrength, N)
                MA = V / alfvenSpeed
                self.cont = C / (N*V**2)**(1/6) * (1 + D*((gamma-1)*MA**2+2)/((gamma+1)*(MA**2-1)))

    def modelR(self, theta, phi):
        if self.modelName == 'J05':
            R = self.cont * J05ReferenceR(theta, phi)
        return R 

    def initMathematicaModel(self):
        '''J05 model by default'''
        if self.modelName == 'J05':
            coeff = 'coeff = {a11 -> 0.45, a22 -> 1, a33 -> 0.8, a12 -> 0.18, a14 -> 46.6, a24 -> -2.2, a34 -> -0.6, a44 -> -618};'
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
        self.wlSession.evaluate(initModelCMD)
        self.wlSession.evaluate(coeff)

    def modelProjection(self, xGSE):
        xGSE_ = 'xGSE = {{x -> {}, y -> {}, z -> {}}};'.format(*xGSE)
        self.wlSession.evaluate(xGSE_)
        self.wlSession.evaluate('contValue = contraction/.coeff/.xGSE')
        numberOfContraction = 0
        for i in range(1, 3):
            if (self.wlSession.evaluate('cont/.contValue[[{}]]'.format(i))) > 0:
                numberOfContraction +=1
                if numberOfContraction == 1:
                    contractionInd = i
                elif numberOfContraction == 2:
                    raise Exception('two contractions > 0')
        cont = self.wlSession.evaluate('cont/.contValue[[{}]]'.format(contractionInd))
        return cont*xGSE


    def modelCurvaturesAndDirections(self, xGSE):
        xGSE_ = 'xGSE = {{x -> {}, y -> {}, z -> {}}};'.format(*xGSE)
        self.wlSession.evaluate(xGSE_)
        self.wlSession.evaluate('contValue = contraction/.coeff/.xGSE')
        numberOfContraction = 0
        for i in range(1, 3):
            if (self.wlSession.evaluate('cont/.contValue[[{}]]'.format(i))) > 0:
                numberOfContraction +=1
                if numberOfContraction == 1:
                    contractionInd = i
                elif numberOfContraction == 2:
                    raise Exception('two contractions > 0')
        self.wlSession.evaluate('zsInXYValue = zsInXY/.coeff/.xGSE[[{1,2}]]'+'/.contValue[[{}]]'.format(contractionInd))
        for i in range(1, 3):
            if np.abs(self.wlSession.evaluate('z/.zsInXYValue[[{}]]'.format(i)) - xGSE[2]) < 10**(-4):
                zInd = i
        curvatures_ = np.array(self.wlSession.evaluate('curvaturesWithDirections[[{}]] /. coeff /. contValue[[{}]] /. xGSE'.format(zInd, contractionInd)))
        curvatures = np.append(curvatures_, 0)
        directions = np.zeros((3, 3))
        directions[:, :2] = np.array(self.wlSession.evaluate('Normalize /@ (directionsSym[[{}]] /. coeff /. contValue[[{}]]) /. xGSE'.format(zInd, contractionInd))).T
        directions[:, 2] = np.array(self.wlSession.evaluate('n[[{}]] /. coeff /. xGSE /. contValue[[{}]]'.format(zInd, contractionInd)))
        return curvatures, directions

    def terminate(self):
        if self.wlSession is not None:
            self.wlSession.terminate()

def J05ReferenceR(theta, phi, a11 = 0.45, a22 = 1, a33 = 0.8, a12 = 0.18, a14 = 46.6, a24 = -2.2, a34 = -0.6, a44 = -618, R0 = 11.8954):
    def mysqrt(x): return np.sqrt((1.+0j)*x)
    aux0=(-4.*(a24*(np.cos((phi+theta)))))+((-4.*(a14*(np.sin((phi-theta))\
    )))+(4.*(a14*(np.sin((phi+theta))))));
    aux1=(((4.*(a24*(np.cos((phi-theta)))))+((8.*(a34*(np.cos(theta))))+\
    aux0))**2);
    aux2=(4.*(a33*(np.cos((2.*theta)))))+((a22*(np.cos(((2.*phi)+(2.*\
    theta)))))+(2.*(a12*(np.sin((2.*phi))))));
    aux3=(-2.*(a11*(np.cos((2.*theta)))))+((-2.*(a22*(np.cos((2.*theta))))\
    )+aux2);
    aux4=(-2.*(a22*(np.cos((2.*phi)))))+((a22*(np.cos(((2.*phi)+(-2.*\
    theta)))))+aux3);
    aux5=((2.*a11)+((2.*a22)+((4.*a33)+((2.*(a11*(np.cos((2.*phi)))))+\
    aux4))))-(a12*(np.sin(((2.*phi)+(2.*theta)))));
    aux6=(aux5-(a12*(np.sin(((2.*phi)+(-2.*theta))))))-(a11*(np.cos(((2.*\
    phi)+(2.*theta)))));
    aux7=mysqrt((aux1+(-32.*(a44*(aux6-(a11*(np.cos(((2.*phi)+(-2.*theta))\
    ))))))));
    aux8=(4.*(a14*(np.sin((phi-theta)))))+((-4.*(a14*(np.sin((phi+theta)))\
    ))+aux7);
    aux9=(-4.*(a24*(np.cos((phi-theta)))))+((-8.*(a34*(np.cos(theta))))+((\
    4.*(a24*(np.cos((phi+theta)))))+aux8));
    aux10=(4.*(a33*(np.cos((2.*theta)))))+((a22*(np.cos(((2.*phi)+(2.*\
    theta)))))+(2.*(a12*(np.sin((2.*phi))))));
    aux11=(-2.*(a11*(np.cos((2.*theta)))))+((-2.*(a22*(np.cos((2.*theta)))\
    ))+aux10);
    aux12=(-2.*(a22*(np.cos((2.*phi)))))+((a22*(np.cos(((2.*phi)+(-2.*\
    theta)))))+aux11);
    aux13=((2.*a11)+((2.*a22)+((4.*a33)+((2.*(a11*(np.cos((2.*phi)))))+\
    aux12))))-(a12*(np.sin(((2.*phi)+(2.*theta)))));
    aux14=(aux13-(a12*(np.sin(((2.*phi)+(-2.*theta))))))-(a11*(np.cos(((2.\
    *phi)+(2.*theta)))));
    output=(0.5*aux9)/(aux14-(a11*(np.cos(((2.*phi)+(-2.*theta))))));
    # the foregoing code is converted from mathematica
    return output/R0


class bowShockAndMagnetopausePositionModels:
    def __init__(self, wlPath=None, modelName=None, **modelParas):
        self.modelName = modelName
        '''
        Parameters:
            wlPath: the path to wolframkernel, for example, '/usr/local/Wolfram/Mathematica/12.0/Executables/WolframKernel'
            modelName: Joy02BS, Joy02MP (doi:10.1029/2001JA009146)
        '''
        modelDict = {'Joy02BS': {'modelType': 'secondOrderSurface'},
                     'Joy02MP': {'modelType': 'secondOrderSurface'},
                     'GruesBeck18MarsBS': {'modelType': 'secondOrderSurface'},
                     'Went11': {'modelType': 'conicSection'},
                     'Kanani10': {'modelType': 'conicSection'},
                     'Shue98': {'modelType': 'conicSection'},
                     }
        self.modelDict = modelDict
        if wlPath is not None:
            self.wlSession = WolframLanguageSession(wlPath)
            self.modelParas = modelParas
#        if self.modelName in ['Joy02BS', 'Joy02MP', 'Went11', 'Kanani10']:
        self.modelType = self.modelDict[self.modelName]['modelType']
        self.initMathematicaModel()

    def initMathematicaModel(self):
        if self.modelName == 'Joy02BS':
            initJoy02Model(self.wlSession, model='BS', **self.modelParas)
        elif self.modelName == 'Joy02MP':
            initJoy02Model(self.wlSession, model='MP', **self.modelParas)
        elif self.modelName == 'Went11':
            initWent11Model(self.wlSession, **self.modelParas)
        elif self.modelName == 'Kanani10':
            initKanani10Model(self.wlSession, **self.modelParas)
        elif self.modelName == 'Shue98':
            initShue98Model(self.wlSession, **self.modelParas)
        elif self.modelName == 'GruesBeck18MarsBS':
            initGruesBeck18MarsBSModel(self.wlSession)

#    def calculatePosition(self, pos, toCal='r', returnR=True, returnXYZ=False):
    def calculatePosition(self, pos, toCal='r'):
        '''
        Parameters:
            pos: position, an array of [..., n]. The last dimension is represent the known parameters for a point.
            toCal: if 'r', then the last dimension of pos represents theta and phi. If 'z', then the last dimension of pos are x and y.
            posPara: a list of two elements. If ['theta', 'phi'], then the last dimension of pos are these components. Other possibility include but not limited to ['x', 'y'], ['r', 'theta']. Make sure the parameter follow the order r>theta>phi, x>y>z
        Return:
            posLastEles: the last element other than posPara.
        '''
        shape = pos.shape
        if len(shape) > 1:
            pos_ = pos.reshape((-1, shape[-1]))
        elif len(shape) == 1:
            pos_ = pos[None, :]
        elif len(shape) == 0:
            pos_ = np.array([[pos]])
        numberOfPoints = len(pos_)
        posLastEles = np.zeros(numberOfPoints)
        if self.modelType == 'secondOrderSurface':
            for ind in range(numberOfPoints):
                if toCal == 'r':
                    theta = pos_[ind, 0]
                    phi = pos_[ind, 1]
                    posLastEle_ = self.calculateRFromThetaPhi(theta=theta, phi=phi)
                    posLastEles[ind] = posLastEle_
                elif toCal == 'z':
                    x = pos_[ind, 0]
                    y = pos_[ind, 1]
                    wlCMD = '''zInXY/.{{x -> {x}, y -> {y}}}
                    '''.format(x=x, y=y)
                    posLastEle_ = np.array(self.wlSession.evaluate(wlCMD))
        elif self.modelType == 'conicSection':
            for ind in range(numberOfPoints):
                if toCal == 'r':
                    if shape[-1] == 1:
                        theta = pos_[ind, 0]
                        wlCMD = '''rInTheta/.{{theta -> {theta}}}
                        '''.format(theta=theta)
                        posLastEle_ = np.array(self.wlSession.evaluate(wlCMD))
                    elif shape[-1] == 2:
                        theta = pos_[ind, 0]
                        phi = pos_[ind, 1]
                        posLastEle_ = self.calculateRFromThetaPhi(theta=theta, phi=phi)
                    posLastEles[ind] = posLastEle_
                elif toCal == 'z':
                    x = pos_[ind, 0]
                    y = pos_[ind, 1]
                    wlCMD = '''z /. NSolve[(eqXYZ == 0)/.{{x->{x}, y->{y}}}, z, Reals]
                    '''.format(x=x, y=y)
                    posLastEle_ = np.array(self.wlSession.evaluate(wlCMD))
                    posLastEles[ind] = np.max(posLastEle_[np.isreal(posLastEle_)])
                else:
                    raise Exception('unknown toCal')
        posLastEles = posLastEles.reshape(shape[:-1])
        return posLastEles

    def calculateRFromThetaPhi(self, theta, phi):
        logging.info('theta: {}'.format(theta))
        logging.info('phi: {}'.format(phi))
        wlCMD = '''r /. NSolve[(eqRTPAtOriginXYZ0 == 0)/.{{\[Theta]->{theta}, \[Phi]->{phi}}}, r, Reals]
        '''.format(theta=theta, phi=phi)
        wlResult = self.wlSession.evaluate(wlCMD)
        if isinstance(wlResult, wolframclient.language.expression.WLSymbol):
            logging.warning('not found posLastEle_, but found:')
            logging.warning(wlResult)
            posLastEle_ = None
        else:
            posLastEle_ = np.array(wlResult).squeeze()
            logging.info('posLastEle_')
            logging.info(posLastEle_)
            posLastEle_ = np.where(posLastEle_ < 10**10, posLastEle_, np.zeros_like(posLastEle_))
            assert np.count_nonzero(posLastEle_ > 0) == 1 
            posLastEle_ = np.max(posLastEle_)
        return posLastEle_


def initJoy02Model(wlSession, model=None, dynamicPressure=None, origin=np.zeros(3)):
    '''
    Purpose: obtain the location of Jovian bow shock and magnetopause using the model doi:10.1029/2001JA009146.
    Parameters:
        wlSession: a wolframe session
        model: 'BS' or 'MP'
        dynamicPressure: upstream solar wind dynamic pressure
        origin: the origin of the coordinate system in the unit of R_J
    '''
    x0, y0, z0 = origin / 120
    if model == 'BS':
        a = -1.107 + 1.591*dynamicPressure**(-1/4)
        b = -0.566 - 0.812*dynamicPressure**(-1/4)
        c = 0.048 - 0.059*dynamicPressure**(-1/4)
        d = 0.077 - 0.038*dynamicPressure
        e = -0.874 - 0.299*dynamicPressure
        f = -0.055 + 0.124*dynamicPressure
    elif model == 'MP':
        a = -0.134 + 0.488*dynamicPressure**(-1/4)
        b = -0.581 - 0.225*dynamicPressure**(-1/4)
        c = -0.186 - 0.016*dynamicPressure**(-1/4)
        d = -0.014 + 0.096*dynamicPressure
        e = -0.814 - 0.811*dynamicPressure
        f = -0.050 + 0.168*dynamicPressure
    initModelCMD = '''eq = {a} + {b} x + {c} x^2 + {d} y + {e} y^2 + {f} x y - z^2 /. {{x -> r Sin[\[Theta]] Cos[\[Phi]], 
   y -> r Sin[\[Theta]] Sin[\[Phi]], z -> r Cos[\[Theta]]}};
        rInThetaPhi = r /. Solve[eq == 0, r];
        eqXYZ = {a} + {b} x + {c} x^2 + {d} y + {e} y^2 + {f} x y - z^2;
        eqXYZAtOriginXYZ0 = eqXYZ /. {{x->x+{x0}, y->y+{y0}, z->z+{z0}}};
        eqRTPAtOriginXYZ0 = eqXYZAtOriginXYZ0 /. {{x -> r Sin[\[Theta]] Cos[\[Phi]], 
   y -> r Sin[\[Theta]] Sin[\[Phi]], z -> r Cos[\[Theta]]}};
        rInThetaPhiAtOriginXYZ0 = r /. Solve[eqRTPAtOriginXYZ0 == 0, r];
        zInXY = z /. Solve[eqXYZ == 0, z]
             '''.format(a=a, b=b, c=c, d=d, e=e, f=f, x0=x0, y0=y0, z0=z0)
    wlSession.evaluate(initModelCMD)


def initWent11Model(wlSession, dynamicPressure=None, origin=np.zeros(3)):
    '''
    Purpose: obtain the location of Kronian bow shock using the model doi:10.1029/2010JA016349. The coordinate system is the aberrated KSM system.
    Parameters:
        wlSession: a wolframe session
        dynamicPressure: upstream solar wind dynamic pressure
        origin: the origin of the coordinate system in the unit of R_J
    '''
    x0, y0, z0 = origin
    c1 = 14.7
    c2 = 5.4
    epsilon = 0.84
    c1Pc2 = c1 * dynamicPressure**(-1/c2)
    initModelCMD = '''rInTheta = (1 + {epsilon}) {c1Pc2} / (1 + {epsilon} Cos[theta]);
                      eqRT = r - (1 + {epsilon}) {c1Pc2} / (1 + {epsilon} Cos[theta]);
                      eqXYZ = eqRT /. {{r -> Sqrt[x^2 + y^2 + z^2], Cos[theta] -> x/Sqrt[x^2 + y^2 + z^2]}};
                      eqXYZAtOriginXYZ0 = eqXYZ /. {{x->x+{x0}, y->y+{y0}, z->z+{z0}}};
                      eqRTPAtOriginXYZ0 = eqXYZAtOriginXYZ0 /. {{x -> r Sin[\[Theta]] Cos[\[Phi]], y -> r Sin[\[Theta]] Sin[\[Phi]], z -> r Cos[\[Theta]]}};
             '''.format(c1Pc2=c1Pc2, epsilon=epsilon, x0=x0, y0=y0, z0=z0)
    wlSession.evaluate(initModelCMD)


def initKanani10Model(wlSession, dynamicPressure=None, origin=np.zeros(3)):
    '''
    Purpose: obtain the location of Kronian magnetopause using the model doi:10.1029/2009JA014262. The coordinate system is KSM system.
    Parameters:
        wlSession: a wolframe session
        dynamicPressure: upstream solar wind dynamic pressure
        origin: the origin of the coordinate system in the unit of R_J
    '''
    x0, y0, z0 = origin
    a1 = 10.3
    a2 = 0.2
    a3 = 0.73
    a4 = 0.4
    r0 = a1 * dynamicPressure**(-a2)
    K = a3 + a4 * dynamicPressure
    initModelCMD = '''rInTheta = {r0} (2 / (1 + Cos[theta]))^{K};
                      eqRT = r - {r0} (2 / (1 + Cos[theta]))^{K};
                      eqXYZ = eqRT /. {{r -> Sqrt[x^2 + y^2 + z^2], Cos[theta] -> x/Sqrt[x^2 + y^2 + z^2]}};
                      eqXYZAtOriginXYZ0 = eqXYZ /. {{x->x+{x0}, y->y+{y0}, z->z+{z0}}};
                      eqRTPAtOriginXYZ0 = eqXYZAtOriginXYZ0 /. {{x -> r Sin[\[Theta]] Cos[\[Phi]], y -> r Sin[\[Theta]] Sin[\[Phi]], z -> r Cos[\[Theta]]}};
             '''.format(r0=r0, K=K, x0=x0, y0=y0, z0=z0)
    wlSession.evaluate(initModelCMD)


def initShue98Model(wlSession, Bz=None, dynamicPressure=None, origin=np.zeros(3)):
    '''
    Purpose: obtain the location of terrestrial magnetopause using the model doi:10.1029/98JA01103. The coordinate system is GSM/GSE system (I am not sure).
    Parameters:
        wlSession: a wolframe session
        Bz: in nT
        dynamicPressure: upstream solar wind dynamic pressure in nPa
        origin: the origin of the coordinate system in the unit of R_J
    '''
    x0, y0, z0 = origin
    r0 = (10.22 + 1.29*np.tanh(0.184 * (Bz + 8.14))) * dynamicPressure**(-1/6.6)
    alpha = (0.58 - 0.007*Bz) * (1+0.024*np.log(dynamicPressure))
    initModelCMD = '''rInTheta = {r0} (2 / (1 + Cos[theta]))^{alpha};
                      eqRT = r - {r0} (2 / (1 + Cos[theta]))^{alpha};
                      eqXYZ = eqRT /. {{r -> Sqrt[x^2 + y^2 + z^2], Cos[theta] -> x/Sqrt[x^2 + y^2 + z^2]}};
                      eqXYZAtOriginXYZ0 = eqXYZ /. {{x->x+{x0}, y->y+{y0}, z->z+{z0}}};
                      eqRTPAtOriginXYZ0 = eqXYZAtOriginXYZ0 /. {{x -> r Sin[\[Theta]] Cos[\[Phi]], y -> r Sin[\[Theta]] Sin[\[Phi]], z -> r Cos[\[Theta]]}};
             '''.format(r0=r0, alpha=alpha, x0=x0, y0=y0, z0=z0)
    wlSession.evaluate(initModelCMD)


def initGruesBeck18MarsBSModel(wlSession, origin=np.zeros(3)):
    '''
    Purpose: obtain the location of Kronian magnetopause using the model doi:10.1029/2018JA025366. The coordinate system is MSO system.
    Parameters:
        wlSession: a wolframe session
        origin: the origin of the coordinate system in the unit of R_J
    '''
    x0, y0, z0 = origin
    Ap = 0.049
    Bp = 0.157
    Cp = 0.153
    Dp = 0.026
    Ep = 0.012
    Fp = 0.051
    Gp = 0.566
    Hp = -0.031
    Ip = 0.019
    initModelCMD = '''eqXYZ = {Ap} x^2 + {Bp} y^2 + {Cp} z^2 + {Dp} x y +{Ep} y z + {Fp} x z + {Gp} x + {Hp} y + {Ip} z - 1;
                      eqRTP = eqXYZ /. {{x -> r Sin[\[Theta]] Cos[\[Phi]], y -> r Sin[\[Theta]] Sin[\[Phi]], z -> r Cos[\[Theta]]}};
                      rInThetaPhi = r /. Solve[eqRTP == 0, r];
                      eqXYZAtOriginXYZ0 = eqXYZ /. {{x->x+{x0}, y->y+{y0}, z->z+{z0}}};
                      eqRTPAtOriginXYZ0 = eqXYZAtOriginXYZ0 /. {{x -> r Sin[\[Theta]] Cos[\[Phi]], y -> r Sin[\[Theta]] Sin[\[Phi]], z -> r Cos[\[Theta]]}};
                      rInThetaPhiAtOriginXYZ0 = r /. Solve[eqRTPAtOriginXYZ0 == 0, r];
                      zInXY = z /. Solve[eqXYZ == 0, z]
             '''.format(Ap=Ap, Bp=Bp, Cp=Cp, Dp=Dp, Ep=Ep, Fp=Fp, Gp=Gp, Hp=Hp, Ip=Ip, x0=x0, y0=y0, z0=z0)
    wlSession.evaluate(initModelCMD)

#    stringList = []
#    for i in range(9):
#        stringList.append('{a}p={a}p'.format(a=chr(ord('A') + i)))
#    string = ', '.join(stringList)
#    print(string)



class shue98MPModel:
    '''This is deprecated, use the mathematica model'''
    def __init__(self, Dp, Bz):
        '''
        see doi:10.1029/98JA01103 for detail
        Parameters:
            Bz: magnetic field in nT
            Dp: dynamic pressure in nPa
        '''
        self.r0 = (10.22 + 1.29*np.tanh(0.184 * (Bz + 8.14))) * Dp**(-1/6.6)
        self.alpha = (0.58 - 0.007*Bz) * (1+0.024*np.log(Dp))

    def get_r(self, theta, phi):
        r = self.r0 * (2/(1+np.cos(theta)))**self.alpha
        return r
