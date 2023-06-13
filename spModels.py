import numpy as np
#import constants as con
import cdflib    # see github.com/MAVENSDC/cdflib
#import tarfile
from datetime import datetime
from datetime import timedelta
import os
import sys
from wolframclient.language import wl
from wolframclient.evaluation import WolframLanguageSession
import dataAnalysisTools as dat


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
