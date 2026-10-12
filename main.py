import numpy as np
import math
import random
import scipy
import sys

##For the inevitable frontend development
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt, QThread, QTimer, QObject, pyqtSignal, QEvent
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QComboBox, QSlider, QLabel, QPushButton, QProgressBar, QCheckBox, QLineEdit, QListWidgetItem, QListWidget, QDialog, QAction, QMenu, QToolBar, QMenuBar
from PyQt5.QtGui import QIntValidator 
#import pyqtgraph as pg ##will use this when i do nice frontend
import matplotlib.pyplot as plt

##Importing defined materials
from materials import *

##defining kinematics modules here


def RelativisticVelocityCalc(p, m):
    return(np.sqrt((p**2)/(p**2 + m**2)))

def Gamma(v):
    return(1/np.sqrt(1-v**2))

def ThresholdKineticEnergy(Qvalue, ExTarget, mTarget, mBeam):
    return((ExTarget - Qvalue)*(mBeam + mTarget)/mTarget)

def ThresholdKineticEnergyInput():
    args = ["Qvalue", "excitedState", "mTarget", "mBeam"]
    for i in range(len(args)):
        aux = input("Please input a " + args[i] + " (MeV):")    
        while True:
            try:
                aux = float(aux)
            except ValueError:
                aux = input("Please input a " + args[i] + " (MeV):")
            except TypeError:
                aux = input("Please input a " + args[i] + " (MeV):")
            else:
                break
            
        args[i] = aux
    return(ThresholdKineticEnergy(*args))

def CSC(theta, gamma):
    return(((1 + gamma**2 + 2*gamma*np.cos(theta)))**(3/2)/abs(1 + gamma*np.cos(theta)))

def AngDist(angles, crossSections):
    normalised = crossSections/sum(crossSections)
    dist = np.random.choice(angles, len(angles), p = normalised) 
    return(dist) #This should return a weighted list of random angles (weighted to the cross section)

def CentreofMassKineticEnergy(Ek, mBeam, mTarget):
    pbeam = np.sqrt(2*mBeam*Ek)
    vbeam = RelativisticVelocityCalc(pbeam, mBeam)
    EBeam = np.sqrt(2*Ek*mBeam + mBeam**2)
    ##V here is the c.o.m velocity for galilean transformations (Will eventually implement to Lorentz transforms)
    V = (vbeam*mBeam)/(mTarget + mBeam)
    Ek_CM = Ek * mTarget/(mBeam+mTarget)
    p_CM = np.sqrt((2*Ek_CM*mBeam*mTarget)/(mBeam+mTarget))
    #vBeam_CM = p_CM/mBeam
    vTarget_CM = -np.sqrt((p_CM**2)/(p_CM**2 + mTarget**2))
    #gammaBeam_CM = 1/np.sqrt(1-vBeam_CM**2)
    #gammaTarget_CM = 1/np.sqrt(1-vTarget_CM**2)
    Ek_beam_CM = p_CM**2/(2*mBeam)
    Ek_target_CM = p_CM**2/(2*mTarget)
    return(Ek_beam_CM, Ek_target_CM, Ek_beam_CM+Ek_target_CM, V)
def TwoBodyDecay(Ek_CM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, ExRecoil, ExEjectile, V): #The two body decay (also valid for general two body scattering) works here since by definition, the total momentum in the CoM frame is 0
    #V is the beam velocity in CoM frame    
    FinalEk_CM = Ek_CM + Qvalue - ExRecoil - ExEjectile
    if FinalEk_CM < 0:
        print("The input energy does not meet the threshold for this reaction. Stopping...")
        sys.exit()
    q = 1 + (Qvalue - ExRecoil - ExEjectile)/Ek_CM
    ## I DONT LIKE q
    print("EK:", FinalEk_CM)
    f = mTarget/mEjectile   
    ## for now just create a basic distribution
    thetaCM = np.arange(0, 181, 0.1)*np.pi/180 #generates thetas from 0-180
    vEjectile = V*np.sqrt(f*q*mRecoil/mBeam)
    ratio = 1/(np.sqrt(f*q*mRecoil/mBeam))
    crossSections = CSC(thetaCM, ratio)
    angularDist = AngDist(thetaCM, crossSections)
    print("Angles", angularDist)
    vEjectile_on_axis = -vEjectile * np.cos(angularDist)#velocity of the ejectile in the CM frame
    vEjectile_off_axis = vEjectile * np.sin(angularDist)
    ##An angle of 180 corresponds to fully forward scattering of the ejectile in the CM frame 
    gammaEjectile_CM = np.sqrt(1/(1-vEjectile**2))
    pEjectile_off_axis = gammaEjectile_CM*mEjectile*vEjectile_off_axis
    pRecoil_off_axis = -pEjectile_off_axis
    pEjectile_on_axis = gammaEjectile_CM*mEjectile*vEjectile_on_axis
    pRecoil_on_axis = -pEjectile_on_axis
    vRecoil_on_axis = pRecoil_on_axis/mRecoil
    vRecoil_off_axis = pRecoil_off_axis/mRecoil ##pRecoil is so low compared to m so we don't need to worry about relativity
    
    EkEjectile = (gammaEjectile_CM - 1) *mEjectile #SALAMI LID DONT FIT
    EkRecoil = FinalEk_CM - EkEjectile
    
    ##from now on it should all be array math
    vEjectile_lab_on_axis = vEjectile_on_axis + V 
    vEjectile_lab_off_axis = vEjectile_off_axis
    vEjectile_lab = np.sqrt(vEjectile_lab_off_axis**2 + vEjectile_lab_on_axis**2)
    ThetaEjectile_lab = np.arctan(vEjectile_lab_off_axis/vEjectile_lab_on_axis)
    print(V, vEjectile, vEjectile_lab_on_axis)
    vRecoil_on_axis_lab = vRecoil_on_axis + V
    vRecoil_off_axis_lab = vRecoil_off_axis
    vRecoil_lab = np.sqrt(vRecoil_off_axis_lab**2 + vRecoil_on_axis_lab**2)
    ThetaRecoil_lab = np.arctan(vRecoil_off_axis_lab/vRecoil_on_axis_lab)
    #print("lab angles", ThetaEjectile_lab*180/np.pi, ThetaRecoil_lab*180/np.pi)
    plt.hist(ThetaEjectile_lab*180/np.pi, bins=50)
    #plt.plot(thetaCM, crossSections)
    plt.show()
    
    gammaEjectile = Gamma(vEjectile_lab)
    gammaRecoil = Gamma(vRecoil_lab)
    EkRecoil_lab = (gammaRecoil - 1)*mRecoil
    EkEjectile_lab = (gammaEjectile - 1)*mEjectile
    print(EkEjectile_lab + EkRecoil_lab) #slightly smaller and I'm not sure why
    #EkEjectile_lab = 1/2 * mEjectile * vEjectile_lab**2
    #EkRecoil_lab = 1/2 * mRecoil * vRecoil_lab**2
    
    return(EkRecoil_lab, EkEjectile_lab, vRecoil_lab, vEjectile_lab)

def TwoBodyCollisionInput():
    args = ["Ek_lab", "Qvalue", "mBeam", "mTarget", "mRecoil", "mEjectile", "exRecoil", "exEjectile"]
    defaults = [800, 5, 86542.44, 1876.139, 87479.78, 938.791, 18, 0]
    default = input("Would you like to use defaults? y/n:").lower()
    if default == "y":
        args = defaults
    else:
        for i in range(len(args)):
            aux = input("Please input a " + args[i] + " (MeV):")
            while True:
                try:
                    aux = float(aux)
                except ValueError:
                    aux = input("Please input a " + args[i] + " (MeV):")
                except TypeError:
                    aux = input("Please input a " + args[i] + " (MeV):")
                else:
                    break
            args[i] = aux
    Ek_beam_CM, Ek_target_CM, Ek_CM, V = CentreofMassKineticEnergy(args[0], args[2], args[3])
    RecoilEk, EjectileEk, vRecoil, vEjectile = TwoBodyDecay(Ek_CM, args[1], args[2], args[3], args[4], args[5], args[6], args[7], V)
    
    return(RecoilEk, EjectileEk, vRecoil, vEjectile)
    
    
choice = input("Would you like to:" + "\n" + "1 Calculate the threshold energy for a reaction" + "\n" + "2 Calculate the reaction" + "\n" + ":")
if choice == "1":
    print(ThresholdKineticEnergyInput())
elif choice == "2":
    print(TwoBodyCollisionInput())
else:
    print("go home jade youre sleepy")