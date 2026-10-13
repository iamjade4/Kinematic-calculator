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
from frontend import MainWindow

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

# def CSC(theta, gamma):
#     return(((1 + gamma**2 + 2*gamma*np.cos(theta)))**(3/2)/abs(1 + gamma*np.cos(theta))) ##from catford's lectures

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

def TwoBodyDecay(Ek, Eex, Q, theta, mParent, mDaughter1, mDaughter2):
    ##Convert to CM frame:
    pLab = np.sqrt(2*mParent*Ek)
    vLab = np.sqrt(pLab**2/(pLab**2 + mParent**2))
    V = vLab #CoM velocity is just the velocity of the parent nucleus before decay
    Ek_CM = Eex + Q
    ##split between daughters using conservation of momentum
    pCM = np.sqrt((2*Ek_CM*mDaughter1*mDaughter2)/(mDaughter1 + mDaughter2))
    ##isotropic release of neutrons (for now have them decay as one)
    thetaCM = np.random.rand(1)*2*np.pi
    vDaughter1_CM = pCM/mDaughter1
    vDaughter2_CM = -pCM/mDaughter2
    vDaughter1OnAxisCM, vDaughter2OnAxisCM = (vDaughter1_CM, vDaughter2_CM) * np.cos(thetaCM)
    vDaughter1OffAxisCM, vDaughter2OffAxisCM = (vDaughter1_CM, vDaughter2_CM) * np.sin(thetaCM)
    #vDaughter2 is already negative, so the angles don't need flipping
    
    vDaughter1OnAxisLab = (vDaughter1OnAxisCM+ V) * np.cos(theta) ##theta is the initial scattering angle from the two body collision. Should be 0 if there was no prior collision
    vDaughter1OffAxisLab = vDaughter1OffAxisCM * np.sin(theta)
    vDaughter2OnAxisLab = (vDaughter2OnAxisCM + V) * np.cos(theta)
    vDaughter2OffAxisLab = vDaughter2OffAxisCM * np.sin(theta)
    
    ##Off axis stays the same
    daughter1Angle = np.arctan(vDaughter1OffAxisLab/vDaughter1OnAxisLab)*180/np.pi
    daughter2Angle = np.arctan(vDaughter2OffAxisLab/vDaughter2OnAxisLab)*180/np.pi
    plt.hist(daughter1Angle, bins=180)
    plt.title("Decay Recoil nucleus angles (deg)")
    plt.show()
    plt.hist(daughter2Angle, bins=180)
    plt.title("Decay Recoil particles angles (deg)")
    plt.show()
    
    #print("Angles:", daughter1Angle, daughter2Angle)
    
def ThreeBodyDecay(Ek, Eex, Q, theta, mParent, mDaughter1, mDaughter2, mDaughter3): #This will be used for 2n emission
    pLab = np.sqrt(2*mParent*Ek)
    vLab = np.sqrt(pLab**2/(pLab**2 + mParent**2))
    V = vLab #CoM velocity is just the velocity of the parent nucleus before decay
    Ek_CM = Eex + Q
    
    

def TwoBodyCollision(Ek_CM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, ExRecoil, ExEjectile, V): #The two body decay (also valid for general two body scattering) works here since by definition, the total momentum in the CoM frame is 0
    #V is the beam velocity in CoM frame    
    FinalEk_CM = Ek_CM + Qvalue - ExRecoil - ExEjectile
    if FinalEk_CM < 0:
        print("The input energy does not meet the threshold for this reaction. Stopping...")
        sys.exit()
    q = 1 + (Qvalue - ExRecoil - ExEjectile)/Ek_CM
    ## I DONT LIKE q
    print("C.o.M EK:", FinalEk_CM)
    f = mTarget/mEjectile   
    ## for now just create a basic distribution
    thetaCM = np.arange(0, 180, 0.01)*np.pi/180 #generates thetas from 0-180
    vEjectile = V*np.sqrt(f*q*mRecoil/mBeam)
    # ratio = 1/(np.sqrt(f*q*mRecoil/mBeam))
    # crossSections = CSC(thetaCM, ratio)
    # angularDist = AngDist(thetaCM, crossSections)
    #print("Angles", angularDist)
    vEjectileOnAxis = -vEjectile * np.cos(thetaCM)#velocity of the ejectile in the CM frame
    vEjectileOffAxis = vEjectile * np.sin(thetaCM)
    ##An angle of 180 corresponds to fully forward scattering of the ejectile in the CM frame 
    gammaEjectile_CM = np.sqrt(1/(1-vEjectile**2))
    pEjectileOffAxis = gammaEjectile_CM*mEjectile*vEjectileOffAxis
    pRecoilOffAxis = -pEjectileOffAxis
    pEjectileOnAxis = gammaEjectile_CM*mEjectile*vEjectileOnAxis
    pRecoilOnAxis = -pEjectileOnAxis
    vRecoilOnAxis = pRecoilOnAxis/mRecoil
    vRecoilOffAxis = pRecoilOffAxis/mRecoil ##pRecoil is so low<
    
    ##from now on it should all be array math
    vEjectile_labOnAxis = vEjectileOnAxis + V 
    vEjectile_labOffAxis = vEjectileOffAxis
    vEjectile_lab = np.sqrt(vEjectile_labOffAxis**2 + vEjectile_labOnAxis**2)
    mask = np.where(vEjectile_labOnAxis < 0, 1, 0)
    print(mask)
    thetaEjectileLab = np.arctan(vEjectile_labOffAxis/vEjectile_labOnAxis) + mask*np.pi
    
    #print(V, vEjectile, vEjectile_labOnAxis)
    vRecoilOnAxislab = vRecoilOnAxis + V
    vRecoilOffAxislab = vRecoilOffAxis
    vRecoil_lab = np.sqrt(vRecoilOffAxislab**2 + vRecoilOnAxislab**2)
    thetaRecoilLab= np.arctan(vRecoilOffAxislab/vRecoilOnAxislab)
    #print("lab angles", ThetaEjectile_lab*180/np.pi, ThetaRecoil_lab*180/np.pi)
    plt.hist(thetaEjectileLab*180/np.pi, bins=180)
    hist = np.histogram(thetaEjectileLab*180/np.pi, bins=180)
    plt.title("Lab angle distribution")
    plt.xlabel("Angle /Deg")
    #window.DisplayHist(hist[0], hist[1][:-1])
    
    maxangle_index = np.argmax(hist[0]) ##finds the bin number with the most counts
    print("Optimal lab angle:", hist[1][maxangle_index])
    plt.show()
    plt.plot(thetaEjectileLab*180/np.pi, thetaCM)
    plt.title("Lab frame angle as a function of COM angle")
    plt.xlabel("Angle /Deg")
    plt.show()
    
    gammaEjectile = Gamma(vEjectile_lab)
    gammaRecoil = Gamma(vRecoil_lab)
    EkRecoil_lab = (gammaRecoil - 1)*mRecoil
    EkEjectile_lab = (gammaEjectile - 1)*mEjectile
    #print(EkEjectile_lab + EkRecoil_lab) #slightly smaller and I'm not sure why
    #EkEjectile_lab = 1/2 * mEjectile * vEjectile_lab**2
    #EkRecoil_lab = 1/2 * mRecoil * vRecoil_lab**2
    return(EkRecoil_lab, EkEjectile_lab, vRecoil_lab, vEjectile_lab, thetaRecoilLab, thetaEjectileLab)

def TwoBodyCollisionInput():
    args = ["Ek_lab", "Qvalue", "mBeam", "mTarget", "mRecoil", "mEjectile", "exRecoil", "exEjectile"]
    defaults = [1600, 5, 86542.44, 1876.139, 87479.78, 938.791, 18, 0]
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
    exRecoil = args[6]
    Ek_beam_CM, Ek_target_CM, Ek_CM, V = CentreofMassKineticEnergy(args[0], args[2], args[3])
    recoilEk, ejectileEk, vRecoil, vEjectile, thetaRecoil, thetaEjectile = TwoBodyCollision(Ek_CM, args[1], args[2], args[3], args[4], args[5], args[6], args[7], V)
    TwoBodyDecay(recoilEk, exRecoil, -16, thetaRecoil, args[4], 85661.69 , 2*939.57)
    #ThreeBodyDecay(recoilEk, exRecoil, -16, thetaRecoil, args[4], 85661.69 , 939.57, 939.57)
    return(recoilEk, ejectileEk, vRecoil, vEjectile)


# if __name__ == "__main__":
#     app = QApplication()
#     window = MainWindow()
#     window.resize(1000,800)
#     window.show()  
#     app.exec()
    
choice = input("Would you like to:" + "\n" + "1 Calculate the threshold energy for a reaction" + "\n" + "2 Calculate the reaction" + "\n" + ":")

if choice == "1":
    print(ThresholdKineticEnergyInput())
elif choice == "2":
    TwoBodyCollisionInput()
else:
    print("go home jade youre sleepy")
