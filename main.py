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

##defining kinematics modules here

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

def CentreofMassKineticEnergy(Ek, mBeam, mTarget):
    pbeam = np.sqrt(2*mBeam*Ek)
    vbeam = np.sqrt(pbeam**2/(mBeam**2 + pbeam**2))
    EBeam = np.sqrt(2*Ek*mBeam + mBeam**2)
    ##V here is the c.o.m velocity for galilean transformations (Will eventually implement to Lorentz transforms)
    V = (vbeam*mBeam)/(mTarget + mBeam)
    print(Ek)
    Ek_CM = Ek * mTarget/(mBeam+mTarget)
    print(Ek_CM)
    p_CM = np.sqrt((2*Ek_CM*mBeam*mTarget)/(mBeam+mTarget))
    #vBeam_CM = p_CM/mBeam
    vTarget_CM = -np.sqrt((p_CM**2)/(p_CM**2 + mTarget**2))
    print(vTarget_CM)
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
    f = mTarget/mEjectile   
    vEjectile = - V*np.sqrt(f*q*mRecoil/mBeam)# *sin(thetaCM)#velocity of the ejectile in the CM frame
    print(vEjectile, V, q, f)
    vEjectile_lab = vEjectile + V 
    gammaEjectile_CM = np.sqrt(1/(1-vEjectile**2))
    EkEjectile = (gammaEjectile_CM - 1) *mEjectile #SALAMI LID DONT FIT
    print(FinalEk_CM)
    print(EkEjectile)
    EkRecoil = FinalEk_CM - EkEjectile
    vRecoil = np.sqrt(EkRecoil*2*mRecoil)/mRecoil
    vRecoil_lab = V + vRecoil
    gammaEjectile = 1/np.sqrt(1-vEjectile_lab**2)
    gammaRecoil = 1/np.sqrt(1-vRecoil_lab**2)
    EkRecoil_lab = (gammaRecoil - 1)*mRecoil
    print(gammaEjectile, gammaRecoil)
    EkEjectile_lab = (gammaEjectile - 1)*mEjectile
    print(EkEjectile_lab + EkRecoil_lab)
    #EkEjectile_lab = 1/2 * mEjectile * vEjectile_lab**2
    #EkRecoil_lab = 1/2 * mRecoil * vRecoil_lab**2
    
    return(EkRecoil_lab, EkEjectile_lab, vRecoil_lab, vEjectile_lab)

def TwoBodyCollisionInput():
    args = ["Ek_lab", "Qvalue", "mBeam", "mTarget", "mRecoil", "mEjectile", "exRecoil", "exEjectile"]
    
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