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

def ThresholdKineticEnergy(Qvalue, ExTarget, targetMass, beamMass):
    return((ExTarget - Qvalue)*(beamMass + targetMass)/targetMass)

def ThresholdKineticEnergyInput():
    args = ["Qvalue", "excitedState", "targetMass", "beamMass"]
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

def CentreofMassKineticEnergy(Ek, beamMass, targetMass):
    pbeam = np.sqrt(2*beamMass*Ek)
    vbeam = pbeam/beamMass
    ##V here is the c.o.m velocity for galilean transformations (Will eventually implement to Lorentz transforms)
    V = beamMass*vbeam/(targetMass + beamMass)
    vbeam_CM = vbeam - V
    vtarget_CM = -V
    
    pbeam_CM = vbeam_CM * beamMass
    ptarget_CM = vtarget_CM * targetMass
    
    Ek_beam_CM = pbeam_CM**2/(2*beamMass)
    Ek_target_CM = ptarget_CM**2/(2*targetMass)
    print(Ek_beam_CM+Ek_target_CM)
    return(Ek_beam_CM, Ek_target_CM, Ek_beam_CM+Ek_target_CM, V)

def TwoBodyDecay(Ek_CM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, ExRecoil, ExEjectile, V): #The two body decay (also valid for general two body scattering) works here since by definition, the total momentum in the CoM frame is 0
    #V is the beam velocity in CoM frame    
    FinalEk_CM = Ek_CM + Qvalue - ExRecoil - ExEjectile
    if FinalEk_CM < 0:
        print("The input energy does not meet the threshold for this reaction. Stopping...")
        sys.exit()
    q = 1 + (Qvalue - ExRecoil - ExEjectile)/Ek_CM
    print(q)
    f = mTarget/mEjectile   
    vEjectile = - V*np.sqrt(f*q*mRecoil/mBeam)# *sin(thetaCM)#velocity of the ejectile in the CM frame
    print(V, vEjectile)
    vEjectile_lab = vEjectile + V 
    
    EkEjectile = 1/2 * vEjectile **2 * mEjectile #SALAMI LID DONT FIT
    print(EkEjectile, vEjectile)
    EkRecoil = FinalEk_CM - EkEjectile
    print(EkRecoil)
    vRecoil = np.sqrt(EkRecoil*2*mRecoil)/mRecoil
    vRecoil_lab = V + vRecoil
    EkEjectile_lab = 1/2 * mEjectile * vEjectile_lab**2
    EkRecoil_lab = 1/2 * mRecoil * vRecoil_lab**2
    
    return(EkRecoil_lab, EkEjectile_lab, vRecoil_lab, vEjectile_lab)

def TwoBodyCollisionInput():
    args = ["Ek_lab", "Qvalue", "mBeam", "mTarget", "mDaughter1", "mDaughter2", "exDaughter1", "exDaughter2"]
    
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
    
    Daughter1Ek, Daughter2Ek, vDaughter1, vDaughter2 = TwoBodyDecay(Ek_CM, args[1], args[2], args[3], args[4], args[5], args[6], args[7], V)
    
    return(Daughter1Ek, Daughter2Ek, vDaughter1, vDaughter2)
    
    
choice = input("Would you like to:" + "\n" + "1 Calculate the threshold energy for a reaction" + "\n" + "2 Calculate the reaction" + "\n" + ":")
if choice == "1":
    print(ThresholdKineticEnergyInput())
elif choice == "2":
    print(TwoBodyCollisionInput())
else:
    print("go home jade youre sleepy")