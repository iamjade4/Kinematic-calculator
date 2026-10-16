##This will be where all the inputs live
import numpy as np
import math
import random
import scipy
import sys

from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt, QThread, QTimer, QObject, pyqtSignal, QEvent
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QComboBox, QSlider, QLabel, QPushButton, QProgressBar, QCheckBox, QLineEdit, QListWidgetItem, QListWidget, QDialog, QAction, QMenu, QToolBar, QMenuBar
from PyQt5.QtGui import QIntValidator 
import pyqtgraph as pg ##will use this when i do nice frontend

from main import ThresholdKineticEnergy, TwoBodyCollision, CentreofMassKineticEnergy, TwoBodyDecay, ThreeBodyDecay

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setWindowTitle("Kinematic calculator")
        self.mainLayout = QHBoxLayout()
        self.leftWidget = QWidget()
        self.inputLayout = QVBoxLayout()
        self.leftWidget.setLayout(self.inputLayout)
        
        self.rightWidget = QWidget()
        self.graphicLayout = QHBoxLayout()
        self.rightWidget.setLayout(self.graphicLayout)
        
        self.choice = 1 ##Default choice
        self.multiAngle = False
        self.colours = ["Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "DarkRed", "DarkGreen", "DarkBlue", "DarkCyan", "DarkMagenta", "White", "darkorchid", "chocolate", "cadetblue", "deeppink", "darkviolet"]
        
        #Starting to define the inputs
        
        self.choiceLabel = QLabel("Would you like to:" + "\n" + "1 Calculate the threshold energy for a reaction" + "\n" + "2 Calculate the reaction" + "\n" + "3 Calculate lab angle distributions" + "\n" + "4 Calculate angle distributions for neutron evapouration")
        self.choiceLabel.setFixedHeight(90)
        self.inputLayout.addWidget(self.choiceLabel) 
        self.choiceInput = QComboBox()
        self.choiceInput.addItems(["1", "2", "3", "4", "5"])
        self.choiceInput.currentTextChanged.connect(self.Choices)
        self.inputLayout.addWidget(self.choiceInput)
        
        self.defaultSwitch = QCheckBox("Toggle defaults")
        self.defaultSwitch.stateChanged.connect(self.Defaults)
        self.inputLayout.addWidget(self.defaultSwitch)
        
        self.multiSwitch = QCheckBox("Toggle Multiangles")
        self.multiSwitch.stateChanged.connect(self.MultiAngle)
        self.inputLayout.addWidget(self.multiSwitch)
        
        self.startButton = QPushButton(text="Start")
        self.startButton.pressed.connect(self.Run)
        self.inputLayout.addWidget(self.startButton)
        
        self.inputLayout.addStretch()
        
        self.mainLayout.addWidget(self.leftWidget)
        self.mainLayout.addWidget(self.rightWidget)
        container = QWidget()
        container.setLayout(self.mainLayout)
        self.setCentralWidget(container)
        
    def Choices(self):
        self.choice = self.choiceInput.currentText()
            
    def Defaults(self, state):
        if state == Qt.Checked:
            self.default = 1
            if self.choice == "2":
                self.defaults = [1600, 5, 86542.44, 1876.139, 87479.78, 938.791, 18, 0, 3]
        else:
            self.default = 2
    
    def MultiAngle(self, state):
        if state == Qt.Checked:
            self.multiAngle = True
        else:
            self.multiAngle = False
    
    def Run(self): ##This happens after the start button gets pressed
        if self.choice == "1":
            self.ThresholdKineticEnergyInput()
            
        elif self.choice == "2":
            
            recoilEk, ejectileEk, vRecoil, vEjectile, EkDecay, angleDecay, thetaEjectile = self.TwoBodyCollisionInput()
            if hasattr(self, 'graph'):
                self.graphicLayout.removeWidget(self.graph)
            
            self.graph = pg.GraphicsLayoutWidget()
            self.graphicLayout.addWidget(self.graph)
            self.plt = self.graph.addPlot(row = 0, col = 0)
            for i in range(len(EkDecay)):
                curve = pg.PlotCurveItem(thetaEjectile*180/np.pi, EkDecay[i], pen=self.colours[i])
                self.plt.addItem(curve)
            self.plt.setTitle("Recoil nucleus energy as a function of Ejectile's angle")
            self.plt.setLabel("left", "Energy /MeV")
            self.plt.setLabel("bottom", "Angle /deg")
        elif self.choice == "3":
            self.LabAngleDistInput()
            
            
        elif self.choice == "4":
            self.ThreeBodyDecaySimInput()
        else:  
            print("Closing...")
    
    def TwoBodyCollisionInput(self): #Used for testing a single lab beam energy and generating a range of angles for both the ejectile in the collision and the recoil nucleus
        self.args = ["Ek_lab", "Qvalue", "mBeam", "mTarget", "mRecoil", "mEjectile", "exRecoil", "exEjectile", "DecayType"]
        self.Defaults(self.defaultSwitch.checkState())
        if self.default == 1:
            self.args = self.defaults
        else:
            for i in range(len(self.args)):
                aux = input("Please input a " + self.args[i] + " (MeV):")
                while True:
                    try:
                        aux = float(aux)
                    except ValueError:
                        aux = input("Please input a " + self.args[i] + " (MeV):")
                    except TypeError:
                        aux = input("Please input a " + self.args[i] + " (MeV):")
                    else:
                        break
                self.args[i] = aux
        EkLab, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile, DecayType = self.args #Giving the args names to make it easier to track them
        EkBeamCm, EkTargetCM, EkCM, V = CentreofMassKineticEnergy(EkLab, mBeam, mTarget) #Calling the module to obtain CoM information for the collision
        #CoMKE only returns a single EKCM, there are no angular distributions yet
        recoilEk, ejectileEk, vRecoil, vEjectile, thetaRecoil, thetaEjectile, bestLabAngle = TwoBodyCollision(EkCM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile, V) #Obtains the post collision particle information
        #TwoBodyCollision returns an array of energies and particles as it uniformly generates CoM angles between 0 and 180
        if DecayType == 2:
            TwoBodyDecay(recoilEk, exRecoil, -16, thetaRecoil, mRecoil, 85661.69 , 2*939.57)
        elif DecayType == 3:
            #fig, ax1 = plt.subplots()
            #ax2 = ax1.twinx()
            if self.multiAngle == True:
                #colour = ["red", "blue", "orange", "fuchsia", "yellow", "black", "cyan", "purple"]
                EkDecay = np.zeros((8, len(recoilEk)))
                angleDecay = np.zeros((8, len(recoilEk)))
                for i in range(8): #Loops through the three body decay over a range of 8 angles to obtain more angular data
                    EkDecaySingle, angleDecaySingle = ThreeBodyDecay(recoilEk, exRecoil, -16, thetaRecoil, mRecoil, 85661.69 , 939.57, 939.57, self.multiAngle, i)                
                    
                    # ax1.plot(thetaEjectile*180/np.pi, EkDecaySingle, label=('CoM Angle = ', i*45), color=colour[i])

                    # ax2.plot(thetaEjectile*180/np.pi, angleDecaySingle, label=('CoM Angle = ', i*45), color=colour[i], linestyle='dashed')
                    EkDecay[i] = EkDecaySingle
                    angleDecay[i] = angleDecaySingle
        return(recoilEk, ejectileEk, vRecoil, vEjectile, EkDecay, angleDecay, thetaEjectile)
        
    

    def ThresholdKineticEnergyInput(self):
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
        return(print(ThresholdKineticEnergy(*args)))

    def LabAngleDistInput(self): #Used for calculating the most common ejectile angle over a range of energies
        self.args = ["EkMin", "EkMax", "Qvalue", "mBeam", "mTarget", "mRecoil", "mEjectile", "exRecoil", "exEjectile"]
        self.defaults = [800, 2400, 5, 86542.44, 1876.139, 87479.78, 938.791, 18, 0]
        if self.default == 1:
            self.args = self.defaults
        else:
            for i in range(len(self.args)):
                aux = input("Please input a " + self.args[i] + " (MeV):")
                while True:
                    try:
                        aux = float(aux)
                    except ValueError:
                        aux = input("Please input a " + self.args[i] + " (MeV):")
                    except TypeError:
                        aux = input("Please input a " + self.args[i] + " (MeV):")
                    else:
                        break
                self.args[i] = aux
        EkMin, EkMax, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile = self.args
        energies = np.arange(EkMin, EkMax, 100)
        labAngles = np.zeros(len(energies))
        if hasattr(self, 'graph'):
            self.graphicLayout.removeWidget(self.graph)
        
        self.graph = pg.GraphicsLayoutWidget()
        self.graphicLayout.addWidget(self.graph)
        self.plt = self.graph.addPlot(row = 0, col = 0)
        for i in range(len(energies)):
            EkBeamCM, EkTargetCm, EkCM, V = CentreofMassKineticEnergy(energies[i], mBeam, mTarget)
            recoilEk, ejectileEk, vRecoil, vEjectile, thetaRecoil, thetaEjectile, bestLabAngle = TwoBodyCollision(EkCM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile, V)
            labAngles[i] = bestLabAngle
            curve = pg.PlotCurveItem(thetaEjectile*180/np.pi, ejectileEk, pen=self.colours[i])
            self.plt.addItem(curve)
        #The outputs are of shape [len(energies), len(angles)] where angles are the uniform distribution in the CoM frame
        self.plt.setTitle("Ejectile energy against lab angle for a range of lab energies")
        self.plt.setLabel("left", "Ejectile energy /MeV")
        self.plt.setLabel("bottom", "Ejectile angle /Deg")
        
    def ThreeBodyDecaySimInput(self): #Going to be used to return a range of energies and momenta for the decay nucleus for a given decay energy
        self.collisionArgs = ["EkMin", "EkMax", "Qvalue", "mBeam", "mTarget", "mRecoil", "mEjectile", "exRecoil", "exEjectile"]
        self.collisionDefaults = [800, 2000, 5, 86542.44, 1876.139, 87479.78, 938.791, 18, 0]
        if self.default == 1:
            self.collisionArgs = self.collisionDefaults
        else:
            print("Input an EkMin of 0 or of the same value as EkMax if you would like to run in single energy mode")
            for i in range(len(self.collisionArgs)):
                aux = input("Please input a " + self.collisionArgs[i] + " (MeV):")
                while True:
                    try:
                        aux = float(aux)
                    except ValueError:
                        aux = input("Please input a " + self.collisionArgs[i] + " (MeV):")
                    except TypeError:
                        aux = input("Please input a " + self.collisionArgs[i] + " (MeV):")
                    else:
                        break
                self.collisionArgs[i] = aux
        EkMin, EkMax, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile = self.collisionArgs
        ## Having two different sets of inputs  
        self.decayArgs = ["mRecoil", "mDaughter1", "mDaughter2"]
        self.decayDefaults = [85661.69, 939.57, 939.57]
        if self.default == 1:
            self.decayArgs = self.decayDefaults
        else:
            for i in range(len(self.decayArgs)):
                aux = input("Please input a" + self.decayArgs[i] + " (MeV):")
                while True:
                    try:
                        aux = float(aux)
                    except ValueError:
                        aux = input("Please input a " + self.decayArgs[i] + " (MeV):")
                    except TypeError:
                        aux = input("Please input a " + self.decayArgs[i] + " (MeV):")
                    else:
                        break
                self.decayArgs[i] = aux
        mRecoil, mDaughter1, mDaughter2 = self.decayArgs
        if hasattr(self, 'graph'):
            self.graphicLayout.removeWidget(self.graph)
        
        self.graph = pg.GraphicsLayoutWidget()
        self.graphicLayout.addWidget(self.graph)
        self.plt = self.graph.addPlot(row = 0, col = 0)
        if EkMin == 0 or EkMin == EkMax: #using one energy
            energies = EkMax
            EkBeamCM, EkTargetCM, EkCM, V = CentreofMassKineticEnergy(energies, mBeam, mTarget)
            recoilEk, ejectileEk, vRecoil, vEjectile, thetaRecoil, thetaEjectile, bestLabAngle = TwoBodyCollision(EkCM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile, V)
            #Multiple energies returned over a uniform CoM angle distribution from TwoBodyCollision
            EDaughterLab, daughterAngle = ThreeBodyDecay(recoilEk, exRecoil, Qvalue, thetaRecoil, mRecoil, mRecoil, mDaughter1, mDaughter2, False, 0)
            #Just a single decay angle chosen
            curve = pg.PlotCurveItem(thetaEjectile*180/np.pi, EDaughterLab) #Plot of the decay nucleus' energy as a function of lab energy
            self.plt.addItem(curve)
        else:
            EStep = 200 # Step between energies for a range of them
            energies = np.arange(EkMin, EkMax, EStep) #Creates an array of beam energies to iterate through

            for i in range(len(energies)):
                EkBeamCM, EkTargetCM, EkCM, V = CentreofMassKineticEnergy(energies[i], mBeam, mTarget)
                recoilEk, ejectileEk, vRecoil, vEjectile, thetaRecoil, thetaEjectile, bestLabAngle = TwoBodyCollision(EkCM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile, V)
                #Gets the kinematic information for the (d,p) etc interaction
                EdaughterLab, thetaDaughterLab = ThreeBodyDecay(recoilEk, exRecoil, Qvalue, thetaRecoil, mRecoil, mRecoil, mDaughter1, mDaughter2, False, 0)
                #Gets the kinematic information for the recoil nucleus' decay based on the outgoing kinematics from the previous collision interaction
                curve = pg.PlotCurveItem(thetaEjectile*180/np.pi, EdaughterLab, pen=self.colours[i])
                self.plt.addItem(curve)
        self.plt.setTitle("Decay Recoil energy as a function of ejectile angle")
        self.plt.setLabel("left", "Decay Recoil Energy /MeV")
        self.plt.setLabel("bottom", "Ejectile angle /deg")
                    
app = QApplication(sys.argv)
window = MainWindow()
window.show()
app.exec()
    