##This will be where all the inputs live
import numpy as np
import math
import random
import scipy
import sys

from PyQt5.QtGui import QPixmap, QFont
from PyQt5.QtCore import Qt, QThread, QTimer, QObject, pyqtSignal, QEvent
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QComboBox, QSlider, QLabel, QPushButton, QProgressBar, QCheckBox, QLineEdit, QListWidgetItem, QListWidget, QDialog, QAction, QMenu, QToolBar, QMenuBar, QGridLayout
from PyQt5.QtGui import QIntValidator 
import pyqtgraph as pg ##will use this when i do nice frontend

from main import ThresholdKineticEnergy, TwoBodyCollision, CentreofMassKineticEnergy, TwoBodyDecay, ThreeBodyDecay, ThreeBodyDecayGenerator

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setWindowTitle("Kinematic calculator")
        self.mainLayout = QHBoxLayout()
        self.leftWidget = QWidget()
        self.inputLayout = QVBoxLayout()
        self.leftWidget.setLayout(self.inputLayout)
        self.showMaximized()
        self.rightWidget = QWidget()
        self.graphicLayout = QVBoxLayout()
        self.rightWidget.setLayout(self.graphicLayout)
        self.font = QFont()
        self.font.setPointSize(20)
        
        self.choice = 1 ##Default choice
        self.default = 0
        self.multiAngle = False
        #self.colours = ["Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "DarkRed", "DarkGreen", "DarkBlue", "DarkCyan", "DarkMagenta", "turquoise", "darkorchid", "chocolate", "cadetblue", "deeppink", "darkviolet"]
        self.colours = ["maroon", "red", "darkorange", "goldenrod", "olive", "forestgreen", "darkgreen", "lightseagreen", "teal", "dodgerblue", "blue", "darkslateblue", "indigo", "purple", "crimson", "black"]
        ## ^Used for when exporting on a white background
        
        #Starting to define the inputs
        
        self.choiceLabel = QLabel("Would you like to:" + "\n" + "1 Calculate the threshold energy for a reaction" + "\n" + "2 Calculate the reaction" + "\n" + "3 Calculate Ejectile lab angle distributions" + "\n" + "4 Calculate angle distributions for neutron evapouration " + "\n" + "5 Calculate the energy variance")
        self.choiceLabel.setFixedHeight(120)
        self.inputLayout.addWidget(self.choiceLabel) 
        self.choiceInput = QComboBox()
        self.choiceInput.addItems(["1", "2", "3", "4", "5", "6"])
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
        
        self.graphButtonBackward = QPushButton(text="Prev")
        self.graphButtonBackward.pressed.connect(self.Backward)
        self.graphButtonForward = QPushButton(text="Next")
        self.graphButtonForward.pressed.connect(self.Forward)

        self.viewEnergyDistribution = QPushButton("View Energy Dsitribution")
        self.viewEnergyDistribution.pressed.connect(self.ViewEnergyDistribution)        
        self.inputLayout.addWidget(self.viewEnergyDistribution)
        self.viewEnergyDistribution.hide()

        self.inputLayout.addStretch()
        
        self.mainLayout.addWidget(self.leftWidget)
        self.mainLayout.addWidget(self.rightWidget)
        container = QWidget()
        container.setLayout(self.mainLayout)
        self.setCentralWidget(container)
        
    def Choices(self):
        self.choice = self.choiceInput.currentText()
    def Backward(self):
        if hasattr(self, 'graphIndex'):
                
            if self.graphIndex == 0:
                print("What are you doing??? Are you crazy? You're going to get us all killed!!")
            else:
                self.graphicLayout.removeWidget(self.graph[self.graphIndex])
                self.graph[self.graphIndex].setParent(None)
                self.graphIndex -= 1
                self.graphicLayout.addWidget(self.graph[self.graphIndex])
        
    def Forward(self):
        if hasattr(self, 'graphIndex'):
            
            if self.graphIndex == (len(self.graph) -1):
                print("What are you doing??? Are you crazy? You're going to get us all killed!!")
            else:
                self.graphicLayout.removeWidget(self.graph[self.graphIndex])
                self.graph[self.graphIndex].setParent(None)
                self.graphIndex += 1
                self.graphicLayout.addWidget(self.graph[self.graphIndex])

    
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
            
    def ViewEnergyDistribution(self): #This will be for changing what plots are visible

        self.plt[self.graphIndex].clear()
        curve = pg.PlotCurveItem(self.thetaEjectile*180/np.pi, self.energyRange[self.graphIndex])
        self.plt[self.graphIndex].addItem(curve)
    def Run(self): ##This happens after the start button gets pressed
        if self.choice == "1":
            self.ThresholdKineticEnergyInput()
            
        elif self.choice == "2":
            
            recoilEk, ejectileEk, vRecoil, vEjectile, EkDecay, angleDecay, thetaEjectile = self.TwoBodyCollisionInput()
            if hasattr(self, 'graph'):
                if hasattr(self, 'graphIndex'):
                    self.graphicLayout.removeWidget(self.graph[self.graphIndex])
                else:
                    self.graphicLayout.removeWidget(self.graph)
            
            self.graph = pg.GraphicsLayoutWidget()
            self.graphicLayout.addWidget(self.graph)
            self.plt = self.graph.addPlot(row = 0, col = 0)
            if self.multiAngle == True:
                self.plt.addLegend()
                for i in range(len(EkDecay)):
                    curve = pg.PlotCurveItem(thetaEjectile*180/np.pi, EkDecay[i], pen={'color':self.colours[i], 'width':3}, name=i*45)
                    self.plt.addItem(curve)
            else:
                curve = pg.PlotCurveItem(thetaEjectile*180/np.pi, EkDecay)
                self.plt.addItem(curve)
            self.plt.setTitle("Recoil nucleus energy as a function of Ejectile's angle")
            self.plt.setLabel("left", "Energy /MeV")
            self.plt.setLabel("bottom", "Angle /deg")
            self.plt.getAxis("bottom").setTickFont(self.font)
            self.plt.getAxis("left").setTickFont(self.font)
        elif self.choice == "3":
            self.LabAngleDistInput()
            
            
        elif self.choice == "4":
            self.ThreeBodyDecaySimInput()
            
        elif self.choice == "5":
            if hasattr(self, 'graph'):
                if hasattr(self, 'graphIndex'):
                    self.graphicLayout.removeWidget(self.graph[self.graphIndex])
                    for i in range(len(self.graph)):
                        self.plt[i].clear()
                        
                else:
                    self.graphicLayout.removeWidget(self.graph)

            self.EnergyRangeInput()
        elif self.choice == "6":
            self.MultiExcitationGenerator()
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
            EkDecay, angleDecay = TwoBodyDecay(recoilEk, exRecoil, -16, thetaRecoil, mRecoil, 85661.69 , 2*939.57)
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
            else:
                i = 0
                EkDecay, angleDecay = ThreeBodyDecay(recoilEk, exRecoil, -16, thetaRecoil, mRecoil, 85661.69 , 939.57, 939.57, self.multiAngle, i)
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
        self.defaults = [1000, 2600, 5, 86542.44, 1876.139, 87479.78, 938.791, 18, 0]
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
        self.plt.addLegend(offset=(0,80))
        for i in range(len(energies)):
            EkBeamCM, EkTargetCm, EkCM, V = CentreofMassKineticEnergy(energies[i], mBeam, mTarget)
            recoilEk, ejectileEk, vRecoil, vEjectile, thetaRecoil, thetaEjectile, bestLabAngle = TwoBodyCollision(EkCM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile, V)
            labAngles[i] = bestLabAngle
            curve = pg.PlotCurveItem(thetaEjectile*180/np.pi, ejectileEk, pen={'color':self.colours[i], 'width':3}, name=energies[i])
            self.plt.addItem(curve)
        #The outputs are of shape [len(energies), len(angles)] where angles are the uniform distribution in the CoM frame

        self.plt.setTitle("Ejectile energy against lab angle for a range of lab energies")
        self.plt.setLabel("left", "Ejectile energy /MeV")
        self.plt.setLabel("bottom", "Ejectile angle /Deg")
        self.plt.getAxis("bottom").setTickFont(self.font)
        self.plt.getAxis("left").setTickFont(self.font)
        
    def ThreeBodyDecaySimInput(self): #Going to be used to return a range of energies and momenta for the decay nucleus for a given decay energy
        self.collisionArgs = ["EkMin", "EkMax", "Qvalue", "mBeam", "mTarget", "mRecoil", "mEjectile", "exRecoil", "exEjectile"]
        self.collisionDefaults = [800, 4000, 5, 86542.44, 1876.139, 87479.78, 938.791, 18, 0]
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
        self.decayDefaults = [85661.69, 939.57, 939.57] #These are for Nb92, neutron, neutron
        if self.default == 1:
            self.decayArgs = self.decayDefaults
        else:
            for i in range(len(self.decayArgs)):
                aux = input("Please input a " + self.decayArgs[i] + " (MeV):")
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
        mDecayRecoil, mDaughter1, mDaughter2 = self.decayArgs
        if hasattr(self, 'graph'):
            if hasattr(self, 'graphIndex'):
                self.graphicLayout.removeWidget(self.graph[self.graphIndex])
                for i in range(len(self.graph)):
                    self.plt[i].clear()
                    
            else:
                self.graphicLayout.removeWidget(self.graph)

        if EkMin == 0 or EkMin == EkMax: #using one energy
            self.graph = pg.GraphicsLayoutWidget()
            self.graphicLayout.addWidget(self.graph)
            self.plt = self.graph.addPlot(row = 0, col = 0)
            energies = EkMax
            EkBeamCM, EkTargetCM, EkCM, V = CentreofMassKineticEnergy(energies, mBeam, mTarget)
            recoilEk, ejectileEk, vRecoil, vEjectile, thetaRecoil, thetaEjectile, bestLabAngle = TwoBodyCollision(EkCM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile, V)
            #Multiple energies returned over a uniform CoM angle distribution from TwoBodyCollision
            EDaughterLab, daughterAngle = ThreeBodyDecay(recoilEk, exRecoil, Qvalue, thetaRecoil, mRecoil, mDecayRecoil, mDaughter1, mDaughter2, False, 0)
            #Just a single decay angle chosen
            curve = pg.PlotCurveItem(thetaEjectile*180/np.pi, EDaughterLab) #Plot of the decay nucleus' energy as a function of lab energy
            self.plt.addItem(curve)
        else:
            EStep = 200 # Step between energies for a range of them
            self.energies = np.arange(EkMin, EkMax, EStep) #Creates an array of beam energies to iterate through
            if self.multiAngle == True:
                
                self.energyRange = np.zeros((len(self.energies), 18000))
                self.plt = [[] for _ in range(len(self.energies))]
                self.graph = [[] for j in range(len(self.energies))]
                for j in range(len(self.energies)): ##Running multi-angle decay multiple times. 
                    self.graph[j] = pg.GraphicsLayoutWidget()
                    self.plt[j] = self.graph[j].addPlot(row = 0, col = 0)

                    EkBeamCM, EkTargetCM, EkCM, V = CentreofMassKineticEnergy(self.energies[j], mBeam, mTarget)
                    recoilEk, ejectileEk, vRecoil, vEjectile, thetaRecoil, self.thetaEjectile, bestLabAngle = TwoBodyCollision(EkCM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile, V)
                    self.EdaughterLab = np.zeros((8, len(recoilEk)))
                    self.angleDecay = np.zeros((8, len(recoilEk)))
                    
                    
                    for i in range(8): #One set of plots for One of the multi energies
                        EkDecaySingle, angleDecaySingle = ThreeBodyDecay(recoilEk, exRecoil, -16, thetaRecoil, mRecoil, 85661.69 , 939.57, 939.57, self.multiAngle, i)                
                    
                    # ax1.plot(thetaEjectile*180/np.pi, EkDecaySingle, label=('CoM Angle = ', i*45), color=colour[i])

                    # ax2.plot(thetaEjectile*180/np.pi, angleDecaySingle, label=('CoM Angle = ', i*45), color=colour[i], linestyle='dashed')
                        self.EdaughterLab[i] = EkDecaySingle
                        self.angleDecay[i] = angleDecaySingle
                        curve = pg.PlotCurveItem(self.thetaEjectile*180/np.pi, self.EdaughterLab[i], pen={'color':self.colours[i], 'width':3}, name=str(i*45)+"Deg")
                        self.plt[j].addLegend()
                        self.plt[j].addItem(curve)
                    firstVals = self.EdaughterLab[:,0]
                    minIndex = np.argmin(firstVals)
                    maxIndex = np.argmax(firstVals)
                    self.energyRange[j] = self.EdaughterLab[maxIndex] - self.EdaughterLab[minIndex]
                    self.plt[j].setTitle("Recoil Nucleus energy as a function of ejectile lab angle for multiple Decay angles. Lab energy of: " +  str(self.energies[j]))
                    self.plt[j].setLabel("left", "Energy /MeV")
                    self.plt[j].setLabel("bottom", "Ejectile Angle /deg")
                    self.plt[j].getAxis("bottom").setTickFont(self.font)
                    self.plt[j].getAxis("left").setTickFont(self.font)
                self.graphIndex = 0
                self.graphicLayout.addWidget(self.graphButtonBackward) ##Adds a button to toggle between graphs
                self.graphicLayout.addWidget(self.graphButtonForward)
                self.graphicLayout.addWidget(self.graph[0]) ##Adds the first graph
                self.viewEnergyDistribution.show()
                
                
                

            else:
                self.graph = pg.GraphicsLayoutWidget()
                self.graphicLayout.addWidget(self.graph)
                self.plt = self.graph.addPlot(row = 0, col = 0)
                for i in range(len(self.energies)):
                    EkBeamCM, EkTargetCM, EkCM, V = CentreofMassKineticEnergy(self.energies[i], mBeam, mTarget)
                    recoilEk, ejectileEk, vRecoil, vEjectile, thetaRecoil, thetaEjectile, bestLabAngle = TwoBodyCollision(EkCM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile, V)
                    #Gets the kinematic information for the (d,p) etc interaction
                    EdaughterLab, thetaDaughterLab = ThreeBodyDecay(recoilEk, exRecoil, Qvalue, thetaRecoil, mRecoil, mDecayRecoil, mDaughter1, mDaughter2, False, 0)
                    #Gets the kinematic information for the recoil nucleus' decay based on the outgoing kinematics from the previous collision interaction
                    curve = pg.PlotCurveItem(thetaEjectile*180/np.pi, EdaughterLab, pen={'color':self.colours[i], 'width':3})
                    self.plt.addItem(curve)
                self.plt.setTitle("Decay Recoil energy as a function of ejectile angle")
                self.plt.setLabel("left", "Decay Recoil Energy /MeV")
                self.plt.setLabel("bottom", "Ejectile angle /deg")
                self.plt.getAxis("bottom").setTickFont(self.font)
                self.plt.getAxis("left").setTickFont(self.font)
    
    def EnergyRangeInput(self, *args): 
##This will be for creating
##1: a set of recoil nucleus energies based on the scattering angle for a single input lab energy and 
##2: creating a range of decay energies based on A: CoM decay angle AND B: Energy distribution. 
##So it will be a 3 step process of energy ranges
        self.collisionArgs = ["Energy", "Qvalue", "mBeam", "mTarget", "mRecoil", "mEjectile", "exRecoil", "exEjectile"]
        self.decayArgs = ["mRecoil", "mDaughter1", "mDaughter2"]
        self.collisionDefaults = [2400, 5, 86542.44, 1876.139, 87479.78, 938.791, 18, 0]
        self.decayDefaults = [85661.69, 939.57, 939.57]
        if self.default == 1:
            if len(args) != 0:
                self.collisionDefaults[6] = args[0] ##allows for MultiExcitationGenerator to parse through excitations
            self.collisionArgs = self.collisionDefaults
            self.decayArgs = self.decayDefaults
        else:
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
            
            for i in range(len(self.decayArgs)):
                aux = input("Please input a " + self.decayArgs[i] + " (MeV):")
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
        EkLab, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile = self.collisionArgs
        mDecayRecoil, mDaughter1, mDaughter2 = self.decayArgs
        EkBeamCM, EkTargetCM, EkCM, V = CentreofMassKineticEnergy(EkLab, mBeam, mTarget) ##Initial EkCM calculation. Array size of 1
        recoilEk, ejectileEk, vRecoil, vEjectile, thetaRecoil, thetaEjectile, bestLabAngle = TwoBodyCollision(EkCM, Qvalue, mBeam, mTarget, mRecoil, mEjectile, exRecoil, exEjectile, V) #List of energies from two body collision, size 18000
        ## This bit will be the weird bit
        EDaughterLab = np.zeros((8, 1000, 18000))
        daughterAngle = np.zeros((8, 1000, 18000))
        EMin = np.zeros((8, 18000))
        EMax = np.zeros((8, 18000))
        self.graph = pg.GraphicsLayoutWidget()
        self.graphicLayout.addWidget(self.graph)
        self.plt = self.graph.addPlot(row = 0, col = 0)
        for i in range(8): ##going through 8 angles
            EDaughterLabSingle, daughterAngleSingle = ThreeBodyDecayGenerator(recoilEk, exRecoil, Qvalue, thetaRecoil, mRecoil, mRecoil, mDaughter1, mDaughter2, True, i)
            ##multiAngle argument is set to true regardless of the state of the checkbox
            ##calls the generator function 8 times, which then calls the three body function 1000 times. Final output of energies will be of shape 8, 1000, 18000
            ##The first dimension is the standard energy from initial scattering angle array
            ##The second dimension is for the CoM decay angle
            ##The Thirs is for a range of CoM decay energies
            EDaughterLab[i] = EDaughterLabSingle
            daughterAngle[i] = daughterAngleSingle
            maxEngIndex = np.argmax(EDaughterLabSingle[:,0])
            minEngIndex = np.argmin(EDaughterLabSingle[:,0])
            EMin[i] = EDaughterLabSingle[minEngIndex,:] ##Im gonna be honest i think im just performing magical operations rn
            EMax[i] = EDaughterLabSingle[maxEngIndex,:]
        lowestEnergyIndex = np.argmin(EMin[:,0])
        largestEnergyIndex = np.argmax(EMax[:,0])
        Minima = pg.PlotCurveItem(thetaEjectile*180/np.pi, EMin[lowestEnergyIndex], pen={'color':self.colours[2], 'width':3}, name="Lowest Energy")
        Maxima =  pg.PlotCurveItem(thetaEjectile*180/np.pi, EMax[largestEnergyIndex], pen={'color':self.colours[5], 'width':3}, name="Largest Energy")       
        median = (EMin[lowestEnergyIndex] + EMax[largestEnergyIndex])/2
        eDiffBeam = EkLab-median
        energyDifference = EMax[largestEnergyIndex] - EMin[lowestEnergyIndex]

        if len(args) == 0:
            difference = pg.PlotCurveItem(thetaEjectile*180/np.pi, eDiffBeam, pen={'color':self.colours[7], 'width':3}, name="Beam energy - Median")       
            Median = pg.PlotCurveItem(thetaEjectile*180/np.pi, median, pen={'color':self.colours[10], 'width':3}, name="Median Energy")       
            
            labEnergy = pg.InfiniteLine(EkLab, angle=0, pen={'width':3})
            self.plt.addLegend(offset=(0,200))
            self.plt.addItem(Minima)
            self.plt.addItem(Maxima)
            self.plt.addItem(Median)
            self.plt.addItem(labEnergy)
            variance = pg.PlotCurveItem(thetaEjectile*180/np.pi, energyDifference, pen={'color':self.colours[12], 'width':3}, name="Energy Variance")
            self.plt.addItem(difference)
            self.plt.addItem(variance)
            self.plt.setTitle("Energy variance for the decaying nucleus as a function of Ejectile angle")
            self.plt.setLabel("left", "Energy difference /MeV")
            self.plt.setLabel("bottom", "Ejectile angle /deg")
            self.plt.getAxis("bottom").setTickFont(self.font)
            self.plt.getAxis("left").setTickFont(self.font)
        
        return(median, eDiffBeam, EkLab, energyDifference)
        
    def MultiExcitationGenerator(self):
        self.args = ["Excitation min", "Excitation Max", "Angle of focus"]
        self.defaults = [16, 24, 150]
        if self.default == 1:
            self.args=self.defaults
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
        excitationEnergies = np.arange(self.args[0], self.args[1], 1)
        median = np.zeros((len(excitationEnergies), 18000))
        eDiffBeam = np.zeros((len(excitationEnergies), 18000))
        EkLab = np.zeros((len(excitationEnergies), 18000))
        energyDifference = np.zeros((len(excitationEnergies),18000))
        for i in range(len(excitationEnergies)):
            medianSingle, eDiffBeamSingle, EkLabSingle, energyDifferenceSingle = self.EnergyRangeInput(excitationEnergies[i])
            median[i] = medianSingle
            eDiffBeam[i] = eDiffBeamSingle
            EkLab[i] = EkLabSingle
            energyDifference[i] = energyDifferenceSingle
            self.graphicLayout.removeWidget(self.graph)
        angleIndex = round(18000 * self.args[2]/180)
        angle = self.args[2]
        self.graph = pg.GraphicsLayoutWidget()
        self.graphicLayout.addWidget(self.graph)
        self.plt = self.graph.addPlot(row = 0, col = 0)
        eDiffBeamData = eDiffBeam[:,angleIndex]
        curve = pg.PlotCurveItem(excitationEnergies, eDiffBeamData, pen={'color':self.colours[8], 'width':3}, name="Beam energy - Median recoil energy")
        data = energyDifference[:,angleIndex]/median[:,angleIndex]
        curve2 = pg.PlotCurveItem(excitationEnergies, data)
        curve3 = pg.PlotCurveItem(excitationEnergies, energyDifference[:,angleIndex], pen={'color':self.colours[10], 'width':3}, name="Energy variance")
        self.plt.addLegend()
        self.plt.addItem(curve)
        #self.plt.addItem(curve2)
        self.plt.addItem(curve3)
        self.plt.setTitle("Decay nucleus energy variance as a function of excitation energy at for a scattering angle of " + str(angle)) 
        self.plt.setLabel("bottom", "Excitation energy /MeV")
        self.plt.setLabel("left", "Energy /MeV")
        self.plt.getAxis("bottom").setTickFont(self.font)
        self.plt.getAxis("left").setTickFont(self.font)
        
        
app = QApplication(sys.argv)
window = MainWindow()
window.show()
app.exec()
    