##frontend FUUUUCK I HATE FRONTEND SO MUCHHHHHH
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt, QThread, QTimer, QObject, pyqtSignal, QEvent
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QComboBox, QSlider, QLabel, QPushButton, QProgressBar, QCheckBox, QLineEdit, QListWidgetItem, QListWidget, QDialog, QAction, QMenu, QToolBar, QMenuBar
from PyQt5.QtGui import QIntValidator 
import pyqtgraph as pg ##will use this when i do nice frontend

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setWindowTitle("Kinematic Calc")
        self.main_layout = QHBoxLayout()
        
    def DisplayHist(self, X, Y):
        self.graph = pg.GraphicsLayoutWidget()
        curve = pg.PlotCurveItem(X, Y, StepMode=True)
        self.main_layout.addWidget(self.graph)
        self.plt = self.graph.addPlot(row=0, column=0)
        self.plt.setXRange(0, 1024, padding = 0)
        self.plt.addItem(curve)
        self.plt.setTitle("Angles")
        self.plt.setLabel("left", "Counts")
        self.plt.setLabel("bottom", "Angle")
