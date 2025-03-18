from PyQt5 import QtCore, QtGui, QtWidgets 
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas)

import columnOutput
import columnInput 
import sys
import columnStart 


class ColumnOutputMain(QtWidgets.QWidget):
    def __init__(self, Results, Graphs):
        super(ColumnOutputMain, self).__init__()
        self.Results = Results
        self.Graphs = Graphs
        self.ui = columnOutput.Ui_Form()
        self.ui.setupUi(self)
        self.displayResults()
        self.show()
        
    
    def displayResults(self):
        #Main Results
        concResults =    self.Results[0]
        reinfStrain =   self.Results[1]
        reinfStress =   self.Results[2]
        reinfUtilRatio = self.Results[3]
        intForce =      self.Results[4]
        RelDiff =       self.Results[5]
        RelDev =        self.Results[6]
        loop =          self.Results[7]
        
        graph1 =        self.Graphs[0]      #Concrete Strain
        graph2 =        self.Graphs[1]      #Concrete Stress
        graph3 =        self.Graphs[2]      #Nx Convergence
        graph4 =        self.Graphs[3]      #Mx Convergence
        graph5 =        self.Graphs[4]      #My Convergence    
        
            #Concrete Results
        concMaxStrain = concResults[0]
        concMaxStress = concResults[1]
        concUtilRatio = concResults[2]
        self.ui.concMaxStrain.setText(str(concMaxStrain))
        self.ui.concMaxStress.setText(str(concMaxStress))
        self.ui.concUtilRatio.setText(str(concUtilRatio))
        
            #Reinforcement Strain
        reinfStrainTens = reinfStrain[0]
        reinfStrainComp = reinfStrain[1]
        self.ui.reinfStrainTens.setText(str(reinfStrainTens))
        self.ui.reinfStrainComp.setText(str(reinfStrainComp))
        
            #Reinforcement Stress
        reinfStressTens = reinfStress[0]
        reinfStressComp = reinfStress[1]
        self.ui.reinfStressTens.setText(str(reinfStressTens))
        self.ui.reinfStressComp.setText(str(reinfStressComp))
        
            #Reinforcement Utility Ratio
        reinfUtilRatioTens = reinfUtilRatio[0]
        reinfUtilRatioComp = reinfUtilRatio[1]
        self.ui.reinfUtilRatioTens.setText(str(reinfUtilRatioTens))
        self.ui.reinfUtilRatioComp.setText(str(reinfUtilRatioComp))
        
            #internal Force
        nx = intForce[0]
        mx = intForce[1]
        my = intForce[2]
        self.ui.Nx.setText(str(nx))
        self.ui.Mx.setText(str(mx))
        self.ui.My.setText(str(my))
        
            #loop
        loop
        self.ui.iteration.setText(str(loop))
        
        
            #raphs
        self.displayGraphs(graph1, self.ui.Graph1VerticalLayout)
        self.displayGraphs(graph2, self.ui.Graph2VerticalLayout)
        self.displayGraphs(graph3, self.ui.Graph3VerticalLayout)
        self.displayGraphs(graph4, self.ui.Graph4VerticalLayout)
        self.displayGraphs(graph5, self.ui.Graph5VerticalLayout)
       
    
       
    def displayGraphs (self, graph, layout):
        
        self.ui.canvas = FigureCanvas(graph)
        layout.addWidget(self.ui.canvas)
        



