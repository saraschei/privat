from PyQt5 import QtCore, QtGui, QtWidgets 
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas)

import beamOutput
import beamInput 
import sys
import beam 


class BeamOutputMain(QtWidgets.QWidget):
    def __init__(self, Results, Graphs):
        super(BeamOutputMain, self).__init__()
        self.Results = Results
        self.Graphs = Graphs
        self.ui = beamOutput.Ui_Form()
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
        #graph5 =        self.Graphs[4]      
        
            #Concrete Results
        concMaxStrain = concResults[0]
        concMaxStress = concResults[1]
        concUtilRatio = concResults[2]
        self.ui.concMaxStrain.setText(str(concMaxStrain))
        self.ui.concMaxStress.setText(str(concMaxStress))
        self.ui.concUtilRatio.setText(str(concUtilRatio))
        
            #Reinforcement Strain
        reinfStrainAsx1 = reinfStrain[0]
        reinfStrainAsx2 = reinfStrain[1]
        self.ui.reinfStrainAsx1.setText(str(reinfStrainAsx1))
        self.ui.reinfStrainAsx2.setText(str(reinfStrainAsx2))
        
            #Reinforcement Stress
        reinfStressAsx1 = reinfStress[0]
        reinfStressAsx2 = reinfStress[1]
        self.ui.reinfStressAsx1.setText(str(reinfStressAsx1))
        self.ui.reinfStressAsx2.setText(str(reinfStressAsx2))
        
            #Reinforcement Utility Ratio
        reinfUtilRatioAsx1 = reinfUtilRatio[0]
        reinfUtilRatioAsx2 = reinfUtilRatio[1]
        self.ui.reinfUtilRatioAsx1.setText(str(reinfUtilRatioAsx1))
        self.ui.reinfUtilRatioAsx2.setText(str(reinfUtilRatioAsx2))
        
            #internal Force
        nx = intForce[0]
        mx = intForce[1]
        self.ui.Nx.setText(str(nx))
        self.ui.Mx.setText(str(mx))
        
            #loop
        loop
        self.ui.iteration.setText(str(loop))
        
        
            #raphs
        self.displayGraphs(graph1, self.ui.Graph1VerticalLayout)
        self.displayGraphs(graph2, self.ui.Graph2VerticalLayout)
        self.displayGraphs(graph3, self.ui.Graph3VerticalLayout)
        self.displayGraphs(graph4, self.ui.Graph4VerticalLayout)
        #self.displayGraphs(graph5, self.ui.Graph5VerticalLayout)
       
    
       
    def displayGraphs (self, graph, layout):
        
        self.ui.canvas = FigureCanvas(graph)
        layout.addWidget(self.ui.canvas)
        



