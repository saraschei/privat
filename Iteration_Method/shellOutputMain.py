from PyQt5 import QtCore, QtGui, QtWidgets 
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas)

import shellOutput
 


class ShellOutputMain(QtWidgets.QWidget):
    def __init__(self, Results, Graphs):
        super(ShellOutputMain, self).__init__()
        self.Results = Results
        self.Graphs = Graphs
        self.ui = shellOutput.Ui_Form()
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
        
        #Graph results
        straincFig1 =        self.Graphs[0]      
        straincFig2 =        self.Graphs[1]      
        stresscFig1 =        self.Graphs[2]      
        stresscFig2 =        self.Graphs[3]      
        #stressFig1 =        self.Graphs[4]      
        #stressFig2 =        self.Graphs[5]      
        forceFig1 =        self.Graphs[4]      
        forceFig2 =        self.Graphs[5]     
        forceFig3 =        self.Graphs[6]      
        momentFig1 =       self.Graphs[7]      
        momentFig2 =       self.Graphs[8]      
        momentFig3 =       self.Graphs[9]      

       
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
        reinfStrainAsy1 = reinfStrain[2]
        reinfStrainAsy2 = reinfStrain[3]
        self.ui.reinfStrainAsx1.setText(str(reinfStrainAsx1))
        self.ui.reinfStrainAsx2.setText(str(reinfStrainAsx2))
        self.ui.reinfStrainAsy1.setText(str(reinfStrainAsy1))
        self.ui.reinfStrainAsy2.setText(str(reinfStrainAsy2))
        
        
            #Reinforcement Stress
        reinfStressAsx1 = reinfStress[0]
        reinfStressAsx2 = reinfStress[1]
        reinfStressAsy1 = reinfStress[2]
        reinfStressAsy2 = reinfStress[3]
        self.ui.reinfStressAsx1.setText(str(reinfStressAsx1))
        self.ui.reinfStressAsx2.setText(str(reinfStressAsx2))
        self.ui.reinfStressAsy1.setText(str(reinfStressAsy1))
        self.ui.reinfStressAsy2.setText(str(reinfStressAsy2))
        
            #Reinforcement Utility Ratio
        reinfUtilRatioAsx1 = reinfUtilRatio[0]
        reinfUtilRatioAsx2 = reinfUtilRatio[1]
        reinfUtilRatioAsy1 = reinfUtilRatio[2]
        reinfUtilRatioAsy2 = reinfUtilRatio[3]
        self.ui.reinfUtilRatioAsx1.setText(str(reinfUtilRatioAsx1))
        self.ui.reinfUtilRatioAsx2.setText(str(reinfUtilRatioAsx2))
        self.ui.reinfUtilRatioAsy1.setText(str(reinfUtilRatioAsy1))
        self.ui.reinfUtilRatioAsy2.setText(str(reinfUtilRatioAsy2))
        
            #internal Force
        nx = intForce[0]
        ny = intForce[1]
        nxy = intForce[2]             
        mx = intForce[3]
        my = intForce[4]
        mxy = intForce[5]
        self.ui.nx.setText(str(nx))
        self.ui.ny.setText(str(ny))
        self.ui.nxy.setText(str(nxy))
        self.ui.mx.setText(str(mx))
        self.ui.my.setText(str(my))
        self.ui.mxy.setText(str(mxy))
        
            #loop
        loop
        self.ui.iteration.setText(str(loop))
        
        
            #raphs
        self.displayGraphs(straincFig1, self.ui.Graph1VerticalLayout)
        self.displayGraphs(straincFig2, self.ui.Graph2VerticalLayout)
        self.displayGraphs(stresscFig1, self.ui.Graph3VerticalLayout)
        self.displayGraphs(stresscFig2, self.ui.Graph4VerticalLayout)
        self.displayGraphs(forceFig1, self.ui.Graph5VerticalLayout)
        self.displayGraphs(forceFig2, self.ui.Graph6VerticalLayout)
        self.displayGraphs(forceFig3, self.ui.Graph7VerticalLayout)
        self.displayGraphs(momentFig1, self.ui.Graph8VerticalLayout)
        self.displayGraphs(momentFig2, self.ui.Graph9VerticalLayout)
        self.displayGraphs(momentFig3, self.ui.Graph10VerticalLayout)
       
       
    def displayGraphs (self, graph, layout):
        
        self.ui.canvas = FigureCanvas(graph)
        layout.addWidget(self.ui.canvas)
              



