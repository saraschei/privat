from PyQt5 import QtCore, QtGui, QtWidgets
import sys
import beamInput
import beam
import beamOutputMain
import noConvergenceDialog
import nonNumericalInputDialog


class BeamInputMain(QtWidgets.QMainWindow):
    def __init__(self):
        super(BeamInputMain, self).__init__()
        self.ui = beamInput.Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.calculateButton.clicked.connect(self.calculate)
        self.show()
    
    def calculate(self):
        try:
            nx = float(self.ui.nx.text())
            mx = float(self.ui.mx.text())
            fck = float(self.ui.fck.text())
            fyk = float(self.ui.fyk.text())
            gamma_c = float(self.ui.gamma_c.text())
            gamma_s = float(self.ui.gamma_s.text())
            concreteModel = int(self.ui.concreteModelComboBox.currentIndex()) + 1   #takes index of combobox and adds 1 
                                                                                    #parab-rec: 1, bilinear: 2
            Asx1 = float(self.ui.Asx1.text())
            Asx2 = float(self.ui.Asx2.text())
            b = float(self.ui.b.text())
            h = float(self.ui.h.text())
            c1 = float(self.ui.c1.text())
            c2 = float(self.ui.c2.text())
            n = int(self.ui.n.text())
            Esx1 = float(self.ui.Esx1.text())
            Esx2 = float(self.ui.Esx2.text())
            epsud = float(self.ui.epsud.text())
            beta = float(self.ui.beta.text())
            maxIt = int(self.ui.maxIt.text())
        
        except ValueError:
            self.nonNumericalInput()
            return

        beamObj = beam.Beam(nx, mx, fck, fyk, gamma_c, gamma_s, concreteModel, \
                        Asx1, Asx2, b, h, c1, c2, n, Esx1, Esx2, epsud, beta, maxIt)   #Create a Beam object
        
        Results = beamObj.iteration()
        loop = Results[7] 
        convergence = Results[11]       #Boolean
        if convergence:
            loop = Results[7] 
            epsciMat = Results[8]
            sigciMat = Results[9]
            sigsjMat = Results[10]
            Graphs = beamObj.plot(loop, epsciMat, sigciMat, sigsjMat)               #Contains 4 Graphs, ConcStrain, ConcStress, Nx, Mx
            
            beamOutputObj = beamOutputMain.BeamOutputMain(Results, Graphs)           #Initialize a BeamOutputMain object
        
            self.ui.window = beamOutputObj
        else:
            self.nonConvergenceDialog(loop)
            
    def nonConvergenceDialog(self,loop):       #method for opening Non Convergence dialog
        dialog = QtWidgets.QDialog()
        dialog.ui = noConvergenceDialog.Ui_Dialog()
        dialog.ui.setupUi(dialog)
        dialog.ui.iterationLabel.setText(str(loop))
        dialog.exec_()
        
        
    def nonNumericalInput(self):
        dialog = QtWidgets.QDialog()
        dialog.ui = nonNumericalInputDialog.Ui_Dialog()
        dialog.ui.setupUi(dialog)
        dialog.exec_()
    
        
        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = BeamInputMain()
    MainWindow.show()
    sys.exit(app.exec_()) 
