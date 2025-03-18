from PyQt5 import QtCore, QtGui, QtWidgets
import sys
import shellInput
import shell
import shellOutputMain
import noConvergenceDialog
import nonNumericalInputDialog


class ShellInputMain(QtWidgets.QMainWindow):
    def __init__(self):
        super(ShellInputMain, self).__init__()
        self.ui = shellInput.Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.calculateButton.clicked.connect(self.calculate)
        self.show()
    
    def calculate(self):
        try:
            nx = float(self.ui.nx.text())
            ny = float(self.ui.ny.text())
            nxy = float(self.ui.nxy.text())
            mx = float(self.ui.mx.text())
            my = float(self.ui.my.text())
            mxy = float(self.ui.mxy.text())
            fck = float(self.ui.fck.text())
            fyk = float(self.ui.fyk.text())
            gamma_c = float(self.ui.gamma_c.text())
            gamma_s = float(self.ui.gamma_s.text())
            concreteModel = int(self.ui.concreteModelComboBox.currentIndex()) + 1   #takes index of combobox and adds 1 
                                                                                    #parab-rec: 1, bilinear: 2
            Asx1 = float(self.ui.Asx1.text())
            Asx2 = float(self.ui.Asx2.text())
            Asy1 = float(self.ui.Asy1.text())
            Asy2 = float(self.ui.Asy2.text())
            h = float(self.ui.h.text())
            c1 = float(self.ui.c1.text())
            c2 = float(self.ui.c2.text())
            n = int(self.ui.n.text())
            Esx1 = float(self.ui.Esx1.text())
            Esx2 = float(self.ui.Esx2.text())
            Esy1 = float(self.ui.Esy1.text())
            Esy2 = float(self.ui.Esy2.text())
            epsud = float(self.ui.epsud.text())
            beta = float(self.ui.beta.text())
            maxIt = int(self.ui.maxIt.text())
            vc = float(self.ui.vc.text())
            
        except ValueError:
            self.nonNumericalInput()
            return
        
        
        shellObj = shell.Shell (nx, ny, nxy, mx, my, mxy, fck, fyk, gamma_c, gamma_s, concreteModel, \
                        Asx1, Asx2, Asy1, Asy2, h, c1, c2, n, Esx1, Esx2, Esy1, Esy2, epsud, beta, maxIt, vc)   #Create a Shell object
        
        Results = shellObj.iteration()
        loop = Results[7] 
        convergence = Results[11]       #Boolean
        
        if convergence == True:
            loop = Results[7] 
            epspciMat = Results[8]
            sigpciMat = Results[9]
            sigpsjMat = Results[10]
            Graphs = shellObj.plot(loop, epspciMat, sigpciMat, sigpsjMat)      #Contains 10 Figures
                                                                                       
            shellOutputObj = shellOutputMain.ShellOutputMain(Results, Graphs)  #Initialize a sellOutputMain object
            
            self.ui.window = shellOutputObj
            
        else:
     
            self.openDialog( loop)
            
            
            
    def openDialog(self, loop):
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
    MainWindow = ShellInputMain()
    MainWindow.show()
    sys.exit(app.exec_()) 
