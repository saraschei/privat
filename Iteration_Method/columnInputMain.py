from PyQt5 import QtCore, QtGui, QtWidgets
import sys
import columnInput
import columnStart
import columnCase1
import columnCase2
import columnCase3
import columnCase4
import columnOutputMain
import noConvergenceDialog
import nonNumericalInputDialog


class ColumnInputMain(QtWidgets.QMainWindow):
    def __init__(self):
        super(ColumnInputMain, self).__init__()
        self.ui = columnInput.Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.calculateButton.clicked.connect(self.calculate)
        self.show()
    
    def calculate(self):
        try:
            nx = float(self.ui.nx.text())
            mx = float(self.ui.mx.text())
            my = float(self.ui.my.text())
            fck = float(self.ui.fck.text())
            fyk = float(self.ui.fyk.text())
            gammac = float(self.ui.gamma_c.text())
            gammas = float(self.ui.gamma_s.text())
            concModel = int(self.ui.concreteModelComboBox.currentIndex()) + 1   #takes index of combobox and adds 1 
                                                                                    #parab-rec: 1, bilinear: 2
            Asx1 = float(self.ui.Asx1.text())
            Asx2 = float(self.ui.Asx2.text())
            Asy1 = float(self.ui.Asy1.text())
            Asy2 = float(self.ui.Asy2.text())
            b = float(self.ui.b.text())
            h = float(self.ui.h.text())
            cx1 = float(self.ui.cx1.text())
            cx2 = float(self.ui.cx2.text())
            cy1 = float(self.ui.cy1.text())
            cy2 = float(self.ui.cy2.text())
            n = int(self.ui.n.text())
            Es = float(self.ui.Es.text())
            epsud = float(self.ui.epsud.text())
            beta = float(self.ui.beta.text())
            maxIt = int(self.ui.maxIt.text())
        
        except ValueError:
            self.nonNumericalInput()
            return


        columnStartObj = columnStart.ColumnStart(mx, my, b, h )   #Create a columnStart object
                         
        iterationCase = columnStartObj.iterationCase()
        ms = iterationCase[0]
        alfa1 = iterationCase[1]
        alfa2 = iterationCase[2]
        alfa3 = iterationCase[3]
        case = iterationCase[4]
        
        if case == 1:
            columnObj = columnCase1.ColumnCase1(nx, mx, my, fck, fyk, gammac, gammas, concModel, Asx1, Asx2, Asy1, Asy2,\
                  b, h, cx1, cx2, cy1, cy2, n, Es, epsud, beta, maxIt, ms, alfa1, alfa2, alfa3)
        elif case == 2: 
            columnObj = columnCase2.ColumnCase2(nx, mx, my, fck, fyk, gammac, gammas, concModel, Asx1, Asx2, Asy1, Asy2,\
                  b, h, cx1, cx2, cy1, cy2, n, Es, epsud, beta, maxIt, ms, alfa1, alfa2, alfa3)
        elif case == 3:
            columnObj = columnCase3.ColumnCase3(nx, mx, my, fck, fyk, gammac, gammas, concModel, Asx1, Asx2, Asy1, Asy2,\
                  b, h, cx1, cx2, cy1, cy2, n, Es, epsud, beta, maxIt, ms, alfa1, alfa2, alfa3)
        else:
            columnObj = columnCase4.ColumnCase4(nx, mx, my, fck, fyk, gammac, gammas, concModel, Asx1, Asx2, Asy1, Asy2,\
                  b, h, cx1, cx2, cy1, cy2, n, Es, epsud, beta, maxIt, ms, alfa1, alfa2, alfa3)
        
        Results = columnObj.iteration()
        loop = Results[7] 
        convergence = Results[11]       #Boolean
        if convergence:
            loop = Results[7] 
            epsciMat = Results[8]
            sigciMat = Results[9]
            sigsjMat = Results[10]
            Graphs = columnObj.plot(loop, epsciMat, sigciMat, sigsjMat)               #Contains 4 Graphs, ConcStrain, ConcStress, Nx, Mx
            
            columnOutputObj = columnOutputMain.ColumnOutputMain(Results, Graphs)           #Initialize a BeamOutputMain object
        
            self.ui.window = columnOutputObj
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
    MainWindow = ColumnInputMain()
    MainWindow.show()
    sys.exit(app.exec_()) 

