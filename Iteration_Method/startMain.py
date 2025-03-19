from PyQt5 import QtCore, QtGui, QtWidgets
import sys
import start
import shellInputMain
import beamInputMain
import columnInputMain


class StartMain(QtWidgets.QMainWindow):
    def __init__(self):
        super(StartMain, self).__init__()
        self.ui = start.Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.runButton.clicked.connect(self.run)
        self.show()
    
    def run(self):
        if self.ui.beamRadioButton.isChecked():
            beamInputObj = beamInputMain.BeamInputMain()
            self.ui.window = beamInputObj
            
        if self.ui.shellRadioButton.isChecked():
            shellInputObj = shellInputMain.ShellInputMain()
            self.ui.window = shellInputObj       

        if self.ui.columnRadioButton.isChecked():
            columnInputObj = columnInputMain.ColumnInputMain()
            self.ui.window = columnInputObj
                                                               
        
        
        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = StartMain()
    MainWindow.show()
    sys.exit(app.exec_()) 
