# -*- coding: utf-8 -*-

import sys
from PyQt5.QtCore import pyqtSlot, pyqtSignal, QThread, Qt, QSize
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QMainWindow, QFileDialog
from PyQt5.QtGui import QPixmap, QIcon, QTextCursor
from PyQt5.uic import loadUi
import numpy as np

import calfem.ui as cfui
import femModel as fm
import gc
import guiHelpMessages as ghm

class SolverThread(QThread):
    """Class to handle model solver calculation in the background"""

    def __init__(self, solver, paramstudy):
        """Class constructor"""
        QThread.__init__(self)
        self.solver = solver
        self.paramstudy = paramstudy

    def __del__(self):
        self.wait()

    def run(self):
        if self.paramstudy:
            self.solver.executeParamStudy()

        else:
            self.solver.execute()

class ImageWindow(QMainWindow):
    """Image Window used to specifically show a preview of the image"""

    def __init__(self):
        """Constructor"""
        super(QMainWindow, self).__init__()
        self.filename = None
        self.ui = loadUi('femgui_image.ui', self)

    def updateImage(self, svg):
        # self.ui.imageLabel.setPixmap(QPixmap(svg).scaled(450, 450, Qt.KeepAspectRatio, transformMode=Qt.SmoothTransformation))
        self.ui.imageLabel.setPixmap(
            QIcon(svg).pixmap(420, 420)) #.scaledToWidth(450, Qt.SmoothTransformation))

class MainWindow(QMainWindow):
    """MainWindow class that manages our main window"""

    def __init__(self):
        """Constructor"""
        super(QMainWindow, self).__init__()

        # --- Store a reference to the application instance in the class
        self.app = app

        # --- Store other necessary reference values
        self.filename = ""
        self.inputData = None
        self.outputData = None
        self.report = None
        self.vis = None
        self.figsOn = [0, 0, 0, 0]
        self.calcDone = False
        self.solver = None
        self.solverThread = None
        self.runparamstudy = 0

        # --- Read interface from the file
        self.ui = loadUi('femgui_final.ui', self)
        self.imageGui = ImageWindow()

        # --- Link controls to event methods
        self.ui.actionNew.triggered.connect(self.onActionNew)
        self.ui.actionOpen.triggered.connect(self.onActionOpen)
        self.ui.actionSave.triggered.connect(self.onActionSave)
        self.ui.actionSaveAs.triggered.connect(self.onActionSaveAs)
        self.ui.actionExit.triggered.connect(self.onActionExit)
        self.ui.actionLoadDXF.triggered.connect(self.onActionLoadDXF)
        self.ui.actionExecute.triggered.connect(self.onActionExecute)
        self.ui.showGeometryButton.clicked.connect(self.onShowGeometry)
        self.ui.showMeshButton.clicked.connect(self.onShowMesh)
        self.ui.showNodalValuesButton.clicked.connect(self.onShowNodalValues)
        self.ui.showElementValuesButton.clicked.connect(self.onShowElementValues)
        self.ui.elSizeSlider.sliderReleased.connect(self.onElSizeSliderMoved)
        self.ui.paramButton.clicked.connect(self.onExecuteParamStudy)
        self.ui.displaydxfButton.clicked.connect(self.onDisplayImage)
        self.ui.clearScreenButton.clicked.connect(self.onScreenClear)
        self.ui.reportEdit.keyReleaseEvent = self.handleKeyRelease
        self.ui.unitsComboBox.currentTextChanged.connect(self.onUnitsChanged)

        # --- Make sure to show the window
        self.ui.show()
        self.ui.raise_()

    def initModel(self):
        """Initialize a new model with default values"""
        self.inputData = fm.InputData()
        self.inputData.load("default.json")        # Load default values
        self.outputData = fm.OutputData()
        self.ui.reportEdit.insertPlainText("\n>>")

    def checkParse(self, number):
        """Checks if the input is a number (vs. an empty string), returns the number if it exists, else
        return -1. If the given input is -1, return an empty string"""
        if type(number) == str:
            if number == "":
                return -1
            elif number[0] == ".":
                return number
            for i in np.linspace(ord('0'), ord('9'), 10):
                if ord(number[0]) == i:
                    return number
            return -1
        elif number == -1:
            return ""
        else:
            return number


    def updateControls(self):
        """Fill in the controls with values from the model"""
        self.ui.d1NameEdit.setText(str(self.inputData.d[0][0]))
        self.ui.d1Edit.setText(str(self.checkParse(self.inputData.d[0][1])))
        self.ui.d2NameEdit.setText(str(self.inputData.d[1][0]))
        self.ui.d2Edit.setText(str(self.checkParse(self.inputData.d[1][1])))
        self.ui.d3NameEdit.setText(str(self.inputData.d[2][0]))
        self.ui.d3Edit.setText(str(self.checkParse(self.inputData.d[2][1])))
        self.ui.tEdit.setText(str(self.inputData.t))
        self.ui.bcEdit.setText(str(self.inputData.bp[1][0]))
        self.ui.bcNameEdit.setText(str(self.inputData.bp[0][0]))
        self.ui.fEdit.setText(str(self.inputData.fp[1][0]))
        self.ui.fNameEdit.setText(str(self.inputData.fp[0][0]))
        self.ui.fAngleEdit.setText(str(self.inputData.fp[2][0]))
        self.ui.aEdit.setText(str(self.checkParse(self.inputData.c['a'])))
        self.ui.bEdit.setText(str(self.checkParse(self.inputData.c['b'])))
        self.ui.aEndEdit.setText(str(self.checkParse(self.inputData.c['aEnd'])))
        self.ui.bEndEdit.setText(str(self.checkParse(self.inputData.c['bEnd'])))
        self.ui.aNameEdit.setText(str(self.inputData.c['aName']))
        self.ui.bNameEdit.setText(str(self.inputData.c['bName']))
        self.ui.elasticEdit.setText(str(self.inputData.E))
        self.ui.poissonEdit.setText(str(self.inputData.v))
        self.ui.paramStep.setValue(int(self.inputData.paramSteps))
        self.ui.elSizeSlider.setValue(int(self.inputData.mp[2] * 1000))
        self.ui.elSizeMaxLabel.setText(str(self.inputData.mp[2]))
        index = self.ui.unitsComboBox.findText(self.inputData.units, Qt.MatchFixedString)
        if index >= 0:
            self.ui.unitsComboBox.setCurrentIndex(index)
        self.onUnitsChanged(self.inputData.units)

    def updateModel(self):
        """Fill in the model with values from the controls"""
        self.inputData.d = [
            [str(self.ui.d1NameEdit.text()), float(self.checkParse(self.ui.d1Edit.text()))],
            [str(self.ui.d2NameEdit.text()), float(self.checkParse(self.ui.d2Edit.text()))],
            [str(self.ui.d3NameEdit.text()), float(self.checkParse(self.ui.d3Edit.text()))]
            ]
        self.inputData.t = float(self.ui.tEdit.text())
        self.inputData.E = float(self.ui.elasticEdit.text())
        self.inputData.v = float(self.ui.poissonEdit.text())
        self.inputData.bp[1][0] = float(self.ui.bcEdit.text())
        self.inputData.bp[0][0] = str(self.ui.bcNameEdit.text())
        self.inputData.fp[1][0] = float(self.ui.fEdit.text())
        self.inputData.fp[0][0] = str(self.ui.fNameEdit.text())
        self.inputData.fp[2][0] = float(self.ui.fAngleEdit.text())
        self.inputData.c['a'] = float(self.checkParse(self.ui.aEdit.text()))
        self.inputData.c['b'] = float(self.checkParse(self.ui.bEdit.text()))
        self.inputData.c['aStart'] = float(self.checkParse(self.ui.aEdit.text()))
        self.inputData.c['bStart'] = float(self.checkParse(self.ui.bEdit.text()))
        self.inputData.c['aEnd'] = float(self.checkParse(self.ui.aEndEdit.text()))
        self.inputData.c['bEnd'] = float(self.checkParse(self.ui.bEndEdit.text()))
        self.inputData.c['paramA'] = self.ui.paramVaryARadio.isChecked()
        self.inputData.c['paramB'] = self.ui.paramVaryBRadio.isChecked()
        self.inputData.c['aName'] = str(self.ui.aNameEdit.text())
        self.inputData.c['bName'] = str(self.ui.bNameEdit.text())
        self.inputData.paramSteps = int(self.ui.paramStep.value())
        self.inputData.refineMesh = self.ui.elSizeRefineRadio.isChecked()
        self.inputData.mp[2] = float(self.ui.elSizeSlider.value()) / 1000.0
        self.inputData.units = str(self.ui.unitsComboBox.currentText())
        self.addInfoText("Updated Model values.")

    def onActionNew(self):
        """Create a new model"""
        print("onActionNew")
        if self.vis != None:
            self.vis.closeAll()

        # --- Reload default file
        print("Creating a new default model.")
        self.addInfoText("Creating a new default model.")
        self.initModel()
        self.updateControls()

    def onActionOpen(self):
        """Open a new model"""
        print("onActionOpen")
        self.filename, _ = QFileDialog.getOpenFileName(self.ui, "Open model", "", "Model File (*.json *.jpg *.bmp)")
        if self.filename != "":
            self.inputData.load(self.filename)
            self.updateControls()
            self.ui.reportEdit.clear()
            print("Successfully loaded file ", self.filename)
            self.addInfoText("Successfully loaded file" + self.filename)
            self.filename = ""                                  # Reset filename
        else:
            self.ui.reportEdit.clear()
            print("Invalid filename specified. Please try again.")
            self.addInfoText("Invalid filename specified. Please try again.")

    def onActionSave(self):
        """Save the current model"""
        print("onActionSave")
        self.updateModel()

        if self.filename == "":
            self.filename, _ = QFileDialog.getSaveFileName(self.ui, "Save model", "", "Model file (*.json)")
        else:
            try:
                self.inputData.save(self.filename)
            except Exception as e:
                import traceback
                print("Invalid filename. Please try using Save As instead.")
                traceback.print_exc()

        if self.filename != "":
            self.inputData.save(self.filename)
            print("Successfully saved file ", self.filename)
            self.addInfoText("Successfully saved file {0}".format(self.filename))
        else:
            print("Invalid filename specified. Please try again.")
            self.addInfoText("Invalid filename specified. Please try again.")


    def onActionSaveAs(self):
        """Save the current model as a new filename"""
        print("onActionSaveAs")
        self.filename, _ = QFileDialog.getSaveFileName(self.ui, "Save model", "", "Model file (*.json)")
        if self.filename != "":
            self.inputData.save(self.filename)
            print("Successfully saved file ", self.filename)
            self.addInfoText("Successfully saved file {0}".format(self.filename))
        else:
            print("Invalid filename specified. Please try again.")
            self.addInfoText("Invalid filename specified. Please try again.")

    def onActionExit(self):
        """Exit the application"""
        print("onActionExit")
        self.close()

    def onActionLoadDXF(self):
        """Loads a specified DXF file"""
        print("onActionLoadDXF")
        self.inputData.dxf_filename, _ = QFileDialog.getOpenFileName(self.ui,
                                                                "Open dxf drawing", "", "DXF Drawing File (*.dxf)")
        if self.inputData.dxf_filename != "":
            self.inputData.dxf.filename = self.inputData.dxf_filename
            self.inputData.dxf.convertDXFtoSVG()                  # convert it to svg to display in gui
            self.imageGui.updateImage(self.inputData.dxf.drawingname)
            self.imageGui.show()

            self.ui.reportEdit.clear()
            print("Successfully loaded DXF file ", self.inputData.dxf_filename)
            self.addInfoText("Successfully loaded file".format(self.inputData.dxf_filename))
            self.filename = ""  # Reset filename
        else:
            self.ui.reportEdit.clear()
            print("Invalid dxf filename specified. Please try again.")
            self.addInfoText("Invalid dxf filename specified. Please try again.")

    def onActionExecute(self):
        """Execute the current FEM model"""
        print("onActionExecute")

        # --- Disable the interface during the calculation
        self.ui.setEnabled(False)

        # --- Update values from the controller
        self.updateModel()

        # --- If running parameter study, update parameter-related values as well
        if self.runparamstudy:
            self.inputData.paramFilename = "paramStudy"
            self.inputData.paramSteps = int(self.ui.paramStep.value())

        # --- Create a solver
        self.solver = fm.Solver(self.inputData, self.outputData)

        # --- Start a thread to run the calculation so that the interface doesn't freeze
        self.solverThread = SolverThread(self.solver, self.runparamstudy)
        if self.runparamstudy:
            self.solverThread.finished.connect(self.onParamStudyFinished)
        else:
            self.solverThread.finished.connect(self.onSolverFinished)
        self.solverThread.start()

    def onElSizeSliderMoved(self):
        """Update the slider value label"""
        print("onElSizeSliderMoved")
        self.ui.elSizeMaxLabel.setText(str(float(self.ui.elSizeSlider.value()) / 1000.0))

    def onSolverFinished(self):
        """Called when the calculation thread ends"""

        # --- Re-enable the user interface
        self.ui.setEnabled(True)
        self.calcDone = True

        # --- Generate results report
        self.report = fm.Report(self.inputData, self.outputData)
        self.ui.reportEdit.clear()
        self.addInfoText(str(self.report))

        # --- Generate figures
        if(self.calcDone):
            self.vis = fm.Visualization(self.inputData, self.outputData, self.figsOn)
            self.vis.show()
            gc.collect()
            self.vis.wait()

    def onParamStudyFinished(self):
        """Called when parameter study finishes"""

        # --- Re-enable the user interface
        self.ui.setEnabled(True)
        self.calcDone = True

        # --- Generate results report
        self.ui.reportEdit.clear()
        self.addInfoText("Parameter Study successfully finished.")
        self.addInfoText("A total of {0} cases were run.".format(self.outputData.paramnum))
        gc.collect()

    def onoroff(self, val):
        if(val):
            return "ON"
        else:
            return "OFF"

    def onShowGeometry(self):
        """Show geometry window"""
        print("onShowGeometry")
        self.figsOn[0] = 1 - self.figsOn[0]
        self.ui.showGeometryButton.setText("Geometry is " + self.onoroff(self.figsOn[0]))
        self.addInfoText("Geometry Figure set to {0}".format(self.onoroff(self.figsOn[0])))

    def onShowMesh(self):
        """Show mesh window"""
        print("onShowMesh")
        self.figsOn[1] = 1 - self.figsOn[1]
        self.ui.showMeshButton.setText("Mesh is " + self.onoroff(self.figsOn[1]))
        self.addInfoText("Meshing set to {0}".format(self.onoroff(self.figsOn[1])))

    def onShowNodalValues(self):
        """Show nodal values"""
        print("onShowNodalValues")
        self.figsOn[2] = 1 - self.figsOn[2]
        self.ui.showNodalValuesButton.setText("Nodal Values is " + self.onoroff(self.figsOn[2]))
        self.addInfoText("Nodal Values Figure set to {0}".format(self.onoroff(self.figsOn[2])))

    def onShowElementValues(self):
        """Show element values"""
        print("onShowElementValues")
        self.figsOn[3] = 1 - self.figsOn[3]
        self.ui.showElementValuesButton.setText("Element Values is " + self.onoroff(self.figsOn[3]))
        self.addInfoText("Element Values Figure set to {0}".format(self.onoroff(self.figsOn[3])))

    def onExecuteParamStudy(self):
        """Toggle Param Study to be run"""
        self.runparamstudy = 1 - self.runparamstudy
        self.ui.paramButton.setText("Param Study is " + self.onoroff(self.runparamstudy))
        self.addInfoText("Run Parameter Study set to {0}".format(self.onoroff(self.runparamstudy)))

    def onDisplayImage(self):
        """Show the dxf image"""
        self.imageGui.updateImage(self.inputData.dxf.drawingname)
        self.imageGui.show()

    def onScreenClear(self):
        """Clear the UI Screen"""
        self.ui.reportEdit.clear()
        self.ui.reportEdit.insertPlainText("\n>>")

    def onUnitsChanged(self, text):
        """Change the units if user specifies"""
        print("onUnitsChanged")
        self.addInfoText("Units changed to {0}".format(text))
        self.units = text
        if text == 'SI':
            self.ui.tLabel.setText("Thickness (m)")
            self.ui.d1ValueLabel.setText("Value (m)")
            self.ui.d2ValueLabel.setText("Value (m)")
            self.ui.d3ValueLabel.setText("Value (m)")
            self.ui.astartendlabel.setText("start, end (m)")
            self.ui.elasticModLabel.setText("Elastic Modulus (GPa)")
            self.ui.fLabel.setText("Force (N/m)")
            self.ui.bcLabel.setText("Boundary Value (m)")
        elif text == 'IMPERIAL':
            self.ui.tLabel.setText("Thickness (in)")
            self.ui.d1ValueLabel.setText("Value (in)")
            self.ui.d2ValueLabel.setText("Value (in)")
            self.ui.d3ValueLabel.setText("Value (in)")
            self.ui.astartendlabel.setText("start, end (in)")
            self.ui.elasticModLabel.setText("Elastic Modulus (psi)")
            self.ui.fLabel.setText("Force (lbf/in)")
            self.ui.bcLabel.setText("Boundary Value (in)")

    def addInfoText(self, text, user=0):
        """Adds info text"""
        self.ui.reportEdit.moveCursor(QTextCursor.End)
        if user:
            # - If user = 1, this is a user-entered text (i.e.: there is an extra \n at the end)
            self.ui.reportEdit.moveCursor(QTextCursor.Up)
        self.ui.reportEdit.moveCursor(QTextCursor.Up)
        self.ui.reportEdit.moveCursor(QTextCursor.EndOfLine)
        self.ui.reportEdit.moveCursor(QTextCursor.Right, mode=QTextCursor.KeepAnchor)
        self.ui.reportEdit.moveCursor(QTextCursor.Right, mode=QTextCursor.KeepAnchor)
        self.ui.reportEdit.moveCursor(QTextCursor.Right, mode=QTextCursor.KeepAnchor)
        self.ui.reportEdit.cut()
        self.ui.reportEdit.moveCursor(QTextCursor.End)
        self.ui.reportEdit.insertPlainText(text)
        self.ui.reportEdit.insertPlainText("\n\n>>")

    def interpretUserInput(self, user_input):
        """Does actions based on received user input"""
        # - First, only read the final line by finding where the '>' was submitted.
        n = len(user_input)
        start_idx = 0
        for i in range(n):
            if (user_input[n - 1 - i] == '>'):
                start_idx = n - i
                break

        read_line = user_input[start_idx:n]
        if read_line == "\n":
            self.addInfoText('', user=1)
        elif read_line in ghm.menu:
            self.addInfoText(ghm.menu[read_line], user=1)
        else:
            self.addInfoText(
                "\nUnable to read input. Please try again.\n\n", user=1
            )

    def handleKeyRelease(self, event):
        """Handles key inputs to report box"""
        if(event.key() == Qt.Key_Return or event.key() == Qt.Key_Enter):
            self.interpretUserInput(self.reportEdit.toPlainText())



if __name__ == '__main__':

    # --- Create application instance
    app = QApplication(sys.argv)

    # --- Create and view main windows
    widget = MainWindow()
    widget.initModel()
    widget.updateControls()
    widget.show()

    # --- Start the event loop
    sys.exit(app.exec_())           # Returns when all application windows have been closed
