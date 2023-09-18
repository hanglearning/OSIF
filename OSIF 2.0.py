# python 3 tkinter import section

from tkinter import *
import tkinter as Tkinter
import tkinter.filedialog as tkFileDialog
import tkinter.messagebox as tkMessageBox

from datetime import datetime
# end python 3 tkinter import section


# python 2 Tkinter import section
'''
from Tkinter import *
import Tkinter
import tkFileDialog
import tkMessageBox
'''
# end python 2 Tkinter import section


import xlrd
import matplotlib

import itertools

matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec  # Allows for custom positioning of sub plots in matlibplot.
import scipy.optimize

# from default python modules
import os
import re
import webbrowser
import sys

# output pretty title, version info and citation prompt.
# print('\n\n\n#########################################################################')
# print('python version: ' + sys.version)
# print('Tkinter version: ' + str(Tkinter.TkVersion))
# print('#########################################################################\n\n')
title_str = "#####     Open Source Impedance Fitter (OSIF) v2.0 Modified Serialization Version    #####"
print('#' * len(title_str))
print(title_str)
print('#' * len(title_str))
print(
    "--------------------------\nWritten by Jason Pfeilsticker for the Hydrogen fuel cell manufacturing group\nat the National Renewable Energy Lab (NREL). Feb, 2018\nV2.0 adapted Oct. 2021\nCode for additional model options provided by Timothy Van cleve (NREL)")
print(
    "--------------------------\nModified by Hang Chen (chenhanginud@gmail.com) for the The Chemours Company in Newark Delaware during the 2023 Univeristy of Delaware DSI Hackathon.")
# print(
#     '\nthis program uses the matplotlib, scipy, and numpy modules and was written in python.\nIf you publish data from this program, please cite them appropriately.\n'
#     'To cite this code specifically, please cite the OSIF github page at: https://github.com/NREL/OSIF\n--------------------------\n\n\n')

VAR_TO_UNIT = {"Rmem:": "[ohm*cm^2]", "Rcl:": "[ohm*cm^2]", "Qdl:": "[F/(cm^2*sec^phi)]", "Phi:": "[ - ]", "Lwire:": "[H*cm^2]", "Theta:": "[ - ]", "Cell Area:": "[cm^2]", "catalyst loading:": "[mg/cm^2]"}
VAR_TO_INIT = {'Rmem:': '0.03', 'Rcl:': '0.1', 'Qdl:': '2.5', 'Phi:': '0.95', 'Lwire:': '2E-5', 'Theta:': '0.95', 'Cell Area:': '50', 'catalyst loading:': '0.1'}
VAR_TO_NAME = {'Rmem:': 'Rmem', 'Rcl:': 'Rcl', 'Qdl:': 'Qdl', 'Phi:': 'Phi', 'Lwire:': 'Lwire', 'Theta:': 'Theta', 'Cell Area:': 'area', 'catalyst loading:': 'loading'}

COST_TO_MODEL = {} # cost_value: [params, eis_model]
RUN_TIME = None

# Main program class which is called on in line 868-ish to run the program.
class OSIF:

    def __init__(self, master):
        master.title("Open Source Impedance Fitter (OSIF) v2.0 UD")
        master.grid()
        buttonFrame = Frame(master, pady=10, )
        InputFrame = Frame(master, padx=10)
        OutputFrame = Frame(master, padx=10)
        CheckboxFrame = Frame(master, padx=10)

        self.plotFrame = Frame(master, bg='blue')
        self.plotFrameToolBar = Frame(master, bg='red')

        Grid.grid_columnconfigure(buttonFrame, 0, weight=1)
        Grid.grid_rowconfigure(buttonFrame, 0, weight=1)
        Grid.grid_columnconfigure(self.plotFrame, 0, weight=1)
        Grid.grid_rowconfigure(self.plotFrame, 0, weight=1)
        Grid.grid_columnconfigure(self.plotFrameToolBar, 0, weight=1)
        Grid.grid_rowconfigure(self.plotFrameToolBar, 0, weight=1)

        buttonFrame.grid(row=0, columnspan=2)
        InputFrame.grid(row=1, column=0, sticky=N, pady=3)
        OutputFrame.grid(row=1, column=1, sticky=N, pady=3)
        CheckboxFrame.grid(row=1, column=1, sticky=N, pady=3)
        self.plotFrame.grid(row=2, pady=1, padx=8, columnspan=5, sticky=N + S + E + W)
        self.plotFrameToolBar.grid(row=3, pady=1, padx=8, columnspan=5, sticky=S + W)

        self.Rmem = Param()
        self.Rcl = Param()
        self.Qdl = Param()
        self.Phi = Param()
        self.Lwire = Param()
        self.Theta = Param()
        self.area = Param()
        self.frequencyRange = Param()
        self.loading = Param()
        self.currentDataDir = Param()
        self.currentFileName = Tkinter.StringVar(master)
        self.model_selection = Tkinter.StringVar(master)
        self.currentFile = NONE
        self.avgResPer = Param()
        self.activeData = Data()

        entryFont = ("Calibri", '12')
        labelFont = ("Calibri", "12")

        # sdPerColumn = 5
        # sdColumn = 4
        # fitValueColumn = 3        
        # nonFitUnitColumn = 2
        
        varNameColumn = 0
        initValueColumn = 1
        incrementColumn = 2
        upperBoundColumn = 3
        unitColumn = 4
        circuitModelColumn = 5

        Label(InputFrame, text="Initial Values", font=labelFont).grid(row=1, column=initValueColumn, sticky=W)
        Label(InputFrame, text="Increment", font=labelFont).grid(row=1, column=incrementColumn, sticky=W)
        Label(InputFrame, text="Upper Bound", font=labelFont).grid(row=1, column=upperBoundColumn, sticky=W) # stopping condition, may not be included as a variable value
        Label(InputFrame, text="Select Circuit Model:", font=labelFont).grid(row=1, column=circuitModelColumn, sticky=W)
        # Label(OutputFrame, text="Fit Values", font=labelFont).grid(row=1, column=4, sticky=W)
        # Label(OutputFrame, text="Estimated SE", font=labelFont).grid(row=1, column=5, sticky=W)
        # Label(OutputFrame, text="SE % of fit value", font=labelFont).grid(row=1, column=6, sticky=W)
        # save Fit, Es SE and SE % to local

        ################################################
        ############ INPUT INITIAL VALUES ##############
        ################################################

        # create columns and prepopulate values
        row_iter = 2
        for var_name, unit in VAR_TO_UNIT.items():
            Label(InputFrame, text=var_name, font=labelFont).grid(row=row_iter, column=varNameColumn, sticky=E)
            Label(InputFrame, text=unit, font=labelFont).grid(row=row_iter, column=unitColumn, sticky=W)
            
            the_var = vars(self)[VAR_TO_NAME[var_name]]

            the_var.IE = Entry(InputFrame, width=10, font=entryFont)
            the_var.IE.grid(row=row_iter, column=initValueColumn)
            the_var.IE.insert(0, VAR_TO_INIT[var_name])

            the_var.Incre = Entry(InputFrame, width=10, font=entryFont)
            the_var.Incre.grid(row=row_iter, column=incrementColumn)
            the_var.Incre.insert(0, 0)

            the_var.UpperBound = Entry(InputFrame, width=10, font=entryFont)
            the_var.UpperBound.grid(row=row_iter, column=upperBoundColumn)
            the_var.UpperBound.insert(0, VAR_TO_INIT[var_name])

            row_iter += 1

        # deal wtih frequency range

        Label(InputFrame, text="Frequency bound:", font=labelFont).grid(row=row_iter, column=varNameColumn, sticky=E)
        Label(InputFrame, text="[Hz]", font=labelFont).grid(row=10, column=unitColumn, sticky=W)

        self.frequencyRange.OE = Entry(InputFrame, width=10, font=entryFont)
        self.frequencyRange.OE.grid(row=row_iter, column=initValueColumn)
        self.frequencyRange.OE.insert(0, "1")

        self.frequencyRange.IE = Entry(InputFrame, width=10, font=entryFont)
        self.frequencyRange.IE.grid(row=row_iter, column=initValueColumn + 2)
        self.frequencyRange.IE.insert(0, "10000")

        



        ################################################
        #################### BUTTONS ###################
        ################################################

        self.simB = Button(buttonFrame, text="Select Data Directory", command=lambda: self.SelectDataDir())
        self.simB.grid(row=0, column=0, sticky=E)

        self.simB = Button(buttonFrame, text="Model Info.", command=self.openModelInfo)
        self.simB.grid(row=0, column=1, sticky=E)

        self.simB = Button(buttonFrame, text="Citation Info.", command=self.openCitationInfo)
        self.simB.grid(row=0, column=2, sticky=W)

        self.fitB = Button(buttonFrame, text="Fit & Serialize", command=self.serialize_fit)
        self.fitB.grid(row=0, column=3, sticky=E)

        # self.simB = Button(buttonFrame, text="Simulate", command=self.PerformSim)
        # self.simB.grid(row=0, column=4, sticky=W)

        # self.simB = Button(OutputFrame, text="Save Model Data", command=self.SaveData)
        # self.simB.grid(row=9, column=sdPerColumn - 2, columnspan=3, sticky=N)

        Label(buttonFrame, text="Current Directory:", font=labelFont).grid(row=1, column=0, sticky=E)
        self.currentDataDir.IE = Entry(buttonFrame, width=60, font=entryFont)
        self.currentDataDir.IE.grid(row=1, column=1, sticky=W, columnspan=2)
        self.currentDataDir.IE.insert(0, "---")
        self.currentDataDir.IE.config(state='readonly')

        Label(buttonFrame, text="Base Data File:", font=labelFont).grid(row=2, column=0, sticky=E)
        self.choices = ['Data Directory not selected']
        self.fileSelectComboBox = OptionMenu(buttonFrame, self.currentFileName, *self.choices)
        self.fileSelectComboBox.grid(row=2, column=1, sticky=EW, columnspan=2)
        self.fileSelectComboBox.config(font=entryFont)

        # Label(InputFrame, text="Select Circuit Model:", font=labelFont).grid(row=3, column=0, sticky=E)
        self.eis_model = ["Transmission Line", "1-D Linear Diffusion", "1-D Spherical Diffusion"]
        # check box for the eis_models - 0,1,2 corresponds to index of self.eis_model
        for i, m in enumerate(self.eis_model):
                            
            vars(self)[f'model_button_{i}'] = IntVar()
            Checkbutton(CheckboxFrame, text = m, 
                        variable = vars(self)[f'model_button_{i}'],
                        offvalue = 0,
                        height = 2,
                        width = 20).grid(row=i, column=0, sticky=E)
        print()
        # self.model_selection.set(self.eis_model[0])
        # self.fileSelectModelBox = OptionMenu(buttonFrame, self.model_selection, *self.eis_model)
        # self.fileSelectModelBox.grid(row=3, column=1, sticky=EW, columnspan=2)
        # self.fileSelectModelBox.config(font=entryFont)
        # self.fileSelectModelBox.pack()

    def openModelInfo(self):
        if self.model_selection.get() == self.eis_model[0]:
            print(
                "--------------------------\nThe model being used in this program is Eq. 2 from the opened Fuller paper.\nThe derivation can be found in their suplimentary information. If you would\nlike to use a different model, contact Jason Pfeilsticker.\n\n")
            webbrowser.open("https://iopscience.iop.org/article/10.1149/2.0361506jes")
        elif self.model_selection.get() == self.eis_model[1]:
            print("1-D linear diffusion model selected.")
            webbrowser.open("https://www.researchgate.net/publication/342833389_Handbook_of_Electrochemical_Impedance_Spectroscopy_DIFFUSION_IMPEDANCES")
        elif self.model_selection.get() == self.eis_model[2]:
            print("1-D spherical diffusion model selected.")
            webbrowser.open("https://www.researchgate.net/publication/342833389_Handbook_of_Electrochemical_Impedance_Spectroscopy_DIFFUSION_IMPEDANCES")





    def openCitationInfo(self):
        print(
            "--------------------------\nThe Program uses the the matplotlib, scipy, and numpy modules and was written in\n"
            "python (Scientific Computing in Python). Please cite them accordingly via the opend website\n"
            "To cite this code specifically please cite the OSIF github repository")
        webbrowser.open("https://www.scipy.org/citing.html")
        webbrowser.open("https://github.com/NREL/OSIF")

    def SelectDataDir(self):
        newDir = tkFileDialog.askdirectory(title="Select EIS data directory") + '/'
        self.currentDataDir.IE.config(state='normal')
        self.currentDataDir.IE.delete(0, END)
        self.currentDataDir.IE.insert(0, newDir)
        self.currentDataDir.IE.config(state='readonly')

        dirList = os.listdir(newDir)
        dirList = [dataFile for dataFile in dirList if
                   ('.txt' in dataFile) | ('.xls' in dataFile) | ('.xlsx' in dataFile)]
        self.fileSelectComboBox.configure(state='normal')  # Enable drop down
        menu = self.fileSelectComboBox.children['menu']

        # Clear the menu.
        menu.delete(0, 'end')
        for file in dirList:
            # Add menu items.
            menu.add_command(label=file, command=lambda v=self.currentFileName, l=file: v.set(l))
        print('Selected Data Directory: ' + self.currentDataDir.IE.get())

    def LoadSElectedFile(self, pass_in_file=""):
        # if nothing is selected, exit this method
        if not pass_in_file:
            pass_in_file = self.currentFileName.get()

        if (len(pass_in_file) == 0) | (pass_in_file is '---'):
            print('attempt to load on null selection')
            return

        # print(
        #     "\n===========================================\n  loading file: " + pass_in_file + '\n===========================================\n\n')
        # clear the variables of any previous data
        self.activeData.rawFrequency = []
        self.activeData.rawzPrime = []
        self.activeData.rawZdoublePrime = []
        self.activeData.rawzMod = []
        self.activeData.rawZExperimentalComplex = []
        self.activeData.rawmodZExperimentalComplex = []

        self.activeData.dataName = pass_in_file

        if ".xlsx" in self.activeData.dataName:
            self.activeData.dataNameNoExt = self.activeData.dataName[0:len(self.activeData.dataName) - 5]
        else:
            self.activeData.dataNameNoExt = self.activeData.dataName[0:len(self.activeData.dataName) - 4]
        # check for different file formats to parse appropriately.
        if pass_in_file.endswith('.txt'):

            self.activeData.dataNameNoExt = self.activeData.dataName[0:len(self.activeData.dataName) - 4]
            self.currentFile = open(self.currentDataDir.IE.get() + pass_in_file)

            # Run through the lines in the data and parse out the variables into the above lists.
            i = 0  # number of lines in the text file (data and non data)
            dataLineString = self.currentFile.readline()
            freqCol = 0
            zPrimeCol = 1
            zDoublePrimeCol = 2
            zModCol = 3
            negateImZ = 1
            if 'Frequency' in dataLineString:
                freqCol = 1
                zPrimeCol = 2
                zDoublePrimeCol = 3
                zModCol = 4
                if ("-Z''" in dataLineString):
                    print("data loaded in has -Z'' instead of Z''; negating -Z'' column.")
                    negateImZ = -1
            k = 0
            
            while dataLineString:

                # For debuging file input regular expressions.

                regExTest = re.match('^\d*.\d*\t\d*.\d*\t\d*.\d*', dataLineString)
                # skip verbose print
                # if regExTest is not None:
                #     print(
                #         '----------------------------------------------------------------------------------------------------------DATA LINE BELOW')
                # if regExTest is None:
                #     print(
                #         '----------------------------------------------------------------------------------------------------------NOT DATA LINE BELOW\n' + dataLineString)

                if (len(dataLineString) > 2) & (dataLineString[0] != '#') & (
                        re.match('^\d*.\d*\t\d*.\d*\t\d*.\d*', dataLineString) is not None) & (dataLineString != ""):
                    lineList = dataLineString.split("\t")
                    length = len(lineList)
                    last = lineList[length - 1]
                    lineList.remove(last)
                    last = last.strip()
                    lineList.insert(length - 1, last)
                    self.activeData.rawFrequency.append(float(lineList[freqCol]))
                    self.activeData.rawzPrime.append(float(lineList[zPrimeCol]))
                    self.activeData.rawZdoublePrime.append(negateImZ * float(lineList[zDoublePrimeCol]))
                    self.activeData.rawzMod.append(float(lineList[zModCol]))
                    if k == 0:
                        # skip verbose print
                        # print('\n\tFrequency,\t\tRe(Z),\t\t\tIm(Z),\t\t\t\t|Z|')
                        k = 1
                    # skip verbose print
                    # print(lineList[freqCol], lineList[zPrimeCol], lineList[zDoublePrimeCol], lineList[zModCol])
                dataLineString = self.currentFile.readline()
            self.currentFile.close()
            i = 0
            for real in self.activeData.rawzPrime:
                self.activeData.rawZExperimentalComplex.append((real + 1j * self.activeData.rawZdoublePrime[i]))
                self.activeData.rawmodZExperimentalComplex.append(abs(self.activeData.rawZExperimentalComplex[i]))
                i += 1

            # Change into a numpy.array type so least_squares can use it.
            self.activeData.rawFrequency = np.array(self.activeData.rawFrequency)


        # Load in default excel spread sheet output format from EIS software
        elif pass_in_file.endswith('.xlsx') | pass_in_file.endswith('.xls'):

            def checkForNegativeImZReturnImZ(dataRowCol):
                numRows = sheet1.col(0).__len__()
                if dataRowCol[0][3].startswith('-'):
                    i = 1
                    while i < numRows:
                        dataRowCol[i][3] = -1 * float(dataRowCol[i][3])
                        i += 1

                    newHeader = dataRowCol[0]
                    newHeader.remove(dataRowCol[0][3])
                    newHeader.insert(3, dataRowCol[0][3].strip("-"))
                    dataRowCol.remove(dataRowCol[0])
                    dataRowCol.insert(0, newHeader)
                return dataRowCol

            def sheetToListRowCol(sheet1):
                returnList = []
                i = 0
                j = 0
                numRows = sheet1.col(0).__len__()
                numCol = sheet1.row(0).__len__()
                while i < numRows:
                    tempRow = []
                    j = 0
                    while j < numCol:
                        tempRow.append(sheet1.cell(i, j).value)
                        j += 1
                    returnList.append(tempRow)
                    i += 1
                returnList = checkForNegativeImZReturnImZ(returnList)
                return returnList

            def getColDataFromData(dataRowCol, colIndex):
                i = 1
                returnCol = []
                numCol = len(dataRowCol)
                while i < numCol:
                    returnCol.append(dataRowCol[i][colIndex])
                    i += 1
                return returnCol

            xlsx = xlrd.open_workbook(self.currentDataDir.IE.get() + pass_in_file)
            sheet1 = xlsx.sheet_by_index(0)
            data = sheetToListRowCol(sheet1)
            xlsx.release_resources()
            del xlsx

            self.activeData.rawFrequency = getColDataFromData(data, 1)
            self.activeData.rawzPrime = getColDataFromData(data, 2)
            self.activeData.rawZdoublePrime = getColDataFromData(data, 3)
            self.activeData.rawzMod = getColDataFromData(data, 4)
            i = 0
            print('\n\tFrequency,\t\tRe(Z),\t\t\tIm(Z),\t\t\t\t|Z|')
            for real in self.activeData.rawzPrime:
                # create things for graphing later
                self.activeData.rawZExperimentalComplex.append((real + 1j * self.activeData.rawZdoublePrime[i]))
                self.activeData.rawmodZExperimentalComplex.append(abs(self.activeData.rawZExperimentalComplex[i]))

                print(self.activeData.rawFrequency[i], self.activeData.rawzPrime[i], self.activeData.rawZdoublePrime[i],
                      self.activeData.rawzMod[i])
                i += 1

            # Change into a numpy.array type so least_squares can use it.
            self.activeData.rawFrequency = np.array(self.activeData.rawFrequency)

        for i in range(len(self.activeData.rawzPrime)):
            self.activeData.rawPhase.append(
                (180 / (np.pi)) * np.arctan(self.activeData.rawZdoublePrime[i] / self.activeData.rawzPrime[i]))

        print(len(self.activeData.rawPhase), len(self.activeData.rawzMod))
        print("===============================\ndone loading file\n===============================")

    def ChopFreq(self):
        tempFreq = []
        self.activeData.frequency = np.array([])
        for freq in self.activeData.rawFrequency:
            if (freq > float(self.frequencyRange.OE.get())) & (freq < float(self.frequencyRange.IE.get())):
                tempFreq.append(freq)

        # chop the data to the frequency range specified in set up
        self.activeData.frequency = np.array(tempFreq)
        minIndex = self.activeData.rawFrequency.tolist().index(self.activeData.frequency[0])
        maxIndex = self.activeData.rawFrequency.tolist().index(
            self.activeData.frequency[self.activeData.frequency.shape[0] - 1])

        self.activeData.zPrime = self.activeData.rawzPrime[minIndex:maxIndex + 1]
        self.activeData.ZdoublePrime = self.activeData.rawZdoublePrime[minIndex:maxIndex + 1]
        self.activeData.zMod = self.activeData.rawzMod[minIndex:maxIndex + 1]
        self.activeData.modZExperimentalComplex = self.activeData.rawmodZExperimentalComplex[minIndex:maxIndex + 1]
        self.activeData.phase = self.activeData.rawPhase[minIndex:maxIndex + 1]

    # function removed
    # def PerformSim(self):
    #     self.LoadSElectedFile()
    #     if len(self.activeData.rawzPrime) == 0:
    #         tkMessageBox.showinfo("Error!", "No data file loaded\nor data is in incorrect format")
    #         return

    #     else:
    #         self.ChopFreq()
    #         ###### /float(self.area.IE.get())      Rmem
    #         params = [float(self.Lwire.IE.get()) / float(self.area.IE.get()),
    #                   float(self.Rmem.IE.get()) / float(self.area.IE.get()),
    #                   float(self.Rcl.IE.get()) / float(self.area.IE.get()),
    #                   float(self.Qdl.IE.get()),
    #                   float(self.Phi.IE.get()),
    #                   float(self.Theta.IE.get())]

    #         self.CreateFigures(params, 'sim')

    #         print("Model Selection Made:", self.model_selection.get())
    #         simResiduals = self.funcCost(params)
    #         # print simResiduals.get()
    #         self.resPercentData = np.sum(simResiduals / self.activeData.zMod * 100) / len(simResiduals)
    #         self.avgResPer.AVGRESPER.config(state='normal')
    #         self.avgResPer.AVGRESPER.delete(0, END)
    #         self.avgResPer.AVGRESPER.insert(0, '%5.4f' % self.resPercentData)
    #         self.avgResPer.AVGRESPER.config(state='readonly')

    
    def serialize_fit(self):
        global RUN_TIME
        RUN_TIME = datetime.now().strftime("%m%d%Y_%H%M%S")
        # for each variable, get lower, inc and upper bound

        var_to_vals = {}
        fit_count = []

        run_combs = ""
        for var_name in VAR_TO_NAME.values():
            
            the_var = vars(self)[var_name]
            lower_bound = float(the_var.IE.get())
            incre = float(the_var.Incre.get())
            upperbound = float(the_var.UpperBound.get())
            if incre == 0:
                var_to_vals[var_name] = [lower_bound]
            else:
                var_to_vals[var_name] = list(np.round(np.arange(start=lower_bound, stop=upperbound + np.finfo(np.float64).eps * 2, step=incre), 10))
            
            run_combs += f"\n{var_name}: {', '.join([str(i) for i in var_to_vals[var_name]])}"
            fit_count.append(len(var_to_vals[var_name]))

        eis_model_count = 0
        for i, m in enumerate(self.eis_model):    
            if vars(self)[f'model_button_{i}'].get():
                eis_model_count += 1
                run_combs += f"\non {m}"
        fit_count.append(eis_model_count)
        
        total_fits = np.product(fit_count)

        # if base data file not selected
        if not self.currentFileName.get() or self.currentFileName.get() == "Data Directory not selected":
            if tkMessageBox.askokcancel("Confirm", "Please select a base data file to continue."):
                return
        
        while total_fits == 0:
            if tkMessageBox.askokcancel("Confirm", "Please select at least one Circuit Model to continue."):
                return
            
        if tkMessageBox.askokcancel("Confirm", f"Confirm the following run combinations:\n {run_combs}\n\n Run total {total_fits} fits on {self.currentFileName.get()} and serialize?\n If OK, look at terminal output for progress."):
            # iterate over all possible combinations of var values
            for _, values in enumerate(list(itertools.product(*var_to_vals.values()))):
                self.PerformFit(values)

        # get the Prameters and Curcuit model that produces the minimum cost
        min_cost = min(COST_TO_MODEL, key=lambda k: int(k))
        min_params, min_model = COST_TO_MODEL[min_cost]

        print(f"The minimum cost value {min_cost} corresponds to params {min_params} with Circuit Model {min_model}. Applying to the rest of the data files.")

        # get all other files under the selected directory
        file_names = [f for f in os.listdir(self.currentDataDir.IE.get()) if os.path.isfile(self.currentDataDir.IE.get()+'/'+f) and not f.startswith('.')]
        # remove the base file
        file_names.remove(self.currentFileName.get())
        
        # fit to other data files
        for other_file in file_names:
            print(f"\nFitting minimum cost model to {other_file}")
            self.PerformFit(min_params, other_file, True, min_model)

        print("Serialization Done.")
        return

        # float(self.Lwire.IE.get()) / float(self.area.IE.get()),
        #               float(self.Rmem.IE.get()) / float(self.area.IE.get()),
        #               float(self.Rcl.IE.get()) / float(self.area.IE.get()),
        #               float(self.Qdl.IE.get()),
        #               float(self.Phi.IE.get()),
        #               float(self.Theta.IE.get())



    def PerformFit(self, values, pass_in_file="", serialization=False, model=None):

        Rmem, Rcl, Qdl, Phi, Lwire, Theta, area, loading = values
        self.LoadSElectedFile(pass_in_file)
        if not pass_in_file:
            pass_in_file = self.currentFileName.get()
        
        print('\n\n\n\n' + 'Sample: ' + pass_in_file + '\n')
        if len(self.activeData.rawzPrime) == 0:
            tkMessageBox.showinfo("Error!", "No data file loaded\nor data is in incorrect format")
            return

        else:
            ### /float(self.area.IE.get())   Rmem
            self.ChopFreq()
            params = [Lwire / area,
                      Rmem / area,
                      Rcl / area,
                      Qdl,
                      Phi,
                      Theta]

            # Perform the fitting using least_squares with the TRF method and max function calls of 10000.
            # iterate over each circuit model and fit if selected
            model_to_identifier = {0: "", 1: "_l", 2: "_s"}
            for i, m in enumerate(self.eis_model):
                if not vars(self)[f'model_button_{i}'].get():
                    # if model not selected, skip
                    continue
                if serialization and m != model:
                    # if serialization model and m is not the model selected based on cost, skip
                    continue
                

                # if (not serialization and vars(self)[f'model_button_{i}'].get()) or ():
                print(f"Now fitting {values} to {m}")
                finalOutput = scipy.optimize.least_squares(getattr(self, f"funcCost{model_to_identifier[i]}"), params,
                                                    bounds=[(0, 0, 0, 0, 0, 0), (1, np.inf, np.inf, np.inf, 1, 1)],
                                                    max_nfev=50000, method='trf', xtol=1e-11,
                                                    ftol=1e-11, gtol=1e-11, verbose=1)
                self.finalParams = finalOutput.x

                # Estimate variance of parameters based on Gauss-Newton approximation of the Hessian of the cost function. See: (https://www8.cs.umu.se/kurser/5DA001/HT07/lectures/lsq-handouts.pdf)
                # basically Covariance matrix = inverse(Jacob^T*Jacob)*meanSquaredError, where Jacob^T*Jacob is the first order estimate for the hessian. The square root of the diagonal elements (c_ii) of Cov are the variances of the parameter b_i

                # sigma squared estimate also called s^2 sometimes = chi^2 reduced
                sigmaSquared = (np.matmul(np.transpose(finalOutput.fun), finalOutput.fun)) / (
                        finalOutput.fun.shape[0] - self.finalParams.size)

                # Plot residuals for error checking. Note: cancels update of the normal plots in the normal UI.
                # plt.figure(100)
                # plt.plot(finalOutput.fun)
                # plt.show()

                Jacob = finalOutput.jac
                estVars = np.matrix.diagonal(
                    np.linalg.inv(Jacob.T.dot(Jacob)) * (np.matmul(np.transpose(finalOutput.fun), finalOutput.fun)) / (
                            finalOutput.fun.shape[0] - self.finalParams.size))
                # estVars = sigmaSquared*np.matrix.diagonal(np.linalg.inv(np.matmul(np.matrix.transpose(finalOutput.jac),finalOutput.jac)))

                # estimated errors in parameters
                self.standardDeviation = np.sqrt(estVars)

                # taking chi^2  = sum((residuals^T)*(residuals)) = sum (r_i)^2
                # print('Sum res^2 = ' + str(np.sum((np.matrix.transpose(finalOutput.fun).__mul__(finalOutput.fun)))))

                self.L2NormOfRes = np.sqrt(np.sum(pow(finalOutput.fun, 2)))

                # print('\nNormalized grad = grad/|grad| = ' + str(finalOutput.grad / np.sqrt(np.sum(pow(finalOutput.grad, 2)))))

                self.resPercentData = np.sum(finalOutput.fun / self.activeData.zMod * 100) / finalOutput.fun.shape[0]

                # print('\nFit to: ' + self.activeData.dataNameNoExt)

                self.percentSigma = self.standardDeviation / self.finalParams * 100

                self.fitOutPutString = '#\n#\n#\t\t\t\t\t\t\t   Fit values\t\t\t~std Error\t\t\t ~std Error %% of value\n#\n#\tRmem  [ohm*cm^2] \t\t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tRcl   [ohm*cm^2] \t\t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f' \
                                    '\n#\tQdl   [F/(sec^phi)]  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tphi   [ ]  \t\t\t\t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tLwire [H*cm^2] \t\t\t  = %.4e\t\t\t%.3e\t\t\t\t%8.2f' \
                                    '\n#\tphi   [ ]  \t\t\t\t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\n#\tQdl/mgpt = %5.6f\n#\tL2 norm of res = %10.8f [ohm*cm^2]' % \
                                    (float(self.finalParams[1]) * float(self.area.IE.get()),
                                        float(self.standardDeviation[1]) * float(self.area.IE.get()), self.percentSigma[1],
                                        float(self.finalParams[2]) * float(self.area.IE.get()),
                                        float(self.standardDeviation[2]) * float(self.area.IE.get()), self.percentSigma[2],
                                        float(self.finalParams[3]), float(self.standardDeviation[3]), self.percentSigma[3],
                                        float(self.finalParams[4]), float(self.standardDeviation[4]), self.percentSigma[4],
                                        float(self.finalParams[0]) * float(self.area.IE.get()),
                                        float(self.standardDeviation[0]) * float(self.area.IE.get()), self.percentSigma[0],
                                        float(self.finalParams[5]), float(self.standardDeviation[5]), self.percentSigma[5],
                                        float(self.finalParams[3]) / (
                                                    float(self.area.IE.get()) * float(self.loading.IE.get())),
                                        float(self.L2NormOfRes) * float(self.area.IE.get()))

                # print(self.fitOutPutString)

                self.realFinalModel = getattr(self, f"funcreal{model_to_identifier[i]}")(self.finalParams)
                self.imagFinalModel = getattr(self, f"funcImg{model_to_identifier[i]}")(self.finalParams)
                self.zModFinalModel = getattr(self, f"funcAbs{model_to_identifier[i]}")(self.finalParams)
                self.phaseFinalModel = getattr(self, f"funcPhase{model_to_identifier[i]}")(self.finalParams)

                # if self.model_selection.get() == self.eis_model[0]:
                #     self.realFinalModel = self.funcreal(self.finalParams)
                #     self.imagFinalModel = self.funcImg(self.finalParams)
                # elif self.model_selection.get() == self.eis_model[1]:
                #     self.realFinalModel = self.funcreal_l(self.finalParams)
                #     self.imagFinalModel = self.funcImg_l(self.finalParams)
                # elif self.model_selection.get() == self.eis_model[2]:
                #     self.realFinalModel = self.funcreal_s(self.finalParams)
                #     self.imagFinalModel = self.funcImg_s(self.finalParams)
                # self.zModFinalModel = self.funcAbs(self.finalParams)
                # self.phaseFinalModel = self.funcPhase(self.finalParams)
                # (180/(np.pi))*np.arctan(self.funcImg(self.finalParams)/self.funcreal(self.finalParams))

                self.AvgRealResPer = np.sum(abs(
                    np.array(abs(self.realFinalModel - self.activeData.zPrime)) / np.array(self.activeData.zPrime)) * 100) / \
                                    self.realFinalModel.shape[0]
                self.AvgImagResPer = np.sum(abs(
                    np.array(abs(self.imagFinalModel - self.activeData.ZdoublePrime)) / np.array(
                        self.activeData.ZdoublePrime)) * 100) / \
                                    self.imagFinalModel.shape[0]

                self.CreateFigures(self.finalParams, 'fit', model_to_identifier[i], m, values, serialization, pass_in_file, m)
        
        
        # print(COST_TO_MODEL)

    def CreateFigures(self, params, fitOrSim, model_identifier, eis_model_name, orig_params, serialization, pass_in_file, m):
        if fitOrSim == 'fit':
            graphLabel = 'Full complex fit: '
        elif fitOrSim == 'sim':
            graphLabel = 'Simulated using: '
        else:
            graphLabel = ''

        plt.close('all')
        # make layout for graphs
        gs0 = gridspec.GridSpec(1, 2)
        gs00 = gridspec.GridSpecFromSubplotSpec(4, 3, subplot_spec=gs0[0])
        gs01 = gridspec.GridSpecFromSubplotSpec(4, 4, subplot_spec=gs0[1])
        f = plt.figure(1, figsize=[8, 3.5], tight_layout='true')

        ###########################################################
        ####          PLOT NYQUIST COMBINED FITTINGS           ####
        ###########################################################

        nyGraph = plt.Subplot(f, gs01[:, :])
        f.add_subplot(nyGraph)
        nyGraph.plot(self.activeData.zPrime, self.activeData.ZdoublePrime, 'bo', ls='--', markersize=2, linewidth=1,
                     label='data: ' + self.activeData.dataNameNoExt)

        x_axis = getattr(self, f"funcreal{model_identifier}")(params)
        y_axis = getattr(self, f"funcImg{model_identifier}")(params)
        
        nyGraph.plot(x_axis, y_axis, 'ro', markersize=2,
                        label='\n%s\nLwire=%.5e\nRmem=%5.8f\nRcl=%5.8f\nQdl=%5.5f\nphi=%5.5f\ntheta=%5.5f' % (
                            graphLabel, params[0] * float(self.area.IE.get()), params[1] * float(self.area.IE.get()),
                            params[2] * float(self.area.IE.get()), params[3], params[4], params[5]))

        plt.gca().invert_yaxis()
        plt.xticks(rotation=20)

        plt.xlabel('Re(Z)')
        plt.ylabel('Im(Z)')
        plt.legend(loc=2, fontsize=6)

        ###########################################################
        ####        PLOT Phase vs w COMBINED FITTINGS          ####
        ###########################################################

        phaseGraph = plt.Subplot(f, gs00[-4, :3])
        f.add_subplot(phaseGraph)
        phaseGraph.plot(self.activeData.frequency, self.activeData.phase, 'bo', ls='--', markersize=2,
                        linewidth=1)
        y_axis = getattr(self, f"funcPhase{model_identifier}")(params)
        phaseGraph.plot(self.activeData.frequency, y_axis, 'ro', markersize=2)
        plt.ylabel('phase')
        # plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')

        ###########################################################
        ####         PLOT |Z| vs w COMBINED FITTINGS           ####
        ###########################################################

        modZgraph = plt.Subplot(f, gs00[-3, :3])
        f.add_subplot(modZgraph)
        modZgraph.plot(self.activeData.frequency, self.activeData.modZExperimentalComplex, 'bo', ls='--', markersize=2,
                       linewidth=1)
        y_axis = getattr(self, f"funcAbs{model_identifier}")(params)
        modZgraph.plot(self.activeData.frequency, y_axis, 'ro', markersize=2)
        plt.ylabel('|Z|')
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')

        ###########################################################
        ####                   PLOT Im(Z)                      ####
        ###########################################################

        imZgraph = plt.Subplot(f, gs00[-2, :3])
        f.add_subplot(imZgraph)
        imZgraph.plot(self.activeData.frequency, self.activeData.ZdoublePrime, 'bo', ls='--', markersize=2, linewidth=1)

        y_axis = getattr(self, f"funcImg{model_identifier}")(params)
        imZgraph.plot(self.activeData.frequency, (y_axis), 'ro', markersize=2)

        plt.ylabel('Im(Z)')
        plt.gca().set_yscale('linear')
        plt.gca().set_xscale('log')
        plt.gca().set_xticks([])

        ###########################################################
        ####                   PLOT Re(Z)                      ####
        ###########################################################

        reZgraph = plt.Subplot(f, gs00[-1, :3])
        f.add_subplot(reZgraph)
        reZgraph.plot(self.activeData.frequency, self.activeData.zPrime, 'bo', ls='--', markersize=2, linewidth=1)

        y_axis = getattr(self, f"funcreal{model_identifier}")(params)
        reZgraph.plot(self.activeData.frequency, y_axis, 'ro', markersize=2)

        plt.xlabel('frequency')
        plt.ylabel('Re(Z)')
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')

        ###########################################################
        ####              Draw figure in Tkinter               ####
        ###########################################################
        # for widget in self.plotFrame.winfo_children():
        #     widget.destroy()
        # for widget in self.plotFrameToolBar.winfo_children():
        #     widget.destroy()

        # dataPlot = FigureCanvasTkAgg(f, master=self.plotFrame)
        # dataPlot.close_event()
        # dataPlot.draw()
        # dataPlot.get_tk_widget().grid(row=0, sticky=N + S + E + W, )
        # # toolbar = NavigationToolbar2TkAgg(dataPlot, self.plotFrameToolBar)
        # # toolbar.update()
        # dataPlot._tkcanvas.grid(row=0, sticky=W + S)

        # save plot
        import secrets
        Rmem, Rcl, Qdl, Phi, Lwire, Theta, area, loading = orig_params
        plt_dir = f"{self.currentDataDir.IE.get()}plots/{RUN_TIME}"
        if serialization:
            plt_dir = f"{plt_dir}/serialized/{eis_model_name}"
        else:
            plt_dir = f"{plt_dir}/initial_fit/{eis_model_name}"
        os.makedirs(plt_dir, exist_ok=True)
        with open(f"{plt_dir}/plot_details_{pass_in_file}.txt", "w+") as f:
            f.write(f"{pass_in_file}\n")
            f.write(f"{m}\n")
            f.write(f"Original params - Rmem: {Rmem}, Rcl: {Rcl}, Qdl: {Qdl}, Phi: {Phi}, Lwire: {Lwire}, Theta: {Theta}, area: {area}, loading: {loading}\n")
            f.write(f"Transformed params {params}\n")
            # calculate cost function
            cost = sum(getattr(self, f"funcCost{model_identifier}")(params))
            COST_TO_MODEL[cost] = [orig_params, eis_model_name]
            f.write(f"Cost value {params}: {cost}\n")
        plt.savefig(f'{plt_dir}/{pass_in_file}.png')
                    
        print('done with plotting')

    def KILLALL(self):
        for widget in self.plotFrame.winfo_children():
            widget.destroy()
        #   for widget in self.plotFrameToolBar.winfo_children():
        #       widget.destroy()
        print("\n\nAll plots killed. \nHave a nice day!")

    def SaveData(self):
        if (len(self.currentFileName.get()) == 0) | (self.currentFileName.get() is '---'):
            print('no data loaded')
            tkMessageBox.showinfo("Error!", "No data file loaded")
            # nothing selected to load
            return

        dataOutFile = open(self.currentDataDir.IE.get() + self.activeData.dataNameNoExt + '_fit.txt', "w+")
        i = 0
        dataOutFile.write('#Fitted model at fitting frequencies:\n#Frequency\t\tRe(Z)\t\t\tIm(Z)\t\t\t|Z|\n')
        for real in self.realFinalModel:
            dataOutFile.write(
                str(self.activeData.frequency[i]) + '\t' + str(self.realFinalModel[i]) + '\t' + str(
                    self.imagFinalModel[i]) + '\t' + str(
                    self.zModFinalModel[i]) + '\n')
            i += 1

        dataOutFile.write(
            '#\n#\n#\t\t\t\t   Fit values\t\t\t~std dev\t\t\t   ~stdDev %% of value\n#\n#\tRmem  [ohm*cm^2] \t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tRcl   [ohm*cm^2] \t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f' \
            '\n#\tQdl   [F/(cm^2*sec^phi)]  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tphi   [ ]  \t\t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\tLwire [H*cm^2] \t\t  = %.4e\t\t\t%.3e\t\t\t\t%8.2f' \
            '\n#\ttheta   [ ]  \t\t  = %5.8f\t\t\t%.3e\t\t\t\t%8.2f\n#\n#\tQdl/mgpt = %5.6f\n#\tL2 norm of res = %10.8f' % (
                float(self.finalParams[1]) * float(self.area.IE.get()),
                float(self.standardDeviation[1]) * float(self.area.IE.get()), self.percentSigma[1],
                float(self.finalParams[2]) * float(self.area.IE.get()),
                float(self.standardDeviation[2]) * float(self.area.IE.get()), self.percentSigma[2],
                float(self.finalParams[3]), float(self.standardDeviation[3]), self.percentSigma[3],
                float(self.finalParams[4]), float(self.standardDeviation[4]), self.percentSigma[4],
                float(self.finalParams[0]) * float(self.area.IE.get()),
                float(self.standardDeviation[0]) * float(self.area.IE.get()), self.percentSigma[0],
                float(self.finalParams[5]), float(self.standardDeviation[5]), self.percentSigma[5],
                float(self.finalParams[3]) / (float(self.area.IE.get()) * float(self.loading.IE.get())),
                float(self.L2NormOfRes)))

        dataOutFile.write('\n#\tAvg. |Z| residual % WRT to data |Z| = ' + str(self.resPercentData))

        dataOutFile.close()
        print("saved data in: " + self.currentDataDir.IE.get() + self.activeData.dataNameNoExt + '_fit.txt')

    # Because built in coth(x) cant deal with complex numbers because exp(x) cant deal with them, but pow(x,y) can.
    def JPcoth(self, x):
        return (pow(np.e, x) + pow(np.e, -x)) / (pow(np.e, x) - pow(np.e, -x))
    
    '''  - Transmission Line
        l - 1-D Linear Diffusion
        s - 1-D Spherical Diffusion
    '''

    # Minimizing this function results in fitting the real and complex parts of the impedance at the same time.
    def funcCost(self, params):
        return np.array(np.sqrt(pow((self.funcreal(params) - self.activeData.zPrime), 2) + pow(
            (self.funcImg(params) - self.activeData.ZdoublePrime), 2)))
    def funcCost_l(self, params):
        return np.array(np.sqrt(pow((self.funcreal_l(params) - self.activeData.zPrime), 2) + pow(
                (self.funcImg_l(params) - self.activeData.ZdoublePrime), 2)))
    def funcCost_s(self, params):
        return np.array(np.sqrt(pow((self.funcreal_s(params) - self.activeData.zPrime), 2) + pow(
                (self.funcImg_s(params) - self.activeData.ZdoublePrime), 2)))

    # Define the functions of the model (equation 2 from "A Physics-Based Impedance Model of Proton Exchange Membrane Fuel
    # Cells Exhibiting Low-Frequency Inductive Loops")

    def funcAbs(self, param):
        return abs(self.funcreal(param) + 1j * self.funcImg(param))
    def funcAbs_l(self, param):
        return abs(self.funcreal_l(param) + 1j * self.funcImg_l(param))
    def funcAbs_s(self, param):
        return abs(self.funcreal_s(param) + 1j * self.funcImg_s(param))


    def funcPhase(self, param):
        return 180 / np.pi * np.arctan(self.funcImg(param) / self.funcreal(param))
    def funcPhase_l(self, param):
        return 180 / np.pi * np.arctan(self.funcImg_l(param) / self.funcreal_l(param))
    def funcPhase_s(self, param):
        return 180 / np.pi * np.arctan(self.funcImg_s(param) / self.funcreal_s(param))

    #Function from Setzler
    def funcreal(self, param):
        return np.real(param[0] * pow((1j * 2 * np.pi * self.activeData.frequency), param[5]) + param[1] + pow(
            (param[2] / (param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4]))), 0.5) * self.JPcoth(
            pow((param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4])), 0.5)))

    def funcImg(self, param):
        return np.imag(param[0] * pow((1j * 2 * np.pi * self.activeData.frequency), param[5]) + param[1] + pow(
            (param[2] / (param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4]))), 0.5) * self.JPcoth(
            pow((param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4])), 0.5)))

    # 1-D linear diffusion model
    def funcreal_l(self, param):
        return np.real(
            param[0] * pow((1j * 2 * np.pi * self.activeData.frequency), param[5]) + param[1] + param[2] * pow(
                (param[2] * (param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4]))),
                -0.5) * self.JPcoth(
                pow((param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4])), 0.5)))

    def funcImg_l(self, param):
        return np.imag(
            param[0] * pow((1j * 2 * np.pi * self.activeData.frequency), param[5]) + param[1] + param[2] * pow(
                (param[2] * (param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4]))),
                -0.5) * self.JPcoth(
                pow((param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4])), 0.5)))

    # 1-D spherical diffusion model
    def funcreal_s(self, param):
        return np.real(param[0] * pow((1j * 2 * np.pi * self.activeData.frequency), param[5]) + param[1] + param[2] /
                       (pow(param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4]), 0.5) *
                        self.JPcoth(
                            pow((param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4])),
                                0.5)) - 1))

    def funcImg_s(self, param):
        return np.imag(param[0] * pow((1j * 2 * np.pi * self.activeData.frequency), param[5]) + param[1] + param[2] /
                       (pow(param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4]), 0.5) *
                        self.JPcoth(
                            pow((param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4])),
                                0.5)) - 1))

        #   pow((-1 + (param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4]), 0.5) *
        #   self.JPcoth( pow((param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4])), 0.5))),
        #   -1))

        # return np.imag(param[0] * pow((1j * 2 * np.pi * self.activeData.frequency), param[5]) + param[1] + param[2] *
        #    pow((-1 + ((param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4]), 0.5)
        #    *self.JPcoth(pow((param[2] * param[3] * pow((1j * 2 * np.pi * self.activeData.frequency), param[4])), 0.5)))),
        #    -1))


class Param():

    def __init__(self):
        self.IE = Entry()
        self.OE = Entry()
        self.OESD = Entry()
        self.OESDP = Entry()
        self.AVGRESPER = Entry()

        self.Incre = Entry()
        self.UpperBound = Entry()


class Data():

    def __init__(self):
        self.dataName = ''
        self.dataNameNoExt = ''

        self.zPrime = []
        self.ZdoublePrime = []
        self.zMod = []
        self.modZExperimentalComplex = []
        self.frequency = np.array([])
        self.phase = []

        self.rawzPrime = []
        self.rawZdoublePrime = []
        self.rawzMod = []
        self.rawmodZExperimentalComplex = []
        self.rawFrequency = []
        self.rawPhase = []


def on_closing():
    if tkMessageBox.askokcancel("Quit", "Do you want to quit?"):
        app.KILLALL()
        root.destroy()
        os._exit(0)


root = Tk()
app = OSIF(root)
root.protocol("WM_DELETE_WINDOW", on_closing)
root.mainloop()