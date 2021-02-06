"""  
University of Manchester, Physics and Astronomy Department
Year 2, Semester 1, Module: Introduction to Programming for Physicists.

                        Extraction of Lambda coefficient and Half life for Rubidium and Strontium
                        from detector data 
                        
Author: Andreas Mastronikolis
Student ID: 10281135
Date: 13/12/2019

Purpose of the program:

This application calculates the lambda coefficient and half life of Rubidium and Strontium
given data from two different detectors. It attempts to analyze the data by plotting a contour
plot of the chi squared distribution and the best fit that describes the said data. Throughout
code, the lambda factor is denoted with the greek letter lambda while half time with the greek
letter tau.

"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as pc
import matplotlib.ticker as mtick
from matplotlib import gridspec
from scipy.optimize import fmin

Data = np.array([])
Data_1 = np.array([])
Data_2 = np.array([])
data_open = False # When this variable is true, the data have been read by the application. Otherwise, the data haven't been read.
# Initiation of Parameters
Chi_Squared_Final = 0 
Reduced_Chi_Squared = 0
Parameters = 2
Deg_Freedom = 0


def Activity_Rb(Time, lambdaf_Sr, lambdaf_Rb):
    """
    Below is the function that measures the activity of Rubidium as a
    function of time
    """
    Avogadros_No = pc.N_A
    Initial_Population = 10**(-6) * Avogadros_No
    Population_Rb_lamdaf = lambdaf_Sr / (lambdaf_Rb - lambdaf_Sr)
    Population_Rb_Exps = np.exp(-1 * lambdaf_Sr * Time) - np.exp(-1 * lambdaf_Rb * Time)
    Population_Rb = (Initial_Population) * (Population_Rb_lamdaf) * (Population_Rb_Exps)
    Activity_Rb = lambdaf_Rb * Population_Rb
    return Activity_Rb

def Chi_Squared(lambda_Vector):
    """

    """
    lambdaf_Sr = lambda_Vector[0]
    lambdaf_Rb = lambda_Vector[1]
    Output = 0
    for i in Data:
        Output += ( (i[1] - Activity_Rb(i[0], lambdaf_Sr, lambdaf_Rb))**2 / i[2]**2 )
    return Output

def Validation_Function(File, Row, Column):
    """
    This function returns True only when the value under
    consideration is not an acceptable value. 
    In order to be deemed such, it should be a strictly 
    positive numerical value, that doesn't deviate more 
    than three standard deviations from the appropriate 
    data set.    
    """
    try:

       if np.isfinite(File[Row, Column]) == False or File[Row, Column] < 0:
           return True
       elif Column == 0:
           if File[Row, Column] > File[Row + 1, Column]:
               return True
           else:
               return False
       elif (Row == 0 and Column == 1):
           Local_Region = np.array([File[Row + 1, Column], File[Row + 2, Column], File[Row + 3, Column]])
           MEAN = np.mean(Local_Region)
           STD = np.std(Local_Region)
           if np.abs(File[Row, Column] - MEAN) > 3*STD:
               return True
           else:
               return False           
       elif (Row == -1 and Column == 1):
           Local_Region = np.array([File[Row - 1, Column], File[Row - 2, Column], File[Row - 3, Column]])
           STD = np.std(Local_Region)
           MEAN = np.mean(Local_Region)
           if np.abs(File[Row, Column] - MEAN) > 3*STD:
               return True
           else:
               return False
       elif Column == 1:
           NDiff = np.abs(File[Row + 1, Column] - File[Row, Column])
           PDiff = np.abs(File[Row - 1, Column] - File[Row, Column])
           if NDiff > 100 and PDiff > 100:
               return True
           else:
               return False
       elif Column == 2:
           
           if File[Row, Column] == 0:
               return True
           else:
               return False     
        
    except:
        return False

"""
The next block of code opens the data files and imports the
raw detector data into the application.

"""
while data_open == False:
    
    Input_No = 1
    while Input_No == 1:
        try:
            Detector_1_Data = str(input('Write the name of the file that contains data from the first detector: '))
            Input_No += 1
        except:
            print('This is not a string')

    while Input_No == 2:
        try:
            Detector_2_Data = str(input('Write the name of the file that contains data from the second detector: '))
            Input_No += 1
        except:
            print('This is not a string')
    while Input_No == 3:
        try:
            File_Type = str(input('Write the type of the two files: '))
            Input_No += 1
            print('Detector 1 data are taken from the file:', Detector_1_Data + '.' + File_Type, 'and detector 2 data are taken from the file:', Detector_2_Data + '.' + File_Type)
        except:
            print('Write a correct file extension')
    try:
        Data_1 = np.genfromtxt(Detector_1_Data + '.' + File_Type, delimiter = ',', skip_header = 1)
        Data_2 = np.genfromtxt(Detector_2_Data + '.' + File_Type, delimiter = ',', skip_header = 1)
        print('The data files have been identified')
        data_open = True
    except:
        print('\nThere seems to be a problem with importing the data. Make sure the data file is in the same directory as the script and you have copied the name correctly')

if data_open:
    Data_1 = np.genfromtxt(Detector_1_Data + '.' + File_Type, delimiter = ',', skip_header = 1)
    Data_2 = np.genfromtxt(Detector_2_Data + '.' + File_Type, delimiter = ',', skip_header = 1)

Data_1_ElementNo = np.prod(Data_1.shape) # This the number of all numerical values accessible for analysis
Data_1_RowNo = len(Data_1) # Number of Columns and Rows in each data file
Data_1_ColNo = len(Data_1[0,:])
Data_2_ElementNo = np.prod(Data_2.shape)
Data_2_RowNo = len(Data_2)
Data_2_ColNo = len(Data_2[0,:])

"""
The next block of code attemps to validate the two imported data files and
concatenates the two files into a single one, which is sorted.
"""
print('_' * 100)
print('\nValidating the imported data files...\n')



""" The following two arrays will collect the rows and
columns of faulty data from detector 1.
"""
Faulty_Row_Number_1 = np.array([])
Faulty_Column_Number_1 = np.array([])

for a in np.arange(Data_1_RowNo):
    for i in np.arange(Data_1_ColNo):
        if Validation_Function(Data_1, a, i) == True:
            Faulty_Row_Number_1 = np.append(Faulty_Row_Number_1, a)
            Faulty_Column_Number_1 = np.append(Faulty_Column_Number_1, i)

""" The following two arrays will collect the rows and
columns of faulty data from detector 2.
"""

Faulty_Row_Number_2 = np.array([])
Faulty_Column_Number_2 = np.array([])

for a in np.arange(Data_2_RowNo):
    for i in np.arange(Data_2_ColNo):
        if Validation_Function(Data_2, a, i) == True:
            Faulty_Row_Number_2 = np.append(Faulty_Row_Number_2, a)
            Faulty_Column_Number_2 = np.append(Faulty_Column_Number_2, i)

Faulty_Data = np.vstack((Data_1, Data_2)) # This array contains every data point from the detectors
Data_1 = np.delete(Data_1, Faulty_Row_Number_1.astype(int), 0) # The rows of elements that don't pass the validation test get deleted
Data_2 = np.delete(Data_2, Faulty_Row_Number_2.astype(int), 0)
Data = np.vstack((Data_1, Data_2)) # This array contains all acceptable data values from both detectors
Data = Data[np.argsort(Data[:,0])] # This line sorts the data in order of increasing time

"""
The following lines convert the units of the in Bq when it comes to
activity and seconds when it comes to time.
"""

Faulty_Data[:,0] = Faulty_Data[:,0] * 3600
Faulty_Data[:,1] = Faulty_Data[:,1] * 10**(12)
Faulty_Data[:,2] = Faulty_Data[:,2] * 10**(12)

Data[:,0] = Data[:,0] * 3600
Data[:,1] = Data[:,1] * 10**(12)
Data[:,2] = Data[:,2] * 10**(12)
Data_Points = len(Data[:,0])
Deg_Freedom = Data_Points - Parameters

"""
Output messages after validating the two files
"""


if len(Faulty_Row_Number_1) == 0 and len(Faulty_Row_Number_2) == 0:
    print('Validation complete! No faulty values were detected on both of the files.')

else:
    print('Validation complete! Both files were processed and problematic values were identified. Below are those values expressed in tuples (a,b), where a represents the rows and b the columns inside the specified file. \n' +
          '_' * 100 + '\n\nFor the file: "' + Detector_1_Data + '.' + File_Type + '": \n')
    
    
    if len(Faulty_Row_Number_1) != 0:
        
        for i in np.arange(0, len(Faulty_Row_Number_1)):
            
            print('{}.'.format(i + 1),'({:g},{:g})'.format(Faulty_Row_Number_1[i] + 1, Faulty_Column_Number_1[i] + 1))
    else:
        print('No faulty points were spotted in this file')        
    
    print('\nFor the file "' + Detector_2_Data + '.' + File_Type + '": \n')
    
    if len(Faulty_Row_Number_2 != 0):
        
        for i in np.arange(0, len(Faulty_Row_Number_2)):
            
            print('{}.'.format(i + 1),'({:g},{:g})'.format(Faulty_Row_Number_2[i] + 1, Faulty_Column_Number_2[i] + 1))
    else:
        print('No faulty points were spotted in this file')
    
    print('\nThus, the number of the data values from for further analysis are {} from "Nuclear_Data_1.csv" \nand {} from "Nuclear_Data_2.csv".'.format(Data_1_ElementNo - len(Faulty_Row_Number_1), Data_2_ElementNo - len(Faulty_Row_Number_2) ))

print('_' * 100)
print('\nAnalyzing the imported data...')


"""
The next block of code executes the minimization of the chi squared distribution

"""

Sr_Interval = np.linspace(0.001, 0.015, 500) # Definition of a sensible ranges in the lambda factors in order to minimize chi squared
Rb_Interval = np.linspace(0.00001, 0.0009, 500)
xx, yy = np.meshgrid(Sr_Interval, Rb_Interval) # We generate a grid with vectors to which we will look for minimal values of chi squared

Fit_Results = fmin(Chi_Squared, [0.005,0.0005], disp = False, full_output = True) # This line executes the minimization of the chi squared distribution
Chi_Squared_Min = Fit_Results[1] # The minimum chi squared and the vector where this is achieved are defined in the two lines below
[Sr_min, Rb_min] = Fit_Results[0]
Chi_Squared_Levels = (Chi_Squared_Min + 1.00, Chi_Squared_Min + 2.30, Chi_Squared_Min + 5.99,
                      Chi_Squared_Min + 9.21) # Definition of confidence levels in chi - squared
hL_Sr = np.log(2) / Sr_min # Definition of the half life of an element
hL_Rb = np.log(2) / Rb_min

New_Sr_Interval = np.linspace(Sr_min - 0.1*Sr_min, Sr_min + 0.1*Sr_min, 500)
New_Rb_Interval = np.linspace(Rb_min - 0.05*Rb_min, Rb_min + 0.05*Rb_min, 500)

N_xx, N_yy = np.meshgrid(New_Sr_Interval, New_Rb_Interval) # This grid is defined such that the contour plot of the chi-squared is closely around its minimal value



""" The following block plots contours of the 
reduced chi squared distribution along with
other important information. """

Parameter_Space = plt.figure('Parameter Space', figsize = (8, 4.8))
Parameter_Space.canvas.manager.window.move(1000,10)
fig = Parameter_Space.add_subplot(111, label = 'PS')
fig.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
fig.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
fig.set_title(r'Contours of $\chi^2$ distribution')
fig.set_ylabel(r'$\lambda$ factor of $^{79}$Rb [$s^{-1}$]')
fig.set_xlabel(r'$\lambda$ factor of $^{79}$Sr [$s^{-1}$]')
color_Filled_Plot = fig.contourf(N_xx, N_yy, Chi_Squared([N_xx,N_yy]), 100, cmap=plt.cm.bone)
pure_Contour = fig.contour(N_xx, N_yy, Chi_Squared([N_xx,N_yy]), levels = Chi_Squared_Levels, colors = 'w')
fig.clabel(pure_Contour, inline = 1, fontsize = 7, colors = 'w')
st_SDeviation = pure_Contour.allsegs[0][0]
fig.scatter(Sr_min, Rb_min, marker = 'x', color = 'w', label = 'minimum')
Parameter_Space.savefig('Parameter Space.png', dpi = 300)

"""
We extract the uncertainties from the results obtained above. The word
'Sigma' in the following block denotes the uncertainty in the value that 
follows it
"""

Maximum_Sr = max(st_SDeviation[:,0])
Minimum_Sr = min(st_SDeviation[:,0])
Maximum_Rb = max(st_SDeviation[:,1])
Minimum_Rb = min(st_SDeviation[:,1])
Sigma_Sr = (Maximum_Sr - Minimum_Sr) / 2
Sigma_Rb = (Maximum_Rb - Minimum_Rb) / 2
Sigma_hL_Sr = (np.log(2) / (Sr_min)**2 ) * Sigma_Sr # Uncertainty on the half life of Sr
Sigma_hL_Rb = (np.log(2) / (Rb_min)**2 ) * Sigma_Rb # Uncertainty on the half life of Rb


""" The next block of code plots the imported data with
errobars. It also includes the best fit of the model, including outliers.
"""

Raw_Data = plt.figure('Raw Detector Data with the Best Fit', figsize = (9.4,8.8))
Sizes = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) # We impose a difference in heights between the residual plot and the main fit plot
Raw_Data.canvas.manager.window.move(10,10) # This line manipulates the starting position of the figure with respect to the screen
Detector_Data = Raw_Data.add_subplot(Sizes[0], label = 'DD')
Detector_Data.grid(True)
Detector_Data.set_title(r'Activity of $^{79}$Rb with respect to time')
Detector_Data.set_ylabel(r'Activity [Bq]')
Detector_Data.set_xlabel(r'Time [s]')
Detector_Data.errorbar(Faulty_Data[:,0], Faulty_Data[:,1], yerr = Faulty_Data[:,2], fmt = '. k', barsabove = False, capsize = 1.5, ecolor = 'red', elinewidth = 0.5, capthick = 1)
Detector_Data.plot(Data[:,0], Activity_Rb(Data[:,0], Sr_min, Rb_min), color = 'blue', linewidth = 2.5, alpha = 0.75)

Residuals = Raw_Data.add_subplot(Sizes[1]) # New subplot that provides information on residuals
Residuals.set_ylabel('Residuals')
Residuals.errorbar(Data[:,0], Activity_Rb(Data[:,0], Sr_min, Rb_min) - Data[:,1], yerr = Data[:,2], fmt = '. k')
Residuals.plot(Data[:,0], 0*Data[:,0])

if 0.5 < Chi_Squared_Min / Deg_Freedom < 2: # Depending on how good the reduced chi squared is, the text includes a fit assessment message
    Residuals.text(0.02, 0.05, r'[$\chi^{2}$' + '= {:4.2f} | '.format(Chi_Squared_Min) + r'Reduced $ \chi^2 = $' + 
               '{:4.2f} | '.format(Chi_Squared_Min / (Deg_Freedom)) 
               + 'Degrees of Freedom: {:}]'.format(Deg_Freedom) + r'$\Rightarrow$' + 'Good Fit', transform=plt.gcf().transFigure)
else:
    Residuals.text(0.02, 0.05, r'[$\chi^{2}$' + '= {:4.2f} | '.format(Chi_Squared_Min) + r'Reduced $ \chi^2 = $' + 
               '{:4.2f} | '.format(Chi_Squared_Min / (Deg_Freedom)) 
               + 'Degrees of Freedom: {:}]'.format(Deg_Freedom) + r'$\Rightarrow$ ' + 'Bad Fit', transform=plt.gcf().transFigure)

Residuals.text(0.02, 0.02, r'[$\lambda_{Sr}$' + r'= ({:5.5f} $\pm$ {:5.5f}) s$^-$$^1$ | '.format(Sr_min, Sigma_Sr) + r'$\lambda_{Rb}$ = ' + 
               r'({:6.6f} $\pm$ {:6.6f}) s$^-$$^1$ | '.format(Rb_min, Sigma_Rb) 
               + r'$\tau_{Sr}$' +' = ({:3.3f} $\pm$ {:3.3f}) min | '.format(hL_Sr / 60, Sigma_hL_Sr / 60) +
               r'$\tau_{Rb}$ = ' + '({:3.3f} $\pm$ {:3.3f}) min]'.format(hL_Rb / 60, Sigma_hL_Rb / 60), transform=plt.gcf().transFigure)
    
Raw_Data.savefig('Raw Detector Data and Fit.png', dpi = 300, orientation = 'landscape') # This line saves the figure on the directory the script is located

"""
Messages after the end of Analysis / Presentation of Results
"""
print('\nThe results of the experiment are presented in the table below with their associated uncertainties.')
print('\nChi Squared: {:4.2f}\nReduced Chi Squared: {:4.2f} (N = {:4.0f})'.format(Chi_Squared_Min, Chi_Squared_Min / (Deg_Freedom), Data_Points ))
print('\nFor the element Sr | Lambda Factor: ({:5.5f} ± {:5.5f}) /s  | Half Life: ({:3.3f} ± {:3.3f}) min'.format(Sr_min,Sigma_Sr, hL_Sr / 60, Sigma_hL_Sr / 60))
print('For the element Rb | Lambda Factor: ({:6.6f} ± {:6.6f}) /s | Half Life: ({:3.3f} ± {:3.3f}) min'.format(Rb_min, Sigma_Rb, hL_Rb / 60, Sigma_hL_Rb / 60))
print('\nEnd of application')
print('_'*100)


