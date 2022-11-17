from easygui import *
import numpy as np
import os
import xraylib
from fractions import Fraction
from scipy.optimize import curve_fit
from PIL import Image
import pandas as pd
import sys
import warnings
warnings.filterwarnings("ignore")
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 15)


def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)
# imges must be in the same folder
image = resource_path('zcompare.png')
image2 = resource_path('dukesim_matfile.png')

msg = "Choose from NIST database or import your own compound:"
choices = ["NIST compounds","My compounds","Exit"]

reply = buttonbox(msg,title = 'Welcome!', image=image, choices=choices)    
if reply == 'NIST compounds':
    nistCompounds = xraylib.GetCompoundDataNISTList()
    mycompounds = multchoicebox("Choose compounds from NIST database:",title = 'List of NIST compounds', choices=nistCompounds)
    if mycompounds == None:
        sys.exit()
elif reply == 'My compounds':

    while True:
        try:
            mypath = multenterbox(msg = 'Enter path to the file.', title = 'Load excel file (.xlsx, .xls)', fields = ['Path to the folder: ', 'File name: '])
            mypath = "".join(mypath)
            mypath = mypath.replace("\\","/")
            if mypath.endswith('csv'):
                mydf = pd.read_csv(mypath, index = True) 
            else:
                mydf = pd.read_excel(mypath, index = True) 
        except:
            if mypath==None:
                sys.exit()
            msgbox('Please insert correct path.')
            continue
        else:
            break
    
    mycompounds = list(mydf.iloc[:,0])
    mydf.index = mycompounds
elif reply == 'Exit':
    sys.exit()

#available methods
spiers = 'Spiers et al. 1946'
glasser = 'Glasser et al. 1947'
hine = 'Hine et al. 1952'
tsa = 'Tsa and Cho 1976 (E < 150 keV)'
pul = 'Puumalainen et al. 1977' 
sirz = 'Champley et al. (SIRZ-2)'

methods = [spiers,glasser,hine,tsa,pul,sirz]
mymethods = multchoicebox("Choose desired methods for Zeff calculation:", choices=methods)
if mymethods==None:
    sys.exit()
if sirz in mymethods:
    Emin = integerbox(msg='For SIRZ-2 method, energy range is needed.\n Enter minimum energy (keV)', lowerbound=1, upperbound=600)
    Emax = integerbox(msg='Enter maximum energy (keV)', lowerbound=Emin+1, upperbound=600)
    Es = range(Emin,Emax)
    def Ze_fun(E,b):
        return ((1-b)*curve1(E) + b*curve2(E))
    Es = [int(E) for E in Es]
if spiers in mymethods:
    myp = enterbox(msg='Enter exponential (p) value for Spiers et al. method:\n Default is p = 2.94, Tsa and Cho (1976) use p = 3.1.',default = 2.94)
    if myp==None:
        sys.exit()
    myp = float(myp)
    p = myp
ptc = 3.1
Zarray = np.zeros((len(mycompounds),len(methods)+1))

for m,material in enumerate(mycompounds):
    if reply == 'NIST compounds':
        #atomic numbers
        Znums = np.array(xraylib.GetCompoundDataNISTByName(material)["Elements"])
        Znums = np.array(Znums)
        #mass fractions
        mfracs = np.array(xraylib.GetCompoundDataNISTByName(material)["massFractions"])
        mfracs = np.array(mfracs).astype('float')
        #density
        rho = float(xraylib.GetCompoundDataNISTByName(material)["density"]) 
    if reply == 'My compounds':
        Znums = np.array(mydf.columns[1:-4])
        mfracs = list(mydf.loc[material][1:-4])
        rho = list(mydf.loc[material])[-4]
        
    #convert the mass fractions into the equivalent of a chemical formula
    #by dividing each element’s fraction-by-mass by its atomic mass
    aw = np.array([xraylib.AtomicWeight(Z_) for Z_ in Znums])

    #atomic fractions
    fracs = mfracs / aw
    #round floats to integer keeping the ratio
    #to get the fraction of the total number of electrons associated with each element
    fractions = [Fraction(val).limit_denominator(100) for val in fracs]
    ratios = np.array([(f.numerator, f.denominator) for f in fractions])
    factor = np.lcm.reduce(ratios[:,1])
    result = [round(v * factor) for v in fracs]
    weights = (result*Znums)/np.sum(result*Znums)
    Z_init = np.round((np.sum(mfracs*Znums**4)/np.sum(mfracs*Znums))**(1/3))
    Zarray[m,6] = rho
    if spiers in mymethods:
        Zs1 = Znums**p
        Zeff1 = np.sum(Zs1*weights)**(1/p)
        Zeff1 = Zeff1
        Zarray[m,0] = Zeff1
    if glasser in mymethods:
        Zeff2 = (np.sum(mfracs*Znums**4)/np.sum(mfracs*Znums))**(1/3)
        Zarray[m,1] = Zeff2
    if hine in mymethods:
        Zeff4 = np.sum((mfracs*Znums)/aw)/np.sum((mfracs)/aw)
        Zarray[m,2] = Zeff4
    if tsa in mymethods:
        Zs5 = Znums**ptc
        Zeff5 = np.sum(Zs5*weights)**(1/ptc)
        Zeff5 = Zeff5
        Zarray[m,3] = Zeff5
    if pul in mymethods:
        ntotal = np.sum(result)
        Zeff8 = np.sum((result/ntotal)*Znums)
        Zarray[m,4] = Zeff8
    if sirz in mymethods:
        ydata = []
        for E in Es:
            lam = 1.239842e-7 / E   
            ydata.append(np.sum(weights*np.array([(xraylib.AtomicWeight(Z)/Z)*xraylib.CS_Total(Z, E) for Z in Znums])))
        
        Zes_water = []
        Res_water = []
        Zs = np.arange(max(Z_init-5,3),min(Z_init+5,98),1)
        Zs = [int(Z) for Z in Zs]
        for ZZ in Zs:

            mus1 = [(xraylib.AtomicWeight(ZZ)/ZZ)*xraylib.CS_Total(ZZ, E) for E in Es]
            crv1 = np.polyfit(Es, mus1, deg = len(Es)-1);
            curve1 = np.poly1d(list(crv1))

            mus2 = [(xraylib.AtomicWeight(ZZ)/ZZ)*xraylib.CS_Total(ZZ+1, E) for E in Es]
            crv2 = np.polyfit(Es, mus2, deg = len(Es)-1);
            curve2 = np.poly1d(list(crv2))


            pars_,covb_  = curve_fit(Ze_fun,Es,ydata)  #,bounds=(0,1)
            b_ = pars_[0]
            if b_ < 0 or b_>1:
                w_b = 0.9   #penal residuals if b_ is not in range 0-1
                residuals = np.sqrt(np.sum((np.array(ydata) - w_b*np.array([Ze_fun(E,b_) for E in Es]))**2))
            else:
                residuals = np.sqrt(np.sum((np.array(ydata) - np.array([Ze_fun(E,b_) for E in Es]))**2))

            #Ze = Z+b_
            Zes_water.append(ZZ+b_)
            Res_water.append(residuals)
        res_min = Res_water.index(min(Res_water))
        winner = Zs[res_min]
        Zarray[m,5] = Zes_water[res_min]

        
df = pd.DataFrame(Zarray, index=mycompounds, columns=methods+[r'Density (g/cm3)'])
df = df.loc[:, (df != 0).any(axis=0)]

msgbox(df)
while True:
    try:
        path = multenterbox(msg = 'Do you want to save your data?', title = 'Save (.csv,.xlsx,.xls)', fields = ['Path to the folder: ', 'File name: '])
        path = "".join(path)
        if path.endswith('csv'):
            df.to_csv(path)
        else:
            df.to_excel(path)
    except:
        if path==None:
            sys.exit()
        msgbox('Please use one of the available extensions: .csv, .xlsx or .xls')
        continue
    else:
        break