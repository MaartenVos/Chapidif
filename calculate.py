import os
import sys
import math
import tkinter as tk
import ctypes
import numpy as np
import time
import constants as cnst
from ctypes import c_double#, POINTER


class calculate:
   # def __init__(self, *args, **kwargs):
    def __init__(self,parent):
        self.MyChapApp=parent
        #self.master=master
        self.initialized = False
        #print (self.master.test)

        # we start with the variables the user is encouraged to change to fit his/her perpose

        # maxoscillators: if you change this you have to recompile the C DLL and change MAXOSC
        # in epslib.cpp to the same value
        self.maxoscillators: int = 250
        # Number of GOS components. if you change this you have to recompile
        # the C DLL after changing MAXGOS to the same value. (changes will
        # cause problems with the layout of the user interface)
        self.maxGOS: int = 10
        # Number of Kaneko components. if you change this you have to recompile
        # the C DLL after changing MAXGOS to the same value. (changes will
        # cause problems with the layout of the user interface)
        self.maxBelkacem: int = 5
        self.maxKaneko: int = 250 # should bethe same as self.maxoscillators:

        # the number of values for which IMFP,stopping and straggling is
        # calculated
        self.NStopping: int = 50  # e.g. 50
        # first energy considered when calculating energy dependence IMFP,
        # stopping and straggling for protons (keV)
        # self.first_proton_energy: float = 5.0  #was 1
        # the same, but for electrons as projectile
        # self.first_electron_energy: float = 0.01
        # the  Beam energy increment factor for which imfp, stopping and
        # straggling are calculated
        self.first_electron_energy: float =round( 0.5 * 0.5**2 * cnst.HARTREE / 1000.0,6)  # E (keV) corresponding to v=0.5 a.u.
        self.first_proton_energy: float =round(0.5 * 0.5**2 * cnst.Mp * cnst.HARTREE / 1000.0,3) # E (keV) corresponding to v=0.5 a.u.

        self.IncrFactor: float = 1.25 # was 1.25
 
        # factor how the width of the 'Penn oscillators' increase with energy
        # loss.
        self.lin_cont_deltaE_Penn: float = 0.025
        self.highest_momentum_considered_stopping: float = 1e6

        # the following parameters affect the quanc8 integration   used to
        # calculate the DIIMFP and hence IMFP, stopping and straggling  and the
        # user may want to check if it fits his/her requirements
        # relates to quanc8 integration. Larger values more precise but take a
        # little longer

        # if the kinematically allowed q_max is smaller then the kinematically allowed limit is used instead
        # same variable but for the SEP, rather than the DIIMFP
        self.max_q_considered_surface: float = 200

        # stepsize for integration of DIIMFP:stepsize= minimum stepsize
        # (self.Stepsize) + lin_cont_deltaE*omega
        # for very small q there are instabilities in the Mermin and MLL
        # algoritm. To avoud this we use a DL instead if omega * q < self.q_transition
        self.q_transition: float = 0.0001
        #same sort of thing, for Kaneko
        self.UMax = 2000.0
        # for very small q.  DL,  Mermin and MLL are the same in the optical
        # limit. 0.03 seems a reasonable value
        self.Compton_k_limit = 10  # was 10
        # the range over which the Comptonprofile is calculated
        self.q_Compton = 60 
        self.omega_ddcs = 100
        self.theta_max= 0.6  # in mrad, max angle theta ddcs color plot
        self.theta_ddcs = 0.5 # in mrad angle for which ddcs is calculated
        # momentum transfer for which the Compton profile is calculated
        # should be much larger than  self.Compton_k_limit

        # end of user-editable variable section
        # ===================================================================================================================
        mypath = os.path.dirname(os.path.realpath(__file__)) + "/"

        if sys.maxsize == 2147483647:
            if os.name == "posix":
                mypath = mypath + "clib/epslib-lin32.so"
            else:  # os.name == 'nt':
                mypath = mypath + "clib/epslib-win32.dll"
        if sys.maxsize == 9223372036854775807:
            if os.name == "posix" and (
                sys.platform == "linux2" or sys.platform == "linux"
            ):
                mypath = mypath + "clib/epslib-lin64.so"
            elif os.name == "posix" and sys.platform == "darwin":
                mypath = mypath + "clib/epslib-mac.so"
            else:  
                mypath = (
                    os.path.dirname(os.path.realpath(__file__))
                    + "\\clib\\epslib-win64.dll"
                )

        self.epslib = ctypes.CDLL(mypath)
        self.epslib.lossfunction_inclGOS_atE.restype = ctypes.c_double
        self.epslib.DSEP.restype = ctypes.c_double
        self.epslib.DDCS_total.restype = ctypes.c_double
        self.epslib.DDCS_longitudinal.restype = ctypes.c_double

        self.DFChoice = tk.IntVar(value=3)
        self.DFChoice.trace_add("write", self.updatedensity)
        self.Projectile = tk.IntVar(value=0) # 1 means proton
  
        self.particle="electron"
        self.ProjectileMass = 1
        self.Projectile.trace_add("write", self.updateProjectileEnergy)
        
        self.Merminize = tk.IntVar(value=1) #0 means Lindhard, 1 means Merminize, 2 Direct approach
        self.RadiativeLosses = tk.IntVar(value=0) # 0 means no radiative losses else Bethe-Heitler approach 
        self.MottCorrection = tk.IntVar(value=0) # 0 means no radiative losses else Bethe-Heitler approach 
        self.MerminText = tk.StringVar(value="Mermin or MLL")
        self.ExchangeCorrection = tk.IntVar(value=0)  # 0 means no exchange correction for electrons
        self.ExchangeCorrectionMethod = tk.IntVar(value=0)  # 0 means no exchange correction for electrons
        self.StoppingUnits=tk.IntVar(value=0)  #0 means eV/Angstrom, 1 = eV / (1E15 atoms/cm2)
        self.DebugMode=tk.IntVar(value=0) # 0 do not print debug information
        self.AddELF = tk.IntVar(value=1) # 1 means add ELF, 0 means add Chi
   
        self.Stopping_calc_quality = tk.IntVar(value=0)
        self.StoppingFactor=1.0  # for the stopping units. 1.0 corresponds to eV/Angstrom. Set in run and plot
        self.Approximations = tk.IntVar(value=0)
        self.ApplySumRuleToGOS = tk.IntVar(value=1) # 1 means we renormalise
        self.FineMeshFactor = tk.IntVar(value = 9)  # KK and sum rules are evaluated on a fine mesh than the plots, stepsize is reduced by this factor

        self.maxEnergyDensityEffect = tk.DoubleVar(value=0.0)  
        self.E0 = 50.0  # beam energy (keV)
        self.ProjectileEnergy=tk.DoubleVar(value=50.0) # copy of beam energy for user interface
        self.ProjectileEnergy.trace_add("write", self.updateProjectileEnergy)
        self.BE_for_exchange=0.0
        
        self.specificweight = tk.DoubleVar(value=2.7)
        self.specificweight.trace_add("write", self.updatedensity)
        self.massunitcell = tk.DoubleVar(value=27.0)

        self.massunitcell.trace_add("write", self.updatedensity)
        mole_per_cm3 = self.specificweight.get() / self.massunitcell.get()
        self.UnitCellDensity = mole_per_cm3 * cnst.NAvogadro / 1.0e24
        self.UnitCellDensityText = tk.StringVar()
   
        self.plasmon_1elec_per_unit_cell = tk.StringVar()
 
        self.Eres           = tk.DoubleVar(value=1.0)
        self.coef1          = tk.DoubleVar(value=0.0)
        self.coef2          = tk.DoubleVar(value=0.0)
        self.coef3          = tk.DoubleVar(value=0.0)
        self.EELS_thickness = tk.DoubleVar(value=1000)  # thickness for transmission EELS
        self.EELS           = tk.IntVar(value=0)  #EELS=1 is TEELS,  0 means REELS

        self.thetaIn        = tk.DoubleVar(value=0.0)
        self.thetaOut       = tk.DoubleVar(value=45.0)
        self.surf_ex_factor = tk.DoubleVar(value=1.0)
        self.fraction_DIIMFP= tk.DoubleVar(value=1.0)
        self.beta_r= 0.0
        self.axis_choice = tk.IntVar(value = 1)

         
        self.NDFPAR = 5 * self.maxoscillators + 4 * self.maxGOS + 4 * self.maxBelkacem + 6 * self.maxKaneko
# length of all the oscillaotor, Gos etc arrays filled in "fill_oscillators()"
         
        self.ParArray = (ctypes.c_double * (self.NDFPAR + 100))() # 100 filled in 'fill_remainder()"
        self.PartIntSum = np.zeros(10)
        self.OOSEnergy = np.zeros(1000)
        self.OOS = np.zeros(1000)
        self.ELF_wide_Energy = (ctypes.c_double * 1300)()  # should match with epslib
        self.ELF_wide = (ctypes.c_double * 1300)()   # should match with epslib
        self.ELF_wideResultArray = (ctypes.c_double * 7)() # for c1,c2,c3, i1, mie and i3
        
        self.N_OOS=1
        self.commentline = ""
        self.LowerELimit=0.0
        self.Stepsize = 0.5  # stepsize
        self.CenterFirstBin = self.LowerELimit +0.5*self.Stepsize # first energy of plot
        self.LastEnergy = 100.0
        self.LastMomentum = 3.0
        self.Stepsize_qplot = 0.03
        self.Energy_qplot = 0.5
        self.NPoints = int((self.LastEnergy - self.CenterFirstBin) / self.Stepsize) + 1
        self.Npointscomparray = 0
        self.q = 0.02
        self.MT_radius = 1.0
        self.precision = 10.0

        self.Overplot = tk.IntVar(value=0)

        self.Eps_bkg = tk.StringVar()  # background dielectric constant  (Drude model)
        self.Eps_bkg.set("1.0")

        self.GOSdenstext = tk.StringVar(value="")
        self.Kanekodenstext = tk.StringVar(value="")
        self.Belkacemtext = tk.StringVar(value="")
        self.myconversiontext1 = tk.StringVar(value="")
        self.myconversiontext2 = tk.StringVar(value="")
        self.myvelocitytext=tk.StringVar(value="")
        self.DF_prop_text = tk.StringVar(value="\n\n")
        self.IntegratePhi = 0
        self.log_choice = tk.IntVar(value=0)
        self.imfp_crosssection_choice = tk.IntVar(value=0)
        self.Energy_Scale_choice = tk.IntVar(value=0)
        self.Kaneko_choice = tk.IntVar(value=0)
        self.Dispersion_choice = tk.IntVar(value=0)
        self.Dispersion_relativistic = tk.IntVar(value=1)
        self.DSEP_choice = tk.IntVar(value=0)
        self.ErrorMessage = tk.StringVar(value="")


        self.weight = 0.0
        self.BohrStraggling = 1.0
        self.StoppingResultArray = (ctypes.c_double * 3)()
        self.xArray = np.zeros(30)
        self.yArray = np.zeros(30)
        self.yArray_relativistic = np.zeros(30)
        self.Nst = 0  # number of impact parameter dependent stopping values calculated (put to 100 in calc.)
        # allocate the variables
        self.Amps = []
        self.sumAi=0.0
        self.Omegas = []
        self.Alphas = []
        self.Gammas = []
        self.Us = []
        for i in range(self.maxoscillators):
            DummyA = tk.StringVar()
            self.Amps.append(DummyA)
            DummyB = tk.StringVar()
            self.Omegas.append(DummyB)
            DummyC = tk.StringVar()
            self.Gammas.append(DummyC)
            DummyD = tk.StringVar()
            self.Alphas.append(DummyD)
            DummyE = tk.StringVar()
            self.Us.append(DummyE)

        self.ConcGOS = []
        self.EdgeGOS = []
        self.nlGOS = []
        self.ZGOS = []
        for _ in range(self.maxGOS):
            DummyA = tk.StringVar()
            self.ConcGOS.append(DummyA)
            DummyB = tk.StringVar()
            self.EdgeGOS.append(DummyB)
            DummyC = tk.StringVar()
            self.nlGOS.append(DummyC)
            DummyD = tk.StringVar()
            self.ZGOS.append(DummyD)
        self.Conc_Belkacem = []
        self.w_Belkacem = []
        self.gamma_Belkacem = []
        for _ in range(self.maxBelkacem):
            DummyA = tk.StringVar()
            self.Conc_Belkacem.append(DummyA)
            DummyB = tk.StringVar()
            self.w_Belkacem.append(DummyB)
            DummyC = tk.StringVar()
            self.gamma_Belkacem.append(DummyC)

        self.N_Kaneko = []  # number of electrons in a shell
        self.Edge_Kaneko = []  # minimum excitation energy of that shell (Archubi&Arista ext.)
        self.Q_Kaneko = []  # characteristic momentum of each shell
        self.width_Kaneko = []  # characteristic width of shell
        self.l_Kaneko = []  # angular momentum
        self.gamma_Kaneko = []  # plasmon position fudge parameter
        for _ in range(self.maxKaneko):
            DummyA = tk.StringVar()
            self.N_Kaneko.append(DummyA)
            DummyB = tk.StringVar()
            self.Edge_Kaneko.append(DummyB)
            DummyC = tk.StringVar()
            self.Q_Kaneko.append(DummyC)
            DummyD = tk.StringVar()
            self.l_Kaneko.append(DummyD)
            DummyE = tk.StringVar()
            self.gamma_Kaneko.append(DummyE)
            DummyF = tk.StringVar()
            self.width_Kaneko.append(DummyF)


        self.SurfExProb = 0.0
        self.PlasmonE=1.0  # dummy value plasmon energy (eV) for 1 elec/unit cell
        self.bulk_eq = True  # dsep eq plot if false
        self.max_eq=0.0
        self.log_range=5
        self.w_p_TPP =tk.DoubleVar(value=15.0)
        self.w_p_Penn=tk.DoubleVar(value=0.0)
        self.l_Kaneko_transform=tk.IntVar(value=0.0)
        self.Q_Kaneko_transform=tk.DoubleVar(value=1.0)
        self.U_MLL_transform=tk.DoubleVar(value=10.0)
        self.N_oscillator_used=tk.IntVar(value= self.maxoscillators)
        self.comp_option_choice=tk.IntVar(value=0)
        self.my_updateProjectileEnergy()
      

        self.init()
        self.initialized = True
        self.SumAi = 0.0
        self.my_updatedensity() 
        self.add_traces()
    
    def add_traces(self):
        for i in range(self.maxoscillators):
            self.Amps[i].trace_add("write", self.updatedensity)
            self.Omegas[i].trace_add("write", self.updatedensity)
            self.Us[i].trace_add("write", self.updatedensity)
        for i in range(self.maxGOS):  
            self.ConcGOS[i].trace_add("write", self.updatedensity) 
        for i in range(self.maxBelkacem): 
            self.Conc_Belkacem[i].trace_add("write", self.updatedensity)
        for i in range(self.maxKaneko):
            self.N_Kaneko[i].trace_add("write", self.updatedensity) 

            
        
    
    def init(self):
        self.initialized = False
        
        for i in range(self.maxoscillators):
            self.Amps[i].set("0.0")
            self.Omegas[i].set("0.0")
            self.Gammas[i].set("0.0")
            self.Alphas[i].set("1.0")
            self.Us[i].set("0.0")
        for i in range(self.maxGOS):
            self.ConcGOS[i].set("0.0")
            self.EdgeGOS[i].set("10.0")
            self.nlGOS[i].set("10")
            self.ZGOS[i].set("6")

        for i in range(self.maxBelkacem):
            self.Conc_Belkacem[i].set("0.0")
            self.w_Belkacem[i].set("10")
            self.gamma_Belkacem[i].set("1.0")
        for i in range(self.maxKaneko):
            self.N_Kaneko[i].set("0.0")
            self.Q_Kaneko[i].set("2.2")
            self.width_Kaneko[i].set("1.0")
            self.Edge_Kaneko[i].set("0.0")
            self.l_Kaneko[i].set("0")
            self.gamma_Kaneko[i].set("1.0")
        self.initialized = True       
       

    def initParArray(self):
        error_code = self.fill_oscillators()
        if error_code == 0:
            error_code = self.fill_remainder()
        else:    
            ic(error_code)    
        return error_code    
            

    def fill_oscillators(self):
        self.ParArray[0] = float(self.Eps_bkg.get()) # not sure if this is still used
        self.ErrorMessage.set("")
        try:
            for i in range(self.maxoscillators):
                self.ParArray[5 * i + 1] = float(self.Amps[i].get())
                self.ParArray[5 * i + 2] = float(self.Omegas[i].get())
                self.ParArray[5 * i + 3] = float(self.Gammas[i].get())
                self.ParArray[5 * i + 4] = float(self.Alphas[i].get())
                self.ParArray[5 * i + 5] = float(self.Us[i].get())
        except ValueError as e:
            ic("error")
            self.ErrorMessage.set("Input error Oscillator no: " + str(i+1) + "\n"  + str(e))   
            return 1   
        offset = 5 * self.maxoscillators
        try: 
            for i in range(self.maxGOS):
                self.ParArray[4 * i + offset + 1] = float(self.ConcGOS[i].get())
                self.ParArray[4 * i + offset + 2] = float(self.ZGOS[i].get())
                self.ParArray[4 * i + offset + 3] = float(self.nlGOS[i].get())
                self.ParArray[4 * i + offset + 4] = float(self.EdgeGOS[i].get())
        except ValueError as e:
            self.ErrorMessage.set("Input error GOS no: " + str(i+1) + "\n"  + str(e))   
            return 1      
        offset += 4 * self.maxGOS
        try:
            for i in range(self.maxBelkacem):
                self.ParArray[4 * i + offset + 1] = float(self.Conc_Belkacem[i].get())
                self.ParArray[4 * i + offset + 2] = float(self.w_Belkacem[i].get())
                self.ParArray[4 * i + offset + 3] = float(self.gamma_Belkacem[i].get())
            offset += 4 * self.maxBelkacem
        except ValueError as e:
            self.ErrorMessage.set("Input error Belkacem no: " + str(i+1) + "\n"  + str(e))      
            return 1              
        try:
            for i in range(self.maxKaneko):
                self.ParArray[6 * i + offset + 1] = float(self.N_Kaneko[i].get())
                self.ParArray[6 * i + offset + 2] = float(self.Q_Kaneko[i].get())
                self.ParArray[6 * i + offset + 3] = float(self.width_Kaneko[i].get())
                self.ParArray[6 * i + offset + 4] = float(self.Edge_Kaneko[i].get())
                self.ParArray[6 * i + offset + 5] = float(self.l_Kaneko[i].get())
                self.ParArray[6 * i + offset + 6] = float(self.gamma_Kaneko[i].get())
            
                
        except ValueError as e:
            self.ErrorMessage.set("Input error Kaneko no: " + str(i+1) + "\n"  + str(e))   
            return 1
        return 0                

    def fill_remainder(self):
        # max_q_considered: the maximum momentum that is used in the integration in terms of q_min
        #  precision:          relates to quanc8 integration. Larger values more precise but take a
        # little longer
        error_code=0   # not really used yet
        if self.Stopping_calc_quality.get() == 0:
            max_q_considered = 3000
            self.lin_cont_deltaE = 0.08
            self.precision = 10.0

        elif self.Stopping_calc_quality.get() == 1:
            max_q_considered = 5000
            self.lin_cont_deltaE = 0.02
            self.precision = 100.0
        elif self.Stopping_calc_quality.get() == 2:
            max_q_considered = 1e20
            self.lin_cont_deltaE = 0.005
            self.precision = 400.0

        self.ParArray[self.NDFPAR + 1] = self.E0 * 1000.0
        # presicion in integration routine for DIIMFP etc.
        self.ParArray[self.NDFPAR + 2] = self.precision
        # maximum of qmax considered in terms of qmin for diimfp
        self.ParArray[self.NDFPAR + 3] = max_q_considered
        # default 500    # maximum of qmax considered in terms of qmin for dsep
        self.ParArray[self.NDFPAR + 4] = self.max_q_considered_surface
        self.ParArray[self.NDFPAR + 5] = self.NPoints
        self.ParArray[self.NDFPAR + 6] = self.CenterFirstBin
        self.ParArray[self.NDFPAR + 7] = self.Stepsize
        # used to be atom  density for GOS
        self.ParArray[self.NDFPAR + 8] = self.UnitCellDensity
        self.ParArray[self.NDFPAR + 9] = self.lin_cont_deltaE
        # lower and upper limit for diimfp integration. partial diimfp and
        # stopping manually overrides these settings
        self.ParArray[self.NDFPAR + 10] = 0
        self.ParArray[self.NDFPAR + 11] = self.highest_momentum_considered_stopping
        self.ParArray[self.NDFPAR + 12] = self.LastMomentum
        self.ParArray[self.NDFPAR + 13] = self.Stepsize_qplot

        self.ParArray[self.NDFPAR + 14] = self.q_transition
        self.ParArray[self.NDFPAR + 15] = float(self.ExchangeCorrection.get())
        self.ParArray[self.NDFPAR + 16] = float(self.Kaneko_choice.get())
        self.ParArray[self.NDFPAR + 17] = float(self.AddELF.get())
        self.ParArray[self.NDFPAR + 18] = float(self.ApplySumRuleToGOS.get())
        self.ParArray[self.NDFPAR + 19] = float(self.Merminize.get())
        self.ParArray[self.NDFPAR + 20] = float(self.Dispersion_choice.get())
        self.ParArray[self.NDFPAR + 21] = self.maxEnergyDensityEffect.get()
        # next parameterts for reels spectrum calculation
        self.ParArray[self.NDFPAR + 22] = self.Eres.get()  # fwhm here
        self.ParArray[self.NDFPAR + 23] = self.coef1.get()
        self.ParArray[self.NDFPAR + 24] = self.coef2.get()
        self.ParArray[self.NDFPAR + 25] = self.coef3.get()
        self.ParArray[self.NDFPAR + 26] = self.thetaIn.get()  # degree here
        self.ParArray[self.NDFPAR + 27] = self.thetaOut.get()
        self.ParArray[self.NDFPAR + 28] = self.surf_ex_factor.get()
        self.ParArray[self.NDFPAR + 29] = float(self.DSEP_choice.get())

        self.ParArray[self.NDFPAR + 30] = self.fraction_DIIMFP.get()
        self.ParArray[self.NDFPAR + 31] = float(self.Projectile.get())
        # 1=proton, 0=electron
        # late arrival dispersion treated  relativistically (yes=1)
        self.ParArray[self.NDFPAR + 32] = float(self.Dispersion_relativistic.get())
        
        self.ParArray[self.NDFPAR + 33] = self.theta_max
        self.ParArray[self.NDFPAR + 34] = float(self.MottCorrection.get())
        self.ParArray[self.NDFPAR + 35] = self.UMax  #controls transition Kaneko DL (a.u.)
        self.ParArray[self.NDFPAR + 36] = self.BE_for_exchange  
        self.ParArray[self.NDFPAR + 37] = float(self.ExchangeCorrectionMethod.get()) #exchange implemented in two ways 
        self.ParArray[self.NDFPAR + 38] = float(self.DebugMode.get()) #controls debugging output
        return 0

    

    def AllocateComparrays(self):

        self.xCompArray = np.zeros(self.Npointscomparray)
        self.yCompArray = np.zeros(self.Npointscomparray)

   


        


    # ============
    def eps1eps2(self):
        self.x_axis = (ctypes.c_double * self.NPoints)()
        self.Result1 = (ctypes.c_double * self.NPoints)()
        self.Result2 = (ctypes.c_double * self.NPoints)()
        self.epslib.GOS_Scaling_init(self.ParArray)
        for i in range(self.NPoints):
            self.x_axis[i] = self.CenterFirstBin + float(i) * self.Stepsize
        self.epslib.Eps1Eps2(self.ParArray, self.Result1, self.Result2,
            ctypes.c_double(self.q), ctypes.c_int(self.DFChoice.get()) )
            
    def eps1eps2_fine(self):  # this one is for sum rules and KK, to make sure there are no significant numerical errors
        self.length_fine = self.NPoints*self.FineMeshFactor.get()
        self.Stepsize_fine= self.Stepsize/self.FineMeshFactor.get()
        self.CenterFirstBin_fine = self.LowerELimit + 0.5 *self.Stepsize_fine # first energy of plot
        self.Result1_fine = (ctypes.c_double * self.length_fine)()
        self.Result2_fine = (ctypes.c_double * self.length_fine)()
        self.Result3_fine = (ctypes.c_double * self.length_fine)()
        self.Result4_fine = (ctypes.c_double * self.length_fine)()
        self.ParArray[self.NDFPAR + 5] = self.length_fine
        self.ParArray[self.NDFPAR + 6] = self.CenterFirstBin_fine
        self.ParArray[self.NDFPAR + 7] = self.Stepsize_fine
        
        self.epslib.GOS_Scaling_init(self.ParArray)
        self.epslib.Eps1Eps2(self.ParArray, self.Result1_fine, self.Result2_fine,
            ctypes.c_double(self.q), ctypes.c_int(self.DFChoice.get()) )            
     
                

    def oneovereps1eps2(self):
        self.eps1eps2()
        for i in range(self.NPoints):
            denominator = (self.Result1[i] * self.Result1[i] + self.Result2[i] * self.Result2[i])
            self.Result1[i] = self.Result1[i] / denominator
            self.Result2[i] = self.Result2[i] / denominator  #this is Im[-1/eps],the loss function, NOT im[1/eps]

    def eps1eps2_q(self):
        Nqstep = int(self.LastMomentum / self.Stepsize_qplot)
        self.x_axis =  (ctypes.c_double * Nqstep)()
        self.Result1 = (ctypes.c_double * Nqstep)()
        self.Result2 = (ctypes.c_double * Nqstep)()
        self.epslib.GOS_Scaling_init(self.ParArray)
        for i in range(Nqstep):
            self.x_axis[i] = 0.01 + float(i) * self.LastMomentum / Nqstep 

        self.epslib.Eps1Eps2_q(
            self.ParArray,
            self.Result1,
            self.Result2,
            ctypes.c_double(self.Energy_qplot),
            ctypes.c_int(self.DFChoice.get()),
        )

    def oneovereps1eps2_q(self):
        self.eps1eps2_q()
        Nqstep = int(self.LastMomentum / self.Stepsize_qplot)
        for i in range(Nqstep):
            denominator = (self.Result1[i] * self.Result1[i] + self.Result2[i] * self.Result2[i])
            self.Result1[i] = self.Result1[i] / denominator
            self.Result2[i] = self.Result2[i] / denominator #this is Im[-1/eps]the loss function NOT im[1/eps]
            
    def oneovereps2(self):

        self.x_axis = (ctypes.c_double * self.NPoints)()
        self.Result1 = (ctypes.c_double * self.NPoints)()
        self.Result2 = (ctypes.c_double * self.NPoints)()
    

        self.epslib.GOS_Scaling_init(self.ParArray)

        for i in range(self.NPoints):
            self.x_axis[i] = self.CenterFirstBin + float(i) * self.Stepsize

        self.epslib.lossfunction_inclGOS(self.ParArray, self.Result2, ctypes.c_double(self.q),  \
            ctypes.c_int(self.DFChoice.get()))        

    # ========================
    
    def eps_kk_test(self):
        self.FineMeshFactor.set(9)
        self.eps1eps2_fine()
        
        self.epslib.Kramers_Kronig_eps1_from_eps2(ctypes.c_double(self.CenterFirstBin_fine),
            ctypes.c_double(self.Stepsize_fine),ctypes.c_int(self.length_fine),
            self.Result2_fine, self.Result3_fine)
        self.epslib.Kramers_Kronig_eps2_from_eps1(ctypes.c_double(self.CenterFirstBin_fine),
            ctypes.c_double(self.Stepsize_fine),ctypes.c_int(self.length_fine),
            self.Result1_fine, self.Result4_fine)  
        #now cast the result on the normal grid we use for plotting  
        self.recast()
        #this wors too, but is a few times  bit slower
        # self.kk_result3=np.zeros(self.length_fine)
        # self.kk_result3=self.kk_numpy_eps1_from_eps2(self.Stepsize_fine, self.Result2_fine[:], 1e-6)
        # self.kk_result4=np.zeros(self.length_fine)
        # self.kk_result4=self.kk_numpy_eps2_from_eps1(self.Stepsize_fine, self.Result1_fine[:], 1e-6)
        # self.recast_numpy()
        
        
     #modified from https://github.com/utf/kramers-kronig/blob/master/kkr.py, attempt to teach me  to use numpy   
    def kk_numpy_eps1_from_eps2(self,de, eps_imag, cshift=1e-3):   

        """Calculate the Kramers-Kronig transformation on imaginary part of dielectric

        Doesn't correct for any artefacts resulting from finite window function.

        Args:
            de (float): Energy grid size at which the imaginary dielectric constant
                is given. The grid is expected to be regularly spaced.
            eps_imag (np.array): A numpy array with dimensions (n, 3, 3), containing
                the imaginary part of the dielectric tensor.
            cshift (float, optional): The implemented method includes a small
                complex shift. A larger value causes a slight smoothing of the
                dielectric function.

        Returns:
            A numpy array with dimensions (n, 3, 3) containing the real part of the
            dielectric function.
        """
        eps_imag = np.array(eps_imag)
        nedos = eps_imag.shape[0]
        cshift = complex(0, cshift)
        w_i = np.arange(0, nedos*de, de, dtype=np.complex_)
        w_i += 0.5*de

        def integration_element( w_r):  # this is an inner function
            factor = w_i / (w_i**2 - w_r**2 + cshift)
            total = np.sum(eps_imag * factor, axis=0)
            return total * (2/math.pi) * de + 1.0

        return np.real([integration_element(w_r) for w_r in w_i[:]])
        
        
     #modified from https://github.com/utf/kramers-kronig/blob/master/kkr.py   
    def kk_numpy_eps2_from_eps1(self,de, eps_real, cshift=1e-3):   

        eps_real = np.array(eps_real)
        nedos = eps_real.shape[0]
        cshift = complex(0, cshift)
        w_i = np.arange(0, nedos*de, de, dtype=np.complex_)
        w_i += 0.5*de

        def integration_element( w_r):  # this is an inner function
            factor = w_r / (w_i**2 - w_r**2 + cshift)
            total = - np.sum((eps_real - 1.0) * factor, axis=0)
            return total * (2/math.pi) * de 

        return np.real([integration_element(w_r) for w_r in w_i[:]])    
    
    def recast_numpy(self):    
        self.x_axis =  (ctypes.c_double * self.NPoints)()
        self.Result1 = (ctypes.c_double * self.NPoints)()
        self.Result2 = (ctypes.c_double * self.NPoints)()   
        self.Result3 = (ctypes.c_double * self.NPoints)()
        self.Result4 = (ctypes.c_double * self.NPoints)() 
        factor=self.FineMeshFactor.get()
        first_i=math.ceil(factor/2.0)-1  
        for i in range(self.NPoints):
            self.x_axis[i] = self.CenterFirstBin + float(i) * self.Stepsize  
            self.Result1[i] = self.Result1_fine[first_i+i*factor]
            self.Result2[i] = self.Result2_fine[first_i+i*factor]
            self.Result3[i] = self.kk_result3[first_i+i*factor]
            self.Result4[i] = self.kk_result4[first_i+i*factor]
        del self.Result1_fine
        del self.Result2_fine
        del self.Result3_fine
        del self.Result4_fine      
    
    
    
    def recast(self):    
        self.x_axis =  (ctypes.c_double * self.NPoints)()
        self.Result1 = (ctypes.c_double * self.NPoints)()
        self.Result2 = (ctypes.c_double * self.NPoints)()   
        self.Result3 = (ctypes.c_double * self.NPoints)()
        self.Result4 = (ctypes.c_double * self.NPoints)() 
        factor=self.FineMeshFactor.get()
        first_i=math.ceil(factor/2.0)-1  
        for i in range(self.NPoints):
            self.x_axis[i] = self.CenterFirstBin + float(i) * self.Stepsize  
            self.Result1[i] = self.Result1_fine[first_i+i*factor]
            self.Result2[i] = self.Result2_fine[first_i+i*factor]
            self.Result3[i] = self.Result3_fine[first_i+i*factor]
            self.Result4[i] = self.Result4_fine[first_i+i*factor]
        del self.Result1_fine
        del self.Result2_fine
        del self.Result3_fine
        del self.Result4_fine 
            
    
    def one_over_eps_kk_test(self):
        self.FineMeshFactor.set(9)
        self.eps1eps2_fine()
        for i in range(self.length_fine):
            denominator = (self.Result1_fine[i] * self.Result1_fine[i] + self.Result2_fine[i] * self.Result2_fine[i])
            self.Result1_fine[i] = self.Result1_fine[i] / denominator
            self.Result2_fine[i] = self.Result2_fine[i] / denominator  #this is Im[-1/eps],the loss function, NOT im[1/eps]
        scratch  = (ctypes.c_double *  self.length_fine)()
        for i in range( self.length_fine):
            scratch[i]= -self.Result2_fine[i]  # now scratch contains im (1/eps)
        self.epslib.Kramers_Kronig_eps1_from_eps2(ctypes.c_double(self.CenterFirstBin_fine),
            ctypes.c_double(self.Stepsize_fine),ctypes.c_int(self.length_fine),
            scratch, self.Result3_fine)
        self.epslib.Kramers_Kronig_eps2_from_eps1(ctypes.c_double(self.CenterFirstBin_fine),
            ctypes.c_double(self.Stepsize_fine),ctypes.c_int(self.length_fine),
            self.Result1_fine, scratch)    
        for i in range( self.length_fine):
            self.Result4_fine[i]= -scratch[i]  # now result4 contains im (-1/eps) (obtained via KK)                
        self.recast()      

    def n_and_k_from_eps1_eps2(self):
        # check this! not sure what it means for q!=0
        # incoming result 1 and result2 is eps1 and eps2
        # transformed into n and k
        # from https://physics.stackexchange.com/questions/91776/\
        # real-and-imaginary-parts-of-dielectric-constant-vs-refractive-index
        # useful for comparisons with https://refractiveindex.info/
        #see also Wooton eq. 3.25,3.26
        self.eps1eps2()
        for i in range(self.NPoints):
            # a=np.longdouble(self.Result1[i])  # problems with precision here for q> 100 or so
            # b=np.longdouble(self.Result2[i])
            # first=np.sqrt(a ** 2 + b ** 2) / 2.0
            # second = a / 2.0
            # realresult=math.sqrt(first + second)
            # imagresult=  math.sqrt(first - second)
            # self.Result1[i]=float(realresult)
            # self.Result2[i] =float(imagresult)
            # the above works fine, here an alternative using pyshon's complex numbers
            epscomplex=complex(self.Result1[i],self.Result2[i])
            tmp=epscomplex**0.5
            self.Result1[i]= tmp.real
            self.Result2[i]= tmp.imag
            
    def n_and_k_kk_test(self):
        self.eps1eps2_fine()
        for i in range(self.length_fine):
            epscomplex=complex(self.Result1_fine[i],self.Result2_fine[i])
            tmp=epscomplex**0.5
            self.Result1_fine[i]= tmp.real
            self.Result2_fine[i]= tmp.imag
        self.epslib.Kramers_Kronig_eps1_from_eps2(ctypes.c_double(self.CenterFirstBin_fine),
            ctypes.c_double(self.Stepsize_fine),ctypes.c_int(self.length_fine),
            self.Result2_fine, self.Result3_fine)
        self.epslib.Kramers_Kronig_eps2_from_eps1(ctypes.c_double(self.CenterFirstBin_fine),
            ctypes.c_double(self.Stepsize_fine),ctypes.c_int(self.length_fine),
            self.Result1_fine, self.Result4_fine)      
        self.recast()  
            
            
            
      
    # ========================

    def calcDIIMFP(self):
        # calculate over the energy range defined in chapidif

        self.x_axis = (ctypes.c_double * self.NPoints)()
        self.Result1 = (ctypes.c_double * self.NPoints)()
        for i in range(self.NPoints):
            self.x_axis[i] = self.CenterFirstBin + float(i) * self.Stepsize

        self.epslib.GOS_Scaling_init(self.ParArray)
        self.epslib.DIIMFP(
            self.ParArray,
            self.Result1,
            ctypes.c_int(self.DFChoice.get()),
            self.StoppingResultArray,
        )

    def sum_rules(self):
        self.FineMeshFactor.set(9)
        self.eps1eps2_fine()
        ImEps  = (ctypes.c_double *  self.length_fine)()
        Ext  = (ctypes.c_double *  self.length_fine)()
        # start with k sum rule put result in result3
        for i in range(self.length_fine):
            ImEps[i]=self.Result2_fine[i]  # we need Im eps later
            epscomplex=complex(self.Result1_fine[i],self.Result2_fine[i])
            tmp=epscomplex**0.5
            Ext[i]= tmp.imag
        # now k in self.Result2_fine
        sumr3 = 0.0
        stepsize_au = self.Stepsize_fine / cnst.HARTREE
        for i in range(self.length_fine):
            E_au = (self.CenterFirstBin_fine + float(i) * self.Stepsize_fine)   / cnst.HARTREE
            sumr3 +=  Ext[i] * E_au * stepsize_au / cnst.PI**2 # number of electrons per a.u.^3
            self.Result3_fine[i] = sumr3/ ( self.UnitCellDensity * cnst.BOHR**3)  # now per atom (or more precise, per unit cell)
     
         #now  make it 1/ eps
        for i in range(self.length_fine):
            denominator = (self.Result1_fine[i] * self.Result1_fine[i] + self.Result2_fine[i] * self.Result2_fine[i])
            self.Result1_fine[i] = self.Result1_fine[i] / denominator
            self.Result2_fine[i] = self.Result2_fine[i] / denominator  #this is Im[-1/eps],the loss function, NOT im[1/eps]

        # integrate  result1 and LossFunction weighted by omega
        sumr1 = 0.0
        sumr2 = 0.0
        sumr4 = 0.0
        for i in range(self.length_fine):
            E_au =(self.CenterFirstBin_fine + float(i) * self.Stepsize_fine)   / cnst.HARTREE
            sumr1 = sumr1 + ImEps[i] * E_au * stepsize_au / 2.0 / cnst.PI**2  # number of electrons per a.u.^3
            sumr2 += self.Result2_fine[i] * E_au * stepsize_au / 2.0 / cnst.PI**2
            sumr4 += (2.0 / cnst.PI) * (self.Result2_fine[i] / E_au) * stepsize_au
            self.Result1_fine[i] = sumr1 / (self.UnitCellDensity * cnst.BOHR**3)
            self.Result2_fine[i] = sumr2 / (self.UnitCellDensity * cnst.BOHR**3)
            self.Result4_fine[i] = sumr4
           # print(sumr1,sumr2,sumr3)
        print("p sum rule: ", self.Result4_fine[self.NPoints - 1])
        # put back on standard grid
        self.recast()
        
    def inertial_rules(self):
        self.FineMeshFactor.set(29)
        self.eps1eps2_fine()
        ReEps        = (ctypes.c_double * self.length_fine)()
        RefracIndex  = (ctypes.c_double * self.length_fine)()
       
       
        for i in range(self.length_fine):
            ReEps[i]=self.Result1_fine[i]  # we need Re eps later
            epscomplex=complex(self.Result1_fine[i],self.Result2_fine[i])
            tmp=epscomplex**0.5
            RefracIndex[i] = tmp.real


        # start with refrac index inertia sum rule put result in result3
        mysum=0.0
        for i in range(self.length_fine):
            cont_this_segment=(RefracIndex[i]-1.0)*self.Stepsize_fine
            mysum +=cont_this_segment
            self.Result3_fine[i]=mysum  # n(0mega)-1 sum rule now in result3_fine will be put in result2 later
            
       
        prefactor = -1.0/(2.0*cnst.PI**2)
        self.limitingvalue =0.0
        for i in range(4):
            E_eV=(self.CenterFirstBin_fine + float(i) * self.Stepsize_fine) 
            limit = self.Result2_fine[i]*E_eV /(4*cnst.PI)
            if limit > self.limitingvalue:
                self.limitingvalue=limit
                
        mysum=0.0
        for i in range(self.length_fine): 
            cont_this_segment=prefactor*(ReEps[i]-1.0)*self.Stepsize_fine  # result1 contains now eps1
            mysum +=cont_this_segment
            self.Result4_fine[i]=mysum  # Re eps-1 sum rule now in result4
            
      
        mysum=0.0     
        for i in range(self.length_fine): 
            denominator = (self.Result1_fine[i] * self.Result1_fine[i] + self.Result2_fine[i] * self.Result2_fine[i])
            Re_one_overeps = self.Result1_fine[i] / denominator
            cont_this_segment=(Re_one_overeps - 1.0)*self.Stepsize_fine  # result1 contains now eps1
            mysum +=cont_this_segment
            self.Result1_fine[i]=mysum  # Re 1/eps -1 sum rule now in result3
            self.Result2_fine[i]=self.Result3_fine[i]  # Only result1 result2 and result 4 and a limiting line are plotted
        self.recast()    
           


    def Mean_Ionisation_Energy(self):
        # make sure we have the latest value,  but always evaluate at zero
        # momentum
        oldq = self.q
        self.q = 0.01  # better avoid 0
        self.oneovereps2()
        self.q = oldq
        # integrate Im -1/eps in  result1 and result 2 weighted by omega and
        # omega log omega
        self.Result3 = (ctypes.c_double * self.NPoints)()  # result1 and result2 are taken care of in oneovereps2
        self.sumBethe = 0.0
        sumr1 = 0.0
        sumr2 = 0.0  # so we never divide by 0
        # also same procedure but now weighted by sqrt(omega) for I' IMFP
        sumr3 = 0.0
        sumr4 =0.0 # so we never divide by 0
        # and for straggling:
        sumr5 = 0.0
        sumr6 = 0.0
        self.Result1[0] = 0.0
        self.Result2[0] = 0.0
        self.Result3[0] = 0.0
        stepsize_au = self.Stepsize / cnst.HARTREE
        for i in range(self.NPoints):
            E_au =  self.x_axis[i] /  cnst.HARTREE
            # NB: C_w=2.0/(cnst.PI*wi_p_au)*self.Result2[i_w]
            # (Penn thing) so we loose one omega (E_au)in the running sum
            sumr1 += (
                2.0 / cnst.PI * self.Result2[i] * math.log(E_au) * stepsize_au
            )  # imfp
            sumr2 += 2.0 / cnst.PI * self.Result2[i] * stepsize_au

            sumr3 += ( 2.0 / cnst.PI * self.Result2[i] * E_au * math.log(E_au) * stepsize_au
            )  # stopping
            sumr4 += 2.0 / cnst.PI * self.Result2[i] * E_au * stepsize_au

            sumr5 += (2.0 / cnst.PI * self.Result2[i] * E_au * E_au * math.log(E_au) * stepsize_au
            )  # straggling
            sumr6 += 2.0 / cnst.PI * self.Result2[i] * E_au * E_au * stepsize_au

            self.sumBethe += self.Result2[i] * E_au * stepsize_au / 2.0 / cnst.PI**2  #seems the same as sumr4 / 4pi
            if sumr2 > 0:
                self.Result1[i] = cnst.HARTREE * math.exp(sumr1 / sumr2)  # in eV
            else:
                self.Result1[i] = self.Result1[i - 1]
            if sumr4 > 0:
                self.Result2[i] = cnst.HARTREE * math.exp(sumr3 / sumr4)
            else:
                self.Result2[i] = self.Result1[i - 1]
            if sumr6 > 0:
                self.Result3[i] = cnst.HARTREE * math.exp(sumr5 / sumr6)
            else:
                self.Result3[i] = self.Result3[i - 1]

        self.I0 = self.Result1[self.NPoints - 1]
        self.MIE = self.Result2[self.NPoints - 1]
        self.Istraggling = self.Result3[self.NPoints - 1]
 
        self.C0 = sumr2 / math.exp(sumr1 / sumr2)
        self.C1 = sumr4 / (math.exp(sumr3 / sumr4)) ** 2 
        self.C2 = sumr6 / (math.exp(sumr5 / sumr6)) ** 3
        n0 = ( self.I0 * self.I0 / (4 * cnst.PI * cnst.HARTREE**2))  # electron density when interpreting I0 as a plamon energy (in e^-per au^3)
        self.C0_from_BetheSum = ( self.sumBethe / n0)  # fraction fo space filled with this oscillator so number of electrons correct.
        n1 = self.MIE * self.MIE / (4 * cnst.PI * cnst.HARTREE**2)
        self.C1_from_BetheSum = self.sumBethe / n1
        n2 = self.Istraggling * self.Istraggling / (4 * cnst.PI * cnst.HARTREE**2)
        self.C2_from_BetheSum = self.sumBethe / n2

        # print("sumbethe", self.sumBethe, "sumr2", sumr4, "sumr4", sumr2, "sumr6", sumr6)
        # print("n_av from sum bethe per angstrom^3", self.sumBethe / (cnst.BOHR**3))
        # print(
            # "c0,c1,c2 from bethe",
            # self.C0_from_BetheSum,
            # self.C1_from_BetheSum,
            # self.C2_from_BetheSum,
        # )
        # print("c0,c1,c2", self.C0, self.C1, self.C2)
        # print("I0,I1,I2", self.I0, self.MIE, self.Istraggling)

    def calccurves(self, smallqonly):
        self.CurvesEnergy = np.zeros(self.NStopping)
        self.CurvesVelocity = np.zeros(self.NStopping)
        self.IMFPEnergy = np.zeros(self.NStopping)
        self.StoppingEnergy = np.zeros(self.NStopping)
        self.CrosssectionEnergy = np.zeros(self.NStopping)
        # assuming 'Bethe dispersion' (constant  at plasmon energy then free electron dispersion)
        self.BetheIMFPEnergy = np.zeros(self.NStopping)
        self.BetheStoppingEnergy = np.zeros(self.NStopping)
        self.BetheStoppingEnergy_Salvat = np.zeros(self.NStopping)
        self.BetheStragglingEnergy = np.zeros(self.NStopping) 
        self.L_0_Bethe=np.zeros(self.NStopping) 
        self.Straggling_Jackson=np.zeros(self.NStopping) 
        self.FCOR_Salvat=np.zeros(self.NStopping) 
        self.DL_IMFPaverage_Energy = np.zeros(self.NStopping)  # DL using single average oscillator
        self.DL_IMFP_sum_Energy = np.zeros(self.NStopping)

        self.DL_StoppingEnergy = np.zeros(self.NStopping)  # DL stopping based on average oscillator
        self.DL_Stopping_sum_Energy = np.zeros(self.NStopping)  # DL stopping based on sum oscillators
        self.DL_StragglingEnergy = np.zeros(self.NStopping)  # DL straggling based on average oscillator
        self.DL_Straggling_sum_Energy = np.zeros(self.NStopping)  # DL straggling based on sum oscillators
        self.DL_Stopping_from_ELF = np.zeros(self.NStopping)
        self.DL_IMFP_from_ELF = np.zeros(self.NStopping)
        self.DL_Straggling_from_ELF = np.zeros(self.NStopping)
        self.StragglingEnergy = np.zeros(self.NStopping)
        self.TPP_IMFPEnergy = np.zeros(self.NStopping)
        self.LinearApprox_lowE = np.zeros(self.NStopping)
        
       
        startime=time.time()
        self.epslib.GOS_Scaling_init(self.ParArray)
        q=0.01

        self.epslib.Loss_wide(self.ParArray,  self.ELF_wide_Energy,  self.ELF_wide,  
            ctypes.c_double(q), ctypes.c_int(self.DFChoice.get()), self.ELF_wideResultArray )
        self.C0 = self.ELF_wideResultArray[0]
        self.C1 = self.ELF_wideResultArray[1]
        self.C2 = self.ELF_wideResultArray[2]
        self.I0 = self.ELF_wideResultArray[3]
        self.MIE = self.ELF_wideResultArray[4]
        self.Istraggling = self.ELF_wideResultArray[5]
        self.sumBethe= self.ELF_wideResultArray[6]
        I1_au =self.MIE / cnst.HARTREE  #mean ionization energy used to separate close and distant collision (if smallqonly is true)

        if  self.particle == "proton":
            CurrentE = self.first_proton_energy
        else:
            CurrentE = self.first_electron_energy
 

        if smallqonly:
            self.ParArray[self.NDFPAR + 11] = math.sqrt(2.0 * I1_au)

        for Ecounter in range(self.NStopping):
            print("calculating step: ", Ecounter," Energy (keV):", CurrentE)
            self.CurvesEnergy[Ecounter] = CurrentE
            self.ParArray[self.NDFPAR + 1] = CurrentE * 1000
            self.CurvesVelocity[Ecounter] = self.velocity_projectile(CurrentE)   # V in a.u. (using relativistic kinematics)
 
            self.epslib.DIIMFP_for_stopping(self.ParArray, ctypes.c_int(self.DFChoice.get()),
                self.StoppingResultArray)
                
            self.IMFPEnergy[Ecounter] = self.StoppingResultArray[0]
            self.CrosssectionEnergy[Ecounter] = 1.0 / (self.UnitCellDensity * self.StoppingResultArray[0])
            self.StoppingEnergy[Ecounter] = self.StoppingResultArray[1]
            self.StragglingEnergy[Ecounter] = self.StoppingResultArray[2]
            CurrentE = CurrentE * self.IncrFactor
        self.ParArray[self.NDFPAR + 11] = self.highest_momentum_considered_stopping
        if self.RadiativeLosses.get()==1 and self.particle == "electron": 
            Z=round(self.Nelec_per_UC)
            print("add estimate radiative losses, Z assumed:",Z)
            for Ecounter in range(self.NStopping):
                TotalE_MeV = (self.CurvesEnergy[Ecounter] + 511.0)/1000.0
                Factor=1.0 +  Z * TotalE_MeV/800.0  # Bethe Heitler estimate, Nikjoo book Interaction Radiation Matter pg 109
                self.StoppingEnergy[Ecounter] = self.StoppingEnergy[Ecounter] * Factor
                
            
        else:
            print(" no radiativelosses added")   
            print(self.particle)
            print (self.RadiativeLosses.get())
        
        finishtime=time.time()
        print( "duration calculation", finishtime - startime, "sec")
        self.Stopping_Linear_V()
        self.Bethe_L_Salvat()
        self.calculate_Straggling_Jackson()
        if self.Projectile.get() == 0:
            self.TanumaPowellPenn()
        if self.Approximations.get() == 1:
            self.calculate_approximations()

    def calculate_approximations(self):
         
        n0 = ( self.I0 * self.I0 / (4 * cnst.PI * cnst.HARTREE**2))  # electron density when interpreting I0 as a plamon energy (in e^-per au^3)
        n1 = self.MIE * self.MIE / (4 * cnst.PI * cnst.HARTREE**2)
        n2 = self.Istraggling * self.Istraggling / (4 * cnst.PI * cnst.HARTREE**2)
        self.sumBethe = n1 * self.C1 
        self.C0_from_BetheSum = ( self.sumBethe / n0)  # fraction of space filled with this oscillator so number of electrons correct.
       
        self.C1_from_BetheSum = self.sumBethe / n1
        
        self.C2_from_BetheSum = self.sumBethe / n2
        
        n1 = self.MIE * self.MIE / (4 * cnst.PI * cnst.HARTREE**2)
        self.C1_from_BetheSum = self.sumBethe / n1
 

        # if self.Projectile.get() == 1:
            # mass = cnst.Mp
        # else:
            # mass = 1
        I1_au = self.MIE / cnst.HARTREE
        I0_au = self.I0 / cnst.HARTREE
        I2_au = self.Istraggling / cnst.HARTREE
        print("self.Istraggling", self.Istraggling)

       
       
     
        
        top0_zerowidth = 0.0
        bottom0_zerowidth = 0.0
        top1_zerowidth = 0.0
        bottom1_zerowidth = 0.0
        top2_zerowidth = 0.0
        bottom2_zerowidth = 0.0
        self.OscillatorsPresent = False
        if self.DFChoice.get() != 1:  # not for ext. Drude 
            for i in range(self.maxoscillators):
                Amp_i = float(self.Amps[i].get())
                if Amp_i > 0:
                    self.OscillatorsPresent = True
                    wi_p_au = float(self.Omegas[i].get()) / cnst.HARTREE
                    top0_zerowidth += Amp_i * wi_p_au * math.log(wi_p_au)
                    bottom0_zerowidth += Amp_i * wi_p_au
                    top1_zerowidth += Amp_i * wi_p_au**2 * math.log(wi_p_au)
                    bottom1_zerowidth += Amp_i * wi_p_au**2
                    top2_zerowidth += Amp_i * wi_p_au**3 * math.log(wi_p_au)
                    bottom2_zerowidth += Amp_i * wi_p_au**3
            if self.OscillatorsPresent:
                I0_zerowidth = math.exp(top0_zerowidth / bottom0_zerowidth)
                C0_zerowidth = bottom0_zerowidth / I0_zerowidth
                I1_zerowidth = math.exp(top1_zerowidth / bottom1_zerowidth)
                C1_zerowidth = bottom1_zerowidth / (I1_zerowidth**2)
                I2_zerowidth = math.exp(top2_zerowidth / bottom2_zerowidth)
                C2_zerowidth = bottom2_zerowidth / (I2_zerowidth**3)
                print("for zero width oscillators:")
                print("C0=", C0_zerowidth, "I0=", I0_zerowidth * cnst.HARTREE, "eV")
                print("C1=", C1_zerowidth, "I1=", I1_zerowidth * cnst.HARTREE, "eV")
                print("C2=", C2_zerowidth, "I2=", I2_zerowidth * cnst.HARTREE, "eV")
                print("end calculate approximations")

        if  self.particle == "proton":
            ConstantA = 1.0/2.0  # notation so we are in line with imfp draft eq. 11
            ConstantB = 1
        else:
            ConstantA = 1
            ConstantB = 2
        for Ecounter in range(self.NStopping):
            E0_au = self.CurvesEnergy[Ecounter] * 1000 / cnst.HARTREE
            velocity = self.CurvesVelocity[Ecounter]
            prefactor = self.C1 * I1_au**2 / velocity**2


    

            # start code for stopping  and IMFP  straggling in high energy limit for the case of
            # 'Bethe disperion'  using the MIE as energy    , (NOT I0 for IMFP, I1 for straggling)
            q_max = velocity / ConstantA #non-relativistic case have to think about this
            
            q_min = I1_au / velocity
            if q_min > q_max:
                q_min = q_max  # this will make logBethe equal 0
                
           
            q_c = np.sqrt(2 * I1_au)
          
            logBethe = np.log(q_max / q_min)
            self.BetheStoppingEnergy[Ecounter] = prefactor * logBethe * cnst.HARTREE / cnst.BOHR
            
                
                
            gamma_r = 1.0 + E0_au/(self.ProjectileMass *(cnst.C)**2)    
            beta_r = np.sqrt((gamma_r * gamma_r -1.0)/(gamma_r*gamma_r))
            print("energy ",E0_au*cnst.HARTREE, "beta^2", beta_r * beta_r,    "np.log(gamma_r * gamma_r)", np.log(gamma_r * gamma_r)  )
     
            if logBethe > 0.0:
                tmp = math.log(q_c / q_min)
                tmp += 0.5 - I1_au / velocity**2

                oneoverlambda = self.C1 * I1_au / velocity**2 * tmp
                if oneoverlambda > 0.0:
                    self.BetheIMFPEnergy[Ecounter] = 1.0 / oneoverlambda  # in a.u.
                else:
                    self.BetheIMFPEnergy[Ecounter] = 1e30
                self.BetheIMFPEnergy[Ecounter] *= cnst.BOHR
            else:
                self.BetheIMFPEnergy[Ecounter] = 1e30
            if (q_c > q_min) and (q_max > q_c):
                logstraggling = np.log(q_c / q_min)
                self.BetheStragglingEnergy[Ecounter] = (
                    self.C1
                    * I1_au**2
                    / velocity**2
                    * (I1_au * logstraggling + (q_max * q_max - q_c * q_c) / 4.0)
                )
                self.BetheStragglingEnergy[Ecounter] *= (
                    cnst.HARTREE * cnst.HARTREE / cnst.BOHR
                )
            else:
                self.BetheStragglingEnergy[Ecounter] = 0
                
           
            # now IMFP and stopping and straggling based on more than one DL oscillator (ie  with gamma 0)
            oneoverlambda_total = 0.0
            stopping = 0.0
            straggling = 0.0
            for i in range(self.maxoscillators):
                Amp_i = float(self.Amps[i].get())
                if Amp_i != 0.0:
                    wi_p_au = float(self.Omegas[i].get()) / cnst.HARTREE
                    if velocity**2 > 2 * ConstantB * wi_p_au:
                        q_min_DL = (
                            velocity - np.sqrt(velocity**2 - 2 * ConstantB * wi_p_au)
                        ) / ConstantB
                        q_max_DL = (
                            velocity + np.sqrt(velocity**2 - 2 * ConstantB * wi_p_au)
                        ) / ConstantB
                    else:
                        q_min_DL = 1
                        q_max_DL = 1

                    logDL = math.log(q_max_DL / q_min_DL)
                    stopping += Amp_i * wi_p_au**2 / velocity**2 * logDL

                    tmp = logDL - 0.5 * math.log(
                        (2 * wi_p_au + q_max_DL**2) / (2 * wi_p_au + q_min_DL**2))
                    oneoverlambda_total += Amp_i * wi_p_au / velocity**2 * tmp

                    tmp = Amp_i * (wi_p_au**3 / velocity**2 * logDL)
                    fudge = (  1.0 ) # fudge here if equivalent oscillator has the wrong density

    
                    straggling += (fudge * Amp_i * 0.25 * wi_p_au**2 / velocity**2
                        * (q_max_DL**2 - q_min_DL**2) + tmp )

            self.DL_Stopping_sum_Energy[Ecounter] = stopping * cnst.HARTREE / cnst.BOHR
            if oneoverlambda_total > 0:
                self.DL_IMFP_sum_Energy[Ecounter] = (
                    1.0 / oneoverlambda_total
                ) * cnst.BOHR
            else:
                self.DL_IMFP_sum_Energy[Ecounter] = 1e20
            self.DL_Straggling_sum_Energy[Ecounter] = (
                straggling * cnst.HARTREE * cnst.HARTREE / cnst.BOHR
            )

            # start code for DL model single oscilator high-energy limit, first stopping
            if velocity**2 > 2 * ConstantB * I1_au:
                q_min_DL = (
                    velocity - np.sqrt(velocity**2 - 2 * ConstantB * I1_au)
                ) / ConstantB
                q_max_DL = (
                    velocity + np.sqrt(velocity**2 - 2 * ConstantB * I1_au)
                ) / ConstantB
            else:
                q_min_DL = 1
                q_max_DL = 1
            logDL = np.log(q_max_DL / q_min_DL)
            prefactor = I1_au**2 / velocity**2
            self.DL_StoppingEnergy[Ecounter] = (
                self.C1 * prefactor * logDL * cnst.HARTREE / cnst.BOHR
            )

            # now start code for IMFP from equivalent oscillator in high E limit
            if velocity**2 > 2 * ConstantB * I0_au:
                q_min_DL_IMFP = (
                    velocity - np.sqrt(velocity**2 - 2 * ConstantB * I0_au)
                ) / ConstantB
                q_max_DL_IMFP = (
                    velocity + np.sqrt(velocity**2 - 2 * ConstantB * I0_au)
                ) / ConstantB
            else:
                q_min_DL_IMFP = 1
                q_max_DL_IMFP = 1

            logDL_IMFP = np.log(q_max_DL_IMFP / q_min_DL_IMFP)

            tmp = logDL_IMFP - 0.5 * math.log(
                (2 * I0_au + q_max_DL_IMFP**2) / (2 * I0_au + q_min_DL_IMFP**2)
            )
            oneoverlambda = self.C0 * I0_au * tmp / velocity**2
            if oneoverlambda > 0.0:
                self.DL_IMFPaverage_Energy[Ecounter] = 1.0 / oneoverlambda  # in a.u.
            else:
                self.DL_IMFPaverage_Energy[Ecounter] = 1e20  # in a.u.
            self.DL_IMFPaverage_Energy[Ecounter] *= cnst.BOHR

            # now start code for straggling from equivalent oscillator in high E limit
            # this does not work for oscillators as then I2_au diverges, then use sum zero width oscillators
            if self.OscillatorsPresent:
                I2_au = I2_zerowidth
                self.C2 = C2_zerowidth
                self.Istraggling = I2_au * cnst.HARTREE

            if velocity**2 > 2 * ConstantB * I2_au:
                q_min_DL_Straggling = (
                    velocity - np.sqrt(velocity**2 - 2 * ConstantB * I2_au)
                ) / ConstantB
                q_max_DL_Straggling = (
                    velocity + np.sqrt(velocity**2 - 2 * ConstantB * I2_au)
                ) / ConstantB
            else:
                q_min_DL_Straggling = 1
                q_max_DL_Straggling = 1
            logDL_Straggling = np.log(q_max_DL_Straggling / q_min_DL_Straggling)
            self.DL_StragglingEnergy[Ecounter] = self.C2 * (
                I2_au**3 / velocity**2 * logDL_Straggling
            )
            # indeed I1_au^2 and C1 not I2_au^2 and C2 in last part, real average density not of average oscillator

            self.DL_StragglingEnergy[Ecounter] += (self.C1 * 0.25 * I1_au**2
                / velocity**2 * (q_max_DL_Straggling**2 - q_min_DL_Straggling**2))

            self.DL_StragglingEnergy[Ecounter] *= (cnst.HARTREE * cnst.HARTREE / cnst.BOHR)
        self.DL_stopping_IMFP_straggling_from_ELF()

    def DL_stopping_IMFP_straggling_from_ELF(self):
        #should work for all models
     
        if self.Projectile.get() == 1:
            ConstantB = 1
        else:
            ConstantB = 2
      
        
        for Ecounter in range(self.NStopping):
            velocity = self.CurvesVelocity[Ecounter]

            self.DL_Stopping_from_ELF[Ecounter] = 0.0
            one_over_lambda = 1e-20  # so we never divide by 0
            self.DL_Straggling_from_ELF[Ecounter] = 0.0
            for i_w in range(1300):  #should match number of steps in epslib.Loss_wide
                wi_p_au = self.ELF_wide_Energy[i_w] / cnst.HARTREE
                if i_w ==0:
                    Stepsize_au =  self.ELF_wide_Energy[0] / cnst.HARTREE
                else:    
                    Stepsize_au = (self.ELF_wide_Energy[i_w] - self.ELF_wide_Energy[i_w-1]) / cnst.HARTREE
                C_w = (2.0 / (cnst.PI * wi_p_au) * self.ELF_wide[i_w])  # G(omega) in Penn's paper for vanishing small gamma
                elec_dens = C_w * wi_p_au * wi_p_au / (4.0 * cnst.PI)
                prefactori = 4 * np.pi / velocity**2 * elec_dens
                if velocity**2 > 2 * ConstantB * wi_p_au:
                    q_min_DL = (
                        velocity - np.sqrt(velocity**2 - 2 * ConstantB * wi_p_au)
                    ) / ConstantB
                    q_max_DL = (
                        velocity + np.sqrt(velocity**2 - 2 * ConstantB * wi_p_au)
                    ) / ConstantB
                else:
                    q_min_DL = 1
                    q_max_DL = 1

                logDL = math.log(q_max_DL / q_min_DL)
                self.DL_Stopping_from_ELF[Ecounter] += (
                    prefactori * logDL * cnst.HARTREE / cnst.BOHR * Stepsize_au
                )
                tmp = logDL - 0.5 * math.log(
                    (2 * wi_p_au + q_max_DL**2) / (2 * wi_p_au + q_min_DL**2)
                )
                one_over_lambda += prefactori / wi_p_au * tmp * Stepsize_au
                self.DL_Straggling_from_ELF[Ecounter] += (
                    prefactori
                    * (wi_p_au * logDL + 0.25 * (q_max_DL**2 - q_min_DL**2))
                    * cnst.HARTREE**2
                    / cnst.BOHR
                    * Stepsize_au
                )

            self.DL_IMFP_from_ELF[Ecounter] = (1.0 / one_over_lambda) * cnst.BOHR
            
    def TanumaPowellPenn(self):
        
        bandgap = 0.0 #bandgap not implementd
        rho = self.specificweight.get()
        if self.w_p_TPP.get() > 0.0:
            beta_TPP_2m = (-0.1 + 0.944 / math.sqrt(self.w_p_TPP.get()**2 + bandgap * bandgap)
                + 0.069 * rho**0.1)  # assume no gap
        # from  TPP SIA 43 689 2011
        gamma_TPP_2m = 0.191 * rho**-0.5
        U = self.w_p_TPP.get()**2 / 829.4
        C_TPP_2m = 1.97 - 0.91 * U
        D_TPP_2m = 53.4 - 20.8 * U
        print("calculated Beta", beta_TPP_2m, "w_p^2*beta",   beta_TPP_2m * self.w_p_TPP.get()**2,
            "gamma_TPP_2m",  gamma_TPP_2m)
                # start code tor TPP IMFP 
       
        for Ecounter in range(self.NStopping):
            if (self.particle=="electron") and (self.w_p_TPP.get() > 0.0):  # electrons
                E_eV = self.CurvesEnergy[Ecounter] * 1000
                tmp = ( beta_TPP_2m * math.log(gamma_TPP_2m * E_eV) - C_TPP_2m / E_eV
                    + D_TPP_2m / E_eV**2)
                self.TPP_IMFPEnergy[Ecounter] = E_eV / (tmp * self.w_p_TPP.get()**2)
            else:
                self.TPP_IMFPEnergy[Ecounter] = 0.0    
                
                
    def calculate_Straggling_Jackson(self):            
     # Jackson Classical electrodynamics eq.13.50, Salvat PRA 2022 eq. 112, only for protons  
        for Ecounter in range(self.NStopping):
            CurrentE= self.CurvesEnergy[Ecounter] # in keV
            E_au=CurrentE*1000.0/cnst.HARTREE
            v=self.velocity_projectile(CurrentE) #subroutine also sets self.gamma_r and self.beta_r
            gamma2=self.gamma_r*self.gamma_r
            beta2=self.beta_r*self.beta_r
            curlyB=2 * cnst.PI/v**2   
          #  straggling_per_UC= curlyB * 2 * self.Nelec_per_UC *(gamma2*beta2*cnst.C**2*(1.0-beta2/2.0))
            straggling_per_UC= curlyB * 2 * self.Nelec_per_UC *(gamma2*v**2*(1.0-beta2/2.0))
            straggling_per_UC *=  cnst.HARTREE**2 * cnst.BOHR**2# now in eV^2/Angstrom^2 per UC
            self.Straggling_Jackson[Ecounter] = straggling_per_UC*self.UnitCellDensity # now eV^2/angstrom
            
    def Stopping_Linear_V(self):
         # a liear approximation
        tmp=self.StoppingEnergy/self.CurvesVelocity   
        Imaxsloop=np.argmax(tmp)
        self.maxsloop=self.StoppingEnergy[Imaxsloop]/self.CurvesVelocity[Imaxsloop]   
        for Ecounter in range(self.NStopping):
            if Ecounter < (Imaxsloop+1):
                self.LinearApprox_lowE[Ecounter]=self.CurvesVelocity[Ecounter]* self.maxsloop
            else:
                self.LinearApprox_lowE[Ecounter]=np.nan
                 
                
            
             
    def Bethe_L_Salvat(self):   
		#subroutine aiming to reproduce high energy limit sBethe from Salvat
     
        I_au=self.MIE / cnst.HARTREE
       
        for Ecounter in range(self.NStopping):
            prefactor =  4 * np.pi / self.CurvesVelocity[Ecounter] ** 2 * self.sumBethe * cnst.HARTREE / cnst.BOHR
            CurrentE= self.CurvesEnergy[Ecounter] # in keV
            E_au=CurrentE*1000.0/cnst.HARTREE
          #  print("I_au", I_au,"E_au",E_au,"sum Bethe", self.sumBethe)
            v=self.velocity_projectile(CurrentE) #subroutine also sets self.gamma_r and self.beta_r
            gamma=self.gamma_r
            gamma2=gamma*gamma
            beta2=self.beta_r*self.beta_r
            if  self.particle == "proton":
                R=1.0/(1.0+(1.0/cnst.Mp**2)+2.0*gamma/cnst.Mp)
                FCOR=math.log(R)+((gamma2-1.0)*R/(gamma*cnst.Mp))**2
                
            else:   
                FCOR=(2.0*gamma2-1.0)/gamma2+0.125*((gamma-1.0)/gamma)**2-(4.0-((gamma-1.0)/gamma)**2)*math.log(2.0)-math.log(gamma+1.0)
            
            self.L_0_Bethe[Ecounter]=math.log(2*v*v/ I_au) +  math.log(gamma2)-beta2
       
            if self.L_0_Bethe[Ecounter] < 0.0:
                self.L_0_Bethe[Ecounter] = 0.0
            self.FCOR_Salvat[Ecounter]= FCOR
            self.BetheStoppingEnergy_Salvat[Ecounter] = prefactor * (self.L_0_Bethe[Ecounter] + 0.5*FCOR)
            if self.BetheStoppingEnergy_Salvat[Ecounter] < 0.0:
                self.BetheStoppingEnergy_Salvat[Ecounter] = 0.0
            
           
               

    def shell_effect(self, smallqonly):
        self.calccurves(smallqonly)

        self.x_axis = (ctypes.c_double * self.NStopping)()  
        self.Result1 = (ctypes.c_double * self.NStopping)()
        self.Result2 = (ctypes.c_double * self.NStopping)()
        self.Result3 = (ctypes.c_double * self.NStopping)()
        self.Result4 = (ctypes.c_double * self.NStopping)()
      
        for Ecounter in range(self.NStopping):  
            prefactor =  4 * np.pi / self.CurvesVelocity[Ecounter] ** 2 * self.sumBethe * cnst.HARTREE / cnst.BOHR
            L_df = self.StoppingEnergy[Ecounter] / prefactor # recover L
            self.Result1[Ecounter] = L_df    # already restricted value if smallqonly == true
            if smallqonly:
                self.Result2[Ecounter] = 0.5 * (self.L_0_Bethe[Ecounter]+  0.5*self.FCOR_Salvat[Ecounter] )
            else:
                self.Result2[Ecounter] = self.L_0_Bethe[Ecounter]+  0.5*self.FCOR_Salvat[Ecounter] 
            if(self.Result2[Ecounter] < 0.0):
                self.Result2[Ecounter] = 0.0
                
            
            self.Result3[Ecounter] =  self.Result1[Ecounter] -self.Result2[Ecounter]   # "shell effect"
            self.Result4[Ecounter] = 0.5 * self.FCOR_Salvat[Ecounter]

    def CalcCompton(self):
        # actual momentum array but plotting routine expect x-axis in this array
        
            
        self.x_axis = (ctypes.c_double * self.NPoints)()
        self.Result1 = (ctypes.c_double * self.NPoints)()
         # k_limit is the largest momentum component of target electrons considered
        if( self.Compton_k_limit > 0.5 * self.q_Compton):
            Compton_k_limit_used = 0.5 * self.q_Compton   # so elower stays >=0
        else:
            Compton_k_limit_used  = self.Compton_k_limit
        if self.Dispersion_relativistic.get():
            Q_recoil_au =self.calc_Qrecoil(self.q_Compton)
        else:
            Q_recoil_au=self.q_Compton**2 / 2.0
                
        E_center = Q_recoil_au * cnst.HARTREE
        width_per_au =math.sqrt(2*Q_recoil_au)
        E_lower = E_center - Compton_k_limit_used * width_per_au * cnst.HARTREE
        E_upper = E_center + Compton_k_limit_used * width_per_au * cnst.HARTREE
 
        E_step = (E_upper - E_lower) / (self.NPoints - 1)

        self.ParArray[self.NDFPAR + 5] = self.NPoints
        self.ParArray[self.NDFPAR + 6] = E_lower
        self.ParArray[self.NDFPAR + 7] = E_step
        self.epslib.GOS_Scaling_init(self.ParArray)
        self.epslib.lossfunction_inclGOS(
            self.ParArray,
            self.Result1,
            ctypes.c_double(self.q_Compton),
            ctypes.c_int(self.DFChoice.get()),
        )
        sumi = 0.0
        E_step_au = E_step / cnst.HARTREE # point separationb in a.u. energy
        Mom_step = (2.0 * Compton_k_limit_used) / (self.NPoints - 1)

        for i in range(self.NPoints):
            # sum rule for DF replacing omega by Q_recoil_au
            self.Result1[i] = self.Result1[i] * Q_recoil_au / (2*cnst.PI**2)
            self.Result1[i] = (self.Result1[i] / cnst.BOHR**3 / self.UnitCellDensity)  # momentum density per atom not per a.u.^3
            self.x_axis[i] = -Compton_k_limit_used+ Mom_step * i
            self.Result1[i] = self.Result1[i] * E_step_au / Mom_step
            # sum rule result, print out if you want to check it
            sumi += self.Result1[i] * Mom_step
        print("sum rule", sumi)    

 

    def colorplot_lossfunction(self):

        current_q = 0.5 * self.Stepsize_qplot  # initial value
        Nqstep = int(self.LastMomentum / self.Stepsize_qplot)
        self.my_image = np.zeros((self.NPoints, Nqstep))
        self.x_axis = (ctypes.c_double * self.NPoints)()

        self.LossFunction = (ctypes.c_double * self.NPoints)()
        self.epslib.GOS_Scaling_init(self.ParArray)
        for i in range(self.NPoints):
            self.x_axis[i] = self.CenterFirstBin + float(i) * self.Stepsize

        for i in range(Nqstep):
            if self.bulk_eq:
                self.epslib.lossfunction_inclGOS(
                    self.ParArray,
                    self.LossFunction,
                    ctypes.c_double(current_q),
                    ctypes.c_int(self.DFChoice.get()))
            else:   
                self.epslib.SurfLossFunc(
                    self.ParArray,
                    self.LossFunction,
                    ctypes.c_double(current_q),
                    ctypes.c_int(self.DFChoice.get()))     

            for j in range(self.NPoints):
                self.my_image[j, i] = self.LossFunction[j]
            current_q = current_q + self.Stepsize_qplot
            
    def Calculate_Integration_limits(self):   
        MC2 = self.ProjectileMass* cnst.C*cnst.C
        for x in range(30):
            q = x * self.LastMomentum / 29.0
            self.xArray[x] = q
            # E0 is in keV, nonreativistic calculation
            p1_min = math.sqrt(2 * self.ProjectileMass * self.E0 * 1000.0 / cnst.HARTREE) - q
 
            omega = self.E0 * 1000.0 - p1_min * p1_min * cnst.HARTREE / (2 * self.ProjectileMass)
            self.yArray[x] = omega
            # and now we try to do the same using relativistic kinematics for the same energy losses
            p1_min_rel = self.p0_rel - q
            # use E^2=(pc)^2+(M_0 C^2), fingers crossed enough precision
            E_square = (p1_min_rel*cnst.C)**2 + (MC2)**2
            E_after = np.sqrt(E_square)
            E_kin_after=E_after - MC2
            omega = self.E0 * 1000.0 - E_kin_after * cnst.HARTREE
            self.yArray_relativistic[x] = omega
 
    def   scale_image(self):   
        self.my_scaled_image=np.copy(self.my_image)
        Nqstep = int(self.LastMomentum / self.Stepsize_qplot)         
        arraymax = np.amax(self.my_scaled_image)
        if self.max_eq > 0.0:
            arraymax = self.max_eq
            for j in range(self.NPoints):
                for i in range(Nqstep):
                    if ( self.my_scaled_image[j, i] > self.max_eq): 
                        self.my_scaled_image[j, i] =  self.max_eq 
            arraymax=self.max_eq       
  
        if self.log_choice.get() == 1:
            myrange=10**self.log_range
            logmin = math.log(arraymax / myrange, 10) 
            for j in range(self.NPoints):
                for i in range(Nqstep):
                    if self.my_scaled_image[j, i] > arraymax / myrange:
                        self.my_scaled_image[j, i] = math.log(self.my_scaled_image[j, i], 10)
                    else:
                        self.my_scaled_image[j, i] = logmin
        return                

        

   

    def colorplot_ddcs(self):

        Nqstep = int(self.LastMomentum / self.Stepsize_qplot)
        self.x_axis = (ctypes.c_double * Nqstep)()
        self.Result1 = (ctypes.c_double * Nqstep)()
        self.Result2 = (ctypes.c_double * Nqstep)()
        
        myarray = (ctypes.c_double * (self.NPoints * Nqstep))()
        self.my_image = np.zeros((self.NPoints, Nqstep))
        self.epslib.GOS_Scaling_init(self.ParArray)
      
        self.stopping_longitudinal=self.epslib.DDCS_longitudinal( self.ParArray, myarray,
            ctypes.c_int(self.DFChoice.get()), ctypes.c_int(self.IntegratePhi) )
        for i in range(self.NPoints):
            for j in range(Nqstep):
                self.my_image[i, j] = myarray[j * self.NPoints + i]
        
    def colorplot_Cerenkov(self):

        Nqstep = int(self.LastMomentum / self.Stepsize_qplot)
  
        
        
        myarray = (ctypes.c_double * (self.NPoints * Nqstep))()
        self.my_image = np.zeros((self.NPoints, Nqstep))
        self.epslib.GOS_Scaling_init(self.ParArray)
      
        self.stopping_total=self.epslib.DDCS_total(self.ParArray, myarray,
            ctypes.c_int(self.DFChoice.get()), ctypes.c_int(self.IntegratePhi) )
 
  
        for i in range(self.NPoints):
            for j in range(Nqstep):
                self.my_image[i, j] = myarray[j * self.NPoints + i]
                
                
    def difference_due_to_Cerenkov(self):  
        Nqstep = int(self.LastMomentum / self.Stepsize_qplot)
       
        self.colorplot_Cerenkov()    
        cerenkov_image= self.my_image
        self.colorplot_ddcs()
        self.my_image=cerenkov_image - self.my_image
        self.stopping_due_to_photons=self.stopping_total - self.stopping_longitudinal
        print(" stopping total", self.stopping_total, "stopping longitudinal",self.stopping_longitudinal )#, "due to photons:",self.stopping_due_to_photons )
        
          
        
        

      
    def DDCS_at_omega(self):
        Nqstep = int(self.LastMomentum / self.Stepsize_qplot)
        self.x_axis = (ctypes.c_double * Nqstep)()
        self.Result1 = (ctypes.c_double * Nqstep)()
        self.Result2 = (ctypes.c_double * Nqstep)()
        self.x_axis[0]= self.theta_max
        for i in range(Nqstep-1):
            self.x_axis[i+1] = self.x_axis[i]*0.95
        print("smallest theta", self.x_axis[Nqstep-1])    
        self.epslib.GOS_Scaling_init(self.ParArray)
 
        self.epslib.DDCS_at_omega(self.ParArray,self.x_axis,self.Result1,self.Result2,
            ctypes.c_int(self.DFChoice.get()),
            ctypes.c_double(self.omega_ddcs))
            
    def DDCS_at_theta(self):
        print("theta", self.theta_ddcs)
        self.x_axis = (ctypes.c_double * self.NPoints)()
        self.Result1 = (ctypes.c_double * self.NPoints)()
        self.Result2 = (ctypes.c_double * self.NPoints)()
        for i in range(self.NPoints):
            self.x_axis[i] = self.CenterFirstBin + float(i) * self.Stepsize
        self.epslib.GOS_Scaling_init(self.ParArray)
        
        self.epslib.DDCS_at_theta(self.ParArray,self.Result1,self.Result2,
            ctypes.c_int(self.DFChoice.get()),
            ctypes.c_double(self.theta_ddcs))            
            
    def DCS(self):

        Nqstep = int(self.LastMomentum / self.Stepsize_qplot)
        self.x_axis = (ctypes.c_double * Nqstep)()
        self.Result1 = (ctypes.c_double * Nqstep)()
        self.Result2 = (ctypes.c_double * Nqstep)()
        self.Result3 = (ctypes.c_double * Nqstep)()
        self.x_axis[0]= self.theta_max
        for i in range(Nqstep-1):
            self.x_axis[i+1] = self.x_axis[i]*0.95
        self.epslib.GOS_Scaling_init(self.ParArray)
        self.epslib.DCS(self.ParArray,self.x_axis,self.Result1,self.Result2,ctypes.c_int(self.DFChoice.get()) ) 
        E0_au= self.E0*1000.0/cnst.HARTREE 
        gamma_r = 1.0 + E0_au/(self.ProjectileMass*(cnst.C)**2)  
        for j in range(Nqstep): # put Rutherford in result3
            current_theta=self.x_axis[j]/1000 # in rad
            q=2.0*np.sin(current_theta/2.0) * self.p0_rel  ## for angles >> larger that theta_0=(|P_1|- |P_0|)/|p_0|
            self.Result3[j] = (4 * self.Nelec_per_UC * gamma_r**2 * self.ProjectileMass**2 / q**4 )/ cnst.BOHR**2 / (4.0 * cnst.PI)  # not sure about the 4 pi   

 
    def surfaceloss(self):
        self.x_axis = (ctypes.c_double * self.NPoints)()
        self.Result1 = (ctypes.c_double * self.NPoints)()
        self.Result2 = (ctypes.c_double * self.NPoints)()
        for i in range(self.NPoints):
            self.x_axis[i] = self.CenterFirstBin + float(i) * self.Stepsize

        self.epslib.SurfLossFunc(
            self.ParArray,
            self.Result1,
            ctypes.c_double(self.q),
            ctypes.c_int(self.DFChoice.get()),
        )

        self.SurfExProb = self.epslib.DSEP(
            self.ParArray, self.Result2, ctypes.c_int(self.DFChoice.get())
        )

    def PartDIIMFP(self):

        self.x_axis = (ctypes.c_double * self.NPoints)()
        self.Result1 = (ctypes.c_double * self.NPoints)()
        for i in range(self.NPoints):
            self.x_axis[i] = self.CenterFirstBin + float(i) * self.Stepsize

        #self.partialresults = (ctypes.c_double * (self.NPoints * 10))()
        self.partialresults =np.zeros(shape=(self.NPoints,10))
        #self.partcross = np.zeros(shape=(10,self.NPoints))
        #self.partcross = (ctypes.c_double * 10)()
        self.epslib.GOS_Scaling_init(self.ParArray)

        q_step = self.LastMomentum / 10.0
        for i in range(10):
            # here we override what was assigned in initParArray as is used by
            # standard diimfp calc, will be put back to normal values in fill_remainder
            self.ParArray[self.NDFPAR + 10] = i * q_step
            self.ParArray[self.NDFPAR + 11] = (i + 1) * q_step
            self.epslib.DIIMFP(
                self.ParArray,
                self.Result1,
                ctypes.c_int(self.DFChoice.get()),
                self.StoppingResultArray,
            )

            if self.weight == 0:
                self.PartIntSum[i] = self.StoppingResultArray[0]
                for j in range(self.NPoints):
                    self.partialresults[j, i] = self.Result1[j]
            else:
                self.PartIntSum[i] = self.StoppingResultArray[1]
                for j in range(self.NPoints):
                    self.partialresults[j, i] = (
                        self.Result1[j] * self.x_axis[j]
                    )
   

    def Penn_from_ELF(self):
        if self.N_oscillator_used.get() > self.maxoscillators:
            self.N_oscillator_used.set(self.maxoscillators)
        self.epslib.GOS_Scaling_init(self.ParArray)
        omega_p=self.w_p_Penn.get()
        if omega_p <=0:  
            omega_p = 1e8   # so the U correction never occur
            
        rangefactor=(self.LastEnergy)/self.CenterFirstBin
        exponent=1.0/self.N_oscillator_used.get()
        incrementfactor=rangefactor**exponent
    
        CurrentStepsize=self.CenterFirstBin*(incrementfactor-1)
        if  CurrentStepsize < 0.5*self.CenterFirstBin:
            CurrentEnergy = self.CenterFirstBin
        else:
            CurrentEnergy = 0.5* CurrentStepsize 
            rangefactor=(self.LastEnergy)/CurrentEnergy
            exponent=1.0/self.N_oscillator_used.get()
            incrementfactor=rangefactor**exponent
            
        self.initialized = False  # avoid check input message

  
        for i in range(self.N_oscillator_used.get()):
            Loss_at_E = self.epslib.lossfunction_inclGOS_atE(self.ParArray,
                ctypes.c_double(CurrentEnergy), ctypes.c_int(self.DFChoice.get()) )
            g_omega_prefactor = 2.0 / (cnst.PI * CurrentEnergy)  # Penn delta function  in the case  gamma=0
            # alternative from nguyen  for gamma !=0 (J of Phys Chem 119 23627 2015) does not work Sum A_i> 1
            # g_omega_prefactor=2.0/(cnst.PI*CurrentEnergy**2*CurrentGamma) * \
            # math.sqrt(2*CurrentEnergy*(CurrentEnergy**2+CurrentGamma**2)*(math.sqrt(CurrentEnergy**2+CurrentGamma**2)-CurrentEnergy))
            
            Amp = g_omega_prefactor * Loss_at_E * CurrentStepsize
            if CurrentEnergy  < omega_p:
                self.Amps[i].set(f"{Amp:.4g}")
                self.Omegas[i].set(f"{CurrentEnergy:.4g}")
                self.Gammas[i].set(f"{1.2*CurrentStepsize:.4g}")
                self.Alphas[i].set("1.0")
                self.Us[i].set("0.0")
            else:
                U=np.sqrt(CurrentEnergy*CurrentEnergy - omega_p*omega_p)
                NewA=Amp*CurrentEnergy*CurrentEnergy/(omega_p*omega_p)
                self.Amps[i].set(f"{NewA:.4g}")
                self.Omegas[i].set(f"{omega_p:.4g}")
                self.Gammas[i].set(f"{1.2*CurrentStepsize:.4g}")
                self.Alphas[i].set("1.0")
                self.Us[i].set(f"{U:.4g}")
                 
            CurrentStepsize *= incrementfactor
            CurrentEnergy = CurrentEnergy + CurrentStepsize
        self.DFChoice.set(value=3)
        for i in range(self.maxGOS):
            self.ConcGOS[i].set("0.0")
        for i in range(self.maxKaneko):
            self.N_Kaneko[i].set("0.0")
        for i in range(self.maxBelkacem):
            self.Conc_Belkacem[i].set("0.0")
        self.initialized = True
        
        
    def convert_DL_to_Kaneko(self): #convert a set of DL (or Mermin)  oscillators to a set of DL oscillators wilth the same optical ELF
        if self.DFChoice.get() == 1:
            print("requires DL or Mermin Oscillators, not Drude-Lindhard oscillators")
            return
        
    
        Q_used =  self.Q_Kaneko_transform.get()
        l_used = self.l_Kaneko_transform.get()
        print ("l used", l_used)
        if l_used == 0:
            doublefact = 1
        elif l_used == 1:
            doublefact = 3
        elif l_used == 2:
            doublefact = 5 * 3
        else:
            doublefact = 7 * 5 * 3
        w_pl_l=math.sqrt(2*(2*l_used+1)*doublefact*Q_used**3*math.exp(l_used)/((2.0*l_used)**l_used*math.sqrt(cnst.PI)))
 
        print("w_pl_l",w_pl_l* cnst.HARTREE)  
     

        for i in range(self.maxKaneko):
            A = float(self.Amps[i].get())
            if abs(A) > 0.0:
                W = float(self.Omegas[i].get()) / cnst.HARTREE
                U = float(self.Us[i].get()) /cnst.HARTREE
                elec_dens = ( A*W*W / (4.0 * cnst.PI) / (cnst.BOHR**3) )
                elec_per_unit_cell = elec_dens / self.UnitCellDensity
                #print("elec_per_unit_cell",elec_per_unit_cell)
                self.N_Kaneko[i].set(f"{elec_per_unit_cell:.5g}")
                self.Q_Kaneko[i].set(f"{Q_used:.4g}")
                self.Edge_Kaneko[i].set(f"{U*cnst.HARTREE:.4g}")
                self.l_Kaneko[i].set(f"{l_used:.4g}")
                width=float(self.Gammas[i].get())
                self.width_Kaneko[i].set(f"{width:.4g}")
                
                current_gamma= W**2 /w_pl_l**2
                self.gamma_Kaneko[i].set(f"{current_gamma:.4g}")  
       
        for i in range(self.maxoscillators):
            self.Amps[i].set("0.0")
        
                
    
  


    def convert_to_MLL(self):
  
        if self.DFChoice.get() == 1:
            print("requires DL or Mermin Oscillators, not Drude-Lindhard oscillators")
            return
        currentU= self.U_MLL_transform.get()             

        for i in range(self.maxoscillators):  # now calculate new values for optimum U
            oldA = float(self.Amps[i].get())
            if oldA > 0.0:
                old_omega = float(self.Omegas[i].get())
                old_U=float(self.Us[i].get())
                old_position=math.sqrt(old_omega**2 + old_U**2)
                if old_position < currentU:
                    print("oscillator", i, "is in the gap and is disgarded")
                    self.Amps[i].set(f"{0.0:.4g}")
                else:
                    new_omega = math.sqrt(old_position**2 - currentU**2)
                NewA = oldA * (old_omega / new_omega) ** 2
                self.Amps[i].set(f"{NewA:.4g}")
                self.Omegas[i].set(f"{new_omega:.5g}")
                self.Us[i].set(f"{currentU:.5g}")
     
        

    def Calc_Os_Strength(self):
        self.oneovereps2() 
        for i in range(self.NPoints):
            self.Result1[i] = (2.0 * self.x_axis[i] / (cnst.PI * self.PlasmonE**2) * self.Result2[i])
    

    
            
    def dyn_struct_factor(self):
        self.FineMeshFactor.set(9)
        self.eps1eps2_fine()
        for i in range(self.length_fine):
            denominator = (self.Result1_fine[i] * self.Result1_fine[i] + self.Result2_fine[i] * self.Result2_fine[i])
            self.Result1_fine[i] = self.Result1_fine[i] / denominator
            self.Result2_fine[i] = self.Result2_fine[i] / denominator  #this is Im[-1/eps],the loss function, NOT im[1/eps]
        elecdens = self.Nelec_per_UC *self.UnitCellDensity * cnst.BOHR**3  # elec/per a.u.^3
        omega_p_square=4.0 * cnst.PI * elecdens #Hartree square
       
        for i in range(self.length_fine):
            self.Result1_fine[i] =  self.q**2 / ( cnst.PI* omega_p_square)  * self.Result2_fine[i]
            self.Result1_fine[i] = self.Result1_fine[i] /cnst.HARTREE
        sum_so_far=0.0
        for i in range(self.length_fine):
            currentE=self.CenterFirstBin_fine+i*self.Stepsize_fine
            sum_so_far+=  self.Result1_fine[i] * currentE * self.Stepsize_fine
            self.Result2_fine[i]= sum_so_far
        self.recast()    
      
    
    def DL_from_OOS(self):
        if self.N_oscillator_used.get() > self.maxoscillators:
            self.N_oscillator_used.set(self.maxoscillators)
        self.PlasmonE= ( np.sqrt(4.0 * cnst.PI * self.UnitCellDensity * cnst.BOHR**3)
                * cnst.HARTREE )  # make sure we have the current one, w_p for one e- per u.c.
       
        omega= self.CenterFirstBin
        rangefactor=(self.LastEnergy)/omega
        exponent=1.0/self.N_oscillator_used.get()
        incrementfactor=rangefactor**exponent
    
        CurrentStepsize=omega*(incrementfactor-1)
        currentBin=0
        for i in range(self.N_oscillator_used.get()):
            while omega >  self.OOSEnergy[currentBin]:
                currentBin +=1
            PosWithinBin =  (omega - self.OOSEnergy[currentBin-1])/ (self.OOSEnergy[currentBin] - self.OOSEnergy[currentBin-1])
           # print("omega,self.OOSEnergy[currentBin-1],self.OOSEnergy[currentBin]",omega,self.OOSEnergy[currentBin-1],self.OOSEnergy[currentBin],"PosWithinBin ",PosWithinBin )
            currentOOS = self.OOS[currentBin-1] + PosWithinBin*(self.OOS[currentBin]-self.OOS[currentBin-1])
        
           
            currentELF=currentOOS/(2.0*omega)*cnst.PI *  self.PlasmonE**2
            g_omega_prefactor = 2.0 / (cnst.PI * omega) 
            Amp = g_omega_prefactor *currentELF * CurrentStepsize
            CurrentStepsize *= incrementfactor    
            self.Omegas[i].set(f"{omega:.4g}")
            self.Amps[i].set(f"{Amp:.5g}")
            self.Gammas[i].set(f"{1.3*CurrentStepsize:.4g}") 
            self.Alphas[i].set("1.0")
            self.Us[i].set("0.0")
            omega=omega + CurrentStepsize 
            if omega > self.OOSEnergy[self.N_OOS - 1]:
                    return

        
             
    def PseudoChargeDensity(self):
        """calculate the distribution of the charge density within the Penn Pseaudo charge picture"""

        # after this call  self.x_axis[] contains omega  and self.Result2[] the
        # loss function
        # sumAi = 0.0
        # for i in range(self.maxoscillators):
            # sumAi += float(self.Amps[i].get())
        self.oneovereps2()
        Energy_step_au = self.Stepsize / cnst.HARTREE  # energy step in a.u.
        vol_so_far = 0.0
        charge_so_far=0.0
        i_overfull = 0.0
        for i in range(self.NPoints):

            # now energy axis in atomic units
            self.x_axis[i] = self.x_axis[i] / cnst.HARTREE
        for i in range(self.NPoints):
            j = self.NPoints - i - 1
            # G(omega from Penn eq. 9b) within the statistical interpretation
            # this is also proportional to the fraction of space filled with
            # this electron density
            GPenn = self.Result2[j] * 2 / (cnst.PI * self.x_axis[j])
            # e density of this electron loss (a.u.)^-3
            elec_dens = self.x_axis[j] ** 2 / (4 * cnst.PI)
            # electron density per angstrom^3
            self.Result1[j] = elec_dens / cnst.BOHR**3
            vol_so_far += GPenn * Energy_step_au
            charge_so_far+=elec_dens* GPenn * Energy_step_au/self.UnitCellDensity / cnst.BOHR**3
            self.x_axis[j] = vol_so_far
            if self.x_axis[j] > 1.0:
                i_overfull += 1

                self.x_axis[j] = 1.0

        self.x_axis[0] = 1.0
        self.Result1[0] = 0.0

    def ConvertToRadialPseudoChargeDensity(self):
        # should be called straight after  calc.PseudoChargeDensity(), only makes sense for pure elements
        self.r_array = (ctypes.c_double * self.NPoints)()
        self.dens_array = (ctypes.c_double * self.NPoints)()
        UC_Volume = 1 / self.UnitCellDensity
        for i in range(self.NPoints):
            current_volume = self.x_axis[i] * UC_Volume
            radius = (current_volume * 3.0 / (4.0 * cnst.PI)) ** (1.0 / 3.0)
            self.r_array[i] = radius
            self.dens_array[i] = self.Result1[i]
            self.x_axis[i] = self.r_array[i]  # so we can plot it by itself
            
            
    def stopping_IMFP_w_p_versus_r(self):
        self.Result3 = (ctypes.c_double * self.NPoints)() # for stoppings result1 for w_p, result2 for lambda
        self.Result4 = (ctypes.c_double * self.NPoints)() # for straggling result1 for w_p, result2 for lambda
        if self.particle == "proton":
            ConstantA = 1.0/2.0  # notation so we are in line with imfp draft eq. 11
            ConstantB = 1
        else:
            ConstantA = 1
            ConstantB = 2
        velocity = self.velocity_projectile(self.E0)  # in a.u.   
        q_max = velocity / ConstantA #non-relativistic case have to think about this
        
        for i in range(self.NPoints):
            current_density = self.Result1[i]
            if(current_density > 0.0):
                PlasmonE= np.sqrt(4.0 * cnst.PI * current_density * cnst.BOHR**3)
                if velocity**2 > 2 * ConstantB * PlasmonE:
                    q_min_DL = (velocity - np.sqrt(velocity**2 - 2 * ConstantB * PlasmonE)
                    ) / ConstantB
                    q_max_DL = (velocity + np.sqrt(velocity**2 - 2 * ConstantB * PlasmonE)
                    ) / ConstantB
                else:
                    q_min_DL = 1
                    q_max_DL = 1
                logDL = np.log(q_max_DL / q_min_DL)    
                tmp = logDL - 0.5 * math.log(
                    (2 * PlasmonE + q_max_DL**2) / (2 * PlasmonE + q_min_DL**2))
                oneoverlambda = PlasmonE * tmp / velocity**2    
                if oneoverlambda > 0.0:
                    Lambda = (1.0 / oneoverlambda)  # in a.u.
                else:
                    Lambda = 0.0 
                stopping =   PlasmonE**2 / velocity**2 * logDL
                tmp =(PlasmonE**3 / velocity**2 * logDL) 
                straggling = 0.25 * PlasmonE**2 / velocity**2 * (q_max_DL**2 - q_min_DL**2) + tmp 
                    
         
    
                    
            else:
                Lambda = 0.0 
                PlasmonE = 0.0 
                stopping=0.0
                straggling = 0.0      
            

            self.Result1[i] = PlasmonE * cnst.HARTREE
            self.Result2[i] = Lambda  * cnst.BOHR
            self.Result3[i] = stopping  * cnst.HARTREE /cnst.BOHR
            self.Result4[i] = straggling  * cnst.HARTREE * cnst.HARTREE / cnst.BOHR
   

    def stopping_versus_r(self):
        """IMFPversus r, to be called straight after  ConvertToRadialPseudoChargeDensity()
        atom in cube with volume Muffin_Tin sphere, i.e. unit cell volume"""
        Cube_length = (1.0 / self.UnitCellDensity) ** (1.0 / 3.0)
        self.MT_radius = (1.0 / self.UnitCellDensity * 3.0 / (4.0 * cnst.PI)) ** (
            1.0 / 3.0
        )
        print("MT radius", self.MT_radius, "side cube", Cube_length)
        # prepare the array for stopping power calculation with one Mermin

        # first store the calculated density versus r values safely
        currentEnergy = self.CenterFirstBin
        currentstepsize = self.Stepsize
        pointswithfixedstepsize = int(
            (self.LastEnergy - self.CenterFirstBin) / self.Stepsize
        )

        for ipoints in range( pointswithfixedstepsize):  # this should match what happens in epslib
  
            if pointswithfixedstepsize > self.max_no_lin_cont_deltaE:
                self.fixedstepsize = False
                self.ParArray[self.NDFPAR + 9] = self.lin_cont_deltaE
                currentstepsize = self.Stepsize + self.lin_cont_deltaE * currentEnergy

            else:
                self.fixedstepsize = True
                self.ParArray[self.NDFPAR + 9] = 0.0

            currentEnergy += currentstepsize
            if currentEnergy > self.LastEnergy:
                break

        diimfparray = (
            ctypes.c_double * self.NPoints
        )()  # calculated but not really used

        self.stop_array = (ctypes.c_double * self.NPoints)()

        for i in range(self.maxoscillators):  # make sure other oscillators in c routine have amplitude 0 
            self.ParArray[4 * i + 1] = 0.0
        for i in range(self.maxGOS):  # make sure all GOS amp are 0
            self.ParArray[4 * i + 4 * self.maxoscillators + 1] = 0.0
        for i in range(self.maxKaneko):  # make sure all GOS amp are 0
            self.ParArray[
                4 * i + 4 * self.maxoscillators + 4 * self.maxKaneko + 1
            ] = 0.0

        # make sure remainder of pararray is filled with correct data
        self.fill_remainder()

        self.ParArray[self.NDFPAR + 5] = (
            ipoints + 1
        )  # put to NPoints by fill remainder but should be ipoints+1 (variable step size)
        for i in range(self.NPoints):
            elec_dens_au = self.dens_array[i] * cnst.BOHR**3
            w_pl = cnst.HARTREE * math.sqrt(4 * cnst.PI * elec_dens_au)  # in eV
            # at distance r_array[i] we have stopping due to DF with plasmon energy w_pl
            self.ParArray[1] = 1.0
            self.ParArray[2] = w_pl
            self.ParArray[3] = (w_pl / 10)  # this is gamma, should not matter may afffect comp. time required
            self.epslib.DIIMFP(
                self.ParArray,
                diimfparray,
                ctypes.c_int(self.DFChoice.get()),
                self.StoppingResultArray,
            )

            self.stop_array[i] = self.StoppingResultArray[1]

    def impact_dep_stop(self): #currently not connected to user interface, newver called
        """impact dependent stopping, to be called straight after  ConvertToRadialPseudoChargeDensity()
        atom in cube with volume Muffin_Tin sphere, i.e. unit cell volume.  It is assumed that the atoms
        are in a simple cubic lattice with nearest neighbor distance Cube_length
        Maximum impact parameter considered is Cube_length/2. Particle impinges perpendicular to cube
        The energy loss for crossing the cube is calculated as a function of the impact parameter.
        Part of the cube is outside the Muffin Tin sphere.  For that part the electron density at
        the edge of the MT shpere is taken.
        This edge density again is determined by the lowest energy considered for the ELF. Take this not too low"""
        # calculate charge density distribution (Penn's pseudo charge density")
        self.PseudoChargeDensity()
        # put it in a radial form
        self.ConvertToRadialPseudoChargeDensity()
        # calculate sthe local stopping for that density (this takes time)
        self.stopping_versus_r()
        # calculate inpact parameter dependent stopping

        Cube_length = (1.0 / self.UnitCellDensity) ** (1.0 / 3.0)
        print("MT radius", self.MT_radius, "side cube", Cube_length)
        Bmax = math.sqrt(self.MT_radius**2 - (Cube_length / 2) ** 2)
        print(Bmax, "BMax for whole trajectory within MT sphere")
        self.Nst = 500
        sumB2 = 0
        sumstop = 0

        self.b_dep_stop = (ctypes.c_double * self.Nst)()
        self.b_param = (ctypes.c_double * self.Nst)()
        steplength = Cube_length / (2 * self.Nst)
        for i in range(self.Nst):
            self.b_dep_stop[i] = 0.0
            self.b_param[i] = (
                i + 0.5
            ) * steplength  # so start at b=Cube_length/(2*Nst) up to almost Cube_length/2
            sumB2 += self.b_param[i] * self.b_param[i]
            for j in range(self.Nst):
                along = (
                    j + 0.5
                ) * steplength  # so from almost  central atom up to cube edge

                r = math.sqrt(along * along + self.b_param[i] * self.b_param[i])
                if r > self.MT_radius:
                    r = self.MT_radius - (
                        r - self.MT_radius
                    )  # folding back when outside mt radius
                for current_index in range(self.NPoints):
                    if r > self.r_array[current_index]:
                        break  # so if we are never outside MT sphere
                current_index = (
                    current_index - 1
                )  # this means if part crystal is empty this part get stopping of zero not the next value
                self.b_dep_stop[i] += (
                    2 * self.stop_array[current_index] * steplength / Cube_length
                )  # the factor of 2 because we integrate only from edge to center
                # not to the other edge
            sumstop += self.b_param[i] * self.b_param[i] * self.b_dep_stop[i]
  
        print("estimate of average stopping:", sumstop / sumB2)
        # now we are going to calculate the contribution of bunching to straggling
        self.straggling_from_bunching = 0.0
        self.straggling_from_bunching_rel_vacuum = 0.0
        for i in range(self.Nst):
            self.straggling_from_bunching += (
                self.b_dep_stop[i] - sumstop / sumB2
            ) ** 2 * self.b_param[i] ** 2
            self.straggling_from_bunching_rel_vacuum += (
                self.b_dep_stop[i]
            ) ** 2 * self.b_param[i] ** 2
        self.straggling_from_bunching = self.straggling_from_bunching / sumB2
        self.straggling_from_bunching_rel_vacuum = (
            self.straggling_from_bunching_rel_vacuum / sumB2
        )
        print(self.straggling_from_bunching, "bunching per unit cell")
        self.straggling_from_bunching = self.straggling_from_bunching / Cube_length
        self.straggling_from_bunching_rel_vacuum = (
            self.straggling_from_bunching_rel_vacuum / Cube_length
        )
        print(self.straggling_from_bunching, "bunching angstrom")
        print(
            self.straggling_from_bunching / self.BohrStraggling,
            "bunching rel to Bohr straggling",
        )
        print(
            self.straggling_from_bunching_rel_vacuum
            / Cube_length
            / self.BohrStraggling,
            "bunching rel to Bohr straggling_rel_vacuum",
        )
        print("self.BohrStraggling", self.BohrStraggling)

    def bunching_versus_E0(self):
        # calculate charge density distribution (Penn's pseudo charge density")
        self.PseudoChargeDensity()
        # put it in a radial form
        self.ConvertToRadialPseudoChargeDensity()
        oldE0 = self.E0
       
        if self.particle  == "proton":
            self.E0 = 2
        else:
            self.E0 = self.first_electron_energy

        StragglingResult1 = (ctypes.c_double * self.NStopping)()
        StragglingResult2 = (ctypes.c_double * self.NStopping)()
        Beam_incr = 1.25
        NBunch = 35
        for Ecounter in range(NBunch):
            print("currently calculating No", Ecounter, "Energy", self.E0)
            self.ParArray[self.NDFPAR + 1] = self.E0 * 1000
            self.impact_dep_stop()
            StragglingResult1[Ecounter] = (
                self.straggling_from_bunching / self.BohrStraggling
            )
            StragglingResult2[Ecounter] = (
                self.straggling_from_bunching_rel_vacuum / self.BohrStraggling
            )
            self.E0 = self.E0 * Beam_incr

        self.E0 = oldE0
        self.x_axis = (ctypes.c_double * NBunch)()
        self.Result1 = (ctypes.c_double * NBunch)()
        self.Result2 = (ctypes.c_double * NBunch)()
        E = self.first_proton_energy
        self.NPoints = NBunch
        for i in range(NBunch):
            self.x_axis[i] = E
            self.Result1[i] = StragglingResult1[i]
            self.Result2[i] = StragglingResult2[i]
            E = E * Beam_incr



    def REELS_spectrum(self):
        self.Start_REELS = -5.0  # hard coded to start at -5
        self.NREELS = (
            int((self.LastEnergy - self.Start_REELS) / self.Stepsize) + 1
        )  # hard coded to start at -5
        self.NDIIMFP = int(self.LastEnergy / self.Stepsize) + 1
        self.x_axis = (ctypes.c_double * self.NREELS)()
        self.Result1 = (ctypes.c_double * self.NREELS)()
        self.NormDIIMFP = (ctypes.c_double * self.NDIIMFP)()
        self.epslib.GOS_Scaling_init(self.ParArray)
        for Ecounter in range(self.NREELS):
            self.x_axis[Ecounter] = -5.0 + Ecounter * self.Stepsize
            self.Result1[Ecounter] = 0.0
        self.epslib.calc_REELS(
            self.ParArray,
            ctypes.c_double(self.Start_REELS),
            ctypes.c_int(self.NREELS),
            ctypes.c_int(self.DFChoice.get()),
            self.Result1,
            ctypes.c_int(self.EELS.get()),
            ctypes.c_double(self.EELS_thickness.get())
            
        )
    def calc_Qrecoil(self, q):  #calculates (relativistically ) the kinetic energy of an initially stationary electron after it's momentum changes to q
        #see appendix J.M. Fernández-Varea et al. / Nucl. Instr. and Meth. in Phys. Res. B 229 (2005) 187–218
        # #high q approach
        # Q= math.sqrt(cnst.C**2*q**2+cnst.C**4)-cnst.C**2
        # #low q approach
        # x=0.5*(cnst.C*q/cnst.C**2)**2
        # Q_low=cnst.C**2*(x-x**2/2+x**3/2)
        # print ("q=",q,"Q high q approach",Q, "Q, low q approach", Q_low, "percentage", 100* (Q_low-Q)/Q, "q^2/2",q*q/2)
        #conclusion of above test: use q^2/2 up to 3 a.u. high-Q approach above that
        if q < 3.0:
            return q**2/2
        else:    
            return math.sqrt(cnst.C**2*q**2+cnst.C**4)-cnst.C**2
        
         
        
        
    def updatedensity(self, var, indx, mode):
        if self.initialized == False:
            return
        
        self.my_updatedensity()   
           
    def my_updatedensity(self):  
        if  not self.initialized:
            return
    
        self.SumAi = 0.0
        SumGOS = 0.0
        SumKaneko = 0.0
        VolumeFractionKaneko=0.0
        SumBelkacem = 0.0
        elec_dens = 0.0
        if self.massunitcell.get() <= 0.0:
            return
        mole_per_cm3 = self.specificweight.get() / self.massunitcell.get()
       
        self.UnitCellDensity = mole_per_cm3 * cnst.NAvogadro / 1.0e24
        self.UnitCellDensityText.set(
            "unit cell  (u.c.) density: "
            + f"{self.UnitCellDensity:.4f}"
            + " per Å³"
        )
        self.PlasmonE= ( np.sqrt(4.0 * cnst.PI * self.UnitCellDensity * cnst.BOHR**3) * cnst.HARTREE)
        self.plasmon_1elec_per_unit_cell.set(
            "ωₚ  for  1 e⁻ per u.c.:  " + f"{self.PlasmonE:.2f}" + " eV"
        )

        self.MT_radius = (1.0 / self.UnitCellDensity * 3.0 / (4.0 * cnst.PI)) ** (1.0 / 3.0)
        try:
            for i in range(self.maxGOS):
                SumGOS += float(self.ConcGOS[i].get())
            for i in range(self.maxBelkacem):
                SumBelkacem += float(self.Conc_Belkacem[i].get())
            for i in range(self.maxKaneko):
                
                if float(self.N_Kaneko[i].get()) > 0.0:
                    currentN = float(self.N_Kaneko[i].get())
                    currentQ = float(self.Q_Kaneko[i].get())
                    currentG = float(self.gamma_Kaneko[i].get())
                    currentl = float(self.l_Kaneko[i].get())
                    SumKaneko += currentN
                    if self.Kaneko_choice.get() == 0:   #modified Kaneko
                        print("modified")
                        if currentl==0:
                            doublefact=1
                        elif currentl==1: 
                            doublefact=3
                        elif currentl==2:
                            doublefact=15
                        elif (currentl==3):
                            doublefact=105
                        #w_pl_l[i]=sqrt(gamma_Kaneko[i]*2*(2*l+1)*doublefact*pow(Q,3)*exp(l)/(pow(2.0*l,l)*sqrt(pi)) )
                        w_p=math.sqrt(currentG*2*(2*currentl+1)*doublefact*currentQ**3 \
                            *math.exp(currentl)/((2.0*currentl)**currentl*math.sqrt(math.pi)) )
                        w_p=w_p*cnst.HARTREE   
                        VolumeFractionKaneko += currentN*self.PlasmonE*self.PlasmonE/(w_p*w_p)  
                        print ("w_p", w_p)
                    else:   # original Kaneko
                        qmean= currentQ*currentN**(1.0/3.0) 
                       # tmp=currentG*qmean**3/math.sqrt(math.pi)
                        w_p = math.sqrt(currentG*qmean**3 / math.sqrt(math.pi) ) * cnst.HARTREE
                        VolumeFractionKaneko += currentN*self.PlasmonE*self.PlasmonE/(w_p*w_p)   
                         
            if self.DFChoice.get() == 1:
                for i in range(self.maxoscillators):
                    A = float(self.Amps[i].get())
                    elec_dens +=  A / cnst.HARTREE**2 / (cnst.BOHR**3) / (4.0 * cnst.PI)
                elec_per_unit_cell = elec_dens / self.UnitCellDensity
            else:
                for i in range(self.maxoscillators):
                    A = float(self.Amps[i].get())
                    if abs(A) > 1e-50:
                        W = float(self.Omegas[i].get()) / cnst.HARTREE
                        U = float(self.Us[i].get()) /cnst.HARTREE
                        
                        if W > 0.0:
                            W_M_square = W**2 + U**2
                            A_M = W*W / W_M_square * A
                        else:
                            A_M=0.0  
                            W_M_square=0.0
                        self.SumAi += A_M
                        # print("A,SumAi", A, SumAi)
                        
                        

                        elec_dens += ( A_M * W_M_square / (4.0 * cnst.PI) / (cnst.BOHR**3) )
        except  ValueError: 
            return
        except Exception as err:
            print(f"Unexpected {err=}, {type(err)=}")
            raise     
        elec_per_unit_cell = elec_dens / self.UnitCellDensity
        OStext= "Oscillatos:  "+  f"{elec_dens:.2f}"  + " e⁻/Å³ or "  + f"{elec_per_unit_cell:.2f}" + " e⁻ per U.C."
        if self.DFChoice.get() != 1:
            OStext  = OStext + "   (\u03A3 c(U)Aₓ ="  + f"{self.SumAi:.3f})" 
           
        self.DF_prop_text.set(OStext)

        density = SumGOS * self.UnitCellDensity
        self.GOSdenstext.set(
            "GOS:            "
            + f"{ density:.3f}"
            + " e⁻/Å³ or "
            + f"{SumGOS:.2f}"
            + " e⁻ per U.C."
        )
        density = SumBelkacem * self.UnitCellDensity
        self.Belkacemtext.set(
            "BelKacem: "
            + f"{ density:.3f}"
            + " e⁻/Å³ or "
            + f"{SumBelkacem:.2f}"
            + " e⁻ per U.C."
        )
        density = SumKaneko * self.UnitCellDensity
        self.Kanekodenstext.set(
            "Kaneko:     "
            + f"{ density:.3f}"
            + " e⁻/Å³ or "
            + f"{SumKaneko:.2f}"
            + " e⁻ per U.C."
            + " volume fraction: " 
            + f"{VolumeFractionKaneko:.3f}"
        )
        self.Nelec_per_UC = elec_dens / self.UnitCellDensity + (SumGOS + SumKaneko + SumBelkacem) 
        self.BohrStraggling = (4 * cnst.PI * self.Nelec_per_UC *self.UnitCellDensity * cnst.BOHR**3 * cnst.HARTREE**2 / cnst.BOHR)


        mytext = ("from eV/ Å to MeV/(mg/cm²):           " + f"{0.1/self.specificweight.get():.3f}")
        self.myconversiontext1.set(mytext)
        mytext = ("from eV/ Å to eV/(10¹⁵ molecules/cm²): "    + f"{self.massunitcell.get()/(self.specificweight.get()*6.022):.3f}")

        
        self.myconversiontext2.set(mytext)
        self.MyChapApp.UpdateColors() 
         
        
    def updateProjectileEnergy(self, var, indx, mode): 
        # this procedure is called  by TK when particle and or energy of projectile has changed (on idle)
        try:
            self.my_updateProjectileEnergy()
        except Exception:
            print("please, check beam energy")
            
    def my_updateProjectileEnergy(self):  #this can be called both by the tk trace facility and from the program      
        self.E0 = self.ProjectileEnergy.get()
        if self.Projectile.get() == 1:
            self.ProjectileMass = cnst.Mp
            self.particle="proton"
        else:
            self.ProjectileMass = 1.0
            self.particle= "electron"
        v=self.velocity_projectile(self.E0)
        self.myvelocitytext.set("p=" + f"{self.p0_rel:.1f}"+"a.u.   v/c ="+f"{self.beta_r:.3f}")         
    
    def velocity_projectile(self, Energy):  # retuns velocity, also sets self.p0_rel, self.beta_rel and self.gamma_r
       E0_au= Energy*1000.0/cnst.HARTREE 
       self.gamma_r = 1.0 + E0_au/(self.ProjectileMass*(cnst.C)**2)
       T = E0_au * (1.0 + self.gamma_r) / (2.0 * self.gamma_r * self.gamma_r)  # Egerton appendix E, T in a.u.
       v = np.sqrt(2.0 * T/ self.ProjectileMass)   # V in a.u.
       self.p0_rel = self.gamma_r * self.ProjectileMass * v
       self.beta_r= v/cnst.C
       return v
                  
      
