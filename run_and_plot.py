import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.ticker import FormatStrFormatter
from tkinter import StringVar, BooleanVar, IntVar, DoubleVar
import numpy as np
import constants as cnst

class run_and_plot:
   # def __init__(self, *args, **kwargs):
    def __init__(self,parent):
        self.MyChapApp=parent
        self.calc = self.MyChapApp.calc
        self.PlotDescription = StringVar(value="")
        self.LiteratureDescription = StringVar(value = "lit. data")
       # plt.rc('font', size=12)
        self.DisplayGrid = BooleanVar(value = False)
        self.CompDataLine = BooleanVar(value = True)
        self.legendfontsize = IntVar(value=12)
        self.axislabelsize = IntVar(value=14)
        self.linewidth  = DoubleVar(value=2)
        self.markersize = IntVar(value=4)
        self.xlabelsize = IntVar(value=12)
        self.ylabelsize = IntVar(value=12)
        self.title_fontsize = IntVar(value=12)
        self.fileformat= StringVar(value="pdf")
        self.screen_width=0
        self.screen_height=0
        self.figurewidth= IntVar(value=540)
        self.figureheight= IntVar(value=400)
        self.legendposition= IntVar(value=0)
        self.start_x, self.start_y, self.dx, self.dy = (0, 0, self.figurewidth.get(),  self.figureheight.get())
 
#===================================================================
    def start_calc_and_plot(self):
        if self.initcalc() != 0:
            print("problems")
            return
        #try:
        routine = getattr(self, self.MyChapApp.plotchoice.get())
        routine()
        
        
        
    def replot(self):
       self.update_plot_settings() 
       plotchoice=self.MyChapApp.plotchoice.get()
       self.calc.max_eq = float( self.MyChapApp.Max_eq.get())
       self.calc.log_range = float( self.MyChapApp.log_range.get())
       if plotchoice == "IMFP_stop_strag":
            if self.calc.axis_choice.get() == 1:
                self.xlabel = self.particle + " energy (keV)"
            else:
                self.xlabel = self.particle +" velocity (a.u.)"
            if self.calc.StoppingUnits.get()==0:
                self.label_curve2 = "stopping (eV/Å) "
                self.calc.StoppingFactor=1.0
            elif self.calc.StoppingUnits.get()==1:
                self.label_curve2 = "stopping (eV / (1E¹⁵ atoms/cm²)) "
                self.calc.StoppingFactor=self.calc.massunitcell.get()/(self.calc.specificweight.get()*6.022)
                print("was here")
            elif self.calc.StoppingUnits.get()==2:
                self.label_curve2 = "stopping (MeV / (mg/cm²))"    
                self.calc.StoppingFactor=0.1/self.calc.specificweight.get()
            self.plot_result3()
       elif plotchoice == "eq_plot"  or  plotchoice == "self_eq_plot"  or plotchoice == "dcs_omega_eq_plot":
           self.calc.scale_image()
           self.colorplot_result() 
       elif plotchoice == "partial_DIIMFP" or plotchoice == "partial_stopping":   
           self.plot_result10()
       elif plotchoice == "diimfp" or plotchoice == "Compton" or plotchoice == "oscillator_strength" \
           or plotchoice == "dyn_struct_factor" or plotchoice == "pseudo_charge_density" or plotchoice == "REELS"\
           or plotchoice =="radial_charge_density" or plotchoice == "DDCS_at_omega" or plotchoice== "bin_enc_peak":
           self.plot_result1()   
       elif plotchoice == "eps_w_kk" or plotchoice == "one_over_eps_w_kk" \
           or plotchoice == "radial_w_p_lambda" or plotchoice == "shell_effect_all":
           self.plot_result4()
       elif plotchoice == "SELF_DSEP" or plotchoice == "S_k_omega_rule":
           self.plot_result2_dif_y_axes()
       elif plotchoice == "sum_rules":
           self.plot_resultSumrules() 
       else:
           self.plot_result2()
                      
            
            
       
#-----------------------------------------------------------------------    
    def initcalc(self):
        
        self.calc.LastEnergy = float(self.MyChapApp.LE.get())
        self.calc.Stepsize = float(self.MyChapApp.step.get())
        self.MyChapApp.calc.LowerELimit =  float(self.MyChapApp.FE.get())
        self.MyChapApp.calc.CenterFirstBin = self.MyChapApp.calc.LowerELimit + 0.5*self.calc.Stepsize
        self.calc.q = float(self.MyChapApp.qused.get())
        self.calc.NPoints = (int((self.MyChapApp.calc.LastEnergy - self.calc.CenterFirstBin) / self.calc.Stepsize) + 1)
   
        # npoints get sometimes changed by calc, number of points to be plotted
        self.calc.Energy_qplot = float(self.MyChapApp.omega_used.get())

        self.calc.Stepsize_qplot = float(self.MyChapApp.stepq.get())
        self.calc.LastMomentum = float(self.MyChapApp.LM.get())

        self.calc.q_Compton= float(self.MyChapApp.Compton_q.get())
        self.calc.omega_ddcs= float(self.MyChapApp.omega_ddcs_entree.get())
        self.calc.theta_max= float(self.MyChapApp.theta_max_entree.get())
        self.calc.theta_ddcs= float(self.MyChapApp.theta_ddcs_entree.get())
        
      #  self.calc.E0 = float(self.E0entry.get())

        self.calc.max_eq = float( self.MyChapApp.Max_eq.get())
        self.calc.log_range = float( self.MyChapApp.log_range.get())
        # so we know what was the last calculation when we save
        self.last_calculation = self.MyChapApp.plotchoice.get()
        self.ylabel = ""
        self.plottitle = ""
        self.label_curve1 = "--"
        self.label_curve2 = "--"
        self.label_curve3 = "--"
        self.label_curve4 = "--"
        if self.calc.Projectile.get() == 1:
            self.particle="proton"
            self.particle_LaTeX =r"H$^{+}$"
        else:
            self.particle="electron" 
            self.particle_LaTeX =r"e$^{-}$"  
        
        self.calc.NStopping = int(self.MyChapApp.NStopping_entree.get()) 
        self.calc.IncrFactor =float(self.MyChapApp.IncrFactor_entree.get())  
        self.calc.first_electron_energy =float(self.MyChapApp.first_E_e_entree.get())     
        self.calc.first_proton_energy =float(self.MyChapApp.first_E_H_entree.get())   
        self.calc.q_transition =float(self.MyChapApp.q_transition_entree.get())  
        self.calc.UMax =float(self.MyChapApp.UMax_entree.get()) 
        self.calc.BE_for_exchange=float(self.MyChapApp.Be_for_exchangeEntry.get()) 
        self.update_plot_settings()
        self.MyChapApp.UpdateColors()
        error_code=self.calc.initParArray()
        return error_code
        
        
    def update_plot_settings(self):    
        plt.rcParams["figure.figsize"] = (self.figurewidth.get()/100.0,self.figureheight.get()/100.0) 
        self.dx=self.figurewidth.get()
        self.dy= self.figureheight.get()
        theme = {'axes.grid':  self.DisplayGrid.get(),
             'grid.linestyle': '--',
             'legend.framealpha': 1,
             'legend.facecolor': 'white',
             'legend.shadow': False,
             'legend.fontsize':self.legendfontsize.get(),
             'legend.title_fontsize': 14,
             'xtick.labelsize': self.xlabelsize.get(),
             'ytick.labelsize': self.ylabelsize.get(),
             'axes.labelsize': self.axislabelsize.get(),
             'axes.titlesize': 20,
             'lines.linewidth': self.linewidth.get(),
             'lines.markersize': self.markersize.get(),
             'figure.dpi': 100,
             'savefig.format':self.fileformat.get()}
    
        matplotlib.rcParams.update(theme) 
        if self.CompDataLine.get():
            self.overplotline='dashed'
        else:
            self.overplotline='none'

    # ================================== start calculations=======================================
    def eps_w_direct(self):
        self.calc.eps1eps2()
        self.xlabel = "ω (eV)"
        self.label_curve1 = "Re [eps]"
        self.label_curve2 = "Im [eps]"
        self.plottitle = "evaluated at q=" + str(self.calc.q) + "a.u."
        self.plot_result2()
        
    def one_over_eps_w_direct(self):  
        self.calc.oneovereps1eps2()
        self.xlabel = "ω (eV)"
        self.label_curve1 = "Re [1/eps]"
        self.label_curve2 = "Im [-1/eps] "
        self.plottitle = "evaluated at q=" + str(self.calc.q) + "a.u."
        self.plot_result2()  
        
    def eps_q_direct(self): 
        self.calc.eps1eps2_q()
        self.xlabel = "q (a.u.)"
        self.label_curve1 = "Re [eps(q)]"
        self.label_curve2 = "Im [eps(q)]"
        self.ylabel = ""
        self.plottitle = "evaluated at ω=" + str(self.calc.Energy_qplot) + "eV"
        self.plot_result2()     
        
    def one_over_eps_q_direct(self):   
        self.calc.oneovereps1eps2_q()
        self.xlabel = "q (a.u.)"
        self.label_curve1 = "Re [1/eps(q)]"
        self.label_curve2 = "Im [-1/eps(q)]"
        self.ylabel = ""
        self.plottitle = "evaluated at ω=" + str(self.calc.Energy_qplot) + "eV"
        self.plot_result2() 
        
    def eps_w_kk(self):   
        self.calc.eps_kk_test()
        self.xlabel = r"$\omega$ (eV)"
        self.label_curve1 = r"Re[$\epsilon$]"
        self.label_curve2 = r"Im [$\epsilon$]"
        self.label_curve3 = r"$1 + \frac{2}{\pi} {\cal P}\int_0^\infty\,\, \frac{\omega'\, {\rm Im} \left[\epsilon (\omega',\, q)\right]}{(\omega')^2-\omega^2}  d\omega' $"
        self.label_curve4 = r"$ -\frac{2 \omega}{\pi}  {\cal P}\int_0^\infty\,\, \frac{ {\rm Re} \left[\epsilon (\omega',\, q)\right]-1}{(\omega')^2-\omega^2}  d\omega'$"
        self.plottitle = "evaluated at q=" + str(self.calc.q) + "a.u."
        self.plot_result4()   
        
    def one_over_eps_w_kk(self):   
        self.calc.one_over_eps_kk_test()
        self.xlabel = r"$\omega$ (eV)"
        self.label_curve1 = r"Re[1/$\epsilon$]"
        self.label_curve2 = r"Im [-1/$\epsilon$]"
        self.label_curve3 = r"$1 + \frac{2}{\pi} {\cal P}\int_0^\infty\,\, \frac{\omega'\, {\rm Im} \left[\frac{1}{\epsilon (\omega',\, q)}\right]}{(\omega')^2-\omega^2}\,  d\omega' $"
        self.label_curve4 = r"$ -\frac{2 \omega}{\pi}  {\cal P}\int_0^\infty \,\, \frac{ {\rm Re} \left[\frac{1}{\epsilon (\omega',\, q)}\right]-1}{(\omega')^2-\omega^2}\,  d\omega'$"
        self.plottitle = "evaluated at q=" + str(self.calc.q) + "a.u."
        self.plot_result4() 
        
    def n_k_kk(self):   
        self.calc.n_and_k_kk_test()
        self.xlabel = r"$\omega$ (eV)"
        self.label_curve1 = r"$n$"
        self.label_curve2 = r"$k$"
        self.label_curve3 = r"$1 + \frac{2}{\pi} {\cal P}\int_0^\infty\,\, \frac{\omega'\, k (\omega',\, q)}{(\omega')^2-\omega^2}\,  d\omega' $"
        self.label_curve4 = r"$ -\frac{2 \omega}{\pi}  {\cal P}\int_0^\infty \,\, \frac{ n(\omega',\, q)-1}{(\omega')^2-\omega^2}\,  d\omega'$"
        self.plottitle = "evaluated at q=" + str(self.calc.q) + "a.u."
        self.plot_result4()    
 
            
    def n_k(self): 
        self.calc.n_and_k_from_eps1_eps2()
        self.xlabel = "ω (eV)"
        self.label_curve1 = "n "
        self.label_curve2 = "k "
        self.plottitle = "evaluated at q=" + str(self.calc.q) + "a.u."
        self.plot_result2()   
        
    def sum_rules(self):  
        self.calc.sum_rules()
        self.xlabel = "ω (eV)"
        self.label_curve1 = (r"$\frac{2}{\pi \Omega_p^2}\int_0^\omega \omega \, {\rm Im} [\epsilon (q,\omega)] d\omega$")
        self.label_curve2 = (r"$\frac{2}{\pi \Omega_p^2}\int_0^\omega \omega \, {\rm Im} [-1/\epsilon (q,\omega)] d\omega$")
        self.label_curve3 = r"$\frac{4}{\pi \Omega_p^2} \int_0^\omega \omega \, k(q,\omega) d\omega$"
        self.label_curve4 = r"$\frac{2}{\pi }\int_0^\omega \frac{1}{\omega}$" + \
            r"$\, {\rm Im} [-1/\epsilon (q,\omega)] d\omega$"
        self.ylabel = "electrons per unit cell "
        self.plottitle = "evaluated at q=" + str(self.calc.q) + "a.u."
        self.plot_resultSumrules()  
        
        
    def inertial_rules(self):  
        self.calc.inertial_rules()
        self.xlabel = "ω (eV)"
        self.label_curve1 = r"$\int_0^\omega  (Re \left[\frac{1}{\epsilon (\omega')}\right]-1) d\omega' $"
        self.label_curve2 = r"$\int_0^\omega  (n(\omega')-1) d\omega' $"
        self.label_curve4 = r"$\frac{-1}{2 \pi^2}\int_0^\omega  (\epsilon_1(\omega')-1) d\omega' $"
        self.label_curve3 = r"$\lim_{\omega \to 0}\,\,  \frac{\omega}{4 \pi}\,  \epsilon_2 (\omega,q) $" 
 
        self.ylabel = r"inertial sum rule (eV$^{-1})$ "
        self.plottitle = "evaluated at q=" + str(self.calc.q) + "a.u."
        self.plot_result2()      
        
    def diimfp(self):    
        self.calc.calcDIIMFP()
        self.xlabel = self.particle + " energy loss ω (eV)"
        self.ylabel =  r"DIIMFP (eV$^{-1} \rm{\AA}^{-1}$)"
        self.label_curve1 = r"DIIMFP, $E_0=${:.1f}".format(
            self.calc.E0)+"keV"
        self.plottitle = r"IMFP: {:.3g} $\AA$, stopping:  {:.3g} eV/$\AA$,".format(float(self.calc.StoppingResultArray[0]), float(self.calc.StoppingResultArray[1]))
        self.plottitle +=  " straggling:  {:.3g} eV$^2/\AA$".format(float(self.calc.StoppingResultArray[2]))
        self.plot_result1()
        
    def IMFP_stop_strag(self): 
        self.calc.calccurves(False)
        if self.calc.axis_choice.get() == 1:
            self.xlabel = self.particle + " energy (keV)"
        else:
            self.xlabel = self.particle +" velocity (a.u.)"
        self.label_curve1 = "IMFP (Å)"
        if self.calc.StoppingUnits.get()==0:
            self.label_curve2 = "stopping (eV/Å) "
            self.calc.StoppingFactor=1.0
        elif self.calc.StoppingUnits.get()==1:
            self.label_curve2 = "stopping (eV / (1E¹⁵ atoms/cm²)) "
            self.calc.StoppingFactor=self.calc.massunitcell.get()/(self.calc.specificweight.get()*6.022)
        elif self.calc.StoppingUnits.get()==2:
            self.label_curve2 = "stopping (MeV / (mg/cm²))"    
            self.calc.StoppingFactor=0.1/self.calc.specificweight.get()
        self.label_curve3 = "straggling (eV²/Å)"
        self.plot_result3()   
        
    def eq_plot(self):  
        self.xlabel = "--"
        self.calc.bulk_eq = True
        if self.calc.log_choice.get() == 0:
            self.plottitle = r"$\rm{Im}\left[ \frac{-1}{\epsilon(\omega,q)}\right]$"
        else:
            self.plottitle = r"$\log\left({\rm{Im}\left[ \frac{-1}{\epsilon(\omega,q)}\right]}\right)$"
        self.calc.colorplot_lossfunction()  
        self.calc.scale_image()  
        self.colorplot_result()  
        
    def self_eq_plot(self):  #self_ means surface energy loss function here!
        self.xlabel = "--"
        self.calc.bulk_eq = False  
        if self.calc.log_choice.get() == 0:
            self.plottitle = r"$\rm{Im} \frac{(\epsilon-1)^2}{\epsilon (\epsilon+1)}$"
        else:
            self.plottitle = r"$\log \rm{Im} \frac{(\epsilon-1)^2}{\epsilon (\epsilon+1)}$"
        self.calc.colorplot_lossfunction()
        self.calc.scale_image()    
        self.colorplot_result()  
        
    def mean_ionization_energy(self):   
        self.calc.Mean_Ionisation_Energy()
        self.xlabel = "ω (eV)"
        self.label_curve1 =  ("$I_0$ =  {:.2f}".format(self.calc.I0)
            + "eV,\n $C_0=$ {:.2E}".format(self.calc.C0))
          #  + ", C from Bethe= {:.2E}".format(self.calc.C0_from_BetheSum))
        self.label_curve2 = ("$I$ =  {:.2f}".format(self.calc.MIE)
            + "eV,\n $C=$ {:.2E}".format(self.calc.C1))
           # + ", C from Bethe= {:.2E}".format(self.calc.C1_from_BetheSum))
        self.label_curve3 = ("$I_1$ =  {:.2f}".format(self.calc.Istraggling)
            + "eV,\n $C_2=$ {:.2E}".format(self.calc.C2))
           # + ", C from Bethe= {:.2E}".format(self.calc.C2_from_BetheSum))
        self.ylabel = "Mean Ionization Energy (eV)"
        self.plottitle = "evaluated at q= 0.01 a.u."
        self.plot_result2() 
    
    def SELF_DSEP(self):
        self.calc.surfaceloss()
        self.xlabel =self.particle + " energy loss(eV)"
        self.label_curve1 = "surf Loss Func"
        self.label_curve2 =  r"DSEP, $\theta = ${:.1f}$^\circ$".format(self.calc.thetaIn.get())+"\n integr.prob.:  {:.3g}".format(
            self.calc.SurfExProb
        ) 
        self.plottitle = "SELF evaluated at q=" + str(self.calc.q) + "a.u."
        self.ylabel  = "Surf. loss function "
        self.ylabel2 = "DSEP at $E_0=${:.1f}".format(
            self.calc.E0)+"keV"
        if self.calc.DSEP_choice.get() == 0:
            self.plottitle += ", from global eps"
        else:
            self.plottitle += ", from  eps per oscillators"
        self.plot_result2_dif_y_axes()
        
    def Compton(self):  
        self.calc.CalcCompton()
        self.xlabel = "k (a.u.)"
        self.label_curve1 = "Compton profile"
        self.plottitle = "evaluated at q=" + str(self.calc.q_Compton) + "a.u."
        self.ylabel = "electrons per a.u. per atom"
        self.plot_result1()  
        
    def partial_DIIMFP(self):
        self.xlabel = "ω (eV)"
        self.ylabel =  r"DIIMFP (eV$^{-1} \rm{\AA}^{-1}$)"
        self.plottitle = "partial DIIMFP, "+self.particle+ ", $E_0=${:.1f}".format(
            self.calc.E0)+"keV"
        self.calc.weight = 0.0
        self.calc.PartDIIMFP()
        self.plot_result10() 
        
    def partial_stopping(self):
        self.xlabel = "ω (eV)"
        self.plottitle = "partial DIIMFP, " +self.particle + " weighted by ω"+self.particle+ ", $E_0=${:.1f}".format(
            self.calc.E0)+"keV"
        self.calc.weight = 1.0
        self.calc.PartDIIMFP()
        self.plot_result10()
        
    def oscillator_strength(self):
        self.calc.Calc_Os_Strength()
        self.xlabel = "ω (eV)"
        self.label_curve1 = r"Oscillator Strength (eV$^{-1}$)"
        self.plottitle = "evaluated at q=" + str(self.calc.q) + "a.u."
        self.ylabel = ""
        self.plot_result1()
        
    def dyn_struct_factor(self):
        self.calc.dyn_struct_factor()
        self.xlabel = "ω (eV)"
        self.label_curve1 = r"$S(k,\omega)$"
        self.plottitle = "evaluated at q=" + str(self.calc.q) + "a.u."
        self.ylabel = r"$S(k,\omega)$, eV$^{-1}$"
        self.plot_result1()
        
        
    def S_k_omega_rule(self):
        self.calc.dyn_struct_factor()
        self.xlabel = "ω (eV)"
        self.label_curve1 = r"$S(k,\omega)$"
        self.label_curve2 = r"$\int_0^\omega S(k,\omega) \omega d \omega$"
        self.plottitle = "evaluated at q=" + str(self.calc.q) + "a.u."
        self.ylabel = r"$S(k,\omega)$, eV$^{-1}$"
        self.ylabel2  = r"sum rule"
        self.plot_result2_dif_y_axes()        
        
    def pseudo_charge_density(self):
        self.calc.PseudoChargeDensity()
        self.xlabel = "fraction of unit cell"
        self.label_curve1 = "minimum charge density"
        self.plottitle = ("calculated from ELF up to " + str(self.calc.LastEnergy) + " eV")
        self.ylabel = "pseudocharge density (e$^-$/ Å$^3$)"
        self.plot_result1()
        
    def radial_charge_density(self):   
        self.calc.PseudoChargeDensity()
        self.calc.ConvertToRadialPseudoChargeDensity()
        self.xlabel = "r (Å)"
        self.label_curve1 = " pseudo charge density (e$^-$/ Å$^3$)"
        self.plottitle = ("from ELF up to " + str(self.calc.LastEnergy)+" eV"
            + ", MT Radius  {:.3f} Å".format(self.calc.MT_radius))
        self.ylabel = "charge density (e$^-$/ Å$^3$)"
        self.plot_result1()    
         
    def radial_w_p_lambda(self):
        self.calc.PseudoChargeDensity()
        self.calc.ConvertToRadialPseudoChargeDensity()
        self.calc.stopping_IMFP_w_p_versus_r()
        self.plottitle = ( self.particle_LaTeX + ", $E_0=${:.1f} keV, ".format(self.calc.E0)+ "from ELF up to " + str(self.calc.LastEnergy)+" eV")
        self.xlabel = "r (Å)"
        self.label_curve1 = r"$\omega_p$ (eV)"
        self.label_curve2 = r" IMFP $\lambda$ ($\AA$)"
        self.label_curve3 = r" stopping (eV/$\AA$)"
        self.label_curve4 = r" straggling  (eV$^2$/$\AA$)"
        self.plot_result4()    
        
    def shell_effect_all(self):
        self.calc.shell_effect(False)
        # if self.calc.axis_choice.get() == 1:
            # self.xlabel = self.particle + " energy (keV)"
        # else:
            # self.xlabel = self.particle + " velocity (a.u.)"
            # for i in range (self.calc.NStopping):
                # self.calc.x_axis[i]= self.calc.CurvesVelocity[i]

        self.plottitle = "shell effect ($I =$ {:.3f} eV)".format(self.calc.MIE)
        self.label_curve1 = "A = L DF"
      
        self.label_curve2 = r"B =  $\max(\ln{\frac{2 v ^2}{I}}+ \ln{\gamma^2}- \beta^2 +0.5f(\gamma),0)$"
        self.label_curve3 = r"shell effect: A - B ($\approx -C/Z_2$)"
        self.label_curve4 =  r"$0.5f(\gamma)$" 
        self.plot_result4()
        
    # def shell_effect_soft(self):
        # self.calc.shell_effect(True)
        # if self.calc.axis_choice.get() == 1:
            # self.xlabel = self.particle + " energy (keV)"
        # else:
            # self.xlabel = self.particle + " velocity (a.u.)"
            # for i in range (self.calc.NStopping):
                # self.calc.x_axis[i]= self.calc.CurvesVelocity[i]

        # self.plottitle = "shell effect, soft col., ($I =$ {:.3f} eV)".format(self.calc.MIE)
        # self.label_curve1 = "L Bethe - L calc"
        # self.label_curve2 = "L Bethe rel. - L calc"
        # self.plot_result2()
        
    def DDCS_at_omega(self): 
        self.calc.DDCS_at_omega()
        self.plottitle = (r"$ \frac{d\sigma}{d\omega d\Omega}$" + "  {:.0f} keV,".format(self.calc.E0)
            +self.particle_LaTeX)
        self.ylabel = r"DDCS  $\AA^2$/sr/eV, per unit cell"
        self.xlabel = r"$\theta$ (mrad)"
        self.label_curve1= r"$\omega = $"+ "  {:.1f} eV, no retardation".format(self.calc.omega_ddcs)
        self.label_curve2= r"$\omega = $"+ "  {:.1f} eV, with retardation".format(self.calc.omega_ddcs)
        self.plot_result2() 
        
    def DDCS_at_theta(self):  
        self.xlabel = r"$\omega$ (eV)"
        self.ylabel = r"DDCS per unit cell  ($\AA^2$/sr/eV)" 
        self.label_curve1 = r"$\theta = $ {:.3g} mrad, no retardation".format(self.calc.theta_ddcs)
        self.label_curve2 = r"$\theta = $ {:.3g} mrad, incl. retardation".format(self.calc.theta_ddcs)
        self.plottitle = (self.particle_LaTeX + " {:.0f} keV,".format(self.calc.E0)) 
        self.calc.DDCS_at_theta()
        self.plot_result2()         
        
    def dcs_omega_eq_plot(self):
        self.calc.colorplot_ddcs()
        self.xlabel = "--"

        if self.calc.log_choice.get() == 0:
            self.plottitle = (
                r"$ \frac{d\sigma}{d\omega d\Omega}$" + r"  ( $\AA^2$/sr/eV, per unit cell )"+ "  {:.1f} keV,".format(self.calc.E0)
                +self.particle_LaTeX
            )
        else:
            self.plottitle = (
               r"$ \log (\frac{d\sigma}{d\omega d\Omega})$" + r"  ( $\AA^2$/sr/eV, per unit cell )  "+ " {:.1f} keV,".format(self.calc.E0)
                + self.particle_LaTeX
            )
        self.calc.scale_image()    
        self.colorplot_result()
        
    def Cerenkov(self):
        self.calc.colorplot_Cerenkov()
        self.xlabel = "--"
 
        if self.calc.log_choice.get() == 0:
            self.plottitle = (
                r"$ \frac{d\sigma}{d\omega d\Omega}$" + r"  ( $\AA^2$/sr/eV, per unit cell )"+ "  {:.1f} keV,".format(self.calc.E0)
                +self.particle_LaTeX
            )
        else:
            self.plottitle = (
               r"$ \log (\frac{d\sigma}{d\omega d\Omega})$" + r"  ( $\AA^2$/sr/eV, per unit cell )  "+ " {:.1f} keV,".format(self.calc.E0)
                + self.particle_LaTeX
            )
        self.calc.scale_image()    
        self.colorplot_result()
        
    def  difference_due_to_Cerenkov(self):
        self.calc.difference_due_to_Cerenkov()
        self.xlabel = "--"
        if self.calc.log_choice.get() == 0:
            self.plottitle = ( self.particle_LaTeX +
                r", $ \frac{d\sigma}{d\omega d\Omega}$ Cerenkov  $\frac{dE}{dx}=$" + "  {:.3g} ".format(self.calc.stopping_due_to_photons ) + r" eV/$\AA$")
        else:
            self.plottitle = (self.particle_LaTeX +
               r"$, \log (\frac{d\sigma}{d\omega d\Omega})$ Cerenkov  $\frac{dE}{dx}=$" + "  {:.3g} ".format(self.calc.stopping_due_to_photons ) + r" eV/$\AA$")
        self.calc.scale_image()    
        self.colorplot_result()   
            
        

        
    def dcs_plot(self):
        self.calc.DCS()
        self.xlabel = r"$\theta$ (mrad)"
        self.ylabel = r"DCS per unit cell $\AA^2$/sr" 
        self.label_curve1 =  r"DCS" 
        self.label_curve2 =  r"DCS incl. retardation"
        self.label_curve3 = "Rutherford" 
        self.plottitle = (
               r"$  (\frac{d\sigma}{ d\Omega})$" + r"  ( $\AA^2$/sr, per unit cell )  "+self.particle_LaTeX + 
               " {:.0f} keV,".format(self.calc.E0))
        self.plot_result2() 
        
       
        
  
        
    def REELS(self):
        if self.particle == "proton":
            return
        self.calc.REELS_spectrum()
        self.xlabel = "ω (eV)"
        self.label_curve1 = "Intensity (eV⁻\u2071)"
        self.plottitle = "REELS spectrum, {:.1f} keV,".format(self.calc.E0)
        self.plot_result1()
            
# ==========================plotting==========================================

    def plot_result1(self):
        
        plt.figure()
        plt.figtext(0.02, 0.95, self.PlotDescription.get(), fontsize=self.legendfontsize.get())
        if(self.title_fontsize.get()> 0): 
            plt.figtext(0.5, 0.9, self.plottitle, horizontalalignment="center", fontsize=self.title_fontsize.get())

        ax = plt.subplot(111)
        ax.tick_params(direction="in", which="both", right=1, top=1)
     
        plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.11)
        plt.plot(self.calc.x_axis, self.calc.Result1,  marker='.', label=self.label_curve1)
        if (self.calc.Overplot.get() == 1) and (self.calc.Npointscomparray > 0):
            plt.plot(self.calc.xCompArray, self.calc.yCompArray, label=self.LiteratureDescription.get(), color="red",marker='.', linestyle=self.overplotline)

        if(self.MyChapApp.LogX.get() == 1):
            plt.xscale("log")
        if (self.MyChapApp.LogY.get() == 1):
            plt.yscale("log")
       
        ymin,ymax = ax.get_ylim()
        if(ymax < 0.001) or (ymax > 10000):
            ax.yaxis.set_major_formatter(FormatStrFormatter('%1.0E')) 
            plt.subplots_adjust(left=0.2)
        xmin,xmax = ax.get_xlim()
        if(xmax < 0.001) or (xmax > 10000):
            ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0E'))     
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        if(self.legendfontsize.get() > 0):
             plt.legend(loc=self.legendposition.get())
        self.plotshow(plt)

      
    

    def plot_result2(self):
        if self.calc.Energy_Scale_choice.get() == 1: 
            ch =self.MyChapApp.plotchoice.get()
            if ch == "eps_w_direct" or ch == "one_over_eps_w_direct" or ch == "Kramers-Kronig tests" \
                or ch == "one_over_eps_w_kk" or ch == "oscillator_strength" or ch == "n_k": 
                    for i in range(self.calc.NPoints):
                        self.calc.x_axis[i] =1239.84/self.calc.x_axis[i]
                    self.xlabel=r"nm$^{-1}$"   
                         

        fig = plt.figure()
        plt.figtext(0.02, 0.95, self.PlotDescription.get(), fontsize=self.legendfontsize.get())
        if(self.title_fontsize.get()> 0): 
            plt.figtext(0.5, 0.90, self.plottitle, horizontalalignment="center",fontsize= self.title_fontsize.get())
        
        ax = fig.add_subplot(111)
        ax.tick_params(direction="in", right=1, top=1)
        plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.11)
        plt.plot(self.calc.x_axis, self.calc.Result1,  marker='.', label=self.label_curve1)
        if(self.MyChapApp.LogX.get() == 1):
            plt.xscale("log")
        if (self.MyChapApp.LogY.get() == 1):
            plt.yscale("log")

        plt.plot(self.calc.x_axis,self.calc.Result2,linestyle="dashed",
                marker='.', label=self.label_curve2)
   
        # if self.MyChapApp.plotchoice.get() == "mean_ionization_energy":
             # plt.plot(self.calc.x_axis,self.calc.Result3,linestyle="dotted", marker='.', label=self.label_curve3)        

        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        if self.MyChapApp.plotchoice.get() == "dcs_plot":
            ax.set_ylim(top=1.1 * np.amax(self.calc.Result2), bottom =0.9 * np.amin(self.calc.Result1))
            plt.plot(self.calc.x_axis,self.calc.Result3,linestyle="dotted", marker='.', label=self.label_curve3)
        if self.MyChapApp.plotchoice.get() ==  "inertial_rules":
            plt.plot(self.calc.x_axis,self.calc.Result4,linestyle="dotted", marker='.', label=self.label_curve4)
            
            plt.plot(
            [self.calc.x_axis[0], self.calc.x_axis[self.calc.NPoints - 1]],
            [self.calc.limitingvalue, self.calc.limitingvalue],
            color="firebrick", label=self.label_curve3, linestyle="dashdot")
            
        if (self.calc.Overplot.get() == 1) and (self.calc.Npointscomparray > 0):
            plt.plot(self.calc.xCompArray, self.calc.yCompArray, label=self.LiteratureDescription.get(), color="red",marker='.',linestyle=self.overplotline)
        if(self.legendfontsize.get() > 0):
            plt.legend(loc=self.legendposition.get())
        self.plotshow(plt)
        
        
    def plot_result4(self):
        ch =self.MyChapApp.plotchoice.get()
        if ch == "shell_effect_all":                 
            if self.calc.axis_choice.get() == 1:
                if self.particle== "electron":
                    self.calc.x_axis = self.calc.CurvesEnergy
                    self.xlabel = "electron energy (keV)"
                else:
                    for i in range(self.calc.NStopping):
                        self.calc.x_axis[i] = self.calc.CurvesEnergy[i]/1000.0
                    self.xlabel = "proton energy (MeV)"
            else:
                self.calc.x_axis = self.calc.CurvesVelocity
                self.xlabel = self.particle + " velocity (a.u.)"
        fig = plt.figure()
        plt.figtext(0.02, 0.95, self.PlotDescription.get(), fontsize=self.legendfontsize.get())
        if(self.title_fontsize.get()> 0): 
            plt.figtext(0.5, 0.90, self.plottitle, horizontalalignment="center", fontsize=self.title_fontsize.get())
      
        ax = fig.add_subplot(111)
        ax.tick_params(direction="in", right=1, top=1)
        plt.plot(self.calc.x_axis, self.calc.Result1, marker='.', label=self.label_curve1)
        plt.plot(self.calc.x_axis,self.calc.Result2,linestyle="dashed", marker='.', label=self.label_curve2)
 
        plt.plot(self.calc.x_axis,self.calc.Result3,linestyle="dashed",label=self.label_curve3,color="brown")
        plt.plot(self.calc.x_axis,self.calc.Result4,linestyle="dotted",label=self.label_curve4, color="black")
            
        if (self.calc.Overplot.get() == 1) and (self.calc.Npointscomparray > 0):
            plt.plot(self.calc.xCompArray, self.calc.yCompArray, label=self.LiteratureDescription.get(), color="red",marker='.', linestyle=self.overplotline)
        if(self.MyChapApp.LogX.get() == 1):
            plt.xscale("log")
        if (self.MyChapApp.LogY.get() == 1):
            plt.yscale("log")
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        if(self.legendfontsize.get() > 0):
            plt.legend(loc=self.legendposition.get())
        self.plotshow(plt)
        
    def plot_result3(self):
        fig, axs = plt.subplots(3, 1, sharex=True, figsize=(10, 10))
        plt.figtext(0.02, 0.95, self.PlotDescription.get(), fontsize=self.legendfontsize.get())

        if(self.title_fontsize.get()> 0): 
            plt.figtext(0.5, 0.90, self.plottitle, horizontalalignment="center", fontsize=self.title_fontsize.get())
        fig.subplots_adjust(hspace=0)
        if self.calc.axis_choice.get() == 1:
            x_axis_array = self.calc.CurvesEnergy
        else:
            x_axis_array = self.calc.CurvesVelocity
            

        if self.calc.imfp_crosssection_choice.get() == 0:
            axs[0].plot(x_axis_array, self.calc.IMFPEnergy, label="IMFP (Å)")
            axs[0].set_ylim(top=1.2 * self.calc.IMFPEnergy[self.calc.NStopping - 1])
            if self.calc.Projectile.get() == 0 and np.amax(self.calc.TPP_IMFPEnergy) > 0.0:
#if np.amax(self.calc.TPP_IMFPEnergy) > 0.0 and :
                axs[0].plot( x_axis_array, self.calc.TPP_IMFPEnergy, color="purple",
                    label=r"TPP-2m $\omega_p$="+ f"{self.calc.w_p_TPP.get():.2f}" + "eV\n$\\rho=$"
                    + f"{self.calc.specificweight.get():.2f}" + " g/cm$^3$",linestyle="dotted")
   
            # assumes maximum IMFP at highest energy
            #       axs[0].plot(x_axis_array, self.calc.BetheIMFPEnergy,
            #             color='firebrick', label="from Bethe dispersion\n $\omega_p$="+
            #             f"{self.calc.w_p_used :.1f}" +"eV", linestyle="--")
            if self.calc.Approximations.get() == 1:
                axs[0].plot(
                    x_axis_array,
                    self.calc.DL_IMFPaverage_Energy,
                    color="orchid",
                    label="Average DL IMFP: $C_{av}$="
                    + f"{self.calc.C0:.3f}"
                    + ", $I'_{av}=$"
                    + f"{self.calc.I0 :.1f}"
                    + "eV",
                    linestyle="-.",
                )

                axs[0].plot(
                    x_axis_array,
                    self.calc.DL_IMFP_sum_Energy,
                    color="black",
                    label="summed DL oscillators",
                    linestyle="dotted",
                )
                axs[0].plot(
                    x_axis_array,
                    self.calc.DL_IMFP_from_ELF,
                    color="limegreen",
                    label="DL IMFP from ELF",
                    linestyle="dotted",
                )
                if np.amax(self.calc.TPP_IMFPEnergy) > 0.0:
                    axs[0].plot(
                        x_axis_array,
                        self.calc.TPP_IMFPEnergy,
                        color="purple",
                        label=r"TPP-2m $\omega_p$="
                        + self.calc.Omegas[0].get()
                        + "eV\n$\\rho=$"
                        + f"{self.calc.specificweight.get():.2f}"
                        + " g/cm$^3$",
                        linestyle="dotted",
                    )
        else:
            axs[0].plot(
                x_axis_array,
                np.reciprocal(self.calc.IMFPEnergy) / self.calc.UnitCellDensity,
                label=r"cross section ($\AA^2$)",
            )
        axs[1].plot(x_axis_array, self.calc.StoppingEnergy*self.calc.StoppingFactor, label=self.label_curve2)
        axs[1].plot(x_axis_array, self.calc.BetheStoppingEnergy_Salvat*self.calc.StoppingFactor, color="firebrick",
            label=r"$\frac{4\pi}{v^2} N Z(L_0 + 0.5*F(\gamma))$" + "\n(" + r"$I=$" + f"{self.calc.MIE :.1f}" + "eV)",
            linestyle="--")
        if self.particle =="proton":
            axs[1].plot(x_axis_array, self.calc.LinearApprox_lowE*self.calc.StoppingFactor, color="limegreen",
                label=r"$-\frac{dE}{dx}=c v$, $c =$"+ f"{self.calc.maxsloop:.1f}" + r"eV/($\AA v_0$)",
                linestyle="--")    
        
        if self.calc.Approximations.get() == 1:
        
           
            axs[1].plot(
                x_axis_array,
                self.calc.DL_StoppingEnergy*self.calc.StoppingFactor,
                color="orchid",
                label="DL Formula from ELF\n$C_1= $"
                + f"{self.calc.C1 :.3f}"
                + r", $\omega_p$="
                + f"{self.calc.MIE :.1f}"
                + "eV",
                linestyle="-.",
            )
            if np.amax(self.calc.DL_Stopping_sum_Energy) > 0.0:
                axs[1].plot(
                    x_axis_array,
                    self.calc.DL_Stopping_sum_Energy*self.calc.StoppingFactor,
                    color="black",
                    label="summed DL oscillators",
                    linestyle="dotted",
                )
            axs[1].plot(
                x_axis_array,
                self.calc.DL_Stopping_from_ELF*self.calc.StoppingFactor,
                color="limegreen",
                label="DL stopping from ELF",
                linestyle="dotted",
            )
        axs[2].plot(x_axis_array, self.calc.StragglingEnergy, label=self.label_curve3)
        if self.particle == "proton":
            scalingfactor=1.0
        else:
            scalingfactor=4.0    
        axs[2].plot(
            [x_axis_array[0], x_axis_array[self.calc.NStopping - 1]],
            [self.calc.BohrStraggling/scalingfactor, self.calc.BohrStraggling/scalingfactor],
            color="darkgreen", label="Bohr Limit", linestyle="dashdot")
        if self.particle == "proton":
            axs[2].plot( x_axis_array, self.calc.Straggling_Jackson, color="firebrick",
                label="Jackson eq. 13.50",
                linestyle="--")
                
            
        if self.calc.Approximations.get() == 1:
            axs[2].plot( x_axis_array, self.calc.BetheStragglingEnergy, color="firebrick",
                label="Bethe Formula\n"+ r"$\omega_p$=" + f"{self.calc.MIE :.1f}" + "eV",
                linestyle="--")
            if (
                np.amax(self.calc.DL_StragglingEnergy) > 0.0
            ) and self.calc.OscillatorsPresent:
                axs[2].plot(
                    x_axis_array,
                    self.calc.DL_StragglingEnergy,
                    color="orchid",
                    label="DL Formula from zero width osc.\n  $C_2$="
                    + f"{self.calc.C2 :.3f}"
                    + r", $\omega_p$="
                    + f"{self.calc.Istraggling :.1f}"
                    + "eV",
                    linestyle="-.",
                )
            elif np.amax(self.calc.DL_StragglingEnergy) > 0.0:
                axs[2].plot(
                    x_axis_array,
                    self.calc.DL_StragglingEnergy,
                    color="orchid",
                    label="DL Formula from ELF\n,  $C_2$="
                    + f"{self.calc.C2 :.3f}"
                    + r", $\omega_p$ ="
                    + f"{self.calc.Istraggling :.1f}"
                    + " eV",
                    linestyle="-.",
                )
            if np.amax(self.calc.DL_Straggling_sum_Energy) > 0.0:
                axs[2].plot(
                    x_axis_array,
                    self.calc.DL_Straggling_sum_Energy,
                    color="black",
                    label="summed DL oscillators",
                    linestyle="dotted",
                )
            axs[2].plot(
                x_axis_array,
                self.calc.DL_Straggling_from_ELF,
                color="limegreen",
                label="DL straggling from ELF",
                linestyle="dotted",
            )

        if self.calc.Overplot.get() == 1:
            axs[self.calc.comp_option_choice.get()].plot(self.calc.xCompArray, self.calc.yCompArray,
                label=self.LiteratureDescription.get(),linestyle=self.overplotline,marker='.', color="red")
          
        plt.xlabel(self.xlabel)
        for i in range(3):
            axs[i].legend()
            if self.calc.axis_choice.get() == 1:
                axs[i].set_xscale("log")
            if (self.MyChapApp.LogY.get() == 1):
                axs[i].set_yscale("log")    
            axs[i].tick_params(
                axis="both", which="both", direction="in", right=1, top=1
            )
            axs[i].set_ylim(bottom=0)
        
        axs[1].set_ylim(top=1.1 * np.nanmax(self.calc.StoppingEnergy)*self.calc.StoppingFactor)
        self.plotshow(plt)
   

    def plot_result10(self):
        ipoints = self.calc.NPoints
        self.my10labels=[]
        plt.figure()
        plt.figtext(0.02, 0.95, self.PlotDescription.get(), fontsize=self.legendfontsize.get())
        if(self.title_fontsize.get()> 0): 
            plt.figtext(0.5, 0.90, self.plottitle, horizontalalignment="center", fontsize=self.title_fontsize.get())
        ax = plt.subplot(111)
        ax.tick_params(direction="in", right=1, top=1)
        plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.11)
        q_step = self.calc.LastMomentum / 10.0
        for i in range(10):
            q_lower = i * q_step
            q_upper = (i + 1) * q_step
            if self.MyChapApp.plotchoice.get() == "partial_DIIMFP":
                self.my10labels.append(f"{q_lower:.1f}" + "< q <" + f"{q_upper:.1f}" + ", λ\u209a="
                    + f"{self.calc.PartIntSum[i]:7.3g}" + "Å")
            else:
                self.my10labels.append(f"{q_lower:.1f}" + "< q <" + f"{q_upper:.1f}" + r", S_p "
                    + f"{self.calc.PartIntSum[i]:7.3g}" + "eV/Å" )
            plt.plot(self.calc.x_axis, self.calc.partialresults[:, i],
                label=self.my10labels[i])
        if(self.MyChapApp.LogX.get() == 1):
            ax = plt.gca()
            xmin,xmax = ax.get_xlim()
           
            plt.xscale("log")
        if (self.MyChapApp.LogY.get() == 1):
            plt.yscale("log")        
        if(self.legendfontsize.get() > 0):
             plt.legend(loc=self.legendposition.get())
             ymin,ymax = ax.get_ylim()
        if(ymax < 0.001) or (ymax > 10000):
            ax.yaxis.set_major_formatter(FormatStrFormatter('%1.0E')) 
            plt.subplots_adjust(left=0.2)
        xmin,xmax = ax.get_xlim()
        if(xmax < 0.001) or (xmax > 10000):
            ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0E'))       
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        self.plotshow(plt)

    def colorplot_result(self):
        fig=plt.figure()
        ax = fig.add_subplot(111)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        
        plt.figtext(0.02, 0.95, self.PlotDescription.get(), fontsize=self.legendfontsize.get())
        if (self.MyChapApp.plotchoice.get() == "eq_plot") or (self.MyChapApp.plotchoice.get() == "self_eq_plot"):
            ratio = self.calc.LastMomentum / (
                self.calc.LastEnergy - self.calc.CenterFirstBin
            )
           # maximum scale and log/linear are taken care of in calc.
            plt.imshow(
                self.calc.my_scaled_image,
                extent=[
                    0,
                    self.calc.LastMomentum,
                    self.calc.LastEnergy,
                    self.calc.CenterFirstBin,
                ],
                aspect=ratio
            )
           
            if (self.MyChapApp.plotchoice.get() == "eq_plot"):        
                plt.xlabel("q (a.u)")
            else:
                plt.xlabel(r"$q_\parallel$"+"(a.u.)")    
            plt.xlim(0, self.calc.LastMomentum)
        else:  # self.MyChapApp.plotchoice.get() == "dcs_omega_eq_plot" or cerenkov
            ratio = self.calc.theta_max / self.calc.LastEnergy
            plt.imshow(self.calc.my_scaled_image,
                extent=[0, self.calc.theta_max, self.calc.LastEnergy, self.calc.CenterFirstBin],
                aspect=ratio)
            
            plt.xlabel(r"$\theta$ (mrad)")
            plt.xlim(0, self.calc.theta_max)
        plt.ylim(self.calc.CenterFirstBin, self.calc.LastEnergy)

        plt.colorbar()
        plt.ylabel("ω (eV)")
        if(self.title_fontsize.get()> 0): 
            plt.figtext(0.5, 0.91, self.plottitle, horizontalalignment="center", fontsize=self.title_fontsize.get())


        if self.calc.E0 < 10000:    
            mylabel = ("upper limit ω\n" + f"{self.calc.E0:.2f}" + " keV " + self.particle_LaTeX)
        else:
            mylabel = ("upper limit ω\n" + f"{self.calc.E0/1000:.2f}" + " MeV " + self.particle_LaTeX)  
        if self.MyChapApp.plotchoice.get() == "eq_plot" and self.calc.E0 > 0.0:
            self.calc.Calculate_Integration_limits()
          #  plt.plot(self.calc.xArray, self.calc.yArray, color="lime", label=mylabel + " non-relativistic")  #uncomment this if you want the non-relativistic boundary as well
           
            plt.plot(self.calc.xArray, self.calc.yArray_relativistic, color="red", label=mylabel)
            if(self.legendfontsize.get() > 0):
                plt.legend(loc=self.legendposition.get())  
        elif  self.MyChapApp.plotchoice.get() == "dcs_omega_eq_plot"  and self.calc.particle == "proton":  # limiting lines for protons

            if self.calc.theta_max > 1000.0/cnst.Mp:
                plt.plot([1000/cnst.Mp, 1000/cnst.Mp], [self.calc.CenterFirstBin,self.calc.LastEnergy ],
                   color="red",linestyle="dashed",  label=r"$\sin(\theta_{\rm max}) =\frac{M_e}{M_p}$")
                   #https://www.physicsforums.com/threads/solving-scattering-angle-problem-when-m-m-or-m-m.28662/ for derivation
                   #or Sigmund 3.1.1
                
            maxloss= 2*(self.calc.beta_r *self.calc.gamma_r* cnst.C)**2 * cnst.HARTREE
            if maxloss < self.calc.LastEnergy :  
                plt.plot([0, self.calc.theta_max], [maxloss,maxloss],
                   color="red",linestyle="dotted", label=r"$2v^2$ * 27.211")
                   
                   
            if(self.legendfontsize.get() > 0  ):
                plt.legend(loc="center left")
        plt.subplots_adjust(bottom=0.14) 
        self.plotshow(plt)
        
        
    def plot_result2_dif_y_axes(self):
        fig = plt.figure()
        plt.figtext(0.02, 0.95, self.PlotDescription.get(), fontsize=self.legendfontsize.get())
        if(self.title_fontsize.get()> 0): 
            plt.figtext(0.5, 0.93, self.plottitle, horizontalalignment="center",fontsize= self.title_fontsize.get())
        ax = fig.add_subplot(111)
        ax.tick_params(direction="in", right=1, top=1)
        plt.plot(self.calc.x_axis, self.calc.Result1, label=self.label_curve1)
        if (self.calc.Overplot.get() == 1) and (self.calc.Npointscomparray > 0):
            plt.plot(self.calc.xCompArray, self.calc.yCompArray, label=self.LiteratureDescription.get(), color="red",marker='.',linestyle=self.overplotline)    
        plt.subplots_adjust(left=0.16, right=0.83, top=0.9, bottom=0.11)
        if(self.legendfontsize.get() > 0):
            ax.legend(loc="lower right")
        ymin,ymax = ax.get_ylim()
        if(ymax < 0.01) or (ymax > 100):
            ax.yaxis.set_major_formatter(FormatStrFormatter('%1.1E')) 
            plt.subplots_adjust(left=0.21)
        xmin,xmax = ax.get_xlim()
        if(xmax < 0.01) or (xmax > 1000):
            ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0E'))    
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        if(self.legendfontsize.get() > 0):
            plt.legend(loc=self.legendposition.get())
        ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

        color = "tab:red"
    
        ax2.plot(self.calc.x_axis, self.calc.Result2, color=color, linestyle="dashed",
            label=self.label_curve2)
        
        if(self.MyChapApp.LogX.get() == 1):
            lower_limit=  int(self.calc.NPoints/10) 
        else:    
            lower_limit= self.calc.NPoints - int(self.calc.NPoints/4)    
        if self.ylabel2  == r"sum rule":    
            ax2.plot(
                [self.calc.x_axis[lower_limit], self.calc.x_axis[self.calc.NPoints - 1]],
                [self.calc.q**2/2.0*cnst.HARTREE, self.calc.q**2/2.0*cnst.HARTREE],
                color="darkgreen", label=r"$q^2/2$", linestyle="dashdot")    
        ax2.tick_params(axis="y", labelcolor=color)

        ax2_ylims = ax2.axes.get_ylim()  # Find y-axis limits set by the plotter
        ax2.set_ylim(top=ax2_ylims[1] * 1.1)  # for esthetics only
        ax2.set_ylabel(self.ylabel2, color=color, rotation=270, labelpad=15)#, labelpad=-15)
        if(ax2_ylims[1] * 1.1 < 0.01) or (ax2_ylims[1] * 1.1 > 1000):
            ax2.yaxis.set_major_formatter(FormatStrFormatter('%1.1E')) 
            plt.subplots_adjust(right=0.8)
        ax2.legend(loc="lower right")
        if(self.MyChapApp.LogX.get() == 1):
            xmin,xmax = ax.get_xlim()
            plt.xscale("log")
        self.plotshow(plt)
        
    
        
    def plot_resultSumrules(self):
        fig = plt.figure()
        plt.figtext(0.02, 0.95, self.PlotDescription.get(), fontsize=self.legendfontsize.get())
        if(self.title_fontsize.get()> 0): 
            plt.figtext(0.5, 0.90, self.plottitle, horizontalalignment="center", fontsize= self.title_fontsize.get())
        ax = fig.add_subplot(111)
        ax.tick_params(direction="in", right=1, top=1)
        plt.plot(self.calc.x_axis, self.calc.Result1, label=self.label_curve1)
        plt.plot(
            self.calc.x_axis, self.calc.Result2, linestyle="dashed", label=self.label_curve2
        )
        plt.plot(
            self.calc.x_axis, self.calc.Result3, linestyle="dotted", label=self.label_curve3
        )
        ax.legend(loc="best")
        plt.xlabel(self.xlabel)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        plt.ylabel(self.ylabel)
        if(self.legendfontsize.get() > 0):
            plt.legend(loc='center right')
        ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

        color = "tab:purple"
        ax2.plot( self.calc.x_axis,  self.calc.Result4, color=color, linestyle="dashed",label=self.label_curve4)
        ax2.tick_params(axis="y", labelcolor=color)
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax2_ylims = ax2.axes.get_ylim()  # Find y-axis limits set by the plotter
        ax2.set_ylim(top=ax2_ylims[1] * 1.1)  # for esthetics only
        ax2.set_ylabel("kk sum rule", color=color, rotation=270, labelpad=10)
        ax2.legend(loc="lower right")
        if(self.MyChapApp.LogX.get() == 1):
       #     xmin,xmax = ax.get_xlim()
            plt.xscale("log")
        self.plotshow(plt)
        
        
    def plotshow(self, plt):
        plt.get_current_fig_manager().window.wm_geometry(f"+{self.start_x}+{self.start_y}") 
        self.start_y += self.dy
        if self.start_y + 480 > self.screen_height:
            self.start_y=0
            self.start_x += self.dx
            if self.start_x + 640 > self.screen_width:
                self.start_x=100
        plt.show()        


    def Close_all_plots(self):
        plt.close("all")
        self.start_x=0
        self.start_y=0
        
        

   
        

# ==============================end plotting======================================      
#from matplotlib documentation
#this seems buggy but useful
        
   ## ***************************************************************************
## * INTERACTIVE KEYMAPS                                                     *
## ***************************************************************************
## Event keys to interact with figures/plots via keyboard.
## See https://matplotlib.org/stable/users/explain/interactive.html for more
## details on interactive navigation.  Customize these settings according to
## your needs. Leave the field(s) empty if you don't need a key-map. (i.e.,
## fullscreen : '')
#keymap.fullscreen: f, ctrl+f   # toggling
#keymap.home: h, r, home        # home or reset mnemonic
#keymap.back: left, c, backspace, MouseButton.BACK  # forward / backward keys
#keymap.forward: right, v, MouseButton.FORWARD      # for quick navigation
#keymap.pan: p                  # pan mnemonic
#keymap.zoom: o                 # zoom mnemonic
#keymap.save: s, ctrl+s         # saving current figure
#keymap.help: f1                # display help about active tools
#keymap.quit: ctrl+w, cmd+w, q  # close the current figure
#keymap.quit_all:               # close all figures
#keymap.grid: g                 # switching on/off major grids in current axes
#keymap.grid_minor: G           # switching on/off minor grids in current axes
#keymap.yscale: l               # toggle scaling of y-axes ('log'/'linear')
#keymap.xscale: k, L            # toggle scaling of x-axes ('log'/'linear')
#keymap.copy: ctrl+c, cmd+c     # copy figure to clipboard     
        
 
        
