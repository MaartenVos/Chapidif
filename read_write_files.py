import re  # for regular expression to extract numbers from string
from tkinter import filedialog
import os
import re
import numpy as np
try:
    from icecream import ic

    ic.configureOutput(includeContext=True)
except ImportError:  # Graceful fallback if IceCream isn't installed.
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a)  # noqa

   # ================================readfile

class read_write_files:
   # def __init__(self, *args, **kwargs):
    def __init__(self,parent):
        self.MyChapApp=parent
        self.calc = self.MyChapApp.calc
        self.runplot  = self.MyChapApp.runplot 
       
        
        
        
        
    def read_df(self):
        self.zero_df()
        filename = filedialog.askopenfilename(
            initialdir=os.getcwd,
            title="Select file",
            filetypes=(("cpd or dat files", ".cpd .dat"), ("all files", "*.*")),
        )
        if not filename:
            return
        fileameonly = os.path.basename(filename)
        self.MyChapApp.fname.set("Comment for/from DF file: " + fileameonly)
    
        myfile = open(filename, "r", errors="ignore")
        myfile.readline()  # skip line
        self.calc.AddELF.set(1)
        line2 = myfile.readline()
        if "Extended" in line2:
            self.calc.DFChoice.set(value=1)
            self.calc.AddELF.set(0)
        elif "Drude Lindhard" in line2:
            self.calc.DFChoice.set(value=2)
        elif " Mermin " in line2 or  "plain Lindhard" in line2:  # keep space around Mermin to avoid being triggered by (Mermin-normalisation) in the Kaneko case
            self.calc.DFChoice.set(value=3)
        if "(non-relativistic)" in line2: 
            self.calc.Dispersion_relativistic.set(0)
        if "(relativistic)" in line2: 
            self.calc.Dispersion_relativistic.set(1)
        if "quadratic dispersion" in line2:    
            self.calc.Dispersion_choice.set(0)
        if "full dispersion" in line2:    
            self.calc.Dispersion_choice.set(1)
        if "add Chi" in line2:
            self.calc.AddELF.set(0)
        elif "add ELF" in line2:
            self.calc.AddELF.set(1)  
        if "(Mermin-renormalisation)"in line2:
            self.calc.Merminize.set(1)
        elif "Direct" in line2:
            self.calc.Merminize.set(2)    
        else:
            self.calc.Merminize.set(0)
                
                  
                
             


        commentline = myfile.readline()  # this is line 3
        if "#" in commentline : commentline  = commentline .replace("#", "")
            
        commentline = commentline.rstrip("\n")
        self.MyChapApp.CommentEntry.delete(0, 'end')
        self.MyChapApp.CommentEntry.insert(0, commentline)
      

        if self.calc.DFChoice.get() == 1:  # these sesinipac files have an extra line
            line4 = myfile.readline()
            result = re.findall(r"\d*\.?\d+", line4)  # from import re
            self.calc.Eps_bkg.set(result[0])
            

        line4 = myfile.readline()
        if line4.find("density") >= 0:
            elements = line4.split(",")
            self.calc.specificweight.set(float(elements[1]))
            self.calc.massunitcell.set(float(elements[3]))
        myfile.readline()  # skip line

        i = 0
        igos = 0
        iKaneko = 0
        iBelkacem = 0
        while 1:
            line = myfile.readline()
            if "#" in line: line = line.replace("#", "")
            
            nums = line.split()
            if len(nums) < 3:
                break
            if not line:
                break
            if "Amplitude" in line:  
                continue
            if "N elec" in line:
                continue
                    
            elif "GOS" in line:
                self.calc.ConcGOS[igos].set(nums[1])
                self.calc.EdgeGOS[igos].set(nums[2])
                self.calc.nlGOS[igos].set(nums[3])
                self.calc.ZGOS[igos].set(nums[4])
                try:
                    self.calc.maxEnergyDensityEffect.set(nums[10])
                except:
                    print("no maximum energy density effect given")    
                igos = igos + 1
                # if "incl. density effect" in line:
                    # self.calc.IncludeGOSDensityEffect.set(1)
                # elif   "no density effect" in line: 
                    # self.calc.IncludeGOSDensityEffect.set(0)
                if "rescaling on" in line:
                    self.calc.ApplySumRuleToGOS.set(1)
                elif  "rescaling off" in line:  
                    self.calc.ApplySumRuleToGOS.set(0)
                    
            elif "Kaneko" in line:
                self.calc.N_Kaneko[iKaneko].set(nums[1])
                self.calc.Q_Kaneko[iKaneko].set(nums[2])
                self.calc.Edge_Kaneko[iKaneko].set(nums[3])
                self.calc.l_Kaneko[iKaneko].set(nums[4])
                self.calc.gamma_Kaneko[iKaneko].set(nums[5])
                iKaneko = iKaneko + 1
                if "modified Kaneko" in line:
                    self.calc.Kaneko_choice.set(0)
                elif "original Kaneko" in line:  
                    self.calc.Kaneko_choice.set(1) 
                        
            elif "Belkacem" in line:
                self.calc.Conc_Belkacem[iBelkacem].set(nums[1])
                self.calc.w_Belkacem[iBelkacem].set(nums[2])
                self.calc.gamma_Belkacem[iBelkacem].set(nums[3])
                iBelkacem = iBelkacem + 1
                
            else:
                self.calc.Amps[i].set(nums[0])
                self.calc.Gammas[i].set(nums[1])
                self.calc.Omegas[i].set(nums[2])
                self.calc.Us[i].set("0.0")
                self.calc.Alphas[i].set("1.0")
                if len(nums) == 4:
                    if self.calc.DFChoice.get() >= 3:
                        self.calc.Us[i].set(nums[3])
                    else:
                        self.calc.Alphas[i].set(nums[3])
                elif len(nums) == 5:
                    self.calc.Alphas[i].set(nums[3])   
                    self.calc.Us[i].set(nums[4])
            i = i + 1
        myfile.close()
 # ================================read oos file
        
    def DF_from_OOS(self):
        #Assumes files in the format the .mat files produced by Salvat's "cbethe"
        self.calc.init()
        self.calc.DF_prop_text.set("\n\n ")
        self.MyChapApp.fname.set("Comment for/from DF file:")
        filename = filedialog.askopenfilename(
            initialdir=os.getcwd,
            title="Select file",
            filetypes=(("oos files", ".mat"), ("all files", "*.*")),
        )
        if not filename:
            return
        fileameonly = os.path.basename(filename)
        self.MyChapApp.fname.set("Comment for/from DF file: " + fileameonly)
        self.MyChapApp.runplot.initcalc()
        

        with open(filename, 'r') as fp:
            for FileLength, line in enumerate(fp):
                pass
        try:
            myfile = open(filename, "r")
        except OSError:
            print("could not open file")
            return 
        Header = True
        HeaderLength = 0
        while  Header:  
            line = myfile.readline()
          #  print(line)
            HeaderLength += 1
            if line.find("Molecular weight") >= 0:
               elements = line.split(" ")
               self.calc.massunitcell.set(float(elements[1]))
            if line.find("Mass density") >= 0:
               elements = line.split(" ")
               self.calc.specificweight.set(float(elements[1]))
            if line.find("Plasma resonance energy") >= 0:
               elements = line.split(" ")
               plasmon_all_electrons=float(elements[1])  
            if line.find("OOS (1/eV") >= 0:
                 Header = False 
        self.calc.N_OOS=FileLength-HeaderLength-1  # there is an empty line at end .mat file. so my file length is off by 1
        for i in range(self.calc.N_OOS):
            line = myfile.readline()
            elements = line.split(" ")
            self.calc.OOSEnergy[i]= float(elements[1])
            self.calc.OOS[i]=float(elements[3])
        myfile.close()
        self.calc.DFChoice.set(2)       
        self.calc.DL_from_OOS()
        return
        
    def write_df_only(self):
        self.write_df()
        
    def zero_df(self):
        self.calc.init()
        self.calc.DF_prop_text.set("\n\n ")
        self.MyChapApp.fname.set("Comment for/from DF file:")  

    def save_calculation(self):  
        self.write_df()  #write dielectric function (can be parsed back into chapidif) 
        self.write_calculation_result() #write calculation details and results (cancurrently not be parsed back into chapidif) 

    # ==========write current configuration and optional calculation result

    def write_df(self):
        # file written so compatible with Sesinipac
        # line one description of file type
        # line 2 model used: Mermin drude etc
        # line 3 comment line
        # line 4 background dielectric constant
        # line 5 amp gamma omega alpha u
        # following lines oscillators
        # lines starting with GOS GOS parameters (not read by sesinipac)

        self.filename = filedialog.asksaveasfilename(
            initialdir=os.getcwd,
            title="Select file",
            filetypes=(("cpd files", "*.cpd"), ("all files", "*.*")),
            defaultextension=".cpd",
        )
        self.now_write()
        
    def now_write(self):
        with open(self.filename, "w", encoding="utf-8") as file:
            fileameonly = os.path.basename(self.filename)
            self.MyChapApp.fname.set("Comment for/from DF file: " + fileameonly)
            # self.copy_values_from_input()  #make sure we have the latest data in the arrays
            file.write("#input file for Chapidif\n")
            model_used="#"
            if self.calc.DFChoice.get() == 1:
                model_used += "Dielectric function model: Extended Drude"
            elif self.calc.DFChoice.get() == 2:
                model_used += "Dielectric function model: Drude Lindhard"
            elif self.calc.DFChoice.get() == 3 and  self.calc.Merminize.get() == 0:
                model_used += "Dielectric function model: plain Lindhard or LL"
            elif self.calc.DFChoice.get() == 3 and  self.calc.Merminize.get() == 1:
                model_used += "Dielectric function model: Mermin MLL"
            elif self.calc.DFChoice.get() == 3 and  self.calc.Merminize.get() == 2:
                model_used += "Dielectric function model: Lindhard-direct"
            if self.calc.AddELF.get() == 0:
                model_used += ", add Chi"
            else:
                model_used += ", add ELF"
            if self.calc.Merminize.get() == 0:    
                model_used +=  ", (without Mermin-renormalisation)"
            elif self.calc.Merminize.get() == 1:
                model_used +=  ", (Mermin-renormalisation)"
            else:     
                model_used +=  ", (Direct)"
            if self.calc.Dispersion_choice.get()== 0:
                model_used +=  ", quadratic dispersion"
            else:
                model_used += ", full dispersion"
            if self.calc.Dispersion_relativistic.get() == 0:
                model_used +=  ", (non-relativistic)"
            else:
                model_used +=  ", (relativistic)"     
            file.write(model_used+"\n")    
            file.write("#"+self.MyChapApp.CommentEntry.get() + "\n")
            if self.calc.DFChoice.get() == 1:
                file.write("#Background dielectric constant:" + self.calc.Eps_bkg.get() + "\n")  #this parameter is not used anymore, sesinipac lagacy

            file.write("#density," + str(self.calc.specificweight.get())
                + ", unit cell mass," + str(self.calc.massunitcell.get()) + "\n")
            if self.calc.DFChoice.get() == 1 or  self.calc.DFChoice.get() == 2:
                file.write("#Osc. Amplitude \tGamma(eV) \tOmega(eV) \talpha \n")
                for i in range(self.calc.maxoscillators):
                    if float(self.calc.Amps[i].get()) != 0.0:
                        file.write("#"+ self.calc.Amps[i].get() + " \t"
                            + self.calc.Gammas[i].get() + " \t"
                            + self.calc.Omegas[i].get() + " \t"
                            + self.calc.Alphas[i].get() + "\t"
                            + self.calc.Us[i].get() + "\n" )
                    i = i + 1
            # elif self.calc.DFChoice.get() == 2:
                 # file.write("#Osc. Amplitude \tGamma(eV) \tOmega(eV)\t U(eV) \n")
                 # for i in range(self.calc.maxoscillators):
                    # if float(self.calc.Amps[i].get()) != 0.0:
                        # file.write( "#" + self.calc.Amps[i].get() + " \t"
                            # + self.calc.Gammas[i].get()  + " \t"
                            # + self.calc.Omegas[i].get()  + " \t"        
                            # + self.calc.Alphas[i].get()  + " \t"
                            # + self.calc.Us[i].get() + "\n" )
                    # i = i + 1
            else:    
                file.write("#Osc. Amplitude \tGamma(eV) \tOmega(eV) \tU(eV) \n")

                for i in range(self.calc.maxoscillators):
                    if float(self.calc.Amps[i].get()) != 0.0:
                        file.write("#" + self.calc.Amps[i].get() + " \t"
                            + self.calc.Gammas[i].get() + " \t"
                            + self.calc.Omegas[i].get() + " \t"
                            + self.calc.Us[i].get() + "\n" )
                    i = i + 1
            file.write("#GOS N elec \tedge(eV) \t 10 * n + l\t  atomic number Z \n")
            for i in range(self.calc.maxGOS):
                if float(self.calc.ConcGOS[i].get()) > 0.0:
                    GOSText = "#GOS \t" + self.calc.ConcGOS[i].get() + "\t" + \
                        self.calc.EdgeGOS[i].get() + " \t" + self.calc.nlGOS[i].get() + " \t" + self.calc.ZGOS[i].get() + "\t"
                    GOSText+= " density effect included up  to: "+ f"{self.calc.maxEnergyDensityEffect.get():.1f}"+ " (eV),\t"
                   
                    if self.calc.ApplySumRuleToGOS.get()== 1:
                        GOSText+= "rescaling on \n"
                    else:
                        GOSText+= "rescaling off\n"
                    file.write( GOSText)
                i = i + 1
            file.write("#Belkacem N elec \tomega(eV) \t Gamma \n")    
            for i in range(self.calc.maxBelkacem):
                if float(self.calc.Conc_Belkacem[i].get()) > 0.0:
                    file.write( "#Belkacem \t"
                        + self.calc.Conc_Belkacem[i].get() + " \t"
                        + self.calc.w_Belkacem[i].get() + " \t"
                        + self.calc.gamma_Belkacem[i].get() + "\n")
                i = i + 1
            file.write("#Kaneko N elec \t Q (a.u.) \t width (eV)\t gap U (eV)\t ang. mom l \t gamma Arista \n")    
            for i in range(self.calc.maxKaneko):
                if float(self.calc.N_Kaneko[i].get()) > 0.0:
                    KanekoText = "#Kaneko \t" + self.calc.N_Kaneko[i].get() + " \t" \
                        + self.calc.Q_Kaneko[i].get() + " \t" +self.calc.width_Kaneko[i].get() + "\t" \
                        + self.calc.Edge_Kaneko[i].get() + " \t" \
                        + self.calc.l_Kaneko[i].get() + " \t" + self.calc.gamma_Kaneko[i].get() + "\t"
                    if self.calc.Kaneko_choice.get() == 0:
                        KanekoText += " modified Kaneko\n"
                    else:
                        KanekoText += " original Kaneko\n"        
                    file.write( KanekoText)
                       
                i = i + 1
            file.write("#\n#==========end description dielectric function==========\n")
            
            
    def write_calculation_result(self):   
        with open(self.filename, "a", encoding="utf-8") as file:
            cmd=self.MyChapApp.plotchoice.get()
            file.write("#"+cmd)
            file.write("#"+self.runplot.plottitle)
          
            if cmd == "REELS":
                file.write("#projectile energy: " + str(self.calc.E0) + "keV\n")
                file.write("#Energy res: " + str(self.calc.Eres) + "eV\n")
                file.write( "#Theta in: "
                    + str(self.calc.thetaIn.get()) + ", Theta Out: "
                    + str(self.calc.thetaOut.get()) + "\n"
                )
                file.write( "#C1: " + str(self.calc.coef1.get())
                    + ", C2: "     + str(self.calc.coef2.get())
                    + ", C3: "     + str(self.calc.coef3.get()) + "\n"
                )
                file.write( "#surface excitation adjustment factor "
                    + str(self.calc.surf_ex_factor.get()) + "\n" )
                file.write("#Fraction DIIMFP in calculated energy range "
                    + str(self.calc.fraction_DIIMFP.get()) + "\n" )
                if self.calc.DSEP_choice.get() == 0:
                    file.write("#Dsep from global eps\n")
                else:
                    file.write("#Dsep per oscillator\n")


            file.write("\n#==========start calculation  output section==========\n")
            # Start Simon's patch 06/12/2020
            if cmd == "eq_plot" or cmd == "self_eq_plot" or cmd == "dcs_omega_eq_plot":
                Nqstep=int(self.calc.LastMomentum / self.calc.Stepsize_qplot)
                if cmd == "dcs_omega_eq_plot":
                    stepsize_b =  self.calc.thetamax.value / Nqstep /1000  # in mrad  
                    file.write("theta [mrad] \t omega (eV) \t DDCS\n")                
                elif cmd =="eq_plot" :
                    stepsize_b =  self.calc.Stepsize_qplot
                    file.write("q [A.U.] \t omega (eV) \t ELF\n")
                else:
                    stepsize_b =  self.calc.Stepsize_qplot
                    file.write("q [A.U.] \t omega (eV) \t SELF\n")    
                current_b = 0.5 *stepsize_b  # initial value
                    
                for i in range( Nqstep):
                    for j in range(self.calc.NPoints):
                        current_E = self.calc.CenterFirstBin + j * self.calc.Stepsize
                        a = "{:.3f} \t{:.4E} \t{:.4E}\n".format(
                            current_b, current_E, self.calc.my_image[j, i]
                        )
                        file.write(a)
                    current_b += stepsize_b
                # end Simon's patch 06/12/2020
                        
                # also write matrix for imagej
                self.save_in_anotherformat()
    
            elif cmd == "mean_ionization_energy":  # MIE calculation
                file.write( self.runplot.xlabel + "\t"
                    + self.runplot.label_curve1      + "\t"
                    + self.runplot.label_curve2      + "\t"
                    + self.runplot.label_curve3      + "\n" )
                for i in range(self.calc.NPoints):
                    a = "{:.3f} \t{:.4E} \t{:.4E}\t{:.4E}\n".format(
                        self.calc.x_axis[i],
                        self.calc.Result2[i],
                        self.calc.Result1[i],
                        self.calc.Result3[i],
                    )  # ordering historical
                    file.write(a)
            elif cmd == "sum_rules":  # sum rule calculations
                file.write(
                    self.runplot.xlabel
                    + "\t1/eps2 sum rule"
                    + "\teps2 sum rule"
                    + "\tk sum rule"
                    + "\tp sum rule"
                    + "\n"
                )
                for i in range(self.calc.NPoints):
                    a = "{:.3f} \t{:.4E} \t{:.4E}\t{:.4E} \t{:.4E}\n".format(
                        self.calc.x_axis[i],
                        self.calc.Result2[i],
                        self.calc.Result1[i],
                        self.calc.Result3[i],
                        self.calc.Result4[i],
                    )  # ordering historical
                    file.write(a)

            elif cmd == "IMFP_stop_strag" and self.calc.calculate_approximations:
                file.write(
                    "#equivalent parameter single oscillator for IMFP: C = "
                    + f"{self.calc.C0 :.5f}"
                    + ", W_p = "
                    + f"{self.calc.I0 :.1f}"
                    + "eV\n"
                )
                file.write(
                    "#equivalent parameter single oscillator for stopping: C = "
                    + f"{self.calc.C1:.5f}"
                    + ", W_p = "
                    + f"{self.calc.MIE :.1f}"
                    + "eV\n"
                )
                if np.amax(self.calc.DL_StragglingEnergy) > 0.0:
                    file.write(
                        "equivalent parameter single oscillator for straggling: C ="
                        + f"{self.calc.C2:.5f}"
                        + ", W_p = "
                        + f"{self.calc.Istraggling:.1f}"
                        + "eV\n"
                    )
                else:
                    file.write("#No meaningful single oscillator for straggling\n")
                file.write(
                    "#cross section per atom (angstrom^2) is (1.0/imfp)/"
                    + f"{self.calc.UnitCellDensity:.3f}"
                    + "\n"
                )

                file.write(
                    self.runplot.xlabel
                    + "\t"
                    + self.runplot.label_curve1
                    + "\t"
                    + self.runplot.label_curve2
                    + "\t"
                    + self.runplot.label_curve3)
                if self.calc.Approximations.get() == 1:    
                    file.write(    
                        "\tIMFP Bethe Disp. \tStop.  Bethe Disp.  \tStrag.  Bethe Disp.  "
                        + "\tDLIMFP average_comp \tDL stopping average_comp \tDL straggling average comp "
                        + "\tDL IMFP sum osc. \tDL stopping sum osc. \tDL strag. sum osc."
                        + "\tDL IMFP from ELF\tDL stopping from ELF \tDL strag. from ELF"
                        + "\t IMFP tpp2m\n"
                    )
                else:
                     file.write("\n")   

                for i in range(self.calc.NStopping):
                    if self.calc.axis_choice.get() == 1:
                        xvalue = self.calc.CurvesEnergy[i]
                    else:
                        xvalue = self.calc.CurvesVelocity[i]
                    if self.calc.Approximations.get() == 1:    
                        a = ("{:.4E}\t"*17).format(
                            xvalue,
                            self.calc.IMFPEnergy[i],
                            self.calc.StoppingEnergy[i]*self.calc.StoppingFactor,
                            self.calc.StragglingEnergy[i],
                            self.calc.BetheIMFPEnergy[i],
                            self.calc.BetheStoppingEnergy[i]*self.calc.StoppingFactor,
                            self.calc.BetheStragglingEnergy[i],
                            self.calc.DL_IMFPaverage_Energy[i],
                            self.calc.DL_StoppingEnergy[i]*self.calc.StoppingFactor,
                            self.calc.DL_StragglingEnergy[i],
                            self.calc.DL_IMFP_sum_Energy[i],
                            self.calc.DL_Stopping_sum_Energy[i]*self.calc.StoppingFactor,
                            self.calc.DL_Straggling_sum_Energy[i],
                            self.calc.DL_IMFP_from_ELF[i],
                            self.calc.DL_Stopping_from_ELF[i]*self.calc.StoppingFactor,
                            self.calc.DL_Straggling_from_ELF[i],
                            self.calc.TPP_IMFPEnergy[i]
                        )
                    else:
                        a = ("{:.4E}\t"*4).format(
                            xvalue,
                            self.calc.IMFPEnergy[i],
                            self.calc.StoppingEnergy[i]*self.calc.StoppingFactor,
                            self.calc.StragglingEnergy[i],
                        )    
                    file.write(a+"\n")

            elif cmd == "IMFP_stop_strag":  # stopping calculated, but no approximations
                file.write( self.runplot.xlabel + "\t"
                    + self.runplot.label_curve1 + "\t"
                    + self.runplot.label_curve2 + "\t"
                    + self.runplot.label_curve3 + "\n" )
                for i in range(self.calc.NStopping):
                    if self.calc.axis_choice.get() == 1:
                        xvalue = self.calc.CurvesEnergy[i]
                    else:
                        xvalue = self.calc.CurvesVelocity[i]

                    a = "{:.3f} \t{:.3f} \t{:.4E} \t{:.4E} \n".format(
                        xvalue,
                        self.calc.IMFPEnergy[i],
                        self.calc.StoppingEnergy[i],
                        self.calc.StragglingEnergy[i],
                    )

                    file.write(a)
            elif cmd == "partial_DIIMFP" or cmd == "partial_stopping": 
                a = self.runplot.xlabel
                for j in range(10):
                            a += "\t"+self.runplot.my10labels[j]
                file.write(a+"\n")  
                for i in range(self.calc.NPoints):
                        a ="{:.3f}".format(self.calc.x_axis[i])
                        for j in range(10):
                            a += "\t{:.4E}".format(self.calc.partialresults[i,j])
                        a += "\n"
                        file.write(a)          

            elif self.runplot.label_curve2 != "--":
                file.write(
                    self.runplot.xlabel + "\t" + self.runplot.label_curve1 + "\t" + self.runplot.label_curve2 + "\n"
                )
                if (self.MyChapApp.plotchoice.get() == "one_over_eps_q_direct") or \
                    (self.MyChapApp.plotchoice.get() == "eps_q_direct"):
                    for i in range(self.calc.NqPoints):
                        a = "{:.3f} \t{:.4E} \t{:.4E}\n".format(
                            self.calc.x_axis[i],
                            self.calc.Result1[i],
                            self.calc.Result2[i],
                        )
                        file.write(a)
                else:
                    for i in range(self.calc.NPoints):
                        a = "{:.3f} \t{:.4E} \t{:.4E}\n".format(
                            self.calc.x_axis[i],
                            self.calc.Result1[i],
                            self.calc.Result2[i],
                        )
                        file.write(a)

            elif self.runplot.label_curve1 != "--":
                file.write(self.runplot.xlabel + "\t" + self.runplot.label_curve1 + "\n")
                for i in range(self.calc.NPoints):
                    a = "{:.4e} \t{:.4E}\n".format(
                        self.calc.x_axis[i], self.calc.Result1[i]
                    )
                    file.write(a)
                    
    def save_in_anotherformat(self):
        filename = filedialog.asksaveasfilename(
            initialdir=os.getcwd,
            title="write file containing matrix for imagej",
            filetypes=(("dat files", "*.dat"), ("all files", "*.*")),
            defaultextension=".txt")
        if filename != '':  # i.e. dialog was not cancelled
            eqmax = np.amax(self.calc.my_image)
            image = (65534 * self.calc.my_image / eqmax)  # so output fits comfortably in 16 bit integer
            myintarray = image.astype(int)
            with open(filename, "w", encoding="utf-8") as file:
                file.seek(0)  # go to the beginning in case the file existed already
                np.savetxt(file, myintarray, delimiter=",", fmt="%i")

    def save_in_EQBetheformat(self):  # not sure where this routine is useful for (currently not called)

        filename = filedialog.asksaveasfilename(
            initialdir=os.getcwd,
            title="write file containing matrix for imagej",
            filetypes=(("dat files", "*.dat"), ("all files", "*.*")),
        )
        with open(filename, "w", encoding="utf-8") as file:
            file.seek(0)  # go to the beginning in case the file existed already
            np.savetxt(file, self.my_image, delimiter=",", fmt="%.4e")
            
    def load_comp_data(self):
        filename = filedialog.askopenfilename(
            initialdir=os.getcwd,
            title="Select file",
            filetypes=(("dat files", "*.dat"), ("txt files", "*.txt"),("all files", "*.*")),
        )
        if not filename:
            return
        count = len(open(filename).readlines())
        
      
        file = open(filename, "r")
        line = file.readline()
        number_comment_lines=1 # first line is not data but always assumed comment line
        for i in range(count - 1): # find the number of lines starting with a "#"
            line = file.readline()
            line.replace(","," ")
            line.replace("\t"," ")
            if line[0] == "#":
                number_comment_lines += 1
        print("number of comment lines",   number_comment_lines)    
        self.calc.Npointscomparray = count - number_comment_lines  
        self.calc.AllocateComparrays()      
           
    
        file.seek(0)   #"rewind"
        line = file.readline()
        print(line)
        idata=0
        for i in range(count - 1):
            line = file.readline()
            if line[0] == "#":
                print("comment line:",line)
            else:    
                line.replace(","," ")  # change a comma delimeted file to a space delimited.
                line.replace("\t"," ") ## change a tab delimeted file to a space delimited.
                nums = line.split()  
                if len(nums) < 2:
                    break
                self.calc.xCompArray[idata] = float(nums[0])
                self.calc.yCompArray[idata] = float(nums[1])
                idata +=1
                print("data:",nums[0],nums[1])
        file.close()
                    

             
