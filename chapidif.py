#!/usr/bin/env python3

import platform
import tkinter as tk
import os
# import tkinter.filedialog
#from tkinter import filedialog
from tkinter import ttk
from tkinter import font as tkFont

# try:
    # from icecream import ic

    # ic.configureOutput(includeContext=True)
# except ImportError:  # Graceful fallback if IceCream isn't installed.
    # ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a)  # noqa

import webbrowser as wb


import constants as cnst
from calculate import calculate
from run_and_plot import run_and_plot
from read_write_files import read_write_files

# this statement is required for the MAC and does not hurt other platforms
#matplotlib.use("TkAgg")
# plt.figure(figsize=(8, 6), dpi=80)

# class from https://gist.github.com/mp035/9f2027c3ef9172264532fcd6262f3b01
# ************************
# Scrollable Frame Class
# ************************


class ScrollFrame(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent)  # create a frame (self)

        self.canvas = tk.Canvas(self, borderwidth=2,width=CentralFrameWidth-30, height=500)  
        # place a frame on the canvas, this frame will hold the child widgets
        self.viewPort = tk.Frame(self.canvas)
        # place a scrollbar on self
        self.vsb = tk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        # attach scrollbar action to scroll of canvas
        self.canvas.configure(yscrollcommand=self.vsb.set)

        # pack scrollbar to right of self
        self.vsb.pack(side="right", fill="y")
        # pack canvas to left of self and expand to fil
        self.canvas.pack(side="left", fill="both", expand=True)
        self.canvas_window = self.canvas.create_window(
            (4, 4), window=self.viewPort, anchor="nw", tags="self.viewPort")

        # bind an event whenever the size of the viewPort frame changes.
        self.viewPort.bind("<Configure>", self.onFrameConfigure)
        # bind an event whenever the size of the viewPort frame changes.
        self.canvas.bind("<Configure>", self.onCanvasConfigure)

        # perform an initial stretch on render, otherwise the
        # scroll region has a tiny border until the first resize
        self.onFrameConfigure(None)

    def onFrameConfigure(self, event):
        """Reset the scroll region to encompass the inner frame"""
        self.canvas.configure(
            scrollregion=self.canvas.bbox("all")
        )  # whenever the size of the frame changes, alter the scroll region respectively.

    def onCanvasConfigure(self, event):
        """Reset the canvas window to encompass inner frame when required"""
        canvas_width = event.width
        # whenever the size of the canvas changes alter the window region respectively.
        self.canvas.itemconfig(self.canvas_window, width=canvas_width)


# class MainWindow(tk.Frame):
class Chapidif(tk.Frame):
    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)


    def UserInterface(self):
        
        small_font = tkFont.Font(family="Verdana", size=9)
        bold_font = tkFont.Font(family="Verdana", size=10, weight=tkFont.BOLD)
        self.master.grid_rowconfigure(0, weight=1)
        self.master.grid_columnconfigure(0, weight=1)
        self.commentline = ""
        self.calc = calculate(self)
        self.runplot = run_and_plot(self)
        self.runplot.screen_width=screen_width
        self.runplot.screen_height=screen_height
        self.readwrite= read_write_files(self)
        
        self.fname = tk.StringVar()
        self.fname.set("Comment for/from DF file:")
        
  
        self.last_calculation = 0
        self.plotchoice = tk.StringVar()
        self.plotchoice.set("one_over_eps_w_direct")
    
        self.BelkacemText = tk.StringVar()
        self.BelkacemText.set("Belkacem parameters (per unit cell)")
        self.KanekoText = tk.StringVar()
        self.KanekoText.set("Kaneko par. (per unit cell) (Mermin cor.)")
        
        self.master.title(
            "Chapidif: CHArged Particle Interactions from DIelectric Function, version 1.3"
        )
      
        # ----------------------dielectric function tab----------------------
        linespacing = 24
        LowerPanelHeight = 160
        UpperPanelHeight = 530
        Boundaryframe = tk.Frame(DFTab, highlightthickness=2,highlightbackground="blue")
        Boundaryframe.place(x=0, y=0, height=UpperPanelHeight, width= CentralFrameWidth )

        
        scrollFrame = ScrollFrame(Boundaryframe)  # add a new scrollable frame.
        scrollFrame.place(x=0,y=0)
        tk.Label(
            scrollFrame.viewPort, text="DF parameters:", justify=tk.LEFT, fg="red"
        ).grid(row=0, column=0, columnspan=3, sticky=tk.W)

        tk.Label(scrollFrame.viewPort, text="ₓ").grid(row=1, column=0)
        tk.Label(scrollFrame.viewPort, text="Aₓ").grid(row=1, column=1)
        tk.Label(scrollFrame.viewPort, text="ωₓ").grid(row=1, column=2)
        tk.Label(scrollFrame.viewPort, text="Γₓ").grid(row=1, column=3)
        tk.Label(scrollFrame.viewPort, text="αₓ").grid(row=1, column=4)
        tk.Label(scrollFrame.viewPort, text="Uₓ").grid(row=1, column=5)
        self.alpha_entrees=[]
        for i in range(self.calc.maxoscillators):
            tk.Label(scrollFrame.viewPort, text=str(i + 1), justify=tk.LEFT, font=small_font
            ).grid(row=i + 2, column=0)
            tk.Entry(scrollFrame.viewPort, width=10, font=small_font, 
                textvariable=self.calc.Amps[i] ).grid(row=i + 2, column=1)
            tk.Entry(scrollFrame.viewPort, width=10, font=small_font,
                textvariable=self.calc.Omegas[i]).grid(row=i + 2, column=2)
            tk.Entry(scrollFrame.viewPort, width=6, font=small_font,
                textvariable=self.calc.Gammas[i]).grid(row=i + 2, column=3)
            b=tk.Entry(scrollFrame.viewPort, width=4, font=small_font,
                textvariable=self.calc.Alphas[i])
            b.grid(row=i + 2, column=4) 
            self.alpha_entrees.append(b)  # keep the object so we can change its color          
            tk.Entry(scrollFrame.viewPort, width=4, font=small_font,
                textvariable=self.calc.Us[i]).grid(row=i + 2, column=5)
#===========dielectric function choice frame
      
        DFchoiceframe = tk.Frame(DFTab, highlightthickness=2,highlightbackground="blue")
        DFchoiceframe.place(x=5 , y=UpperPanelHeight+5, height=LowerPanelHeight-30+linespacing, width=200)
        
        ypos = 0
        tk.Label(DFchoiceframe, text="dielectric function model:", justify=tk.LEFT, fg="red"
        ).place(relx=0, y=ypos)
        
        ypos += linespacing
        tk.Frame(DFchoiceframe, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98)    
       
        ypos += 1
        tk.Radiobutton(DFchoiceframe, text="Extended Drude", variable=self.calc.DFChoice, value=1, highlightthickness=0.
            ).place(relx=0, y=ypos)
        
        ypos += linespacing
        tk.Radiobutton(DFchoiceframe, text="Drude-Lindhard (+U)", variable=self.calc.DFChoice, value=2,
            ).place(relx=0, y=ypos)
            
        ypos += linespacing
        tk.Radiobutton( DFchoiceframe,  variable=self.calc.DFChoice, value=3).place(relx=0, y=ypos)
        self.MerminLabel=tk.Label(DFchoiceframe, textvariable=self.calc.MerminText)
        self.MerminLabel.place(relx=0.15, y=ypos)
        # ypos += linespacing
        # tk.Radiobutton(DFchoiceframe, text="'Mermin-direct'", variable=self.calc.DFChoice, value=4,
            # ).place(relx=0, y=ypos)
        
        ypos += linespacing
        tk.Frame(DFchoiceframe, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98) 
        
        ypos += 1
        tk.Radiobutton(DFchoiceframe, text="add Chi", variable=self.calc.AddELF, value=0
            ).place(relx=0, y=ypos)
        tk.Radiobutton(DFchoiceframe, text="add ELF", variable=self.calc.AddELF, value=1
            ).place(relx=0.5, y=ypos)
     
 #==============gos frame========================================
        GOSFrameHeight=350
        self.GOSframe = tk.Frame(DFTab, highlightthickness=2,highlightbackground="blue")
        self.GOSframe.place(x=CentralFrameWidth+10, y=0, height=GOSFrameHeight, width= CentralFrameWidth )

        tk.Label(
            self.GOSframe,
            text="Hydrogenic GOS parameters (per unit cell):",
            justify=tk.LEFT,
            fg="red",
        ).grid(columnspan=5, row=0, column=0, sticky=tk.W)
        tk.Label(self.GOSframe, text="No").grid(row=1, column=0)
        tk.Label(self.GOSframe, text="N elec.").grid(row=1, column=1)
        tk.Label(self.GOSframe, text="Uₓ (eV)").grid(row=1, column=2)
        tk.Label(self.GOSframe, text="10*n+l").grid(row=1, column=3)
        tk.Label(self.GOSframe, text="Z").grid(row=1, column=4)

        for i in range(self.calc.maxGOS):  # Rows
            tk.Label(self.GOSframe, text=str(i + 1), justify=tk.LEFT, font=small_font).grid(row=i + 2, column=0)
            tk.Entry(self.GOSframe, width=7, font=small_font, textvariable=self.calc.ConcGOS[i]).grid(row=i + 2, column=1)
            tk.Entry(self.GOSframe, width=7, font=small_font, textvariable=self.calc.EdgeGOS[i]).grid(row=i + 2, column=2)
            tk.Entry(self.GOSframe, width=7, font=small_font, textvariable=self.calc.nlGOS[i]  ).grid(row=i + 2, column=3)
            tk.Entry(self.GOSframe, width=7, font=small_font, textvariable=self.calc.ZGOS[i]   ).grid(row=i + 2, column=4)

        tk.Radiobutton(self.GOSframe, text="GOS rescaling off", variable=self.calc.ApplySumRuleToGOS, value=0,
            ).grid(row=self.calc.maxGOS + 3, column=0, columnspan=3)
        tk.Radiobutton(self.GOSframe, text="GOS rescaling on",  variable=self.calc.ApplySumRuleToGOS, value=1,
            ).grid(row=self.calc.maxGOS + 3, column=3, columnspan=3)
        tk.Label(self.GOSframe, text="Apply Density correction up to (eV):",
            ).grid(row=self.calc.maxGOS + 4, column=0, columnspan=4)
        tk.Entry(self.GOSframe, width=5, font=small_font, textvariable=self.calc.maxEnergyDensityEffect,
            ).grid(row= self.calc.maxGOS + 4, column=4)

 #==============Belkacem frame========================================
        BelKacemFrameHeight = UpperPanelHeight -  GOSFrameHeight-5
        self.Belkacemframe = tk.Frame(DFTab, highlightthickness=2,highlightbackground="blue")
        self.Belkacemframe.place(x=CentralFrameWidth+10,y=GOSFrameHeight+5, height=BelKacemFrameHeight, width= CentralFrameWidth)

        self.BelkacemLabel = tk.Label(self.Belkacemframe, textvariable=self.BelkacemText, justify=tk.LEFT, fg="red")
        self.BelkacemLabel.grid(columnspan=4, row=0, column=0)
        tk.Label(self.Belkacemframe, text="No").grid(row=1, column=0)
        tk.Label(self.Belkacemframe, text="N elec.").grid(row=1, column=1)
        tk.Label(self.Belkacemframe, text="ωₓ").grid(row=1, column=2)
        tk.Label(self.Belkacemframe, text="Γₓ").grid(row=1, column=3)
     

        for i in range(self.calc.maxBelkacem):  # Rows
            tk.Label(self.Belkacemframe, text=str(i + 1), justify=tk.LEFT, font=small_font
                ).grid(row=i + 2, column=0)
            tk.Entry(self.Belkacemframe, width=6, font=small_font, textvariable=self.calc.Conc_Belkacem[i],
                ).grid(row=i + 2, column=1)
            tk.Entry(self.Belkacemframe, width=6, font=small_font, textvariable=self.calc.w_Belkacem[i],
                ).grid(row=i + 2, column=2)
            tk.Entry(self.Belkacemframe, width=6, font=small_font, textvariable=self.calc.gamma_Belkacem[i],
                ).grid(row=i + 2, column=3)
#=====================dispersion frame
    
        Dispersionframe = tk.Frame(DFTab, highlightthickness=2,highlightbackground="blue")
        Dispersionframe.place(
            x=225, y=UpperPanelHeight+5, height=LowerPanelHeight, width=400)
        ypos = 0
        tk.Label(
            Dispersionframe,
            text="RPA (Lindhard, Kaneko):",
            justify=tk.LEFT).place(relx=0.0, y=ypos)
        # tk.Checkbutton(Dispersionframe, command= self.update_Merminize,
            # variable=self.calc.Merminize,onvalue=1, offvalue=0,height=1, width=1, borderwidth=0,
            # ).place(relx=0.88, y=ypos)  
        tk.Radiobutton(Dispersionframe, text="Plain",command= self.update_Merminize,
            variable=self.calc.Merminize, value=0, borderwidth=0,pady=0).place(relx=0.5, y=ypos)  
        ypos += linespacing 
        tk.Radiobutton(Dispersionframe, text="Mermin corrected",command= self.update_Merminize,
            variable=self.calc.Merminize, value=1, borderwidth=0,pady=0).place(x=0.0, y=ypos)  
        tk.Radiobutton(Dispersionframe, text="Direct", command= self.update_Merminize,
            variable=self.calc.Merminize, value=2, borderwidth=0,pady=0).place(relx=0.5, y=ypos)      
              
        ypos += linespacing+2
        tk.Frame(Dispersionframe, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98) 
        ypos += 2
        tk.Label(
            Dispersionframe,
            text="dispersion (ext. Drude or Drude-Lindhard):",
            justify=tk.LEFT,
        ).place(relx=0.0, y=ypos)
        ypos += linespacing
        tk.Radiobutton(Dispersionframe, text="ω(q) = ω(0) + αQ",
            variable=self.calc.Dispersion_choice, value=0, borderwidth=0,pady=0).place(x=0.0, y=ypos)

        tk.Radiobutton(
            Dispersionframe, text="ω(q)² = ω(0)² + 2/3 α(v_f)² Q + Q²",
            variable=self.calc.Dispersion_choice, value=1, borderwidth=0,pady=0).place(x=148, y=ypos)
        ypos += linespacing+3
        tk.Frame(Dispersionframe, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98) 
        ypos += 4
        tk.Label(Dispersionframe, text="non-relativistic or relativistic dispersion:",
            justify=tk.LEFT).place(relx=0.0, y=ypos)
        ypos += linespacing
        tk.Radiobutton(
            Dispersionframe, text="Q = q²/2",
            variable=self.calc.Dispersion_relativistic,
            value=0).place(relx=0.0, y=ypos)

        tk.Radiobutton(
            Dispersionframe,
            text="Q = sqrt(c²q²+m²c⁴) - mc²",
            variable=self.calc.Dispersion_relativistic,
            value=1).place(x=148, y=ypos)
        
     
 ###################Kaneko frame     
        secondBoundaryframe = tk.Frame(DFTab, highlightthickness=2,highlightbackground="blue")
        secondBoundaryframe.place(x=2*(CentralFrameWidth+10), y=0, height=UpperPanelHeight, width= CentralFrameWidth )

        self.Kanekoframe = ScrollFrame(secondBoundaryframe)
        self.Kanekoframe.place(x=0,y=0)
        tk.Label(self.Kanekoframe.viewPort,textvariable= self.KanekoText,justify=tk.LEFT,fg="red",
            ).grid(columnspan=7, row=0, column=0)
        tk.Radiobutton(self.Kanekoframe.viewPort,text="Modiied \nKaneko", variable=self.calc.Kaneko_choice, value=0,
            ).grid(row=1, column=0, columnspan=3)
        tk.Radiobutton(self.Kanekoframe.viewPort,text="original \nKaneko", variable=self.calc.Kaneko_choice, value=1,
            ).grid(row=1, column=3, columnspan=3)
        tk.Label(self.Kanekoframe.viewPort, text="").grid(row=5, column=0)
        tk.Label(self.Kanekoframe.viewPort, text="Nₓ").grid(row=5, column=1)
        tk.Label(self.Kanekoframe.viewPort, text="Qₓ").grid(row=5, column=2)
        tk.Label(self.Kanekoframe.viewPort, text="Γₓ").grid(row=5, column=3)
        tk.Label(self.Kanekoframe.viewPort, text="lₓ").grid(row=5, column=4)
        tk.Label(self.Kanekoframe.viewPort, text="Uₓ").grid(row=5, column=5)
        tk.Label(self.Kanekoframe.viewPort, text="γₓ").grid(row=5, column=6)
        tk.Label(self.Kanekoframe.viewPort, text="ₓ").grid(row=6, column=0)
        tk.Label(self.Kanekoframe.viewPort, text="elec.").grid(row=6, column=1)
        tk.Label(self.Kanekoframe.viewPort, text="(a.u.)").grid(row=6, column=2)
        tk.Label(self.Kanekoframe.viewPort, text="(eV)").grid(row=6, column=3)
        tk.Label(self.Kanekoframe.viewPort, text="(<4)").grid(row=6, column=4)
        tk.Label(self.Kanekoframe.viewPort, text="(eV)").grid(row=6, column=5)


        for i in range(self.calc.maxKaneko):  # Rows
            tk.Label(self.Kanekoframe.viewPort, text=str(i + 1), justify=tk.LEFT, font=small_font,
                ).grid(row=i + 10, column=0)
            tk.Entry(self.Kanekoframe.viewPort, width=6, font=small_font, textvariable=self.calc.N_Kaneko[i],
                ).grid(row=i + 10, column=1)
            tk.Entry(self.Kanekoframe.viewPort, width=6, font=small_font, textvariable=self.calc.Q_Kaneko[i],
                ).grid(row=i + 10, column=2)
            tk.Entry(self.Kanekoframe.viewPort, width=6, font=small_font, textvariable=self.calc.width_Kaneko[i],
                ).grid(row=i + 10, column=3)
            tk.Entry(self.Kanekoframe.viewPort, width=2, font=small_font, textvariable=self.calc.l_Kaneko[i],
                ).grid(row=i + 10, column=4)
            tk.Entry(self.Kanekoframe.viewPort, width=6, font=small_font, textvariable=self.calc.Edge_Kaneko[i],
                ).grid(row=i + 10, column=5)
            tk.Entry(self.Kanekoframe.viewPort, width=6, font=small_font, textvariable=self.calc.gamma_Kaneko[i],
                ).grid(row=i + 10, column=6)
 
#=====================DF properties frame 
      
        OscillatorFrame = tk.Frame(DFTab, highlightthickness=2,highlightbackground="blue")
        OscillatorFrame.place( x=640,  y=UpperPanelHeight+5, height=LowerPanelHeight+20, width=430  )
        
        ypos = 0
        tk.Label(OscillatorFrame, text="Target:", justify=tk.LEFT, fg="red").place(relx=0, y=ypos)
        tk.Label(OscillatorFrame, text="Density (g/cm³):", justify=tk.LEFT).place(relx=0.15, y=ypos)
        tk.Entry(OscillatorFrame, width=5, textvariable=self.calc.specificweight
            ).place( relx=0.45, y=ypos)
        tk.Label(OscillatorFrame, text="Molar Mass (g):").place(relx=0.59, y=ypos)
        tk.Entry(OscillatorFrame, width=5, textvariable=self.calc.massunitcell).place(
            relx=0.87, y=ypos)
            
        ypos += linespacing+3   
        tk.Frame(OscillatorFrame, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98)
        
        ypos += 3    
        tk.Label(OscillatorFrame, textvariable=self.calc.UnitCellDensityText).place(x=0, y=ypos)
        
        ypos += linespacing    
        tk.Label( OscillatorFrame, textvariable=self.calc.plasmon_1elec_per_unit_cell).place(x=0, y=ypos)
        
        ypos += linespacing
        self.DFLabel=tk.Label(OscillatorFrame, textvariable=self.calc.DF_prop_text, justify=tk.LEFT )
        self.DFLabel.place(x=0, y=ypos) 
        
        ypos += linespacing
        tk.Label(OscillatorFrame, textvariable=self.calc.GOSdenstext  , justify=tk.LEFT).place(x=0, y=ypos)
        
        ypos += linespacing
        tk.Label(OscillatorFrame, textvariable=self.calc.Belkacemtext , justify=tk.LEFT).place(x=0, y=ypos)
        
        ypos += linespacing
        tk.Label(OscillatorFrame, textvariable=self.calc.Kanekodenstext,justify=tk.LEFT).place(x=0, y=ypos)
        
        self.buttonreadDF = tk.Button( DFTab, text="read DF from file",  command=self.readwrite.read_df)
        self.buttonreadDF.place(x=10, y=UpperPanelHeight + LowerPanelHeight )
        self.buttonwriteDF = tk.Button(DFTab, text="write DF to file",   command=self.readwrite.write_df_only ) 
        self.buttonwriteDF.place(x=10, y=UpperPanelHeight + LowerPanelHeight+38)
        self.buttonzeroDF = tk.Button(DFTab, text="zero DF ",   command=self.readwrite.zero_df ) 
        self.buttonzeroDF.place(x=140, y=UpperPanelHeight + LowerPanelHeight+38)
        
        tk.Label(DFTab, textvariable=self.fname).place(x=240, rely=0.92)#y=UpperPanelHeight + LowerPanelHeight + 20 )
       
        self.CommentEntry = tk.Entry(DFTab, width=60, font=small_font)
        self.CommentEntry.place(x=240,rely=0.95)  # y=UpperPanelHeight + LowerPanelHeight +20 + linespacing)
        tk.Button(DFTab, text="Calculate and Plot",command=self.runplot.start_calc_and_plot, fg="red", 
            bg="blue", font=bold_font).place(relx=.64, rely=0.95)
        tk.Button(DFTab, text="replot", command=self.runplot.replot).place(relx=.78, rely=0.95)
        tk.Button(DFTab, text="close plots", command=self.runplot.Close_all_plots).place(relx=0.85,rely=0.95) 
        tk.Button(DFTab, text="help", command=self.Help).place(relx=0.94,rely=0.95)
      
            
             
# ====================begin  Action Tab===============
        CurrentHeight = 0
        ypos = 0
        linespacing = 25
        PlotchoiceframeHeight = 700
        FirstColWidth = 430
        SecondColWidth =330
        ThirdColWidth = 350

        Plotchoiceframe = tk.Frame(ActionTab,  highlightthickness=2,highlightbackground="blue")
        Plotchoiceframe.place(x=5,y=5, height=PlotchoiceframeHeight, width= FirstColWidth)

        tk.Label(Plotchoiceframe, text="PLOT CHOICE:", justify=tk.LEFT, fg="red").place(relx=0,y=ypos)
        ypos += linespacing + 2

        tk.Label(Plotchoiceframe, text="Dielectric function", justify=tk.LEFT, fg="red").place(relx=0,y=ypos)
        ypos += linespacing
        tk.Radiobutton(
            Plotchoiceframe, text="Re[ϵ(ω)],Im[ϵ(ω)]",
            variable=self.plotchoice, value="eps_w_direct", highlightthickness=0).place(relx=0,y=ypos)
        ypos += linespacing
        tk.Radiobutton(Plotchoiceframe, text="Re[1/ϵ(ω)], Im[-1/ϵ(ω)] ",
            variable=self.plotchoice,value="one_over_eps_w_direct",).place(relx=0,y=ypos)
        ypos += linespacing    
        tk.Radiobutton(Plotchoiceframe,text="Re[ϵ(q)],Im[ϵ(q)]",
            variable=self.plotchoice,value="eps_q_direct").place(relx=0,y=ypos)
        ypos += linespacing
        tk.Radiobutton(Plotchoiceframe,text="Re[1/ϵ(q)], Im[-1/ϵ(q)] ",
            variable=self.plotchoice,value="one_over_eps_q_direct").place(relx=0,y=ypos)
        ypos += linespacing
        tk.Radiobutton(Plotchoiceframe, text="color plot  Im[-1/ϵ(q,ω)]",variable=self.plotchoice,
            value="eq_plot").place(relx=0,y=ypos)
        ypos += linespacing      
        tk.Radiobutton(Plotchoiceframe, text="n(ω), k(ω)", variable=self.plotchoice, value="n_k").place(relx=0,y=ypos)    
        ypos += linespacing + 2

        tk.Label(Plotchoiceframe, text="Kramers-Kronig tests", justify=tk.LEFT, fg="red").place(relx=0,y=ypos)
        ypos += linespacing-1
        tk.Radiobutton(
            Plotchoiceframe,text="Re[ϵ(ω)],Im[ϵ(ω)]",variable=self.plotchoice,
            value="eps_w_kk").place(relx=0,y=ypos)
        ypos += linespacing
        tk.Radiobutton(Plotchoiceframe,text="Re[1/ϵ(ω)], Im[-1/ϵ(ω)] ",
            variable=self.plotchoice,value="one_over_eps_w_kk").place(relx=0,y=ypos)
        ypos += linespacing
        tk.Radiobutton(Plotchoiceframe,text="n(ω), k(ω)",
            variable=self.plotchoice,value="n_k_kk").place(relx=0,y=ypos)
            #============================================================
        ypos += linespacing
      # ttk.Separator(Plotchoiceframe, orient=tk.HORIZONTAL).place(relx=0,y=ypos,relwidth=0.5)
        ypos += 2
        tk.Label(Plotchoiceframe, text="Others:", justify=tk.LEFT, fg="red").place(relx=0,y=ypos)
        ypos += linespacing-1
 
        tk.Radiobutton(Plotchoiceframe, text="Oscillator Strength", variable=self.plotchoice,
            value="oscillator_strength").place(relx=0,y=ypos)
        ypos += linespacing      
        tk.Radiobutton(Plotchoiceframe, text="dyn. struct. factor S(k,ω)", variable=self.plotchoice,
            value="dyn_struct_factor").place(relx=0.,y=ypos)    
        ypos += linespacing
       
        tk.Radiobutton(Plotchoiceframe, text="sum rules", variable=self.plotchoice, value="sum_rules"
            ).place(relx=0.0,y=ypos)
        ypos += linespacing      
        tk.Radiobutton(Plotchoiceframe, text="inertial sum rules", variable=self.plotchoice, value="inertial_rules"
            ).place(relx=0.,y=ypos)
        ypos += linespacing  
        tk.Radiobutton(Plotchoiceframe, text="S(k,ω) sum rule", variable=self.plotchoice, value="S_k_omega_rule"
            ).place(relx=0.,y=ypos)    
        ypos += linespacing 
        tk.Radiobutton(Plotchoiceframe, text="mean ionization energy", variable=self.plotchoice,
            value="mean_ionization_energy").place(relx=0,y=ypos)
        ypos += linespacing  
        tk.Radiobutton(Plotchoiceframe, text="Compton profile", variable=self.plotchoice, 
            value="Compton").place(x=0,y=ypos)
        ypos += linespacing  
        tk.Label(Plotchoiceframe, text="(at q = ", justify=tk.LEFT).place(relx=.05,y=ypos)     
        self.Compton_q = tk.Entry(Plotchoiceframe, width=5)
        self.Compton_q.place(relx=.17,y=ypos+2)
        self.Compton_q.insert(0, self.calc.q_Compton) 
        tk.Label(Plotchoiceframe, text="a.u.)", justify=tk.LEFT).place(relx=.30,y=ypos) 
        ypos += linespacing 
         
        ypos += 2  
        tk.Label(Plotchoiceframe, text="Local density approximation:", justify=tk.LEFT, fg="red").place(relx=0.0,y=ypos)
        ypos += linespacing 
        tk.Radiobutton(Plotchoiceframe, text="pseudo charge-density", variable=self.plotchoice, 
            value="pseudo_charge_density").place(relx=0,y=ypos)
        ypos += linespacing   
        tk.Radiobutton(Plotchoiceframe, text="radial pseudo charge-dens.", variable=self.plotchoice, 
            value="radial_charge_density").place(relx=0,y=ypos)
        ypos += linespacing 
        tk.Radiobutton(Plotchoiceframe, text="ωₚ(r), λ(r)", variable=self.plotchoice, 
            value= "radial_w_p_lambda").place(relx=0,y=ypos)
            
#===============new column            
        ypos = 0+2*linespacing     
 
        tk.Label(Plotchoiceframe, text="Quantities involving", justify=tk.LEFT, fg="red").place(relx=0.5,y=ypos)
        ypos += linespacing-4
        tk.Label(Plotchoiceframe, text="momentum integration", justify=tk.LEFT, fg="red").place(relx=0.5,y=ypos)
        ypos += linespacing
        tk.Radiobutton(Plotchoiceframe, text="diimfp", variable=self.plotchoice, value="diimfp"
            ).place(relx=0.5,y=ypos)
       
        tk.Radiobutton(Plotchoiceframe, text="partial diimfp", variable=self.plotchoice, 
            value="partial_DIIMFP").place(relx=0.7,y=ypos)
        ypos += linespacing
        # tk.Radiobutton(Plotchoiceframe, text="partial diimfp weighted by  ω", variable=self.plotchoice, 
            # value="partial_stopping").place(relx=0.5,y=ypos) 
        tk.Radiobutton(Plotchoiceframe, text="partial diimfp weighted by  ω", variable=self.plotchoice, 
             value="partial_stopping").place(relx=0.5,y=ypos) 
        ypos += linespacing  
        tk.Radiobutton(Plotchoiceframe, text="energy dependence of \nimfp,stopping,straggling", variable=self.plotchoice,
            value="IMFP_stop_strag").place(relx=0.5,y=ypos)
        ypos += 2*linespacing        

        tk.Radiobutton(Plotchoiceframe, text="L (stopping number)", variable=self.plotchoice,
            value="shell_effect_all").place(relx=0.5,y=ypos)
        ypos += linespacing
        # tk.Radiobutton(Plotchoiceframe, text="Shell effect soft collisions", variable=self.plotchoice,
            # value="shell_effect_soft").place(relx=0.5,y=ypos)
        # ypos += linespacing
      # ttk.Separator(Plotchoiceframe, orient=tk.HORIZONTAL).place(relx=0,y=ypos,relwidth=0.5)
        ypos +=  linespacing
        tk.Label(Plotchoiceframe, text="Theta plots:", justify=tk.LEFT, fg="red").place(relx=0.52,y=ypos)
        ypos += linespacing  
        tk.Label(Plotchoiceframe, text="( up to θ = ", justify=tk.LEFT).place(relx=.52,y=ypos)
        self.theta_max_entree = tk.Entry(Plotchoiceframe, width=4)
        self.theta_max_entree.place(relx=.7,y=ypos-2)
        self.theta_max_entree.insert(0, self.calc.theta_max)   
        tk.Label(Plotchoiceframe, text="mrad)", justify=tk.LEFT).place(relx=.770,y=ypos)
        ypos += linespacing     
       
        tk.Radiobutton(Plotchoiceframe, text="dσ/(dωdΩ) at ω=", variable=self.plotchoice, 
            value="DDCS_at_omega").place(relx=0.5,y=ypos)
        self.omega_ddcs_entree = tk.Entry(Plotchoiceframe, width=4)
        self.omega_ddcs_entree.place(relx=.84,y=ypos-2)
        self.omega_ddcs_entree.insert(0, self.calc.omega_ddcs)
        tk.Label(Plotchoiceframe, text="eV", justify=tk.LEFT).place(relx=.92,y=ypos)   
        ypos += linespacing             
        tk.Radiobutton(Plotchoiceframe, text="dσ/(dωdΩ) at θ=", variable=self.plotchoice, 
            value="DDCS_at_theta").place(relx=0.5,y=ypos)
        self.theta_ddcs_entree = tk.Entry(Plotchoiceframe, width=4)
        self.theta_ddcs_entree.place(relx=.82,y=ypos-2)
        self.theta_ddcs_entree.insert(0, self.calc.theta_ddcs)
        tk.Label(Plotchoiceframe, text="mrad", justify=tk.LEFT).place(relx=.9,y=ypos)       

       
           
        ypos += linespacing        
        tk.Radiobutton(Plotchoiceframe, text="dcs  dσ/(dΩ)",
            variable=self.plotchoice, value="dcs_plot").place(relx=0.5,y=ypos)  
        ypos += linespacing  
        tk.Radiobutton(Plotchoiceframe, text="2d-plot  dσ/(dωdΩ), \nno retardation",
            variable=self.plotchoice, value="dcs_omega_eq_plot").place(relx=0.5,y=ypos)  
        ypos +=2*linespacing       
        tk.Radiobutton(Plotchoiceframe, text="2d-plot  dσ/(dωdΩ), \nincl. retardation",
            variable=self.plotchoice, value="Cerenkov").place(relx=0.5,y=ypos)  
        ypos +=2*linespacing       
        tk.Radiobutton(Plotchoiceframe, text="retarded - not-retarded",
            variable=self.plotchoice, value="difference_due_to_Cerenkov").place(relx=0.5,y=ypos)      
        ypos += 2*linespacing+10  
        tk.Label(Plotchoiceframe, text="Electron specific:", justify=tk.LEFT, fg="red").place(relx=0.52,y=ypos)
        ypos += linespacing
        tk.Radiobutton(Plotchoiceframe, text="SELF, DSEP", variable=self.plotchoice, 
            value="SELF_DSEP").place(relx=0.52,y=ypos)
        ypos += linespacing
        tk.Radiobutton(Plotchoiceframe, text="color plot SELF", variable=self.plotchoice, 
            value="self_eq_plot").place(relx=0.52,y=ypos)
        ypos += linespacing
        tk.Radiobutton(Plotchoiceframe, text="(R)EELS spectrum", variable=self.plotchoice, 
            value="REELS").place(relx=0.52,y=ypos)
        ypos += 2*linespacing
        # tk.Label(Plotchoiceframe, text="Proton specific:", justify=tk.LEFT, fg="red").place(relx=0.52,y=ypos)
        # ypos += linespacing
        # tk.Radiobutton(Plotchoiceframe, text="DDCS at 0 degrees",
            # variable=self.plotchoice, value="bin_enc_peak").place(relx=0.5,y=ypos)      

################### plot parameters        
        PlotParframeheight= 320
        PlotParframe = tk.Frame(ActionTab,  highlightthickness=2,highlightbackground="blue")
        PlotParframe.place(x=  FirstColWidth+10,y=5, height=PlotParframeheight, width=SecondColWidth)
      
        ypos=0
        tk.Label(PlotParframe, text= "Energy plots:", fg='red').place(relx=0,y=ypos)
        
        ypos += linespacing
        tk.Label(PlotParframe, text="from ").place(x=0,y=ypos)
        self.FE = tk.Entry(PlotParframe, width=5)
        self.FE.place(relx=0.13,y=ypos)
        self.FE.insert(0, self.calc.LowerELimit)
        tk.Label(PlotParframe, text="eV  to", justify=tk.LEFT).place(relx=.29,y=ypos)
        self.LE = tk.Entry(PlotParframe, width=5)
        self.LE.place(relx=0.43,y=ypos)
        self.LE.insert(0, self.calc.LastEnergy)
        tk.Label(PlotParframe, text="eV,  step ").place(relx=.58,y=ypos)
        # tk.Label(PlotParframe, text="step").place(relx=0,y=ypos)
        self.step = tk.Entry(PlotParframe, width=4)
        self.step.place(relx=.78,y=ypos)
        self.step.insert(0, self.calc.Stepsize)
        tk.Label(PlotParframe, text="eV").place(relx=.92,y=ypos)
        
        ypos += linespacing
        tk.Label(PlotParframe, text="at q =").place(x=0,y=ypos)
        self.qused = tk.Entry(PlotParframe, width=5)
        self.qused.place(relx=.13,y=ypos)
        self.qused.insert(0, self.calc.q)
        tk.Label(PlotParframe, text="a.u.",).place(relx=.29,y=ypos)
        
        ypos += linespacing
        tk.Label(PlotParframe, text="Momentum plots:", justify=tk.LEFT, fg="red"
        ).place(relx=0.0,y=ypos)
        
        ypos += linespacing
        tk.Label(PlotParframe, text="from 0.0 a.u. to q =  ", justify=tk.LEFT).place(x=0.0,y=ypos)
        self.LM = tk.Entry(PlotParframe, width=5)
        self.LM.place(relx=.42,y=ypos)
        self.LM.insert(0, self.calc.LastMomentum)
        tk.Label(PlotParframe, text="a.u,  step", justify=tk.LEFT).place(relx=.56,y=ypos)
        self.stepq = tk.Entry(PlotParframe, width=4)
        self.stepq.place(relx=.76,y=ypos)
        self.stepq.insert(0, self.calc.Stepsize_qplot)
        tk.Label(PlotParframe, text="a.u.",).place(relx=.90,y=ypos)
        
        ypos += linespacing
        tk.Label(PlotParframe, text="at ω =", justify=tk.LEFT).place(x=0,y=ypos)
        self.omega_used = tk.Entry(PlotParframe, width=5)
        self.omega_used.place(x=42,y=ypos)
        self.omega_used.insert(0, self.calc.Energy_qplot)
        tk.Label(PlotParframe, text="eV", justify=tk.LEFT).place(x=70,y=ypos)
        ypos += linespacing+4
        tk.Frame( PlotParframe, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98)
        
        ypos += +4
        tk.Label( PlotParframe, text="energy scale:").place(x=0,y=ypos)

        tk.Radiobutton(PlotParframe, text="eV", variable=self.calc.Energy_Scale_choice,
            value=0).place(x=90,y=ypos)
        tk.Radiobutton(PlotParframe, text="nm⁻¹",variable=self.calc.Energy_Scale_choice,
            value=1,).place(relx=0.63,y=ypos)
            
        ypos += linespacing    
        self.LogX=tk.IntVar()
        self.LogY=tk.IntVar()
        tk.Checkbutton(PlotParframe, text = "log scale x", variable = self.LogX).place(x=0,y=ypos)
        tk.Checkbutton(PlotParframe, text = "log scale y", variable = self.LogY).place(relx=0.63,y=ypos)
        
        ypos += linespacing
        tk.Checkbutton(PlotParframe, text = "log  color plots, range:", 
            variable = self.calc.log_choice).place(relx=0,y=ypos)
        self.log_range= tk.Entry(PlotParframe, width=2)
        self.log_range.place(relx=.65,y=ypos+3)
        self.log_range.insert(0, self.calc.log_range)
        ypos += linespacing
        tk.Label(PlotParframe, text="max. color plot (autoscale=0)").place(x=0,y=ypos)
        self.Max_eq= tk.Entry(PlotParframe, width=4)
        self.Max_eq.place(relx=.65,y=ypos)
        self.Max_eq.insert(0, self.calc.max_eq)
        ypos += linespacing+3
        tk.Frame( PlotParframe, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98)
        ypos += 2
         
        # tk.Label( PlotParframe, text="plot IMFP or Crosssection/Unit Cell:").place(relx=0,y=ypos)
        # ypos += linespacing 
        # tk.Radiobutton(PlotParframe, text="IMFP", variable=self.calc.imfp_crosssection_choice,
            # value=0).place(relx=0.,y=ypos)
        # tk.Radiobutton(PlotParframe, text="Crosssection/U.C.",variable=self.calc.imfp_crosssection_choice,
            # value=1,).place(x=80,y=ypos)
 
        REELSFrameHeight = 220
        REELSFrame = tk.Frame(ActionTab, highlightthickness=2,highlightbackground="blue")
        REELSFrame.place(x= FirstColWidth+10, y=PlotParframeheight + 10
            , height=REELSFrameHeight,width=SecondColWidth)
            
        ypos = 0
        tk.Label(REELSFrame, text="EELS parameters:", justify=tk.LEFT, fg="red").place(
            x=0, y=ypos)
        tk.Label(REELSFrame, text="Energy res. (eV):", justify=tk.LEFT).place(
            relx=.45, y=ypos)
        tk.Entry(REELSFrame, width=4,textvariable=self.calc.Eres).place(relx= .85, y=ypos)


        ypos = ypos + linespacing
        tk.Label(REELSFrame, text="θ₀:").place(relx=0., y=ypos)
        tk.Entry(REELSFrame, width=4, textvariable=self.calc.thetaIn ).place(relx=0.08, y=ypos)
        tk.Label(REELSFrame, text="θ₁:").place(relx=0.3, y=ypos)
        tk.Entry(REELSFrame, width=4, textvariable=self.calc.thetaOut).place(relx=0.37, y=ypos)

        ypos = ypos + linespacing
        tk.Label(
            REELSFrame, text="Fraction DIIMFP in energy loss range:", justify=tk.LEFT
        ).place(x=0, y=ypos)

        tk.Entry(REELSFrame, width=4, textvariable=self.calc.fraction_DIIMFP).place(relx=.8, y=ypos)

        ypos += linespacing
        tk.Label( REELSFrame, text="surf. plasmon scaling factor:", justify=tk.LEFT).place(x=0, y=ypos)
        tk.Entry(REELSFrame, width=4, textvariable= self.calc.surf_ex_factor).place(relx=.8, y=ypos)
        
        ypos = ypos + linespacing
        tk.Radiobutton(REELSFrame, text="DSEP from total eps", variable=self.calc.DSEP_choice,
            value=0, padx=0).place(x=0, y=ypos)

        tk.Radiobutton(REELSFrame, text="DSEP per oscil.", variable=self.calc.DSEP_choice,
            value=1, padx=0).place(relx=.55, y=ypos)
        
        ypos += linespacing+4
        tk.Frame( REELSFrame, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98)
        ypos += 4
        tk.Radiobutton(REELSFrame, text="Reflection:", fg='red', variable=self.calc.EELS,
            value=0, padx=0).place(x=0, y=ypos)
        tk.Label( REELSFrame,text="Part.Int(i)=1+c₁ i + c₂ i² +c₃ i³",
        ).place(relx=.4, y=ypos)
        
        ypos += linespacing
        tk.Label(REELSFrame, text="c₁:").place(x=0, y=ypos)
        tk.Entry(REELSFrame, width=6, textvariable=self.calc.coef1).place(relx=0.08, y=ypos)
        tk.Label(REELSFrame, text="c₂:").place(relx=0.29, y=ypos)
        tk.Entry(REELSFrame, width=6, textvariable=self.calc.coef2).place(relx=0.36, y=ypos)
        tk.Label(REELSFrame, text="c₃:").place(relx=0.6, y=ypos)
        tk.Entry(REELSFrame, width=6, textvariable=self.calc.coef3).place(relx=0.68, y=ypos)

        ypos += linespacing+4
        tk.Frame( REELSFrame, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98)
        ypos += 4
        tk.Radiobutton(REELSFrame, text="Transmission:", fg='red', variable=self.calc.EELS,
            value=1, padx=0).place(x=0, y=ypos)
        tk.Label(REELSFrame, text="Thickness (Å)").place(relx=0.48, y=ypos)  
        tk.Entry(REELSFrame, width=6, textvariable=self.calc.EELS_thickness).place(relx=.8,y=ypos) 

        PlotEtcFrameHeight = 150
        PlotEtcFrame = tk.Frame(ActionTab, highlightthickness=2,highlightbackground="blue")
        PlotEtcFrame.place(x= FirstColWidth+10, y=PlotParframeheight+REELSFrameHeight+15
            , height=PlotEtcFrameHeight,width=SecondColWidth)
        ypos = 3
        self.buttonOpenCompData = tk.Button(
            PlotEtcFrame,text="Load comp. data",command=self.readwrite.load_comp_data)
        self.buttonOpenCompData.place(relx=.02, y=ypos)
        C1 = tk.Checkbutton(PlotEtcFrame, text="plot comp. data",  variable=self.calc.Overplot,
            onvalue=1, offvalue=0, height=1, width=12, pady=0, borderwidth = 0 )
        C1.place(relx=.5, y=ypos+5 )
        
        ypos += 1.5*linespacing
        tk.Label(PlotEtcFrame, text="legend comp. data", justify=tk.LEFT).place(x=0,y=ypos)
        tk.Entry(PlotEtcFrame, font=small_font, textvariable=self.runplot.LiteratureDescription, width=20).place(relx=.41, y=ypos)
        ypos += linespacing
        tk.Label(PlotEtcFrame, text="plot with:", justify=tk.LEFT).place(x=0,y=ypos)
        tk.Radiobutton(PlotEtcFrame, text="imfp", variable=self.calc.comp_option_choice, value=0
        ).place(relx=0.19, y=ypos)
        tk.Radiobutton(PlotEtcFrame, text="stopping", variable=self.calc.comp_option_choice, value=1
        ).place(relx=0.40, y=ypos)
        tk.Radiobutton(PlotEtcFrame, text="straggling", variable=self.calc.comp_option_choice, value=2
        ).place(relx=0.68, y=ypos)
         
        ypos += linespacing
        tk.Label(PlotEtcFrame, text="description of this calculation:", justify=tk.LEFT, fg="red").place(x=0,y=ypos)
        ypos += linespacing    
        tk.Entry(PlotEtcFrame, font=small_font, textvariable=self.runplot.PlotDescription, width=35).place(relx=.01, y=ypos)



        addvaluesframeHeight = 120
        currentHeight = 5
        addvaluesframe = tk.Frame(ActionTab,highlightthickness=2,highlightbackground="blue")
        addvaluesframe.place(x=FirstColWidth+SecondColWidth+15, y=currentHeight,height=addvaluesframeHeight
            ,width=ThirdColWidth  )
        currentHeight += addvaluesframeHeight +5
        ypos = 0
        tk.Label(addvaluesframe, text="Projectile:", fg="red").place(relx=0, y=ypos)
        tk.Radiobutton(
            addvaluesframe, text="e⁻", variable=self.calc.Projectile, value=0
        ).place(relx=0.5, y=ypos)
        tk.Radiobutton(
            addvaluesframe, text="H⁺", variable=self.calc.Projectile, value=1).place(relx=0.7, y=ypos)
        ypos += linespacing
        tk.Label(addvaluesframe, text="E₀ (keV):").place(relx=0.0, y=ypos)
        self.E0entry = tk.Entry(addvaluesframe,textvariable=self.calc.ProjectileEnergy, width=7)
        self.E0entry.place(relx=0.20, y=ypos)
        tk.Label(addvaluesframe, textvariable=self.calc.myvelocitytext).place(relx=0.43, y=ypos)
        ypos += linespacing+4
        tk.Frame(addvaluesframe, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98)
        ypos += +4
        tk.Label(addvaluesframe, text="H⁺: incl. Mott correction").place(relx=0.0, y=ypos)
       
        tk.Checkbutton(addvaluesframe,variable=self.calc.MottCorrection,onvalue=1, offvalue=0,height=1, 
            width=1, borderwidth=0).place(relx=0.45, y=ypos)  
      
       
        ypos += linespacing    
        tk.Label(addvaluesframe, text="e⁻: incl. exchange eff.").place(relx=0.0, y=ypos)
        tk.Checkbutton(addvaluesframe,variable=self.calc.ExchangeCorrection, onvalue=1, offvalue=0,height=1, 
            width=1, borderwidth=0).place(relx=0.4, y=ypos)  
        tk.Label(addvaluesframe, text="B.E. for exch.:").place(relx=0.5, y=ypos)
        self.Be_for_exchangeEntry = tk.Entry(addvaluesframe,width=5)
        self.Be_for_exchangeEntry.insert(0, self.calc.BE_for_exchange) 
        self.Be_for_exchangeEntry.place(relx=0.78, y=ypos) 
        tk.Label(addvaluesframe, text="eV").place(relx=0.93, y=ypos)
      #  ypos += linespacing
        # tk.Radiobutton(
            # addvaluesframe, text="None", variable=self.calc.ExchangeCorrection, value=0
        # ).place(relx=0.336, y=ypos)
        # tk.Radiobutton(
            # addvaluesframe,
            # text="valence only",
            # variable=self.calc.ExchangeCorrection,
            # value=1,
        # ).place(relx=0.52, y=ypos)
        # tk.Radiobutton(
            # addvaluesframe, text="all", variable=self.calc.ExchangeCorrection, value=2
        # ).place(relx=0.85, y=ypos)
   
        StoppingCalcFrameHeight = 250
        StoppingCalcFrame = tk.Frame(ActionTab, highlightthickness=2,highlightbackground="blue")
        StoppingCalcFrame.place(x= FirstColWidth+ SecondColWidth+15, y= currentHeight,height=StoppingCalcFrameHeight
            ,width=ThirdColWidth)
        currentHeight += StoppingCalcFrameHeight+5    
       
        ypos = 0
        tk.Label(StoppingCalcFrame, text="DIIMFP integration:").place(x=0, y=ypos)
       # ypos += linespacing
        tk.Radiobutton(StoppingCalcFrame, text="rough"  ,variable=self.calc.Stopping_calc_quality, value=0,
            ).place(relx=0.35, y=ypos)
        tk.Radiobutton(StoppingCalcFrame, text="medium" ,variable=self.calc.Stopping_calc_quality, value=1,
            ).place(relx=0.57, y=ypos)
        tk.Radiobutton(StoppingCalcFrame, text="fine",variable=self.calc.Stopping_calc_quality, value=2,
            ).place(relx=0.82, y=ypos)
        ypos += linespacing 
        tk.Label(StoppingCalcFrame, text="plot:").place(x=0, y=ypos)   
        tk.Radiobutton(StoppingCalcFrame, text="IMFP", variable=self.calc.imfp_crosssection_choice,
            value=0).place(relx=0.25,y=ypos)
        tk.Radiobutton(StoppingCalcFrame, text="Crosssection/U.C.",variable=self.calc.imfp_crosssection_choice,
            value=1,).place(relx=0.5,y=ypos)
        ypos += linespacing
        tk.Label(StoppingCalcFrame, text="IMPF, Stop., Strag. versus:").place(x=0, y=ypos)
        tk.Radiobutton(StoppingCalcFrame, text="velocity", variable=self.calc.axis_choice, value=2,
            ).place(relx=0.50, y=ypos)
        tk.Radiobutton(StoppingCalcFrame, text="Energy", variable=self.calc.axis_choice, value=1,
            ).place(relx=0.76, y=ypos)    
            
        ypos += linespacing+4
        tk.Frame( StoppingCalcFrame, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98)
        
        ypos += 4
        tk.Label(StoppingCalcFrame, text="stopping conversion:", justify=tk.LEFT, fg="red",
            ).place(x=0, y=ypos)
            
        ypos += linespacing
        tk.Label(StoppingCalcFrame, textvariable=self.calc.myconversiontext1).place(x=0, y=ypos)
        
        ypos += linespacing
        tk.Label(StoppingCalcFrame, textvariable=self.calc.myconversiontext2).place(x=0, y=ypos)
        ypos += linespacing
        tk.Label(StoppingCalcFrame, text= "e⁻: add estimate radiative stopping").place(relx=0.0, y=ypos)
        tk.Checkbutton(StoppingCalcFrame,variable=self.calc.RadiativeLosses,onvalue=1, offvalue=0,height=1, 
            width=1, borderwidth=0).place(relx=0.85, y=ypos)  
        
        ypos += linespacing+4
        tk.Frame( StoppingCalcFrame, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98)
        
        ypos += 4
        tk.Label(StoppingCalcFrame, text="plot  approximations:").place( x=0, y=ypos)
        tk.Radiobutton(StoppingCalcFrame, text="yes", variable=self.calc.Approximations, value=1,
            ).place(relx=0.6, y=ypos)
        tk.Radiobutton(StoppingCalcFrame, text="no", variable=self.calc.Approximations, value=0,
            ).place(relx=0.8,y=ypos)
            
        ypos += linespacing
        tk.Label(StoppingCalcFrame, text="ωₚ  for TPP (eV):").place(relx=0, y=ypos)
        tk.Entry(StoppingCalcFrame, width=4, textvariable= self.calc.w_p_TPP).place(relx=0.7, y=ypos)

        
        
          
        TransformationFrameHeight=175
        TransformationFrame = tk.Frame(ActionTab, highlightthickness=2,highlightbackground="blue")
        TransformationFrame.place(x=FirstColWidth+ SecondColWidth+15, 
            y=  currentHeight , height=TransformationFrameHeight, width=ThirdColWidth)
        currentHeight+= TransformationFrameHeight+5  
        
        ypos=0
        tk.Label(TransformationFrame, text="DF transformations:", fg="red").place(relx=0.0, y=ypos)
        tk.Label(TransformationFrame, text="(using", ).place(relx=0.4, y=ypos)
        tk.Entry(TransformationFrame, textvariable = self.calc.N_oscillator_used, width=3).place(relx=0.53, y=ypos)   
        tk.Label(TransformationFrame, text="Oscillators)", ).place(relx=0.62, y=ypos)
        ypos += linespacing 
      
        tk.Button(TransformationFrame,  text="DF from OOS",command=self.readwrite.DF_from_OOS, pady=0).place(x=5, y=ypos)
        ypos += 1.2*linespacing 
        tk.Button( TransformationFrame, text="'Penn'",command=self.Convert_to_Penn, pady=0).place(x=5, y=ypos)
        tk.Label(TransformationFrame, text="max. value ωₚ (eV) Penn:" ).place(relx=0.3, y=ypos)
      
        tk.Entry(TransformationFrame, textvariable = self.calc.w_p_Penn, width=4).place(relx=0.75, y=ypos-3)
    
        ypos += 1.2*linespacing
        tk.Frame(TransformationFrame, bg="blue", height=1, bd=0).place(relx=0.01, y=ypos, relwidth=.98)
        ypos += .5*linespacing
        tk.Button(TransformationFrame, text="DL/Mermin to 'MLL'", command=self.Convert_to_MLL, pady=0).place(x=5, y=ypos)
        tk.Label(TransformationFrame, text="with U=" ).place(relx=0.5, y=ypos)
        tk.Entry(TransformationFrame, textvariable = self.calc.U_MLL_transform, width=5).place(relx=0.65, y=ypos-3) 
        tk.Label(TransformationFrame, text="eV" ).place(relx=0.8, y=ypos)
        ypos += 1.5*linespacing
        tk.Button(TransformationFrame, text=" DL/Mermin to Kaneko",  command=self.Pennify_Kaneko, pady=0).place(x=5, y=ypos)
       
        tk.Label(TransformationFrame, text=" l =" ).place(relx=0.5, y=ypos)
        tk.Entry(TransformationFrame, textvariable = self.calc.l_Kaneko_transform, width=1).place(relx=0.56, y=ypos - 3)
        tk.Label(TransformationFrame, text=", Q =" ).place(relx=0.6, y=ypos)
        tk.Entry(TransformationFrame, textvariable = self.calc.Q_Kaneko_transform, width=4).place(relx=0.7, y=ypos - 3)
        tk.Label(TransformationFrame, text="a.u." ).place(relx=0.8, y=ypos)
   
        ypos = 4     
        LastFrameHeight=140
        LastFrame = tk.Frame(ActionTab, highlightthickness=2,highlightbackground="blue")
        LastFrame.place(x=FirstColWidth+SecondColWidth+15, y= currentHeight, 
            height=LastFrameHeight, width=ThirdColWidth)
            
        tk.Button(LastFrame, text="save last calculation",command=self.readwrite.save_calculation).place(relx=0.01, y=ypos)  
      #  tk.Button(LastFrame, text="replot calculation",command=self.runplot.replot).place(relx=0.6, y=ypos)        
            
        ypos += 1.5*linespacing    
        tk.Button(LastFrame, text="close plots", command=self.runplot.Close_all_plots).place(relx=0.01, y=ypos) 
    #    tk.Button(LastFrame, text="help", command=self.Help).place(relx=0.6, y=ypos) 
        ypos += linespacing
        self.ErrorLabel=tk.Label(LastFrame, textvariable=self.calc.ErrorMessage)
        self.ErrorLabel.config(fg='red', font ='bold')    
        self.ErrorLabel.place(relx=0, y=ypos)
        #tk.Button(ActionTab, text="replot calculation",command=self.runplot.replot).place(relx=0.7, rely=0.95)  
        tk.Button(ActionTab, text="Calculate and Plot",command=self.runplot.start_calc_and_plot, fg="red", 
            bg="blue", font=bold_font).place(relx=.64, rely=0.95)
        tk.Button(ActionTab, text="replot", command=self.runplot.replot).place(relx=.78, rely=0.95)
        tk.Button(ActionTab, text="close plots", command=self.runplot.Close_all_plots).place(relx=0.85,rely=0.95) 
        tk.Button(ActionTab, text="help", command=self.Help).place(relx=0.94,rely=0.95)
#====================begin  DefaultsTab===============
        CurrentHeight = 0
        ypos = 0
        linespacing = 25
        DefaultsframeHeight = 700
        FirstColWidth = 390
        SecondColWidth =390
        #ThirdColWidth = 300 

        Defaultsframe = tk.Frame(DefaultsTab,  highlightthickness=2,highlightbackground="blue")
        Defaultsframe.place(x=5,y=5, height=DefaultsframeHeight, width= FirstColWidth)

        tk.Label( Defaultsframe, text="Change some default calculation paramaeters here:", justify=tk.LEFT, fg="red").place(relx=0,y=ypos)
        ypos += linespacing
        tk.Label( Defaultsframe, text="NStopping", justify=tk.LEFT, ).place(relx=0,y=ypos)
        self.NStopping_entree = tk.Entry( Defaultsframe, width=4)
        self.NStopping_entree.place(relx=.68,y=ypos)
        self.NStopping_entree.insert(0, self.calc.NStopping)  
        ypos += linespacing
        tk.Label( Defaultsframe, text="Energy increment factor (> 1)", justify=tk.LEFT, ).place(relx=0,y=ypos)
        self.IncrFactor_entree = tk.Entry( Defaultsframe, width=4)
        self.IncrFactor_entree.place(relx=.68,y=ypos)
        self.IncrFactor_entree.insert(0, self.calc.IncrFactor)  
        ypos += linespacing
        tk.Label( Defaultsframe, text="first energy e^- (keV)", justify=tk.LEFT, ).place(relx=0,y=ypos)
        self.first_E_e_entree = tk.Entry( Defaultsframe, width=7)
        self.first_E_e_entree.place(relx=.68,y=ypos)
        self.first_E_e_entree.insert(0, self.calc.first_electron_energy)  
        ypos += linespacing
        tk.Label( Defaultsframe, text="first energy H^+ (keV)", justify=tk.LEFT, ).place(relx=0,y=ypos)
        self.first_E_H_entree = tk.Entry( Defaultsframe, width=7)
        self.first_E_H_entree.place(relx=.68,y=ypos)
        self.first_E_H_entree.insert(0, self.calc.first_proton_energy) 
        ypos += 5*linespacing
        tk.Label( Defaultsframe, text="Replace  Mermin by DL  q*ω <", justify=tk.LEFT, ).place(relx=0,y=ypos)
        self.q_transition_entree = tk.Entry( Defaultsframe, width=7)
        self.q_transition_entree.place(relx=.68,y=ypos)
        self.q_transition_entree.insert(0, self.calc.q_transition) 
        ypos += linespacing
        tk.Label( Defaultsframe, text="value UMax (= omega/qQ, in a.u.):", justify=tk.LEFT, ).place(relx=0,y=ypos)
        self.UMax_entree = tk.Entry( Defaultsframe, width=7)
        self.UMax_entree.place(relx=.68,y=ypos)
        self.UMax_entree.insert(0, self.calc.UMax) 
        
        ypos += 2*linespacing
        tk.Label( Defaultsframe, text="electrons calculate exchange used, as in:", justify=tk.LEFT, ).place(relx=0,y=ypos)
        ypos += linespacing
        tk.Radiobutton( Defaultsframe, text="Ashley", variable=self.calc.ExchangeCorrectionMethod, value=0).place(relx=.0, y=ypos)
        tk.Radiobutton(Defaultsframe, text="Salvat (SBethe)", variable=self.calc.ExchangeCorrectionMethod, value=1).place(relx=.5, y=ypos)
        ypos += 2*linespacing
        tk.Label( Defaultsframe, text="Stopping units", justify=tk.LEFT, ).place(relx=0,y=ypos)
        ypos += linespacing
        tk.Radiobutton( Defaultsframe, text="eV/Å", variable=self.calc.StoppingUnits, value=0).place(relx=.0, y=ypos)
        tk.Radiobutton(Defaultsframe, text="eV / (10¹⁵ atoms/cm²)", variable=self.calc.StoppingUnits, value=1).place(relx=.17, y=ypos)
        tk.Radiobutton(Defaultsframe, text="MeV / (mg/cm²)", variable=self.calc.StoppingUnits, value=2).place(relx=.63, y=ypos)
       
        
        
        ypos=0
        PlotParamFrame = tk.Frame(DefaultsTab,  highlightthickness=2,highlightbackground="blue")
        PlotParamFrame.place(x=FirstColWidth+ 15,y=5, height=DefaultsframeHeight, width= SecondColWidth)
        tk.Label( PlotParamFrame, text="Modify plot appearance:", justify=tk.LEFT, fg="red").place(relx=0,y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="figure width:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Entry(PlotParamFrame, width=3,  textvariable=self.runplot.figurewidth).place(relx=0.4,y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="figure height:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Entry(PlotParamFrame, width=3,  textvariable=self.runplot.figureheight).place(relx=0.4,y=ypos)
        
        ypos += linespacing
        tk.Label( PlotParamFrame, text="Display grid:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Radiobutton( PlotParamFrame, text="yes", variable=self.runplot.DisplayGrid, value=True).place(relx=.3, y=ypos)
        tk.Radiobutton( PlotParamFrame, text="No", variable=self.runplot.DisplayGrid, value=False).place(relx=.5, y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="Title font size:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Entry(PlotParamFrame, width=3,  textvariable=self.runplot.title_fontsize).place(relx=0.4,y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="Legend font size:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Entry(PlotParamFrame, width=3,  textvariable=self.runplot.legendfontsize).place(relx=0.4,y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="Legend position:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Radiobutton( PlotParamFrame, text="best", variable=self.runplot.legendposition, value=0).place(relx=.25, y=ypos)
        tk.Radiobutton( PlotParamFrame, text="upper right", variable=self.runplot.legendposition, value=1).place(relx=.5, y=ypos)
        tk.Radiobutton( PlotParamFrame, text="upper left", variable=self.runplot.legendposition, value=2).place(relx=.75, y=ypos)
        ypos += linespacing
        tk.Radiobutton( PlotParamFrame, text="lower left", variable=self.runplot.legendposition, value=3).place(relx=.0, y=ypos)
        tk.Radiobutton( PlotParamFrame, text="lower right", variable=self.runplot.legendposition, value=4).place(relx=.25, y=ypos)
        tk.Radiobutton( PlotParamFrame, text="right", variable=self.runplot.legendposition, value=5).place(relx=.5, y=ypos)
        tk.Radiobutton( PlotParamFrame, text="left", variable=self.runplot.legendposition, value=6).place(relx=.75, y=ypos)
        ypos += linespacing
        tk.Radiobutton( PlotParamFrame, text="lower center",variable=self.runplot.legendposition, value=8).place(relx=.0, y=ypos)
        tk.Radiobutton( PlotParamFrame, text="upper center", variable=self.runplot.legendposition, value=9).place(relx=.25, y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="axis label font size:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Entry(PlotParamFrame, width=3,  textvariable=self.runplot.axislabelsize).place(relx=0.4,y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="x label font size:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Entry(PlotParamFrame, width=3,  textvariable=self.runplot.xlabelsize).place(relx=0.4,y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="y label font size:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Entry(PlotParamFrame, width=3,  textvariable=self.runplot.ylabelsize).place(relx=0.4,y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="marker size:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Entry(PlotParamFrame, width=3,  textvariable=self.runplot.markersize).place(relx=0.4,y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="line width:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Entry(PlotParamFrame, width=3,  textvariable=self.runplot.linewidth).place(relx=0.4,y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="Comparison Data:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Radiobutton( PlotParamFrame, text="Line + Marker", variable=self.runplot.CompDataLine, value=True).place(relx=.3, y=ypos)
        tk.Radiobutton( PlotParamFrame, text="Marker only", variable=self.runplot.CompDataLine, value=False).place(relx=.65, y=ypos)
        ypos += linespacing
        tk.Label( PlotParamFrame, text="default file format figures:", justify=tk.LEFT).place(relx=0,y=ypos)
        tk.Radiobutton( PlotParamFrame, text="pdf", variable=self.runplot.fileformat, value="pdf").place(relx=.45, y=ypos)
        tk.Radiobutton( PlotParamFrame, text="png", variable=self.runplot.fileformat, value="png").place(relx=.65, y=ypos)
        # ypos += linespacing
        # tk.Button(PlotParamFrame, text="replot calculation",command=self.runplot.replot).place(relx=0,y=ypos)  
        ypos += 3*linespacing
        T = tk.Text(PlotParamFrame, height = 15, width = 43) #height and width in characters and line units
      
        T.place(relx=0.05,y=ypos)
        xpos=FirstColWidth + SecondColWidth+50
        ypos=25
        tk.Label(DefaultsTab, text="debug mode", justify=tk.LEFT, fg="red").place(x=xpos, y=ypos)
        tk.Radiobutton(DefaultsTab, text="off", variable=self.calc.DebugMode, value=0).place(x=xpos+100,  y=ypos)
        tk.Radiobutton(DefaultsTab, text="on", variable=self.calc.DebugMode, value=1).place(x=xpos+150, y=ypos)
       
        my_text="""keymap from matplotlib documentation
toggle fullscreen:       f, ctrl+f    
home (reset):            h, r, home        
back:                    left arrow, c
forward:                 right arrow , v
pan:                     p 
zoom:                    o 
toggle major grid:       g  
toggle minor grid:       G   
toggle log-linear y:     l      
toggle log-linear x:     L,k          
save:                    s, ctrl+s
close:                   ctrl+w """

        T.insert(tk.END, my_text)  
        T.config(state=tk.DISABLED)  
        tk.Button(DefaultsTab, text="Calculate and Plot",command=self.runplot.start_calc_and_plot, fg="red", 
            bg="blue", font=bold_font).place(relx=.64, rely=0.95) 
        tk.Button(DefaultsTab, text="replot", command=self.runplot.replot).place(relx=.78, rely=0.95)
        tk.Button(DefaultsTab, text="close plots", command=self.runplot.Close_all_plots).place(relx=0.85,rely=0.95) 
        tk.Button(DefaultsTab, text="help", command=self.Help).place(relx=0.94,rely=0.95) 
           
                    

        
    def UpdateColors(self):
        try:
            if self.calc.SumAi > 1.0:
                self.DFLabel.config(fg='red')
                self.alpha_entrees[2].config(fg='red')
            else:
                self.DFLabel.config(fg='black')
                self.alpha_entrees[2].config(fg='black')
            self.calc.Stepsize =    float(self.step.get())  
            self.calc.CenterFirstBin = float(self.FE.get())  +  0.5* self.calc.Stepsize
            if self.calc.CenterFirstBin <= 0.0:
                self.FE.config(fg='red')
            else:    
                self.FE.config(fg='black')
           
            if float(self.FE.get()) != 0.0:
                self.FE.config(fg= 'orange')
            else:
                self.FE.config(fg= 'black')    
            for i in range(self.calc.maxoscillators):  
                if float(self.calc.Alphas[i].get()) != 1.0:
                    self.alpha_entrees[i].config(fg='orange')
                else:
                    self.alpha_entrees[i].config(fg='black')
        except ValueError:
            return
        except AttributeError:  # this means calc not yet created
            return
           
        except Exception as err:
            print(f"Unexpected {err=}, {type(err)=}")
            #raise    
                
           
        
    def update_Merminize(self):
        if self.calc.Merminize.get() == 0:   
             self.calc.MerminText.set("Lindhard or LL")
             self.KanekoText.set("Kaneko par. (per unit cell)")    
        elif self.calc.Merminize.get() == 1:   
             self.calc.MerminText.set("Mermin or MLL")
             self.KanekoText.set("Kaneko par. (per unit cell) (Mermin cor.)")
        elif self.calc.Merminize.get() == 2:   
             self.calc.MerminText.set("RPA Direct")
             self.KanekoText.set("Kaneko par. (per unit cell) (Direct method")
             
            
    def Help(self):
        print(os.path.dirname(os.path.abspath(__file__)))
        a=os.path.join(os.getcwd(), "chapidif_manual_plus_background.pdf")
        wb.open_new("file:" + a)

    # -----------------------------end right panel-------------------------------------

    # ========================end user interface=======================================

 
    def Convert_to_Penn(self):
        self.runplot.initcalc()
        self.calc.Penn_from_ELF()

    def Pennify_Kaneko(self):
        self.runplot.initcalc()
        self.calc.Stepsize = float(self.step.get())
        self.calc.CenterFirstBin = float(self.FE.get())+ 0.5*self.calc.Stepsize 
        #self.calc.Pennified_Kaneko_from_ELF()
        self.calc.convert_DL_to_Kaneko()

    def Convert_to_MLL(self):
        self.runplot.initcalc()
        self.calc.convert_to_MLL()

        
    def Cleanup_and_exit(self):
        self.runplot.Close_all_plots()
        root.destroy()
   
    
#================================================================================================

  
if __name__ == "__main__":
    root = tk.Tk()  # create a program derived from tk
    root.title("Chapidif")
    MainFrameHeight = 800
    MainFrameWidth = 1130
    LeftFrameWidth = 280
    CentralFrameWidth = 350 # also width for Rightframe
    tabControl = ttk.Notebook(root, height=MainFrameHeight+15)
    DFTab = ttk.Frame(tabControl)
    ActionTab = ttk.Frame(tabControl)
    DefaultsTab = ttk.Frame(tabControl)
    tabControl.add(DFTab, text='Define dielectric Funtion')
    tabControl.add(ActionTab, text='Perform calculation')
    tabControl.add(DefaultsTab, text='change defaults')
 
    RightFrameWidth = MainFrameWidth - LeftFrameWidth - CentralFrameWidth

    screen_width = root.winfo_screenwidth() 
    screen_height = root.winfo_screenheight()    
    left_offset=screen_width - MainFrameWidth-25
    root.geometry(f"{MainFrameWidth}x{MainFrameHeight}+{left_offset}+0")
    tabControl.pack(expand = 1, fill ="both")
     
   
    MyUserInterface = Chapidif()  # create the class that provides the user interface, 
                          #input output etc in it, a class is created that does the calculations
    MyUserInterface.pack(side="top", fill="both", expand=True)  # put my userinterface in the main frame
   
    MyUserInterface.UserInterface()  # calls the user interface
    root.protocol("WM_DELETE_WINDOW",MyUserInterface.Cleanup_and_exit)    
    root.mainloop()  # makes  the program interact with user input
    
    
      
       # atexit.register(Chapidif.runplot.Close_all_plots)    
