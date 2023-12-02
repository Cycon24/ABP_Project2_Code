import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dataClasses import CompressorStageData
from parameters import CompressorBaseParameters
from freeVortexCompressor import freeVortexCompressorMeanLine
from tools import getWDF
import warnings

'''
Class for a:
    Compressor with a Constant Axial Velocity Constant Mean Line

'''
class ConstantAxialVelocityConstantMeanLineCompressor(CompressorBaseParameters):
    root:CompressorStageData = None
    mean:CompressorStageData = None
    tip:CompressorStageData  = None

    def __init__(self, C_am, r_m, n_inf_c, m_dot, Delta_T_o, N, number) -> None:
        super().__init__()

        # setting class parameters
        self.C_am       = C_am
        self.r_m        = r_m
        self.n_inf_c    = n_inf_c
        self.U_m        = 2*np.pi*N*r_m # computing the tangential velocity at the mean line
        self.m_dot      = m_dot
        self.N          = N
        self.number     = number
        self.Delta_T_o  = Delta_T_o
    
    def calculate(self, T_o1, P_o1, Lambda=None, C_w1m=0):
        self.T_o1 = T_o1
        self.P_o1 = P_o1
        # Currently assuming Constant Axial Velocity

        # gets work done factor through interpolation
        wdf = getWDF(self.number)

        # calculate mean line props
        # if no degree of reaction is provided, find it
        if Lambda == None:
            # Solve for whirl difference and whirl vel
            Delta_C_wm = self.c_p*self.Delta_T_o/ (wdf*self.U_m)
            C_w2m = Delta_C_wm + C_w1m

            # Solve for Gas angles
            beta_1m  = np.arctan((self.U_m - C_w1m) / self.C_am)
            beta_2m  = np.arctan((self.U_m - C_w2m) / self.C_am)
            alpha_1m = np.arctan(C_w1m/self.C_am) # Should be 0 for first stage
            alpha_2m = np.arctan(C_w2m/self.C_am)
            
            # Calculate Lambda
            Lambda = 1 - (C_w1m + C_w2m)/(2*self.U_m)
        else: 
            # We already know lamda
            # Equations to help solve for betas
            eq1= (self.c_p*self.Delta_T_o) / (wdf*self.U_m*self.C_am) # = tan(B1) - tan(B2)
            eq2= 2*Lambda*self.U_m/self.C_am # = tan(B1) + tan(B2)
            
            # Solve for gas angles
            beta_1m  = np.arctan((eq1+eq2)/2)
            beta_2m  = np.arctan(eq2 - np.tan(beta_1m))
            alpha_1m = np.arctan(self.U_m/self.C_am  - np.tan(beta_1m))
            alpha_2m = np.arctan(self.U_m/self.C_am  - np.tan(beta_2m))

            # Solve for whirl velocities
            C_w1m = self.C_am*np.tan(alpha_1m)
            C_w2m = self.C_am*np.tan(alpha_2m)
            
        # The rest is the same whether first, last, or intermediate stage    
        # Solve for absolute gas velocities
        C_1m = np.sqrt(C_w1m**2 + self.C_am**2)
        C_2m = np.sqrt(C_w2m**2 + self.C_am**2)

        # get static props at mean line
        T_1m   = T_o1 - (C_1m**2)/(2*self.c_p)
        P_1m   = P_o1*(T_1m/T_o1)**(self.gamma/(self.gamma-1))
        rho_1m = P_1m/(self.R*T_1m)

        # Solve for blade dimensions
        self.h   = self.m_dot/(2*np.pi * rho_1m * self.C_am * self.r_m)
        self.r_t = self.r_m + self.h/2
        self.r_r = self.r_m - self.h/2

        self.A = np.pi*(self.r_t**2 - self.r_r**2)

        # create a free vortex blade and send it the properties at the mean line
        obj = freeVortexCompressorMeanLine(
            beta_1m, beta_2m, alpha_1m, alpha_2m, self.C_am, Lambda, self.r_m, T_o1, P_o1, self.N, self.Delta_T_o, wdf, self.n_inf_c
        )

        try:
            print("  Stage", self.number, end="\r")
            # recompute at the mean, to standardize the data collection
            self.mean = obj.calculate(self.r_m, "mean")

            # compute root and tip properties
            self.tip  = obj.calculate(self.r_t, "tip")
            self.root = obj.calculate(self.r_r, "root")

        except ValueError as e:
            # the "obj.calculate" calls will throw an error if the design criteria are violated
            print(f"While running {self.number}")
            print("  "+str(e))
            exit()

        # # Calcualte exit stag properties
        P_o3 = P_o1*(1 + self.n_inf_c*self.Delta_T_o/T_o1)**(self.gamma/(self.gamma - 1))
        T_o3 = T_o1 + self.Delta_T_o

        self.T_o3 = T_o3
        self.P_o3 = P_o3

        # returns stag props to update next stage
        return T_o3, P_o3

    # gets the data from the root mean and tip of the blade
    def getData(self):
        mean_data = self.mean.toDataFrame()
        tip_data  = self.tip.toDataFrame()
        root_data = self.root.toDataFrame()

        return pd.concat([root_data, mean_data, tip_data], keys=3*[f"Stage {self.number}"])


class DummyConstantAxialVelocityConstantMeanLineCompressor(CompressorBaseParameters):
    root:CompressorStageData = None
    mean:CompressorStageData = None
    tip:CompressorStageData  = None

    def __init__(self, C_am, r_m, m_dot, T_o1, P_o1) -> None:
        super().__init__()       

        self.C_am = C_am
        self.r_m = r_m
        self.m_dot = m_dot
        self.Delta_T_o = 0

        # calculate parameters
        self.calculate(T_o1, P_o1)

    def calculate(self, T_o1, P_o1):
        C_w1m = 0
        # The rest should be the same wether first, last, or intermediate stage    
        # Solve for absolute gas velocities
        C_1m = np.sqrt(C_w1m**2 + self.C_am**2)

        # get static props
        T_1m   = T_o1 - (C_1m**2)/(2*self.c_p)
        P_1m   = P_o1*(T_1m/T_o1)**(self.gamma/(self.gamma-1))
        rho_1m = P_1m/(self.R*T_1m)
        
        # Solve for blade dimensions
        self.h = self.m_dot/(2*np.pi * rho_1m * self.C_am * self.r_m)
        self.r_t = self.r_m + self.h/2
        self.r_r = self.r_m - self.h/2
        self.A = np.pi*(self.r_t**2 - self.r_r**2)

        # just setting to what is required in the combustor
        self.root = CompressorStageData({
            "alpha_1" : 0, 
            "C_1" : self.C_am, 
            "C_w1" : 0, 
            "r" : self.r_r 
        }, "root")
        self.tip = CompressorStageData({
            "alpha_1" : 0, 
            "C_1" : self.C_am, 
            "C_w1" : 0, 
            "r" : self.r_t
        }, "tip")
        self.mean = CompressorStageData({
            "alpha_1" : 0,
            "C_1" : self.C_am,
            "C_w1" : 0,
            "r" : self.r_m
        }, "mean")

    # returns the data
    def getData(self):
        mean_data = self.mean.toDataFrame()
        tip_data  = self.tip.toDataFrame()
        root_data = self.root.toDataFrame()

        return pd.concat([root_data, mean_data, tip_data], keys=3*[f"Dummy"])



############
#
# The rest is just data output, not important
#
############



def outputData(stage_data, stages):
    out = pd.concat(stage_data)
    for_formatting = CompressorStageData({}, "")
    with open("Compressor_Data.tex", 'w') as f:
        one = out[["alpha_1", "alpha_2", "beta_1", "beta_2", "C_w1", "C_w2", "C_1",  "C_2", "U"]]
        two = out[["V_w1",     "V_w2",      "V_1",    "V_2",   'M_1', 'M_2', 'Lambda', 'deHaller', 'r']]

        names = for_formatting.getFormattedColumns()
        one = one.rename(columns=names)
        two = two.rename(columns=names)
        
        f.write(
            "\\begin{center}\n" + one.to_latex(float_format='%.3f') + "\\end{center}\n\n" +
            "\\begin{center}\n" + two.to_latex(float_format='%.3f') + "\\end{center}\n\n"
        )

    out_T = out.T

    A = []
    h = []
    r_m = []
    names = []
    stags = []
    for i, stage in enumerate(stages):
        stags.append(pd.DataFrame({"$T_{o1}$" : stage.T_o1, "$P_{o1}$" : stage.P_o1, "$T_{o3}$" : stage.T_o3, "$P_{o3}$" : stage.P_o3, "$P_{o3}/P_{o1}$" : stage.P_o3/stage.P_o1}, index=[f"Stage {i+1}"]))

        A.append(stage.A)
        h.append(stage.h)
        r_m.append(stage.r_m)
        names.append(f"Stage {stage.number}")
    
    annulus_data = pd.DataFrame({"$A$" : A, "$h$" : h, "$r_m$" : r_m}, index=names)
    with open("Compressor_Annulus_Data.tex", 'w') as f:
        f.write(
            "\\begin{center}\n" + (annulus_data.T).to_latex(float_format='%.3f') + "\\end{center}\n\n"
        )
    
    stag_data = pd.concat(stags)

    with open("Compressor_Stagnation_Data.tex", "w") as f:
        f.write("\\begin{center}\n" + (stag_data).to_latex(float_format='%.3f') + "\\end{center}\n\n")

    parameters = CompressorBaseParameters()

    plt.rcParams["text.usetex"] = True

    plt.figure()
    for name in ["root", "mean", "tip"]:
        count = 1
        mach_data = []
        deHallar_data = []
        Lambda_data = []
        for item in out_T:
            if item[1] == name:
                mach_data.append(out_T[item[0]][name]["M_1"])
                count+=1
        plt.plot(list(range(1, count)), mach_data, label=name[0].upper() + name[1:])

    plt.plot(list(range(1, count)), (count-1)*[parameters.Mach_max], "--", label="Limit")
    plt.title(f"Mach Number along the Stages")
    plt.legend()
    plt.ylabel("Mach number")
    plt.xlabel("Stage")
    plt.grid(True)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(f"Compressor_Mach.png")

    plt.figure()
    for name in ["root", "mean", "tip"]:
        count = 1
        deHallar_data = []
        for item in out_T:
            if item[1] == name:
                deHallar_data.append(out_T[item[0]][name]["deHaller"])
                count+=1
                # print(out_T[item[0]][name]["M_1"])
        plt.plot(list(range(1, count)), deHallar_data, label=name[0].upper() + name[1:])

    plt.plot(list(range(1, count)), (count-1)*[parameters.de_Haller_min], "--", label="Limit")
    plt.title(f"de-Haller Number along the Stages")
    plt.legend()
    plt.ylabel("de-Haller Number")
    plt.xlabel("Stage")
    plt.grid(True)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(f"Compressor_deHaller.png")

    plt.figure()
    for name in ["root", "mean", "tip"]:
        count = 1
        Lambda_data = []
        for item in out_T:
            if item[1] == name:
                Lambda_data.append(out_T[item[0]][name]["Lambda"])
                count+=1
        plt.plot(list(range(1, count)), Lambda_data, label=name[0].upper() + name[1:])

    plt.title(f"Degree of Reaction $\\Lambda$ along the Stages")
    plt.legend()
    plt.ylabel("Degree of Reaction $\\Lambda$")
    plt.xlabel("Stage")
    plt.grid(True)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(f"Compressor_Lambda.png")


def drawCompressorGeometry(stage_data):
    warnings.filterwarnings("ignore")
    out = pd.concat(stage_data)
    out_T = out.T

    plt.rcParams["text.usetex"] = True

    r_rs = []
    for item in out_T:
        r_rs.append(out_T[item[0]]["root"]["r"])

    r_ts = []
    for item in out_T:
        r_ts.append(out_T[item[0]]["tip"]["r"])
    
    r_ms = []
    for item in out_T:
        r_ms.append(out_T[item[0]]["mean"]["r"])
    
    plt.figure()
    x = np.linspace(0, 0.6, len(r_rs))
    plt.plot(x, r_rs, "k")
    plt.plot(x, r_ts, "k")
    plt.plot(x, -np.array(r_rs), "k")
    plt.plot(x, -np.array(r_ts), "k")
    plt.plot(x, r_ms, "--b", label="Mean Line")
    plt.plot(x, -np.array(r_ms), "--b")

    plt.title(f"Representative Geometry of the Compressor")
    plt.legend()
    plt.ylabel("Radius (m)")
    plt.xlabel("Position")
    plt.grid(True)
    plt.minorticks_on()
    plt.gca().set_aspect("equal")
    plt.tight_layout()
    plt.savefig(f"Compressor_Geometry.png")
    warnings.filterwarnings("error")
