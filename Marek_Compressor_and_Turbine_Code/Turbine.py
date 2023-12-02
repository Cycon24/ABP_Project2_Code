import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings

from parameters import TurbineBaseParameters
from freeVortexTurbine import freeVortexCompressorMeanLine
from dataClasses import TurbineStageData


class ConstantAxialVelocityConstantMeanLineTurbine(TurbineBaseParameters):
    def __init__(self, C_am, r_m, m_dot, Delta_T_o, Lambda_m, phi_m, psi_m, N, number) -> None:
        super().__init__()
        # setting class parameters
        self.Lambda_m   = Lambda_m
        self.phi_m      = phi_m
        self.psi_m      = psi_m
        self.C_am       = C_am
        self.r_m        = r_m
        self.U_m        = 2*np.pi*N*r_m # computing the tangential velocity at the mean line
        self.m_dot      = m_dot
        self.N          = N
        self.number     = number
        self.Delta_T_o  = Delta_T_o
    
    def calculate(self, T_o1, P_o1):
        self.T_o1 = T_o1
        self.P_o1 = P_o1
        beta_3  = np.arctan(1/(2*self.phi_m)*(self.psi_m/2 + 2*self.Lambda_m))
        alpha_3 = np.arctan(np.tan(beta_3) - 1/self.phi_m) # rad

        C_w3m   = self.C_am*np.tan(alpha_3)

        C_3     = np.sqrt(self.C_am**2 + C_w3m**2)

        T_o3    = T_o1 - self.Delta_T_o
        PR      = (1 - self.Delta_T_o/(self.n_inf_t*T_o1))**(self.gamma/(self.gamma - 1))
        P_o3    = PR*P_o1

        T_3     = T_o3 - C_3**2 / (2*self.c_p) # K
        P_3     = P_o3*(T_3/T_o3)**(self.gamma/(self.gamma - 1))
        rho_3   = P_3/(self.R*T_3) # kg/m^3
        self.h  = self.m_dot / (2*np.pi*rho_3*self.C_am*self.r_m) # m

        self.r_t  = self.r_m + self.h/2 # m
        self.r_r  = self.r_m - self.h/2 # m
        self.A = np.pi*(self.r_t**2 - self.r_r**2)
        # print(self.r_r, self.r_t)
        # exit()

        obj = freeVortexCompressorMeanLine(
            self.C_am, self.Lambda_m, self.psi_m, self.phi_m, self.r_m, T_o1, self.N, self.Delta_T_o, P_o1
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

        self.T_o3 = T_o3
        self.P_o3 = P_o3

        return T_o3, P_o3

    # gets the data from the root mean and tip of the blade
    def getData(self):
        mean_data = self.mean.toDataFrame()
        tip_data  = self.tip.toDataFrame()
        root_data = self.root.toDataFrame()

        return pd.concat([root_data, mean_data, tip_data], keys=3*[f"Stage {self.number}"])


class DummyConstantAxialVelocityConstantMeanLineTurbine(TurbineBaseParameters):
    root:TurbineStageData = None
    mean:TurbineStageData = None
    tip:TurbineStageData  = None

    def __init__(self, C_am, r_m, m_dot, T_o1, P_o1) -> None:
        super().__init__()       

        self.C_am = C_am
        self.r_m = r_m
        self.m_dot = m_dot
        self.Delta_T_o = 0

        self.calculate(T_o1, P_o1)
    
    def calculate(self, T_o1, P_o1):
        self.T_o1 = T_o1
        self.P_o1 = P_o1
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

        self.root = TurbineStageData({
            "alpha_3" : 0, 
            "C_3" : self.C_am, 
            "C_w3" : 0, 
            "r" : self.r_r 
        }, "root")
        self.tip = TurbineStageData({
            "alpha_3" : 0, 
            "C_3" : self.C_am, 
            "C_w3" : 0, 
            "r" : self.r_t
        }, "tip")
        self.mean = TurbineStageData({
            "alpha_3" : 0,
            "C_3" : self.C_am,
            "C_w3" : 0,
            "r" : self.r_m
        }, "mean")
    
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
    for_formatting = TurbineStageData({}, "")
    with open("Turbine_Data.tex", 'w') as f:
        one = out[["beta_2", "beta_3", "alpha_2", "alpha_3", "C_w2", "C_w3", "C_2", "C_3"]] 
        two = out[["V_w2", "V_w3",      "V_2",      "V_3",    "U", "Lambda", "phi", "psi", "r"]]

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
    with open("Turbine_Annulus_Data.tex", 'w') as f:
        f.write(
            "\\begin{center}\n" + (annulus_data.T).to_latex(float_format='%.3f') + "\\end{center}\n\n"
        )
    
    stag_data = pd.concat(stags)

    with open("Turbine_Stagnation_Data.tex", "w") as f:
        f.write("\\begin{center}\n" + (stag_data).to_latex(float_format='%.3f') + "\\end{center}\n\n")


    parameters = TurbineBaseParameters()

    plt.rcParams["text.usetex"] = True

    plt.figure()
    for name in ["root", "mean", "tip"]:
        count = 1
        phi_data = []
        for item in out_T:
            if item[1] == name:
                phi_data.append(out_T[item[0]][name]["phi"])
                count+=1
        plt.plot(list(range(1, count)), phi_data, label=name[0].upper() + name[1:])

    plt.plot(list(range(1, count)), (count-1)*[parameters.phi_min], "--", label="Limit")
    plt.title(f"Flow Coefficient ($\\phi$) along the Stages")
    plt.legend()
    plt.ylabel("Flow Coefficient ($\\phi$)")
    plt.xlabel("Stage")
    plt.grid(True)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(f"Turbine_phi.png")

    plt.figure()
    for name in ["root", "mean", "tip"]:
        count = 1
        psi_data = []
        for item in out_T:
            if item[1] == name:
                psi_data.append(out_T[item[0]][name]["psi"])
                count+=1
                # print(out_T[item[0]][name]["M_1"])
        plt.plot(list(range(1, count)), psi_data, label=name[0].upper() + name[1:])

    plt.plot(list(range(1, count)), (count-1)*[parameters.psi_max], "--", label="Limit")
    plt.title(f"Blade Loading Coefficient ($\\psi$) along the Stages")
    plt.legend()
    plt.ylabel("Blade Loading Coefficient ($\\psi$)")
    plt.xlabel("Stage")
    plt.grid(True)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(f"Turbine_psi.png")

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
    plt.savefig(f"Turbine_Lambda.png")


def drawTurbineGeometry(stage_data):
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
    x = np.linspace(0, 0.1, len(r_rs))
    plt.plot(x, r_rs, "k")
    plt.plot(x, r_ts, "k")
    plt.plot(x, -np.array(r_rs), "k")
    plt.plot(x, -np.array(r_ts), "k")
    plt.plot(x, r_ms, "--b", label="Mean Line")
    plt.plot(x, -np.array(r_ms), "--b")

    plt.title(f"Representative Geometry of the Turbine")
    plt.legend()
    plt.ylabel("Radius (m)")
    plt.xlabel("Position")
    plt.xlim([-0.1, 0.2])
    plt.grid(True)
    plt.minorticks_on()
    plt.gca().set_aspect("equal")
    plt.tight_layout()
    plt.savefig(f"Turbine_Geometry.png")
    warnings.filterwarnings("error")






# Cw2t        = Cw2m * rm/rt # m/s
# C2t         = np.sqrt(Cw2t**2 + Ca**2) # m/s
# Ut2         = Um * rt/rm # m/s
# V2t         = np.sqrt((Cw2t - Ut2)**2 + Ca**2) # m/s
# beta2t      = np.arccos(Ca / V2t) # rad
# alpha2t     = np.arctan(Cw2t / Ca) # rad
# T2t         = To1 - C2t**2 / (2*cpg) # K
# M2t         = V2t / np.sqrt(gamma_g * R * T2t)