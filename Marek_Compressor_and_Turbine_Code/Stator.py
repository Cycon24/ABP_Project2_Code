import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from parameters import CompressorBaseParameters
from dataClasses import StatorStageData
from Compressor import ConstantAxialVelocityConstantMeanLineCompressor as Compressor

'''
Class for a:
    Stator with a Constant Axial Velocity Constant Mean Line

'''
class ConstantAxialVelocityConstantMeanLineStator(CompressorBaseParameters):
    def __init__(self, number, stage_before:Compressor, stage_after:Compressor) -> None:
        super().__init__()

        self.number     = number
        print("  Stage", self.number, end="\r")

        """
        This basically just sets the parameters of this stator stage by
        getting the parameters of the compressor stage before and after. 
        It then checks to make sure it does not violate deHaller criteria
        using 

            C_3/C_2 >= 0.65
        
        """
        
        self.C_3m       = stage_after.mean["C_1"]
        self.alpha_3m   = np.deg2rad(stage_after.mean["alpha_1"])
        self.C_2m       = stage_before.mean["C_2"]
        self.alpha_2m   = np.deg2rad(stage_before.mean["alpha_2"])
        
        r_t1            = stage_before.r_t
        self.r_t3       = stage_after.r_t
        
        r_r1            = stage_before.r_r
        self.r_r3       = stage_after.r_r

        self.r_t2       = (r_t1 + self.r_t3)
        self.r_r2       = (r_r1 + self.r_r3)

        self.h_3        = stage_after.h
        self.h_2        = self.r_t2 - self.r_r2
                
        self.C_w2m      = self.C_2m*np.sin(self.alpha_2m)
        self.C_w3m      = self.C_3m*np.sin(self.alpha_3m)

        if abs(self.C_w2m - stage_before.mean["C_w2"]) > 1e-12:
            raise ValueError(f"Did not get matching C_w2ms: calculated = {self.C_w2m},  given = {stage_before.mean['C_w2']}")
        if abs(self.C_w3m - stage_after.mean["C_w1"]) > 1e-12:
            raise ValueError(f"Did not get matching C_w3ms: calculated = {self.C_w3m},  given = {stage_before.mean['C_w1']}")
        
        self.C_3r       = stage_after.root["C_1"]
        self.alpha_3r   = np.deg2rad(stage_after.root["alpha_1"])
        self.C_2r       = stage_before.root["C_2"]
        self.alpha_2r   = np.deg2rad(stage_before.root["alpha_2"])
        self.h_3        = stage_after.h
        self.h_2        = (stage_before.h + self.h_3)/2
                
        self.C_w2r      = self.C_2r*np.sin(self.alpha_2r)
        self.C_w3r      = self.C_3r*np.sin(self.alpha_3r)

        if abs(self.C_w2r - stage_before.root["C_w2"]) > 1e-12:
            raise ValueError(f"Did not get matching C_w2rs: calculated = {self.C_w2r},  given = {stage_before.root['C_w2']}")
        if abs(self.C_w3r - stage_after.root["C_w1"]) > 1e-12:
            raise ValueError(f"Did not get matching C_w3rs: calculated = {self.C_w3r},  given = {stage_before.root['C_w1']}")

        self.C_3t       = stage_after.tip["C_1"]
        self.alpha_3t   = np.deg2rad(stage_after.tip["alpha_1"])
        self.C_2t       = stage_before.tip["C_2"]
        self.alpha_2t   = np.deg2rad(stage_before.tip["alpha_2"])
        self.h_3        = stage_after.h
        self.h_2        = (stage_before.h + self.h_3)/2
                
        self.C_w2t      = self.C_2t*np.sin(self.alpha_2t)
        self.C_w3t      = self.C_3t*np.sin(self.alpha_3t)

        if abs(self.C_w2t - stage_before.tip["C_w2"]) > 1e-12:
            raise ValueError(f"Did not get matching C_w2ts: calculated = {self.C_w2t},  given = {stage_before.tip['C_w2']}")
        if abs(self.C_w3t - stage_after.tip["C_w1"]) > 1e-12:
            raise ValueError(f"Did not get matching C_w3ts: calculated = {self.C_w3t},  given = {stage_before.tip['C_w1']}")
        
        # Check de Haller Criteria
        deHaller_m = self.C_3m/self.C_2m
        if deHaller_m  < self.de_Haller_min:
            raise ValueError('de Haller check failed at {}: C3/C2 = {:.4f}'.format("mean", deHaller_m))

        # Check de Haller Criteria
        deHaller_r = self.C_3r/self.C_2r
        if deHaller_r  < self.de_Haller_min:
            raise ValueError('de Haller check failed at {}: C3/C2 = {:.4f}'.format("root", deHaller_r))
        
        # Check de Haller Criteria
        deHaller_t = self.C_3t/self.C_2t
        if deHaller_t  < self.de_Haller_min:
            raise ValueError('de Haller check failed at {}: C3/C2 = {:.4f}'.format("tip", deHaller_t))


        data = {
            "alpha_2" : np.rad2deg(self.alpha_2m),
            "alpha_3" : np.rad2deg(self.alpha_3m),
            "C_2" : self.C_2m,
            "C_3" : self.C_3m,
            "deHaller" : deHaller_m,
            # "r" : self.r_m
        }
        self.mean = StatorStageData(data, "mean")

        data = {
            "alpha_2" : np.rad2deg(self.alpha_2r),
            "alpha_3" : np.rad2deg(self.alpha_3r),
            "C_2" : self.C_2r,
            "C_3" : self.C_3r,
            "deHaller" : deHaller_r,
            # "r" : self.r_r
        }
        self.root = StatorStageData(data, "root")

        data = {
            "alpha_2" : np.rad2deg(self.alpha_2t),
            "alpha_3" : np.rad2deg(self.alpha_3t),
            "C_2" : self.C_2t,
            "C_3" : self.C_3t,
            "deHaller" : deHaller_t,
            # "r" : self.r_t
        }
        self.tip = StatorStageData(data, "tip")
    
    def getData(self):
        mean_data = self.mean.toDataFrame()
        tip_data  = self.tip.toDataFrame()
        root_data = self.root.toDataFrame()

        return pd.concat([root_data, mean_data, tip_data], keys=3*[f"Stator {self.number}"])



############
#
# The rest is just data output, not important
#
############



def outputData(stage_data):
    out = pd.concat(stage_data)
    for_formatting = StatorStageData({}, "")
    with open("Stator_Data.tex", 'w') as f:
        # angles = out[["alpha_2", "alpha_3", "beta_1", "beta_2"]]

        names = for_formatting.getFormattedColumns()
        _out = out.rename(columns=names)
        #angles = angles.rename(columns=names)
        #absolute_velocities = absolute_velocities.rename(columns=names)
        #relative_velocities = relative_velocities.rename(columns=names)
        #misc = misc.rename(columns=names)
        
        f.write(
            "\\begin{center}\n" + _out.to_latex(float_format='%.3f') + "\\end{center}\n\n"
        )

    out_T = out.T

    parameters = CompressorBaseParameters()

    plt.rcParams["text.usetex"] = True   

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
    plt.title(f"de-Haller Number along the Stator Stages")
    plt.legend()
    plt.ylabel("de-Haller Number")
    plt.xlabel("Stage")
    plt.grid(True)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(f"Stator_deHaller.png")