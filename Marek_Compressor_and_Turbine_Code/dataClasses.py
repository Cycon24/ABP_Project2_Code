from collections import OrderedDict
import pandas as pd


############
#
# This file is just data output, not important
#
############


'''
Just holds information about a compressor stage, not important to the structure of the code

'''
class CompressorStageData:
    parameters = [val.strip() for val in "beta_1, beta_2, alpha_1, alpha_2, C_w1, C_w2, C_1, C_2, V_w1, V_w2, V_1, V_2, U, M_1, M_2, Lambda, deHaller, r".split(",")]
    max_len = max([len(name) for name in parameters])
    
    def __init__(self, input_data, name) -> None:
        self.name = name
        self.data = OrderedDict()
        for param in self.parameters:
            if param in input_data:
                self.data[param] = input_data[param]
            else:
                self.data[param] = None
    
    def __getitem__(self, key):
        return self.data[key]
    
    def getFormattedColumns(self):
        names = {}
        specials = ["beta", "alpha", "Lambda", "deHaller"]
        for param in self.parameters:
            if "_" in param:
                temp = param[:param.find("_")+1] + "{" + param[param.find("_")+1:] + "}"
                names[param] = f"$\\{temp}$" if param[:param.find("_")] in specials else f"${temp}$"
            else:
                names[param] = f"$\\{param}$" if param in specials else f"${param}$"

        names["Ca"] = "$C_a$"
        return names

    def toDataFrame(self):
        return pd.DataFrame(self.data, index=[self.name])

    def print(self):
        print(2*" " + self.name)
        for param in self.data:
            if self.data[param] is not None:
                print(4*" " + param + (self.max_len-len(param))*" " + ": {:.4f}".format(self.data[param]))
            else:
                print(4*" " + param + (self.max_len-len(param))*" " +": None")


'''
Just holds information about a compressor stage, not important to the structure of the code

'''
class StatorStageData:
    parameters = [val.strip() for val in "alpha_2, alpha_3, C_2, C_3, deHaller".split(",")]
    max_len = max([len(name) for name in parameters])
    
    def __init__(self, input_data, name) -> None:
        self.name = name
        self.data = OrderedDict()
        for param in self.parameters:
            if param in input_data:
                self.data[param] = input_data[param]
            else:
                self.data[param] = None
    
    def __getitem__(self, key):
        return self.data[key]
    
    def getFormattedColumns(self):
        names = {}
        specials = ["beta", "alpha", "Lambda", "deHaller"]
        for param in self.parameters:
            if "_" in param:
                temp = param[:param.find("_")+1] + "{" + param[param.find("_")+1:] + "}"
                names[param] = f"$\\{temp}$" if param[:param.find("_")] in specials else f"${temp}$"
            else:
                names[param] = f"$\\{param}$" if param in specials else f"${param}$"

        # names["Ca"] = "$C_a$"
        return names

    def toDataFrame(self):
        return pd.DataFrame(self.data, index=[self.name])

    def print(self):
        print(2*" " + self.name)
        for param in self.data:
            if self.data[param] is not None:
                print(4*" " + param + (self.max_len-len(param))*" " + ": {:.4f}".format(self.data[param]))
            else:
                print(4*" " + param + (self.max_len-len(param))*" " +": None")



'''
Just holds information about a turbine stage, not important to the structure of the code

'''
class TurbineStageData:
    parameters = [val.strip() for val in "beta_2, beta_3, alpha_2, alpha_3, C_w2, C_w3, C_2, C_3, V_w2, V_w3, V_2, V_3, U, Lambda, phi, psi, r".split(",")]
    max_len = max([len(name) for name in parameters])
    
    def __init__(self, input_data, name) -> None:
        self.name = name
        self.data = OrderedDict()
        for param in self.parameters:
            if param in input_data:
                self.data[param] = input_data[param]
            else:
                self.data[param] = None
    
    def __getitem__(self, key):
        return self.data[key]
    
    def getFormattedColumns(self):
        names = {}
        specials = ["beta", "alpha", "Lambda", "deHaller", "phi", "psi"]
        for param in self.parameters:
            if "_" in param:
                temp = param[:param.find("_")+1] + "{" + param[param.find("_")+1:] + "}"
                names[param] = f"$\\{temp}$" if param[:param.find("_")] in specials else f"${temp}$"
            else:
                names[param] = f"$\\{param}$" if param in specials else f"${param}$"

        # names["Ca"] = "$C_a$"
        return names

    def toDataFrame(self):
        return pd.DataFrame(self.data, index=[self.name])

    def print(self):
        print(2*" " + self.name)
        for param in self.data:
            if self.data[param] is not None:
                print(4*" " + param + (self.max_len-len(param))*" " + ": {:.4f}".format(self.data[param]))
            else:
                print(4*" " + param + (self.max_len-len(param))*" " +": None")
