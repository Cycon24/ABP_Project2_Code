import numpy as np
import warnings
import matplotlib.pyplot as plt
import pandas as pd

# Function for interpolating WDF based on the stage number
def getWDF(stage_number):
    # Provided data
    data = """
    0.9529037215940281 0.9755790975957844
    1.220331540234934 0.9668349983532769
    1.577999780436933 0.9569875946865737
    1.936107146777914 0.947864749149193
    2.3543747941596234 0.9380063673290152
    2.8025030189922067 0.9274179383027774
    3.282907015040072 0.9200845317817543
    3.8234712921286658 0.9120155889779339
    4.243715007135801 0.9054177187397079
    4.724777692392141 0.8991711494126687
    5.26600065868921 0.8921890438028324
    5.898781424964323 0.886277308156768
    6.531562191239436 0.8803655725107036
    7.073663409814472 0.8748325831595125
    7.798001976067631 0.8699912174772203
    8.613678779229339 0.8658579426940388
    9.339554286968934 0.8635525304643759
    10.004171698320349 0.8601712591942035
    10.700186628609075 0.8585958941706004
    11.516961247118239 0.8562740147107255
    12.243495444066316 0.8550554396750466
    13.000109781534753 0.853469096497969
    13.84784279284225 0.8522285651553407
    14.6042375672412 0.8502799429136019
    15.270172357009557 0.8490723460313974
    15.936765835986392 0.8489515863431769
    16.69338017345483 0.8473652431660994
    17.541332747831824 0.8464869908881325
    18.42002415193765 0.846327807662751
    19.116478208365358 0.8454770007684707
    19.873751235042278 0.844977494785377
    20.57042485453948 0.8444889669557579
    """

    # Splitting the data into two arrays
    lines = data.strip().split('\n')
    stage_data, wdf_data = zip(*(map(float, line.split()) for line in lines))

    # Using numpy's interp function for linear interpolation
    wdf_interpolated = np.interp(stage_number, stage_data, wdf_data)
    return wdf_interpolated


def drawTotalGeometry(compressor_stage_data, turbine_stage_data):
    warnings.filterwarnings("ignore")
    compressor_data = pd.concat(compressor_stage_data)
    compressor_data = compressor_data.T

    plt.rcParams["text.usetex"] = True

    r_rs = []
    for item in compressor_data:
        r_rs.append(compressor_data[item[0]]["root"]["r"])

    r_ts = []
    for item in compressor_data:
        r_ts.append(compressor_data[item[0]]["tip"]["r"])
    
    r_ms = []
    for item in compressor_data:
        r_ms.append(compressor_data[item[0]]["mean"]["r"])
    
    plt.figure()
    x = np.linspace(0, 0.6, len(r_rs))
    plt.plot(x, r_rs, "k")
    plt.plot(x, r_ts, "k")
    plt.plot(x, -np.array(r_rs), "k")
    plt.plot(x, -np.array(r_ts), "k")
    plt.plot(x, r_ms, "--b", label="Compressor Mean Line")
    plt.plot(x, -np.array(r_ms), "--b")

    turbine_data = pd.concat(turbine_stage_data)
    turbine_data = turbine_data.T

    r_rs = []
    for item in turbine_data:
        r_rs.append(turbine_data[item[0]]["root"]["r"])

    r_ts = []
    for item in turbine_data:
        r_ts.append(turbine_data[item[0]]["tip"]["r"])
    
    r_ms = []
    for item in turbine_data:
        r_ms.append(turbine_data[item[0]]["mean"]["r"])

    x = np.linspace(0.9, 1, len(r_rs))
    plt.plot(x, r_rs, "k")
    plt.plot(x, r_ts, "k")
    plt.plot(x, -np.array(r_rs), "k")
    plt.plot(x, -np.array(r_ts), "k")
    plt.plot(x, r_ms, ".-b", label="Turbine Mean Line")
    plt.plot(x, -np.array(r_ms), ".-b")

    plt.title(f"Representative Geometry of the Turbine and Compressor")
    plt.legend()
    plt.ylabel("Radius (m)")
    plt.xlabel("Position")
    plt.grid(True)
    plt.minorticks_on()
    plt.gca().set_aspect("equal")
    plt.tight_layout()
    plt.savefig(f"Total_Geometry.png")
    warnings.filterwarnings("error")