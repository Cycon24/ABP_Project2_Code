import numpy as np

blade_type = ... # straight_backed/different, string

M3 = ... # Mach @ Rotor Exit
beta2r = ... # beta2 @ root, deg
beta2t = ... # beta2 @ tip, deg

"""
The blade incidence is independnet of Mach
"""

theta2r = beta2r - 5 # Degrees
theta2t = beta2t + 5 # Degrees

"""
The blade outlet is dependent of Mach
"""
o = ... # opening between consecutive blade
s = ... # pitch 
e = ... # radius of curfvature of aft suction side of blade

if blade_type == "straight_backed":
    theta3 = np.arccos(o/s) # Should be same as beta3, make sure to play with o & s to get theta3 = beta3 for straight_backed.
elif blade_type == "different":
    if M3 < 0.5:
        theta3 = 4 * (s/e) # set e such that theta3 = beta3 at mean. (Not so sure)
    elif M3 == 1:
        f_s_over_e = 0.0541 * (s/e) / (1 - 1.49 * (s/e) + 0.742 * (s/e)**2)
        theta3 = np.arccos(o/s) + f_s_over_e * np.arcsin(o/s)
    elif 0.5 <= M3 <= 1:
        theta3 = np.arccos(o/s)
        
# theta3 = np.rad2deg(theta3) # Activate this line if theta3 seems too small.