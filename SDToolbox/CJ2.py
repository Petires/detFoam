'''
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
'''

from SDToolbox import *
from numpy import *
import math

def CJspeed2(P1, T1, q, mech):

    """

    CJspeed2
    Calculates CJ detonation velocity and CJ state

    FUNCTION
    SYNTAX
    [cj_speed,gas] = CJspeed(P1,T1,q,mech)

    INPUT
    P1 = initial pressure (Pa)
    T1 = initial temperature (K)
    q = string of reactant species mole fractions
    mech = cti file containing mechanism data (i.e. 'gri30.cti')

    OUTPUT
    cj_speed = CJ detonation speed (m/s)
    gas = gas object at CJ state

    """
    
    gas2 = Solution(mech)
    gas1 = Solution(mech)
    gas  = Solution(mech)

    #INTIAL CONDITIONS
    gas.TPX  = T1, P1, q;
    gas1.TPX = T1, P1, q;
    gas2.TPX = T1, P1, q; 
    
    #INITIALIZE ERROR VALUES & CHANGE VALUES
    ERRFT = 1.0*10**-4;  ERRFV = 1.0*10**-4;

    r1 = gas1.density; V1 = 1/r1;
    P1 = gas1.P; T1 = gas1.T;
    i = 0;
    #PRELIMINARY GUESS
    Vg = V1/10; rg = 1/Vg; 
    
    gas.TD = T1,rg; 
    gas.equilibrate('UV')
    Tg = gas.T; 
    gas2.TDX = Tg, rg, gas.X
    
    #SAVE STATE
    V = Vg; r = rg;
    T = Tg;
    deltaT = 1000; deltaV = 1000; cj_speed = 0;
    #START LOOP
    while(abs(deltaT) > ERRFT*T or abs(deltaV) > ERRFV*V):
        i = i + 1
        if i == 500:
            print "CJ speed 2 calc did not converge"
            return gas
        
        #CALCULATE FH & FP FOR GUESS 1
        [FH,FP,cj_speed] = FHFP_CJ2(gas,gas1,gas2)


        #TEMPERATURE PERTURBATION
        DT = T*0.01; Tper = T + DT;
        Vper = V; Rper = 1/Vper;
        
        gas.TD = Tper, Rper
        gas.equilibrate('TV',2)
        gas2.TDX = Tper, Rper, gas.X

        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX,cj_speed] = FHFP_CJ2(gas,gas1,gas2)
        #ELEMENTS OF JACOBIAN
        DFHDT = (FHX-FH)/DT; DFPDT = (FPX-FP)/DT;

        #VOLUME PERTURBATION
        DV = 0.01*V; Vper = V + DV;
        Tper = T; Rper = 1/Vper;
        
        gas.TD = Tper, Rper
        gas.equilibrate('TV',2)
        gas2.TDX = Tper, Rper, gas.X
        
        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX,cj_speed] = FHFP_CJ2(gas,gas1,gas2)
        #ELEMENTS OF JACOBIAN
        DFHDV = (FHX-FH)/DV; DFPDV = (FPX-FP)/DV;

        #INVERT MATRIX
        J = DFHDT*DFPDV - DFPDT*DFHDV
        b = [DFPDV, -DFHDV, -DFPDT, DFHDT]
        a = [-FH, -FP]
        deltaT = (b[0]*a[0]+b[1]*a[1])/J; deltaV = (b[2]*a[0]+b[3]*a[1])/J;

        #CHECK & LIMIT CHANGE VALUES
        #TEMPERATURE
        DTM = 0.2*T
        if abs(deltaT) > DTM:
            deltaT = DTM*deltaT/abs(deltaT)
        #VOLUME
        V2X = V + deltaV
        if V2X > V1:
            DVM = 0.5*(V1 - V)
        else:
            DVM = 0.2*V
        if abs(deltaV) > DVM:
            deltaV = DVM*deltaV/abs(deltaV)
        #MAKE THE CHANGES
        T = T + deltaT; V = V + deltaV; r = 1/V;
        gas.TD = T, r
        gas.equilibrate('TV',2)
        gas2.TDX = T, r, gas.X

    [FH,FP,cj_speed] = FHFP_CJ2(gas,gas1,gas2)
    
    return [gas,cj_speed]


def FHFP_CJ2(gas,gas1,gas2):

    """

    FHFP_CJ2
    Uses the momentum and energy conservation equations and the equilibrium sound speed to calculate error in current pressure and enthalpy guesses.  In this case, state 2 is in equilibrium.
    
    FUNCTION
    SYNTAX
    [FH,FP,cj_speed] = FHFP_CJ2(gas,gas1,gas2)

    INPUT
    gas = working gas object
    gas1 = gas object at initial state 
    gas2 = dummy gas object (for calculating numerical derivatives)
    
    OUTPUT
    FH,FP = error in enthalpy and pressure
    cj_speed = CJ detonation speed (m/s)
    
    """
    
    P1 = gas1.P
    H1 = gas1.enthalpy_mass
    r1 = gas1.density
    P2 = gas.P
    H2 = gas.enthalpy_mass
    r2 = gas.density
     
    speeds = equilSoundSpeeds(gas2)
    w2s=(speeds[0])**2
    w1s = w2s*(r2/r1)**2
    FH = H2 + 0.5*w2s - (H1 + 0.5*w1s)
    FP = P2 + r2*w2s - (P1 + r1*w1s)
    return [FH, FP, sqrt(w1s)]

def equilSoundSpeeds(gas):

    """

    equilSoundSpeeds
    Calculates equilibrium and frozen sound speeds. For the equilibrium sound speed, the gas is equilibrated holding entropy and specific volume constant.

    FUNCTION
    SYNTAX
    [aequil,afrozen] = equilSoundSpeeds(gas)

    INPUT
    gas = working gas object (modified inside function)

    OUTPUT
    aequil = equilibrium sound speed (m/s)
    afrozen = frozen sound speed (m/s)

    """

    # set the gas to equilibrium at its current T and P
    gas.equilibrate('TP')

    # save properties
    s0 = gas.entropy_mass
    p0 = gas.P
    r0 = gas.density

    # perturb the density
    r1 = r0*1.0001
 
    # set the gas to a state with the same entropy and composition but
    # the perturbed density
    gas.SV = s0, 1.0/r1

    # save the pressure for this case for the frozen sound speed
    pfrozen = gas.P
        
    # now equilibrate the gas holding S and V constant
    gas.equilibrate("SV")
    
    p1 = gas.P

    # equilibrium sound speed
    aequil = math.sqrt((p1 - p0)/(r1 - r0));

    # frozen sound speed
    afrozen = math.sqrt((pfrozen - p0)/(r1 - r0));    
    return (aequil, afrozen)

# test program
if __name__ == "__main__":

    gas = GRI30()
    gas.set(X = 'CH4:0.1.1, O2:2.0, N2:3.76')
    for n in range(27):
        temp = 300.0 + n*100.0
        gas.set(T = temp, P = OneAtm)
        print temp, equilSoundSpeeds(gas)



