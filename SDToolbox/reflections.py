'''
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
'''

from numpy     import *
from SDToolbox import *
import math


def reflected_fr(gas1,gas2,gas3,UI):
    """

    reflected_fr
    Calculates frozen post-reflected-shock state assumming u1 = 0

    FUNCTION
    SYNTAX
    [p3,UR,gas3] = reflected_fr(gas1,gas2,gas3,UI)

    INPUT:
    gas1 = gas object at initial state
    gas2 = gas object at post-incident-shock state (already computed)
    gas3 = working gas object
    UI = incident shock speed (m/s)

    OUTPUT:
    p3 = post-reflected-shock pressure (Pa)
    UR = reflected shock speed (m/s)
    gas3 = gas object at frozen post-reflected-shock state

    """
        
    p2 = gas2.P;
    p1 = gas1.P;
    rho2 = gas2.density;
    v2=1/rho2;
    rho1 = gas1.density;
    v1=1/rho1;
    gamma = gas1.cp_mole/gas1.cv_mole;
    T1 = gas1.T;
    T2 = gas2.T;

    u2 = sqrt((p2-p1)*(v1-v2)); #particle velocity
    w2 = UI - u2; # velocity in shock fixed frame

    #BASIC PRELIMINARY GUESS
    v3 = .2/rho2
    r3 = 1/v3;
    p3 = p2 + rho2*(UI**2)*(1-v3/v2);
    T3 = T2*p3*v3/(p2*v2);

    gas3.TPX = T3, p3, gas2.X; 
    gas3 = PostReflectedShock_fr(u2, gas2,gas3);
    p3 = gas3.P;
    UR = (p3-p2)/u2/rho2-u2;
    
    return [p3,UR,gas3]


def reflected_eq(gas1,gas2,gas3,UI):
    """

    reflected_eq
    Calculates equilibrium post-reflected-shock state assumming u1 = 0

    FUNCTION
    SYNTAX:
    [p3,UR,gas3] = reflected_eq(gas1,gas2,gas3,UI)

    INPUT:
    gas1 = gas object at initial state
    gas2 = gas object at post-incident-shock state (already computed)
    gas3 = working gas object
    UI = incident shock speed (m/s)

    OUTPUT:
    p3 = post-reflected-shock pressure (Pa)
    UR = reflected shock speed (m/s)
    gas3 = gas object at equilibrium post-reflected-shock state

    """
    
    p2 = gas2.P
    p1 = gas1.P
    rho2 = gas2.density
    v2=1/rho2;
    rho1 = gas1.density
    v1=1/rho1;
    gamma2 = gas2.cp_mole/gas2.cv_mole
    T1 = gas1.T;
    T2 = gas2.T


    u2 = sqrt((p2-p1)*(v1-v2)); #particle velocity
    w2 = UI - u2; # velocity in shock fixed frame

    #BASIC PRELIMINARY GUESS
    v3 = 0.2/rho2;
    r3 = 1/v3;
    p3 = p2 + rho2*(UI**2)*(1-v3/v2);
    T3 = T2*p3*v3/(p2*v2);

    gas3.TPX = T3, p3, gas2.X
    gas3 = PostReflectedShock_eq(u2, gas2,gas3);
    p3 = gas3.P;
    UR = (p3-p2)/u2/rho2-u2;    
    
    return [p3,UR,gas3]


def PostReflectedShock_fr(u2, gas2,gas3):
    """

    PostReflectedShock_fr
    Calculates frozen post-rReflected-shock state for a specified shock velocity

    FUNCTION
    SYNTAX
    [gas3] = PostReflectedShock_fr(u2,gas2,gas3)

    INPUT
    u2 = current post-incident-shock lab frame particle speed
    gas2 = gas object at post-incident-shock state (already computed)
    gas3 = working gas object

    OUTPUT
    gas3 = gas object at frozen post-reflected-shock state

    """

    #INITIALIZE ERROR VALUES
    ERRFT = 1.0*10**-4;
    ERRFV = 1.0*10**-4;

    #CALCULATE POST-REFLECTED SHOCK STATE
    P2 = gas2.P
    H2 = gas2.enthalpy_mass
    T2 = gas2.T
    r2 = gas2.density
    V2 = 1/r2
    
    print 'State 2'
    print 'Post of Incident Shock'
    print "T2 = %s" % (T2)
    print "V2 = %s" % (V2)
    print "P2 = %s" % (P2)
    print "H2 = %s" % (H2)

    j = 0;
    deltaT = 1000;
    deltaV = 1000;

    ##################################################################################################
    #PRELIMINARY GUESS
    P = gas3.P;
    H = gas3.enthalpy_mass;
    T = gas3.T;
    r = gas3.density;
    V = 1/r;
    ##################################################################################################
    #START LOOP

    flag = 0;

    while ((abs(deltaT) > ERRFT*T) or (abs(deltaV) > ERRFV*V)):
        j = j + 1;
        
        if j == 500:
            print "Calculation did not converge for U = %s" % (u2)
            flag = 1;
            return
                
        
        ##############################################################################################
        #CALCULATE FH & FP FOR GUESS 1
        [FH,FP] = FHFP_reflected_fr(u2,gas3,gas2);
        ##############################################################################################
        #TEMPERATURE PERTURBATION
        DT = T*0.02;
        Tper = T + DT;
        Vper = V;
        Rper = 1/Vper;        
        [Pper, Hper] = state(gas3,Rper,Tper);
        
        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX] = FHFP_reflected_fr(u2,gas3,gas2);

        #ELEMENTS OF JACOBIAN
        DFHDT = (FHX-FH)/DT;
        DFPDT = (FPX-FP)/DT;
        ##############################################################################################
        #VOLUME PERTURBATION
        DV = 0.02*V;
        Vper = V + DV;
        Tper = T;
        Rper = 1/Vper;
        
        [Pper, Hper] = state(gas3,Rper,Tper);
        
        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX] = FHFP_reflected_fr(u2,gas3,gas2);
       
        #ELEMENTS OF JACOBIAN
        DFHDV = (FHX-FH)/DV;
        DFPDV = (FPX-FP)/DV;
        ##############################################################################################
        #INVERT MATRIX
        J = DFHDT*DFPDV - DFPDT*DFHDV
        b = [DFPDV, -DFHDV, -DFPDT, DFHDT]
        a = [-FH, -FP]
        deltaT = (b[0]*a[0]+b[1]*a[1])/J; deltaV = (b[2]*a[0]+b[3]*a[1])/J;

        
        ############################
        #CHECK & LIMIT CHANGE VALUES
        ############################
        #TEMPERATURE
        #VOLUME
        DTM = 0.2*T;
        if (abs(deltaT) > DTM):
            deltaT = DTM*deltaT/abs(deltaT);
         
        ############################
        #VOLUME
        V3X = V + deltaV;
        if V3X > V2:
            DVM = 0.5*(V2 - V);
        else:
            DVM = 0.2*V;
         
        if abs(deltaV) > DVM:
            deltaV = DVM*deltaV/abs(deltaV);
         
        ##############################################################################################
        #MAKE THE CHANGES
        T = T + deltaT;
        V = V + deltaV;
        r = 1/V;
        [P, H] = state(gas3,r,T);

     

    T3 = T
    V3 = V
    P3 = P
    H3 = H
    r3 = r
    print 'State 3'
    print 'Post of Reflected Shock'
    print "T3 = %s" % (T3)
    print "V3 = %s" % (V3)
    print "P3 = %s" % (P3)
    print "H3 = %s" % (H3)

    return gas3


def PostReflectedShock_eq(u2, gas2,gas3):
    """

    PostReflectedShock_eq
    Calculates equilibrium post-reflected-shock state for a specified shock velocity

    FUNCTION
    SYNTAX
    [gas3] = PostReflectedShock_fr(u2,gas2,gas3)

    INPUT
    u2 = current post-incident-shock lab frame particle speed
    gas2 = gas object at post-incident-shock state (already computed)
    gas3 = working gas object

    OUTPUT
    gas3 = gas object at equilibrium post-reflected-shock state

    """

    #INITIALIZE ERROR VALUES
    ERRFT = 1.0*10**-4;
    ERRFV = 1.0*10**-4;

    #CALCULATE POST-REFLECTED SHOCK STATE
    P2 = gas2.P
    H2 = gas2.enthalpy_mass
    T2 = gas2.T
    r2 = gas2.density
    V2 = 1/r2
    
    print 'State 2'
    print 'Post of Incident Shock'
    print "T2 = %s" % (T2)
    print "V2 = %s" % (V2)
    print "P2 = %s" % (P2)
    print "H2 = %s" % (H2)
    
    j = 0;
    deltaT = 1000;
    deltaV = 1000;

    ##################################################################################################
    #PRELIMINARY GUESS
    P = gas3.P;
    H = gas3.enthalpy_mass;
    T = gas3.T;
    r = gas3.density;
    V = 1/r;
    [P, H] = eq_state(gas3,r,T);
    ##################################################################################################
    #START LOOP

    flag = 0;

    while (abs(deltaT) > ERRFT*T or abs(deltaV) > ERRFV*V):
        j = j + 1;
        if j == 500:
            print "Calculation did not converge for U = %s" % (u2)
            flag = 1;
            return gas
                
        ##############################################################################################
        #CALCULATE FH & FP FOR GUESS 1
        [FH,FP] = FHFP_reflected_fr(u2,gas3,gas2);
        ##############################################################################################
        #TEMPERATURE PERTURBATION
        DT = T*0.02;
        Tper = T + DT;
        Vper = V;
        Rper = 1/Vper;
        [Pper, Hper] = eq_state(gas3,Rper,Tper);
        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX] = FHFP_reflected_fr(u2,gas3,gas2);

        #ELEMENTS OF JACOBIAN
        DFHDT = (FHX-FH)/DT;
        DFPDT = (FPX-FP)/DT;
        ##############################################################################################

        #VOLUME PERTURBATION
        DV = 0.02*V
        Vper = V + DV
        Tper = T
        Rper = 1/Vper
        [Pper, Hper] = eq_state(gas3,Rper,Tper)
        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX] = FHFP_reflected_fr(u2,gas3,gas2)
        #ELEMENTS OF JACOBIAN
        DFHDV = (FHX-FH)/DV
        DFPDV = (FPX-FP)/DV
        ##############################################################################################

        #INVERT MATRIX
        J = DFHDT*DFPDV - DFPDT*DFHDV
        b = [DFPDV, -DFHDV, -DFPDT, DFHDT]
        a = [-FH, -FP]
        deltaT = (b[0]*a[0]+b[1]*a[1])/J; deltaV = (b[2]*a[0]+b[3]*a[1])/J;

        
        ############################
        #CHECK & LIMIT CHANGE VALUES
        ############################
        #TEMPERATURE
        #VOLUME
        DTM = 0.2*T;
        if (abs(deltaT) > DTM):
            deltaT = DTM*deltaT/abs(deltaT);
         
        ############################
        #VOLUME
        V3X = V + deltaV;
        if V3X > V2:
            DVM = 0.5*(V2 - V);
        else:
            DVM = 0.2*V;
         
        if abs(deltaV) > DVM:
            deltaV = DVM*deltaV/abs(deltaV);
         
        ##############################################################################################
        #MAKE THE CHANGES
        T = T + deltaT;
        V = V + deltaV;
        r = 1/V;
        [P, H] = eq_state(gas3,r,T);

    T3 = T
    V3 = V
    P3 = P
    H3 = H
    r3 = r
    print 'State 3'
    print 'Post of Reflected Shock'
    print "T3 = %s" % (T3)
    print "V3 = %s" % (V3)
    print "P3 = %s" % (P3)
    print "H3 = %s" % (H3)
    
    return gas3


def FHFP_reflected_fr(u2,gas3,gas2):
    """

    FHFP_reflected_fr
    Uses the momentum and energy conservation equations to calculate error in current pressure and enthalpy guesses. In this case, state 3 is frozen.

    FUNCTION
    SYNTAX
    FHFP_reflected_fr(u2,gas3,gas2)

    INPUT
    u2 = current post-incident-shock lab frame particle speed
    gas3 = working gas object
    gas2 = gas object at post-incident-shock state (already computed)

    OUTPUT
    FH,FP = error in enthalpy and pressure

    """
    P2 = gas2.P;
    H2 = gas2.enthalpy_mass;
    r2 = gas2.density;
    
    P3 = gas3.P;
    H3 = gas3.enthalpy_mass;
    r3 = gas3.density;
    
    FH = H3 - H2 - 0.5*(u2**2)*((r3/r2)+1)/(r3/r2-1);
    FP = P3 - P2 - r3*(u2**2)/(r3/r2-1);

    return [FH,FP]

