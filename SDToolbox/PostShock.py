'''
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
'''

from numpy     import *
from cantera   import *
from SDToolbox import *

def FHFP_CJ(gas,gas1,w1):
    """

    FHFP_CJ
    Uses the momentum and energy conservation equations to calculate error in current pressure and enthalpy guesses.  In this case, state 2 is in equilibrium.

    FUNCTION
    SYNTAX
    [FH,FP] = FHFP_CJ(gas,gas1,w1)
    
    INPUT
    gas = working gas object
    gas1 = gas object at initial state
    w1 = current guess for initial velocity (m/s)

    OUTPUT
    FH,FP = error in enthalpy and pressure

    """

    P1 = gas1.P
    H1 = gas1.enthalpy_mass
    r1 = gas1.density
    P2 = gas.P
    H2 = gas.enthalpy_mass
    r2 = gas.density
    
    w2 = w1*(r1/r2)
    w1s = w1**2
    w2s = w2**2
    FH = H2 + 0.5*w2s - (H1 + 0.5*w1s)
    FP = P2 + r2*w2s - (P1 + r1*w1s)
    return [FH, FP]

def FHFP_fr(U1,gas,gas1):
    """

    FHFP_fr
    Uses the momentum and energy conservation equations to calculate error in current pressure and enthalpy guesses.  In this case, state 2 is frozen.

    FUNCTION
    SYNTAX
    [FH,FP] = FHFP_fr(U1,gas,gas1)
    
    INPUT
    U1 = shock speed (m/s)
    gas = working gas object
    gas1 = gas object at initial state

    OUTPUT
    FH,FP = error in enthalpy and pressure

    """


    P1 = gas1.P
    H1 = gas1.enthalpy_mass
    r1 = gas1.density
    P2 = gas.P
    H2 = gas.enthalpy_mass
    r2 = gas.density
    
    w1s = U1**2;
    w2s = w1s*(r1/r2)**2;
    FH = H2 + 0.5*w2s - (H1 + 0.5*w1s);
    FP = P2 + r2*w2s - (P1 + r1*w1s);
    return [FH, FP]

def CJ_calc(gas, gas1, ERRFT, ERRFV, x):
    """

    CJ_calc
    Calculates the Chapman-Jouguet wave speed using Reynolds' iterative method.

    FUNCTION
    SYNTAX
    [gas,w1] = CJ_calc(gas,gas1,ERRFT,ERRFV,x)

    INPUT
    gas = working gas object
    gas1 = gas object at initial state
    ERRFT,ERRFV = error tolerances for iteration
    x = density ratio

    OUTPUT
    gas = gas object at equilibrium state
    w1 = initial velocity to yield prescribed density ratio

    """
    T = 2000; r1 = gas1.density
    V1 = 1/r1; P1 = gas1.P
    i = 0; DT = 1000; DV = 1000; DP = 1000;
    #PRELIMINARY GUESS
    V = V1/x; r = 1/V; w1 = 2000;
    [P, H] = eq_state(gas,r,T)
    #START LOOP
    while (abs(DT) > ERRFT*T or abs(DW) > ERRFV*w1):
        i = i + 1
	if i == 500:
            'i = 500'
            return
    	#CALCULATE FH & FP FOR GUESS 1
	[FH,FP] = FHFP_CJ(gas,gas1,w1)

	#TEMPERATURE PERTURBATION
	DT = T*0.02; Tper = T + DT;
	Vper = V; Rper = 1/Vper;
	Wper = w1;
	[Pper, Hper] = eq_state(gas,Rper,Tper)
	#CALCULATE FHX & FPX FOR "IO" STATE
	[FHX,FPX] = FHFP_CJ(gas,gas1,Wper)
	#ELEMENTS OF JACOBIAN
	DFHDT = (FHX-FH)/DT; DFPDT = (FPX-FP)/DT;

	#VELOCITY PERTURBATION
	DW = 0.02*w1; Wper = w1 + DW;
	Tper = T; Rper = 1/V;
	[Pper, Hper] = eq_state(gas,Rper,Tper)
	#CALCULATE FHX & FPX FOR "IO" STATE
	[FHX,FPX] = FHFP_CJ(gas,gas1,Wper)
	#ELEMENTS OF JACOBIAN
	DFHDW = (FHX-FH)/DW; DFPDW = (FPX-FP)/DW;

	#INVERT MATRIX
	J = DFHDT*DFPDW - DFPDT*DFHDW
	b = [DFPDW, -DFHDW, -DFPDT, DFHDT]
	a = [-FH, -FP]
	DT = (b[0]*a[0]+b[1]*a[1])/J; DW = (b[2]*a[0]+b[3]*a[1])/J;
	
	#CHECK & LIMIT CHANGE VALUES
	#VOLUME
	DTM = 0.2*T
	if abs(DT) > DTM:
            DT = DTM*DT/abs(DT)
        #MAKE THE CHANGES
	T = T + DT; w1 = w1 + DW;
	[P, H] = eq_state(gas,r,T)
    return [gas, w1]

def CJspeed(P1, T1, q, mech, plt_num):
    """

    CJspeed
    Calculates CJ detonation velocity

    FUNCTION
    SYNTAX
    [cj_speed,R2] = CJspeed(P1,T1,q,mech,plt_num)

    INPUT
    P1 = initial pressure (Pa)
    T1 = initial temperature (K)
    q = string of reactant species mole fractions
    mech = cti file containing mechanism data (i.e. 'gri30.cti')
    plt_num = unused

    OUTPUT
    cj_speed = CJ detonation speed (m/s)
    R2 = R-squared value of LSQ curve fit

    """
    #DECLARATIONS
    numsteps = 20; maxv = 2.0; minv = 1.5; 


    w1 = zeros(numsteps+1,float)
    rr = zeros(numsteps+1,float)

    gas1 = Solution(mech)
    gas  = Solution(mech)
    #INTIAL CONDITIONS
    gas.TPX  = T1, P1, q
    gas1.TPX = T1, P1, q
    
    #INITIALIZE ERROR VALUES & CHANGE VALUES
    ERRFT = 1.0*10**-4;  ERRFV = 1.0*10**-4;

    i = 1;
    T1 = gas1.T; P1 = gas1.P;
    
    counter = 1; R2 = 0.0; cj_speed = 0.0
    a = 0.0; b = 0.0; c = 0.0; dnew = 0.0
    while (counter <= 4) or (R2 < 0.99999):
	step = (maxv-minv)/float(numsteps)
	i = 0
	x = minv
	while x <= maxv:
            gas.TPX = T1, P1, q
            [gas, temp] = CJ_calc(gas, gas1, ERRFT, ERRFV, x)
            w1[i] = temp
            rr[i] = gas.density/gas1.density
            i = i + 1; x = x + step;
        [a,b,c,R2,SSE,SST] = LSQ_CJspeed(rr,w1)

        dnew = -b/(2.0*a)
        minv = dnew - dnew*0.001
        maxv = dnew + dnew*0.001
        counter = counter + 1;

    cj_speed = a*dnew**2 + b*dnew + c
    return [cj_speed,R2]

def PostShock_fr(U1, P1, T1, q, mech):
    """
     
    PostShock_fr
    Calculates frozen post-shock state for a specified shock velocity

    FUNCTION
    SYNTAX
    [gas] = PostShock_fr(U1,P1,T1,q,mech)
    
    INPUT
    U1 = shock speed (m/s)
    P1 = initial pressure (Pa)
    T1 = initial temperature (K)
    q = string of reactant species mole fractions
    mech = cti file containing mechanism data (i.e. 'gri30.cti')
    
    OUTPUT
    gas = gas object at frozen post-shock state
    
    """
    gas1 = Solution(mech)
    gas  = Solution(mech)
    #INTIAL CONDITIONS
    gas.TPX = T1, P1, q;
    gas1.TPX= T1, P1, q;
    #INITIALIZE ERROR VALUES
    ERRFT = 1.0*10**-4; ERRFV = 1.0*10**-4;
    #CALCULATES POST-SHOCK STATE
    gas = shk_calc(U1, gas, gas1, ERRFT, ERRFV)
    return gas

def PostShock_eq(U1, P1, T1, q, mech):
    """

    PostShock_eq
    Calculates equilibrium post-shock state for a specified shock velocity

    FUNCTION
    SYNTAX
    [gas] = PostShock_eq(U1,P1,T1,q,mech)

    INPUT
    U1 = shock speed (m/s)
    P1 = initial pressure (Pa)
    T1 = initial temperature (K)
    q = string of reactant species mole fractions
    mech = cti file containing mechanism data (i.e. 'gri30.cti')

    OUTPUT
    gas = gas object at equilibrium post-shock state

    """
    gas1 = Solution(mech);
    gas  = Solution(mech);
    #INTIAL CONDITIONS
    gas.TPX = T1, P1, q
    gas1.TPX= T1, P1, q
    #INITIALIZE ERROR VALUES
    ERRFT = 1.0*10**-4; ERRFV = 1.0*10**-4;
    #CALCULATES POST-SHOCK STATE
    gas = shk_eq_calc(U1, gas, gas1, ERRFT, ERRFV)
    return gas

def shk_calc(U1, gas, gas1, ERRFT, ERRFV):
    """

    shk_calc
    Calculates frozen post-shock state using Reynolds' iterative method

    FUNCTION
    SYNTAX
    [gas] = shk_calc(U1,gas,gas1,ERRFT,ERRFV)
    
    INPUT
    U1 = shock speed (m/s)
    gas = working gas object
    gas1 = gas object at initial state
    ERRFT,ERRFV = error tolerances for iteration

    OUTPUT
    gas = gas object at frozen post-shock state

    """
    r1 = gas1.density; V1 = 1/r1;
    P1 = gas1.P; T1 = gas1.T;
    i = 0;
    deltaT = 1000; deltaV = 1000;
    #PRELIMINARY GUESS
    Vg = V1/5; rg = 1/Vg;
    Pg = P1 + r1*(U1**2)*(1-Vg/V1); Tg = T1*Pg*Vg/(P1*V1);
    [Pg, Hg] = state(gas,rg,Tg)
    #SAVE STATE
    V = Vg; r = rg; P = Pg;
    T = Tg; H = Hg;
    #START LOOP
    while(abs(deltaT) > ERRFT*T or abs(deltaV) > ERRFV*V):
        i = i + 1
        if i == 500:
            print "shk_calc did not converge for U = %s" % (U1)
            return gas
        #CALCULATE FH & FP FOR GUESS 1
        [FH,FP] = FHFP_fr(U1,gas,gas1)

        #TEMPERATURE PERTURBATION
        DT = T*0.02; Tper = T + DT;
        Vper = V; Rper = 1/Vper;
        [Pper, Hper] = state(gas,Rper,Tper)
        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX] = FHFP_fr(U1,gas,gas1)
        #ELEMENTS OF JACOBIAN
        DFHDT = (FHX-FH)/DT; DFPDT = (FPX-FP)/DT;

        #VOLUME PERTURBATION
        DV = 0.02*V; Vper = V + DV;
        Tper = T; Rper = 1/Vper;
        [Pper, Hper] = state(gas,Rper,Tper)
        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX] = FHFP_fr(U1,gas,gas1)
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
        [P, H] = state(gas,r,T)

    T2 = T; V2 = V; P2 = P;
    H2 = H; r2 = r;
    return gas

def shk_eq_calc(U1, gas, gas1, ERRFT, ERRFV):
    """

    shk_eq_calc
    Calculates equilibrium post-shock state using Reynolds' iterative method

    FUNCTION
    SYNTAX
    [gas] = shk_calc(U1,gas,gas1,ERRFT,ERRFV)
    
    INPUT
    U1 = shock speed (m/s)
    gas = working gas object
    gas1 = gas object at initial state
    ERRFT,ERRFV = error tolerances for iteration

    OUTPUT
    gas = gas object at equilibrium post-shock state

    """
    r1 = gas1.density; V1 = 1/r1;
    P1 = gas1.P; T1 = gas1.T;
    i = 0;
    deltaT = 1000; deltaV = 1000;
    #PRELIMINARY GUESS
    V = V1/5; r = 1/V;
    P = P1 + r1*(U1**2)*(1-V/V1); T = T1*P*V/(P1*V1);
    [P, H] = eq_state(gas,r,T)
    #START LOOP
    while(abs(deltaT) > ERRFT*T or abs(deltaV) > ERRFV*V):
        i = i + 1
        if i == 500:
            print "shk_calc did not converge for U = %s" % (U1)
            return gas
        #CALCULATE FH & FP FOR GUESS 1
        [FH,FP] = FHFP_fr(U1,gas,gas1)

        #TEMPERATURE PERTURBATION
        DT = T*0.02; Tper = T + DT;
        Vper = V; Rper = 1/Vper;
        [Pper, Hper] = eq_state(gas,Rper,Tper)
        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX] = FHFP_fr(U1,gas,gas1)
        #ELEMENTS OF JACOBIAN
        DFHDT = (FHX-FH)/DT; DFPDT = (FPX-FP)/DT;

        #VOLUME PERTURBATION
        DV = 0.02*V; Vper = V + DV;
        Tper = T; Rper = 1/Vper;
        [Pper, Hper] = eq_state(gas,Rper,Tper)
        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX] = FHFP_fr(U1,gas,gas1)
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
        [P, H] = eq_state(gas,r,T)

    T2 = T; V2 = V; P2 = P;
    H2 = H; r2 = r;
    return gas

