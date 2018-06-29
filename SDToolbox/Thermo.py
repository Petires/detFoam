'''
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
'''

from numpy   import *
from cantera import *

def eq_state(gas, r1, T1):
    """

    eq_state
    Calculates equilibrium state given T & rho

    FUNCTION
    SYNTAX
    [P,H] = eq_state(gas,r1,T1)

    INPUT
    gas = working gas object
    r1,T1 = desired density and temperature

    OUTPUT
    P,H = equilibrium pressure and enthlapy at given temperature and density

    """
    gas.TD = T1, r1
    gas.equilibrate('TV')
    P = gas.P
    H = gas.enthalpy_mass
    return [P, H]

def state(gas,r1,T1):
    """

    state
    Calculates frozen state given T & rho

    FUNCTION
    SYNTAX
    [P,H] = state(gas,r1,T1)

    INPUT
    gas = working gas object
    r1,T1 = desired density and temperature

    OUTPUT
    P,H = pressure and enthlapy

    """
    gas.TD = T1, r1
    P = gas.P;
    H = gas.enthalpy_mass
    return [P, H]

def LSQ_CJspeed(x,y):
    """

    LSQ_CJspeed
    Determines least squares fit of parabolic data.

    FUNCTION
    SYNTAX
    [a,b,c,R2,SSE,SST] = LSQ_CJspeed(x,y)

    INPUT
    x = independent data points
    y = dependent data points

    OUTPUT
    a,b,c = coefficients of quadratic function (ax^2 + bx + c = 0)
    R2 = R-squared value
    SSE = sum of squares due to error
    SST = total sum of squares

    """
    #Calculate Sums
    k = 0
    X = 0.0; X2 = 0.0; X3 = 0.0; X4 = 0.0;
    Y = 0.0; Y1 = 0.0; Y2 = 0.0;
    a = 0.0; b = 0.0; c = 0.0; R2 = 0.0
    n = size(x)

    while k < n:
        X = X + x[k]
        X2 = X2 + x[k]**2
        X3 = X3 + x[k]**3
        X4 = X4 + x[k]**4
        Y = Y + y[k]
        Y1 = Y1 + y[k]*x[k]
        Y2 = Y2 + y[k]*x[k]**2
        k= k + 1
    m = float(Y)/float(n)

    den = (X3*float(n) - X2*X)
    temp = (den*(X*X2-X3*float(n))+X2*X2*(X*X-float(n)*X2)-X4*float(n)*(X*X-X2*float(n)))
    temp2 = (den*(Y*X2-Y2*float(n)) + (Y1*float(n)-Y*X)*(X4*float(n)-X2*X2))

    b = temp2/temp
    a = 1.0/den*(float(n)*Y1 - Y*X - b*(X2*float(n)-X*X))
    c = 1/float(n)*(Y - a*X2 - b*X)

    k= 0; SSE = 0.0; SST = 0.0;

    f = zeros(size(x),float)
    
    while k < size(x):
        f[k] = a*x[k]**2 + b*x[k] + c
        SSE = SSE + (y[k] - f[k])**2
        SST = SST + (y[k] - m)**2
        k = k + 1
    R2 = 1 - SSE/SST

    return [a,b,c,R2,SSE,SST]

def hug_fr(x,vb,h1,P1,v1,gas):
    """

    hug_fr
    Algebraic expressions of frozen (reactant) Hugoniot pressure and enthalpy. Passed to root solver 'fsolve'.

    FUNCTION
    SYNTAX
    fval = fsolve(hug_fr,Ta,args=(vb,h1,P1,v1,gas))

    INPUT
    Ta = initial guess for frozen Hugoniot temperature (K)
    vb = desired frozen Hugoniot specific volume (m^3/kg)
    h1 = enthalpy at state 1 (J/kg)
    P1 = pressure at state 1 (Pa)
    v1 = specific volume at state 1 (m^3/kg)
    gas = working gas object

    OUTPUT
    fval = frozen Hugoniot temperature corresponding to vb (K)

    """
    gas.TD = x, 1.0/vb
    hb1 = gas.enthalpy_mass
    Pb = gas_constant*x/(gas.mean_molecular_weight*vb)
    
    hb2 = h1 + 0.5*(Pb-P1)*(vb+v1)
    return hb2-hb1
    
def hug_eq(x,vb,h1,P1,v1,gas):
    """ 

    hug_eq
    Algebraic expressions of equilibrium (product) Hugoniot pressure and enthalpy. Passed to root solver 'fsolve'.

    FUNCTION
    SYNTAX
    fval = fsolve(hug_eq,Ta,args=(vb,h1,P1,v1,gas))

    INPUT
    Ta = initial guess for equilibrium Hugoniot temperature (K)
    vb = desired equilibrium Hugoniot specific volume (m^3/kg)
    h1 = enthalpy at state 1 (J/kg)
    P1 = pressure at state 1 (Pa)
    v1 = specific volume at state 1 (m^3/kg)
    gas = working gas object

    OUTPUT
    fval = equilibrium Hugoniot temperature corresponding to vb (K)

    """
    gas.TD = x, 1.0/vb
    gas.equilibrate('TV')
    hb1 = gas.enthalpy_mass
    Pb = gas_constant*x/(gas.mean_molecular_weight*vb)
    hb2 = h1 + 0.5*(Pb-P1)*(vb+v1)
    return hb2-hb1
    
