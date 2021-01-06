import numpy as np #Recommended
import random

def deviation_func(exp,model,deviationType,moleFrac = False):
    if moleFrac and exp > 0.5:
        exp = 1 - exp
        model = 1 - model
    if deviationType == 'ARD':
            deviation = np.abs((model - exp) / exp) * 100
    if deviationType == 'RD':
            deviation = (model - exp) / exp * 100
    elif deviationType == 'AD':
        deviation = np.abs(model - exp)
        
    return deviation


def PBubble_comparison(Case,expT,expP,expComposition,deviationType,Pini = 1.0):
    """
    deviation = PBubble_comparison(expT,expP,expMoles,deviationType)
    Inputs:
    expT: Experimental temperature (K)
    expP: Experimental bubble pressure (bar)
    expComposition: Experimental feed composition (mole)
    deviationType: Type of deviation. Must be a string. Either "AD" (absolute deviation) or "ARD" (absolute relative deviation)
    Absolute deviation = abs(model - exp_data)
    Absolute relative deviation = abs( (model - exp_data) / exp_data)
    
    Outputs:
    deviation: Deviation in either (%) (if deviationType == "ARD") or (bar) (if deviationType == "AD")
    """
    if not isinstance(deviationType,str):
        raise SyntaxError('deviationType must be a string')    
   
    if not isinstance(expT,(int,float,list,np.ndarray)):
        raise SyntaxError('expT must be numeric')
        
    if not isinstance(expP,(int,float,list,np.ndarray)):
        raise SyntaxError('expP must be numeric')
        
    if isinstance(expT,(float,int)):
        expT = [expT]
        
    if isinstance(expP,(float,int)):
        expP = [expP]
    
    for element in expT:
        if not isinstance(element,(float,int)):
            raise SyntaxError('Each element of expT must be numeric')
    
    for element in expP:
        if not isinstance(element,(float,int)):
            raise SyntaxError('Each element of expP must be numeric')
    
    if len(expP) != len(expT):
        raise SyntaxError('expT and expP must be of same length')
    
    deviation = np.zeros(len(expP))
    P = np.zeros(len(expP))
    
    for i in range(0,len(expT)):
        [P[i], LnK, ierr] = Case.PBubble(expT[i], expComposition, Pini)
        deviation[i] = deviation_func(expP[i],P[i],deviationType)
    
    if len(deviation) == 1:
        deviation = deviation[0]
    return deviation

def TBubble_comparison(Case,expT,expP,expComposition,deviationType,Tini = 400):
    """
    deviation = PBubble_comparison(expT,expP,expMoles,deviationType)
    Inputs:
    expT: Experimental bubble temperature (K)
    expP: Experimental pressure (bar)
    expComposition: Experimental feed composition (mole)
    deviationType: Type of deviation. Must be a string. Either "AD" (absolute deviation) or "ARD" (absolute relative deviation)
    Absolute deviation = abs(model - exp_data)
    Absolute relative deviation = abs( (model - exp_data) / exp_data)
    
    Outputs:
    deviation: Deviation in either (%) (if deviationType == "ARD") or (K) (if deviationType == "AD")
    """
    if not isinstance(deviationType,str):
        raise SyntaxError('deviationType must be a string')    
   
    if not isinstance(expT,(int,float,list,np.ndarray)):
        raise SyntaxError('expT must be numeric')
        
    if not isinstance(expP,(int,float,list,np.ndarray)):
        raise SyntaxError('expP must be numeric')
        
    if isinstance(expT,(float,int)):
        expT = [expT]
        
    if isinstance(expP,(float,int)):
        expP = [expP]
    
    for element in expT:
        if not isinstance(element,(float,int)):
            raise SyntaxError('Each element of expT must be numeric')
    
    for element in expP:
        if not isinstance(element,(float,int)):
            raise SyntaxError('Each element of expP must be numeric')
    
    if len(expP) != len(expT):
        raise SyntaxError('expT and expP must be of same length')
    
    deviation = np.zeros(len(expP))
    T = np.zeros(len(expP))
    
    for i in range(0,len(expT)):
        [T[i], LnK, ierr] = Case.TBubble(expP[i], expComposition, Tini)
        deviation[i] = deviation_func(expT[i],T[i],deviationType)
    
    if len(deviation) == 1:
        deviation = deviation[0]
    
    return deviation


def LiqRho_comparison(Case, expT, expRho, expComposition, deviationType, Pini=1):
    """
    deviation = LiqRho_comparison(expT,expP,expMoles,deviationType)
    Inputs:
    expT: Experimental temperature (K)
    expRho: Experimental liquid density (mol/L)
    expComposition: Experimental feed composition (mole)
    deviationType: Type of deviation. Must be a string. Either "AD" (absolute deviation) or "ARD" (absolute relative deviation)
    Absolute deviation = abs(model - exp_data)
    Absolute relative deviation = abs( (model - exp_data) / exp_data)

    Outputs:
    deviation: Deviation in either (%) (if deviationType == "ARD") or (K) (if deviationType == "AD")
    """
    if not isinstance(deviationType, str):
        raise SyntaxError('deviationType must be a string')

    if not isinstance(expT, (int, float, list,np.ndarray)):
        raise SyntaxError('expT must be numeric')

    if not isinstance(expRho, (int, float, list,np.ndarray)):
        raise SyntaxError('expRho must be numeric')

    if isinstance(expT, (float, int)):
        expT = [expT]

    if isinstance(expRho, (float, int)):
        expRho = [expRho]

    for element in expT:
        if not isinstance(element, (float, int)):
            raise SyntaxError('Each element of expT must be numeric')

    for element in expRho:
        if not isinstance(element, (float, int)):
            raise SyntaxError('Each element of expRho must be numeric')

    if len(expRho) != len(expT):
        raise SyntaxError('expT and expRho must be of same length')

    deviation = np.zeros(len(expRho))
    Rho = np.zeros(len(expRho))

    for i in range(0, len(expT)):
        Rho[i] = Case.LiqRho(expT[i], expComposition, Pini)
        deviation[i] = deviation_func(expRho[i], Rho[i], deviationType)

    if len(deviation) == 1:
        deviation = deviation[0]

    return deviation


def PDew_comparison(Case,expT,expP,expComposition,deviationType,Pini = 1.0):
    """
    deviation = PDew_comparison(expT,expP,expMoles,deviationType)
    Inputs:
    expT: Experimental temperature (K)
    expP: Experimental dew pressure (bar)
    expComposition: Experimental feed composition (mole)
    deviationType: Type of deviation. Must be a string. Either "AD" (absolute deviation) or "ARD" (absolute relative deviation)
    Absolute deviation = abs(model - exp_data)
    Absolute relative deviation = abs( (model - exp_data) / exp_data)
    
    Outputs:
    deviation: Deviation in either (%) (if deviationType == "ARD") or (bar) (if deviationType == "AD")
    """
    if not isinstance(deviationType,str):
        raise SyntaxError('deviationType must be a string')    
   
    if not isinstance(expT,(int,float,list,np.ndarray)):
        raise SyntaxError('expT must be numeric')
        
    if not isinstance(expP,(int,float,list,np.ndarray)):
        raise SyntaxError('expP must be numeric')
        
    if isinstance(expT,(float,int)):
        expT = [expT]
        
    if isinstance(expP,(float,int)):
        expP = [expP]
    
    for element in expT:
        if not isinstance(element,(float,int)):
            raise SyntaxError('Each element of expT must be numeric')
    
    for element in expP:
        if not isinstance(element,(float,int)):
            raise SyntaxError('Each element of expP must be numeric')
    
    if len(expP) != len(expT):
        raise SyntaxError('expT and expP must be of same length')
    
    deviation = np.zeros(len(expP))
    P = np.zeros(len(expP))
    
    for i in range(0,len(expT)):
        [P[i], LnK, ierr] = Case.PDew(expT[i], expComposition, Pini)
        deviation[i] = deviation_func(expP[i],P[i],deviationType)
    
    if len(deviation) == 1:
        deviation = deviation[0]
    
    return deviation

def TDew_comparison(Case,expT,expP,expComposition,deviationType,Tini = 400):
    """
    deviation = TDew_comparison(expT,expP,expMoles,deviationType)
    Inputs:
    expT: Experimental dew temperature (K)
    expP: Experimental pressure (bar)
    expComposition: Experimental feed composition (mole)
    deviationType: Type of deviation. Must be a string. Either "AD" (absolute deviation) or "ARD" (absolute relative deviation)
    Absolute deviation = abs(model - exp_data)
    Absolute relative deviation = abs( (model - exp_data) / exp_data)
    
    Outputs:
    deviation: Deviation in either (%) (if deviationType == "ARD") or (K) (if deviationType == "AD")
    """
    if not isinstance(deviationType,str):
        raise SyntaxError('deviationType must be a string')    
   
    if not isinstance(expT,(int,float,list,np.ndarray)):
        raise SyntaxError('expT must be numeric')
        
    if not isinstance(expP,(int,float,list,np.ndarray)):
        raise SyntaxError('expP must be numeric')
        
    if isinstance(expT,(float,int)):
        expT = [expT]
        
    if isinstance(expP,(float,int)):
        expP = [expP]
    
    for element in expT:
        if not isinstance(element,(float,int)):
            raise SyntaxError('Each element of expT must be numeric')
    
    for element in expP:
        if not isinstance(element,(float,int)):
            raise SyntaxError('Each element of expP must be numeric')
    
    if len(expP) != len(expT):
        raise SyntaxError('expT and expP must be of same length')
    
    deviation = np.zeros(len(expP))
    T = np.zeros(len(expP))
    
    for i in range(0,len(expT)):
        [T[i], LnK, ierr] = Case.TDew(expP[i], expComposition, Tini)
        deviation[i] = deviation_func(expT[i],T[i],deviationType)
    
    if len(deviation) == 1:
        deviation = deviation[0]
    
    return deviation

def BinaryVLE_Comparison(Case, expT, expP, expX, expY, deviationType, expZ = 1):
    N = len(expT)
    
    if isinstance(expZ, (list)):
        automaticFeedComp = False
    else:
        automaticFeedComp = True

    deviation_x = []
    deviation_y = []

    for n in range(0, N):
        T = expT[n]
        P = expP[n]
        y = expY[n]
        x = expX[n]

        if automaticFeedComp == False:
            z = expZ[n]
            nfas, PhaseFrac, PhaseComp, PhaseType, ierr = Case.PTFlash(T,P,z)
            if nfas != 2:
                deviation_x.append(None)
                deviation_y.append(None)
            elif not ((PhaseType[0] == -1 and PhaseType[1] == 1) or (PhaseType[0] == 1 and PhaseType[1] == -1)):
                deviation_x.append(None)
                deviation_y.append(None)
            else:
                if PhaseType[0] == 1:
                    deviation_x.append(deviation_func(x,PhaseComp[0][0],deviationType,True))
                    deviation_y.append(deviation_func(y,PhaseComp[1][0],deviationType,True))
                else:
                    deviation_x.append(deviation_func(x,PhaseComp[1][0],deviationType,True))
                    deviation_y.append(deviation_func(y,PhaseComp[0][0],deviationType,True))
        else:
            k = 0.5
            noScenario = True
            for i in range(0,100):
                z = [k,1-k]
                k = random.random()
                nfas, PhaseFrac, PhaseComp, PhaseType, ierr = Case.PTFlash(T,P,z)
                if nfas == 2:
                    noScenario = False
                    break
            if not noScenario:
                if not ((PhaseType[0] == -1 and PhaseType[1] == 1) or (PhaseType[0] == 1 and PhaseType[1] == -1)):
                    deviation_x.append(None)
                    deviation_y.append(None)
                else:
                    if PhaseType[0] == 1:
                        deviation_x.append(deviation_func(x,PhaseComp[0][0],deviationType,True))
                        deviation_y.append(deviation_func(y,PhaseComp[1][0],deviationType,True))
                    else:
                        deviation_x.append(deviation_func(x,PhaseComp[1][0],deviationType,True))
                        deviation_y.append(deviation_func(y,PhaseComp[0][0],deviationType,True))
            else:
                deviation_x.append(None)
                deviation_y.append(None)

    return np.array(deviation_x), np.array(deviation_y)
        







