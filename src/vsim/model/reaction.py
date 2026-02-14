from typing import Optional
from dataclasses import asdict

from vsim.model.params import ModelParams, SCALED_KEYS
from vsim.model.lattice import get_neighbors

def generate_reactions(cells, dX0: float = 0.0, params: Optional[ModelParams] = None, param_scale:float = 4.1):
    if params is None:
        raise ValueError("params must be provided to generate_reactions()")
    
    # write reactions
    S = [];
    for i in range(len(cells)):
        
        #create the neigbor sum term
        near_idx = get_neighbors(cells,i)
        #first neighbor sum block (for R2)
        ss2 = '*('
        count = 0
        for n in near_idx:
            if count+1 < len(near_idx):
                ss2 = ss2+'I_{}+DA_{}+'.format(n,n)
            else:
                ss2 = ss2+'I_{}+DA_{})'.format(n,n)
            count = count+1
        #second neighbor sum block (for R3)
        ss2_1 = '*('
        count = 0
        for n in near_idx:
            if count+1 < len(near_idx):
                ss2_1 = ss2_1+'k_A*I_{}+k_da*DA_{}+'.format(n,n)
            else:
                ss2_1 = ss2_1+'k_A*I_{}+k_da*DA_{})'.format(n,n)
            count = count+1    
        
        # first reaction
        S.append('UI_{} => I_{}; k_I*UI_{}*(1/(1+k_dXI*dX_{}))'.format(i,i,i,i)+ss2)
        # second reaction
        ss1 = 'UI_{} => DB_{}; k_B*UI_{}'.format(i,i,i)
        S.append(ss1+ss2+'+k_XB*X_{}*UI_{}'.format(i,i)+ss2) #+ss2
        # third reaction
        ss3 = 'I_{} => DA_{}; I_{}'.format(i,i,i)
        S.append(ss3+ss2_1+'+k_XA*X_{}*I_{}'.format(i,i))
        ss3_1 = 'DA_{} => DE_{} ; k_E*DA_{}'.format(i,i,i)
        S.append(ss3_1)
        # fourth reaction
        S.append('=> X_{}; k_Xmax*dX_{}/(k_Xhalf+dX_{})'.format(i,i,i))
        # fifth reaction
        S.append('X_{} => ; k_Xdeg*X_{}'.format(i,i))
        
        # six reaction: death diffused factor generation: happens when 
        # cell is infected
        S.append('=> dX_{}; k_fX*(I_{})'.format(i,i))
        # last set of reactions: diffusion of dX to neighbor
        for n in near_idx:
            # conversion rate is proportional to conc gradient
            # if dX_n < dX_i, this rate is zero
            # bc there's a reverse equation taking care of that
            S.append('dX_{} => dX_{};'.format(i,n) + 
                     'piecewise(k_d*(dX_{}-dX_{}),(dX_{}-dX_{}) > 0, 0)'.format(i,n,i,n))
    
    # write species value
    C = [];
    for i in range(len(cells)):
        # seed all cells to be uninfected for now
        C.append(f'UI_{i}=1')
        C.append(f'I_{i}=0')
        C.append(f'DA_{i}=0')
        C.append(f'DB_{i}=0')
        C.append(f'DE_{i}=0')
        C.append(f'X_{i}=0')
        C.append(f'dX_{i}={dX0}')

    
    # write rate constant values, placeholder for now
    K = []
    for k, v in asdict(params).items():
        if k in SCALED_KEYS:
            K.append(f"{k} = {v / param_scale}")
        else:
            K.append(f"{k} = {v}")
    # combine everything into the tellurium script
    # at specified time t_start, set the center cell to be infected
    
    r = ("""
        {}
        
        {}
        
        {}
        
        """.format(';'.join(S),';'.join(C),';'.join(K)))
    
    return r
