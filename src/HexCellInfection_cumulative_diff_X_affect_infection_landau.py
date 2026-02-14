# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
# modeling cells being infected as points in a hexagonal coordinate system
# user specfies the number of lines and the number of cells on each line
# more lines and more cells equals more realistic-looking circle

# first convert inputs into a dataframe with cells coordinate
# use that to populate a tellurium script for stochastic Gillespie simulation


import numpy as np 
from numpy import asarray as ar,exp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import tellurium as te
from joblib import Parallel, delayed
import argparse
from pathlib import Path
import json
from dataclasses import dataclass, asdict, replace
from typing import Dict, Any, Optional
############# FUNCTION ZOO ####################################################################
###############################################################################################
@dataclass
class ModelParams:
    k_I: float = 0.05
    k_dXI: float = 1.0
    k_A: float = 0.008
    k_B: float = 0.0001
    k_XA: float = 0.0
    k_XB: float = 0.03
    k_da: float = 0.01
    k_E: float = 0.006
    k_Xmax: float = 0.4
    k_Xhalf: float = 1.0
    k_Xdeg: float = 0.001
    k_fX: float = 0.2
    k_d: float = 1.0

SCALED_KEYS = {
    "k_I", "k_A", "k_B", "k_XA", "k_XB", "k_da", "k_E", "k_Xmax", "k_Xhalf", "k_Xdeg", "k_fX"
}

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pretreat", action="store_true",
                   help="Enable pretreatment (uniform initial dX and pre-run before infection).")
    p.add_argument("--t_start", type=float, default=0.0,
                   help="Pretreatment duration (time units). Ignored if --pretreat is not set.")
    p.add_argument("--dX0", type=float, default=10.0,
                   help="Uniform initial dX value for pretreatment (default 10). Ignored if --pretreat is not set.")
    p.add_argument("--runs", type=int, default=40)
    p.add_argument("--cell_iter", type=int, default=10000)
    p.add_argument("--n_jobs", type=int, default=20)
    p.add_argument("--out_dir", type=str, default="./",
               help="Output directory for CSV files.")
    p.add_argument("--param_file", type=str, default=None,
               help="JSON file with parameter overrides (e.g., rigoe.json).")
    p.add_argument("--param_scale", type=float, default=4.1,
               help="Scale factor applied to kinetic parameters in SCALED_KEYS (default 4.1).")
    return p.parse_args()

def load_params_json(path: str) -> Dict[str, Any]:
    with open(path, "r") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("Parameter file must be a JSON object (dict).")
    return data

def override_params(base: ModelParams, updates: Dict[str, Any]) -> ModelParams:
    allowed = set(base.__dataclass_fields__.keys())
    unknown = set(updates.keys()) - allowed
    if unknown:
        raise KeyError(f"Unknown parameter keys in file: {sorted(unknown)}")
    return replace(base, **updates)
            
#%% get neighbors idx
def get_neighbors(cells,cell_idx):
    
    cell_coor_now = cells[cell_idx];
    
    X = np.transpose(np.array([0,-1,-1,0,1,1]))
    Y = np.transpose(np.array([1,1,0,-1,-1,0]))
    Z = np.transpose(np.array([-1,0,1,1,0,-1]))
    
    neighbors_og = np.column_stack((X,Y,Z))
    neighbors_now = neighbors_og + cell_coor_now
    
    neighbors_idx = []
    for vec in neighbors_now.tolist():
        idx = np.where(np.all(cells==np.array(vec),axis=1))
        if idx[0].size > 0:
            # only append cells included in the cell list
            neighbors_idx.append(idx[0][0])
            
    return neighbors_idx

#%% generate cells with hexagonal coordinates
def generate_cells(iter_num):
    cell_num = iter_num;
    
    X = np.transpose(np.array([0,0,-1,-1,0,1,1]))
    Y = np.transpose(np.array([0,1,1,0,-1,-1,0]))
    Z = np.transpose(np.array([0,-1,0,1,1,0,-1]))
    layer = np.transpose(np.array([1,1,1,1,1,1,1]))
    
    M = np.column_stack((layer,X,Y,Z))
    M_noL = np.column_stack((X,Y,Z))
    current_layer = 1
    idx_chosen = 0
    vec_chosen_box = np.array([], dtype=np.int64).reshape(0,3)
    
    #create the model neighbor matrix
    M_fixed = M[1:,1:]
    
    for i in range(cell_num):
        # get the array with the current layer only
        subM = M[M[:,0]==current_layer,:]
        # also get the array without the L
        subM_noL = subM[:,1:]
        # reset chosen cell to spawn neighbors if chosen cell id exceeds number of
        # cells in layer; also reset chosen cell data container
        if len(subM) < idx_chosen +1: 
            # first add the vector chosen box back to M
            L = np.array(([current_layer+1]*len(vec_chosen_box)))
            vec_chosen_box = np.column_stack((L,vec_chosen_box))
            M = np.concatenate((M,vec_chosen_box))
            # reset idx_chose and vec_chosen_box
            idx_chosen = 0
            vec_chosen_box = np.array([], dtype=np.int64).reshape(0,3)
            # update current_layer
            current_layer = current_layer+1
            # get the new subM and subM_noL
            subM = M[M[:,0]==current_layer,:]
            subM_noL = subM[:,1:]
            
        # get x and y and z of neighbors of chosen cell
        vec_chosen = subM_noL[idx_chosen,:]
        if np.array_equal(vec_chosen,np.array([0,0,0])):
            idx_chosen = idx_chosen + 1
            continue
        # generate neighbors of chosen cells
        neighbors_chosen = M_fixed + vec_chosen
        # remove duplicated cells
        neighbors_chosen_unq = [i for i in neighbors_chosen.tolist() if i not in M_noL.tolist()]
        neighbors_chosen_unq = np.array(neighbors_chosen_unq)       
        # add this into the chosen box and retains only unique ones
        vec_chosen_box = np.concatenate((vec_chosen_box,neighbors_chosen_unq),axis=0)
        vec_chosen_box = np.unique(vec_chosen_box,axis=0)
        # update idx_chosen
        idx_chosen = idx_chosen + 1
    
    # only retain unique id for cells 
    M =  M[:,1:] 
    M = np.unique(M,axis=0)    
    M = M.tolist()
    
    return M

#%% let's generate tellurium script
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

def getCellInfo(c,cells,result_array,t_step,names):

    IDX_ = c
    X_ = cells[c][0]
    Y_ = cells[c][1]
    Z_ = cells[c][2]
    H_ = cells[c][0]
    V_ = 2. * (np.sin(np.radians(60)) * (cells[c][1] - cells[c][2]) /3.)
    T_ = result_array[t_step,0]
    T_step = t_step
    
    # find the status of the cell at this time point    
    Df_= result_array[t_step,1+names.index('dX_{}'.format(c))]
    Deathf = result_array[t_step,1+names.index('X_{}'.format(c))]
    UI = result_array[t_step,1+names.index('UI_{}'.format(c))]
    if UI ==1:
        S_ = 0
        return (IDX_,X_,Y_,Z_,H_,V_,T_,T_step,Df_,S_,Deathf)
    I = result_array[t_step,1+names.index('I_{}'.format(c))]
    if I ==1:
        S_ = 1
        return (IDX_,X_,Y_,Z_,H_,V_,T_,T_step,Df_,S_,Deathf)
    DA = result_array[t_step,1+names.index('DA_{}'.format(c))]
    if DA ==1:
        S_ = 2
        return (IDX_,X_,Y_,Z_,H_,V_,T_,T_step,Df_,S_,Deathf)
    DB = result_array[t_step,1+names.index('DB_{}'.format(c))]
    if DB == 1:
        S_ = 3
        return (IDX_,X_,Y_,Z_,H_,V_,T_,T_step,Df_,S_,Deathf)
    DE = result_array[t_step,1+names.index('DE_{}'.format(c))]
    if DE ==1:
        S_ = 4
        return (IDX_,X_,Y_,Z_,H_,V_,T_,T_step,Df_,S_,Deathf)
    
#################################################################
#################################################################
#%% MAIN SCRIPT
args = parse_args()
out_dir = Path(args.out_dir)
out_dir.mkdir(parents=True, exist_ok=True)

pretreat = args.pretreat
t_start = args.t_start if pretreat else 0.0
dX0 = args.dX0 if pretreat else 0.0

base_params = ModelParams()  # Control default

if args.param_file is not None:
    updates = load_params_json(args.param_file)
    for k, v in updates.items():
        if not isinstance(v, (int, float)):
            raise TypeError(f"{k} must be a number, got {type(v)}")
    
    base_params = override_params(base_params, updates)
    
for run in range( args.runs):
    # generate some cells
    cells = generate_cells(args.cell_iter)
    print('cell num: {}'.format(len(cells)))
    #print('S1')
    # generate reaction object
    S = generate_reactions(cells, dX0=dX0, params = base_params, param_scale = args.param_scale)
    #print('S2')
    # set parameters; feel free to change them!
    r = te.loada(S)
    #print('S3')
    r.integrator = 'gillespie'
    r.integrator.variable_step_size = False
    r.setSeed(run)
    
    # Start the simulation
    if t_start > 0:
        _ = r.simulate(0,t_start, 2)
    # after t_start, infect the cell in the middle and then simulate again
    center = cells.index([0,0,0])
    r[f"UI_{center}"] = 0
    r[f"I_{center}"] = 1
    t_end = t_start + 500
    result = r.simulate(t_start,t_end,100)
    print('Simulaton done')
    
     
        
    #%% compile result into a dataframe
    IDX_ = []
    X_ =[]
    Y_ = []
    Z_ = []
    H_ = []
    V_ = []
    T_ = []
    T_step = [];
    S_ = []
    Df_ = []
    Deathf = []
    
    
    df = pd.DataFrame()
    out_list = []
    names = r.getFloatingSpeciesIds()
    result_array = np.array(result) # need to convert NamedArray to array for pickable obj
    for t_step in range(len(result)):
        
        sub_out_list = Parallel(n_jobs=args.n_jobs)(delayed(getCellInfo)(c,cells,result_array,t_step,names) 
                                      for c in range(len(cells)) ) 
        out_list = out_list + sub_out_list    

    df['IDX'] = [out[0] for out in out_list]
    df['X'] = [out[1] for out in out_list]
    df['Y'] = [out[2] for out in out_list]
    df['Z'] = [out[3] for out in out_list]
    df['H'] = [out[4] for out in out_list]
    df['V'] = [out[5] for out in out_list]
    df['T'] = [out[6] for out in out_list]
    df['T'] = df['T'] - t_start
    df['T_step'] = [out[7] for out in out_list]
    df['Df'] = [out[8] for out in out_list]
    df['S']  = [out[9] for out in out_list]
    df['Deathf'] = [out[10] for out in out_list]
    print('done with run {}'.format(run+1))
    out_path = out_dir / f"run{run}_pretreat{int(pretreat)}_dX0{dX0}_tstart{t_start}.csv"
    df.to_csv(out_path, index=None)
    used_param_path = out_dir / f"run{run}_params_used.json"
    with open(used_param_path, "w") as f:
        json.dump(asdict(base_params), f, indent=2, sort_keys=True)

    

