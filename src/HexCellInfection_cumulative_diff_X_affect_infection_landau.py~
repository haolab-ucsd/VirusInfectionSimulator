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
from scipy import asarray as ar,exp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import tellurium as te
from joblib import Parallel, delayed

############# FUNCTION ZOO ####################################################################
###############################################################################################


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
def generate_reactions(cells):
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
        if i == cells.index([0,0,0]):
            C.append('UI_{}=1'.format(i))
            C.append('I_{}=0'.format(i))
            C.append('DA_{}=0'.format(i))
            C.append('DB_{}=0'.format(i))
            C.append('DE_{}=0'.format(i))
            C.append('X_{}=0'.format(i))
            C.append('dX_{}=0'.format(i))
            
        else:    
            C.append('UI_{}=1'.format(i))
            C.append('I_{}=0'.format(i))
            C.append('DA_{}=0'.format(i))
            C.append('DB_{}=0'.format(i))
            C.append('DE_{}=0'.format(i))
            C.append('X_{}=0'.format(i))
            C.append('dX_{}=0'.format(i))

    
    # write rate constant values, placeholder for now
    K = []
    K.append('k_I = 0.4')
    K.append('k_dXI = 1')
    K.append('k_A = 0.2')
    K.append('k_B = 0.01')
    K.append('k_XA = 0.05')
    K.append('k_XB = 0.01')
    K.append('k_da = 2')
    K.append('k_E = 0.1')
    K.append('k_Xmax = 0.001')
    K.append('k_Xhalf = 0.01')
    K.append('k_Xdeg = 0.0002')
    K.append('k_fX = 0.0001')
    K.append('k_d = 0.0001')
    
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

t_start = 0
for run in range(40):
    # generate some cells
    cells = generate_cells(10000)
    print('cell num: {}'.format(len(cells)))
    #print('S1')
    # generate reaction object
    S = generate_reactions(cells)
    #print('S2')
    # set parameters; feel free to change them!
    r = te.loada(S)
    #print('S3')
    r.integrator = 'gillespie'
    r.integrator.variable_step_size = False
    r.setSeed(run)
    r.k_I = 0.05/4.1 # infection rate
    r.k_dXI = 1 # attenuating factor for diffusion molecule dX on the reduction of infection rate  
    r.k_A  = 0.008/4.1 # death type A rate - neighbor induced was 0.0006 0.005
    r.k_B = 0.0001/4.1 # death type B rate - neighbor induced
    r.k_XA  = 0.000/4.1 # death type A rate - death factor induced
    r.k_XB  = 0.03/4.1 # death type B rate - death factor induced was 0.03
    r.k_da = 0.01/4.1 # fold increase in death type A induction strength
               # wrt being next to a infected dead cell vs a live infected one 
    r.k_E  = 0.006/4.1 # defunct rate - rate where a dead cell no longer contribute to anythin
    r.k_Xmax = 0.4/4.1 # max death factor producing rate - death factor is a protein being expressed in the cell
    r.k_Xhalf = 1.0/4.1 # k1/2 of the effect of diffusing molecule on death factor expression - death factor is a protein being expressed in the cell
    r.k_Xdeg = 0.001/4.1 # death factor degradation rate
    r.k_fX = 0.2/4.1 # production rate of the diffusing molecule was 0.1
    r.k_d = 1 # diffusion rate of the diffusion molecule to the neighbor
    
    # Start the simulation
    if t_start > 0:
        result = r.simulate(0,t_start,t_start)
    # after t_start, infect the cell in the middle and then simulate again
    exec('r.UI_{} = 0'.format(cells.index([0,0,0])))
    exec('r.I_{} = 1'.format(cells.index([0,0,0])))
    result =r.simulate(0,500,100)
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
        
        sub_out_list = Parallel(n_jobs=20)(delayed(getCellInfo)(c,cells,result_array,t_step,names) 
                                      for c in range(len(cells)) ) 
        out_list = out_list + sub_out_list    

    df['IDX'] = [out[0] for out in out_list]
    df['X'] = [out[1] for out in out_list]
    df['Y'] = [out[2] for out in out_list]
    df['Z'] = [out[3] for out in out_list]
    df['H'] = [out[4] for out in out_list]
    df['V'] = [out[5] for out in out_list]
    df['T'] = [out[6] for out in out_list]
    df['T_step'] = [out[7] for out in out_list]
    df['Df'] = [out[8] for out in out_list]
    df['S']  = [out[9] for out in out_list]
    df['Deathf'] = [out[10] for out in out_list]
    print('done with run {}'.format(run+1))      
    df.to_csv('/home/phn013/VirusHost/data/csvs/run{}_tstart{}.csv'.format(run,t_start),index=None)        
   
    

