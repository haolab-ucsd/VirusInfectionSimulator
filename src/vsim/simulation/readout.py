import numpy as np
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

