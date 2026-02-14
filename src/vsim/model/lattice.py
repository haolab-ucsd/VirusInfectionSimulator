import numpy as np
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
