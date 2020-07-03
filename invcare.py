import numpy as np
def invcare(lamb,chimax):

    
    tol=0.00000000001;

    invec = np.zeros((chimax,chimax))
    for i in range(chimax):
    
     if (lamb[i,i]<tol): 
        invec[i,i]=0
     else:
        invec[i,i]=1/lamb[i,i]
          

     return invec   
