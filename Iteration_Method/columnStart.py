
import numpy as np
from matplotlib.figure import Figure


class ColumnStart():
    
    def __init__(self, mx, my, b, h ):
                  
        
        # Input
       
        self.__mx = mx * 1e6 / b             #kNm * 10^6 /b =  Nmm/mm
        
        self.__my = my * 1e6 / h             #kNm * 10^6 /h =  Nmm/mm
        
        self.__h = h
        
        self.__b = b
        

    
    def iterationCase(self):
        
       # Value of Ms
        
        if self.__mx>0 :
            ms = np.sqrt(self.__mx**2 + self.__my**2)
        elif self.__mx<0:
            ms = -np.sqrt(self.__mx**2 + self.__my**2)
        elif self.__mx==0 and self.__my>=0:
            ms = np.sqrt(self.__mx**2 + self.__my**2)
        else:
            ms = -np.sqrt(self.__mx**2 + self.__my**2)
            
        # Value of angle alfa1
        
        if self.__mx ==0 and self.__my==0:
            alfa1 = 0
        elif self.__mx==0:
            alfa1 = np.pi/2
        else:
            alfa1 = np.arctan(self.__my / self.__mx)
            
        alfa2 = np.pi/2 - abs(alfa1)            #Always positive
            
        # Value of angle alfa3
        
        alfa3 = np.arctan(self.__h / self.__b)
        
        
        # Choice of cases
        
        if alfa1 >= 0:
            if abs(alfa1) <= alfa3:
                case = 1
            else:
                case = 2
        else:
            if abs(alfa1) <= alfa3:
                case = 3
            else:
                case = 4
                
        
        
        return ms, alfa1, alfa2, alfa3, case
            
        
    