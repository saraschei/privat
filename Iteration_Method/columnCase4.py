
import numpy as np
from matplotlib.figure import Figure


class ColumnCase4():
    
    def __init__(self, nx, mx, my, fck, fyk, gammac, gammas, concModel, Asx1, Asx2, Asy1, Asy2,\
                  b, h, cx1, cx2, cy1, cy2, n, Es, epsud, beta, maxIt, ms, alfa1, alfa2, alfa3):
        
        # Input
    
        self.__nx =  nx * 1000 * np.cos(alfa2) / h            #N/mm
        
        self.__mx = mx * 1e6 / b             #kNm * 10^6 /b =  Nmm/mm
        
        self.__my = my * 1e6 / h             #kNm * 10^6 /h =  Nmm/mm
        
        self.__ms = ms
        self.__alfa1 = abs(alfa1)
        self.__alfa1Orig = alfa1                #value of alfa1 with its original sign
        self.__alfa2 = alfa2
        self.__alfa3 = alfa3
        
        
        self.__fck = fck
        self.__fyk = fyk
        self.__gammac = gammac
        self.__gammas = gammas
        self.__concModel = concModel       #Concrete Model, 1 or 2
        
        #Reinforcement: 1: lower rf, 2: upper rf
        self.__Asx1 = Asx1 / b      #mm^2/mm
        self.__Asx2 = Asx2 / b      #mm^2/mm
        self.__Asy1 = Asy1 / h      #mm^2/mm
        self.__Asy2 = Asy2 / h      #mm^2/mm
        
        
        self.__b = b        #mm
        self.__h = h        #mm
        self.__cx1 = cx1      #lower cover x-direction
        self.__cx2 = cx2      #upper cover x-direction
        self.__cy1 = cy1      #lower cover y-direction
        self.__cy2 = cy2      #upper cover y-direction
        self.__n = n        #number of layers
        
        
        self.__Es = Es
        self.__epsud = epsud    #reinforcement design ultimate strain
        
        self.__beta = beta
        self.__maxIt = maxIt
    
    
        #Basic Calculations based on the input
        
            #Matrix diplay of input (All terms are in N and mm)
        self.R = np.array([[self.__nx], [self.__ms]])
        self.ROrig = np.array([[nx], [mx], [my]])      # original forces
        self.Asx = np.array([self.__Asx1, self.__Asx2])
        self.Asy = np.array([self.__Asy1, self.__Asy2])
        
            #Empty lists for plotting forces
        self.nxSList = []
        self.mxSList = []
        self.mySList = []
        self.loopList = []
        
            #Material Charcteristics
        self.fcd = 0.85*fck/gammac
        self.fcm = fck + 8
        self.Ecm = 22000*((self.fcm/10)**0.3)
        self.Ecd = self.Ecm
        self.fyd = fyk/gammas
        
            #Geometry data
        self.bc = b/2 - h/2 * np.tan(alfa2)
        self.bz = (b-self.bc) * np.cos(alfa2)
        self.deltah = 2*self.bz/n
        
        self.convergence = True     #Initialize convergence as True
        
    def iteration(self):
        
        fck = self.__fck
        concModel= self.__concModel
        b = self.__b
        h = self.__h
        cx1 = self.__cx1
        cx2 = self.__cx2
        cy1 = self.__cy1
        cy2 = self.__cy2
        n = self.__n
        Es = self.__Es
        epsud = self.__epsud
        beta = self.__beta
        maxIt = self.__maxIt
        alfa1 = self.__alfa1
        alfa2 = self.__alfa2
        
        KcNew = 0
        KsxNew = 0
        
        
        #Functions
        def transformTo6x6(mat1, mat2, mat3, mat4):     #from 4 3x3 Matrix to one 6x6 Matrix
            row1 = np.hstack((mat1[0], mat2[0]))
            row2 = np.hstack((mat1[1], mat2[1]))
            row3 = np.hstack((mat1[2], mat2[2]))
            row4 = np.hstack((mat3[0], mat4[0]))
            row5 = np.hstack((mat3[1], mat4[1]))
            row6 = np.hstack((mat3[2], mat4[2]))
            
            mat6x6 = np.vstack((row1, row2, row3, row4, row5, row6))
            return mat6x6
        
        
        
        R = self.R
        Asx = self.Asx
        Asy = self.Asy
        
        
        #Geometry and reinforcement positions
        
        bc = self.bc
        bz = self.bz
        deltah = self.deltah
        
        if alfa2==0:                                 #layer horizontal width
            deltahVert = np.inf                        
        else: 
            deltahVert = deltah/np.sin(alfa2)
            
        
        #Reinforcement positions
            #Asx1
        if cx1 < deltahVert:
            Asx1Start = 0
        else:
            Asx1Start = int(cx1/deltahVert)
        
        Asx1Finish = Asx1Start + int(b * np.cos(alfa2) / deltah)
        if Asx1Finish==n:
            Asx1Finish = n-1    
        
            #Asx2
        Asx2Start = int((h-cx2)/deltahVert)
        
        if cx2 > deltahVert:
            Asx2Finish = Asx2Start + int(b * np.cos(alfa2)/deltah)
        else:
            Asx2Finish = n-1
        if Asx2Finish == n:
            Asx2Finish = n-1
            
            
            #Asy2
        if cy2 < deltah/np.cos(alfa2):
            Asy2Start = 0
        else:
            Asy2Start = int(cy2*np.cos(alfa2)/deltah)
            
        Asy2Finish = Asy2Start + int(h * np.sin(alfa2)/deltah)
        if Asy2Finish == n:
            Asy2Finish = n-1
        
        
            #Asy1
        Asy1Start = int((b-cy1)*np.cos(alfa2)/deltah)
        
        if cy1 > deltah/np.cos(alfa2):
            Asy1Finish = Asy1Start + int(h * np.sin(alfa2) / deltah)
        else:
            Asy1Finish = n-1
        if Asy1Finish==n:
            Asy1Finish = n-1
        
        
        
        #Convergence values (differences are taken in consideration)
        RN = R[0][0]        # Axial Force
        RM = R[1][0]        # Moment
        
        if RN !=0:
            RNAbs = abs(RN)
            betaN = RNAbs * beta  
            
            if RM != 0:
                RMAbs = abs(RM)
                betaM = RMAbs/1000 * beta   
                
            else:       #RM==zeros 
                betaM = betaN       
                
        else:           #RN == zeros
            RMAbs = abs(RM)
            betaM = RMAbs/1000 * beta       
            betaN = betaM                   
        
        
        fcd = self.fcd
        Ecd = self.Ecd
        
            #Parabola-rectangle diagram
        
        if concModel == 1:
            if fck<=50:
                nc = 2
                epsc2 = 2/1000
                epscu2 = 3.5/1000
            else:
                nc = 1.4+23.4*((90-fck)/100)**4
                epsc2 = (2+0.085*(fck-50)**0.53)/1000
                epscu2 = (2.6+35*((90-fck)/100)**4)/1000
                
            #Bilinear stress-strain relation
        elif concModel == 2:
            if fck<=50:
                nc = 2
                epsc3 = 1.75/1000
                epscu3 = 3.5/1000
            else:
                nc = 1.4+23.4*((90-fck)/100)**4
                epsc3 = (1.75+0.55*((fck-50)/40))/1000
                epscu3 = (2.6+35*((90-fck)/100)**4)/1000
        
        #In case of pure compression
        if RM == 0 and RN < 0:
            if concModel == 1:
                epscu2 = epsc2
            else:
                epscu3 = epsc3
        
        #Reinforcement Characteristics
       
        fyd = self.fyd
        
        #Distance between layers and the middle plane of the shell
        zc = np.zeros(n)        # 1Xn matrix
        for i in range(n):
            zc[i] = (2*bz/n)*i + (2*bz/n)/2 - bz
        
        
            
        
        
            
        #Main loop
        
        start = 0                       #to separate the first loop from the rest
        loop = 0                        #number of loops
        loopList = self.loopList        #list of iterations
        devMax = 2 * beta               #maximum relative deviation
        diffNMax = 2 * betaN            #maximum difference for N
        diffMMax = 2* betaM             #maximum difference for M
        
        #Values for internal Forces
        SMat = np.zeros((2,maxIt))
        
        while devMax > beta or diffNMax > betaN or diffMMax > betaM:
            
            loopList.append(loop)
            
            #Concrete stiffness Matrix
            Cc0 = Ecd
            
            if start ==0:
                Kc = np.zeros((2,2))
                for i in range(n):
                    if abs(zc[i])< bc*np.cos(alfa2):
                        KcTemp = (h/np.cos(alfa2)) *  np.array([[Cc0, -zc[i]*Cc0], [-zc[i]*Cc0, ((zc[i])**2)*Cc0]])  #temporary Kc
                    else:
                        if alfa2==0:
                             KcTemp = ((bz-abs(zc[i]))*h) *  np.array([[Cc0, -zc[i]*Cc0], [-zc[i]*Cc0, ((zc[i])**2)*Cc0]])  #temporary Kc
                        else:
                            KcTemp = ((bz-abs(zc[i]))*(np.tan(alfa2)+1/np.tan(alfa2))) *  np.array([[Cc0, -zc[i]*Cc0], [-zc[i]*Cc0, ((zc[i])**2)*Cc0]])  #temporary Kc
                    Kc = Kc + KcTemp
                Kc = deltah * np.cos(alfa2)/h * Kc      
     
            else:
                Kc = KcNew
                
                
            #Reinforcement stiffness Matrix
            
            Ksx1 = np.zeros((2,2))
            Ksx2 = np.zeros((2,2))
            Ksy1 = np.zeros((2,2))
            Ksy2 = np.zeros((2,2))
           
            
            if start ==0:
                Csx10 = Es
                Csx20 = Es
                Csy10 = Es
                Csy20 = Es
                
                #Asx1
                for i in range(Asx1Start, Asx1Finish+1):
                    
                    Ksx0Temp = np.array([[Csx10, -zc[i]*Csx10], [-zc[i]*Csx10, (zc[i]**2)*Csx10]])
                    Ksx10 = Asx[0]/(Asx1Finish-Asx1Start+1) * Ksx0Temp
                   
                    Ksx1 = Ksx1 + Ksx10
               
                #Asx2
                for i in range(Asx2Start, Asx2Finish+1):
                    
                    Ksx0Temp = np.array([[Csx20, -zc[i]*Csx20], [-zc[i]*Csx20, (zc[i]**2)*Csx20]])
                    Ksx20 = Asx[1]/(Asx2Finish-Asx2Start+1) * Ksx0Temp
                   
                    Ksx2 = Ksx2 + Ksx20
                
                #Asy1
                for i in range(Asy1Start, Asy1Finish+1):
                    
                    Ksy0Temp = np.array([[Csy10, -zc[i]*Csy10], [-zc[i]*Csy10, (zc[i]**2)*Csy10]])
                    Ksy10 = Asy[0]/(Asy1Finish-Asy1Start+1) * Ksy0Temp
                   
                    Ksy1 = Ksy1 + Ksy10
               
                #Asy2
                for i in range(Asy2Start, Asy2Finish+1):
                    
                    Ksy0Temp = np.array([[Csy20, -zc[i]*Csy20], [-zc[i]*Csy20, (zc[i]**2)*Csy20]])
                    Ksy20 = Asy[1]/(Asy2Finish-Asy2Start+1) * Ksy0Temp
                   
                    Ksy2 = Ksy2 + Ksy20
                   
                Ksx = Ksx1 + Ksx2 + Ksy1 + Ksy2
                   
            else:
                Ksx = KsxNew
               
            Ks = Ksx 
            
            #Total stiffness matrix
            K = Kc + Ks
            
            
            #Strains and curvatures in the middle plane of the shell
            
            cond = np.linalg.cond(K, p=1)
            
            
            try:
                #epst = np.matmul((np.linalg.inv(K)), R)
                epst = np.linalg.solve(K,R)
            except np.linalg.LinAlgError as err:
                if 'Singular matrix' in str(err):
                    epst = np.matmul(np.linalg.pinv(K, hermitian= True), R)
            
            
            if cond > 10:     
                  epst = np.matmul(np.linalg.pinv(K, hermitian= True), R)
            
            
            
            #In-plane strains and stresses in each concrete layer
            
            sigciMat = np.zeros((1,n))  #Matrix for the collection of stresses in gobal direction 
            epsciMat = np.zeros((1,n))  #Matrix for the collection of strains in global direction
            for i in range(n):
                Aci = np.array([[1, -zc[i]]])
                epsci = np.matmul(Aci, epst)        #matrix 3X1
                epsciMat[:, i:i+1] = epsci
                
                
                
                
                #Principal stresses ( In principal directions)
                sigci = np.zeros((1,1))
                if concModel == 1:
                    
                    if (epsci[0][0])<0:
                        if epsci[0][0]>-epsc2:
                            sigci[0][0] = -fcd*(1-(1-(-epsci[0][0]/epsc2))**nc)
                        elif epsci[0][0]<=-epsc2 and epsci[0][0]>=-epscu2:
                            sigci[0][0] = -fcd
                        else:
                            sigci[0][0] = 0
                    else:
                        sigci[0][0] = 0
                elif concModel == 2:
                    
                    if (epsci[0][0])<0:
                        if epsci[0][0]>-epsc3:
                            sigci[0][0] = fcd*epsci[0][0]/epsc3
                        elif epsci[0][0]<=-epsc3 and epsci[0][0]>-epscu3:
                            sigci[0][0] = -fcd
                        else:
                            sigci[0][0] = 0
                    else:
                        sigci[0][0] = 0
                        
               
                
                #Collection of principal stresses in one Matrix 1xn
                sigciMat[:,i:i+1] = sigci
            
            epsciMax = np.amin(epsciMat)
            sigciMax = np.amin(sigciMat)
            
            if concModel == 1:
                concUtilRatio = abs(epsciMax/epscu2)      #Utilization ratio for concrete
            else:
                concUtilRatio = abs(epsciMax/epscu3)
            
            
            
            #Concrete contribution to internal force vector
            
            Sc = np.zeros((2,1))        #main iternal concrete vector
            
            for i in range(n):
        
                if abs(zc[i]) < bc*np.cos(alfa2):
                    Ai = deltah * h/np.cos(alfa2)
                else:
                    if alfa2==0:
                        Ai = deltah * (bz-abs(zc[i]))*h
                    else:
                        Ai = deltah * ((bz-abs(zc[i]))*(np.tan(alfa2)+1/np.tan(alfa2)))
                                                                  #change
                sc = np.zeros((2,1))    #internal concrete vector for each layer
                sc[0][0] = Ai * sigciMat[0][i]
                sc[1][0] = Ai * (-zc[i]) * sigciMat[0][i] 
                
                
                Sc = Sc + sc
                
            Sc = Sc * np.cos(alfa2)/h
        
            
            #Strains and stresses in reinforcement
            
            sigsj = np.zeros((1,1))
            sigsjMat = np.zeros((1,n))  #Stress Matrix  
            epssjMat = np.zeros((1,n)) #Strain Matrix 
            
            for j in range(n):          
                if zc[j] == 0:
                    Asj = np.zeros((1,2))
                else:
                    Asj = np.array([[1, -zc[j]]])
                    
                epssj = np.matmul(Asj,epst)     #3x1 matrix
                       
                
                epssjMat[:,j:j+1] = epssj
                
                
                #Stresses in principal directions
                
                if abs(epssj[0][0]) <= fyd/Es:
                    sigsj[0][0] =  (epssj[0][0]) * Es
                elif abs(epssj[0][0]) <= epsud:
                    sigsj[0][0] = np.sign(epssj[0][0]) * fyd
                else:
                    sigsj[0][0] = 0
                
                # no need to specify sigpsj[2][0] since it is already 0
                
               
                sigsjMat[:, j:j+1] = sigsj
                
            
            
            #Relevant reinforcement List
            
            reinfList = []
            strainList = []
            stressList = []
            if Asx[0] != 0:
                reinfList.append(Asx1Start)
                reinfList.append(Asx1Finish)
                strainList.append(epssjMat[0][Asx1Start])
                strainList.append(epssjMat[0][Asx1Finish])
                stressList.append(sigsjMat[0][Asx1Start])
                stressList.append(sigsjMat[0][Asx1Finish])
            if Asx[1] != 0:
                reinfList.append(Asx2Start)
                reinfList.append(Asx2Finish)
                strainList.append(epssjMat[0][Asx2Start])
                strainList.append(epssjMat[0][Asx2Finish])
                stressList.append(sigsjMat[0][Asx2Start])
                stressList.append(sigsjMat[0][Asx2Finish])
            if Asy[0] != 0:
                reinfList.append(Asy1Start)
                reinfList.append(Asy1Finish)
                strainList.append(epssjMat[0][Asy1Start])
                strainList.append(epssjMat[0][Asy1Finish])
                stressList.append(sigsjMat[0][Asy1Start])
                stressList.append(sigsjMat[0][Asy1Finish])
            if Asy[1] != 0:
                reinfList.append(Asy2Start)
                reinfList.append(Asy2Finish)
                strainList.append(epssjMat[0][Asy2Start])
                strainList.append(epssjMat[0][Asy2Finish])
                stressList.append(sigsjMat[0][Asy2Start])
                stressList.append(sigsjMat[0][Asy2Finish])
            
            #Maximum and Minimum values for reinforcement
            
            epssjMin = min(strainList)
             
            epssjMax = max(strainList)
                
            sigsjMin = min(stressList)
            
            sigsjMax = max(stressList)
            
            
            epssMax = max(abs(epssjMin), abs(epssjMax))         # Strain Highest absolute value
            sigsMax = max(abs(sigsjMin), abs(sigsjMax))         # Stress Highest absolute value
            
                #Maximum reinforcement tension and compression
            if epssjMax<=0 and epssjMin<0:
                epssjTens = 0.0
                epssjComp = epssjMin
            elif epssjMax>0 and epssjMin>=0:
                epssjTens = epssjMax
                epssjComp = 0.0
            elif epssjMax>0 and epssjMin<0:
                epssjTens= epssjMax
                epssjComp = epssjMin
            else:
                epssjTens= 0.0
                epssjComp = 0.0
                
            if sigsjMax<=0 and sigsjMin<0:
                sigsjTens = 0.0
                sigsjComp = sigsjMin
            elif sigsjMax>0 and sigsjMin>=0:
                sigsjTens = sigsjMax
                sigsjComp = 0.0
            elif sigsjMax>0 and sigsjMin<0:
                sigsjTens= sigsjMax
                sigsjComp = sigsjMin
            else:
                sigsjTens= 0.0
                sigsjComp = 0.0
                
                
            #Reinforcement contribution to internal force vector
            
            Ss = np.zeros((2,1))
            
            #Asx1
            SsAsx1= np.zeros((2,1))
            for j in range (Asx1Start, Asx1Finish+1):
                ss = np.zeros((2,1))
                ss[0][0] = Asx[0]/(Asx1Finish-Asx1Start+1) * sigsjMat[0][j]
                ss[1][0] = -zc[j] *Asx[0]/(Asx1Finish-Asx1Start+1) * sigsjMat[0][j]
                
                SsAsx1 = SsAsx1 + ss
            
            #Asx2
            SsAsx2= np.zeros((2,1))
            for j in range (Asx2Start, Asx2Finish+1):
                ss = np.zeros((2,1))
                ss[0][0] = Asx[1]/(Asx2Finish-Asx2Start+1) * sigsjMat[0][j]
                ss[1][0] = -zc[j] *Asx[1]/(Asx2Finish-Asx2Start+1) * sigsjMat[0][j]
                
                SsAsx2 = SsAsx2 + ss
                
            #Asy1
            SsAsy1= np.zeros((2,1))
            for j in range (Asy1Start, Asy1Finish+1):
                ss = np.zeros((2,1))
                ss[0][0] = Asy[0]/(Asy1Finish-Asy1Start+1) * sigsjMat[0][j]
                ss[1][0] = -zc[j] *Asy[0]/(Asy1Finish-Asy1Start+1) * sigsjMat[0][j]
                
                SsAsy1 = SsAsy1 + ss    
                
            #Asy2
            SsAsy2= np.zeros((2,1))
            for j in range (Asy2Start, Asy2Finish+1):
                ss = np.zeros((2,1))
                ss[0][0] = Asy[1]/(Asy2Finish-Asy2Start+1) * sigsjMat[0][j]
                ss[1][0] = -zc[j] *Asy[1]/(Asy2Finish-Asy2Start+1) * sigsjMat[0][j]
                
                SsAsy2 = SsAsy2 + ss 
                
            Ss = SsAsx1 + SsAsx2 + SsAsy1 + SsAsy2
                
                
            #Internal Force Vector
            
            S = Sc + Ss
            
            
            SMat[:,loop:loop+1] = S
            
            #Addition of values of S to the internal forces list
            
            self.nxSList.append(S[0][0] * h/(1000 * np.cos(alfa2)))   #converted to kN
            self.mxSList.append(S[1][0] * np.cos(self.__alfa1Orig) * b/1e6)    #converted to kNm
            self.mySList.append(S[1][0] * np.sin(self.__alfa1Orig) * h/1e6)    #converted to kNm
            
            #Calculation of the maximum difference between S and R
            
            Diff = R - S
            RelDev = np.zeros((2,1))    #Matrix containing the relative deviation
            RelDiff = np.zeros((2,1))   #Matrix containing the differences    
            for d in range (2):
                if R[d][0] == 0:
                    RelDev[d][0] = 0
                    RelDiff[d][0] = abs(Diff[d][0])
                else:
                    RelDev[d][0] = abs(Diff[d][0]/R[d][0])
            
                    
            devMax = np.amax(RelDev)    #Maximum realtive deviation
            diffNMax = RelDiff[0:1,0:1]  #Maximum difference for the values where values of N in R are 0
            diffMMax = RelDiff[1:2,0:1]  #Maximum difference for the values where values of M in R are 0
            
            
            
            # New Young's Modulus and Material Matrix for Concrete
            
            CcgiMat = np.zeros((1,n)) #Matrix containing all Material Matrices (for every concrete layer)
            for i in range(n):
                if epsciMat[0][i] == 0:
                    E11 = 0
                else:
                    E11 = (sigciMat[0][i])/(epsciMat[0][i])
                    
                
                
                #Material Matrix in the principal direction
                Ccpi = E11
            
                #Material Matrix in the global direction (Ccgi)
                Ccgi = Ccpi
                CcgiMat[:,i:i+1] = Ccgi     
                
            
            
                
            
            # New Young's Modulus and Material Matrix for the Reinforcement
            
            CsgxMat = np.zeros((1,n))   #Matrix collecting Ex
            for j in range(n):
                
                if epssjMat[0][j] ==0:
                    EsxNew = 0
                else:
                    EsxNew = sigsjMat[0][j]/epssjMat[0][j]
                
                    
             
                
                    #Young's Modulus global direction  
        
                Csgx = np.array([[EsxNew]])
                
                
                CsgxMat[:,j:j+1]= Csgx
                
            
                
            #New stiffness Matrix for concrete
            
            KcNew = np.zeros((2,2))
            for i in range(n):
                if abs(zc[i])< bc*np.cos(alfa2):
                    KcNewTemp =(h/np.cos(alfa2)) * np.array([[CcgiMat[0,i], -zc[i]*CcgiMat[0,i]], [-zc[i]*CcgiMat[0,i], ((zc[i])**2)*CcgiMat[0,i]]])  #temporary KcNew
                else:
                    if alfa2==0:
                        KcNewTemp = ((bz-abs(zc[i]))*h) *  np.array([[CcgiMat[0,i], -zc[i]*CcgiMat[0,i]], [-zc[i]*CcgiMat[0,i], ((zc[i])**2)*CcgiMat[0,i]]])
                    else:
                        KcNewTemp = ((bz-abs(zc[i]))*(np.tan(alfa2)+1/np.tan(alfa2))) *  np.array([[CcgiMat[0,i], -zc[i]*CcgiMat[0,i]], [-zc[i]*CcgiMat[0,i], ((zc[i])**2)*CcgiMat[0,i]]])
                KcNew = KcNew + KcNewTemp
                
            KcNew = deltah * np.cos(alfa2)/h * KcNew    
             
            
            #New stiffness Matrix for the reinforcement
            
            KsxNew = np.zeros((2,2))  
            #Asx1
            KsAsx1= np.zeros((2,2))
            for j in range (Asx1Start, Asx1Finish+1):
                Csx = CsgxMat[0,j]
                KsxTemp = np.array([[Csx, -zc[j]*Csx], [-zc[j]*Csx, (zc[j]**2)*Csx]])
                Ksx = Asx[0]/(Asx1Finish-Asx1Start+1) *  KsxTemp
                
                KsAsx1 = KsAsx1 + Ksx
                
            #Asx2
            KsAsx2= np.zeros((2,2))
            for j in range (Asx2Start, Asx2Finish+1):
                Csx = CsgxMat[0,j]
                KsxTemp = np.array([[Csx, -zc[j]*Csx], [-zc[j]*Csx, (zc[j]**2)*Csx]])
                Ksx = Asx[1]/(Asx2Finish-Asx2Start+1) *  KsxTemp
                
                KsAsx2 = KsAsx2 + Ksx
                
            #Asy1
            KsAsy1= np.zeros((2,2))
            for j in range (Asy1Start, Asy1Finish+1):
                Csx = CsgxMat[0,j]
                KsyTemp = np.array([[Csx, -zc[j]*Csx], [-zc[j]*Csx, (zc[j]**2)*Csx]])
                Ksy = Asy[0]/(Asy1Finish-Asy1Start+1) *  KsyTemp
                
                KsAsy1 = KsAsy1 + Ksy
                
            #Asy2
            KsAsy2= np.zeros((2,2))
            for j in range (Asy2Start, Asy2Finish+1):
                Csx = CsgxMat[0,j]
                KsyTemp = np.array([[Csx, -zc[j]*Csx], [-zc[j]*Csx, (zc[j]**2)*Csx]])
                Ksy = Asy[1]/(Asy2Finish-Asy2Start+1) *  KsyTemp
                
                KsAsy2 = KsAsy2 + Ksy
                
                 
            KsxNew = KsAsx1 + KsAsx2 + KsAsy1 + KsAsy2
             
        
            loop = loop + 1
            start = 1
            if loop >= maxIt:
                self.convergence = False
                break
            
            if np.all((K==0)):
                self.convergence = False
                break
             
        #End of While Loop
       
        if concUtilRatio > 1:
            self.convergence = False
            
        
            #Convergence result
        convergence = self.convergence
        
            #Internal Original Forces convert to N= KN and M= kNm
        SOrig = np.zeros((3,1))
        SOrig[0][0] = S[0][0] * h / (np.cos(alfa2) * 1000)
        SOrig[1][0] = S[1][0] * np.cos(self.__alfa1Orig)* b/1e6
        SOrig[2][0] = S[1][0] * np.sin(self.__alfa1Orig)* h/1e6
        
            #Utilization Ratio for steel
        reinfUtilRatioTens = abs(epssjTens/epsud)
        reinfUtilRatioComp = abs(epssjComp/epsud)
       
        
        
            #Results appropriately rounded
            

        
        concResult = [np.round(epsciMax *1000, 3), np.round(sigciMax, 2), np.round(concUtilRatio , 2)]
                #reinfStrain (reinforcemnt strain) is represented as per mille
        reinfStrain = [np.round(epssjTens*1000 , 3), np.round(epssjComp*1000, 3)]
        reinfStress = [np.round(sigsjTens, 2), np.round(sigsjComp, 2)]
        reinfUtilRatio = [np.round(reinfUtilRatioTens, 2), np.round(reinfUtilRatioComp, 2)]
        intForce = np.round(SOrig, 2).flatten()
        RelDiff = RelDiff.flatten()
        RelDev = RelDev.flatten()
        loop
                    
        return concResult, reinfStrain, reinfStress, reinfUtilRatio, intForce, RelDiff, RelDev, loop, epsciMat, sigciMat, sigsjMat, convergence  
        
        
        
    #Plotting 
    def plot(self, loop, epsciMat, sigciMat, sigsjMat ): 
    
            #Concrete
            
                #Strain (displayed in graph as per mille)
        
        strainpc1 = (1000 * epsciMat[0,:]).flatten()
        layers = list(range(1,self.__n+1))
        
    
        straincFig = Figure(figsize=(7,7), constrained_layout= True)      #Figure 1: concrete strain
        straincFig.suptitle('Concrete strain (\u03B5c) distribution', fontsize=10)
        chart1 = straincFig.add_subplot()
        #chart1.set_title('Concrete strain (\u03B5c) distribution', fontsize=10)
        chart1.set_xlabel('Strain \u2030 (\u03B5c)', fontsize=8)
        chart1.set_ylabel('Layers', fontsize=8)
        chart1.tick_params(axis = 'x', labelsize=8)
        chart1.tick_params(axis = 'y', labelsize=8)
        chart1.set_ylim([1,self.__n])
        chart1.plot(strainpc1, layers)
        
    
                #Stress
        
        stresspc1 = (sigciMat[0,:]).flatten()
       
        stresscFig = Figure(figsize=(7,7), constrained_layout= True)      #Figure 2: concrete stress
        stresscFig.suptitle('Concrete stress (\u03C3c) distribution', fontsize=10)
        chart2 = stresscFig.add_subplot()
        #chart2.set_title('Concrete stress (\u03C3c) distribution', fontsize=10)
        chart2.set_xlabel('Stress (\u03C3c)', fontsize=8)
        chart2.set_ylabel('Layers', fontsize=8)
        chart2.tick_params(axis = 'x', labelsize=8)
        chart2.tick_params(axis = 'y', labelsize=8)
        chart2.set_ylim([1,self.__n])
        chart2.plot(stresspc1, layers)
            
        
            #Convergence

        
        
                #Force (N)
        
        forceFig = Figure(figsize=(7,7), constrained_layout=True)       #Figure 3: N convergence
        forceFig.suptitle('N Convergence', fontsize=10)
        nxPlot = forceFig.add_subplot()
        #nxPlot.set_title('N Convergence', fontsize=10)
        nxPlot.set_xlabel('Iteration', fontsize=8)
        nxPlot.set_ylabel('N (kN)', fontsize=8)
        nxPlot.tick_params(axis = 'x', labelsize=8)
        nxPlot.tick_params(axis = 'y', labelsize=8)
        nxPlot.hlines(self.ROrig[0][0], 0, loop)
        nxPlot.plot(self.loopList, self.nxSList, 'r')
        
        
                #Moment (Mx)
             
        momentFig1 = Figure(figsize=(7,7), constrained_layout= True)     #Figure 4: Mx convergence
        momentFig1.suptitle('Mx Convergence', fontsize=10)
        mxPlot = momentFig1.add_subplot()
        #mxPlot.set_title('M Convergence', fontsize=10)
        mxPlot.set_xlabel('Iteration', fontsize=8)
        mxPlot.set_ylabel('Mx (kNm)', fontsize=8)
        mxPlot.tick_params(axis = 'x', labelsize=8)
        mxPlot.tick_params(axis = 'y', labelsize=8)
        mxPlot.hlines(self.ROrig[1][0], 0, loop)
        mxPlot.plot(self.loopList, self.mxSList, 'r')
        
        
                #Moment (My)
       
        momentFig2 = Figure(figsize=(7,7), constrained_layout= True)     #Figure 5: My convergence
        momentFig2.suptitle('My Convergence', fontsize=10)
        myPlot = momentFig2.add_subplot()
        #mxPlot.set_title('M Convergence', fontsize=10)
        myPlot.set_xlabel('Iteration', fontsize=8)
        myPlot.set_ylabel('My (kNm)', fontsize=8)
        myPlot.tick_params(axis = 'x', labelsize=8)
        myPlot.tick_params(axis = 'y', labelsize=8)
        myPlot.hlines(self.ROrig[2][0], 0, loop)
        myPlot.plot(self.loopList, self.mySList, 'r')
        
        return straincFig, stresscFig, forceFig, momentFig1, momentFig2
        


        
                
        
                                 
    
    