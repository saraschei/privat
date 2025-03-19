
import numpy as np
from matplotlib.figure import Figure


class Beam():
    
    def __init__(self, nx, mx, fck, fyk, gammac, gammas, concModel, Asx1, Asx2,\
                  b, h, c1, c2, n, Esx1, Esx2, epsud, beta, maxIt):
        
        # Input
    
        self.__nx =  nx * 1000 / b            #N/mm
        
        self.__mx = mx * 1e6 / b             #kNm * 10^6 /b =  Nmm/mm
        
        
        self.__fck = fck
        self.__fyk = fyk
        self.__gammac = gammac
        self.__gammas = gammas
        self.__concModel = concModel       #Concrete Model, 1 or 2
        
        #Reinforcement: 1: lower rf, 2: upper rf
        self.__Asx1 = Asx1 / b      #mm^2/mm
        self.__Asx2 = Asx2 / b      #mm^2/mm
        
        
        self.__b = b        #mm
        self.__h = h        #mm
        self.__c1 = c1      #lower cover
        self.__c2 = c2      #upper cover
        self.__n = n        #number of layers
        
        # Reinforcement Young's Modulus (N/mm2) 
        self.__Esx1 = Esx1
        self.__Esx2 = Esx2
        self.__epsud = epsud    #reinforcement design ultimate strain
        
        self.__beta = beta
        self.__maxIt = maxIt
    
    
        #Basic Calculations based on the input
        
            #Matrix diplay of input (All terms are in N and mm)
        self.R = np.array([[self.__nx], [self.__mx]])
        self.Esx = np.array([Esx1, Esx2])
        self.Asx = np.array([self.__Asx1, self.__Asx2])
        
            #Empty lists for plotting forces
        self.nxSList = []
        self.mxSList = []
        self.loopList = []
        
            #Material Charcteristics
        self.fcd = 0.85*fck/gammac
        self.fcm = fck + 8
        self.Ecm = 22000*((self.fcm/10)**0.3)
        self.Ecd = self.Ecm
        self.fyd = fyk/gammas
        
        self.convergence = True     #Initialize convergence as True
        
    def iteration(self):
        
        fck = self.__fck
        concModel= self.__concModel
        b = self.__b
        h = self.__h
        c1 = self.__c1
        c2 = self.__c2
        n = self.__n
        epsud = self.__epsud
        beta = self.__beta
        maxIt = self.__maxIt
    
        

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
        Esx = self.Esx
        Asx = self.Asx
        
        
        
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
            zc[i] = (h/n)*i + (h/n)/2 - h/2
        
        zs = np.zeros(2)        #1X2 matrix
        if c1 ==0:
            zs[0] = 0
        else:
            zs[0] = c1-(h/2)
            
        if c2 ==0:
            zs[1] = 0
        else:
            zs[1] = (h/2)-c2
            
        
        
            
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
                    KcTemp = np.array([[Cc0, -zc[i]*Cc0], [-zc[i]*Cc0, ((zc[i])**2)*Cc0]])  #temporary Kc
                    Kc = Kc + KcTemp
                Kc = (h/n)*Kc   
            else:
                Kc = KcNew
                
                
            #Reinforcement stiffness Matrix
            
            Ksx = np.zeros((2,2))
           
            
            if start ==0:
                
                for i in range(2):
                    Csx0 = Esx[i]
                   
                    Ksx0Temp = np.array([[Csx0, -zs[i]*Csx0], [-zs[i]*Csx0, (zs[i]**2)*Csx0]])
                    Ksx0 = Asx[i] * Ksx0Temp
                   
                    Ksx = Ksx + Ksx0
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
                utilratioc = abs(epsciMax/epscu2)      #Utilization ratio for concrete
            else:
                utilratioc = abs(epsciMax/epscu3)
            
            
            
            
            #Concrete contribution to internal force vector
            
            Sc = np.zeros((2,1))        #main iternal concrete vector
            
            for i in range(n):
                sc = np.zeros((2,1))    #internal concrete vector for each layer
                sc[0][0] = (h/n) * sigciMat[0][i]
                sc[1][0] = (h/n) * (-zc[i]) * sigciMat[0][i] 
                
                
                Sc = Sc + sc
        
            
            #Strains and stresses in reinforcement
            
            sigsj = np.zeros((1,1))
            sigsjMat = np.zeros((1,2))  #Stress Matrix  
            epssjMat = np.zeros((1,2)) #Strain Matrix 
            
            for j in range(2):          #j represents upper and lower reinforcement
                if zs[j] == 0:
                    Asj = np.zeros((1,2))
                else:
                    Asj = np.array([[1, -zs[j]]])
                    
                epssj = np.matmul(Asj,epst)     #3x1 matrix
                       
                
                epssjMat[:,j:j+1] = epssj
                
                
                #Stresses in principal directions
                
                if abs(epssj[0][0]) <= fyd/Esx[j]:
                    sigsj[0][0] =  (epssj[0][0]) * Esx[j]
                elif abs(epssj[0][0]) <= epsud:
                    sigsj[0][0] = np.sign(epssj[0][0]) * fyd
                else:
                    sigsj[0][0] = 0
                
                # no need to specify sigpsj[2][0] since it is already 0
                
               
                sigsjMat[:, j:j+1] = sigsj
                
            
                #Maximum and Minimum values for reinforcement (unused in GUI results)
            epssjMin = np.amin(epssjMat)
            epssjMax = np.amax(epssjMat)
            sigsjMin = np.amin(sigsjMat)
            sigsjMax = np.amax(sigsjMat)
            
                
                
            #Reinforcement contribution to internal force vector
            
            Ss = np.zeros((2,1))
            
            for j in range (2):
                ss = np.zeros((2,1))
                ss[0][0] = Asx[j] * sigsjMat[0][j]
                ss[1][0] = -zs[j] *Asx[j] * sigsjMat[0][j]
                
                Ss = Ss + ss
                
                
            #Internal Force Vector
            
            S = Sc + Ss
            
            
            SMat[:,loop:loop+1] = S
            
            #Addition of values of S to the internal forces list
            
            self.nxSList.append(S[0][0] * b/1000)   #converted to kN
            self.mxSList.append(S[1][0] * b/1e6)    #converted to kNm
            
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
            
            CsgxMat = np.zeros((1,2))   #Matrix collecting Ex
            for j in range(2):
                
                if epssjMat[0][j] ==0:
                    EsxNew = 0
                else:
                    EsxNew = sigsjMat[0][j]/epssjMat[0][j]
                
                    
             
                
                    #Young's Modulus global direction x 
        
                Csgx = np.array([[EsxNew]])
                
                
                CsgxMat[:,j:j+1]= Csgx
                
            
                
            #New stiffness Matrix for concrete
            
            KcNew = np.zeros((2,2))
            for i in range(n):
                KcNewTemp = np.array([[CcgiMat[0,i], -zc[i]*CcgiMat[0,i]], [-zc[i]*CcgiMat[0,i], ((zc[i])**2)*CcgiMat[0,i]]])  #temporary KcNew
                KcNew = KcNew + KcNewTemp
            KcNew = (h/n)*KcNew  
             
            
            #New stiffness Matrix for the reinforcement
            
            KsxNew = np.zeros((2,2))  
            for j in range(2):
                Csx = CsgxMat[0,j]
               
                KsxTemp = np.array([[Csx, -zs[j]*Csx], [-zs[j]*Csx, (zs[j]**2)*Csx]])
                Ksx = Asx[j] * KsxTemp
                 
                KsxNew = KsxNew + Ksx
             
        
            loop = loop + 1
            start = 1
            if loop >= maxIt:
                self.convergence = False
                break
            
            if np.all((K==0)):
                self.convergence = False
                break
        #End of While Loop
        
        
            #Convergence result
        convergence = self.convergence
        
            #Internal forces convert to N= KN and M= kNm
        S[0][0] = S[0][0] * b/1000
        S[1][0] = S[1][0] * b/1e6
        
            #Utilization Ratio for steel
        reinfUtilRatioAsx1 = abs(epssjMat[0][0]/epsud)
        reinfUtilRatioAsx2 = abs(epssjMat[0][1]/epsud)
        
        
            #Results appropriately rounded
            
                #epsciMax (max strain) is represente as per mille
        concResult = [(epsciMax *1000).round(3), sigciMax.round(2), utilratioc.round(2)]
                #reinfStrain (reinforcemnt strain) is represented as per mille
        reinfStrain = (epssjMat.flatten() *1000).round(3)
        reinfStress = (sigsjMat.flatten()).round(2)
        reinfUtilRatio = [reinfUtilRatioAsx1.round(2),reinfUtilRatioAsx2.round(2)]
        intForce = S.round(2).flatten()
        RelDiff = RelDiff.flatten()
        RelDev = RelDev.flatten()
        loop
                    
        return concResult, reinfStrain, reinfStress, reinfUtilRatio, intForce, RelDiff, RelDev, loop, epsciMat, sigciMat, sigsjMat, convergence  
        
        
        
    #Plotting 
    def plot(self, loop, epsciMat, sigciMat, sigsjMat ): 
    
            #Concrete
            
                #Strain (displayed in graph as per mille)
        
        strainpc1 = ((1000 * epsciMat[0,:]).round(3)).flatten()
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

                #R values in kN and kNm
        R = np.array([[self.R[0][0] * self.__b/1000], [self.R[1][0] * self.__b/1e6]])
        
                #Force (N)
        
        forceFig = Figure(figsize=(7,7), constrained_layout=True)       #Figure 3: N convergence
        forceFig.suptitle('N Convergence', fontsize=10)
        nxPlot = forceFig.add_subplot()
        #nxPlot.set_title('N Convergence', fontsize=10)
        nxPlot.set_xlabel('Iteration', fontsize=8)
        nxPlot.set_ylabel('N (kN)', fontsize=8)
        nxPlot.tick_params(axis = 'x', labelsize=8)
        nxPlot.tick_params(axis = 'y', labelsize=8)
        nxPlot.hlines(R[0][0], 0, loop)
        nxPlot.plot(self.loopList, self.nxSList, 'r')
        
        
                #Moment (M)
             
        momentFig = Figure(figsize=(7,7), constrained_layout= True)     #Figure 4: M convergence
        momentFig.suptitle('M Convergence', fontsize=10)
        mxPlot = momentFig.add_subplot()
        #mxPlot.set_title('M Convergence', fontsize=10)
        mxPlot.set_xlabel('Iteration', fontsize=8)
        mxPlot.set_ylabel('M (kNm)', fontsize=8)
        mxPlot.tick_params(axis = 'x', labelsize=8)
        mxPlot.tick_params(axis = 'y', labelsize=8)
        mxPlot.hlines(R[1][0], 0, loop)
        mxPlot.plot(self.loopList, self.mxSList, 'r')
       
        
        return straincFig, stresscFig, forceFig, momentFig
        


        
                
        
                                 
    
    