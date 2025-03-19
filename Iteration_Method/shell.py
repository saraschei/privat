
import numpy as np
from matplotlib.figure import Figure


class Shell():
    
    def __init__(self, nx, ny, nxy, mx, my, mxy, fck, fyk, gammac, gammas, concModel, Asx1, Asx2,\
                 Asy1, Asy2, h, c1, c2, n, Esx1, Esx2, Esy1, Esy2, epsud, beta, maxIt, vc):
        
        # Input
    
        self.__nx = nx              #kN/m  = N/mm
        self.__ny = ny              #kN/m  = N/mm
        self.__nxy = nxy            #kN/m  = N/mm
        
        self.__mx = mx * 1000        #kNm * 1000  =  Nmm/mm
        self.__my = my * 1000        #kNm * 1000  =  Nmm/mm
        self.__mxy = mxy * 1000      #kNm * 1000  =  Nmm/mm
        
        
        self.__fck = fck
        self.__fyk = fyk
        self.__gammac = gammac
        self.__gammas = gammas
        self.__concModel = concModel       #Concrete Model, 1 or 2
        
        #Reinforcement: 1: lower rf, 2: upper rf
        self.__Asx1 = Asx1 / 1000      #mm^2/mm
        self.__Asx2 = Asx2 / 1000      #mm^2/mm
        self.__Asy1 = Asy1 / 1000      #mm^2/mm
        self.__Asy2 = Asy2 / 1000      #mm^2/mm
        self.__reinfAngle1 = 0
        self.__reinfAngle2 = 0
        
        
        self.__h = h        #mm
        self.__c1 = c1      #lower cover
        self.__c2 = c2      #upper cover
        self.__n = n        #number of layers
        
        # Reinforcement Young's Modulus (N/mm2) 
        self.__Esx1 = Esx1
        self.__Esx2 = Esx2
        self.__Esy1 = Esy1
        self.__Esy2 = Esy2
        self.__epsud = epsud    #reinforcement design ultimate strain
        
        self.__beta = beta
        self.__maxIt = maxIt
        self.__vc = vc          #poisson's ratio for concrete 
    

        #Basic Calculations based on the input
        
            #Matrix diplay of input (All terms are in N and mm)
        self.R = np.array([[self.__nx], [self.__ny], [self.__nxy], [self.__mx], [self.__my], [self.__mxy]])
        self.Esx = np.array([Esx1, Esx2])
        self.Esy = np.array([Esy1, Esy2])
        self.Asx = np.array([self.__Asx1, self.__Asx2])
        self.Asy = np.array([self.__Asy1, self.__Asy2])
        self.alfa = np.array([self.__reinfAngle1, self.__reinfAngle2])
        
            #Empty lists for plotting forces
        self.nxSList = []
        self.nySList = []
        self.nxySList = []
        self.mxSList = []
        self.mySList = []
        self.mxySList = []
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
        h = self.__h
        c1 = self.__c1
        c2 = self.__c2
        n = self.__n
        epsud = self.__epsud
        beta = self.__beta
        maxIt = self.__maxIt
        vc = self.__vc
        
        #Initialization of updated concrete and reinforcement stiffness matrices
        KcNew = 0
        KsxNew = 0
        KsyNew = 0
        
        
        #Functions
        def transformTo6x6(mat1, mat2, mat3, mat4): #from 4 3x3 Matrix to one 6x6 Matrix
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
        Esy = self.Esy
        Asx = self.Asx
        Asy = self.Asy
        alfa = self.alfa
        
        
        #Convergence values (differences are taken in consideration)
        RN = R[0:3, 0:1]
        RM = R[3:6, 0:1] 
        
        if np.any((RN != 0)):
            RNAbs = abs(RN)
            RNMin = min([i for i in RNAbs if i > 0])
           
            betaN = RNMin * beta
            
            if np.any((RM != 0)):
                RMAbs = abs(RM)
                RMMin = min([i for i in RMAbs if i > 0])
                betaM = RMMin/1000 * beta
                
            else:       #RM==zeros 
                betaM = betaN
                
        else:           #RN==zeros
            RMAbs = abs(RM)
            RMMin = min([i for i in RMAbs if i > 0])
            betaM = RMMin/1000 * beta
            betaN = betaM
        
        
        #Material characteristics variable names 
        fcd = self.fcd
        fcm = self.fcm
        Ecm = self.Ecm
        Ecd = self.Ecd
        fyd = self.fyd
        
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
        
        #In case of pure compression(no moment and at least on negative value for axial forces)
        if np.all((RM == 0)) and np.any((RN < 0)):
            if concModel == 1:
                epscu2 = epsc2
            else:
                epscu3 = epsc3
        
        
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
        
        start = 0           #to separate the first loop from the rest
        loop = 0            #number of loops
        loopList = self.loopList        #list of iterations
        devMax = 2 * beta               #maximum relative deviation
        diffNMax = 2 * betaN            #maximum N difference
        diffMMax = 2* betaM             #maximum M difference
        
        while devMax > beta or diffNMax > betaN or diffMMax > betaM:
            
            loopList.append(loop)
            
            #Concrete stiffness Matrix
            Cc0 = (Ecd/(1-vc**2)) * np.array([[1, vc, 0], [vc, 1, 0], [0, 0, (1-vc)/2]])
            
            if start ==0:
                Kc = np.zeros((6,6))
                for i in range(n):
                    KcTemp = transformTo6x6(Cc0, -zc[i]*Cc0, -zc[i]*Cc0, ((zc[i])**2)*Cc0)  #temporary Kc
                    Kc = Kc + KcTemp
                Kc = (h/n)*Kc   
            else:
                Kc = KcNew

                
            #Reinforcement stiffness Matrix
            
            Ksx = np.zeros((6,6))
            Ksy = np.zeros((6,6))
            
            if start ==0:
                
                for i in range(2):
                    Csx0 = np.array([[Esx[i], 0, 0], [0, 0, 0], [0, 0, 0]])
                    Csy0 = np.array([[0, 0, 0], [0, Esy[i], 0], [0, 0, 0]])
                    Ksx0Temp = transformTo6x6(Csx0, -zs[i]*Csx0, -zs[i]*Csx0, (zs[i]**2)*Csx0)
                    Ksx0 = Asx[i] * Ksx0Temp
                    Ksy0Temp = transformTo6x6(Csy0, -zs[i]*Csy0, -zs[i]*Csy0, (zs[i]**2)*Csy0)
                    Ksy0 = Asy[i] * Ksy0Temp
                    
                    Ksx = Ksx + Ksx0
                    Ksy = Ksy + Ksy0

            else:
                Ksx = KsxNew
                Ksy = KsyNew
                
            Ks = Ksx + Ksy
            
            #Total stiffness matrix
            K = Kc + Ks
            
            
            #Strains and curvatures in the middle plane of the shell
            
            cond = np.linalg.cond(K, p=1)
            
            try:
                epst = np.linalg.solve(K,R)
            except np.linalg.LinAlgError as err:
                if 'Singular matrix' in str(err):
                    epst = np.matmul(np.linalg.pinv(K,  hermitian= True), R)
            
            
            if cond > 10:
                epst = np.matmul(np.linalg.pinv(K, hermitian= True), R)
        
            
            #In-plane strains and stresses in each concrete layer
            
            sigciMat = np.zeros((3,n))  #Matrix for the collection of stresses in gobal direction 
            sigpciMat = np.zeros((3,n)) #Matrix for the collection of stresses in principal direction
            epsciMat = np.zeros((3,n))  #Matrix for the collection of strains in global direction
            epspciMat = np.zeros((3,n)) #Matrix for the collection of strains in principal direction
            TepsciMat = np.zeros((3, 3*n))
            thetaci = np.zeros(n)           #list of theta angles for every concrete layer
            for i in range(n):
                Aci = np.array([[1, 0, 0, -zc[i], 0, 0], [0, 1, 0, 0, -zc[i], 0], [0, 0, 1, 0, 0, -zc[i]]])
                epsci = np.matmul(Aci, epst) #matrix 3X1
                epsciMat[:,i:i+1] = epsci
                
                epscxi = epsci[0][0]
                epscyi = epsci[1][0]
                gammacxyi = epsci[2][0]
                
                
                
                if (epscxi-epscyi)==0:
                    thetaci[i] = 0
                elif gammacxyi == 0:
                    thetaci[i] = 0
                else:
                    thetaci[i] = 0.5 * np.arctan(gammacxyi/(epscxi-epscyi))
                
                cc = np.cos(thetaci[i])
                sc = np.sin(thetaci[i])
                Tepsci = np.array([[cc**2, sc**2, sc*cc], [sc**2, cc**2, -sc*cc], [-2*sc*cc, 2*sc*cc, cc**2-sc**2]])
                TepsciMat[:,3*i:3*i+3] = Tepsci  
                
                #Principal strains (epspci) 
                epspci = np.matmul(Tepsci, epsci)   # Matrix 3X1
                epspciMat[:,i:i+1] = epspci
            
                
                #Principal stresses ( In principal directions)
                sigpci = np.zeros((3,1))
                if concModel == 1:
                    for j in range(2):
                        if (epspci[j][0])<0:
                            if epspci[j][0]>-epsc2:
                                sigpci[j][0] = -fcd*(1-(1-(-epspci[j][0]/epsc2))**nc)
                            elif epspci[j][0]<=-epsc2 and epspci[j][0]>=-epscu2:
                                sigpci[j][0] = -fcd
                            else:
                                sigpci[j][0] = 0
                        else:
                            sigpci[j][0] = 0
                elif concModel == 2:
                    for j in range(2):
                        if (epspci[j][0])<0:
                            if epspci[j][0]>-epsc3:
                                sigpci[j][0] = fcd*epspci[j][0]/epsc3
                            elif epspci[j][0]<=-epsc3 and epspci[j][0]>=-epscu3:
                                sigpci[j][0] = -fcd
                            else:
                                sigpci[j][0] = 0
                        else:
                            sigpci[j][0] = 0
                        
                
                #Collection of principal stresses in one Matrix 3xn
                sigpciMat[:,i:i+1] = sigpci
                
                # Conversion of principal stresses to global xy directions
                TepsciTransp= np.transpose(Tepsci)
                sigci = np.matmul(TepsciTransp, sigpci)
                sigciMat[:,i:i+1] = sigci
            
            epspciMax = np.amin(epspciMat)
            sigpciMax = np.amin(sigpciMat)
            if concModel == 1:
                utilratioc = abs(epspciMax/epscu2)      #Utilization ratio for concrete
            else:
                utilratioc = abs(epspciMax/epscu3)
            
            
            #Concrete contribution to internal force vector
            
            Sc = np.zeros((6,1))        #main iternal concrete vector
            
            for i in range(n):
                sc = np.zeros((6,1))    #internal concrete vector for each layer
                sc[0][0] = (h/n) * sigciMat[0][i]
                sc[1][0] = (h/n) * sigciMat[1][i]
                sc[2][0] = (h/n) * sigciMat[2][i]
                sc[3][0] = (h/n) * (-zc[i]) * sigciMat[0][i] 
                sc[4][0] = (h/n) * (-zc[i]) * sigciMat[1][i]
                sc[5][0] = (h/n) * (-zc[i]) * sigciMat[2][i]
                
                Sc = Sc + sc
                
            
            #Strains and stresses in reinforcement
            
            sigpsj = np.zeros((3,1))
            sigpsjMat = np.zeros((3,2)) #Stress Matrix in principal direction
            sigsjMat = np.zeros((3,2))  #Stress Matrix in global direction
            epspsjMat = np.zeros((3,2)) #Strain Matrix in principal direction
            TepssjMat = np.zeros((3,6)) #Matrix to contain Tranformation Matrices for alfa
            
            for j in range(2):          #j represents upper and lower reinforcement
                if zs[j] == 0:
                    Asj = np.zeros((3,6))
                else:
                    Asj = np.array([[1, 0, 0, -zs[j], 0, 0], [0, 1, 0, 0, -zs[j], 0], [0, 0, 1, 0, 0, -zs[j]]])
                    
                epssj = np.matmul(Asj,epst)     #3x1
                
                alfasj =  alfa[j]*np.pi / 180   #convert from degrees to radians
                cs = np.cos(alfasj)
                ss = np.sin(alfasj)
                Tepssj = np.array([[cs**2, ss**2, ss*cs], [ss**2, cs**2, -ss*cs], [-2*ss*cs, 2*ss*cs, cs**2-ss**2]])
                TepssjMat[:,3*j:3*j+3] = Tepssj     
                
                epspsj = np.matmul(Tepssj, epssj)   #3x1
                epspsjMat[:,j:j+1] = epspsj
                
                
                #Stresses in principal directions
                
                if abs(epspsj[0][0]) <= fyd/Esx[j]:
                    sigpsj[0][0] =  (epspsj[0][0]) * Esx[j]
                elif abs(epspsj[0][0]) <= epsud:
                    sigpsj[0][0] = np.sign(epspsj[0][0]) * fyd
                else:
                    sigpsj[0][0] = 0
                    
                if abs(epspsj[1][0]) <= fyd/Esy[j]:
                    sigpsj[1][0] =  (epspsj[1][0]) * Esy[j]
                elif abs(epspsj[1][0]) <= epsud:
                    sigpsj[1][0] = np.sign(epspsj[1][0]) * fyd
                else: 
                      sigpsj[1][0] = 0
                
                # no need to specify sigpsj[2][0] since it is already 0
                
                TepssjTransp = np.transpose(Tepssj)
                sigsj = np.matmul(TepssjTransp, sigpsj)
                sigsjMat[:, j:j+1] = sigsj
                
                sigpsjMat[:,j:j+1] = sigpsj
                
                #Maximum and minimum values for reinforcement (unused in GUI results)
            epspsjMin = np.amin(epspsjMat)
            epspsjMax = np.amax(epspsjMat)
            sigpsjMin = np.amin(sigpsjMat)
            sigpsjMax = np.amax(sigpsjMat)

                
            #Reinforcement contribution to internal force vector
            
            Ss = np.zeros((6,1))
            
            for j in range (2):
                ss = np.zeros((6,1))
                ss[0][0] = Asx[j] * sigsjMat[0][j]
                ss[1][0] = Asy[j] * sigsjMat[1][j]
                ss[2][0] = 0
                ss[3][0] = -zs[j] *Asx[j] * sigsjMat[0][j]
                ss[4][0] = -zs[j] *Asy[j] * sigsjMat[1][j]
                ss[5][0] = 0
                
                Ss = Ss + ss
                
            #Internal Force Vector
            
            S = Sc + Ss
            
            #Addition of values of S to the internal forces list
            
            self.nxSList.append(S[0][0])
            self.nySList.append(S[1][0])
            self.nxySList.append(S[2][0])
            self.mxSList.append(S[3][0]/1000)
            self.mySList.append(S[4][0]/1000)
            self.mxySList.append(S[5][0]/1000)
            
            #Calculation of the maximum difference between S and R
            
            Diff = R - S
            RelDev = np.zeros((6,1))    #Matrix containing the relative deviation
            RelDiff = np.zeros((6,1))   #Matrix containing the differences    
            for d in range (6):
                if R[d][0] == 0:
                    RelDev[d][0] = 0
                    RelDiff[d][0] = abs(Diff[d][0])
                else:
                    RelDev[d][0] = abs(Diff[d][0]/R[d][0])
                    
            devMax = np.amax(RelDev)    #Maximum realtive deviation
            diffNMax = np.amax(RelDiff[0:3,0:1])  #Maximum difference for the values where values of N in R are 0
            diffMMax = np.amax(RelDiff[3:6,0:1])  #Maximum difference for the values where values of M in R are 0
            
            
            # New Young's Modulus and Material Matrix for Concrete
            
            CcgiMat = np.zeros((3,3*n)) #Matrix containing all Material Matrices (for every concrete layer)
            for i in range(n):
                if epspciMat[0][i] == 0:
                    E11 = 0
                else:
                    E11 = (sigpciMat[0][i])/(epspciMat[0][i])
                    
                if epspciMat[1][i] == 0:
                    E22 = 0
                else:
                    E22 = (sigpciMat[1][i])/(epspciMat[1][i])
                    
                E12 = (E11+E22)/2
                
                #Material Matrix in the principal direction
                Ccpi = (1/(1-vc**2)) * np.array([[E11, vc*E12, 0], [vc*E12, E22, 0], [0, 0, (1-vc)*E12/2]])
            
                #Material Matrix in the global direction (Ccgi)
                Tepsci = TepsciMat[:,3*i:3*i+3]    
                TepsciTransp = np.transpose(Tepsci)
                Ccgi = np.matmul((np.matmul(TepsciTransp,Ccpi)),Tepsci)
                CcgiMat[:,3*i:3*i+3] = Ccgi     #recent change
                
            
            # New Young's Modulus and Material Matrix for the Reinforcement
            
            CsgxMat = np.zeros((3,6))   #Matrix collecting Ex
            CsgyMat = np.zeros((3,6))   #Matrix collecting Ey
            for j in range(2):
                
                if epspsjMat[0][j] ==0:
                    EsxNew = 0
                else:
                    EsxNew = sigpsjMat[0][j]/epspsjMat[0][j]
                
                if epspsjMat[1][j] ==0:
                    EsyNew = 0
                else:
                    EsyNew = sigpsjMat[1][j]/epspsjMat[1][j]
                    
                    #Young's Modulus principal direction x and y   
                Cspx = np.array([[EsxNew, 0, 0], [0, 0, 0], [0, 0, 0]])
                Cspy = np.array([[0, 0, 0], [0, EsyNew, 0], [0, 0, 0]])
                
                    #Young's Modulus global direction x and y (to be used for New Stiffness Matrix)
                Tepssj = TepssjMat[:,3*j:3*j+3]     
                TepssjTransp = np.transpose(Tepssj)
                Csgx = np.matmul((np.matmul(TepssjTransp,Cspx)),Tepssj)
                Csgy = np.matmul((np.matmul(TepssjTransp,Cspy)),Tepssj)
                
                CsgxMat[:,3*j:3*j+3]= Csgx
                CsgyMat[:,3*j:3*j+3]= Csgy
                
                
            #New stiffness Matrix for concrete
            
            KcNew = np.zeros((6,6))
            for i in range(n):
                KcNewTemp = transformTo6x6(CcgiMat[:,3*i:3*i+3], -zc[i]*CcgiMat[:,3*i:3*i+3], -zc[i]*CcgiMat[:,3*i:3*i+3], ((zc[i])**2)*CcgiMat[:,3*i:3*i+3])  #temporary KcNew
                KcNew = KcNew + KcNewTemp
            KcNew = (h/n)*KcNew  
            
                
            
            #New stiffness Matrix for the reinforcement
            
            KsxNew = np.zeros((6,6))
            KsyNew = np.zeros((6,6))
                
            for j in range(2):
                Csx = CsgxMat[:,3*j:3*j+3]
                Csy = CsgyMat[:,3*j:3*j+3]
                KsxTemp = transformTo6x6(Csx, -zs[j]*Csx, -zs[j]*Csx, (zs[j]**2)*Csx)
                Ksx = Asx[j] * KsxTemp
                KsyTemp = transformTo6x6(Csy, -zs[j]*Csy, -zs[j]*Csy, (zs[j]**2)*Csy)
                Ksy = Asy[j] * KsyTemp
                 
                KsxNew = KsxNew + Ksx
                KsyNew = KsyNew + Ksy
        
            loop = loop + 1
            start = 1
            if loop >= maxIt:
                self.convergence = False
                break
            
            if np.all((K==0)):
                self.convergence = False
                break
        
            #Convergence result
        convergence = self.convergence
        
    
            #Internal forces convert moment values to  M= kNm
        S[3][0] = S[3][0] /1000
        S[4][0] = S[4][0] /1000
        S[5][0] = S[5][0] /1000
        
        #Utilization Ratio for steel
        reinfUtilRatioAsx1 = abs( epspsjMat[0][0]/epsud)
        reinfUtilRatioAsx2 = abs( epspsjMat[0][1]/epsud)
        reinfUtilRatioAsy1 = abs (epspsjMat[1][0]/epsud)
        reinfUtilRatioAsy2 = abs (epspsjMat[1][1]/epsud)    
        
        
            #Results appropriately rounded
            
                 #epsciMax (max strain) is represente as per mille
        concResult = [(epspciMax*1000).round(3), sigpciMax.round(2), utilratioc.round(2)]
                #reinfStrain (reinforcemnt strain) is represented as per mille
        reinfStrain = (epspsjMat.flatten()*1000).round(3)
        reinfStress = (sigpsjMat.flatten()).round(2)
        reinfUtilRatio =  [reinfUtilRatioAsx1.round(2),reinfUtilRatioAsx2.round(2), reinfUtilRatioAsy1.round(2),reinfUtilRatioAsy2.round(2)]
        intForce = S.round(2).flatten()
        RelDiff = RelDiff.flatten()
        RelDev = RelDev.flatten()
        loop
                    
        return concResult, reinfStrain, reinfStress, reinfUtilRatio, intForce, RelDiff, RelDev, loop, epspciMat, sigpciMat, sigpsjMat, convergence  
        
    
    
    
    
    #Plotting 
    def plot(self, loop, epspciMat, sigpciMat, sigpsjMat ): 
    
            #Concrete
                #Strain (displayed in graph as per mille)
        
        strainpc1 = ((1000 * epspciMat[0,:]).round(3)).flatten() 
        strainpc2 = ((1000 * epspciMat[1,:]).round(3)).flatten()
        
                #Principal strains for a coesive strain graphs
        strainpc1New = np.zeros(self.__n)           
        strainpc2New = np.zeros(self.__n)

        for i in range (self.__n):
            if strainpc1[i]< strainpc2[i]:
                strainpc1New[i]=strainpc1[i]
                strainpc2New[i]=strainpc2[i]
            else:
                strainpc1New[i]=strainpc2[i]
                strainpc2New[i]=strainpc1[i]
                
        
        layers = list(range(1,self.__n+1))
        
        straincFig1 = Figure(figsize=(7,7), constrained_layout=True)      #Figure 1: concrete strain 1
        straincFig1.suptitle('Concrete principal strain (\u03B51) distribution', fontsize=10)
        chart1 = straincFig1.add_subplot()
        #chart1.set_title('Concrete principal strain (\u03B51) distribution', fontsize=10)
        chart1.set_xlabel('Strain \u2030 (\u03B51)', fontsize=8)
        chart1.set_ylabel('Layers', fontsize=8)
        chart1.tick_params(axis = 'x', labelsize=8)
        chart1.tick_params(axis = 'y', labelsize=8)
        chart1.set_ylim([1,self.__n])
        chart1.plot(strainpc1New, layers)
        
        straincFig2 = Figure(figsize=(7,7), constrained_layout=True)      #Figure 2 : concrete strain 2
        straincFig2.suptitle('Concrete principal strain (\u03B52) distribution', fontsize=10)
        chart2 = straincFig2.add_subplot()
        #chart2.set_title('Concrete principal strain (\u03B52) distribution', fontsize=10)
        chart2.set_xlabel('Strain \u2030 (\u03B52)', fontsize=8)
        chart2.set_ylabel('Layers', fontsize=8)
        chart2.tick_params(axis = 'x', labelsize=8)
        chart2.tick_params(axis = 'y', labelsize=8)
        chart2.set_ylim([1,self.__n])
        chart2.plot(strainpc2New, layers)
        
                #Stress
        
        stresspc1 = (sigpciMat[0,:]).flatten()
        stresspc2 = (sigpciMat[1,:]).flatten()
        
            #Principal stresses for a coesive graph
        stresspc1New = np.zeros(self.__n)           
        stresspc2New = np.zeros(self.__n)

        for i in range (self.__n):
            if stresspc1[i]< stresspc2[i]:
                stresspc1New[i]=stresspc1[i]
                stresspc2New[i]=stresspc2[i]
            else:
                stresspc1New[i]=stresspc2[i]
                stresspc2New[i]=stresspc1[i]
                
        
        stresscFig1 = Figure(figsize=(7,7), constrained_layout=True)      #Figure 3: concrete stress 1
        stresscFig1.suptitle('Concrete principal stress (\u03C31) distribution', fontsize=10)
        chart3 = stresscFig1.add_subplot()
        #chart3.set_title('Concrete principal stress (\u03C31) distribution', fontsize=10)
        chart3.set_xlabel('Stress (\u03C31)', fontsize=8)
        chart3.set_ylabel('Layers', fontsize=8)
        chart3.tick_params(axis = 'x', labelsize=8)
        chart3.tick_params(axis = 'y', labelsize=8)
        chart3.set_ylim([1,self.__n])
        chart3.plot(stresspc1New, layers)
        
        stresscFig2 = Figure(figsize=(7,7), constrained_layout=True)        #Figure 4: concrete stress 2
        stresscFig2.suptitle('Concrete principal stress (\u03C32) distribution', fontsize=10)
        chart4 = stresscFig2.add_subplot()
        #chart4.set_title('Concrete principal stress (\u03C32) distribution', fontsize=10)
        chart4.set_xlabel('Stress (\u03C32)', fontsize=8)
        chart4.set_ylabel('Layers', fontsize=8)
        chart4.tick_params(axis = 'x', labelsize=8)
        chart4.tick_params(axis = 'y', labelsize=8)
        chart4.set_ylim([1,self.__n])
        chart4.plot(stresspc2New, layers)
        
            
        
            #Reinforcement stress  (Not part of the final return of the plot method)
        
        stressx = (sigpsjMat[0,:]).flatten()
        stressy = (sigpsjMat[1,:]).flatten()
        reinfH1 = self.__c1                            #Reinforcement Height layer 1
        reinfH2 = self.__h - self.__c2                 #Reinforcement Height layer 2
        
        
        
                #x-direction
        
        stressFig1 = Figure(figsize=(7,7), constrained_layout=True)   #Figure 5: reinforcement stress x
        chart5 = stressFig1.add_subplot()
        chart5.set_title('Steel stress (\u03C3s) distribution in x-direction', fontsize=10)
        chart5.set_xlabel('Stress (\u03C3s)', fontsize=8)
        chart5.set_ylabel('Thickness (mm)', fontsize=8)
        chart5.tick_params(axis = 'x', labelsize=8)
        chart5.tick_params(axis = 'y', labelsize=8)
        chart5.set_ylim([0,self.__h+1]) 
        chart5.set_xlim([-self.fyd,self.fyd])
        chart5.hlines(reinfH1, 0, stressx[0])
        chart5.hlines(reinfH2, 0, stressx[1])
        
                #y-direction
        
        stressFig2 = Figure(figsize=(7,7), constrained_layout=True)   #Figure 6: reinforcement stress y
        chart6 = stressFig2.add_subplot()
        chart6.set_title('Steel stress (\u03C3s) distribution in y-direction', fontsize=10)
        chart6.set_xlabel('Stress (\u03C3s)', fontsize=8)
        chart6.set_ylabel('Thickness (mm)', fontsize=8)
        chart6.tick_params(axis = 'x', labelsize=8)
        chart6.tick_params(axis = 'y', labelsize=8)
        chart6.set_ylim([0,self.__h+1])
        chart6.set_xlim([-self.fyd, self.fyd])
        chart6.hlines(reinfH1, 0, stressy[0])
        chart6.hlines(reinfH2, 0, stressy[1])
        
        
            #Convergence
            
                #R values in kN and kNm
        R = np.array([[self.R[0][0]], [self.R[1][0]], [self.R[2][0]], [self.R[3][0]/1000], [self.R[4][0]/1000], [self.R[5][0]/1000]])
            
                #Forces
        
                    #Nx
        
        forceFig1 = Figure(figsize=(7,7), constrained_layout=True)     #Figure 7: Nx
        forceFig1.suptitle('Nx Convergence', fontsize=10)
        nxPlot = forceFig1.add_subplot()
        #nxPlot.set_title('Nx Convergence', fontsize=10)
        nxPlot.set_xlabel('Iteration', fontsize=8)
        nxPlot.set_ylabel('Nx (kN)', fontsize=8)
        nxPlot.tick_params(axis = 'x', labelsize=8)
        nxPlot.tick_params(axis = 'y', labelsize=8)
        nxPlot.hlines(R[0][0], 0, loop)
        nxPlot.plot(self.loopList, self.nxSList, 'r')
        
                    #Ny
        
        forceFig2 = Figure(figsize=(7,7), constrained_layout=True)     #Figure 8: Ny
        forceFig2.suptitle('Ny Convergence', fontsize=10)
        nyPlot = forceFig2.add_subplot()
        #nyPlot.set_title('Ny Convergence', fontsize=10)
        nyPlot.set_xlabel('Iteration', fontsize=8)
        nyPlot.set_ylabel('Ny (kN)', fontsize=8)
        nyPlot.tick_params(axis = 'x', labelsize=8)
        nyPlot.tick_params(axis = 'y', labelsize=8)
        nyPlot.hlines(R[1][0], 0, loop)
        nyPlot.plot(self.loopList, self.nySList, 'r')
        
                    #Nxy
        
        forceFig3 = Figure(figsize=(7,7), constrained_layout=True)     #Figure 9: Nxy
        forceFig3.suptitle('Nxy Convergence', fontsize=10)
        nxyPlot = forceFig3.add_subplot()
        #nxyPlot.set_title('Nxy Convergence', fontsize=10)
        nxyPlot.set_xlabel('Iteration', fontsize=8)
        nxyPlot.set_ylabel('Nxy (kN)', fontsize=8)
        nxyPlot.tick_params(axis = 'x', labelsize=8)
        nxyPlot.tick_params(axis = 'y', labelsize=8)
        nxyPlot.hlines(R[2][0], 0, loop)
        nxyPlot.plot(self.loopList, self.nxySList, 'r')
        
        
                #Moments
                
        
        
                      #Mx
        
        momentFig1 = Figure(figsize=(7,7), constrained_layout=True)     #Figure 10: Mx
        momentFig1.suptitle('Mx Convergence', fontsize=10)
        mxPlot = momentFig1.add_subplot()
        #mxPlot.set_title('Mx Convergence', fontsize=10)
        mxPlot.set_xlabel('Iteration', fontsize=8)
        mxPlot.set_ylabel('Mx (kNm)', fontsize=8)
        mxPlot.tick_params(axis = 'x', labelsize=8)
        mxPlot.tick_params(axis = 'y', labelsize=8)
        mxPlot.hlines(R[3][0], 0, loop)
        mxPlot.plot(self.loopList, self.mxSList, 'r')
        
                    #My
        
        momentFig2 = Figure(figsize=(7,7), constrained_layout=True)     #Figure 11: My
        momentFig2.suptitle('My Convergence', fontsize=10)
        myPlot = momentFig2.add_subplot()
        #myPlot.set_title('My Convergence', fontsize=10)
        myPlot.set_xlabel('Iteration', fontsize=8)
        myPlot.set_ylabel('My (kNm)', fontsize=8)
        myPlot.tick_params(axis = 'x', labelsize=8)
        myPlot.tick_params(axis = 'y', labelsize=8)
        myPlot.hlines(R[4][0], 0, loop)
        myPlot.plot(self.loopList, self.mySList, 'r')
                    #Mxy
        
        momentFig3 = Figure(figsize=(7,7), constrained_layout=True)     #Figure 12: Mxy
        momentFig3.suptitle('Mxy Convergence', fontsize=10)
        mxyPlot = momentFig3.add_subplot()
        #mxyPlot.set_title('Mxy Convergence', fontsize=10)
        mxyPlot.set_xlabel('Iteration', fontsize=8)
        mxyPlot.set_ylabel('Mxy (kNm)', fontsize=8)
        mxyPlot.tick_params(axis = 'x', labelsize=8)
        mxyPlot.tick_params(axis = 'y', labelsize=8)
        mxyPlot.hlines(R[5][0], 0, loop)
        mxyPlot.plot(self.loopList, self.mxySList, 'r')
       
        
        return straincFig1, straincFig2, stresscFig1, stresscFig2, forceFig1, forceFig2, forceFig3, momentFig1, momentFig2, momentFig3
            

