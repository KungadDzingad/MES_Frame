import math


PARAM = 25

#-----------------------------------------------------------------------#


class GlobalData:
    def __init__(self, fileName: str):
        try:
            with open(fileName, 'r') as file:
                self.H = float(file.readline())
                self.W = float(file.readline())
                self.nH = int(file.readline())
                self.nW = int(file.readline())
                self.param = float(file.readline())
                self.pc = int(file.readline())
                self.c = int(file.readline())
                self.ro = int(file.readline())
                self.alfa = int(file.readline())
                self.t = int(file.readline())
                self.inittT = int(file.readline())
                self.time = int(file.readline())
                self.step = int(file.readline())
        except FileNotFoundError:
            raise

        self.nE = (self.nH - 1) * (self.nW - 1)
        self.nN = self.nH * self.nW
        self.deltaX = self.W / (self.nW - 1)
        self.deltaY = self.H / (self.nH - 1)


#------------------------------------------------------------------------#


class Node:

    def __init__(self, x: float, y: float, idf: int, bc):
        self.x = x
        self.y = y
        self.id = idf
        self.bc = bc  #etykieta ustalająca czy węzeł jest na brzegu

    def __str__(self):
        return "Id: {}, x: {}, y: {}".format(self.id, self.x, self.y) 


#-----------------------------------------------------------------------#


class Element:
    def __init__(self, nodes, nr,alfa,t):
        if len(nodes) != 4:
            print("Wrong Node number")
            raise Exception("Wrong number of Nodes")

        self.nodes = nodes
        self.nr = nr
        self.alfa = alfa
        self.t = t

    
    def printElem(self):
        print("element nr: {}".format(self.nr))
        for node in self.nodes:
            print("\t", node)



    def calculateHCLocal(self, param, nr,c,ro):
        e4 =  Elem4(nr, self, param,c,ro,self.alfa,self.t)
        self.h = e4.h   
        self.c = e4.calculateC()
        self.pl = e4.pl
        

#--------------------------------------------------------------------------#


class Frame:
    def __init__(self, data: GlobalData):
        self.data = data

        self.nodes = []
        for i in range(data.nW):
            for j in range(data.nH):
                x = i * data.deltaX
                y = j * data.deltaY
                bc = 0
                if (x == data.W or x == 0 ) or (y==data.H or y==0):
                    bc = 1
                self.nodes.append(Node(x,y, i * data.nH + j,bc))

        self.elements = []
        for i in range(data.nW -1 ):
            for j in range(data.nH -1):
                id1 = j + i* data.nH
                id2 = id1 + data.nH
                id3 = id2 + 1
                id4 = id1 + 1
                nodesOfElement = [self.nodes[id1], self.nodes[id2], self.nodes[id3], self.nodes[id4]]
                e = Element(nodesOfElement, (i * (data.nH-1) + j),data.alfa,data.t)
                e.calculateHCLocal(data.param,data.pc,self.data.c,self.data.ro)
                self.elements.append(e)


    def printElements(self):
        for e in self.elements:
            e.printElem() 


    def calculateHGlobal(self):
        h = [[0 for _ in range(len(self.nodes))] for _ in range(len(self.nodes))]

        for elem in self.elements:
            for i in range(len(elem.nodes)):
                for j in range(len(elem.nodes)):
                    h[elem.nodes[i].id][elem.nodes[j].id] += elem.h[i][j]
        return h


    def calculateCGlobal(self):
        c = [[0 for _ in range(len(self.nodes))] for _ in range(len(self.nodes))]

        for elem in self.elements:
            for i in range(len(elem.nodes)):
                for j in range(len(elem.nodes)):
                    c[elem.nodes[i].id][elem.nodes[j].id] += elem.c[i][j]
        return c

    
    def calculatePlGlobal(self):
        pl = [0 for _ in range(len(self.nodes))]
        for elem in self.elements:
            for i in range(len(elem.nodes)):
                pl[elem.nodes[i].id] += elem.pl[i]


        return pl

        
#---------------------------------------------------------------------------------------#


class DNCalculator:
    def __init__(self, nr):
        self.nr = nr

        self.fN = [
            lambda eta,ksi: 0.25*(1-ksi)*(1-eta),
            lambda eta,ksi: 0.25*(1+ksi)*(1-eta),
            lambda eta,ksi: 0.25*(1+ksi)*(1+eta),
            lambda eta,ksi: 0.25*(1-ksi)*(1+eta)    
        ]

        self.fKsi = [
            lambda eta: -0.25 * (1 - eta),
            lambda eta: 0.25 * (1 - eta),
            lambda eta: 0.25 * (1 + eta),
            lambda eta: -0.25 * (1 + eta)
        ]

        self.fEta = [
            lambda ksi: -0.25 * (1 - ksi),
            lambda ksi: -0.25 * (1 + ksi),
            lambda ksi: 0.25 * (1 + ksi),
            lambda ksi: 0.25 * (1 - ksi)
        ]


        if nr == 2:
            a = 1 / math.sqrt(3)
            self.ksi = [-a, a, a, -a]
            self.eta = [-a, -a, a, a]
            self.wKsi = [1,1,1,1]
            self.wEta = [1,1,1,1]

        elif nr == 3:
            self.ksi = [-(math.sqrt(3/5)), 0, (math.sqrt(3/5)), -(math.sqrt(3/5)), 0, (math.sqrt(3/5)), -(math.sqrt(3/5)), 0, (math.sqrt(3/5))]
            self.eta = [-(math.sqrt(3/5)), -(math.sqrt(3/5)), -(math.sqrt(3/5)), 0, 0, 0, (math.sqrt(3/5)), (math.sqrt(3/5)), (math.sqrt(3/5))]
            self.wKsi = [5 / 9, 8 / 9, 5 / 9, 5 / 9, 8 / 9, 5 / 9, 5 / 9, 8 / 9, 5 / 9]
            self.wEta = [5 / 9, 5 / 9, 5 / 9, 8 / 9, 8 / 9, 8 / 9, 5 / 9, 5 / 9, 5 / 9]

        elif nr == 4:
            x = math.sqrt((3/7)-(2/7)*math.sqrt(6/5))
            wx = (18 + math.sqrt(30))/36
            y = math.sqrt((3/7) + (2/7)*math.sqrt(6/5))
            wy = (18 - math.sqrt(30))/36

            self.ksi = [-y,-x,x,y,-y,-x,x,y,-y,-x,x,y,-y,-x,x,y]
            self.eta = [-y,-y,-y,-y,-x,-x,-x,-x,x,x,x,x,y,y,y,y]
            self.wKsi = [wy,wx,wx,wy,wy,wx,wx,wy,wy,wx,wx,wy,wy,wx,wx,wy]
            self.wEta = [wy,wy,wy,wy,wx,wx,wx,wx,wx,wx,wx,wx,wy,wy,wy,wy]


    def calculatedNdKsi(self):
        Ns = []
        for i in range(self.nr**2):
            N = []
            for j in range(4):
                N.append(self.fKsi[j](self.eta[i]))
            Ns.append(N)

        return Ns


    def calculatedNdEta(self):
        Ns = []
        for i in range(self.nr **2):
            N = []
            for j in range(4):
                N.append(self.fEta[j](self.ksi[i]))
            Ns.append(N)

        return Ns
     

    def calculateN(self):
        N = []
        for i in range(self.nr**2):
            n = []
            for j in range(4):
                n.append(self.fN[j](self.eta[i],self.ksi[i]))
            N.append(n)
        return N


    def calulatePcBc(self,i,j):
        h = []
        H = [[0 for _ in range(4)] for _ in range(4)]

        if i == 0 and j == 1:
            for m in range(self.nr):
                s=[]
                for n in self.fN:
                    s.append(n(-1,self.ksi[m])) #s.append(n(-1,self.ksi[i]))
                h.append(s)

        elif i == 1 and j == 2:
            for m in range(self.nr):
                s=[]
                for n in self.fN:
                    s.append(n(self.ksi[m],1))
                h.append(s)
            
        elif i==2 and j ==3:
            for m in range(self.nr):
                s=[]
                for n in self.fN:
                    s.append(n(1,self.ksi[m]))
                h.append(s)

        elif i ==3 and j ==0:
            for m in range(self.nr):
                s=[]
                for n in self.fN:
                    s.append(n(self.ksi[m],-1))
                h.append(s)

        h2 = []
        for i in range(len(h)):
            h22 = []
            for j in range(4):
                h222 = []
                for k in range(4):
                    h222.append(h[i][j] * h[i][k])
                h22.append(h222)
            h2.append(h22)

       
        for i in range(len(h)):
            for j in range(4):
                for k in range(4):
                    H[j][k] += h2[i][j][k] * self.wKsi[i]

        return H


    def calculatePLL(self,i,j,t):
        pl = []

        if i == 0 and j == 1:
            for m in range(self.nr):
                s=[]
                for n in self.fN:
                    s.append(n(-1,self.ksi[m])* self.wKsi[m]) 
                pl.append(s)

        elif i == 1 and j == 2:
            for m in range(self.nr):
                s=[]
                for n in self.fN:
                    s.append(n(self.ksi[m],1)* self.wKsi[m])
                pl.append(s)
            
        elif i==2 and j ==3:
            for m in range(self.nr):
                s=[]
                for n in self.fN:
                    s.append(n(1,self.ksi[m]) * self.wKsi[m])
                pl.append(s)

        elif i ==3 and j ==0:
            for m in range(self.nr):
                s=[]
                for n in self.fN:
                    s.append(n(self.ksi[m],-1)* self.wKsi[m])
                pl.append(s)
        
        return pl


#----------------------------------------------------------------------------------#


def det2x2(matrix):
    return matrix[0][0] * matrix[1][1] - (matrix[0][1] * matrix[1][0])


#----------------------------------------------------------------------------------#


class Elem4:
    def __init__(self, nr, element:Element, parameter, c, ro,alfa,t):
        self.calculator = DNCalculator(nr)
        self.element = element
        self.c = c
        self.ro = ro
        self.parameter = parameter
        self.alfa = alfa
        self.t = t
        self.dNdKsi = self.calculator.calculatedNdKsi()
        self.dNdEta = self.calculator.calculatedNdEta()
        self.N = self.calculator.calculateN()
        self.det =  []
        self.jacobians = self.calculateJacobians()
        self.pl = [0,0,0,0]
        for j in self.jacobians:
            self.det.append(det2x2(j))
        self.h = self.calculateH()
        self.hbcAndPl()

    def calculateJacobians(self): 
        jac = []
        for i in range(self.calculator.nr**2):
            xKsi, yEta, xEta, yKsi = 0,0,0,0
            for j in range(4):
                xKsi += self.dNdKsi[i][j] * self.element.nodes[j].x
                yEta += self.dNdEta[i][j] * self.element.nodes[j].y
                xEta += self.dNdEta[i][j] * self.element.nodes[j].x
                yKsi += self.dNdKsi[i][j] * self.element.nodes[j].y
    
            jac.append([[xKsi,xEta] , [yKsi, yEta]])
        return jac


    def calculateN(self):
        
        dNdX, dNdY = [], []
        for i,jac in enumerate(self.jacobians):
            nx, ny = [], []
            for k in range(4):
                x = (1/self.det[i])*(self.dNdKsi[i][k] * jac[1][1] + self.dNdEta[i][k] * (-jac[0][1]))
                y = (1/self.det[i])*(self.dNdKsi[i][k] * (-jac[1][0]) + self.dNdEta[i][k] * jac[0][0])
                nx.append(x)
                ny.append(y)
            dNdX.append(nx)
            dNdY.append(ny)
        
        return dNdX, dNdY


    def calculateHPoint(self):
        dNdX, dNdY = self.calculateN()

        H = []
        for i in range(len(dNdX)):
            h=[]
            for j in range(len(dNdX[i])):
                s =[]
                for k in range(len(dNdX[i])):
                    x = dNdX[i][j] * dNdX[i][k] * self.det[i] * self.parameter
                    y = dNdY[i][j] * dNdY[i][k] * self.det[i] * self.parameter
                    s.append((x+y))
                h.append(s)
            H.append(h)
        return H


    def calculateH(self):
        hL = self.calculateHPoint()
    
        hs = [[0 for _ in range(len(self.element.nodes))] for _ in range(len(self.element.nodes))]

        for i in range(self.calculator.nr**2):
            for j, hss in enumerate(hs):
                for k, h in enumerate(hss):    
                    hs[j][k] +=  (hL[i][j][k] * self.calculator.wKsi[i] * self.calculator.wEta[i])
        return hs

    
    def calculateCPoint(self):
        CPoints = []

        for i in range(len(self.N)):
            cp = []
            for j in range(len(self.N[i])):
                c = []
                for k in range(len(self.N[i])):
                    
                    c.append(self.c  * self.ro * self.N[i][j] * self.N[i][k] * self.det[i])
                cp.append(c)
            CPoints.append(cp)
        return CPoints


    def calculateC(self):
        cL = self.calculateCPoint()

        C = [[0 for _ in range(len(self.element.nodes))] for _ in range(len(self.element.nodes))]

        for i in range(self.calculator.nr**2):
            for j in range(len(C)):
                for k in range(len(C[j])):
                    C[j][k] += cL[i][j][k] * self.calculator.wKsi[i] * self.calculator.wEta[i]    
        
        return C


    def hbcAndPl(self):
        if self.element.nodes[0].bc == 1 and self.element.nodes[1].bc == 1:
            l = math.fabs(self.element.nodes[1].x - self.element.nodes[0].x)/2 
            self.addHbcToH(0,1,l)
            self.addToPl(0,1,l)

        if self.element.nodes[1].bc == 1 and self.element.nodes[2].bc == 1:
            l = math.fabs(self.element.nodes[2].y - self.element.nodes[1].y)/2
            self.addHbcToH(1,2,l)
            self.addToPl(1,2,l)

        if self.element.nodes[2].bc == 1 and self.element.nodes[3].bc == 1:
            l = math.fabs(self.element.nodes[2].x - self.element.nodes[3].x)/2
            self.addHbcToH(2,3,l)
            self.addToPl(2,3,l)

        if self.element.nodes[3].bc == 1 and self.element.nodes[0].bc == 1:
            l = math.fabs(self.element.nodes[0].y - self.element.nodes[3].y)/2
            self.addHbcToH(3,0,l)
            self.addToPl(3,0,l)


    def addHbcToH(self,i,j,l):
        hbc = self.calculator.calulatePcBc(i,j)
        for m in range(4):
            for n in range(4):
                self.h[m][n] += hbc[m][n] * l * self.alfa
                
    
    def addToPl(self,i,j,l):
        pll = self.calculator.calculatePLL(i,j,self.t)
        for m in range(self.calculator.nr):
            for n in range(4):
                self.pl[n] -= pll[m][n] * self.alfa * l * self.t


#---------------------------------------------- ---------------------------------------------------------#
