from frame import *
from gauss import gaussElimination

def uklad(data : GlobalData):
    iter = int(data.time/data.step)
    temp = [data.inittT for _ in range(data.nN)] 
    frame = Frame(data)
    hg = frame.calculateHGlobal()
    cg = frame.calculateCGlobal()
    pg = frame.calculatePlGlobal()
    for iteracja in range(iter):
        
        matrix = [[0 for _ in range(len(hg) + 1)] for _ in range(len(hg))]
        for i in range(len(hg)):
            for j in range(len(hg[i])):
                matrix[i][j] = hg[i][j] + cg[i][j]/data.step
        
    

        for i in range(len(pg)):
            cxt = 0 
            for j in range(len(cg[i])):
                cxt += (cg[i][j]/data.step) * temp[j]
            matrix[i][len(hg)] = -pg[i] + cxt
    
        tempNew = gaussElimination(matrix)
        #print(iteracja + 1,":\tmin: {}, max: {}".format(min(tempNew),max(tempNew)))
        print(iteracja+1," ",min(tempNew)," ",max(tempNew))
        temp = tempNew
            



if __name__ == '__main__':
    data = GlobalData("data.txt")
    frame = Frame(data)
    #frame.printElements()

    # for elem in frame.elements:
    #     for h in elem.h:
    #         print(h)
    #     print()

    # for elem in frame.elements:
    #     for c in elem.c:
    #         print(c)
    #     print()

    # for h in frame.calculateHGlobal():
    #     print(h)

    # print()
    # print(frame.calculatePlGlobal())

    # for elem in frame.elements:
    #     print(elem.pl)
    # print()

    # e = Element([Node(0,0,1), Node(4,0,2), Node(4,6,3), Node(0,6,4)],0)
    # el4 = Elem4(2,e,30)

    # for h in el4.calculateH():
    #     print(h)

    uklad(data)

