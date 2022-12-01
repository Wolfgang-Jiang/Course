import math
import numpy as np
import matplotlib.pyplot as plt
import sympy


#basic parameters of concrete and steel
b = 150     #mm
h = 280     #mm
tas=34.2    #mm top steel a
h0 =251.4
bas=h-h0    #mm base steel a
As1= 339  #mm2
As= 339  #mm2
fc = 37.9   #N/mm2
ft= 4      #N/mm2
Ec = 32500  #N/mm2
Es = 200000 #N/mm2
fy = 521    #N/mm2

Nmax = fc*b*h+fy*(As1+As)
"""N here!!!!!!!!!"""
N=1900000   #N
print(Nmax)
#calculation1
"""no crack state calculation:"""
aE = Es/Ec                                      #const
Aso = b*h +(aE-1)*As1+(aE-1)*As            #mm2
ch = (0.5*h*h*b+(aE-1)*(As1*tas+As*h0))/Aso   #centroid   #mm
Ix = (b*h**3/12 + b*h*(h/2-ch)**2 +(aE-1)*(As1*(tas-ch)**2+As*(h0-ch)**2) ) #mm^4
#M crack
Mcr = Ix/(h0-ch)*(ft-N*((h0-ch)*(ch-h/2)/Ix-1/(Aso)))   #Nmm #This is Mcr as calculated from 1

#phi1, Mphi1 are two linspace
#calculation 2
#Max N?

#The following calculate MAXIMUM moment:
x = sympy.symbols('x')
epsstl = 0.0033*(x-h0)/x
Fn= 0.798*fc*b*x-0.667*fc*b*x*(1-(h/x))**2*(3.3-2.732*(1-h/x))
res = sympy.solve([Fn+fy*As1+epsstl*Es*As-N],[x])
for i in range(3):
    result = str(res[i])[1:-2]
    try:
        r = float(result)
    except:
        pass
#r is now the resulting x
print(res)
print(r)
yc = (0.318-0.222*(1-h/r)**2*(1+2*h/r)*(3.3-2.723*(1-h/r))) / (0.798-0.667*(1-(h/r))**2*(3.3-2.732*(1-h/r))) *r
Mu = (0.798*fc*b*r-0.667*fc*b*r*(1-(h/r))**2*(3.3-2.732*(1-h/r))) *(h0-yc) +fy*As1*(h0-tas)-N*(h0-h/2)
print(yc)
print(Mu/1000000)
#Nu = Fn +fy*As1 + epsstl*Es*As

if Mcr<Mu:
    phicr = Mcr/Ec/Ix
    phi_1 = np.linspace(0,phicr,100)
    Mphi = Ec*Ix*phi_1/1000000
    #plt.plot(phi_1,Mphi,color="red",linewidth=2)


    """phi_2 = np.linspace(phicr,2*phicr,100)
    Mcr2 =np.linspace(Mcr/1000000,Mcr/1000000,100)
    plt.plot(phi_2,Mcr2,color="red",linewidth="2")
    """
    #cracked section
    sect = 40
    epstop = np.linspace(0.0001,0.0033,sect)
    eps0 = 0.002
    r = epstop/eps0
    alpha = np.zeros((sect))
    beta = np.zeros((sect))
    for i in range(sect):
        if r[i] < 1:
            alpha[i] = 1 - (1 - (1 - r[i]) ** 3) / (3 * r[i])
            beta[i] = (1 / 2 - (1 + (-1 - 3 * r[i]) * (1 - r[i]) ** 3) / (12 * r[i] ** 2)) / alpha[i]
        else:
            alpha[i] = 1 - 1 / (3 * r[i])
            beta[i] = (1 / 2 - 1 / 12 * r[i] ** 2) / alpha[i]

    x = np.ones(sect)
    M = np.zeros(sect)
    Phi = np.zeros(sect)
    for i in range(sect):
        xi = sympy.symbols("x") #xi is the variable x, if viable, then return to x list.
        #first case:
        def first_case(i,xi):
            res0 = sympy.solve([xi**2*b*fc*alpha[i]+(Es*epstop[i]*(As+As1)-N)*xi - (Es*epstop[i]*(As*h0+As1*tas))],[xi])
            sol1 = float(str(res0[0])[1:-2])
            sol2 = float(str(res0[1])[1:-2])
            if 0<sol1<h:
                sol = sol1
            else:
                sol=sol2
            if 0 < sol < h: #set global variables?
                strain_steel_com = epstop[i]*(sol-tas)/sol
                strain_steel_ten = epstop[i]*(h0-sol)/sol
                if strain_steel_ten < fy/Es and strain_steel_com < fy/Es:   #数量级！
                    #print(strain_steel_com/(fy/Es),strain_steel_ten/(fy/Es))
                    Fs1 = Es*As1*epstop[i]*(sol-tas)/sol
                    Fs = Es*As*epstop[i]*(h0-sol)/sol
                    M[i] = (h/2-tas)*Fs1+(h0-h/2)*Fs +(h/2-sol*(1-beta[i]))*alpha[i]*fc*b*sol
                    print("1 case: " + " Fs: "+str(Fs) + " Fs1: "+ str(Fs1) +"x value:" + str(sol) + " M = " + str(M[i]))

                    return sol
                else:
                    return second_case(i,xi)
            else:
                return second_case(i,xi)
                #function returns the solved x for scenario 1.
        def second_case(i,x):
            res = sympy.solve([x**2*b*fc*alpha[i]+(Es*epstop[i]*As1-fy*As-N)*x -Es*epstop[i]*As1*tas],[x])
            sol1 = float(str(res[0])[1:-2])
            sol2 = float(str(res[1])[1:-2])
            if 0<sol1<h:
                sol=sol1
            elif 0<sol2<h:
                sol=sol2
            else:
                return third_case(i,x)
            strain_steel_com = epstop[i] * (sol - tas) / sol
            strain_steel_ten = epstop[i] * (h0 - sol) / sol
            if strain_steel_ten > fy/Es and strain_steel_com < fy/Es:
                Fs = fy * As
                Fs1 = Es * As1 * epstop[i] * (sol - tas) / sol
                M[i] = (h/2-tas)*Fs1+(h0-h/2)*Fs +(h/2-sol*(1-beta[i]))*alpha[i]*fc*b*sol
                print("1 case: " + " Fs: "+ str(Fs) + " Fs1: " + str(Fs1) + "x value:" + str(sol) + " M = " + str(M[i]))
                #print("case 2:"+str(epstop[i]/sol))
                return sol
            else:
                return third_case(i,x)

        def third_case(i,x):
            res = sympy.solve([x**2*b*fc*alpha[i]+(Es*epstop[i]*As+fy*As1-N)*x -Es*epstop[i]*As*h0],[x])
            sol1 = float(str(res[0])[1:-2])
            sol2 = float(str(res[1])[1:-2])
            if 0<sol1<h:
                sol=sol1
            elif 0<sol2<h:
                sol=sol2
            else:
                return fourth_case(i,x)
            strain_steel_com = epstop[i] * (sol - tas) / sol
            strain_steel_ten = epstop[i] * (h0 - sol) / sol
            #print(strain_steel_ten /( fy / Es))
            if strain_steel_ten < fy / Es and strain_steel_com > fy / Es:
                Fs = Es * As * epstop[i] * (h0 - sol) / sol
                Fs1 = fy * As1
                M[i] = (h / 2 - tas) * Fs1 + (h0 - h / 2) * Fs + (h / 2 - sol * (1 - beta[i])) * alpha[i] * fc * b * sol
                print("3 case: "+str(M[i])+" X = "+ str(sol))

                return sol
            else:
                return fourth_case(i, x)
        def fourth_case(i,x):
            sol = (N+fy*As-fy*As1)/(alpha[i]*b*fc)
            if 0<sol<h:
                strain_steel_com = epstop[i] * (sol - tas) / sol
                strain_steel_ten = epstop[i] * (h0 - sol) / sol
                if strain_steel_ten > fy / Es and strain_steel_com > fy / Es:
                    Fs = fy * As
                    Fs1 = As1 * fy
                    M[i] = (h/2-tas)*Fs1+(h0-h/2)*Fs +(h/2-sol*(1-beta[i]))*alpha[i]*fc*b*sol
                    print("4 case: " + str(M[i])+ " X = "+str(sol))
                    return sol
                else:
                    return 0
            else:
                return 0

        x[i] = first_case(i,xi)
        if x[i]!=0:
            Phi[i] = epstop[i]/x[i]
        else:Phi[i]=phicr
    #recalibration
    realM = np.zeros(sect)
    boo = False
    for i in range(sect):
        if M[i] <=Mcr:
            if boo==False:
                realM[i] = max(M[i]/1000000,Mcr/1000000)
            else:
                realM[i] = M[i]/1000000
        else:
            realM[i] = M[i]/1000000
            boo = True

    print(realM)

    realPhi = np.zeros(sect)
    for i in range(sect):
        realPhi[i] =max(Phi[i],phicr)
    print(h)
    #concatenation:
    for i in range(sect):
        if realPhi[i] > phicr:
            h=realPhi[i]
            break
    overphi = np.linspace(phicr,h,20)
    overM = np.linspace(Mcr/1000000,Mcr/1000000, 20)
    realPhi = np.concatenate([phi_1,overphi,realPhi],axis=0 )
    realM = np.concatenate([Mphi,overM, realM],axis=0)
    plt.plot(realPhi,realM,color="red",linewidth="2")

    print(phicr)
    #plt.plot(Phi,M,color="blue")


    #plt.plot(Phi,M/1000000,color="blue")

else:
    phiu = Mu / Ec / Ix
    realPhi = np.linspace(0, phiu, 100)
    realM = Ec * Ix * realPhi / 1000000
    plt.plot(realPhi, realM, color="red", linewidth=2)
plt.title("Curvature-Moment relationship")
plt.xlabel("curvature")
plt.ylabel("Moment/kN·m")
plt.text(max(realPhi)*3/7,max(realM)*1/3,"Normal Force N: {0}  kN\nMaximal Moment: {1:.2f}  kN.m\nMaximum Curvature = {1:.2f} e-6".format(N/1000,max(realM),max(realPhi)))
plt.show()

