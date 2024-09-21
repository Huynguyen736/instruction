
#init var
a = 5.0 #cm
Cx = 90 #cm
Cy = 90 #cm

#Static
N = 1290.08 #P (Tons)
Mx = 8.36 #M2 (T.m)
My = 63.14 #M3 (T.m)

#Static dict
EBT = {
    "B35":2100
    }

Es = {
    "400V": 190000
}

RnBT = {
    "B35":15
    }

RsT = {
    "400V":2300
}

gammab3 = {
    "full":1,
    "notfull":0.85
}

gammab5 = {
    "full":1,
    "notfull":0.85
}

Station = [457] #column 4

#calc var
l0 = max(Station)
lambX = l0/(0.288*Cx)
lambY = l0/(0.288*Cy)
lamb = max(lambX, lambY)
Eax = max(l0/0.7/600, Cx/30)
Eay = max(l0/0.7/600, Cy/30)
Jx = Cy*(Cx**3)/12 # cm4
Jy = Cx*(Cy**3)/12 # cm4
Nthx = 2.5*EBT["B35"]*10*Jx/l0**2 #KG
Nthy = 2.5*EBT["B35"]*10*Jy/l0**2 #KG

if lambX < 28:  muyX = 1
else:   muyX = 1/(1-(N*1000/Nthx))

if lambY < 28:  muyY = 1
else:   muyY = 1/(1-(N*1000/Nthy))

Mx1 = muyX*Mx
My1 = muyY*My

if Mx1/Cx < My1/Cy: h = Cy
else: h = Cx

if My1/Cy < Mx1/Cx: b = Cx
else: b = Cy

if Mx1/Cx < My1/Cy: M1 = My1
else: M1 = Mx1

if My1/Cy < Mx1/Cx: M2 = Mx1
else: M2 = My1

x1 = N*1000/(RnBT["B35"]*10*b)

if x1 > (h-a):  m0 = 0.4
else: m0 = 1-0.6*x1/(h-a)

if Mx1/Cx < My1/Cy: Ea = Eay+0.2*Eax
else: Ea = Eax+0.2*Eay

e1 = (M1+m0*M2*h/b)/N*100
e0 = max(Ea, e1)

exl = e0/(h-a)

if exl>0.3:
    xicmae = False
else: xicmae = 1/((0.5-exl)*(2+exl))

if lamb > 14: phi = 1.028 - 0.0000288*lamb**2-0.0016*lamb
else: phi = 1

if exl > 0.3: phie = False
else: phie = phi+(1-phi)*exl/0.3

if exl > 0.3:
    AstLTRB = False
else:   AstLTRB = ((lamb*N*10**4/phie-RnBT["B35"]*gammab3["full"]*gammab5["full"]*h*b*100)/(RsT["400V"]-RnBT["B35"]*gammab3["full"]*gammab5["full"]))/100

if exl > 0.3:
    if x1 > Es["400V"]*(h-a): x = Es["400V"]/(1+50*(e0/h)**2)*(h-a)
    else: x = False
else: x = False

if exl > 0.3:
    if x1 > Es["400V"]*(h-a): AstLTB = (N*10**4*(exl+h/2-a)*10-RnBT["B35"])*gammab3["full"]*gammab5["full"]*b*x*(h-a-x/2)*10**3/(0.4*RsT["400V"]*(h-2*a)*10/100)
    AstLTB = False
else:   AstLTB = False

if exl > 0.3:
    if x1 > Es["400V"]*(h-a): AstL = False
    else: AstL = N*10**4*((e0+h/2-a)*10+x1/2*10-(h-a)*10)/(0.4*RsT["400V"]*(h-2*a)*10)/100
else: AstL = False

Asct = 0.01*b*(h-a)/100

# calc res

#show direction
if Mx1/(Cx/100) < My1/(Cy/100): print("Phương tính toán: Y") 
else: print("Phương tính toán: X")

#show lech tam
if (xicmae != False or phi != False or phie != False):  
    print("Lệch tâm rất bé")
    As = AstLTRB
elif (x != False): 
    print("Lệch tâm bé")
    As = AstLTB
else: 
    print("Lệch tâm lớn")
    As = AstL

#show As
print(f"As: {As}")

#show Asct
print(f"Asct: {Asct}")