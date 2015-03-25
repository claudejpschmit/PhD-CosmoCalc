from math import *
def bjl(l, x):

    if (l < 0):
        print "Can not evaluate Spherical Bessel Function with index l<0"
        return -1
    ln2=0.6931471805599453094
    onemln2=0.30685281944005469058277
    pid2=1.5707963267948966192313217
    pid4=0.78539816339744830961566084582
    rootpi12 = 21.269446210866192327578
    gamma1 =   2.6789385347077476336556 #!/* Gamma function of 1/3 */
    gamma2 =   1.3541179394264004169452 #!/* Gamma function of 2/3 */
    pi=3.141592653589793238463

    ax = abs(x)
    ax2 = ax**2
    if (l<7):
        if (l==0):
            if (ax < 0.1):
                jl = 1.0 - ax2/6.0 * (1.0 - ax2/20.0)
            else:
                jl = sin(ax)/ax

        elif (l == 1):
            if (ax < 0.2):
                jl = ax/3.0*(1.0 - ax2/10.0 * (1.0 - ax2/28.0))
            else:
                jl = (sin(ax)/ax - cos(ax))/ax

        elif (l == 2):
            if (ax < 0.3):
                jl = ax2/15.0 * (1.0 - ax2/14.0 * (1.0-ax2/36.0))
            else:
                jl = (-3.0 * cos(ax)/ax - sin(ax) * (1.0 - 3.0/ax2))/ax
        elif (l == 3):
            if (ax < 0.4):
                jl = ax*ax2/105.0 * (1.0 - ax2/18.0*(1.0 - ax2/44.0))
            else:
                jl = (cos(ax)*(1.0-15.0/ax2)-sin(ax) * (6.0-15.0/ax2)/ax)/ax
        elif (l == 4):
            if (ax < 0.6):
                jl = ax2**2/945.0 * (1.0-ax2/22.0 * (1.0 - ax2/52.0))
            else:
                jl = (sin(ax)*(1.0-(45.0-105.0/ax2)/ax2)+cos(ax)*(10.0-105.0/ax2)/ax)/ax
        elif (l == 5):
            if (ax < 1.0):
                jl = ax2**2 * 2 * ax/10395.0*(1.0 - ax2/26.0 * (1.0 - ax2/60.0))
            else:
                jl = (sin(ax) * (15.0 - (420.0 - 945.0/ax2)/ax2)/ax - cos(ax)*(1.0 - (105.0-945.0/ax2)/ax2))/ax
        else:
            if (ax < 1.0):
                jl = ax2**3/135135.0 * (1.0 - ax2/30.0*(1.0-ax2/68.0))
            else:
                jl = (sin(ax) * (-1.0 + (210.0 - (4725.0 - 10395.0/ax2)/ax2)/ax2)+ cos(ax) * (-21.0 + (1260.0-10395.0/ax2)/ax2)/ax)/ax
    else:
        nu = l + 0.5
        nu2 = nu**2
        if (ax < 1e-40):
            jl = 0.0
        elif ((ax2/l)<0.5):
            jl = exp(l * log(ax/nu) - ln2 + nu * onemln2 - (1.0 - (1.0 - 3.5/nu2)/nu2/30.0)/12.0/nu)/nu * (1.0 - ax2/(4.0*nu+4.0)*(1.0-ax2/(8.0*nu + 16.0)*(1.0-ax2/(12.0*nu + 36.0))))
        elif ((double(l)**2/ax)<0.5):
            beta = ax - pid2*(l+1)
            jl = (cos(beta) * (1.0-(nu2 - 0.25)*(nu2-2.25)/8.0/ax2*(1.0-(nu2-6.25)*(nu2-12.25)/48.0/ax2)) - sin(beta)*(nu2-0.25)/2.0/ax*(1.0-(nu2-2.25)*(nu2-6.25)/24.0/ax2*(1.0-(nu2-12.25)*(nu2-20.25)/80.0/ax2)))/ax
        else:
            l3=nu**0.325
            if (ax < (nu -1.31*l3)):
                cosb = nu/ax
                sx = sqrt(nu2-ax2)
                cotb = nu/sx
                secb = ax/nu
                beta = log(cosb+sx/ax)
                cot3b = cotb**3
                cot6b = cot3b**2
                sec2b = secb**2
                expterm=( (2.0+3.0*sec2b)*cot3b/24.0 - ( (4.0+sec2b)*sec2b*cot6b/16.0  + ((16.0-(1512.0+(3654.0+375.0*sec2b)*sec2b)*sec2b)*cot3b/5760.0 + (32.0+(288.0+(232.0+13.0*sec2b)*sec2b)*sec2b)*sec2b*cot6b/128.0/nu)*cot6b/nu)/NU)/NU

                jl = sqrt(cotb*cosb)/(2.0*nu)*exp(-nu*beta+nu/cotb-expterm)

#          /**************** Region 2: x >> l ****************/

            elif (ax > (nu + 1.48 * l3)):
                COSB=nu/ax
                SX=sqrt(ax2-nu2)
                COTB=nu/sx
                SECB=ax/nu
                BETA=acos(COSB)
                COT3B=COTB**3
                COT6B=COT3B**2
                SEC2B=SECB**2
                TRIGARG=nu/COTB-nu*BETA-pid4-((2.0+3.0*SEC2B)*COT3B/24.0+(16.0-(1512.0+(3654.0+375.0*SEC2B)*SEC2B)*SEC2B)*COT3B*COT6B/5760.0/nu2)/nu
                EXPTERM=( (4.0+SEC2B)*SEC2B*COT6B/16.0-(32.0+(288.0+(232.0+13.0*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B**2/128.0/nu2)/nu2

                jl=sqrt(COTB*COSB)/nu*exp(-EXPTERM)*cos(TRIGARG)

#          /***************** Region 3: x near l ****************/

            else:

                BETA=ax-nu
                BETA2=BETA**2
                SX=6.0/ax
                SX2=SX**2
                SECB=SX**0.3333333333333333
                SEC2B=SECB**2
                jl=( gamma1*SECB + BETA*gamma2*SEC2B -(BETA2/18.0-1.0/45.0)*BETA*SX*SECB*gamma1 -((BETA2-1.0)*BETA2/36.0+1.0/420.0)*SX*SEC2B*gamma2 +(((BETA2/1620.0-7.0/3240.0)*BETA2+1.0/648.0)*BETA2-1.0/8100.0)*SX2*SECB*gamma1 +(((BETA2/4536.0-1.0/810.0)*BETA2+19.0/11340.0)*BETA2-13.0/28350.0)*BETA*SX2*SEC2B*gamma2 -((((BETA2/349920.0-1.0/29160.0)*BETA2+71.0/583200.0)*BETA2-121.0/874800.0)* BETA2+7939.0/224532000.0)*BETA*SX2*SX*SECB*gamma1)*sqrt(SX)/rootpi12



    if ((x < 0) and (l%2 != 0)):
        jl=-jl

    return jl
