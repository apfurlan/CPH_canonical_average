#!/usr/bin/env python

import matplotlib,os
matplotlib.use('TkAgg') 
from pylab import *
from math import *
from numpy import *

rc('text', usetex=True)
rc('font', **{'family': 'serif','serif':['Computer Modern']})
rc('xtick', labelsize=22)
rc('ytick', labelsize=22)

subplots_adjust(bottom=0.11,left=0.13,right=0.965,top=0.97,
                wspace=0.32,hspace=0.)



L=int(raw_input("L ="))
Npart=int((L*L)*(3./4.))

Nmus=13

output_file='avgcan_L_'+str(L)+'_N_'+str(Npart)+'.dat'
ftmd=open(output_file,'w')
ftmd.write("# count    temp \t \t <E>  \t  \t  Cv \t \t phi \t \t xhi(phi) \t \tU4 \t \t <mu> \n")
ftmd.write("#"+"-"*135+"\n")


fTemp=6.8 ; iTemp=0.10 ; delTemp=0.001
Ntemps=int((fTemp-iTemp)/delTemp)

filename='alg_DOS_L_'+str(L)+'_N_'+str(Npart)+'.dat'
#filename='alg_ran_0_DOS_L_'+str(L)+'_N_'+str(Npart)+'.dat'

if os.path.exists(filename): 

    data=loadtxt(filename)
    En  =data[:,1]
    lngE=data[:,2]
    H   =data[:,3]
    op  =data[:,4]
    op2 =data[:,5]
    op3 =data[:,6]
    op4 =data[:,7]
#    mu  =data[:,8]

    Cv    =zeros((Ntemps), dtype=float64)
    avEn  =zeros((Ntemps), dtype=float64)
    avEn2 =zeros((Ntemps), dtype=float64)
    avop  =zeros((Ntemps), dtype=float64)
    avop2 =zeros((Ntemps), dtype=float64)
    avop4 =zeros((Ntemps), dtype=float64)
    xhiop =zeros((Ntemps), dtype=float64)
    binU4 =zeros((Ntemps), dtype=float64)
#    avMu  =zeros((Ntemps), dtype=float64)
    
#    filename2='alg_pbmu_L_'+str(L)+'_N_'+str(Npart)+'.dat'
#    if os.path.exists(filename2):

#        Enbins=185
#        Nmus=13

#        pdmuE=zeros((Enbins,Nmus),float64)
        
#        data2=loadtxt(filename2)
#        for i in range(Enbins):
#            for j in range(Nmus):
#                pdmuE[i][j]=data2[i*Nmus+j,3]
    
    
    for j in range(Ntemps):
        temp=iTemp+j*delTemp
        beta=1./temp
        beta2=beta*beta

        Zi=0.
        A=zeros(len(En), dtype=float64)
        P=zeros(len(En), dtype=float64)
        M=zeros(len(En), dtype=float64)
        
        for i in range(len(En)):
            A[i]=(lngE[i]-beta*En[i])
            
 #           if os.path.exists(filename2):
 #               for k in range(Nmus):
 #                   M[i]=M[i]+pdmuE[i][k]*exp(-beta*(Nmus-6))
            

        maxA=max(A)
        Zi=sum(exp(A))

        for i in range(len(En)):
            P[i]=exp(A[i]-maxA)
            
        Zi=sum(P)
        for i in range(len(En)):
            avEn[j] =avEn[j] + En[i]*P[i]   
            avEn2[j]=avEn2[j]+ (En[i]*En[i])*P[i]
            avop[j] =avop[j] + op[i] *P[i]
            avop2[j]=avop2[j]+ op2[i]*P[i]
            avop4[j]=avop4[j]+ op4[i]*P[i]
#            avMu[j]=avMu[j] + M[i]*P[i]
            

        avEn[j]=avEn[j]/Zi
        avEn2[j]=avEn2[j]/Zi
        Cv[j]=beta2*(avEn2[j]-avEn[j]*avEn[j])

        avop[j]= avop[j] /Zi
        avop2[j]=avop2[j]/Zi
        avop4[j]=avop4[j]/Zi

        xhiop[j]=beta*(L*L)*(avop2[j] - avop[j]*avop[j] )
        binU4[j]=1.-(avop4[j]/(3*avop2[j]*avop2[j]))
        
#        avMu[j]=avMu[j]/Zi
        #avMu[j]=avMu[j]

        ftmd.writelines([str('%4i' %j )       ,'       ', 
                         str('%.5f' %temp )   ,'       ', 
                         str('%.5e' %avEn[j]) ,'       ',
                         str('%.5e' %Cv[j])   ,'       ',
                         str('%.5e' %avop[j]) ,'       ',
                         str('%.5e' %xhiop[j]),'       ',
                         str('%.5e' %binU4[j]),'\n'])#'       ',
#                         str('%.5e' %avMu[j]),'\n'])
