clear;
clc;
close all;
% Param�tres d'entr�e n�cessaires au calcul de l'affaiblissement
hs=0.068; 
thetadeg=38.35; 
thetarad=(thetadeg*2*pi)/360; 
phi=36.83; 
f=17.9; 
Re=8500; 
%�tape1 : calcul de la hauteur de la pluie en km
hR=ituRP839(3,36.83); 
%�tape2 : calcul de la longueur du trajet oblique en km
Ls=(hR-hs)/sin(thetarad); 
%�tape3 : calcul de la projection horizontale du trajet oblique en km
LG=Ls*cos(thetarad); 
%�tape4 : taux de pr�cipitation d�pass� pour 0,01% d'une ann�e moyenne en mm/h
R001=42; 
%�tape5 : calcul de l'affaiblissement lin�ique en dB/km 
[k,alpha]=ituRP838('H',thetadeg,phi,f); 
gammaR=k*(R001)^(alpha) 
%�tape6 : calcul du facteur de r�duction horizontale pour 0,01% du temps
r001=1/(1+0.78*sqrt((LG*gammaR)/f)-0.38*(1-exp(-2*LG))) 
%�tape7 : calcul du facteur d'ajustement vertical pour 0,01% du temps
zetarad=atan((hR-hs)/(LG*r001))  
zetadeg=zetarad*180/pi
if zetadeg>thetadeg 
    LR=(LG*r001)/cos(thetarad)
else
    LR=(hR-hs)/sin(thetarad)
end
if  abs(phi)<36
    chi=36-abs(phi)
else
    chi=0
end
nu001=1/(1+sqrt(sin(thetarad))*(31*(1-exp(-(thetadeg/(1+chi))))*((sqrt(LR*gammaR))/(f^2))-0.45)) 
%�tape8 : calcul de la longueur effective du trajet
LE=LR*nu001  
%�tape9 : calcul de l'affaiblissement pr�vu d�pass� pour 0,01% d'une ann�e moyenne
A001=gammaR*LE  
%�tape10 : calcul de l'affaiblissement pour d'autres pourcentages d'une ann�e moyenne
p=0.03 
if p>=0.01 || abs(phi)>=36
    beta=0
elseif p<0.01 && abs(phi)<36 && thetadeg>=25
    beta=-0.005*(abs(phi)-36)
else 
    beta=-0.005*(abs(phi)-36)+1.8-4.25*sin(thetarad)
end
Ap=A001*(p/0.01)^(-(0.655+0.033*log(p)-0.045*log(A001)-beta*(1-p)*sin(thetarad))) 
ApLin=10^(Ap/10);    
      


