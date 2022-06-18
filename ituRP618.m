clear;clc;close all;

longitude = 3.01;
latitude = 36.77;
hs=0.184;
thetadeg=38.35;
thetarad=(thetadeg*2*pi)/360;
phi=latitude;
f=17300000000;
Re=8500;
%étape1
hR=ituRP839(longitude,latitude);
%étape2
Ls=(hR-hs)/sin(thetarad);
%étape3
LG=Ls*cos(thetarad);
%étape4
R001=42;
%étape5
[k,alpha] = ituRP838('H',thetadeg,phi,f/(1e9));
gammaR=k*(R001)^(alpha);
%étape6
r001=1/(1+0.78*sqrt((LG*gammaR)/(f/1e9))-0.38*(1-exp(-2*LG)));
%étape7
zetarad=atan((hR-hs)/(LG*r001));
zetadeg = zetarad*180/pi;
if zetadeg>thetadeg 
    LR=(LG*r001)/cos(thetarad);
else
    LR=(hR-hs)/sin(thetarad);
end
if  abs(phi)<36
    chi=36-abs(phi);
else
    chi=0;
end
nu001=1/(1+sqrt(sin(thetarad))*(31*(1-exp(-(thetadeg/(1+chi))))*((sqrt(LR*gammaR))/((f/1e9)^2))-0.45));
%étape8
LE=LR*nu001;
%étape9
A001=gammaR*LE
%etape 10:
p=0.01;
if p>=1 || abs(phi)>=36
        beta=0;
    elseif p<1 && abs(phi)<36 && thetadeg>=25
        beta=-0.005*(abs(phi)-36);
    else 
        beta=-0.005*(abs(phi)-36)+1.8-4.25*sin(thetarad);
end
Ap=A001*(p/0.01)^(-(0.655+0.033*log(p)-0.045*log(A001)-beta*(1-p)*sin(thetarad)))

%%
%étape10
index=1;
for p=0.001:0.001:5
    if p>=1 || abs(phi)>=36
        beta=0;
    elseif p<1 && abs(phi)<36 && thetadeg>=25
        beta=-0.005*(abs(phi)-36);
    else 
        beta=-0.005*(abs(phi)-36)+1.8-4.25*sin(thetarad);
    end
    Ap=A001*(p/0.01)^(-(0.655+0.033*log(p)-0.045*log(A001)-beta*(1-p)*sin(thetarad)));
    ys(index)=Ap;
    xs(index)=p;
    index=index+1;
end

valeurs_satmaster = [99.999,26.75;
                     99.98 , 9.83;
                     99.9  ,4.49;
                     99.75 ,2.75;
                     99.5  ,1.73;
                     99    ,1.09;
                     98    ,0.66;
                     97    ,0.49;
                     96    ,0.39;
                     95    ,0.33];
valeurs_satmaster(:,1) = 100-valeurs_satmaster(:,1); 
                    
plot(xs,ys,'-o',valeurs_satmaster(:,1),valeurs_satmaster(:,2),'-r*')
