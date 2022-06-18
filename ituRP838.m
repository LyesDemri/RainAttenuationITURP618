function [k,alpha] = ituRP838(polarisation,elevation,latitude,f)

ajkH=[-5.33980,-0.35351,-0.23789,-0.94158];
bjkH=[-0.10008, 1.26970, 0.86036, 0.64552];
cjkH=[ 1.13098, 0.45400, 0.15354, 0.16817];
mkkH=-0.18961;
ckkH= 0.71147;

ajkV=[-3.80595,-3.44965,-0.39902, 0.50167];
bjkV=[ 0.56934,-0.22911, 0.73042, 1.07319];
cjkV=[ 0.81061, 0.51059, 0.11899, 0.27195];
mkkV=-0.16398;
ckkV= 0.63297;

ajalphaH=[-0.14318,0.29591,0.32177,-5.37610, 16.1721];
bjalphaH=[ 1.82442,0.77564,0.63773,-0.96230,-3.29980];
cjalphaH=[-0.55187,0.19822,0.13164, 1.47828, 3.43990];
malphaH=0.67849;
calphaH=-1.95537;

ajalphaV=[-0.07771,0.56727,-0.20238,-48.2991  ,48.5833];
bjalphaV=[ 2.33840,0.95545, 1.14520,  0.791669, 0.791459];
cjalphaV=[-0.76284,0.54039, 0.26809,  0.116226, 0.116479];
malphaV=-0.053739;
calphaV=0.83433;

thetadeg = elevation;
phi=latitude;
taudeg = thetadeg + phi;
taurad = taudeg * 2 * pi / 360;


if strcmp(polarisation,'H')
    k=10.^(sum(ajkH.*exp(-(((log10(f)-bjkH)./cjkH).^2)))+mkkH*log10(f)+ckkH);
    alpha=(sum(ajalphaH.*exp(-(((log10(f)-bjalphaH)./cjalphaH).^2)))+malphaH*log10(f)+calphaH);
elseif strcmp(polarisation,'V')
    k=10.^(sum(ajkV.*exp(-(((log10(f)-bjkV)./cjkV).^2)))+mkkV*log10(f)+ckkV);
    alpha=(sum(ajalphaV.*exp(-(((log10(f)-bjalphaV)./cjalphaV).^2)))+malphaV*log10(f)+calphaV);
elseif strcmp(polarisation,'linear')
    kH=10.^(sum(ajkH.*exp(-(((log10(f)-bjkH)./cjkH).^2)))+mkkH*log10(f)+ckkH);
    alphaH=(sum(ajalphaH.*exp(-(((log10(f)-bjalphaH)./cjalphaH).^2)))+malphaH*log10(f)+calphaH);
    kV=10.^(sum(ajkV.*exp(-(((log10(f)-bjkV)./cjkV).^2)))+mkkV*log10(f)+ckkV);
    alphaV=(sum(ajalphaV.*exp(-(((log10(f)-bjalphaV)./cjalphaV).^2)))+malphaV*log10(f)+calphaV);
    k=(kH+kV+(kH-kV)*((cos(thetarad))^2)*cos(2*taurad))/2;
    alpha=(kH*alphaH+kV*alphaV+(kH*alphaH-kV*alphaV)*((cos(thetarad))^2)*cos(2*taurad))/(2*k);
elseif strcmp(polarisation,'circular')
    taudeg=45;
    taurad=taudeg*2*pi/360;
    kH=10.^(sum(ajkH.*exp(-(((log10(f)-bjkH)./cjkH).^2)))+mkkH*log10(f)+ckkH);
    alphaH=(sum(ajalphaH.*exp(-(((log10(f)-bjalphaH)./cjalphaH).^2)))+malphaH*log10(f)+calphaH);
    kV=10.^(sum(ajkV.*exp(-(((log10(f)-bjkV)./cjkV).^2)))+mkkV*log10(f)+ckkV);
    alphaV=(sum(ajalphaV.*exp(-(((log10(f)-bjalphaV)./cjalphaV).^2)))+malphaV*log10(f)+calphaV);
    k=(kH+kV+(kH-kV)*((cos(thetarad))^2)*cos(2*taurad))/2;
    alpha=(kH*alphaH+kV*alphaV+(kH*alphaH-kV*alphaV)*((cos(thetarad))^2)*cos(2*taurad))/(2*k);
else
    error('Invalid value for polarisation')
end

