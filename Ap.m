function msfuntmpl_basic(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2010 The MathWorks, Inc.

%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C-Mex counterpart: mdlInitializeSizes
%%
function setup(block)

% Register number of ports
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
% block.InputPort(1).Dimensions        = 3;
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Complex';
block.InputPort(1).DirectFeedthrough = true;

% Override output port properties
% block.OutputPort(1).Dimensions       = 3;
block.OutputPort(1).DatatypeID  = 0; % double
block.OutputPort(1).Complexity  = 'Complex';

% Register parameters
block.NumDialogPrms     = 8;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [-1 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
%% The MATLAB S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------

block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup

%%
%% PostPropagationSetup:
%%   Functionality    : Setup work areas and state variables. Can
%%                      also register run-time methods here
%%   Required         : No
%%   C-Mex counterpart: mdlSetWorkWidths
%%
function DoPostPropSetup(block)
block.NumDworks = 1;
  
  block.Dwork(1).Name            = 'x1';
  block.Dwork(1).Dimensions      = 2;
  block.Dwork(1).DatatypeID      = 0;      % double
  block.Dwork(1).Complexity      = 'Real'; % real
  block.Dwork(1).UsedAsDiscState = true;


%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is 
%%                      present in an enabled subsystem configured to reset 
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C-MEX counterpart: mdlInitializeConditions
%%
function InitializeConditions(block)

%end InitializeConditions


%%
%% Start:
%%   Functionality    : Called once at start of model execution. If you
%%                      have states that should be initialized once, this 
%%                      is the place to do it.
%%   Required         : No
%%   C-MEX counterpart: mdlStart
%%
function Start(block)

block.Dwork(1).Data = [0,0];

%end Start

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%
function Outputs(block)
hs=block.DialogPrm(1).Data;
thetadeg=block.DialogPrm(2).Data;
thetarad=(thetadeg*2*pi)/360;
phi=block.DialogPrm(3).Data;
f=block.DialogPrm(4).Data;
Re=block.DialogPrm(5).Data;
%étape1
% hR=block.DialogPrm(6).Data;
hR=ituRP839(3,36.83);
%étape2
Ls=(hR-hs)/sin(thetarad);
%étape3
LG=Ls*cos(thetarad);
%étape4
R001=block.DialogPrm(6).Data;
%étape5
polarisation=block.DialogPrm(7).Data;
[k,alpha]=ituRP838(polarisation,thetadeg,phi,f);
gammaR=k*(R001)^(alpha);
%étape6
r001=1/(1+0.78*sqrt((LG*gammaR)/f)-0.38*(1-exp(-2*LG)));
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
nu001=1/(1+sqrt(sin(thetarad))*(31*(1-exp(-(thetadeg/(1+chi))))*((sqrt(LR*gammaR))/(f^2))-0.45));
%étape8
LE=LR*nu001;
%étape9
A001=gammaR*LE
%étape10
p=block.DialogPrm(8).Data;
if p>=1 || abs(phi)>=36
    beta=0;
elseif p<1 && abs(phi)<36 && thetadeg>=25
    beta=-0.005*(abs(phi)-36);
else 
    beta=-0.005*(abs(phi)-36)+1.8-4.25*sin(thetarad);
end
Ap=A001*(p/0.01)^(-(0.655+0.033*log(p)-0.045*log(A001)-beta*(1-p)*sin(thetarad)));
ApLin=10^(Ap/10); % calcul de l'affaiblissement en liénaire
disp(['Ap=' num2str(Ap) 'dB,ApLin=' num2str(ApLin)])
block.OutputPort(1).Data = block.InputPort(1).Data/ApLin;


%end Outputs

%%
%% Update:
%%   Functionality    : Called to update discrete states
%%                      during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlUpdate
%%
function Update(block)

% block.Dwork(1).Data = block.InputPort(1).Data;

%end Update

%%
%% Derivatives:
%%   Functionality    : Called to update derivatives of
%%                      continuous states during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlDerivatives
%%
function Derivatives(block)

%end Derivatives

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C-MEX counterpart: mdlTerminate
%%
function Terminate(block)

%end Terminate

