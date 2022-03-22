clc;
clear all;
close all;
%% Final Project
%  Treston Brown
%%
w_rps = logspace(-3,3,1000);
R2D = (180/pi);
time = 30;
amp_g = 32.2 *20;
%==========================================================================>
%==========================================================================>
%% Plant Data
sel_state_alpha = 1;
sel_state_q     = 4;
sel_input_dele  = 2;
ol = [0;0];
% Coupled Dynamics
A = [-1.57        0      0     1     0; 
         0     -.50    .17     0    -1; 
    -21.13  -2876.7  -2.10  -.14  -.05;
    -82.92   -11.22   -.01  -.57     0;
      -.19   -11.86   -.01     0  -.57];
  
B = [     0    -.10       0;
       -.07       0     .11;
    -1234.7  -30.49 -1803.2;
      -4.82 -119.65      -7;
      14.84     .27 -150.58];
  
  eigs_ol = eig(A);
  disp('Open Loop EigenValues');
  disp(eigs_ol);

  %==========================================================================>
  %==========================================================================>
  %% Define Linearized Flight Dynamics
  
  %AeroDynamic Data
  Mach    = 2.5;
  c_fps   = 968.1;            %at 40000 feet
  V       = Mach*c_fps;
  Z_alpha = A(sel_state_alpha,sel_state_alpha)*V;
  Z_dele  = B(sel_state_alpha,sel_input_dele)*V;
  
  % Short Period Approximation

  % Short Period Dynamics
  Asp = [-1.57   1.00;
      -82.92 -0.57];
  Bsp = [-0.10;
         -119.65];
  Csp = [Z_alpha 0; 1 0; 0 1];
  Dsp = [Z_dele; 0; 0];

  sys_sp_plant = ss(Asp,Bsp,Csp,Dsp);
  
  Eigs_sp = eig(Asp);
  disp('Short Period EigenValues');
  disp(Eigs_sp);

  figure(1)
  plot(real(eig(A))   ,imag(eig(A)),'or','markersize',10); hold on
  plot(real(eig(Asp)) ,imag(eig(Asp)),'color', [0 1 .2],'linestyle','none','marker','x','markersize',10)
  set(gcf,'color','w');
  set(gca,'fontsize',14);
  xlabel('Re','fontsize',14);
  ylabel('Im','fontsize',14);
  title('Eigenvalues(OL Dynamics/Short period Dynamics)','fontsize',12);
  legend({'Open Loop \lambda','Short Period \lambda'},'Location','northeast','Orientation','vertical');
  axis([-2 0 -30 30])
  grid on
  hold off
  %==========================================================================>
  %==========================================================================>
  %% Define Actuator model
  wa_rps = 35*(2*pi);  %natural freq, rad/s
  za = .71;            %damping ratio, n/d
  
  Aact = [0 1; -wa_rps^2  -2*za*wa_rps]; nAact = size(Aact,1);
  Bact = [0 ; wa_rps^2];  mBact = size(Bact,2);
  Cact = [1 0];
  Dact = 0;
  sys_actuator =  ss(Aact,Bact,Cact,Dact,'StateName', {'dele','dele_dot'},'InputName',{'dele_cmd'},'OutputName',{'dele'}); %Actuator State Space
  
  %==========================================================================>
  %==========================================================================>
  %% Combined Short Period and Actuator Models (Analysis Model)
  sys_plant = sys_sp_plant * sys_actuator;
  Ap = sys_plant.a; nAsp_act = size(Ap,1);
  Bp = sys_plant.b;
  Cp = sys_plant.c;
  Dp = sys_plant.d;
  
  sel_sys_ActSP_output_Az     = 1;
  sel_sys_ActSP_output_alpha  = 2;
  sel_sys_ActSP_output_q      = 3;
  
  %==========================================================================>
  %==========================================================================>
  %% LQR Control Design
  %regulation matrices
  C_reg = [Z_alpha 0]; nC_reg = size(C_reg,1);
  D_reg = Z_dele;
  
  %RS system
  % Design model Matrices
  Aw = [0 C_reg; ol Asp];
  Bw = [Z_dele; Bsp]; [nBw,mBw] = size(Bw);
  
  % Design Matrices Eigenvalues
  eigs_AW = eig(Aw);
  disp('Design Matrices Eigen Values')
  disp(eigs_AW)
  
  % Determine if the system is controllable
  controllable = ctrb(Aw,Bw);
  control_rank = rank(controllable);
  [m,n] = size(controllable);
  if(control_rank == n)
      disp('System is controllable')
  else
      disp('System not controllable')
  end
  
  %LQR setup
  Q = zeros(3);
  R = 1;
  t = linspace(0, 2, 50); %time scale
  %Penalty Vectors
  q_eI  = logspace(-6,-3,100);
  N = length(q_eI);

for ii = 1 : length(q_eI)
    
    wq(ii) = ii; % Iteration Number
    Q(1,1) = q_eI(ii);
    Q(2,2) = 0;
    Q(3,3) = 0;
    % Use LQR system
    [Kw,~,~] = lqr(Aw,Bw,Q,R); % LQR solver
    
    %Dynamic Controller Matrices
    Ac  = [0];
    Bc1 = [1 0 0];
    Bc2 = [-1];
    Cc  = -Kw(1);
    Dc1 = [0 -Kw(2:3)];
    Dc2 = [0];
    sys_control =  ss(Ac,Bc1,Cc,Dc1);
      
    % Closed Loop System Dynamics
    I = eye(size(Dc1*Dp));
    Z = inv(I - Dc1*Dp); 
    Acl_act = [Ap+Bp*Z*Dc1*Cp Bp*Z*Cc; Bc1*(eye(3)+Dp*Z*Dc1)*Cp Ac+Bc1*Dp*Z*Cc]; 
    Bcl_act = [Bp*Z*Dc2; Bc2+Bc1*Dp*Z*Dc2];
    Ccl_act = [(eye(3)+Dp*Z*Dc1)*Cp Dp*Z*Cc]; 
    Dcl_act = [Dp*Z*Dc2]; 
    
    %Z       = eye(1) - Dc1 * Dsp;
    %Zinv    = Z\eye(1);
    %Acl_act = [Ac+Bc1*Dp*Zinv*Cc   Bc1*Cp+Bc1*Dp*Zinv*Dc1*Cp;
    %           Bp*Zinv*Cc              Ap+Bp*Zinv*Dc1*Cp];
    %Bcl_act = [Bc2; zeros(nAsp_act,1)];
    %Ccl_act = [0 Z_alpha 0 0 0 ] - Z_dele*[ Kw 0 0];
    %Dcl_act = 0;
    
    % Produce closed loop ss system
    sys_cl_act = ss(Acl_act,Bcl_act,Ccl_act,Dcl_act,'StateName',{'Az','alpha','q','dele','dele_dot'},'InputName',{'Az_cmd'},'OutputName',{'Az'});
    
    % Step Response of CL System
    opt = stepDataOptions('StepAmplitude', amp_g); % Add in 20g
    [ydesign,tdesign,xdesign] = step(sys_cl_act,opt);
    S = stepinfo(sys_cl_act);
         
    % Max fin Displacement and Max Fin Rate 
    dmax(ii)  =  max(abs(xdesign(:,3)))*(180/pi); % deg
    ddmax(ii) =  max(abs(xdesign(:,4)))*(180/pi); % 

    % Az Data
    OS(ii)           = S(1,1).Overshoot;
    US(ii)           = S(1,1).Undershoot;
    risetime(ii)     = S.RiseTime;
    settlingtime(ii) = S.SettlingTime;
    
    % Loop Gain
    A_Lu   = [Ap zeros(4,1); Bc1*Cp Ac];
    B_Lu   = [Bp; Bc1*Dp];
    C_Lu   = -[Dc1*Cp Cc]; % Change sign for loop gain
    D_Lu   = -[Dc1*Dp]; % Change signe for loop gain
    sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);
    

    % Margins and Gain Crossover Frequency
    S_Lu       = allmargin(sys_Lu);
    MGM        = max(S_Lu.GainMargin); % Maximum Gain Margin
    GM(ii)     = 20.*log10(MGM); %dB
    PM(ii)     = S_Lu.PhaseMargin;
    PM_rad     = deg2rad(PM);
    w_rpsp(ii) = S_Lu.PMFrequency;
    w_cg(ii)   = w_rpsp(ii)/2/pi; %Hz
  
    % Sensitivity and Co-sensitivity
    % Return Difference at input
    sImAinv = ((w_rpsp(ii))*eye(max(size(A_Lu)))-A_Lu)\eye(max(size(A_Lu)));
    % Compute loop gain at input
    Lu_eval = C_Lu*sImAinv*B_Lu + D_Lu;
    % Compute return difference
    RD = svd(eye(max(size(Lu_eval))) + Lu_eval);
    RDdesired(ii) = min(RD);
    
    %sensitivity
    sensitivity           = svd(1/(1+Lu_eval));
    sensitivity_dB(ii)    = abs(20*log10(sensitivity));
    cosensitivity         = svd(Lu_eval/(1+Lu_eval));
    cosensitivity_dB(ii)  = abs(20*log10(cosensitivity));
    
    StabilityRobustness   = svd(1+inv(Lu_eval));
    sig_minL(ii)          = RDdesired(ii);
    sig_minLinv(ii)       = min(StabilityRobustness);
     
    %% Properties
    [Wn,Z,P]         = damp(sys_cl_act);
    damping_min(ii)  = min(Z);
    Eigs_Cl(:,ii)    = P;    
    
    figure(2)
    hold on
    plot(tdesign,ydesign(:,1))
    title('Step Response of Az')
    xlabel('Time (s)')
    ylabel('Az (ft/s^2)')
    grid on

end
RD_min = min(sig_minL);
StabilityRobustness_min = min(sig_minLinv);
%==========================================================================>
%==========================================================================>
%% Design Point
q_eI_datapoint = 36;
% Design Data to Select Design Point
design_data_project = table(wq',q_eI',GM',PM',dmax',ddmax',w_cg',US',OS',sensitivity_dB',cosensitivity_dB');
design_data_project.Properties.VariableNames = {'Iteration';'q_eI (Penalty Value)';'Gain Margin';'Phase Margin';'fin displacement';'fin rate';'Crossover Frequency';'% undershoot';'% overshoot';'sensitivity dB';'cosensitivity dB'};
filename = 'Project_Data.xlsx';
writetable(design_data_project,filename);

Q(1,1) = q_eI(q_eI_datapoint);
Q(2,2) = 0;
Q(3,3) = 0;

% Use LQR system
Kw = lqr(Aw,Bw,Q,R); % LQR solver
eigs_CL = eig(Aw-Bw*Kw); % Eigenvalues of the closed loop system
eigs_CL_Realmin(ii) = min(real(eigs_CL));

%Dynamic Controller Matrices
Ac  = [0];
Bc1 = [1 0 0];
Bc2 = [-1];
Cc  = [-Kw(1)];
Dc1 = [0 -Kw(2:3)];
Dc2 = [0];
sys_control =  ss(Ac,Bc1,Cc,Dc1);

% Closed Loop System Dynamics
I = eye(size(Dc1*Dp));
Z = inv(I - Dc1*Dp); 
Acl_act = [Ap+Bp*Z*Dc1*Cp Bp*Z*Cc; Bc1*(eye(3)+Dp*Z*Dc1)*Cp Ac+Bc1*Dp*Z*Cc]; 
Bcl_act = [Bp*Z*Dc2; Bc2+Bc1*Dp*Z*Dc2];
Ccl_act = [(eye(3)+Dp*Z*Dc1)*Cp Dp*Z*Cc]; 
Dcl_act = [Dp*Z*Dc2]; 

% Produce closed loop ss system
sys_cl_act = ss(Acl_act,Bcl_act,Ccl_act,Dcl_act,'StateName',{'Az','alpha','q','dele','dele_dot'},'InputName',{'Az_cmd'},'OutputName',{'Az'});

% Step Response of CL System
opt = stepDataOptions('StepAmplitude', amp_g); % Add in 20g
[yd,td,xd] = step(sys_cl_act,opt);
S = stepinfo(sys_cl_act);

% Loop Gain
A_Lu   = [Ap zeros(4,1); Bc1*Cp Ac];
B_Lu   = [Bp; Bc1*Dp];
C_Lu   = -[Dc1*Cp Cc]; % Change sign for loop gain
D_Lu   = -[Dc1*Dp]; % Change signe for loop gain
sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu);

%% Design Charts

figure(3)
plot(-abs(real(eig(Aw))),imag(eig(Aw)),'or','markersize',10); hold on 
set(gcf,'color','w');
set(gca,'fontsize',14);
xlabel('Re','fontsize',14);
ylabel('Im','fontsize',14);
title('Root Locus of ServoMechanism: q_{11} Range 10^{-6} to 10^{-3}','fontsize',12);
grid on
%for i = 1:length(q_eI)
 %   Q(1,1) = q_eI(i);
 %   [Ks, ~, ~] = lqr(Aw, Bw, Q, R);
  %  A_cl = Aw - Bw*Ks;
  %  lamdaAcl = eig(A_cl);
  %  plot(real(lamdaAcl), imag(lamdaAcl), '.b', 'markersize', 12);
  %  axis([-25 0 -25 25])
  %  pause(0.05)
%end

figure(4)
hold on
plot(td,yd(:,1))
title('Step Response of Az')
xlabel('Time (s)')
ylabel('Az (ft/s^2)')
grid on
hold off

figure(5)
subplot(2,1,1)
plot(w_cg,sensitivity_dB); hold on
plot(w_cg(q_eI_datapoint),sensitivity_dB(q_eI_datapoint),'-o');
xlabel('Frequency (Hz)')
ylabel('Sensitivity (dB)')
grid on

subplot(2,1,2)
plot(w_cg,cosensitivity_dB); hold on
plot(w_cg(q_eI_datapoint),cosensitivity_dB(q_eI_datapoint),'-o');
xlabel('Frequency (Hz)')
ylabel('Cosensitivity (dB)')
grid on

figure(6)
margin(sys_Lu)
grid on

figure(7)
hold on
plot(w_cg,damping_min); hold on 
plot(w_cg(q_eI_datapoint),damping_min(q_eI_datapoint),'-o')
ylabel('Damping Co_eff')
xlabel('Crossover Frequency (Hz)')
title('Damping Minumum vs Crossover Freq.');
grid on
hold off

figure(8)
hold on
plot(w_cg,risetime,w_cg,settlingtime); hold on 
plot(w_cg(q_eI_datapoint),risetime(q_eI_datapoint),'-o',w_cg(q_eI_datapoint),settlingtime(q_eI_datapoint),'-o')
ylabel('Rise Time, Settling Time (sec)')
xlabel('Crossover Frequency (Hz)')
legend({'Rise Time','Settling Time'},'Location','northeast','Orientation','vertical') % Creatiing a legend to explain the plot
title('Rise time,Settling Time vs Crossover Freq.');
grid on
hold off

figure(9)
subplot(2,2,1)
plot(w_cg,OS);hold on
plot(w_cg(q_eI_datapoint),OS(q_eI_datapoint),'-o')
xlabel('Crossover Frequency (Hz)')
ylabel('% Overshoot')

grid on
hold off

subplot(2,2,2)
plot(w_cg,US);hold on
plot(w_cg(q_eI_datapoint),US(q_eI_datapoint),'-o')
xlabel('Crossover Frequency (Hz)')
ylabel('% Undershoot')
grid on
hold off 

subplot(2,2,3)
plot(w_cg,dmax);hold on
plot(w_cg(q_eI_datapoint),dmax(q_eI_datapoint),'-o')
title('Fin Displacement')
ylabel('Fin Displacement(deg)')
xlabel('Crossover Frequency (Hz)')
grid on

subplot(2,2,4)
plot(w_cg,ddmax);hold on
plot(w_cg(q_eI_datapoint),ddmax(q_eI_datapoint),'-o')
title('Fin Rate')
ylabel('Fin Rate (deg/s)')
xlabel('Crossover Frequency (Hz)')
grid on

% Nyquist
figure(10)
[Re_Lu, Im_Lu] = nyquist(sys_Lu,w_rps);
Re_Lu = squeeze(Re_Lu);
Im_Lu = squeeze(Im_Lu);
plot(Re_Lu, Im_Lu)
hold on
plot(Re_Lu, -Im_Lu)
grid on
%this adds the circle
theta = linspace(0,2*pi,1000);
xplot = (.85)*cos(theta);
yplot = (.85)*sin(theta);
plot(xplot-1,yplot)
hold on
plot(-1,0,'-x')
xlabel('Re', 'fontsize', 14);
ylabel('Im', 'fontsize', 14);
title('Nyquist', 'fontsize', 14);
axis equal
axis([-8 8 -11 11])
grid on 
set(gca,'fontsize', 14);

% Nyquist zoomed
figure(11)
[Re_Lu, Im_Lu] = nyquist(sys_Lu,w_rps);
Re_Lu = squeeze(Re_Lu);
Im_Lu = squeeze(Im_Lu);
plot(Re_Lu, Im_Lu)
hold on
plot(Re_Lu, -Im_Lu)
grid on
%this adds the circle
theta = linspace(0,2*pi,1000);
xplot = (.86)*cos(theta);
yplot = (.86)*sin(theta);
plot(xplot-1,yplot)
hold on
plot(-1,0,'-x')
xlabel('Re', 'fontsize', 14);
ylabel('Im', 'fontsize', 14);
title('Nyquist', 'fontsize', 14);
axis equal
axis([-2 1 -2 2])
grid on 
set(gca,'fontsize', 14);

figure(12)
subplot(2,1,1)
plot(w_cg,GM);hold on
plot(w_cg(q_eI_datapoint),GM(q_eI_datapoint),'-o')
title('Crossover Frequency Vs Gain Margin')
ylabel('Gain Margin (dB)')
xlabel('Crossover Frequency (Hz)')
grid on

subplot(2,1,2)
plot(w_cg,PM);hold on
plot(w_cg(q_eI_datapoint),PM(q_eI_datapoint),'-o')
title('Crossover Frequency Vs Phase Margin')
ylabel('Phase Margin (deg)')
xlabel('Crossover Frequency (Hz)')
grid on

figure(13)
plot(w_cg,q_eI); hold on
plot(w_cg(q_eI_datapoint),q_eI(q_eI_datapoint),'-o');
title('Loop Gain Crossover Freq < Actuator Natural Freq');
xlabel('Crossover Freq (Hz)')
ylabel('q_eI design points')
grid on

figure(14)
subplot(2,1,1)
semilogx(w_cg, sig_minL); hold on
plot(w_cg(q_eI_datapoint),sig_minL(q_eI_datapoint),'-o')
%xlim([10^-3 10^3])
%ylim([-20 100])
title('Sigma_min of Return Difference vs Frequency')
xlabel('Crossover Frequency (Hz)')
ylabel('Sigma (I+L,Return Difference) dB')
grid on
hold off

subplot(2,1,2)
semilogx(w_cg, sig_minLinv); hold on
plot(w_cg(q_eI_datapoint),sig_minLinv(q_eI_datapoint),'-o')
%xlim([10^-3 10^3])
%ylim([-20 50])
title('Sigma_min of Stability Robustness vs Frequency')
xlabel('Crossover Frequency (Hz)')
ylabel('Sigma (I+L^-1,Stability Robustness) dB')
grid on
hold off
sImAinv = ((w_rpsp(q_eI_datapoint))*eye(max(size(A_Lu)))-A_Lu)\eye(max(size(A_Lu)));
Lu_eval_LQR  = C_Lu*sImAinv*B_Lu + D_Lu;
%==========================================================================>
%==========================================================================>
%% Leuenberger Observer
Cw_meas  = [1 0 0 ; 0 0 1];
nCw_meas = size(Cw_meas,1);
Qo   =  diag([1e4 1e4 1e5]); 
Ro   =  [1 0; 0 1];
ol   = [0;0];

%Define Observer Model
Awob = [0 C_reg ;ol Asp];  nAw = size(Awob,1);
Bwob = [D_reg;Bsp];
nBsp = size(Bsp,1);
Bcmd = [-1;0;0];
Cwob = [1 0 0];
Dwob = 0;

% Observer Integrators
Ai  = 0; nAi = size(Ai,1);
Bi1 = [1 0];
Bi2 = -1;
Ci  = [1;0]; %Error
Di1 = [0 0; 0 1];
Di2 = [0;0];

%solve dual ARE
L = lqr(Awob',Cw_meas',Qo,Ro);
L = L';

% Building the state space model of controller,observer, integrator model
Ac  = [Ai 0 0 0; L*Ci (Awob - Bwob*(Kw) - L*Cw_meas)];
Bc1 = [Bi1;L*Di1];
Bc2 = [Bi2; L*Di2+Bcmd]; nBci2 = size(Bc2,2);
Cc  = [0 -Kw];
Dc1 = zeros(1,size(Bc1,2));
Dc2 = 0;
sys_obs_comp = ss(Ac,Bc1,Cc,Dc1);

%observer loop gain
A_Lu   = [Ap zeros(4,4); Bc1*Cp([1,3],:)  Ac];
B_Lu   = [Bp; Bc1*Dp([1,3],:) ];
C_Lu   = -[Dc1*Cp([1,3],:)  Cc]; % Change sign for loop gain
D_Lu   = -[Dc1*Dp([1,3],:) ]; % Change signe for loop gain
sys_Lu_ob = ss(A_Lu,B_Lu,C_Lu,D_Lu);

sImAinv_OB = ((w_rpsp(q_eI_datapoint))*eye(max(size(A_Lu)))-A_Lu)\eye(max(size(A_Lu)));
Lu_eval_OB  = C_Lu*sImAinv_OB*B_Lu + D_Lu;

% Closing the Loop
H = Dc1*Dp([1,3],:) ; nH = size(H,1);
IZ = inv(eye(nH) - H);
PL = Dp([1,3],:)*IZ*Dc1; nPL = size(PL,1);
IPL = eye(nPL) - PL;
A = [Ap+Bp*IZ*Dc1*Cp([1,3],:)    Bp*IZ*Cc;     Bc1*IPL*Cp([1,3],:)   Ac+Bc1*Dp([1,3],:)*IZ*Cc];
B = [Bp*IZ*Dc2 ; Bc2+Bc1*Dp([1,3],:)*IZ*Dc2];
C = [(IPL)*Cp([1,3],:)  Dp([1,3],:)*IZ*Cc];
D = Dp([1,3],:)*IZ*Dc2;

sys_combined = ss(A,B,C,D);

opt = stepDataOptions('StepAmplitude', amp_g); % Add in 20g
[yob,tob,xob] = step(sys_combined,opt);

figure(15)
plot(tob,yob(:,1))
title('Observer Step Response')
xlabel('Time (s)')
ylabel('Az (m/s^2)')
grid on

% Eigenvalue of the closed loop with the Observer
eigs_observer = eig(A);
disp('Closed Loop Observer Eigen Values')
disp(eigs_observer)

figure(16)
title('Closed Loop Observer Eigen Values')
plot(real(eigs_observer),imag(eigs_observer),'or','markersize',10)
xlabel('Real Numbers')
ylabel('Imaginary Numbers')
title('Closed Loop Observer Eigen Values')
grid on

%==========================================================================>
%==========================================================================>
%% Observer Design charts
figure(17)
hold on
plot(tob,yob(:,1))
plot(td,yd(:,1))
title('Obv+Ctrl+Az Step Response')
xlabel('Time (s)')
ylabel('Az (m/s^2)')
xlim([0 1])
legend({'Compensator','LQR'},'Location','northeast','Orientation','vertical') % Creatiing a legend to explain the plot
grid on
hold off

figure(18)
subplot(3,1,1)
hold on
plot(tob ,xob(:,1))
plot(td,xd(:,1))
title('State Estimated Vs. Actual')
ylabel('e|Az')
xlabel('Time (s)')
xlim([0 1])
grid on
hold off

subplot(3,1,2)
hold on
plot(tob,xob(:,2))
plot(td,xd(:,2))
title('State Estimated Vs. Actual')
ylabel('\alpha')
xlabel('Time (s)')
xlim([0 1])
grid on
hold off

subplot(3,1,3)
hold on
plot(tob,xob(:,3))
plot(td,xd(:,3))
title('State Estimated Vs. Actual')
ylabel('q')
xlabel('Time (s)')
xlim([0 1])
grid on
hold off

figure(19)
margin(sys_Lu_ob)
grid on

figure(20)
bodemag(sys_Lu,sys_Lu_ob)
title('Compensator Loop gain within.25Hz Of LQR Loop gain at plant input')
legend({'Compensator','LQR'},'Location','northeast','Orientation','vertical')
grid on

%==========================================================================>
%==========================================================================>
%% Guidance Law
% Define time vector
t0 = 0;  
tf = 8.86;

% Parameters
nT = 3*32.2;
HE_rad = -20*pi/180;

linespecs = {'b-', 'b-.','b-x','b-sq','b-+','b-^'};

% Initial Conditions
beta_rad = 0;
RT1      = 20000;
RT2      = 40000;
RM1      = 0;
RM2      = 40000;
VM       = V;
VT       = 300;
VT1      = -VT*cos(beta_rad);
VT2      =  VT*sin(beta_rad);

% relative positions and velocities
RTM1 = RT1 - RM1;
RTM2 = RT2 - RM2;

% relative distance
RTM = sqrt(RTM1^2 + RTM2^2);

% line of sight angle and time derivative
lambda = atan2( RTM2, RTM1 );

% missile lead angle
L = asin( VT*sin( beta_rad + lambda )/VM );

% missile velocity components
VM1  = VM*cos(lambda + L + HE_rad);
VM2  = VM*sin(lambda + L + HE_rad);

% Pointers to states [beta, RT1, RT2, RM1, RM2, VT1, VT2, VM1, VM2]
sel_beta = 1;
sel_RT1  = 2;
sel_RT2  = 3;
sel_RM1  = 4;
sel_RM2  = 5;
sel_VT1  = 6;
sel_VT2  = 7;
sel_VM1  = 8;
sel_VM2  = 9;

% Initial condition vector
y0 = [beta_rad, RT1, RT2, RM1, RM2, VT1, VT2, VM1, VM2, 0, 0, 0, 0, 0, 0, 0, 0]';

options = odeset('abstol', 1e-3, 'reltol', 1e-3);

% Integrate nonlinear 2-D engagement situation
[t,y] = ode45(@nlinpronav_project, [t0, tf], y0, options, HE_rad, 5, nT, A, B, C, D);

% target and missile velocity magnitudes
VT = sqrt( y(:,sel_VT1).^2 + y(:,sel_VT2).^2 );

% relative positions and velocities
RTM1 = y(:,sel_RT1) - y(:,sel_RM1);
RTM2 = y(:,sel_RT2) - y(:,sel_RM2);
VTM1 = y(:,sel_VT1) - y(:,sel_VM1);
VTM2 = y(:,sel_VT2) - y(:,sel_VM2);

% relative distance
RTM = sqrt(RTM1.^2 + RTM2.^2);
miss = min(RTM)

% line of sight angle and time derivative
lambda     = atan2( RTM2, RTM1 );
lambda_dot = (RTM1.*VTM2 - RTM2.*VTM1)./RTM.^2;

% closing velocity
VC = -(RTM1.*VTM1 + RTM2.*VTM2)./RTM;

% Compute acc commands
nc = 5*VC.*lambda_dot;

%achieved acceleration
Az = Z_alpha*y(:,10) + Z_dele*y(:,12);

%Guidance Plots
% Missile and target positions in inertial coordinate system
figure(21)
plot(y(:,sel_RM1), y(:,sel_RM2), 'linewidth', 2); hold on
xlabel('Downrange [ft]', 'fontsize', 14);
ylabel('Crossrange [ft]', 'fontsize', 14);
set(gca, 'fontsize', 14);
set(gcf, 'color', 'w');
grid on
h = legend('N'' = 5');
set(h,'location','best');
plot(y(:,sel_RT1), y(:,sel_RT2), 'r--', 'linewidth', 2); 
plot(y(1,sel_RM1), y(1,sel_RM2), 'ob', 'linewidth', 2);
plot(y(1,sel_RT1), y(1,sel_RT2), 'or', 'linewidth', 2);

figure(22)
plot(t, nc./32.2, 'linewidth', 2); hold on
plot(t, Az./32.2, 'linewidth', 2);
xlabel('Time [s]', 'fontsize', 14);
ylabel('Acceleration [G]', 'fontsize', 14);
title('3 G Target Maneuver','fontsize',14);
set(gca, 'fontsize', 14, 'xlim', [0 10]);
set(gcf, 'color', 'w');
grid on
h = legend({'Commanded Accel.','Acheieved Accel.'},'Location','northeast','Orientation','vertical')
set(h,'location','best');

figure(23)
subplot(2,1,1)
plot(t,y(:,12)*(180/pi));hold on
title('Fin Displacement')
ylabel('Fin Displacement(deg)')
xlabel('Time (sec)')
grid on

subplot(2,1,2)
plot(t,y(:,13)*(180/pi));hold on
title('Fin Rate')
ylabel('Fin Rate (deg/s)')
xlabel('Time (sec)')
grid on

figure(24)
subplot(3,1,1)
hold on
plot(t,y(:,14));hold on
plot(t,y(:,15));
title('State Estimated Vs. Actual')
ylabel('e|Az')
xlabel('Time (s)')
xlim([0 1])
grid on

subplot(3,1,2)
hold on
plot(t,y(:,10));hold on
plot(t,y(:,16));
title('State Estimated Vs. Actual')
ylabel('\alpha')
xlabel('Time (s)')
xlim([0 1])
grid on

subplot(3,1,3)
hold on
plot(t,y(:,11));hold on
plot(t,y(:,17));
title('State Estimated Vs. Actual')
ylabel('q')
xlabel('Time (s)')
xlim([0 1])
grid on

figure(25)
hold on
plot(t,RTM); hold on
title('Miss Distance')
xlabel('Time (s)')
ylabel('Distance [Feet]')
%xlim([0 1])
grid on
hold off

