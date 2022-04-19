clc;clear;close
%% Linear Kalman Filter

% Variable Init
Sig_w = 5; %Process Noise
Sig_v = 5; %Sensor Noise
A = 1; B = 1; C = 2; D = 0; % Simple State-Space Model

%Initial State
x = 0; xhat = 0; u = 0;
SigX = 0;

%Variable Store
iteration = 50;
xstore = zeros(1,iteration+1);
xstore(1,1) = x; % initial value
xhatstore = zeros(1,iteration);
SigXstore = zeros(1^2,iteration);

% For Later Use: Steady-state KF
ustore = zeros(1,iteration+1);
ustore(1,1) = u;
zstore = zeros(1,iteration);

for k=1:iteration
    %true value simulation
    u_prev = u;
    u = sin(k/(2*pi))+randn/3;
    w = chol(Sig_w)'*randn;
    v = chol(Sig_v)*randn;
    z = C*x+D*u+v;
    x = A*x+B*u+w;
    
    %Kalman Filter Prediction 1 (state estimate prediction)
    xhat = A*xhat+B*u_prev;
    %Kalman Filter Prediction 2 (covariance update)
    SigX = A*SigX*A'+Sig_w;
    %Kalman Filter Prediction 3 (output prediction)
    zhat = C*xhat+D*u;

    %Kalman Filter Update 1 (Kalman Gain)
    L = SigX*C'/(C*SigX*C'+Sig_v);
    %Kalman Filter Update 2 (state estimate update)
    xhat = xhat + L*(z-zhat);
    %Kalman Filter Update 3 (Error covariance update)
    SigX = SigX-L*C*SigX;

    xstore(1,k+1) = x;
    xhatstore(1,k) = xhat;
    SigXstore(1,k) = SigX;

    %for SS KF Comparison
    ustore(1,k+1) = u;
    zstore(1,k)=z;
end

%% Steady-state Kalman Filter (infinite horizon)

%Obtain Steady-State Kalman Gain
hamiltonian = [A^(-1)' A^(-1)'*C'*Sig_v^(-1)*C; Sig_w*A^(-1)' A+Sig_w*A^(-1)'*C'*Sig_w^(-1)*C];
[evector,evalue] = eig(hamiltonian);
psi12 = evector(1,2);
psi22 = evector(2,2);
sigX_ss_minus = psi22*psi12^(-1);
L_ss = sigX_ss_minus*C'*(C*sigX_ss_minus*C'+Sig_v)^(-1);
sigX_ss_plus = sigX_ss_minus-L_ss*C*sigX_ss_minus;

%initialization
x2 = 0; xhat2 = 0;
xhatstore2 = zeros(1,iteration);

for t=1:iteration
    u2_prev = ustore(1,t);
    u2 = ustore(1,t+1);
    z2 = zstore(1,t);
    
    %Kalman Filter Prediction 1 (state estimate prediction)
    xhat2 = A*xhat2+B*u2_prev;
    %Kalman Filter Prediction 3 (output prediction)
    zhat2 = C*xhat2+D*u2;

    %Kalman Filter Update 2 (state estimate update)
    xhat2 = xhat2 + L_ss*(z2-zhat2);


    xhatstore2(1,t) = xhat2;

end

figure(1)
hold on
plot(0:iteration-1,xstore(1:iteration)','k-')
plot(0:iteration-1,xhatstore','r-')
plot(0:iteration-1,xhatstore'+3*sqrt(SigXstore'), 'b-.', ...
    0:iteration-1,xhatstore'-3*sqrt(SigXstore'), 'b-.')
grid on
legend('true x','estimate x (Linear KF)','error bounds')
hold off

figure(2)
hold on
plot(0:iteration-1,xstore(1:iteration)','k-')
plot(0:iteration-1,xhatstore2','r-')
plot(0:iteration-1,xhatstore2'+3*sqrt(sigX_ss_plus'), 'b-.', ...
    0:iteration-1,xhatstore2'-3*sqrt(sigX_ss_plus'), 'b-.')
grid on
legend('true x','estimate x (SS KF)','error bounds')
hold off
