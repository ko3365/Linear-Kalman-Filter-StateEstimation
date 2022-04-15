clc;clear;close
% Variable Init
Sig_w = 1; %Process Noise
Sig_v = 1; %Sensor Noise
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
end
figure(1)
hold on
plot(0:iteration-1,xstore(1:iteration)','k-')
plot(0:iteration-1,xhatstore','r-')
plot(0:iteration-1,xhatstore'+3*sqrt(SigXstore'), 'b-.', ...
    0:iteration-1,xhatstore'-3*sqrt(SigXstore'), 'b-.')
grid on
legend('true x','estimate x','error bounds')
hold off


