clc;clear;close all;
%% MIMO Linear Kalman Filter

% Variable Init
SigmaW = [0.5 0;0 1]; %Process Noise
SigmaV = [.5 0;0 .5]; %Sensor Noise
A = [.1 .2; .5 .2]; B = [2 0;1 2]; C = [.2 0;-.1 .3]; D = 0; % Simple State-Space Model

%Initial State
x = [0;0]; xhat =[0;0]; u = [0;0];
SigX = zeros(length(A));

%Variable Store
iteration = 40;
xstore = zeros(length(x),iteration+1); 
xstore(:,1) = x;
xhatstore = zeros(length(xhat),iteration);
SigXstore = zeros(length(xhat)^2,iteration);

%For Later Use: Comparing the result with Robust KF
ustore = zeros(length(u),iteration+1); 
ustore(:,1)=u;
zstore = zeros(length(C*x),iteration);

for k = 1:iteration
    % True value simulation
    u_prev = u;
    u = [0.5*randn(1) + cos(k/pi); 0.3*randn(1)-sin(k/pi)]; % for example...
    w = chol(SigmaW)'*randn(size(x,1),size(x,2));
    v = chol(SigmaV)'*randn(size(C*x,1),size(C*x,2));
    z = C*x + D*u + v; % z is based on present x and u
    x = A*x + B*u + w; % future x is based on present u

    % Kalman Filter Prediction
    xhat = A*xhat + B*u_prev;
    SigX = A*SigX*A' + SigmaW;
    zhat = C*xhat + D*u;

    % Kalman Filter Correction
    L = SigX*C'/(C*SigX*C' + SigmaV);
    xhat = xhat + L*(z - zhat);
    SigX = SigX - L*C*SigX;

    % Store for the simulation plot
    xstore(:,k+1) = x; xhatstore(:,k) = xhat;
    SigXstore(:,k) = SigX(:);

    % Store to apply the same input and output to the robust KF
    ustore(:,k+1)=u;
    zstore(:,k)=z;
end

figure(1)
subplot(2,1,1);
hold on
plot(0:iteration-1,xstore(1,1:iteration)','k-')
plot(0:iteration-1,xhatstore(1,:)','r-')
plot(0:iteration-1,xhatstore(1,:)'+3*sqrt(SigXstore(1,:)'), 'b--', ...
    0:iteration-1,xhatstore(1,:)'-3*sqrt(SigXstore(1,:)'), 'b--')
grid on
legend('true x1','estimate x1','bounds')
ylim([-10,10])
hold off

subplot(2,1,2)
hold on
plot(0:iteration-1,xstore(2,1:iteration)','k-')
plot(0:iteration-1,xhatstore(2,:)','r-')
plot(0:iteration-1,xhatstore(2,:)'+3*sqrt(SigXstore(4,:)'), 'b--', ...
    0:iteration-1,xhatstore(2,:)'-3*sqrt(SigXstore(4,:)'), 'b--')
grid on
legend('true x2','estimate x2','bounds')
ylim([-10,10])
hold off

%% Robust Linear KF (Sequential Measurement Processing + Joseph Form Covariance Update)

%initialization
xhat2 = [0;0];
SigX2 = zeros(length(A));
xhatstore2 = zeros(length(xhat2),iteration);
SigXstore2 = zeros(length(xhat)^2,iteration);
for t=1:iteration
    u2_prev = ustore(:,t);
    u2 = ustore(:,t+1);
    z2 = zstore(:,t);

    % Kalman Filter Prediction
    xhat2= A*xhat2 + B*u2_prev;
    SigX2 = A*SigX2*A' + SigmaW;
    zhat2 = C*xhat2 + D*u2;
    
    % Kalman Filter Correction (Using the Sequential Measurement Processing!)
    for m=1:length(C*xhat2)
        C_m = C(m,:)';
        zvar = C_m'*SigX2*C_m+SigmaV(m,m);
        L2 = SigX2*C_m/zvar;
        xhat2 = xhat2+L2*(z2(m)-C_m'*xhat2);
        % (Joseph Form Covariance Update)
        SigX2 = (eye(length(xhat2))-L2*C_m')*SigX2*(eye(length(xhat2))-L2*C_m')'+L2*zvar*L2';
    end
    

    xhatstore2(:,t) = xhat2;
    SigXstore2(:,t) = SigX2(:);
end

figure(2)
subplot(2,1,1);
hold on
plot(0:iteration-1,xstore(1,1:iteration)','k-')
plot(0:iteration-1,xhatstore2(1,:)','r-')
plot(0:iteration-1,xhatstore2(1,:)'+3*sqrt(SigXstore2(1,:)'), 'b--', ...
    0:iteration-1,xhatstore2(1,:)'-3*sqrt(SigXstore2(1,:)'), 'b--')
grid on
legend('true x1','estimate x1 (robust KF)','bounds')
ylim([-10,10])
hold off

subplot(2,1,2)
hold on
plot(0:iteration-1,xstore(2,1:iteration)','k-')
plot(0:iteration-1,xhatstore2(2,:)','r-')
plot(0:iteration-1,xhatstore2(2,:)'+3*sqrt(SigXstore2(4,:)'), 'b--', ...
    0:iteration-1,xhatstore2(2,:)'-3*sqrt(SigXstore2(4,:)'), 'b--')
grid on
legend('true x2','estimate x2 (robust KF)','bounds')
ylim([-10,10])
hold off
