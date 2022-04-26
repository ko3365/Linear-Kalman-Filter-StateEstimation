# 1. Linear Kalman Filter Derivation
Kalman Filter is an estimation method that uses a series of measurments observed over time. It implements minimum mean square error that minimize the difference between the true state and state estimate:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\hat{x}^{MMSE}=\text{min}_{\hat{x}}(\mathbb{E}[||x-\hat{x}||^2\&space;|\&space;\mathbb{Z}])&space;}">
Z is the series of measurements. Thus, given all the measurement, we want to find an estimate of the state to minimize the mean square error.
Setting the derivative to 0, the optimal state estimation can be obtained: 

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\frac{d}{d\hat{x}}\mathbb{E}[x^Tx-2x^T\hat{x}&plus;\hat{x}^T\hat{x}|\mathbb{Z}]=0}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\therefore\hat{x}=\mathbb{E}[\&space;x\&space;|\&space;\mathbb{Z}\&space;]}">

### Deriving Optimal Estimate
First, define the prediction error and solve for the state estimation:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\tilde{x}_k^-=x_k-\hat{x_k}^-&space;\&space;\text{where}\&space;\hat{x}_k^-=\mathbb{E}[x_k\&space;|\&space;\mathbb{Z}_{k-1}]}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\mathbb{E}[\tilde{x}_k^-&space;|&space;\mathbb{Z}_k]=\mathbb{E}[x_k|\mathbb{Z}_k]-\mathbb{E}[\hat{x}_k^-|\mathbb{Z}_k]}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\hat{x}^&plus;=\hat{x}_k^-&plus;\mathbb{E}[\tilde{x}_k^-|z_k]}">

Since, all densities are assumed to be Gaussian, we can find the mean of the conditional PDF of estimation error and measurement to obtain the expected value:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}f(\tilde{x}_k^-|z_k)=\frac{f(\tilde{x}_k^-,z_k)}{f(z_k)}}">
The mean of the conditional pdf and kalman gain L is: 

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\mathbb{E}[\tilde{x}_k^-|z_k]=\mathbb{E}[\tilde{x}_k^-]&plus;\Sigma_{\tilde{x}\tilde{z},k}^{-1}\Sigma_{\tilde{z},k}^{-1}(z-\mathbb{E}[z])=\Sigma_{\tilde{x}\tilde{z},k}^{-1}\Sigma_{\tilde{z},k}^{-1}\tilde{z}_k}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}L_k&space;=&space;\Sigma_{\tilde{x}\tilde{z},k}^{-1}\Sigma_{\tilde{z},k}^{-1}}&space;">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\hat{x}_k^&plus;=\hat{x}_k^-&plus;L_k\tilde{z}_k}">

### Prediction and Correction
Thus, we have the equations required to obtain the optimal estimation: (Prediction + Correction)

#### Prediction 
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\hat{x}_k^-=\mathbb{E}[f_{k-1}(x_{k-1},u_{k-1},w_{k-1})|\mathbb{Z}_{k-1}]}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\Sigma_{\tilde{x},k}^-=\mathbb{E}[(\tilde{x}_k^-)(\tilde{x}_k^-)^T]}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\hat{z}_k=\mathbb{E}[y_k(x_k,u_k,v_k)|\mathbb{Z}_{k-1}]}">

#### Correction
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}L_k=\Sigma_{\tilde{x}\tilde{z},k}^-\Sigma_{\tilde{z},k}^{-1}}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\hat{x}_k^&plus;=\hat{x}_k^-&plus;L_k(z_k-\hat{z}_k)}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\Sigma_{\tilde{x},k}^&plus;=\Sigma_{\tilde{x},k}^--L_k\Sigma_{\tilde{z},k}L_k^T}">

# 2. Applying to an Arbitrary System (SISO and Infinite Horizon):
Now, let's apply the Kalman Filter to an arbitrary system A=1, B=1, C=2, D=0 with process and measurement noise random variance 5. 
```Matlab
Sig_w = 5; %Process Noise
Sig_v = 5; %Sensor Noise
A = 1; B = 1; C = 2; D = 0; % Simple State-Space Model
```
Also, arbitrarily define the input (deterministic)
```Matlab
u = sin(k/(2*pi))+randn/3;
```
The iteration steps of prediction and correction to obtain the optimal estimation
```Matlab
for k=1:iteration
    %Kalman Filter Prediction
    xhat = A*xhat+B*u_prev;
    SigX = A*SigX*A'+Sig_w;
    zhat = C*xhat+D*u;

    %Kalman Filter Correction 
    L = SigX*C'/(C*SigX*C'+Sig_v);
    xhat = xhat + L*(z-zhat);
    SigX = SigX-L*C*SigX;
end
```


The real state is unknown and we estimate the state using the Kalman Filter derived from above. The figure below shows the comparison between the real state and estimated state:
<p align="center">
  <img 
    width="600"
    src="images/1.linear_kf.png"
  >
</p>

## Steady State Kalman Filter (Infinite Horizon)
It is possible for the Kalman Filter recursion to converge to a unique steady state solution when the system is non-time varying, {A,C} observable, and {A,Process noise variance} controllerable. 

### Hamiltonian Approach. 
If we define our state prediction variance to be:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\Sigma_{\tilde{x},k}^-=S_kZ_k^{-1}}">
We can derive the following relationship:
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\begin{bmatrix}Z_{k&plus;1}\\S_{k&plus;1}\end{bmatrix}=\mathcal{H}\begin{bmatrix}Z_k\\S_k\end{bmatrix}}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\mathcal{H}=\begin{bmatrix}A^{-T}&A^{-T}C^T\Sigma_{\tilde{v}}^{-1}C\\\Sigma_{\tilde{w}}A^{-T}&A&plus;\Sigma_{\tilde{w}}A^{-T}C^T\Sigma_{\tilde{v}}^{-1}C\end{bmatrix}}">
Perform the diagonal transformation of the Hamiltonian Matrix and define Y:
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\mathcal{H}=\Psi\begin{bmatrix}\Lambda^{-1}&0\\0&\Lambda\end{bmatrix}\Psi^{-1}}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\begin{bmatrix}Y_{1,k}\\Y_{2,k}\end{bmatrix}=\Psi^{-1}\begin{bmatrix}Z_k\\S_k\end{bmatrix}}">
If the equation is from k steps after time zero and assuming k is arbitrarily high number,

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\begin{bmatrix}Z_k\\S_k\end{bmatrix}=\begin{bmatrix}\Psi_{11}&\Psi_{12}\\\Psi_{21}&\Psi_{22}\end{bmatrix}\begin{bmatrix}\Lambda^{-k}Y_{1,0}\\\Lambda^kY_{2,0}\end{bmatrix}=\begin{bmatrix}\Psi_{11}&\Psi_{12}\\\Psi_{21}&\Psi_{22}\end{bmatrix}\begin{bmatrix}0\\\Lambda^kY_{2,0}\end{bmatrix}}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\therefore\Sigma_{\tilde{x}}^-=S_kZ_k^{-1}=\Psi_{22}\Psi_{12}^{-1}}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}L=\Sigma_{\tilde{x}}^-C^T(C\Sigma_{\tilde{x}}^-C^T&plus;\Sigma_{\tilde{v}})^{-1}}">
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\Sigma_{\tilde{x}}^&plus;=\Sigma_{\tilde{x}}^--LC\Sigma_{\tilde{x}}^-}">

To achieve steady state kalman gain and state variance:
```Matlab
hamiltonian = [A^(-1)' A^(-1)'*C'*Sig_v^(-1)*C; Sig_w*A^(-1)' A+Sig_w*A^(-1)'*C'*Sig_w^(-1)*C];
[evector,evalue] = eig(hamiltonian);
psi12 = evector(1,2);
psi22 = evector(2,2);
sigX_ss_minus = psi22*psi12^(-1);
L_ss = sigX_ss_minus*C'*(C*sigX_ss_minus*C'+Sig_v)^(-1);
sigX_ss_plus = sigX_ss_minus-L_ss*C*sigX_ss_minus;
```


Using the steady-state approach, the kalman gain and the state variance are not a function of step k. It is computationally simple with a slight penalty on optimality. The comparison between estimated state and real state is shown in the figure below (if we compare with the Linear KF approach above, the result is almost the same even if we used the constant Kalman gain:

<p align="center">
  <img 
    width="600"
    src="images/2.linear_kf_ss.png"
  >
</p>

# 3. Applying to an Arbitrary System (MIMO and Robust Model):
If we have two inputs and two outputs with two states, the A,B,C matrix is now a 2x2 matrix:
```Matlab
SigmaW = [0.5 0;0 1]; %Process Noise
SigmaV = [.5 0;0 .5]; %Sensor Noise
A = [.1 .2; .5 .2]; B = [2 0;1 2]; C = [.2 0;-.1 .3]; D = 0; % Simple State-Space Model
```
The iteration steps of prediction and correction to obtain the optimal estimation are shown above, and the figure below shows the comparison between the two real states and estimated states:

<p align="center">
  <img 
    width="600"
    src="images/MIMO_plot.png"
  >
</p>

## Robust Kalman Filter (Joseph Form Covariance Update & Sequential Processing of Measurement)
During the covariance update, we can lose the postive definiteness of matrix due to the rounding error due to the subtraction:

### Joseph Form Covariance Update
<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\Sigma_{\tilde{x},k}^&plus;=\Sigma_{\tilde{x},k}^--L_kC_k\Sigma_{\tilde{x},k}^-}">
We can convert this to the Joseph form so that the subtraction occurs in the "square" term, and thus, the postive definiteness of covariance matrix is guaranteed:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}\Sigma_{\tilde{x},k}^&plus;=[I-L_kC_k]\Sigma_{\tilde{x},k}^-[I-L_kC_k]^T&plus;L_K\Sigma_{\tilde{v}}L_k^T}">

### Sequential Processing of Measurement
If there exists large number of sensors, the Kalman gain calculation can be computationally intensive due to the matrix inverse operation _O_(m<sup>3</sup>), where m is the number of sensors. If we break the m measurements into single measurement, we can improve the efficiency improves to _O_(m<sup>2</sup>). The measurments can be splitted as shown below:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}z_k=\begin{bmatrix}z_{k:1}\\\vdots\\z_{k:m}\end{bmatrix}=\begin{bmatrix}C_{k:1}^Tx_k&plus;v_{k:1}\\\vdots\\C_{k:m}^Tx_k&plus;v_{k:m}\end{bmatrix}}">

And now we need to perform m steps of prediction correction:

<img src="https://latex.codecogs.com/svg.image?\large&space;{\color{Gray}&space;L_{k:i}&space;=&space;\frac{\Sigma_{\tilde{x},k:i-1}^&plus;C_{k:i}}{\sigma_{\tilde{z},k:i}^2}&space;&space;}">




## References
[1] Optimal State Estimation, D. Simon, 2006 

[2] State Estimation: Applied Kalman Filter, G. Plett (Lecture Notes)
