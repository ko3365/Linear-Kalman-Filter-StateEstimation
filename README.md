# Linear Kalman Filter
## Part 1. Linear Kalman Filter Derivation
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

### Applying to an Arbitrary System:
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
The real state is unknown and we estimate the state using the Kalman Filter derived from above. The figure below shows the comparison between the real state and estimated state:
<p align="center">
  <img 
    width="600"
    src="images/1.linear_kf.png"
  >
</p>

## Part 2. Steady State Kalman Filter (Infinite Horizon)
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

Using the steady-state approach, the kalman gain and the state variance are not a function of step k. It is computationally simple with a slight penalty on optimality. The comparison between estimated state and real state is shown in the figure below (if we compare with the Linear KF approach above, the result is almost the same even if we used the constant Kalman gain:

<p align="center">
  <img 
    width="600"
    src="images/2.linear_kf_ss.png"
  >
</p>

