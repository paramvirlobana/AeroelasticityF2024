% MECH 6481 Project
% P-method

clc;
clear
close all;

% Inputs for validation:
a = -0.2;
e = -0.1;
mu = 20;
r_squared = 0.24;
sigma = 0.4;
N = 100001; % Change as wanted for higher accuracy
 
% Setup
xtheta = e-a;
V_vec = linspace(0,4,N); 
lambda = zeros(4,N); 

% Roots for polynomial 
root1 = zeros(1,N); 
root2 = zeros(1,N); 
root3 = zeros(1,N); 
root4 = zeros(1,N); 

% P-METHOD
% ---- zeta_h and zeta_theta changing ----
zeta_h = 0;
zeta_theta = 0;

for i = 2:N
    V = V_vec( i ); % reuse V
    
    p = [r_squared-xtheta^2 ... %P4
        ((2*(r_squared)/V)*zeta_theta+(r_squared*2*(sigma/V)*zeta_h)) ... %P3
        (r_squared/V^2)-(2/mu)*(a+0.5)+(2*zeta_h*sigma/V)*2*(r_squared/V)*zeta_theta+r_squared*((sigma^2)/(V^2))-xtheta*2/mu ... %P2
        ((2*(sigma/V)*zeta_h)*((r_squared/V^2)-2/mu*(a+0.5)))+(sigma^2/V^2)*2*(r_squared/V)*zeta_theta ... %P1
        (sigma/V)^2*(r_squared/V^2-(2/mu)*(a+0.5))];%P0

    root = roots(p); % roots function
    root1(i) = root(1);
    root2(i) = root(2);
    root3(i) = root(3); 
    root4(i) = root(4);
end

% =================== PLOTTING ======================

% Dimensionless frequency
figure;
plot(V_vec,V_vec.*imag(root1),'ro','MarkerSize',10); 
hold on;
plot(V_vec,V_vec.*imag(root2),'bd','MarkerSize',8);
plot(V_vec,V_vec.*imag(root3),'k*','MarkerSize',6);
plot(V_vec,V_vec.*imag(root4),'gs','MarkerSize',4);
xlabel('V');
ylabel('\Omega/\omega_\theta');
title(sprintf('Dimensionless flutter speed versus Reduced Speed for \\zeta_h = %.2f and \\zeta_\\theta = %.2f', zeta_h, zeta_theta));

grid on;
hold off;

% Dimensionless damping
figure;
plot(V_vec,V_vec.*real(root1),'ro','MarkerSize',10);
hold on;
plot(V_vec,V_vec.*real(root2),'bd','MarkerSize',8);
plot(V_vec,V_vec.*real(root3),'k*','MarkerSize',6);
plot(V_vec,V_vec.*real(root4),'gs','MarkerSize',4);
xlabel('V');
ylabel('\Gamma/\omega_\theta');
title(sprintf('Dimensionless damping versus Reduced Speed for \\zeta_h = %.2f and \\zeta_\\theta = %.2f', zeta_h, zeta_theta));
grid on;
hold off;