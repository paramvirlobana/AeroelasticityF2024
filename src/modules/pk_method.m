% MECH 6481 Project Fall 2024
% Team #1 P-K Method

clc;
clear all;
close all;

% Inputs for validation
a = -0.25;
e = -0.1;
mu = 20;
r_squared = 6/25;
sigma = 0.4;
x_theta = e - a;

% Setup
V_vec = linspace(0.1, 10, 1000);
N = length(V_vec);
converged_roots = zeros(N, 4);
k1 = zeros(N, 4);
gamma = zeros(N, 4);
k = 5;                  % Initial guess
err = 1;        
iteration = 0;

% P-K Loop
for m = 1:N
    for n = 1:4
        while err >= 0.0001
            V = V_vec(m);
            C_func = ((0.01365 + (0.2808 * k * 1i) - ((k.^ 2) / 2))) / (0.01365 + (0.3455 * k * 1i) - (k.^2));
            f11 = ((sigma^ 2) / (V^ 2)) - ((k.^ 2) / (mu)) + ((2 * 1i * k.* C_func) / (mu));
            f12 = ((k.*(1i + (a.*k))) + ((2 + ((1i.*k).*(1 -(2.*a)))).*C_func)) / (mu);
            f21 = ((a * (k.^2)) - ((1i * k) * (1 + (2 * a)).*C_func)) / (mu);
            f22 = (((8.*mu.*r_squared) / (V.^2)) + ((4.*1i) * (1 + (2.*a)) * ((2.*1i) - (k.*(1 - (2.*a)))).*C_func) - (k.*(k - (4.*1i) ...
                + ((8.*a).*(1i + (a.*k)))))) / (8.*mu);
            p = [r_squared - (x_theta.^2), 0, f22 + (r_squared.*f11) - (x_theta.*f21) - (x_theta.*f12), 0, (f11.*f22) - (f12.*f21)];
            polynomial_roots = roots(p);
            [~,sorted_indices] = sort(imag(polynomial_roots));
            sorted_roots = polynomial_roots(sorted_indices);
            k_new = imag(sorted_roots(n));
            gamma = real(sorted_roots(n)) / k_new;
            err = abs(k - k_new);
            k = k_new;

        end
        converged_roots(m,n) = sorted_roots(n);
        k1(m,n) = k_new;
        gamma(m,n) = gamma;
        iteration = 0;
        err = 1;
    end
end

converged_roots1=converged_roots(:,1)';
converged_roots2=converged_roots(:,2)';
converged_roots3=converged_roots(:,3)';
converged_roots4=converged_roots(:,4)';

% Plotting
figure;
plot(V_vec, V_vec.* real(converged_roots1));
hold on;
plot(V_vec, V_vec.* real(converged_roots2));
plot(V_vec, V_vec.* real(converged_roots3));
plot(V_vec, V_vec.* real(converged_roots4));
ylabel('Γ/ωθ');
xlabel('V');
legend('Real Converged root 1','Real Converged root 2', 'Real Convergerd root 3','Real Converged root 4')

figure;
plot(V_vec, V_vec.* imag(converged_roots1));
hold on;
plot(V_vec, V_vec.* imag(converged_roots2));
plot(V_vec, V_vec.* imag(converged_roots3));
plot(V_vec, V_vec.* imag(converged_roots4));
ylabel('Ω/ωθ');
xlabel('V');
legend('Imaginary Converged root 1','Imaginary Converged root 2', 'Imaginary Converged root 3','Imaginary Converged root 4')