% ********************** quadrotor control with mpc ********************* %
close all; clearvars;

%% model and simulation parameters

% drone parameters
m = 4; % mass of drone
g = 9.8; % gravity constant
% drone inertia in body frame
Ix = 0.033; 
Iy = 0.033;
Iz = 0.066;
l = 0.5; % length of drone arm
a =10^-6; % rotors thrust coefficient T = a*omega^2
b =2*10^-7; % rotors counter-torque coefficient t = (+-)b*omega^2

% simulation parameters
tf = 100; % final time
ts = 0.05; % sampling time
N = 100; % prediction horizon
Nf = tf/ts+1; % simulation steps number
t = linspace(0,tf,Nf);

% system matrices
A = zeros(14,14);
B = zeros(14,4);
C = zeros(4,14);
for i=1:13
    A(i,i+1) = 1;
end
A(4,5) = 0;
A(8,9) = 0;
A(12,13) = 0;
B(4,1) = 1;
B(8,2) = 1;
B(12,3) = 1;
B(14,4) = 1;
C(1,1) = 1;
C(2,5) = 1;
C(3,9) = 1;
C(4,13) = 1;
% discrete system matrices
sysd = c2d(ss(A,B,C,[]),0.05,'zoh');
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
% state observer
Q = 0.01*(Bd'*Bd);
R = 0.0001*eye(4);
[~,L,~] = kalman(sysd,Q,R,[]);
Ad = sparse(Ad);
Bd = sparse(Bd);
Cd = sparse(Cd);
% mapping matrices
Phi_map = zeros(14*N,14);
Psi_map = zeros(14*N,4*N);
for i=1:N
    Phi_map(1+14*(i-1):14*i,:) = Ad^i;
    for j=1:N
        Psi_map(1+14*(i-1):14*i,1+4*(j-1):4*j) = (i>=j)*Ad^(i-j)*Bd;
    end
end
Phi_map = sparse(Phi_map);
Psi_map = sparse(Psi_map);
Cd_map = kron(speye(N),Cd);

%% Initialization

% trajectory to track
Yt = [30*sin(0.2*t);30*cos(0.2*t);-20-0.5*t;pi/4*ones(1,Nf)];

% bernstein polynomials
n = 4; % degree of polynomials
P = zeros(1,n);
Pm = zeros(4*N,4*n);
for i=1:N
    for j=1:n
        P(j) = nchoosek(n-1,j-1)*((i-1)/(N-1))^(j-1)*((N-i)/(N-1))^(n-j);
    end
    Pm(1+4*(i-1):4*i,:) = blkdiag(P,P,P,P);
end

% initial state and control
X0 = [-10 0 0 0 -5 0 0 0 2.5 0 0 0 0 0]';
Xold = X0;
Xold_true = X0;
p0 = zeros(4*n,1);

% state and control history
Xlist = zeros(14,Nf);
Xlist(:,1) = X0;
Xlist_true = zeros(14,Nf);
Xlist_true(:,1) = X0;
Vlist = zeros(4,Nf);

% linear inequality constraints matrices
Ax = zeros(N,14*N);
Ay = zeros(N,14*N);
Az = zeros(N,14*N);
Apsi = zeros(N,4*N);
for i=1:N
    Ax(i,3+14*(i-1)) = 1;
    Ay(i,7+14*(i-1)) = 1;
    Az(i,11+14*(i-1)) = 1;
    Apsi(i,4*i) = 1;
end
Ax = sparse(Ax);
Ay = sparse(Ay);
Az = sparse(Az);
Apsi = sparse(Apsi);
Alin_ineq = [Ax*Psi_map*Pm;-Ax*Psi_map*Pm;Ay*Psi_map*Pm;-Ay*Psi_map*Pm;...
    Az*Psi_map*Pm;-Az*Psi_map*Pm;Apsi*Pm;-Apsi*Pm];


%% simulation loop

tic
for k = 1:Nf-1
    if rem(k,10)==0
        disp(k);
    end

    % reducing prediction horizon and depending matrices near end
    if Nf-k<N
        N = Nf-k;
        Pm = Pm(1:4*N,:);
        Phi_map = Phi_map(1:14*N,:);
        Psi_map = Psi_map(1:14*N,1:4*N);
        Cd_map = kron(speye(N),Cd);
        Ax = Ax(1:N,1:14*N);
        Ay = Ay(1:N,1:14*N);
        Az = Az(1:N,1:14*N);
        Apsi = Apsi(1:N,1:4*N);
        Alin_ineq = [Ax*Psi_map*Pm;-Ax*Psi_map*Pm;Ay*Psi_map*Pm;...
            -Ay*Psi_map*Pm;Az*Psi_map*Pm;-Az*Psi_map*Pm;Apsi*Pm;-Apsi*Pm];
    end

    % solving mpc optimization
    axy_lim = 2;
    az_lim = 1.5;
    apsi_lim = 0.1;
    blin_ineq = -[Ax*Phi_map*Xold-axy_lim*ones(N,1);...
        -Ax*Phi_map*Xold-axy_lim*ones(N,1);...
        Ay*Phi_map*Xold-axy_lim*ones(N,1);...
        -Ay*Phi_map*Xold-axy_lim*ones(N,1);...
        Az*Phi_map*Xold-az_lim*ones(N,1);...
        -Az*Phi_map*Xold-az_lim*ones(N,1);...
        -apsi_lim*ones(N,1);-apsi_lim*ones(N,1)];
    fun = @(x) cost_function(Phi_map,Psi_map,Cd_map,x,Yt(:,k+1:k+N),Xold,Pm);
    HessFunc = @(x,lambda) HessianJ(x,lambda,Pm,Psi_map,Cd_map);
    options = optimoptions("fmincon","Algorithm","interior-point",...
        "SpecifyObjectiveGradient",true,"HessianFcn",HessFunc,"Display","none");
    p = fmincon(fun,p0,Alin_ineq,blin_ineq,[],[],[],[],[],options);
    p0 = p;
    % apply first control and update states and history lists
    v = Pm*p;
    vk = v(1:4);
    Vlist(:,k) = vk;
    Xnew_true = Ad*Xold_true+Bd*vk;
    Y_mes = Cd*Xnew_true+randn(4,1)*0.0001;
    Xold_true = Xnew_true;
    Xlist_true(:,k+1) = Xnew_true;
    Xnew = Ad*Xold+Bd*vk+L*(Y_mes-Cd*Xold);
    Xold = Xnew;
    Xlist(:,k+1) = Xnew;
end
toc

%% control inputs and remaining states from flat output

theta = atan((Xlist(3,:).*cos(Xlist(13,:))+Xlist(7,:).*sin(Xlist(13,:)))...
    ./(Xlist(11,:)-g));
phi = atan(sin(theta).*(Xlist(7,:).*cos(Xlist(13,:))-Xlist(3,:).*sin(Xlist(13,:)))...
    ./(Xlist(3,:).*cos(Xlist(13,:))+Xlist(7,:).*sin(Xlist(13,:))));
phi(1) = 0;
U = -m*sqrt(Xlist(3,:).^2+Xlist(7,:).^2+(Xlist(11,:)-g).^2);
theta_dot = zeros(1,Nf);
phi_dot = zeros(1,Nf);
theta_sec = zeros(1,Nf);
phi_sec = zeros(1,Nf);
for k=2:Nf
    theta_dot(k) = (theta(k)-theta(k-1))/ts;
    phi_dot(k) = (phi(k)-phi(k-1))/ts;
    theta_sec(k) = (theta_dot(k)-theta_dot(k-1))/ts;
    phi_sec(k) = (phi_dot(k)-phi_dot(k-1))/ts;
end
% computes angular rates and accelerations
rate_p = phi_dot-sin(theta).*Xlist(14,:);
rate_q = cos(phi).*theta_dot + sin(phi).*cos(theta).*Xlist(14,:);
rate_r = -sin(phi).*theta_dot + cos(phi).*cos(theta).*Xlist(14,:);
acc_p = phi_sec-sin(theta).*phi_sec - theta_dot.*cos(theta).*Xlist(14,:);
acc_q = cos(phi).*theta_sec + sin(phi).*cos(theta).*Vlist(4,:) - phi_dot.*sin(phi).*theta_dot ...
    + phi_dot.*cos(phi).*cos(theta).*Xlist(14,:) - theta_dot.*sin(theta).*sin(phi).*theta_dot;
acc_r = -sin(theta).*theta_sec + cos(phi).*cos(theta).*Vlist(4,:) - phi_dot.*cos(phi).*theta_dot ...
    -phi_dot.*sin(phi).*cos(theta).*Xlist(14,:) - theta_dot.*sin(theta).*cos(phi).*Xlist(14,:);
tau_phi = acc_p*Ix-(Iy-Iz)*rate_r.*rate_q;
tau_theta = acc_q*Iy-(Iz-Ix)*rate_r.*rate_p;
tau_psi = acc_r*Iz;

% MIX-RPM: compute motors rotation speeds from control inputs
mix_mat = [a a a a;...
    l*a*cosd(45) -l*a*cosd(45) -l*a*cosd(45) l*a*cosd(45);...
    l*a*sind(45) l*a*sin(45) -l*a*sin(45) -l*a*sin(45); ...
    -b b -b b];
omega = sqrt(mix_mat\[-U;tau_phi;tau_theta;tau_psi]);

%% plot results

figure();
plot(t,Yt(1,:),'--b',t,Xlist(1,:),'-r');
title('x plot');
figure();
plot(t,Yt(2,:),'--b',t,Xlist(5,:),'-r');
title('y plot');
figure();
plot(t,Yt(3,:),'--b',t,Xlist(9,:),'-r');
title('z plot');
figure();
plot(t,Yt(4,:),'--b',t,Xlist(13,:),'-r');
title('psi plot');
figure();
plot(t,phi);
title('phi plot');
figure();
plot(t,theta);
title('theta plot');
figure();
plot(t,-U,'-r',t,33.2*ones(1,Nf),'--k',t,46.6*ones(1,Nf),'--k');
title('U plot');
figure();
plot(t,tau_phi);
title('tau phi plot');
figure();
plot(t,tau_theta);
title('tau theta plot');
figure();
plot(t,tau_psi);
title('tau psi plot');
figure();
plot(t,omega);
title('motors speeds');
