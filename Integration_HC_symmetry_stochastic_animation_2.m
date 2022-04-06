%This program simulate four degrees of freedom half-car system response for
%harmonic excitation
close all


%% System inputs 
MRf=0.6;   %Motion ratio front wheel
MRr=0.6;   %Motion ratio rear wheel
m1f=20;    %Half front unsprung mass
m1r=20;    %Half rear unsprung mass
m2=(316/2)-m1f-m1r+140;  %Half sprung mass
kt=200000;     %Tyre vertical stiffnes
ksr=38000;      %Rear suspension spring
kr=ksr*(MRr^2);  %Rear suspension vertical stiffness
ksf=32000;      %Front suspension spring
kf=ksf*(MRf^2);  %Front suspension vertical stiffness
cdr=1500;       %Rear damper
cr=cdr*(MRr^2);  %Rear vertical damping
cdf=2000;       %Front damper
cf=cdf*(MRf^2);  %Front vertical damping
lr=1.189;    %Distance from front axle to GC %1189/1000;
lf=1.411;    %Distance from rear axle to GC %1411/1000
L=lr+lf;         %Wheel base
J=m2*lf*lr;      %Car moment of inertia  

%Critical damping
ckr=2*sqrt(kr*m2);
cpass=0.2*ckr;
g=9.81;

%% Stochastic excitation

ISOk = 7; % Values For ISO Road A-B Roughness Classification, from 3 to 9
N  = 2500; % Number of data points
RoadL  = 250;  % Length Of Road Profile (m)
B  = L/N ; % Sampling Interval (m)
dn = 1/RoadL;  % Frequency Band
n0 = 0.1;        % Spatial Frequency (cycles/m)
n  = dn : dn : N*dn; % Spatial Frequency Band
phi =  2*pi*rand(size(n)); % Random Phase Angle
Amp1 = sqrt(dn)*(2^ISOk)*(1e-3)*(n0./n); % Amplitude for Road  Class A-B
x = 0:B:RoadL-B;                          % Abscissa Variable from 0 to L
v = 10;
% tt=L/v;
ttf = x/v;
t0=0;
tf = x(end)/v;
delay=2.6/v;
ttr = ttf-delay; %Delay zmienione na '-' dle próby !!!!!
hf = zeros(size(x));
hr = zeros(size(x));
for i=1:length(ttf)
    hf(i) = sum(Amp1.*cos(2*pi*n*ttf(i)+ phi));
    hr(i) = sum(Amp1.*cos(2*pi*n*ttr(i)+ phi));
end

s_time=RoadL/v;
tspan=[t0,tf];
step=s_time/N*2*2;

%% Initial conditions 
y1f0=-(m2+m1f)*10*lr/(kt*L); y1r0=-(m2+m1r)*10*lf/(kt*L); y1f0dot=0; y1r0dot=0; y20=-0.1; y20dot=0; theta0=0; theta0dot=0;
IC=[y1f0, y1f0dot, y1r0, y1r0dot, y20, y20dot, theta0, theta0dot];

%% Call ODE solver
%%%%%% Ta funkcja u¿ywa interpolacji wyników do animacji %%%%%%%% [0:step:tf-step]
options=odeset('MaxStep',1/100);
[time, state_values] = ode45(@(t,s) HC_symmetry_stochastic_interpolation(s,t,m1f,m1r,m2,kt,kf,kr,cr,cf,J,g,hf,hr,ttf,ttr,lr,lf,L),[0:step:tf-step],IC,options);
%% Extract individaul values
n_of_del=floor(1*length(time));
time = wkeep(time,n_of_del,'r');
y1f = state_values(:,1);        %Front wheel displacement
y1f = wkeep(y1f,n_of_del,'r');
y1fdot = state_values(:,2);     %Front wheel velocity
y1fdot = wkeep(y1fdot,n_of_del,'r');
y1r = state_values(:,3);        %Rear wheel displacement
y1r = wkeep(y1r,n_of_del,'r');
y1rdot = state_values(:,4);     %Rear wheel velocity
y1rdot = wkeep(y1rdot,n_of_del,'r');
y2 = state_values(:,5);         %Car vertical displacement
y2 = wkeep(y2,n_of_del,'r');
y2dot  = state_values(:,6);     %Car vertical velocity
y2dot = wkeep(y2dot,n_of_del,'r');
theta  = state_values(:,7);     %Car pitch angle
theta = wkeep(theta,n_of_del,'r');
thetadot  = state_values(:,8);  %Car pitch angular velocity
thetadot = wkeep(thetadot,n_of_del,'r');

y2f=y2+theta*lf;    %Front suspension mounting displacement
y2fdot=y2dot+thetadot*lf;  %Front suspension mounting vel
y2r=y2-theta*lr;    %Rear suspension mounting displacement
y2rdot=y2dot-thetadot*lr;    %Rear suspension mounting vel
RSf=y2f-y1f;        %Rattle space front
RSr=y2r-y1r;        %Rattle space rear

%% Ground
Front_exct=zeros(size(time));
Rear_exct=zeros(size(time));

for i=1:length(time);
Front_exct(i)=sum(Amp1.*cos(2*pi*n*time(i) + phi));
Rear_exct(i)=sum(Amp1.*cos(2*pi*n*(time(i) - delay)+ phi));
end

%% Accelerations
y1fdotdot = (-kf*(y1f-y2-lf*theta)-kt*(y1f-Front_exct)-cf*(y1fdot-y2dot-lf*thetadot)-m1f*g)/m1f;
y1rdotdot = (-kr*(y1r-y2+lr*theta)-kt*(y1r-Rear_exct)-cr*(y1rdot+lr*thetadot-y2dot)-m1r*g)/m1r;
y2dotdot = (-kr*(y2-lr*theta-y1r)-kf*(y2-y1f+lf*theta)-cr*(y2dot-lr*thetadot-y1rdot)-cf*(y2dot-y1fdot+lf*thetadot)-m2*g)/m2;
thetadotdot = (+lr*kr*(y2-lr*theta-y1r)+lr*cr*(y2dot-lr*thetadot-y1rdot)-lf*kf*(y2+lf*theta-y1f)-lf*cf*(y2dot+lf*thetadot-y1fdot))/J;
y2fdotdot=y2dotdot+thetadotdot*lf;    %Rear suspension mounting vel
y2rdotdot=y2dotdot-thetadotdot*lr;    %Rear suspension mounting vel

%% Animation 

timeall=time;
x_road=time*v;
roadall=time*v;
InitBody=y2(1);
zeros13=zeros(1,13);
% y1r=[zeros13;y1r];
% y1r=y1r(1:end-13,:);

    for i=1:length(time)
        Plotting_solar_stochastic(RSr(i),RSf(i),InitBody,v,x_road(i),roadall,lr,lf,theta(i),y2(i),y1f(i),y1r(i),time(i),timeall,Front_exct);
        
        drawnow
        pause(0.1);
        F(i)=getframe(gcf);
    end

video = VideoWriter('Animation','MPEG-4');%MPEG-4
open(video);  %%Blad bo zmieniaja sie axis w czasie symulacji
writeVideo(video,F);
close(video)

%% Result ploting

figure
subplot(3,2,1);
plot(time,Rear_exct);
hold on;
plot(time,Front_exct);
ylabel('Ground at rear wheel [m]');
xlabel('Time [s]');

subplot(3,2,2);
plot(time,Front_exct);
ylabel('Ground at front wheel [m]');
xlabel('Time [s]');
% axis([ 15 20 -inf inf]);

subplot(3,2,3);
plot(time,y1r);
hold on;
plot(time,y1f);
ylabel('Rear wheel displ [m]');
xlabel('Time [s]');
%axis([3.5 7 -0.02 0.02])

subplot(3,2,4);
plot(time,y1f);
hold on;
plot(time,0);
ylabel('Front wheel displ [m]');
xlabel('Time [s]');
% axis([15 20 -inf inf])

subplot(3,2,5);
plot(time,y2);
ylabel('Car vertical displ [m]');
xlabel('Time [s]');

subplot(3,2,6);
plot(time,theta);
ylabel('Car pitch angle [rad]');
xlabel('Time [s]');

figure
plot(time,RSf);
ylabel('Front relative displ [m]');
xlabel('Time [s]');

figure
plot(time,RSr);
ylabel('Front relative displ [m]');
xlabel('Time [s]');

% figure
% plot(time,y1f);
% ylabel('Front wheel displ [m]');
% xlabel('Time [s]');
% 
% figure
% plot(time,y1r);
% ylabel('Rear wheel displ [m]');
% xlabel('Time [s]');
% 
% figure
% plot(time,y2);
% ylabel('Car vertical displ [m]');
% xlabel('Time [s]');
% 
% figure
% plot(time,theta);
% ylabel('Theta [rad]');
% xlabel('Time [s]');

% y2max=max(abs(y2))
% Amplitude=max(abs(y2))-min(abs(y2))