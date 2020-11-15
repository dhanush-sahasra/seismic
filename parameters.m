%% Analysis of Seismic Parameters for Ridgecrest earthquake

%% Accelogram Data
a = importdata('accnew.txt');
a = reshape(a',[], 1).*10^-2; %m/sec^2

v = importdata('velnew.txt');
v = reshape(v',[], 1).*10^-2; %m/sec^2

% Time step (delta T)
dt=0.02; % sec
NumPoints=numel(a);
% Time vector
t = zeros(1, numel(a));
t(1) = 0;
for i = 2:numel(a)
    t(i) = t(i-1) + dt;
end
%% 1. PGA
pga = max(abs(a)); %m/sec^2

%% 2. PGV
pgv = max(abs(v)); %m/sec

%% 3. PGA/PGV
peak_ratio = pga/pgv;

%% 4. Arias Intensity
% Values of Arias Intensity
ariasTH = 1/9.81*cumsum(a.^2)*pi*dt/2;
% Total Arias Intensity at the end of the ground motion
arias = ariasTH(end);

%% 5. RMS acceleration
rmsaTH = sqrt((cumsum(a.^2)*dt)/t(end));
rmsa = rmsaTH(end); %m/sec^2

%% 6. Strong motion duration of Trifunac/Brady
% CUMULATIVE ENERGY
% Values of cumulative energy
EcumTH = cumsum(a.^2)*dt;
% Total cumulative energy at the end of the ground motion
Ecum = EcumTH(end);
% time history of the normalized cumulative energy
Ecum_Norm = EcumTH/Ecum;

timed = t(EcumTH>=0.05*Ecum & EcumTH<=0.95*Ecum);
% starting and ending points of the significant duration
t_5_95 = [timed(1),timed(end)]; %D5-95 parameter
% significant duration
Td = timed(end)-timed(1); %sec

%% 7. Seismic power
H_1 = cumsum(a(t<=timed(1)).^2)*dt;
H_2 = cumsum(a(t<=timed(end)).^2)*dt;
P = (H_2(end) - H_1(end))/Td;

%% 8. Spectral Intensity of Housner
g = 9.81;
Ag=a(:,1);          % Acceleration Vector
zet=2;              % Damping Ratio (%)
endp=5;             % End Period of Spectra (sec)

[T, Spa, Spv, Sd, Sv, Sa]=SPECTRA(dt, Ag, zet, g, endp);
Spv1 = cumsum(Spv(t<=2.5));
Spv2 = cumsum(Spv(t<=0.1));
SI_H =  Spv1(end)- Spv2(end);

%% 9. Spectral Intensity of Kappos
Tn = 1.2; %sec == funamental period of the structure
tk = 0.2*Tn; 
LL = Tn - tk; %lower limit of integration
UL = Tn + tk; %upper limit of integration
Spv3 = cumsum(Spv(t<=UL));
Spv4 = cumsum(Spv(t<=LL));
SI_K = Spv3(end) - Spv4(end); 

%% 10. Spectral Intensity of Martinez Rueda
Th = 1.5; %sec == hardening period
Ty = 1.2; %sec == yield period
Spv5 = cumsum(Spv(t<=Th));
Spv6 = cumsum(Spv(t<=Ty));
SI_MR = (Spv5(end) - Spv6(end))/(Th - Ty);

%% 11. Effective Peak Acceleration
g = 9.81;
Ag=a(:,1);          % Acceleration Vector
zet=5;              % Damping Ratio (%)
endp=5;             % End Period of Spectra (sec)

[T, Spa2, Spv2, Sd2, Sv2, Sa2]=SPECTRA(dt, Ag, zet, g, endp);

SA = mean(Sa2(6:26));
EPA = SA/2.5;

%% 12. EPA max
SA_sw = zeros(length(Sa2),1);
EPA_sw = zeros(length(Sa2),1);

for i = 2:230
    SA_sw(i) = mean(Sa2(i:i+20));
    EPA_sw(i) = SA_sw(i)/2.5;
end

EPA_max = max(abs(EPA_sw));

%% 13. Seismic Energy Input
g = 9.81;
Ag=a(:,1);          % Acceleration Vector
zet=5;              % Damping Ratio (%)
endp=t(end);             % End Period of Spectra (sec)

[T3, Spa3, Spv3, Sd3, Sv3, Sa3]=SPECTRA(dt, Ag, zet, g, endp);

% E_inpTH = zeros(length(Sa3),1);
% for i =2:t(end)
xggt = a(1:15055,1);
xgt = v(1:15055,1);
E_inpTH = dt*cumsum((Sa3+xggt).*xgt);
% end
E_inp = E_inpTH(end);

%% 14. Cumulative Energy Input
CAVTH = cumsum(abs(a))*dt;
CAV = CAVTH(end);

%% 15. Seismic damage potential of Araya/Sarogoni
b = zeros(numel(a),1);
for i = 2:numel(a)
    if a(i) <= 0 && a(i-1) >= 0
        b(i) = 1;
    else
        if a(i) >= 0 && a(i-1) <= 0
            b(i) = 1;
        end
    end
end
v0TH = cumsum(b);
v0 = v0TH(end)./t(end);

DP_AS = arias/(v0^2);

%% 16. Central Period
c = zeros(numel(a),1);
for i = 2:numel(a)
    if a(i) <= 0 && a(i-1) >= 0
        c(i) = 1;
    end
end
npTH = cumsum(c);
np = npTH(end)./t(end);
CP = 1/np;

%% 17. Spectral Displacement
figure(1)
plot(T,Sd);
xlabel('Period');
ylabel('Sd');
title('Spectral Displacement');

%% 18. Spectral Velocity
figure(2)
plot(T,Sv);
xlabel('Period');
ylabel('Sv');
title('Spectral Velocity');

%% 19. Spectral Acceleration
figure(3)
plot(T,Sa);
xlabel('Period');
ylabel('Sa');
title('Spectral Acceleration');

%% 20. Intensity of Fajifar/Vidic/Fishinger
I_FVF = pgv*(Td^0.25);