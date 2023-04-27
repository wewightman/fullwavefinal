%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WREN WIGHTMAN (adapted from scripts by GIANMARCO PINTON)
% WRITTEN: MAR 29, 2023
% Launch Fullwave code, easy matlab wrapper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare folder and workspace
clear all
close all

fwpath = '/hpc/group/ultrasound/wew12/modules/fullwave_bmme890/';
addpath(fwpath);
%%
command = sprintf("rm *.dat\nrm runit\ncp %s/try6_nomex runit", fwpath);
system(command);
%% Set up problem
% read in parameter file and map
mapparams = readjson('matparams.json');

% read in material parameters
matkey = readjson('matkey.json');

% read in imaging paramters
imparams = readjson('impar.json');

% matparams = [...
%     % speed of sound (m/s), density, attenuation @ 1MHz, Non-Linearity,
%     % fractional standard deviation
%     1465,  985, 0.40, 8.5, 0.010;...  % fatty tissue
%     1580, 1050, 0.74, 6.6, 0.001;...  % skeletal muscle
%     1613, 1120, 1.57, 0.0, 0.015;...  % connective tissue
%     ];

%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = mapparams.c0;    % speed of sound (m/s)
omega0 = mapparams.f0*2*pi;
lambda = mapparams.lam;
nx = mapparams.Nx;  % number of lateral elements
nz = mapparams.Nz;  % number of depth elements
duration = 80e-6;  % duration of simulation (s)
p0 = 1e5; % pressure in Pa

matmap = readmap('map.bin', nx, nz);
%%
theta = imparams.theta;
seed = 0;
rng(seed);

%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = mapparams.ppw;           % number of points per spatial wavelength
cfl = 0.4;         % Courant-Friedrichs-Levi condition

%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nT = round(duration*c0/lambda*ppw/cfl);
dx = lambda/ppw;
dz = lambda/ppw;
dT = dx/c0*cfl;

wx = dx*nx;
wz = dz*nz;

%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = ones(nx, nz)*1540;   % speed of sound map (m/s)
rhomap = ones(nx, nz)*1000; % density map (kg/m^3)
Amap = ones(nx, nz)*0.0;    % attenuation map (dB/MHz/cm)
boveramap = -2*ones(nx, nz);    % nonlinearity map 

%% Import material models
keys = ["skin", "fat", "con", "musc", "inter"];
for key = keys
    mask = matmap == matkey.(key);  % get mask number for this material
    mparams = matkey.params.(key);  % get parameters for this material
    N = sum(mask, "all");
    cmap(mask)      = normrnd(mparams(1), mparams(1)*mparams(5), N,1);
    rhomap(mask)    = normrnd(mparams(2), mparams(2)*mparams(5), N,1);
    Amap(mask)      = normrnd(mparams(3), mparams(3)*mparams(5), N,1);
    boveramap(mask) = normrnd(mparams(4), mparams(4)*mparams(5), N,1);
end

%%
%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nx, nz); 
inmap(:,1:3) = 1;
incoords = mapToCoords(inmap); % note zero indexing for compiled code

%%
%%% Generate initial conditions based on input coordinates %%%%%%

t = (0:nT-1)*dT;                    % time vector of sim

dele = mapparams.probe.pitch;

xz = [dx, dz] .* incoords;          % spatial coordinates
xz = [dele*(floor(xz(:,1)/dele)+1/2), xz(:,2)];
nsteer = [sin(theta), cos(theta)];  % normal vector of steered plane wave
tau = (nsteer * xz')'/c0;           % calculate the travel time to source pixels
tau = tau - min(tau);

tau0 = 2*pi/omega0;
sigma = 0.8*tau0;

tprime = t - tau0 - tau;
icmat = p0 * exp(-((tprime)/sigma).^2).*cos(omega0*tprime);

%%
%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
outmap = zeros(nx,nz);
outmap(:, 4) = 1;
outcoords = mapToCoords(outmap);

%% Run the simulation
%%% Launch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
launchTotalFullWave2(c0,omega0,wx,wz,duration,p0,ppw,cfl,cmap',rhomap',Amap',boveramap',incoords,outcoords,icmat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
!./runit
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dispaly results
ncoordsout=size(outcoords,1);
nRun=sizeOfFile('genout.dat')/4/ncoordsout;
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));

%% convert to nele channel binary file
dele = mapparams.probe.pitch;
wele = mapparams.probe.width;
nele = mapparams.probe.noElements;

channels = zeros(size(genout,1), nele);
for iele = 1:nele
    indmin = max(1, round(ppw*dele*(iele-1)/lambda));
    indmax = min(round(indmin+ppw*wele/lambda), size(genout,2));
    channels(:,iele) = sum(genout(:,indmin:indmax), 2);
end


%% save channel data
fileID = fopen('channels.bin','w');
fwrite(fileID,channels,'double');
fclose(fileID);

chanparams.nT = nT;
chanparams.nele = nele;
chanparams.dele = dele;
chanparams.dT = dT;
chanparams.theta = theta;

fileID = fopen('channels.json','w');
fwrite(fileID,jsonencode(chanparams), 'char');
fclose(fileID);