%WETLAND3P.  A 3-point dynamic model for the moprhological evolution of a backbarrier basin with marshes and mudflats
%Copyright (C) 2014, Giulio Mariotti
%Developer can be contacted by <email> and <paper mail>
%This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.function dX=funMARSH(X,PAR)

clear;close all;clc

%just to avoid long time close to dfo~0
OPT=odeset('AbsTol',10^-6,'RelTol',10^-6,'Events',@POOLstopp5);
OPTfzero=optimset('Algorithm','Levenberg-Marquardt','TolFun',10^-28,'TolX',10^-28,'MaxFunEvals',10000);

rhos=1000;rhom=1000;%bulk densities [kg/m3]
B=10000;%basin width [m]
P=12.5*3600*1; % tidal period [s]
ws=0.5*10^-3; %settlign velocity[m/s]
tcr=0.1; %tau cr mudflat [Pa]
wind=6; %reference wind speed [m/s]
amp=1.4/2; % tidal amplitude [m]
Ba=2; %marsh progradation coeff [-]
Be=0.16/(365*24*3600); %marsh erosion coeff [m/yr/ (W/m)]
lamda=0.0001; % mudflat erodability coeff [-]
dist=10; %reference dist from marsh bank [m].
Co=50/1000; % ref conc. [kg/m3]
RSLR=5*(10^-3)./(3600*24*365); %relative sea level rise [m/s]

%marsh parameters
BMax=2.500;% kg/m2
nuGp=.0138;
por=1000/2650;
chiref=0.15;

%time
years=100000;TN=1*years+1;to=linspace(1,3600*24*365*years,TN);





%%%%%%%%%%%%%
%OPTION 1, single plot
%%%%%%%%%%%%%%%%

%initial conditions
dfo=3;%initial mudfflat depth
dmo=0.01; %initial marsh depth
bfo=B/2; %initial mudflat width

PAR=[rhos P B ws tcr Co wind Ba Be amp RSLR NaN rhom lamda dist bfo BMax nuGp por chiref];
[t,X]=ode23s(@(t,X) funMARSH(X,PAR),to,[bfo dfo dmo],OPT);bf=X(:,1);df=X(:,2);dm=X(:,3);bm=B-bf;

%plot
Dmax=.7167*2*amp-.483; 
bf=bf./1000;
axis([0 B./1000 -7.5 1.5])
plot(bf,-df,'-b',bf,-dm,'-g');hold on; 
plot(bf(end),-df(end),'ob',bf(end),-dm(end),'dg',[0 B],0*[1 1],'k',[0 B],-Dmax.*[1 1],'k');hold on; 
plot(bf(1),-df(1),'xr',bf(1),-dm(1),'xr')
xlabel('b_f [%]')
ylabel('d_m, d_f [m]')
axis([0 B/1000 -4.5 0.1])




%%%%%%%%%%%%%
%OPTION 2, multi plot. Recreate Figure 2 in Mariotti and Carr (2014) WRR
%%%%%%%%%%%%%%%%
% CO=[10 50 100]./1000;
% RSLR=[ 5 1 5].*(10^-3)./(3600*24*365);
% in=[.1 .5 .9]; 
%  
% figure
% for i=1:3;
% for j=1:3;bfo=B.*in(j);
%         
% Dmax=.7167*2*amp-.483;                      
% PAR=[rhos P B ws tcr CO(i) wind Ba Be amp RSLR(i) NaN rhom lamda dist bfo BMax nuGp por chiref];
% deq=fsolve(@(X) eq1(X,PAR),2,OPTfzero);
% if deq>0.5;dfo=deq;else;dfo=0.5;end
% dmo=Dmax/2;
% [t,X]=ode23s(@(t,X) funMARSH(X,PAR),to,[bfo dfo dmo],OPT);bf=X(:,1);df=X(:,2);dm=X(:,3);bm=B-bf; 
% 
% subplot(3,1,i)
% bf=bf./1000;
% axis([0 B./1000 -7.5 1.5])
% hold on; 
% plot(bf,-df,'-b',bf,-dm,'-g');hold on; 
% plot(bf(end),-df(end),'ob',bf(end),-dm(end),'dg',[0 B],0*[1 1],'k',[0 B],-Dmax.*[1 1],'k');hold on; 
% plot(bf(1),-df(1),'xr',bf(1),-dm(1),'xr')
% xlabel('b_f [%]')
% ylabel('d_m, d_f [m]')
% axis([0 B/1000 -4.5 0.1])
% 
% end
% end
