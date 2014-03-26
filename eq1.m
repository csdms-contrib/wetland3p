function er=eq1(X,PAR)   
rhos=PAR(1);
P=PAR(2);
B=PAR(3);
ws=PAR(4);
tcr=PAR(5);
Co=PAR(6);
wind=PAR(7);
Ba=PAR(8);
Be=PAR(9);
amp=PAR(10);
RSLR=PAR(11);
rhom=PAR(13);
lamda=PAR(14);
dist=PAR(15);

%%%%%%%%%%%%
fetch=PAR(16);%X(1);
df=X(1);
dm=0;%X(3);
bm=B-fetch;

%%%SALT MARSH GROWTH%%%%%%%%%%%%%%%%%%%%%%%%%%%5
BMax=2.500;% kg/m2
Dmin=0;
Dmax=.7167*2*amp-.483;
delT=0;
sigB=0.06;%degrees C-1=
AA=.25.*(-Dmin-Dmax)*(Dmax-3*Dmin);
Bpeak=BMax*(dm-Dmax)*(dm-Dmin)/AA;
if (Bpeak<=1e-3);Bpeak=0;end

Bfrac=(Bpeak/BMax);

nuGp=.0138;
AMC=(180.)*Bpeak*(nuGp)/(365*24*3600);%
%%%%%%%%%5
por=1000/2650;
chiref=0.15;%
Rref=AMC*chiref;
po=1000; %dens of org
FFm=(1/por)*(Rref/po);

%%%%%%%%%%%%%%%%%%%%%%%%%
%average df
fac=min(1,df/(2*amp));
fac2=min(1,dm/(2*amp));
Df=(df+(df-fac*2*amp))/2;
Dm=(dm+(dm-fac2*2*amp))/2; 

tw=  wavetau(fetch,wind,Df,B);
if Dm >1e-4  
    
    tw2=wavetauBmod(fetch,wind,Dm,Df,Bfrac,B);
else
    tw2=0;
end

tau=max((tw-tcr)/tcr,0)*lamda;
tau2=max((tw2-tcr)/tcr,0)*lamda;
Cr=rhos*tau/(1+tau);
Cm=rhos*tau2/(1+tau2);
hb=dm+(df-dm)*(1-exp(-dist*0.1/df)); %scarp height according to profile
W=waveTRNS(amp,df,Df,dist,wind,fetch,dm,hb);

E=(Be*W/(hb-dm)-Ba*Cr*ws/rhom);%*(hb-dm)/hb
Fm=(Cr-Cm)*min(2*amp,dm)/P/rhom; %flux from tidal flat to marsh
Fc=(Cr-Co)*(fac*2*amp)/P/rhom;

%[Bfrac tw tw2 dm Dmax Dm]



er=-E*(df-dm)/fetch +Fm*bm/fetch +(Cr-Co)*(fac*2*amp)/P/rhom   +RSLR; %*1000/(1000+fetch)  +(Cr-Co);  
%%%VERDE%%