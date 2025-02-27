Case Study word document : 

<a href="/case_study.doc" target="_blank">Case Study</a>

Sinsyn_man word document 

<a href="/sinsyn_man.doc" target="_blank">Sinsyn_Man</a>

``` data title="Excieee1.dat"

.397
.09
0.01
.02
1.0
1.0
400
.06
1.0
6
-6

```


``` data title="excst1.dat"

1
10
400
.025
6
-6
```


``` data title="hydrogov.dat"

5
.2
.04
.3
.05
0
1
1.2
.1
```


``` data title="machine.dat"
0.00327
1.5845
1.7572
1.04
0.4245
6.66
0.44
3.542
0
50
1
```


``` data title="pss.dat"
10
12
.048
.032
0
0
6
-6
```


``` octave title="Object C source code.m"

echo on
% This program is associated with singen.mdl which simulates the dynamics of a 
% SINGLE MACHINE INFINITE BUS(SMIB)SYSTEM. The Synchronous generator is powered by
% hydraulic/steam turbine, with excitation system(IEEE Type-1/ Static) and PSS. For simulation 
% of SMIB system, type " singen " at the MATLAB command window.The Simulink model for SMIB
% system is opened. Double click  on the " more info " button for further instructions.
echo off

clear all
clc

disp 'enter the time at which simulation has to end'
t_sim=input('simulation time''\n');

load machine.dat

ra=machine(1,1);
xq=machine(2,1);
xd=machine(3,1);
xq1=machine(4,1);
xd1=machine(5,1);
tdo1=machine(6,1);
tqo1=machine(7,1);
h=machine(8,1);
d=machine(9,1);
f=machine(10,1);
eb=machine(11,1);
wref=2*pi*f;
wb=2*pi*f;
w=wb;
%Network data
disp 'enter h and z parameters'
h1=input('real part of h parameter''\n');
h2=input('imaginary part of h parameter''\n');
zr=input('real part of z parameter''\n');
zi=input('imaginary part of z parameter''\n');

h12=(h1+sqrt(-1)*h2);
z11=(zr+sqrt(-1)*zi);

%Excitation system

 exsys=input('if excitation system is included enter 1 or else 0''\n');
 if exsys==1
 type=input('enter type of excitation system-1for IEEE type-1,2 for static excitation system''\n');
 
 	if type==2
   avr_chk_2=1;avr_chk_1=0;
   exc_choice=2;%choice of static excitation system

   load excst1.dat
       
   tc_exsys=excst1(1,1);
   tb_exsys=excst1(2,1);
   ka2=excst1(3,1);
   ta2=excst1(4,1);
   efdmax=excst1(5,1);
   efdmin=excst1(6,1);
   %Dummy values for ieee type-1  excitation system
   a_exsys=0.3970;   b_exsys= 0.0900;   tr_exsys=0.0100;   ta1=0.0200;   tf_exsys=1.0000;   te_exsys=1.0000;
   ka1=400.0;   kf_exsys=0.06;   ke_exsys=1.0;   vrmax=1;vrmin=-1;
   
	else
   avr_chk_1=1;avr_chk_2=0;
   exc_choice=1; % choice of ieee type-1 excitation system

   load excieee1.dat
   
   a_exsys=excieee1(1,1);
   b_exsys=excieee1(2,1);
   tr_exsys=excieee1(3,1);
   ta1=excieee1(4,1);
   tf_exsys=excieee1(5,1);
   te_exsys=excieee1(6,1);
   ka1=excieee1(7,1);
   kf_exsys=excieee1(8,1);
   ke_exsys=excieee1(9,1);
   vrmax=excieee1(10,1);
   vrmin=excieee1(11,1);
   % Dummy value  for static excitation system
   tc_exsys=1.0000;tb_exsys=10;ka2=400;ta2=0.025;efdmax=1;efdmin=-1;
   end
   
else 
   avr_chk_1=0;avr_chk_2=0;
   exc_choice=3;% excitation system is not used 
   tc_exsys=1.0000;tb_exsys=10;ka2=400;ta2=0.025;efdmax=1;efdmin=-1;
   a_exsys=0.3970;   b_exsys= 0.0900;   tr_exsys=0.0100;   ta1=0.0200;   tf_exsys=1.0000;   te_exsys=1.0000;
   ka1=400.0;   kf_exsys=0.06;   ke_exsys=1.0;   vrmax=1;vrmin=-1;
end
if exsys~=0 
pss=input('If pss is included enter1 or else 0''\n');
if pss==1
   pss_chk=1;
   pss_choice=1;
	load pss.dat
  
   ks=pss(1,1);
   tw=pss(2,1);
   t1=pss(3,1);
   t2=pss(4,1);
   t3=pss(5,1);
   t4=pss(6,1);
   vsmax=pss(7,1);
   vsmin=pss(8,1);
else
   pss_chk=0;
   pss_choice=2;
   tw=0;ks=0;t1=1;t2=1;t3=1;t4=1;vsmax=1;vsmin=-1;pss_ch=0;
   %vs=0;
end
%  Extra added Codes
else
      pss_choice=2;
      tw=0;ks=0;t1=1;t2=1;t3=1;t4=1;vsmax=1;vsmin=-1;pss_ch=0;
end
turb=input('enter 1 if turbine is included else 0''\n');
if turb==1
   type=input('1 for hydroturbine 2 for steam turbine''\n');
   if type==1
   hydro_chk=1;steam_chk=0;
   gov_choice=1;%choice of hydro governor
	load hydrogov.dat
   
   Tr_ht=hydrogov(1,1);
   Tg_ht=hydrogov(2,1);
   Tp_ht=hydrogov(3,1);
   del_ht=hydrogov(4,1);
   sig_ht=hydrogov(5,1);
   Ta_ht=(1/sig_ht)*Tr_ht*Tg_ht;
   Tb_ht=(1/sig_ht)*(((del_ht+sig_ht)*Tr_ht)+Tg_ht);
   t1_ht=(Tb_ht/2)+sqrt((Tb_ht*Tb_ht/4)-Ta_ht);
   t3_ht=(Tb_ht/2)-sqrt((Tb_ht*Tb_ht/4)-Ta_ht);
   t2_ht=hydrogov(6,1);
   k_ht=1/sig_ht;
   tw_ht=hydrogov(7,1);
   p_hydromax=hydrogov(8,1);
   p_hydromin=hydrogov(9,1);
   % Dummy values for steam turbine
    tch_st=0.3; trh_st=8; tco_st=0.4; fhp_st=0.3; fip_st=0.3; p_steammin=p_hydromin;p_steammax=p_hydromax;
    flp_st=0.4; t1_st=0.3; t2_st=0; t3_st=0.1; k_st=20; p_up=p_hydromax;p_down=p_hydromin;

   else
   hydro_chk=0;steam_chk=1;
   gov_choice=2;%choice of steam governor
	load steamgov.dat
	   
   tch_st=steamgov(1,1);
   trh_st=steamgov(2,1);
   tco_st=steamgov(3,1);
   fhp_st=steamgov(4,1);
   fip_st=steamgov(5,1);
   flp_st=steamgov(6,1);
   t1_st=steamgov(7,1);
   t2_st=steamgov(8,1);
   t3_st=steamgov(9,1);
   sig_st=steamgov(10,1);
   k_st=1/sig_st;
   p_up=steamgov(11,1);
   p_down=steamgov(12,1);
   p_steammax=steamgov(13,1);
   p_steammin=steamgov(14,1);
   %Dummy values for hydro turbine
   t1_ht=38.4803; t2_ht=0; t3_ht=0.5197; tw_ht=1; k_ht=20; p_hydromax=p_steammax;p_hydromin=p_steammin;     
   end
else
   hydro_chk=0;steam_chk=0;
   gov_choice=3; % governor is not considered
 
    tch_st=0.3; trh_st=8; tco_st=0.4; fhp_st=0.3; fip_st=0.3; p_steammin=0.1;p_steammax=1.1;
    flp_st=0.4; t1_st=0.3; t2_st=0; t3_st=0.1; k_st=20; p_up=1.1;p_down=0.1;
    t1_ht=38.4803; t2_ht=0; t3_ht=0.5197; tw_ht=1; k_ht=20; p_hydromax=1.1;p_hydromin=0.1;   
end

   
eb=input('infinite bus voltage''\n');
   
delbo=0;
ebc=eb*(cos(delbo)+sqrt(-1)*sin(delbo));
EbQ=real(ebc);
EbD=imag(ebc);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
VtQ=1;VtD=0;IQ=0.6;ID=0;% initial conditions for newton's method of network solution

disp 'enter 1 for PV bus 2 for PQ bus'
pvpqch=input('choice''\n');

if pvpqch==1
   disp 'enter P'
   pg=input('Pg''\n');
   vg=input('Vg''\n');
   Vtm=vg;
   val(1,1)=VtQ-zr*IQ+zi*ID-h1*EbQ+h2*EbD;
   val(2,1)=VtD-zi*IQ-zr*ID-h1*EbD-h2*EbQ;
   val(3,1)=Vtm^2-VtQ^2-VtD^2;
   val(4,1)=pg-VtD*ID-VtQ*IQ;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
while max(abs(val))>1e-14
J(1,1)=1;
J(1,3)=-zr;
J(1,4)=zi;

J(2,2)=1;
J(2,3)=-zi;
J(2,4)=-zr;

J(3,1)=-2*VtQ;
J(3,2)=-2*VtD;

J(4,1)=-IQ;
J(4,2)=-ID;
J(4,3)=-VtQ;
J(4,4)=-VtD;

update=-J\val;
VtQ=VtQ+update(1);
VtD=VtD+update(2);
IQ=IQ+update(3);
ID=ID+update(4);

val(1,1)=VtQ-zr*IQ+zi*ID-h1*EbQ+h2*EbD;
val(2,1)=VtD-zi*IQ-zr*ID-h1*EbD-h2*EbQ;
val(3,1)=Vtm^2-VtQ^2-VtD^2;
val(4,1)=pg-VtD*ID-VtQ*IQ;

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

else
   disp 'Enter Pg and Qg'
   pg=input('Pg''\n');
   qg=input('Qg''\n');
   val(1,1)=VtQ-zr*IQ+zi*ID-h1*EbQ+h2*EbD;
   val(2,1)=VtD-zi*IQ-zr*ID-h1*EbD-h2*EbQ;
   val(3,1)=qg-VtD*IQ+VtQ*ID;
   val(4,1)=pg-VtD*ID-VtQ*IQ;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
while max(abs(val))>1e-14
J(1,1)=1;
J(1,3)=-zr;
J(1,4)=zi;

J(2,2)=1;
J(2,3)=-zi;
J(2,4)=-zr;

J(3,1)=ID;
J(3,2)=-IQ;
J(3,3)=-VtD;
J(3,4)=VtQ;

J(4,1)=-IQ;
J(4,2)=-ID;
J(4,3)=-VtQ;
J(4,4)=-VtD;

update=-J\val;
VtQ=VtQ+update(1);
VtD=VtD+update(2);
IQ=IQ+update(3);
ID=ID+update(4);

val(1,1)=VtQ-zr*IQ+zi*ID-h1*EbQ+h2*EbD;
val(2,1)=VtD-zi*IQ-zr*ID-h1*EbD-h2*EbQ;
val(3,1)=qg-VtD*IQ+VtQ*ID;
val(4,1)=pg-VtD*ID-VtQ*IQ;

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
vgc=VtQ+sqrt(-1)*VtD;
vg=abs(vgc);
thetao=angle(vgc);
qg=VtD*IQ-VtQ*ID;


V1=vg;
V2=abs(h12*ebc);
%theta2=angle(h12*ebc);
%thetao=acos((V1*V1*zr/(zr*zr+zi*zi)-pg)*sqrt(zr*zr+zi*zi)/(V1*V2))+theta2-atan2(zi,zr);
vtoc=vg*(cos(thetao)+sqrt(-1)*sin(thetao));
iaoc=(vtoc-h12*ebc)/z11;
iao=abs(iaoc);
fio=angle(iaoc);
eq1=vtoc+(ra+sqrt(-1)*xq)*iaoc;
eqo=abs(eq1);
deltao=angle(eq1);
ido=-iao*sin(deltao-fio);
iqo=iao*cos(deltao-fio);
vdo=-vg*sin(deltao-thetao);
vqo=vg*cos(deltao-thetao);
efdo=eqo-(xd-xq)*ido;
eqo1=efdo+(xd-xd1)*ido;
edo1=-(xq-xq1)*iqo;
tmo=eqo1*iqo+edo1*ido+(xd1-xq1)*ido*iqo;
smo=0.0;
po=tmo;
vs=0;


if avr_chk_2==0
   if avr_chk_1~=0
	 vrefo=vg+(1/ka1)*(ke_exsys*efdo+a_exsys*exp(b_exsys*efdo));% for ieee type-1 excitation system.
   end
end

if avr_chk_1==0
   if avr_chk_2~=0
   vrefo=vg+(efdo/ka2);% for static excitation system
   end
end

if avr_chk_2==0
   if avr_chk_1==0
   vrefo=0;
   end
end
```


### Data File: steamgov.dat

Here is the content of the `steamgov.dat` file:

```plaintext
% .3
% 8
% .4
% .3
% .3
% .4
% .3
% 0
% .1
% .05
% .1
% -.1
% 0.9
% .1


.3
8
.4
.3
.3
.4
.3
0
.1
.05
.1
-.1
1.2
.1
```