```c title="AC4A_Exciter.c"


%----------------------Alternator supplied Controlled Rectifier: AC4A type-

load exc_AC4A.dat

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[mm order_AC4A]=sort([setxor(1:nb,exc_AC4A(:,1));exc_AC4A(:,1)]);
Tr_AC4A=exc_AC4A(:,2)';
kA_AC4A=exc_AC4A(:,3)';
TA_AC4A=exc_AC4A(:,4)';
TC_AC4A=exc_AC4A(:,5)';
TB_AC4A=exc_AC4A(:,6)';
VImax_AC4A=exc_AC4A(:,7)';
VImin_AC4A=exc_AC4A(:,8)';

VRmin_AC4A=exc_AC4A(:,9);
VRmax_AC4A=exc_AC4A(:,10);
KC_AC4A=exc_AC4A(:,11);

Efd0_AC4A=EFD0(exc_AC4A(:,1));

Vref_AC4A=Efd0_AC4A./kA_AC4A + lfl(exc_AC4A(:,1),2)';
Tgr0_AC4A=(Vref_AC4A - lfl(exc_AC4A(:,1),2)').*(1-TC_AC4A./TB_AC4A);
```


```c title="B_bus_From.c"
%---------------ybus--------------------------
Y=sparse(zeros(nb,nb));
Yd=sparse(zeros(nb,nb));

for i=1:nline
incr=1/(nt(i,3)+j*nt(i,4));
Y((nt(i,1)),(nt(i,2)))=Y((nt(i,1)),(nt(i,2)))-incr;
Y((nt(i,2)),(nt(i,1)))=Y((nt(i,2)),(nt(i,1)))-incr;
Y((nt(i,1)),(nt(i,1)))=Y((nt(i,1)),(nt(i,1)))+incr+j*nt(i,5)/2;
Y((nt(i,2)),(nt(i,2)))=Y((nt(i,2)),(nt(i,2)))+incr+j*nt(i,5)/2;
end

for i=nline+1:(ntrans+nline)
incr=1/(nt(i,3)+j*nt(i,4));  
incr1=incr/nt(i,5);
incr2=(1-nt(i,5))*incr/(nt(i,5)*nt(i,5));
incr3=(nt(i,5)-1)*incr/nt(i,5);
Y((nt(i,1)),(nt(i,2)))=Y((nt(i,1)),(nt(i,2)))-incr1;
Y((nt(i,2)),(nt(i,1)))=Y((nt(i,2)),(nt(i,1)))-incr1;
Y((nt(i,1)),(nt(i,1)))=Y((nt(i,1)),(nt(i,1)))+incr1+incr2;
Y((nt(i,2)),(nt(i,2)))=Y((nt(i,2)),(nt(i,2)))+incr1+incr3; 
end

for i=1:nshunt
incr=shunt(i,2)+j*shunt(i,3);   
Y((shunt(i,1)),(shunt(i,1)))=Y((shunt(i,1)),(shunt(i,1)))+incr;
end

%------Preparation of Y matrix to obtain Bd matrix----
for i=1:(nline + ntrans) 
   incr=1/(j*nt(i,4));
   Yd((nt(i,1)),(nt(i,2)))=Yd((nt(i,1)),(nt(i,2)))-incr;
   Yd((nt(i,2)),(nt(i,1)))=Yd((nt(i,2)),(nt(i,1)))-incr;
   Yd((nt(i,1)),(nt(i,1)))=Yd((nt(i,1)),(nt(i,1)))+incr;
   Yd((nt(i,2)),(nt(i,2)))=Yd((nt(i,2)),(nt(i,2)))+incr;
end
  Bd = -imag(Yd);
  Bd(sl,:) = [];
  Bd(:,sl) = [];
%-------------------------------------------------------------------------%
%------------------------------formation of Bdd matrix--------------------%      
  Bdd_bus = -imag(Y);
%-------------------------------------------------------------------------%
```


```data title="busno.dat"
100		%Slack bus number.
0.001  	%Load flow convergence tolerance.
145  		%No. of buses.
401		%No. of Lines.
52		%No. of transformers.
49		%No. of PV buses.
0		%Q-bit (please set this bit to zero only).
64		%No. of loads. 
97		%No. of shunts.
1.014		%Slack bus voltage magnitude.
50          %Nominal system frequency in Hz.
```


```c title="DC1A_Exciter.c"
%------------------------IEEE type-1 Exciters-------------------
load exc_DC1A.dat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mm order_DC1A]=sort([setxor(1:nb,exc_DC1A(:,1)) exc_DC1A(:,1)]);
Tr_DC1A=exc_DC1A(:,2)';
KA_DC1A=exc_DC1A(:,3)';
TA_DC1A=exc_DC1A(:,4)';
TC_DC1A=exc_DC1A(:,5)';
TB_DC1A=exc_DC1A(:,6)';
VRmax_DC1A=exc_DC1A(:,7)';
VRmin_DC1A=exc_DC1A(:,8)';
KE_DC1A=exc_DC1A(:,9)';
TE_DC1A=exc_DC1A(:,10)';
E1=exc_DC1A(:,11)';
SE1=exc_DC1A(:,12)';
E2=exc_DC1A(:,13)';
SE2=exc_DC1A(:,14)';
KF_DC1A=exc_DC1A(:,15)';
TF_DC1A=exc_DC1A(:,16)';

Efd0_DC1A=EFD0(exc_DC1A(:,1));

%calculation of saturation function coef., A and B in  SE=A*exp(B*Efd)
B_DC1A=(1./(E1-E2)).*log(SE1./SE2);
A_DC1A=SE1./(exp(B_DC1A.*E1));

VR0_DC1A=(A_DC1A.*exp(B_DC1A.*Efd0_DC1A)+KE_DC1A).*Efd0_DC1A;
Vref_DC1A=VR0_DC1A./KA_DC1A + lfl(exc_DC1A(:,1),2)';
Fst0_DC1A=Efd0_DC1A.*KF_DC1A;
Tgr0_DC1A=(Vref_DC1A - lfl(exc_DC1A(:,1),2)').*(1-TC_DC1A./TB_DC1A);
```


```data title="delPw_pss.dat"
60 10 10 10 10 0.01 10 54 1 0 0.1 0.1 0.05 0.1 0.05 10 0.1 -0.1
```

```c title="delPw_pss_settings.c"

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        FORMATION OF A MATRIX with PSS

TAGmod=sparse([zeros(ngen) zeros(ngen) zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag21)  zeros(ngen)  diag(Ag23)    diag(Ag24)   diag(Ag25)   diag(Ag26)    zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)]);    
              
 TAG=TAGmod(ii,ii);
 clear TAGmod;
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 TBGmod=sparse([zeros(ngen)      zeros(ngen)
              diag(Bg21rect)   diag(Bg22rect)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen) ]);
           
           
  TBG=TBGmod(ii,kk);
  clear TBGmod; 
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 if(npss~=0)
   TAGp=([TAG+sparse(zeros(size(EG*PS1'*DPS*PS1*FG))) sparse(zeros(size(EG*PS1'*CPS))); BPS*PS1*FG APS]);
   TBGp=([TBG; sparse(zeros(size(APS,1),2*ngen))]);
   TCGp=CG;
   TEGp=([EG; sparse(zeros(size(APS,1),ngen))]); %padded zeros to EG suitable.
   TFGp=([FG  sparse(zeros(ngen,size(APS,1)))]); %padded zeros to FG suitable.
   TEATd=-TAGp-TBGp*PG'*YDQdi*PG*TCGp;
else
   TEGp=EG;
   TFGp=FG; 
  TEATd=-TAG-TBG*PG'*YDQdi*PG*CG;
 end
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 H2_diag = 2*diag(H);
 TEATD=H2_diag*TEATd([2:16:16*ngen],:);

APSS=[];
BPSS=[];
CPSS=[];
DPSS=[];

for k=1:size(delPw_pss,1)
[N1,D1]=series([Tw1_delPw(k) 0],[Tw1_delPw(k) 1],[Tw2_delPw(k) 0],[Tw2_delPw(k) 1]);
[N1,D1]=series(N1,D1,[0 1],[T6_delPw(k) 1]);
[APW,BPW,CPW,DPW]=tf2ss(N1,D1);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[N2,D2]=series([Tw3_delPw(k) 0],[Tw3_delPw(k) 1],[Tw4_delPw(k) 0],[Tw4_delPw(k) 1]);
[N2,D2]=series(N2,D2,[0 T7_delPw(k)/2/H_delPw(k)],[T7_delPw(k) 1]);
[APT,BPT,CPT,DPT]=tf2ss(N2,D2);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[N3,D3]=series([T8_delPw(k) 1],[(T9_delPw(k)*T9_delPw(k)) (2*T9_delPw(k)) 1],[T8_delPw(k) 1],[(T9_delPw(k)*T9_delPw(k)) (2*T9_delPw(k)) 1]);
[N3,D3]=series(N3,D3,[T8_delPw(k) 1],[(T9_delPw(k)*T9_delPw(k)) (2*T9_delPw(k)) 1]);
[N3,D3]=series(N3,D3,[T8_delPw(k) 1],[(T9_delPw(k)*T9_delPw(k)) (2*T9_delPw(k)) 1]);
[APF,BPF,CPF,DPF]=tf2ss(N3,D3);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[N4,D4]=series([T1_delPw(k) 1],[T2_delPw(k) 1],[Ks1_delPw(k)],[1]);
[N4,D4]=series(N4,D4,[T3_delPw(k) 1],[T4_delPw(k) 1]);
[APC,BPC,CPC,DPC]=tf2ss(N4,D4);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
APSSd=[APW    zeros(size(APW,1),size(APT,1)) zeros(size(APW,1),size(APF,1))  zeros(size(APW,1),size(APC,1))
    zeros(size(APT,1),size(APW,1))    APT    zeros(size(APT,1),size(APF,1))  zeros(size(APT,1),size(APC,1))
    (BPF*CPW)      (BPF*Ks3_delPw(k)*CPT)                APF                    zeros(size(APF,1),size(APC,1))
    (BPC*DPF*CPW)  (BPC*(DPF*Ks3_delPw(k)*CPT-CPT))    (BPC*CPF)                APC];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
BPSSd=[     BPW         zeros(size(APW,1),1)
       zeros(size(APT,1),1)     BPT
       (BPF*DPW)       (BPF*Ks3_delPw(k)*DPT)
       (BPC*DPF*DPW)   (BPC*(DPF*Ks3_delPw(k)*DPT-DPT))];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
CPSSd=[DPC*DPF*CPW   (DPC*(DPF*Ks3_delPw(k)*CPT-CPT))  DPC*CPF    CPC];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DPSSd=[DPC*DPF*DPW   (DPC*(DPF*Ks3_delPw(k)*DPT-DPT))];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
[APSS,BPSS,CPSS,DPSS]=append(APSS,BPSS,CPSS,DPSS,APSSd,BPSSd,CPSSd,DPSSd);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
APSS=sparse(APSS);
BPSS=sparse(BPSS);
CPSS=sparse(CPSS);
DPSS=sparse(DPSS);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for k=1:size(delPw_pss,1)
    for p=1:ngen
        if (delPw_pss(k,1)==gen(p,1))
            temp(k)=p;
        end
    end
end

PS2=[];
PS2m=[];
for q=1:size(delPw_pss,1)
PS21=zeros(1,ngen);
PS21(1,temp(q))=1;
PS2m1=[PS21 zeros(1,ngen);zeros(1,ngen) PS21];
PS2=[PS2;PS21];
PS2m=[PS2m;PS2m1];
end

FGd=[TFGp;TEATD];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AG=([AG+TEGp*PS2'*DPSS*PS2m*FGd TEGp*PS2'*CPSS; BPSS*PS2m*FGd APSS]);
BG=([BG; sparse(zeros(size(APSS,1),2*ngen))]);
CG=([CG  sparse(zeros(size(APSS,1),2*ngen)')]);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear C1 D1;
clear A1 B1;
clear temp;
```


```data title="exc_AC4A.dat"
 93 0.02 185.00 0.020 1.0 1.0 1000 -1000 -2.00 8.89 0 
104 0.02 253.00 0.015 1.0 1.0 1000 -1000 -7.00 8.86 0
105 0.02  54.63 0.468 1.0 1.0 1000 -1000  0.00 7.38 0
106 0.02  54.63 0.468 1.0 1.0 1000 -1000  0.00 7.38 0
110 0.02 185.00 0.020 1.0 1.0 1000 -1000 -2.00 8.89 0
111 0.02 253.00 0.015 1.0 1.0 1000 -1000 -7.00 8.86 0
```


```data title="exc_DC1A.dat"
60  0.02 200 0.02 1 10 6.0  -6.0  -0.0485 0.250 3.5461 0.0800 4.7281 0.260 0.0400 1.000
```


```data title="exc_static.dat"
60 200 0.02 -6 6
```


```c title="exciter_setting.c"

%%%%%%%%%%%%%%%%%%%% Exciter settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


AVR_dash=~(AVR);
AVR_sm=AVR_dash(gen(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
static_sm=zeros(1,nb);
static_sm(exc_static(:,1))=-1./TA_static;
static_sm=static_sm(gen(:,1));

SL_static_dash=~(SL_static);
SL_static_sm=SL_static_dash(gen(:,1));
static_dEfd=SL_static_sm.*static_sm.*AVR_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static_bm=zeros(1,nb);
static_bm(exc_static(:,1))=kA_static./TA_static;
static_bm=static_bm(gen(:,1));

static_bdEfd=SL_static_sm.*static_bm.*AVR_sm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%DC1A_EXCITER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Ag77%%%%%%%%%%%%%%%%5555
DC1A_Ag77=zeros(1,nb);
DC1A_Ag77(exc_DC1A(:,1))=(-(KE_DC1A+A_DC1A.*exp(B_DC1A.*Efd0_DC1A).*(B_DC1A.*Efd0_DC1A+1)))./TE_DC1A;
DC1A_Ag77=DC1A_Ag77(gen(:,1));

SL_DC1A_dash=~(SL_DC1A);
SL_DC1A_sm=SL_DC1A_dash(gen(:,1));

DC1A_dAg77=SL_DC1A_sm.*DC1A_Ag77.*AVR_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag78~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DC1A_Ag78=zeros(1,nb);
DC1A_Ag78(exc_DC1A(:,1))=1./TE_DC1A;
DC1A_Ag78=DC1A_Ag78(gen(:,1));

DC1A_dAg78=SL_DC1A_sm.*DC1A_Ag78.*AVR_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag87~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DC1A_Ag87=zeros(1,nb);
DC1A_Ag87(exc_DC1A(:,1))=-(KA_DC1A.*TC_DC1A.*KF_DC1A)./(TA_DC1A.*TB_DC1A.*TF_DC1A);
DC1A_Ag87=DC1A_Ag87(gen(:,1));
DC1A_dAg87=SL_DC1A_sm.*DC1A_Ag87.*AVR_sm;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag88~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DC1A_Ag88=zeros(1,nb);
DC1A_Ag88(exc_DC1A(:,1))=-1./TA_DC1A;
DC1A_Ag88=DC1A_Ag88(gen(:,1));
DC1A_dAg88=SL_DC1A_sm.*DC1A_Ag88.*AVR_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag89~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DC1A_Ag89=zeros(1,nb);
DC1A_Ag89(exc_DC1A(:,1))=KA_DC1A./TA_DC1A;
DC1A_Ag89=DC1A_Ag89(gen(:,1));
DC1A_dAg89=SL_DC1A_sm.*DC1A_Ag89.*AVR_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag810~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DC1A_Ag810=zeros(1,nb);
DC1A_Ag810(exc_DC1A(:,1))=(KA_DC1A.*TC_DC1A)./(TA_DC1A.*TB_DC1A);
DC1A_Ag810=DC1A_Ag810(gen(:,1));
DC1A_dAg810=SL_DC1A_sm.*DC1A_Ag810.*AVR_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag97~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DC1A_Ag97=zeros(1,nb);
DC1A_Ag97(exc_DC1A(:,1))=-((1-TC_DC1A./TB_DC1A).*(KF_DC1A./TF_DC1A))./TB_DC1A;
DC1A_Ag97=DC1A_Ag97(gen(:,1));
DC1A_dAg97=SL_DC1A_sm.*DC1A_Ag97.*AVR_sm;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag99~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DC1A_Ag99=zeros(1,nb);
DC1A_Ag99(exc_DC1A(:,1))=-1./TB_DC1A;
DC1A_Ag99=DC1A_Ag99(gen(:,1));
DC1A_dAg99=SL_DC1A_sm.*DC1A_Ag99.*AVR_sm;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag910~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DC1A_Ag910=zeros(1,nb);
DC1A_Ag910(exc_DC1A(:,1))=(1-TC_DC1A./TB_DC1A)./TB_DC1A;
DC1A_Ag910=DC1A_Ag910(gen(:,1));
DC1A_dAg910=SL_DC1A_sm.*DC1A_Ag910.*AVR_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`Ag107~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DC1A_Ag107=zeros(1,nb);
DC1A_Ag107(exc_DC1A(:,1))=(KF_DC1A./TF_DC1A)./TF_DC1A;
DC1A_Ag107=DC1A_Ag107(gen(:,1));
DC1A_dAg107=SL_DC1A_sm.*DC1A_Ag107.*AVR_sm;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag1010~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DC1A_Ag1010=zeros(1,nb);
DC1A_Ag1010(exc_DC1A(:,1))=-1./TF_DC1A;
DC1A_Ag1010=DC1A_Ag1010(gen(:,1));
DC1A_dAg1010=SL_DC1A_sm.*DC1A_Ag1010.*AVR_sm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Bg82~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DC1A_Bg82=zeros(1,nb);
DC1A_Bg82(exc_DC1A(:,1))=-(KA_DC1A.*TC_DC1A)./(TA_DC1A.*TB_DC1A);
DC1A_Bg82=DC1A_Bg82(gen(:,1));
DC1A_dBg82=SL_DC1A_sm.*DC1A_Bg82.*AVR_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Bg92~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DC1A_Bg92=zeros(1,nb);
DC1A_Bg92(exc_DC1A(:,1))=-(1-TC_DC1A./TB_DC1A)./TB_DC1A;
DC1A_Bg92=DC1A_Bg92(gen(:,1));
DC1A_dBg92=SL_DC1A_sm.*DC1A_Bg92.*AVR_sm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Eg8vr~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DC1A_Eg8vr=zeros(1,nb);
DC1A_Eg8vr(exc_DC1A(:,1))=(KA_DC1A.*TC_DC1A)./(TA_DC1A.*TB_DC1A);
DC1A_Eg8vr=DC1A_Eg8vr(gen(:,1));
DC1A_dEg8vr=SL_DC1A_sm.*DC1A_Eg8vr.*AVR_sm;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Eg9xb~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DC1A_Eg9xb=zeros(1,nb);
DC1A_Eg9xb(exc_DC1A(:,1))=(1-TC_DC1A./TB_DC1A)./TB_DC1A;
DC1A_Eg9xb=DC1A_Eg9xb(gen(:,1));
DC1A_dEg9xb=SL_DC1A_sm.*DC1A_Eg9xb.*AVR_sm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AC4A EXCITER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag77~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AC4A_Ag77=zeros(1,nb);
AC4A_Ag77(exc_AC4A(:,1))=-1./TA_AC4A;
AC4A_Ag77=AC4A_Ag77(gen(:,1));

SL_AC4A_dash=~(SL_AC4A);
SL_AC4A_sm=SL_AC4A_dash(gen(:,1));

AC4A_dAg77=SL_AC4A_sm.*AC4A_Ag77.*AVR_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag79~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AC4A_Ag79=zeros(1,nb);
AC4A_Ag79(exc_AC4A(:,1))=kA_AC4A./TA_AC4A;
AC4A_Ag79=AC4A_Ag79(gen(:,1));
AC4A_dAg79=SL_AC4A_sm.*AC4A_Ag79.*AVR_sm;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Ag99~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AC4A_Ag99=zeros(1,nb);
AC4A_Ag99(exc_AC4A(:,1))=-1./TB_AC4A;
AC4A_Ag99=AC4A_Ag99(gen(:,1));
AC4A_dAg99=SL_AC4A_sm.*AC4A_Ag99.*AVR_sm;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Bg72~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AC4A_Bg72=zeros(1,nb);
AC4A_Bg72(exc_AC4A(:,1))=-(kA_AC4A.*TC_AC4A)./(TA_AC4A.*TB_AC4A);
AC4A_Bg72=AC4A_Bg72(gen(:,1));
AC4A_dBg72=SL_AC4A_sm.*AC4A_Bg72.*AVR_sm;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Bg92~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AC4A_Bg92=zeros(1,nb);
AC4A_Bg92(exc_AC4A(:,1))=-(1-TC_AC4A./TB_AC4A)./TB_AC4A;
AC4A_Bg92=AC4A_Bg92(gen(:,1));
AC4A_dBg92=SL_AC4A_sm.*AC4A_Bg92.*AVR_sm;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Eg7~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AC4A_Eg7Efd=zeros(1,nb);
AC4A_Eg7Efd(exc_AC4A(:,1))=(kA_AC4A.*TC_AC4A)./(TA_AC4A.*TB_AC4A);
AC4A_Eg7Efd=AC4A_Eg7Efd(gen(:,1));
AC4A_dEg7Efd=SL_AC4A_sm.*AC4A_Eg7Efd.*AVR_sm;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Eg9~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AC4A_Eg9xb=zeros(1,nb);
AC4A_Eg9xb(exc_AC4A(:,1))=(1-TC_AC4A./TB_AC4A)./TB_AC4A;
AC4A_Eg9xb=AC4A_Eg9xb(gen(:,1));
AC4A_dEg9xb=SL_AC4A_sm.*AC4A_Eg9xb.*AVR_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```


```c title="fdlf_jacob_form.c"
%~~~~~~~~~~~~~Form Jacobian and solve for bus voltages and angles~~~~~~~~~~%
Vmag=sparse(ones(nb,1));
Vang=sparse(zeros(nb,1));
Vmag(pv_data(:,1))=pv_data(:,2);  %Replaces the voltages of P-V buses in Vmag 
Vmag(sl) = Vsl;
Vang(sl) = 0;
kp = 0;
kq = 0;
fid=fopen('report.dat', 'w');
%------------------Iteration begins from here onwards---------------------%
for k=1:200
   Vbus=Vmag.*(cos(Vang)+ j*sin(Vang));
%----------------------------------bus powers calculations----------------%
   S=Vbus.*(conj(Y*Vbus));
   Pc=real(S);
   Qc=imag(S);
%---------------------finding non slack and non pv buses------------------%
   pv_num=pv_data(:,1);
   pv_sl_num=[pv_num;sl];
   num_no_sl_pv =[1:nb]'; 
   num_no_sl_pv(pv_sl_num,:)=[];
 %------------------------bus power mis matches -------------------------%
   delP=Psp-Pc;
   delP(sl,:)=[];
   delQ=Qsp-Qc;
   delQ(pv_sl_num,:)=[];
   if (max(abs(delP))<= tole)
     if (max(abs(delQ))<= tole)
      fprintf(fid,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
      fprintf(fid,'Converged Real power mismatch iteration No = %i, Max.mismatch = %10.8f\n',kp, full(max(abs(delP))));
      fprintf(fid,'Converged Reactive power mismatch iteration No = %i, Max. mismatch = %10.8f\n',kq, full(max(abs(delQ))));
      convergence_bit=1;
      S(sl)= S(sl)+(sl_load(1,2)+j*sl_load(1,3));
      fprintf(fid,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
      break;
     end
   else
     kp=kp+1;
     if(k~=1)
      del_Pm=sparse(zeros(nb,1)); 
      del_Pm(sort([pv_num;num_no_sl_pv]))=delP;
      Bno=[1:nb]';
      vio_busP = Bno(abs(del_Pm)==max(abs(del_Pm)));
      fprintf(fid,'Iter. No = %i Max. Real power mismatch at bus = %i, Max. mismatch = %10.8f\n',(kp-1),full(vio_busP),full(max(abs(del_Pm))));
     end
  
    %---------------------Updation of bus angles at each bus----------------% 
     del_Vang = Bd\(delP./Vmag(num_nosl));
   
     Vang(num_nosl)=Vang(num_nosl)+ del_Vang;
  end
  %--------------------Updation of bus voltages at each bus-----------------% 
  Vbus = Vmag.*(cos(Vang)+ j*sin(Vang));
  S=Vbus.*(conj(Y*Vbus));
  Qc = imag(S);
  delQ=Qsp-Qc;
  delQ(pv_sl_num,:)=[];
     
  %-----------Perform the following if Q-limit is accounted--------------%
  pv_num_Uvio=sparse(zeros(nb,1));
  pv_num_Lvio=sparse(zeros(nb,1));
     
  if (Q_bit~=0 & max(abs(delQ))<=0.1)
 
    Qc(pvpq_buses) = Qc(pvpq_buses) + Qload(pvpq_buses); 
       
  %----------------------checking Qlimit voilation-----------------------%
    if (sum(Qc(pv_num)>QUlim(pv_num))>=1)
     pv_num_Uvio(pv_num)=[Qc(pv_num)>QUlim(pv_num)];
        
     pv_num=[pv_num(Qc(pv_num)<QUlim(pv_num))];
     pv_sl_num=[pv_num;sl];
     num_no_sl_pv = [1:nb]';
     num_no_sl_pv(pv_sl_num,:)=[];
    end
        
    if (sum(Qc(pv_num)<QLlim(pv_num))>=1)
     pv_num_Lvio(pv_num) =[Qc(pv_num)<QLlim(pv_num)];
          
     pv_num=[pv_num(Qc(pv_num)>QLlim(pv_num))];
     pv_sl_num=[pv_num;sl];
     num_no_sl_pv = [1:nb]';
     num_no_sl_pv(pv_sl_num,:)=[];
    end
        
    Qsp = sparse(zeros(nb,1));
    Qsp(pq_data(:,1)) = -pq_data(:,3);
    Qsp = Qsp + (pv_num_Uvio.*QUlim)+ (pv_num_Lvio.*QLlim);
         
    Vmag(pv_num)=pv_Vmag(pv_num);  %Again Replace the voltages of P-V buses in Vmag 
       
   %------------------------Re-calculation of bus powers --------------------%
    Vbus=Vmag.*(cos(Vang)+ j*sin(Vang));
    S=Vbus.*(conj(Y*Vbus));
    Qc = imag(S);
  end       %----------------Q-limit part ends here----------------------%
      
  delP=Psp-real(S);
  delP(sl,:)=[];
  delQ=Qsp-Qc;
  delQ(pv_sl_num,:)=[];
   
  if (max(abs(delQ))<= tole)
   if (max(abs(delP))<= tole)
    fprintf(fid,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    fprintf(fid,'Converged Real power mismatch iteration No = %i, Max.mismatch = %10.8f\n',kp, full(max(abs(delP))));
    fprintf(fid,'Converged Reactive power mismatch iteration No = %i, Max. mismatch = %10.8f\n',kq, full(max(abs(delQ))));
    convergence_bit=1;
    S(sl)= S(sl)+(sl_load(1,2)+j*sl_load(1,3));
    fprintf(fid,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    break;
   end
  else
   kq = kq+1; 
   if(k~=1)
    del_Qm=sparse(zeros(nb,1));
    del_Qm(num_no_sl_pv)=delQ;
    Bno=[1:nb]';
    vio_busQ=Bno(abs(del_Qm)==max(abs(del_Qm)));
    fprintf(fid,'Iter. No = %i Max. Reactive power mismatch at bus = %i, Max. mismatch = %10.8f\n',(kq-1),full(vio_busQ),full(max(abs(del_Qm))));
   end 
        
   Bdd = Bdd_bus;
   Bdd(:,pv_sl_num)=[];
   Bdd(pv_sl_num,:)=[];
   
   del_Vmag = Bdd\(delQ./Vmag(num_no_sl_pv));
   Vmag(num_no_sl_pv)=Vmag(num_no_sl_pv) + del_Vmag;
     
   if (sum(pv_num_Uvio + pv_num_Lvio)~=0)
    pv_Vmag(pv_data(:,1)) = Vmag(pv_data(:,1));
   end
  end
end

fclose(fid);
```



```c title="fdlf_loadflow.c"
clear all;
t0=clock;
convergence_bit=0;    % this bit will be set to 1 if convergence is reached.
%----------------------------------Load the system data--------------------------%
load busno.dat;
load nt.dat;
load pvpq.dat;
%--------------------------------------------------------------------------------%
sl=busno(1);
tole=busno(2);
nb=busno(3);
nline=busno(4);
ntrans=busno(5);
npv=busno(6);
Q_bit=busno(7);
npq=busno(8);                                                                                                                                                                          
nshunt=busno(9);
Vsl=busno(10);
if (nshunt~=0)
  load shunt.dat;
end   
if (Q_bit~=0)
  load Qlim_data.dat;
end 

pv_data=pvpq(1:npv,:);      %P-V bus data
load_data=pvpq(npv+1:npv+npq,:); % load data
pq_data = load_data((load_data(:,1)~=sl),:);  %P-Q bus data
pv=sparse(zeros(nb,1));
pq=sparse(zeros(nb,1));
pv(pv_data(:,1))=pv_data(:,1);
pq(pq_data(:,1))=pq_data(:,1);
pvpq_buses=pq(pv==pq);   
pvpq_buses((pvpq_buses==0),:)=[]; %pvpq bus nos


sl_load=sparse(zeros(1,3));
if(sum(pvpq(:,1)==sl)==1)
  sl_load=pvpq((pvpq(:,1)==sl),:);   %slack bus load data
end


num_nosl = [1:nb]';
num_nosl(sl,:) = [];                                                      
%------------------------call B_bus_form.m  to construct the Bd and Bdd matrices-----%
B_bus_form;       
%--------------------------------------------------------------------------%
Psp = sparse(zeros(nb,1));
Psp(pv_data(:,1)) = pv_data(:,3);
Psp(pq_data(:,1)) = Psp(pq_data(:,1))-pq_data(:,2);
Qsp = sparse(zeros(nb,1));
Qsp(pq_data(:,1)) = -pq_data(:,3);

pv_Vmag =sparse(zeros(nb,1));
pv_Vmag(pv_data(:,1))= pv_data(:,2);

Qload =sparse(zeros(nb,1));
Qload(pq_data(:,1))= pq_data(:,3);

if (Q_bit~=0)
  QLlim = sparse(zeros(nb,1)); 
  QUlim = sparse(zeros(nb,1));
  QLlim(pv_data(:,1))=-9999;
  QUlim(pv_data(:,1))=9999;
  QUlim(Qlim_data(:,1)) = Qlim_data(:,3); 
  QLlim(Qlim_data(:,1)) = Qlim_data(:,2);          
end
  
%-----------------call fdlf_jacob_form.m to obtain the bus voltages----------%
  
  fdlf_jacob_form;
  
%-------------------perform the power flow calculations----------------%
  
  powerflow;     
  
%-------------------prepare the lfl.dat and report.dat--------------------------%
  lfl_result
  
if (convergence_bit==1)
 fprintf('\n\n Please wait... \n\n')    
 disp('Solution Converged !!!' );
else
 disp('Convergence is not reached !!!' );
end
 
elapsed_time = etime(clock,t0)
```


```c title="freq_responce.c"
%frequency_response;
TAGmod=sparse([zeros(ngen) zeros(ngen) zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag21)  zeros(ngen)  diag(Ag23)    diag(Ag24)   diag(Ag25)   diag(Ag26)    zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)]);    
              
 TAG=TAGmod(ii,ii);
 clear TAGmod;
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 TBGmod=sparse([zeros(ngen)      zeros(ngen)
              diag(Bg21rect)   diag(Bg22rect)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen) ]);
           
           
  TBG=TBGmod(ii,kk);
  clear TBGmod; 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(npss~=0)
   TAGp=([TAG+sparse(zeros(size(EG*PS1'*DPS*PS1*FG))) sparse(zeros(size(EG*PS1'*CPS))); BPS*PS1*FG APS]);
   TBGp=([TBG; sparse(zeros(size(APS,1),2*ngen))]);
   TEGp=([EG; sparse(zeros(size(APS,1),ngen))]);%padded zeros to EG suitable.
  else
  TAGp=TAG;
  TBGp=TBG;
  TEGp=EG;
end
 if(npss1~=0)        
    TAGp=([TAGp+sparse(zeros(size(TEGp*PS2'*DPSS*PS2m*FGd))) sparse(zeros(size(TEGp*PS2'*CPSS))); BPSS*PS2m*FGd APSS]); 
    TBGp=([TBGp; sparse(zeros(size(APSS,1),2*ngen))]);
    TEGp=([TEGp; sparse(zeros(size(APSS,1),ngen))]);%padded zeros to EG suitable.
 end
  if(npss2~=0)        
    TAGp=([TAGp+sparse(zeros(size(TEGp*PS3'*DPSS2*PS3*TEATD2))) sparse(zeros(size(TEGp*PS3'*CPSS2))); BPSS2*PS3*TEATD2 APSS2]); 
    TBGp=([TBGp; sparse(zeros(size(APSS2,1),2*ngen))]);
    TEGp=([TEGp; sparse(zeros(size(APSS2,1),ngen))]);%padded zeros to EG suitable.
 end
 
   if((npss+npss1+npss2)~=0)
     TA=Ap;% Ap represents the overall state matrix is obtained in the small_sig programme with PSS
     else
     TA=Anp;% Anp represents the overall state matrix is obtained in the small_sig programme without PSS
    end 
   TCGp=CG;
   TEATd=-TAGp-TBGp*PG'*YDQdi*PG*TCGp;
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TEAT=TEATd([2:16:16*ngen],:);
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mtpfreq=1;
  while(mtpfreq==1)
  rr=input('Enter the generator number for which you want to obtaine the frequency response:  ');
  for gg=1:1:ngen
     if(gen(gg,1)==rr)
        te=gg;
     end
  end
  TEATR=TEAT(te,:)*2*H(te);%taking the required machine Te data.
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  EGT=TEGp(:,te);% EG formed in the main programme is used.
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  [V2,Deig]=eig(full(TA));
  V2=sparse(V2);
  W2=sparse(V2\eye(size(V2)));
  for i=1:470
    w(i)=0.05*i;
    smDinv = 1./(j*w(i)-diag(Deig)+0.0001);
    smDinv = sparse(diag(smDinv));
    Xwd=(V2*smDinv*W2)*EGT;
    Te_Vs(i)=TEATR*Xwd;
   end
 delete(get(0,'children'))
 plot(w/pi/2,abs(Te_Vs))
 grid
 mtpfreq=input('Enter 1 to repeate frequency plot for another machine otherwise 0 :    ');
end
clear Xwd;
clear Te_Vs;
```


```data title="gen.dat"
 93 0.09842 0.02400 0.02400  8.50 0.05 0.09673 0.03655 0.03655 1.24 0.04000  115.0366 0.0 
104 0.10160 0.01220 0.01220 10.00 0.05 0.09820 0.01440 0.01440 1.50 0.04000   73.8528 0.0 
105 0.11440 0.02080 0.02080  6.61 0.05 0.10920 0.03149 0.03149 1.50 0.04000   84.3915 0.0 
106 0.17165 0.03118 0.03118  6.61 0.05 0.16377 0.04720 0.04720 1.50 0.04000   56.2610 0.0 
110 0.09842 0.02400 0.02400  8.50 0.05 0.09673 0.03655 0.03655 1.24 0.04000  115.0500 0.0 
111 0.10160 0.01220 0.01220 10.00 0.05 0.09820 0.01440 0.01440 1.50 0.04000   73.8528 0.0 
 60 2.86840 0.4769  0.4769  10000 0.05 0.4769  0.4769  0.4769  0.40 0.04000    1.4100 0.0 
 67 0.12780 0.0213  0.0213  10000 0.05 0.0213  0.0213  0.0213  0.40 0.04000   52.1796 0.0 
 79 0.77520 0.1292  0.1292  10000 0.05 0.1292  0.1292  0.1292  0.40 0.04000    6.6500 0.0 
 80 3.98880 0.6648  0.6648  10000 0.05 0.6648  0.6648  0.6648  0.40 0.04000    1.2857 0.0 
 82 3.17460 0.5291  0.5291  10000 0.05 0.5291  0.5291  0.5291  0.40 0.04000    2.1150 0.0 
 89 0.35100 0.0585  0.0585  10000 0.05 0.0585  0.0585  0.0585  0.40 0.04000   20.5602 0.0 
 90 9.60000 1.6000  1.6000  10000 0.05 1.6000  1.6000  1.6000  0.40 0.04000    0.7628 0.0 
 91 2.23080 0.3718  0.3718  10000 0.05 0.3718  0.3718  0.3718  0.40 0.04000    1.6848 0.0 
 94 0.50340 0.0839  0.0839  10000 0.05 0.0839  0.0839  0.0839  0.40 0.04000   17.3424 0.0 
 95 0.97140 0.1619  0.1619  10000 0.05 0.1619  0.1619  0.1619  0.40 0.04000    5.4662 0.0 
 96 2.89440 0.4824  0.4824  10000 0.05 0.4824  0.4824  0.4824  0.40 0.04000    2.1216 0.0 
 97 1.27500 0.2125  0.2125  10000 0.05 0.2125  0.2125  0.2125  0.40 0.04000    5.4912 0.0 
 98 1.70000 0.0795  0.0795  10000 0.05 0.0795  0.0795  0.0795  0.40 0.04000   13.9600 0.0 
 99 0.47700 0.1146  0.1146  10000 0.05 0.1146  0.1146  0.1146  0.40 0.04000   17.1080 0.0 
100 0.83760 0.1386  0.1386  10000 0.05 0.1386  0.1386  0.1386  0.40 0.04000    7.5600 0.0 
101 0.55440 0.0924  0.0924  10000 0.05 0.0924  0.0924  0.0924  0.40 0.04000   12.2844 0.0 
102 0.08100 0.0135  0.0135  10000 0.05 0.0135  0.0135  0.0135  0.40 0.04000   78.4366 0.0 
103 0.63780 0.1063  0.1063  10000 0.05 0.1063  0.1063  0.1063  0.40 0.04000    8.1600 0.0 
108 0.14880 0.0248  0.0248  10000 0.05 0.0248  0.0248  0.0248  0.40 0.04000   30.4320 0.0 
109 1.21740 0.2029  0.2029  10000 0.05 0.2029  0.2029  0.2029  0.40 0.04000    2.6622 0.0 
112 0.55440 0.0924  0.0924  10000 0.05 0.0924  0.0924  0.0924  0.40 0.04000   12.2844 0.0 
115 0.01440 0.0024  0.0024  10000 0.05 0.0024  0.0024  0.0024  0.40 0.04000   97.3300 0.0 
116 0.01320 0.0022  0.0022  10000 0.05 0.0022  0.0022  0.0022  0.40 0.04000  105.5000 0.0 
117 0.01020 0.0017  0.0017  10000 0.05 0.0017  0.0017  0.0017  0.40 0.04000  102.1600 0.0 
118 0.00840 0.0014  0.0014  10000 0.05 0.0014  0.0014  0.0014  0.40 0.04000  162.7400 0.0 
119 0.00120 0.0002  0.0002  10000 0.05 0.0002  0.0002  0.0002  0.40 0.04000  348.2200 0.0 
121 0.01020 0.0017  0.0017  10000 0.05 0.0017  0.0017  0.0017  0.40 0.04000  116.5400 0.0 
122 0.05340 0.0089  0.0089  10000 0.05 0.0089  0.0089  0.0089  0.40 0.04000   39.2400 0.0 
124 0.01020 0.0017  0.0017  10000 0.05 0.0017  0.0017  0.0017  0.40 0.04000  116.8600 0.0 
128 0.00060 0.0001  0.0001  10000 0.05 0.0001  0.0001  0.0001  0.40 0.04000  503.8700 0.0 
130 0.00600 0.0010  0.0010  10000 0.05 0.0010  0.0010  0.0010  0.40 0.04000  230.9000 0.0 
131 0.00060 0.0001  0.0001  10000 0.05 0.0001  0.0001  0.0001  0.40 0.04000 1101.7200 0.0 
132 0.00960 0.0016  0.0016  10000 0.05 0.0016  0.0016  0.0016  0.40 0.04000  120.3500 0.0 
134 0.00180 0.0003  0.0003  10000 0.05 0.0003  0.0003  0.0003  0.40 0.04000  802.1200 0.0 
135 0.00480 0.0008  0.0008  10000 0.05 0.0008  0.0008  0.0008  0.40 0.04000  232.6300 0.0 
136 0.00060 0.0001  0.0001  10000 0.05 0.0001  0.0001  0.0001  0.40 0.04000 2018.1700 0.0 
137 0.00240 0.0004  0.0004  10000 0.05 0.0004  0.0004  0.0004  0.40 0.04000  469.3200 0.0 
139 0.00060 0.0001  0.0001  10000 0.05 0.0001  0.0001  0.0001  0.40 0.04000 2210.2000 0.0 
140 0.00180 0.0003  0.0003  10000 0.05 0.0003  0.0003  0.0003  0.40 0.04000  899.1900 0.0 
141 0.00060 0.0001  0.0001  10000 0.05 0.0001  0.0001  0.0001  0.40 0.04000 1474.2200 0.0 
142 0.00180 0.0003  0.0003  10000 0.05 0.0003  0.0003  0.0003  0.40 0.04000  950.8000 0.0 
143 0.01380 0.0023  0.0023  10000 0.05 0.0023  0.0023  0.0023  0.40 0.04000  204.3000 0.0 
144 0.00240 0.0004  0.0004  10000 0.05 0.0004  0.0004  0.0004  0.40 0.04000  443.2200 0.0 
145 0.01080 0.0018  0.0018  10000 0.05 0.0018  0.0018  0.0018  0.40 0.04000  518.0800 0.0 
```


```c title="genamt.c"

% Formation of the generator matrices


%****************************************************************************************************
%****************************************************************************************************
%             FORMATION OF THE `A' MATRIX BEGINS
%****************************************************************************************************
%****************************************************************************************************

dmth0=delta0-thetag0;
C1=(xdd-xddd)./xdd;
C2=(xd-xdd).*xddd./xd./xdd;
C3=(xqd-xqdd)./xqd;
C4=(xq-xqd).*xqdd./xq./xqd;
C5=(xddd-xqdd)./xddd./xqdd;

Ag12=wB*ones(1,ngen);

Ag21=(C2.*Vg0.*sif0.*cos(dmth0)./xddd + C1.*Vg0.*sih0.*cos(dmth0)./xddd +C4.*Vg0.*sig0.*sin(dmth0)./xqdd);
Ag21=(Ag21+C3.*Vg0.*sik0.*sin(dmth0)./xqdd + C5.*Vg0.*Vg0.*cos(2.0*(dmth0)))*(-0.5)./H;
Ag22=-(Dm./H)*(0.5);

Ag23=-C2.*Vg0.*sin(dmth0)*(0.5)./H./xddd;

Ag24=-C1.*Vg0.*sin(dmth0)*(0.5)./H./xddd;

Ag25=C4.*Vg0.*cos(dmth0)*(0.5)./H./xqdd;

Ag26=C3.*Vg0.*cos(dmth0)*(0.5)./H./xqdd;

Ag31=-Vg0.*sin(dmth0)./Tdd;

Ag33=-ones(1,ngen)./Tdd;

Ag37=(xdd./(xd-xdd))./Tdd;

Ag41=-Vg0.*sin(dmth0)./Tddd;

Ag44=-ones(1,ngen)./Tddd;

Ag51=Vg0.*cos(dmth0)./Tqd;

Ag55=-ones(1,ngen)./Tqd;

Ag61=Vg0.*cos(dmth0)./Tqdd;

Ag66=-ones(1,ngen)./Tqdd;

Ag77=static_dEfd+DC1A_dAg77+AC4A_dAg77;

Ag78=DC1A_dAg78;

Ag79=AC4A_dAg79;

Ag87=DC1A_dAg87;

Ag88=DC1A_dAg88;

Ag89=DC1A_dAg89;

Ag810=DC1A_dAg810;

Ag97=DC1A_dAg97;

Ag99=DC1A_dAg99+AC4A_dAg99;

Ag910=DC1A_dAg910;

Ag107=DC1A_dAg107;

Ag1010=DC1A_dAg1010;

Ag211=RHST_dAg211;

Ag212=RHST_dAg212;

Ag213=RHST_dAg213;

Ag215=HYDRO_dAg215;

Ag216=HYDRO_dAg216;

Ag1111=RHST_dAg1111;

Ag1115=RHST_dAg1115;

Ag1211=RHST_dAg1211;

Ag1212=RHST_dAg1212;

Ag1312=RHST_dAg1312;

Ag1313=RHST_dAg1313;

Ag142=RHST_dAg142+HYDRO_dAg142;

Ag1414=RHST_dAg1414+HYDRO_dAg1414;

Ag152=RHST_dAg152+HYDRO_dAg152;

Ag1514=RHST_dAg1514+HYDRO_dAg1514;

Ag1515=RHST_dAg1515+HYDRO_dAg1515;

Ag1615=HYDRO_dAg1615;

Ag1616=HYDRO_dAg1616;




AGmod=sparse([zeros(ngen) diag(Ag12)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag21)  diag(Ag22)  diag(Ag23)    diag(Ag24)   diag(Ag25)   diag(Ag26)    zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  diag(Ag211)  diag(Ag212)  diag(Ag213)  zeros(ngen)   diag(Ag215)  diag(Ag216)
              diag(Ag31)  zeros(ngen) diag(Ag33)    zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag37)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag41)  zeros(ngen) zeros(ngen)   diag(Ag44)   zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag51)  zeros(ngen) zeros(ngen)   zeros(ngen)  diag(Ag55)   zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag61)  zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  diag(Ag66)    zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag77)  diag(Ag78)    diag(Ag79)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag87)  diag(Ag88)    diag(Ag89)   diag(Ag810)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag97)  zeros(ngen)   diag(Ag99)   diag(Ag910)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag107) zeros(ngen)   zeros(ngen)  diag(Ag1010) zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  diag(Ag1111) zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag1115) zeros(ngen)                                        
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  diag(Ag1211) diag(Ag1212) zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)   zeros(ngen) diag(Ag1312) diag(Ag1313) zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) diag(Ag142) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)  zeros(ngen)  diag(Ag1414)  zeros(ngen)  zeros(ngen)
              zeros(ngen) diag(Ag152) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)  zeros(ngen)  diag(Ag1514)  diag(Ag1515) zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag1615) diag(Ag1616)]); 
             
           
           
   i=[1:ngen
  (ngen+1):(2*ngen)
  (2*ngen+1):(3*ngen)
  (3*ngen+1):(4*ngen)       
  (4*ngen+1):(5*ngen)
  (5*ngen+1):(6*ngen)
  (6*ngen+1):(7*ngen)
  (7*ngen+1):(8*ngen)
  (8*ngen+1):(9*ngen)
  (9*ngen+1):(10*ngen)
  (10*ngen+1):(11*ngen)
  (11*ngen+1):(12*ngen)
  (12*ngen+1):(13*ngen)
  (13*ngen+1):(14*ngen)
  (14*ngen+1):(15*ngen)
  (15*ngen+1):(16*ngen)];

ii=i(:);

AG=AGmod(ii,ii);

clear AGmod;
%clear Ag21 Ag22 Ag23 Ag24 Ag25 Ag26
clear Ag12  Ag31 Ag33 Ag37 Ag41 Ag44  Ag51 Ag55 Ag61 Ag66 Ag77 Ag78 Ag87 Ag88 Ag89 Ag810 Ag97 Ag99 Ag910 Ag107 Ag1010;
clear Ag211 Ag212 Ag213 Ag214 Ag215 Ag216 Ag1111 Ag1115 Ag1211 Ag1212 Ag1312 Ag1313 Ag142 Ag1414 Ag152 Ag1514 Ag1515 Ag1615 Ag1616; 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bg21=(C2.*sif0.*cos(dmth0)./xddd + C1.*sih0.*cos(dmth0)./xddd +C4.*sig0.*sin(dmth0)./xqdd);

Bg21=(Bg21+C3.*sik0.*sin(dmth0)./xqdd + C5.*Vg0.*cos(2.0*(dmth0)))*(0.5)./H;

Bg22=(C2.*sif0.*sin(dmth0)./xddd + C1.*sih0.*sin(dmth0)./xddd -C4.*sig0.*cos(dmth0)./xqdd);

Bg22=(Bg22-C3.*sik0.*cos(dmth0)./xqdd + C5.*Vg0.*sin(2.0*(dmth0)))*(-0.5)./H;

Bg31=sin(dmth0)./Tdd;

Bg32=cos(dmth0)./Tdd;

Bg41=sin(dmth0)./Tddd;

Bg42=cos(dmth0)./Tddd;

Bg51=-cos(dmth0)./Tqd;

Bg52=sin(dmth0)./Tqd;

Bg61=-cos(dmth0)./Tqdd;

Bg62=sin(dmth0)./Tqdd;

Bg72=-static_bdEfd+AC4A_dBg72;

Bg82=DC1A_dBg82;

Bg92=DC1A_dBg92+AC4A_dBg92;

Bg21rect=(Bg21.*(-VDg0)+Bg22.*VQg0)./Vg0;
Bg22rect=(Bg21.*VQg0+Bg22.*VDg0)./Vg0;

Bg31rect=(Bg31.*(-VDg0)+Bg32.*VQg0)./Vg0;
Bg32rect=(Bg31.*VQg0+Bg32.*VDg0)./Vg0;

Bg41rect=(Bg41.*(-VDg0)+Bg42.*VQg0)./Vg0;
Bg42rect=(Bg41.*VQg0+Bg42.*VDg0)./Vg0;

Bg51rect=(Bg51.*(-VDg0)+Bg52.*VQg0)./Vg0;
Bg52rect=(Bg51.*VQg0+Bg52.*VDg0)./Vg0;

Bg61rect=(Bg61.*(-VDg0)+Bg62.*VQg0)./Vg0;
Bg62rect=(Bg61.*VQg0+Bg62.*VDg0)./Vg0;

Bg71rect=(Bg72.*VQg0)./Vg0;
Bg72rect=(Bg72.*VDg0)./Vg0;


Bg81rect=(Bg82.*VQg0)./Vg0;
Bg82rect=(Bg82.*VDg0)./Vg0;


Bg91rect=(Bg92.*VQg0)./Vg0;
Bg92rect=(Bg92.*VDg0)./Vg0;



BGmod=sparse([zeros(ngen)      zeros(ngen)
              diag(Bg21rect)   diag(Bg22rect)
              diag(Bg31rect)   diag(Bg32rect)
              diag(Bg41rect)   diag(Bg42rect)
              diag(Bg51rect)   diag(Bg52rect)
              diag(Bg61rect)   diag(Bg62rect)
              diag(Bg71rect)   diag(Bg72rect)
              diag(Bg81rect)   diag(Bg82rect)
              diag(Bg91rect)   diag(Bg92rect)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen) ]);
        
       

BG=BGmod(ii,kk);

clear BGmod;
%clear Bg21 Bg22
clear  Bg31 Bg32 Bg41 Bg42 Bg51 Bg52 Bg61 Bg62 Bg72 Bg82 Bg92;
%clear Bg21rect Bg22rect
clear  Bg31rect Bg32rect Bg41rect Bg42rect Bg51rect Bg52rect Bg61rect Bg62rect Bg72rect;
clear Bg81rect Bg82rect Bg91rect Bg92rect;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Input matrix used to interface PSS -see pssmat.m
EG=sparse(zeros(ngen*16,ngen));
EG([7:16:ngen*16],[1:ngen])=diag(static_bdEfd+AC4A_dEg7Efd);
EG([8:16:ngen*16],[1:ngen])=diag(DC1A_dEg8vr);
EG([9:16:ngen*16],[1:ngen])=diag(DC1A_dEg9xb+AC4A_dEg9xb);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%This matrix is useful to extract 'del s_m' from the 
%x=[deldelta,dels_m,delSif,delSih,delSig,delSik,delEfd,delvr,delxb,delxf,delx1,delx2,delx3,dely1,delPgv,delz1 ]'
FG=sparse(zeros(ngen,ngen*16));
FG([1:ngen],(16*[1:ngen]-14))=eye(ngen);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
Cg11=iQg0+Vg0.*(sin(delta0).*cos(dmth0)./xqdd-sin(dmth0).*cos(delta0)./xddd);    
 
Cg13=-C2.*cos(delta0)./xddd;

Cg14=-C1.*cos(delta0)./xddd;

Cg15=-C4.*sin(delta0)./xqdd;

Cg16=-C3.*sin(delta0)./xqdd;

Cg21=Vg0.*(sin(delta0).*sin(dmth0)./xddd+cos(delta0).*cos(dmth0)./xqdd)-iDg0; 
     
Cg23=C2.*sin(delta0)./xddd;

Cg24=C1.*sin(delta0)./xddd;

Cg25=-C4.*cos(delta0)./xqdd;

Cg26=-C3.*cos(delta0)./xqdd;


CGmod=sparse([diag(Cg11) zeros(ngen) diag(Cg13) diag(Cg14)  diag(Cg15)  diag(Cg16)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen)
              diag(Cg21) zeros(ngen) diag(Cg23) diag(Cg24)  diag(Cg25)  diag(Cg26)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen)]);
       
CG=CGmod(kk,ii);

      clear CGmod;
clear Cg11 Cg21 Cg13 Cg14 Cg15 Cg16 Cg23 Cg24  Cg25 Cg26;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Dg11=cos(delta0).*sin(dmth0)./xddd-sin(delta0).*cos(dmth0)./xqdd;

Dg12=cos(delta0).*cos(dmth0)./xddd+sin(delta0).*sin(dmth0)./xqdd;

Dg21=-sin(delta0).*sin(dmth0)./xddd-cos(delta0).*cos(dmth0)./xqdd;

Dg22=-sin(delta0).*cos(dmth0)./xddd+cos(delta0).*sin(dmth0)./xqdd;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dg11rect=(Dg11.*(-VDg0)+Dg12.*VQg0)./Vg0;
Dg12rect=(Dg11.*VQg0+Dg12.*VDg0)./Vg0;
Dg21rect=(Dg21.*(-VDg0)+Dg22.*VQg0)./Vg0;
Dg22rect=(Dg21.*VQg0+Dg22.*VDg0)./Vg0;

DGmod=sparse([diag(Dg11rect) diag(Dg12rect); diag(Dg21rect) diag(Dg22rect)]);
 
DG=DGmod(kk,kk);
YG=-DG;
clear DGmod DG;
clear Dg11rect Dg22rect Dg21rect Dg12rect;
clear Dg11 Dg22 Dg21 Dg12;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pack;

% Now back to the main program
```


```c title="genmat1.c"

% Formation of the generator matrices


%****************************************************************************************************
%****************************************************************************************************
%             FORMATION OF THE `A' MATRIX BEGINS
%****************************************************************************************************
%****************************************************************************************************

dmth0=delta0-thetag0;
C1=(xdd-xddd)./xdd;
C2=(xd-xdd).*xddd./xd./xdd;
C3=(xqd-xqdd)./xqd;
C4=(xq-xqd).*xqdd./xq./xqd;
C5=(xddd-xqdd)./xddd./xqdd;

Ag12=wB*ones(1,ngen);

Ag21=(C2.*Vg0.*sif0.*cos(dmth0)./xddd + C1.*Vg0.*sih0.*cos(dmth0)./xddd +C4.*Vg0.*sig0.*sin(dmth0)./xqdd);
Ag21=(Ag21+C3.*Vg0.*sik0.*sin(dmth0)./xqdd + C5.*Vg0.*Vg0.*cos(2.0*(dmth0)))*(-0.5)./H;
Ag22=-(Dm./H)*(0.5);

Ag23=-C2.*Vg0.*sin(dmth0)*(0.5)./H./xddd;

Ag24=-C1.*Vg0.*sin(dmth0)*(0.5)./H./xddd;

Ag25=C4.*Vg0.*cos(dmth0)*(0.5)./H./xqdd;

Ag26=C3.*Vg0.*cos(dmth0)*(0.5)./H./xqdd;

Ag31=-Vg0.*sin(dmth0)./Tdd;

Ag33=-ones(1,ngen)./Tdd;

Ag37=(xdd./(xd-xdd))./Tdd;

Ag41=-Vg0.*sin(dmth0)./Tddd;

Ag44=-ones(1,ngen)./Tddd;

Ag51=Vg0.*cos(dmth0)./Tqd;

Ag55=-ones(1,ngen)./Tqd;

Ag61=Vg0.*cos(dmth0)./Tqdd;

Ag66=-ones(1,ngen)./Tqdd;

Ag77=static_dEfd+DC1A_dAg77+AC4A_dAg77;

Ag78=DC1A_dAg78;

Ag79=AC4A_dAg79;

Ag87=DC1A_dAg87;

Ag88=DC1A_dAg88;

Ag89=DC1A_dAg89;

Ag810=DC1A_dAg810;

Ag97=DC1A_dAg97;

Ag99=DC1A_dAg99+AC4A_dAg99;

Ag910=DC1A_dAg910;

Ag107=DC1A_dAg107;

Ag1010=DC1A_dAg1010;

Ag211=RHST_dAg211;

Ag212=RHST_dAg212;

Ag213=RHST_dAg213;

Ag215=HYDRO_dAg215;

Ag216=HYDRO_dAg216;

Ag1111=RHST_dAg1111;

Ag1115=RHST_dAg1115;

Ag1211=RHST_dAg1211;

Ag1212=RHST_dAg1212;

Ag1312=RHST_dAg1312;

Ag1313=RHST_dAg1313;

Ag142=RHST_dAg142+HYDRO_dAg142;

Ag1414=RHST_dAg1414+HYDRO_dAg1414;

Ag152=RHST_dAg152+HYDRO_dAg152;

Ag1514=RHST_dAg1514+HYDRO_dAg1514;

Ag1515=RHST_dAg1515+HYDRO_dAg1515;

Ag1615=HYDRO_dAg1615;

Ag1616=HYDRO_dAg1616;




AGmod=sparse([zeros(ngen) diag(Ag12)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag21)  diag(Ag22)  diag(Ag23)    diag(Ag24)   diag(Ag25)   diag(Ag26)    zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  diag(Ag211)  diag(Ag212)  diag(Ag213)  zeros(ngen)   diag(Ag215)  diag(Ag216)
              diag(Ag31)  zeros(ngen) diag(Ag33)    zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag37)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag41)  zeros(ngen) zeros(ngen)   diag(Ag44)   zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag51)  zeros(ngen) zeros(ngen)   zeros(ngen)  diag(Ag55)   zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag61)  zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  diag(Ag66)    zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag77)  diag(Ag78)    diag(Ag79)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag87)  diag(Ag88)    diag(Ag89)   diag(Ag810)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag97)  zeros(ngen)   diag(Ag99)   diag(Ag910)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag107) zeros(ngen)   zeros(ngen)  diag(Ag1010) zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  diag(Ag1111) zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag1115) zeros(ngen)                                        
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  diag(Ag1211) diag(Ag1212) zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)   zeros(ngen) diag(Ag1312) diag(Ag1313) zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) diag(Ag142) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)  zeros(ngen)  diag(Ag1414)  zeros(ngen)  zeros(ngen)
              zeros(ngen) diag(Ag152) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)  zeros(ngen)  diag(Ag1514)  diag(Ag1515) zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)   zeros(ngen) zeros(ngen)  zeros(ngen)  zeros(ngen)   diag(Ag1615) diag(Ag1616)]); 
             
           
           
   i=[1:ngen
  (ngen+1):(2*ngen)
  (2*ngen+1):(3*ngen)
  (3*ngen+1):(4*ngen)       
  (4*ngen+1):(5*ngen)
  (5*ngen+1):(6*ngen)
  (6*ngen+1):(7*ngen)
  (7*ngen+1):(8*ngen)
  (8*ngen+1):(9*ngen)
  (9*ngen+1):(10*ngen)
  (10*ngen+1):(11*ngen)
  (11*ngen+1):(12*ngen)
  (12*ngen+1):(13*ngen)
  (13*ngen+1):(14*ngen)
  (14*ngen+1):(15*ngen)
  (15*ngen+1):(16*ngen)];

ii=i(:);

AG=AGmod(ii,ii);

clear AGmod;
%clear Ag21 Ag22 Ag23 Ag24 Ag25 Ag26
clear Ag12  Ag31 Ag33 Ag37 Ag41 Ag44  Ag51 Ag55 Ag61 Ag66 Ag77 Ag78 Ag87 Ag88 Ag89 Ag810 Ag97 Ag99 Ag910 Ag107 Ag1010;
clear Ag211 Ag212 Ag213 Ag214 Ag215 Ag216 Ag1111 Ag1115 Ag1211 Ag1212 Ag1312 Ag1313 Ag142 Ag1414 Ag152 Ag1514 Ag1515 Ag1615 Ag1616; 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bg21=(C2.*sif0.*cos(dmth0)./xddd + C1.*sih0.*cos(dmth0)./xddd +C4.*sig0.*sin(dmth0)./xqdd);

Bg21=(Bg21+C3.*sik0.*sin(dmth0)./xqdd + C5.*Vg0.*cos(2.0*(dmth0)))*(0.5)./H;

Bg22=(C2.*sif0.*sin(dmth0)./xddd + C1.*sih0.*sin(dmth0)./xddd -C4.*sig0.*cos(dmth0)./xqdd);

Bg22=(Bg22-C3.*sik0.*cos(dmth0)./xqdd + C5.*Vg0.*sin(2.0*(dmth0)))*(-0.5)./H;

Bg31=sin(dmth0)./Tdd;

Bg32=cos(dmth0)./Tdd;

Bg41=sin(dmth0)./Tddd;

Bg42=cos(dmth0)./Tddd;

Bg51=-cos(dmth0)./Tqd;

Bg52=sin(dmth0)./Tqd;

Bg61=-cos(dmth0)./Tqdd;

Bg62=sin(dmth0)./Tqdd;

Bg72=-static_bdEfd+AC4A_dBg72;

Bg82=DC1A_dBg82;

Bg92=DC1A_dBg92+AC4A_dBg92;

Bg21rect=(Bg21.*(-VDg0)+Bg22.*VQg0)./Vg0;
Bg22rect=(Bg21.*VQg0+Bg22.*VDg0)./Vg0;

Bg31rect=(Bg31.*(-VDg0)+Bg32.*VQg0)./Vg0;
Bg32rect=(Bg31.*VQg0+Bg32.*VDg0)./Vg0;

Bg41rect=(Bg41.*(-VDg0)+Bg42.*VQg0)./Vg0;
Bg42rect=(Bg41.*VQg0+Bg42.*VDg0)./Vg0;

Bg51rect=(Bg51.*(-VDg0)+Bg52.*VQg0)./Vg0;
Bg52rect=(Bg51.*VQg0+Bg52.*VDg0)./Vg0;

Bg61rect=(Bg61.*(-VDg0)+Bg62.*VQg0)./Vg0;
Bg62rect=(Bg61.*VQg0+Bg62.*VDg0)./Vg0;

Bg71rect=(Bg72.*VQg0)./Vg0;
Bg72rect=(Bg72.*VDg0)./Vg0;


Bg81rect=(Bg82.*VQg0)./Vg0;
Bg82rect=(Bg82.*VDg0)./Vg0;


Bg91rect=(Bg92.*VQg0)./Vg0;
Bg92rect=(Bg92.*VDg0)./Vg0;



BGmod=sparse([zeros(ngen)      zeros(ngen)
              diag(Bg21rect)   diag(Bg22rect)
              diag(Bg31rect)   diag(Bg32rect)
              diag(Bg41rect)   diag(Bg42rect)
              diag(Bg51rect)   diag(Bg52rect)
              diag(Bg61rect)   diag(Bg62rect)
              diag(Bg71rect)   diag(Bg72rect)
              diag(Bg81rect)   diag(Bg82rect)
              diag(Bg91rect)   diag(Bg92rect)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen) ]);
        
       

BG=BGmod(ii,kk);

clear BGmod;
%clear Bg21 Bg22
clear  Bg31 Bg32 Bg41 Bg42 Bg51 Bg52 Bg61 Bg62 Bg72 Bg82 Bg92;
%clear Bg21rect Bg22rect
clear  Bg31rect Bg32rect Bg41rect Bg42rect Bg51rect Bg52rect Bg61rect Bg62rect Bg72rect;
clear Bg81rect Bg82rect Bg91rect Bg92rect;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Input matrix used to interface PSS -see pssmat.m
EG=sparse(zeros(ngen*16,ngen));
EG([7:16:ngen*16],[1:ngen])=diag(static_bdEfd+AC4A_dEg7Efd);
EG([8:16:ngen*16],[1:ngen])=diag(DC1A_dEg8vr);
EG([9:16:ngen*16],[1:ngen])=diag(DC1A_dEg9xb+AC4A_dEg9xb);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%This matrix is useful to extract 'del s_m' from the 
%x=[deldelta,dels_m,delSif,delSih,delSig,delSik,delEfd,delvr,delxb,delxf,delx1,delx2,delx3,dely1,delPgv,delz1 ]'
FG=sparse(zeros(ngen,ngen*16));
FG([1:ngen],(16*[1:ngen]-14))=eye(ngen);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
Cg11=iQg0+Vg0.*(sin(delta0).*cos(dmth0)./xqdd-sin(dmth0).*cos(delta0)./xddd);    
 
Cg13=-C2.*cos(delta0)./xddd;

Cg14=-C1.*cos(delta0)./xddd;

Cg15=-C4.*sin(delta0)./xqdd;

Cg16=-C3.*sin(delta0)./xqdd;

Cg21=Vg0.*(sin(delta0).*sin(dmth0)./xddd+cos(delta0).*cos(dmth0)./xqdd)-iDg0; 
     
Cg23=C2.*sin(delta0)./xddd;

Cg24=C1.*sin(delta0)./xddd;

Cg25=-C4.*cos(delta0)./xqdd;

Cg26=-C3.*cos(delta0)./xqdd;


CGmod=sparse([diag(Cg11) zeros(ngen) diag(Cg13) diag(Cg14)  diag(Cg15)  diag(Cg16)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen)
              diag(Cg21) zeros(ngen) diag(Cg23) diag(Cg24)  diag(Cg25)  diag(Cg26)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen) zeros(ngen)]);
       
CG=CGmod(kk,ii);

      clear CGmod;
clear Cg11 Cg21 Cg13 Cg14 Cg15 Cg16 Cg23 Cg24  Cg25 Cg26;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Dg11=cos(delta0).*sin(dmth0)./xddd-sin(delta0).*cos(dmth0)./xqdd;

Dg12=cos(delta0).*cos(dmth0)./xddd+sin(delta0).*sin(dmth0)./xqdd;

Dg21=-sin(delta0).*sin(dmth0)./xddd-cos(delta0).*cos(dmth0)./xqdd;

Dg22=-sin(delta0).*cos(dmth0)./xddd+cos(delta0).*sin(dmth0)./xqdd;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dg11rect=(Dg11.*(-VDg0)+Dg12.*VQg0)./Vg0;
Dg12rect=(Dg11.*VQg0+Dg12.*VDg0)./Vg0;
Dg21rect=(Dg21.*(-VDg0)+Dg22.*VQg0)./Vg0;
Dg22rect=(Dg21.*VQg0+Dg22.*VDg0)./Vg0;

DGmod=sparse([diag(Dg11rect) diag(Dg12rect); diag(Dg21rect) diag(Dg22rect)]);
 
DG=DGmod(kk,kk);
YG=-DG;
clear DGmod DG;
clear Dg11rect Dg22rect Dg21rect Dg12rect;
clear Dg11 Dg22 Dg21 Dg12;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%pack;

% Now back to the main program
```


```c title="hydro_turbine.c"
%speed Governer system with Hydro Turbine.
load turb_hydro.dat
[mm order_hydro]=sort([setxor(1:nb,turb_hydro(:,1)) turb_hydro(:,1)]);
T_wht=turb_hydro(:,2)';
TG_ht=turb_hydro(:,3)';
% sig_ht -Regulation on system base= R(on m/c base)*100/P_base of m/c
sig_ht=turb_hydro(:,4)'; 
T2_ht=turb_hydro(:,5)';
Pgv_ht_max=turb_hydro(:,6)';
Pgv_ht_min=turb_hydro(:,7)';

TR_ht= 5*T_wht;
beta_ht=1.25*T_wht./H_a(turb_hydro(:,1));
TA_ht=(1./sig_ht).*TR_ht.*TG_ht;
TB_ht=(1./sig_ht).*((beta_ht+sig_ht).*TR_ht + TG_ht);
T1_ht=TB_ht/2+ sqrt(TB_ht.^2/4-TA_ht);
T3_ht=TB_ht/2-sqrt(TB_ht.^2/4-TA_ht);
K_ht=1./sig_ht;

Tm0_hydr0=TmA0(turb_hydro(:,1));
Pgv_ht0=Tm0_hydr0;

Tgr_ht0=(1-(-T_wht)./(0.5*T_wht)).*Pgv_ht0;
```


```c title="initcond.c"

%%%%%%%%%%%%%% This is called from the main program%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%% LOAD FLOW DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load lfl.dat;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%%%%%%%%%%%%%%%%%%INITIAL CONDITIONS FOR GENERATOR %%%%%%%%%%%%%%%%%%%                  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

xd=gen(:,2)';
xdd=gen(:,3)';
xddd=gen(:,4)';
Td0d=gen(:,5)';
Td0dd=gen(:,6)';
xq=gen(:,7)';
xqd=gen(:,8)';
xqdd=gen(:,9)';
Tq0d=gen(:,10)';
Tq0dd=gen(:,11)';
H=gen(:,12)';
Dm=gen(:,13)';
Vg0=lfl(gen(:,1),2)';
thetag0=pi/180*lfl(gen(:,1),3)';

Vg0bar=(Vg0.*cos(thetag0))+j*Vg0.*sin(thetag0);
Ig0bar=conj((lfl(gen(:,1),4)'+j*lfl(gen(:,1),5)')./Vg0bar);

Eqbar=(Vg0bar+(j*xq).*Ig0bar);
delta0=angle(Eqbar);
Eq=abs(Eqbar);

VQg0=real(Vg0bar);
VDg0=imag(Vg0bar);

Ig0bardq=Ig0bar.*(cos(delta0)-j*sin(delta0));
Vg0bardq=Vg0bar.*(cos(delta0)-j*sin(delta0));

iq0=real(Ig0bardq);
id0=imag(Ig0bardq);

iDg0=imag(Ig0bar);
iQg0=real(Ig0bar);

vqg0=real(Vg0bardq);
vdg0=imag(Vg0bardq);

siq0=-vdg0;
sid0=vqg0;
Efd0=Eq-(xd-xq).*id0;

% d-axis states
sih0=sid0;
sif0=sid0 + (xdd./(xd-xdd)).*Efd0;

% q-axis states
sig0=siq0;
sik0=siq0;

% dummy coil
Edummydd0=-(xqdd-xddd).*iq0;
Tdummy=0.01;

%~~~~~~~~~~~~~~~~~~~~~~field current~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
xfl=(xd.*xdd)./(xd-xdd);
ifd0= (sif0-sid0 )./xfl;
IFD0=ifd0.*xd;

%~~~~~~~~~~~~~~~~ To determine Tdd and Tddd from Td0d and Td0dd~~~~~~~~~%
   a=(1- xd./xdd + xd./xddd);
   b=-(Td0d + Td0dd);
   c=(xddd./xdd).*Td0d.*Td0dd;
   Tddd1= (-b + sqrt(b.*b - 4*a.*c))./(2*a);
   Tddd2= (-b - sqrt(b.*b - 4*a.*c))./(2*a);
   Tddd= min(Tddd1,Tddd2);
   Tdd = Td0d.*Td0dd.*(xddd./xd)./Tddd;

   
%~~~~~~~~~~~~~~~~~~~~~To determine Tqd and Tqdd from Tq0d and Tq0dd~~~~~~~~~%
   a=(1- xq./xqd + xq./xqdd);
   b=-(Tq0d + Tq0dd);
   c=(xqdd./xqd).*Tq0d.*Tq0dd;
   Tqdd1= (-b + sqrt(b.*b - 4*a.*c))./(2*a);
   Tqdd2= (-b - sqrt(b.*b - 4*a.*c))./(2*a);
   Tqdd= min(Tqdd1,Tqdd2);
   Tqd = Tq0d.*Tq0dd.*(xqdd./xq)./Tqdd;

%---------------------------------------------------------------------------%
Eqdd0=((xdd-xddd)./xdd).*sih0 + ((xd-xdd)./xd).*(xddd./xdd).*sif0;
Eddd0=-(((xqd-xqdd)./xqd).*sik0 + ((xq-xqd)./xq).*(xqdd./xqd).*sig0);

Tm0=Eqdd0.*iq0+Eddd0.*id0+id0.*iq0.*(xddd-xqdd);
%----------------------------------------------------------------------------%
clear iq0,iQg0;
clear id0,iDg0;
clear Ig0bardq;
clear Vg0bardq;


[mm order_gen]=sort([setxor(1:nb,gen(:,1)); gen(:,1)]);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Exciter Initial Condition calculations
EFD01=[zeros(1, nb-ngen) Efd0];
EFD0=EFD01(order_gen);
IFD01=[zeros(1, nb-ngen) IFD0];
IFD0_a=IFD01(order_gen);
Vg0bar1=[zeros(1, nb-ngen) Vg0bar];
Vg0bar_a=Vg0bar1(order_gen);
Ig0bar1=[zeros(1, nb-ngen) Ig0bar];
Ig0bar_a=Ig0bar1(order_gen);

%----------------------Single time constant Static type exciters-----~~~--

static_exciter

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%----------------------- IEEE DC1A -Type%Exciters------------------------

DC1A_exciter

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%----------------IEEE AC4A -Type exciter~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AC4A_exciter

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%-----------Main Selector for excitation systems.

%AVR = 1's (DISABLED); =0's (ENABLED)
%AVR=ones(1,nb);
AVR=zeros(1,nb);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SL_static(1:nb)=1;

SL_DC1A(1:nb)=1;

SL_AC4A(1:nb)=1;

%~~~~~~~~~~~~~Indicate the generator number on which a specfic type of ~~%
%~~~~~~~~~~~~~~~~~~~~exciter is to be enabled, otherwise null~~~~~~~~~~~~%

%load ng_static.dat
ng_static=[1,2,3,4];

%load ng_DC1A.dat
ng_DC1A=[];

%load ng_AC4A.dat
ng_AC4A=[];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SL_static(ng_static)=0;

SL_DC1A(ng_DC1A)=0;

SL_AC4A(ng_AC4A)=0;


%-------------------------Load Modelling----------------------------------

%~~~~~~~~~~~~~~Load varables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
VL0=lfl(ld(:,1),2)';
theL0=lfl(ld(:,1),3)'*pi/180;
PL0=ld(:,2)';
QL0=ld(:,3)';

%%%%%%%%%%%%%%%%% set p1, p2 and p3, r1, r2 and r3  in load_zip_model.m

load_zip_model

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%------Turbine and Speed Governers- Initial condition calculations.------
H1=[zeros(1, nb-ngen) H];
H_a=H1(order_gen);

Tm01=[zeros(1, nb-ngen) Tm0];
TmA0=Tm01(order_gen);

%speed Governer system with Hydro Turbine.

hydro_turbine

%Speed governer with reheat type turbine

Reheat_turbine

% ---------------------Main Selector for Turbine models---------------%
%TURB = 1's (DISABLED); =0's (ENABLED)
TURB=ones(1,nb);
%TURB=zeros(1,nb);
%TURB(4)=1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SL_HYDRO(1:nb)=1; 

SL_RHST(1:nb)=1;

%-----Indicate the generator number on which a specfic type of --------%
%------- speed-governor-turbine is to be enabled, otherwise null-------%

%load ng_hydro.dat
ng_hydro=[1,2,3,4];

%load ng_rht.dat
ng_rht=[];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SL_HYDRO(ng_hydro)=0;

SL_RHST(ng_rht)=0;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%---Main Selector for PSS----
%PSS = 1's (DISABLED); =0's (ENABLED)
PSS=ones(1,nb);
%PSS=zeros(1,nb);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SL_slip_pss(1:nb)=1;

SL_delPw_pss(1:nb)=1;

SL_power_pss(1:nb)=1;
%-----Indicate the generator number on which a specfic type of 
%------- PSS is to be enabled, otherwise nullmatrix---

%load ng_slip_pss.dat
ng_slip_pss=[1];

%load ng_delPw_pss.dat
ng_delPw_pss=[];

%load ng_power_pss.dat
ng_power_pss=[];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SL_slip_pss(ng_slip_pss)=0;

SL_delPw_pss(ng_delPw_pss)=0;

SL_power_pss(ng_power_pss)=0;

%---------------PSS models----------------------
%------Slip signal based PSS----m


pss_slip_signal


%----Delta Power-Omega signal based PSS----
pss_delPw_signal


%----Power signal based PSS----
pss_power_signal
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%---------------------------------------------------------------------
busang=lfl(:,3)*pi/180;
Vpre=lfl(:,2).*(cos(busang)+j*sin(busang));
%---------------------------------------------------------------------
```


```data title="Id.dat"
  34     0.4505    0.4656     
  35     0.4919    0.2753    
  51     0.5845    0.2844     
  58     0.7630   -0.1080     
  66     1.0220    0.2670     
  68     0.00     -0.0741     
  70     0.00      0.5663     
  71     0.00     -0.2120     
  74     0.8190    0.4370     
  78     0.8900    0.2680     
  79     0.0910    0.0300   
  80     0.1710    0.0500    
  81     0.8220   -0.9310     
  82     0.0210    0.0110   
  84     0.2430    0.0820     
  85     0.2740    0.0030     
  88     0.6900    0.2090     
  89     0.0060    0.0020   
  90     0.0460    0.0150    
  92     0.00      0.3102     
  93     1.0040    0.7320  
  94     0.1540    0.0760   
  95     0.0670    0.0220    
  99     0.1046    0.0523   
 101     0.1780    0.0450   
 102     0.3760    0.0920  
 104     0.3020    0.0760  
 105     0.9600    1.6740  
 106     0.6400    0.1600  
 107    -0.1750   -0.1280   
 110     1.0040    0.7320   
 111     0.6040   11.6600 
 112     0.1860    0.0460  
 115     6.8350    1.8470  
 116     7.9260    3.1550  
 117     4.8530    0.7140  
 118     6.5190    3.2840  
 119    20.9400   37.7400  
 120    -4.0800    1.7510  
 121     2.3770   -0.1730  
 122     0.2920    0.0700  
 123    -0.8400   -0.1900  
 124     0.9410    7.8030  
 125    -7.1200   -3.1900  
 126    -3.3300   -1.6000  
 127    -5.4600   -0.7200   
 128    40.7500    7.0350 
 129    -4.8200   -1.2200 
 130    43.2800    9.4430 
 131   218.4000   43.2000 
 132     4.9190    1.1020 
 133    -0.8300   -0.3630 
 134   223.0900   74.0200 
 135    42.9800   12.6400 
 136   529.5100  135.5200 
 137   129.4600   26.0800 
 138    -3.6300   -1.8800    
 139   577.1800  139.3600 
 140   247.7500   66.7600 
 141   327.9900  113.6100 
 142  177.3700    39.3400 
 143   46.7200    17.0900 
 144   96.0200    22.0300 
 145   91.7300    15.5500 
```


```data title="lfl.dat"
  1 	     1.080998 	   -5.029944	    0.000000	    0.000000	 0.000000	 0.000000 
  2 	     1.080853 	   -5.096552	    0.000000	    0.000000	 0.000000	 0.000000 
  3 	     1.101531 	   -4.732669	    0.000000	    0.000000	 0.000000	 0.000000 
  4 	     1.101531 	   -4.732669	    0.000000	    0.000000	 0.000000	 0.000000 
  5 	     1.101814 	   -4.732477	    0.000000	    0.000000	 0.000000	 0.000000 
  6 	     1.043300 	   -8.545121	    0.000000	    0.000000	 0.000000	 0.000000 
  7 	     1.076301 	    2.501586	    0.000000	    0.000000	 0.000000	 0.000000 
  8 	     1.113652 	    0.440080	    0.000000	    0.000000	 0.000000	 0.000000 
  9 	     1.039591 	   -8.760132	    0.000000	    0.000000	 0.000000	 0.000000 
 10 	     1.039568 	   -8.760931	    0.000000	    0.000000	 0.000000	 0.000000 
 11 	     1.093670 	  -11.365751	    0.000000	    0.000000	 0.000000	 0.000000 
 12 	     1.038865 	   -9.477825	    0.000000	    0.000000	 0.000000	 0.000000 
 13 	     1.098189 	  -12.138574	    0.000000	    0.000000	 0.000000	 0.000000 
 14 	     1.038524 	   -9.885644	    0.000000	    0.000000	 0.000000	 0.000000 
 15 	     1.068304 	  -10.519430	    0.000000	    0.000000	 0.000000	 0.000000 
 16 	     1.068638 	  -10.572840	    0.000000	    0.000000	 0.000000	 0.000000 
 17 	     1.001236 	  -10.148720	    0.000000	    0.000000	 0.000000	 0.000000 
 18 	     1.074645 	  -11.587839	    0.000000	    0.000000	 0.000000	 0.000000 
 19 	     1.070821 	  -11.667798	    0.000000	    0.000000	 0.000000	 0.000000 
 20 	     1.113056 	  -11.666558	    0.000000	    0.000000	 0.000000	 0.000000 
 21 	     1.108556 	  -11.949572	    0.000000	    0.000000	 0.000000	 0.000000 
 22 	     1.031053 	   -4.591935	    0.000000	    0.000000	 0.000000	 0.000000 
 23 	     1.097860 	   -6.217242	    0.000000	    0.000000	 0.000000	 0.000000 
 24 	     1.027216 	    1.596959	    0.000000	    0.000000	 0.000000	 0.000000 
 25 	     1.037961 	  -10.579214	    0.000000	    0.000000	 0.000000	 0.000000 
 26 	     1.089437 	  -12.082964	    0.000000	    0.000000	 0.000000	 0.000000 
 27 	     1.038853 	  -13.782173	    0.000000	    0.000000	 0.000000	 0.000000 
 28 	     1.076213 	  -15.987627	    0.000000	    0.000000	 0.000000	 0.000000 
 29 	     1.074609 	  -16.150450	    0.000000	    0.000000	 0.000000	 0.000000 
 30 	     1.073074 	   -6.057738	    0.000000	    0.000000	 0.000000	 0.000000 
 31 	     1.090534 	  -12.522370	    0.000000	    0.000000	 0.000000	 0.000000 
 32 	     1.093709 	  -11.373291	    0.000000	    0.000000	 0.000000	 0.000000 
 33 	     1.139245 	   -4.770083	    0.000000	    0.000000	 0.000000	 0.000000 
 34 	     1.138711 	   -4.710749	    0.000000	    0.000000	 0.450500	 0.465600 
 35 	     1.139002 	   -4.788903	    0.000000	    0.000000	 0.491900	 0.275300 
 36 	     1.138542 	   -4.523437	    0.000000	    0.000000	 0.000000	 0.000000 
 37 	     1.123516 	   -6.940806	    0.000000	    0.000000	 0.000000	 0.000000 
 38 	     1.130577 	   -6.001223	    0.000000	    0.000000	 0.000000	 0.000000 
 39 	     1.126959 	   -8.624708	    0.000000	    0.000000	 0.000000	 0.000000 
 40 	     1.126939 	   -8.628814	    0.000000	    0.000000	 0.000000	 0.000000 
 41 	     1.118840 	  -11.140431	    0.000000	    0.000000	 0.000000	 0.000000 
 42 	     1.118813 	  -11.154647	    0.000000	    0.000000	 0.000000	 0.000000 
 43 	     1.118947 	  -11.112560	    0.000000	    0.000000	 0.000000	 0.000000 
 44 	     1.118920 	  -11.126533	    0.000000	    0.000000	 0.000000	 0.000000 
 45 	     1.117292 	  -12.123368	    0.000000	    0.000000	 0.000000	 0.000000 
 46 	     1.117299 	  -12.118313	    0.000000	    0.000000	 0.000000	 0.000000 
 47 	     1.127532 	   -7.433984	    0.000000	    0.000000	 0.000000	 0.000000 
 48 	     1.127833 	   -7.414530	    0.000000	    0.000000	 0.000000	 0.000000 
 49 	     1.127903 	   -7.405086	    0.000000	    0.000000	 0.000000	 0.000000 
 50 	     1.127608 	   -7.423941	    0.000000	    0.000000	 0.000000	 0.000000 
 51 	     1.112399 	  -10.867749	    0.000000	    0.000000	 0.584500	 0.284400 
 52 	     1.111781 	  -11.840752	    0.000000	    0.000000	 0.000000	 0.000000 
 53 	     1.111781 	  -11.842287	    0.000000	    0.000000	 0.000000	 0.000000 
 54 	     1.113112 	  -12.494277	    0.000000	    0.000000	 0.000000	 0.000000 
 55 	     1.113120 	  -12.494304	    0.000000	    0.000000	 0.000000	 0.000000 
 56 	     1.107153 	  -10.650914	    0.000000	    0.000000	 0.000000	 0.000000 
 57 	     1.107201 	  -10.652361	    0.000000	    0.000000	 0.000000	 0.000000 
 58 	     1.106654 	  -10.471996	    0.000000	    0.000000	 0.763000	-0.108000 
 59 	     1.116455 	  -11.550636	    0.000000	    0.000000	 0.000000	 0.000000 
 60 	     1.137000 	   -7.075880	    0.510000	    0.329164	 0.000000	 0.000000 
 61 	     1.114429 	  -12.599427	    0.000000	    0.000000	 0.000000	 0.000000 
 62 	     1.056600 	  -15.177319	    0.000000	    0.000000	 0.000000	 0.000000 
 63 	     1.110920 	  -14.688829	    0.000000	    0.000000	 0.000000	 0.000000 
 64 	     1.098000 	   -9.996327	    0.000000	    0.000000	 0.000000	 0.000000 
 65 	     1.097975 	   -9.998114	    0.000000	    0.000000	 0.000000	 0.000000 
 66 	     1.112887 	    0.609349	    0.000000	    0.000000	 1.022000	 0.267000 
 67 	     1.090000 	   -6.366411	   14.860000	    2.852061	 0.000000	 0.000000 
 68 	     1.208597 	  -31.694046	    0.000000	    0.000000	 0.000000	-0.074100 
 69 	     1.096753 	  -11.125235	    0.000000	    0.000000	 0.000000	 0.000000 
 70 	     0.999829 	  -14.875302	    0.000000	    0.000000	 0.000000	 0.566300 
 71 	     1.027483 	  -14.968905	    0.000000	    0.000000	 0.000000	-0.212000 
 72 	     1.100723 	  -11.902767	    0.000000	    0.000000	 0.000000	 0.000000 
 73 	     1.097523 	  -11.768184	    0.000000	    0.000000	 0.000000	 0.000000 
 74 	     1.097266 	  -12.169410	    0.000000	    0.000000	 0.819000	 0.437000 
 75 	     1.117872 	  -15.896103	    0.000000	    0.000000	 0.000000	 0.000000 
 76 	     1.020878 	    4.828444	    0.000000	    0.000000	 0.000000	 0.000000 
 77 	     0.987963 	    6.013543	    0.000000	    0.000000	 0.000000	 0.000000 
 78 	     1.073979 	   -5.896450	    0.000000	    0.000000	 0.890000	 0.268000 
 79 	     1.052000 	  -10.218608	    2.502000	   -0.159510	 0.091000	 0.030000 
 80 	     1.069000 	   -8.918395	    0.470000	   -0.150649	 0.171000	 0.050000 
 81 	     1.130389 	  -26.573264	    0.000000	    0.000000	 0.822000	-0.931000 
 82 	     0.975000 	  -19.372126	    0.700000	    0.171541	 0.021000	 0.011000 
 83 	     1.098486 	   -6.088671	    0.000000	    0.000000	 0.000000	 0.000000 
 84 	     1.115578 	  -10.147721	    0.000000	    0.000000	 0.243000	 0.082000 
 85 	     1.116493 	  -13.754570	    0.000000	    0.000000	 0.274000	 0.003000 
 86 	     1.056689 	  -14.718246	    0.000000	    0.000000	 0.000000	 0.000000 
 87 	     1.065151 	   -7.879955	    0.000000	    0.000000	 0.000000	 0.000000 
 88 	     1.109420 	   -9.053508	    0.000000	    0.000000	 0.690000	 0.209000 
 89 	     1.066000 	    2.975973	    6.730000	    1.363931	 0.006000	 0.002000 
 90 	     0.950000 	   -8.062662	    0.220000	   -0.038652	 0.046000	 0.015000 
 91 	     1.000000 	   -9.985229	    0.640000	   -0.015385	 0.000000	 0.000000 
 92 	     0.956120 	  -13.459889	    0.000000	    0.000000	 0.000000	 0.310200 
 93 	     1.000000 	   -2.626796	    7.000000	    3.738139	 1.004000	 0.732000 
 94 	     1.020000 	   -1.450779	    3.000000	    0.190526	 0.154000	 0.076000 
 95 	     0.920000 	   18.172369	    1.310000	    0.101211	 0.067000	 0.022000 
 96 	     1.000000 	   -9.685438	    0.600000	    0.211089	 0.000000	 0.000000 
 97 	     0.967000 	   -5.053482	    1.400000	    0.456264	 0.000000	 0.000000 
 98 	     0.970000 	    4.479341	    4.260000	   -0.327207	 0.000000	 0.000000 
 99 	     1.000000 	    0.388998	    2.000000	   -0.083408	 0.104600	 0.052300 
100 	     1.014000 	    0.000000	    1.701011	    0.587239	 0.000000	 0.000000 
101 	     1.039000 	   -6.798184	    3.109000	    1.486615	 0.178000	 0.045000 
102 	     1.019000 	   -5.473517	   20.400000	    4.889027	 0.376000	 0.092000 
103 	     1.000000 	    0.807694	    1.350000	    0.049585	 0.000000	 0.000000 
104 	     1.005900 	   12.967922	   20.000000	    4.999287	 0.302000	 0.076000 
105 	     1.007000 	   -3.503571	   16.200000	    3.883489	 0.960000	 1.674000 
106 	     1.005000 	   -3.458237	   10.800000	    2.093641	 0.640000	 0.160000 
107 	     1.021077 	  -14.282204	    0.000000	    0.000000	-0.175000	-0.128000 
108 	     1.014000 	  -14.739737	    8.000000	    0.772871	 0.000000	 0.000000 
109 	     0.915000 	  -19.163968	    0.520000	   -0.155512	 0.000000	 0.000000 
110 	     1.000000 	   -2.015779	    7.000000	    5.198349	 1.004000	 0.732000 
111 	     1.000000 	    7.262588	   20.000000	    5.637547	 0.604000	11.660000 
112 	     1.037000 	   -6.972694	    3.000000	    1.401114	 0.186000	 0.046000 
113 	     0.977971 	   -5.096552	    0.000000	    0.000000	 0.000000	 0.000000 
114 	     0.977971 	   -5.096552	    0.000000	    0.000000	 0.000000	 0.000000 
115 	     1.049000 	  -16.320584	   24.930000	    1.427262	 6.835000	 1.847000 
116 	     1.043000 	  -17.571642	   27.130000	    6.318378	 7.926000	 3.155000 
117 	     1.030000 	  -16.033238	   26.270000	    2.585555	 4.853000	 0.714000 
118 	     1.010000 	  -18.502919	   42.200000	    6.603727	 6.519000	 3.284000 
119 	     1.013000 	  -60.121676	   89.540000	   47.485070	20.940000	37.740000 
120 	     1.033103 	  -52.314276	    0.000000	    0.000000	-4.080000	 1.751000 
121 	     1.046000 	  -20.905301	   29.970000	   -1.601911	 2.377000	-0.173000 
122 	     1.000000 	   -3.498467	   10.090000	    1.740455	 0.292000	 0.070000 
123 	     1.017115 	  -33.831875	    0.000000	    0.000000	-0.840000	-0.190000 
124 	     1.000000 	   -2.594616	   30.050000	    5.692040	 0.941000	 7.803000 
125 	     1.008381 	  -33.303513	    0.000000	    0.000000	-7.120000	-3.190000 
126 	     1.052383 	  -74.609654	    0.000000	    0.000000	-3.330000	-1.600000 
127 	     1.006960 	  -37.108854	    0.000000	    0.000000	-5.460000	-0.720000 
128 	     1.025000 	  -40.415264	  129.630000	   26.108349	40.750000	 7.035000 
129 	     0.980186 	  -73.782875	    0.000000	    0.000000	-4.820000	-1.220000 
130 	     1.057000 	  -52.574928	   59.370000	   18.349495	43.280000	 9.443000 
131 	     1.042000 	  -25.026583	  283.000000	   74.730319	218.400000	43.200000 
132 	     1.042000 	   -7.951271	   30.950000	    6.334228	 4.919000	 1.102000 
133 	     1.092218 	  -12.308207	    0.000000	    0.000000	-0.830000	-0.363000 
134 	     1.044000 	  -11.531059	  206.260000	   74.021400	223.090000	74.020000 
135 	     1.107000 	   28.334386	   59.820000	   15.648412	42.980000	12.640000 
136 	     1.083000 	    3.677557	  519.500000	  144.535002	529.510000	135.520000 
137 	     1.064000 	  -73.438423	  120.680000	   34.506890	129.460000	26.080000 
138 	     1.113793 	   11.301156	    0.000000	    0.000000	-3.630000	-1.880000 
139 	     1.040000 	  -11.267381	  568.340000	  158.495909	577.180000	139.360000 
140 	     1.050000 	  -26.874223	  231.230000	   67.104727	247.750000	66.760000 
141 	     1.053000 	   -9.830327	  379.110000	  116.694923	327.990000	113.610000 
142 	     1.155000 	  -11.441815	  244.490000	   54.961354	177.370000	39.340000 
143 	     1.031000 	  -14.373714	   52.540000	   21.586259	46.720000	17.090000 
144 	     0.997000 	   -9.287162	  113.970000	   26.868461	96.020000	22.030000 
145 	     1.052000 	    4.309572	  141.186200	   29.871238	91.730000	15.550000 
```


```data title="lfl_result.c"


Bno=[1:nb]';
PG=sparse(zeros(nb,1));
QG=sparse(zeros(nb,1));
PL=sparse(zeros(nb,1));
QL=sparse(zeros(nb,1));
QG_pv = Qc(pv_data(:,1))+ Qload(pv_data(:,1));

Pc=real(S);
Qc=imag(S);
PG([sl;pv_data(:,1)])=[Pc(sl); pv_data(:,3)];
QG([sl;pv_data(:,1)])=[Qc(sl); QG_pv];
PL([sl;pq_data(:,1)])=[sl_load(1,2);pq_data(:,2)];
QL([sl;pq_data(:,1)])=[sl_load(1,3);pq_data(:,3)];

fid=fopen('report.dat', 'a');
fprintf(fid,' Loadflow results:\n'); 
fprintf(fid,'---------------------------------------------------------------------------------------\n');
fprintf(fid,'\nBus No     VbO          thetaO           PGO             QGO         PLO         QLO \n');
fprintf(fid,'---------------------------------------------------------------------------------------\n');
fprintf(fid,'%3i \t %12.6f \t%12.6f\t%12.6f\t%12.6f\t%9.6f\t%9.6f \n',[Bno  full(Vmag) full(Vang*180/pi) full(PG)  full(QG)  full(PL) full(QL)]');
fprintf(fid,'---------------------------------------------------------------------------------------\n');
fprintf(fid,'Line flows:\n');
fprintf(fid,'----------------------------------------------------------------------------------\n');
fprintf(fid,'                       Line flows                                Line flows \n');
fprintf(fid,'                   _____________________                   _______________________\n');
fprintf(fid,' From   To         P-flow         Q-flow    From   To      P-flow           Q-flow\n');
fprintf(fid,'----------------------------------------------------------------------------------\n');
fprintf(fid,'%3i \t %3i\t %12.4f\t %12.4f\t %3i\t %4i\t %12.4f\t %12.4f\n',[nt(:,1) nt(:,2) full(real(Spq)) full(imag(Spq)) nt(:,2) nt(:,1) full(real(Sqp)) full(imag(Sqp))]');
fprintf(fid,'----------------------------------------------------------------------------------\n');
fprintf(fid,'Total real power losses in the system = %9.6f\n',full(real(S_line_loss)));
fprintf(fid,'Total reactive power losses in the system = %9.6f\n',full(imag(S_line_loss)));

if (convergence_bit==1)
 fid1=fopen('lfl.dat','w');
 fprintf(fid1,'%3i \t %12.6f \t%12.6f\t%12.6f\t%12.6f\t%9.6f\t%9.6f \n',[Bno  full(Vmag) full(Vang*180/pi) full(PG)  full(QG)  full(PL) full(QL)]');
 fclose(fid1);
else
 fprintf(fid,'Convergence is not reached !!!' );
end
fclose(fid);
```


```c title="load_zip_model.c"
%-------------------Load Modelling------------------

%Real power component: (p1+p2+p3=1)
p1=0;
p2=0;
p3=1;
%Reactive power component: (r1+r2+r3 =1)
r1=0;
r2=0;
r3=1;

% constant Power type
kcpr=r1*QL0;
kcp=p1*PL0;
% constant current type
kccr=r2*QL0./VL0;
kcc= p2*PL0./VL0;
% constant impedance type
kcir=r3*QL0./(VL0.^2);
kci=p3*PL0./(VL0.^2);

TL=0.01;
% Selector for  loads%

[mm order_load]=sort([setxor(1:nb,ld(:,1));ld(:,1)]);
```


```data title="nt.dat"
   1    2    0.000030  0.000800   0.06320         
   1    2    0.000030  0.000800   0.06320
   1    6    0.001940  0.020900   2.37920   
   2    6    0.001940  0.020900   2.37920         
   3   33    0.000200  0.022100   0.00000         
   4   33    0.000200  0.022100   0.00000         
   5   33    0.000200  0.021900   0.00000         
   6    7    0.001290  0.013900   1.46520         
   6    9    0.000160  0.001700   0.17520         
   6   10    0.000160  0.001700   0.17520         
   6   12    0.000200  0.002100   0.87760         
   6   12    0.000200  0.002100   0.87760         
   8   66    0.000200  0.029900   0.00000         
   8   66    0.000200  0.022100   0.00000         
  11   69    0.000200  0.026200   0.00000        
  12   14    0.000960  0.009100   0.85560        
  12   14    0.000960  0.009100   0.85560        
  12   25    0.000510  0.005500   0.62500     
  12   25    0.000510  0.005500   0.62500     
  13   72    0.000200  0.026000   0.00000        
  13   72    0.000300  0.026200   0.00000     
  13   72    0.000200  0.026000   0.00000        
  14   17    0.003390  0.036700   3.45820     
  14   17    0.003520  0.036700   3.45160     
  15   58    0.000200  0.025500   0.00000       
  16   58    0.000200  0.022000   0.00000       
  17   22    0.002280  0.027600   2.62040     
  18   59    0.000200  0.029800   0.00000       
  19   59    0.000000  0.062900   0.00000      
  20   59    0.000000  0.063800   0.00000    
  21   59    0.000200  0.032900   0.00000    
  22   24    0.001730  0.020800   1.96480 
  23   83    0.000400  0.059500   0.00000    
  23   83    0.000300  0.059700   0.00000   
  25   27    0.002300  0.026600   3.05080 
  25   27    0.002300  0.026600   3.05080 
  26   73    0.000300  0.026700   0.00000   
  28   75    0.000200  0.029000   0.00000   
  29   75    0.000200  0.026900   0.00000   
  30   78    0.000000  0.033500   0.00000    
  31   74    0.000300  0.027900   0.00000   
  32   69    0.000200  0.026500   0.00000   
  33   34    0.000060  0.000900   0.00060  
  33   35    0.000060  0.000900   0.00060  
  33   37    0.009960  0.070700   0.11160   
  33   38    0.009950  0.069300   0.11100 
  33   39    0.008500  0.069900   0.10060 
  33   40    0.008490  0.069800   0.10040 
  33   49    0.005600  0.049300   0.07780 
  33   50    0.005600  0.049300   0.07780 
  34   36    0.000250  0.002200   0.00060  
  37   88    0.003100  0.165100   0.00000   
  38   88    0.003100  0.163800   0.00000   
  39   43    0.006020  0.049500   0.07120 
  39   84    0.007220  0.278600   0.00000    
  40   44    0.006030  0.049600   0.07140 
  40   84    0.007290  0.275600   0.00000   
  41   42    0.000500  0.151400   0.00000
  41   43    0.000010  0.000900   0.00060    
  42   44    0.000010  0.000900   0.00060    
  43   46    0.006180  0.050800   0.07320 
  44   45    0.006180  0.050800   0.07320 
  45   61    0.004450  0.036600   0.05260 
  45   85    0.000000  0.260000   0.00000    
  46   61    0.004450  0.036600   0.05260 
  46   85    0.000000  0.259200   0.00000    
  47   48   -0.010000  0.230600   0.00000    
  47   50    0.000010  0.000900   0.00060  
  47   87    0.083100  0.401000   0.00000   
  48   49    0.000010  0.000900   0.00060    
  48   87    0.099800  0.436000   0.00000    
  49   51    0.008980  0.079000   0.12480 
  50   51    0.008980  0.079000   0.12480 
  51   52    0.002900  0.027900   0.04660  
  51   53    0.002900  0.027900   0.04660 
  51   56    0.007590  0.048300   0.07120 
  51   57    0.007590  0.048300   0.07120 
  52   53   -0.006700  0.391100   0.00000   
  52   54    0.004700  0.029300   0.04620    
  53   55    0.004700  0.029300   0.04620    
  54   55   -0.055300  0.928900   0.00000    
  54   61    0.001410  0.008700   0.01380  
  55   61    0.001410  0.008700   0.01380  
  56   57   -0.009000  0.389500   0.00000   
  56   58    0.001900  0.012000   0.01780   
  57   58    0.001900  0.012000   0.01780    
  58   59    0.667400  2.217500   0.00000   
  58   72    0.030200  0.236400   0.00000   
  58   87    0.086300  0.390600   0.00000   
  58   98    0.013100  0.176500   0.00000   
  58  100    0.119300  1.269000   0.00000    
  58  103    0.841600  5.538300   0.00000  
  59   60   -0.180300  5.965900   0.00000   
  59   72    0.861300  3.048500   0.00000   
  59   79    0.009900  0.264400   0.00000    
  59   80    0.287600  2.389800   0.00000   
  59   89    0.342100  9.057100   0.00000   
  59   92   -0.007000  0.567800   0.00000     
  59   94    0.704100  5.988500   0.00000   
  59   98    0.106000  0.584500   0.00000    
  59  100    0.018300  0.201600   0.00000   
  59  103    0.036800  0.334100   0.00000    
  59  107    0.037200  0.883400   0.00000    
  60  135   -1.831000  9.796400   0.00000   
  60   79   -0.037500  1.106800   0.00000   
  60   80    0.065500  2.644100   0.00000   
  60   90   -0.020100  1.513500   0.00000  
  60   92   -0.264000  3.713900   0.00000    
  60   94    0.001200  0.077500   0.00000  
  60   95   -0.085500  0.992600   0.00000    
  60  138   -0.363900  1.793600   0.00000   
  61   63    0.008120  0.078200   0.13180 
  61   63    0.008120  0.078200   0.13180 
  61   64    0.002420  0.031800   0.05680 
  61   65    0.002420  0.031800   0.05680 
  62   86    0.003600  0.050100   0.00000   
  62   86    0.001300  0.083800   0.00000   
  63   64    0.014700  0.282500   0.00000    
  63   65    0.014700  0.281300   0.00000    
  63   66    0.005600  0.090000   0.00000   
  63   67    0.032100  0.278500   0.00000   
  63   69    0.010700  0.157100   0.00000   
  63  102    0.010600  0.158300   0.00000  
  63  102    0.010600  0.157600   0.00000   
  63  102    0.010700  0.160400   0.00000    
  63  102    0.010400  0.154200   0.00000  
  63  116   -0.389700  6.858800   0.00000     
  63  117    0.003000  0.056000   0.00000   
  63  118   -0.012500  0.242500   0.00000   
  63  124   -0.126500  2.022000   0.00000   
  64   65    0.001300  0.167400   0.00000    
  64   66    0.003900  0.068400   0.00000       
  64   67    0.023300  0.212000   0.00000   
  64   69    0.007500  0.119600   0.00000    
  64   97   -0.433600  8.292300   0.00000     
  64  124   -0.104100  1.537500   0.00000    
  65   66    0.003900  0.068200   0.00000   
  65   67    0.023300  0.211100   0.00000   
  65   69    0.007500  0.119100   0.00000   
  65   97   -0.429200  8.258200   0.00000   
  65  124   -0.103200  1.531200   0.00000     
  66   67    0.008100  0.067500   0.00000    
  66   68   -2.473000  2.472000   0.00000   
  66   69    0.002800  0.038100   0.00000   
  66   97   -0.111900  2.643200   0.00000     
  66  111    0.000000  0.026400   0.00000   
  66  111    0.000570  0.026600   0.00000   
  66  111    0.000000  0.027300   0.00000  
  66  111    0.000570  0.026400   0.00000    
  66  124   -0.028300  0.490200   0.00000    
  67   68   -3.443000  3.717200   0.00000    
  67   69    0.006100  0.055000   0.00000   
  67   97    0.006300  0.116600   0.00000    
  67  119   -0.221300  9.391800   0.00000   
  67  120   -0.003400  1.784700   0.00000    
  67  121    0.008200  1.170000   0.00000    
  67  122   -0.004700  0.447300   0.00000   
  67  124    0.000300  0.006500   0.00000   
  67  125    0.006200  0.251900   0.00000     
  67  132   -0.319400  4.356600   0.00000     
  68   69   -0.692000  0.698400   0.00000       
  69   70    0.008500  0.333300   0.00000   
  69   71    0.007500  0.312000   0.00000       
  69   72    0.001300  0.010000   0.00000   
  69   73    0.009800  0.074700   0.00000    
  69   74    0.013500  0.074100   0.00000     
  69   97   -0.067400  1.584900   0.00000    
  69  101    0.017400  0.218800   0.00000        
  69  112    0.017500  0.220100   0.00000     
  69  124   -0.026700  0.398600   0.00000     
  70   71   -0.489100  2.661300   0.00000      
  70   72   -0.006200  0.121600   0.00000    
  70   73   -0.042400  0.912500   0.00000       
  70   74    0.003200  0.913800   0.00000   
  70  101   -0.124800  1.040900   0.00000        
  70  112   -0.125700  1.047100   0.00000   
  71   72   -0.006000  0.113800   0.00000   
  71   73   -0.040900  0.854100   0.00000   
  71   74    0.001800  0.855300   0.00000   
  71  101   -0.159200  1.230300   0.00000   
  71  112   -0.160300  1.237700   0.00000    
  72   73    0.001500  0.027500   0.00000   
  72   74    0.002800  0.027400   0.00000    
  72   98    0.013800  0.241700   0.00000   
  72  100    0.133700  1.738400   0.00000     
  72  101    0.000200  0.080200   0.00000   
  72  103    1.022400  7.594500   0.00000     
  72  112    0.000200  0.080600   0.00000   
  73   74   -0.000700  0.039300   0.00000     
  73   75    0.014700  0.258100   0.00000      
  73   81   -0.012200  0.306800   0.00000     
  73   82    0.003600  2.016900   0.00000     
  73   91    0.027100  0.573200   0.00000   
  73   96    0.024500  0.480500   0.00000   
  73  101    0.004400  0.601400   0.00000  
  73  105    0.000700  0.032500   0.00000   
  73  105    0.000700  0.032500   0.00000   
  73  105    0.000600  0.029500   0.00000     
  73  108   -0.018200  0.583200   0.00000   
  73  109    0.052400  3.005900   0.00000   
  73  112    0.004300  0.605000   0.00000    
  73  121   -0.026800  1.765300   0.00000    
  74   75    0.021500  0.327700   0.00000    
  74   81   -0.033300  0.463100   0.00000   
  74   82   -0.009800  1.985900   0.00000   
  74   91    0.041300  0.751100   0.00000    
  74   96    0.435000  7.690100   0.00000    
  74  101    0.034400  0.600500   0.00000   
  74  106    0.003000  0.033500   0.00000    
  74  106    0.000500  0.032800   0.00000   
  74  108   -0.018700  0.454400   0.00000   
  74  109    0.100400  3.469700   0.00000   
  74  112    0.034500  0.604200   0.00000    
  74  121   -0.034800  1.375700   0.00000     
  75   82    0.077700  1.125000   0.00000    
  75   91   -0.225500  3.144200   0.00000    
  75   96   -0.451600  4.631000   0.00000   
  75  108    0.004200  0.104900   0.00000    
  75  109    0.104600  1.446500   0.00000   
  75  121    0.017800  0.317200   0.00000    
  76   77    0.000200  0.016000   0.00000   
  76   89    0.001100  0.022100   0.00000     
  79   80    0.044000  0.099100   0.00000   
  79   90    0.050600  2.471000   0.00000   
  79   92    0.001700  0.303200   0.00000     
  79   94    0.127500  1.119500   0.00000   
  79   95    0.305000  6.415400   0.00000    
  79  107    0.078600  1.414000   0.00000   
  80   90    0.465800  5.875600   0.00000  
  80   92    0.119200  1.505300   0.00000    
  80   94    0.460000  2.647500   0.00000    
  82   91   -0.234900  2.418800   0.00000    
  82  108   -0.074200  0.727800   0.00000   
  82  109   -0.007100  0.263400   0.00000   
  82  121   -0.189200  2.205400   0.00000    
  83   89    0.058200  0.385500   0.00000    
  89  103   -1.073000  4.143300   0.00000    
  90   92   -0.138000  8.295900   0.00000   
  90   94    0.068900  1.071700   0.00000   
  91   96   -0.122400  4.246300   0.00000  
  91  108   -0.107800  0.699400   0.00000    
  91  109   -0.269900  4.263400   0.00000   
  91  121   -0.292400  2.121000   0.00000     
  92   94    0.288300  3.771700   0.00000   
  92  107    0.017600  3.022700   0.00000    
  94   95    0.053400  0.996000   0.00000    
  94  138   -0.112500  1.838500   0.00000   
  95  138   -0.073200  0.638900   0.00000    
  96  108   -0.821500  6.114300   0.00000      
  97  124   -0.379300  1.955700   0.00000      
  98  100   -0.006300  0.326900   0.00000     
  98  103    0.054400  1.435800   0.00000    
 100  103   -0.024900  0.489100   0.00000     
 101  112   -0.013800  0.361000   0.00000   
 102  117   -0.000300  0.019000   0.00000      
 102  118   -0.026700  0.322200   0.00000      
 108  109   -0.082500  1.271300   0.00000    
 108  121   -0.000900  0.043100   0.00000     
 109  121   -0.188100  3.849900   0.00000    
 115  116    0.000800  0.029100   0.00000    
 115  117   -0.009200  0.222200   0.00000     
 115  118   -0.004400  0.067700   0.00000     
 115  143   -0.101700  0.492400   0.00000     
 116  117    0.001910  0.028800   0.00000   
 116  118   -0.001000  0.044000   0.00000     
 116  143   -0.218700  1.289600   0.00000     
 117  118    0.000800  0.008100   0.00000     
 117  143   -0.083400  0.685400   0.00000   
 118  131   -0.892500  6.238500   0.00000     
 118  132   -0.696700  8.143000   0.00000   
 118  143   -0.001100  0.023100   0.00000     
 119  120    0.001000  0.023600   0.00000        
 119  121   -0.011000  0.290100   0.00000   
 119  122   -0.601300  5.894100   0.00000   
 119  124   -0.261800  3.394000   0.00000    
 119  125   -0.008200  0.259500   0.00000   
 119  126    0.001530  0.017900   0.00000    
 119  127   -0.117200  1.393200   0.00000    
 119  128   -0.005400  0.051600   0.00000     
 119  129    0.003400  0.064200   0.00000     
 119  130   -0.002200  0.016300   0.00000     
 119  131   -0.004400  0.024200   0.00000      
 119  132   -0.413700  2.402700   0.00000       
 119  144   -0.851100  3.835800   0.00000     
 120  121    0.000900  0.077900   0.00000       
 120  122   -0.061000  0.930500   0.00000     
 120  123   -0.046600  0.501100   0.00000      
 120  124   -0.025900  0.472200   0.00000     
 120  125   -0.000200  0.055500   0.00000    
 120  127    0.002000  0.181800   0.00000    
 120  128   -0.002900  0.074300   0.00000       
 120  129   -0.022900  0.491100   0.00000     
 120  130   -0.167400  1.067500   0.00000     
 120  131   -0.068700  0.451600   0.00000      
 120  132   -0.025500  0.456600   0.00000      
 121  122   -0.010800  0.483000   0.00000       
 121  123   -0.171200  1.948200   0.00000     
 121  124   -0.006000  0.349400   0.00000    
 121  125    0.000000  0.012400   0.00000   
 121  127   -0.020400  0.833800   0.00000    
 121  128   -0.027800  0.309500   0.00000    
 121  129   -0.454500  4.254000   0.00000       
 121  131   -0.218300  1.506600   0.00000        
 121  132   -0.130800  1.381500   0.00000       
 122  123   -0.584000  4.860900   0.00000       
 122  124   -0.000900  0.055200   0.00000      
 122  125   -0.006900  0.158300   0.00000       
 122  131   -0.243300  1.935000   0.00000      
 122  132   -0.018700  0.257200   0.00000      
 122  133   -0.098000  0.982100   0.00000       
 122  143   -0.031200  0.488800   0.00000      
 123  124   -0.223000  1.967000   0.00000     
 123  125   -0.082100  0.606200   0.00000      
 123  131   -0.178300  1.253500   0.00000       
 123  132   -0.135500  1.204100   0.00000      
 124  125   -0.001700  0.094900   0.00000        
 124  128   -1.153000  8.251300   0.00000      
 124  131   -0.106200  0.818500   0.00000    
 124  132   -0.009400  0.161200   0.00000      
 124  133   -0.034200  1.179800   0.00000     
 124  143   -0.007800  0.760700   0.00000    
 125  127   -0.079100  0.985100   0.00000     
 125  128   -0.062000  0.599100   0.00000     
 125  129   -0.421700  3.970200   0.00000      
 125  130   -1.974000  8.485400   0.00000      
 125  131   -0.125100  0.693900   0.00000      
 125  132   -0.053600  0.508600   0.00000       
 127  128   -0.002600  0.124000   0.00000      
 127  129   -0.039200  1.108200   0.00000       
 128  129   -0.001000  0.020700   0.00000      
 128  130   -1.100000  2.992400   0.00000     
 128  131   -1.559000  4.086900   0.00000   
 130  131   -0.002700  0.015400   0.00000     
 130  132   -0.650900  3.031000   0.00000      
 130  144   -0.753200  3.066400   0.00000      
 131  132   -0.003200  0.041100   0.00000   
 131  133   -1.077000  5.528500   0.00000   
 131  143   -0.058800  0.405500   0.00000     
 131  144   -0.002200  0.015100   0.00000     
 132  133   -0.091600  0.822900   0.00000      
 132  143   -0.004900  0.096500   0.00000     
 132  144   -0.110800  0.982700   0.00000       
 133  143   -0.360000  2.630900   0.00000    
 134  131   -0.404200  0.914400   0.00000      
 134  136   -0.069800  0.642800   0.00000       
 134  139   -0.035300  0.166000   0.00000    
 134  141   -0.023000  0.117900   0.00000     
 134  142   -0.026300  0.116700   0.00000        
 134  144   -0.014500  0.043500   0.00000      
 134  145   -0.003400  0.021600   0.00000     
 135   95   -0.344800  3.484500   0.00000       
 135  136   -0.003100  0.017800   0.00000      
 135  138   -0.008400  0.172900   0.00000    
 135  141   -0.129000  0.699300   0.00000    
 136  115   -0.012000  0.085500   0.00000      
 136  116   -1.200000  4.265500   0.00000      
 136  117   -2.969000  9.087500   0.00000       
 136  118   -0.574900  1.620600   0.00000       
 136  138   -0.158100  0.548500   0.00000       
 136  139   -0.005900  0.029300   0.00000      
 136  140   -2.403000  9.378000   0.00000   
 136  141   -0.002600  0.017500   0.00000       
 136  142   -0.046700  0.170900   0.00000      
 136  143   -1.762000  3.454900   0.00000      
 136  145   -0.004900  0.053900   0.00000      
 137  139   -0.018300  0.093600   0.00000      
 137  140   -2.229000  8.022800   0.00000        
 137  145   -0.085200  0.407100   0.00000      
 139  140   -0.005400  0.023900   0.00000      
 139  141   -0.008300  0.046000   0.00000     
 139  142   -0.310200  1.267000   0.00000      
 139  145   -0.000900  0.008000   0.00000       
 140  145   -0.108800  0.480000   0.00000    
 141  115   -0.000700  0.013100   0.00000   
 141  116   -0.156800  0.744800   0.00000     
 141  117   -0.370200  1.382000   0.00000      
 141  118   -0.041400  0.143900   0.00000        
 141  131   -0.233100  0.812900   0.00000      
 141  132   -1.628000  7.093600   0.00000      
 141  142   -0.001800  0.010500   0.00000       
 141  143   -0.070200  0.177800   0.00000   
 141  144   -0.075600  0.244100   0.00000     
 141  145   -0.003800  0.035800   0.00000      
 142  115   -0.016600  0.156300   0.00000      
 142  116   -0.691600  2.630200   0.00000      
 142  117   -0.559600  2.228400   0.00000    
 142  118   -0.018500  0.103700   0.00000    
 142  119   -0.274200  1.861100   0.00000        
 142  120   -0.604300  7.353000   0.00000    
 142  122   -0.258900  2.173200   0.00000    
 142  124   -0.173600  2.134700   0.00000      
 142  125   -1.090000  8.616000   0.00000     
 142  130   -0.360800  1.861800   0.00000    
 142  131   -0.001300  0.015700   0.00000     
 142  132   -0.005500  0.081000   0.00000    
 142  133   -1.636000  9.172500   0.00000    
 142  143   -0.003800  0.018700   0.00000     
 142  144   -0.002000  0.022900   0.00000     
 142  145   -0.073800  0.438000   0.00000      
 143  144   -0.486300  2.328200   0.00000    
 144  145   -0.383500  1.205200   0.00000
   1    3   -0.009000 -0.171800   0.9350   
   1    4   -0.009000 -0.171800   0.9350    
   1    5   -0.008900 -0.169700   0.9350  
   1   33    0.000100  0.006000   0.9350    
   1   93    0.000200  0.013800   1.1036    
   1   93    0.000200  0.013800   1.1036  
   2  113    0.000000  0.014800   1.1052    
   2  114    0.000180  0.014500   1.1052   
   7    8   -0.011200 -0.151600   0.9716    
   7   66    0.000150  0.009700   0.9716    
   7  104    0.000360  0.019000   1.1052    
   7  104    0.000410  0.017400   1.1052  
   9   11   -0.021700 -0.306200   0.9166   
   9   69    0.000400  0.018800   0.9166    
  10   32   -0.027000 -0.304100   0.9166    
  10   69    0.000400  0.018700   0.9166   
  12   13   -0.022300 -0.309900   0.9166    
  12   13   -0.023700 -0.316000   0.9166    
  12   13   -0.023700 -0.316000   0.9166  
  12   72    0.000300  0.018900   0.9166    
  12   72    0.000300  0.019000   0.9166    
  12   72    0.000300  0.019000   0.9166   
  14   15   -0.041500 -0.399600   0.9164    
  14   16   -0.010000 -0.166900   0.9164    
  14   58    0.000200  0.009700   0.9164    
  17   18   -0.318100 -1.315000   0.8708    
  17   19    0.000000 -0.847000   0.8634    
  17   20    0.000000 -0.867600   0.8634    
  17   21   -0.009500 -0.161500   0.8708   
  17   59    0.000100  0.007100   0.8708 
  22   23    0.000000 -0.378700   0.9322  
  22   30    0.000000 -0.306600   0.9532   
  22   78    0.000000  0.026800   0.9532    
  22   83    0.000000  0.034900   0.9322   
  24   76    0.000200  0.008800   0.9898   
  24   77   -0.002300 -0.060300   0.9898   
  25   26   -0.006000 -0.137500   0.9166  
  25   31   -0.008200 -0.164800   0.9166   
  25   73    0.000300  0.017200   0.9166   
  25   74    0.000400  0.017900   0.9166  
  27   28   -0.115300 -0.745300   0.9074    
  27   29   -0.016300 -0.261800   0.9074    
  27   75    0.000160  0.010000   0.9074   
  33  110    0.000240  0.015700   1.1800    
  33  110    0.000230  0.015600   1.1800   
  36   99    0.000800  0.045500   1.1291    
  37   87    0.000930  0.044200   1.0500 
  61   62   -0.036200 -0.260800   1.0500   
  61   62   -0.047200 -0.543800   1.0500 
  61   86    0.001320  0.032000   1.0500   
  61   86    0.001100  0.037000   1.0500   
  61   86    0.001100  0.037000   1.0500  
```


```c title="pmat.c"
% This is called from the main program

% Form the P matrices
k=[1:ngen;(ngen+1):2*ngen];
kk=k(:);

PGmod1=sparse(zeros(nb,ngen));
PGmod1((gen(:,1)),[1:ngen])=diag(ones(1,ngen));
PGmod=[PGmod1 zeros(nb,ngen);zeros(nb,ngen) PGmod1];
PG=PGmod(ll,kk);

clear PGmod1;
clear PGmod;

m=[1:nload;(nload+1):2*nload];
mmm=m(:);

PLmod1=sparse(zeros(nb,nload));
PLmod1((ld(:,1)),[1:nload])=diag(ones(1,nload));
PLmod=[PLmod1 zeros(nb,nload);zeros(nb,nload) PLmod1];
PL=PLmod(ll,mmm);

clear PLmod1;
clear PLmod;

%pack;
```


```data title="power_pss.dat"
60 10 0.05 0.03 0.1 -0.1
```


```c title="power_pss_settings.c"

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        FORMATION OF A MATRIX with PSS

TAGmod=sparse([zeros(ngen) zeros(ngen) zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag21)  zeros(ngen)  diag(Ag23)    diag(Ag24)   diag(Ag25)   diag(Ag26)    zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)]);    
              
 TAG=TAGmod(ii,ii);
 clear TAGmod;
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 TBGmod=sparse([zeros(ngen)      zeros(ngen)
              diag(Bg21rect)   diag(Bg22rect)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen) ]);
           
           
  TBG=TBGmod(ii,kk);
  clear TBGmod; 
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TCGp2=CG;
  
 if(npss~=0)
   TAGp2=([TAG+sparse(zeros(size(EG*PS1'*DPS*PS1*FG))) sparse(zeros(size(EG*PS1'*CPS))); BPS*PS1*FG APS]);
   TBGp2=([TBG; sparse(zeros(size(APS,1),2*ngen))]);
   TEGp2=([EG; sparse(zeros(size(APS,1),ngen))]); %padded zeros to EG suitable.
 else
   TAGp2=TAG;
   TBGp2=TBG;
   TEGp2=EG;
 end
 
if(npss1~=0)  
TAGp2=([TAGp2+sparse(zeros(size(TEGp2*PS2'*DPSS*PS2m*FGd))) sparse(zeros(size(TEGp2*PS2'*CPSS))); BPSS*PS2m*FGd APSS]);
TBGp2=([TBGp2; sparse(zeros(size(APSS,1),2*ngen))]);
TEGp2=([TEGp2; sparse(zeros(size(APSS,1),ngen))]); %padded zeros to EG suitable.
end


TEATd=-TAGp2-TBGp2*PG'*YDQdi*PG*TCGp2; 
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 H2_diag = 2*diag(H);
 TEATD2=H2_diag*TEATd([2:16:16*ngen],:);

APSS2=[];
BPSS2=[];
CPSS2=[];
DPSS2=[];

for k=1:size(power_pss,1)
[NP1,DP1]=series([-Kpower_pss(k)],[Td1_power(k) 1],[Tw_power(k) 0],[Tw_power(k) 1]);
[APSS2d,BPSS2d,CPSS2d,DPSS2d]=tf2ss(NP1,DP1);
[APSS2,BPSS2,CPSS2,DPSS2]=append(APSS2,BPSS2,CPSS2,DPSS2,APSS2d,BPSS2d,CPSS2d,DPSS2d);
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
APSS2=sparse(APSS2);
BPSS2=sparse(BPSS2);
CPSS2=sparse(CPSS2);
DPSS2=sparse(DPSS2);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for k=1:size(power_pss,1)
    for p=1:ngen
        if (power_pss(k,1)==gen(p,1))
            temp(k)=p;
        end
    end
end

PS3=sparse(zeros(size(power_pss,1),ngen));
PS3(1:(size(power_pss,1)),temp)=(eye(size(power_pss,1)));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AG=([AG+TEGp2*PS3'*DPSS2*PS3*TEATD2 TEGp2*PS3'*CPSS2; BPSS2*PS3*TEATD2 APSS2]);
BG=([BG; sparse(zeros(size(APSS2,1),2*ngen))]);
CG=([CG  sparse(zeros(size(APSS2,1),2*ngen)')]);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear C1 D1;
clear A1 B1;
clear temp;
```


```c title="powerflow.c"
%------------------------------------
S_sh = sparse(zeros(nshunt,1));
for i=1:nline
   Ipq(i)=(Vbus(nt(i,1)) - Vbus(nt(i,2)))/(nt(i,3)+j*nt(i,4)) + j*Vbus(nt(i,1))*(nt(i,5)/2);
   Iqp(i)=(Vbus(nt(i,2)) - Vbus(nt(i,1)))/(nt(i,3)+j*nt(i,4)) + j*Vbus(nt(i,2))*(nt(i,5)/2);
end
for i=nline+1:nline+ntrans
   incr=1/(nt(i,3)+j*nt(i,4));
   incr1=incr/nt(i,5);
   incr2=(1-nt(i,5))*incr/(nt(i,5)*nt(i,5));
   incr3=(nt(i,5)-1)*incr/nt(i,5);
   Ipq(i)=(Vbus(nt(i,1)) - Vbus(nt(i,2)))*incr1 + Vbus(nt(i,1))*incr2 ;
   Iqp(i)=(Vbus(nt(i,2)) - Vbus(nt(i,1)))*incr1 + Vbus(nt(i,2))*incr3;  
end
if (nshunt~=0)
 for i=1:nshunt
  incr=shunt(i,2)+j*shunt(i,3);   
  Ish(i)= Vbus(shunt(i,1))*incr;
end
 S_sh = Vbus(shunt(:,1)).*Ish';
end
Spq = Vbus(nt(:,1)).* Ipq';
Sqp = Vbus(nt(:,2)).* Iqp';

S_line_loss = sum(Spq+Sqp) + sum(S_sh);
```


```c title="primemover_settings.c"
%------------------------------------
S_sh = sparse(zeros(nshunt,1));
for i=1:nline
   Ipq(i)=(Vbus(nt(i,1)) - Vbus(nt(i,2)))/(nt(i,3)+j*nt(i,4)) + j*Vbus(nt(i,1))*(nt(i,5)/2);
   Iqp(i)=(Vbus(nt(i,2)) - Vbus(nt(i,1)))/(nt(i,3)+j*nt(i,4)) + j*Vbus(nt(i,2))*(nt(i,5)/2);
end
for i=nline+1:nline+ntrans
   incr=1/(nt(i,3)+j*nt(i,4));
   incr1=incr/nt(i,5);
   incr2=(1-nt(i,5))*incr/(nt(i,5)*nt(i,5));
   incr3=(nt(i,5)-1)*incr/nt(i,5);
   Ipq(i)=(Vbus(nt(i,1)) - Vbus(nt(i,2)))*incr1 + Vbus(nt(i,1))*incr2 ;
   Iqp(i)=(Vbus(nt(i,2)) - Vbus(nt(i,1)))*incr1 + Vbus(nt(i,2))*incr3;  
end
if (nshunt~=0)
 for i=1:nshunt
  incr=shunt(i,2)+j*shunt(i,3);   
  Ish(i)= Vbus(shunt(i,1))*incr;
end
 S_sh = Vbus(shunt(:,1)).*Ish';
end
Spq = Vbus(nt(:,1)).* Ipq';
Sqp = Vbus(nt(:,2)).* Iqp';

S_line_loss = sum(Spq+Sqp) + sum(S_sh);
```


```c title="pss_delPw_signal.c"

%----Delta Power-Omega signal based PSS----
load delPw_pss.dat
%Modfy the delPw_pss to suit the enabled PSS
PSS_dash=~(PSS);
PSS_sm=PSS_dash(gen(:,1));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SL_delPw_pss_dash=~(SL_delPw_pss);
SL_delPw_pss_dash=SL_delPw_pss_dash(gen(:,1));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PS2d=SL_delPw_pss_dash.*PSS_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
npss1=sum(PS2d);
if (npss1~=0&npss1~=size(delPw_pss,1))
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
delPw_pss_temp=delPw_pss;
clear delPw_pss;
k=1;
for i=1:ngen
   if(PS2d(i)~=0)
      x(k)=gen(i,1);
      k=k+1;
   end
end
m=1;
for i=1:size(delPw_pss_temp,1)
   for np=1:npss1 
      if(delPw_pss_temp(i)==x(np))
         delPw_pss(m,:)=delPw_pss_temp(i,:);
         m=m+1;
      end
   end
end
%-------------------------------------------------------------------------------------
clear delPw_pss_temp;
end
[mm order_delPw_pss]=sort([setxor(1:nb,delPw_pss(:,1)) delPw_pss(:,1)]);

Tm0_delPw=TmA0(delPw_pss(:,1)');
Tw1_delPw=delPw_pss(:,2)';
Tw2_delPw=delPw_pss(:,3)';
Tw3_delPw=delPw_pss(:,4)';
Tw4_delPw=delPw_pss(:,5)';
T6_delPw=delPw_pss(:,6)';
T7_delPw=delPw_pss(:,7)';
H_delPw=delPw_pss(:,8)';
Ks3_delPw=delPw_pss(:,9)';
T8_delPw=delPw_pss(:,10)';
T9_delPw=delPw_pss(:,11)';
T1_delPw=delPw_pss(:,12)';
T2_delPw=delPw_pss(:,13)';
T3_delPw=delPw_pss(:,14)';
T4_delPw=delPw_pss(:,15)';
Ks1_delPw=delPw_pss(:,16)';
VS_delPw_max=delPw_pss(:,17)';
VS_delPw_min=delPw_pss(:,18)';
```


```c title="pss_design.c"
%pss_selection;
TAGmod=sparse([zeros(ngen) zeros(ngen) zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              diag(Ag21)  diag(Ag22)  diag(Ag23)    diag(Ag24)   diag(Ag25)   diag(Ag26)    zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)
              zeros(ngen) zeros(ngen) zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)  zeros(ngen)   zeros(ngen)  zeros(ngen)]);    
              
 TAG=TAGmod(ii,ii);
 clear TAGmod;
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 TBGmod=sparse([zeros(ngen)      zeros(ngen)
              diag(Bg21rect)   diag(Bg22rect)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen)
              zeros(ngen)      zeros(ngen) ]);
           
           
  TBG=TBGmod(ii,kk);
  clear TBGmod; 
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TEATd=-TAG-TBG*PG'*YDQdi*PG*CGnp;
  TEAT=TEATd([2:16:16*ngen],:);
  TEAT(:,[1:16:16*ngen])=[]; % eliminating delta for all machines
  TEAT(:,[1:15:15*ngen])=[]; % eliminating Sm for all machines
  rr=input('Enter the generator number for which you want to obtain the angle of GEPS(s):  ');
  for gg=1:1:ngen
     if(gen(gg,1)==rr)
        te=gg;
     end
  end
  TEATR=TEAT(te,:)*2*H(te);%taking the required machine Te data.
  
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TA=Anp; % Anp represents the overall state matrix is obtained in the small_sig programme without PSS. 
  
  TA([1:16:16*ngen],:)=[]; % eliminating delta (in row-wise) for all machines 
  TA([1:15:15*ngen],:)=[]; % eliminating Sm (in row-wise)for all machines
  TA(:,[1:16:16*ngen])=[]; % eliminating delta (in column-wise)for all machines
  TA(:,[1:15:15*ngen])=[]; % eliminating Sm (in column-wise)for all machines
  
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  EGT=EG(:,te);      % EG formed in the main programme is used.
  EGT([1:16:16*ngen],:)=[]; % eliminating delta for all machines
  EGT([1:15:15*ngen],:)=[]; % eliminating Sm for all machines
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for i=1:470
    w(i)=0.05*i;
    Xwd(:,i)=(j*w(i)*eye((16-2)*ngen)-TA)\EGT;
  end
  Te_Vs=TEATR*Xwd;
  
  PHASE1=angle(Te_Vs)*180/pi;
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 disp('Enter 1  :To design single input PSS - Slip signal ')
 disp('      2  :To design double input PSS ')
 disp('      3  :To design single input PSS -Power signal ')
 psstype = input('Enter your choice : ');
 if(psstype==1)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mpssd=1;
while(mpssd==1)
 delete(get(0,'children'))
    figure(1)
 plot( w/2/pi, PHASE1);
W_m=input('Enter the center frequency f_m for the PSS (in Hz) :  ' );
Phi_m=input('Enter the amount of phase lead required (in Degrees):  ');
Kslip_pssf=input('Enter the PSS gain Ks:  ' );
alpha=(1-sin(Phi_m*pi/180))/(1+sin(Phi_m*pi/180));
T1_slipf=1/(sqrt(alpha)*2*pi*W_m);
T2_slipf=alpha*T1_slipf;
if(1/alpha<10)
    fprintf('\n The ratio of T1 to T2 = %5.4f is less than 10\n', 1/alpha);
else
    fprintf('\n The ratio of T1 to T2 = %5.4f is greater than 10\n', 1/alpha);
end
 disp('Enter 1 : Only compensator')
 disp('      2 : Washout only')
 disp('      3 : Washout and measuring ckt.')
 disp('      4 : Washout, measuring ckt. and torsional filter')
 config = input('Enter your choice : ');
 
  if (config==1) 
    TFPSS=Kslip_pssf*((1+(j*w).*T1_slipf)./(1+(j*w).*T2_slipf));
  end

  if (config==2)
     Tw_slipf=input('Enter the value of Tw (in s) for the wash-out circuit [1 - 20]s :  ' );
     TFPSS=Kslip_pssf*(((j*w).*Tw_slipf)./(1+(j*w).*Tw_slipf)).*((1+(j*w).*T1_slipf)./(1+(j*w).*T2_slipf));
  end
   
  if (config==3)
   Tw_slipf=input('Enter the value of Tw (in s) for the wash-out circuit [1 - 20]s :  ' );
   TR_slipf=input('Enter the value of TR (in s) for the measurement circuit [0.02-0.05]s :  ' );
   TFPSS=Kslip_pssf*(1./(1+(j*w).*TR_slipf)).*(((j*w).*Tw_slipf)./(1+(j*w).*Tw_slipf)).*((1+(j*w).*T1_slipf)./(1+(j*w).*T2_slipf));
  end
 
  if (config==4)
   Tw_slipf=input('Enter the value of Tw (in s) for the wash-out circuit [1 - 20]s :  ' );
   TR_slipf=input('Enter the value of TR (in s) for the measurement circuit [0.02-0.05]s :  ' );  
   a0 = input('Enter the value of a0 [570] :');
   a1 = input('Enter the value of a1 [35] :');
   TFPSS=Kslip_pssf*(1./(1+(j*w).*TR_slipf)).*(((j*w).*Tw_slipf)./(1+(j*w).*Tw_slipf)).*((1+(j*w).*T1_slipf)./(1+(j*w).*T2_slipf)).*(a0./(((j*w).*(j*w))+a1.*(j*w)+570));
  end 
 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PHASE2=angle(TFPSS)*180/pi;
PHASEC=PHASE1+PHASE2;
figure(2)
plot(w/2/pi, PHASEC, 'r', w/2/pi, PHASE1, 'b', w/2/pi, PHASE2, 'k');
grid on
clear PHASEC  PHASE2;
mpssd=input('Enter 0 if you are satisfied with the design, otherwise 1 to re-design the pss:  ');
%-------------------------------------------------------------------

end%while loop end if you are satisfied...psstype==1
if (config==2)
fprintf('\n Update the slip_pss.dat for the machine %d \n',rr);
disp('---------------------------------------------------------')
disp('gen.no   Ks      Tw        T1             T2 ')
disp('---------------------------------------------------------')
fprintf('%d \t %d \t %d \t  %5.5f \t %5.5f \n\n', rr, Kslip_pssf, Tw_slipf, T1_slipf, T2_slipf)
disp('---------------------------------------------------------')
disp('Use typical value for TR (0.02 s), a0 =570, and a1 = 35')
disp('---------------------------------------------------------')
end

if (config==3)
fprintf('\n Update the slip_pss.dat for the machine %d \n',rr);
disp('---------------------------------------------------------')
disp('gen.no   TR              Ks      Tw        T1             T2 ')
disp('---------------------------------------------------------')
fprintf('%d \t %5.5f \t %d \t %d \t  %5.5f \t %5.5f \n\n', rr, TR_slipf, Kslip_pssf, Tw_slipf, T1_slipf, T2_slipf)
disp('---------------------------------------------------------')
disp('Use typical value for a0 =570, and a1 = 35')
disp('---------------------------------------------------------')
end

if (config==4)
fprintf('\n Update the slip_pss.dat for the machine %d \n',rr);
disp('-------------------------------------------------------------------------')
disp('gen.no   TR      Ks      Tw        T1         T2       a0       a1')
disp('-------------------------------------------------------------------------')
fprintf('%d \t %5.3f \t %d \t %d \t  %8.5f  %8.5f  %5.2f  %5.2f \n\n', rr, TR_slipf, Kslip_pssf, Tw_slipf, T1_slipf, T2_slipf,a0,a1)
disp('-------------------------------------------------------------------------')
end

end %if(psstype==1) loop end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(psstype==2)   
  mpssd=1;
while(mpssd==1)
 delete(get(0,'children'))
 figure(1)
 plot( w/2/pi, PHASE1);
 disp('Now you are designing compensator Gc1(s)..............')
 W_m=input('Enter the center frequency f_m for  Gc1 (in Hz) :  ' );
 Phi_m=input('Enter the amount of phase lead required (in Degrees):  ');
 KdelPw_pssf=input('Enter the PSS gain Ks:  ' );
 alpha=(1-sin(Phi_m*pi/180))/(1+sin(Phi_m*pi/180));
 T1_delPwf=1/(sqrt(alpha)*2*pi*W_m);
 T2_delPwf=alpha*T1_delPwf;
if(1/alpha<10)
    fprintf('\n The ratio of T1 to T2 = %5.4f is less than 10\n', 1/alpha);
else
    fprintf('\n The ratio of T1 to T2 = %5.4f is greater than 10\n', 1/alpha);
 end
 
disp('Enter y : if compensator Gc2(s) is same as the compensator Gc1(s)')
disp('      n : if compensator Gc2(s) is different from the compensator Gc1(s) ')
compensat3= input('Enter your choice (as a character input y or n): '); 
if ( compensat3=='y') 
 T3_delPwf=T1_delPwf;
 T4_delPwf=T2_delPwf;
else
 W_m=input('Enter the center frequency f_m for Gc2 (in Hz) :  ' );
 Phi_m=input('Enter the amount of phase lead required (in Degrees):  ');
 alpha=(1-sin(Phi_m*pi/180))/(1+sin(Phi_m*pi/180));
 T3_delPwf=1/(sqrt(alpha)*2*pi*W_m);
 T4_delPwf=alpha*T3_delPwf;
if(1/alpha<10)
    fprintf('\n The ratio of T3 to T4 = %5.4f is less than 10\n', 1/alpha);
else
    fprintf('\n The ratio of T3 to T4 = %5.4f is greater than 10\n', 1/alpha);
 end
end % if( compensat3=='y') loop end

 disp('Enter 1 : Compensators Gc1(s)*Gc2(s) only')
 disp('      2 : Compensators with Washout only')
 disp('      3 : All Blocks')
 config2 = input('Enter your choice : ');
 if (config2==1) 
    TFPSS=KdelPw_pssf*((1+(j*w).*T1_delPwf)./(1+(j*w).*T2_delPwf)).*((1+(j*w).*T3_delPwf)./(1+(j*w).*T4_delPwf));
  end 
 if(config2==2)
  Tw_delPwf=input('Enter the value of Tw (in s) for the wash-out circuit [1 - 20]s :  ' );
  TFPSS=KdelPw_pssf.*(((j*w).*Tw_delPwf)./(1+(j*w).*Tw_delPwf)).*((1+(j*w).*T1_delPwf)./(1+(j*w).*T2_delPwf)).*((1+(j*w).*T3_delPwf)./(1+(j*w).*T4_delPwf));
end 
if (config2==3)
    Ks3_delPw=1;
    Tw1_delPwf=input('Enter the value of Tw1 (in s) for the wash-out circuit-1 [1 - 20]s :  ' );
    Tw2_delPwf=input('Enter the value of Tw2 (in s) for the wash-out circuit-2 [1 - 20]s :  ' );
    Tw3_delPwf=input('Enter the value of Tw3 (in s) for the wash-out circuit-3 [1 - 20]s :  ' );
    Tw4_delPwf=input('Enter the value of Tw4 (in s) for the wash-out circuit-4 [1 - 20]s :  ' );
    T6_delPwf=input('Enter the value of T6 (in s) for the slip-path delay [0.01 - 0.05]s :  ' );
    T7_delPwf=input('Enter the value of T7 (in s) for the equivalent integrator [1 - 10]s :  ' );
    T8_delPwf=input('Enter the value of T8 (in s) for the Filter circuit [0 - 0.01]s :  ' );
    T9_delPwf=input('Enter the value of T9 (in s) for the Filter circuit [0.1 - 0.2]s :  ' );
   TFPSS1=(((j*w).*Tw1_delPwf)./(1+(j*w).*Tw1_delPwf)).*(((j*w).*Tw2_delPwf)./(1+(j*w).*Tw2_delPwf)).*(1./(1+(j*w).*T6_delPwf)).*((1+j*w.*T8_delPwf)./(1+(j*w).*T9_delPwf).^2).^4;
   TFPSS2=(((j*w).*Tw3_delPwf)./(1+(j*w).*Tw3_delPwf)).*(((j*w).*Tw4_delPwf)./(1+(j*w).*Tw4_delPwf)).*(((-j*w).*T7_delPwf)./(1+(j*w).*T7_delPwf)).*((((1+j*w.*T8_delPwf)./(1+(j*w).*T9_delPwf).^2).^4).*Ks3_delPw-1);
   TFPSS3=KdelPw_pssf.*((1+(j*w).*T1_delPwf)./(1+(j*w).*T2_delPwf)).*((1+(j*w).*T3_delPwf)./(1+(j*w).*T4_delPwf));
   TFPSS=(TFPSS1+ TFPSS2).*TFPSS3;
end
PHASE2=angle(TFPSS)*180/pi;
PHASEC=PHASE1+PHASE2;
figure(2)
plot(w/2/pi, PHASEC, 'r', w/2/pi, PHASE1, 'b', w/2/pi, PHASE2, 'k');
grid on
clear PHASEC  PHASE2;
mpssd=input('Enter 0 if you are satisfied with the design, otherwise 1 to re-design the pss:  ');
%-------------------------------------------------------------------

end %while loop end if you are satisfied...psstype==2


fprintf('\n Update the delPW_pss.dat for the machine %d \n',rr);
disp('----------------------------------------------------------------------')
disp('gen.no    Ks        T1             T2               T3           T4 ')
disp('-----------------------------------------------------------------------')
fprintf('%d \t %d \t %5.5f \t %5.5f \t  %5.5f \t %5.5f \n\n', rr, KdelPw_pssf, T1_delPwf, T2_delPwf, T3_delPwf, T4_delPwf)
disp('------------------------------------------------------------------------')
end %if(psstype==2) loop end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(psstype==3)
mpssd=1;
while(mpssd==1)
 delete(get(0,'children'))
    figure(1)
    plot( w/2/pi, PHASE1);
    disp('Enter 1 : Plain power type PSS')
    disp('      2 : With measurement delay and washout time constant')
    config_power = input('Enter your choice : ');
    
    if(config_power==1)
       Kpower_pssf=input('Enter the PSS gain Ks:  ' );
       PHASEC=PHASE1+90;
       figure(2)
       plot(w/2/pi, PHASEC, 'r', w/2/pi, PHASE1, 'b');
       Tw_powerf=10;Td_powerf=0.05;
    end
    
 if(config_power==2)
  Kpower_pssf=input('Enter the PSS gain Ks:  ' );
  Tw_powerf=input('Enter the value of Tw (in s) for the wash-out circuit [1 - 20]s :  ' );
  Td_powerf=input('Enter the value of TR (in s) for the measurement delay [0.01 - 0.05]s :  ' );
  TFPOWER= Kpower_pssf.*(((j*w).*Tw_powerf)./(1+(j*w).*Tw_powerf)).*(1./(1+(j*w).*Td_powerf));
  PHASE2=angle(TFPOWER)*180/pi+90; 
  PHASEC=PHASE1+PHASE2;
  figure(2)
  plot(w/2/pi, PHASEC, 'r', w/2/pi, PHASE1, 'b', w/2/pi, PHASE2, 'k');
 end

grid on
clear PHASEC;
mpssd=input('Enter 0 if you are satisfied with the design, otherwise 1 to re-design the pss:  ');
%-------------------------------------------------------------------

end %while loop end if you are satisfied...psstype==3
fprintf('\n Update the power_pss.dat for the machine %d \n',rr);
disp('----------------------------------------------------------------------')
disp('gen.no     Ks             TR            Tw    ')
disp('-----------------------------------------------------------------------')
fprintf('%d \t %5.5f \t %5.5f \t %5.5f  \n\n', rr, Kpower_pssf, Td_powerf,Tw_powerf)
disp('------------------------------------------------------------------------')
end % if end of psstype==3

disp('Press any key to obtain the amplitude response of GEPS(iw)')
pause
figure(3)
plot(w/2/pi,abs(Te_Vs))
   
clear PHASE1;
clear Xwd;
clear Te_Vs;
```


```c title="pss_power_signal.c"

%----Power signal based PSS----
load power_pss.dat
%Modfy the power_pss to suit the enabled PSS
PSS_dash=~(PSS);
PSS_sm=PSS_dash(gen(:,1));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SL_power_pss_dash=~(SL_power_pss);
SL_power_pss_dash=SL_power_pss_dash(gen(:,1));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PS3d=SL_power_pss_dash.*PSS_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
npss2=sum(PS3d);
if (npss2~=0&npss2~=size(power_pss,1))
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
power_pss_temp=power_pss;
clear power_pss;
k=1;
for i=1:ngen
   if(PS3d(i)~=0)
      x(k)=gen(i,1);
      k=k+1;
   end
end
m=1;
for i=1:size(power_pss_temp,1)
   for np=1:npss2 
      if(power_pss_temp(i)==x(np))
         power_pss(m,:)=power_pss_temp(i,:);
         m=m+1;
      end
   end
end
%-------------------------------------------------------------------------------------
clear power_pss_temp;
end
[mm order_power_pss]=sort([setxor(1:nb,power_pss(:,1)) power_pss(:,1)]);
Tm0_power=TmA0(power_pss(:,1)');
Tw_power=power_pss(:,2)';
Td1_power=power_pss(:,3)';
Kpower_pss=power_pss(:,4)';
VS_power_max=power_pss(:,5)';
VS_power_min=power_pss(:,6)';
```


```c title="pss_selection.c"
 YES2=input('Enter 1 to use EIG function, otherwise 0 to use EIGS function:  ');
if(YES2)
   [V1 D]=eig(full(AT));
   V1=sparse(V1);
   W=V1\eye(size(V1));
   PF=V1.*conj(W');
   D=diag(D);
   n_eig=1;
   for k=1:size(AT,1)-3
    if (sum(abs(D(k+[0:3])))>0.1)
     n_eig = n_eig+1;
    end
   end
   D=diag(D(1:n_eig))-0.000001;
else
   n_eig=10;
   f_eig=1;
   fprintf('You  are scanning %d eigenvalues around %5.3f Hz using EIGS function.\n', n_eig,f_eig);
   disp('Please press a key: ')
   pause
   options.tol=1e-12;
   options.disp =0;
   options.maxit = 25;
   [V1 D Flag1]= eigs(AT,n_eig,j*(f_eig*2*pi),options);
   [W1 D Flag2]= eigs((AT)',n_eig,j*(f_eig*2*pi),options);
   VW=conj(V1)'*W1;
   PF=(V1.*W1)/VW;
   EIG=diag(D)-0.000001;
   D=diag(EIG);
   if ((Flag1==1)|(Flag2==1))
      disp('----------------------------NOTE---------------------------------------');
      disp('Determination of eigenvalues with EIGS is NOT accurate. Please use EIG function');
      end
end

SL_lamb=(1:n_eig)';
disp('---------------------------------------------------------');
disp(   'SL_number         Eigenvalue       dampingfactor     frequency(Hz)')
[SL_lamb diag(D) -real(diag(D))./sqrt(real(diag(D)).*real(diag(D))+imag(diag(D)).*imag(diag(D)))  abs(imag(diag(D)))./(2*pi) ]   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('---------------------------------------------------------');
disp('---------------------------------------------------------');
EIG=diag(D);
nodp=0; % to count number of swing modes with damping factor less than 0.05  
disp('swing modes for which dampingfactor is less than 0.05')  
disp('---------------------------------------------------------'); 
disp(   'SL_number          Eigenvalue    dampingfactor       frequency(Hz)')

for dp=1:1:n_eig-1
   if ((abs(imag(EIG(dp)))>1e-3)&(abs(imag(EIG(dp)))~=abs(imag(EIG(dp+1))))) %To filter out Non_Oscillatory modes 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
PF_lamb=PF(:,dp);
PF_lamb_nr = PF_lamb/max(PF_lamb);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
slip_PF_lamb=PF_lamb(1:16:ngen*16);
slip_PF_nr=slip_PF_lamb/max(slip_PF_lamb);
[m_sort_slip  m_ind_slip]=sort(abs(slip_PF_nr));
rev_ind_slip=[ngen:-1:1];
n_ind_slip=m_ind_slip(rev_ind_slip);
n_sort_slip=m_sort_slip(rev_ind_slip);
gen_num=gen(n_ind_slip,1);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
abs_sort=0.05;% Sorting Slips for which Normalised slip_PF_nr  is greater than 0.05
k=1;
r_sort_slip(k)=0;
for i=1:ngen
   if(abs(n_sort_slip(i))>abs_sort)
   r_sort_slip(k)=n_sort_slip(i);
      k=k+1;
   end
end
r_ind_slip=n_ind_slip(1:size(r_sort_slip,2));
par_gen_numd=gen_num(1:size(r_sort_slip,2));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
V=V1(:,dp);
slip_V=V(1:16:ngen*16);
%-----To filter out a mode in which participation of slip is very small----
splip_PF_norm = sparse(zeros(1,nb));
splip_PF_norm(gen(:,1)) = PF_lamb_nr(1:16:ngen*16);
splip_PF_norm_filt = splip_PF_norm(par_gen_numd);

slip_Vpart = 0.3;     %  use this variable to decide the slip part. level
if(sum(abs(splip_PF_norm_filt))> slip_Vpart)
n_slip_Vd=slip_V(r_ind_slip);
n_slip_Vd_nr=n_slip_Vd/max(n_slip_Vd);
[m_sort_slipV  m_ind_slipV]=sort(abs(n_slip_Vd_nr));
rev_ind_slipV=[size(n_slip_Vd,1):-1:1];
n_ind_slipV=m_ind_slipV(rev_ind_slipV);
par_gen_num=par_gen_numd(n_ind_slipV);
n_slip_V=n_slip_Vd(n_ind_slipV);
slip_V_angle=angle(n_slip_V)*180/pi;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
swing_ang=90; %Angle Calculations..
p1=1;
p2=1;
Group1(p1)=par_gen_num(1);
Group2(p2)=0;
a=slip_V_angle(1);
for i=2:size(r_sort_slip,2)
   b=slip_V_angle(i);
   if((abs(a-b)>swing_ang)&(abs(a-b)<(360-swing_ang )))
      Group2(p2)=par_gen_num(i);
      p2=p2+1;
   else
      p1=p1+1;
      Group1(p1)=par_gen_num(i);
   end
end
damp_fac_lim = 0.05;
if(sum(Group2)~=0) % Only swing modes, Filteringout Non_Swingmodes
     if (-real(EIG(dp))./sqrt(real(EIG(dp)).*real(EIG(dp))+imag(EIG(dp)).*imag(EIG(dp)))<damp_fac_lim) 
         %swing modes to which dp is less than 0.05
      [dp   EIG(dp)   -real(EIG(dp))./sqrt(real(EIG(dp)).*real(EIG(dp))+imag(EIG(dp)).*imag(EIG(dp)))   abs(imag(EIG(dp)))./(2*pi)]
      nodp=nodp+1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X1=('');
X3=('');
for i=1:size(r_sort_slip,2)
X1=strvcat(X1,'Slip-');
X2d=par_gen_numd(i,1);
X2=num2str(X2d);
X3=strvcat(X3,X2);
end

Xt=strcat(X1,X3);

slip_PF_sortd=slip_PF_nr(r_ind_slip);

if(npss~=0)
    disp('PSS is Enabled......')
else
    disp('PSS is Disabled.....')
end
disp('---------------------------------------------------------');
disp('State variable     Mag(slip-PF-nr)       angle(slip-PF-nr) in deg.');
disp('---------------------------------------------------------');
 
for i=1:size(r_sort_slip,2)
   fprintf('%s \t\t %12.4f \t\t %10.2f \n',Xt(i,:),abs(slip_PF_sortd(i)),angle(slip_PF_sortd(i))*180/pi);
end 
      %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear Group1 Group2;
clear slip_V_angle;
clear n_slip_Vd;
clear n_slip_V;
clear par_gen_numd;
clear par_gen_num;
clear r_ind_slip;
clear rev_ind;
clear n_ind;
clear p1 p2;
clear r_sort_slip k;
clear gen_num;
clear rev_ind_slipV;
clear Xt;
clear slip_PF_sortd;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

     end  %swing modes to which dp is less than 0.05 loop end

     end  % (Group2(p2)~=0) if loop end..
   
  end % (sum(abs(splip_PF_norm_filt))> slip_Vpart) if loop end
  
end %(imag(D(dp,dp))>1e-3)
end % for loop end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp('---------------------------------------------------------')
if(nodp==0)
    disp('NO swing modes contains damping factor less than 0.05')
else
rpt=1;
while(rpt==1)                                                                                       
disp('Please press a key to obtain the angle of GEPS(s) for the selected machine ')
pause
pss_design %PSS Design...........
rpt=input('Enter 1 if you want to design PSS  for another machine, otherwise 0: ') ;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end  % while loop end............
  %else 
  %disp('All Eigenvalues contained dampingfactor of greater than 0.05,no need of design pss')
end % if loop end  if nodp loop

clear Group1 Group2;
clear slip_V_angle;
clear n_slip_Vd;
clear n_slip_V;
clear par_gen_numd;
clear par_gen_num;
clear r_ind_slip;
clear rev_ind;
clear n_ind;
clear p1 p2;
clear r_sort_slip k;
clear gen_num;
clear rev_ind_slipV;
clear Xt;
clear slip_PF_sortd;
```


```c title="pss_slip_signal.c"
%------Slip signal based PSS----
load slip_pss.dat
%Modfy the slip_pss to suit the enabled PSS
PSS_dash=~(PSS);
PSS_sm=PSS_dash(gen(:,1));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SL_slip_pss_dash=~(SL_slip_pss);
SL_slip_pss_dash=SL_slip_pss_dash(gen(:,1));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PS1d=SL_slip_pss_dash.*PSS_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
npss=sum(PS1d);
if (npss~=0&npss~=size(slip_pss,1))
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
slip_pss_temp=slip_pss;
clear slip_pss;
k=1;
for i=1:ngen
   if(PS1d(i)~=0)
      x(k)=gen(i,1);
      k=k+1;
   end
end
m=1;
for i=1:size(slip_pss_temp,1)
   for np=1:npss 
      if(slip_pss_temp(i)==x(np))
         slip_pss(m,:)=slip_pss_temp(i,:);
         m=m+1;
      end
   end
end
%---------------------------------------------------------
clear slip_pss_temp;
end


[mm order_slip_pss]=sort([setxor(1:nb,slip_pss(:,1));slip_pss(:,1)]);
Kslip_pss=slip_pss(:,2)';
Td1_slip=slip_pss(:,3)';
Tw_slip=slip_pss(:,4)';
T1_slip=slip_pss(:,5)';
T2_slip=slip_pss(:,6)';
VS_slip_max=slip_pss(:,7)';
VS_slip_min=slip_pss(:,8)';
a0_slip=slip_pss(:,9)';
a1_slip=slip_pss(:,10)';
TRF_slip=slip_pss(:,11)';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Globalization of all parameters..
Kslip_pss_g=sparse(zeros(1,nb));
Kslip_pss_g(slip_pss(:,1))=Kslip_pss;
Td1_slip_g=sparse(zeros(1,nb));
Td1_slip_g(slip_pss(:,1))=Td1_slip;
Tw_slip_g=sparse(zeros(1,nb));
Tw_slip_g(slip_pss(:,1))=Tw_slip;
T1_slip_g=sparse(zeros(1,nb));
T1_slip_g(slip_pss(:,1))=T1_slip;
T2_slip_g=sparse(zeros(1,nb));
T2_slip_g(slip_pss(:,1))=T2_slip;
a0_slip_g=sparse(zeros(1,nb));
a0_slip_g(slip_pss(:,1))=a0_slip;
a1_slip_g=sparse(zeros(1,nb));
a1_slip_g(slip_pss(:,1))=a1_slip;
TRF_slip_g=sparse(zeros(1,nb));
TRF_slip_g(slip_pss(:,1))=TRF_slip;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xpss_nt=[];
xpss_t=[];
gnt=1;
gt=1;
for i=1:size(slip_pss,1)
   if(slip_pss(i,11)==1)
      xpss_nt(gnt)=slip_pss(i,1);
      gnt=gnt+1;
   else
      xpss_t(gt)=slip_pss(i,1);
      gt=gt+1;
   end
end

%-To enable or disable measurement time delay--set it '0' for Enable and '1' to disable---
Tmd_slip_nt = 1;  %--- Without Torsional Filter 
Tmd_slip_t = 0;   %--- With Torsional Filter 

Tmd_slip = zeros(1,nb);
Tmd_slip(xpss_nt) = Tmd_slip_nt;
Tmd_slip(xpss_t) = Tmd_slip_t;
Tmd_slip = Tmd_slip(slip_pss(:,1));
```


```data title="pvpq.dat"
  60  1.137        0.51
  67  1.0900       14.86     
  91  1.0000       0.64      
  96  1.0000       0.60     
  97  0.9670       1.40     
  98  0.9700       4.26       
 103  1.0000       1.3500     
 108  1.0140       8.0000     
 109  0.9150       0.5200     
  79  1.0520       2.5020   
  80  1.0690       0.4700    
  82  0.9750       0.7000    
  89  1.0660       6.7300   
  90  0.9500       0.2200   
  93  1.0000       7.0000   
  94  1.0200       3.0000   
  95  0.9200       1.3100   
  99  1.0000       2.0000   
 101  1.0390       3.1090  
 102  1.0190       20.4000  
 104  1.0059       20.0000  
 105  1.0070       16.2000  
 106  1.0050       10.8000  
 110  1.0000        7.0000  
 111  1.0000       20.0000  
 112  1.0370        3.0000  
 115  1.0490       24.9300  
 116  1.0430       27.1300  
 117  1.0300       26.2700  
 118  1.0100       42.2000   
 119  1.0130       89.5400 
 121  1.0460       29.9700 
 122  1.0000       10.0900   
 124  1.0000       30.0500   
 128  1.0250      129.6300  
 130  1.0570       59.3700 
 131  1.0420      283.0000 
 132  1.0420       30.9500   
 134  1.0440      206.2600 
 135  1.1070       59.8200 
 136  1.0830      519.5000  
 137  1.0640      120.6800 
 139  1.0400      568.3400 
 140  1.0500      231.2300   
 141  1.0530      379.1100  
 142  1.1550      244.4900  
 143  1.0310       52.5400   
 144  0.9970      113.9700  
 145  1.0520      141.1862  
  34  0.4505        0.4656     
  35  0.4919        0.2753     
  51  0.5845        0.2844    
  58  0.7630       -0.1080      
  66  1.0220        0.2670    
  68  0.00         -0.0741    
  70  0.00          0.5663      
  71  0.00         -0.2120     
  74  0.8190        0.4370     
  78  0.8900        0.2680   
  81  0.8220       -0.9310  
  84  0.2430        0.0820       
  85  0.2740        0.0030    
  88  0.6900        0.2090     
  92  0.00          0.3102     
 107 -0.1750       -0.1280  
 120 -4.0800        1.7510  
 123 -0.8400       -0.1900     
 125 -7.1200       -3.1900     
 126 -3.3300       -1.6000      
 127 -5.4600       -0.7200     
 129 -4.8200       -1.2200     
 133 -0.8300       -0.3630     
 138 -3.6300       -1.8800  
  79  0.0910        0.0300    
  80  0.1710        0.0500       
  82  0.0210        0.0110       
  89  0.0060        0.0020    
  90  0.0460        0.0150       
  93  1.0040        0.7320      
  94  0.1540        0.0760    
  95  0.0670        0.0220     
  99  0.1046        0.0523    
 101  0.1780        0.0450    
 102  0.3760        0.0920    
 104  0.3020       0.0760   
 105  0.9600       1.6740   
 106  0.6400       0.1600   
 110  1.0040       0.7320     
 111  0.6040      11.6600   
 112  0.1860       0.0460    
 115  6.8350       1.8470   
 116  7.9260       3.1550  
 117  4.8530       0.7140  
 118  6.5190       3.2840 
 119 20.9400      37.7400    
 121  2.3770      -0.1730   
 122  0.2920       0.0700     
 124  0.9410       7.8030    
 128 40.7500       7.0350    
 130 43.2800       9.4430    
 131 218.4000     43.2000  
 132  4.9190       1.1020     
 134 223.0900     74.0200  
 135  42.9800     12.6400   
 136 529.5100    135.5200   
 137 129.4600     26.0800  
 139 577.1800    139.3600   
 140 247.7500     66.7600   
 141 327.9900    113.6100   
 142 177.3700     39.3400    
 143  46.7200     17.0900    
 144  96.0200     22.0300    
 145  91.7300     15.5500    
```


```c title="r_eig_plot.c"
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X5=('');
    for i=1:size(n_slip_V,1)
         X5=strvcat(X5,'SG-');
     end
  XGN=num2str(par_gen_num);     
  XGEN=strcat(X5,XGN);
   
  %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  col=['r','b','m','k','g','c'];
  num_r_eig=6;
  nplts=fix(size(n_slip_V,1)/num_r_eig);
  delete(get(0,'children'))
  for p=1:nplts+1
     figure(p)
     cl=1;
     for k=1:num_r_eig
        if ((size(n_slip_V,1)-num_r_eig*(p-1))>=k)
        compass(n_slip_V(k+num_r_eig*(p-1)),col(cl));
        gtext({XGEN(k+num_r_eig*(p-1),:)});
            if(cl==1)
                hold
              end
             cl=cl+1;
          end
       end
    end  
```


```c title="Reheat_turbine.c"
%Speed governer with reheat type turbine

load turb_rhst.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mm order_rhst]=sort([setxor(1:nb,turb_rhst(:,1)) turb_rhst(:,1)]);
T1_rht=turb_rhst(:,2)';
T2_rht=turb_rhst(:,3)';
T3_rht=turb_rhst(:,4)';
sig_rht=turb_rhst(:,5)';
Pgv_rht_max=turb_rhst(:,6)';
Pgv_rht_min=turb_rhst(:,7)';
TCH_rht=turb_rhst(:,8)';
TRH_rht=turb_rhst(:,9)';
TCO_rht=turb_rhst(:,10)';
FHP_rht=turb_rhst(:,11)';
FIP_rht=turb_rhst(:,12)';
FLP_rht=turb_rhst(:,13)';

K_rht=1./sig_rht;
Pgv_rht0=TmA0(turb_rhst(:,1));
HP_rht0=Pgv_rht0;
IP_rht0=Pgv_rht0;
LP_rht0=Pgv_rht0;
```


```data title="Report.dat"
Iter. No = 1 Max. Real power mismatch at bus = 131, Max. mismatch = 6.62271308
Iter. No = 1 Max. Reactive power mismatch at bus = 125, Max. mismatch = 2.11343348
Iter. No = 2 Max. Real power mismatch at bus = 124, Max. mismatch = 1.29028696
Iter. No = 2 Max. Reactive power mismatch at bus = 120, Max. mismatch = 0.34304918
Iter. No = 3 Max. Real power mismatch at bus = 17, Max. mismatch = 0.53093485
Iter. No = 3 Max. Reactive power mismatch at bus = 120, Max. mismatch = 0.02379278
Iter. No = 4 Max. Real power mismatch at bus = 17, Max. mismatch = 0.12275449
Iter. No = 4 Max. Reactive power mismatch at bus = 120, Max. mismatch = 0.00559679
Iter. No = 5 Max. Real power mismatch at bus = 139, Max. mismatch = 0.02691164
Iter. No = 5 Max. Reactive power mismatch at bus = 59, Max. mismatch = 0.00375592
Iter. No = 6 Max. Real power mismatch at bus = 139, Max. mismatch = 0.01063162
Iter. No = 7 Max. Real power mismatch at bus = 139, Max. mismatch = 0.00422150
Iter. No = 6 Max. Reactive power mismatch at bus = 59, Max. mismatch = 0.00128057
Iter. No = 8 Max. Real power mismatch at bus = 139, Max. mismatch = 0.00167527
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Converged Real power mismatch iteration No = 9, Max.mismatch = 0.00066652
Converged Reactive power mismatch iteration No = 7, Max. mismatch = 0.00011477
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Loadflow results:
---------------------------------------------------------------------------------------

Bus No     VbO          thetaO           PGO             QGO         PLO         QLO 
---------------------------------------------------------------------------------------
  1 	     1.080998 	   -5.029944	    0.000000	    0.000000	 0.000000	 0.000000 
  2 	     1.080853 	   -5.096552	    0.000000	    0.000000	 0.000000	 0.000000 
  3 	     1.101531 	   -4.732669	    0.000000	    0.000000	 0.000000	 0.000000 
  4 	     1.101531 	   -4.732669	    0.000000	    0.000000	 0.000000	 0.000000 
  5 	     1.101814 	   -4.732477	    0.000000	    0.000000	 0.000000	 0.000000 
  6 	     1.043300 	   -8.545121	    0.000000	    0.000000	 0.000000	 0.000000 
  7 	     1.076301 	    2.501586	    0.000000	    0.000000	 0.000000	 0.000000 
  8 	     1.113652 	    0.440080	    0.000000	    0.000000	 0.000000	 0.000000 
  9 	     1.039591 	   -8.760132	    0.000000	    0.000000	 0.000000	 0.000000 
 10 	     1.039568 	   -8.760931	    0.000000	    0.000000	 0.000000	 0.000000 
 11 	     1.093670 	  -11.365751	    0.000000	    0.000000	 0.000000	 0.000000 
 12 	     1.038865 	   -9.477825	    0.000000	    0.000000	 0.000000	 0.000000 
 13 	     1.098189 	  -12.138574	    0.000000	    0.000000	 0.000000	 0.000000 
 14 	     1.038524 	   -9.885644	    0.000000	    0.000000	 0.000000	 0.000000 
 15 	     1.068304 	  -10.519430	    0.000000	    0.000000	 0.000000	 0.000000 
 16 	     1.068638 	  -10.572840	    0.000000	    0.000000	 0.000000	 0.000000 
 17 	     1.001236 	  -10.148720	    0.000000	    0.000000	 0.000000	 0.000000 
 18 	     1.074645 	  -11.587839	    0.000000	    0.000000	 0.000000	 0.000000 
 19 	     1.070821 	  -11.667798	    0.000000	    0.000000	 0.000000	 0.000000 
 20 	     1.113056 	  -11.666558	    0.000000	    0.000000	 0.000000	 0.000000 
 21 	     1.108556 	  -11.949572	    0.000000	    0.000000	 0.000000	 0.000000 
 22 	     1.031053 	   -4.591935	    0.000000	    0.000000	 0.000000	 0.000000 
 23 	     1.097860 	   -6.217242	    0.000000	    0.000000	 0.000000	 0.000000 
 24 	     1.027216 	    1.596959	    0.000000	    0.000000	 0.000000	 0.000000 
 25 	     1.037961 	  -10.579214	    0.000000	    0.000000	 0.000000	 0.000000 
 26 	     1.089437 	  -12.082964	    0.000000	    0.000000	 0.000000	 0.000000 
 27 	     1.038853 	  -13.782173	    0.000000	    0.000000	 0.000000	 0.000000 
 28 	     1.076213 	  -15.987627	    0.000000	    0.000000	 0.000000	 0.000000 
 29 	     1.074609 	  -16.150450	    0.000000	    0.000000	 0.000000	 0.000000 
 30 	     1.073074 	   -6.057738	    0.000000	    0.000000	 0.000000	 0.000000 
 31 	     1.090534 	  -12.522370	    0.000000	    0.000000	 0.000000	 0.000000 
 32 	     1.093709 	  -11.373291	    0.000000	    0.000000	 0.000000	 0.000000 
 33 	     1.139245 	   -4.770083	    0.000000	    0.000000	 0.000000	 0.000000 
 34 	     1.138711 	   -4.710749	    0.000000	    0.000000	 0.450500	 0.465600 
 35 	     1.139002 	   -4.788903	    0.000000	    0.000000	 0.491900	 0.275300 
 36 	     1.138542 	   -4.523437	    0.000000	    0.000000	 0.000000	 0.000000 
 37 	     1.123516 	   -6.940806	    0.000000	    0.000000	 0.000000	 0.000000 
 38 	     1.130577 	   -6.001223	    0.000000	    0.000000	 0.000000	 0.000000 
 39 	     1.126959 	   -8.624708	    0.000000	    0.000000	 0.000000	 0.000000 
 40 	     1.126939 	   -8.628814	    0.000000	    0.000000	 0.000000	 0.000000 
 41 	     1.118840 	  -11.140431	    0.000000	    0.000000	 0.000000	 0.000000 
 42 	     1.118813 	  -11.154647	    0.000000	    0.000000	 0.000000	 0.000000 
 43 	     1.118947 	  -11.112560	    0.000000	    0.000000	 0.000000	 0.000000 
 44 	     1.118920 	  -11.126533	    0.000000	    0.000000	 0.000000	 0.000000 
 45 	     1.117292 	  -12.123368	    0.000000	    0.000000	 0.000000	 0.000000 
 46 	     1.117299 	  -12.118313	    0.000000	    0.000000	 0.000000	 0.000000 
 47 	     1.127532 	   -7.433984	    0.000000	    0.000000	 0.000000	 0.000000 
 48 	     1.127833 	   -7.414530	    0.000000	    0.000000	 0.000000	 0.000000 
 49 	     1.127903 	   -7.405086	    0.000000	    0.000000	 0.000000	 0.000000 
 50 	     1.127608 	   -7.423941	    0.000000	    0.000000	 0.000000	 0.000000 
 51 	     1.112399 	  -10.867749	    0.000000	    0.000000	 0.584500	 0.284400 
 52 	     1.111781 	  -11.840752	    0.000000	    0.000000	 0.000000	 0.000000 
 53 	     1.111781 	  -11.842287	    0.000000	    0.000000	 0.000000	 0.000000 
 54 	     1.113112 	  -12.494277	    0.000000	    0.000000	 0.000000	 0.000000 
 55 	     1.113120 	  -12.494304	    0.000000	    0.000000	 0.000000	 0.000000 
 56 	     1.107153 	  -10.650914	    0.000000	    0.000000	 0.000000	 0.000000 
 57 	     1.107201 	  -10.652361	    0.000000	    0.000000	 0.000000	 0.000000 
 58 	     1.106654 	  -10.471996	    0.000000	    0.000000	 0.763000	-0.108000 
 59 	     1.116455 	  -11.550636	    0.000000	    0.000000	 0.000000	 0.000000 
 60 	     1.137000 	   -7.075880	    0.510000	    0.329164	 0.000000	 0.000000 
 61 	     1.114429 	  -12.599427	    0.000000	    0.000000	 0.000000	 0.000000 
 62 	     1.056600 	  -15.177319	    0.000000	    0.000000	 0.000000	 0.000000 
 63 	     1.110920 	  -14.688829	    0.000000	    0.000000	 0.000000	 0.000000 
 64 	     1.098000 	   -9.996327	    0.000000	    0.000000	 0.000000	 0.000000 
 65 	     1.097975 	   -9.998114	    0.000000	    0.000000	 0.000000	 0.000000 
 66 	     1.112887 	    0.609349	    0.000000	    0.000000	 1.022000	 0.267000 
 67 	     1.090000 	   -6.366411	   14.860000	    2.852061	 0.000000	 0.000000 
 68 	     1.208597 	  -31.694046	    0.000000	    0.000000	 0.000000	-0.074100 
 69 	     1.096753 	  -11.125235	    0.000000	    0.000000	 0.000000	 0.000000 
 70 	     0.999829 	  -14.875302	    0.000000	    0.000000	 0.000000	 0.566300 
 71 	     1.027483 	  -14.968905	    0.000000	    0.000000	 0.000000	-0.212000 
 72 	     1.100723 	  -11.902767	    0.000000	    0.000000	 0.000000	 0.000000 
 73 	     1.097523 	  -11.768184	    0.000000	    0.000000	 0.000000	 0.000000 
 74 	     1.097266 	  -12.169410	    0.000000	    0.000000	 0.819000	 0.437000 
 75 	     1.117872 	  -15.896103	    0.000000	    0.000000	 0.000000	 0.000000 
 76 	     1.020878 	    4.828444	    0.000000	    0.000000	 0.000000	 0.000000 
 77 	     0.987963 	    6.013543	    0.000000	    0.000000	 0.000000	 0.000000 
 78 	     1.073979 	   -5.896450	    0.000000	    0.000000	 0.890000	 0.268000 
 79 	     1.052000 	  -10.218608	    2.502000	   -0.159510	 0.091000	 0.030000 
 80 	     1.069000 	   -8.918395	    0.470000	   -0.150649	 0.171000	 0.050000 
 81 	     1.130389 	  -26.573264	    0.000000	    0.000000	 0.822000	-0.931000 
 82 	     0.975000 	  -19.372126	    0.700000	    0.171541	 0.021000	 0.011000 
 83 	     1.098486 	   -6.088671	    0.000000	    0.000000	 0.000000	 0.000000 
 84 	     1.115578 	  -10.147721	    0.000000	    0.000000	 0.243000	 0.082000 
 85 	     1.116493 	  -13.754570	    0.000000	    0.000000	 0.274000	 0.003000 
 86 	     1.056689 	  -14.718246	    0.000000	    0.000000	 0.000000	 0.000000 
 87 	     1.065151 	   -7.879955	    0.000000	    0.000000	 0.000000	 0.000000 
 88 	     1.109420 	   -9.053508	    0.000000	    0.000000	 0.690000	 0.209000 
 89 	     1.066000 	    2.975973	    6.730000	    1.363931	 0.006000	 0.002000 
 90 	     0.950000 	   -8.062662	    0.220000	   -0.038652	 0.046000	 0.015000 
 91 	     1.000000 	   -9.985229	    0.640000	   -0.015385	 0.000000	 0.000000 
 92 	     0.956120 	  -13.459889	    0.000000	    0.000000	 0.000000	 0.310200 
 93 	     1.000000 	   -2.626796	    7.000000	    3.738139	 1.004000	 0.732000 
 94 	     1.020000 	   -1.450779	    3.000000	    0.190526	 0.154000	 0.076000 
 95 	     0.920000 	   18.172369	    1.310000	    0.101211	 0.067000	 0.022000 
 96 	     1.000000 	   -9.685438	    0.600000	    0.211089	 0.000000	 0.000000 
 97 	     0.967000 	   -5.053482	    1.400000	    0.456264	 0.000000	 0.000000 
 98 	     0.970000 	    4.479341	    4.260000	   -0.327207	 0.000000	 0.000000 
 99 	     1.000000 	    0.388998	    2.000000	   -0.083408	 0.104600	 0.052300 
100 	     1.014000 	    0.000000	    1.701011	    0.587239	 0.000000	 0.000000 
101 	     1.039000 	   -6.798184	    3.109000	    1.486615	 0.178000	 0.045000 
102 	     1.019000 	   -5.473517	   20.400000	    4.889027	 0.376000	 0.092000 
103 	     1.000000 	    0.807694	    1.350000	    0.049585	 0.000000	 0.000000 
104 	     1.005900 	   12.967922	   20.000000	    4.999287	 0.302000	 0.076000 
105 	     1.007000 	   -3.503571	   16.200000	    3.883489	 0.960000	 1.674000 
106 	     1.005000 	   -3.458237	   10.800000	    2.093641	 0.640000	 0.160000 
107 	     1.021077 	  -14.282204	    0.000000	    0.000000	-0.175000	-0.128000 
108 	     1.014000 	  -14.739737	    8.000000	    0.772871	 0.000000	 0.000000 
109 	     0.915000 	  -19.163968	    0.520000	   -0.155512	 0.000000	 0.000000 
110 	     1.000000 	   -2.015779	    7.000000	    5.198349	 1.004000	 0.732000 
111 	     1.000000 	    7.262588	   20.000000	    5.637547	 0.604000	11.660000 
112 	     1.037000 	   -6.972694	    3.000000	    1.401114	 0.186000	 0.046000 
113 	     0.977971 	   -5.096552	    0.000000	    0.000000	 0.000000	 0.000000 
114 	     0.977971 	   -5.096552	    0.000000	    0.000000	 0.000000	 0.000000 
115 	     1.049000 	  -16.320584	   24.930000	    1.427262	 6.835000	 1.847000 
116 	     1.043000 	  -17.571642	   27.130000	    6.318378	 7.926000	 3.155000 
117 	     1.030000 	  -16.033238	   26.270000	    2.585555	 4.853000	 0.714000 
118 	     1.010000 	  -18.502919	   42.200000	    6.603727	 6.519000	 3.284000 
119 	     1.013000 	  -60.121676	   89.540000	   47.485070	20.940000	37.740000 
120 	     1.033103 	  -52.314276	    0.000000	    0.000000	-4.080000	 1.751000 
121 	     1.046000 	  -20.905301	   29.970000	   -1.601911	 2.377000	-0.173000 
122 	     1.000000 	   -3.498467	   10.090000	    1.740455	 0.292000	 0.070000 
123 	     1.017115 	  -33.831875	    0.000000	    0.000000	-0.840000	-0.190000 
124 	     1.000000 	   -2.594616	   30.050000	    5.692040	 0.941000	 7.803000 
125 	     1.008381 	  -33.303513	    0.000000	    0.000000	-7.120000	-3.190000 
126 	     1.052383 	  -74.609654	    0.000000	    0.000000	-3.330000	-1.600000 
127 	     1.006960 	  -37.108854	    0.000000	    0.000000	-5.460000	-0.720000 
128 	     1.025000 	  -40.415264	  129.630000	   26.108349	40.750000	 7.035000 
129 	     0.980186 	  -73.782875	    0.000000	    0.000000	-4.820000	-1.220000 
130 	     1.057000 	  -52.574928	   59.370000	   18.349495	43.280000	 9.443000 
131 	     1.042000 	  -25.026583	  283.000000	   74.730319	218.400000	43.200000 
132 	     1.042000 	   -7.951271	   30.950000	    6.334228	 4.919000	 1.102000 
133 	     1.092218 	  -12.308207	    0.000000	    0.000000	-0.830000	-0.363000 
134 	     1.044000 	  -11.531059	  206.260000	   74.021400	223.090000	74.020000 
135 	     1.107000 	   28.334386	   59.820000	   15.648412	42.980000	12.640000 
136 	     1.083000 	    3.677557	  519.500000	  144.535002	529.510000	135.520000 
137 	     1.064000 	  -73.438423	  120.680000	   34.506890	129.460000	26.080000 
138 	     1.113793 	   11.301156	    0.000000	    0.000000	-3.630000	-1.880000 
139 	     1.040000 	  -11.267381	  568.340000	  158.495909	577.180000	139.360000 
140 	     1.050000 	  -26.874223	  231.230000	   67.104727	247.750000	66.760000 
141 	     1.053000 	   -9.830327	  379.110000	  116.694923	327.990000	113.610000 
142 	     1.155000 	  -11.441815	  244.490000	   54.961354	177.370000	39.340000 
143 	     1.031000 	  -14.373714	   52.540000	   21.586259	46.720000	17.090000 
144 	     0.997000 	   -9.287162	  113.970000	   26.868461	96.020000	22.030000 
145 	     1.052000 	    4.309572	  141.186200	   29.871238	91.730000	15.550000 
---------------------------------------------------------------------------------------
Line flows:
----------------------------------------------------------------------------------
                       Line flows                                Line flows 
                   _____________________                   _______________________
 From   To         P-flow         Q-flow    From   To      P-flow           Q-flow
----------------------------------------------------------------------------------
  1 	   2	       1.7028	       0.0951	   2	    1	      -1.7028	      -0.1670
  1 	   2	       1.7028	       0.0951	   2	    1	      -1.7028	      -0.1670
  1 	   6	       3.4691	       0.3392	   6	    1	      -3.4441	      -2.7555
  2 	   6	       3.4055	       0.3340	   6	    2	      -3.3813	      -2.7579
  3 	  33	       0.0201	      -1.8800	  33	    3	      -0.0195	       1.9444
  4 	  33	       0.0201	      -1.8800	  33	    4	      -0.0195	       1.9444
  5 	  33	       0.0204	      -1.8834	  33	    5	      -0.0198	       1.9474
  6 	   7	     -15.4370	      -0.3449	   7	    6	      15.7197	       1.7446
  6 	   9	       2.5859	       1.9417	   9	    6	      -2.5843	      -2.1148
  6 	  10	       2.5960	       1.9550	  10	    6	      -2.5944	      -2.1279
  6 	  12	       8.5403	       0.9808	  12	    6	      -8.5265	      -1.7872
  6 	  12	       8.5403	       0.9808	  12	    6	      -8.5265	      -1.7872
  8 	  66	      -0.1223	       0.0295	  66	    8	       0.1223	      -0.0291
  8 	  66	      -0.1653	       0.0403	  66	    8	       0.1653	      -0.0398
 11 	  69	      -0.1932	      -0.1268	  69	   11	       0.1932	       0.1280
 12 	  14	       0.8390	      -0.5083	  14	   12	      -0.8383	      -0.4088
 12 	  14	       0.8390	      -0.5083	  14	   12	      -0.8383	      -0.4088
 12 	  25	       3.7554	      -0.4785	  25	   12	      -3.7487	      -0.1234
 12 	  25	       3.7554	      -0.4785	  25	   12	      -3.7487	      -0.1234
 13 	  72	      -0.1922	      -0.1052	  72	   13	       0.1922	       0.1062
 13 	  72	      -0.1911	      -0.1037	  72	   13	       0.1911	       0.1047
 13 	  72	      -0.1922	      -0.1052	  72	   13	       0.1922	       0.1062
 14 	  17	       0.2257	      -0.8303	  17	   14	      -0.2221	      -2.7299
 14 	  17	       0.2292	      -0.8278	  17	   14	      -0.2256	      -2.7254
 15 	  58	      -0.0510	      -1.6062	  58	   15	       0.0514	       1.6639
 16 	  58	      -0.1114	      -1.8455	  58	   16	       0.1120	       1.9114
 17 	  22	      -3.6716	      -1.9160	  22	   17	       3.7031	      -0.4091
 18 	  59	      -0.0363	      -1.5075	  59	   18	       0.0367	       1.5661
 19 	  59	      -0.0389	      -0.7768	  59	   19	       0.0389	       0.8100
 20 	  59	      -0.0394	      -0.0592	  59	   20	       0.0394	       0.0595
 21 	  59	      -0.2635	      -0.2636	  59	   21	       0.2636	       0.2673
 22 	  24	      -5.4115	      -0.1073	  24	   22	       5.4605	      -1.3835
 23 	  83	      -0.0456	      -0.0112	  83	   23	       0.0456	       0.0113
 23 	  83	      -0.0454	      -0.0112	  83	   23	       0.0454	       0.0113
 25 	  27	       2.2506	      -1.8095	  27	   25	      -2.2397	      -1.3544
 25 	  27	       2.2506	      -1.8095	  27	   25	      -2.2397	      -1.3544
 26 	  73	      -0.2497	      -0.3264	  73	   26	       0.2497	       0.3302
 28 	  75	      -0.0769	      -1.5454	  75	   28	       0.0773	       1.6053
 29 	  75	      -0.2111	      -1.7263	  75	   29	       0.2116	       1.7967
 30 	  78	      -0.0968	      -0.0289	  78	   30	       0.0968	       0.0292
 31 	  74	      -0.2670	      -0.2595	  74	   31	       0.2670	       0.2627
 32 	  69	      -0.1969	      -0.1238	  69	   32	       0.1969	       0.1250
 33 	  34	      -1.4412	       0.7721	  34	   33	       1.4413	      -0.7711
 33 	  35	       0.4919	       0.2747	  35	   33	      -0.4919	      -0.2753
 33 	  37	       0.7092	       0.0941	  37	   33	      -0.7051	      -0.2081
 33 	  38	       0.4119	       0.0156	  38	   33	      -0.4106	      -0.1491
 33 	  39	       1.2457	       0.0250	  39	   33	      -1.2355	      -0.0702
 33 	  40	       1.2488	       0.0255	  40	   33	      -1.2386	      -0.0701
 33 	  49	       1.2155	       0.1011	  49	   33	      -1.2090	      -0.1441
 33 	  50	       1.2244	       0.1073	  50	   33	      -1.2178	      -0.1493
 34 	  36	      -1.8918	       0.3055	  36	   34	       1.8925	      -0.3000
 37 	  88	       0.2801	       0.0958	  88	   37	      -0.2799	      -0.0843
 38 	  88	       0.4106	       0.1491	  88	   38	      -0.4101	      -0.1247
 39 	  43	       1.1144	       0.0257	  43	   39	      -1.1085	      -0.0669
 39 	  84	       0.1211	       0.0445	  84	   39	      -0.1210	      -0.0408
 40 	  44	       1.1165	       0.0253	  44	   40	      -1.1105	      -0.0664
 40 	  84	       0.1221	       0.0448	  84	   40	      -0.1220	      -0.0412
 41 	  42	       0.0021	       0.0002	  42	   41	      -0.0021	      -0.0002
 41 	  43	      -0.6780	      -0.1254	  43	   41	       0.6780	       0.1250
 42 	  44	      -0.6839	      -0.1250	  44	   42	       0.6839	       0.1246
 43 	  46	       0.4305	      -0.0581	  46	   43	      -0.4296	      -0.0259
 44 	  45	       0.4266	      -0.0582	  45	   44	      -0.4257	      -0.0260
 45 	  61	       0.2892	       0.0206	  61	   45	      -0.2889	      -0.0835
 45 	  85	       0.1366	       0.0054	  85	   45	      -0.1366	      -0.0015
 46 	  61	       0.2921	       0.0205	  61	   46	      -0.2918	      -0.0834
 46 	  85	       0.1374	       0.0054	  85	   46	      -0.1374	      -0.0015
 47 	  48	      -0.0018	      -0.0016	  48	   47	       0.0018	       0.0016
 47 	  50	      -0.2487	      -0.0938	  50	   47	       0.2487	       0.0931
 47 	  87	       0.0572	       0.1636	  87	   47	      -0.0553	      -0.1542
 48 	  49	      -0.2339	      -0.0847	  49	   48	       0.2339	       0.0840
 48 	  87	       0.0566	       0.1493	  87	   48	      -0.0546	      -0.1406
 49 	  51	       0.9751	       0.0601	  51	   49	      -0.9683	      -0.1565
 50 	  51	       0.9692	       0.0562	  51	   50	      -0.9624	      -0.1533
 51 	  52	       0.7479	      -0.0755	  52	   51	      -0.7466	       0.0305
 51 	  53	       0.7491	      -0.0756	  53	   51	      -0.7477	       0.0307
 51 	  56	      -0.0756	       0.0888	  56	   51	       0.0758	      -0.1756
 51 	  57	      -0.0752	       0.0877	  57	   51	       0.0753	      -0.1745
 52 	  53	       0.0001	       0.0000	  53	   52	      -0.0001	      -0.0000
 52 	  54	       0.4622	      -0.1504	  54	   52	      -0.4613	       0.0987
 53 	  55	       0.4611	      -0.1506	  55	   53	      -0.4602	       0.0988
 54 	  55	       0.0000	      -0.0000	  55	   54	      -0.0000	       0.0000
 54 	  61	       0.2284	      -0.2139	  61	   54	      -0.2283	       0.1974
 55 	  61	       0.2285	      -0.2128	  61	   55	      -0.2284	       0.1964
 56 	  57	       0.0001	      -0.0001	  57	   56	      -0.0001	       0.0001
 56 	  58	      -0.3038	       0.0838	  58	   56	       0.3040	      -0.1046
 57 	  58	      -0.3057	       0.0885	  58	   57	       0.3058	      -0.1093
 58 	  59	       0.0083	      -0.0073	  59	   58	      -0.0082	       0.0075
 58 	  72	       0.1303	       0.0127	  72	   58	      -0.1299	      -0.0094
 58 	  87	      -0.1047	       0.1438	  87	   58	       0.1069	      -0.1337
 58 	  98	      -1.4821	       1.1727	  98	   58	       1.5203	      -0.6580
 58 	 100	      -0.1504	       0.1097	 100	   58	       0.1538	      -0.0738
 58 	 103	      -0.0345	       0.0304	 103	   58	       0.0359	      -0.0209
 59 	  60	      -0.0165	      -0.0037	  60	   59	       0.0164	       0.0051
 59 	  72	       0.0038	       0.0047	  72	   59	      -0.0038	      -0.0046
 59 	  79	      -0.0929	       0.2768	  79	   59	       0.0936	      -0.2588
 59 	  80	      -0.0199	       0.0251	  80	   59	       0.0202	      -0.0231
 59 	  89	      -0.0325	       0.0116	  89	   59	       0.0328	      -0.0030
 59 	  92	       0.0587	       0.3170	  92	   59	      -0.0593	      -0.2697
 59 	  94	      -0.0305	       0.0245	  94	   59	       0.0313	      -0.0172
 59 	  98	      -0.4336	       0.4304	  98	   59	       0.4653	      -0.2554
 59 	 100	      -1.0539	       0.7768	 100	   59	       1.0791	      -0.4995
 59 	 103	      -0.6559	       0.5388	 103	   59	       0.6771	      -0.3457
 59 	 107	       0.0665	       0.1192	 107	   59	      -0.0660	      -0.1060
 60 	 135	      -0.0769	       0.0129	 135	   60	       0.0683	       0.0331
 60 	  79	       0.0562	       0.0908	  79	   60	      -0.0565	      -0.0811
 60 	  80	       0.0155	       0.0291	  80	   60	      -0.0154	      -0.0269
 60 	  90	       0.0104	       0.1407	  90	   60	      -0.0107	      -0.1174
 60 	  92	       0.0283	       0.0592	  92	   60	      -0.0292	      -0.0468
 60 	  94	      -1.4388	       1.8108	  94	   60	       1.4437	      -1.4902
 60 	  95	      -0.4761	       0.3082	  95	   60	       0.4548	      -0.0613
 60 	 138	      -0.2237	       0.0053	 138	   60	       0.2096	       0.0641
 61 	  63	       0.5773	      -0.0812	  63	   61	      -0.5751	      -0.0609
 61 	  63	       0.5773	      -0.0812	  63	   61	      -0.5751	      -0.0609
 61 	  64	      -1.6910	       0.7089	  64	   61	       1.6976	      -0.6910
 61 	  65	      -1.6897	       0.7096	  65	   61	       1.6963	      -0.6918
 62 	  86	      -0.1777	       0.0116	  86	   62	       0.1778	      -0.0102
 62 	  86	      -0.1067	       0.0010	  86	   62	       0.1067	      -0.0001
 63 	  64	      -0.3489	       0.0834	  64	   63	       0.3504	      -0.0540
 63 	  65	      -0.3502	       0.0839	  65	   63	       0.3518	      -0.0544
 63 	  66	      -3.5818	       0.6853	  66	   63	       3.6421	       0.2845
 63 	  67	      -0.6064	       0.1991	  67	   63	       0.6170	      -0.1072
 63 	  69	      -0.4720	       0.1473	  69	   63	       0.4741	      -0.1162
 63 	 102	      -1.0910	       0.8104	 102	   63	       1.1068	      -0.5735
 63 	 102	      -1.0955	       0.8143	 102	   63	       1.1115	      -0.5764
 63 	 102	      -1.0769	       0.7996	 102	   63	       1.0925	      -0.5657
 63 	 102	      -1.1195	       0.8325	 102	   63	       1.1359	      -0.5893
 63 	 116	       0.0078	       0.0117	 116	   63	      -0.0079	      -0.0106
 63 	 117	       0.5641	       1.5807	 117	   63	      -0.5572	      -1.4529
 63 	 118	       0.2827	       0.4871	 118	   63	      -0.2859	      -0.4248
 63 	 124	      -0.1192	       0.0657	 124	   63	       0.1173	      -0.0353
 64 	  65	       0.0002	       0.0002	  65	   64	      -0.0002	      -0.0002
 64 	  66	      -3.2736	       0.2529	  66	   64	       3.3084	       0.3588
 64 	  67	      -0.3474	       0.0909	  67	   64	       0.3499	      -0.0683
 64 	  69	       0.1984	       0.0010	  69	   64	      -0.1982	       0.0029
 64 	  97	      -0.0119	       0.0172	  97	   64	       0.0118	      -0.0142
 64 	 124	      -0.0967	       0.0694	 124	   64	       0.0955	      -0.0513
 65 	  66	      -3.2836	       0.2539	  66	   65	       3.3187	       0.3597
 65 	  67	      -0.3490	       0.0914	  67	   65	       0.3515	      -0.0686
 65 	  69	       0.1989	       0.0007	  69	   65	      -0.1987	       0.0032
 65 	  97	      -0.0120	       0.0173	  97	   65	       0.0118	      -0.0142
 65 	 124	      -0.0971	       0.0697	 124	   65	       0.0959	      -0.0515
 66 	  67	       2.2120	       0.2449	  67	   66	      -2.1796	       0.0250
 66 	  68	       0.1248	       0.1659	  68	   66	      -0.2108	      -0.0799
 66 	  69	       6.5638	       0.6584	  69	   66	      -6.4654	       0.6803
 66 	  97	       0.0374	       0.0650	  97	   66	      -0.0379	      -0.0530
 66 	 111	      -4.8841	       5.0426	 111	   66	       4.8841	      -3.9921
 66 	 111	      -4.7379	       5.1062	 111	   66	       4.7603	      -4.0641
 66 	 111	      -4.7230	       4.8764	 111	   66	       4.7230	      -3.8605
 66 	 111	      -4.7730	       5.1457	 111	   66	       4.7956	      -4.0957
 66 	 124	       0.1115	       0.2663	 124	   66	      -0.1134	      -0.2333
 67 	  68	       0.0820	       0.0752	  68	   67	      -0.1178	      -0.0365
 67 	  69	       1.7749	      -0.2558	  69	   67	      -1.7584	       0.4046
 67 	  97	      -0.1444	       1.1600	  97	   67	       0.1517	      -1.0259
 67 	 119	       0.0934	       0.0592	 119	   67	      -0.0957	       0.0375
 67 	 120	       0.4530	       0.2279	 120	   67	      -0.4538	       0.1584
 67 	 121	       0.2451	       0.0705	 121	   67	      -0.2447	      -0.0064
 67 	 122	      -0.1242	       0.2211	 122	   67	       0.1240	      -0.1969
 67 	 124	     -10.2960	      15.9307	 124	   67	      10.3868	     -13.9623
 67 	 125	       1.9958	       0.7775	 125	   67	      -1.9719	       0.1952
 67 	 132	       0.0063	       0.0126	 132	   67	      -0.0063	      -0.0118
 68 	  69	      -0.4937	      -0.1747	  69	   68	       0.3638	       0.3058
 69 	  70	       0.2233	       0.3203	  70	   69	      -0.2223	      -0.2780
 69 	  71	       0.2480	       0.2457	  71	   69	      -0.2473	      -0.2141
 69 	  72	       1.5567	      -0.6267	  72	   69	      -1.5537	       0.6501
 69 	  73	       0.1764	      -0.0334	  73	   69	      -0.1762	       0.0354
 69 	  74	       0.2856	      -0.0569	  74	   69	      -0.2846	       0.0621
 69 	  97	      -0.0746	       0.0904	  97	   69	       0.0739	      -0.0723
 69 	 101	      -0.3664	       0.3335	 101	   69	       0.3700	      -0.2888
 69 	 112	      -0.3472	       0.3389	 112	   69	       0.3507	      -0.2958
 69 	 124	      -0.4261	       0.2681	 124	   69	       0.4205	      -0.1841
 70 	  71	       0.0025	      -0.0099	  71	   70	      -0.0025	       0.0102
 70 	  72	      -0.4265	      -0.8392	  72	   70	       0.4211	       0.9469
 70 	  73	      -0.0602	      -0.1081	  73	   70	       0.0595	       0.1220
 70 	  74	      -0.0570	      -0.1051	  74	   70	       0.0571	       0.1181
 70 	 101	      -0.1350	      -0.0439	 101	   70	       0.1324	       0.0649
 70 	 112	      -0.1311	      -0.0418	 112	   70	       0.1287	       0.0617
 71 	  72	      -0.4961	      -0.6732	  72	   71	       0.4921	       0.7486
 71 	  73	      -0.0696	      -0.0855	  73	   71	       0.0692	       0.0954
 71 	  74	      -0.0646	      -0.0821	  74	   71	       0.0646	       0.0910
 71 	 101	      -0.1212	      -0.0165	 101	   71	       0.1189	       0.0339
 71 	 112	      -0.1178	      -0.0148	 112	   71	       0.1157	       0.0313
 72 	  73	      -0.0959	       0.1334	  73	   72	       0.0959	      -0.1328
 72 	  74	       0.2171	       0.1172	  74	   72	      -0.2170	      -0.1158
 72 	  98	      -1.1978	       0.8431	  98	   72	       1.2222	      -0.4151
 72 	 100	      -0.1264	       0.0784	 100	   72	       0.1288	      -0.0467
 72 	 101	      -1.2665	       0.9068	 101	   72	       1.2669	      -0.7462
 72 	 103	      -0.0289	       0.0220	 103	   72	       0.0300	      -0.0138
 72 	 112	      -1.2148	       0.9257	 112	   72	       1.2152	      -0.7705
 73 	  74	       0.2144	       0.0117	  74	   73	      -0.2144	      -0.0102
 73 	  75	       0.3369	      -0.0934	  75	   73	      -0.3354	       0.1196
 73 	  81	       1.0310	       0.0577	  81	   73	      -1.0418	       0.2139
 73 	  82	       0.0703	       0.0712	  82	   73	      -0.0703	      -0.0544
 73 	  91	      -0.0506	       0.1900	  91	   73	       0.0515	      -0.1716
 73 	  96	      -0.0714	       0.2279	  96	   73	       0.0726	      -0.2052
 73 	 101	      -0.1634	       0.1151	 101	   73	       0.1636	      -0.0952
 73 	 105	      -4.8126	       3.5138	 105	   73	       4.8332	      -2.5558
 73 	 105	      -4.8126	       3.5138	 105	   73	       4.8332	      -2.5558
 73 	 105	      -5.3067	       3.8648	 105	   73	       5.3282	      -2.8093
 73 	 108	       0.0938	       0.1627	 108	   73	      -0.0944	      -0.1456
 73 	 109	       0.0442	       0.0687	 109	   73	      -0.0439	      -0.0520
 73 	 112	      -0.1564	       0.1175	 112	   73	       0.1566	      -0.0983
 73 	 121	       0.1026	       0.0418	 121	   73	      -0.1029	      -0.0238
 74 	  75	       0.2383	      -0.0767	  75	   74	      -0.2371	       0.0938
 74 	  81	       0.6624	       0.0533	  81	   74	      -0.6746	       0.1165
 74 	  82	       0.0672	       0.0721	  82	   74	      -0.0673	      -0.0561
 74 	  91	      -0.0477	       0.1458	  91	   74	       0.0485	      -0.1311
 74 	  96	      -0.0054	       0.0143	  96	   74	       0.0055	      -0.0128
 74 	 101	      -0.1706	       0.1246	 101	   74	       0.1719	      -0.1023
 74 	 106	      -4.6437	       3.8177	 106	   74	       4.7337	      -2.8122
 74 	 106	      -5.0378	       3.5512	 106	   74	       5.0536	      -2.5163
 74 	 108	       0.1013	       0.2077	 108	   74	      -0.1021	      -0.1875
 74 	 109	       0.0369	       0.0587	 109	   74	      -0.0365	      -0.0449
 74 	 112	      -0.1634	       0.1265	 112	   74	       0.1646	      -0.1051
 74 	 121	       0.1254	       0.0537	 121	   74	      -0.1259	      -0.0325
 75 	  82	       0.0683	       0.1390	  82	   75	      -0.0669	      -0.1174
 75 	  91	      -0.0396	       0.0410	  91	   75	       0.0390	      -0.0328
 75 	  96	      -0.0288	       0.0271	  96	   75	       0.0282	      -0.0213
 75 	 108	      -0.1734	       1.1161	 108	   75	       0.1777	      -1.0090
 75 	 109	       0.0515	       0.1542	 109	   75	      -0.0492	      -0.1236
 75 	 121	       0.3358	       0.2485	 121	   75	      -0.3333	      -0.2042
 76 	  77	      -1.2771	       2.1296	  77	   76	       1.2783	      -2.0349
 76 	  89	       1.4857	      -2.1325	  89	   76	      -1.4785	       2.2758
 79 	  80	      -0.2809	      -0.0528	  80	   79	       0.2842	       0.0601
 79 	  90	      -0.0143	       0.0440	  90	   79	       0.0144	      -0.0392
 79 	  92	       0.1895	       0.3369	  92	   79	      -0.1892	      -0.2960
 79 	  94	      -0.1396	       0.0572	  94	   79	       0.1422	      -0.0342
 79 	  95	      -0.0697	       0.0431	  95	   79	       0.0715	      -0.0042
 79 	 107	       0.0550	       0.0219	 107	   79	      -0.0548	      -0.0174
 80 	  90	      -0.0009	       0.0217	  90	   80	       0.0011	      -0.0193
 80 	  92	       0.0599	       0.0776	  92	   80	      -0.0589	      -0.0649
 80 	  94	      -0.0480	       0.0316	  94	   80	       0.0494	      -0.0240
 82 	  91	      -0.0647	      -0.0110	  91	   82	       0.0636	       0.0219
 82 	 108	      -0.1038	      -0.0584	 108	   82	       0.1026	       0.0692
 82 	 109	      -0.0183	       0.2216	 109	   82	       0.0179	      -0.2079
 82 	 121	       0.0149	      -0.0299	 121	   82	      -0.0152	       0.0325
 83 	  89	      -0.4486	       0.1982	  89	   83	       0.4602	      -0.1214
 89 	 103	       0.0050	       0.0184	 103	   89	      -0.0053	      -0.0171
 90 	  92	       0.0103	      -0.0000	  92	   90	      -0.0103	       0.0010
 90 	  94	      -0.1073	      -0.0491	  94	   90	       0.1083	       0.0657
 91 	  96	      -0.0012	      -0.0000	  96	   91	       0.0012	       0.0000
 91 	 108	       0.1196	       0.0034	 108	   91	      -0.1212	       0.0066
 91 	 109	       0.0327	       0.0248	 109	   91	      -0.0331	      -0.0176
 91 	 121	       0.0934	       0.0001	 121	   91	      -0.0960	       0.0184
 92 	  94	      -0.0543	      -0.0064	  94	   92	       0.0552	       0.0187
 92 	 107	       0.0045	      -0.0205	 107	   92	      -0.0045	       0.0220
 94 	  95	      -0.3071	       0.1736	  95	   94	       0.3135	      -0.0545
 94 	 138	      -0.1336	      -0.0450	 138	   94	       0.1315	       0.0801
 95 	 138	       0.2196	      -0.2424	 138	   95	      -0.2289	       0.3231
 96 	 108	       0.0146	       0.0003	 108	   96	      -0.0147	       0.0010
 97 	 124	      -0.0175	      -0.0193	 124	   97	       0.0172	       0.0207
 98 	 100	       0.2372	      -0.1168	 100	   98	      -0.2377	       0.1411
 98 	 103	       0.0425	      -0.0205	 103	   98	      -0.0424	       0.0239
100 	 103	      -0.0306	       0.0277	 103	  100	       0.0306	      -0.0269
101 	 112	       0.0089	       0.0061	 112	  101	      -0.0089	      -0.0061
102 	 117	      10.1154	       0.5053	 117	  102	     -10.1451	       1.3717
102 	 118	       0.7061	       0.1692	 118	  102	      -0.7197	      -0.0056
108 	 109	       0.0508	       0.0844	 109	  108	      -0.0516	      -0.0724
108 	 121	       2.6546	      -0.5551	 121	  108	      -2.6611	       0.8634
109 	 121	       0.0090	      -0.0306	 121	  109	      -0.0093	       0.0353
115 	 116	       0.8265	       0.2025	 116	  115	      -0.8259	      -0.1834
115 	 117	      -0.0281	       0.0886	 117	  115	       0.0280	      -0.0869
115 	 118	       0.5536	       0.6516	 118	  115	      -0.5565	      -0.6067
115 	 143	      -0.0794	       0.0232	 143	  115	       0.0788	      -0.0201
116 	 117	      -0.9651	       0.5482	 117	  116	       0.9672	      -0.5156
116 	 118	       0.3711	       0.7938	 118	  116	      -0.3718	      -0.7628
116 	 143	      -0.0470	       0.0030	 143	  116	       0.0466	      -0.0004
117 	 118	       5.7412	       2.0955	 118	  117	      -5.7130	      -1.8103
117 	 143	      -0.0441	      -0.0062	 143	  117	       0.0440	       0.0075
118 	 131	       0.0194	      -0.0013	 131	  118	      -0.0197	       0.0036
118 	 132	      -0.0233	      -0.0038	 132	  118	       0.0230	       0.0082
118 	 143	      -3.2005	      -0.9536	 143	  118	       3.1885	       1.2061
119 	 120	      -6.0322	      -0.1963	 120	  119	       6.0677	       1.0340
119 	 121	      -2.3328	       0.6190	 121	  119	       2.2703	       1.0278
119 	 122	      -0.1501	       0.0642	 122	  119	       0.1345	       0.0888
119 	 124	      -0.2612	       0.1220	 124	  119	       0.2400	       0.1529
119 	 125	      -1.7881	       0.3849	 125	  119	       1.7614	       0.4611
119 	 126	      14.7632	      -1.5967	 126	  119	     -14.4344	       5.4431
119 	 127	      -0.2895	       0.0383	 127	  119	       0.2797	       0.0774
119 	 128	      -6.8094	       0.2303	 128	  119	       6.5652	       2.1039
119 	 129	       3.6930	       0.7597	 129	  119	      -3.6459	       0.1296
119 	 130	      -8.1860	      -3.2703	 130	  119	       8.0194	       4.5046
119 	 131	     -25.4568	       2.0873	 131	  119	      22.6594	      13.2983
119 	 132	      -0.3634	       0.0951	 132	  119	       0.3065	       0.2352
119 	 144	      -0.2160	       0.0533	 144	  119	       0.1749	       0.1317
120 	 121	      -7.2068	       1.9449	 121	  120	       7.2538	       2.1220
120 	 122	      -0.8592	       0.3596	 122	  120	       0.8096	       0.3967
120 	 123	      -0.6721	       0.0786	 123	  120	       0.6521	       0.1363
120 	 124	      -1.7103	       0.7520	 124	  120	       1.6256	       0.7924
120 	 125	      -6.1197	       1.4619	 125	  120	       6.1123	       0.5966
120 	 127	      -1.4968	       0.3654	 127	  120	       1.5012	       0.0390
120 	 128	      -2.9505	       0.3038	 128	  120	       2.9266	       0.3087
120 	 129	       0.7412	       0.2889	 129	  120	      -0.7548	       0.0023
120 	 130	       0.0081	      -0.0218	 130	  120	      -0.0082	       0.0224
120 	 131	      -1.1045	       0.0769	 131	  120	       1.0256	       0.4418
120 	 132	      -1.6796	       0.5582	 132	  120	       1.6048	       0.7820
121 	 122	      -0.6520	       0.1842	 122	  121	       0.6474	       0.0184
121 	 123	       0.1187	       0.0398	 123	  121	      -0.1211	      -0.0119
121 	 124	      -0.9452	       0.2731	 124	  121	       0.9399	       0.0361
121 	 125	      18.2632	       5.1570	 125	  121	     -18.2632	      -1.0755
121 	 127	       0.3499	       0.1077	 127	  121	      -0.3524	      -0.0056
121 	 128	       1.1236	       0.3708	 128	  121	      -1.1592	       0.0252
121 	 129	       0.1782	       0.1308	 129	  121	      -0.1985	       0.0592
121 	 131	       0.0503	       0.0119	 131	  121	      -0.0508	      -0.0083
121 	 132	      -0.1775	       0.0063	 132	  121	       0.1737	       0.0335
122 	 123	       0.1012	       0.0373	 123	  122	      -0.1080	       0.0193
122 	 124	      -0.2857	      -0.0024	 124	  122	       0.2857	       0.0069
122 	 125	       3.1259	       0.9259	 125	  122	      -3.1992	       0.7566
122 	 131	       0.1926	       0.0401	 131	  122	      -0.2020	       0.0348
122 	 132	       0.3238	      -0.1275	 132	  122	      -0.3261	       0.1587
122 	 133	       0.1766	      -0.0632	 133	  122	      -0.1801	       0.0977
122 	 143	       0.3980	      -0.0001	 143	  122	      -0.4029	       0.0776
123 	 124	      -0.2741	       0.0527	 124	  123	       0.2573	       0.0954
123 	 125	      -0.0173	       0.0124	 125	  123	       0.0172	      -0.0121
123 	 131	      -0.1254	      -0.0281	 131	  123	       0.1226	       0.0481
123 	 132	      -0.3869	       0.0237	 132	  123	       0.3672	       0.1511
124 	 125	       5.3995	       1.4984	 125	  124	      -5.4528	       1.4814
124 	 128	       0.0716	       0.0331	 128	  124	      -0.0787	       0.0182
124 	 131	       0.4720	       0.1063	 131	  124	      -0.4969	       0.0853
124 	 132	       0.6149	      -0.1965	 132	  124	      -0.6188	       0.2636
124 	 133	       0.1579	      -0.0603	 133	  124	      -0.1589	       0.0940
124 	 143	       0.2768	      -0.0094	 143	  124	      -0.2774	       0.0677
125 	 127	       0.0677	       0.0092	 127	  125	      -0.0680	      -0.0046
125 	 128	       0.2128	       0.0073	 128	  125	      -0.2156	       0.0194
125 	 129	       0.1528	       0.0830	 129	  125	      -0.1653	       0.0351
125 	 130	       0.0391	       0.0103	 130	  125	      -0.0422	       0.0033
125 	 131	      -0.2053	      -0.0701	 131	  125	       0.1996	       0.1022
125 	 132	      -0.8887	       0.0387	 132	  125	       0.8470	       0.3571
127 	 128	       0.4826	      -0.1225	 128	  127	      -0.4833	       0.1528
127 	 129	       0.5242	       0.2192	 129	  127	      -0.5367	       0.1336
128 	 129	      26.1404	      11.4825	 129	  128	     -26.9163	       4.5784
128 	 130	       0.0681	       0.0222	 130	  128	      -0.0735	      -0.0076
128 	 131	      -0.0622	      -0.0186	 131	  128	       0.0560	       0.0351
130 	 131	     -33.6453	       3.2395	 131	  130	      30.8843	      12.5086
130 	 132	      -0.2666	       0.0527	 132	  130	       0.2236	       0.1476
130 	 144	      -0.2487	       0.0531	 144	  130	       0.2051	       0.1244
131 	 132	      -7.8003	       0.5572	 132	  131	       7.6201	       1.7578
131 	 133	      -0.0428	      -0.0128	 133	  131	       0.0409	       0.0229
131 	 143	      -0.4902	       0.0029	 143	  131	       0.4772	       0.0869
131 	 144	     -19.0859	       2.9042	 144	  131	      18.3307	       2.2791
132 	 133	       0.1103	      -0.0473	 133	  132	      -0.1116	       0.0582
132 	 143	       1.2325	       0.2512	 143	  132	      -1.2397	      -0.1106
132 	 144	       0.0190	       0.0501	 144	  132	      -0.0193	      -0.0475
133 	 143	       0.0117	       0.0273	 143	  133	      -0.0120	      -0.0253
134 	 131	       0.2193	       0.1321	 131	  134	      -0.2436	      -0.0771
134 	 136	      -0.4559	      -0.0512	 136	  134	       0.4424	       0.1753
134 	 139	      -0.0339	       0.0180	 139	  134	       0.0339	      -0.0178
134 	 141	      -0.2524	      -0.1248	 141	  134	       0.2507	       0.1334
134 	 142	       0.1977	      -0.9485	 142	  134	      -0.2203	       1.0490
134 	 144	      -1.1871	       0.7507	 144	  134	       1.1608	      -0.6719
134 	 145	     -13.7809	      -0.6250	 145	  134	      13.1872	       4.3963
135 	  95	       0.0448	       0.0684	  95	  135	      -0.0467	      -0.0494
135 	 136	      25.9810	      12.1584	 136	  135	     -28.0625	      -0.2065
135 	 138	       2.0709	       0.3699	 138	  135	      -2.1013	       0.2545
135 	 141	       0.9173	       0.6110	 141	  135	      -1.0452	       0.0822
136 	 115	       4.2868	       1.8335	 115	  136	      -4.5092	      -0.2488
136 	 116	       0.0816	       0.0511	 116	  136	      -0.0911	      -0.0174
136 	 117	       0.0334	       0.0244	 117	  136	      -0.0378	      -0.0112
136 	 118	       0.1952	       0.1680	 118	  136	      -0.2277	      -0.0763
136 	 138	      -0.2584	      -0.1158	 138	  136	       0.2476	       0.1533
136 	 139	       8.9680	       4.6955	 139	  136	      -9.4835	      -2.1356
136 	 140	       0.0529	       0.0342	 140	  136	      -0.0610	      -0.0025
136 	 141	      14.3607	       5.7928	 141	  136	     -14.8922	      -2.2151
136 	 142	       1.8280	       0.2966	 142	  136	      -1.9646	       0.2031
136 	 143	       0.0664	       0.0661	 143	  136	      -0.0796	      -0.0402
136 	 145	      -0.2875	       0.5980	 145	  136	       0.2857	      -0.5778
137 	 139	     -11.3084	       4.3651	 139	  137	       8.9332	       7.7831
137 	 140	      -0.1056	       0.0160	 140	  137	       0.0831	       0.0648
137 	 145	      -3.0147	       1.5665	 145	  137	       2.1461	       2.5842
139 	 140	      11.4267	       3.8312	 140	  139	     -12.1518	      -0.6217
139 	 141	      -0.5282	      -0.3817	 141	  139	       0.5249	       0.3998
139 	 142	       0.0245	      -0.0884	 142	  139	      -0.0269	       0.0982
139 	 145	     -36.6502	      -0.6600	 145	  139	      35.5322	      10.5985
140 	 145	      -1.2041	       0.0552	 145	  140	       1.0607	       0.5773
141 	 115	       9.4580	       1.3673	 115	  141	      -9.5157	      -0.2884
141 	 116	       0.1846	       0.0664	 116	  141	      -0.1901	      -0.0406
141 	 117	       0.0736	       0.0418	 117	  141	      -0.0760	      -0.0329
141 	 118	       0.9232	       0.6648	 118	  141	      -0.9715	      -0.4968
141 	 131	       0.3106	       0.1505	 131	  141	      -0.3357	      -0.0632
141 	 132	      -0.0052	       0.0005	 132	  141	       0.0052	      -0.0004
141 	 142	       4.8603	      -9.3501	 142	  141	      -5.0405	      10.4017
141 	 143	       0.3674	       0.2945	 143	  141	      -0.3814	      -0.2590
141 	 144	      -0.1055	       0.2091	 144	  141	       0.1018	      -0.1970
141 	 145	      -7.5763	       0.1627	 145	  141	       7.3795	       1.6914
142 	 115	       0.5667	       0.8716	 115	  142	      -0.5801	      -0.7449
142 	 116	       0.0330	       0.0605	 116	  142	      -0.0355	      -0.0511
142 	 117	       0.0245	       0.0727	 117	  142	      -0.0270	      -0.0628
142 	 118	       1.0462	       1.8870	 118	  142	      -1.1108	      -1.5251
142 	 119	       0.4186	       0.3634	 119	  142	      -0.4818	       0.0653
142 	 120	       0.1007	       0.0670	 120	  142	      -0.1073	       0.0136
142 	 122	      -0.0827	       0.0776	 122	  142	       0.0802	      -0.0567
142 	 124	      -0.0900	       0.0830	 124	  142	       0.0880	      -0.0590
142 	 125	       0.0459	       0.0352	 125	  142	      -0.0486	      -0.0136
142 	 130	       0.3741	       0.2951	 130	  142	      -0.4356	       0.0218
142 	 131	      17.0228	      11.8672	 131	  142	     -17.4424	      -6.7994
142 	 132	      -1.0112	       1.5702	 132	  142	       0.9969	      -1.3584
142 	 133	       0.0006	       0.0080	 133	  142	      -0.0007	      -0.0076
142 	 143	       1.6171	       8.0708	 143	  142	      -1.8101	      -7.1210
142 	 144	      -2.5701	       7.7801	 144	  142	       2.4694	      -6.6276
142 	 145	      -0.7938	       0.2420	 145	  142	       0.7557	      -0.0159
143 	 144	      -0.0409	       0.0083	 144	  143	       0.0401	      -0.0045
144 	 145	      -0.1797	      -0.0783	 145	  144	       0.1649	       0.1249
  1 	   3	       0.0191	      -0.3687	   3	    1	      -0.0201	       0.3511
  1 	   4	       0.0191	      -0.3687	   4	    1	      -0.0201	       0.3511
  1 	   5	       0.0195	      -0.3713	   5	    1	      -0.0204	       0.3537
  1 	  33	      -0.9410	       3.2748	  33	    1	       0.9419	      -3.2227
  1 	  93	      -2.9958	      -1.3479	  93	    1	       2.9980	       1.5031
  1 	  93	      -2.9958	      -1.3479	  93	    1	       2.9980	       1.5031
  2 	 113	       0.0000	      -0.0000	 113	    2	       0.0000	       0.0000
  2 	 114	       0.0000	       0.0000	 114	    2	      -0.0000	      -0.0000
  7 	   8	      -0.2884	       0.0591	   8	    7	       0.2876	      -0.0698
  7 	  66	       4.1876	      -0.5809	  66	    7	      -4.1854	       0.7221
  7 	 104	      -9.3774	      -0.6072	 104	    7	       9.4109	       2.3762
  7 	 104	     -10.2416	      -0.6157	 104	    7	      10.2871	       2.5471
  9 	  11	      -0.1941	      -0.1405	  11	    9	       0.1932	       0.1268
  9 	  69	       2.7785	       2.2553	  69	    9	      -2.7745	      -2.0681
 10 	  32	      -0.1981	      -0.1375	  32	   10	       0.1969	       0.1238
 10 	  69	       2.7925	       2.2654	  69	   10	      -2.7885	      -2.0774
 12 	  13	      -0.1950	      -0.1190	  13	   12	       0.1941	       0.1064
 12 	  13	      -0.1916	      -0.1161	  13	   12	       0.1906	       0.1038
 12 	  13	      -0.1916	      -0.1161	  13	   12	       0.1906	       0.1038
 12 	  72	       2.8242	       1.9732	  72	   12	      -2.8214	      -1.7986
 12 	  72	       2.8091	       1.9631	  72	   12	      -2.8064	      -1.7893
 12 	  72	       2.8091	       1.9631	  72	   12	      -2.8064	      -1.7893
 14 	  15	      -0.0521	      -0.1790	  15	   14	       0.0510	       0.1682
 14 	  16	      -0.1129	      -0.4326	  16	   14	       0.1114	       0.4066
 14 	  58	       1.3868	       3.0873	  58	   14	      -1.3850	      -3.0008
 17 	  18	      -0.0374	      -0.0570	  18	   17	       0.0363	       0.0523
 17 	  19	      -0.0389	      -0.1221	  19	   17	       0.0389	       0.1118
 17 	  20	      -0.0394	      -0.0628	  20	   17	       0.0394	       0.0592
 17 	  21	      -0.2646	      -0.2819	  21	   17	       0.2635	       0.2636
 17 	  59	       4.4993	       5.3889	  59	   17	      -4.4956	      -5.1242
 22 	  23	      -0.0909	      -0.0252	  23	   22	       0.0909	       0.0224
 22 	  30	      -0.0968	      -0.0316	  30	   22	       0.0968	       0.0289
 22 	  78	       0.9868	       0.3218	  78	   22	      -0.9868	      -0.2972
 22 	  83	       0.9093	       0.2513	  83	   22	      -0.9093	      -0.2260
 24 	  76	      -6.7335	       2.3403	  76	   24	       6.7429	      -1.9251
 24 	  77	       1.2729	      -0.9568	  77	   24	      -1.2783	       0.8148
 25 	  26	      -0.2506	      -0.3460	  26	   25	       0.2497	       0.3264
 25 	  31	      -0.2680	      -0.2787	  31	   25	       0.2670	       0.2595
 25 	  73	       1.5392	       2.2851	  73	   25	      -1.5374	      -2.1833
 25 	  74	       1.9756	       2.2054	  74	   25	      -1.9729	      -2.0831
 27 	  28	      -0.0783	      -0.0946	  28	   27	       0.0769	       0.0860
 27 	  29	      -0.2127	      -0.2980	  29	   27	       0.2111	       0.2712
 27 	  75	       4.7704	       3.1014	  75	   27	      -4.7665	      -2.8544
 33 	 110	      -2.9857	      -2.0072	 110	   33	       2.9890	       2.2252
 33 	 110	      -3.0037	      -2.0217	 110	   33	       3.0070	       2.2411
 36 	  99	      -1.8925	       0.3000	  99	   36	       1.8954	      -0.1357
 37 	  87	       0.4250	       0.1123	  87	   37	      -0.4248	      -0.1048
 61 	  62	      -0.1930	       0.0031	  62	   61	       0.1918	      -0.0117
 61 	  62	      -0.0930	      -0.0033	  62	   61	       0.0927	      -0.0009
 61 	  86	       1.3010	       0.1253	  86	   61	      -1.2989	      -0.0767
 61 	  86	       1.1243	       0.1213	  86	   61	      -1.1230	      -0.0793
 61 	  86	       1.1243	       0.1213	  86	   61	      -1.1230	      -0.0793
----------------------------------------------------------------------------------
Total real power losses in the system = 748.026022
Total reactive power losses in the system = 192.552024
```



```data title="shunt.dat"
   3   0.0000 -1.2600    
   4   0.0000 -1.2600   
   5   0.0000 -1.2600    
  15   0.0000 -1.2600    
  16   0.0000 -1.2600    
  17   0.0000 -2.5000    
  18   0.0000 -1.2600    
  19   0.0000 -0.5800    
  28   0.0000 -1.2600    
  29   0.0000 -1.2600    
  41   0.5400 -0.1000    
  42   0.5480 -0.1000    
  47   0.1520  0.0537    
  48   0.1380  0.0520    
  52   0.2300 -0.0970    
  53   0.2320 -0.0970    
  54   0.1880 -0.0930    
  55   0.1870 -0.0920    
  56   0.1860 -0.0750    
  57   0.1880 -0.0700   
  58   1.2100  1.4000   
  59   5.0800  0.0863    
  60   2.0100  1.6500    
  61   0.0000  1.4900   
  63   8.2300  5.2500    
  64   1.2300 -0.2600    
  65   1.2400 -0.2600    
  66   2.1600 18.9700    
  67  18.2100 12.8000    
  68   0.5630 -0.2500   
  69   9.7100 -1.4400    
  70   1.0300 -0.8600    
  71   1.0600 -1.2200   
  72  10.1900  0.2460    
  73  12.2200  8.4100    
  74   8.5700  5.7400    
  75   3.8700  1.9900    
  76  -6.6700 -1.8500    
  77   0.0000 -1.2500    
  79  2.3800  0.2710   
  80  -0.0008  0.2990    
  81   0.7000 -0.4700    
  82   1.1100 -0.2800   
  83   1.0500 -0.0042    
  86   2.9200 -0.2200    
  87   0.3770 -0.4700    
  89   6.7800  0.7110    
  90   0.2950 -0.1900    
  91   0.1930 -0.2700    
  92   0.4340 -0.4300    
  94   1.4000 -1.4100    
  95   0.2720 -0.5800    
  96   0.4780 -0.4500    
  97   1.2900 -1.7700    
  98   0.8210 -1.2100    
 100  0.5910  -1.0100   
 101   0.6470 -2.3800    
 102   4.5800 -6.1900    
 103   0.6240 -0.4500    
 105   0.2420 -9.9900    
 106   0.3690 -7.1900    
 107   0.2880 -0.2200    
 108   5.2000 -2.4400    
 109  0.8450 -0.4700   
 111  0.2330 -9.9900   
 112   0.6430 -2.3600    
 115  28.5600  0.0941    
 116  19.3000 -1.9500   
 117  24.1000 -0.5400    
 118  47.8800 -9.7900    
 119  99.9900 -9.9900    
 120  20.2600  8.5500    
 121   3.0600 10.5500    
 122   3.9700 -0.6100    
 123  1.1800   0.0140    
 124  7.6600  -9.9900    
 125  28.3100 -0.7100    
 126  16.0400  3.4700    
 127   3.0500 -0.5100   
 128  52.5200 -4.7200    
 129  38.5500  3.8700    
 130 38.3000  -0.9100    
 131  51.4500 -7.8000    
 132  12.3900 -2.2600   
 133   1.0300 -0.0590    
 134  -1.4100 -0.7800    
 135  -9.9900  8.3500    
 136  -9.9900  3.7500    
 137   4.9900 -2.1900    
 138   4.3300 -0.8100    
 139  16.0900 -9.9900    
 140  -2.8900 -0.7700    
 141  52.1200 -9.9900    
 142  43.2300 22.1000    
 143  5.8600 -9.9900    
 144  -4.3600 -9.9900    
 145  -9.9900  4.5700   
```


```dat title="slip_pss.dat"
104 15 0.05 10 0.09189 0.03063 0.1 -0.1  570  35  1
111 20 0.05 10 0.09189 0.03063 0.1 -0.1  570  35  1
 93 05 0.05 10 0.09189 0.03063 0.1 -0.1  570  35  1
110 05 0.05 10 0.09189 0.03063 0.1 -0.1  570  35  1
```


```c title="slip_pss_settings.c"
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        FORMATION OF A MATRIX with PSS


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
APS1=[];
BPS1=[];
CPS1=[];
DPS1=[];
APS2=[];
BPS2=[];
CPS2=[];
DPS2=[];
APS=[];
BPS=[];
CPS=[];
DPS=[];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(sum(xpss_nt)~=0)
 for k=1:size(xpss_nt,2)
  [NPSS1,DPSS1]=series([Tw_slip_g(xpss_nt(k))*Kslip_pss_g(xpss_nt(k)) 0],[Tw_slip_g(xpss_nt(k)) 1],[T1_slip_g(xpss_nt(k)) 1],[T2_slip_g(xpss_nt(k)) 1]);
  if (Tmd_slip_nt==0)
      [NPSS1,DPSS1]=series(NPSS1,DPSS1,[1],[Td1_slip_g(xpss_nt(k)) 1]);
  end    
  [A1,B1,C1,D1]=tf2ss(NPSS1,DPSS1);
  [APS1,BPS1,CPS1,DPS1]=append(APS1,BPS1,CPS1,DPS1,A1,B1,C1,D1);
  end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(sum(xpss_t)~=0)
 for k=1:size(xpss_t,2)
  [NPSS2,DPSS2]=series([Tw_slip_g(xpss_t(k))*Kslip_pss_g(xpss_t(k)) 0],[Tw_slip_g(xpss_t(k)) 1],[T1_slip_g(xpss_t(k)) 1],[T2_slip_g(xpss_t(k)) 1]);
  [NPSS2,DPSS2]=series(NPSS2,DPSS2,[a0_slip_g(xpss_t(k))],[1 a1_slip_g(xpss_t(k)) a0_slip_g(xpss_t(k))]);
  if (Tmd_slip_t==0)
      [NPSS2,DPSS2]=series(NPSS2,DPSS2,[1],[Td1_slip_g(xpss_t(k)) 1]);
  end 
  [A2,B2,C2,D2]=tf2ss(NPSS2,DPSS2);
  [APS2,BPS2,CPS2,DPS2]=append(APS2,BPS2,CPS2,DPS2,A2,B2,C2,D2);
  end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[APS,BPS,CPS,DPS]=append(APS1,BPS1,CPS1,DPS1,APS2,BPS2,CPS2,DPS2);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
APS=sparse(APS);
BPS=sparse(BPS);
CPS=sparse(CPS);
DPS=sparse(DPS);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SL_x_pss=[xpss_nt xpss_t]';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for k=1:size(slip_pss,1)
    for p=1:ngen
        if (SL_x_pss(k,1)==gen(p,1))
            temp(k)=p;
        end
    end
end

PS1=sparse(zeros(size(slip_pss,1),ngen));
PS1(1:(size(slip_pss,1)),temp)=(eye(size(slip_pss,1)));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PS1 matrix is used to adjust the dimension of DPS 
AG=([AG+EG*PS1'*DPS*PS1*FG EG*PS1'*CPS; BPS*PS1*FG APS]);
BG=([BG; sparse(zeros(size(APS,1),2*ngen))]);
CG=([CG  sparse(zeros(size(APS,1),2*ngen)')]);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear C1 D1;
clear A1 B1;
%clear APS BPS;
%clear CPS DPS;
clear temp;
%end
```


```c title="small_sig.c"




% A MATLAB program for the small signal analysis of a multimachine power
% system. 2.2 Machine model is used with a STATIC,DC1A,AC4A EXCITERS;HYDRO,RHT STEAM TURBINES;SLIP PSS ;



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all;

load busno.dat;

load nt.dat;

load gen.dat;

load ld.dat;

	if (busno(9)~=0)
load shunt.dat;
end

nb=busno(3);
nline=busno(4);
ntrans=busno(5);
ngen=size(gen,1);
nshunt=busno(9);
nload=size(ld,1);
wB=busno(11)*2*pi


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Calculate the initial conditions
initcond
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~
% Formation of YDQ and YBUS matrix
yform
%~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Formation of the PG and  PL matrices
pmat
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Selection of type of exciters: Single time constant static, DC1A, AC4A
exciter_setting
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Selection of type of turbines: Hydro and Reheat steam turbine
primemover_setting
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Formation of Generator matrices:AG, BG, CG and DG (YG)
genmat1
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~

YDQdash=PG*YG*PG';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Formation of YL matrices for Static loads 
statld
YDQdash=YDQdash+ PL*YL*PL';
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% YDQ  matrix is formed in yform.m
YDQdash=YDQdash+YDQ;
YDQdi=YDQdash\eye(size(YDQdash));

% Anp -state matrix without PSS.
%------------------------------------
Anp = sparse(AG+BG*PG'*YDQdi*PG*CG);
CGnp=CG;
%Updation of AG, BG, CG matrices to include slip signal-based PSS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AT = Anp;

if (npss~=0)
   slip_pss_settings
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (npss1~=0)
   delPw_pss_settings
end
if (npss2~=0)
   power_pss_settings
end

if ((npss+npss1+npss2)~=0)
   AT=sparse(AG+BG*PG'*YDQdi*PG*CG);
   Ap=AT;
end
%Unreduced A matrix for the system

% Formation of Reduced A matrix
AT(1,:)=[];
AT(:,1)=[];
aa=sparse(zeros(size(AT)*[1 0]',size(AT)*[0 1]'));
aa(16:16:16*ngen-1,1)=ones(ngen-1,1);
AT=(AT-wB*aa);

% The Order of occurrences of states X for a genarator is decided by its
% row position in gen.dat
```


```c title="static_Exciter.c"
%----------------------Static - type exciters--------------------------

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load exc_static.dat
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[mm order_static]=sort([setxor(1:nb,exc_static(:,1)) exc_static(:,1)]);

kA_static=exc_static(:,2)';
TA_static=exc_static(:,3)';

Efd_min_static=exc_static(:,4)';
Efd_max_static=exc_static(:,5)';

Efd0_static=EFD0(exc_static(:,1));

Vref_static=Efd0_static./kA_static + lfl(exc_static(:,1),2)';
```


```c title="statId.c"

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~               
%            Formation of YL (STATIC LOADS)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mp=0;
mi=1;
mz=2;
np=0;
ni=1;
nz=2;

mk=mp*p1+mi*p2+mz*p3;
nk=np*r1+ni*r2+nz*r3;

% VL0 and theL0 are obtained in initcond.m
VDl0=VL0.*sin(theL0);
VQl0=VL0.*cos(theL0);


BDQ=(QL0.*((nk-2).*VQl0.*VQl0./VL0./VL0+1)./VL0./VL0)-(PL0.*((mk-2).*VQl0.*VDl0./VL0./VL0)./VL0./VL0);

GDD=(PL0.*((mk-2).*VDl0.*VDl0./VL0./VL0+1)./VL0./VL0)-(QL0.*((nk-2).*VQl0.*VDl0./VL0./VL0)./VL0./VL0);

GQQ=(PL0.*((mk-2).*VQl0.*VQl0./VL0./VL0+1)./VL0./VL0)+(QL0.*((nk-2).*VQl0.*VDl0./VL0./VL0)./VL0./VL0);

BQD=(QL0.*((nk-2).*VDl0.*VDl0./VL0./VL0+1)./VL0./VL0)+(PL0.*((mk-2).*VQl0.*VDl0./VL0./VL0)./VL0./VL0);



YLmod=sparse([diag(-BDQ) diag(GDD);diag(GQQ) diag(BQD)]);


YL=YLmod(mmm,mmm);

clear YLmod;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

```c title="trace.c"
%~~~~~~~~~~~~~~~~~~~~~eig_test_pss_final~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
format short g
YES2=input('Enter 1 to display  ALL eigenvalues (EIG), otherwise 0 (using EIGS):  ');
if(YES2)
   [V1 D]=eig(full(AT));
   V1=sparse(V1);
   W=V1\eye(size(V1));
   PF=V1.*conj(W');
   D=diag(D);
   n_eig=1;
   for k=1:size(AT,1)-3
    if (sum(abs(D(k+[0:3])))>0.1)
     n_eig = n_eig+1;
    end
   end
   D=diag(D(1:n_eig))-0.000001;
else
   n_eig=10;
   f_eig=1;
   fprintf('You  are scanning %d eigenvalues around %5.3f Hz .....\n', n_eig,f_eig);
   disp('Please press a key: ')
   pause
   options.tol=1e-12;
   options.disp =0;
   options.maxit = 25;
   [V1 D Flag1]= eigs(AT,n_eig,j*(f_eig*2*pi),options);
   [W1 D Flag2]= eigs((AT)',n_eig,j*(f_eig*2*pi),options);
   VW=conj(V1)'*W1;
   PF=(V1.*W1)/VW;
   EIG=diag(D)-0.000001;
   D=diag(EIG);
   if ((Flag1==1)|(Flag2==1))
      disp('----------------------------NOTE---------------------------------------');
      disp('Determination of eigenvalues with EIGS is NOT accurate. Please use EIG function');
      end
end

SL_lamb=(1:n_eig)';
disp('-------------------------------------------------------------------');
disp(   'SL_number         Eigenvalue       dampingfactor     frequency(Hz)')
disp('--------------------------------------------------------------------');
[SL_lamb diag(D) -real(diag(D))./sqrt(real(diag(D)).*real(diag(D))+imag(diag(D)).*imag(diag(D)))  abs(imag(diag(D)))./(2*pi) ]   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('---------------------------------------------------------------------');
EIG=diag(D);
rpt=1;
while (rpt==1)
SL_num=input('Enter the serial number of the eigenvalue for which you want to obtain the P.factor:  ');                                                                                              
PF_lamb=PF(:,SL_num);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PF_lamb_nr = PF_lamb/max(PF_lamb);
[m_sort m_ind]=sort(abs(PF_lamb_nr));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (npss~=0)
  if (Tmd_slip_nt==0)
     nm_st_pss_nt = 3;  % if you are using input measuring time delay 
  else
      nm_st_pss_nt = 2; 
  end
  
  if (Tmd_slip_t==0)
     nm_st_pss_t = 5;  % if you are using input measuring time delay 
  else
      nm_st_pss_t = 4; 
   end
end %-----End of npss if loop-------------------


 rev_ind=[size(AT,1):-1:1];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_ind=m_ind(rev_ind);
n_ind=n_ind+1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X1=('');
X3=('');
for i=1:ngen
X1=strvcat(X1,'Delta-','Slip-','Field-','DampH-','DampG-','DampK-','Efd-','VR_DC-','xB_DC_AC-','xF_DC-','x1_st_tu-','x2_st_tu-','x3_st_tu-','y1_gv-','PG-','z-');
X2d=gen(i,1)*ones(16,1);
X2=num2str(X2d);
X3=strvcat(X3,X2);
end

Xt=strcat(X1,X3);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(npss~=0)
   XPt_nt =('');
   XPt_t=('');
   XP_nt=('');
   X4_nt=('');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(sum(xpss_nt)~=0)
for i=1:size(xpss_nt,2);
    
  if (Tmd_slip_nt==0)
     XP_nt=strvcat(XP_nt,'xpss_1-','xpss_2-','xpss_3-');  % if you are using input measuring time delay 
  else
      XP_nt=strvcat(XP_nt,'xpss_1-','xpss_2-'); 
  end

   XPd_nt=xpss_nt(i)*ones(nm_st_pss_nt,1);
   XP2_nt=num2str(XPd_nt);
   X4_nt=strvcat(X4_nt,XP2_nt);
end
XPt_nt=strcat(XP_nt,X4_nt);
end
     %-------------------------------------------------
if(sum(xpss_t)~=0)
   XP_t=('');
   X4_t=(''); 
 for i=1:size(xpss_t,2);
     
  if (Tmd_slip_t==0)
     XP_t=strvcat(XP_t,'xpss_1-','xpss_2-','xpss_3-','xpss_4-','xpss_5-');  % if you are using input measuring time delay 
   else
      XP_t=strvcat(XP_t,'xpss_1-','xpss_2-','xpss_3-','xpss_4-');
   end
  
   XPd_t=xpss_t(i)*ones(nm_st_pss_t,1);
   XP2_t=num2str(XPd_t);
   X4_t=strvcat(X4_t,XP2_t);
  end
XPt_t=strcat(XP_t,X4_t);
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 XPt=strvcat(XPt_nt,XPt_t);
 Xt=strvcat(Xt,XPt);
end   %-----End of npss if loop-------
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(npss1~=0)
   XdelPW=('');
   X5=('');
for i=1:size(delPw_pss,1);
   XdelPW=strvcat(XdelPW,'xdelPw_1-','xdelPw_2-','xdelPw_3-','xdelPw_4-','xdelPw_5-','xdelPw_6-','xdelPw_7-','xdelPw_8-','xdelPw_9-','xdelPw_10-','xdelPw_11-','xdelPw_12-','xdelPw_13-','xdelPw_14-','xdelPw_15-','xdelPw_16-');
   XdelPWd=delPw_pss(i,1)*ones(16,1);
   XdelPWd2=num2str(XdelPWd);
   X5=strvcat(X5,XdelPWd2);
end
  XPt=strcat(XdelPW,X5);
  Xt=strvcat(Xt,XPt);
end  %-----End of npss1 if loop-------

if(npss2~=0)
   Xpower=('');
   X6=('');
for i=1:size(power_pss,1);
   Xpower=strvcat(Xpower,'xpower_1-','xpower_2-');
   Xpowerd=power_pss(i,1)*ones(2,1);
   Xpowerd2=num2str(Xpowerd);
   X6=strvcat(X6,Xpowerd2);
end
  XPt=strcat(Xpower,X6);
  Xt=strvcat(Xt,XPt);
end  %-----End of npss2 if loop-------
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Yt=Xt(n_ind,:);
PF_lamb_sort=PF_lamb_nr(n_ind-1);

disp('----------------------------------------------------------------------------------');
disp('      State variable      Mag(Norm PF)    ang(Norm PF)deg.  Mag(PF)    ang(PF)deg.');
disp('----------------------------------------------------------------------------------');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   for i=1:(size(AT,1))
      if (abs(PF_lamb_sort(i))>=0.1)
      fprintf('\t %s %15.4f %13.2f %15.4f %13.2f\n',Yt(i,:),abs(PF_lamb_sort(i)),angle(PF_lamb_sort(i))*180/pi,abs(PF_lamb_sort(i)*full(max(PF_lamb))),angle(PF_lamb_sort(i)*full(max(PF_lamb)))*180/pi);
       end
    end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 nonosc=EIG(SL_num);
if(abs(imag(nonosc))<1e-3)
      disp('-----------------------------------------------------------------------------');
      disp('It is a NON-OSCILLATORY MODE...')
      disp('-----------------------------------------------------------------------------');
else
 slip_PF_lamb=PF_lamb(1:16:ngen*16);
 slip_PF_nr=slip_PF_lamb/max(slip_PF_lamb);
[m_sort_slip  m_ind_slip]=sort(abs(slip_PF_nr));
rev_ind_slip=[ngen:-1:1];
n_ind_slip=m_ind_slip(rev_ind_slip);
n_sort_slip=m_sort_slip(rev_ind_slip);
gen_num=gen(n_ind_slip,1);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
abs_sort=0.05;
k=1;
r_sort_slip(k)=0;
for i=1:ngen
   if(abs(n_sort_slip(i))>abs_sort)
      r_sort_slip(k)=n_sort_slip(i);
      k=k+1;
   end
end
r_ind_slip=n_ind_slip(1:size(r_sort_slip,2));
par_gen_numd=gen_num(1:size(r_sort_slip,2));
   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
V=V1(:,SL_num);
slip_V=V(1:16:ngen*16);
%-----To filter out a mode in which participation of slip is very small----
splip_PF_norm = sparse(zeros(1,nb));
splip_PF_norm(gen(:,1)) = PF_lamb_nr(1:16:ngen*16);
splip_PF_norm_filt = splip_PF_norm(par_gen_numd);
slip_Vpart = 0.3;     %  use this variable to decide the slip part. level
if(sum(abs(splip_PF_norm_filt))< slip_Vpart)
    disp('---------------------------------------------------------------------------------');
    disp('You have chosen a NON-SWING MODE with very low slip participation...')
    disp('---------------------------------------------------------------------------------');
  else
n_slip_Vd=slip_V(r_ind_slip);
n_slip_Vd_nr=n_slip_Vd/max(n_slip_Vd);
[m_sort_slipV  m_ind_slipV]=sort(abs(n_slip_Vd_nr));
rev_ind_slipV=[size(n_slip_Vd,1):-1:1];
n_ind_slipV=m_ind_slipV(rev_ind_slipV);
par_gen_num=par_gen_numd(n_ind_slipV);
n_slip_V=n_slip_Vd(n_ind_slipV);
slip_V_angle=angle(n_slip_V)*180/pi;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
swing_ang=90;
p1=1;
p2=1;
Group1(p1)=par_gen_num(1);
Group2(p2)=0;
a=slip_V_angle(1);
for i=2:size(r_sort_slip,2)
   b=slip_V_angle(i);
   if((abs(a-b)>swing_ang)&(abs(a-b)<(360-swing_ang )))
      Group2(p2)=par_gen_num(i);
      p2=p2+1;
   else
      p1=p1+1;
      Group1(p1)=par_gen_num(i);
   end
end
%------------------------------------------------------------------
disp('----------------------------------------------------------------------------------');
if(sum(Group2)~=0)
   disp('      You have chosen a SWING-MODE       ')
else
   disp('      You have chosen a NON SWING-MODE   ')
end
%-------------------------------------------------------------------


disp('----------------------------------------------------------------------------------');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp('The generator(s) in group-1 is(are) ...       ')
Group1
disp('The generator(s) in group-2 is(are) ...      ')
Group2
disp('----------------------------------------------------------------------------------');

%--------------------------------------------------------------------------

TRACE=input('Enter 1 if you want to plot the compass plot, otherwise 0:   ');
if(TRACE)
  disp('NOTE: Use mouse click on the plot to identify the generator')
  disp('Press any key')
  pause
  r_eig_plot
end
end %It is non swing mode of less participation if loop end
end %It is Non-oscillatory mode if loop end.
rpt=input('Enter 1 if you want to repeat  for another eigenvalue, otherwise 0: ') ;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear Group1 Group2;
clear slip_V_angle;
clear n_slip_Vd;
clear n_slip_V;
clear par_gen_numd;
clear par_gen_num;
clear r_ind_slip;
clear rev_ind;
clear n_ind;
clear p1 p2;
clear r_sort_slip k;
clear gen_num;
clear rev_ind_slipV;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end
```


```data title="TRANS.tbl"
F AC4A_exciter.m                                                                                                                                                                                                  	AC4A_exciter.m
F B_bus_form.m                                                                                                                                                                                                    	B_bus_form.m
F DC1A_exciter.m                                                                                                                                                                                                  	DC1A_exciter.m
F Reheat_turbine.m                                                                                                                                                                                                	Reheat_turbine.m
F busno.dat                                                                                                                                                                                                       	busno.dat
F delPw_pss.dat                                                                                                                                                                                                   	delPw_pss.dat
F delPw_pss_settings.m                                                                                                                                                                                            	delPw_pss_settings.m
F exc_AC4A.dat                                                                                                                                                                                                    	exc_AC4A.dat
F exc_DC1A.dat                                                                                                                                                                                                    	exc_DC1A.dat
F exc_static.dat                                                                                                                                                                                                  	exc_static.dat
F exciter_setting.m                                                                                                                                                                                               	exciter_setting.m
F fdlf_jacob_form.m                                                                                                                                                                                               	fdlf_jacob_form.m
F fdlf_loadflow.m                                                                                                                                                                                                 	fdlf_loadflow.m
F freq_response.m                                                                                                                                                                                                 	freq_response.m
F gen.dat                                                                                                                                                                                                         	gen.dat
F genmat.m                                                                                                                                                                                                        	genmat.m
F hydro_turbine.m                                                                                                                                                                                                 	hydro_turbine.m
F initcond.m                                                                                                                                                                                                      	initcond.m
F ld.dat                                                                                                                                                                                                          	ld.dat
F lfl.dat                                                                                                                                                                                                         	lfl.dat
F lfl_result.m                                                                                                                                                                                                    	lfl_result.m
F load_zip_model.m                                                                                                                                                                                                	load_zip_model.m
F nt.dat                                                                                                                                                                                                          	nt.dat
F pmat.m                                                                                                                                                                                                          	pmat.m
F power_pss.dat                                                                                                                                                                                                   	power_pss.dat
F power_pss_settings.m                                                                                                                                                                                            	power_pss_settings.m
F powerflow.m                                                                                                                                                                                                     	powerflow.m
F primemover_setting.m                                                                                                                                                                                            	primemover_setting.m
F pss_delPw_signal.m                                                                                                                                                                                              	pss_delPw_signal.m
F pss_design.m                                                                                                                                                                                                    	pss_design.m
F pss_power_signal.m                                                                                                                                                                                              	pss_power_signal.m
F pss_selection.m                                                                                                                                                                                                 	pss_selection.m
F pss_slip_signal.m                                                                                                                                                                                               	pss_slip_signal.m
F pvpq.dat                                                                                                                                                                                                        	pvpq.dat
F r_eig_plot.m                                                                                                                                                                                                    	r_eig_plot.m
F report.dat                                                                                                                                                                                                      	report.dat
F shunt.dat                                                                                                                                                                                                       	shunt.dat
F slip_pss.dat                                                                                                                                                                                                    	slip_pss.dat
F slip_pss_settings.m                                                                                                                                                                                             	slip_pss_settings.m
F small_sig.m                                                                                                                                                                                                     	small_sig.m
F static_exciter.m                                                                                                                                                                                                	static_exciter.m
F statld.m                                                                                                                                                                                                        	statld.m
F trace_mode.m                                                                                                                                                                                                    	trace_mode.m
F transtability.mdl                                                                                                                                                                                               	transtability.mdl
F turb_hydro.dat                                                                                                                                                                                                  	turb_hydro.dat
F turb_rhst.dat                                                                                                                                                                                                   	turb_rhst.dat
F yform.m                                                                                                                                                                                                         	yform.m
```


```data title="turb_hydro.dat"
60 1 0.2 0.05 0 1.1 0.1
```


```data title="turb_rhst.dat"
60 0.2 0 0.1 0.05 1.1 0.1 0.3 10 0.4 0.3 0.3 0.4
```


```c title="yform.c"
% Formation of Ybus
% This is called by the main program 

Y=sparse(zeros(nb,nb));

   for i=1:nline
incr=1/(nt(i,3)+j*nt(i,4));
Y((nt(i,1)),(nt(i,2)))=Y((nt(i,1)),(nt(i,2)))-incr;
Y((nt(i,2)),(nt(i,1)))=Y((nt(i,2)),(nt(i,1)))-incr;
Y((nt(i,1)),(nt(i,1)))=Y((nt(i,1)),(nt(i,1)))+incr+j*nt(i,5)/2;
Y((nt(i,2)),(nt(i,2)))=Y((nt(i,2)),(nt(i,2)))+incr+j*nt(i,5)/2;
   end

 for i=nline+1:(ntrans+nline)
incr=1/(nt(i,3)+j*nt(i,4));  
incr1=incr/nt(i,5);
incr2=(1-nt(i,5))*incr/(nt(i,5)*nt(i,5));
incr3=(nt(i,5)-1)*incr/nt(i,5);
Y((nt(i,1)),(nt(i,2)))=Y((nt(i,1)),(nt(i,2)))-incr1;
Y((nt(i,2)),(nt(i,1)))=Y((nt(i,2)),(nt(i,1)))-incr1;
Y((nt(i,1)),(nt(i,1)))=Y((nt(i,1)),(nt(i,1)))+incr1+incr2;
Y((nt(i,2)),(nt(i,2)))=Y((nt(i,2)),(nt(i,2)))+incr1+incr3;   
end

   for i=1:nshunt
incr=shunt(i,2)+j*shunt(i,3);   
Y((shunt(i,1)),(shunt(i,1)))=Y((shunt(i,1)),(shunt(i,1)))+incr;
   end

  
YDQmod=[imag(Y) real(Y);real(Y) -imag(Y)];

l=[1:nb;(nb+1):2*nb];
ll=l(:);
YDQ=YDQmod(ll,ll);

clear YDQmod;
clear incr;

% The following portion is for running transient stability programme

kcir1=QL0./VL0.^2;
kci1=PL0./VL0.^2;

for i=1:nload
incr=kci1(i)-j*kcir1(i);  
Y((ld(i,1)),(ld(i,1)))=Y((ld(i,1)),(ld(i,1)))+incr;
end

for i=1:ngen
incr=1/(j*gen(i,4));
Y((gen(i,1)),(gen(i,1)))=Y((gen(i,1)),(gen(i,1)))+incr;
end
YES=input('Enter 1 if you want to run transtability programme for network disturbances, otherwise 0:  ');
if(YES)
display('If NO action to be taken, PRESS ENTER for any/every prompt.')

Tfault=input('Fault initiation time (s), Tfault= ');
Tclear=input('Fault Duration,(s) Tclear= ');
fbus=input('Faulted Bus: ');
Line=input('Line(s) to be tripped, [ , ]= ');

if isempty(Tclear)
      Tclear=0;
   end
   if isempty(Tfault)
      Tfault=1000;
   end
   
   Yf=Y;
   if ~isempty(fbus)
      Yf(fbus,fbus)=Yf(fbus,fbus)+100000;
   end

Ypf=Y;

% Tripping of lines to clear fault
if ~isempty(Line)
for tline = Line
   if (tline>nline)
% For transformer tripping
  incr=1/(nt(tline,3)+j*nt(tline,4));  
  incr1=incr/nt(tline,5);
  incr2=(1-nt(tline,5))*incr/(nt(tline,5)*nt(tline,5));
  incr3=(nt(tline,5)-1)*incr/nt(tline,5);
  Ypf((nt(tline,1)),(nt(tline,2)))=Ypf((nt(tline,1)),(nt(tline,2)))+incr1;
  Ypf((nt(tline,2)),(nt(tline,1)))=Ypf((nt(tline,2)),(nt(tline,1)))+incr1;
  Ypf((nt(tline,1)),(nt(tline,1)))=Ypf((nt(tline,1)),(nt(tline,1)))-incr1-incr2;
  Ypf((nt(tline,2)),(nt(tline,2)))=Ypf((nt(tline,2)),(nt(tline,2)))-incr1-incr3;   
else
%For line tripping   
  incr=1/(nt(tline,3)+j*nt(tline,4));
  Ypf(nt(tline,1),nt(tline,1))=Ypf(nt(tline,1),nt(tline,1))-incr-j*nt(tline,5)/2;
  Ypf(nt(tline,2),nt(tline,2))=Ypf(nt(tline,2),nt(tline,2))-incr-j*nt(tline,5)/2;
  Ypf(nt(tline,1),nt(tline,2))=Ypf(nt(tline,1),nt(tline,2))+incr;
  Ypf(nt(tline,2),nt(tline,1))=Ypf(nt(tline,2),nt(tline,1))+incr;  
 end
end   
end
purt_vref=0;
else 
purt_vref= input('Enter 1 if you want to run transtability programme for perturbation of VREF, otherwise 0:  ');
if(purt_vref)
 Tclear = 0;
 Tfault = 1000;
 Yf=Y;
 Ypf=Y;
Gvref = input('Enter the generator number whose Vref needs to be perturbed: ');

if(isempty(ng_static))
   ng_static=0;
end

if (sum(ng_static==Gvref)*(~AVR(Gvref))~=0)
    for gg=1:1:size(exc_static,1)
     if(exc_static(gg,1)==Gvref)
        Gn_vref=gg;
     end
    end
  Vref_static(Gn_vref) = Vref_static(Gn_vref) + 0.01;
end

if(isempty(ng_DC1A))
   ng_DC1A=0;
end

if (sum(ng_DC1A==Gvref)*(~AVR(Gvref))~=0)
    for gg=1:1:size(exc_DC1A,1)
     if(exc_DC1A(gg,1)==Gvref)
        Gn_vref=gg;
     end
    end
  Vref_DC1A(Gn_vref) = Vref_DC1A(Gn_vref) + 0.01;
end

if(isempty(ng_AC4A))
   ng_AC4A=0;
end

if (sum(ng_AC4A==Gvref)*(~AVR(Gvref))~=0)
    for gg=1:1:size(exc_AC4A,1)
     if(exc_AC4A(gg,1)==Gvref)
        Gn_vref=gg;
     end
    end
  Vref_AC4A(Gn_vref) = Vref_AC4A(Gn_vref) + 0.01;
end
end %end of purt_verf
end  %end of YES


if((YES+purt_vref)==0)
purt_Tm= input('Enter 1 if you want to run transtability programme for ramping of Tm, otherwise 0:  ');
if(purt_Tm)
 Tclear = 0;
 Tfault = 1000;
 Yf=Y;
 Ypf=Y;
 gtm = input('Enter the generator number whose Tm needs to be ramped-up/down: ');
if(TURB(gtm)==0) 
   disp('The Turbine-governor on the chosen generator must be disabled using TURB selector')
end
else
   gtm=gen(1,1); % adummy number
   Tclear = 0;
   Tfault = 1000;
   Yf=Y;
   Ypf=Y;
end %end of purt_Tm

else
   gtm=gen(1,1); % adummy number
   purt_Tm=0;
end %end 

%------------------------------------------------------------
%For ramping of Tm with a rate of Tm/s.
   [mm order_gen_tm]=sort([setxor(1:nb,gtm),gtm]);
   t_low=1; %time instant of ramping low
   t_up=11; %time instant of ramping up
 %----------------------------------------------------------  
%Frequency Dependancy of Load parameters:: Not accounted.
%please do not tamper this.
kpf=1.5*0;
kqf=2*0;
if((YES+purt_vref+purt_Tm)==0)
   disp('----You can run the transient stability programme without any disturbance----')
end
```