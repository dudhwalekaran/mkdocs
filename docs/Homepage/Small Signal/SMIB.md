```c title="AC4A_Exciter.c"


%----------------------Alternator supplied Controlled Rectifier: AC4A type-

load exc_AC4A.dat

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[mm order_AC4A]=sort([setxor(1:nb,exc_AC4A(:,1)),exc_AC4A(:,1)']);
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


```c title="B_bus_form.c"
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
 1  
0.001 
3
2   
0  
1
0 
1
0
1.0 
60 
```

```c title="DC1A_Exciter.c"
%------------------------IEEE type-1 Exciters-------------------
load exc_DC1A.dat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mm order_DC1A]=sort([setxor(1:nb,exc_DC1A(:,1)),exc_DC1A(:,1)']);
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
1 10 10 10 10 0.01 10 5 1 0 0.1 0.05926 0.03489 0.05926 0.03489 10 0.1 -0.01 
3 10 10 10 10 0.01 10 1000 1 0 0.1 0.06322 0.04452 0.06322 0.04452 10 0.1 -0.1 
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
1 0.02 200 0.015 1.0 12 10 -10 -4.53 5.64 0 
```


```data title="exc_DC1A.dat"
1  0.02 20 0.06 1 10 6.0  -6.0  -0.0485 0.250 3.5461 0.0800 4.7281 0.260 0.0400 1.000
```


```data title="exc_static.dat"
1 200 0.05 -6 6
```


```c title="Exciter_Settings.c"

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


```c title="fdlf_jacob_form.c"%~~~~~~~~~~~~~Form Jacobian and solve for bus voltages and angles~~~~~~~~~~%
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


```c title="fdlf_loadflow.c"clear all;
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


```data title="freq_response.dat"
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
1  1.6     0.32     0.32    6.00   0.05  1.55     1.55    1.55     .44    0.04   5.000   0     
3  0.0006  0.0001   0.0001  10000  0.05  0.0001   0.0001  0.0001   .44    0.04   10000   0   
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
pack;

% Now back to the main program
```


```c title="Hydro_turbine.c"
%speed Governer system with Hydro Turbine.
load turb_hydro.dat
[mm order_hydro]=sort([setxor(1:nb,turb_hydro(:,1)),turb_hydro(:,1)']);
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


[mm order_gen]=sort([setxor(1:nb,gen(:,1)),gen(:,1)']);

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
AVR(3)=1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SL_static(1:nb)=1;

SL_DC1A(1:nb)=1;

SL_AC4A(1:nb)=1;

%~~~~~~~~~~~~~Indicate the generator number on which a specfic type of ~~%
%~~~~~~~~~~~~~~~~~~~~exciter is to be enabled, otherwise null~~~~~~~~~~~~%

%load ng_static.dat
ng_static=[1];

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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SL_HYDRO(1:nb)=1; 

SL_RHST(1:nb)=1;

%-----Indicate the generator number on which a specfic type of --------%
%------- speed-governor-turbine is to be enabled, otherwise null-------%

%load ng_hydro.dat
ng_hydro=[];

%load ng_rht.dat
ng_rht=[];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SL_HYDRO(ng_hydro)=0;

SL_RHST(ng_rht)=0;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%---Main Selector for PSS----
%PSS = 1's (DISABLED); =0's (ENABLED)
%PSS=ones(1,nb);
PSS=zeros(1,nb);
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
3  1.0  0
```


```data title="LFL.dat"
  1 	     1.000000 	    23.573715	    0.999911	    0.208206	 0.000000	 0.000000 
  2 	     0.979002 	    11.786900	    0.000000	    0.000000	 0.000000	 0.000000 
  3 	     1.000000 	    0.0000000 	    0.000000	    0.208206	 1.000000	 0.000000 
```


```c title="LFL_Result.c"


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


```c title="Load_zip_model.c"
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

[mm order_load]=sort([setxor(1:nb,ld(:,1)),ld(:,1)']);
```


```data title="nt.dat"
1  2  0  0.2  0 
2  3  0  0.2  0 
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
1 10 0.05 0.1 0.1 -0.1
```


```c title="power_pss_setting.c"

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


```c title="primemover_setting.c"
%%%%%%%%%%%%%%%%%%%% Primemover_settings %%%%%%%%%%%%%%%%%%%%%%%  
TURB_dash=~(TURB);
TURB_sm=TURB_dash(gen(:,1));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%%%%%%%%%%%%%%%%%%REHEAT TURBINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_rht_tur = zeros(1,nb);
H_rht_tur(gen(:,1)) = H;
H_rht_tur = H_rht_tur(turb_rhst(:,1));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
SL_RHST_dash=~(SL_RHST);
SL_RHST_sm=SL_RHST_dash(gen(:,1));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
RHST_Ag211=zeros(1,nb);
RHST_Ag211(turb_rhst(:,1))=(FHP_rht)*(0.5)./H_rht_tur;
RHST_Ag211=RHST_Ag211(gen(:,1));
RHST_dAg211=SL_RHST_sm.*RHST_Ag211.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag212=zeros(1,nb);
RHST_Ag212(turb_rhst(:,1))=(FIP_rht)*(0.5)./H_rht_tur;
RHST_Ag212=RHST_Ag212(gen(:,1));
RHST_dAg212=SL_RHST_sm.*RHST_Ag212.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag213=zeros(1,nb);
RHST_Ag213(turb_rhst(:,1))=(FLP_rht)*(0.5)./H_rht_tur;
RHST_Ag213=RHST_Ag213(gen(:,1));
RHST_dAg213=SL_RHST_sm.*RHST_Ag213.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag1111=zeros(1,nb);
RHST_Ag1111(turb_rhst(:,1))=-1./TCH_rht;
RHST_Ag1111=RHST_Ag1111(gen(:,1));
RHST_dAg1111=SL_RHST_sm.*RHST_Ag1111.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag1115=zeros(1,nb);
RHST_Ag1115(turb_rhst(:,1))=1./TCH_rht;
RHST_Ag1115=RHST_Ag1115(gen(:,1));
RHST_dAg1115=SL_RHST_sm.*RHST_Ag1115.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag1211=zeros(1,nb);
RHST_Ag1211(turb_rhst(:,1))=1./TRH_rht;
RHST_Ag1211=RHST_Ag1211(gen(:,1));
RHST_dAg1211=SL_RHST_sm.*RHST_Ag1211.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag1212=zeros(1,nb);
RHST_Ag1212(turb_rhst(:,1))=-1./TRH_rht;
RHST_Ag1212=RHST_Ag1212(gen(:,1));
RHST_dAg1212=SL_RHST_sm.*RHST_Ag1212.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag1312=zeros(1,nb);
RHST_Ag1312(turb_rhst(:,1))=1./TCO_rht;
RHST_Ag1312=RHST_Ag1312(gen(:,1));
RHST_dAg1312=SL_RHST_sm.*RHST_Ag1312.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag1313=zeros(1,nb);
RHST_Ag1313(turb_rhst(:,1))=-1./TCO_rht;
RHST_Ag1313=RHST_Ag1313(gen(:,1));
RHST_dAg1313=SL_RHST_sm.*RHST_Ag1313.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag142=zeros(1,nb);
RHST_Ag142(turb_rhst(:,1))=K_rht.*(1-T2_rht./T1_rht)./T1_rht;
RHST_Ag142=RHST_Ag142(gen(:,1));
RHST_dAg142=SL_RHST_sm.*RHST_Ag142.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag1414=zeros(1,nb);
RHST_Ag1414(turb_rhst(:,1))=-1./T1_rht;
RHST_Ag1414=RHST_Ag1414(gen(:,1));
RHST_dAg1414=SL_RHST_sm.*RHST_Ag1414.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag152=zeros(1,nb);
RHST_Ag152(turb_rhst(:,1))=-K_rht.*T2_rht./T1_rht./T3_rht;
RHST_Ag152=RHST_Ag152(gen(:,1));
RHST_dAg152=SL_RHST_sm.*RHST_Ag152.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag1514=zeros(1,nb);
RHST_Ag1514(turb_rhst(:,1))=-1./T3_rht;
RHST_Ag1514=RHST_Ag1514(gen(:,1));
RHST_dAg1514=SL_RHST_sm.*RHST_Ag1514.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RHST_Ag1515=zeros(1,nb);
RHST_Ag1515(turb_rhst(:,1))=-1./T3_rht;
RHST_Ag1515=RHST_Ag1515(gen(:,1));
RHST_dAg1515=SL_RHST_sm.*RHST_Ag1515.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%%%%%%%%%%%%%%%%HYDRO TURBINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_hyd_tur = zeros(1,nb);
H_hyd_tur(gen(:,1)) = H;
H_hyd_tur = H_hyd_tur(turb_hydro(:,1));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
SL_HYDRO_dash=~(SL_HYDRO);
SL_HYDRO_sm=SL_HYDRO_dash(gen(:,1));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HYDRO_Ag215=zeros(1,nb);
HYDRO_Ag215(turb_hydro(:,1))=-1./H_hyd_tur;
HYDRO_Ag215=HYDRO_Ag215(gen(:,1));
HYDRO_dAg215=SL_HYDRO_sm.*HYDRO_Ag215.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HYDRO_Ag216=zeros(1,nb);
HYDRO_Ag216(turb_hydro(:,1))=(0.5)./H_hyd_tur;
HYDRO_Ag216=HYDRO_Ag216(gen(:,1));
HYDRO_dAg216=SL_HYDRO_sm.*HYDRO_Ag216.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HYDRO_Ag142=zeros(1,nb);
HYDRO_Ag142(turb_hydro(:,1))=K_ht.*(1-T2_ht./T1_ht)./T1_ht;
HYDRO_Ag142=HYDRO_Ag142(gen(:,1));
HYDRO_dAg142=SL_HYDRO_sm.*HYDRO_Ag142.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HYDRO_Ag1414=zeros(1,nb);
HYDRO_Ag1414(turb_hydro(:,1))=-1./T1_ht;
HYDRO_Ag1414=HYDRO_Ag1414(gen(:,1));
HYDRO_dAg1414=SL_HYDRO_sm.*HYDRO_Ag1414.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HYDRO_Ag152=zeros(1,nb);
HYDRO_Ag152(turb_hydro(:,1))=-K_ht.*T2_ht./T1_ht./T3_ht;
HYDRO_Ag152=HYDRO_Ag152(gen(:,1));
HYDRO_dAg152=SL_HYDRO_sm.*HYDRO_Ag152.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HYDRO_Ag1514=zeros(1,nb);
HYDRO_Ag1514(turb_hydro(:,1))=-1./T3_ht;
HYDRO_Ag1514=HYDRO_Ag1514(gen(:,1));
HYDRO_dAg1514=SL_HYDRO_sm.*HYDRO_Ag1514.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HYDRO_Ag1515=zeros(1,nb);
HYDRO_Ag1515(turb_hydro(:,1))=-1./T3_ht;
HYDRO_Ag1515=HYDRO_Ag1515(gen(:,1));
HYDRO_dAg1515=SL_HYDRO_sm.*HYDRO_Ag1515.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HYDRO_Ag1615=zeros(1,nb);
HYDRO_Ag1615(turb_hydro(:,1))=6./T_wht;
HYDRO_Ag1615=HYDRO_Ag1615(gen(:,1));
HYDRO_dAg1615=SL_HYDRO_sm.*HYDRO_Ag1615.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HYDRO_Ag1616=zeros(1,nb);
HYDRO_Ag1616(turb_hydro(:,1))=-2./T_wht;
HYDRO_Ag1616=HYDRO_Ag1616(gen(:,1));
HYDRO_dAg1616=SL_HYDRO_sm.*HYDRO_Ag1616.*TURB_sm;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
[mm order_delPw_pss]=sort([setxor(1:nb,delPw_pss(:,1)),delPw_pss(:,1)']);

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
[mm order_power_pss]=sort([setxor(1:nb,power_pss(:,1)),power_pss(:,1)']);
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


[mm order_slip_pss]=sort([setxor(1:nb,slip_pss(:,1)),slip_pss(:,1)']);
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


```data title="pvpq.data"
3  1.0   0
3  1.0   0
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


```c title="Reheat_turbine.c"%Speed governer with reheat type turbine

load turb_rhst.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mm order_rhst]=sort([setxor(1:nb,turb_rhst(:,1)),turb_rhst(:,1)']);
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


```data title="report.dat"
Iter. No = 1 Max. Real power mismatch at bus = 3, Max. mismatch = 0.02645414
Iter. No = 1 Max. Reactive power mismatch at bus = 2, Max. mismatch = 0.01043612
Iter. No = 2 Max. Real power mismatch at bus = 3, Max. mismatch = 0.00214322
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Converged Real power mismatch iteration No = 3, Max.mismatch = 0.00008915
Converged Reactive power mismatch iteration No = 2, Max. mismatch = 0.00085633
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Loadflow results:
---------------------------------------------------------------------------------------

Bus No     VbO          thetaO           PGO             QGO         PLO         QLO 
---------------------------------------------------------------------------------------
  1 	     1.000000 	    0.000000	    0.999911	    0.208206	 0.000000	 0.000000 
  2 	     0.979002 	  -11.786858	    0.000000	    0.000000	 0.000000	 0.000000 
  3 	     1.000000 	  -23.573715	    0.000000	    0.208206	 1.000000	 0.000000 
---------------------------------------------------------------------------------------
Line flows:
----------------------------------------------------------------------------------
                       Line flows                                Line flows 
                   _____________________                   _______________________
 From   To         P-flow         Q-flow    From   To      P-flow           Q-flow
----------------------------------------------------------------------------------
  1 	   2	       0.9999	       0.2082	   2	    1	      -0.9999	       0.0004
  2 	   3	       0.9999	       0.0004	   3	    2	      -0.9999	       0.2082
----------------------------------------------------------------------------------
Total real power losses in the system =  0.000000
Total reactive power losses in the system =  0.417269

```

```data title="result_lfl.dat"
  1 	     1.000000 	    0.000000	    0.999911	    0.208206	 0.000000	 0.000000 
  2 	     0.979002 	  -11.786858	    0.000000	    0.000000	 0.000000	 0.000000 
  3 	     1.000000 	  -23.573715	    0.000000	    0.208206	 1.000000	 0.000000 
```

```data title="slip_pss.dat"
1 16  0.02  2 0.08 0.027 0.05 -0.05 570 35 1
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

[mm order_static]=sort([setxor(1:nb,exc_static(:,1)),exc_static(:,1)']);

kA_static=exc_static(:,2)';
TA_static=exc_static(:,3)';

Efd_min_static=exc_static(:,4)';
Efd_max_static=exc_static(:,5)';

Efd0_static=EFD0(exc_static(:,1));

Vref_static=Efd0_static./kA_static + lfl(exc_static(:,1),2)';
```


```c title="StatId.c"

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


```c title="Trace_Mode.c"
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
F Result_lfl.dat                                                                                                                                                                                                  	Result_lfl.dat
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


```data title="Turb_hydro.dat"
1 1 0.2 0.05 0 1.1 0.1
```


```data title="turb_rhst.dat"
1 0.2 0 0.1 0.05 1.1 0.1 0.3 10 0.4 0.3 0.3 0.4
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