``` octave title="B_bus_form.m"
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

``` octave title="decdata.m"
%--------------------------Load the HVDC link data-------------------------%
load hvdc.dat;
load spdata.dat;

Vdor = sparse(zeros(nl,1));
Vdoi = sparse(zeros(nl,1));
Vdr = sparse(zeros(nl,1));
Vdi = sparse(zeros(nl,1));
Alpha = sparse(zeros(nl,1));
Gama = sparse(zeros(nl,1));
Phi_r = sparse(zeros(nl,1));
Phi_i = sparse(zeros(nl,1));
Pr =  sparse(zeros(nl,1));
Pi =  sparse(zeros(nl,1));
Qr =  sparse(zeros(nl,1));
Qi =  sparse(zeros(nl,1));
Id =  sparse(zeros(nl,1));


[Lno LIn] = sort(spdata(:,1));
Link_dat = spdata(LIn,:);

Re_bus  = Link_dat(:,2);  %List of Rectifier bus nos
In_bus  = Link_dat(:,3);  %List of Inverter bus nos 
link    = [Re_bus In_bus]; 

Re_Ran_data = hvdc((1:2:(2*nl-1)),:);
In_Ran_data = hvdc((2:2:2*nl),:);

for i=1:nl
 Re_data(i,:) = Re_Ran_data(Re_Ran_data(:,1)==i,:); 
 In_data(i,:) = In_Ran_data(In_Ran_data(:,1)==i,:); 
end

Br = Re_data(:,3);
Tr_max = Re_data(:,4);
Tr_min = Re_data(:,5);
Tr = Re_data(:,6);
Xcr = Re_data(:,7); 
Al_min = (pi/180)*Re_data(:,8); 
Tpr = Re_data(:,9); 
Bfr = Re_data(:,10);

Bi = In_data(:,3);
Ti_max = In_data(:,4);
Ti_min = In_data(:,5);
Ti = In_data(:,6);
Xci = In_data(:,7);
Ga_min = (pi/180)*In_data(:,8);
Tpi = In_data(:,9);
Bfi = In_data(:,10);

Vdis = Link_dat(:,4);
Idcs = Link_dat(:,5);
Pdcs = Link_dat(:,6);
Im   = Link_dat(:,7);
Rdc  = Link_dat(:,8);
Al_Lr = (pi/180)*Link_dat(:,9);
Al_Up = (pi/180)*Link_dat(:,10); 
cn = Link_dat(:,11); 
pn = Link_dat(:,12); 
Vb = Link_dat(:,13); 
Vdrn = Link_dat(:,14);

Rcr = 3/pi*Xcr.*Br;
Rci = 3/pi*Xci.*Bi;
kr = (3*sqrt(2)/pi)*Br;  
ki = (3*sqrt(2)/pi)*Bi; 

del_Tr = (Tr_max - Tr_min)./Tpr; 
del_Ti = (Ti_max - Ti_min)./Tpi; 

if(nshunt==0)
  shunt=[];
end
nshunt = nshunt + 2*nl;
shunt = [shunt;Re_bus zeros(nl,1) Bfr; In_bus zeros(nl,1) Bfi];

Idc =  sparse(zeros(nl,1));

Isp_link = Link_dat(cn==1,1);
Psp_link = Link_dat(pn==1,1);

Idc(Isp_link) = Idcs(Isp_link);
Vdr(Isp_link) = Vdrn(Isp_link);
Vdi = Vdis;
```

``` octave title="dcflow.m"
%~~~~~~~~~~~~~~~Form Real & Reactive power loads from HVDC Link side~~~~~~~~~~~~~~~~~%

if (sum(Psp_link)~=0)
 SQ_term = sqrt(Vdi(Psp_link).*Vdi(Psp_link) + 4*Pdcs(Psp_link).*Rdc(Psp_link));
 Idc(Psp_link) = (-Vdi(Psp_link) + SQ_term)./(2*Rdc(Psp_link));
 Vdr(Psp_link) = Pdcs(Psp_link)./Idc(Psp_link);
end

Vdor_Ms = ((kr.*Tr).*Vmag(Re_bus));
Al_Ms = acos((Vdr+Rcr.*Idc)./Vdor_Ms); 

%-----------------------Identifiying the Mode of operation--------------------------%

mod_con = Al_Ms >=Al_min;
mod1_link = find(mod_con==1);
mod2_link = find(mod_con==0);

%---------------------------Entering in Mode 1--------------------------------------%

mode_bit= sparse(zeros(nl,1));

if (sum(mod1_link)~=0)
 mode_bit(mod1_link)=1;  
 Gama(mod1_link) = Ga_min(mod1_link);
 Id(mod1_link) = Idc(mod1_link); 
 Vdi(mod1_link) = Vdis(mod1_link);
 
 mod1_rbus =Re_bus(mod1_link); % List of rectifier buses on mode-1
 mod1_ibus =In_bus(mod1_link); % List of Inverter buses on mode-1
 
%------------------Calculating Inverter side Parameters in mode-1--------------------%

 Ti(mod1_link) = (Vdi(mod1_link) + Rci(mod1_link).*Id(mod1_link))./((ki(mod1_link).*Vmag(mod1_ibus)).*cos(Ga_min(mod1_link)));
 Vdoi(mod1_link) = (ki(mod1_link).*Ti(mod1_link)).*Vmag(mod1_ibus);
 
%-----------Perform the following if Inverter T/f Tap limits are accounted-----------%

 if ((sum(Ti(mod1_link)>Ti_max(mod1_link))~=0)|(sum(Ti(mod1_link)<Ti_min(mod1_link))~=0))
  Ti_Uvio_All = sparse(zeros(nl,1));   
  Ti_Lvio_All = sparse(zeros(nl,1));
  mod1_sp_link = sparse(zeros(nl,1));   
  Ti_Uvio_con =find(Ti>Ti_max);
  Ti_Lvio_con =find(Ti<Ti_min);
  
  Ti_Uvio_All(Ti_Uvio_con) = Ti_Uvio_con;
  Ti_Lvio_All(Ti_Lvio_con) = Ti_Lvio_con;
  mod1_sp_link(mod1_link) = mod1_link;
  
  Ti_Uvio_int =mod1_sp_link(mod1_sp_link==Ti_Uvio_All) ;
  Ti_Lvio_int =mod1_sp_link(mod1_sp_link==Ti_Lvio_All) ;
    
  if (sum(Ti_Uvio_int)~=0)
   Ti_Uvio_link = Ti_Uvio_int(Ti_Uvio_int~=0);
  else 
   Ti_Uvio_link=[];
  end
  
  if (sum(Ti_Lvio_int)~=0)
    Ti_Lvio_link = Ti_Lvio_int(Ti_Lvio_int~=0);
  else
     Ti_Lvio_link=[];
  end
  
  if (sum(Ti_Uvio_link)~=0)
   Ti(Ti_Uvio_link) = Ti_max(Ti_Uvio_link) ;
  end
  if (sum(Ti_Lvio_link)~=0)
   Ti(Ti_Lvio_link) = Ti_min(Ti_Lvio_link) ;
  end
  Ti_vio_link = [Ti_Uvio_link;Ti_Lvio_link] ;
  Vdoi(mod1_link) = (ki(mod1_link).*Ti(mod1_link)).*Vmag(mod1_ibus);
  Vdi(Ti_vio_link) = Vdoi(Ti_vio_link).*cos(Gama(Ti_vio_link))- Rci(Ti_vio_link).*Id(Ti_vio_link);
 end%--------------- End of Inverter T/f tap limits setting--------------------------%
 
 Phi_i(mod1_link) = acos(Vdi(mod1_link)./Vdoi(mod1_link));
 Pi(mod1_link) = Vdi(mod1_link).*Id(mod1_link);
 Qi(mod1_link) = Pi(mod1_link).*tan(Phi_i(mod1_link));
 Vdr(mod1_link) = Vdi(mod1_link) + Rdc(mod1_link).*Id(mod1_link);
 Vdor(mod1_link) = (kr(mod1_link).*Tr(mod1_link)).*Vmag(mod1_rbus);
 Alpha(mod1_link) = acos((Vdr(mod1_link)+Rcr(mod1_link).*Id(mod1_link))./Vdor(mod1_link));
 
 
%------------Perform the following if Alpha-limits are Accounted-----------------------% 
 
 Al_Uvio_link = find(Alpha(mod1_link)>Al_Up(mod1_link));
 Al_Lvio_link = find(Alpha(mod1_link)<Al_Lr(mod1_link));
 Re_Uvio_bus = Re_bus(Al_Uvio_link);
 Re_Lvio_bus = Re_bus(Al_Lvio_link); 
 
 if(sum(Al_Uvio_link)~=0)
  while (Alpha(Al_Uvio_link) > Al_Up(Al_Uvio_link))
   Tr(Al_Uvio_link) = Tr(Al_Uvio_link) - del_Tr(Al_Uvio_link);
   Vdor(Al_Uvio_link) = (kr(Al_Uvio_link).*Tr(Al_Uvio_link)).*Vmag(Re_Uvio_bus);
   Alpha(Al_Uvio_link) = acos((Vdr(Al_Uvio_link)+Rcr(Al_Uvio_link).*Id(Al_Uvio_link))./Vdor(Al_Uvio_link));
  end
 end
 
 if(sum(Al_Lvio_link)~=0)
  while (Alpha(Al_Lvio_link)< Al_Lr(Al_Lvio_link))
   Tr(Al_Lvio_link) = Tr(Al_Lvio_link) + del_Tr(Al_Lvio_link);
   Vdor(Al_Lvio_link) = (kr(Al_Lvio_link).*Tr(Al_Lvio_link)).*Vmag(Re_Lvio_bus);
   Alpha(Al_Lvio_link) = acos((Vdr(Al_Lvio_link)+Rcr(Al_Lvio_link).*Id(Al_Lvio_link))./Vdor(Al_Lvio_link));
  end
 end %--------Alpha-Limits parts ends here--------------------------------------------%
 
 %-----------Perform the following if Rectifier T/f Tap limits are accounted----------%

 if ((sum(Tr(mod1_link)>Tr_max(mod1_link))~=0)|(sum(Tr(mod1_link)<Tr_min(mod1_link))~=0))
  Tr_Uvio_All = sparse(zeros(nl,1));   
  Tr_Lvio_All = sparse(zeros(nl,1));
  mod1_sp_link = sparse(zeros(nl,1));   
  
  Tr_Uvio_con =find(Tr>Tr_max);
  Tr_Lvio_con =find(Tr<Tr_min);
  
  Tr_Uvio_All(Tr_Uvio_con) = Tr_Uvio_con;
  Tr_Lvio_All(Tr_Lvio_con) = Tr_Lvio_con;
  mod1_sp_link(mod1_link) = mod1_link;
  
  Tr_Uvio_int =mod1_sp_link(mod1_sp_link==Tr_Uvio_All) ;
  Tr_Lvio_int =mod1_sp_link(mod1_sp_link==Tr_Lvio_All) ;
  
  if (sum(Tr_Uvio_int)~=0)
   Tr_Uvio_link = Tr_Uvio_int(Tr_Uvio_int~=0);
  else
   Tr_Uvio_link=[];
  end   
  if (sum(Tr_Lvio_int)~=0)
   Tr_Lvio_link = Tr_Lvio_int(Tr_Lvio_int~=0);
  else
   Tr_Lvio_link=[];
  end
  
  if (sum(Tr_Uvio_link)~=0)
   Tr(Tr_Uvio_link) = Tr_max(Tr_Uvio_link) ;
  end
  
  if (sum(Tr_Lvio_link)~=0)
   Tr(Tr_Lvio_link) = Tr_min(Tr_Lvio_link) ;
  end
  Tr_vio_link = [Tr_Uvio_link;Tr_Lvio_link];
  Vdor(mod1_link) = (kr(mod1_link).*Tr(mod1_link)).*Vmag(mod1_rbus);
  Alpha(Tr_vio_link) = acos((Vdr(Tr_vio_link)+Rcr(Tr_vio_link).*Id(Tr_vio_link))./Vdor(Tr_vio_link));
 end %--------------------------Alpha-Limits parts ends here-------------------------%
 Phi_r(mod1_link) = acos(Vdr(mod1_link)./Vdor(mod1_link));
 Pr(mod1_link) = Vdr(mod1_link).*Id(mod1_link);
 Qr(mod1_link) = Pr(mod1_link).*tan(Phi_r(mod1_link));
end
  
%------------------------------- Entering in Mode 2----------------------------------%

if(sum(mod2_link)~=0)
 mode_bit(mod2_link)=2;
 Alpha(mod2_link) = Al_min(mod2_link);
 
 mod2_sp_link = sparse(zeros(nl,1));
 Psp_sp_link = sparse(zeros(nl,1)); 
 mod2_sp_link(mod2_link)  = mod2_link;
 Psp_sp_link(Psp_link) = Psp_link;
 
 mod2_Psp_con = mod2_sp_link(mod2_sp_link==Psp_sp_link);
 
 if (sum(mod2_Psp_con)~=0)
   mod2_Psp_link = mod2_Psp_con(mod2_Psp_con~=0);  % list of mode2 Pspecified link nos
 else
   mod2_Psp_link=[];
 end
 
 mod2_Isp_con = mod2_sp_link(mod2_sp_link~=Psp_sp_link);% list of mode2 Ispecified link nos
 if (sum(mod2_Isp_con)~=0)
   mod2_Isp_link = mod2_Isp_con(mod2_Isp_con~=0);  % list of mode2 Pspecified link nos
 else
   mod2_Isp_link=[];
 end
 
 mod2_rbus =Re_bus(mod2_link); % List of rectifier buses on mode-2
 mod2_ibus =In_bus(mod2_link); % List of Inverter buses on mode-2
 
 Vdor(mod2_link) = (kr(mod2_link).*Tr(mod2_link)).*Vmag(mod2_rbus);
 
 if (sum(mod2_Isp_link)~=0)
  Id(mod2_Isp_link) = Idc(mod2_Isp_link)-Im(mod2_Isp_link);  
  Vdr(mod2_Isp_link) = Vdor(mod2_Isp_link).*cos(Alpha(mod2_Isp_link))- Rcr(mod2_Isp_link).*Id(mod2_Isp_link) ;
  Pr(mod2_Isp_link) = Vdr(mod2_Isp_link).*Id(mod2_Isp_link);
 end
 
 if (sum(mod2_Psp_link)~=0)
   Pr(mod2_Psp_link) = Pdcs(mod2_Psp_link);  
   Vdor_cp = Vdor(mod2_Psp_link).*cos(Alpha(mod2_Psp_link));
   Id(mod2_Psp_link) =  (Vdor_cp - sqrt(Vdor_cp.*Vdor_cp - 4*Pr(mod2_Psp_link).*Rcr(mod2_Psp_link)))./(2*Rcr(mod2_Psp_link));
   Vdr(mod2_Psp_link) =  Pr(mod2_Psp_link)./Id(mod2_Psp_link);
 end
  
 Phi_r(mod2_link) = acos(Vdr(mod2_link)./Vdor(mod2_link));
 Qr(mod2_link) = Pr(mod2_link).*tan(Phi_r(mod2_link));
  
 Vdi(mod2_link) = Vdr(mod2_link)- Rdc(mod2_link).*Id(mod2_link);
 Vdoi(mod2_link) = (ki(mod2_link).*Ti(mod2_link)).*Vmag(mod2_ibus);
 Gama(mod2_link) = acos((Vdi(mod2_link)+ Rci(mod2_link).*Id(mod2_link))./Vdoi(mod2_link));
  
%------------Perform the following if Gama-limits are Accounted----------------------% 
 
 Ga_vio_link = find(Gama(mod2_link) < Ga_min(mod2_link));
 Ga_vio_bus = In_bus(Ga_vio_link);
 
 if(sum(Ga_vio_link)~=0)
  while (Gama(mod2_link) < Ga_min(mod2_link))  
   Ti(Ga_vio_link) = Ti(Ga_vio_link) - del_Ti(Ga_vio_link);
   Vdoi(Ga_vio_link) = (ki(Ga_vio_link).*Ti(Ga_vio_link)).*Vmag(Ga_vio_bus);
   Gama(Ga_vio_link) = acos((Vdi(Ga_vio_link)+ Rci(Ga_vio_link).*Id(Ga_vio_link))./Vdoi(Ga_vio_link));
  end
 end   %--------------------------Gama-Limit parts ends here-------------------------%
 
 
%-----------Perform the following if Inverter T/f Tap limits are accounted-----------%
 
 if ((sum(Ti(mod2_link)>Ti_max(mod2_link))~=0)|(sum(Ti(mod2_link)<Ti_min(mod2_link))~=0))
  Ti_Uvio_All = sparse(zeros(nl,1));   
  Ti_Lvio_All = sparse(zeros(nl,1));
  mod2_sp_link = sparse(zeros(nl,1));   
  
  Ti_Uvio_con =find(Ti>Ti_max);
  Ti_Lvio_con =find(Ti<Ti_min);
  
  Ti_Uvio_All(Ti_Uvio_con) = Ti_Uvio_con;
  Ti_Lvio_All(Ti_Lvio_con) = Ti_Lvio_con;
  mod2_sp_link(mod2_link) = mod2_link;
  
  Ti_Uvio_int =mod2_sp_link(mod2_sp_link==Ti_Uvio_All) ;
  Ti_Lvio_int =mod2_sp_link(mod2_sp_link==Ti_Lvio_All) ;
  
  if (sum(Ti_Uvio_int)~=0)
   Ti_Uvio_link = Ti_Uvio_int(Ti_Uvio_int~=0);
  else 
   Ti_Uvio_link=[];
  end
  
  if (sum(Ti_Lvio_int)~=0)
    Ti_Lvio_link = Ti_Lvio_int(Ti_Lvio_int~=0);
  else
     Ti_Lvio_link=[];
  end
  
  if (sum(Ti_Uvio_link)~=0)
   Ti(Ti_Uvio_link) = Ti_max(Ti_Uvio_link) ;
  end
  if (sum(Ti_Lvio_link)~=0)
   Ti(Ti_Lvio_link) = Ti_min(Ti_Lvio_link) ;
  end
  Ti_vio_link = [Ti_Uvio_link;Ti_Lvio_link] ;
  Vdoi(mod2_link) = (ki(mod2_link).*Ti(mod2_link)).*Vmag(mod2_ibus);
  Gama(Ti_vio_link) = acos((Vdi(Ti_vio_link)+ Rci(Ti_vio_link).*Id(Ti_vio_link))./Vdoi(Ti_vio_link));
 end%---------- End of Inverter T/f tap limits setting -----------------------------%
 
 Vdi_vio_link = find(Vdi>Vdoi);
 if (sum(Vdi_vio_link)~=0)
   Vdoi(Vdi_vio_link) = Vdi(Vdi_vio_link);
 end
 Phi_i(mod2_link) = acos(Vdi(mod2_link)./Vdoi(mod2_link)) ;
 Pi(mod2_link) = Vdi(mod2_link).*Id(mod2_link);
 Qi(mod2_link) = Pi(mod2_link).*tan(Phi_i(mod2_link));
end
```

``` octave title="powerflow.m"
%-----------------------------Calculating Power flows--------------------------------%
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

if (nl~=0)
 S_line_loss = sum(Spq+Sqp) + sum(S_sh) + sum(Id.^2.*Rdc);
else 
 S_line_loss = sum(Spq+Sqp) + sum(S_sh);
end


%--------------------------------------------------------------------------------------%          
```

``` octave title="fdlf_jacob_form.m"
%~~~~~~~~~~~~~Form Jacobian and solve for bus voltages and angles~~~~~~~~~~%
Vmag=sparse(ones(nb,1));
Vang=sparse(zeros(nb,1));
Vmag(pv_data(:,1))=pv_data(:,2);  %Replaces the voltages of P-V buses in Vmag 
Vmag(sl) = Vsl;
Vang(sl) = 0;
Vo=Vmag;
kp = 0;
kq = 0;
fid=fopen('report.dat', 'w');
fprintf(fid,'                   Detailed Report of Load flow (FDLF method)\n' ); 
fprintf(fid,'                   ----------------------------\n' ); 

%------------------Iteration begins from here onwards---------------------%
for k=1:300
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
%---------------------Load modelling starts here-------------------------%
   
   Pdc =  sparse(zeros(nb,1));
   Qdc =  sparse(zeros(nb,1));
   if (nl~=0)
    dcflow;
    Pdc([Re_bus;In_bus]) = [Pr;-Pi];
    Qdc([Re_bus;In_bus]) = [Qr;Qi];
   end 
   
   if (Load_bit==0&nl~=0)
    Psp = sparse(zeros(nb,1));
    Psp(pv_data(:,1)) = pv_data(:,3);
    Psp(pq_data(:,1)) = Psp(pq_data(:,1))-pq_data(:,2);
    Qsp = sparse(zeros(nb,1));
    Qsp(pq_data(:,1)) = -pq_data(:,3);
    Qload =sparse(zeros(nb,1));
    Qload(pq_data(:,1))= pq_data(:,3);
   
    Psp([Re_bus;In_bus]) = Psp([Re_bus;In_bus])- Pdc([Re_bus;In_bus]); 
    Qsp([Re_bus;In_bus]) = Qsp([Re_bus;In_bus])- Qdc([Re_bus;In_bus]); 
    Qload([Re_bus;In_bus]) = Qload([Re_bus;In_bus]) + Qdc([Re_bus;In_bus]);
   end

   if(Load_bit~=0)
    PL = sparse(zeros(nb,1));
    PL_cp = a1.*PLo;
    PL_cc = a2.*PLo.*(Vmag./Vo);
    PL_ci = a3.*PLo.*(Vmag./Vo).^2;
    PL = PL_cp + PL_cc + PL_ci; 
   
    QL = sparse(zeros(nb,1));
    QL_cp = b1.*QLo;
    QL_cc = b2.*QLo.*(Vmag./Vo);
    QL_ci = b3.*QLo.*(Vmag./Vo).^2;
    QL = QL_cp + QL_cc + QL_ci; 
   
    Psp = sparse(zeros(nb,1));
    Psp(pv_data(:,1)) = pv_data(:,3);
    Psp(pq_data(:,1)) = Psp(pq_data(:,1))-PL(pq_data(:,1));
    Psp([Re_bus;In_bus])= Psp([Re_bus;In_bus])- Pdc([Re_bus;In_bus]);

    Qsp = sparse(zeros(nb,1));
    Qsp(pq_data(:,1)) = -QL(pq_data(:,1));
    Qsp([Re_bus;In_bus]) =Qsp([Re_bus;In_bus]) - Qdc([Re_bus;In_bus]);
   end    %-----------Load modelling ends here---------%
   
   if (Freq_bit~= 0)
    if(Load_bit==0) 
     PL = PLo;
     QL = QLo;
    end
    
    Psp_sl  = PG_sl - PL(sl);
    delP_D = Pc(sl) - Psp_sl;
    delF = delF + (-delP_D*R);
     
    PL = PL*(1 + Kpf*delF);
    QL = QL*(1 + Kqf*delF);
     
    Psp = sparse(zeros(nb,1));
    Psp(pv_data(:,1)) = pv_data(:,3);
    Psp(pq_data(:,1)) = Psp(pq_data(:,1))-PL(pq_data(:,1));
    Psp([Re_bus;In_bus])= Psp([Re_bus;In_bus])- Pdc([Re_bus;In_bus]);

    Qsp = sparse(zeros(nb,1));
    Qsp(pq_data(:,1)) = -QL(pq_data(:,1));
    Qsp([Re_bus;In_bus]) =Qsp([Re_bus;In_bus]) - Qdc([Re_bus;In_bus]);
   end
   
 %------------------------bus power mis matches -------------------------%  
   delP=Psp-Pc;
   delP(sl,:)=[];
   delQ=Qsp-Qc;
   delQ(pv_sl_num,:)=[];
   
   if (Freq_bit~= 0)
    Psp_sl  = PG_sl - PL(sl);
    delP_sl = Pc(sl) - Psp_sl;
    nosl_check = (max(abs(delP))<= tole);
    sl_check =  (abs(delP_sl)<= tole);
    del_P_sl = [delP;delP_sl];
   else
    nosl_check = (max(abs(delP))<= tole);
    sl_check = 1;
    del_P_sl = delP;
   end

   if (nosl_check~=0 & sl_check~=0)
     if (max(abs(delQ))<= tole)
      fprintf(fid,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
      fprintf(fid,'Converged Real power mismatch iteration No = %i, Max.mismatch = %10.8f\n',kp, full(max(abs(del_P_sl))));
      fprintf(fid,'Converged Reactive power mismatch iteration No = %i, Max. mismatch = %10.8f\n',kq, full(max(abs(delQ))));
      convergence_bit=1;
      if(sum(pvpq(:,1)==sl)==1&(Load_bit~=0|Freq_bit~=0))
       S(sl)= S(sl)+PL(sl)+j*QL(sl)+Pdc(sl)+j*Qdc(sl);
      else
       S(sl)= S(sl)+(sl_load(1,2)+j*sl_load(1,3))+Pdc(sl)+j*Qdc(sl);
      end
      fprintf(fid,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
      break;
     end
   else
    kp=kp+1;
    if(k~=1)
     del_Pm=sparse(zeros(nb,1)); 
     del_Pm(sort([pv_num;num_no_sl_pv]))=delP;
     
     if (Freq_bit~= 0)
       del_Pm(sl) = delP_sl;
     end
      
     Bno=[1:nb]';
     vio_busP = Bno(abs(del_Pm)==max(abs(del_Pm)));
     fprintf(fid,'Iter. No = %i Max. Real power mismatch at bus = %i, Max. mismatch = %10.8f\n',(kp-1),full(vio_busP),full(max(abs(del_Pm))));
    end
    
 %---------------------Updation of bus angles at each bus------------------% 
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
    
    Qstart_cri = 0.1;
    if (Q_bit~=0 & max(abs(delQ))<=Qstart_cri)
     if(Load_bit~=0 | Freq_bit~= 0)  
      Qc(pvpq_buses) = Qc(pvpq_buses) + QL(pvpq_buses); 
      Qc(pvpq_buses) = Qc(pvpq_buses) + Qdc(pvpq_buses);  
     else
      Qc(pvpq_buses) = Qc(pvpq_buses) + Qload(pvpq_buses); 
     end
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
     if(Load_bit~=0 | Freq_bit~= 0)
      Qsp(pq_data(:,1)) = -QL(pq_data(:,1));
      Qsp([Re_bus;In_bus]) = Qsp([Re_bus;In_bus])- Qdc([Re_bus;In_bus]); 
     else  
      Qsp(pq_data(:,1)) = -pq_data(:,3);
      Qsp([Re_bus;In_bus]) = Qsp([Re_bus;In_bus])- Qdc([Re_bus;In_bus]); 
     end 
     
     Vmag(pv_num)=pv_Vmag(pv_num); 
     if(Load_bit~=0)
      PL_cp(pvpq_buses) = a1(pvpq_buses).*PLo(pvpq_buses);
      PL_cc(pvpq_buses) = a2(pvpq_buses).*PLo(pvpq_buses).*(Vmag(pvpq_buses)./Vo(pvpq_buses));
      PL_ci(pvpq_buses) = a3(pvpq_buses).*PLo(pvpq_buses).*(Vmag(pvpq_buses)./Vo(pvpq_buses)).^2;
      PL = PL_cp + PL_cc + PL_ci; 
      
      QL_cp(pvpq_buses) = b1(pvpq_buses).*QLo(pvpq_buses);
      QL_cc(pvpq_buses) = b2(pvpq_buses).*QLo(pvpq_buses).*(Vmag(pvpq_buses)./Vo(pvpq_buses));
      QL_ci(pvpq_buses) = b3(pvpq_buses).*QLo(pvpq_buses).*(Vmag(pvpq_buses)./Vo(pvpq_buses)).^2;
      QL = QL_cp + QL_cc + QL_ci; 
      
      if (Freq_bit~= 0)
       PL = PL*(1 + Kpf*delF);
       QL = QL*(1 + Kqf*delF);
      end
               
      
      Psp = sparse(zeros(nb,1));
      Psp(pv_data(:,1)) = pv_data(:,3);
      Psp(pq_data(:,1)) = Psp(pq_data(:,1))-PL(pq_data(:,1));
      Psp([Re_bus;In_bus])= Psp([Re_bus;In_bus])- Pdc([Re_bus;In_bus]);

      Qsp = sparse(zeros(nb,1));
      Qsp(pq_data(:,1)) = -QL(pq_data(:,1));
      Qsp([Re_bus;In_bus]) = Qsp([Re_bus;In_bus])- Qdc([Re_bus;In_bus]);
     end 
         
     Qsp = Qsp + (pv_num_Uvio.*QUlim)+ (pv_num_Lvio.*QLlim);  
     
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
    if (Freq_bit~= 0)
     Psp_sl  = PG_sl - PL(sl);
     delP_sl = Pc(sl) - Psp_sl;
     nosl_check = (max(abs(delP))<= tole);
     sl_check =  (abs(delP_sl)<= tole);
     del_P_sl = [delP;delP_sl];
    else
     nosl_check = (max(abs(delP))<= tole);
     sl_check = 1;
     del_P_sl = delP;
    end
    if (nosl_check~=0 & sl_check~=0)
     fprintf(fid,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
     fprintf(fid,'Converged Real power mismatch iteration No = %i, Max.mismatch = %10.8f\n',kp, full(max(abs(del_P_sl))));
     fprintf(fid,'Converged Reactive power mismatch iteration No = %i, Max. mismatch = %10.8f\n',kq, full(max(abs(delQ))));
     convergence_bit=1;
     if(sum(pvpq(:,1)==sl)==1&(Load_bit~=0|Freq_bit~=0))
      S(sl)= S(sl)+PL(sl)+j*QL(sl)+Pdc(sl)+j*Qdc(sl);
     else
      S(sl)= S(sl)+(sl_load(1,2)+j*sl_load(1,3))+Pdc(sl)+j*Qdc(sl);
     end
     
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

``` octave title="lfl_results.m"
%--------------------Preparation of Load flow Results------------------------------%
Bno=[1:nb]';
PG=sparse(zeros(nb,1));
QG=sparse(zeros(nb,1));

Pc=real(S);
Qc=imag(S);
if(Load_bit~=0 | Freq_bit~=0)
  QG_pv = Qc(pv_data(:,1))+ QL(pv_data(:,1))+ Qdc(pv_data(:,1));
else
 PL = sparse(zeros(nb,1));
 QL = sparse(zeros(nb,1));
 PL([sl;pq_data(:,1)])=[sl_load(1,2);pq_data(:,2)];
 QL([sl;pq_data(:,1)])=[sl_load(1,3);pq_data(:,3)];
 QG_pv = Qc(pv_data(:,1))+ Qload(pv_data(:,1));
end

PG([sl;pv_data(:,1)])=[Pc(sl); pv_data(:,3)];
QG([sl;pv_data(:,1)])=[Qc(sl); QG_pv];


if (nl~=0)
 com_ani = acos(2*Vdi./Vdoi-cos(Gama))-Gama;
 com_anr = acos(2*Vdr./Vdor-cos(Alpha))-Alpha;
end

fid=fopen('report.dat', 'a');
if (Q_bit~=0)
   fprintf(fid,' Q-limits at PV-buses accounted\n ' );
   fprintf(fid,'------------------------------------------------------------\n');

end
if (Load_bit~=0)
   fprintf(fid,' Voltage-dependent load models are accounted\n  ' ); 
   fprintf(fid,'------------------------------------------------------------\n');
end

if(Freq_bit~=0)
 F = 50*(1 + delF);
 fprintf(fid,' Frequency-dependent load models are accounted\n  ' ); 
 fprintf(fid,'------------------------------------------------------------\n');
 fprintf(fid,' System Nominal Frequency in Hz = %i \n',50); 
 fprintf(fid,' System New Frequency in Hz = %3.4f\n',F); 
 fprintf(fid,'-----------------------------------\n');
end
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
fprintf(fid,'----------------------------------------------------------------------------------\n');

%-------------------------------Hvdc Results-------------------------%
if (nl~=0)
 for i=1:nl
  fprintf(fid,' Results of The D.C.Link:%3i\t\n',i); 
  fprintf(fid,'----------------------------------------------------------\n');
  fprintf(fid,'\nParameter                     Rectifier      Inverter \n');
  fprintf(fid,'----------------------------------------------------------\n');
  fprintf(fid,'Bus no %30i\t%13i\n',[Re_bus(i) In_bus(i)]);
  fprintf(fid,'D.C.Voltage(pu) %27.6f\t%12.6f\n',[Vdr(i) Vdi(i)]);
  fprintf(fid,'Transformer tap Position(pu) %14.6f\t%12.6f\n',[Tr(i) Ti(i)]);
  fprintf(fid,'Control Angles(Deg) %23.6f\t%12.6f\n',[Alpha(i)*180/pi Gama(i)*180/pi]);
  fprintf(fid,'Commutation overlap Angles(Deg) %11.6f\t%12.6f\n',[com_anr(i)*180/pi com_ani(i)*180/pi]);
  fprintf(fid,'Real Power flow(pu) %23.6f\t%12.6f\n',[Pr(i) Pi(i)]);
  fprintf(fid,'Reactive power consumption(pu) %12.6f\t%12.6f\n',[Qr(i) Qi(i)]);
  fprintf(fid,'Power factor %30.6f\t%12.6f\n',[cos(Phi_r(i)) cos(Phi_i(i)) ]);
  fprintf(fid,'Current in the D.C.Link(pu) %15.6f\t\n',[Id(i)]);
  fprintf(fid,'----------------------------------------------------------\n');
  fprintf(fid,'Voltage Base in KV = %5.2f\t\n',Vb(i));
 fprintf(fid,'----------------------------------------------------------\n');
 end
 end
%---------------------------------------------------------------------------------%

if (convergence_bit==1)
 fid1=fopen('lfl.dat','w');
 fprintf(fid1,'%3i \t %12.6f \t%12.6f\t%12.6f\t%12.6f\t%9.6f\t%9.6f \n',[Bno  full(Vmag) full(Vang*180/pi) full(PG)  full(QG)  full(PL) full(QL)]');
 fclose(fid1);
 if (nl~=0) 
  fid2=fopen('hvdc_res.dat','w');
  fprintf(fid2,'%3i \t %9.6f \t%12.6f\t%12.6f\t%9.6f\t%9.6f\t%9.6f\t%5.2f\t \n',[[Re_bus;In_bus] [full(Vdr); full(Vdi)] [full(Tr);full(Ti)] [full(Pr);full(Pi)] [full(Qr);full(Qi)]  [full(Alpha*180/pi); full(Gama*180/pi)]  [cos(full(Phi_r)); cos(full(Phi_i))] [Vb(1:nl,1);Vb(1:nl,1)]]');
  fclose(fid2);
 end
else
 fprintf(fid,'Convergence is not reached !!!' );
end
fclose(fid);
```

``` octave title="fdlf_loadflow.m"
clear all;
clc
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
Load_bit=busno(8);
Freq_bit = busno(9);
npq=busno(10);                                                                                                                                                                          
nshunt=busno(11);
Vsl=busno(12);
nl = busno(13);

if (nshunt~=0)
  load shunt.dat;
end   
if (Q_bit~=0)
  load Qlim_data.dat;
end 
if (Load_bit~=0)
  load load_model.dat;
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

%-----------------Call dcdata.m to read Hvdc link input data----------------------%

if (nl~=0)
  dcdata;
else
   Re_bus = 1;
   In_bus = 2; % dummy  bus numbers
end   

%------------------------call B_bus_form.m  to construct the Bd and Bdd matrices-----%
B_bus_form;       
%--------------------------------------------------------------------------%
if (Load_bit~=0)
 PLo = sparse(zeros(nb,1));
 PLo(load_data(:,1)) = load_data(:,2);
 a1 = sparse(zeros(nb,1));
 a2 = sparse(zeros(nb,1));
 a3 = sparse(zeros(nb,1));
 a1(load_data(:,1)) = 1;
 a1(load_model(:,1)) = load_model(:,2);
 a2(load_model(:,1)) = load_model(:,3);
 a3(load_model(:,1)) = load_model(:,4);

 QLo = sparse(zeros(nb,1));
 QLo(load_data(:,1)) = load_data(:,3);
 b1 = sparse(zeros(nb,1));
 b2 = sparse(zeros(nb,1));
 b3 = sparse(zeros(nb,1));
 b1(load_data(:,1)) = 1;
 b1(load_model(:,1)) = load_model(:,5);
 b2(load_model(:,1)) = load_model(:,6);
 b3(load_model(:,1)) = load_model(:,7);
else
 Psp = sparse(zeros(nb,1));
 Psp(pv_data(:,1)) = pv_data(:,3);
 Psp(pq_data(:,1)) = Psp(pq_data(:,1))-pq_data(:,2);
 Qsp = sparse(zeros(nb,1));
 Qsp(pq_data(:,1)) = -pq_data(:,3);
 Qload =sparse(zeros(nb,1));
 Qload(pq_data(:,1))= pq_data(:,3);
end

if (Freq_bit~=0)
 load freq_model.dat;
  
 R = freq_model(1);
 Kpf = freq_model(2);
 Kqf = freq_model(3);
 PG_sl = freq_model(4);
 delF = 0;
 
 if (Load_bit==0)
  PLo = sparse(zeros(nb,1));
  QLo = sparse(zeros(nb,1));
  PLo(load_data(:,1)) = load_data(:,2);
  QLo(load_data(:,1)) = load_data(:,3);
 end
end 

pv_Vmag =sparse(zeros(nb,1));
pv_Vmag(pv_data(:,1))= pv_data(:,2);

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
disp('See file: report.dat for details' );

```