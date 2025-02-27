``` octave title="y_bus_from.m"
% -------------------------------------------------------------%
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
```


``` octave title="powerflow.m"
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

``` octave title="loadflow.m"
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
%------------------------call ybus_form.m  to construct the YBUS-----------%
ybus_form;       
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
%-----------------call jacob_form.m to obtain the bus voltages----------%

jacob_form

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

``` octave title="jacob_form.m"
%~~~~~~~~~~~~~~~~~~~~Form Jacobian and solve for bus voltages and angles~~~~~~~~~~%
G=real(Y);
B=imag(Y);
theta_pq = sparse(zeros(nb,nb));
VpVq = sparse(zeros(nb,nb));
Vmag=sparse(ones(nb,1));
Vang=sparse(zeros(nb,1));
Vmag(pv_data(:,1))=pv_data(:,2);  %Replaces the voltages of P-V buses in Vmag 

Vmag(sl) = Vsl;
Vang(sl) = 0;
fid=fopen('report.dat', 'w');
fprintf(fid,'                   Detailed Report of Load flow (NR method)\n' ); 
fprintf(fid,'                   ----------------------------\n' );
%------------------Iteration begins from here onwards---------------------%
for k=1:10
   Vbus=Vmag.*(cos(Vang)+ j*sin(Vang));
%----------------------------------bus powers calculations------------------------%
   S=Vbus.*(conj(Y*Vbus));
   Pc=real(S);
   Qc=imag(S);
 %---------------------------------finding non slack and non pv buses------------------%
   pv_num=pv_data(:,1);
   pv_sl_num=[pv_num;sl];
   num_no_sl_pv =[1:nb]'; 
   num_no_sl_pv(pv_sl_num,:)=[];
   
   %-----------Perform the following if Q-limit is accounted--------------%
   start_cri=2;
   if (Q_bit~=0 &k>start_cri)
     Qc(pvpq_buses) = Qc(pvpq_buses) + Qload(pvpq_buses); 
     
 %-----------------------checking Qlimit voilation-----------------------------%
     pv_num_Uvio=sparse(zeros(nb,1));
     pv_num_Lvio=sparse(zeros(nb,1));
      
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
     
 %----------------------------------Re-calculation of bus powers --------------------%
     Vbus=Vmag.*(cos(Vang)+ j*sin(Vang));
     S=Vbus.*(conj(Y*Vbus));
     Pc=real(S);
     Qc=imag(S);
   end       %-------Q-limit part ends here---------------
   
  %------------------------bus power mis matches -------------%
   delP=Psp-Pc;
   delP(sl,:)=[];
   delQ=Qsp-Qc;
   delQ(pv_sl_num,:)=[];
   del_PQ=[delP;delQ]; 
   
   if (max(abs(del_PQ))<= tole)
     convergence_bit=1;
     S(sl)= S(sl)+(sl_load(1,2)+j*sl_load(1,3));
     fprintf(fid,'__________________________________________________________________________\n');
     fprintf(fid,'Converged loadflow iteration No = %i, Max. mismatch = %10.8f\n',k-1, full(max(abs(del_PQ))));
     fprintf(fid,'----------------------------------\n'); 
     break;
   else
      del_Pm=sparse(zeros(nb,1)); 
      del_Qm=sparse(zeros(nb,1));
      del_Pm(sort([pv_num;num_no_sl_pv]))=delP;
      del_Qm(num_no_sl_pv)=delQ;
      Bno=[1:nb]';
      
      if (max(abs(del_Pm))>=max(abs(del_Qm)))
         vio_busP = Bno(abs(del_Pm)==max(abs(del_Pm)));
         if (k~=1)
           fprintf(fid,'Iter. No = %i Max. Real power mismatch at bus = %i, Max. mismatch = %10.8f\n',(k-1),full(vio_busP),full(max(abs(del_Pm))));
         end  
      else
         vio_busP=Bno(abs(del_Qm)==max(abs(del_Qm)));
         if (k~=1)
           fprintf(fid,'Iter. No = %i Max. Reactive power mismatch at bus = %i, Max. mismatch = %10.8f\n',(k-1),full(vio_busP),full(max(abs(del_Qm))));
         end    
      end
   end
      
 %--------------------------------calculation of Jacobian elements-----------%
   for i=1:nline+ntrans
     theta_pq(nt(i,1),nt(i,2))=Vang(nt(i,1))- Vang(nt(i,2));
     theta_pq(nt(i,2),nt(i,1))=Vang(nt(i,2))- Vang(nt(i,1));
     VpVq(nt(i,1),nt(i,2)) = Vmag(nt(i,1))* Vmag(nt(i,2));
     VpVq(nt(i,2),nt(i,1)) = VpVq(nt(i,1),nt(i,2));
   end
   
%---------------------------------------off diagonal------------------------%
   Hpq=VpVq.*(G.*sin(theta_pq)-B.*cos(theta_pq));
   Npq=VpVq.*(G.*cos(theta_pq)+B.*sin(theta_pq));
   
%---------------------------diagonal---------------------------------------%
   Hpp=-Qc-(diag(B).*(Vmag.*Vmag));
   Lpp=-Hpp-2*diag(B).*(Vmag.*Vmag);
   Npp=Pc+diag(G).*(Vmag.*Vmag);
   Mpp=Npp-2*diag(G).*(Vmag.*Vmag);
%----------------------------------H matrix--------------------------------%
   H=Hpq+sparse(diag(Hpp));
%---------------------------------L martrix--------------------------------%
   L=Hpq+sparse(diag(Lpp));
%---------------------------------N matrix---------------------------------%
   N=Npq+sparse(diag(Npp));
%---------------------------------M matrix---------------------------------%
   M=(-Npq)+sparse(diag(Mpp));
%---------------------------formation of jacobian matrix-------------------%
   H(sl,:)=[];
   H(:,sl)=[];
   L(:,pv_sl_num)=[];
   L(pv_sl_num,:)=[];
   N(:,pv_sl_num)=[];
   N(sl,:)=[];
   M(pv_sl_num,:)=[];
   M(:,sl)=[];
   J=[H N;M L];
%----------------------Calculation of correction vectors----------------%
   delx=J\del_PQ;
%----------------------Updation of bus angles at each bus----------------% 
   Vang(num_nosl)=Vang(num_nosl)+ delx(1:nb-1);
   
%----------------------Updation of bus voltages at each bus----------------% 
   
   del_Vmag=delx(nb:size(delx,1)).*Vmag(num_no_sl_pv);
   Vmag(num_no_sl_pv)=Vmag(num_no_sl_pv) + del_Vmag;
end
fclose(fid);
```

``` octave title="lfl_results.m"


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
if (Q_bit~=0)
   fprintf(fid,' Q-limits at PV-buses accounted\n ' );
   fprintf(fid,'------------------------------------------------------------\n');

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

if (convergence_bit==1)
 fid1=fopen('lfl.dat','w');
 fprintf(fid1,'%3i \t %12.6f \t%12.6f\t%12.6f\t%12.6f\t%9.6f\t%9.6f \n',[Bno  full(Vmag) full(Vang*180/pi) full(PG)  full(QG)  full(PL) full(QL)]');
 fclose(fid1);
else
 fprintf(fid,'Convergence is not reached !!!' );
end
fclose(fid);
```
