clear all 
%close all
ct=500;
snrdb = [20 ]; %dbm
Mx = [2 :2 : 10]; 
 
al=3;
alphax = [0.1 0.01 0.001 0.0001];
sigma2_dbm= -94;%+10*log10(BW)+Nf; %Thermal noise in dBm -90+10+10
sigman=10^((sigma2_dbm-30)/10);
disc=3;

R0=0.5;%the target data rate is 2 bits/s/Hz   
eps0 = 2^R0-1;

P0 = 10^((snrdb-30)/10);    
for aict = 1: length(alphax)%length(Mx)%
    M = 8;%Mx(aict); 
    N=16; K=N;
    tauvec = 0.01*ones(K,1);%tau =0.01
    alpha1 = alphax(aict);
    sum1 = 0;
    for ict = 1 : ct

     %BD users, backward and forward links      
     locm = sign(rand(M,2)-0.5) .* disc.*rand(M,2);% - [disc*ones(M,1) zeros(M,1)]; %location  of BD [-1 1; -1.5 2; -4 -1; -2 3];%
     locd = sign(rand(K,2)-0.5) .* disc.*rand(K,2);% + [disc*ones(K,1) zeros(K,1)]; %downlink users' locations  [3 0; 3 1; 3 -1]; 
     dist = max(1, sqrt(sum(locm.^2,2))); %distance to the base station
     Hmk_hfading = complex(sqrt(0.5)*randn(K,M),sqrt(0.5)*randn(K,M));% N subcarriers with M users, forward channel
     Hmk =   Hmk_hfading  * diag(1./dist.^al); % forward channel for BD-BS K by M
     Fmk_hfading = complex(sqrt(0.5)*randn(K,M),sqrt(0.5)*randn(K,M));% N subcarriers with M users, backward channel
     Fmk =   Fmk_hfading  * diag(1./dist.^al); % backward channel for BD-BS K by M

     %BD-DU channels
     for m =1 : M
        for k = 1 : K
            distUD(m,k) = max(1,sqrt(sum((locd(k,:)-locm(m,:)).^2)));
        end
     end
     gmk_fading = complex(sqrt(0.5)*randn(M,K),sqrt(0.5)*randn(M,K));
     Gmk = gmk_fading ./ (distUD.^al); %composite channels for BD-U0, M by K

     %downlink users 
     gfading = complex(sqrt(0.5)*randn(K,K),sqrt(0.5)*randn(K,K));% K by K, each column is a user's channel
     distd = max(1, sqrt(sum(locd.^2,2))); %distance to the base station
     Gk = gfading * diag(1./distd.^al); % channel for BS-downlink

     %self interference 
     hsi =  complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1)) ;
     hsi2 = abs(hsi)^2;

     % signal
     x0=   complex(sqrt(0.5)*randn(K,1),sqrt(0.5)*randn(K,1));
     Dx = diag(x0);

     %define the linear constraints     
     A1 = [abs(Gmk.').^2.*abs(Hmk).^2 ];   
     a1 = tauvec;%P0*tauk.'-sigman*ones(K,1); %the first constraint A1*y <= a1;  
     A2 = [-eye(M)];
     a2 = zeros(M,1); % the second constraint y>=0 
     A3 = [eye(M) ];
     a3 = [ones(M,1)]; % ym<=1 
     Aall = [A1; A2; A3]; 
     aall = [a1; a2; a3];
     
     x00=rand(M,1); %random initial
     %find barH
     breveH = 1/sqrt(alpha1*P0*hsi2+sigman)* Hmk'.*Fmk';%M by K 
     barH = sqrt(P0)*Dx*breveH';

     options = optimoptions('fmincon','Display', 'off','MaxFunctionEvaluations', 300000); %display off
     x0 = fmincon(@(x) 1,x00,Aall, aall,[],[],[],[],[],options);%a feasible starting point    
     if min(x0)<=0 | max(Aall*x0 - aall)>0 %for some reasons, a negative x0 is very damaging
         f=0; f2=0; fr=0; f2r=0;
     else
         %Sum capacity solution         
         y = fmincon(@(y) 100*myFSMAC(y,barH,M,K),x0,Aall, aall,[],[],[],[],[],options);  
 
        %Sum capacity  
        f2 = real(log2(det(eye(K)+barH*diag(y)*barH')));%real(log2(det(eye(N)+y(1)*hms0(1)*H(:,1)*H(:,1)' + y(2)*hms0(2)*H(:,2)*H(:,2)') ));
 
        %random Sum capacity
        f2r = real(log2(det(eye(K)+barH*diag(x0)*barH')));%real(log2(det(eye(N)+x0(1)*hms0(1)*H(:,1)*H(:,1)' + x0(2)*hms0(2)*H(:,2)*H(:,2)') ));%capacity random

  
     end
     
     %OMA
     rateomax2(ict)=0;%sum capacity
     for m = 1 : M
         etaomamx = min(tauvec./A1(:,m));;%K by 1
         etaomam = max(0,min(1,etaomamx));
         %if etaomam<=1 & etaomam>0 
             rateomax2(ict) = rateomax2(ict) + real(log2(1+etaomam*barH(:,m)'*barH(:,m)));%sum capacity
         %end
     end
     
    ratex2(ict) = f2;%capacity
    ratex4(ict) = f2r; %capacity random
    end
    rate2(aict) = mean(ratex2);%capacity
    rate4(aict) = mean(ratex4);%capacity random
    rateoma2(aict) = mean(rateomax2)/M; %oma based on sum capacity
end
%plot( Mx, rateoma2,Mx,rate3,Mx,rate1, Mx,rate4,Mx,rate2)
semilogx( alphax, rateoma2/K,  alphax,rate4/K,alphax,rate2/K)
  
 
 %define the objective function for mac
 function [f] = myFSMAC(y,barH,M,K) 
      sumf=0;
      for m =1 : M
          sumf = sumf + y(m)*barH(:,m)*barH(:,m)';
      end
      f = -real(log2(det(eye(K)+sumf)));%- real(log2(det(eye(N)+y(1)*hms0(1)*H(:,1)*H(:,1)' + y(2)*hms0(2)*H(:,2)*H(:,2)') ));
 end
 
 
 
 