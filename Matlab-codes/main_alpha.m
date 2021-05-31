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
 

P0 = 10^((snrdb-30)/10);    
for aict = 1: length(alphax)%length(Mx)%
    M = 10;%Mx(aict); number of BDs
    N=M; K=M; %numbers of antennas and downlink users
    tauvec = 0.01*ones(K,1);%tau =0.01
    alpha1 = alphax(aict);
    sum1 = 0;
    for ict = 1 : ct

     %BD users
     hfading = complex(sqrt(0.5)*randn(N,M),sqrt(0.5)*randn(N,M));% N by M, each column is a user's channel
     locm = sign(rand(M,2)-0.5) .* disc.*rand(M,2);% - [disc*ones(M,1) zeros(M,1)]; %location  of BD [-1 1; -1.5 2; -4 -1; -2 3];%
     locd = sign(rand(K,2)-0.5) .* disc.*rand(K,2);% + [disc*ones(K,1) zeros(K,1)]; %downlink users' locations  [3 0; 3 1; 3 -1]; 
     dist = max(1, sqrt(sum(locm.^2,2))); %distance to the base station
     Hm =   hfading  * diag(1./dist.^al); % channel for BD-BS

     %BD-DU channels
     for m =1 : M
        for k = 1 : K
            distUD(m,k) = max(1,sqrt(sum((locd(k,:)-locm(m,:)).^2)));
        end
     end
     gmk_fading = complex(sqrt(0.5)*randn(M,K),sqrt(0.5)*randn(M,K));
     Gmk = gmk_fading ./ (distUD.^al); %composite channels for BD-U0

     %downlink users
     %Ax = complex(sqrt(0.5)*randn(max(N,K),max(N,K)),sqrt(0.5)*randn(max(N,K),max(N,K)));%this matrix is used to get orthogonal beams
     %[U,V,d] = svd(Ax);
     gfading = complex(sqrt(0.5)*randn(N,K),sqrt(0.5)*randn(N,K));%U(1:N,1:K);%complex(sqrt(0.5)*randn(N,K),sqrt(0.5)*randn(N,K));% N by K, each column is a user's channel
     distd = max(1, sqrt(sum(locd.^2,2))); %distance to the base station
     Gk = gfading * diag(1./distd.^al); % channel for BS-downlink

     %self interference 
     Hsi =  complex(sqrt(0.5)*randn(N,N),sqrt(0.5)*randn(N,N)) ;
     Csi = Hsi*Hsi';

     % signal
     W = Gk./sqrt(sum(abs(Gk).^2))/sqrt(K); %each column, wk, is a normalized version of each user's channel vector
     x0=   complex(sqrt(0.5)*randn(K,1),sqrt(0.5)*randn(K,1));
     s0 =   W * x0; 
     hms0 = abs(Hm.'*s0).^2;

     %QR for BD channel matrix
     [Q, R] = qr(sqrtm(inv(sigman*eye(N)+P0*alpha1*Csi))*Hm);
     hmw = sum(abs(Hm.'*W).^2,2);%M rows for M BD 

     %define the linear constraints     
     A1 = [abs(Gmk.').^2*diag(hmw) ]; %[P0*abs(Gmk.').^2*diag(hmw) ]; 
     a1 = tauvec;%P0*tauk.'-sigman*ones(K,1); %the first constraint A1*y <= a1;  
     A2 = [-eye(M)];
     a2 = zeros(M,1); % the second constraint y>=0 
     A3 = [eye(M) ];
     a3 = [ones(M,1)]; % ym<=P0*hml 
     Aall = [A1; A2; A3]; 
     aall = [a1; a2; a3];
     
     x00=rand(M,1); %random initial
     H = sqrtm(inv(sigman*eye(N)+P0*alpha1*Csi))*Hm;

     options = optimoptions('fmincon','Display', 'off','MaxFunctionEvaluations', 300000); %display off
     x0 = fmincon(@(x) 1,x00,Aall, aall,[],[],[],[],[],options);%a feasible starting point    
     if min(x0)<=0 | max(Aall*x0 - aall)>0 %for some reasons, a negative x0 is very damaging
         f=0; f2=0; fr=0; f2r=0;
     else
         % QR solution
         x = fmincon(@(x) 100*myFS(x,R,sigman,M,hmw,P0),x0,Aall, aall,[],[],[],[],[],options);%x = fmincon(@(x) myFS(x,R,sigman,M),x0,Aall, aall,[],[],[],[],[],options);  

         %Sum capacity solution         
         y = fmincon(@(y) 100*myFSMAC(y,H,hms0,N,M),x0,Aall, aall,[],[],[],[],[],options);  

        %QR sum rate
        rr = diag(R(1:M,1:M)).^2; %get those rii
        sumtemp= -exp(1./rr./hmw./x/P0).*ei( -1./rr./hmw./x/P0 );
          for m = 1 :M
              vm = 1./rr(m)./hmw(m)./x(m)/P0;
              if vm>=700
                  sumtemp(m) = 1/vm;%-exp(vm)*ei(-vm);%some need approximation due to the resolution of Matlab
              end
          end
          f = log2(exp(1))*sum(sumtemp);
        %Sum capacity  
        f2 = real(log2(det(eye(N)+H*diag(y.*hms0)*H')));%real(log2(det(eye(N)+y(1)*hms0(1)*H(:,1)*H(:,1)' + y(2)*hms0(2)*H(:,2)*H(:,2)') ));

        %random QR
        sumtemp= -exp(1./rr./hmw./x0/P0).*ei( -1./rr./hmw./x0/P0 );
          for m = 1 :M
              vm = 1./rr(m)./hmw(m)./x0(m)/P0;
              if vm>=700
                  sumtemp(m) = 1/vm;%-exp(vm)*ei(-vm);%some need approximation due to the resolution of Matlab
              end
          end
        fr = log2(exp(1))*sum(sumtemp);
        %random Sum capacity
        f2r = real(log2(det(eye(N)+H*diag(x0.*hms0)*H')));%real(log2(det(eye(N)+x0(1)*hms0(1)*H(:,1)*H(:,1)' + x0(2)*hms0(2)*H(:,2)*H(:,2)') ));%capacity random

     end
     
     %OMA
     rateomax2(ict)=0;%sum capacity
     for m = 1 : M
         etaomamx = min(tauvec./hmw(m)./abs(Gmk(m,:).').^2);;%min((tauk.'-sigman/P0)./hmw(m)./abs(Gmk(m,:).').^2);
         etaomam = max(0,min(1,etaomamx));
         %if etaomam<=1 & etaomam>0 
             rateomax2(ict) = rateomax2(ict) + real(log2(1+etaomam*hms0(m)*H(:,m)'*H(:,m)));%sum capacity
         %end
     end
    
     if f<fr
         dfd=0;
     end
    ratex1(ict) = f;%qr
    ratex2(ict) = f2;%capacity
    ratex3(ict) = fr;%qr random
    ratex4(ict) = f2r; %capacity random
    end
    rate1(aict) = mean(ratex1);%qr
    rate2(aict) = mean(ratex2);%capacity
    rate3(aict) = mean(ratex3);%qr random
    rate4(aict) = mean(ratex4);%capacity random
    %rateoma(aict) = mean(rateomax)/M; %oma based on QR
    rateoma2(aict) = mean(rateomax2)/M; %oma based on sum capacity
end
%plot( Mx, rateoma2,Mx,rate3,Mx,rate1, Mx,rate4,Mx,rate2)
semilogx( alphax, rateoma2,alphax,rate3,alphax,rate1, alphax,rate4,alphax,rate2)
 
 %define the objective function for QR
 function [f] = myFS(x,R,sigman,M,hmw,P0) 
      rr = diag(R(1:M,1:M)).^2; %get those rii
      sumtemp= -exp(1./rr./hmw./x/P0).*ei( -1./rr./hmw./x/P0 );
      for m = 1 :M
          vm = 1./rr(m)./hmw(m)./x(m)/P0;
          if vm>=700
              sumtemp(m) = 1/vm;%-exp(vm)*ei(-vm);%some need approximation due to the resolution of Matlab
          end
      end
      f = -sum(sumtemp); 
 end
 
 
 %define the objective function for QR
 function [f] = myFSMAC(y,H,hms0,N,M) 
      sumf=0;
      for m =1 : M
          sumf = sumf + y(m)*hms0(m)*H(:,m)*H(:,m)';
      end
      f = -real(log2(det(eye(N)+sumf)));%- real(log2(det(eye(N)+y(1)*hms0(1)*H(:,1)*H(:,1)' + y(2)*hms0(2)*H(:,2)*H(:,2)') ));
 end
 
 
 
 