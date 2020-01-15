clear all;
clc;

%% parameter init
fs = 16000;
c = 340;

MicNum = 6;
%PMic = zeros(3, MicNum);
%PMic(1, :) = 0.026 * [0:MicNum-1];
% figure,plot(PMic(1,:), PMic(2,:), 'o', 'LineWidth', 1.5); grid on
PMic=[[0.0298 0 0]' [0.0149 -0.0258 0]' [-0.0149 -0.02598 0]' [-0.0298 0 0]' [-0.0149 0.0258 0]' [0.0149 0.0258 0]'];
%ststic analyse

FrameLen = 512;
FrameInc = FrameLen/2;
freDis = fs/FrameLen;
freMin = 100;
freMax = 4000;
freBand = freDis;
FreMin = floor(freMin/fs*FrameLen);
FreMax = floor(freMax/fs*FrameLen);
FreDis = floor(freBand/fs*FrameLen);
FreUse = FreMin:FreDis:FreMax;
FreLen = length(FreUse);
hammWin = hamming(FrameLen);

thetaDis = 2;
thetaDisp = (0:thetaDis:360);
theta = thetaDisp / 180 * pi;
thetaLen = length(theta);

steeringA = zeros(FreLen, MicNum, thetaLen);
sourceR = 1.0;
sourcePos = zeros(3,1);
sourcePos(3) = sourceR * sin(20/180*pi);
for n = 1 : thetaLen
    sourcePos(1) = sourceR * cos(20/180*pi) * cos(theta(n));
    sourcePos(2) = sourceR * cos(20/180*pi) * sin(theta(n));
    for m = 1 : MicNum
        sourceMicDis = sqrt(sum((sourcePos - PMic(:,m)).^2));
        diffDis = sourceR - sourceMicDis;
        fre = FreUse'/FrameLen*fs;
        steeringA(:, m, n) = exp(j*2*pi*fre*diffDis/c);
    end
end

filepath = 'C:\data\??so????\songruitao 1-25\';
filename = 'alg_dir_info_0';
InfoLoca = load([filepath filename]);

%% load micArray signal
% read micArray signal
%filepath = 'D:\Project\20180614_toXiaotan\LineArray\huawei\Bach2_20180604\aud_rec\';
micNo = 0:5;
waveName = 'rec_mic_0_0.pcm';
pathInfo = dir([filepath waveName]);
L = pathInfo.bytes/2 + 4*fs;
sigOri = zeros(L, MicNum);
for m = 1 : MicNum
    waveName = ['rec_mic_' num2str(micNo(m)) '_0.pcm'];
    pathInfo = dir([filepath waveName]);
    L = pathInfo.bytes/2;
    fid = fopen([filepath waveName],'rb');
    sigOri(1:L,m) = fread(fid,[L,1],'int16');
    fclose(fid);
end
SigLen = size(sigOri,1);
% SigLen = 10*fs;
% FrameNum = floor((SigLen-FrameLen)/FrameInc);
TimeStart = 43;
TimeEnd = 44;
FrameStart = floor(TimeStart*fs/FrameInc);
FrameEnd = floor(TimeEnd*fs/FrameInc);
FrameNum = FrameEnd - FrameStart + 1;
frameDisp = (FrameStart:FrameEnd)*FrameInc/fs;


%static analyse
L1 = 1;%floor(TimeStart*fs/FrameInc);             %begin frame number
L2 = FrameNum;%floor(TimeEnd*fs/FrameInc);             %max frame number
tt = (L1:L2)*FrameInc/fs;                   % 8ms step
sL1 = floor(TimeStart*fs);
sL2 = floor(TimeEnd*fs);
stt = (sL1:sL2)/fs;                   % each sample point time
tt = (L1:L2)*FrameInc/fs+TimeStart;

fstart=128*FrameStart;
fend=128*FrameEnd;

ll=128;
SL1 = floor(TimeStart*fs/ll);             %begin frame number
SL2 = floor(TimeEnd*fs/ll);             %max frame number
tt1 = (SL1:SL2)*ll/fs;                   % 8ms step





%% localization process
tN = 30;
select_num=80;
dataFFT = zeros(FrameLen, MicNum, tN);

tic
fbfEnergy = zeros(FrameNum, thetaLen);
mvdrEnergy = zeros(FrameNum, thetaLen);
fbfEnergyNorm = zeros(FrameNum, thetaLen);
mvdrEnergyNorm = zeros(FrameNum, thetaLen);
fbfEnergyNorm2 = zeros(FrameNum, thetaLen);
mvdrEnergyNorm2 = zeros(FrameNum, thetaLen);

fbfselect_Energy = zeros(FrameNum, thetaLen);
musicD = zeros(FrameNum, 4);



fbfTemp_select = zeros(select_num, 1);

hwt = waitbar(0, 'localization process');
for p = 1 : FrameNum
    waitbar(p/FrameNum, hwt);
    data = sigOri((p-1+FrameStart)*FrameInc + 1 : (p-1+FrameStart)*FrameInc + FrameLen, :);
    for m = 1 : MicNum
        data(:,m) = data(:,m).*hammWin;
    end
    dataFFT(:, :, tN) = fft(data);
    
    if  p >= tN 
        for ff=1:FreLen
            mfre=FreUse(ff);
            coxy=sum(dataFFT(mfre,1,:).*conj(dataFFT(mfre,4,:)));
            coxx=sum(dataFFT(mfre,1,:).*conj(dataFFT(mfre,1,:)));
            coyy=sum(dataFFT(mfre,4,:).*conj(dataFFT(mfre,4,:)));
            cofre=(coxy*coxy)/(coxx*coyy);
            cofre1(ff)=abs(cofre);
        end
            [aa,bb]=sort(cofre1,'descend');
            flag=1;
    end
        
    if flag==1
         for n = 1 : thetaLen
            %mvdrTemp = zeros(FreLen, 1);
            X = zeros(MicNum, tN);
             for fn = 1:FreLen %: select_num                                    %????????????????
                    %mfre = FreUse(bb(fn));
                    %mf= mfre/FrameLen*fs;
                    %A=exp(i*2*pi*mf*(0:MicNum-1)*0.026*cos(theta(n))/c);
                    A=steeringA(fn,:,n);
                    for m = 1 : tN
                        X(:, m) = dataFFT(mfre, :, m)';
                    end
                    Y = X ./ abs(X);
                    fbfTemp_select(fn) = abs(A * Y * Y' * A');          %??????????
                    
                    for m = 1 : tN
                        X(:, m) = dataFFT(mfre, :, m)';
                    end
                    %music ????
%                     Rxx = X * X'/tN;       %????????????????????
%                     [EV,D] = eig(Rxx);
%                     EVA=diag(D)';
%                     [EVA,I]=sort(EVA);
%                     EVA=fliplr(EVA);
%                     EV=fliplr(EV(:,I));
%                     En=EV(:,2+1:4);
%                     mvdrTemp(fn) = abs(1/ (A * En* En' * A'));
                    %mvdr ????
                     Corr = X * X' + 0.0000001*diag(ones(MicNum,1));
                     mvdrTemp(fn) = 1 / abs(A * inv(Corr) * A');
                    
                    
             end
              fbfselect_energy(p, n) = mean(fbfTemp_select.^2);
              mvdrEnergy(p, n) = mean(mvdrTemp.^2); 
         end
         mvdrEnergyNorm(p, :) = mvdrEnergy(p, :) / max(mvdrEnergy(p, :));
         mvdrEnergyNorm2(p, :) = mvdrEnergyNorm(p, :) >=1;
         
         fbfselect_EnergyNorm(p, :) = fbfselect_energy(p, :) / max(fbfselect_energy(p, :));
         fbfselect_EnergyNorm2(p, :) = fbfselect_EnergyNorm(p, :) >=1;
    end
    
    for m = 1 : tN-1
        dataFFT(:, :, m) = dataFFT(:, :, m+1);
    end
    if(p==1000)
        p=p;
    end
end
close(hwt)
toc

figure
subplot(121)
imagesc(thetaDisp, frameDisp, mvdrEnergy)
xlabel('angle'); ylabel('time')
subplot(122)
imagesc(thetaDisp, frameDisp, fbfselect_EnergyNorm)
xlabel('angle'); ylabel('time')

figure
subplot(121)
imagesc(thetaDisp, frameDisp, mvdrEnergyNorm)
xlabel('angle'); ylabel('time')
subplot(122)
imagesc(thetaDisp, frameDisp, fbfselect_energy)
xlabel('angle'); ylabel('time')

figure
subplot(121)
imagesc(thetaDisp, frameDisp, mvdrEnergyNorm2)
xlabel('angle'); ylabel('time')
subplot(122)
imagesc(thetaDisp, frameDisp, fbfselect_EnergyNorm2)
xlabel('angle'); ylabel('time')

%angel analyse
histNum=zeros(FrameNum,thetaLen);
locationAngel=zeros(FrameNum,1);

for n = 1 : FrameNum
       
            amount=sum(fbfselect_EnergyNorm2(n,:));
            [value,location]=sort(fbfselect_EnergyNorm2(n,:),'descend');
            for i=1:amount
                tick=location(i);
                histNum(n,tick) = histNum(n,tick) + 1;
                tick1 = tick - 1;
                if tick1 < 1
                    tick1 = tick1 + thetaLen;
                    histNum(n,tick1) = histNum(n,tick1) + 1;
                else
                    histNum(n,tick1) = histNum(n,tick1) + 1;
                end
                tick2 = tick + 1;
                if tick2 > thetaLen
                    tick2 = tick2 - thetaLen;
                    histNum(n,tick2) = histNum(n,tick2) + 1;
                else 
                    histNum(n,tick2) = histNum(n,tick2) + 1;
                end    
            end
      
        [maxa,maxb]=max(histNum(n,:));
        if maxa>0
            locationAngel(n,1)=maxb*2;    %??????????
        end
end

   for i=L1:L2
       if locationAngel(i,1)>180 
            locationAngel(i,1)=locationAngel(i,1)-360;
       end
   end

 figure
   subplot(211)
   plot(stt, sigOri(sL1:sL2,1)/max(abs(sigOri(sL1:sL2,1)))*100);
   hold on, plot(tt, locationAngel(L1:L2, 1),'.c');         % 1st LocaDetectPos
   grid on
   locationAngel1=locationAngel;
   histNum=zeros(FrameNum,thetaLen);
   locationAngel=zeros(FrameNum,1);

  
for n = 1 : FrameNum
       
            amount=sum(mvdrEnergyNorm2(n,:));
            [value,location]=sort(mvdrEnergyNorm2(n,:),'descend');
            for i=1:amount
                tick=location(i);
                histNum(n,tick) = histNum(n,tick) + 1;
                tick1 = tick - 1;
                if tick1 < 1
                    tick1 = tick1 + thetaLen;
                    histNum(n,tick1) = histNum(n,tick1) + 1;
                else
                    histNum(n,tick1) = histNum(n,tick1) + 1;
                end
                tick2 = tick + 1;
                if tick2 > thetaLen
                    tick2 = tick2 - thetaLen;
                    histNum(n,tick2) = histNum(n,tick2) + 1;
                else 
                    histNum(n,tick2) = histNum(n,tick2) + 1;
                end    
            end
      
        [maxa,maxb]=max(histNum(n,:));
        if maxa>0
            locationAngel(n,1)=maxb*2;    %??????????
        end
end

   for i=L1:L2
       if locationAngel(i,1)>180 
            locationAngel(i,1)=locationAngel(i,1)-360;
       end
   end
  
   hold on, plot(tt, locationAngel(L1:L2, 1),'.g');         % 1st LocaDetectPos
   for i=L1:L2
       if locationAngel1(i,1)>180 
            locationAngel1(i,1)=locationAngel1(i,1)-360;
       end
   end
   
   subplot(212)
   plot(stt, sigOri(sL1:sL2,1)/max(abs(sigOri(sL1:sL2,1)))*100);
   hold on, plot(tt, locationAngel1(L1:L2, 1),'.r');         % 1st LocaDetectPos
   
   hold on, plot(tt, InfoLoca(SL1:2:SL2+1, 1),'.c');         % 1st LocaDetectPos
  % hold on, plot(tt, -InfoLoca(SL1:2:SL2, 2),'.c');         % 2nd LocaDetectPos
   grid on
    axis([min(stt) max(stt) -200 200])
