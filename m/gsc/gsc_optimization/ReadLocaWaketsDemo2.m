clear
clc

ll = 128;
fs = 16000;

% read wakeup pcm file
filepath = '/Users/daiyi02/Downloads/AEC_100db_5m_d1/aud_rec/';
% read wakeup pcm file
waveName = 'rec_mic_0_0.pcm';
pathInfo = dir([filepath waveName]);
L = pathInfo.bytes/2;
fid = fopen([filepath waveName],'rb');
sigMic = fread(fid,[L,1],'int16');
fclose(fid);

frameLen = floor(L/128);
WakFlag = zeros(frameLen, 11);
sigFlag = zeros(frameLen, 1);
wakAngle = zeros(frameLen, 1);

% read wakeup time info
filename = 'alg_wake_info_0';
fidTxtR = fopen([filepath filename],'r+','n','utf-8');
WakNum = 0;
WakNum1=0;
WakNum2=0;
WakInfo = [];
tickRepFlag = 0;
temp=0;
startFlag=0;
while 1
    content = fgetl(fidTxtR);
    contentLen = length(content);
    if contentLen < 2
        fclose(fidTxtR);
        break;
    end
    if contentLen >= 12
        if(sum(content(1:11)=='wake_status')==11)
            temp=str2double(content(13:end));
            if (temp==1) 
                WakNum=WakNum+1;
                WakInfo(WakNum, 11) = temp;
                startFlag = 1;
            elseif temp==20
                WakNum1=WakNum1+1;
            else
                WakNum2=WakNum2+1;
            end
            
        end        
    end

    if startFlag == 1
        while 1
            content = fgetl(fidTxtR);
            contentLen = length(content);
            if contentLen < 2
                fclose(fidTxtR);
                break;
            end
            if(content(1)=='-')
                break;
            end
            
            if(sum(content(1:15)=='wake_timetick_l')==15)
                WakInfo(WakNum, 1) = 65535 * tickRepFlag + str2double(content(17:end));
                if (WakNum > 1)
                    if (WakInfo(WakNum, 1) < WakInfo(WakNum-1, 1))
                        tickRepFlag = tickRepFlag + 1;
                        WakInfo(WakNum, 1) = 65535 + WakInfo(WakNum, 1);
                    end
                end                
            end

            if(contentLen >= 22 && sum(content(1:22)=='wakeword_start_time[0]')==22)
                WakInfo(WakNum, 2) = str2double(content(24:end));
                n1 = WakInfo(WakNum, 1) - WakInfo(WakNum, 2);
            end
            
            if(contentLen >= 22 && sum(content(1:22)=='wakeword_start_time[1]')==22)
                WakInfo(WakNum, 4) = str2double(content(24:end));
            end
            
            if(contentLen >= 22 && sum(content(1:22)=='wakeword_start_time[2]')==22)
                WakInfo(WakNum, 6) = str2double(content(24:end));
            end
            
            if(contentLen >= 22 && sum(content(1:22)=='wakeword_start_time[3]')==22)
                WakInfo(WakNum, 8) = str2double(content(24:end));
            end
            
            if(contentLen >= 20 && sum(content(1:20)=='wakeword_end_time[0]')==20)
                WakInfo(WakNum, 3) = str2double(content(22:end));
            end
            
            if(contentLen >= 20 && sum(content(1:20)=='wakeword_end_time[1]')==20)
                WakInfo(WakNum, 5) = str2double(content(22:end));
            end
            
            if(contentLen >= 20 && sum(content(1:20)=='wakeword_end_time[2]')==20)
                WakInfo(WakNum, 7) = str2double(content(22:end));
            end

            if(contentLen >= 20 && sum(content(1:20)=='wakeword_end_time[3]')==20)
                WakInfo(WakNum, 9) = str2double(content(22:end));
                n2 = WakInfo(WakNum, 1) - WakInfo(WakNum, 9);
            end

            if(contentLen >= 16 && sum(content(1:16)=='wakeup_direction')==16)
                WakInfo(WakNum, 10) = str2double(content(18:end));
                WakFlag(WakInfo(WakNum, 1), 1) = 1;
                WakFlag(WakInfo(WakNum, 1), 2:11) = WakInfo(WakNum, 2:11);
                nn = (n1:n2);       %????????????????
                sigFlag(nn) = 1;
                wakAngle(nn) = WakInfo(WakNum, 10);
%                 wakAngle(WakInfo(WakNum, 1)) = WakInfo(WakNum, 10);
            end
        end
        startFlag = 0;
    end
end

WakNum = size(WakInfo, 1);
disp(['WakNum = ' num2str(WakNum)]);

filename = 'Debug_WakeupTime2.txt';
fidTxtW = fopen([filepath filename],'w+','n','utf-8');
frameLenTemp = max(frameLen, size(WakFlag, 1));
for n = 1 : frameLenTemp
    for ii = 1 : 11
        fprintf(fidTxtW, '%f ', WakFlag(n,ii));
    end
    fprintf(fidTxtW, '\n');
end
fclose(fidTxtW);

filename = 'Debug_WakeupFlag.txt';
fidTxtW = fopen([filepath filename],'w+','n','utf-8');
frameLenTemp = size(sigFlag, 1);
for n = 1 : frameLenTemp
    fprintf(fidTxtW, '%d\n', sigFlag(n));
end
fclose(fidTxtW);

t1 = 0.1;
t2 = min(frameLen, size(WakFlag, 1))*8/1000;
% t1 = 20;
% t2 = 200;
L1 = floor(t1*fs/ll);
L2 = floor(t2*fs/ll);
tt = (L1:L2)*ll/fs;
sL1 = floor(t1*fs);
sL2 = floor(t2*fs);
stt = (sL1:sL2)/fs;


figure
plot(stt, sigMic(sL1:sL2)/max(abs(sigMic(sL1:sL2)))*100); 
% hold on, plot(tt, sigFlag(L1:L2) * 70,'g');
hold on, plot(tt, wakAngle(L1:L2),'r');
grid on; xlabel('time(s)'); axis([min(stt) max(stt) -180 180])