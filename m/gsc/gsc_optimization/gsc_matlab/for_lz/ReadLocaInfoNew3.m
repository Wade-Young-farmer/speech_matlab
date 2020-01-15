% clear
% clc

ll = 128;
fs = 16000;
refectFlag = 0;

filepath = './aj_185_1m_d2/matlab/construct/';
% waveName = 'rec_recog_0.pcm';
% waveName = 'opt_outsig_asr.pcm';
waveName = 'wak_asr_0.pcm';
% waveName = 'rec_mic_0_0.pcm';
% waveName = 'test_asr_a.raw';
pathInfo = dir([filepath waveName]);
L = pathInfo.bytes/2;
fid = fopen([filepath waveName],'rb');
sigMic = fread(fid,[L,1],'int16');
fclose(fid);

% waveName = 'opt_outsig_asr.pcm';
pathInfo = dir([filepath waveName]);
L = pathInfo.bytes/2;
fid = fopen([filepath waveName],'rb');
sigMicOri = fread(fid,[L,1],'int16');
fclose(fid);

L = min(length(sigMic), length(sigMicOri));

% filename = 'alg_dir_info_0';
% filename = 'opt_locainfo_track.txt';
filename = 'wak_locainfo_0.txt';
% filename = 'test_float_a.txt';
InfoLoca = load([filepath filename]);

t1 = 0.1;
t2 = min(length(InfoLoca(:,1))*8/1000, L/fs) - 1;
% t2 = floor(t2/10)*10;
%t1 = 121;
%t2 = 180;
L1 = floor(t1*fs/ll);
L2 = floor(t2*fs/ll);
tt = (L1:L2)*ll/fs;
sL1 = floor(t1*fs);
sL2 = floor(t2*fs);
stt = (sL1:sL2)/fs;

wakNum = 0;
correct_loc = -90;
correct_loc_num = 0;
DisturbInfo = zeros(10,12);
for n = L1 : L2-1
    if (InfoLoca(n, 7) == 0 && InfoLoca(n+1, 7) > 0 && InfoLoca(n+1, 7) < 10)
        wakNum = wakNum + 1;
        DisturbInfo(wakNum, :) = InfoLoca(n+1, :);
%         DisturbInfo(wakNum, 7) = (n+1) * 0.008;
        if (InfoLoca(n+1, 4) > correct_loc - 30 && InfoLoca(n+1, 4) < correct_loc + 30)
            correct_loc_num = correct_loc_num + 1;
        end
    end
    if InfoLoca(n, 7) > 10
        InfoLoca(n, 7) = -1;
    end
end
disp (['wakNum:' num2str(wakNum) ' locaNum:' num2str(correct_loc_num)]);

ratio1 = floor(InfoLoca(L1:L2, 9)) / 100;
ratio2 = InfoLoca(L1:L2, 9) - floor(InfoLoca(L1:L2, 9));

if refectFlag == 0

figure
subplot(211);
plot(stt, sigMic(sL1:sL2)/max(abs(sigMic(sL1:sL2)))*100); grid on
hold on, plot(tt, InfoLoca(L1:L2, 1),'.c');         % 1st LocaDetectPos
hold on, plot(tt, InfoLoca(L1:L2, 2),'.g');         % 2nd LocaDetectPos
hold on, plot(tt, InfoLoca(L1:L2, 4),'.r');         % LocaWakeup
hold on, plot(tt, InfoLoca(L1:L2, 7)*20,'m','LineWidth',2);        % wakeup flag
hold on, plot(tt, InfoLoca(L1:L2, 8),'.k');         % disturb pos
% hold on, plot(tt, (InfoLoca(L1:L2, 8).*(InfoLoca(L1:L2, 10)>0.2)),'.k');         % disturb pos
% hold on, plot(tt, InfoLoca(L1:L2, 10),'.g');        % wakeup energy
hold on, plot(tt, InfoLoca(L1:L2, 5)*20,'c');        % DoubleTalkState
hold on, plot(tt, ratio1*100,'r');        % disturb ratio 1
hold on, plot(tt, ratio2*100,'g');        % disturb ratio 2
% hold on, plot(tt, InfoLoca(L1:L2, 11)*50,'k');        % asr flag
% hold on, plot(tt, InfoLoca(L1:L2, 11),'.m');        % 2nd disturb pos
% hold on, plot(tt, InfoLoca(L1:L2, 6)*20,'m');        % disturb vad result
grid on; xlabel('time(s) 阵列输出'); axis([min(stt) max(stt) -180 180])
subplot(212);
plot(stt, sigMicOri(sL1:sL2)/max(abs(sigMicOri(sL1:sL2)))*100); grid on
hold on, plot(tt, InfoLoca(L1:L2, 4),'.r');         % LocaWakeup
hold on, plot(tt, InfoLoca(L1:L2, 7)*40,'m','LineWidth',2);        % wakeup flag
% hold on, plot(tt, InfoLoca(L1:L2, 10),'.g');        % history wakeup loca
hold on, plot(tt, InfoLoca(L1:L2, 10)*100,'g');        % wakeup energy
hold on, plot(tt, InfoLoca(L1:L2, 11),'.k');        % second disturb loca
% hold on, plot(tt, InfoLoca(L1:L2, 12),'c');        % disturb energy
grid on; xlabel('time(s) 原始信号'); axis([min(stt) max(stt) -180 180])

else

figure
subplot(211);
plot(stt, sigMic(sL1:sL2)/max(abs(sigMic(sL1:sL2)))*100); grid on
hold on, plot(tt, InfoLoca(L1:L2, 1),'.c');         % 1st LocaDetectPos
hold on, plot(tt, InfoLoca(L1:L2, 2),'.g');         % 2nd LocaDetectPos
hold on, plot(tt, InfoLoca(L1:L2, 4),'.r');         % LocaWakeup
hold on, plot(tt, InfoLoca(L1:L2, 7)*20,'m','LineWidth',2);        % wakeup flag
hold on, plot(tt, InfoLoca(L1:L2, 8),'.k');         % disturb pos
hold on, plot(tt, InfoLoca(L1:L2, 5)*20,'c');        % DoubleTalkState
hold on, plot(tt, ratio1*100,'r');        % disturb ratio 1
hold on, plot(tt, ratio2*100,'g');        % disturb ratio 2
grid on; xlabel('time(s) 阵列输出'); axis([min(stt) max(stt) -180 180])
subplot(212);
plot(stt, sigMicOri(sL1:sL2)/max(abs(sigMicOri(sL1:sL2)))*100); grid on
hold on, plot(tt, InfoLoca(L1:L2, 4),'.r');         % LocaWakeup
hold on, plot(tt, InfoLoca(L1:L2, 7)*40,'m','LineWidth',2);        % wakeup flag
hold on, plot(tt, InfoLoca(L1:L2, 10),'.g');        % wakloca_reflect1
hold on, plot(tt, InfoLoca(L1:L2, 11),'.m');        % wakloca_reflect2
hold on, plot(tt, InfoLoca(L1:L2, 12),'.c');        % Disturb_MostlikelyPos
grid on; xlabel('time(s) 原始信号'); axis([min(stt) max(stt) -180 180])

end