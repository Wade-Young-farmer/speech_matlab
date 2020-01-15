clear
clc

fs = 16000;
ChNum = 8;
MicNum = 6;

filepath = './aj_185_3m_d4/';
waveName = 'test_a.raw';
pathInfo = dir([filepath waveName]);
L = pathInfo.bytes / 2 / ChNum;
fid = fopen([filepath waveName], 'rb');
sigOri = fread(fid, [ChNum, L], 'int16');
fclose(fid);

sigMic = zeros(MicNum, L);
sigMic(1,:) = sigOri(1,:);
sigMic(2,:) = sigOri(2,:);
sigMic(3,:) = sigOri(5,:);
sigMic(4,:) = sigOri(6,:);
sigMic(5,:) = sigOri(7,:);
sigMic(6,:) = sigOri(8,:);
sigRef = sigOri(3,:);
% sigRef = zeros(1, L);

% sumTemp = mean(abs(sigOri(4, end-fs:end)));
% if(sumTemp > 200)
%     disp('录音通道错误：3、4通道有数据')
%     return;
% end


% waveName = 'test_asr_wakeup_a.raw';
% pathInfo = dir([filepath waveName]);
% L2 = pathInfo.bytes / 2 / 2;
% fid = fopen([filepath waveName], 'rb');
% sigOut = fread(fid, [2, L2], 'int16');
% fclose(fid);
% 
% waveName = 'test_asr_a.pcm';
% fid = fopen([filepath waveName], 'wb');
% fwrite(fid, sigOut(1,:), 'int16');
% fclose(fid);
% 
% waveName = 'test_wakeup_a.pcm';
% fid = fopen([filepath waveName], 'wb');
% fwrite(fid, sigOut(2,:), 'int16');
% fclose(fid);



% frameLen = 128;
% frameNum = floor(L / frameLen);
% MicEnergy = zeros(MicNum, frameNum);
% MicRatio = zeros(MicNum, frameNum) + 1;
% for p = 2 : frameNum
%     for m = 1 : MicNum
%         MicEnergy(m,p) = 0.7*MicEnergy(m,p-1) + 0.3*mean(sigMic(m,(p-1)*frameLen+1:p*frameLen).^2);
%         MicRatio(m,p) = MicEnergy(1,p) / MicEnergy(m,p);
%     end
% end
% % fs = 16000;
% % t1 = 70;
% % t2 = 71;
% % sL1 = floor(t1*fs/frameLen);
% % sL2 = floor(t2*fs/frameLen);
% % tt = (sL1:sL2)*frameLen/fs;
% % figure,plot(tt,MicRatio(:,sL1:sL2)');grid on
% figure,plot(MicRatio');grid on
% 
% for p = 1 : frameNum
%     for m = 1 : MicNum
%         sigMic(m,(p-1)*frameLen+1:p*frameLen) = sigMic(m,(p-1)*frameLen+1:p*frameLen) * MicRatio(m,p);
%     end
% end
% 
for m = 1 : MicNum
    waveName = ['rec_mic_' num2str(m-1) '.pcm'];
    fid = fopen([filepath waveName], 'wb');
    fwrite(fid, sigMic(m,:), 'int16');
    fclose(fid);
end

waveName = 'rec_spk_l_0.pcm';
fid = fopen([filepath waveName], 'wb');
fwrite(fid, sigRef, 'int16');
fclose(fid);