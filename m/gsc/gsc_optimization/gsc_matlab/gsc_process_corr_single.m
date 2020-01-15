close all
clc
clear
% % aj_185_1m_d2
% moments = [47.00811006267665, 53.05853704568263, 61.989149017868954, 68.16328923800685, 68.95018946214208, 74.24663327843685, 74.24663327843685, 80.53549708102877, 80.53549708102877, 86.94643588631027, 88.31370314267308, 93.42855459955203, 93.7614739251477, 99.63296021292591, 99.63296021292591, 105.52537042255975, 106.52412839934676, 111.51791828328183, 112.03242996829333, 117.3288737845881, 118.54413091604155, 123.44712462026871, 124.50641338352767, 129.71206102011453, 131.07400371573317, 136.55203989144377, 136.55203989144377, 141.90901449421048, 142.87250581150647, 149.47036153694793, 150.34805794079108, 155.2813170382542, 157.52095613771598, 163.7556271443258, 176.64466753941454, 180.39757630067484, 187.17702438553212, 193.13930685301824, 194.72739170366913, 199.47905844171643, 200.5383472049754, 206.22824113333778, 207.43885686277656, 213.52220090320657, 214.97258974385028, 221.38885310987595, 223.32583827697803, 228.8341398459246, 238.31848885027318, 243.00962480184856, 245.37032547425423, 252.66428524412305, 254.69206659093305, 261.4412492825544, 264.8179958150178, 271.5066477201672, 274.61666109016176, 280.7302705238277, 284.483179285088, 290.6573195052259, 291.942750755613, 297.9352986163351, 300.1446723225609, 305.5319123185636, 306.7727934412384, 311.4639293928137, 527.851643118506, 533.7590499176224, 536.3858011386621, 542.4086143926202, 547.3418734900832, 553.3041559575694, 553.576544496693, 559.4480307844713, 561.0304605905634, 566.8414160918696, 567.9612356416005, 572.6826369864118, 574.4077644008621, 579.3107581050892, 580.4305776548201, 584.9144975878546, 586.3672364631811, 591.6334148862398, 593.600665446578, 599.9563980261316, 604.5613791923305, 611.1785209261867, 612.9944445203449, 618.1395613704598, 619.1486157938202, 624.4044977703055, 625.0201020815985, 630.0138919655335, 631.1942423017363, 635.9459090397836, 637.2170555556943, 643.3306649893602, 644.1074140976335, 649.131469334462, 650.7052697700951, 656.3648982597755];
% aj_185_3m_d2
moments = [34.11367839929327, 40.31982116262949, 46.70184651440034, 52.87242286609165, 73.47363836627233, 78.44860301056615, 86.46550516329854, 91.8424871525858, 95.83753451845811, 101.79241644117343, 105.23620830316104, 110.56293812432412, 113.47756387552656, 118.92992411700008, 121.52063785881046, 127.02325026840816, 134.23443639422797, 138.9581401979009, 140.89340261805188, 146.3960150276496, 148.55685825698933, 153.6574533215936, 157.37611376278292, 162.90385225644272, 172.15537966992122, 177.50723557514638, 180.19572656979003, 185.07018687783548, 202.0823084731743, 208.28845123651055, 210.4209489098776, 215.2451570497989, 220.3206260303411, 226.2503818689943, 232.24515727531718, 236.91860891086594, 240.2352520070618, 245.96399917321833, 247.52181638506792, 252.3208984409271, 254.66283438039704, 260.09006853780846, 262.7283073643279, 268.2309197739256, 271.54756287012145, 277.5962005777167, 279.12889170550415, 284.45562152666724, 288.04865154754606, 293.42563353683335, 294.8666408779042, 301.85169224716526, 304.1632919808775, 311.1232172660764, 317.41581390493013, 323.82296534076306, 332.83278919147045, 338.61178852575114, 344.16465310347303, 349.49138292463607, 351.8883088496221, 357.86831685639953, 360.10453833792553, 365.3056377387781, 366.7634229186832, 371.91427015141164, 375.08015674323497, 381.1606690862607, 389.4020246586262, 395.8275536271529, 398.61654895804486, 404.52117871263596, 406.2548785129202, 411.73236483845574, 415.86403689959144, 421.7435405701205, 424.78379674163335, 429.7101902600407, 433.00170727217443, 438.22793275708915, 439.6601195486283, 444.68533636104627, 445.71550580759197, 451.86817693089347, 455.7878460445795, 461.7678540513569, 463.7025625241378, 470.7442852158406, 472.7041197726836, 478.5333712750884, 479.7142972260067, 485.5184226443494, 485.99581824152915, 490.96756341391244, 497.67622785849045, 502.6009403346601, 505.51556608586253, 511.8249402688624, 512.553596706663, 517.9557047800123];
% % aj_185_3m_d4
% moments = [4.560112072201626, 9.736085388992155, 16.721136758253163, 20.791562376311735, 20.791562376311735, 25.113248834991207, 30.888754725388615, 36.1903584624896, 39.80851456743055, 44.45684011891719, 44.45684011891719, 49.101467397082935, 51.36281496267103, 55.78500575759887, 55.78500575759887, 60.05644004815416, 60.05644004815416, 66.11880028537189, 174.2972386052191, 180.1013640235619, 180.55363353667954, 186.43313720720863, 186.7849023840779, 191.25734534712993, 194.11423882909497, 199.0640773893267, 204.35421412120476, 209.8065743626783, 210.76136555703775, 215.2086824360277, 216.26397796663548, 221.4201005642282, 227.34985640288144, 232.22431671092693, 238.94813377141608, 244.29998967664127, 255.16014326208372, 260.3612426629364, 260.8637643441782, 266.08998982909293, 276.3842039465566, 281.1581599183537, 281.66516221081173, 287.67029630165126, 287.87130497414796, 294.02719556936006, 302.58979359470214, 308.8461885261626, 309.1225754508456, 314.02216184295315, 317.99946676189256, 323.42670091930404, 323.42670091930404, 328.6026742360946, 330.66301312918597, 336.6022363131309, 339.692744652768, 345.39636573486246, 346.774806914394, 352.629184500861, 353.2070844342891, 357.930788237962, 357.930788237962, 364.061552749112, 364.061552749112, 369.0432862953421, 370.65135567531587, 376.00321158054106, 376.2042202530378, 382.9882629498021, 384.02581628035705, 390.2319590436933, 390.3827155480658, 396.11146271422234, 396.3627235548432, 402.28878112017554, 402.5902941289206, 408.2436630428909, 408.9723194806915, 414.5503101424755, 414.8266970671585, 420.55544423331503, 422.8757160649109, 428.227571970136, 428.227571970136, 433.27791486661613, 434.2075799769135, 440.41372274024974, 442.3872163759319, 447.81445053334335, 454.4226106416731, 458.4176580075454, 459.0567707445814, 463.679970212006, 463.679970212006, 469.3333391259763, 475.26309496462954, 480.0873031045508];
% % aj_185_5m_d3
% moments = [48.081830074557296, 54.6085793993898, 54.6085793993898, 60.7789375, 61.887602277222726, 67.67808880873814, 88.42845358404482, 94.17308486595026, 94.93682197279317, 101.31236651687318, 102.46356250000001, 108.50228958570645, 109.1996147702152, 115.2098937414573, 117.39239431816316, 123.83435078457734, 124.36564616325067, 130.7411907073307, 132.06033629369236, 138.10382122610156, 138.70152852710908, 144.44615980901455, 145.06798021086738, 150.5137578422691, 166.78131786001953, 171.96144780208456, 187.61118267308538, 194.6704131111915, 203.1711391699649, 209.34744794704244, 223.43487834016784, 228.88065597156955, 237.87947144784923, 244.96499046488086, 249.45475000000002, 255.73063926525754, 256.66570375245954, 262.01186350036, 263.76227329801185, 269.95577779369313, 269.95577779369313, 276.377177587383, 276.377177587383, 282.51000206693135, 282.51000206693135, 288.85352612588474, 288.85352612588474, 295.8063125, 302.7794731433355, 308.922575959246, 310.4077510928829, 317.67985658847425, 319.2548264902313, 326.2944902576531, 326.2944902576531, 332.86518031091964, 333.59571145659544, 339.2075188938326, 339.2075188938326, 347.0510431115059, 347.0510431115059, 353.4335050378276, 354.5293017563414, 360.4399628440823, 360.4399628440823, 366.45715919756594, 367.58616187724675, 373.59644084848895, 374.0854374056733, 379.79686272641175, 381.02548328959386, 386.3716430374943, 394.81860646237647, 399.60026487043655, 400.26438409377823, 405.8281608033251, 406.95716348300596, 412.3697351532406, 413.16118750000004, 419.4277801153064, 420.88884240665806, 426.401031960394, 427.6869715855883, 433.6972505568305, 444.7846806325158, 450.2968701862517, 453.1020625, 458.6315664391897, 459.2213635117557, 464.53431729848916, 465.2854375, 471.8064227940805, 471.8064227940805, 477.05296465847977, 478.14876137699355, 483.32889131905864, 496.441885025936, 501.89458003957924, 501.89458003957924, 506.941886136976];

%% parameter define
fs = 16000;
c = 340;
MicNum = 6;

%% load micArray signal
filepath = './aj_185_3m_d2/';
wavename = 'rec_mic_0.pcm';
pathInfo = dir([filepath wavename]);
SigLen = pathInfo.bytes/2;
sig_mic = zeros(MicNum, SigLen);
for m = 1 : MicNum
    wavename = ['rec_mic_' num2str(m-1) '.pcm'];
    pathInfo = dir([filepath wavename]);
    L = pathInfo.bytes/2;
    if L > SigLen
        L = SigLen;
    end
    fid = fopen([filepath wavename],'rb');
    sig_mic(m, 1:L) = fread(fid, [1,L], 'int16');
    fclose(fid);
end

%sig_mic = sig_mic(:, 60*fs + 1 : 120*fs);
%% gsc ctrl process
tic
SigLen = size(sig_mic, 2);
% FrameNum = floor((SigLen-FrameLen)/FrameInc);
hwt = waitbar(0, 'ATF process');

% ATF estimate
% atf estimate parameter define
Frame_len_atf = 512;
Frame_inc_atf = Frame_len_atf / 2;
FreLen_atf = Frame_len_atf / 2 + 1;
win_atf = hanning(Frame_len_atf, 'periodic');
alpha = 0.9;

FrameNum_atf = floor((SigLen-Frame_len_atf)/Frame_inc_atf);

moments_frame = floor((moments * fs - Frame_len_atf)/Frame_inc_atf);
moment_frame_index = 1;
moment_frame_size = size(moments_frame, 2);

% moment_length = 0;
moment_length_index = 0;
atf_moments_plot_1 = zeros(FreLen_atf, moment_frame_size/2);
atf_moments_plot_2 = zeros(FreLen_atf, moment_frame_size/2);
atf_moments_plot_3 = zeros(FreLen_atf, moment_frame_size/2);
atf_moments_plot_4 = zeros(FreLen_atf, moment_frame_size/2);
atf_moments_plot_5 = zeros(FreLen_atf, moment_frame_size/2);
atf_moments_plot_6 = zeros(FreLen_atf, moment_frame_size/2);

voice1 = enframe(sig_mic(1, :), win_atf, Frame_inc_atf);
voice2 = enframe(sig_mic(2, :), win_atf, Frame_inc_atf);
voice3 = enframe(sig_mic(3, :), win_atf, Frame_inc_atf);
voice4 = enframe(sig_mic(4, :), win_atf, Frame_inc_atf);
voice5 = enframe(sig_mic(5, :), win_atf, Frame_inc_atf);
voice6 = enframe(sig_mic(6, :), win_atf, Frame_inc_atf);

atf_length = 20;
h_mat_plot_1 = ones(FreLen_atf, FrameNum_atf);
h_mat_plot_2 = ones(FreLen_atf, FrameNum_atf);
h_mat_plot_3 = ones(FreLen_atf, FrameNum_atf);
h_mat_plot_4 = ones(FreLen_atf, FrameNum_atf);
h_mat_plot_5 = ones(FreLen_atf, FrameNum_atf);
h_mat_plot_6 = ones(FreLen_atf, FrameNum_atf);

h_mat_plot_1_corr = zeros(FreLen_atf, 1);
h_mat_plot_2_corr = zeros(FreLen_atf, 1);
h_mat_plot_3_corr = zeros(FreLen_atf, 1);
h_mat_plot_4_corr = zeros(FreLen_atf, 1);
h_mat_plot_5_corr = zeros(FreLen_atf, 1);
h_mat_plot_6_corr = zeros(FreLen_atf, 1);

corr_threshold = 0.8;
corr_smooth_factor = 0.955;
is_first_value = ones(FreLen_atf, 1);


sig_mic_reconstruct = zeros(SigLen, MicNum);
% sig_mic_reconstruct = sig_mic';

%只能存储某一刻的ATF值
h_mat = ones(FreLen_atf, MicNum);
h_mat_temp = ones(FreLen_atf, MicNum);

query_energy = zeros(MicNum, FreLen_atf);
query_energy_max = zeros(MicNum, FreLen_atf);
query_energy_id = zeros(MicNum, FreLen_atf);

coxy = zeros(MicNum, FreLen_atf);
coxx = zeros(MicNum, FreLen_atf);
coyy = zeros(MicNum, FreLen_atf);

G = zeros(FreLen_atf, MicNum - 1);
G_auxiliary = zeros(FreLen_atf, MicNum - 1);
P_est = zeros(1, FreLen_atf);

for p = atf_length : FrameNum_atf
    waitbar(p/FrameNum_atf)
    lp = (p - 1) * Frame_inc_atf + 1:(p - 1) * Frame_inc_atf + Frame_len_atf;
    if p < moments_frame(moment_frame_index)
        sig_mic_reconstruct(lp, 1) = sig_mic_reconstruct(lp, 1) + voice1(p, :)';
        sig_mic_reconstruct(lp, 2) = sig_mic_reconstruct(lp, 2) + voice2(p, :)';
        sig_mic_reconstruct(lp, 3) = sig_mic_reconstruct(lp, 3) + voice3(p, :)';
        sig_mic_reconstruct(lp, 4) = sig_mic_reconstruct(lp, 4) + voice4(p, :)';
        sig_mic_reconstruct(lp, 5) = sig_mic_reconstruct(lp, 5) + voice5(p, :)';
        sig_mic_reconstruct(lp, 6) = sig_mic_reconstruct(lp, 6) + voice6(p, :)';
        continue;
    elseif p > moments_frame(moment_frame_index + 1)
%         return;
        moment_length_index = moment_length_index + 1;
            
        atf_moments_plot_1(:, moment_length_index) = h_mat_plot_1_corr;
        atf_moments_plot_2(:, moment_length_index) = h_mat_plot_2_corr;
        atf_moments_plot_3(:, moment_length_index) = h_mat_plot_3_corr;
        atf_moments_plot_4(:, moment_length_index) = h_mat_plot_4_corr;
        atf_moments_plot_5(:, moment_length_index) = h_mat_plot_5_corr;
        atf_moments_plot_6(:, moment_length_index) = h_mat_plot_6_corr;
        
%         moment_length = 0;
        is_first_value = ones(FreLen_atf, 1);
        h_mat_plot_1_corr = ones(FreLen_atf, 1);
        h_mat_plot_2_corr = ones(FreLen_atf, 1);
        h_mat_plot_3_corr = ones(FreLen_atf, 1);
        h_mat_plot_4_corr = ones(FreLen_atf, 1);
        h_mat_plot_5_corr = ones(FreLen_atf, 1);
        h_mat_plot_6_corr = ones(FreLen_atf, 1);
        
        moment_frame_index = moment_frame_index + 2;
        if moment_frame_index > moment_frame_size
            break;
        end
        continue;
    else
%         moment_length = moment_length + 1;
%         lp = (p - 1) * Frame_inc_atf + 1:(p - 1) * Frame_inc_atf + Frame_len_atf;
        
        % 相关性的长度为2
%         corr_length = p - 1: p;
          
        coxy = zeros(MicNum, FreLen_atf);
        coxx = zeros(MicNum, FreLen_atf);
        coyy = zeros(MicNum, FreLen_atf);
          
%         for i = corr_length
%             fv1 = fft(voice1(i, :));
%             fv2 = fft(voice2(i, :));
%             fv3 = fft(voice3(i, :));
%             fv4 = fft(voice4(i, :));
%             fv5 = fft(voice5(i, :));
%             fv6 = fft(voice6(i, :));
%         
%             fv1 = fv1(1: FreLen_atf);
%             fv2 = fv2(1: FreLen_atf);
%             fv3 = fv3(1: FreLen_atf);
%             fv4 = fv4(1: FreLen_atf);
%             fv5 = fv5(1: FreLen_atf);
%             fv6 = fv6(1: FreLen_atf);
%             
%             coxx(1, :) = coxx(1, :) + (fv1 .* conj(fv1));
%             coxx(2, :) = coxx(2, :) + (fv2 .* conj(fv2));
%             coxx(3, :) = coxx(3, :) + (fv3 .* conj(fv3));
%             coxx(4, :) = coxx(4, :) + (fv4 .* conj(fv4));
%             coxx(5, :) = coxx(5, :) + (fv5 .* conj(fv5));
%             coxx(6, :) = coxx(6, :) + (fv6 .* conj(fv6));
%             
%             coyy(1, :) = coyy(1, :) + (fv1 .* conj(fv1));
%             coyy(2, :) = coyy(2, :) + (fv2 .* conj(fv1));
%             coyy(3, :) = coyy(3, :) + (fv3 .* conj(fv1));
%             coyy(4, :) = coyy(4, :) + (fv4 .* conj(fv1));
%             coyy(5, :) = coyy(5, :) + (fv5 .* conj(fv1));
%             coyy(6, :) = coyy(6, :) + (fv6 .* conj(fv1));
%         end
%           
%         coxy(1, :) = coyy(1, :)  ./ sqrt(coxx(1, :) .* coxx(1, :));
%         coxy(2, :) = coyy(2, :)  ./ sqrt(coxx(2, :) .* coxx(1, :));
%         coxy(3, :) = coyy(3, :)  ./ sqrt(coxx(3, :) .* coxx(1, :));
%         coxy(4, :) = coyy(4, :)  ./ sqrt(coxx(4, :) .* coxx(1, :));
%         coxy(5, :) = coyy(5, :)  ./ sqrt(coxx(5, :) .* coxx(1, :));
%         coxy(6, :) = coyy(6, :)  ./ sqrt(coxx(6, :) .* coxx(1, :));
        % 相关性计算结束
          
        atf_frame_length = p - atf_length + 1: p;
    
        fai1 = zeros(FreLen_atf, MicNum);
        fai3 = zeros(FreLen_atf, MicNum);
    
        out_signal = zeros(Frame_len_atf, MicNum);
    
        for i = atf_frame_length
            fv1 = fft(voice1(i, :));
            fv2 = fft(voice2(i, :));
            fv3 = fft(voice3(i, :));
            fv4 = fft(voice4(i, :));
            fv5 = fft(voice5(i, :));
            fv6 = fft(voice6(i, :));
       
            fv = [fv1;fv2;fv3;fv4;fv5;fv6];
            
            fv = fv(:, 1:FreLen_atf); 
            
            coxx(1, :) = coxx(1, :) + (fv(1, :) .* conj(fv(1, :)));
            coxx(2, :) = coxx(2, :) + (fv(2, :) .* conj(fv(2, :)));
            coxx(3, :) = coxx(3, :) + (fv(3, :) .* conj(fv(3, :)));
            coxx(4, :) = coxx(4, :) + (fv(4, :) .* conj(fv(4, :)));
            coxx(5, :) = coxx(5, :) + (fv(5, :) .* conj(fv(5, :)));
            coxx(6, :) = coxx(6, :) + (fv(6, :) .* conj(fv(6, :)));
            
            coyy(1, :) = coyy(1, :) + (fv(1, :) .* conj(fv(1, :)));
            coyy(2, :) = coyy(2, :) + (fv(2, :) .* conj(fv(1, :)));
            coyy(3, :) = coyy(3, :) + (fv(3, :) .* conj(fv(1, :)));
            coyy(4, :) = coyy(4, :) + (fv(4, :) .* conj(fv(1, :)));
            coyy(5, :) = coyy(5, :) + (fv(5, :) .* conj(fv(1, :)));
            coyy(6, :) = coyy(6, :) + (fv(6, :) .* conj(fv(1, :)));
            
            fv = fv.';
            
            for j = 1 : MicNum
                fai3(:,j) = fai3(:, j) + conj(fv(:, 1)) .* fv(:, j);
                fai1(:,j) = fai1(:, j) + conj(fv(:, 1)) .* fv(:, 1) .* conj(fv(:, 1)) .* fv(:, j);
            end
        end

        coxy(1, :) = coyy(1, :)  ./ sqrt(coxx(1, :) .* coxx(1, :));
        coxy(2, :) = coyy(2, :)  ./ sqrt(coxx(2, :) .* coxx(1, :));
        coxy(3, :) = coyy(3, :)  ./ sqrt(coxx(3, :) .* coxx(1, :));
        coxy(4, :) = coyy(4, :)  ./ sqrt(coxx(4, :) .* coxx(1, :));
        coxy(5, :) = coyy(5, :)  ./ sqrt(coxx(5, :) .* coxx(1, :));
        coxy(6, :) = coyy(6, :)  ./ sqrt(coxx(6, :) .* coxx(1, :));
        %% 计算atf 矩阵
    
        %% 一次最小二乘解
        for j = 1: MicNum
            h_mat(:, j) = (fai1(:, j) - fai3(:, 1) .* fai3(:, j)) ./ (fai1(:, 1) - fai3(:, 1) .* fai3(:, 1));
        end
        
        for i = 1: FreLen_atf
            % 选定2, 3, 4, 5, 6号mic的相关性取并集作为参考
            if abs(coxy(2, i)) > corr_threshold || abs(coxy(3, i)) > corr_threshold || abs(coxy(4, i)) > corr_threshold || abs(coxy(5, i)) > corr_threshold || abs(coxy(6, i)) > corr_threshold
                if is_first_value(i) == 1
                    is_first_value(i) = 0;
                    h_mat_plot_1_corr(i) = h_mat(i, 1);
                    h_mat_plot_2_corr(i) = h_mat(i, 2);
                    h_mat_plot_3_corr(i) = h_mat(i, 3);
                    h_mat_plot_4_corr(i) = h_mat(i, 4);
                    h_mat_plot_5_corr(i) = h_mat(i, 5);
                    h_mat_plot_6_corr(i) = h_mat(i, 6);
                else
                    h_mat_plot_1_corr(i) = corr_smooth_factor * h_mat_plot_1_corr(i) + (1 - corr_smooth_factor) * h_mat(i, 1);
                    h_mat_plot_2_corr(i) = corr_smooth_factor * h_mat_plot_2_corr(i) + (1 - corr_smooth_factor) * h_mat(i, 2);
                    h_mat_plot_3_corr(i) = corr_smooth_factor * h_mat_plot_3_corr(i) + (1 - corr_smooth_factor) * h_mat(i, 3);
                    h_mat_plot_4_corr(i) = corr_smooth_factor * h_mat_plot_4_corr(i) + (1 - corr_smooth_factor) * h_mat(i, 4);
                    h_mat_plot_5_corr(i) = corr_smooth_factor * h_mat_plot_5_corr(i) + (1 - corr_smooth_factor) * h_mat(i, 5);
                    h_mat_plot_6_corr(i) = corr_smooth_factor * h_mat_plot_6_corr(i) + (1 - corr_smooth_factor) * h_mat(i, 6);
                end 
            end
        end
        
        a = 2;
        a = a+ 1;
        
%         h_mat_plot_1(:, p) = h_mat(:, 1);
%         h_mat_plot_2(:, p) = h_mat(:, 2);
%         h_mat_plot_3(:, p) = h_mat(:, 3);
%         h_mat_plot_4(:, p) = h_mat(:, 4);
%         h_mat_plot_5(:, p) = h_mat(:, 5);
%         h_mat_plot_6(:, p) = h_mat(:, 6);
    end
end
toc
close(hwt)

moment_length_index = 0;
moment_frame_index = 1;

tic
hwt_2 = waitbar(0, 'reconstruct process');
for p = atf_length : FrameNum_atf
    waitbar(p/FrameNum_atf)
    lp = (p - 1) * Frame_inc_atf + 1:(p - 1) * Frame_inc_atf + Frame_len_atf;
    if p < moments_frame(moment_frame_index)
        continue;
    elseif p > moments_frame(moment_frame_index + 1)
        moment_length_index = moment_length_index + 1;
        moment_frame_index = moment_frame_index + 2;
        
        if moment_frame_index > moment_frame_size
            break;
        end
        continue;
    else
        out_signal = zeros(Frame_len_atf, MicNum);
        fv1 = fft(voice1(p, :));
        fv1 = fv1.';
        fv1 = fv1(1:FreLen_atf, :);

        out_signal(1:FreLen_atf, 1) = atf_moments_plot_1(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 1) = conj(out_signal(Frame_len_atf/2: -1: 2, 1));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 1)));
        sig_mic_reconstruct(lp, 1) = sig_mic_reconstruct(lp, 1) + out_signal_temp;
        
        out_signal(1:FreLen_atf, 2) = atf_moments_plot_2(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 2) = conj(out_signal(Frame_len_atf/2: -1: 2, 2));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 2)));
        sig_mic_reconstruct(lp, 2) = sig_mic_reconstruct(lp, 2) + out_signal_temp;
            
        out_signal(1:FreLen_atf, 3) = atf_moments_plot_3(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 3) = conj(out_signal(Frame_len_atf/2: -1: 2, 3));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 3)));
        sig_mic_reconstruct(lp, 3) = sig_mic_reconstruct(lp, 3) + out_signal_temp;
            
        out_signal(1:FreLen_atf, 4) = atf_moments_plot_4(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 4) = conj(out_signal(Frame_len_atf/2: -1: 2, 4));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 4)));
        sig_mic_reconstruct(lp, 4) = sig_mic_reconstruct(lp, 4) + out_signal_temp;
            
        out_signal(1:FreLen_atf, 5) = atf_moments_plot_5(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 5) = conj(out_signal(Frame_len_atf/2: -1: 2, 5));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 5)));
        sig_mic_reconstruct(lp, 5) = sig_mic_reconstruct(lp, 5) + out_signal_temp;
            
        out_signal(1:FreLen_atf, 6) = atf_moments_plot_6(:, moment_length_index + 1) .* fv1;
        out_signal(Frame_len_atf/2 + 2:Frame_len_atf, 6) = conj(out_signal(Frame_len_atf/2: -1: 2, 6));
        out_signal_temp = real(ifft(out_signal(1:Frame_len_atf, 6)));
        sig_mic_reconstruct(lp, 6) = sig_mic_reconstruct(lp, 6) + out_signal_temp;    
    end
end

toc
close(hwt_2)

%return;
%% save result for atf
result_folder = './aj_185_3m_d2/matlab/corr/';
for m = 1 : MicNum
%     wavname_ori = ['rec_mic_ori_' num2str(m) '.wav'];
    wavname_construct = ['rec_mic_construct_' num2str(m) '.wav'];
%     filetemp = [result_folder wavname_ori];
%     audiowrite(filetemp, sig_mic(m, :)'/32768, fs);
    filetemp = [result_folder wavname_construct];
    audiowrite(filetemp, sig_mic_reconstruct(:, m)/32768, fs);
end

return;