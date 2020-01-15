function [st,en]=VAD(x, fs)
x=double(x);
x=x/max(abs(x));
framelen= floor(fs*40/1000);%
frameinc= floor(fs*10/1000);%
y=enframe(x,framelen,frameinc);
%计算短时间能量 
amp=sum(abs(y),2);    
%开始端点检测
tmp1=enframe(x(1:length(x)-1),framelen,frameinc);
tmp2=enframe(x(2:length(x)),framelen,frameinc);
signs=(tmp1.*tmp2)<0;    %
diffs=(tmp1-tmp2)>0.01;
zcr=sum(signs.*diffs,2);  %    
zcr=[zcr;zcr(end)];
zcr_yu=0.2*mean(zcr);%
yuzhi=0.2*mean(amp);%
minspeech=10;
count=0;%
start=[];
tail=[];
N=length(amp);
flag=0;%


biaozhi2=0;
biaozhi3=0;
%%%%%%%%%%%%
for n=1:N 
    if amp(n)<yuzhi || (zcr(n)<zcr_yu)
        continue;
    end
    kaitou=n;
    break;
end
for n=N:-1:1
    if amp(n)<yuzhi || (zcr(n)<zcr_yu)
        continue;
    end
    jiewei=n;
    break;
end
noise=[amp(1:kaitou-1);amp(jiewei+1:end)];
noise_mean=mean(noise);
noise_var = std(noise);
speech_mean=mean(amp(kaitou:jiewei));%noise_mean
yuzhi1= 0.3*speech_mean; %
yuzhi2= max(0.3*speech_mean , (noise_mean + noise_var)*1.3);%



noise_zcr=[zcr(1:kaitou-1);zcr(jiewei+1:end)];
noise_zcr_mean = mean(noise_zcr);
noise_zcr_std = std(noise_zcr);
speech_zcr_mean=mean(zcr(kaitou:jiewei));
zcr_yu1 = 0.3*speech_zcr_mean;
zcr_yu2 = max(0.3*speech_zcr_mean , (noise_zcr_mean+noise_zcr_std)*0.3);
%%%%%%%%%%%%%%

st = [];
en = [];
bstart_state = 0;
bend_state = 0;
segment = 0;
unvoice = 0;
voice_min_len = 7;% 最短语音长度70ms  
unvoice_min_len = 5;%结束段最小持续50ms
st_candicate = 0;
en_candicate = 0;
for  i = 1:N    

    if( (amp(i) >= yuzhi2 && zcr(i) >= zcr_yu1) && ~bstart_state ) %find  start
            bstart_state = 1;  %     
            if(~st_candicate)
                st_candicate = i;
            end
            segment = segment + 1;            
            
    elseif( (amp(i) >= yuzhi2 || zcr(i) >= zcr_yu1) && bstart_state )%
        if(unvoice >= unvoice_min_len) %
            st = [st; st_candicate];
            en_candicate = en_candicate + unvoice_min_len - 1;
            en = [en; en_candicate];
       
            bstart_state = 0;
            bend_state = 1;
            unvoice = 0;
            segment = 0;
            st_candicate = 0;
        else %
            unvoice = 0;  %
            bend_state = 0;
            segment = segment + 1;       
        end
    elseif( (amp(i) < yuzhi2 && zcr(i) < zcr_yu1)  && bstart_state )     %
         if segment >= voice_min_len  %  
             unvoice = unvoice + 1; %
             if ~bend_state  %
                 en_candicate = i;
             end
             bend_state = 1; %
         else %              
             bstart_state = 0; %
             segment = 0;
             bend_state = 0;
             unvoice = 0;      
             st_candicate = 0;
             %prepare_for_start = 0;
             en_candicate = 0;
         end         
            
    elseif((amp(i) >= yuzhi2 || zcr(i) >= zcr_yu1) && ~bstart_state) %
          if ~ st_candicate
              st_candicate = i;
          end                   
    else%
        st_candicate = 0;
        continue;
    end
    
end


if(unvoice >= unvoice_min_len )%
    st = [st; st_candicate];              
    en = [en; (en_candicate + unvoice_min_len -1)];
    segment = 0;
end
if( segment >= voice_min_len) %
    st = [st; st_candicate];  
    en = [en; N];
end

figure(1);
subplot(3,1,1)
plot(x);  %原始语音波形
axis([1,length(x),-1,1])
ylabel('speech');
xlabel('样本点');

for  k=1:length(st)
line([st(k)*frameinc,st(k)*frameinc],[-1,1],'linestyle',':','color','blue','LineWidth',2);
line([en(k)*frameinc,en(k)*frameinc],[-1,1],'linestyle',':','color','red','LineWidth',2);
end

subplot(3,1,2)
plot(amp);%原始语音能量
axis([1,length(amp),0,max(amp)])
ylabel('energy');
xlabel('帧数');
line([1,N],[yuzhi1,yuzhi1],'color','yellow','LineWidth',2);%由语音能量
line([1,N],[yuzhi2,yuzhi2],'color','red','LineWidth',2);%由噪声平均能量和语音能量比较而得

for  k=1:length(st)
line([st(k),st(k)],[min(amp),max(amp)],'linestyle',':','color','blue','LineWidth',2);
line([en(k),en(k)],[min(amp),max(amp)],'linestyle',':','color','red','LineWidth',2);
end

subplot(3,1,3)
plot(zcr);%原始语音过零率
axis([1,length(zcr),0,max(zcr)])
ylabel('zcr');
xlabel('帧数');
line([1,N],[zcr_yu1,zcr_yu1],'color','yellow','LineWidth',2);%由语音能量
line([1,N],[zcr_yu2,zcr_yu2],'color','red','LineWidth',2);%语音加噪音

for  k=1:length(st)
line([st(k),st(k)],[min(zcr),max(zcr)],'linestyle',':','color','blue','LineWidth',2);
line([en(k),en(k)],[min(zcr),max(zcr)],'linestyle',':','color','red','LineWidth',2);

end

end