clear all
close all
speed_sound=343; %m/s

Mic=[5920,2580,2140;
     5920,2180,2140;
     4200,0,2140;
     3800,0,2140;
     2200,0,2140;
     1800,0,2140;
     0,625,2140;
     0,1025,2140;
     0,3750,2140;
     0,4150,2140;
     1800,4750,2140;
     2200,4750,2140;
     3800,4750,2140;
     4200,4750,2140]; %in mm
 
Nch=length(Mic);
pairs=[1 2;3 4;5 6;7 8;9 10;11 12;13 14]; %%7 microphone pairs
Npairs=length(pairs);

NominalSource=[1350 3650]; %Nominal source position

%%Load raw audio files
Fs=44100;
fin=fopen('speech_14ch.raw','r');
data=fread(fin,Inf,'short');
fclose(fin);
s=reshape(data,Nch,length(data)/Nch);


%%STFT Parameters
Winlen=4096;
Overlap=2;
Nfft=4096;
Nframes=floor(length(s)/(Winlen/Overlap)-1);
S=zeros(Nch,Nfft);

TauMax=52; %%Maximum time delay for a microphone distance of 40 cm at
fr=1;

for k = 1:Winlen/Overlap:length(s)-Winlen
    
    for m=1:Nch
        z=s(m,k:k+Winlen-1).*blackman(Winlen).';
        S(m,:)=fft(z);
    end
    
    for p=1:Npairs
        
        CC = S(pairs(p,1),:).*conj(S(pairs(p,2),:));
        cc=real(ifft(CC./abs(CC))).';        
        
        gcc(fr,p,:)=[flipud(cc(1:TauMax+1,:));flipud(cc(end-TauMax+1:end,:))];        
        
        [~,im]=max(gcc(fr,p,:));
        TDOA(fr,p)=im-TauMax-1;
        Peak(fr,p)=gcc(fr,p,im);
        
    end
    fr=fr+1;
end

%%For each microphone pair:
%% --> Compare the behaviour of the GCC-PHAT for each microphone pair
%% --> Compare the TDOA estimaton with the reference TDOA


%%Acoustic Map%%
xstep=100;
ystep=100;
Max_X = 5920; %Room boundaries
Max_Y = 4750;

%Define grid
x=1:xstep:Max_X;
y=1:ystep:Max_Y;

%Pre-compute delays for each point of the search grid
for ix=1:length(x)
    for iy=1:length(y)
        for p=1:Np
            %delay(ix,iy,p)= %% Compute the integer time difference of arrival
        end
    end
end
   

for fr=1:size(gcc,1)    
    map=zeros(length(x),length(y));
    for ix=1:length(x)
        for iy=1:length(y)
            map(ix,iy)=0;
            for p=1:Np
                map(ix,iy)=map(ix,iy)+gcc(fr,p,delay(ix,iy,p));%delay(ix,iy,p)= %% Compute the time difference of arrival
            end
        end
    end
    imagesc(map)
    %%Estimate the source position by maximizing map
end

    






