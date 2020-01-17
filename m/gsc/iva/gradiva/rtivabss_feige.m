function [y, W] = rtivabss_feige(x, nfft, eta, nsou,theta)
%    [y, W] = ivabss(x, nfft, maxiter, tol, eta, nsou)
%     y : separated signals (nsou x N)
%     W : unmixing matrices (nsou x nmic x nfft/2+1)
%     x : observation signals (nmic x N),
%           where nsou is # of sources, nmic is # of mics, and N is # of time frames
%     nfft : # of fft points (default =1024)
%     eta : learning rate (default =0.1)
%     nsou : # of sources (default =nmic)

[nmic, nn] = size(x);
if ~exist('nfft','var')|isempty(nfft), nfft = 1024; end
if ~exist('eta','var')|isempty(eta), eta = 0.1; end
if ~exist('nsou','var')|isempty(nsou), nsou = nmic; end

if(nmic~= nsou)
    nmic = nsou;
end
win = 2*hanning(nfft,'periodic')/nfft;
nol = fix(2*nfft/4);
for l=1:nmic,
    X(l,:,:) = conj(stft(x(l,:)', nfft, win, nol)');
end
X1 = X;

fs = 16000;
c = 342;
d = 0.08;
target = theta*pi/180;
LL = length(target);
MicNum = nsou;
H = zeros(MicNum,LL,nfft/2+1);
delta = d/(MicNum-1);

PMic = zeros(3, MicNum);
for i = 1:MicNum
    PMic(1, i) = -d/2+delta * (i-1);
end
freMin = floor(0 / fs * nfft)+1;

for f = 1 : nfft/2+1
    Steering = zeros(MicNum,LL);
    for ll = 1:LL
        fre = (f + freMin-1) / nfft * fs;
        omiga = 2 * pi * fre;
        d0 = zeros(1,MicNum);
        for i = 1:MicNum
            d0(i) = exp(-1j * omiga*(i-1)*delta* cos(target(ll))/c);
        end 
        Steering(:,ll) =  d0.';     
        
    end
    H(:,:,f) = Steering;

end



% clear x;

N = size(X,2);
nfreq = size(X,3);

epsi = 1e-10;
pObj = Inf;

% Meomory allocations
Wp = zeros(nsou,nsou,nfreq);
W = zeros(nsou,nsou,nfreq);
dWp = zeros(size(Wp));
Q = zeros(nsou,nsou,nfreq);
Xp = zeros(nsou,N,nfreq);
S = zeros(nsou,N,nfreq);
Ssq = zeros(nsou,N);
Ssq1 = zeros(nsou,N);
Ssq2 = zeros(nsou,N);
norm_err = zeros(N,nfreq);
eta_k=zeros(N,nfreq);

% Execute PCA and initialize
dmatrix= zeros(1,N);    
    
   
for k=1:nfreq,   
    Wp(:,:,k) = eye(nsou);
end
Xp= X;  
norm_X =  mean(abs(Xp).^2,1);
gama = zeros(nsou,nsou,nfreq);
N
for n=1:1:N
    n
     for k=1:nfreq,
         if (n ==1)
            norm_X(1,n,k) = 0.5*0 + 0.5*norm_X(1,n,k);
        else
            norm_X(1,n,k) = 0.5*norm_X(1,n-1,k) + 0.5*norm_X(1,n,k);
        end
     end
    
    for iter = 1:1
    dlw = 0;
    for k=1:nfreq,
        S(:,n,k) = Wp(:,:,k)*Xp(:,n,k);
    end
    Ssq(:,n) = sum(abs(S(:,n,:)).^2,3).^(1/2);
    
     Ssq1(:,n) = (Ssq(:,n)+epsi).^-1;

    for k=1:nfreq,
        % Calculate multivariate score function and gradients
        
        Phi = Ssq1(:,n).*S(:,n,k);
        
        tmp = Phi*S(:,n,k)';
        gama(:,:,k) = tmp;%0.5* gama(:,:,k) + 0.5*
        err = (eye(nsou).*gama(:,:,k) - gama(:,:,k));
         
         Phi1 = 2*Phi - Phi.^3;% Ssq2(:,n).*S(:,n,k).*S(:,n,k).*S(:,n,k);
         norm_ww1 = 2*err*Phi1*Xp(:,n,k)';
        
        dWp(:,:,k) = err*Wp(:,:,k);
        n11 = real(trace(norm_ww1*dWp(:,:,k)'));

        u = norm(err,'fro')^2/(2*n11+epsi);
       
        dWp(:,:,k) = dWp(:,:,k) * eta/ sqrt(norm_X(1,n,k));
        %%gss
        Egc = Wp(:,:,k)*H(:,:,k) - eye(nsou);%%diag(diag(Wp(:,:,k)*H(:,:,k)));%%
        Jgc = Egc*H(:,:,k)';
        ugc = norm(Egc,'fro')/(2*norm(Jgc,'fro'));
        dgss = ugc * Jgc;
        
        
        
        dWp(:,:,k) = dWp(:,:,k) - dgss;
        aaa = dWp(:,:,k);

    end
    
    % Update unmixing matrices
    Wp = Wp + dWp;
    end
    for k=1:nfreq,
         W(:,:,k) = Wp(:,:,k);%*Q(:,:,k);
        for isrc = 1:nsou
            norm_w = norm(W(isrc,:,k));
            if(norm_w  < epsi)
                norm_w = norm_w + epsi;
            end
            
            W(isrc,:,k) = W(isrc,:,k)/(norm_w); 
        end
%         W(:,:,k) = diag(diag(pinv(W(:,:,k))))*W(:,:,k);
        S(:,n,k) = W(:,:,k)*X(:,n,k);
    end
%     Wp(:,:,:) = W(:,:,:);
end


% Re-synthesize the obtained source signals
for k=1:nsou,
    
    y(k,:) = istft(conj(squeeze(S(k,:,:))'), nn, win, nol)';
end