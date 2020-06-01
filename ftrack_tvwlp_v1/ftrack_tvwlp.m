function [Fi,Ak] = ftrack_tvwlp(s,fs,lptype,nwin,nshift,p,q,npeaks,PREEMP,fint,PLOT_FLAG)

%%%%%%%%%%%%%%%%%%%%%%
% Usage: [Fi,Ak] = ftrack_tvwlp(s,fs,lptype,nwin,nshift,p,q,npeaks,PREEMP,fint,PLOT_FLAG)
%
%   Computes the continuous time-varying formant tracks at every sample
% shift. Uses the GLOAT gci detection tool kit (http://tcts.fpms.ac.be/~drugman/Toolbox/GLOAT.zip).
%
% INPUTS:
%   s       - speech signal (will be always converted to 8kHz)
%   fs      - sampling rate 
%   lptype  - one of {'tvwlp_l2','tvwlp_l1','tvlp_l2','tvlp_l1'} [default 'tvwlp_l2']
%   nwin    - window or block size for TVLP analysis in samples (@ given fs)
%                [default 1600] @ 8kHz (200ms)
%   nshift  - window shift for TVLP [default 1600 or 200ms]
%   p       - LP or TVLP order [default 8]
%   q       - polynomail order [default 3]
%   npeaks  - number of formants to track [default 3]
%   PREEMP  - preemphasis factor [default 0.97]
%   fint    - formant interval in samples fs/fint gives the output formant rate  [default 80] (10 ms)
%   PLOT_FLAG - plot flag [default 0]
%
% OUTPUTS:
%   Fi      - Formant tracks [npeaks x Ns] where Ns is signal length
%   Ak      - Time-varying LPCs [p+1 x Ns]
%
% Eg: [s,fs]=audioread('test-dr8-mjln0-sx279.wav');
%     [Fi,Ak] = ftrack_tvwlp(s,fs); 
%     OR
%     [Fi,Ak] = ftrack_tvwlp(s,fs,'tvwlp_l2',3200,3200,8,3,3,0.5,1);
%
% Author: D.Gowda, 24 Oct, 2016
%%%%%%%%%%%%%%%%%%%%%
addpath ./GLOAT/

    fs_ref=8000;
    n1ms=floor(fs_ref/1000);

    if(~exist('lptype','var'))
        lptype='tvwlp_l2';
    end
    if(~exist('nwin','var'))
        nwin=200*n1ms;
    end
    if(~exist('nshift','var'))
        nshift=200*n1ms;
    end
    if(~exist('p','var'))
        p=8;
    end      
    if(~exist('q','var'))
        q=3;
    end
    if(~exist('npeaks','var'))
        npeaks=3;
    end      
    if(~exist('PREEMP','var'))
        PREEMP=0.97;
    end  
    if(~exist('PLOT_FLAG','var'))
        PLOT_FLAG = 0;
    end
    
    if(fs ~= fs_ref)
        s=resample(s,fs_ref,fs);
        nwin=floor(nwin*fs_ref/fs);
        nshift=floor(nshift*fs_ref/fs);
        fs=fs_ref;
    end

    if(exist('PREEMP','var'))
        shat=filter([1 -PREEMP],1,s);
    else
        shat=s;
    end
    Ns=length(shat);

    switch lower(lptype)
        case {'tvwlp_l2','tvwlp_l1'}
            %%% QCP weight function computation
            %%% GCI detection using SEDREAMS
            wmin = 0.00001; % Minimum value of weighting function
            DQ = 0.7; % Duration Quotient (relative to fundamental period)
            PQ = 0.05; % Position Quotient (relative to fundamental period)
            Nramp = 3; % Length of linear ramp (in samples)

            f0min=60;
            f0max=600;
            % Pitch tracking using a method based on the Summation of the Residual Harmonics
            [f0,VUV] = SRH_PitchTracking(s,fs,f0min,f0max);
            f0_tmp=f0.*VUV;
            pos= (f0_tmp~=0);
            f0_tmp=f0_tmp(pos);
            f0mean=mean(f0_tmp);
            % Oscillating Moment-based Polarity Detection
            %pol = OMPD_PolarityDetection(s,fs,f0mean,VUV);
            % Speech Event Detection using the Residual Excitation And a Mean-based Signal
            %[gc] = SEDREAMS_GCIDetection(pol*s,fs,f0mean);
            [gc] = SEDREAMS_GCIDetection(s,fs,f0mean);
            w = qcp_wt(s,p,DQ,PQ,wmin,Nramp,gc,fs);
        otherwise
            w=ones(size(s)); % just dummy as of now no use
    end

    Fi=zeros(npeaks,floor((nwin-nshift)/2));
    Ak=zeros(p+1,floor((nwin-nshift)/2));
    for j=0:nshift:Ns-nwin;

        if(j<=Ns-nwin-nshift)            
            x=shat(j+[1:nwin]);%.*hamming(nwin);
            wj=w(j+[1:nwin]);
        else % last window
            x=shat(j+1:end);%take all samples
            wj=w(j+1:end);
        end

        switch lower(lptype)
            case {'tvlp_l2'}
                [aki]=tvlp_l2(x,p,q);
            case {'tvwlp_l2'}
                [aki]=tvwlp_l2(x,p,q,wj);
            case {'tvlp_l1'}
                [aki]=tvlp_l1(x,p,q);
            case {'tvwlp_l1'}
                [aki]=tvwlp_l1(x,p,q,wj);
            otherwise
                disp('Unknown method.')
        end

        Nx=length(x);
        [fi,ak]=tvlptoformants_akitofi(aki,Nx,npeaks,fs);

        if(j==0)
            Fi=fi(1:floor((nwin-nshift)/2),1:npeaks);
        end
        
        Fi = [Fi;fi(floor((nwin-nshift)/2)+[1:nshift],1:npeaks)];
        Ak = [Ak ak(:,floor((nwin-nshift)/2)+[1:nshift])];
        
        if(j>Ns-nwin-nshift)
            Fi=[Fi;fi(nshift+floor((nwin-nshift)/2)+1:end,1:npeaks)];
        end
    end
    
    Fi(end+1:Ns,:)=0;
    Ak(:,end+1:Ns)=0;   

    Fi=Fi(fint:fint:end,:)';
    Fi=medfilt1(Fi,5,[],2);
    
    if(PLOT_FLAG)
        n1ms=floor(fs/1000);
        [S,f,t]=spectrogram(s,20*n1ms,15*n1ms,2^10,fs);
        figure;imagesc(t,f/1000,20*log10(abs(S)));axis xy;colorbar;
        xlabel('Time (s)');ylabel('Freq (kHz)');
        hold on;plot((fint:fint:Ns)/fs,Fi'/1000,'.k');
    end
    
return;
    
function [ w ] = qcp_wt( x, p, DQ, PQ, d, Nramp, gci_ins, fs )
%   Create a AME weight function for frame x for LPC order p
%   DQ = duration quotient (from 0 to 1)
%   PQ = Position Quotient (from 0 to 1)
%   d = minimum value of the weight function
%   Nramp = length of the linear ramp (in samples)
%   gci_ins = Glottal Closure Instants of the frame x

N = length(x);
%Nramp = 1;
if Nramp > 0
UPramp = linspace(d,1,2+Nramp);
UPramp = UPramp(2:end-1);
DOWNramp = UPramp(end:-1:1);
end

if DQ+PQ > 1
    DQ = 1-PQ;
end

w = d.*ones(1,N+p);

for i = 1:length(gci_ins)-1
   T = gci_ins(i+1)-gci_ins(i);
   T1 = round(DQ*T);
   T2 = round(PQ*T);
   while T1+T2 > T
       T1 = T1-1;
   end
   w(gci_ins(i)+T2:gci_ins(i)+T2+T1-1) = 1;
   if Nramp > 0
       w(gci_ins(i)+T2:gci_ins(i)+T2+Nramp-1) = UPramp;
       if gci_ins(i)+T2+T1-Nramp > 0
           w(gci_ins(i)+T2+T1-Nramp:gci_ins(i)+T2+T1-1) = DOWNramp;
       end
   end
end

Nend = N-(T2+gci_ins(i+1));

if T2+gci_ins(i+1) < N
    if T1+T2 < Nend
        w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+T1-1) = 1;
        if Nramp > 0
            w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+Nramp-1) = UPramp;
            w(gci_ins(i+1)+T2+T1-Nramp:gci_ins(i+1)+T2+T1-1) = DOWNramp;
        end
    else
        T1 = Nend-T2;
                w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+T1-1) = 1;
        if Nramp > 0
            w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+Nramp-1) = UPramp;
        end
    end
end

return;
