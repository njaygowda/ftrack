function [gci,MeanBasedSignal] = SEDREAMS_GCIDetection(wave,Fs,F0mean)

% USAGE:
% [gci,MeanBasedSignal] = SEDREAMS_GCIDetection(wave,Fs,F0mean)
% 
% INPUTS:
%     - wave is the speech signal
%     - Fs is the sampling frequency (Hz)
%     - F0mean is the average pitch value (Hz) for the considered speaker
%    
% OUPUTS:
%     - gci are the Glottal Closure Instant locations (in samples)
%     - MeanBasedSignal is the mean-based signal used by the SEDREAMS algorithm 
%     
% For more information, please read:
% 
% T. Drugman, M. Thomas, J. Gudnason, P. Naylor, T. Dutoit,
% "Detection of Glottal Closure Instants from Speech Signals: a Quantitative Review",
% IEEE Transactions on Audio, Speech and Language Processing, Accepted for publication.
%
% and 
%
% T.Drugman, T.Dutoit, "Glottal Closure and Opening Instant Detection from Speech Signals",
% Interspeech09, Brighton, U.K, 2009
%
% Please refer to these works in your publication if you use this code.
% 
% Code written by Thomas Drugman in TCTS Lab, University of Mons, Belgium.
% 
% Copyright (C) 2000-2011 Thomas Drugman - TCTS Lab
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Matlab version used for developing this code: 7.0
%
% http://tcts.fpms.ac.be/~drugman/      thomas.drugman@umons.ac.be


%display('GCI detection with the SEDREAMS algorithm')
%pause(0.00001)

res = GetLPCresidual(wave,round(25/1000*Fs),round(5/1000*Fs),round(Fs/1000)+2);
Res = res;

% Calculation of the mean-based signal
MeanBasedSignal=zeros(1,length(wave));
T0mean=round(Fs/F0mean);

halfL=round((1.6*T0mean)/2);
Blackwin=blackman(2*halfL+1);

for m=halfL+1:length(wave)-halfL
    vec=wave(m-halfL:m+halfL);
        
    vec=vec(:).*Blackwin(:);
    MeanBasedSignal(m)=mean(vec);
end

% Remove the low-frequency contents of the mean-based signal
Ws = 30/(Fs/2);
Wp = 50/(Fs/2);
Rp = 3; Rs = 60;
[n,Wp] = ellipord(Wp,Ws,Rp,Rs);
[b,a] = ellip(n,Rp,Rs,Wp,'high');

MeanBasedSignal=filtfilt(b,a,MeanBasedSignal);
MeanBasedSignal=MeanBasedSignal/max(abs(MeanBasedSignal));


% Detect the minima and maxima of the mean-based signal
PotMaxis=[];
PotMinis=[];

for m=2:length(MeanBasedSignal)-1    
    if (MeanBasedSignal(m)>MeanBasedSignal(m-1))&&(MeanBasedSignal(m)>MeanBasedSignal(m+1))
        PotMaxis=[PotMaxis m];
    elseif (MeanBasedSignal(m)<MeanBasedSignal(m-1))&&(MeanBasedSignal(m)<MeanBasedSignal(m+1))
        PotMinis=[PotMinis m];
    end
end

while PotMaxis(1)<PotMinis(1)
    PotMaxis(1)=[];
end

while PotMinis(end)>PotMaxis(end)
    PotMinis(end)=[];
end

Minis=PotMinis;
Maxis=PotMaxis;


% Determine the median position of GCIs within the cycle
res=res/max(abs(res));
%figure;plot(res);pause;
%Posis=find(res>0.4);
%Posis=find(abs(res)>0.4);
Posis=find(diff(diff([0 abs(res)])>0)==-1);Posis=Posis(abs(res(Posis))>0.4);
% figure;plot(res);
% figure;plot(abs(res));hold on;stem(Posis,ones(size(Posis)),'r');
% pause;
for k=1:length(Posis)

    Dists=abs(Minis-Posis(k));
    [mini,pos]=min(Dists);
    interv=Maxis(pos)-Minis(pos);

    RelPosis(k)=(Posis(k)-Minis(pos))/interv;
end
%figure;plot(find(res>0.4),RelPosis,'.');pause;
RatioGCI=median(RelPosis);

% Detect GCIs from the residual signal using the presence intervals derived
% from the mean-based signal
gci=zeros(1,length(Minis));

Ind=1;
for k=1:length(Minis)
    interv=Maxis(k)-Minis(k);    
    alpha=RatioGCI-0.25;
    start=Minis(k)+round(alpha*interv);    
    alpha=RatioGCI+0.35;
    stop=Minis(k)+round(alpha*interv);

    if start<1
        start=1;
    end
    if stop>length(res)
        stop=length(res);
    end

    vec=res(start:stop);
%     if(isempty(vec))
%         [Minis(k) Maxis(k) interv RatioGCI start stop]
%         figure;plot(res)
%         pause;
%     end
    [maxi,posi]=max(vec);
    gci(Ind)=start+posi-1;
    Ind=Ind+1;
end
