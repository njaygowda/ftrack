function [gci,MeanBasedSignal] = SEDREAMS_GCIDetection_ComputationalPerformanceOptimized(wave,Fs,F0mean)

% USAGE:
% [gci,MeanBasedSignal] = SEDREAMS_GCIDetection_ComputationalPerformanceOptimized(wave,Fs,F0mean)
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


display('GCI detection with the SEDREAMS algorithm')
pause(0.00001)

res = GetLPCresidual(wave,round(25/1000*Fs),round(5/1000*Fs),round(Fs/1000)+2);

%%%%%%%%%%%%%%%

% Obtain the mean-based signal
MeanBasedSignal=zeros(1,length(wave));
T0mean=round(Fs/F0mean);

halfL=round((1.7*T0mean)/2);
BlackWin=blackman(2*halfL+1);

Maxis=[];
Minis=[];
PosMeanBasedSignal=zeros(1,length(wave));
StepExp=4;
Step=2^StepExp;

for m=halfL+1:Step:length(wave)-halfL
    vec=wave(m-halfL:m+halfL);
    
    vec=vec.*BlackWin;
    MeanBasedSignal(m)=mean(vec);
    PosMeanBasedSignal(m)=1;
end

posis=find(PosMeanBasedSignal==1);
for tmp=2:length(posis)-1
    if (MeanBasedSignal(posis(tmp))>MeanBasedSignal(posis(tmp-1))) && (MeanBasedSignal(posis(tmp))>MeanBasedSignal(posis(tmp+1)))
        Maxis=[Maxis posis(tmp)];
    end
    
    if (MeanBasedSignal(posis(tmp))<MeanBasedSignal(posis(tmp-1))) && (MeanBasedSignal(posis(tmp))<MeanBasedSignal(posis(tmp+1)))
        Minis=[Minis posis(tmp)];
    end
end


for StepExp=3:-1:0
    Step=2^StepExp;
    
    for tmp=1:length(Maxis)
        m=Maxis(tmp)-Step;
        vec=wave(m-halfL:m+halfL);
        vec=vec.*BlackWin;
        MeanBasedSignal(m)=mean(vec);
        PosMeanBasedSignal(m)=1;
        
        m=Maxis(tmp)+Step;
        vec=wave(m-halfL:m+halfL);
        vec=vec.*BlackWin;
        MeanBasedSignal(m)=mean(vec);
        PosMeanBasedSignal(m)=1;
        
        VecTmp=[MeanBasedSignal(Maxis(tmp)-Step) MeanBasedSignal(Maxis(tmp)) MeanBasedSignal(Maxis(tmp)+Step)];
        [ManTmp,PosTmp]=max(VecTmp);
        Maxis(tmp)=Maxis(tmp)+(PosTmp-2)*Step;
    end
    
    for tmp=1:length(Minis)
        m=Minis(tmp)-Step;
        vec=wave(m-halfL:m+halfL);
        vec=vec.*BlackWin;
        MeanBasedSignal(m)=mean(vec);
        PosMeanBasedSignal(m)=1;
        
        m=Minis(tmp)+Step;
        vec=wave(m-halfL:m+halfL);
        vec=vec.*BlackWin;
        MeanBasedSignal(m)=mean(vec);
        PosMeanBasedSignal(m)=1;
        
        VecTmp=[MeanBasedSignal(Minis(tmp)-Step) MeanBasedSignal(Minis(tmp)) MeanBasedSignal(Minis(tmp)+Step)];
        [MinTmp,PosTmp]=min(VecTmp);
        Minis(tmp)=Minis(tmp)+(PosTmp-2)*Step;
    end
    
end


MeanBasedSignal=MeanBasedSignal/max(abs(MeanBasedSignal));

if Maxis(1)<Minis(1)
    Maxis(1)=[];
end

if Minis(end)>Maxis(end)
    Minis(end)=[];
end


% Determine the median position of GCIs within the cycle
res=res/max(abs(res));
Posis=find(res>0.4);
for k=1:length(Posis)

    Dists=abs(Minis-Posis(k));
    [mini,pos]=min(Dists);
    interv=Maxis(pos)-Minis(pos);

    RelPosis(k)=(Posis(k)-Minis(pos))/interv;
end

RatioGCI=median(RelPosis);

% Detect GCIs from the residual wave using the presence intervals derived
% from the mean-based wave
gci=[];

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
    [maxi,posi]=max(vec);
    gci=[gci start+posi-1];
end

