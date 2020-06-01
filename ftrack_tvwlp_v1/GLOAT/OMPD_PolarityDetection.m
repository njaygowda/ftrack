function [Polarity] = OMPD_PolarityDetection(wave,Fs,F0mean,VUVDecisions)

% USAGE:
% [Polarity] = OMPD_PolarityDetection(wave,Fs,F0mean,VUVDecisions)
% 
% INPUTS:
%     - wave is the speech signal
%     - Fs is the sampling frequency (Hz)
%     - F0mean is the average pitch value (Hz) for the considered speaker
%     - VUVDecisions is the vector of voiced-unvoiced VUVDecisions (with an hopsize of 10ms)
%    
% OUPUT:
%     - Polarity is the speech polarity for the input signal
%     
% For more information, please read:
% 
% T.Drugman, T.Dutoit, "Oscillating Statistical Moments for Speech Polarity Detection",
% Non-Linear Speech Processing Workshop (NOLISP11), Las Palmas, Gran Canaria, Spain, 2011
%  
% Please refer to this work in your publication if you use this code.
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


display('Polarity Detection with OMPD Method')
pause(0.0001)

if Fs>16000
    wave=resample(wave,16000,Fs);
    Fs=16000;
end

T0mean=round(Fs/F0mean);

% Calculation of the oscillating moments
EvenMoment=zeros(1,length(wave));
OddMoment=zeros(1,length(wave));

halfL=round((1.5*T0mean)/2);
halfL2=round((2.5*T0mean)/2);

BlackWin=blackman(2*halfL+1);
BlackWin2=blackman(2*halfL2+1);

wave2=wave.^2;

Step=round(Fs/4000); %These moments are actually downsampled to 4kHz, in order to reduce computation time

for k=halfL2+1:Step:length(wave)-halfL2
    vec=wave(k-halfL:k+halfL);
    vec=vec.*BlackWin;
    vec2=wave2(k-halfL2:k+halfL2);
    vec2=vec2.*BlackWin2;

    EvenMoment(k)=mean(vec2);    
    OddMoment(k)=mean(vec);     
end

pos=find(EvenMoment==0);
EvenMoment(pos)=[];
OddMoment(pos)=[];

% Remove the low-frequency contents of these oscillating moments
Fs2=Fs/Step;

Ws = 40/(Fs2/2);
Wp = 50/(Fs2/2);
Rp = 3; Rs = 60;
[n,Wp] = ellipord(Wp,Ws,Rp,Rs);
[b,a] = ellip(n,Rp,Rs,Wp,'high');

OddMoment=filtfilt(b,a,OddMoment);
EvenMoment=filtfilt(b,a,EvenMoment);


% For each frame, estimate the phase shift between the two oscillating
% moments
L=round(25/1000*Fs2);
Shift=round(10/1000*Fs2);
Start=1;
Stop=L;
Ind=1;
IndTot=1;
Win=hanning(L)';
Nframes=ceil((length(OddMoment))/Shift);
Phi=zeros(1,Nframes);
while Stop<length(OddMoment)
    
    if VUVDecisions(IndTot)~=0
        Seg1=EvenMoment(Start:Stop);
        Seg1=Seg1.*Win;
        Seg2=OddMoment(Start:Stop);
        Seg2=Seg2.*Win;
        
        Corr=xcorr(Seg1,Seg2);
        Corr2=xcorr(Seg2,Seg2);
        
        [maxi,posi]=max(Corr);
        posi=posi-L;
        
        Corr2(1:L)=[];
        
        pos2=2;
        while (Corr2(pos2)<Corr2(pos2-1))||(Corr2(pos2)<Corr2(pos2+1))
            pos2=pos2+1;
            if pos2==length(Corr2)
                break
            end
        end
        
        Phi(Ind)=posi/pos2;
        Ind=Ind+1;
    end
    
    IndTot=IndTot+1;
    Start=Start+Shift;
    Stop=Stop+Shift;
end

Phi(Ind:end) = [];

Phi=mod(Phi,1);
Phi=Phi-0.5;


% Inspect the phase shifts to determine the polarity
pos=find((Phi>-0.12)&(Phi<0.38));

if length(pos)>0.5*length(Phi)
    Polarity=1;
    display(['Positive polarity with a probability of ' num2str(100*length(pos)/length(Phi)) '%'])
    pause(0.00001)
    
else
    Polarity=-1;
    display(['Negative polarity with a probability of ' num2str(100*(length(Phi)-length(pos))/length(Phi)) '%'])
    pause(0.00001)
end
