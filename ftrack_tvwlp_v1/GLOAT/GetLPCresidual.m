function [res] = GetLPCresidual(wave,L,shift,order,gci,type,t0)

% %%%
%  
% Use: [res] = GetLPCresidual(wave,L,shift,order,gci,type,t0)
% 
% 
% L=window length (samples) (typ.25ms)
% shift=window shift (samples) (typ.5ms)
% order= LPC order
% gci=gci position (samples)
% type=vector of voicing decisions (=0 if Unvoiced, =1 if Voiced)
% t0=vector of period values (in samples)
% 
% Written by Thomas Drugman, TCTS Lab.
% 
% %%%

if nargin<5
    doPS=0;
else
    doPS=1;
end

if doPS==0
    start=1;
    stop=start+L;   
    
    HannWin = hanning(L+1);
    
    res=zeros(1,length(wave));
    while stop<length(wave)

        segment=wave(start:stop);
% size(segment)
% size(HannWin)
        segment=segment(:).*HannWin(:);
        
        [A,e]=lpc(segment,order);
        
        inv=filter(A,1,segment);
        
        inv=inv*sqrt(sum(segment.^2)/sum(inv.^2));

        res(start:stop)=res(start:stop)+inv';

        start=start+shift;
        stop=stop+shift;
    end

    res=res/max(abs(res));

else
    
    Ltot=length(wave);
    [begin,ending,frametype] = CompleteFramingWithType(gci,16000,Ltot,shift,L,type,t0);
    
    res=zeros(1,length(wave));
    for k=1:length(begin)
       
        start=begin(k);
        stop=ending(k);
        
        segment=wave(start:stop);
        segment=segment.*hanning(length(segment)); % length(segment) is changing f(Pitch) --> cannot be optimized

        [A,e]=lpc(segment,order);

        inv=filter(A,1,segment);

        res(start:stop)=res(start:stop)+inv';
    end
        
    res=res/max(abs(res));
    
end
