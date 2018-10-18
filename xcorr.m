function [ c,d ] = xcorrC186( a,type )
 
if strcmp(type,'coeff')
corrLength=length(a)+length(a)-1;

c=fftshift(ifft(fft(a,corrLength).*conj(fft(a,corrLength))));
c=c/(c(length(a)));
d=-length(a):length(a);
end
end

