function [deriv] = fftDeriv (vec, dx)

if mod(length(vec),2) ~= 0
    error('Length(vec) must be even')
end
if isrow(vec)
    vec = vec';
end

N_orig = length(vec);

%N_buff = ceil(0.5*N_orig);
N_pow2 = 2^ceil(log((N_orig))/log(2D0));
N_buff = (N_pow2 - N_orig) / 2;
extrap1 = vec(1) - (1:N_buff)*(vec(2)-vec(1));
extrap2 = vec(end) + (1:N_buff)*(vec(end)-vec(end-1));
%vec = [flipud(extrap1'); vec; ones(N_buff,1)*vec(end)];
vec = [flipud(extrap1'); vec; extrap2'];
%figure;plot(1:length(vec),vec)

vec = [vec; flipud(vec)];

%a = 0;
%b = 4*pi;
%N = 4096;
%dx = (b-a)/N;
%t = a + dx*(0:N-1);
%f = sin(t);
N = length(vec);
fftx = fft(vec);
k = 2*pi*[0:N/2-1, 0, -N/2+1:-1]/(dx*N);
deriv = ifft(i*k'.*fftx);
%figure;plot(1:N,deriv);
deriv = deriv(N_buff+(1:N_orig));

endfunction
