clc;
clear;
%% Task-1
t= [0:0.1:10*pi];    %Times at which to sample the sine function
signal=sin(t);       %Original signal, a sine wave
samples = signal;
maxsig=max(samples);
minsig = min(samples);
interv=(maxsig-minsig)/(2^8-1); %interval length for 256 levels resolution
u=minsig-interv;
partition = [minsig:interv:maxsig]; 
codebook = [u:interv:maxsig]; 
[index, quantized]= quantiz(samples,partition,codebook);

figure
plot(t,signal,'-',t,quantized,'-');
legend('Original signal','Quantized signal');


%% Task 2 - packetization
k = 127;
packet_row = ceil(length(index)/k);
sequence = zeros(1,(k*packet_row));
sequence(1,1:length(index)) = index;
packets = zeros(k,k);
L = length(sequence);
j=1;
%Packetizing
for i=1:k:L
    packets(j,:) = sequence(1,i:i+k-1);
    j=j+1;
end

%% Task - 3
%RS Encoding
m= 8;%bits per symbol
n= (2^m)-1;%codeword length 
k= 127;%message length
msgwords=gf(packets, m);%represent data by using a Galois array
codes = rsenc(msgwords,n,k);    %perform RS encoder
codewords = codes.x; %extract rows of codewords from the GF array  

%% RS Decoding
[decmsg,cnumerr] = rsdec(codes,n,k);
isequal(decmsg,msgwords)

%% Task 2 - depacketization
temp = decmsg.x;
for i1=1:k
    de_packets(:,i1) = temp(:,i1);
end
len = length(index);
re_seq = zeros(1,len);
 j2=1;
for i2=1:packet_row-1
    re_seq(j2:j2+k-1) = de_packets(i2,:);
    j2=j2+k;
end
re_seq(j2:j2+k-((k*packet_row)-len)-1) = de_packets(packet_row,1:k-((k*packet_row)-len));

%% Converting quantization indices to quantized values
indices = uint8(re_seq);
sqdec = dsp.ScalarQuantizerDecoder;
sqdec.CodebookSource = 'Input port';
qout = sqdec(uint8(index),codebook);


figure
plot(signal);
hold on
plot(qout);
legend('original signal','quantized signal');
title("Quantized signal after RS decoding and depacketization");
