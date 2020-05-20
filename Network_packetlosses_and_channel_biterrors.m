clc;
clear;
%% Load 2D image
I = imread("lena.bmp");
I_gray = mat2gray(I);
figure;
imshow(I_gray);

%% Set compression ratio
Rcompression=0.5;
N1 = round(0.5*16*16);
Nc = 256-N1;
%% Block based 2D DCT
[m,n] = size(I);
for i =1:16:m
    for j=1:16:n
    dct_trimg(i:i+15,j:j+15)=dct2(I_gray(i:i+15,j:j+15));
    end
end

%% zigzag scanning DCT coeficients and inverse zigzag
cnt = 1;
for i = 1:16:256
    for j = 1:16:256 
        temp = zigzag(dct_trimg(i:i+15,j:j+15));
        zig_mat(cnt:cnt +Nc -1 ) = temp(1:Nc);
        cnt = cnt +Nc;
    end
end
signal = zig_mat;

%% Block Quantization
samples = signal;
maxsig = max(samples);
minsig = min(samples);
interv=(maxsig-minsig)/(2^8); %interval length for 256 levels resolution
u1=minsig-interv/2;
u2=minsig+0.001;
partition = [u2:interv:maxsig]; 
codebook = [u1:interv:maxsig]; 
[index, quantized]= quantiz(samples,partition,codebook);

%% Block Packetizer
k = 127;
index=uint8(index);
packet_row = ceil(length(index)/k);
seq_1d = zeros(1,(k*packet_row));
seq_1d(1,1:length(index)) = index;
packets = zeros(k,k);
L = length(seq_1d);
j=1;
for i=1:k:L
    packets(j,:) = seq_1d(1,i:i+k-1);
    j=j+1;
end

%% Block RS Encoding
m= 8;%bits per symbol
n= (2^m)-1;% codeword length
k= 127;%message length
msgwords=gf(packets, m);
codes =rsenc(msgwords, n, k);
codewords = codes.x; 

%% Introducing Noice
choice = input('Set Noice: enter 1 for packet error and 2 for channel bit error');
switch choice
    case 1
        disp("Packet error losses");
          rand_indices = randsample(1:length(codes),8);
          for ind=1:length(rand_indices)
               epacket=zeros(1, length(codewords(ind,:)));
               errorpacket=gf(epacket, m);
               codes(rand_indices(ind),:) =errorpacket;
          end
    case 2
        disp("Channel bit error")
        noise= (1 +randi(254,n,(2^m)-1)).*randerr(255,n, 60);
        noise=gf(noise, m);
        noise(256:259,:)=0;
        codes=codes+noise;
    otherwise
        disp("wrong choice");
end

%% Block RS DEcoding
[decmsg,cnumerr] = rsdec(codes,n,k);

%% Block Depacketization
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

%% Block Quantization indices to quantized values
sqdec = dsp.ScalarQuantizerDecoder;
sqdec.CodebookSource = 'Input port';
reconstructed_sig = sqdec(uint8(re_seq),codebook);

%% Inverse zigzag
cnt = 1;
temp = zeros(16,16);
for i = 1:16:256
    for j = 1:16:256
        temp(1:Nc) = reconstructed_sig(cnt:cnt +Nc -1 );
        izig_mat(i:i+15,j:j+15) = izigzag(temp,16,16);
        cnt = cnt + Nc;
    end
end
%% inverce IDCT
for i =1:16:256
    for j=1:16:256
      idct_img(i:i+15,j:j+15)= idct2(izig_mat(i:i+15,j:j+15));
    end
end
figure
subplot(1,2,1)
imshow(I_gray)
title("Original Image");
subplot(1,2,2)
imshow(idct_img);
title("Reconstructed Image");


%%
peaksnr = psnr(idct_img,I_gray)
for i =1:16:256
    for j=1:16:256
    SSIM(i,j)=ssim(I_gray(i:i+15,j:j+15),idct_img(i:i+15,j:j+15));
    end
end
MSSIM=sum(SSIM,'all')/(16*16)
