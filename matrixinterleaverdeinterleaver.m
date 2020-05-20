clc;
clear;

%% Load 2D image
I = imread("lena.bmp");
I_gray = mat2gray(I);
figure;
imshow(I_gray);

%% Set compression ratio
Rcompression=0.5;
N1 = round(Rcompression*16*16);
Nc = 256-N1;
%% Block based 2D DCT
[m, n] = size(I);
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
        cnt = cnt+Nc;
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
j=1;
for i=1:k:length(seq_1d)
    packets(j,:) = seq_1d(1,i:i+k-1);
    j=j+1;
end

%-----Creating square matrices-----
packets_sq1 = packets(1:255,:);
packets_sq2 = zeros(255,127);
packets_sq2(1:4,:) = packets(256:259,:);

%% Block RS Encoding
m= 8;%bits per symbol
n= (2^m)-1;
k= 127;
msgwords1=gf(packets_sq1, m);
codes1 =rsenc(msgwords1, n, k); 
codewords1 = codes1.x;

msgwords2=gf(packets_sq2, m);
codes2 =rsenc(msgwords2, n, k); 
codewords2 = codes2.x;

%% Matrix Interleaving
codes1_tr = codewords1';
codes1_tr  = gf(codes1_tr ,8);
codes2_tr = codewords2';
codes2_tr  = gf(codes2_tr ,8);

%% loss packet
choice = input('Set Noice: enter 1 for packet error and 2 for channel bit error ');
switch choice
    case 1
        disp("Packet error losses");
          rand_indices = randsample(1:length(codes1), 7);
          for ind=1:length(rand_indices)
               epacket=zeros(1, length(codes1(ind,:)));
               errorpacket=gf(epacket, m);
               codes1_tr(rand_indices(ind),:) =errorpacket;
          end
          rand_indices = randsample(1:length(codes2), 1);
          for ind=1:length(rand_indices)
               epacket=zeros(1, length(codes2(ind,:)));
               errorpacket=gf(epacket, m);
               codes2_tr(rand_indices(ind),:) =errorpacket;
          end
    case 2
        disp("Channel bit error")
        noise= (1 +randi(254,n,(2^m)-1)).*randerr(255,n, 60);
        noise=gf(noise, m);
        codes1_tr=codes1_tr+noise;
        codes2_tr=codes2_tr+noise;      
    otherwise
        disp("wrong choice");
end
%% Matrix de-interleaving
  temp = codes1_tr.x';
  codes1 = gf(temp,8);
  temp = codes2_tr.x';
  codes2 = gf(temp,8);
%% Block RS Decoding
[decmsg1,cnumerr] = rsdec(codes1,n,k);
[decmsg2,cnumerr1] = rsdec(codes2,n,k);

%% Block Depacketization
temp1 = decmsg1.x
temp2 = decmsg2.x
for i1=1:k
    de_packets1(:,i1) = temp1(:,i1);
end
for i1=1:k
    de_packets2(:,i1) = temp2(:,i1);
end
re_seq = zeros(1,length(index));
 j2=1;
for i2=1:255
    re_seq(j2:j2+k-1) = de_packets1(i2,:);
    j2=j2+k;
end
for i2=1:3
    re_seq(j2:j2+k-1) = de_packets2(i2,:);
    j2=j2+k;
end
re_seq(j2:j2+1) = de_packets2(4,1:2);

%% Block Quantization indices to quantized values
sqdec = dsp.ScalarQuantizerDecoder;
sqdec.CodebookSource = 'Input port';
reconstructed_sig = sqdec(uint8(re_seq),codebook);
%% Inverse DCT
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
title("Reconstructed Image with 3% packet loss");

%%
peaksnr = psnr(idct_img,I_gray)
for i =1:16:256
    for j=1:16:256
    SSIM(i,j)=ssim(I_gray(i:i+15,j:j+15),idct_img(i:i+15,j:j+15));
    end
end
MSSIM=sum(SSIM,'all')/(16*16)
