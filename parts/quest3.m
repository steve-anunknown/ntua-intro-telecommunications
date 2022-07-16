%3ο ερώτημα
AM = 319051;
%a
%random_bits = randi([0 1], 36,1);
random_bits = [1 1 0 1 0 1 1 0 0 0 0 0 1 1 1 1 0 1 0 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 1 1];
Tb = 0.25;
A = 1;
temp = int32(AM);
sum_temp = 0;
decimal_place = 10;
while not(temp==0)
    sum_temp = sum_temp + mod(temp,10);
    temp = idivide(temp,10);
end
if mod(sum_temp,2)==0
    fc = 1;
else
    fc = 2;
end

random_bits_qpsk = 2*random_bits(1:2:end-1) + random_bits(2:2:end);   % συνένωση δύο υπακολουθιών
random_bits_qpsk_gray = dec2graydec(random_bits_qpsk);   % απόδοση δεικτών gray
random_bits_8psk = 4*random_bits(1:3:end-2) + 2*random_bits(2:3:end-1) + random_bits(3:3:end);    %συνένωση τριών υπακολουθιών
random_bits_8psk_gray = dec2graydec(random_bits_8psk);  % απόδοση δεικτών gray

dt = Tb / 250;
t_bpsk = 0:dt:length(random_bits)*Tb-dt;
bpsk = A*sin(2*pi*fc*t_bpsk + pi*random_bits(1+fix(t_bpsk/Tb)));

t_qpsk = 0:dt:length(random_bits_qpsk)*Tb-dt;
qpsk = A*sin(2*pi*fc*t_qpsk + pi/4*random_bits_qpsk_gray(1+fix(t_qpsk/Tb)));

t_8psk = 0:dt:length(random_bits_8psk)*Tb-dt;
psk8 = A*sin(2*pi*fc*t_8psk + pi/8*random_bits_8psk_gray(1+fix(t_8psk/Tb)));

figure;
subplot(3,1,1);
plot(t_bpsk,bpsk)
title('Διαμόρφωση με BPSK')
subplot(3,1,2)
plot(t_qpsk,qpsk)
title('Διαμόρφωση με QPSK')
ylabel('Πλάτος (V)')
subplot(3,1,3)
plot(t_8psk,psk8)
title('Διαμόρφωση με 8-PSK')
xlabel('Χρόνος (sec)')

%b
AM = 319051;
A = mod(AM,1000);
while A > 9
    A = mod(A,10) + fix(A/10);
end
t_bpam = t_bpsk;
bpam = A*(2*random_bits(1+fix(t_bpam/Tb))-1);    % κάνουμε τον παλμό από (1)-(0) -> (1)-(-1)
Ebit = A^2*Tb;
figure;
plot(t_bpam,bpam);
title('Διαμόρφωση B-PAM');
xlabel('Χρόνος (sec)');
ylabel('Πλάτος (volt)');
grid on;

%c
scatterplot(bpam);
title('Διάγραμμα Αστερισμού B-PAM Διαμόρφωσης');
grid on;

%d
Eb_N0_db = [5 15];
snr_db = Eb_N0_db;
snr_lin = 10.^(snr_db./10);
N0_lin = Ebit ./ snr_lin;
%N0_db = 10*log10(N0_lin);
numSamples = length(bpam);
wgn = sqrt(N0_lin/4).*(randn(numSamples,1)+1i*randn(numSamples,1));
bpam_noised = bpam'+wgn;

figure;
subplot(3,1,1)
plot(t_bpam,bpam);
title('Διαμόρφωση B-PAM');
xlabel('Χρόνος (sec)');
ylabel('Πλάτος (volt)');
grid on;
subplot(3,1,2)
plot(t_bpam,real(bpam_noised(:,1)));
title('Διαμόρφωση B-PAM με AWGN_{b}/N_{0}=5db');
xlabel('Χρόνος (sec)');
ylabel('Πλάτος (volt)');
subplot(3,1,3)
plot(t_bpam,real(bpam_noised(:,2)));
title('Διαμόρφωση B-PAM με AWGN E_{b}/N_{0}=15db');
xlabel('Χρόνος (sec)');
ylabel('Πλάτος (volt)');

%e
maxx = max(abs(real(bpam_noised(:,1))));
maxy = max(abs(imag(bpam_noised(:,1))));
figure;
scatter(real(bpam_noised(:,1)), imag(wgn(:,1)),2,'filled')
xlim([-maxx maxx]);
ylim([-maxy maxy]);
title('Διάγραμμα Αστερισμού B-PAM με AWGN E_{b}/N_{0}=5db');
figure;
maxx = max(abs(real(bpam_noised(:,2))));
maxy = max(abs(imag(bpam_noised(:,2))));
scatter(real(bpam_noised(:,2)), imag(wgn(:,2)),2,'filled')
xlim([-maxx maxx]);
ylim([-maxy maxy]);
title('Διάγραμμα Αστερισμού B-PAM με AWGN E_{b}/N_{0}=15db');

random_bits = randi([0 1],50,1);
dt = 5e-5;
t_bpam = 0:dt:length(random_bits)*Tb-dt;
bpam = A*(2*random_bits(1+fix(t_bpam/Tb))-1);
 
Eb_N0_db = 0:15;
Eb_N0_lin = 10.^(Eb_N0_db./10);
snr_db = Eb_N0_db;
snr_lin = 10.^(snr_db./10);
N0_lin = Ebit ./ snr_lin;
N0_db = 10*log10(N0_lin);
 
numSamples = length(bpam);
noise_factor =sqrt(N0_lin*2);
wgn = noise_factor.*(randn(numSamples,1)+1i*randn(numSamples,1));
bpam_noised = bpam+real(wgn);
 
bit_error_rate = (bpam.*bpam_noised)<0;
bit_error_rate = sum(bit_error_rate)/length(bpam);
theor_ber = qfunc(sqrt(2*Eb_N0_lin));
figure()
hold on;
plot(snr_db,bit_error_rate+eps,'LineWidth',2)
plot(snr_db,qfunc(sqrt(2*10.^(snr_db/10))),'LineWidth',2)
hold off;
set(gca,'yscale','log')
grid on;















Eb_N0_db = [5 15];
snr_db = Eb_N0_db;
snr_lin = 10.^(snr_db./10);
N0_lin = Ebit ./ snr_lin;
% N0_db = 10*log10(N0_lin);

% N0 = Ebit./Eb_N0;   %db
% N0_lin = 10.^(N0./10);
% snr = Eb_N0;   %db
% nr = Ebit./10.^(N0./10) + 3;
numSamples = length(bpam);
wgn = sqrt(N0_lin/2).*(randn(numSamples,1)+1i*randn(numSamples,1));
wgn2_real = sqrt(N0_lin).*randn(numSamples,1);
wgn2_imag = sqrt(N0_lin).*randn(numSamples,1);
wgn2 = wgn2_real+1i*wgn2_imag;
 
 
figure;
subplot(4,1,1)
plot(t_bpam,awgn(bpam,snr(1),'measured'));
title('awgn(bpam,snr(1),"measured")');
axis tight;
subplot(4,1,2)
plot(t_bpam,bpam'+real(wgn(:,1)));
axis tight;
title('bpam. + real(wgn(:,1))');
subplot(4,1,3)
plot(t_bpam,awgn(bpam,snr(1),'measured'));
title('awgn(bpam,snr(1),"measured")');
axis tight;
subplot(4,1,4)
plot(t_bpam,bpam'+ real(wgn2));
title('bpam. + wgn2 real');
axis tight;

% subplot(4,1,4)
% plot(t_bpam,bpam'+wgn(:,1));
% title('bpam. + wgn(:,1)');

% figure;
% subplot(3,1,1)
% plot(t_bpam,bpam);
% title('Διαμόρφωση B-PAM');
% xlabel('Χρόνος (sec)');
% ylabel('Πλάτος (volt)');
% grid on;
% subplot(3,1,2)
% plot(t_bpam,bpam'+real(wgn(:,1)));
% title('Διαμόρφωση B-PAM με Λευκό Προσθετικό Θόρυβο Gauss E_{b}/N_{0}=5db');
% xlabel('Χρόνος (sec)');
% ylabel('Πλάτος (volt)');
% subplot(3,1,3)
% plot(t_bpam,bpam'+real(wgn(:,2)));
% title('Διαμόρφωση B-PAM με Λευκό Προσθετικό Θόρυβο Gauss E_{b}/N_{0}=15db');
% xlabel('Χρόνος (sec)');
% ylabel('Πλάτος (volt)');


% scatterplot(awgn(bpam,snr(1),'measured'));
% title('Διάγραμμα Αστερισμού B-PAM με Λευκό Προσθετικό Θόρυβο Gauss E_{b}/N_{0}=5db');
% scatterplot(awgn(bpam,snr(2),'measured'));
% title('Διάγραμμα Αστερισμού B-PAM με Λευκό Προσθετικό Θόρυβο Gauss E_{b}/N_{0}=15db');
%--------------------------------------------------------------------------
%}





function [res] = dec2graydec(number)
    [dim1, dim2] = size(number);
    max_num = ceil(log2(max(number,[],'all')));
    x = repmat('0',1,max_num);
    y = repmat('0',1,max_num);
    res = zeros(dim1,dim2);
    for i = 1:dim1
        for j = 1:dim2
            in1 = number(i,j);
            in2 = fix(number(i,j)/2);
            bin1 = dec2bin(in1);
            bin2 = dec2bin(in2);
            x(end-length(bin1)+1:end) = bin1;
            y(end-length(bin2)+1:end) = bin2;
            resi = xor(x-'0', y-'0');
            resi = bin2dec(char(resi + '0'));
            res(i,j) = resi;
        end
    end
end
