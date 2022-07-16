%4ο ερώτημα
AM = 319051;

%α
%--------------------------------------------------------------------------
%random_bits = randi([0 1], 36 ,1);
random_bits = [1 1 0 1 0 1 1 0 0 0 0 0 1 1 1 1 0 1 0 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 1 1];

qpsk_mod = comm.QPSKModulator('PhaseOffset',pi/4,'BitInput',true,'SymbolMapping','Gray');
qpsk_modData = qpsk_mod(random_bits');
constDiagram1 = comm.ConstellationDiagram('Title', 'Διάγραμμα Αστερισμού QPSK');
%constDiagram1(qpsk_modData)
%--------------------------------------------------------------------------

%β
%--------------------------------------------------------------------------
sum_of_AM_dig = mod(AM,1000);
while sum_of_AM_dig > 9
    sum_of_AM_dig = mod(sum_of_AM_dig,10) + fix(sum_of_AM_dig/10);
end
Tb = 0.25;
Ebit = sum_of_AM_dig^2*Tb;

Eb_N0_db = [5 15];
snr_db = Eb_N0_db;
snr_lin = 10.^(snr_db./10);
N0_lin = Ebit ./ snr_lin;
%N0_db = 10*log10(N0_lin);
samples = length(qpsk_modData);
X = sqrt(N0_lin/4).*randn(samples,1);
Y = sqrt(N0_lin/4).*randn(samples,1);
Z = X + 1i*Y;

qpsk_with_wgn = qpsk_modData+Z;
max_qpsk_wgn_x = max(abs(real(qpsk_with_wgn)));
max_qpsk_wgn_y = max(abs(imag(qpsk_with_wgn)));
constDiagram2 = comm.ConstellationDiagram('Title','Διάγραμμα Αστερισμού QPSK Eb/No=5db' ...
    ,'ShowTrajectory',false ...
    ,'XLimits', [-max_qpsk_wgn_x(1) max_qpsk_wgn_x(1)] ...
    ,'YLimits', [-max_qpsk_wgn_y(1) max_qpsk_wgn_y(1)]);
%constDiagram2(qpsk_with_wgn(:,1))

constDiagram3 = comm.ConstellationDiagram('Title','Διάγραμμα Αστερισμού QPSK Eb/No=15db' ...
    ,'ShowTrajectory',false ...
    ,'XLimits', [-max_qpsk_wgn_x(2) max_qpsk_wgn_x(2)] ...
    ,'YLimits', [-max_qpsk_wgn_y(2) max_qpsk_wgn_y(2)]);
%constDiagram3(qpsk_with_wgn(:,2))
%--------------------------------------------------------------------------


%γ
%--------------------------------------------------------------------------
num_of_bits = 10^6;
random_bits = randi([0 1], num_of_bits,1);

qpsk_mod = comm.QPSKModulator('BitInput', true);
qpsk_demod = comm.QPSKDemodulator('BitOutput',true);
qpsk_modData = qpsk_mod(random_bits);
BER = zeros(16,1);

Ebit = sum_of_AM_dig^2*Tb/2;
Eb_N0_db = 0:15;
Eb_N0_lin = 10.^(Eb_N0_db./10);
snr_db = Eb_N0_db;
snr_lin = 10.^(snr_db./10);
N0_lin = Ebit ./ snr_lin;
N0_db = 10*log10(N0_lin);

index = 1;
for eb_no = Eb_N0_db
    noise_channel = comm.AWGNChannel('EbNo',eb_no,'BitsPerSymbol',2);
    qpsk_modDataNoised = noise_channel(qpsk_modData);
    qpsk_demodData = qpsk_demod(qpsk_modDataNoised);
    errorRate = comm.ErrorRate;
    errorStats = errorRate(random_bits,qpsk_demodData);
    BER(index) = errorStats(1);
    index = index + 1;
end
theor_ber = qfunc(sqrt(2*Eb_N0_lin));

% figure()
% hold on;
% plot(snr_db,BER+eps,'*')
% plot(snr_db,theor_ber,'LineWidth',2)
% hold off;
% set(gca,'yscale','log')
% grid on;
%--------------------------------------------------------------------------

%δ1
%--------------------------------------------------------------------------
sum_of_AM_dig = mod(AM,1000);
while sum_of_AM_dig > 9
    sum_of_AM_dig = mod(sum_of_AM_dig,10) + fix(sum_of_AM_dig/10);
end
text_arr = [];
if mod(sum_of_AM_dig,2)==0
    fileName = 'rice_even.txt';
else
    fileName = 'rice_odd.txt';
end
fileID = fopen(fileName,'r');
formatSpec = '%c';
text_arr = fscanf(fileID,formatSpec);
fclose(fileID);

ascii_text = double(text_arr);
[binary_vector,binary_string] = text2bin(ascii_text);

bits = 8;
Mmax = 2^bits-1;
Mmin = 0;
levels = 2^bits;
D = (Mmax-Mmin)/levels-1;
partition = linspace(Mmin+D/2,Mmax-D/2,levels-1);
codebook = linspace(Mmin,Mmax,levels);
[index,quantized_text] = quantiz(ascii_text,partition,codebook);
plot(quantized_text)
%--------------------------------------------------------------------------

%δ1
%--------------------------------------------------------------------------
A=1;
qpskModulator = comm.QPSKModulator('BitInput',true);
qpskModulatedData = qpskModulator(binary_vector');
constDiagram4 = comm.ConstellationDiagram('Title', 'Διάγραμμα Αστερισμού QPSK Κβαντισμένου Σήματος');
% constDiagram4(qpskModulatedData);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
qpskDemodulator = comm.QPSKDemodulator('BitOutput',true);
constDiagram5 = comm.ConstellationDiagram('Title', ...
    'Διάγραμμα Αστερισμού QPSK Κβαντισμένου Σήματος Es/N0 = 5db' ...
    ,'SamplesPerSymbol',2);
constDiagram6 = comm.ConstellationDiagram('Title', ...
    'Διάγραμμα Αστερισμού QPSK Κβαντισμένου Σήματος Es/N0 = 15db' ...
    ,'SamplesPerSymbol',2);

Es_N0 = [5 15];
Es_N0_lin = 10.^(Es_N0./10);
BER = [0 0];
theor_BER = qfunc(sqrt(2*Es_N0_lin));

%πρόσθεση θορύβου
noise_channel1 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Es/No)' ...
    ,'EsNo',Es_N0(1));
qpskModulatedDataNoised1 = noise_channel1(qpskModulatedData);

noise_channel2 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Es/No)' ...
    ,'EsNo',Es_N0(2));
qpskModulatedDataNoised2 = noise_channel2(qpskModulatedData);

%αποδιαμόρφωση και ανίχνευση σφαλμάτων
qpskDemodulatedData1 = qpskDemodulator(qpskModulatedDataNoised1);
errorRate = comm.ErrorRate;
errorStats = errorRate(binary_vector',qpskDemodulatedData1);
BER(1) = errorStats(1);

qpskDemodulatedData2 = qpskDemodulator(qpskModulatedDataNoised2);
errorRate = comm.ErrorRate;
errorStats = errorRate(binary_vector',qpskDemodulatedData2);
BER(2) = errorStats(1);

%διαγράμματα αστερισμού τεθορυβεμένων σημάτων
constDiagram5(qpskModulatedDataNoised1);
constDiagram6(qpskModulatedDataNoised2);

%--------------------------------------------------------------------------

text_recovered_1 = bin2text(qpskDemodulatedData1);
text_recovered_2 = bin2text(qpskDemodulatedData2);

fileID = fopen('recovered_text_test_1.txt', 'w');
fwrite(fileID, text_recovered_1);
fclose(fileID);
fileID = fopen('recovered_text_test_2.txt', 'w');
fwrite(fileID, text_recovered_2);
fclose(fileID);

function text = bin2text(binVS)
    btxt  = reshape(binVS,[8, length(binVS)/8])';
    if length(class(btxt))== 6
        text  = char(bin2dec(char(btxt+48)))';
    else
        text  = char(bin2dec(btxt))';
    end
end

function [binV, binS] = text2bin(text)
    binS = dec2bin(text,8);
    binS = binS';
    binS = binS(:)';
    binV = binS-48;
end