%5ο ερώτημα
%α
%--------------------------------------------------------------------------
AM = 319051;
sum_of_AM_dig = mod(AM,1000);
while sum_of_AM_dig > 9
    sum_of_AM_dig = mod(sum_of_AM_dig,10) + fix(sum_of_AM_dig/10);
end
if mod(sum_of_AM_dig,2)==0
    fileName = 'soundfile2_lab2.wav';
else
    fileName = 'soundfile1_lab2.wav';
end
[sound_sampled Fs] = audioread(fileName);
time_sound = 0:1/Fs:(length(sound_sampled)-1)/Fs;
%--------------------------------------------------------------------------
%β
%--------------------------------------------------------------------------
bits = 8;
Mmax = max(sound_sampled);
Mmin = min(sound_sampled);
levels = 2^bits;
D = (Mmax-Mmin)/levels-1;
partition = linspace(Mmin+D/2,Mmax-D/2,levels-1);
codebook = linspace(Mmin,Mmax,levels);
[index,sound_quantized] = quantiz(sound_sampled,partition,codebook);
figure;
plot(time_sound, sound_quantized)
%--------------------------------------------------------------------------
%γ
%--------------------------------------------------------------------------
sound_bits = dec2bin(index);
sound_bits = sound_bits';
sound_bits_string = sound_bits(:)';
sound_bits_vector = sound_bits_string-48;
qpskModulatorSound = comm.QPSKModulator('BitInput',true);
qpskModulatedSound = qpskModulatorSound(sound_bits_vector');
constDiagramSound = comm.ConstellationDiagram('Title' ...
    ,'Διάγραμμα Αστερισμού QPSK Κβαντισμένου Ηχητικού Σήματος ' ...
    ,'SamplesPerSymbol',2 ...
    ,'SymbolsToDisplaySource', 'Property' ...
    ,'SymbolsToDisplay', 50000);
constDiagramSound(qpskModulatedSound)
%--------------------------------------------------------------------------

%δ
%--------------------------------------------------------------------------
qpskDemodulatorSound = comm.QPSKDemodulator('BitOutput',true);
constDiagramSound2 = comm.ConstellationDiagram('Title', ...
    'Διάγραμμα Αστερισμού QPSK Κβαντισμένου Ηχητικού Σήματος Es/N0 = 5db' ...
    ,'SamplesPerSymbol',2 ...
    ,'SymbolsToDisplaySource', 'Property' ...
    ,'SymbolsToDisplay', 10^5);
constDiagramSound3 = comm.ConstellationDiagram('Title', ...
    'Διάγραμμα Αστερισμού QPSK Κβαντισμένου Σήματος Es/N0 = 15db' ...
    ,'SamplesPerSymbol',2 ...
    ,'SymbolsToDisplaySource', 'Property' ...
    ,'SymbolsToDisplay', 10^5);

Es_N0 = [4 14];
Es_N0_lin = 10.^(Es_N0./10);
BER = [0 0];
theor_BER = qfunc(sqrt(2*Es_N0_lin));

%πρόσθεση θορύβου
noise_channel_sound1 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Es/No)' ...
    ,'EsNo',Es_N0(1));
qpskModulatedSoundNoised1 = noise_channel_sound1(qpskModulatedSound);

noise_channel_sound2 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Es/No)' ...
    ,'EsNo',Es_N0(2));
qpskModulatedSoundNoised2 = noise_channel_sound2(qpskModulatedSound);

%αποδιαμόρφωση και ανίχνευση σφαλμάτων
qpskDemodulatedSound1 = qpskDemodulatorSound(qpskModulatedSoundNoised1);
errorRate = comm.ErrorRate;
errorStats = errorRate(sound_bits_vector',qpskDemodulatedSound1);
BER(1) = errorStats(1);

qpskDemodulatedSound2 = qpskDemodulatorSound(qpskModulatedSoundNoised2);
errorRate = comm.ErrorRate;
errorStats = errorRate(sound_bits_vector',qpskDemodulatedSound2);
BER(2) = errorStats(1);

%διαγράμματα αστερισμού τεθορυβεμένων σημάτων
% constDiagramSound2(qpskModulatedSoundNoised1);
% constDiagramSound3(qpskModulatedSoundNoised2);

recovered_indexes_1 = bin2dec(char( ...
    reshape(qpskDemodulatedSound1,[bits, length(qpskDemodulatedSound1)/bits])'+48))';
recovered_indexes_2 = bin2dec(char( ...
    reshape(qpskDemodulatedSound2,[bits, length(qpskDemodulatedSound2)/bits])'+48))';

recover1_sound = 2*Mmax*(recovered_indexes_1/(2^bits-1)-0.5);
recover2_sound = 2*Mmax*(recovered_indexes_2/(2^bits-1)-0.5);

audiowrite('recovered_sound1_test.wav',recover1_sound,Fs);
audiowrite('recovered_sound2_test.wav',recover2_sound,Fs);


