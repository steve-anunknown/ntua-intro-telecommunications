AM = 319051;

%--------------------------------------------------------------------------
%1ο ερώτημα
%α
A = 4;  %   πλάτος τριγωνικού παλμού
N = 4;  %   επιθυμητό πλήθος περιόδων
fm = mod(AM,1000);
while fm > 9
    fm = mod(fm,10) + fix(fm/10);
end
fm = fm * 1000;  %   συχνότητα όπως προκύπτει από τον ΑΜ
fs1 = 30*fm;    %   πρώτη συχνότητα δειγματοληψίας
fs2 = 50*fm;    %   δεύτερη συχνότητα δειγματοληψίας
t1 = 0:1/fs1:N/fm;  %   πρώτος χρονικός άξονας
t2 = 0:1/fs2:N/fm;  %   δεύτερος χρονικός άξονας
trig1 = A*sawtooth(2*pi*fm*t1,1/2); %   δειγματοληψία συχνότητας fs1 από τριγωνική κυματομορφή συχνότητας fs, πλάτους Α
trig2 = A*sawtooth(2*pi*fm*t2,1/2); %   δειγματοληψία συχνότητας fs2 από τριγωνική κυματομορφή συχνότητας fs, πλάτους Α

plot(t1,trig1,'o')  %   γράφημα δειγμάτων πρώτης δειγματοληψίας
grid on;
title('Δείγματα trig_1(t) κατόπιν δειγματοληψίας με συχνότητα f_{s1}=30f_m');
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');

figure;
plot(t2,trig2,'x')  %   γράφημα δειγμάτων δεύτερης δειγματοληψίας
grid on;
title('Δείγματα trig_2(t) κατόπιν δειγματοληψίας με συχνότητα f_{s2}=50f_m');
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');

figure;
plot(t1,trig1,'o',t2,trig2,'x')
grid on;
title('Κοινό διάγραμμα δειγμάτων trig_1(t), trig_2(t)'); %  κοινό γράφημα δειγμάτων αμφοτέρων των δειγματοληψιών
legend('Τριγωνική κυματομορφή συχνότητας f_{s1}=30f_m', 'Τριγωνική κυματομορφή συχνότητας f_{s2}=50f_m')
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');
%--------------------------------------------------------------------------


%β
%--------------------------------------------------------------------------
fs = 4*fm;
t = 0:1/fs:4/fm;
trig = A*sawtooth(2*pi*fm*t,1/2);
figure;
plot(t,trig,'o');
grid on;
title('Δείγματα trig(t) κατόπιν δειγματοληψίας με συχνότητα f_s=4f_m');
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');
%--------------------------------------------------------------------------


%γ1
%--------------------------------------------------------------------------
Az = 1; %   πλάτος κυματομορφής
z1 = Az*sin(2*pi*fm*t1);    %   δείγματα κατόπιν δειγματοληψίας συχνότητας fs1
z2 = Az*sin(2*pi*fm*t2);    %   δείγματα κατόπιν δειγματοληψίας συχνότητας fs2

figure;
plot(t1,z1,'o')  %  γράφημα δειγμάτων πρώτης δειγματοληψίας
grid on;
title('Δείγματα z_1(t) κατόπιν δειγματοληψίας με συχνότητα f_{s1}=30f_m');
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');

figure;
plot(t2,z2,'x')  %  γράφημα δειγμάτων δεύτερης δειγματοληψίας
grid on;
title('Δείγματα z_2(t) κατόπιν δειγματοληψίας με συχνότητα f_{s2}=50f_m');
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');

figure;
plot(t1,z1,'o',t2,z2,'x')
grid on;
title('Κοινό διάγραμμα δειγμάτων z_1(t), z_2(t)'); %    κοινό γράφημα δειγμάτων αμφοτέρων των δειγματοληψιών
legend('Ημίτονο συχνότητας f_{s1}=30f_m', 'Ημίτονο συχνότητας f_{s2}=50f_m')
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');

t = 0:1/fs:4/fm;
z = Az*sin(2*pi*fm*t);
figure;
plot(t,z,'o');
grid on;
title('Δείγματα κατόπιν z(t) δειγματοληψίας με συχνότητα f_s=4f_m');
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');
%--------------------------------------------------------------------------

%γ2
%--------------------------------------------------------------------------
L = 10^3;
T = 1/gcd(fm,fm+L);
tq1 = 0:1/fs1:fix(T*fm)/fm; %   εύρος μίας περιόδου αντί τεσσάρων με δείγματα ανά fs1
tq2 = 0:1/fs2:fix(T*fm)/fm; %   εύρος μίας περιόδου αντί τεσσάρων με δείγματα ανά fs2
zq1 = Az*sin(2*pi*fm*tq1);  %   δείγματα της Ζ (με συχνότητα fs1) στο κατάλληλο χρονικό εύρος
zq2 = Az*sin(2*pi*fm*tq2);  %   δείγματα της Ζ (με συχνότητα fs2) στο κατάλληλο χρονικό εύρος
q1 = zq1 + Az*sin(2*pi*(fm+L)*tq1);
q2 = zq2 + Az*sin(2*pi*(fm+L)*tq2);

figure;
plot(tq1,q1,'o')  % γράφημα δειγμάτων πρώτης δειγματοληψίας
grid on;
title('Δείγματα q_1(t) κατόπιν δειγματοληψίας με συχνότητα f_{s1}=30f_m');
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');

figure;
plot(tq2,q2,'x')  % γράφημα δειγμάτων δεύτερης δειγματοληψίας
grid on;
title('Δείγματα q_2(t) κατόπιν δειγματοληψίας με συχνότητα f_{s2}=50f_m');
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');

figure;
plot(tq1,q1,'o',tq2,q2,'x')
grid on;
title('Κοινό διάγραμμα δειγμάτων q_1(t), q_2(t)'); %    κοινό γράφημα δειγμάτων αμφοτέρων των δειγματοληψιών
legend('Ημίτονο συχνότητας f_{s1}=30f_m', 'Ημίτονο συχνότητας f_{s2}=50f_m')
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');

t = 0:1/fs:fix(T*fm)/fm;
q = Az*sin(2*pi*fm*t) + Az*sin(2*pi*(fm+L)*t);
figure;
plot(t,q,'o');
grid on;
title('Δείγματα q(t) κατόπιν δειγματοληψίας με συχνότητα f_s=4f_m');
ylabel('Πλάτος (Volt)');
xlabel('Χρόνος (sec)');
%--------------------------------------------------------------------------

%2ο ερώτημα
%--------------------------------------------------------------------------
y = trig1;  %   αποθήκευση των δειγμάτων σε μεταβλητή σύμφωνη με την εκφώνηση της ασκήσεως
if mod(fm/1000,2) == 0
    num_of_bits = 4;
else
    num_of_bits = 5;
end

levels = 2^num_of_bits; % επίπεδα κβάντισης
Mmax = max(y);  % μέγιστη τιμή κβαντιστέου σήματος
D = 2*Mmax/levels;  % διάστημα λήψεως αποφάσεως
partition = linspace(-Mmax+D/2,Mmax-D/2,levels-1);
codebook = linspace(-Mmax,Mmax,levels);
[index,quants] = quantiz(y,partition,codebook);
gray_code = graycode(num_of_bits);

figure;
plot(t1,quants,'.')
title('Έξοδος κβαντιστή mid-riser με είσοδο τριγωνική κυματομορφή')
xlabel('Χρόνος (sec)')
ylabel('Πλάτος (Volt)')
levels_linspace = linspace(-Mmax,Mmax,levels);
set(gca,'YTick', levels_linspace)
set(gca,'YTickLabel', gray_code)
grid on

error = y - quants; % συνολικά σφάλματα
error1 = y(1:10) - quants(1:10);    % σφάλματα για πρώτα δέκα δείγματα
error2 = y(1:20) - quants(1:20);    % σφάλματα για πρώτα είκοσι δείγματα

var_error = var(error); % διασπορά σφαλμάτων
var_error1 = var(error1);
var_error2 = var(error2);

stand_deviations = [sqrt(var_error), sqrt(var_error1), sqrt(var_error2)];   % τυπική απόκλιση σφαλμάτων

SNR = ( quants.^2 * ones(size(y,2),1) ) / ( error.^2 * ones(size(y,2),1) );
SNR1 = ( quants(1:10).^2 * ones(10,1) ) / ( error1.^2 * ones(10,1) );
SNR2 = ( quants(1:20).^2 * ones(20,1) ) / ( error2.^2 * ones(20,1) );

F = max(quants)/rms(quants);
theor_SNR = (3*2^(2*num_of_bits))/F^2;

bit_dur = 2*10^-3;  %   διάρκεια bit
ampl_volt = fm / 1000;  %   πλάτος τάσης ίσο με την συχνότητα σε kHz
bit_array = gray_code(index(1:31)+1);   %   πίνακας περιέχων τις στάθμες ανά σημείο σε κώδικα Gray
bit_stream = zeros(1,4*num_of_bits*length(bit_array)+1); %    ροή μετάδοσης
time = zeros(1,4*num_of_bits*length(bit_array)+1);   %   άξονας χρόνου
index = 2;
for bit_group = bit_array   %   για κάθε στοιχείο στον πίνακα
    bit_group = convertStringsToChars(bit_group);   %   μετατρέπουμε τις συμβολοσειρές σε διακριτούς χαρακτήρες
    for bit = 1:num_of_bits % προσπελάτουμε κάθε χαρακτήρα ανά συμβολοσειρά
        if bit_group(bit) == '1'    %   εάν το ψηφίο είναι άσσος, τότε διαμορφώνουμε αναλόγως τον παλμό
            bit_stream(index-1:index+2) = [ampl_volt ampl_volt 0 0];
        else    %   ομοίως εάν είναι μηδέν
            bit_stream(index-1:index+2) = [-ampl_volt -ampl_volt 0 0];
        end
        if index==2
            tend = 0;
        else
            tend = time(index-2);
        end
        time(index-1:index+2) = [tend tend+bit_dur/2 tend+bit_dur/2 tend+bit_dur];  %   ανανεώνουμε τον χρονικό άξονα
        index = index+4;
    end
end
figure;
plot(time, bit_stream)
title('Ροή Μετάδοσης με Κωδικοποίηση Γραμμής POLAR RZ');
xlabel('Χρόνος (sec)');
ylabel('Πλάτος (Volt)');
%--------------------------------------------------------------------------


%3ο ερώτημα
%α
%--------------------------------------------------------------------------
%random_bits = randi([0 1], 36,1);
random_bits = [1 1 0 1 0 1 1 0 0 0 0 0 1 1 1 1 0 1 0 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 1 1];
Tb = 0.25;
A = 1;
temp = int32(AM);   % αναγκαίο για την χρήση της ακέραιας διαίρεσης παρακάτων
sum_temp = 0;
decimal_place = 10;
while not(temp==0)
    sum_temp = sum_temp + mod(temp,10);
    temp = idivide(temp,10);    % ακέραια διαίρεση
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

t_qpsk = 0:dt:length(random_bits_qpsk)*2*Tb-dt;
qpsk = A*sin(2*pi*fc*t_qpsk + pi/4*random_bits_qpsk_gray(1+fix(t_qpsk/(2*Tb))));

t_8psk = 0:dt:length(random_bits_8psk)*3*Tb-dt;
psk8 = A*sin(2*pi*fc*t_8psk + pi/8*random_bits_8psk_gray(1+fix(t_8psk/(3*Tb))));

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
%--------------------------------------------------------------------------

%β
%--------------------------------------------------------------------------
AM = 319051;
A = mod(AM,1000);
while A > 9
    A = mod(A,10) + fix(A/10);
end

t_bpam = t_bpsk;
bpam = A*(2*random_bits(1+fix(t_bpam/Tb))-1);    % κάνουμε τον παλμό από (1)-(0) -> (1)-(-1)
Ebit = A^2*Tb;  % ενέργεια ψηφίου 

figure;
plot(t_bpam,bpam);
title('Διαμόρφωση B-PAM');
xlabel('Χρόνος (sec)');
ylabel('Πλάτος (volt)');
grid on;
%--------------------------------------------------------------------------

%γ
%--------------------------------------------------------------------------
scatterplot(bpam);  % διάγραμμα αστερισμού
title('Διάγραμμα Αστερισμού B-PAM Διαμόρφωσης');
grid on;
%--------------------------------------------------------------------------

%δ
%--------------------------------------------------------------------------
Eb_N0_db = [5 15];
snr_db = Eb_N0_db;
snr_lin = 10.^(snr_db./10);
N0_lin = Ebit ./ snr_lin;
% N0_db = 10*log10(N0_lin);
numSamples = length(bpam);
wgn = sqrt(N0_lin/4).*(randn(numSamples,1)+1i*randn(numSamples,1)); % ιδιόφτιαχτος θόρυβος, αργότερα θα χρησιμοποιηθούν αντικείμενα
bpam_noised = bpam'+wgn;

figure;
subplot(3,1,1)
plot(t_bpam,bpam);
title('Διαμόρφωση B-PAM');
xlabel('Χρόνος (sec)');
ylabel('Πλάτος (volt)');
grid on;
subplot(3,1,2)
plot(t_bpam,real(bpam_noised(:,1)));    % χρήση μόνο του πραγματικού του μέρους
title('Διαμόρφωση B-PAM με AWGN E_{b}/N_{0}=5db');
xlabel('Χρόνος (sec)');
ylabel('Πλάτος (volt)');
subplot(3,1,3)
plot(t_bpam,real(bpam_noised(:,2)));    % χρήση μόνο του πραγματικού του μέρους
title('Διαμόρφωση B-PAM με AWGN E_{b}/N_{0}=15db');
xlabel('Χρόνος (sec)');
ylabel('Πλάτος (volt)');
%--------------------------------------------------------------------------

%ε
%--------------------------------------------------------------------------
maxx = max(abs(real(bpam_noised(:,1))));
maxy = max(abs(imag(bpam_noised(:,1))));
figure;
scatter(real(bpam_noised(:,1)), imag(wgn(:,1)),2,'filled')  % διάγραμμα αστερισμού
xlim([-maxx maxx]);
ylim([-maxy maxy]);
title('Διάγραμμα Αστερισμού B-PAM με AWGN E_{b}/N_{0}=5db');
figure;
maxx = max(abs(real(bpam_noised(:,2))));
maxy = max(abs(imag(bpam_noised(:,2))));
scatter(real(bpam_noised(:,2)), imag(wgn(:,2)),2,'filled')  % διάγραμμα αστερισμού
xlim([-maxx maxx]);
ylim([-maxy maxy]);
title('Διάγραμμα Αστερισμού B-PAM με AWGN E_{b}/N_{0}=15db');

random_bits = randi([0 1],1000,1);  % παραγωγή χιλίων τυχαίων δυαδικών ψηφίων
dt = (5e-5);
t_bpam = 0:dt:length(random_bits)*Tb-dt;
bpam = A*(2*random_bits(1+fix(t_bpam/Tb))-1);
 
Eb_N0_db = 0:15;
Eb_N0_lin = 10.^(Eb_N0_db./10);
snr_db = Eb_N0_db;
snr_lin = 10.^(snr_db./10);
N0_lin = Ebit ./ snr_lin;
%N0_db = 10*log10(N0_lin);
 
numSamples = length(bpam);
noise_factor = sqrt(N0_lin*2);  % καθορίζων παράγοντας εντάσεως θορύβου
wgn = noise_factor.*(randn(numSamples,1)+1i*randn(numSamples,1));
bpam_noised = bpam+real(wgn);   % προσθήκη του πραγματικού μέρους του θορύβου
 
bit_error_rate = (bpam.*bpam_noised)<0; % ουσιαστικά πολλαπλασιάζουμε τα αθορύβητα και τα θορυβημένα σήματα
                                        % και λογίζουμε ως σφάλματα τα ετερόσημα αποτελέσματα
                                        % αποθηκεύοντές τα ως άσσους σε
                                        % έναν πίνακα
bit_error_rate = sum(bit_error_rate)/length(bpam);  % προσθέτοντες τους άσσους, μετράμε τα σφάλματα
theor_ber = qfunc(sqrt(2*Eb_N0_lin));
figure()
hold on;
plot(snr_db,bit_error_rate+eps,'LineWidth',2)
plot(snr_db,qfunc(sqrt(2*10.^(snr_db/10))),'LineWidth',2)
title('Θεωρητική Πιθανότητα Εσφαλμένου Ψηφίου B-PAM και Πειραματική')
hold off;
set(gca,'yscale','log')
grid on;
%--------------------------------------------------------------------------


%4ο ερώτημα
%--------------------------------------------------------------------------
AM = 319051;

%α
%--------------------------------------------------------------------------
%random_bits = randi([0 1], 36 ,1);
random_bits = [1 1 0 1 0 1 1 0 0 0 0 0 1 1 1 1 0 1 0 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 1 1];
reference_constellation = A/2.*[sqrt(2)+1i*sqrt(2), -sqrt(2)+1i*sqrt(2),-sqrt(2)-1i*sqrt(2),sqrt(2)-1i*sqrt(2)];

qpsk_mod = comm.QPSKModulator('PhaseOffset',pi/4,'BitInput',true,'SymbolMapping','Gray');
qpsk_modData = A*qpsk_mod(random_bits');
constDiagram1 = comm.ConstellationDiagram('Title', 'Διάγραμμα Αστερισμού QPSK' ...
    ,'ReferenceConstellation',reference_constellation ...
    ,'XLimits',[-A, A] ...
    ,'YLabel',[-A, A]);
constDiagram1(qpsk_modData)

%--------------------------------------------------------------------------

%β
%--------------------------------------------------------------------------
sum_of_AM_dig = mod(AM,1000);
while sum_of_AM_dig > 9
    sum_of_AM_dig = mod(sum_of_AM_dig,10) + fix(sum_of_AM_dig/10);
end
Tb = 0.25;
%Ebit = sum_of_AM_dig^2*Tb/2;

Eb_N0_db = [5 15];
snr_db = Eb_N0_db + 10*log10(2);
%snr_lin = 10.^(snr_db./10);
%N0_lin = Ebit ./ snr_lin;
%N0_db = 10*log10(N0_lin);
% samples = length(qpsk_modData);
% X = sqrt(N0_lin/4).*randn(samples,1);
% Y = sqrt(N0_lin/4).*randn(samples,1);
% Z = X + 1i*Y;
% qpsk_with_wgn = qpsk_modData+Z;

noise_channel_01 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)' ...
        ,'SNR',snr_db(1),'SignalPower', A^2);
noise_channel_02 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)' ...
        ,'SNR',snr_db(2),'SignalPower', A^2);
qpsk_with_wgn_01 = noise_channel_01(qpsk_modData);
qpsk_with_wgn_02 = noise_channel_02(qpsk_modData);

max_qpsk_wgn_x = [max(abs(real(qpsk_with_wgn_01))), max(abs(real(qpsk_with_wgn_02)))];  % υπολογισμών των περιθωριών
max_qpsk_wgn_y = [max(abs(imag(qpsk_with_wgn_01))), max(abs(imag(qpsk_with_wgn_02)))];  % των διαγραμμάτων αστερισμού

constDiagram2 = comm.ConstellationDiagram('Title','Διάγραμμα Αστερισμού QPSK Eb/No=5db' ...
    ,'ReferenceConstellation', reference_constellation ...
    ,'XLimits',sqrt(2)*[-max_qpsk_wgn_x(1), max_qpsk_wgn_x(1)] ...
    ,'YLabel',sqrt(2)*[-max_qpsk_wgn_y(1), max_qpsk_wgn_y(1)]);
constDiagram2(qpsk_with_wgn_01)

constDiagram3 = comm.ConstellationDiagram('Title','Διάγραμμα Αστερισμού QPSK Eb/No=15db' ...
    ,'ReferenceConstellation', reference_constellation ...
    ,'XLimits',sqrt(2)*[-max_qpsk_wgn_x(2), max_qpsk_wgn_x(2)] ...
    ,'YLabel',sqrt(2)*[-max_qpsk_wgn_y(2), max_qpsk_wgn_y(2)]);
constDiagram3(qpsk_with_wgn_02)

% scatterplot(qpsk_with_wgn(:,1))
% scatterplot(qpsk_with_wgn(:,2))
%--------------------------------------------------------------------------


%γ
%--------------------------------------------------------------------------
num_of_bits = 10^6;
random_bits = randi([0 1], num_of_bits,1);  % ένα εκατομμύριο τυχαία δυαδικά ψηφία

qpsk_mod = comm.QPSKModulator('BitInput', true);    
qpsk_demod = comm.QPSKDemodulator('BitOutput',true);
qpsk_modData = A*qpsk_mod(random_bits); % πλάτους Α
BER = zeros(16,1);

Ebit = sum_of_AM_dig^2*Tb/2;
Eb_N0_db = 0:15;
Eb_N0_lin = 10.^(Eb_N0_db./10);
snr_db = Eb_N0_db + 10*log10(2);
snr_lin = 10.^(snr_db./10);
N0_lin = Ebit ./ snr_lin;
N0_db = 10*log10(N0_lin);

index = 1;
for snr = snr_db
    noise_channel = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)' ...
        ,'SNR',snr,'SignalPower', A^2);
    qpsk_modDataNoised = noise_channel(qpsk_modData);
    qpsk_demodData = qpsk_demod(qpsk_modDataNoised);
    errorRate = comm.ErrorRate; % αντικείμενο για τον υπολογισμό σφάλματος
    errorStats = errorRate(random_bits,qpsk_demodData);
    BER(index) = errorStats(1); % το πρώτο στοιχείο του πίνακα errorStats δίνει το BER
    index = index + 1;
end
theor_ber_qpsk = qfunc(sqrt(2.*Eb_N0_lin));

figure;
hold on;
plot(Eb_N0_db,BER+eps,'*')
plot(Eb_N0_db,theor_ber,'LineWidth',2)
title('Θεωρητική Πιθανότητα Εσφαλμένου Ψηφίου QPSK και Πειραματική')
hold off;
set(gca,'yscale','log')
grid on;

bpsk_mod = comm.BPSKModulator;  % υλοποίηση του αντιστοίχου συστήματος BPSK διαμόρφωσης
bpsk_demod = comm.BPSKDemodulator;

bpsk_modData = A*bpsk_mod(random_bits);
index = 1;
snr_db_bpsk = snr_db - 10*log10(2);
for snr = snr_db_bpsk
    noise_channel = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)' ...
        ,'SNR',snr,'SignalPower', A^2);
    bpsk_modDataNoised = noise_channel(bpsk_modData);
    bpsk_demodData = bpsk_demod(bpsk_modDataNoised);
    errorRate = comm.ErrorRate;
    errorStats = errorRate(random_bits,bpsk_demodData);
    BER(index) = errorStats(1);
    index = index + 1;
end

theor_ber_bpsk = qfunc(sqrt(2.*Eb_N0_lin));

hold on;
plot(Eb_N0_db,BER+eps,'+')
plot(Eb_N0_db,theor_ber_bpsk,'LineWidth',2)
set(gca,'yscale','log')
grid on;
%--------------------------------------------------------------------------

%δ1,2
%--------------------------------------------------------------------------
sum_of_AM_dig = mod(AM,1000);
while sum_of_AM_dig > 9
    sum_of_AM_dig = mod(sum_of_AM_dig,10) + fix(sum_of_AM_dig/10);
end
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
[binary_vector,binary_string] = text2bin(ascii_text);   % συνάρτηση για την μετατροπή κειμένου σε ψηφία
                                                        % επιστρέφει τα
                                                        % ψηφία σε μορφή
                                                        % διανύσματος
                                                        % ακεραίων ή σε
                                                        % μορφή
                                                        % συμβολοσειράς
bits = 8;
Mmax = max(ascii_text);
Mmin = 0;
levels = 2^bits;
D = (Mmax-Mmin)/levels-1;
partition = linspace(Mmin+D/2,Mmax-D/2,levels-1);
codebook = linspace(Mmin,Mmax,levels);
[~ , quantized_text] = quantiz(ascii_text,partition,codebook); % η μεταβλητή index αντικαταστάθηκε από ~ διότι θα έμενε άχρηστη
figure;
plot(quantized_text)
xlabel('Χαρακτήρες');
ylabel('Επίπεδο Κβάντισης')
title('Κβαντισμένοι Χαρακτήρες');
%--------------------------------------------------------------------------

%δ3
%--------------------------------------------------------------------------
qpskModulator = comm.QPSKModulator('BitInput',true);
qpskModulatedData = qpskModulator(binary_vector');
constDiagram4 = comm.ConstellationDiagram('Title', 'Διάγραμμα Αστερισμού QPSK Κβαντισμένου Σήματος');
constDiagram4(qpskModulatedData);
%--------------------------------------------------------------------------

%δ4,5,6
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

%δ7
%--------------------------------------------------------------------------
text_recovered_1 = bin2text(qpskDemodulatedData1);  % συνάρτηση μετατροπής κειμένου από δυαδικά ψηφία
text_recovered_2 = bin2text(qpskDemodulatedData2);

fileID = fopen('recovered_text1.txt', 'w');
fwrite(fileID, text_recovered_1);
fclose(fileID);
fileID = fopen('recovered_text2.txt', 'w');
fwrite(fileID, text_recovered_2);
fclose(fileID);
%--------------------------------------------------------------------------


%5
%α
%--------------------------------------------------------------------------
if mod(sum_of_AM_dig,2)==0
    fileName = 'soundfile2_lab2.wav';
else
    fileName = 'soundfile1_lab2.wav';
end
[sound_sampled, Fs] = audioread(fileName);
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
plot(time_sound, sound_sampled);
title('Ηχητικό Σήμα');
xlabel('Χρόνος (sec)');
ylabel('Πλάτος');
axis tight;
figure;
plot(time_sound, sound_quantized)
title('Κβαντισμένο Ηχητικό Σήμα');
xlabel('Χρόνος (sec)');
ylabel('Στάθμες Κβάντισης');
axis tight;
%--------------------------------------------------------------------------

%γ
%--------------------------------------------------------------------------
sound_bits = dec2bin(index);    % μετατροπή των δεκαδικών σταθμών σε δυαδική αναπαράσταση
sound_bits = sound_bits';   % αναστροφή
sound_bits_string = sound_bits(:)'; % συνένωση στηλών και αναστροφή σε γραμμή
sound_bits_vector = sound_bits_string-48;   % τέχνασμα για την δημιούργια πίνακος ψηφίων τύπου int αντί τύπου char
qpskModulatorSound = comm.QPSKModulator('BitInput',true);
qpskModulatedSound = qpskModulatorSound(sound_bits_vector');
constDiagramSound = comm.ConstellationDiagram('Title' ...
    ,'Διάγραμμα Αστερισμού QPSK Κβαντισμένου Ηχητικού Σήματος ' ...
    ,'SamplesPerSymbol',2 ...
    ,'SymbolsToDisplaySource', 'Property' ...
    ,'SymbolsToDisplay', 50000);    % χρήση 50.000 ψηφίων προς αναπαράσταση
                                    % διότι είναι αδύνατον να
                                    % αναπαρασταθούν όλα τα 2.000.000
constDiagramSound(qpskModulatedSound)
%--------------------------------------------------------------------------

%δ
%--------------------------------------------------------------------------
qpskDemodulatorSound = comm.QPSKDemodulator('BitOutput',true);
constDiagramSound2 = comm.ConstellationDiagram('Title', ...
    'Διάγραμμα Αστερισμού QPSK Κβαντισμένου Ηχητικού Σήματος Es/N0 = 4db' ...
    ,'SamplesPerSymbol',2 ...
    ,'SymbolsToDisplaySource', 'Property' ...
    ,'SymbolsToDisplay', 10^5);
constDiagramSound3 = comm.ConstellationDiagram('Title', ...
    'Διάγραμμα Αστερισμού QPSK Κβαντισμένου Σήματος Es/N0 = 14db' ...
    ,'SamplesPerSymbol',2 ...
    ,'SymbolsToDisplaySource', 'Property' ...
    ,'SymbolsToDisplay', 10^5); % ομοίως, χρήση περιορισμένων αλλά αρκετών 
                                % ψηφίων προς αναπαράσταση

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
constDiagramSound2(qpskModulatedSoundNoised1);
constDiagramSound3(qpskModulatedSoundNoised2);

recovered_indexes_1 = bin2dec(char( ...
    reshape(qpskDemodulatedSound1,[bits, length(qpskDemodulatedSound1)/bits])'+48))';
recovered_indexes_2 = bin2dec(char( ...
    reshape(qpskDemodulatedSound2,[bits, length(qpskDemodulatedSound2)/bits])'+48))';

recover1_sound = 2*Mmax*(recovered_indexes_1/(2^bits-1)-0.5);
recover2_sound = 2*Mmax*(recovered_indexes_2/(2^bits-1)-0.5);

audiowrite('recovered_sound1.wav',recover1_sound,Fs);
audiowrite('recovered_sound2.wav',recover2_sound,Fs);
