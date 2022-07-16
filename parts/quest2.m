%2ο ερώτημα
y = trig1;  %   αποθήκευση των δειγμάτων σε μεταβλητή σύμφωνη με την εκφώνηση της ασκήσεως
if mod(fm/1000,2) == 0
    num_of_bits = 4;
else
    num_of_bits = 5;
end
levels = 2^num_of_bits;
Mmax = max(y);
D = 2*Mmax/levels;
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

error = y - quants;
error1 = y(1:10) - quants(1:10);
error2 = y(1:20) - quants(1:20);

SNR = ( quants.^2 * ones(size(y,2),1) ) / ( error.^2 * ones(size(y,2),1) );
SNR1 = ( quants(1:10).^2 * ones(10,1) ) / ( error1.^2 * ones(10,1) );
SNR2 = ( quants(1:20).^2 * ones(20,1) ) / ( error2.^2 * ones(20,1) );

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

function [arr] = graycode(nbits)
    arr = strings(1,2^nbits);
    arr(1:2) = ["0" "1"];
    i = 2;
    while(i<2^nbits)
        for j = i:-1:1
            arr(2*i-j+1) = arr(j);
        end
        for j = 1:i 
            arr(j) = "0" + arr(j);
        end
        for j = i+1:2*i
            arr(j) = "1" + arr(j);
        end
        i = i * 2;
    end
end