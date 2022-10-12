close all; clear all; clc;


# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
To = ;

a0 = ;

terms = 1000;

tbottom = ;
ttop = ;
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


t = [tbottom:(ttop - tbottom) / 1000:ttop];
w0 = (2 * pi) / (To);
fx = a0 + (0 .* t);
tx = a0 + (0 .* t);
for n = 1:terms

  w = n * w0;

  
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  an = ;

  bn = ;
  
  tf = ;
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  fx += an .* cos(w .* t) + bn .* sin(w .* t);
  
  tx += abs(tf) .* (an .* cos(w .* t + angle(tf)) + bn .* sin(w .* t + angle(tf)));

end


subplot(2,1,1);
plot(t, fx, '-r');
xlim([tbottom ttop]);
ylim([(min(fx) - abs(min(fx) * 0.2)) (max(fx) + abs(max(fx) * 0.2))]);
grid minor on;
xlabel("Time");
ylabel("Amplitude");


subplot(2,1,2);
plot(t, tx, '-b');
xlim([tbottom ttop]);
ylim([(min(tx) - abs(min(tx) * 0.2)) (max(tx) + abs(max(tx) * 0.2))]);
grid minor on;
xlabel("Time");
ylabel("Transfer Amplitude");

