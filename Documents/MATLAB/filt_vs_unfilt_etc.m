f1 = 3.215e-9;freq = 3.2150e-8;Rp = 1 ;Rs = 5 ;falloff = .3e-8 ;
f2=3.215e-8;
nyq = freq/2.0;
h_Wp = f1 / nyq;
h_Ws = (f1+falloff)/nyq;
[high_n,high_Wn] = buttord(h_Wp,h_Ws,Rp,Rs);
[high_b,high_a] = butter(high_n,high_Wn,'high');

hpf_am_n34 = filtfilt(high_b,high_a,double(n34_ind));

n34_mv = movingvar(n34_ind,30);

hpf_mv = movingvar(hpf_am_n34,30);

for i=1:499-30
    n34_rv(i) = var(n34_ind(i:i+29));
    hpf_rv(i) = var(hpf_am_n34(i:i+29));
end

plot(n34_mv(15:end),'k'); hold on;
plot(hpf_mv(15:end),'r');
plot(n34_rv,'b');
plot(hpf_rv,'g'); hold off;

legend('unfiltered+biased','filtered+unbiased','unfiltered+unbiased','filtered+biased')