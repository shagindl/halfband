BUF_SIZE = 2048;
freq = 1000;

buf = test_sine(BUF_SIZE, freq, 100*freq);

orig_sum = 0;
for j = 1:4:BUF_SIZE
    if(j > BUF_SIZE / 3)
        orig_sum = orig_sum + pow2(buf(j));
    end
end

plot(buf)

function buf = test_sine(n, F, Fs)
    for i = 1:n
        t = (1/Fs) * i;
        buf(i) = sin( 2*pi*F * t);
    end
end

