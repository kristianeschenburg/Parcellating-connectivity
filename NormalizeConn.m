function D = NormalizeConn(D)
D = cast(D, 'double');

off_diags = true(size(D));
for i = 1:size(D,1)
    off_diags(i,i) = false;
end
off_diags = off_diags(:);

D = D - mean(D(off_diags));
D = D./std(D(off_diags));

D = cast(D, 'single');
disp(num2str(mean(D(off_diags)))); disp(num2str(std(D(off_diags))));
end