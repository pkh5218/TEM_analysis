function S2 = compute_S2(T_ave,p)
N_p = length(p);
Rd = 287;
cp = 1004;
S2 = T_ave;
dTdp = compute_dfdp(T_ave,p);
for j = 1:N_p
    S2(:,j) = -(Rd/p(j))*(dTdp(:,j)-(Rd/cp)*T_ave(:,j)/p(j));
end
end