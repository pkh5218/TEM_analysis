function [eumpsi, tempsi] = compute_TEM_ERA5(v,dthdp,vt_ave,p,cosphi)

g = 9.8;
a = 6371000;
v_ave = squeeze(mean(v,[1,4]));
v_star = v_ave - compute_dfdp(vt_ave./dthdp,p);
dp = compute_df(p);
dv_EUM = 2*pi*a.*cosphi/g.*v_ave;
dv_TEM = 2*pi*a.*cosphi/g.*v_star;
eumpsi = zeros(size(v,2),size(v,3));
tempsi = zeros(size(v,2),size(v,3));
for j = 2:size(v,3)
    eumpsi(:,j) = eumpsi(:,j-1)+dv_EUM(:,j)*dp(j);
    tempsi(:,j) = tempsi(:,j-1)+dv_TEM(:,j)*dp(j);
end


