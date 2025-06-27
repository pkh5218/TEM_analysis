function var_3days = compute_var_3days(var,m,day,N_period,...
    p_ini, p_start, p_end, N_new_levels)

[N_x,~,~,N_day] = size(var);

if m == 1
   dd = 0;
else
   dd = 1;
end

if m == 1 && day == 1
    var_back = [];
    for j = 1:N_x
        var_mid(j,:,:) = compute_int(squeeze(var(j,:,:,day+dd)),p_ini, p_start, p_end, N_new_levels);
    end
    for j = 1:N_x
        var_for(j,:,:) = compute_int(squeeze(var(j,:,:,day+dd+1)),p_ini, p_start, p_end, N_new_levels);
    end
else
if m == N_period && day == N_day
    var_for = [];
    for j = 1:N_x
        var_mid(j,:,:) = compute_int(squeeze(var(j,:,:,day+dd)),p_ini, p_start, p_end, N_new_levels);
    end
    for j = 1:N_x
        var_back(j,:,:) = compute_int(squeeze(var(j,:,:,day+dd-1)),p_ini, p_start, p_end, N_new_levels);
    end
else
for j = 1:N_x
    var_mid(j,:,:) = compute_int(squeeze(var(j,:,:,day+dd)),p_ini, p_start, p_end, N_new_levels);
end
for j = 1:N_x
    var_for(j,:,:) = compute_int(squeeze(var(j,:,:,day+dd+1)),p_ini, p_start, p_end, N_new_levels);
end
for j = 1:N_x
    var_back(j,:,:) = compute_int(squeeze(var(j,:,:,day+dd-1)),p_ini, p_start, p_end, N_new_levels);
end
end
end

var_3days = cat(4,var_back,var_mid,var_for);
