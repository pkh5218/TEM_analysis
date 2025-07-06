function var_log = plot_figure_zmlog(lat,lev,var,mymap,coeff,xlim_value)
var_p = var/1e9;
var_p(find(var<0))=0;

for j = 1:size(var_p,1)
    for k = 1:size(var_p,2)
        var_p_log(j,k) = log2(var_p(j,k));
        if var_p_log(j,k)<=1
            var_p_log(j,k) = var_p(j,k)/2;
        end
    end
end
var_p_log(isnan(var_p_log==1))=0;

var_n = -var/1e9;
var_n(find(var>0))=0;

for j = 1:size(var_n,1)
    for k = 1:size(var_n,2)
        var_n_log(j,k) = log2(var_n(j,k));
        if var_n_log(j,k)<=1
            var_n_log(j,k) = var_n(j,k)/2;
        end
    end
end
var_n_log(isnan(var_p_log==1))=0;
var_n_log = -var_n_log;
var_n(isnan(var_n==1))=0;

var_log = var_p_log+var_n_log;

plot_figure_zm(lat,lev,var_log,mymap,coeff,xlim_value)