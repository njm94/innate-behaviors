function new_dat = get_on_time(data)
new_dat = zeros(size(data));
d = diff(data);
t_on = find(d>0)+ 1;

new_dat(t_on) = 1;

end