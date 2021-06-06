function[cell] = Ackley_function(~)
format long g

number_of_solution = 100;
generations = 100;
nofs = number_of_solution;
son = nofs * 2;
nos = number_of_solution/2;
sum_aptd = zeros(1,generations);
ipg = zeros(1,generations);
gp = zeros(1,generations);
cell_num = zeros(son,5);
cell_num = num2cell(cell_num);

[cell] =  poblation_worm_float_Ackley(nofs);
[cell] = objetive_function_float_Ackley(cell);

for k = 1:nofs
    cell_num(k,:) = cell(k,:);
end

for i = 1:generations
    for j = 1:nos
        [pred1,pred2] = sel_pad_best_float_Ackley(cell);
        [desc1,desc2] = two_point_crossover_float_Ackley(pred1,pred2);

        
        cell_num{j+nofs,1} = desc1;
        cell_num{j+nos*3,1} = desc2;
    end    
   
    [cell_num] = objetive_function_float_Ackley(cell_num);
    [cell_num] = biology_competition_float_Ackley(cell_num);
   
    b = mod(i,10);
    if b == 0
        disp('heuristic')
        [cell_num] = heuristic_mutation_float_Ackley(cell_num);
    end
   
    [sum_apt] = apt_for_generation_float_Ackley(cell_num);
    sum_aptd(i) =  sum_apt;
    ipg(i) = cell_num{1,4};
    gp(i) = i;
    cell = cell_num;

    
end
%disp(sum_aptd)
disp(cell{1,4})
disp(cell{1,3})

[z] =  plot3d(cell);


%plot(gp,ipg)
%xlabel('Generación')
%ylabel('f(x)')
%plot(gp,sum_aptd)
cell = cell_num;

end









%poblation_worm_float_Ackley.m
function [cell] =  poblation_worm_float_Ackley(nofs)
format long g

%number_of_posibles_solution = nofs;
%nofs = 100;
chromosome_bin = zeros(1,82);
%chromosome_dec = zeros(1,2);
%chromosome_np = zeros(1,2);
cell = zeros(100,5);
cell = num2cell(cell);


for i = 1:nofs
    
    n = 16;         % number bits for integer part of your number      
    m = 25;         % number bits for fraction part of your number
    
    k = rand(1);
    p = randi([0 30], 1, 1);
    a = k*p;

   % binary number
    %d2b = fix(rem(a*pow2(-(n-1):m),2)); 
   
    % the inverse transformation
    %b2d = d2b*pow2(n-1:-1:-m).';
    
    %generete a x1 in binary and float decimal
    k = rand(1);
    p = randi([0 30], 1, 1);
    a = k*p;
    
    x1bin = fix(rem(a*pow2(-(n-1):m),2)); 
    %x1dec = x1bin*pow2(n-1:-1:-m).';
    
    %generete a x2 in binary and float decimal
    k = rand(1);
    p = randi([0 30], 1, 1);
    a = k*p;
    
    x2bin = fix(rem(a*pow2(-(n-1):m),2)); 
    %x2dec = x2bin*pow2(n-1:-1:-m).';
    
    for k = 1: length(x1bin)
       chromosome_bin(k) = x1bin(k);
       chromosome_bin(k+41) = x2bin(k); 
    end
    
    %{
    chromosome_np(1) = 1;
    chromosome_np(2) = 1;

    if x1dec <= 15
        posneg = randi([0 1], 1, 1);
        chromosome_np(1) = posneg;
        if posneg == 0
            x1dec = x1dec * -1;
        end
    end    
     if x2dec <= 15
        posneg = randi([0 1], 1, 1);
        chromosome_np(2) = posneg;
        if posneg == 0
            x2dec = x2dec * -1;
        end
    end 
       
    chromosome_dec(1) = x1dec;
    chromosome_dec(2) = x2dec;
    %}
    cell{i,1} = chromosome_bin;
    %cell{i,2} = chromosome_np;
    %cell{i,3} = chromosome_dec;

end

end









%objetive_function_float_Ackley.m
function[cell] = objetive_function_float_Ackley(cell)
format long g

a = 20;
b = .2;
d = 2;
c = 2*pi;
n = 16;         % number bits for integer part of your number      
m = 25;         % number bits for fraction part of your number
chromosome_dec = zeros(1,2);
chromosome_np = zeros(1,2);
x1bin = zeros(1,41);
x2bin = zeros(1,41);
len = length(cell);

for i = 1:len
    chromosome_bin = cell{i,1};
    
    for k = 1:41
        x1bin(k) = chromosome_bin(k);
        x2bin(k) = chromosome_bin(41+k);
    end
    
    x1dec = x1bin*pow2(n-1:-1:-m).';
    x2dec = x2bin*pow2(n-1:-1:-m).';
    
    %if d_x1 == 31
    %    d_x1 = 30;
    %    chromosome(1,2) = 0;
        
    %end
    %if d_x2 == 31
    %    d_x2 = 30;
    %    chromosome(1,8) = 0;
    %end
    
    
    
    chromosome_np(1) = 1;
    chromosome_np(2) = 1;
    if x1dec <= 15
        posneg = randi([0 1], 1, 1);
        chromosome_np(1) = posneg;
        if posneg == 0
            x1dec = x1dec * -1;
        end
    end    
     if x2dec <= 15
        posneg = randi([0 1], 1, 1);
        chromosome_np(2) = posneg;
        if posneg == 0
            x2dec = x2dec * -1;
        end
     end
     
     
    chromosome_dec(1) = x1dec;
    chromosome_dec(2) = x2dec;
    
    x1 = x1dec;
    x2 = x2dec;
    
    fx = a + exp(1) - a*exp(-b*sqrt( (1/d)* ((x1.^2) + (x2.^2)))) - exp((1/d) * (cos(c*x1) + cos(c*x2)) );

    fitness = 1/fx;
    
    cell{i,2} = chromosome_np;
    cell{i,3} = chromosome_dec;
    cell{i,4} = fx;
    cell{i,5} = fitness;
    
end

[cell] = ord_insertion_float_Ackley(cell);


end









%two_point_crossover_float_Ackley.m
function[desc1,desc2] = two_point_crossover_float_Ackley(pred1,pred2)

temp = zeros(1,8);

nir = randi([12,33],1,1);
for i = 1:8
    temp(i) = pred1(nir+i);
    pred1(nir+i) = pred2(nir+i);
    pred2(nir+i) = temp(i);
     
end

desc1 = pred1;
desc2 = pred2;

end









ord_insertion_float_Ackley.m
function[celda_ord] = ord_insertion_float_Ackley(celda_num)

for ls = 1:length(celda_num)
    d = ls;
    while((d > 1) && (celda_num{d,5}) > (celda_num{d-1,5}))
        % aptitud.
        var_temp = celda_num{d,5};
        celda_num{d,5} = celda_num{d-1,5};
        celda_num{d-1,5} = var_temp;
        
        var_temp1 = celda_num{d,4};
        celda_num{d,4} = celda_num{d-1,4};
        celda_num{d-1,4} = var_temp1;
        
        %distance. 
        var_temp1 = celda_num{d,3};
        celda_num{d,3} = celda_num{d-1,3};
        celda_num{d-1,3} = var_temp1;
        %gene. 
        var_temp2 = celda_num{d,2};
        celda_num{d,2} = celda_num{d-1,2};
        celda_num{d-1,2} = var_temp2;
        
        var_temp2 = celda_num{d,1};
        celda_num{d,1} = celda_num{d-1,1};
        celda_num{d-1,1} = var_temp2;
        
        d = d-1;
        
    end
end
celda_ord= celda_num;
end









%sel_pad_best_float_Ackley.m
function [padre1,padre2] =  sel_pad_best_float_Ackley(Ce_cmp)    
%global selection_metodo 
padre1 = zeros(1,82);
padre2 = zeros(1,82);
for rl = 1:2
    if rl ==1
       [index] = my_own_RWS_best_float_Ackley(Ce_cmp);
       padre1 = Ce_cmp{index,1};
    end
    if rl ==2
        [index] = my_own_RWS_best_float_Ackley(Ce_cmp);
        padre2 = Ce_cmp{index,1};
    end
end
 
end









%my_own_RWS_best_float_Ackley.m
function [index] = my_own_RWS_best_float_Ackley(Ce_cmp)
% generamos la probabilidad de que sean seleccionados, esta aumenta 
%dependiendo de su fitness
%creamos un valor prioridad.
len_v = length(Ce_cmp);
vec_prio = zeros(len_v,1);
for le = 1:len_v
    vec_prio(le,1) = Ce_cmp{le,5};
end
%[1] = previous_probability + (fitness / sum_of_fitness) = 0.0 + (1 / 10) = 0.1
%previous_probability = 0.1
vec_prio = flipud(vec_prio); %invertimos los valores 
sum_vec_prio = sum(vec_prio);
prob_selec = vec_prio/sum_vec_prio;   %Generamos todo la matriz con los resultados (fitness / sum_of_fitness)
%Generamos la probabilidad de ser seleccionado la suma de pre_pro + (fit/sum)
proba = zeros(1,len_v);
prev_proba = 0;
for km = 1:len_v
    proba(km) = prev_proba + prob_selec(km,1);
    prev_proba = proba(km);
end
%Escogemos al asar el numero en index que necesitamos.
%xbp = flipud(proba);
num_rand = rand;
%disp(num_rand)
for ksr = 1:len_v
    if num_rand < proba(ksr)
        index = ksr;
        return
    end
end

end









%biology_competition_float_Ackley.m
function[more_stronge] = biology_competition_float_Ackley(cell_pred_desc)
format long g

[cell_out] =  delete_repeated_float_Ackley(cell_pred_desc);

len_bc = length(cell_out);
len_bc_two = len_bc/2;
cell_str = zeros(len_bc_two,5);
cell_str = num2cell(cell_str);

for str = 1:len_bc_two
    cell_str(str,:) = cell_out(str,:);
end    


more_stronge = cell_str;



end









%delete_repeated_float_Ackley.m
function[cell_nor] =  delete_repeated_float_Ackley(cell_act)
format long g


cookie_loc = zeros(18,1);

len_dt = length(cell_act);
%disp(len_dt)
for dtl = 1:len_dt
    for dtla = 1:len_dt
        exist_cookie = any(cookie_loc(:) == dtla);
        if cell_act{dtl,5} == cell_act{dtla,5} & dtl ~= dtla & exist_cookie == 0
            cookie_loc(dtl) = dtl;
            gen_loc_rep = cell_act{dtla,1};
            %disp(gen_loc_rep)
            [gen_mut] = scramble_met_per_one_Ackley(gen_loc_rep);
            cell_act{dtla,1} = gen_mut;    
        end
    end  
     %disp(dtl)

end
[cell_act] = objetive_function_float_Ackley(cell_act);

%[cell_act] = make_dist_apt(cell_act,dista_pid);
%disp(cell_act)
cell_nor = cell_act;
end









%scramble_met_per_one_float_Ackley.m
function[anser] = scramble_met_per_one_float_Ackley(gen_loc_rep)

%gen_loc_rep = randi([0,1],1,82);
%disp(gen_loc_rep)

nsm = randi([12,74],1,1);
nsn = nsm + 8;
nso = 82 - nsn;
mat_rep_one = zeros(1,8);
mat_rep_two = zeros(1,74);

for i = 1:8
    mat_rep_one(i) = gen_loc_rep(nsm+i);
end

for j = 1:nsm
    mat_rep_two(j) = gen_loc_rep(j);
end    

for k = 1:nso 
    mat_rep_two(8+k) = gen_loc_rep(nsn + k);
end



nsm = randi([12,74],1,1);
nsn = nsm + 8;
nso = 82 - nsn;

for ii = 1:nsm
    gen_loc_rep(ii) = mat_rep_two(ii);
end
for jj = 1:8
    gen_loc_rep(nsm+jj) = mat_rep_one(jj);
end
for kk = 1:nso
    gen_loc_rep(nsn + kk) = mat_rep_two(nsm+kk);
end    

anser = gen_loc_rep;


end









%heuristic_mutation_float_Ackley.m
function[cell_muted] = heuristic_mutation_float_Ackley(cell_for_mut)
%in this function acept the cell and back 
% a cell with the 10 with mutation.

len_hm = length(cell_for_mut);
len_por = len_hm/20;
len_ph = len_hm - len_por;

for hma = 1:len_por
    one_man = cell_for_mut{len_ph+hma,1};
    [one_man_mut] = permu_loc_float_Ackley(one_man);
    cell_for_mut(len_ph+hma,:) = one_man_mut(1,:);
end

[cell_for_mut] = ord_insertion_float_Ackley(cell_for_mut);

%[cell_for_mut] = ord_insertion_Ackley(cell_for_mut);
%[cell_for_mut] = ord_insertion(cell_for_mut);

cell_muted = cell_for_mut;


end









%permu_loc_float_Ackley.m
function[one_man_mut] = permu_loc_float_Ackley(one_man)

%one_man = randi([0,1],1,82);
one_man_per = zeros(1,6);
rni = randi([12,76],1,1);

for i = 1:6
    one_man_per(i) = one_man(rni+i);
end

one_man_per = perms(one_man_per);
cell = zeros(length(one_man_per),5);
cell = num2cell(cell);

for j = 1:length(one_man_per)
    for k = 1:6
        one_man(rni+k) = one_man_per(j,k);
    end    
    cell{j,1} = one_man;
end
 
[cell] = objetive_function_float_Ackley(cell);
one_man_mut = cell(1,:);


end  









%apt_for_generation_float_Ackley.m
function[sum_apt] = apt_for_generation_float_Ackley(cell_act)
len_cel = length(cell_act); 
mat_apt = zeros(1,len_cel);
for ju = 1:len_cel
    mat_apt(ju) = cell_act{ju,5};
end
sum_apt = sum(mat_apt);

end









%plot3d.m
function[z] =  plot3d(cell)
format long g

len = length(cell);
x1 = zeros(1,len);
x2 = zeros(1,len);
fx = zeros(1,len);

for i = 1:len
    x1x2 = cell{i,3};
    x1(i) = x1x2(1);
    x2(i) = x1x2(2);
    fx(i) = cell{i,4}; 
end    
set(fsurf(@(x,y) ackleyfcn([x,y]),[-15 30 -15 30]),'AdaptiveMeshDensity',0,'MeshDensity',60)

hold on

plot3(x1,x2,fx,'*')

z= 1;
end








%ackleyfcn.m
function z = ackleyfcn(xx)
format long g

% Ackley's function
% Search domain: [-15,30]
% Global minimum: f(x) = 0 | x = (0,...,0)
%disp(xx)

a = 20;
b = .2;
d = 2;
c = 2*pi;

xx = max(-32,min(32,xx));

z = a + exp(1) - a*exp(-b*sqrt(1/d*sum(xx.^2,2))) - exp(1/d*sum(cos(c*xx),2));
%z = a + exp(1)- a*exp(-b*sqrt( (1/d)* ((xx.^2) + (xx.^2)))) - exp((1/d) * (cos(c*xx) + cos(c - xx)) );

hold on


end
