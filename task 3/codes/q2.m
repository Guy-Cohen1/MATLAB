%-------------------------------Q2---------------------------------------
%--------------------------part b------------------------------------
n = [2,3,4];
theta_sample = linspace(0,pi,41); 
for j = 1:length(n)
    rand_theta = linspace(0,pi,n(j));
    for i = 1:length(rand_theta)
        y1(i) = potential2(rand_theta(i),'Q2');
    end
    [a(j),b(j),c(j)] = find_mekadmim(rand_theta,y1);
    least_squeres_b(j,:) = a(j)+b(j)*sin(theta_sample)+c(j)*cos(theta_sample);
end
for i = 1:length(theta_sample)
    real_phi_b(i) = potential2(theta_sample(i),'Q2');
end
graph_b(theta_sample,real_phi_b,least_squeres_b)
%--------------------------part c------------------------------------
rand_theta = linspace(0,pi,41);
r0 = 10;
r = r0*2.^(0:-1:-8);
n = 4; 
theta_sample = linspace(0,pi,n);
for j = 1:length(r)
    for i = 1:length(theta_sample)
        y1(i) = potential2_c(theta_sample(i),r(j));
    end
    [a(j),b(j),c(j)] = find_mekadmim(theta_sample, y1);
    least_squeres_c(j,:) = a(j)+b(j)*sin(rand_theta)+c(j)*cos(rand_theta);
    for i = 1:length(rand_theta)
        real_phi_c(j,i) = potential2_c(rand_theta(i),r(j));
    end
end
relerrc = sqrt(sum((least_squeres_c-real_phi_c).^2,2))./sqrt(sum(real_phi_c.^2,2));
graph_c(r,relerrc)

%--------------------------part d------------------------------------
n = 2.^(2:18);
rand_theta = linspace(0,pi,41);
for j = 1:length(n)
    theta_sample = linspace(0,pi,n(j)); 
    for i = 1:length(theta_sample)
        y(i) = potential2(theta_sample(i),'Q2');
    end
    yerr = (1+(rand(1,n(j))-0.5)*10^-1).*y;
    [a(j),b(j),c(j)] = find_mekadmim(theta_sample, yerr);
    least_squeres_d(j,:) = a(j)+b(j)*sin(rand_theta)+c(j)*cos(rand_theta);
end
for i = 1:length(rand_theta)
    real_phi_d(i) = potential2(rand_theta(i),'Q2');
end
rel_err_M = sqrt(sum((least_squeres_d-real_phi_d).^2,2))./sqrt(sum(real_phi_d.^2,2));
graph_d1(n,rel_err_M)

%--------------------------part d, different delta------------------------
n = 2.^(2:18);
rand_theta = linspace(0,pi,41);
for j = 1:length(n)
    theta_sample = linspace(0,pi,n(j)); 
    for i = 1:length(theta_sample)
        y_2(i) = potential2(theta_sample(i),'Q2');
    end
    y2err = (1+(rand(1,n(j))-0.5)*10^-4).*y_2;
    [a(j),b(j),c(j)] = find_mekadmim(theta_sample, y2err);
    least_squeres_d(j,:) = a(j)+b(j)*sin(rand_theta)+c(j)*cos(rand_theta);
end
for i = 1:length(rand_theta)
    real_phi_d(i) = potential2(rand_theta(i),'Q2');
end
rel2_err_M = sqrt(sum((least_squeres_d-real_phi_d).^2,2))./sqrt(sum(real_phi_d.^2,2));
graph_d2(n,rel2_err_M)

%---------------------------functions-----------------------------
function [phi] = potential2(x,part) 
    if part == 'A'
        r = 0.05;
    elseif part == 'C'
        r = 4*10^(-3);
    elseif part == 'Q2'
        r = 0.1;
    end
    q_minus = -sum([2 0 7 8 8 1 0 0]) - sum([3 1 8 8 8 0 0]);  %formula
    q_plus = sum([0 7 8 8 1 0 0 4]) + sum([1 8 8 8 0 0 5]);  %formula
    phi = (q_plus/(4*pi*radius(x,'+',r))) + (q_minus/(4*pi*radius(x,'-',r)));
end


function [a,b,c] = find_mekadmim(teta,y)
    f0 = ones(length(teta),1);
    f1 = sin(teta');
    f2 = cos(teta');
    F = [f0 f1 f2];
    vec = (inv(F'*F))*F'*y';
    a = vec(1);
    b = vec(2);
    c = vec(3);
end


function graph_b(theta_sample,real_phi_b,least_squeres_b)
    figure(7)
    p = plot(theta_sample, real_phi_b, theta_sample, least_squeres_b,'-*');
    p(1).LineWidth = 2;
    p(1).Color = 'k';
    title("Q2 partb- \phi(\teta) LS")
    xlabel("\teta")
    ylabel("\phi(\theta)")
    legend("real", "2","3","4")
    grid on
end

function [val] = radius(x,siman,r)  
    delta = 5 * 10^(-3);%initialization
    if siman == '+'
        val = sqrt((r*cos(x))^2+(r*sin(x) - delta/2)^2); %formula  
    elseif siman == '-'
        val = sqrt((r*cos(x))^2+(r*sin(x) + delta/2)^2);%formula
    end
end    

function [phi] = potential2_c(x,r) 
    q_minus = -sum([2 0 7 8 8 1 0 0]) - sum([3 1 8 8 8 0 0]);  %formula
    q_plus = sum([0 7 8 8 1 0 0 4]) + sum([1 8 8 8 0 0 5]);  %formula
    phi = (q_plus/(4*pi*radius(x,'+',r))) + (q_minus/(4*pi*radius(x,'-',r)));
end


function graph_c(r,relerrc)
    figure(8)
    loglog(r', relerrc, "-o")
    title("Q2 part c Rel Error")
    xlabel("Radius")
    ylabel("relative error")
    grid on
end


function graph_d1(n,rel_err_M)
    figure(9)
    loglog(n', rel_err_M, '-+')
    title("Q2 part d-Relative Error")
    xlabel("n+1")
    ylabel("rel error")
    grid on
end


function graph_d2(n,rel2_err_M)
    figure(10)
    loglog(n', rel2_err_M, '-+')
    title("Q2 part D-Rel err with different delta")
    xlabel("n+1")
    ylabel("rel error")
    grid on
end



