
%------------------Q1------------------------------------------------

%comment: part a is coming in the functions part!

%------------------------------part b-----------------------------------
n = [2,3,4,5]; 
theta_sample = linspace(0,pi,41); 
for j = 1:length(n)
    rand_theta = linspace(0,pi,n(j)); 
    for i = 1:length(theta_sample)
        app_phi(j,i) = Lagrange_Int(theta_sample(i),rand_theta,'B');
    end
end
for i = 1:length(theta_sample)  
   re_phi(i) = potential(theta_sample(i),'B');
end
graph_b(re_phi,theta_sample,app_phi) %lets print a graph
%------------------------------part c-----------------------------------
list_err = 2:2:20;
for i = 1:length(list_err)
    relerrc(i) = RelError(rand_theta,list_err(i),'B');
end
graph_c(list_err,relerrc)  %lets print a graph
%------------------------------part d-----------------------------------
n = [3,7,11,15]; 
theta_sample_d = linspace(0,pi,41); 
for j = 1:length(n)
    rand_theta = linspace(0,pi,n(j)); 
    for i = 1:length(theta_sample_d)
        app_phi_d(j,i) = Lagrange_Int(theta_sample_d(i),rand_theta,'D');
    end
end
for i = 1:length(theta_sample_d)  
   re_phi_d(i) = potential(theta_sample_d(i),'D');
end
graph_d_1(theta_sample,re_phi_d,app_phi_d) %lets print a graph
for i = 1:length(list_err)
    relerrd(i) = RelError(rand_theta,list_err(i),'D');
end
graph_d_2(list_err,relerrd) %lets print a graph
%------------------------------part e-----------------------------------
n = [3,7,11,15]; 
theta_sample_e = linspace(0,pi,41); 
re_phi_e = 0;
for j = 1:length(n) 
    for i = 1:length(theta_sample_e)
        app_phi_e(j,i) = Lagrange_Int(theta_sample_e(i),Chebyshev(n(j)-1,0,pi),'D');
    end
end
for i = 1:length(theta_sample_e)  
   re_phi_e(i) = potential(theta_sample_e(i),'D');
end

graph_e_1(theta_sample_e,re_phi_e,app_phi_e)

for i = 1:length(list_err)
    relerre(i) = RelError_e(theta_sample_e,list_err(i),'D');
end
graph_e_2(list_err,relerre,relerrd)






%--------------------------------------part a - Q1---------------------------
function [LN] = Lagrange_Int(x,teta,part)
    sum = 0; %initialization
    for i = 1:length(teta)
        mone = 1;%initialization
        mehane = 1;%initialization
        for j = 1:length(teta)
            if (j~=i)
                mone = mone * (x-teta(j)); 
                mehane = mehane * (teta(i)-teta(j));  
            end
        end
        sum = sum + potential(teta(i),part) * (mone/mehane);  %formula
    end
    LN = sum;
end



function [val] = radius(x,siman,r)  
    delta = 5 * 10^(-3);%initialization
    if siman == '+'
        val = sqrt((r*cos(x))^2+(r*sin(x) - delta/2)^2); %formula  
    elseif siman == '-'
        val = sqrt((r*cos(x))^2+(r*sin(x) + delta/2)^2);%formula
    end
end    


function [phi] = potential(x,part) 
    if part == 'B' %for part b
        r = 0.05;%initialization
    elseif part == 'D' %for part d
        r = 4*10^(-3);%initialization
    end
    q_minus = -sum([2 0 7 8 8 1 0 0]) - sum([3 1 8 8 8 0 0]);  %formula
    q_plus = sum([0 7 8 8 1 0 0 4]) + sum([1 8 8 8 0 0 5]);  %formula
    phi = (q_plus/(4*pi*radius(x,'+',r))) + (q_minus/(4*pi*radius(x,'-',r))); %formula
end



%------------------functions-------------------------------
function graph_b(re_phi,theta_sample,app_phi)
    figure(1)
    p = plot(theta_sample, re_phi, theta_sample, app_phi,'-+');
    p(1).LineWidth = 2;
    p(1).Color = 'b';
    legend('Real Phi','Approximate Phi2: 2 samples','Approximate Phi3: 3 samples','Approximate Phi4: 4 samples','Approximate Phi5: 5 samples','Location','northeast')
    title('Q1 Section B- \phi(\teta)')
    xlabel("\teta")
    ylabel("\phi(\teta)")
    grid on
end


function [relerr] = RelError(rand_theta,n,part)
    theta_sample = linspace(0,pi,n);
    mone = 0;%initialization
    mehane = 0;%initialization
    for i = 1:length(rand_theta)
        mone = mone + (Lagrange_Int(rand_theta(i),theta_sample,part) - potential(rand_theta(i),part))^2;
        mehane = mehane + (potential(rand_theta(i),part))^2;
    end
    relerr = sqrt(mone / mehane);
end

function [relerr] = RelError_e(rand_theta,n,part)
    mone = 0;%initialization
    mehane = 0;%initialization
    for i = 1:length(rand_theta)
        mone = mone + (Lagrange_Int(rand_theta(i), Chebyshev(n-1,0,pi),part) - potential(rand_theta(i),part))^2;
        mehane = mehane + (potential(rand_theta(i),part))^2;
    end
    relerr = sqrt(mone / mehane);
end

function graph_c(list_err,relerrc)
    figure(2)
    semilogy(list_err,relerrc,'-o')
    title('Q1 part c-Rel Error')
    xlabel("n+1")
    ylabel("Rel Error")
    grid on
end



function graph_d_1(theta_sample,re_phi_d,app_phi_d)
    figure(3)
    p = plot(theta_sample, re_phi_d, theta_sample, app_phi_d,'-+');
    p(1).LineWidth = 2;
    p(1).Color = 'b';
    legend('Real Phi','Approximate Phi3: 3 samples','Approximate Phi7: 7 samples','Approximate phi11: 11 samples','Approximate Phi15: 15 samples','Location','northeast')
    title('Q1 part d - \phi(\teta) Changed Radius')
    xlabel("\teta")
    ylabel("\phi(\teta)")
    grid on
end


function graph_d_2(list_err,relerrd)
    figure(4)
    semilogy(list_err,relerrd,'-+')
    title('Q1 part d  - Rel Error ner r')
    xlabel("n+1")
    ylabel("Rel Error")
    grid on
end

function [roots] = Chebyshev(n, start, End)
    n = n + 1;%initialization
    teta = [];%initialization
    for i = 1:n
        x = cos(pi*(2*i-1)/(2*n));%formula
        t = ((End-start)*x+End+start)/2;%formula
        teta = [teta t];
    end
    roots = teta;
end

function graph_e_1(theta_sample_e,re_phi_e,app_phi_e)
    figure(5)
    p = plot(theta_sample_e, re_phi_e, theta_sample_e, app_phi_e,'-*');
    p(1).LineWidth = 2;
    p(1).Color = 'b';
    legend('Real Phi','Approximate Phi3: 3 samples','Approximate Phi7: 7 samples','Approximate phi11: 11 samples','Approximate Phi15: 15 samples','Location','northeast')
    title('Q1 part e- \phi(\teta) Chebyshev')
    xlabel("\theta")
    ylabel("\phi(\theta)")
    grid on
end

function graph_e_2(list_err,relerre,relerrd)
    figure(6)
    semilogy(list_err,relerre,'-+',list_err,relerrd,'-+')
    legend('Rel error e, chebyshev','Rel error d: without chebyshev')
    title('Q1 part e Rel Error')
    xlabel("n+1")
    ylabel("Relative Error")
    grid on
end







