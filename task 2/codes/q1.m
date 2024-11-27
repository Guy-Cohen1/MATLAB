a = 1;
b = 5;
I1 = 318888005; %initialization 
I2 = 207881004; %initialization 
tollerence = 10^(-12); %initialization 
s = 3^(1/4); 
x0 = a + (b-a)*(I1/(I1+I2)); %initialization 
[X_n_q1,error_list,Xn_minus_xn_minus_one_list] = Newton_Raphson(x0,tollerence,s); %lets run NR
graph_q1(error_list); %print graph
num_of_iter = [1:length(error_list)]'; 
X_n = X_n_q1(:,1:length(X_n_q1)-1)'; 
difference = Xn_minus_xn_minus_one_list'; 
Err = error_list'; 
Table = table(num_of_iter,X_n, difference, Err); 
disp(Table); 
fprintf('It took %d iteration to converge\n',length(num_of_iter));





%------------------------------functions-----------------------------------
function [X_n,error_list,Xn_minus_xn_minus_one_list] = Newton_Raphson(x0,tollerence,s)
    X_n = zeros;%initialization 
    error_list = zeros; %initialization 
    Xn_minus_xn_minus_one_list = zeros; %initialization 
    X_n(1) = x0; 
    X_n(2) = x0 - h_q1(x0);
    iter = 2;
    difference = abs(X_n(iter)-X_n(iter-1)); %|xn-xn-1|
    Xn_minus_xn_minus_one_list(1) = difference;
    error_list(1) = abs(X_n(1) - s); %|xn - s| 
    while difference >= tollerence    %10^-12                 
        X_n(iter+1) = X_n(iter) - h_q1(X_n(iter)); %formula
        error_list(iter) = abs(X_n(iter) - s);     
        iter=iter+1;    
        difference = abs(X_n(iter)-X_n(iter-1));  
        Xn_minus_xn_minus_one_list(iter-1) = difference;       
    end
end

function f = func_q1(x) % f(x)
        f = x^4-3;
end

function f_derive = deriv_func_q1(x) % f'(x) 
        f_derive = 4*x^3;    
end


function h = h_q1(x)  %  f(x) / f'(x)
        h = func_q1(x) / deriv_func_q1(x);     
end


function graph_q1(error_list)  
    figure(1);
    x_label = log(error_list(1:end-1)); %log(epsilon_(n-1))
    y_label = log(error_list(2:end));  %log(epsilon_(n))
    plot(x_label,y_label,'-o');
    title('q1 - Newton Raphson solution');
    xlabel('log(epsilon_n-1))'); 
    ylabel('log(epsilon_n)');
    grid on;     
end














