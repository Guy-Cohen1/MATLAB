

s_4 = 1.895494267034;
x0 = pi/2;
tollerence = 10^(-12);  
[X_n_A,error_list_A,Xn_minus_xn_minus_one_list_A] = Fixed_Point_q4_part_A(x0,tollerence,s_4); 
graph_q4_part_A(error_list_A)
n = [1:length(error_list_A)]'; 
X_n = X_n_A(:,1:length(X_n_A)-1)'; 
difference = Xn_minus_xn_minus_one_list_A'; 
Err = error_list_A'; 
T4_A = table(n,X_n, difference, Err); 
disp(T4_A); 
fprintf(' %d iter to converge\n',length(n));
%---------------------------Question 4: Part B-----------------------------
[X_n_4_B,error_list_B,Xn_minus_xn_minus_one_list_B] = Newton_Raphson_q4(x0,tollerence,s_4); 
graph_q4_part_B(error_list_B)
n = [1:length(error_list_B)]'; 
X_n = X_n_4_B(:,1:length(X_n_4_B)-1)'; 
difference = Xn_minus_xn_minus_one_list_B'; 
Err = error_list_B'; 
T4_B = table(n,X_n, difference, Err); 
disp(T4_B); 
fprintf(' %d iter to converge\n',length(n));
%---------------------------Question 4: Part D-----------------------------
s_4_D = 0;
x0 = 1;
[X_n_4_D,error_list_D,Xn_minus_xn_minus_one_list_D] = Fixed_Point_q4_part_D(x0,tollerence,s_4_D); 
graph_q4_part_D(error_list_D)
n = [1:length(error_list_D)]'; 
X_n = X_n_4_D(:,1:length(X_n_4_D)-1)'; 
X_n_Diff = Xn_minus_xn_minus_one_list_D'; 
Error_n = error_list_D'; T4_D = table(n,X_n, X_n_Diff, Error_n); 
disp(T4_D); 
fprintf('It took %d iteration to converge',length(n));




function [X_n_A,error_list_A,Xn_minus_xn_minus_one_list_A] = Fixed_Point_q4_part_A(x0,tollerence,s)%part a
    iter = 2;
    error_list_A = zeros;
    Xn_minus_xn_minus_one_list_A = zeros;
    X_n_A(1) = x0;
    X_n_A(2) = g_A(X_n_A(1));
    while abs(X_n_A(iter) - X_n_A(iter-1)) >= tollerence
        X_n_A(iter+1) = g_A(X_n_A(iter));
        error_list_A(iter) = abs(X_n_A(iter) - s);     
        iter=iter+1;    
        difference = abs(X_n_A(iter)-X_n_A(iter-1));  
        Xn_minus_xn_minus_one_list_A(iter-1) = difference;
    end
end

function [X_n_D,error_list_D,Xn_minus_xn_minus_one_list_D] = Fixed_Point_q4_part_D(x0,tollerence,s) %part 
    iter = 2;
    error_list_D = zeros;
    Xn_minus_xn_minus_one_list_D = zeros;
    X_n_D(1) = x0;
    X_n_D(2) = g_D(X_n_D(1));
    while abs(X_n_D(iter) - X_n_D(iter-1)) >= tollerence
        X_n_D(iter+1) = g_D(X_n_D(iter));
        error_list_D(iter) = abs(X_n_D(iter) - s);     
        iter=iter+1;    
        Xn_Diff = abs(X_n_D(iter)-X_n_D(iter-1));  
        Xn_minus_xn_minus_one_list_D(iter-1) = Xn_Diff;
    end
end

function g_A = g_A(x) %g(x) function
    g_A = 2*sin(x);
end

function g_D = g_D(x) %g(x) function
    g_D = asin(x/2);
end

function [X_n,error_list,Xn_minus_xn_minus_one_list] = Newton_Raphson_q4(x0,tollerence,s) %part b
    X_n = zeros;%initialization 
    error_list = zeros; %initialization 
    Xn_minus_xn_minus_one_list = zeros; %initialization 
    X_n(1) = x0; 
    X_n(2) = x0 - h_q4(x0);
    iter = 2;
    difference = abs(X_n(iter)-X_n(iter-1)); %|xn-xn-1|
    Xn_minus_xn_minus_one_list(1) = difference;
    error_list(1) = abs(X_n(1) - s); %|xn - s| 
    while difference >= tollerence    %10^-12                 
        X_n(iter+1) = X_n(iter) - h_q4(X_n(iter)); %formula
        error_list(iter) = abs(X_n(iter) - s);     
        iter=iter+1;    
        difference = abs(X_n(iter)-X_n(iter-1));  
        Xn_minus_xn_minus_one_list(iter-1) = difference;       
    end
end

function f = func_q4(x) % f(x)
        f = x - 2*sin(x);
end

function f_derive = deriv_func_q4(x) % f'(x) 
        f_derive = 1- 2*cos(x);    
end


function h = h_q4(x)  %  f(x) / f'(x)
        h = func_q4(x) / deriv_func_q4(x);     
end

function graph_q4_part_A(error_list_A)  
    figure(6)
    x_axis = log(error_list_A(1:end-1)); 
    y_axis = log(error_list_A(2:end)); 
    plot(x_axis,y_axis,'-*');
    title('Q4 Part A, Fixed Point');
    xlabel('log(epsilon_n-1))'); 
    ylabel('log(epsilon_n)');
    grid on;   
end



function graph_q4_part_B(error_list_B)  
    figure(7)
    x_axis = log(error_list_B(1:end-1)); 
    y_axis = log(error_list_B(2:end)); 
    plot(x_axis,y_axis,'-*');
    title('Q4 Part B, Newton Raphson');
    xlabel('log(epsilon_n-1))'); 
    ylabel('log(epsilon_n)');
    grid on; 
end



function graph_q4_part_D(error_list_D)  
    figure(8)
    x_axis = log(error_list_D(1:end-1)); 
    y_axis = log(error_list_D(2:end)); 
    plot(x_axis,y_axis,'-o');
    title('Q4 Part D, Fixed Point');
    xlabel('log(epsilon_n-1))'); 
    ylabel('log(epsilon_n)');
    grid on;
end



