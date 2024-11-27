a = 1;
b = 5;
I1 = 318888005; %initialization 
I2 = 207881004; %initialization 
tollerence = 10^(-12); %initialization 
s = 3^(1/4); 
x0 = a + (b-a)*(I1/(I1+I2)); %initialization 
x1 = x0 + (b-x0)*(I1/(I1+I2)); %formula
[X_n_q2,error_list,Xn_minus_xn_minus_one_list] = Secant_method(x0,x1,tollerence,s); %lets do Secant method 
graph_q2(error_list); %lets print a graph
n = [1:length(error_list)]'; 
X_n = X_n_q2(:,1:length(X_n_q2)-1)'; 
difference = Xn_minus_xn_minus_one_list'; 
Err = error_list'; 
T2 = table(n,X_n, difference, Err); 
disp(T2); 
fprintf(' %d iter to converge\n',length(n));












%-----------------------------------functions---------------------------------------
function [X_n,error_list,Xn_minus_xn_minus_one_list] = Secant_method(x0, x1, tollerence, s)
    num_of_iter = 2; %initialization 
    error_list = zeros; %initialization 
    Xn_minus_xn_minus_one_list = zeros;%initialization 
    X_n(1) = x0; %given
    X_n(2) = x1; %given
    difference = abs(X_n(num_of_iter)-X_n(num_of_iter-1)); %|Xn - Xn-1|
    Xn_minus_xn_minus_one_list(1) = difference;
    error_list(1) = abs(X_n(1) - s); %|xn - s| 
    while difference > tollerence
        X_n(num_of_iter+1) = Secant_formula(num_of_iter, X_n);
        error_list(num_of_iter) = abs(X_n(num_of_iter) - s);  % |Xn -s|
        num_of_iter = num_of_iter + 1;    
        difference = abs(X_n(num_of_iter)-X_n(num_of_iter-1));  %|Xn - Xn-1|
        Xn_minus_xn_minus_one_list(num_of_iter-1) = difference;
    end
end



function Secant = Secant_formula(num_of_iter, X_n) %define function
        Secant = X_n(num_of_iter) - func_q2(X_n(num_of_iter))*(X_n(num_of_iter)-X_n(num_of_iter-1))/(func_q2(X_n(num_of_iter))-func_q2(X_n(num_of_iter-1)));
end



function fx = func_q2(x) %given
    fx = x^4-3; 
end


function graph_q2(error_list)  
    figure(2)
    x_label = log(error_list(1:end-1)); 
    y_label = log(error_list(2:end)); 
    plot(x_label,y_label,'-o');
    title('q2: Secant-M');
    xlabel('log(epsilon_n-1))'); 
    ylabel('log(epsilon_n)');
    grid on; 
end




