
x0 = 5; %initialization 
s = 2;  %initialization 
tollerence = 10^-12; %initialization 
[X_n_q3,error_list,Xn_minus_xn_minus_one_list] = Newton_Raphson_q3(x0,tollerence,s); %lets run NR
graph_q3_part_a(error_list); %print graph
num_of_iter = [1:length(error_list)]'; 
X_n = X_n_q3(:,1:length(X_n_q3)-1)'; 
difference = Xn_minus_xn_minus_one_list'; 
Err = error_list'; 
Table3_A = table(num_of_iter,X_n, difference, Err); 
disp(Table3_A); 
fprintf('It took %d iteration to converge\n',length(num_of_iter));

[X_n_q3_B,error_list_2,Xn_minus_xn_minus_one_list_2] = Newton_Raphson_q3_part_b(x0,tollerence,s); %Run Newton Raphson function
graph_q3_partb(error_list_2); %lets print a graph
n = [1:length(error_list_2)]'; 
X_n = X_n_q3_B(:,1:length(X_n_q3_B)-1)';
difference = Xn_minus_xn_minus_one_list_2'; 
error_list_2 = error_list_2'; 
Table3_B = table(n,X_n, difference, error_list_2);
fprintf(''); 
disp(Table3_B); 
fprintf(' %d iteration to converge\n',length(n));

mul = 3;
[X_n_q3_C,error_list_3,Xn_minus_xn_minus_one_list_3] = Newton_Raphson_q3_part_c(x0,mul,tollerence,s); %Run Newton Raphson function
graph_q3_partc(error_list_3);
n = [1:length(error_list_3)]'; 
X_n = X_n_q3_C(:,1:length(X_n_q3_C)-1)'; 
difference = Xn_minus_xn_minus_one_list_3'; 
Error_n = error_list_3'; 
Table3_C = table(n,X_n, difference, Error_n);
fprintf('%60s\n','<strong> NR2 </strong>'); 
disp(Table3_C); 
fprintf('It took %d iteration to converge\n',length(n));


%---------------------------------functions-------------------------------------------
function [X_n,error_list,Xn_minus_xn_minus_one_list] = Newton_Raphson_q3_part_c(x0,mul,tollerence,s)
    iter = 2; %Number of iterations.
    error_list = zeros; %initialization 
    Xn_minus_xn_minus_one_list = zeros;%initialization 
    X_n(1) = x0;
    X_n(2) = x0 - mul*u(x0);
    difference = abs(X_n(iter)-X_n(iter-1)); %|X_n - X_(n-1)|
    Xn_minus_xn_minus_one_list(1) = difference;
    error_list(1) = abs(X_n(1) - s);
    while (abs(X_n(iter)-X_n(iter-1))>tollerence) && (abs(X_n(iter)-s)<(abs(X_n(iter-1)-s))) % condition for converging
        X_n(iter+1) = X_n(iter) - mul*u(X_n(iter));
        error_list(iter) = abs(X_n(iter) - s);  % |X_n -s|   
        iter = iter+1;    
        difference = abs(X_n(iter)-X_n(iter-1));  
        Xn_minus_xn_minus_one_list(iter-1) = difference; 
    end
end



function [X_n,error_list,Xn_minus_xn_minus_one_list] = Newton_Raphson_q3_part_b(x0,tollerence,s)
    iter = 2; %Number of iterations.
    error_list = zeros; %initialization 
    Xn_minus_xn_minus_one_list = zeros;%initialization 
    X_n(1) = x0;
    X_n(2) = x0 - g(x0);
    difference = abs(X_n(iter)-X_n(iter-1)); %|X_n - X_(n-1)|
    Xn_minus_xn_minus_one_list(1) = difference;
    error_list(1) = abs(X_n(1) - s);
    while (abs(X_n(iter)-X_n(iter-1))>tollerence) && (abs(X_n(iter)-s)<(abs(X_n(iter-1)-s))) % condition for converging
        X_n(iter+1) = X_n(iter) - g(X_n(iter));
        error_list(iter) = abs(X_n(iter) - s);  % |X_n -s|   
        iter = iter+1;    
        difference = abs(X_n(iter)-X_n(iter-1));  
        Xn_minus_xn_minus_one_list(iter-1) = difference; 
    end
end


function [X_n,error_list,Xn_minus_xn_minus_one_list] = Newton_Raphson_q3(x0,tollerence,s) %for part a
    X_n = zeros;%initialization 
    error_list = zeros; %initialization 
    Xn_minus_xn_minus_one_list = zeros; %initialization 
    X_n(1) = x0; 
    X_n(2) = x0 - h_q3(x0);
    iter = 2;
    difference = abs(X_n(iter)-X_n(iter-1)); %|xn-xn-1|
    Xn_minus_xn_minus_one_list(1) = difference;
    error_list(1) = abs(X_n(1) - s); %|xn - s| 
    while difference >= tollerence    %10^-12                 
        X_n(iter+1) = X_n(iter) - h_q3(X_n(iter)); %formula
        error_list(iter) = abs(X_n(iter) - s);     
        iter=iter+1;    
        difference = abs(X_n(iter)-X_n(iter-1));  
        Xn_minus_xn_minus_one_list(iter-1) = difference;       
    end
end



function fx = f(x)  %f(x)
    fx = x^5-6*x^4+14*x^3-20*x^2+24*x-16;
end

function f_tag_x = f_tag_x(x) %f'(x)
    f_tag_x = 5*x^4-24*x^3+42*x^2-40*x+24;
end

function f_tag_tag_x = f_tag_tag_x(x) %f''(x)
    f_tag_tag_x = 20*x^3-72*x^2+84*x-40;
end

function y = u(x) % u = f(x) / f'(x) ; wanted to use the letter u for this
    y = f(x) / f_tag_x(x);
end

function y = u_tag(x) % u'(x) = 1 - (f(x) * f''(x))/(f'(x) * f'(x))
    y =  1-(f(x)*f_tag_tag_x(x))/(f_tag_x(x)^2);
end


function [y] = g(x)  %choose h depents in the section
        y = u(x)/u_tag(x);
    end
    



function f = func_q3(x) % f(x)for part a
        f = x^5-6*x^4+14*x^3-20*x^2+24*x-16 ;
end

function f_derive = deriv_func_q3(x) % for part a
        f_derive = 5*x^4-24*x^3+42*x^2-40*x+24;    
end


function h = h_q3(x)  %  f(x) / f'(x) for part a
        h = func_q3(x) / deriv_func_q3(x);     
end


function graph_q3_part_a(error_list)  
    figure(3);
    x_label = log(error_list(1:end-1)); %log(epsilon_(n-1))
    y_label = log(error_list(2:end));  %log(epsilon_(n))
    plot(x_label,y_label,'-o');
    title('q3 part a - Newton Raphson solution');
    xlabel('log(epsilon_n-1))'); 
    ylabel('log(epsilon_n)');
    grid on;     
end



function graph_q3_partb(error_list_2)
    figure(4)
    x_label = log(error_list_2(1:end-1)); %log(epsilon_(n-1))
    y_label = log(error_list_2(2:end)); %log(epsilon_(n))
    plot(x_label,y_label,'-o');
    title('q3 part b - Newton Raphson multiple solution');
    xlabel('log(\bf\epsilon_n_-_1))');
    ylabel('log(\bf\epsilon_n)');
    grid on;
end

function graph_q3_partc(error_list_3)
    figure(5)
    x_label = log(error_list_3(1:end-1)); %log(epsilon_(n-1))
    y_label = log(error_list_3(2:end)); %log(epsilon_(n))
    plot(x_label,y_label,'-o');
    title('q3 part c - Newton Raphson mul solution, q = 3 ');
    xlabel('log(\bf\epsilon_n_-_1))');
    ylabel('log(\bf\epsilon_n)');
end



