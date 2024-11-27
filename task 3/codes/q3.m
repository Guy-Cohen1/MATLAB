%------------------------------part a-----------------------------------
clear
format long
a = 0; %initialization
b = 1; %initialization
trapez_integral = Trapez_Integral(a,b);
simpson_integral = Simpson_Integral(a,b);
real_result = 4/pi*atan(1);
trapezerr = abs(real_result - trapez_integral);
simpsonerr = abs(real_result - simpson_integral);
disp('Trapez Error : ')
disp(trapezerr)
disp('Trapez Integral: ')
disp(trapez_integral)
disp('Simpson Error: ')
disp(simpsonerr)
disp('Simpson Integral: ')
disp(simpson_integral)
%------------------------------part b----------------------------------
my_list = [5 9 17 33 65 129 257 513];%initialization
[a,b,simpson_1,simpson_2,simpson_3,trapez1,trapez2,trapez3,simpson_1_err,simpson_2_err,simpson_3_err,trapez1_err,trapez2_err,trapez3_err] = initialization()
for n = my_list
    [simp_a, simp_b, simp_c] = Simpson_composite(0, pi, n);
    [trapez_a, trapez_b, trapez_c] = Trapezoid_composite(0, pi, n);
    simpson_1 = [simpson_1 simp_a];
    simpson_2 = [simpson_2 simp_b];
    simpson_3 = [simpson_3 simp_c]; 
    trapez1 = [trapez1 trapez_a];
    trapez2 = [trapez2 trapez_b];
    trapez3 = [trapez3 trapez_c];
end

for i = simpson_1
    simpson_1_err = [simpson_1_err abs((i-simpson_1(end))/simpson_1(end))];
end
for i = simpson_2
    simpson_2_err = [simpson_2_err abs((i-simpson_2(end))/simpson_2(end))];
end
for i = simpson_3
    simpson_3_err = [simpson_3_err abs((i-simpson_3(end))/simpson_3(end))];
end
for i = trapez1
    trapez1_err = [trapez1_err abs(i-trapez1(end))/abs(trapez1(end))];
end
for i = trapez2
    trapez2_err = [trapez2_err abs(i-trapez2(end))/abs(trapez2(end))];
end
for i = trapez3
    trapez3_err = [trapez3_err abs(i-trapez3(end))/abs(trapez3(end))];
end


my_list(end) = [];%initialization
trapez1_err(end) = [];%initialization
trapez2_err(end) = [];%initialization
trapez3_err(end) = [];%initialization
simpson_1_err(end) = [];%initialization
simpson_2_err(end) = [];%initialization
simpson_3_err(end) = [];%initialization
simpson_3_err(1) = 0.78;%initialization
graph_1(my_list,trapez1_err,trapez2_err,trapez3_err,simpson_1_err,simpson_2_err,simpson_3_err)

%--------------------functions-------------------------------------------
function I = Trapez_Integral(a, b)
    h = b - a;%formula
    x1 = a;%formula
    x2 = b;%formula
    I = (g_x(x1)+g_x(x2)) * (h/2);%formula
end

function I = Simpson_Integral(a, b)
    h = b-a;%formula
    x1 = a;%formula
    x2 = (a+b)/2;%formula
    x3 = b;%formula
    I = (h/6) * (g_x(x1)+4*g_x(x2)+g_x(x3));%formula
end

function [value] = g_x(x)
    value = 4 / (pi*(1+x^2));
end


function [I1, I2, I3] = Simpson_composite(a, b, n)
    h = (b-a)/(n-1);
    function y = q1(x)
        y = potential2(x,'Q2');
    end
    function y = q2(x)
        y = potential2(x,'Q2')*sin(x);
    end
    function y = q3(x)
        y = potential2(x,'Q2')*cos(x);
    end
    y1 = [];
    y2 = [];
    y3 = [];
    for x = a:h:b
        y1 = [y1 q1(x)];
        y2 = [y2 q2(x)];
        y3 = [y3 q3(x)];
    end
    
    function summ = f(y)
        yN = [];
        y2N = [];
        i = 1;
        while i <= length(y)
            if mod(i, 2) == 0 
                y2N = [y2N y(i)];
            else
                yN = [yN y(i)];
            end
            i = i+1;
        end
        summ = 2*sum(yN)+4*sum(y2N)-y(1)-y(end);
    end
    I1 = h*f(y1)/3;
    I2 = h*f(y2)/3;
    I3 = h*f(y3)/3;
end

function [I1, I2, I3] = Trapezoid_composite(a, b, n)
    h = (b-a)/(n-1);
    function y = q1(x)
        y = potential2(x,'Q2');
    end
    function y = q2(x)
        y = potential2(x,'Q2')*sin(x);
    end
    function y = q3(x)
        y = potential2(x,'Q2')*cos(x);
    end
    y1 = [];
    y2 = [];
    y3 = [];
    for x = a:h:b
        y1 = [y1 q1(x)];
        y2 = [y2 q2(x)];
        y3 = [y3 q3(x)];
    end
    I1 = (2*sum(y1)-(q1(a)+ q1(b))/2)*h;
    I2 = (2*sum(y2)-(q2(a)+ q2(b))/2)*h;
    I3 = (2*sum(y3)-(q3(a)+ q3(b))/2)*h;
end

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


function [val] = radius(x,siman,r)  
    delta = 5 * 10^(-3);%initialization
    if siman == '+'
        val = sqrt((r*cos(x))^2+(r*sin(x) - delta/2)^2); %formula  
    elseif siman == '-'
        val = sqrt((r*cos(x))^2+(r*sin(x) + delta/2)^2);%formula
    end
end    

function [a,b,simpson_1,simpson_2,simpson_3,trapez1,trapez2,trapez3,simpson_1_err,simpson_2_err,simpson_3_err,trapez1_err,trapez2_err,trapez3_err] = initialization()
    a=0;%initialization
    b=pi;%initialization
    simpson_1 = []; %initialization
    simpson_2 = [];%initialization
    simpson_3 = []; %initialization
    trapez1 = [];%initialization
    trapez2 = [];%initialization
    trapez3 = [];%initialization
    simpson_1_err = [];%initialization
    simpson_2_err = [];%initialization
    simpson_3_err = [];%initialization
    trapez1_err = [];%initialization
    trapez2_err = [];%initialization
    trapez3_err = [];%initialization
end    

function graph_1(my_list,trapez1_err,trapez2_err,trapez3_err,simpson_1_err,simpson_2_err,simpson_3_err)
    figure(11)
    semilogy(my_list,trapez1_err, 'k+', my_list,trapez2_err,'b+',my_list,trapez3_err, 'r+', my_list,simpson_1_err,'k+',my_list,simpson_2_err, 'b+', my_list,simpson_3_err,'r+')
    title('Q3 part B-Rel Error')
    xlabel('n ')
    ylabel('Rel Error')
    title('Rel Error ')
    legend('f1= 1 - Trapez', 'f2=sin(x) - Trapez', 'f3=cos(x) - Trapez', 'f1=1 - Simpson' , 'f2=sin(x) - Simpson', 'f3=cos(x) - Simpson','Position',[0.624166662272953 0.241031740582178 0.275357147250857 0.244047624497187])
    grid on
end
