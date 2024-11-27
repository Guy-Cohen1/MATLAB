%-----------------------part a and b :)-----------------------------------%
A_determinata = zeros(5,1);%initialization
Rel_error_q = zeros(5,1);%initialization
q = [2;0;7;8;8;1;0;0;4;3;1;8;8;8;8;0;0;5];%initialization
M = 18;%initialization
rho = 1 ;%initialization
h_list = [1/5,1/2,2,5,10];%initialization
h_graph = h_list.*((pi * rho) ./ M);%initialization
for i=1:5
     h = h_graph(i);
     A = matrix_q3(h,18);
     v = A * q; 
     A_determinata(i) = abs(det(A)); %list of final determinanats
     A_transpose = transpose(A);
     q_aprx = inv(A_transpose * A) * A_transpose * v; %formula
     Rel_error_q(i) = (norm(q_aprx - q)) / norm(q); %formula, list of final relative errors
end

%lets print a graph :)
figure(7);
lg_3 = loglog(h_graph,A_determinata,"*-",h_graph,Rel_error_q,"*-");
lg_3(1).LineWidth = 2;
lg_3(2).LineWidth = 2;
title('Least Squares');
xlabel('h');
legend('|det(A)|','|q"-q|/|q|','Location','southwest');
grid on;

%---------------------------------------functions----------------%
function A = matrix_q3(h,M) %lets build the matrix for q3
rho = 1 ;%initialization
A = zeros(M,M);%initialization
for m = 1:M
    for n = 1:M
        R_m_n = sqrt((h+rho*sin((m*pi)/M)-rho*sin((n*pi)/M)).^2+(rho*cos((m*pi)/M)-rho*cos((n*pi)/M)).^2);
        A(m,n) = 1 ./ (4*pi*R_m_n) ;
    end
end
end