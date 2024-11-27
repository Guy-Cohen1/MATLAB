%-----------------------part a-----------------------------------%
q = [2;0;7;8;8;1;0;0;4;3;1;8;8;8;8;0;0;5];
list = [1,2,5,10,20,50];
qError_B = zeros(1,6)'; %Relative Error for q for section B
qError_C = zeros(1,6)'; %Relative Error for q for section C
qError_D = zeros(1,6)'; %Relative Error for q for section D
qcon_B = zeros(1,6)'; %Relative Error for q for section B
M = 18;
rho = 1;
h = pi*rho/M;
A = matrix_q1(h,M);
v = A*q;
[L,U,P] = lu(A);
k = cond(A); 
qnorm = norm(q,2);  %norm 2
vnorm = norm(v,2); %norm 2
Anorm_fro = norm(A, 'fro'); % Frobinius norm 
%-----------------------part b-----------------------------------%
for i = 1:6
    M = 18;
    rho = 1;
    h = list(i)*pi*rho/M;
    A = matrix_q1(h,M);
    v = A*q;
    [L,U,P] = lu(A);
    k = cond(A);
    qcon_B(i)=k;
    qnorm = norm(q,2);  %norm 2
    vnorm = norm(v,2); %norm 2
    Anorm_fro = norm(A, 'fro'); % Frobinius norm 
    vector_y = Ly_b(L,P*v); 
    vector_x = Ux_y(U,vector_y); 
    q_aprx = vector_x';
    q_aprx_norm = norm(q_aprx-q,2);
    q_rel_error = q_aprx_norm / qnorm;
    qError_B(i) = q_rel_error;
end
%-----------------------part c-----------------------------------%
for i = 1:6
    M = 18;
    rho = 1;
    h = list(i)*pi*rho/M;
    A_C = matrix_q1(h,M);
    deltav = vnorm .* 10.^(-3); %the new error
    newv = v + deltav;
    [L5,U5,P5] = lu(A_C);
    vector_y_new = Ly_b(L5,P5*newv);
    vector_x_new = Ux_y(U5,vector_y_new);
    q_aprx_c = vector_x_new';%like and common part b
    q_aprx_norm_c = norm(q_aprx_c-q,2);
    q_rel_error_c = q_aprx_norm_c / qnorm;
    qError_C(i) = q_rel_error_c;
end
%-----------------------part d-----------------------------------%
for i = 1:6
    M = 18;
    rho = 1;
    h = list(i)*pi*rho/M;
    A = matrix_q1(h,M);
    deltaA = Anorm_fro .* 10^(-3);
    newA = deltaA + A;
    [L1,U1,P1] = lu(newA);
    vector_y_d = Ly_b(L1,P1*v);
    vector_x_d = Ux_y(U1,vector_y_d);
    q_aprx_d = vector_x_d';%like and common part b,c
    q_aprx_norm_d = norm(q_aprx_d-q,2);
    q_rel_error_d = q_aprx_norm_d / qnorm;
    qError_D(i) = q_rel_error_d; 
end
%-----------------------part e-----------------------------------%
q = [2;0;7;8;8;1;0;0;4;3;1;8;8;8;8;0;0;5];
M = 18;
rho = 1;
for i=1:6
h_e = list(i)*pi*rho/M;
A_e = matrix_q1(h_e,M);
v_e = A_e*q;
[L_e,U_e,P_e] = lu(A_e);
k = cond(A_e); 
qnorm = norm(q,2);  %norm 2
vnorm = norm(v,2); %norm 2
Anorm_fro = norm(A, 'fro'); % Frobinius norm
end
h_graph = list.*((pi * rho) ./ M);
figure(1);
lg = loglog(h_graph,qcon_B,"-",h_graph,qError_B,"-",h_graph,qError_C,"-",h_graph,qError_D,"-");
lg(1).LineWidth = 1.5; %Change width of the line in the graph
lg(2).LineWidth = 1.5;
lg(3).LineWidth = 1.5;
lg(4).LineWidth = 1.5;
title("Relative error and Condition Number of h");
ylabel('Condition Number and Relative error');
xlabel('h(k)');
legend('k(A)[h]','Relative Err','Relative Err v','Relative Err A','Location','northwest');
grid on
%---------------------------functions-------------------------%
function A = matrix_q1(h,M)
rho = 1;
A = zeros(M,M);
for m = 1:M
    for n = 1:M
        R_m_n = sqrt((h+rho*sin(((m*pi)/M))-rho*sin(((n*pi)/M))).^2+(rho*cos((m*pi)/M)-rho*cos((n*pi)/M)).^2);
        A(m,n) = 1 ./ (4*pi*R_m_n) ;
    end
end
end

function y = Ly_b(L,b) 
%lets find vector y
y = zeros(1,length(L));
y(1) = b(1) / L(1,1);
my_sum = 0;
for j=2:length(L)
    for i=1:j-1
        my_sum = my_sum + L(j,i) .* y(i);
    end
    y(j) = (b(j) - my_sum) ./ L(j,j);
    my_sum = 0;
end
end

function x = Ux_y(U,y)  
%lets find vector y
x = zeros(1,length(U));
x(length(U)) = y(length(U)) / U(length(U),length(U));
my_sum = 0;
for i=length(U)-1:-1:1
    for j=length(U):-1:i+1
        my_sum = my_sum + U(i,j).*x(j);
    end
    x(i) = (y(i) - my_sum) ./ U(i,i);
    my_sum = 0;
end
end



