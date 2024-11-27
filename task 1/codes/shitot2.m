%-------------------------------part a and b----------------------------%
q = [2;0;7;8;8;1;0;0;4;3;1;8;8;8;8;0;0;5];
rho = 1;
M = 18;
list = [5,2,1];
Rel_dist_between_iter_A = zeros;%initialization
Rel_dist_between_iter_B = zeros;%initialization
Rel_dist_between_iter_C = zeros;%initialization
Rel_err_with_q_A = zeros;%initialization
Rel_err_with_q_B = zeros;%initialization
Rel_err_with_q_C = zeros;%initialization
num_of_iter_A=zeros;%initialization
num_of_iter_B=zeros;%initialization
num_of_iter_C=zeros;%initialization

for i = 1:3
    h = pi ./ (list(i).*M);
    A = matrix_q2(h,M);
    v = A * q;
    if i == 1
        [Rel_dist_between_iter_A,Rel_err_with_q_A,num_of_iter_A,qkA] = Gauss_Seidel(A,v,q); 
    end
    if i == 2
        [Rel_dist_between_iter_B,Rel_err_with_q_B,num_of_iter_B,qkB] = Gauss_Seidel(A,v,q); 
    end
    if i == 3 
        [Rel_dist_between_iter_C,Rel_err_with_q_C,num_of_iter_C,qkC] = Gauss_Seidel(A,v,q); 
    end          
end

figure(2); %plot for h = pi/5*M
plotA = semilogy(num_of_iter_A,Rel_dist_between_iter_A,num_of_iter_A,Rel_err_with_q_A,'-*');
comment2 = sprintf('Took %d iterations',length(num_of_iter_A));
comment_2_2 = sprintf('Final Relative error of q^k and q is %d',Rel_err_with_q_A(end));
plotA(1).LineWidth = 2;
plotA(2).LineWidth = 2;
text(2,10^-4,comment2,'Color','blue','FontSize',12);
text(2,10^-4.5,comment_2_2,'Color','blue','FontSize',12);
title('Gauss-Seidel- h=pi/5M');
xlabel('Iter');
legend('Rel Destination between two iters','Rel error between iter and q');

figure(3); %plot for h = pi/2*M
plotB = semilogy(num_of_iter_B,Rel_dist_between_iter_B,num_of_iter_B,Rel_err_with_q_B,'-*');
comment_3 = sprintf('Took %d iterations to converge',length(num_of_iter_B));
comment_3_2 = sprintf('Final Relative error of q^k and q is %d',Rel_err_with_q_B(end));
plotB(1).LineWidth = 2;
plotB(2).LineWidth = 2;
text(2,10^-4,comment_3,'Color','blue','FontSize',14);
text(2,10^-4.5,comment_3_2,'Color','black','FontSize',14);
title('Gauss-Seidel- h=pi/2M');
legend('Rel Destination between two iters','Rel error between iter and q');

figure(4); %plot for h = pi/M
plotC = semilogy(num_of_iter_C,Rel_dist_between_iter_C,num_of_iter_C,Rel_err_with_q_C,'-*');
comment_4 = sprintf('Solution Doesnt Converge :(');
text(2,10^-4,comment_4,'Color','blue','FontSize',14);
plotC(1).LineWidth = 2;
plotC(2).LineWidth = 2;
title('Gauss-Seidel- h=pi/M');
legend('Rel Destination between two iters','Rel error between iter and q');



%------------------------------part C----------------------------------%
h = pi ./ (5.*M);
A = matrix_q2(h,M);
v = A * q;
Rel_dist_between_iter = zeros;%initialization
Rel_err_with_q = zeros;%initialization
num_of_iter=zeros;%initialization
[Rel_dist_between_iter,Rel_err_with_q,num_of_iter,qkD,norm_G] = Jacobi(A,v,q); %h=pi/5M
figure(5);
plot_D = semilogy(num_of_iter,Rel_dist_between_iter,num_of_iter,Rel_err_with_q,'-*');
plot_D(1).LineWidth = 2;
plot_D(2).LineWidth = 2;
comment_5 = sprintf('Solution Doesnt Converge :(');
text(10,10^6,comment_5,'Color','green','FontSize',14);
title('Jacobi- h=pi/5M');
legend('Rel Destination between two iters','Rel error between iter and q');

%---------------------------------------part D--------------------------%
h = pi ./ (5.*M);
A_1 = matrix_q2_part_D(h,M);
v_1= A_1 * q;
Rel_dist_between_iter = zeros;%initialization
Rel_err_with_q = zeros;%initialization
num_of_iter=zeros;%initialization
[Rel_dist_between_iter,Rel_err_with_q,num_of_iter,qkE] = Jacobi(A_1,v_1,q); %h=pi/5M
figure(6);
plot_E = semilogy(num_of_iter,Rel_dist_between_iter,num_of_iter,Rel_err_with_q,'-*');
plot_E(1).LineWidth = 2;
plot_E(2).LineWidth = 2;
comment = sprintf('Took %d iterations',length(num_of_iter));
comment_7 = sprintf('Final Relative error of q^k and q is %d',Rel_err_with_q(end));
comment_5 = sprintf('Solution Converge :)');
text(1.5,10^-4.1,comment,'Color','green','FontSize',14);
text(1.5,10^-4.4,comment_7,'Color','blue','FontSize',12);
text(1.5,10^-4.7,comment_5,'Color','blue','FontSize',12);
title('Jacobi- h=pi/5M , part D');
legend('Rel Destination between two iters','Rel error between iter and q');

%---------------------------functions-------------------------%
function A = matrix_q2(h,M)
rho = 1;
A = zeros(M,M);
for m = 1:M
    for n = 1:M
        R_m_n = sqrt((h+rho*sin(((m*pi)/M))-rho*sin(((n*pi)/M))).^2+(rho*cos((m*pi)/M)-rho*cos((n*pi)/M)).^2);
        A(m,n) = 1 ./ (4*pi*R_m_n) ;
    end
end
end

function A = matrix_q2_part_D(h,M)
rho = 1;
A = zeros(M,M);
for m = 1:M
    for n = 1:M
        r_mn = (h+rho*sin(((m*pi)/M))-rho*sin(((n*pi)/M))).^2+(rho*cos((m*pi)/M)-rho*cos((n*pi)/M)).^2;
        A(m,n) = 1 ./ (4*pi*r_mn) ;
    end
end
end

function [Rel_dist_between_iter,Rel_err_with_q,num_of_iter,qk_out] = Gauss_Seidel(A,v,q)
    Max_Error = 10 ^ (-3);
    D = diag(diag(A)); %find D
    L = tril(A,-1); %find L
    Q = L + D; 
    U_minus = Q - A; %find -U
    Q_minus_one = inv(Q);%Q^(-1)
    G = Q_minus_one * U_minus  ; % by definition
    C =  Q_minus_one * v;   % by definition
    qk = C;  %iter number one
    Rel_dist_between_iter = zeros; %initialization
    Rel_err_with_q = zeros;%initialization
    iter = 1;   
    num_of_iter = zeros;%initialization
    error = max(abs(q-qk)); %infinity norm
    counter = 500; %Limit Iterations
    while abs(error) > Max_Error && iter <= counter
       num_of_iter(iter) = iter;
       qk_minus1 = qk;   %initialization
       qk = G*(qk_minus1) + C; 
       error = norm(q-qk,'inf');
        Rel_err_with_q(iter) = norm(qk - q,'inf') / norm(q,'inf');
       Rel_dist_between_iter(iter) = norm(qk - qk_minus1, 'inf') / norm(qk_minus1, 'inf');
       iter = iter + 1;
    end
    qk_out = qk;
end

function [Rel_dist_between_iter,Rel_err_with_q,num_of_iter,qk_out,norm_G] = Jacobi(A,v,q)   
    Max_Error = 10 ^ (-3); %maximum erorr
    D = diag(diag(A)); %lets find D
    Q = D ;
    I = eye(18);
    Q_minus_one = inv(Q); %Q^(-1)
    G = I - (Q_minus_one * A);
    norm_G = norm(G,'inf'); %lets check if converge
    C = Q_minus_one * v; % by definition
    qk = C;  %iter number one
    Rel_dist_between_iter = zeros; %initialization
    Rel_err_with_q = zeros; %initialization
    iter = 1;   
    num_of_iter = zeros;  %initialization
    error = max(abs(q-qk));  %infinity norm
    counter = 500; %Limit Iterations
    while abs(error) > Max_Error && iter <=counter  
       qk_minus1 = qk;
       qk = G*(qk_minus1) + C; %initialization
       error = norm(q-qk,'inf');
       Rel_err_with_q(iter) = norm(qk - q,'inf') / norm(q,'inf');
       Rel_dist_between_iter(iter) = norm(qk - qk_minus1, 'inf') / norm(qk_minus1, 'inf');
       num_of_iter(iter) = iter;
       iter = iter + 1;
    end
    qk_out = qk;
end