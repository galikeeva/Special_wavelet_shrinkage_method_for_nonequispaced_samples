clear
clc
J = 10;
N = 2^J;
eps = 0.01;
for i = 1:N
    z(i) = normrnd(0,1);
    % H^(-1) для определения промежутков времени
    t(i) = (i/N)^2;
    %t(i) = exp(i/N)-1;
end;

H = t.^(0.5);
%H = log(t+1);
for i = 1:N
    if i < 256 || i > 768
        f(i) = -0.5;
    else
        f(i) = 0.5;
    end;
end;
%f = (t.*(1-t)).^0.5.*sin((2.1*pi)./(t+0.05));
%f = exp(t);
%f = -t.^2;
%f = sin(t*(2*pi));
%функция с шумом
y = f+0.1*z;

%приближение масштабирующей и вейвлет функций
[phi,psi,xval] = wavefun('sym8',J);

h = 2*t;
%h = exp(t);

% step 1
%входная матрица Р
for i = 1:N
    for k = 1:N
        % индикатор ближайшей точки t1 на сетке
        [~,idx] = min(abs(t - k/N));
        % индикатор ближайшей точки в разложении phi
        [~,idx1] = min(abs(xval - (N*H(idx) - i)));
        P(k,i) = phi(idx1);
    end;
end;

% step 2
% вейвлет-разложение
[coef, l] = wavedec(P*((N^(-0.5))*(y')), J, 'sym8');



%step 3
%пороговая обработка
theta = coef;
theta1 = coef;
for lev = 2:length(l) - 1
    curr_theta = sum(l(1:lev-1));
    for i = curr_theta:curr_theta + l(lev) - 1
        %рассматриваемый промежуток времени
        tmp = 2^(-(J-lev+2))*(i - curr_theta + 1+N);
        %t_2 = exp(tmp) - 1;
        t_2 = (tmp)^2;
        %супремум производной на отрезке
        %h = exp(t_2);
        h = 2*(t_2);
        lambda = eps * (2* h *(1/N)*log2(N)).^0.5;
        theta(i) = theta(i)*(abs(theta(i)) > lambda);
        ch = (abs(theta(i)) > lambda);
        theta1(i) = sign(theta1(i)) * max(0, (abs(theta1(i)) - lambda));
        %строгая оценка
    end;
end;

%step 4
%обратное вейвлет-преобразование
f1 = waverec(theta, l, 'sym8');
f3 = waverec(theta1, l, 'sym8');

%step 5
% матрица обратного преобразования Р2
for i = 1:N
    for k = 1:N
        [~,idx2] = min(abs(xval - (N*t(k) - i)));
        P2(k,i) = N^0.5*phi(idx2);
    end;
end;

%задействуем Р2
f2 = P2*f1;
f4 = P2*f3;

[c, levels] = wavedec((N^(-0.5))*y', J, 'sym8');
cl = wdenoise(c,10,'DenoisingMethod', 'UniversalThreshold');
cl2 = wdenoise(c,10,'DenoisingMethod', 'Minimax');
f_ch = waverec(N^0.5*cl, levels, 'sym8');
f_ch2 = waverec(N^0.5*cl2, levels, 'sym8');

ind = t<=1;
t = t(ind);
y = y(ind);f = f(ind);f2 = f2(ind);f4 = f4(ind);f_ch = f_ch(ind);f_ch2 = f_ch2(ind);

%графики
plot(t, y, 'LineWidth', 1); hold on; plot(t, f, 'LineWidth', 1); hold on;
plot(t, f2, 'LineWidth', 1);hold on;plot(t, f4, 'LineWidth', 1);hold on;
plot(t, f_ch, 'LineWidth', 1);hold on; plot(t, f_ch2, 'LineWidth', 1);
legend('Зашумленная функция',"Оригинальная функция",'Специальный метод, hard-threshold','Специальный метод, soft-threshold', "UniversalThreshold",'Minimax');