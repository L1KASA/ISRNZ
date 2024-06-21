%16.	Получить продольное распределение интенсивности светового поля (вдоль оси z при r=0) в фокусе апланатического объектива при фокусировке радиально-поляризованного света
%l = 1 моды R-TEM01 (параметр β считать равным 1) с длиной волны 0,532 мкм в воздухе (n = 1).
%Вычисления произвести для значения числовой апертуры NA равного 0,9.
%Повторить вычисления для параметра β равного 0,5, 1,5, 2,0. Оценить влияние параметра в уменьшение размеров фокусного пятна.

% Очистка рабочего пространства и закрытие графических окон
close all; % Закрытие всех открытых графических окон
clear;     % Очистка рабочего пространства от переменных

%% Параметры
% Длина волны
lambda = 0.532; % мкм
n = 1;
% Заданные значения числовой апертуры и параметров бета
NA = 0.9;
Betas = [0.5, 1.5, 2.0];

% Константы
A = 1;
A_ = 1;
N = 100; % количество дроблений интервала

R = 1;

%%
% Волновое число
k = 2*pi/lambda; % Вычисляем волновое число k
x = linspace(-R, R, N); % Создаем вектор x с равномерно распределенными значениями от -R до R
y = linspace(-R, R, N); % Создаем вектор y с равномерно распределенными значениями от -R до R
[X, Y] = meshgrid(x, y); % Создаем сетку X и Y с использованием x и y
[~, r] = cart2pol(X, Y); % Преобразуем координаты из декартовых в полярные и получаем радиус r
z = linspace(-R, R, N); % Создаем вектор z с равномерно распределенными значениями от -R до R

T_func = @(tetha) cos(tetha).^0.5;  % Функция апотизации зрачка для апланатического объектива
l_func = @(A, tetha, alpha, beta)  A*sin(tetha).*... % Начальное распределение напряженности электрического поля
    exp(-(beta^2*sin(tetha).^2)/(sin(alpha)^2));

% Функция для E_r
fun = @(tetha, alpha, beta) l_func(A_, tetha, alpha, beta) * ...
    T_func(tetha) * sin(2*tetha) * ...
    exp(1i*k*z*cos(tetha)) * besselj(1, k*r*sin(tetha));
% Функция для E_z
fun2 = @(tetha, alpha, beta) l_func(A_, tetha, alpha, beta) *...
    T_func(tetha) * (sin(tetha))^2 * ...
    exp(1i*k*z*cos(tetha)) * besselj(0, k*r*sin(tetha));

I = zeros(N, numel(Betas)); % Массив для сохранения интенсивности светового поля для каждого значения beta
E_r = zeros(N, numel(Betas)); % Массив для сохранения компоненты E_r
E_z = zeros(N, numel(Betas)); % Массив для сохранения компоненты E_z

for i = 1:numel(Betas)
    beta = Betas(i);

    alpha = asin(NA / n); % Вычисляем максимальный угол

    % Интегрируем обе функции
    integral_r = integral(@(tetha)fun(tetha, alpha, beta),  0, alpha, 'ArrayValued', true); % Интегрируем функцию fun по углу tetha от 0 до alpha
    integral_z = integral(@(tetha)fun2(tetha, alpha, beta), 0, alpha, 'ArrayValued', true); % Интегрируем функцию fun2 по углу tetha от 0 до alpha
    E_r(:, i) = A * integral_r; % Вычисляем E_r
    E_z(:, i) = A * 1i * 2 * integral_z; % Вычисляем E_z
    % Находим интенсивности
    I_r = abs(E_r(:, i)).^2; % Вычисляем интенсивность I_r
    I_z = abs(E_z(:, i)).^2; % Вычисляем интенсивность I_z
    I(:, i) = I_r + I_z; % Записываем суммарную интенсивность в массив I
end

% Визуализация интенсивности для разных значений beta
figure;
for i = 1:numel(Betas)
    subplot(2, numel(Betas), i);
    plot(z, real(E_r(:, i)), 'LineWidth', 2); % Строим график вещественной части E_r по z для каждого значения beta
    xlabel('z, мкм');
    ylabel('Re(E_r), a.u.');
    title(sprintf('Re(E_r) для \\beta = %.1f', Betas(i))); % Заголовок графика
    grid on;

    subplot(2, numel(Betas), numel(Betas) + i);
    plot(z, real(E_z(:, i)), 'LineWidth', 2); % Строим график вещественной части E_z по z для каждого значения beta
    xlabel('z, мкм');
    ylabel('Re(E_z), a.u.');
    title(sprintf('Re(E_z) для \\beta = %.1f', Betas(i))); % Заголовок графика
    grid on;
end

% Визуализация интенсивности для разных значений beta
figure;
for i = 1:numel(Betas)
    subplot(1, numel(Betas), i);
    plot(z, I(:, i), 'LineWidth', 2); % Строим график интенсивности по z для каждого значения beta
    xlabel('z, мкм');
    ylabel('I, a.u.');
    title(sprintf('Интенсивность для \\beta = %.1f', Betas(i))); % Заголовок графика
    grid on;

    % Определение FWHM и отметка на графике
    I_focused = I(:, i); % Интенсивность светового поля для текущего значения параметра beta
    [max_intensity, max_pos] = max(I_focused); % Находим максимальное значение интенсивности и его позицию
    half_max_intensity = max_intensity / 2; % Находим половину максимальной интенсивности
    left_half = find(I_focused(1:max_pos) < half_max_intensity, 1, 'last'); % Левый край
    right_half = find(I_focused(max_pos:end) < half_max_intensity, 1) + max_pos - 1; % Правый край
    FWHM = z(right_half) - z(left_half); % Ширина фокусного пятна на полувысоте
    % Отмечаем значение FWHM на графике
    hold on;
    plot([z(left_half), z(right_half)], [half_max_intensity, half_max_intensity], 'r--', 'LineWidth', 1.5);
    text(z(right_half), max_intensity/2, sprintf('FWHM=%.2f', FWHM), 'HorizontalAlignment', 'left');
    hold off;
end

% Оценка влияния параметра beta на уменьшение размеров фокусного пятна
FWHM_values = zeros(size(Betas)); % Создаем массив для сохранения значений FWHM

for i = 1:numel(Betas)
    I_focused = I(:, i); % Интенсивность светового поля для текущего значения параметра beta
    % Находим максимальное значение интенсивности и его позицию
    [max_intensity, max_pos] = max(I_focused);

    % Находим полуширину фокусного пятна на полувысоте
    half_max_intensity = max_intensity / 2;
    left_half = find(I_focused(1:max_pos) < half_max_intensity, 1, 'last'); % Левый край
    right_half = find(I_focused(max_pos:end) < half_max_intensity, 1) + max_pos - 1; % Правый край
    FWHM_values(i) = z(right_half) - z(left_half); % Ширина фокусного пятна на полувысоте
end

% График зависимости FWHM от параметра beta
figure;
plot(Betas, FWHM_values, 'o-', 'LineWidth', 2);
xlabel('Параметр \beta');
ylabel('Ширина FWHM, мкм');
title('Зависимость ширины фокусного пятна на полувысоте от параметра \beta');
grid on;

