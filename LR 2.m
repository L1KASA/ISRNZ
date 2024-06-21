% 4.	Рассмотреть с помощью интеграла Рэлея-Зоммерфельда первого рода дифракцию плоской монохроматической волны на прямоугольной апертуре с размерами 5λ×3λ.
%Получить распределение интенсивности вдоль оси x на расстояниях от экрана равных 10λ и 50λ (рисунок 2).
%Сравнить результаты полученные при дифракции на прямоугольной апертуре с аналитическим решением в дальней зоне (sinc).

% Очистка рабочего пространства и закрытие графических окон
close all; % Закрытие всех открытых графических окон
clear;     % Очистка рабочего пространства от переменных

% Параметры сцены
lambda = 1;            % Длина волны
a = 5 * lambda;        % Ширина сцены
b = 3 * lambda;        % Высота сцены
k = 2 * pi / lambda;   % Волновое число
grid_size = 10;        % Размер сетки
N = 1000;              % Количество точек в сетке
x = linspace(-a/2, a/2, N); % Координаты по оси x
y = linspace(-b/2, b/2, N); % Координаты по оси y
[X, Y] = meshgrid(x, y);    % Создание сетки координат
z = [10 50] * lambda;       % Расстояния до экрана

% Начальное поле
initial_wave = @(z) exp(1i*2*pi*(z/lambda)); % Функция начальной волны

% Функция интенсивности (аналитическое)
intensity_function = @(x, y, z) (sin(k*a.*X/(2*z))./(k*a.*X/(2*z))).^2 .* (sin(k*b.*Y/(2*z))./(k*b.*Y/(2*z))).^2;

% Расчет интеграла
for z_index = 1:length(z)
    wave_field = zeros(2*grid_size, 2*grid_size); % Создание массива для волнового поля
    for i = 1:(2*grid_size+1)
        for j = 1:(2*grid_size+1)
            X_val = i-(grid_size+1); % Координата x текущей ячейки
            Y_val = j-(grid_size+1); % Координата y текущей ячейки
            % Функция подынтегрального выражения
            integrand = @(x,y) initial_wave(z(z_index)).*(exp(1i.*k.*(sqrt((X_val-x).^2+(Y_val-y).^2+z(z_index).^2)))./((X_val-x).^2+(Y_val-y).^2+z(z_index).^2)).*(1+(1i./(k.*sqrt((X_val-x).^2+(Y_val-y).^2+z(z_index).^2))));
            % Вычисление интеграла по двумерной области
            wave_field(i,j) = integral2(integrand,-2.5*lambda,2.5*lambda, -1.5*lambda,1.5*lambda);
        end
    end
    % Дальнейшая обработка волнового поля (коррекции фазы волнового поля и масштабирование амплитуды)
    wave_field = -(1i*k*z(z_index))/(2*pi).*wave_field;
    % Вычисление интенсивности
    intensity = intensity_function(X, Y, z(z_index));
    intensity = intensity / max(intensity(:)); % Нормировка интенсивности

    % Визуализация результатов
    % Визуализация волнового поля
    figure;
    imagesc(-grid_size:grid_size, -grid_size:grid_size, wave_field.*conj(wave_field));
    xlabel('x');
    ylabel('y');
    title(['z = ' num2str(z(z_index)) 'λ']); % Заголовок с расстоянием до экрана
    colorbar;

    % Визуализация интенсивности
    figure;
    imagesc(-grid_size:grid_size, -grid_size:grid_size, intensity);
    xlabel('x');
    ylabel('y');
    title(['z = ' num2str(z(z_index)) 'λ']); % Заголовок с расстоянием до экрана
    colorbar;
end

