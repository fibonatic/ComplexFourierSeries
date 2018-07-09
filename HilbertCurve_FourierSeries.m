clearvars
close all

%% Hilbert curve

T1 = [0  1;  1  0] / 2;
T2 = [1  0;  0  1] / 2;
T3 = [1  0;  0  1] / 2;
T4 = [0 -1; -1  0] / 2;

H = [0; 0];
for k = 1 : 5
    H = [T1 * H + [-1; -1] / 2 ...
         T2 * H + [-1;  1] / 2 ...
         T3 * H + [ 1;  1] / 2 ...
         T4 * H + [ 1; -1] / 2];
end
% A Hilbert curve is not a closed, so it's mirrored to turn it into a loop.
x = [H(2,:)+1 -flip(H(2,:))-1 H(2,1)+1];
y = [H(1,:) flip(H(1,:)) H(1,1)];

%% Fourier series

% http://www.wolframalpha.com/input/?i=(a%2Bb*T)*exp(-i*n*T)dT
cn = @(a,b,n,T) (b - (a + b * T(2)) * n * 1i) / n^2 * exp(1i * n * T(2)) ...
              - (b - (a + b * T(1)) * n * 1i) / n^2 * exp(1i * n * T(1));
M = length(x);
P = complex(x, y);
T = linspace(0, 2 * pi, M);     % Integral is evaluated between these angles
N = 999;                        % Order of the Fourier series
C = zeros(2, N);
for m = 2 : M
    a = (T(m) * P(m-1) - T(m-1) * P(m)) / T(2); % Linear interpolation of 
    b = (P(m) - P(m-1)) / T(2);                 % Hilbert curve: P(T)=a+b*T
    for n = 1 : N
        C(1,n) = C(1,n) + cn(a, b,  n, T(m-1:m)) / 2 / pi;
        C(2,n) = C(2,n) + cn(a, b, -n, T(m-1:m)) / 2 / pi;
    end
end

K = 10 * M;
Ni = 17;                        % Order Fourier series used in epicycles
t = linspace(0, 2 * pi, K);
Y = zeros(1, K);
for n = 1 : 2 : Ni
    Y = Y + C(1,n) * exp(-1i * n * t) + C(2,n) * exp(1i * n * t);
end

%% Animation

tt = linspace(0, 2 * pi, 100);
cc = complex(cos(tt), sin(tt));
col = [0.23 0.35 0.43];
Col = [1 1 0.3];
k = round(linspace(1, K, K / 30));

% Draw epicycles
cL = [1 0.8 0.8];
cR = [0.8 1 0.8];

fig = figure('Position', [10 50 [1920 1080]*5/6]);
fig.InvertHardcopy = 'off';
for ki = 1 : K/30
    plot(real(Y(1:k(ki))), imag(Y(1:k(ki))), 'Color', Col, 'linewidth', 2)
    hold on
    set(gca, 'Color', 'k')
    pp = 0;
    for n = 1 : 2 : Ni
        p1 = C(2,n) * exp( 1i * n * t(k(ki)));
        plot(real(pp+[0 p1]), imag(pp+[0 p1]), '.-', 'Color', cL)
        plot(real(pp+p1*cc), imag(pp+p1*cc), 'Color', cL, 'markersize', 6)
        pp = pp + p1;
        p2 = C(1,n) * exp(-1i * n * t(k(ki)));
        plot(real(pp+[0 p2]), imag(pp+[0 p2]), '.-', 'Color', cR)
        plot(real(pp+p2*cc), imag(pp+p2*cc), 'Color', cR, 'markersize', 6)
        pp = pp + p2;
    end
    hold off
    axis equal
    set(gca,'XColor','none')
    set(gca,'YColor','none')
    axis([-3 3 -1.6875 1.6875])
    
    name = ['animationTemp2\frame_' num2str(ki,'%0.3d')];
    print(name, '-painters', '-dpng', '-r180')
    A = imread([name '.png']);
    imwrite(A((1:1080)+270, (1:1920)+593, :), [name '.png'])
    % One might have to add a pause inbetween print, imread or imwrite to
    % ensure that the previous command is done with running the file.
end

% Interpolate Fourier order
x0 = 10;
a = (N * x0 - Ni) / (N - Ni);
b = (1 - x0) / (N - Ni);
xn = @(n) round(a + b * n);

fO = 0;
for n = Ni+2 : 2 : N
    yn = C(1,n) * exp(-1i * n * t) + C(2,n) * exp(1i * n * t);
    X = xn(n);
    xx = linspace(0, 1, X);
    for xi = 1 : X
        plot(x, y, 'Color', col * sqrt((n-Ni) / (N-Ni)), 'linewidth', 2)
        hold on
        yy = Y + yn * xx(xi);
        plot(real(yy), imag(yy), 'Color', Col, 'linewidth', 2)
        set(gca, 'Color', 'k')
        hold off
        axis equal
        set(gca,'XColor','none')
        set(gca,'YColor','none')
        axis([-3 3 -1.6875 1.6875])
        text(-0.25, 1.1, ['Order = ' num2str(n - 1 + xx(xi), '%.1f')], ...
            'Color', 'w', 'FontSize', 14)
        
        fO = fO + 1;
        name = ['animationTemp\frame_' num2str(fO,'%0.4d')];
        print(name, '-painters', '-dpng', '-r180')
        A = imread([name '.png']);
        imwrite(A((1:1080)+270, (1:1920)+593, :), [name '.png'])
    end
    Y = Y + yn;
end

%% Write video

outputVideo = VideoWriter('Hilbert_Fourier.mp4', 'MPEG-4');
open(outputVideo)

for ii = 1 : K/30
    name = ['epicycles\frame_' num2str(ii,'%0.3d')];
    A = imread([name '.png']);
    writeVideo(outputVideo, A)
end

name = ['interpolation\frame_' num2str(1,'%0.4d')];
B = imread([name '.png']);
w = linspace(0, 1, 32);
for ii = 2 : 31
    writeVideo(outputVideo, (1 - w(ii)) * A + w(ii) * B)
end

for ii = 1 : fO
    name = ['interpolation\frame_' num2str(ii,'%0.4d')];
    A = imread([name '.png']);
    writeVideo(outputVideo, A)
end

close(outputVideo)
