function [magnitude, integphase] = func_integphase190511(phase)

% 【入力】
%  phase -> ベクトル. n個のphase(radian)を並べたベクトル. (全てが単位ベクトルを想定)

% 【出力】
%  magnitude -> 合成したベクトルの大きさ
%　integphase -> 合成したベクトルのphase.

numelmol = length(phase);

X = sum(cos(phase));
Y = sum(sin(phase));

magnitude = norm([X,Y])/numelmol;
integphase = acos(X/(numelmol*magnitude)) * sign(Y);
