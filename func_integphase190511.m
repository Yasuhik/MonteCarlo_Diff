function [magnitude, integphase] = func_integphase190511(phase)

% �y���́z
%  phase -> �x�N�g��. n��phase(radian)����ׂ��x�N�g��. (�S�Ă��P�ʃx�N�g����z��)

% �y�o�́z
%  magnitude -> ���������x�N�g���̑傫��
%�@integphase -> ���������x�N�g����phase.

numelmol = length(phase);

X = sum(cos(phase));
Y = sum(sin(phase));

magnitude = norm([X,Y])/numelmol;
integphase = acos(X/(numelmol*magnitude)) * sign(Y);
