function [endpoint,centerpoint,newstate] = endpoint2_140805(Rie,center,startpoint,theta,flight,trans,state)

%cp�͂��̎��_�ł�Centerpoint
%state�͌��ݒn���זE�����זE�O��
%R�͊e�~(�זE)�̔��a. center�͂����̒��S���W
sp = startpoint;
cp = sp;
A = flight;
th1 = theta;
rotate = 0;
hit = zeros(length(Rie),4);

while A > 0.000000005
    cp = 0;
    if state == 0 % �זE�O
        R = Rie(:,2);
        for k = 1:length(R)
            xx = sp(1)-center(k,1); yy = sp(2)-center(k,2); % �ʓ|�Ȃ̂ł�������e�~�̒��S�����_�ɂ���. xx,yy�͂��̏ꍇ�̃X�^�[�g�n�_
            dist = abs(xx*tan(th1)-yy)/sqrt(tan(th1)^2 + 1); % �����ƌ��_�̋���
            if dist < R(k) %�זE�O�Œ����ƌ��_�̋���=���a�Ȃ�ڐ��Ȃ̂ŋC�ɂ��Ȃ��Ă悢.
                al = 1 + tan(th1)^2; be = -tan(th1)*(-yy + tan(th1)*xx); ga = (tan(th1)*xx - yy)^2 - R(k)^2; 
                % �~�ƒ����̌�_�����߂邽�߂�2�����c�̊e����: (al)+x^2 + 2*(be)*x + ga = 0 
                t = sqrt(be^2 - (al*ga));
                xpoint = [ (-be+t)/al , (-be-t)/al ]; ypoint = tan(th1).* xpoint + yy - tan(th1)*xx;
                xmov = xpoint - xx;
                ymov = ypoint - yy;
                if isempty(min(abs(xmov(abs(xmov)>0.000001)))) == 0
                    sel = abs(xmov) == min(abs(xmov(abs(xmov)>0.000001)));
                    if sum(sel) ~= 1
                        sel = abs(ymov) == min(abs(ymov(abs(ymov)>0.000001)));
                    end
                    if sum(sel) ~= 1
                        A = 0;
                    else
                        hit(k,2) = sqrt(xmov(sel)^2 + ymov(sel)^2); % �Ԃ���܂łɏ�����򋗗�
                        hit(k,1) = (xmov(sel)*cos(th1) >= 0 & ymov(sel)*sin(th1) >=0) & (hit(k,2)<=A);
                        % Particle�̔��ł�������͂Ԃ���������ǂ���,�򋗗��͈̔͂łԂ��邩�̔���
                        hit(k,3) = acos(xpoint(sel)/R(k)); % ��_���~���猩��X���W���牽�x�ɂ��邩(theta2)
                        if ypoint(sel) < 0; hit(k,3) = -hit(k,3); end % acos�͂���Ȃ���1,2�ی��̒l��Ԃ��̂�
                            hit(k,4) = th1 - hit(k,3); % �~�ւ̓��ˊp(phi).
                        if cos(hit(k,4))<0; hit(k,4) = hit(k,4)-pi; end
                    end
                else
                    A =0;
                end
            end
        end
        
        if sum(hit(:,1)) == 0
            ep = [sp(1)+A*cos(th1), sp(2)+A*sin(th1)];
            cp = cp + [ (ep(1) + sp(1))/2 , (ep(2) + sp(2))/2 ] * (A/flight);
            sp = ep;
            A = 0;
        else
            target = find(hit(:,2) == min(nonzeros(hit(:,2).*hit(:,1)))); %�{���ɂԂ������̂͂ǂ̍זE�Ȃ̂�
            A = A - hit(target,2);
            if rand() > trans %���߂��Ȃ�
                th1 = -hit(target,4)+hit(target,3);
                sp2 = [R(target)*cos(hit(target,3)) + center(target,1) , R(target)*sin(hit(target,3)) + center(target,2)];
                cp = cp +  [(sp2(1)+sp(1))/2 , (sp2(2)+sp(2))/2] * (hit(target,2)/flight);
                sp = sp2;
            else %���߂���
                state = target;
                sp2 = [R(target)*cos(hit(target,3)) + center(target,1) , R(target)*sin(hit(target,3)) + center(target,2)];
                cp = cp + [(sp2(1)+sp(1))/2 , (sp2(2)+sp(2))/2] * (hit(target,2)/flight);
                sp = sp2;
            end
        end
        
    elseif state > 0 % �זE��. state�̒l�͂ǂ̍זE�̒��ɂ��邩������
        R = Rie(:,1);
        RR = R(state);
        C = center(state,:);
        sp2 = sp - C; %���݂���~�̒��S�����_�Ɏ����Ă���

        xx = sp2(1); yy = sp2(2);
        dist = sqrt((xx + A*(cos(th1)))^2 + (yy + A*(sin(th1)))^2);
        if dist <= RR
            ep = [ xx + A*cos(th1) , yy+A*sin(th1)];
            sp2 = rotatevec(ep,-rotate)+C;
            cp = cp + (A/flight).*(sp2+sp)/2;
            sp = sp2;
            A = 0;
            rotate = 0;
            
        else
            al = 1 + tan(th1)^2; be = -tan(th1)*(-yy + tan(th1)*xx); ga = (tan(th1)*xx - yy)^2 - RR^2;
            t = sqrt(be^2 - (al*ga));
        
            xpoint = [ (-be+t)/al , (-be-t)/al ];
            ypoint = tan(th1)* xpoint + yy - tan(th1)*xx;
            xmov = xpoint - xx;
            ymov = ypoint - yy;
            
            sel = xmov * cos(th1) >=0 & ymov * sin(th1) >= 0;
            if sum(sel) ~= 1
                sel = xmov.^2+ymov.^2 <= A^2 & xmov.^2 + ymov.^2 > 0.0000000001; % �ǂ�����I�ׂȂ��悤�ȂƂ��͂��蓾�Ȃ������܂�����
            end
            if sum(sel) ~= 1
                A = 0;
            else
                th2 = acos(xpoint(sel)/RR); %�Ԃ������ꏊ���~�̒��S(�ړ����Ă���̂Ō��_)���牽�x�̏ꏊ�ɂ��邩
                if ypoint(sel) <0
                    th2 = -th2;
                end
            
                r = sqrt(xmov(sel)^2 + ymov(sel)^2); %������������
                A = A-r;      
                sp2 = rotatevec(sp2, -th2); %�~�ɂԂ������_��x��(��)�̈ʒu�܂ŉ�]���Ď����Ă���悤�ɃX�^�[�g�n�_����](-th2��]).
                rotate = rotate - th2;
                phi = sign(sp2(2))*( pi - atan(abs(sp2(2))/abs(RR-sp2(1))) ); % �ǂɓ����Ă���p�x(��]���W)
        
                if rand() > trans %���߂��Ȃ�
                    th1 = -phi - rotate;
                    sp2 = rotatevec([RR,0],-rotate)+C;
                    cp = cp + (r/flight)*(sp+sp2)/2;
                    sp = sp2;
                    rotate = 0;
                else
                
                sp2 = rotatevec([RR,0],-rotate)+C;
                cp = cp + (r/flight)*(sp+sp2)/2;
                sp = sp2;
                th1 = phi + pi - rotate;

                
                rotate = 0;
                state = 0;
                end
            end
        end
    elseif state < 0
        % �זE�ڒ����l������Ƃ��̂��߉���
    end
end

endpoint = sp;
centerpoint = cp;
newstate = state;

        
function [p2] = rotatevec(p1,theta)
    fun = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    p2 = (fun * p1(:))';
       
        
        