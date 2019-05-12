function [endpoint,centerpoint,newstate] = endpoint2_140805(Rie,center,startpoint,theta,flight,trans,state)

%cpはその時点でのCenterpoint
%stateは現在地が細胞内か細胞外か
%Rは各円(細胞)の半径. centerはそれらの中心座標
sp = startpoint;
cp = sp;
A = flight;
th1 = theta;
rotate = 0;
hit = zeros(length(Rie),4);

while A > 0.000000005
    cp = 0;
    if state == 0 % 細胞外
        R = Rie(:,2);
        for k = 1:length(R)
            xx = sp(1)-center(k,1); yy = sp(2)-center(k,2); % 面倒なのでいったん各円の中心を原点にする. xx,yyはその場合のスタート地点
            dist = abs(xx*tan(th1)-yy)/sqrt(tan(th1)^2 + 1); % 直線と原点の距離
            if dist < R(k) %細胞外で直線と原点の距離=半径なら接線なので気にしなくてよい.
                al = 1 + tan(th1)^2; be = -tan(th1)*(-yy + tan(th1)*xx); ga = (tan(th1)*xx - yy)^2 - R(k)^2; 
                % 円と直線の交点を求めるための2次式…の各成分: (al)+x^2 + 2*(be)*x + ga = 0 
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
                        hit(k,2) = sqrt(xmov(sel)^2 + ymov(sel)^2); % ぶつかるまでに消費した飛距離
                        hit(k,1) = (xmov(sel)*cos(th1) >= 0 & ymov(sel)*sin(th1) >=0) & (hit(k,2)<=A);
                        % Particleの飛んでいる方向はぶつかる方向かどうか,飛距離の範囲でぶつかるかの判定
                        hit(k,3) = acos(xpoint(sel)/R(k)); % 交点が円から見てX座標から何度にあるか(theta2)
                        if ypoint(sel) < 0; hit(k,3) = -hit(k,3); end % acosはからなず第1,2象限の値を返すので
                            hit(k,4) = th1 - hit(k,3); % 円への入射角(phi).
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
            target = find(hit(:,2) == min(nonzeros(hit(:,2).*hit(:,1)))); %本当にぶつかったのはどの細胞なのか
            A = A - hit(target,2);
            if rand() > trans %透過しない
                th1 = -hit(target,4)+hit(target,3);
                sp2 = [R(target)*cos(hit(target,3)) + center(target,1) , R(target)*sin(hit(target,3)) + center(target,2)];
                cp = cp +  [(sp2(1)+sp(1))/2 , (sp2(2)+sp(2))/2] * (hit(target,2)/flight);
                sp = sp2;
            else %透過する
                state = target;
                sp2 = [R(target)*cos(hit(target,3)) + center(target,1) , R(target)*sin(hit(target,3)) + center(target,2)];
                cp = cp + [(sp2(1)+sp(1))/2 , (sp2(2)+sp(2))/2] * (hit(target,2)/flight);
                sp = sp2;
            end
        end
        
    elseif state > 0 % 細胞内. stateの値はどの細胞の中にいるかを示す
        R = Rie(:,1);
        RR = R(state);
        C = center(state,:);
        sp2 = sp - C; %現在いる円の中心を原点に持ってくる

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
                sel = xmov.^2+ymov.^2 <= A^2 & xmov.^2 + ymov.^2 > 0.0000000001; % どちらも選べないようなときはあり得ない方をまず消す
            end
            if sum(sel) ~= 1
                A = 0;
            else
                th2 = acos(xpoint(sel)/RR); %ぶつかった場所が円の中心(移動しているので原点)から何度の場所にあるか
                if ypoint(sel) <0
                    th2 = -th2;
                end
            
                r = sqrt(xmov(sel)^2 + ymov(sel)^2); %今回消費した距離
                A = A-r;      
                sp2 = rotatevec(sp2, -th2); %円にぶつかった点をx軸(正)の位置まで回転して持ってくるようにスタート地点を回転(-th2回転).
                rotate = rotate - th2;
                phi = sign(sp2(2))*( pi - atan(abs(sp2(2))/abs(RR-sp2(1))) ); % 壁に入ってくる角度(回転座標)
        
                if rand() > trans %透過しない
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
        % 細胞接着を考慮するときのため温存
    end
end

endpoint = sp;
centerpoint = cp;
newstate = state;

        
function [p2] = rotatevec(p1,theta)
    fun = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    p2 = (fun * p1(:))';
       
        
        