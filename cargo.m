load data.mat
global c_f = 0.47;
global mass = 100;
global S = 0;
global Wind = importdata("out.csv");
global XWdata = 1:1600;
global YWdata = 1:1600;
global Fdata = 1:30000;
global Rdata;

function Fgenerate()
  global c_f;
  global mass;
  global S;
  global Fdata;
  global Rdata;
  svec(1) = 0;
  for i = 1:rows(Rdata)
    v = fix(Rdata(i, 1) * 100 + 1);
    Fdata(v) = Rdata(i, 2);
    if Rdata(i, 1) != 0 && Rdata(i, 2) != 0
      svec(i) = 2 * Rdata(i, 2) / (1.29 * c_f * Rdata(i, 1) * Rdata(i, 1));
    endif
  endfor
  if svec(1) == 0
    svec(1) = svec(2);
  endif
  S = sum(svec) / rows(Rdata);
  svec(rows(Rdata) + 1) = svec(rows(Rdata));
  Rdata(rows(Rdata) + 1, 1) = 300;
  for u = 1:(rows(Rdata) - 1)
    first = Rdata(u, 1) * 100 + 1;
    second = Rdata(u + 1, 1) * 100;
    fvalue = svec(u);
    svalue = svec(u + 1);
    for  i = first:second
      if Fdata(i) == i
        sectional = fvalue + (((svalue - fvalue) / ( -first + second)) * (i - first)); 
        v = (i - 1) / 100;
        Fdata(i) = 0.5  * c_f *  1.29 * v * v * sectional;
      endif
    endfor
  endfor
end

function Wgenerate()
  global Wind;
  global XWdata = 1:1600;
  global YWdata = 1:1600;
  Wind(rows(Wind) + 1, :) = Wind(rows(Wind),:);
  Wind(end, 1) = 1600; 
  for i = 1:rows(Wind)
    v = Wind(i, 1) + 1;
    XWdata(v) = Wind(i, 2);
    YWdata(v) = Wind(i, 3);
  endfor
  for u = 1:(rows(Wind) - 1)
    first = Wind(u, 1) + 1;
    second = Wind(u + 1, 1);
    xfvalue = Wind(u, 2);
    yfvalue = Wind(u, 3);
    for  i = first:second
      if XWdata(i) == i
        XWdata(i) = xfvalue;
      endif
      if YWdata(i) == i
        YWdata(i) = yfvalue;
      endif
    endfor
  endfor
end

function dYdt = f(Y, t)
  global c_f;
  global mass;
  global S;
  global Fdata;
  global XWdata;
  global YWdata;
  g = 9.81; % ”скорение свободного падени€
  vLength = sqrt(sum((Y(4:6)).*(Y(4:6))));
  R = Fdata(fix(vLength * 100));
  first = 0;
  second = 0;
  if Y(2) > 0 && Y(2) < 1600
    first = XWdata(fix(Y(2))  + 1);
    second = YWdata(fix(Y(2))  + 1);
  endif
  dYdt(1) = Y(4) + first; % dx/dt = vx
  dYdt(2) = Y(5); % dy/dt = vy
  dYdt(3) = Y(6) + second;  % dz/dt = vz
  dYdt(4) = - (R / mass) * Y(4) / vLength;    % dvx/dt = 0
  dYdt(5) = - (R / mass) * Y(5) / vLength - g;   % dvy/dt = -g
  dYdt(6) = - (R / mass) * Y(6) / vLength;   % dvz/dt = 0
endfunction


Fgenerate();
Wgenerate();
%-------------------------------------------------------------------------------
% ѕараметры стрельбы
%-------------------------------------------------------------------------------

% начальна€ скорость снар€да
v0 = 400;
% угол наклона ствола пушки к горизонту
vector = [10, 0, 10];
vector /= norm(vector, 1); 
%-------------------------------------------------------------------------------
% Ќачальные услови€
%-------------------------------------------------------------------------------
x0 = 0;
y0 = 1000;
z0 = 0;
vx0 = v0 * vector(1);
vy0 = v0 * vector(2);
vz0 = v0 * vector(3);
Y0 = [x0; y0; z0; vx0; vy0; vz0];

%-------------------------------------------------------------------------------
% ѕараметры временного интервала
%-------------------------------------------------------------------------------

% начальный момент времени
t0 = 0;
% конечный момент времени
tend = 37.0;
% шаг выдачи решени€
deltaT = 0.1;

% ћассив интересующих нас моментов времени
t = [t0:deltaT:tend];

% »нтрегрируем систему уравнений движени€
Y = lsode("f", Y0, t);
axis([0 3000 0 3000 0 3000]);

% –исуем траекторию полета снар€да
plot3(Y(:,1), Y(:,3), Y(:,2));
xlabel("x");
ylabel("z");
zlabel("y");