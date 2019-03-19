
global c_f;
global mass;
global Wind;
global XWdata;
global YWdata;
global Fdata;
global Rdata;


c_f = 0.47;
Wind = csvread("W.csv");
XWdata = 1:1600;
YWdata = 1:1600;
Fdata = 1:200000;
mass = input("Input mass (kg): ");
Rdata = csvread("F.csv");

Fgenerate();
Wgenerate();

% начальная скорость снаряда
v0 = input("Input start speed (m/s): ");
vector = [10, 0, 10];
vector = vector / norm(vector, 1);
vLength = sqrt(sum(vector.*vector)); 
x0 = 0;
y0 = input("Input start height (m): ");
z0 = 0;
vx0 = v0 * vector(1) / vLength;
vy0 = v0 * vector(2) / vLength;
vz0 = v0 * vector(3) / vLength;
Y0 = [x0; y0; z0; vx0; vy0; vz0];
finish(1) = input("Input finish x coordinate (m): ");
finish(2) = 0;
finish(3) = input("Input finish z coordinate (m): ");
% начальный момент времени
t0 = 0;
% конечный момент времени
tend = 37.0;
% шаг выдачи решения
deltaT = 0.004;

t = t0:deltaT:tend;

% Интрегрируем систему уравнений движения
[tim ,Y] = ode45(@f, t, Y0);
result = zeros(0, 6);
ploter = zeros(0, 3);
for i = 1:length(Y)
  if Y(i, 2) > 0
    result = [result; Y(i, :)]; 
    if length(Y) < 500 || mod(i  - 1, fix(length(Y) / 500)) == 0
      ploter = [ploter; Y(i,1:3)];
    end
  end
end
vec = result(end, 1:3);
vec(2) = 0;
vec = finish - vec;
result(:,1) = result(:,1) + vec(1);
result(:,3) = result(:,3) + vec(3);
ploter(:,1) = ploter(:,1) + vec(1);
ploter(:,3) = ploter(:,3) +  vec(3);
csvwrite("result.csv", result);
plot3(ploter(:,1), ploter(:,3), ploter(:,2));
vector
startpoint = result(1, 1:3)
grid on;
view([-45 45 45]);

xlabel("x");
ylabel("z");
zlabel("y");
clear
function Fgenerate()
  global c_f;
  global mass;
  global Fdata;
  global Rdata;
  svec = 1:(length(Rdata)  + 1);
  svec(1) = 0;
  for i = 1:length(Rdata)
    v = fix(Rdata(i, 1) * 100 + 1);
    Fdata(v) = Rdata(i, 2);
    if Rdata(i, 1) ~= 0 && Rdata(i, 2) ~= 0
      svec(i) = 2 * Rdata(i, 2) / (1.29 * c_f * Rdata(i, 1) * Rdata(i, 1));
    end
  end
  svec(end) = svec(end - 1);
  Rdata(end + 1, 1) = 2000;
  for u = 1:(length(Rdata) - 1)
    first = Rdata(u, 1) * 100 + 1;
    second = Rdata(u + 1, 1) * 100;
    fvalue = svec(u);
    svalue = svec(u + 1);
    for  i = first:second
      if Fdata(i) == i
        sectional = fvalue + (((svalue - fvalue) / ( -first + second)) * (i - first)); 
        v = (i - 1) / 100;
        Fdata(i) = 0.5  * c_f *  1.29 * v * v * sectional;
      end
    end
  end
end

function Wgenerate()
  global Wind;
  global XWdata;
  global YWdata;
  Wind(length(Wind) + 1, :) = Wind(length(Wind),:);
  Wind(end, 1) = 1600; 
  for i = 1:length(Wind)
    v = Wind(i, 1) + 1;
    XWdata(v) = Wind(i, 2);
    YWdata(v) = Wind(i, 3);
  end
  for u = 1:(length(Wind) - 1)
    first = Wind(u, 1) + 1;
    second = Wind(u + 1, 1);
    xfvalue = Wind(u, 2);
    yfvalue = Wind(u, 3);
    for  i = first:second
      if XWdata(i) == i
        XWdata(i) = xfvalue;
      end
      if YWdata(i) == i
        YWdata(i) = yfvalue;
      end
    end
  end
end

function [dYdt] = f(t, Y)
  global mass;
  global Fdata;
  global XWdata;
  global YWdata;
  g = 9.81; % Ускорение свободного падения
  vLength = sqrt(sum((Y(4:6)).*(Y(4:6))));
  R = Fdata(fix(vLength * 100));
  first = 0;
  second = 0;
  dYdt = zeros(6, 1);
  if Y(2) > 0 && Y(2) < 1600
    first = XWdata(fix(Y(2))  + 1);;
    second = YWdata(fix(Y(2))  + 1);
  end
  dYdt(1) = Y(4) + first; 
  dYdt(2) = Y(5); 
  dYdt(3) = Y(6) + second;  
  dYdt(4) = - (R / mass) * Y(4) / vLength;    
  dYdt(5) = - (R / mass) * Y(5) / vLength - g;   
  dYdt(6) = - (R / mass) * Y(6) / vLength;
  
end
  