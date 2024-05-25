function [P_mag, P_pha] = PFC_DrawTransducer3D(transducer_object_list, linestyle)

if(nargin == 1)
    linestyle = 'none';
else
end

transducer_object = PFC_TransducerArray2Transducer(transducer_object_list);

P = transducer_object.Pressure;

Tri_Points = transducer_object.Tri_Points;
Tri_ConnectivityList = transducer_object.Tri_ConnectivityList;

source_coordinate = this_TriCenter(Tri_Points, Tri_ConnectivityList);
source_nv = this_TriNV(Tri_Points, Tri_ConnectivityList);
source_ds = this_TriArea(Tri_Points, Tri_ConnectivityList);


P_mag = abs(P);
P_pha = angle(P);

P_alpha = P_mag./max(P_mag);


max_pos = max(Tri_Points, [], 1);
min_pos = min(Tri_Points, [], 1);
max_pos_dis = max(max_pos - min_pos);
cen_pos = 0.5.*(max_pos + min_pos);


%% Draw Pressure
fig1 = figure;
trisurf(Tri_ConnectivityList,...
    Tri_Points(:, 1), Tri_Points(:, 2), Tri_Points(:, 3),...
    P_mag, 'LineStyle', linestyle);

colormap('parula');
colorbar;
axis equal;
title('Pressure Magnitude');
xlabel('x');
ylabel('y');
zlabel('z');

xlim(0.5.*[-max_pos_dis , max_pos_dis] + cen_pos(1));
ylim(0.5.*[-max_pos_dis , max_pos_dis] + cen_pos(2));
zlim(0.5.*[-max_pos_dis , max_pos_dis] + cen_pos(3));

cursorMode = datacursormode(fig1);
cursorMode.Enable = 'on';
cursorMode.UpdateFcn = @this_displayMagnitude;


%% Draw Phase
fig2 = figure;
trisurf(Tri_ConnectivityList,...
    Tri_Points(:, 1), Tri_Points(:, 2), Tri_Points(:, 3),...
    P_pha, 'LineStyle', linestyle, 'FaceAlpha', 'flat', 'AlphaDataMapping', 'none', 'FaceVertexAlphaData', P_alpha);

colormap('hsv');
ax = gca;
ax.CLim = [-pi, pi];
c = colorbar;
c.Limits = [-1, 1].*pi;
c.Ticks = (-1:0.5:1).*pi;
c.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};

axis equal;
title('Pressure Phase');
xlabel('x');
ylabel('y');
zlabel('z');

xlim(0.5.*[-max_pos_dis , max_pos_dis] + cen_pos(1));
ylim(0.5.*[-max_pos_dis , max_pos_dis] + cen_pos(2));
zlim(0.5.*[-max_pos_dis , max_pos_dis] + cen_pos(3));

cursorMode = datacursormode(fig2);
cursorMode.Enable = 'on';
cursorMode.UpdateFcn = @this_displayPhase;


%%
fig3 = figure; %#ok<NASGU>
trisurf(Tri_ConnectivityList,...
    Tri_Points(:, 1), Tri_Points(:, 2), Tri_Points(:, 3),...
    source_ds, 'LineStyle', linestyle);

hold on;
quiver3(source_coordinate(:, 1), source_coordinate(:, 2), source_coordinate(:, 3),...
    source_nv(:, 1), source_nv(:, 2), source_nv(:, 3));

axis equal;
title('Normal Vectors');
xlabel('x');
ylabel('y');
zlabel('z');

xlim(0.5.*[-max_pos_dis , max_pos_dis] + cen_pos(1));
ylim(0.5.*[-max_pos_dis , max_pos_dis] + cen_pos(2));
zlim(0.5.*[-max_pos_dis , max_pos_dis] + cen_pos(3));

end


function txt = this_displayMagnitude(~,info)
disp('a');
assignin('base', 'a', info);
x = info.Position(1);
y = info.Position(2);
z = info.Position(3);

[~,col] = ind2sub(size(info.Target.XData),info.Target.Children.DataIndex);

Pmag = info.Target.CData(col);

txt = {['X ', num2str(x)], ['Y ', num2str(y)], ['Z ', num2str(z)], ['P_{mag} ', num2str(Pmag)]};

end


function txt = this_displayPhase(~,info)

x = info.Position(1);
y = info.Position(2);
z = info.Position(3);

[~,col] = ind2sub(size(info.Target.XData),info.Target.Children.DataIndex);

Ppha = info.Target.CData(col);
Ppha_rad = Ppha/pi;
Ppha_deg = Ppha/pi*180;
txt = {['X ', num2str(x)], ['Y ', num2str(y)], ['Z ', num2str(z)],...
    ['P_{pha} ', num2str(Ppha_rad), '\pi (', num2str(Ppha_deg), '\circ)']};

end


function Center_Point = this_TriCenter(Points, ConnectivityList)
P = Points;
C = ConnectivityList;

P1 = P(C(:, 1), :);
P2 = P(C(:, 2), :);
P3 = P(C(:, 3), :);

Center_Point = (P1 + P2 + P3)./3;
end

function nv = this_TriNV(Points, ConnectivityList)
P = Points;
C = ConnectivityList;

P1 = P(C(:, 1), :);
P2 = P(C(:, 2), :);
P3 = P(C(:, 3), :);

v1 = P2 - P1;
v2 =  P3 - P2;

C = cross(v1, v2);
C_length = sqrt(sum(C.^2, 2));
C_length_matrix = C_length*ones(1, 3);
nv = C./C_length_matrix;
end

function A = this_TriArea(Points, ConnectivityList)
P = Points;
C = ConnectivityList;

P1 = P(C(:, 1), :);
P2 = P(C(:, 2), :);
P3 = P(C(:, 3), :);

a = sqrt(sum((P1 - P2).^2, 2));
b = sqrt(sum((P2 - P3).^2, 2));
c = sqrt(sum((P3 - P1).^2, 2));

s = (a + b + c)/2;

A = sqrt(s.*(s - a).*(s - b).*(s - c));
end

