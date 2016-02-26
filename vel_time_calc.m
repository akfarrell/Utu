function delay2 = vel_time_calc(sta_s)

% ----------- load in velocity structure --------------
vel_s = load('Heather_vel_things/3D_1.75_Utur_0.1km.mat') %structure of velocities
vel_deps = vel_s.zmod;
p_vels = 1000*vel_s.apvel; %p-wave velocity in m/s
elevs = sta_s.elev; %elevations, positive. Sea level = 0, deeper => negative
vel_deps_corr = (-vel_deps+3.5)*1000; %depths corrected for sea level, in m
dist = sta_s.derp;

for s = 1:numel(sta_s.map_dist)
    %dist = sta_s.map_dist(s); %maximum map distance from wave to sta, in m
    %depth = dist*0.5*sind(sta_s.aoi*2) %determine max depth, in m
    depth = dist(s)*cosd(sta_s.aoi);
    depth_val = elevs(1)-depth; %depth value, in m
    p_vels_extrap = p_vels;
    p_vels_extrap(76:329) = p_vels(75); %use 5712.2447 for values that would be in magma body
    %23 km deep below sea level is index 231

    dep_index = -ceil((depth_val-3500)/100-1); %index of upper vel_deps_corr that matches depth_val
    if dep_index <=2
        
        dist_in_top_layer1 = (elevs(s)-vel_deps_corr(1))/cosd(sta_s.aoi);
        dist_in_top_layer = (elevs(1)-vel_deps_corr(1))/cosd(sta_s.aoi)+(elevs(s)-elevs(1))/cosd(sta_s.aoi);
        if depth <= dist_in_top_layer
            dist_in_top_layer = dist(s);
            s
        end
        dist_in_top_layer - dist_in_top_layer1;
        total_time(s) = dist_in_top_layer/p_vels(1);
        total_time_extrap(s) = dist_in_top_layer/p_vels_extrap(1);
        %s
    else
        %max_vel = p_vels(dep_index);
        last_full_vel = dep_index-1;
        dist_in_vel = 100/cosd(sta_s.aoi);
        for i = 2:last_full_vel %calculates travel time through most of layers, but not top or last
            time(i) = dist_in_vel/p_vels(i);
            time2(i) = dist_in_vel/p_vels_extrap(i);
        end
        time_full(s) = sum(time);
        time_full2(s) = sum(time2);
        
        dist_in_top_layer1 = (elevs(s)-vel_deps_corr(1))/cosd(sta_s.aoi);
        dist_in_top_layer = (elevs(1)-vel_deps_corr(1))/cosd(sta_s.aoi)+(elevs(s)-elevs(1))/cosd(sta_s.aoi);
        if depth <= dist_in_top_layer
            dist_in_top_layer = dist(s);
        end
        dist_in_top_layer - dist_in_top_layer1;
        time_in_top_layer = dist_in_top_layer/p_vels(1);
        time_in_top_layer2 = dist_in_top_layer/p_vels_extrap(1);
        dist_in_last_layer = -(depth_val-vel_deps_corr(dep_index));
        time_in_last_layer = dist_in_last_layer/p_vels(dep_index);
        time_in_last_layer2 = dist_in_last_layer/p_vels_extrap(dep_index);
        total_time(s) = time_full(s)+time_in_last_layer+time_in_top_layer;
        total_time2(s) = time_full2(s)+time_in_last_layer2+time_in_top_layer2;
    end
    dist_in_top_layer
    depth
    clear time
    clear time2
    clear dist_in_top_layer
    clear depth
end
%total_time
delay2 = sta_s.time_vals_ref-total_time;
%delay_extrap = sta_s.time_vals_ref-total_time2;