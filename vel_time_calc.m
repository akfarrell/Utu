%function delay2 = vel_time_calc(sta_s)

% ----------- load in velocity structure --------------
vel_s = load('Heather_vel_things/3D_1.75_Utur_0.1km.mat') %structure of velocities
vel_deps = vel_s.zmod;
p_vels = 1000*vel_s.apvel; %p-wave velocity in m/s
elevs = -sta_s.elev; %elevations, negative. Sea level = 0, deeper => positive
vel_deps_corr = (vel_deps-3.5)*1000; %depths corrected for sea level, in m

for s = 1:numel(sta_s.map_dist)
    dist = sta_s.map_dist(s); %maximum map distance from wave to sta, in m
    depth = dist*tand(sta_s.aoi); %determine max depth, in m
    depth_val = elevs(1)+depth %depth value, in m

    dep_index = ceil((depth_val+3500)/100+1) %index of upper vel_deps_corr that matches depth_val
    if dep_index <=0
        total_time(1) = 0;
        if s>1
            dist_in_top_layer = (vel_deps_corr(1)-elevs(s))/cosd(sta_s.aoi);
            total_time(s) = dist_in_top_layer/p_vels(1);
        end
    else
        max_vel = p_vels(dep_index);

        last_full_vel = dep_index-1;
        dist_in_vel = 100/cosd(sta_s.aoi);
        for i = 2:last_full_vel %calculates travel time through most of layers, but not top or last
            time(i) = dist_in_vel/p_vels(i);
        end
        time_full(s) = sum(time);
        
        dist_in_top_layer = (vel_deps_corr(1)-elevs(s))/cosd(sta_s.aoi);
        time_in_top_layer = dist_in_top_layer/p_vels(1);
        dist_in_last_layer = depth_val-vel_deps_corr(dep_index-1);
        time_in_last_layer = dist_in_last_layer/p_vels(dep_index);
        total_time(s) = time_full(s)+time_in_last_layer+time_in_top_layer;
    end
    clear time
end
%total_time
delay2 = sta_s.time_vals_ref-total_time