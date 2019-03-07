function [turb_coords, fname_turb, fname_wr] = getTurbLocYAML(file_name)
    %-- Function to pull the turbine locations from the example .yaml file
    
    fname_turb = 'no such file';
    fname_wr = 'no such file';
    
    TurbLocStruct = ReadYaml(file_name);    % Pull our .yaml file into a struct
    % Shortcut to reused .yaml path
    defs = TurbLocStruct(1).definitions(1);
    coords = defs.position(1).items(1);
    turb_coords.x = cell2mat(coords.xc(:));
    turb_coords.y = cell2mat(coords.yc(:));
   % Read Turbine locations into an array
   
   ref_list_turbs = defs.wind_plant(1).properties(1).layout(1);
   ref_list_wr = defs.plant_energy(1).properties(1).wind_resource_selection(1).properties(1);
   
   % Get size of lists
   num_turb_files = numel(ref_list_turbs.items);
   num_wr_files = numel(ref_list_wr.items);
   
   % Iterate until first reference not lead by a '#'
   for i = 1:num_turb_files
        ref = ref_list_turbs.items(i);  % Rip the reference list
        ref = ref{1}(1).x0x24ref;       % Take string out of cell and structs
        if ~contains(ref, '#')          % Search for first file w/o '#'
            fname_turb = ref;           % Note filename for other .yaml files
        end
   end
   for i = 1:num_wr_files
      ref = ref_list_wr.items(i);       % Rip the reference list
        ref = ref{1}(1).x0x24ref;       % Take string out of cell and structs
        if ~contains(ref, '#')          % Search for first file w/o '#'
            fname_wr = ref;             % Note filename for other .yaml files
        end
   end
end