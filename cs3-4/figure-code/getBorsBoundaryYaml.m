function [IIIa] = getBorsBoundaryYaml(file_name)
    %-- Function to pull the relevant data for the given wind rose .yaml file (cs1&2)

    BoundaryStruct = ReadYaml(file_name);    % Pull our .yaml file into a struct
    % Shortcuts for the .yaml paths
    bound = BoundaryStruct(1).boundaries(1);
    % Pull the data we need
    IIIa = cell2mat(bound.IIIa(:,:));
end