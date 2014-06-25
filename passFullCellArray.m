function passFullCellArray(number, info)

% Describe the INFO input.
[dim1,dim2] = size(info);  typ = class(info);
fprintf('\nInput 2 is a %d-by-%d %s array.\n', dim1, dim2, typ);

% Show the result.
% fprintf('\nThe %dth element of the array contains:\n', number)
% info{number, :};

end

