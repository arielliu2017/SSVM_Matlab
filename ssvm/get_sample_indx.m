function sample_indx = get_sample_indx(all_sample_ids, coord_ids, indx, varargin)

assert(isnumeric(indx) || strcmp(indx, ':'), 'INDICES OR : EXPECTED');
sample_size = size(all_sample_ids, 1);

EPSILON = 1e-10;

if(strcmp(indx, ':'))
    sample_indx = (1:sample_size)';
else
    % Default search from cell string
    if(nargin == 3)
        coord_id = coord_ids{indx};
        cmp_cell_array =  cellfun(@(x) strcmp(x, coord_id),  all_sample_ids);
        sample_indx = find(cmp_cell_array == 1);
    else
        % Otherwise search from integer array
        assert(nargin == 4, 'UNEXPECTED NUMBER OF ARGUMENTS');
        
        assert(strcmp('INT', varargin), 'EXPECTED INT');
        coord_id = coord_ids(indx);
        % May be floating point like time
        sample_indx = find(abs(all_sample_ids -coord_id) < EPSILON);
    end
end

end