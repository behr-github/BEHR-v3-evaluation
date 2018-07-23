classdef RunningAverage < handle
    %RUNNINGAVERAGE A class build to simplify carrying out running averages
    %   Instantiate an instance of RunningAverage and add new data
    %   (weighted or not) with the addData method. Calculate the weighted
    %   average of the data given so far with getWeightedAverage.
    
    properties(SetAccess = protected)
        values = [];
        weights = [];
    end
    
    methods
        function obj = RunningAverage(varargin)
            %RunningAverage Instantiate a new running average instance
            %   RunningAverage() initializes its values and weights to
            %   empty arrays, which will be populated the first time
            %   addData() is called.
            %
            %   RunningAverage(SHAPE) will initialize values and weights as
            %   NaN arrays with size SHAPE. All data added via addData must
            %   have that size.
            
            if nargin > 0
                obj.values = nan(varargin{1});
                obj.weights = nan(varargin{1});
            end
        end
        
        function addData(obj, values, weights)
            %addData Add additional data to the running average
            %   addData(VALUES) adds the data VALUES to the running
            %   average, with weights of 1. VALUES must be a numeric array
            %   the same size as the existing values. Any NaNs in the array
            %   are given a weight of 0. VALUES is added to the existing
            %   values with NANADD(), which treats NaNs as 0.
            %
            %   addData(VALUES, WEIGHTS) uses the WEIGHTS given to weight
            %   the data in the average, instead of assuming weights of 1.
            %   Any NaN in VALUES still has its weight set to 0.
            
            if nargin < 3
                weights = ones(size(values));
            end
            
            E = JLLErrors;
            if ~isnumeric(values) || ~isnumeric(weights)
                E.badinput('VALUES and WEIGHTS must both be numeric')
            elseif ~isequal(size(values), size(weights))
                E.badinput('VALUES and WEIGHTS must be equal in size')
            elseif any(isnan(weights(:)) & ~isnan(values(:)))
                E.badinput('There is a NaN weight for a non-NaN value');
            end
            
            % By default, give zero weight to NaN values as they should not
            % really contribute to the average.
            weights(isnan(values)) = 0;
            
            if isempty(obj.values)
                obj.values = values .* weights;
                obj.weights = weights;
            else
                if ~isequal(size(values), size(obj.values))
                    E.badinput('Size of input VALUES is different than existing values (%s)', mat2str(size(obj.values)));
                elseif ~isequal(size(weights), size(obj.weights))
                    E.badinput('Size of input WEIGHTS is different than existing weights (%s)', mat2str(size(obj.weights)));
                end
                obj.values = nanadd(obj.values, values .* weights);
                obj.weights = nanadd(obj.weights, weights);
            end
        end
        
        function avg = getWeightedAverage(obj)
            % getWeightedAverage Returns the current weighted average of
            % the data.
            avg = obj.values ./ obj.weights;
        end
    end
end

