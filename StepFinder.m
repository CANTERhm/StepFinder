classdef StepFinder < afm
    % STEPFINDER find steps in afm force vs. time curves using the
    % MSF-Algorithm from J. Opfer et al. , 2012 
    %   
    % Copyright (C) 2019  Julian Blaser
    %
    % This file is part of StepFinder.
    %
    % StepFinder is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    % StepFinder is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
    
    properties
        path % file path to valid afm-txt-curve
        x_data % data-vector of the indipendend variable
        y_data % data-vector of the dependend variable
        sigma % sigma of the gaussian kernel used for data smoothing
        w % half window_width
        N % number of data-points within x/y-data
    end
    
    methods
        function obj = StepFinder(varargin)
            %STEPFINDER Construct an instance of this class
            
            % input parser
            p = inputParser;
            
            addRequired(p, 'x_data');
            addRequired(p, 'y_data');
            addParameter(p, 'window_width', 30);
            addParameter(p, 'sigma', 2);
            
            parse(p, varargin{:});
            
            x_data = p.Results.x_data;
            y_data = p.Results.y_data;
            window_width = p.Results.window_width;
            sigma = p.Results.sigma;
            
            % call to superclass constructor
            
            % validation of input parameter
            if isempty(x_data)  || isempty(y_data)
                ME = MException('StepFinder:invalidInput', 'specify valid vectors for x/y_data');
                throw(ME);
            else
                if length(x_data) ~= length(y_data)
                    ME = MException('StepFinder:invalidInput', 'the vectors of "x_data" and "y_data" must have the same length');+
                    throw(ME);
                end
                if ~isvector(x_data) || ~isvector(y_data)
                    ME = MException('StepFinder:invalidInput', 'x/y_data inputs have to be double vectors of the size Nx1 or 1xN');
                    throw(ME);
                else
                    if ~isnumeric(x_data) || isnumeric(y_data)
                        ME = MException('StepFinder:invalidInput', 'x/y_data have to be numerical vectors');
                        throw(ME);
                    end
                end
            end
            
            % make sure thad x/y_data are column vectors
            [x_data, y_data] = prepareCurveData(x_data, y_data);
            
            % construct StepFinder-object
            obj.x_data = x_data;
            obj.y_data = y_data;
            obj.sigma = sigma;
            obj.w = round(window_width/2);
            
        end
        
        function obj = SmoothData(obj)
            % SMOOTHDATA data smoothing by convolution of "y_data" with
            % a gaussian kernel and the standard deviation "sigma"
            
            
        end
    end % normal methods
    
end % classdef

function GaussFunction(sigma)

