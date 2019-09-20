classdef StepFinder
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
        y_conv % smoothed y_data
        theta = []% vector of RSS values
        window_width % width of the moving window
        step_refinement
        numpnts % number of data-points within x/y-data
        smoothing_sigma % standard deviation of gaussion smoothing 
        smoothing_window_width % window withd of gaussian smoothing
        step_results
    end % properties
       
    methods
        function obj = StepFinder(varargin)
            %STEPFINDER Construct an instance of this class
            
            % input parser
            p = inputParser;
            
            addRequired(p, 'x_data');
            addRequired(p, 'y_data');
            addParameter(p, 'window_width', 30);
            addParameter(p, 'step_refinement', 1);
            addParameter(p, 'smoothing_sigma', 5);
            addParameter(p, 'smoothing_window_width',50);
            
            parse(p, varargin{:});
            
            x_data = p.Results.x_data;
            y_data = p.Results.y_data;
            window_width = p.Results.window_width;
            step_refinement = p.Results.step_refinement;
            smoothing_sigma = p.Results.smoothing_sigma;
            smoothing_window_width = p.Results.smoothing_window_width;
            
            % validation of input parameter
            if isempty(x_data)  || isempty(y_data)
                ME = MException('StepFinder:invalidInput', 'specify valid vectors for x/y_data');
                throw(ME);
            else
                if length(x_data) ~= length(y_data)
                    ME = MException('StepFinder:invalidInput', 'the vectors of "x_data" and "y_data" must have the same length');
                    throw(ME);
                end
                if ~isvector(x_data) || ~isvector(y_data)
                    ME = MException('StepFinder:invalidInput', 'x/y_data inputs have to be double vectors of the size Nx1 or 1xN');
                    throw(ME);
                else
                    if ~isnumeric(x_data) || ~isnumeric(y_data)
                        ME = MException('StepFinder:invalidInput', 'x/y_data have to be numerical vectors');
                        throw(ME);
                    end
                end
            end
            if mod(window_width, 2)
                ME = MException('StepFinder:invalidInput', 'input for "window_width" must be even');
                throw(ME);
            end
            
            % make sure thad x/y_data are column vectors
            [x_data, y_data] = prepareCurveData(x_data, y_data);
            
            % construct StepFinder-object
            obj.x_data = x_data;
            obj.y_data = y_data;
            obj.smoothing_sigma = smoothing_sigma;
            obj.smoothing_window_width = smoothing_window_width;
            obj.window_width = window_width;
            obj.step_refinement = step_refinement;
            obj.numpnts = length(x_data);
            
        end
        
        function obj = SmoothData(obj, varargin)
            % SMOOTHDATA data smoothing by convolution of "y_data" with
            % a gaussian kernel 
            
            % input parser
            p = inputParser;
            
            addRequired(p, 'obj');
            addOptional(p, 'sigma', []);
            addOptional(p, 'width', []);
            
            parse(p, obj, varargin{:});
            
            obj = p.Results.obj;
            sigma = p.Results.sigma;
            width = p.Results.width;
            
            % check input parameters
            if ~isempty(sigma)
                obj.smoothing_sigma = sigma;
            end
            if ~isempty(width)
                obj.smooting_window_width = width;
            end
            
            % smooth data
            kernel = StepFinder.gaussian_kernel(obj.window_width, obj.smoothing_sigma);
            obj.y_conv = conv(obj.y_data, kernel, 'same');
            
        end % SmoothData
        
        function obj = StepSearch(obj, varargin)
            % STEPSEARCH find steps wihtin the given data using the
            % MSF-Algorithm
            
            % input parser
            p = inputParser;
            
            addRequired(p, 'obj');
            addOptional(p, 'window_width', []);
            
            parse(p, obj, varargin{:});
            
            obj = p.Results.obj;
            width = p.Results.window_width;
            
            % check input parameters
            if ~isempty(width)
                w = round(width/2);
                obj.window_width = width;
                win_width = width;
            else
                w = round(obj.window_width/2);
                win_width = obj.window_width;
            end
            N = obj.numpnts;
            t = zeros(length(obj.x_data),1);
            
            % calculate theta
            f = NaN(win_width-1, N);
            g = NaN(win_width-1, N);
            ydat = NaN(win_width-1, N);
            RSS_g = NaN(1, N);
            RSS_f = NaN(1, N);
            indices = NaN(2, N);
            
            m = NaN(1, N);
            m_0 = NaN(1, N);
            t_l = NaN(1, N);
            t_r = NaN(1, N);
            t_0 = NaN(1, N);
            
            for i = 1+w:win_width/obj.step_refinement:N
                
                % calculate needed parameters
                xvec = obj.x_data;
                yvec = obj.y_conv;

                
                a = i-w;
                b = i+w-1;
                if a < 1
                    w = i-1;
                    a = i-w;
                    b = i+w_1;
                end
                if b > N
                    w = N-i+1;
                    a = i-w;
                    b = N-i+1;
                end

                xfit = xvec(a:b);
                yfit = yvec(a:b);
                
                aux_ones = ones(round((b-a)/2),1);
                aux_zeros = zeros(round((b-a)/2),1);
                c_left = vertcat(aux_ones, aux_zeros);
                c_right = vertcat(aux_zeros, aux_ones);
                
                coeffs_piecewise = [xfit c_left c_right]\yfit;
                m_i = coeffs_piecewise(1);
                t_i_l = coeffs_piecewise(2);
                t_i_r = coeffs_piecewise(3);
                
                coeffs_global = [xfit vertcat(aux_ones, aux_ones)]\yfit;
                m_i_0 = coeffs_global(1);
                t_i_0 = coeffs_global(2);

                m(1,i) = m_i;
                m_0(1,i) = m_i_0;
                t_l(1,i) = t_i_l;
                t_r(1,i) = t_i_r;
                t_0(1,i) = t_i_0;
                
                a = i-w;
                b = i+w-1;
                c = win_width-1;
                if a < 1
                    a = 1;
                end
                if b > N
                    b = N;
                end
                
                f_i = zeros(c,1);
                for j = 1:c
                    if j < w
                        try
                            f_i(j) = m_i*xvec(j+a)+t_i_l;
                        catch
                        end
                    else
                        try
                            f_i(j) = m_i*xvec(j+a)+t_i_r;
                        catch 
                        end
                    end
                end
                
                g_i = zeros(c,1);
                for j = 1:c
                    try
                        g_i(j) = m_i_0*xvec(j+a)+t_i_0;
                    catch
                    end
                end
                
                y_i = zeros(c,1);
                for j = 1:c
                    try
                        y_i(j) = yvec(j+a);
                    catch
                    end
                end
                
                % calculation of theta
                RSS_gi = StepFinder.msfSum((g_i-y_i).^2);
                RSS_fi = StepFinder.msfSum((f_i-y_i).^2);
                t(i) = (RSS_gi-RSS_fi)*(t_i_r-t_i_l);
                
                RSS_g(1,i) = RSS_gi;
                RSS_f(1,i) = RSS_fi;
                g(:,i) = g_i;
                f(:,i) = f_i;
                ydat(:,i) = y_i;
                indices(1, i) = a;
                indices(2, i) = b;
                
                disp(num2str(i))
                
            end
            
            mask = isnan(m);
            obj.step_results.m = m(~mask);
            
            mask = isnan(m_0);
            obj.step_results.m_0 = m_0(~mask);
            
            mask = isnan(t_l);
            obj.step_results.t_l = t_l(~mask);
            
            mask = isnan(t_r);
            obj.step_results.t_r = t_r(~mask);
            
            mask = isnan(t_0);
            obj.step_results.t_0 = t_0(~mask);
            
            mask = isnan(g);
            g(:,mask(1,:)) = [];
            obj.step_results.g = g;
            
            mask = isnan(f);
            f(:,mask(1,:)) = [];
            obj.step_results.f = f;
            
            mask = isnan(ydat);
            ydat(:,mask(1,:)) = [];
            obj.step_results.y = ydat;
            
            mask = isnan(indices);
            indices(:,mask(1,:)) = [];
            obj.step_results.indices = indices;
            
            mask = isnan(RSS_g);
            obj.step_results.RSS_g = RSS_g(~mask);
            
            mask = isnan(RSS_f);
            obj.step_results.RSS_f = RSS_f(~mask);
            
            % calculate theta
            obj.theta = t./max(t);
            
        end % StepSearch
        
    end % normal methods
    
    methods 
        
        function obj = set.smoothing_sigma(obj, sigma)
            obj.smoothing_sigma = sigma;
        end % set.smoothing_sigma
        
        function obj = set.smoothing_window_width(obj, window_width)
            obj.smoothing_window_width = window_width;
        end % set.smoothing_window_withd
        
%         function obj = set.window_width(obj, width)
%             
%             if mod(width, 2)
%                 ME = MException('StepFinder:ivalidInput', 'input for "window_width" must be even');
%                 throw(ME);
%             end
%             obj.window_width = width;
%             
%         end
        
        function val = get.y_conv(obj)
            if isempty(obj.y_conv)
                val = obj.y_data;
            else
                val = obj.y_conv;
            end
        end % get.y_conv
        
    end % setter/getter methods
    
    methods % getter for dependent properties
        
    end
    
    methods(Static)

        function kernel = gaussian_kernel(width, sigma)
            
            vec = linspace(-(width-1)/2, (width-1)/2, width);
            gvec = exp(-vec.^2./(2*sigma)^2);
            
            kernel = gvec./sum(gvec);
        end % gaussian_kernel
        
        function s = msfSum(inVec1, varargin)

            % input parser
            p = inputParser;

            validVec = @(x) isvector(x) && isnumeric(x);

            addRequired(p, 'inVec1', validVec);
            addOptional(p, 'inVec2', []);
            addParameter(p, 'start_index', 1);
            addParameter(p, 'stop_index', []);

            parse(p, inVec1, varargin{:});

            inVec1 = p.Results.inVec1;
            inVec2 = p.Results.inVec2;
            start_index = p.Results.start_index;
            stop_index = p.Results.stop_index;

            % check input parameter
            if ~isempty(inVec2)
                if length(inVec1) ~= length(inVec2)
                    ME = MException('smfSum:invalidInput', 'input vectors are not the same size');
                    throw(ME);
                end
            end
            if isempty(stop_index)
                stop_index = length(inVec1);
            end

            % calculate sum
            vec1 = inVec1(start_index:stop_index);
            if ~isempty(inVec2)
                vec2 = inVec2(start_index:stop_index);
            else
                vec2 = ones(length(vec1), 1);
            end

            warning('off');
            [vec1, vec2] = prepareCurveData(vec1, vec2);
            warning('on');
            s = vec1'*vec2;

        end % sum
        
    end % static methods
    
end % classdef

%% helper functions
% function s = msfSum(inVec1, varargin)
% 
%     % input parser
%     p = inputParser;
%     
%     validVec = @(x) isvector(x) && isnumeric(x);
%     
%     addRequired(p, 'inVec1', validVec);
%     addOptional(p, 'inVec2', []);
%     addParameter(p, 'start_index', 1);
%     addParameter(p, 'stop_index', []);
%     
%     parse(p, inVec1, varargin{:});
%     
%     inVec1 = p.Results.inVec1;
%     inVec2 = p.Results.inVec2;
%     start_index = p.Results.start_index;
%     stop_index = p.Results.stop_index;
%     
%     % check input parameter
%     if ~isempty(inVec2)
%         if length(inVec1) ~= length(inVec2)
%             ME = MException('smfSum:invalidInput', 'input vectors are not the same size');
%             throw(ME);
%         end
%     end
%     if isempty(stop_index)
%         stop_index = length(inVec1);
%     end
%     
%     % calculate sum
%     vec1 = inVec1(start_index:stop_index);
%     if ~isempty(inVec2)
%         vec2 = inVec2(start_index:stop_index);
%     else
%         vec2 = ones(length(vec1), 1);
%     end
%     
%     warning('off');
%     [vec1, vec2] = prepareCurveData(vec1, vec2);
%     warning('on');
%     s = vec1'*vec2;
% 
% end % sum
